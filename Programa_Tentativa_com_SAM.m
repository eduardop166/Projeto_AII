% PROJETO 2 - CC segmentation (SAM com prompts automáticos + rápido)

% carregar modelo SAM
model = segmentAnythingModel;

% carregar volume (.hdr + .img)
[hdr, path] = uigetfile('*.hdr');
V = double(analyze75read(fullfile(path, hdr)));

% escolher 5 slices onde aparece o CC
z = round(linspace(62, 72, 5));

% guardar máscaras
masks = cell(1, numel(z));

% parâmetros de pré-processamento e speed-up
sigma = 1;
scale = 0.5;

figure;
for i = 1:numel(z)

    % extrair slice e corrigir orientação
    S = squeeze(V(:,:,z(i)));
    S = rot90(S,2);

    % normalizar e suavizar
    p = prctile(S(:), [2 98]);
    Sn = (S - p(1)) / (p(2)-p(1) + eps);
    Sn = min(max(Sn,0),1);
    Sn = imgaussfilt(Sn, sigma);

    % definir ROI central
    [H,W] = size(Sn);
    x1 = round(0.35*W); x2 = round(0.75*W);
    y1 = round(0.28*H); y2 = round(0.62*H);
    R  = Sn(y1:y2, x1:x2);
    [hR,wR] = size(R);

    % prompts automáticos na ROI (1 positivo = máximo, negativos = bordas)
    [~, idxMax] = max(R(:));
    [py, px] = ind2sub(size(R), idxMax);
    pointPrompt = [px py];

    backgroundPoints = [
        1 1;
        round(wR/2) 1;
        wR 1;
        1 round(hR/2);
        wR round(hR/2);
        1 hR;
        round(wR/2) hR;
        wR hR
    ];

    pointPromptG = [pointPrompt(1) + (x1-1), pointPrompt(2) + (y1-1)];
    backgroundPointsG = [backgroundPoints(:,1) + (x1-1), backgroundPoints(:,2) + (y1-1)];

    % reduzir ROI para acelerar SAM (ajusta prompts e depois volta ao tamanho original)
    Rsmall = imresize(R, scale, "bilinear");
    R8 = uint8(255 .* rescale(Rsmall));
    imageSize = size(R8);

    pointPromptS = round(pointPrompt * scale);
    pointPromptS(1) = max(1, min(imageSize(2), pointPromptS(1)));
    pointPromptS(2) = max(1, min(imageSize(1), pointPromptS(2)));

    backgroundPointsS = round(backgroundPoints * scale);
    backgroundPointsS(:,1) = max(1, min(imageSize(2), backgroundPointsS(:,1)));
    backgroundPointsS(:,2) = max(1, min(imageSize(1), backgroundPointsS(:,2)));

    % embeddings + segmentação SAM
    embeddings = extractEmbeddings(model, R8);

    [candMasks, scores] = segmentObjectsFromEmbeddings(model, embeddings, imageSize, ...
        ForegroundPoints=pointPromptS, BackgroundPoints=backgroundPointsS);

    if ndims(candMasks) == 3
        [~,k] = max(scores);
        Msmall = candMasks(:,:,k);
    else
        Msmall = candMasks;
    end

    Mroi = imresize(Msmall, [hR wR], "nearest");

    % colocar máscara no tamanho total da slice
    maskFull = false(H,W);
    maskFull(y1:y2, x1:x2) = Mroi;
    masks{i} = rot90(maskFull,2);

    % mostrar resultado
    subplot(1,5,i);
    imagesc(Sn); axis image off; colormap gray; hold on;
    rectangle('Position',[x1 y1 (x2-x1+1) (y2-y1+1)], 'EdgeColor','y', 'LineWidth',1);
    contour(maskFull, [0.5 0.5], 'r', 'LineWidth', 1);
    plot(pointPromptG(1), pointPromptG(2), 'g*');
    plot(backgroundPointsG(:,1), backgroundPointsG(:,2), 'r*');
    title(sprintf('z=%d', z(i)));
end
