% PROJETO 2 - CC segmentation (baseline Active Contour)

% escolher ficheiro hdr (tem que estar o .img tambem na pasta )
[hdr, path] = uigetfile('*.hdr');
V = double(analyze75read(fullfile(path, hdr)));

% slices em que aparece o CC (escolhidas manualmente)
z = round(linspace(61, 71, 10));

% guardar máscaras
masks = cell(1, numel(z));

% parâmetros 
sigma = 1.2;        % suavização gaussiana (0.8-2.0)
percSeed = 92;      % percentil para seed (88-95)
itersFirst = 70;   % iterações primeira slice
itersNext  = 100;    % iterações seguintes (porque usa seed anterior)
minArea = 120;      % remove blobs pequenos

figure;
for i = 1:numel(z)

    % obter slice e rodar 
    S = squeeze(V(:,:,z(i)));
    S = rot90(S,2);   % 180 graus

    % pré-processamento: normalizar + suavizar
    p = prctile(S(:), [2 98]);
    Sn = (S - p(1)) / (p(2)-p(1) + eps);
    Sn = min(max(Sn,0),1);
    Sn = imgaussfilt(Sn, sigma);

    [H,W] = size(Sn);

    % ROI central
    % (CC costuma estar na zona central, um pouco acima do meio)
    x1 = round(0.30*W); x2 = round(0.70*W);
    y1 = round(0.30*H); y2 = round(0.70*H);
    R  = Sn(y1:y2, x1:x2);

    % seed
    if i == 1
        % seed automática na primeira slice
        t = prctile(R(:), percSeed);  % CC mais claro
        init = R > t;
        init = imopen(init, strel('disk',2));
        init = imclose(init, strel('disk',5));
        init = imfill(init,'holes');
        init = bwareaopen(init, minArea);
    else
        % seed = máscara anterior
        initPrev = masks{i-1};
        init = initPrev(y1:y2, x1:x2);

        % pequena dilatação para dar "margem" ao contorno
        init = imdilate(init, strel('disk',2));
    end

    % Se o init ficar vazio, volta a seed automática
    if nnz(init) < 50
        t = prctile(R(:), percSeed);
        init = R > t;
        init = imopen(init, strel('disk',2));
        init = imclose(init, strel('disk',5));
        init = imfill(init,'holes');
        init = bwareaopen(init, minArea);
    end

    % Active Contour (Chan–Vese)
    if i == 1
        bw = activecontour(R, init, itersFirst, "Chan-Vese");
    else
        bw = activecontour(R, init, itersNext, "Chan-Vese");
    end

    % pós-processamento: manter componente mais central
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw, minArea);

    CC = bwconncomp(bw);
    if CC.NumObjects > 0
        stats = regionprops(CC,'Area','Centroid');
        cent = vertcat(stats.Centroid);
        area = vertcat(stats.Area);

        c0 = [size(R,2)/2, size(R,1)/2];
        d  = hypot(cent(:,1)-c0(1), cent(:,2)-c0(2));

        [~,k] = max(area .* (1./(1+d))); % grande e central
        bw2 = false(size(bw));
        bw2(CC.PixelIdxList{k}) = true;
        bw = bw2;
    end

    % colocar no tamanho total da slice
    maskFull = false(H,W);
    maskFull(y1:y2, x1:x2) = bw;
    masks{i} = maskFull;

    % Mostrar
    subplot(2,5,i);
    imagesc(Sn); axis image off; colormap gray; hold on;
    contour(maskFull, [0.5 0.5], 'r', 'LineWidth', 1);
    title(sprintf('Slice %d (z=%d)', i, z(i)));
    drawnow
end