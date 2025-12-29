% PROJETO 2 - CC segmentation (baseline Active Contour)

% escolher ficheiro hdr (tem que estar o .img tambem na pasta )
[hdr, path] = uigetfile('*.hdr');
V = double(analyze75read(fullfile(path, hdr)));

% slices em que aparece o CC (escolhidas manualmente)
z = round(linspace(64, 74, 10));

% guardar máscaras
masksSeed  = cell(1, numel(z));   % máscara para seed (do active contour)
masksFinal = cell(1, numel(z));   % máscara final (depois do recorte + binarização)


% parâmetros 
sigma = 1.3;        % suavização gaussiana (0.8-2.0)
percSeed = 92;      % percentil para seed (88-95)
itersFirst = 70;   % iterações primeira slice
itersNext  = 100;    % iterações seguintes (porque usa seed anterior)
minArea = 120;      % remove blobs pequenos

figAC   = figure('Name','(1) Sn + Active Contour mask');
figMask = figure('Name','(2) Sn recortada pela máscara');
figBin  = figure('Name','(3) Depois da binarização');

for i = 1:numel(z)

    % obter slice e rodar 
    S = squeeze(V(:,:,z(i)));
    S = rot90(S,2);   % 180 graus

    % pré-processamento: normalizar + suavizar
    p = prctile(S(:), [2 98]);
    Sn = (S - p(1)) / (p(2)-p(1) + eps);
    Sn = min(max(Sn,0),1);
    
    % melhor que gaussian forte: suaviza mas preserva bordas
    Sn = imbilatfilt(Sn, 0.06, 3);   % ajusta: (0.04–0.10) e (2–5)

    Sn = imgaussfilt(Sn, 0.6);
    
    [H,W] = size(Sn);


    % ROI (dinâmica a partir da slice anterior)
    if i > 1 && ~isempty(masksSeed{i-1}) && nnz(masksSeed{i-1}) > 0
        prevFull = masksSeed{i-1};
        [rr, cc] = find(prevFull);
        pad = 15; % margem em pixels
    
        y1 = max(min(rr)-pad, 1);  y2 = min(max(rr)+pad, H);
        x1 = max(min(cc)-pad, 1);  x2 = min(max(cc)+pad, W);
    else
        % ROI central na 1ª slice (ou se a anterior falhar)
        x1 = round(0.30*W); x2 = round(0.70*W);
        y1 = round(0.30*H); y2 = round(0.70*H);
    end
    
    R = Sn(y1:y2, x1:x2);

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
        initPrev = masksSeed{i-1};
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
        bw = activecontour(R, init, itersNext, "Chan-Vese","SmoothFactor", 1.8, "ContractionBias", 0.15);
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
    
    % guardar máscara para seed (active contour)
    masksSeed{i} = maskFull;
    
    % recortar a máscara e binarizar
    
    % 1) recorte (fora da máscara -> 0)
    SnMasked = Sn;
    SnMasked(~maskFull) = 0;
    
    % 2) threshold (Otsu calculado só dentro da máscara)
    pix = SnMasked(maskFull);
    thr= 0.7;
    
    bwFinal = (SnMasked > thr) & maskFull;
    
    % 3) limpeza leve
    bwFinal = bwareaopen(bwFinal, minArea);
    bwFinal = imfill(bwFinal,'holes');
    
    % guardar máscara final
    masksFinal{i} = bwFinal;


    % FIGURA 1: Sn + mask AC (como antes)
    figure(figAC);
    subplot(2,5,i);
    imagesc(Sn); axis image off; colormap gray; hold on;
    contour(maskFull, [0.5 0.5], 'r', 'LineWidth', 1);
    title(sprintf('Slice %d (z=%d)', i, z(i)));
    hold off;
    
    % FIGURA 2: Sn recortada pela máscara (fora = 0)
    figure(figMask);
    subplot(2,5,i);
    imagesc(SnMasked); axis image off; colormap gray;
    title(sprintf('Recorte z=%d', z(i)));

    % FIGURA 3: Depois da binarização (máscara final)

    figure(figBin);
    subplot(2,5,i);
    imagesc(masksFinal{i}); axis image off; colormap gray;
    title(sprintf('Binarizado z=%d', z(i)));
    
    drawnow

end