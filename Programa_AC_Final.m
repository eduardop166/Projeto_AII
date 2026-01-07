% PROJETO 2 - CC segmentation (baseline Active Contour) - VERSÃO MELHORADA

[hdr, path] = uigetfile('*.hdr');
V = double(analyze75read(fullfile(path, hdr)));

fprintf('A detetar slices com Corpus Callosum...\n');

[H, W, numSlices] = size(V);

scores = zeros(numSlices, 1);
for zz = 1:numSlices
    S = squeeze(V(:,:,zz));
    S = rot90(S, 2);
    
    p = prctile(S(:), [2 98]);
    Sn = (S - p(1)) / (p(2) - p(1) + eps);
    Sn = min(max(Sn, 0), 1);
    
    y1 = round(0.25*H); y2 = round(0.35*H);
    x1 = round(0.35*W); x2 = round(0.70*W);
    
    roi = Sn(y1:y2, x1:x2);
    
    brightPixels = roi(roi > prctile(roi(:), 75));
    scores(zz) = mean(brightPixels);
end

scoresSmooth = movmean(scores, 5);

zMin = 60;
zMax = 80;
win  = 8;

scoresSmoothWin = scoresSmooth;
scoresSmoothWin(1:zMin-1) = -Inf;
scoresSmoothWin(zMax+1:end) = -Inf;

validStarts = zMin:(zMax - win + 1);

winScore = zeros(size(validStarts));
for k = 1:numel(validStarts)
    s = validStarts(k);
    winScore(k) = mean(scoresSmoothWin(s:s+win-1));
end

[~, bestIdx] = min(winScore);      
startBest = validStarts(bestIdx);     
z = startBest:(startBest + win - 1);

% PASTAS para resultados
outRoot = fullfile(pwd, ['RESULTADOS_CC_' datestr(now,'yyyymmdd_HHMMSS')]);
%mkdir(outRoot);

dirGraph = fullfile(outRoot, 'GRAFICO');
dirRaw   = fullfile(outRoot, 'FIG_RAW');
dirAC    = fullfile(outRoot, 'FIG_AC');
dirMask  = fullfile(outRoot, 'FIG_MASK');
dirBin   = fullfile(outRoot, 'FIG_BIN');

%mkdir(dirGraph); mkdir(dirRaw); mkdir(dirAC); mkdir(dirMask); mkdir(dirBin);

% ARRAYS para resultados (uma linha por slice selecionado)
nSel = numel(z);
Area = nan(nSel,1);
Ecc  = nan(nSel,1);
BBx  = nan(nSel,1);
BBy  = nan(nSel,1);
BBw  = nan(nSel,1);
BBh  = nan(nSel,1);

fprintf('Slices selecionados: %s\n', mat2str(z));

figScore = figure('Name', 'Deteção de Slices');
plot(1:numSlices, scoresSmooth, 'b-', 'LineWidth', 1.5); hold on;
plot(z, scoresSmooth(z), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Slice number'); ylabel('Score (intensidade)');
title('Deteção automática de slices com CC');
legend('Score suavizado', 'Slices selecionados');
grid on;
%exportgraphics(figScore, fullfile(dirGraph,'grafico_scores.png'), 'Resolution', 300);
%savefig(figScore, fullfile(dirGraph,'grafico_scores.fig')); % opcional


% Mostrar as 9 slices com o ROI onde se calcula o score
figRaw = figure('Name','9 slices onde se calcula o score');

% ROI usada no score


for i = 1:numel(z)
    Sraw = squeeze(V(:,:,z(i)));
    Sraw = rot90(Sraw, 2);  % igual ao teu fluxo

    subplot(2,4,i);
    imagesc(Sraw); axis image off; colormap gray;
    hold on;
    %rectangle('Position',[x1, y1, (x2-x1), (y2-y1)], ...
    %          'EdgeColor','r','LineWidth',2);
    title(sprintf('z=%d', z(i)));
    putSliceLabel(z(i));
    %exportgraphics(gca, fullfile(dirRaw, sprintf('raw_z%03d.png', z(i))), 'Resolution', 300);
    hold off;
end
%exportgraphics(figRaw, fullfile(dirRaw,'FIG_RAW_mosaico.png'), 'Resolution', 300);
%savefig(figRaw, fullfile(dirRaw,'FIG_RAW_mosaico.fig')); % opcional



masksSeed  = cell(1, numel(z));
masksFinal = cell(1, numel(z));

% parâmetros 
sigma = 1.3;        % suavização gaussiana (0.8-2.0)
percSeed = 92;      % percentil para seed (88-95)
itersFirst = 110;   % iterações primeira slice
itersNext  = 110;    % iterações seguintes (porque usa seed anterior)
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


    % ROI para CC

    x1 = round(0.34*W); x2 = round(0.75*W);
    y1 = round(0.28*H); y2 = round(0.60*H);
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

        % dilatação para dar "margem" ao contorno
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
        bw = activecontour(R, init, itersNext, "Chan-vese","SmoothFactor", 1.8, ContractionBias=0.2);
    end

    % pós-processamento: manter componente mais central
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw, minArea);

    % colocar no tamanho total da slice
    maskFull = false(H,W);
    maskFull(y1:y2, x1:x2) = bw;
    
    % guardar máscara para seed (active contour)
    masksSeed{i} = maskFull;
 
    % recortar a máscara e binarizar

    % 1) recorte (fora da máscara -> 0)
    SnMasked = Sn;
    SnMasked(~maskFull) = 0;
    
    % --- 1.5) esticar histograma (contraste) SÓ dentro da máscara ---
    pix = SnMasked(maskFull);
    
    % proteção: se a máscara for pequena / vazia
    if numel(pix) > 50
        % limites robustos (evita ser afetado por outliers)
        lo = prctile(pix, 2);
        hi = prctile(pix, 98);
    
        % aplicar estiramento (clipping + remapeamento)
        SnCE = SnMasked;  % CE = contrast enhanced
        SnCE(maskFull) = imadjust(pix, [lo hi], [0 1]);
    else
        SnCE = SnMasked;
    end
    
    % 2) threshold depois do contraste 
    pix2 = SnCE(maskFull);
    thr = 0.30;               
    
    bwFinal = (SnCE > thr) & maskFull;

    
    % 3) limpeza leve
    bwFinal = bwareaopen(bwFinal, minArea);
    bwFinal = imfill(bwFinal,'holes');
    % manter só o maior objeto
    CCf = bwconncomp(bwFinal);
    if CCf.NumObjects > 0
        statsF = regionprops(CCf,'Area');
        areasF = [statsF.Area];
        [~, kF] = max(areasF);
    
        bwKeep = false(size(bwFinal));
        bwKeep(CCf.PixelIdxList{kF}) = true;
        st = regionprops(bwFinal, 'Area', 'Eccentricity', 'BoundingBox');
        if ~isempty(st)
            Area(i) = st.Area;
            Ecc(i)  = st.Eccentricity;
            bb      = st.BoundingBox;   % [x y w h]
            BBx(i)  = bb(1);
            BBy(i)  = bb(2);
            BBw(i)  = bb(3);
            BBh(i)  = bb(4);
        end
        bwFinal = bwKeep;
    end

    
    % guardar máscara final
    masksFinal{i} = bwFinal;

    % FIGURA 1: Sn + mask AC
    figure(figAC);
    subplot(2,4,i);
    imagesc(Sn); axis image off; colormap gray; hold on;
    contour(maskFull, [0.5 0.5], 'r', 'LineWidth', 1);
    rectangle('Position',[x1, y1, (x2-x1), (y2-y1)], 'EdgeColor','b','LineWidth',2);
    putSliceLabel(z(i));
    title(sprintf('Slice %d (z=%d)', i, z(i)));
    hold off;
    
    % guardar imagem
    %exportgraphics(gca, fullfile(dirAC, sprintf('ac_z%03d.png', z(i))), 'Resolution', 300);


    % FIGURA 2: Sn recortada pela máscara
    figure(figMask);
    subplot(2,4,i);
    imagesc(SnCE); axis image off; colormap gray; hold on;
    putSliceLabel(z(i));
    title(sprintf('Slice z=%d', z(i)));
    hold off;
    
    %exportgraphics(gca, fullfile(dirMask, sprintf('mask_z%03d.png', z(i))), 'Resolution', 300);
    
    % FIGURA 3: Depois da binarização
    figure(figBin);
    subplot(2,4,i);
    imagesc(masksFinal{i}); axis image off; colormap gray; hold on;
    putSliceLabel(z(i));
    title(sprintf('Binarizado z=%d', z(i)));
    hold off;
    
    %exportgraphics(gca, fullfile(dirBin, sprintf('bin_z%03d.png', z(i))), 'Resolution', 300);


end
% exportgraphics(figAC,   fullfile(dirAC,  'FIG_AC_mosaico.png'),   'Resolution', 300);
% exportgraphics(figMask, fullfile(dirMask,'FIG_MASK_mosaico.png'), 'Resolution', 300);
% exportgraphics(figBin,  fullfile(dirBin, 'FIG_BIN_mosaico.png'),  'Resolution', 300);
% 
% savefig(figAC,   fullfile(dirAC,  'FIG_AC_mosaico.fig')); 
% savefig(figMask, fullfile(dirMask,'FIG_MASK_mosaico.fig')); 
% savefig(figBin,  fullfile(dirBin, 'FIG_BIN_mosaico.fig'));  
% 
% % TABELA
% 
% SliceTxt = string(z(:));
% T = table(SliceTxt, Area, Ecc, BBx, BBy, BBw, BBh, ...
%     'VariableNames', {'Slice','Area','Excentricidade','BB_x','BB_y','BB_largura','BB_altura'});
% 
% meanRow = table("Media", mean(Area,'omitnan'), mean(Ecc,'omitnan'), ...
%     mean(BBx,'omitnan'), mean(BBy,'omitnan'), mean(BBw,'omitnan'), mean(BBh,'omitnan'), ...
%     'VariableNames', T.Properties.VariableNames);
% 
% T2 = [T; meanRow];
% writetable(T2, fullfile(outRoot,'metricas_CC.xlsx'));


%% Funcoes auxiliares

function putSliceLabel(z)
    % escreve "z=XX" no canto superior esquerdo
    text(0.02, 0.98, sprintf('z=%d', z), ...
        'Units','normalized','Color','w','FontSize',12,'FontWeight','bold', ...
        'BackgroundColor','k','Margin',2,'Interpreter','none', ...
        'VerticalAlignment','top','HorizontalAlignment','left');
end
