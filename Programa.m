% PROJETO 2

% escolher ficheiro hr ( tem que estar o .img tambem na pasta ) 
[hdr, path] = uigetfile('*.hdr');
V = analyze75read(fullfile(path, hdr));

% 5 slices entre 60 e 80
z = round(linspace(55, 75, 10));

figure;
for i = 1:10
    S = squeeze(V(:,:,z(i)));   % RODA 180 graus (dataset ta ao contrario)
    S = rot90(S,1);             
    S = rot90(S,1);  
    subplot(2,5,i)
    imagesc(S)
    axis image off
    colormap gray
    title(sprintf('Axial %d', z(i)))
end
