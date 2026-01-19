close all; clear; clc;

jcfolder='preprocessing/4T1_CP/';
scans=1:11;
bscfolders=46:62;

bscdata4T1CP=struct();

for j=1:length(bscfolders)
    disp(['CP ',num2str(bscfolders(j))]);

for i=scans
    disp(['   scan ',num2str(i)])

    % Build folder path
        folderPath = fullfile(jcfolder, ['A', sprintf('%03d', bscfolders(j))]);

    % Find files using dir
    bscFile = dir(fullfile(folderPath, ['*bscFULL-', num2str(i), '.mat']));
    roiFile = dir(fullfile(folderPath, ['*ROI_', num2str(i), '.mat']));

    % Safety check
    if isempty(bscFile) || isempty(roiFile)
        warning('Missing files for CP %d scan %d', bscfolders(j), i);
        continue
    end

    % Load files
    load(fullfile(folderPath, bscFile(1).name));
    load(fullfile(folderPath, roiFile(1).name));

    % disp(bscFile(1).name)
    % disp(roiFile(1).name)
    % 
    % continue

    pos_l=ROI_positions.left;
    pos_r=ROI_positions.right;
    pos_b=ROI_positions.bottom;
    pos_t=ROI_positions.top;

    if ~isempty(pos_l)
        pos_meanx=pos_l+(pos_r-pos_l)/2;
        pos_meany=pos_t+(pos_b-pos_t)/2;
    
        nb_posx=(max(pos_l)-min(pos_l))./min(abs(diff(pos_l)));
        nb_posxint=int16(nb_posx);
        nb_posyint=length(pos_meany([1; 1+find(diff(pos_meany))]));
        % UIB_grille=zeros(nb_posxint,nb_posyint,4);
        UIB_grille=zeros(nb_posxint,nb_posyint);
        UIB_grille_bsc=zeros(nb_posxint,nb_posyint,128);
    
        for ii=1:length(pos_meanx)
            if pos_l(ii)>0
                indf=int16(pos_l(ii)./min(abs(diff(pos_l))));
                indc=round((pos_b(ii)-pos_b(1))./(max(diff(pos_b))+eps)+1);
                % UIB_grille(indf,indc,:)=UIB(:,ii);
                UIB_grille(indf,indc)=ii;
                UIB_grille_bsc(indf,indc,:)=BSC(:,ii);
            end
        end

        % remove single zero rows
        zr = all(all(UIB_grille_bsc==0,3),2); % zero rows
        % UIB_grille(zr & [true; ~zr(1:end-1)] & [~zr(2:end); true],:,:) = [];
        UIB_grille_bsc(zr & [true; ~zr(1:end-1)] & [~zr(2:end); true],:,:) = [];

        bscdata4T1CP.(['CP',num2str(bscfolders(j))]).(['scan',num2str(i)])=UIB_grille_bsc;
    
    
        % figure(j); subplot(4,3,i);
        % if isvector(UIB_grille)
        %     imagesc(UIB_grille>0);
        %     set(gca, 'YDir', 'normal');
        % else
        %     surf((UIB_grille>0)*1,'EdgeColor','none','LineStyle','none');
        % end
        % view(90,90);        
    end
end
end

save('data/bscdata4T1CP','bscdata4T1CP','-v7.3');