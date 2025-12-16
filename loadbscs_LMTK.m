close all; clear; clc;

jcfolder='LMTK_CP/';
scans=1:11;
bscfolders={
    'A130_40_48dB_';
    };

bscdataLMTKCP=struct();

for j=1:length(bscfolders)
    disp(['CP ',num2str(j)]);
% for j=1
for i=scans
    disp(['   scan ',num2str(i)])

    filebsc=[jcfolder,bscfolders{j},'bscFULL-',num2str(i),'.mat'];
    fileroi=[jcfolder,bscfolders{j},'ROI_',num2str(i),'.mat'];
    load(fileroi);
    load(filebsc);

    % % load('inversion_mean_CP1_fminconC1_fminsearch')
    % load('inversion_mean_CP1_fminconC1')
    % % load('inversion_mean_CP1_fminconC7') % up to C7

    % Params_SFMi=eval(['Params_SFM',num2str(i)]);
    % UIB=Params_SFMi;

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

        bscdataLMTKCP.(['CP',num2str(j)]).(['scan',num2str(i)])=UIB_grille_bsc;
    
    
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

save('bscdataLMTKCP','bscdataLMTKCP');