close all; clear; clc;

jcfolder='JC_CP_data/JC_CP_';
cps=[1,2,3,4,5,6,10,11,12,13];
bscfolders={
    'JC_CP1_2_40_40dB';
    'JC_CP2_40_40dB';
    'JC_CP3_40_44dB';
    'JC_CP4_40_44dB';
    'JC_CP5_40_44dB';
    'JC_CP6_40_44dB';
    'JC_CP10_40_44dB';
    'JC_CP11_40_44dB';
    'JC_CP12_40_45dB';
    'JC_CP13_40_45dB';
    };

bscdataJCCP=struct();

for j=1:length(bscfolders)
    disp(['CP ',num2str(j)]);
% for j=1
for i=1:11
    disp(['   scan ',num2str(i)])

    fileroi=[jcfolder,num2str(cps(j)),'/BSC/',bscfolders{j},'_ROI_',num2str(i),'.mat'];
    filebsc=[jcfolder,num2str(cps(j)),'/BSC/',bscfolders{j},'/',bscfolders{j},'_bscFULL-',num2str(i),'.mat'];
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

        bscdataJCCP.(['CP',num2str(j)]).(['scan',num2str(i)])=UIB_grille_bsc;
    
    
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

save('bscdataJCCP','bscdataJCCP');