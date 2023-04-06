clear
close all

dataFolder = 'W:\05_Topography\Data';
saveFolder = 'W:\05_Topography\Outputs';
mkdir(saveFolder)
toSave = true;

imFiles = dir(fullfile(dataFolder,'Images','*.mat'));
featureMaps_All = [];   % To save all the feature maps across all the images
patientNums = [];
fprintf('\n')
for imNum = 1:numel(imFiles)
    ptNum = imFiles(imNum).name(8:10);
    fprintf('(%03d/%03d) P%s',imNum,numel(imFiles),ptNum)
    patientNums = [patientNums; ptNum];

    % Reading the image:
    load(fullfile(dataFolder,'Images',imFiles(imNum).name));
    fprintf('    (%1.1f x %1.1f x %1.1f)\n', pixDim(1), pixDim(2), pixDim(3))

    % Reading the Segmentation Labels:
    load(fullfile(dataFolder,'Labels',sprintf('Label_0%s_placenta.mat',ptNum)))
    if exist('mrLabel','var')
        pLabel = logical(mrLabel);   %Placenta Label
    elseif exist('plLabel','var')
        pLabel = logical(plLabel);   %Placenta Label
    end
    load(fullfile(dataFolder,'Labels',sprintf('Label_0%s_Uterus.mat',ptNum)),'utLabel')
    uLabel = logical(utLabel);   %Uterus Label
    clear mrLabel plLabel utLabel
    
    %% Up-sampling the image and the labels (isotropic):
    pSlices = find(sum(sum(pLabel,1),2));
    uSlices = find(sum(sum(uLabel,1),2));
    if pSlices(1) > 1
        [pCx1,pCy1] = ind2sub(size(pLabel(:,:,pSlices(1))), find(pLabel(:,:,pSlices(1))));
        pC1 = round(mean([pCx1,pCy1])); % Placenta center at the first slice
        pLabel(pC1(1),pC1(2),pSlices(1)-1) = 1;
    end
    if pSlices(end) < size(pLabel,3)
        [pCx2,pCy2] = ind2sub(size(pLabel(:,:,pSlices(end))), find(pLabel(:,:,pSlices(end))));
        pC2 = round(mean([pCx2,pCy2])); % Placenta center at the last slice
        pLabel(pC2(1),pC2(2),pSlices(end)+1) = 1;
    end
    
    if uSlices(1) > 1
        [uCx1,uCy1] = ind2sub(size(uLabel(:,:,uSlices(1))), find(uLabel(:,:,uSlices(1))));
        uC1 = round(mean([uCx1,uCy1])); % Uterus center at the first slice
        uLabel(uC1(1),uC1(2),uSlices(1)-1) = 1;
    end
    if uSlices(end) < size(uLabel,3)
        [uCx2,uCy2] = ind2sub(size(uLabel(:,:,uSlices(end))), find(uLabel(:,:,uSlices(end))));
        uC2 = round(mean([uCx2,uCy2])); % Uterus center at the last slice
        uLabel(uC2(1),uC2(2),uSlices(end)+1) = 1;
    end
    orig_z_size = size(mrImage,3);
    
    upCoe = 1;
    [xq,yq,zq] = meshgrid(1:upCoe:size(mrImage,2), 1:upCoe*pixDim(1)/pixDim(2):size(mrImage,1), 1:upCoe*pixDim(1)/pixDim(3):size(mrImage,3));
    mrImage = interp3(double(mrImage),xq,yq,zq,'cubic');
    pLabel = level_set_Interp3(double(pLabel),xq,yq,zq,'cubic');
    uLabel = level_set_Interp3(double(uLabel),xq,yq,zq,'cubic');
    
    pixDimIso = upCoe * pixDim(1);

    [originX, originY, originZ] = ind2sub(size(uLabel),find(uLabel));
    uCenter = round(mean([originX, originY, originZ]));  %Point of Observation (uterus center of mass)
    [pcentX, pcentY, pcentZ] = ind2sub(size(pLabel),find(pLabel));
    pCenter = round(mean([pcentX, pcentY, pcentZ]));  %Placenta center of mass
    
    % center placenta
    trans_vec = [128, 128, 128]-pCenter;
    mrImage = imtranslate(mrImage, [trans_vec(2) trans_vec(1) trans_vec(3)], 'FillValues', 0);
    pLabel = imtranslate(pLabel, [trans_vec(2) trans_vec(1) trans_vec(3)], 'FillValues', 0);
    
    % pad to 256x256x256
    if size(mrImage, 3) < 256
        pad_amt = ceil((256 - size(mrImage, 3))/2);
        mrImage = padarray(mrImage, [0 0 pad_amt], 0, 'both');
        pLabel = padarray(pLabel, [0 0 pad_amt], 0, 'both');
    end
    
    % resize to 256x256x256
    if size(mrImage, 3) > 256
        crop_amt = ceil((size(mrImage,3)-256)/2);
        mrImage = mrImage(:, :, crop_amt:crop_amt+255);
        pLabel = pLabel(:, :, crop_amt:crop_amt+255);
    end

    [dists_sag_left, dists_sag_right] = get_dists(pLabel, 1); % distances along x axis from patient left side
    [dists_cor_ant, dists_cor_post] = get_dists(pLabel, 2); % distances along y axis from patient anterior
    
    intensity_sag_left = zeros(size(dists_sag_left));
    intensity_sag_right = zeros(size(dists_sag_left));
    for y = 1:256
        for z = 1:256
            if dists_sag_left(y, z) > 1
                intensity_sag_left(y, z) = mrImage(dists_sag_left(y,z), y, z);
            end
            if dists_sag_right(y, z) > 1
                intensity_sag_right(y, z) = mrImage(dists_sag_right(y,z), y, z);
            end
        end
    end
    
    intensity_cor_ant = zeros(size(dists_cor_ant));
    intensity_cor_post = zeros(size(dists_cor_ant));
    for x = 1:256
        for z = 1:256
            if dists_cor_ant(x, z) > 1
                intensity_cor_ant(x, z) = mrImage(dists_cor_ant(x,z), x, z);
            end
            if dists_cor_post(x, z) > 1
                intensity_cor_post(x, z) = mrImage(dists_cor_post(x,z), x, z);
            end
        end
    end
    
    mod_dists_sag_left = pixDim(1).*dists_sag_left.';
    mod_dists_sag_right = pixDim(1).*dists_sag_right.';
    thickness_sag = mod_dists_sag_right - mod_dists_sag_left;
    mod_dists_cor_ant = pixDim(2).*dists_cor_ant.';
    mod_dists_cor_post = pixDim(2).*dists_cor_post.';
    thickness_cor = mod_dists_cor_post - mod_dists_cor_ant;
    intensity_sag_left = intensity_sag_left.';
    intensity_sag_right =  intensity_sag_right.';
    intensity_cor_ant = intensity_cor_ant.';
    intensity_cor_post = intensity_cor_post.';
    
    
    %% Display and save
    imSaveFolder = strcat(saveFolder, '\', imFiles(imNum).name(1:end-4));
    if toSave
        if ~exist(imSaveFolder, 'dir')
            mkdir(imSaveFolder);
        end
    end
    
    figure('Position', [1, 41, 1920, 1080]),
    imagesc(mod_dists_sag_left)
    ax = gca;
    ax.XTick = 0:20:pixDim(1)*size(mrImage,1);
    ax.XTickLabel = num2cell(0:20:pixDim(1)*256);
    ax.YTick = 0:20:pixDim(3)*orig_z_size;
    ax.YTickLabel = num2cell(0:20:pixDim(3)*256);
    ax.FontSize = 14;
    colorbar
    colormap hot
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Ant.','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Post.','Color','Blue','FontSize',16)
    text(256,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Sagittal left side topography map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Map_Left.png']))
    end
    
    figure('Position', [1, 41, 1920, 1080]),
    imagesc(mod_dists_sag_right)
    ax = gca;
    ax.XTick = 0:20:pixDim(1)*size(mrImage,1);
    ax.XTickLabel = num2cell(0:20:pixDim(1)*256);
    ax.YTick = 0:20:pixDim(3)*orig_z_size;
    ax.YTickLabel = num2cell(0:20:pixDim(3)*256);
    ax.FontSize = 14;
    colorbar
    colormap hot
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Ant.','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Post.','Color','Blue','FontSize',16)
    text(256,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Sagittal right side  topography map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Map_Right.png']))
    end
    
    figure('Position', [1, 41, 1920, 1080]),
    imagesc(thickness_sag)
    ax = gca;
    ax.XTick = 0:20:pixDim(1)*size(mrImage,1);
    ax.XTickLabel = num2cell(0:20:pixDim(1)*256);
    ax.YTick = 0:20:pixDim(3)*orig_z_size;
    ax.YTickLabel = num2cell(0:20:pixDim(3)*256);
    ax.FontSize = 14;
    colorbar
    colormap hot
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Ant.','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Post.','Color','Blue','FontSize',16)
    text(256,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Sagittal thickness (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Thickness_Sag.png']))
    end
    
    figure('Position', [1, 41, 1920, 1080]),
    imagesc(mod_dists_cor_ant)
    ax = gca;
    ax.XTick = 0:20:pixDim(1)*size(mrImage,1);
    ax.XTickLabel = num2cell(0:20:pixDim(1)*256);
    ax.YTick = 0:20:pixDim(3)*orig_z_size;
    ax.YTickLabel = num2cell(0:20:pixDim(3)*256);
    ax.FontSize = 14;
    colorbar
    colormap hot
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Right','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Left','Color','Blue','FontSize',16)
    text(256,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Coronal anterior side topography map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Map_Ant.png']))
    end
    
    figure('Position', [1, 41, 1920, 1080]),
    imagesc(mod_dists_cor_post)
    ax = gca;
    ax.XTick = 0:20:pixDim(1)*size(mrImage,1);
    ax.XTickLabel = num2cell(0:20:pixDim(1)*256);
    ax.YTick = 0:20:pixDim(3)*orig_z_size;
    ax.YTickLabel = num2cell(0:20:pixDim(3)*256);
    ax.FontSize = 14;
    colorbar
    colormap hot
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Right','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Left','Color','Blue','FontSize',16)
    text(256,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Coronal posterior side topography map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Map_Post.png']))
    end
    
    figure('Position', [1, 41, 1920, 1080]),
    imagesc(thickness_cor)
    ax = gca;
    ax.XTick = 0:20:pixDim(1)*size(mrImage,1);
    ax.XTickLabel = num2cell(0:20:pixDim(1)*256);
    ax.YTick = 0:20:pixDim(3)*orig_z_size;
    ax.YTickLabel = num2cell(0:20:pixDim(3)*256);
    ax.FontSize = 14;
    colorbar
    colormap hot
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Right','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Left','Color','Blue','FontSize',16)
    text(256,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Coronal thickness (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Thickness_Cor.png']))
    end
    
    figure,
    imshow(intensity_sag_left, [])
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Ant.','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Post.','Color','Blue','FontSize',16)
    xlabel(['Sagittal left side intensity map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Intensity_Left.png']))
    end
    
    figure,
    imshow(intensity_sag_right, [])
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Ant.','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Post.','Color','Blue','FontSize',16)
    xlabel(['Sagittal right side intensity map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Intensity_Right.png']))
    end
    
    figure,
    imshow(intensity_cor_ant, [])
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Right','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Left','Color','Blue','FontSize',16)
    xlabel(['Coronal anterior side intensity map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Intensity_Ant.png']))
    end
    
    figure,
    imshow(intensity_cor_post, [])
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    text(123,5,'Sup.','Color','Blue','FontSize',16)
    text(5,123,'Right','Color','Blue','FontSize',16)
    text(123,251,'Inf.','Color','Blue','FontSize',16)
    text(240,123,'Left','Color','Blue','FontSize',16)
    xlabel(['Coronal posterior side intensity map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(imSaveFolder,[imFiles(imNum).name(1:end-4),'_Intensity_Post.png']))
    end
    
    
    
    % Feature maps order: (1)topography fetal side, (2)topography maternal side,
    %                     (3)thickness, (4)intensity fetal side,
    %                     (5)intensity maternal side, (6)average intensity
    %                     fetal side, (7)average intensity maternal side,
    %                     (8)intensity std fetal side, (9)intensity std
    %                     materal side, (10)entropy fetal side, (11)entropy
    %                     maternal side
%     featureMaps_All = cat(4,featureMaps_All,cat(3,topoMapF,topoMapM,thicknessMap,...
%                                                   intensityMapF,intensityMapM,...
%                                                   intAvMapF,intAvMapM,...
%                                                   intStdMapF,intStdMapM,...
%                                                   entropyMapF,entropyMapM));
    
    %pause(1)
    close all
end
% fMapOrder = {'(1)topography fetal side', '(2)topography maternal side', '(3)thickness', '(4)intensity fetal side',...
%                          '(5)intensity maternal side', '(6)average intensity fetal side', '(7)average intensity maternal side',...
%                          '(8)intensidt std fetal side', '(9)intensity std materal side', '(10)entropy fetal side',...
%                          '(11)entropy maternal side'};
% if toSave == true
%     save(fullfile(saveFolder,'FeatureMaps.mat'),'patientNums','featureMaps_All','patchSize_mm','fMapOrder')
% end