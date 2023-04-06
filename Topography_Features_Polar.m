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
        pLabel = single(mrLabel);   %Placenta Label
    elseif exist('plLabel','var')
        pLabel = single(plLabel);   %Placenta Label
    end
    load(fullfile(dataFolder,'Labels',sprintf('Label_0%s_Uterus.mat',ptNum)),'utLabel')
    uLabel = single(utLabel);   %Uterus Label
    clear mrLabel plLabel utLabel
    
    pLabel = logical(pLabel);
    uLabel = logical(uLabel);
    
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
    
    upCoe = 1;
    [xq,yq,zq] = meshgrid(1:upCoe:size(mrImage,2), 1:upCoe*pixDim(1)/pixDim(2):size(mrImage,1), 1:upCoe*pixDim(1)/pixDim(3):size(mrImage,3));
    mrImage = interp3(double(mrImage),xq,yq,zq,'cubic');
    pLabel = level_set_Interp3(double(pLabel),xq,yq,zq,'cubic');
    uLabel = level_set_Interp3(double(uLabel),xq,yq,zq,'cubic');
    
    %showCross(mrImage);
    pixDimIso = upCoe * pixDim(1);

    [originX, originY, originZ] = ind2sub(size(uLabel),find(uLabel));
    uCenter = round(mean([originX, originY, originZ]));  %Point of Observation (uterus center of mass)
    [pcentX, pcentY, pcentZ] = ind2sub(size(pLabel),find(pLabel));
    pCenter = round(mean([pcentX, pcentY, pcentZ]));  %Placenta center of mass
    bDist = 1/3;
    
    % For origin outside
%     Origin = uCenter;
%     v_u2p = pCenter - Origin;
%     Origin = round(pCenter + (1+bDist)*v_u2p);
%     fprintf(['New calc origin: ',mat2str(Origin),'\n'])
    
    
    % For origin inside uterus
    Origin = uCenter;
    Origin = round(((1+bDist)*Origin)-(bDist*pCenter));
    fprintf(['Original calc origin: ',mat2str(Origin)])
%     while pLabel(Origin(1),Origin(2),Origin(3))
%         Origin = round(((1+bDist)*Origin)-(bDist*pCenter));
%         fprintf('  (Adjusted point of observation!)')
%     end
    
    mrBorder = pLabel - (1-(bwdist(1-pLabel) <= 1));  % 3D Border
    % Border Points:
    [bX,bY,bZ] = ind2sub(size(mrBorder),find(mrBorder)); % Cartesian border
    [bAz, bEl, bR] = cart2sph(bX-Origin(1), bY-Origin(2), bZ-Origin(3)); % Polar border
    [pcAz, pcEl, pcR] = cart2sph(pCenter(1)-Origin(1), pCenter(2)-Origin(2), pCenter(3)-Origin(3)); % Placenta polar

    %%
    bAz_sort = sort(bAz);
    bAz_gap = bAz_sort(1:end-1)-bAz_sort(2:end);
    [gap_dist,gap_index] = max(abs(bAz_gap));
    if abs(gap_dist) > 0.2  % Any large gap between azimuth angles of the surface points?
        bAz(bAz <= bAz_sort(gap_index)) = bAz(bAz <= bAz_sort(gap_index)) + 2*pi;
        pcAz(pcAz <= bAz_sort(gap_index)) = pcAz(pcAz <= bAz_sort(gap_index)) + 2*pi;
    end

%     cAz = mean([min(bAz),max(bAz)]); cEl = mean([min(bEl),max(bEl)]);   % Center point of Placenta
    cAz = pcAz; cEl = mean([min(bEl),max(bEl)]);   % Center point of Placenta
    bAz = bAz - cAz;
    bAz(bAz > pi) = bAz(bAz > pi) - 2*pi;
    bAz(bAz < -pi) = 2*pi + bAz(bAz < -pi);
    bEl = bEl - cEl;
    bAz(bEl > pi/2) = pi + bAz(bEl > pi/2);
    bAz(bEl < -pi/2) = pi + bAz(bEl < -pi/2);
    bEl(bEl > pi/2) = pi - bEl(bEl > pi/2);
    bEl(bEl < -pi/2) = -pi - bEl(bEl < -pi/2);
    bAz(bAz > pi) = bAz(bAz > pi) - 2*pi;
    bAz(bAz < -pi) = 2*pi + bAz(bAz < -pi);

    linePos = [0 90 180 270] - (cAz*180/pi);
    linePos(linePos < 0) = linePos(linePos < 0) + 360;

    longestD = ceil(max(bR)*1.1);
    pixInt = mrImage(logical(mrBorder));   % border pixel intensities

    topoMapF = ones(181,361) * longestD;
    topoMapM = ones(181,361) * longestD;
    intensityMapF = nan(181,361);
    intensityMapM = nan(181,361);
    for c = 1:numel(bAz)
        X = round((bEl(c)*180/pi+91));
        Y = round((bAz(c)*180/pi+181));
        if topoMapF(X,Y) ==  longestD
            topoMapF(X,Y) = bR(c);
            intensityMapF(X,Y) = pixInt(c);
        elseif topoMapF(X,Y) > bR(c)
            topoMapF(X,Y) = bR(c);
            intensityMapF(X,Y) = pixInt(c);
        end        
        if topoMapM(X,Y) ==  longestD
            topoMapM(X,Y) = bR(c);
            intensityMapM(X,Y) = pixInt(c);
        elseif topoMapM(X,Y) < bR(c)
            topoMapM(X,Y) = bR(c);
            intensityMapM(X,Y) = pixInt(c);
        end        
    end
    pixDimIsoHlf = 0.5 * pixDimIso;
    topoMapF = topoMapF * pixDimIsoHlf;
    topoMapM = topoMapM * pixDimIsoHlf;
    thicknessMap = topoMapM-topoMapF;
%     topoMapF_t = flip(flip(topoMapF,1),2);
%     topoMapF_t = imtranslate(topoMapF_t,[181, 90]) + imtranslate(topoMapF_t,[-180, 90]) + imtranslate(topoMapF_t,[181, -91]) + imtranslate(topoMapF_t,[-180,-91]);
%     topoMap360 = topoMapM + topoMapF_t;
    

%% Feature extraction for Fetal side:    
    [centEl,centAz] = ind2sub(size(topoMapF),find(topoMapF <longestD*pixDimIsoHlf));
    centR = topoMapF(topoMapF < longestD*pixDimIsoHlf);
    centPsph = [centAz,centEl,centR];   % Surface points in spherical coordinate system
    centEl = ((centEl-91)*pi/180)+cEl;
    centAz = ((centAz-181)*pi/180)+cAz;
    [centX,centY,centZ] = sph2cart(centAz,centEl,centR);
    centPcart = round([centX,centY,centZ]+repmat(Origin,numel(centX),1)); % Surface points in cartesian coordinate system
    
    patchSize_mm = 15;  % patch dim. = (15 x 15 x 15) mm
    pRad = round((patchSize_mm/2)/(upCoe*pixDim(1)))-1;  
     
    entropyMapF = nan(181,361);
    intAvMapF = nan(181,361);
    intStdMapF = nan(181,361);
    for pNum = 1:size(centPcart,1)
        Patch = im2double(mrImage(max(centPcart(pNum,1)-pRad,1):min(centPcart(pNum,1)+pRad,size(mrImage,1)),...
                                  max(centPcart(pNum,2)-pRad,1):min(centPcart(pNum,2)+pRad,size(mrImage,2)),...
                                  max(centPcart(pNum,3)-pRad,1):min(centPcart(pNum,3)+pRad,size(mrImage,3))));
        entropyMapF(centPsph(pNum,2),centPsph(pNum,1)) = entropy(Patch);
        intAvMapF(centPsph(pNum,2),centPsph(pNum,1)) = mean(Patch(:));
        intStdMapF(centPsph(pNum,2),centPsph(pNum,1)) = std(Patch(:));
    end

%% Feature extraction for Maternal side:    
    [centEl,centAz] = ind2sub(size(topoMapM),find(topoMapM <longestD*pixDimIsoHlf));
    centR = topoMapM(topoMapM < longestD*pixDimIsoHlf);
    centPsph = [centAz,centEl,centR];   % Surface points in spherical coordinate system
    centEl = ((centEl-91)*pi/180)+cEl;
    centAz = ((centAz-181)*pi/180)+cAz;
    [centX,centY,centZ] = sph2cart(centAz,centEl,centR);
    centPcart = round([centX,centY,centZ]+repmat(Origin,numel(centX),1)); % Surface points in crtesian coordinate sytem
    
    entropyMapM = nan(181,361);
    intAvMapM = nan(181,361);
    intStdMapM = nan(181,361);
    for pNum = 1:size(centPcart,1)
        Patch = im2double(mrImage(max(centPcart(pNum,1)-pRad,1):min(centPcart(pNum,1)+pRad,size(mrImage,1)),...
                                  max(centPcart(pNum,2)-pRad,1):min(centPcart(pNum,2)+pRad,size(mrImage,2)),...
                                  max(centPcart(pNum,3)-pRad,1):min(centPcart(pNum,3)+pRad,size(mrImage,3))));
        entropyMapM(centPsph(pNum,2),centPsph(pNum,1)) = entropy(Patch);
        intAvMapM(centPsph(pNum,2),centPsph(pNum,1)) = mean(Patch(:));
        intStdMapM(centPsph(pNum,2),centPsph(pNum,1)) = std(Patch(:));
    end
    
% %% Feature extraction for 360 view:    
%     [centEl,centAz] = ind2sub(size(topoMapM),find(topoMapM <longestD*pixDimIsoHlf));
%     centR = topoMapM(topoMapM < longestD*pixDimIsoHlf);
%     centPsph = [centAz,centEl,centR];   % Surface points in spherical coordinate system
%     centEl = ((centEl-91)*pi/180)+cEl;
%     centAz = ((centAz-181)*pi/180)+cAz;
%     [centX,centY,centZ] = sph2cart(centAz,centEl,centR);
%     centPcart = round([centX,centY,centZ]+repmat(Origin,numel(centX),1)); % Surface points in crtesian coordinate sytem
%     
%     entropyMapM = nan(181,361);
%     intAvMapM = nan(181,361);
%     intStdMapM = nan(181,361);
%     for pNum = 1:size(centPcart,1)
%         Patch = im2double(mrImage(max(centPcart(pNum,1)-pRad,1):min(centPcart(pNum,1)+pRad,size(mrImage,1)),...
%                                   max(centPcart(pNum,2)-pRad,1):min(centPcart(pNum,2)+pRad,size(mrImage,2)),...
%                                   max(centPcart(pNum,3)-pRad,1):min(centPcart(pNum,3)+pRad,size(mrImage,3))));
%         entropyMapM(centPsph(pNum,2),centPsph(pNum,1)) = entropy(Patch);
%         intAvMapM(centPsph(pNum,2),centPsph(pNum,1)) = mean(Patch(:));
%         intStdMapM(centPsph(pNum,2),centPsph(pNum,1)) = std(Patch(:));
%     end
    
%% Displaying and saving:
%     [x,y,z] = ndgrid(1:size(pLabel, 1), 1:size(pLabel, 2), 1:size(pLabel, 3)); %get coordinates of all points
%     xx = x(pLabel == 1); %keep only coordinates for A == 1 
%     yy = y(pLabel == 1); %these 3 lines also reshape the 3d array
%     zz = z(pLabel == 1); %into column vectors
%     figure, subplot(1,2,1,'align'), plot3(xx,yy,zz,'.','markersize',1)
%     hold on;
%     plot3(Origin(1),Origin(2),Origin(3),'+')
%     plot3(pCenter(1),pCenter(2),pCenter(3),'O')
%     view([0 0 5]);
%     legend('Placenta center of mass','Point of observation')
%     title('Axial view of placenta label')
%     
%     subplot(1,2,2,'align'), plot3(xx,yy,zz,'.','markersize',1)
%     hold on;
%     plot3(Origin(1),Origin(2),Origin(3),'+')
%     plot3(pCenter(1),pCenter(2),pCenter(3),'O')
%     view([0 5 0]);
%     legend('Placenta center of mass','Point of observation')
%     title('Saggital view of placenta label')
%     saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Origin.png']))

    figure('Position', [1, 41, 1920, 1080]),
    imagesc(flip(topoMapF,1))
    ax = gca;
    ax.XTick = 1:180:361;
    ax.XTickLabel = num2cell(-180:180:180);
    ax.YTick = 1:180:181;
    ax.YTickLabel = num2cell(90:-180:-90);
    ax.FontSize = 14;
    colorbar
    colormap(flip(colormap(hot),1))
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    text(367,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Fetal side topography (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Topography_Fetal.png']))
    end

    figure('Position', [1, 41, 1920, 1080]),
    imagesc(flip(topoMapM,1))
    ax = gca;
    ax.XTick = 1:180:361;
    ax.XTickLabel = num2cell(-180:180:180);
    ax.YTick = 1:180:181;
    ax.YTickLabel = num2cell(90:-180:-90);
    ax.FontSize = 14;
    colormap(flip(colormap(hot),1))
    colorbar
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    text(367,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Maternal side topography (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Topography_Maternal.png']))
    end
    
%     figure('Position', [1, 41, 1920, 1080]),
%     imagesc(flip(topoMap360,1))
%     ax = gca;
%     ax.XTick = 1:180:361;
%     ax.XTickLabel = num2cell(-180:180:180);
%     ax.YTick = 1:180:181;
%     ax.YTickLabel = num2cell(90:-180:-90);
%     ax.FontSize = 14;
%     colorbar
%     colormap(flip(colormap(hot),1))
%     line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
%     line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
%     line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
%     line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
%     hold on
%     text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
%     text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
%     text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
%     text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
%     text(367,-5,'mm','Color','Black','FontSize',16)
%     xlabel(['Fetal side topography (Case ',imFiles(imNum).name(end-7:end-4), ')'])
%     saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Topography_360.png']))

    figure('Position', [1, 41, 1920, 1080]),
    imagesc(flip(thicknessMap,1))
    ax = gca;
    ax.XTick = 1:180:361;
    ax.XTickLabel = num2cell(-180:180:180);
    ax.YTick = 1:180:181;
    ax.YTickLabel = num2cell(90:-180:-90);
    ax.FontSize = 14;
    colorbar
    colormap hot
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    text(367,-5,'mm','Color','Black','FontSize',16)
    xlabel(['Placenta Thickness (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Thickness.png']))
    end

    figure
    imshow(flip(intensityMapF,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Fetal side intensity map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Intensity_Fetal.png']))
    end

    figure
    imshow(flip(intensityMapM,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Maternal side intensity map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Intensity_Maternal.png']))
    end

    figure
    imshow(flip(entropyMapF,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Fetal side entropy map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Entropy_Fetal.png']))
    end

    figure
    imshow(flip(entropyMapM,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Maternal side entropy map (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_Entropy_Maternal.png']))
    end

    figure
    imshow(flip(intAvMapF,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Fetal side average of intensities (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_LocalAverageIntensity_Fetal.png']))
    end

    figure
    imshow(flip(intAvMapM,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Maternal side average of intensities (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_LocalAverageIntensity_Maternal.png']))
    end

    figure
    imshow(flip(intStdMapF,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Fetal side side std of intensities (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_LocalStdIntensity_Fetal.png']))
    end

    figure
    imshow(flip(intStdMapM,1),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = gca;
    ax.XTick = 1:90:361;
    ax.XTickLabel = num2cell(-180:90:180);
    ax.YTick = 1:90:181;
    ax.YTickLabel = num2cell(90:-90:-90);
    ax.FontSize = 16;
    line([linePos(1) linePos(1)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(2) linePos(2)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(3) linePos(3)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    line([linePos(4) linePos(4)],[1 180],'Color','Blue','LineStyle','--','LineWidth',2)
    hold on
    text(linePos(1)-5,-5,'Left','Color','Blue','FontSize',16)
    text(linePos(2)-5,-5,'Ant.','Color','Blue','FontSize',16)
    text(linePos(3)-5,-5,'Right','Color','Blue','FontSize',16)
    text(linePos(4)-5,-5,'Post.','Color','Blue','FontSize',16)
    xlabel(['Maternal side std of intensities (Case ',imFiles(imNum).name(end-7:end-4), ')'])
    if toSave == true
        saveas(gcf, fullfile(saveFolder,[imFiles(imNum).name(1:end-4),'_LocalStdIntensity_Maternal.png']))
    end
    
    % Feature maps order: (1)topography fetal side, (2)topography maternal side,
    %                     (3)thickness, (4)intensity fetal side,
    %                     (5)intensity maternal side, (6)average intensity
    %                     fetal side, (7)average intensity maternal side,
    %                     (8)intensidt std fetal side, (9)intensity std
    %                     materal side, (10)entropy fetal side, (11)entropy
    %                     maternal side
    featureMaps_All = cat(4,featureMaps_All,cat(3,topoMapF,topoMapM,thicknessMap,...
                                                  intensityMapF,intensityMapM,...
                                                  intAvMapF,intAvMapM,...
                                                  intStdMapF,intStdMapM,...
                                                  entropyMapF,entropyMapM));
    
    pause(1)
    close all
end
fMapOrder = {'(1)topography fetal side', '(2)topography maternal side', '(3)thickness', '(4)intensity fetal side',...
                         '(5)intensity maternal side', '(6)average intensity fetal side', '(7)average intensity maternal side',...
                         '(8)intensidt std fetal side', '(9)intensity std materal side', '(10)entropy fetal side',...
                         '(11)entropy maternal side'};
if toSave == true
    save(fullfile(saveFolder,'FeatureMaps.mat'),'patientNums','featureMaps_All','patchSize_mm','fMapOrder')
end