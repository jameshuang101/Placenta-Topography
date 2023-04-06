clear
close all

dataFolder = 'W:\05_Topography\Data';
saveFolder = 'W:\05_Topography\Outputs';
mkdir(saveFolder)
toSave = true;
toReg = true;

imFiles = dir(fullfile(dataFolder,'Images','*.mat'));
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
    
    if toReg
        [mrImage, pLabel, uLabel] = simple_reg(mrImage, pLabel, uLabel, 100);
    end
    
    if toSave
        save(fullfile(saveFolder,imFiles(imNum).name,'Reg.mat'), mrImage);
        save(fullfile(saveFolder,sprintf('Label_0%s_placenta_Reg.mat',ptNum)), pLabel);
        save(fullfile(saveFolder,sprintf('Label_0%s_uterus_Reg.mat',ptNum)), uLabel);
    end
end