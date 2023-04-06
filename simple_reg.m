function [mrImage, pLabel, uLabel] = simple_reg(mrImage, pLabel, uLabel, itr)
    
    oddIm = mrImage(:,:,1:2:end); % Odd slices of the MR image
    evenIm = mrImage(:,:,2:2:end); % Even slices of the MR image
    oddU = uLabel(:,:,1:2:end); 
    evenU = uLabel(:,:,2:2:end);
    oddP = pLabel(:,:,1:2:end);
    evenP = pLabel(:,:,2:2:end);
    sharpRad = 1.5;
    sharpAmt = 1.5;
    sharpThrs = 0.7;
    
    % Determine if odd or even sequence is at the time of greatest
    % inhalation
    [originX, originY, originZ] = ind2sub(size(oddU),find(oddU));
    oddCenter = mean([originX, originY, originZ]);  %Point of Observation (uterus center of mass)
    [originX, originY, originZ] = ind2sub(size(evenU),find(evenU));
    evenCenter = mean([originX, originY, originZ]);  %Point of Observation (uterus center of mass)
    clear originX originY originZ
    
    % If odd slices have greatest inhalation then register odd slices, else
    % register even slices
    if oddCenter(2) < evenCenter(2)
        % Make it so that there is 1 less even slice than odd slices
        if size(evenIm,3) == size(oddIm,3)
            evenIm = evenIm(:,:,1:end-1);
%             evenU = evenU(:,:,1:end-1);
%             evenP = evenP(:,:,1:end-1);
        end
        
        % Interpolate between odd slices to create target even slices
        [xq,yq,zq] = meshgrid(1:size(oddIm,2), 1:size(oddIm,1), 1:0.5:size(oddIm,3));
        evenTarget = interp3(double(oddIm), xq,yq,zq, 'spline');
        evenTarget = evenTarget(:,:,2:2:end);
        for n = 1:size(evenTarget,3)
            evenTarget(:,:,n) = imsharpen(evenTarget(:,:,n),'Radius',sharpRad,'Amount',sharpAmt,'Threshold',sharpThrs);
        end
        oddP_Int = level_set_Interp3(double(oddP),xq,yq,zq,'cubic');
        oddU_Int = level_set_Interp3(double(oddU),xq,yq,zq,'cubic');
%         pLabel(:,:,1:size(oddP_Int,3)) = oddP_Int;
%         uLabel(:,:,1:size(oddU_Int,3)) = oddU_Int;
        for s = 1:size(oddP_Int,3)
            if sum(uint8(oddP_Int(:,:,s)),'all') > 0
                pLabel(:,:,s) = oddP_Int(:,:,s);
            end
            if sum(uint8(oddU_Int(:,:,s)),'all') > 0
                uLabel(:,:,s) = oddU_Int(:,:,s);
            end
        end
        
        % Register even slices to interpolated target even slices
        evenIm = double(evenIm);
        [D, evenImReg] = imregdemons(evenIm, evenTarget, itr);
        evenpLabel = imwarp(double(evenP), D);
        evenuLabel = imwarp(double(evenU), D);
        for i = 2:2:2*size(evenIm,3)
           mrImage(:,:,i) = evenImReg(:,:,i/2);
           pLabel(:,:,i) = evenpLabel(:,:,i/2);
           uLabel(:,:,i) = evenuLabel(:,:,i/2);
        end
        
        
    else 
        % Make it so that there is 1 less odd slice than even slices
        if size(oddIm,3) == size(evenIm,3)
            oddIm = oddIm(:,:,2:end);
%             oddU = oddU(:,:,2:end);
%             oddP = oddP(:,:,2:end);
        elseif size(oddIm,3) > size(evenIm,3)
            oddIm = oddIm(:,:,2:end-1);
%             oddU = oddU(:,:,2:end-1);
%             oddP = oddP(:,:,2:end-1);
        end
        
        % Interpolate between even slices to create target odd slices
        [xq,yq,zq] = meshgrid(1:size(evenIm,2), 1:size(evenIm,1), 1:0.5:size(evenIm,3));
        oddTarget = interp3(double(evenIm), xq,yq,zq, 'spline');
        oddTarget = oddTarget(:,:,2:2:end);
        for n = 1:size(oddTarget,3)
            oddTarget(:,:,n) = imsharpen(oddTarget(:,:,n),'Radius',sharpRad,'Amount',sharpAmt,'Threshold',sharpThrs);
        end
        evenP_Int = level_set_Interp3(double(evenP),xq,yq,zq,'cubic');
        evenU_Int = level_set_Interp3(double(evenU),xq,yq,zq,'cubic');
        pLabel(:,:,2:size(evenP_Int,3)+1) = evenP_Int;
        uLabel(:,:,2:size(evenU_Int,3)+1) = evenU_Int;
        for s = 1:size(evenP_Int,3)
            if sum(uint8(evenP_Int(:,:,s)),'all') > 0
                pLabel(:,:,s) = evenP_Int(:,:,s);
            end
            if sum(uint8(evenU_Int(:,:,s)),'all') > 0
                uLabel(:,:,s) = evenU_Int(:,:,s);
            end
        end
        
        % Register odd slices to interpolated target odd slices
        oddIm = double(oddIm);
        [D, oddImReg] = imregdemons(oddIm, oddTarget, itr);
%         oddpLabel = imbinarize(imwarp(double(oddP), D), 0.5);
%         odduLabel = imbinarize(imwarp(double(oddU), D), 0.5);
        for i = 3:2:size(mrImage,3)-1
           mrImage(:,:,i) = oddImReg(:,:,floor(i/2));
%            pLabel(:,:,i) = oddpLabel(:,:,floor(i/2));
%            uLabel(:,:,i) = odduLabel(:,:,floor(i/2));
        end
    end

end