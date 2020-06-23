%Corepressive consortium majority wins analysis code
%Razan Alnahhas - Matthew Bennett Lab 2020

clearvars

for t=10    %trap number
    trap = sprintf('%02d',t);

    name = '08-20-18-toggle-dimeric-hallway-only-xylxy';
    fill = 13;  %frame where trap is filled with cells
    total = 181; %total number of frames

    %select only trap to crop images
    img = imread(strcat(name,trap,'c1t001.tif'));
    figure(1);
    imshow(img, [0 4095]);
    title('Select Trap')
    rect = getrect;
    close

    %crop dimensions
    x1 = round(rect(1,1));
    x2 = round(rect(1,1) + rect(1,3));
    y1 = round(rect(1,2));
    y2 = round(rect(1,2) + rect(1,4));

    %select region in trap without cells for thresholds
    figure(2);
    imshow(img, [0 4095]);
    title('Select empty region')
    rect = getrect;
    close

    xmin = round(rect(1,1));
    xmax = round(rect(1,1) + rect(1,3));
    ymin = round(rect(1,2));
    ymax = round(rect(1,2) + rect(1,4));
    
    %calculate background 'darkness' from empty region selected
    backgroundphase = img(ymin:ymax,xmin:xmax);
    meanphase = mean2(backgroundphase);
    stdphase = std2(backgroundphase);

    %threshold darkness for cells
    thresholdPHASE = meanphase-stdphase
    
    %data save value (trap fill frame becomes frame 1)
    ii=1;
    
    
    for i=fill:total

        data(ii,1) = 6*(ii-1); %imaged ever 6 minutes, trap fill image is t=0 min
        num = sprintf('%03d',i);
        
        %create mask of cells
        phase = imread(strcat(name,trap,'c1t',num,'.tif')); %open phase image
        cropph = phase(y1:y2, x1:x2);   %crop to only trap from above
        mask = cropph<thresholdPHASE;   %threshold for darkness of cells
        mask16 = uint16(mask);  %convert T/F to 16-bit to matrix multiply below
        
        yfp = imread(strcat(name,trap,'c2t',num,'.tif'));   %open yfp image
        cropy = yfp(y1:y2, x1:x2);  %crop to only trap from above
        cropymask = cropy.*mask16;  %multiply to cell mask

        cfp = imread(strcat(name,trap,'c3t',num,'.tif'));   %open cfp image
        cropc = cfp(y1:y2, x1:x2);  %crop to only trap from above
        cropcmask = cropc.*mask16;  %multiply to cell mask
        
        ratio = cropymask./cropcmask;   %calculate ratio of yfp:cfp

        yellow = ratio>0.5; %if y:c>0.5 the cell is yellow
        cyan = ratio<=0.5 & mask16==1;  %if y:c<=0.5 the cell is cyan (if also cell on mask)
        y16 = uint16(yellow);   %convert T/F to 16-bit
        c16 = uint16(cyan); %convert T/F to 16-bit
         
         %apply strain masks to FPs
         truey = cropy.*y16;
         truec = cropc.*c16;
         
        %remove zeroes from mask before averaging fluorescence
         truey = truey(truey>0);
         truec = truec(truec>0);
        
%         %uncomment below if want to save masks of ID as yfp or cfp to
%         %compare to images
%         imwrite(y16, strcat('yellow',trap,'t',num,'var.png'));
%         imwrite(c16, strcat('cyan',trap,'t',num,'var.png'));

        data(ii,2) = sum(sum(yellow2)); %save total number yellow pixels over time
        data(ii,3) = sum(sum(cyan2)); %save total number cyan pixels over time
        data(ii,4) = mean2(truey); %save average yellow fluorescence over time
        data(ii,5) = mean2(truec); %save average cyan fluorescence over time
        data(ii,6) =  sum(sum(yellow2))+sum(sum(cyan2)); %save total number pixels over time

        disp(i);    %dsiplay current frame while code is running to monitor progress
        
        ii=ii+1;    %increase frame count

    end

    %save data file for trap
    file = strcat(name,trap,'data.mat');
    save(file, 'data');

    %calculate strain ratios
    yratio(:) = data(:,2)./data(:,6);
    cratio(:) = data(:,3)./data(:,6);

    %graph fluor & ratio over time
    yyaxis left
    plot(data(:,1),yratio(:),'g');
    hold on
    plot(data(:,1),cratio(:),'b');
    xlabel('Time(min)', 'FontSize', 12);
    ylabel('Strain Ratio', 'FontSize', 12);
    yyaxis right
    plot(data(:,1),data(:,4),'.','MarkerEdgeColor','g');
    plot(data(:,1),data(:,5),'.','MarkerEdgeColor','b');
    ylabel('Average Fluorescence', 'FontSize', 12);
    title(strcat('Trap ', trap));
    hold off

end
