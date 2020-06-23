%Corepressive consortium minority wins analysis code
%Razan Alnahhas - Matthew Bennett Lab 2020

clearvars

for t=12    %trap number
    trap = sprintf('%02d',t);

    name = '8-13-19-minority-hallwayxy';
    fill = 22;  %frame where trap is filled with cells
    total = 133; %total number of frames

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
        cropc = cfp(y1:y2, x1:x2); %crop to only trap from above
        cropcmask = cropc.*mask16;  %multiply to cell mask
        
        
        if i<=21 %frames where CFP is on - change for each experiment/trap
             cyan = cropcmask>200;  %cfp threshold to be a cyan cell - empirically determined for all experimetns/traps
             cyan2 = bwareaopen(cyan,20);   %remove pixels smaller than a cell
             yellow = cyan2<1 & mask16==1;  %whatever is not cyan cell is a yellow cell (if is a cell on mask too)
             yellow2 = bwareaopen(yellow,20);   %remove pixels smaller than a cell
             y16 = uint16(yellow2); %convert T/F to 16-bit
             c16 = uint16(cyan2);   %convert T/F to 16-bit
        elseif i>21 & i<30 %frames where both on/dim - change for each experiment/trap
            yellow = cropymask>120; %yfp threshold to be a yellow cell - empirically determined for all experimetns/traps
            yellow2 = bwareaopen(yellow,20);    %remove pixels smaller than a cell
            cyan = yellow2<1 & mask16==1;   %whatever is not yellow cell is a cyan cell (if is a cell on mask too)
            cyan2 = bwareaopen(cyan,20);    %remove pixels smaller than a cell
            y16 = uint16(yellow2);  %convert T/F to 16-bit
            c16 = uint16(cyan2);    %convert T/F to 16-bit
         else    %frames where YFP is on - change for each experiment/trap
            yellow = cropymask>300; %yfp threshold to be a yellow cell - empirically determined for all experimetns/traps
            yellow2 = bwareaopen(yellow,20);    %remove pixels smaller than a cell
            cyan = yellow2<1 & mask16==1;   %whatever is not yellow cell is a cyan cell (if is a cell on mask too)
            cyan2 = bwareaopen(cyan,20);    %remove pixels smaller than a cell
            y16 = uint16(yellow2);  %convert T/F to 16-bit
            c16 = uint16(cyan2);    %convert T/F to 16-bit
         end
         
         %apply strain masks to FPs
         truey = cropy.*y16;
         truec = cropc.*c16;
        %remove zeroes from mask before averaging fluorescence
         truey = truey(truey>0);
         truec = truec(truec>0);
         

%         %uncomment below if want to save masks of ID as yfp or cfp to
%         %compare to images
%         imwrite(y16, strcat('yellow',trap,'t',num,'sw30.png'));
%         imwrite(c16, strcat('cyan',trap,'t',num,'sw30.png'));


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
