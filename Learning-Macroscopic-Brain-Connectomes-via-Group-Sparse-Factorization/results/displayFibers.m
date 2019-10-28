function displayFibers(fibers)
    numFascicles = length(fibers);
    colourIndex = 1;
    for i = 1:numFascicles
        fascicle = fibers{i};
        hold on
        hsvColour(:,1) = colourIndex/numFascicles;
        hsvColour(:,2) = 1;
        hsvColour(:,3) = 1;
        scatter3(fascicle(1,:), fascicle(2,:), fascicle(3,:), 10, hsv2rgb(hsvColour));
        colourIndex = colourIndex + 1;
    end
end