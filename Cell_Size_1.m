
%% ========================================================================
% CHLAMYDOMONAS CELL SIZE ANALYSIS - MANUAL CENTER AND EDGE
%
% LEICA LIF CHANNELS:
%   Channel 1 = GFP
%   Channel 2 = Grayscale
%   Channel 3 = Chlorophyll
%
% WORKFLOW:
% 1. Select Leica LIF file
% 2. Select series (all selected by default)
% 3. Preview chlorophyll Z-stack for each series
% 4. Select Z-slices
% 5. Overlay chlorophyll (magenta) on grayscale
% 6. Manually click:
%       First click  = CELL CENTER
%       Second click = OUTER CELL EDGE
% 7. Create a circular cell mask from these points
% 8. Review:
%       Accept circle
%       Adjust ellipse
%       Trace boundary manually
%       Retry center and edge
%       Skip
% 9. Measure major axis, minor axis, equivalent diameter
% 10. Select largest-area accepted slice per series
% 11. Generate grouped bar graph for all series
%
% Requires:
%   Bio-Formats MATLAB toolbox
%   Image Processing Toolbox
%
% NOTE:
% Measurements are in physical units from LIF calibration.
% Circle and ellipse are user-verified estimates of cell boundary.
%% ========================================================================

clear;
clc;
close all;

%% ====================== USER SETTINGS ============================

CH_GRAY = 2;
CH_CHL  = 3;

% Leave empty to read from LIF metadata
manualPixelX = [];
manualPixelY = [];

% Display
chlorophyllOverlayStrength = 0.55;
previewColumns = 6;
exportDPI = 300;

% Minimum acceptable cell mask area
minimumCellPixels = 30;

% Save image formats
savePNG = true;
saveFIG = true;

%% ====================== SELECT LIF FILE ==========================

[fileName,filePath] = uigetfile( ...
    '*.lif','Select Leica LIF file');

if isequal(fileName,0)
    error('No LIF file selected.');
end

if exist('bfGetReader','file') ~= 2
    error('Bio-Formats toolbox not found.');
end

fullFileName = fullfile(filePath,fileName);

reader = bfGetReader(fullFileName);

readerCleanup = onCleanup(@() reader.close());

numberSeries = reader.getSeriesCount();

%% ====================== SELECT SERIES ============================

seriesNames = cell(numberSeries,1);

for s = 1:numberSeries

    reader.setSeries(s-1);

    seriesNames{s} = sprintf( ...
        'Series %d | %d x %d | Z=%d | T=%d | C=%d', ...
        s,reader.getSizeX(),reader.getSizeY(), ...
        reader.getSizeZ(),reader.getSizeT(), ...
        reader.getSizeC());

end

[selectedSeries,ok] = listdlg( ...
    'PromptString','Select series to analyze:', ...
    'ListString',seriesNames, ...
    'SelectionMode','multiple', ...
    'InitialValue',1:numberSeries, ...
    'ListSize',[600 400]);

if ~ok || isempty(selectedSeries)
    error('No series selected.');
end

%% ====================== OUTPUT DIRECTORY =========================

[~,baseName] = fileparts(fileName);

outputDirectory = fullfile(filePath, ...
    [baseName '_Manual_CellSize']);

if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

allResults = table();
centralResults = table();

%% =================================================================
% MAIN LOOP: ALL SELECTED SERIES
%% =================================================================

for s = selectedSeries(:)'

    reader.setSeries(s-1);

    W = reader.getSizeX();
    H = reader.getSizeY();
    nZ = reader.getSizeZ();
    nT = reader.getSizeT();
    nC = reader.getSizeC();

    if nC < CH_CHL
        warning('Series %d has fewer than 3 channels.',s);
        continue;
    end

    fprintf('\n========================================\n');
    fprintf('SERIES %d\n',s);
    fprintf('========================================\n');

    %% ================= PIXEL CALIBRATION ==========================

    if ~isempty(manualPixelX) && ~isempty(manualPixelY)

        px = manualPixelX;
        py = manualPixelY;

    else

        try

            store = reader.getMetadataStore();

            umUnit = ome.units.UNITS.MICROMETER;

            px = double(store.getPixelsPhysicalSizeX( ...
                s-1).value(umUnit).doubleValue());

            py = double(store.getPixelsPhysicalSizeY( ...
                s-1).value(umUnit).doubleValue());

        catch

            answer = inputdlg( ...
                {'Pixel size X (um/pixel):', ...
                 'Pixel size Y (um/pixel):'}, ...
                sprintf('Series %d calibration',s), ...
                [1 45],{'',''});

            if isempty(answer)
                error('Calibration is required.');
            end

            px = str2double(answer{1});
            py = str2double(answer{2});

        end
    end

    assert(isfinite(px) && px>0 && ...
           isfinite(py) && py>0, ...
           'Invalid pixel calibration.');

    fprintf('Pixel X: %.6f um/pixel\n',px);
    fprintf('Pixel Y: %.6f um/pixel\n',py);

    %% ================= SELECT TIME POINT ==========================

    if nT > 1

        timeNames = arrayfun( ...
            @(k)sprintf('T=%d',k), ...
            1:nT,'UniformOutput',false);

        [selectedTime,ok] = listdlg( ...
            'PromptString','Select time point:', ...
            'ListString',timeNames, ...
            'SelectionMode','single');

        if ~ok
            continue;
        end

    else
        selectedTime = 1;
    end

    %% ================= LOAD Z STACKS ==============================

    grayStack = zeros(H,W,nZ,'single');
    chlStack  = zeros(H,W,nZ,'single');

    for z = 1:nZ

        idxGray = reader.getIndex( ...
            z-1,CH_GRAY-1,selectedTime-1)+1;

        idxChl = reader.getIndex( ...
            z-1,CH_CHL-1,selectedTime-1)+1;

        grayStack(:,:,z) = single( ...
            bfGetPlane(reader,idxGray));

        chlStack(:,:,z) = single( ...
            bfGetPlane(reader,idxChl));

    end

    %% ================= PREVIEW CHLOROPHYLL Z STACK =================

    stackMin = double(min(chlStack(:)));
    stackMax = double(max(chlStack(:)));

    if stackMax<=stackMin
        stackMax=stackMin+1;
    end

    hPreview = figure( ...
        'Color','w', ...
        'Name',sprintf('Series %d Chlorophyll Z-stack',s), ...
        'Position',[60 60 1300 850]);

    tl = tiledlayout( ...
        ceil(nZ/previewColumns),previewColumns, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for z = 1:nZ

        nexttile;

        C = (double(chlStack(:,:,z))-stackMin)/ ...
            (stackMax-stackMin);

        C = max(0,min(1,C));

        imshow(cat(3,C,zeros(H,W),C));

        title(sprintf('Z=%d',z));

    end

    title(tl,sprintf('Series %d: Chlorophyll',s));

    %% ================= SELECT Z SLICES ============================

    zNames = arrayfun( ...
        @(k)sprintf('Z=%d',k), ...
        1:nZ,'UniformOutput',false);

    [selectedZ,ok] = listdlg( ...
        'PromptString','Select Z slices for measurement:', ...
        'ListString',zNames, ...
        'SelectionMode','multiple', ...
        'InitialValue',1:nZ, ...
        'ListSize',[350 450]);

    if isgraphics(hPreview)
        close(hPreview);
    end

    if ~ok || isempty(selectedZ)
        continue;
    end

    %% ================= SERIES OUTPUT FOLDER =======================

    seriesDirectory = fullfile(outputDirectory, ...
        sprintf('Series_%03d_T%03d',s,selectedTime));

    if ~exist(seriesDirectory,'dir')
        mkdir(seriesDirectory);
    end

    seriesResults = table();

    %% =============================================================
    % PROCESS EACH SELECTED Z SLICE
    %% =============================================================

    for z = selectedZ(:)'

        fprintf('\nSeries %d | Z=%d\n',s,z);

        gray = normalizeMicroscopy(grayStack(:,:,z));
        chl  = normalizeMicroscopy(chlStack(:,:,z));

        RGB = overlayChlorophyll( ...
            gray,chl,chlorophyllOverlayStrength);

        accepted = false;
        skipped = false;
        currentMask = false(H,W);
        measurementMethod = "Manual center and edge";
        centerPixel = [NaN NaN];
        edgePixel = [NaN NaN];

        %% =========================================================
        % CLICK CENTER AND OUTER EDGE ON OVERLAY
        %% =========================================================

        while ~accepted && ~skipped

            hSelect = figure( ...
                'Color','w', ...
                'Name',sprintf('Series %d Z=%d Selection',s,z), ...
                'Position',[100 80 1200 650]);

            ax = axes(hSelect);

            imshow(RGB,'Parent',ax);
            hold(ax,'on');

            title(ax,{ ...
                sprintf('Series %d | Z=%d',s,z), ...
                'CLICK 1: CELL CENTER', ...
                'CLICK 2: OUTER CELL EDGE'});

            [xClick,yClick] = ginput(2);

            if isgraphics(hSelect)
                close(hSelect);
            end

            if numel(xClick)<2
                skipped=true;
                break;
            end

            centerPixel = [xClick(1),yClick(1)];
            edgePixel = [xClick(2),yClick(2)];

            %% ---- CALCULATE PHYSICAL RADIUS -----------------------

            deltaX_um = ...
                (edgePixel(1)-centerPixel(1))*px;

            deltaY_um = ...
                (edgePixel(2)-centerPixel(2))*py;

            radius_um = hypot(deltaX_um,deltaY_um);

            if radius_um<=0
                warning('Invalid radius. Please select again.');
                continue;
            end

            %% ---- INITIAL CIRCLE -------------------------------

            % The circle is physically circular.
            % If pixel dimensions differ, it appears
            % elliptical in pixel coordinates.

            radiusX_pixels = radius_um/px;
            radiusY_pixels = radius_um/py;

            ellipsePosition = [ ...
                centerPixel(1)-radiusX_pixels, ...
                centerPixel(2)-radiusY_pixels, ...
                2*radiusX_pixels, ...
                2*radiusY_pixels];

            currentMask = ellipseMaskFromPosition( ...
                ellipsePosition,0,H,W);

            measurementMethod = "Manual center-edge circle";

            %% =====================================================
            % REVIEW CIRCLE AND OPTIONALLY ADJUST ELLIPSE
            %% =====================================================

            reviewDone = false;
            retryPoints = false;

            while ~reviewDone && ~skipped

                hReview = figure( ...
                    'Color','w', ...
                    'Name',sprintf('Review S%d Z%d',s,z), ...
                    'Position',[80 80 1350 600]);

                tiledlayout(1,3, ...
                    'Padding','compact', ...
                    'TileSpacing','compact');

                %% OVERLAY

                nexttile;
                imshow(RGB);
                hold on;

                visboundaries(currentMask, ...
                    'Color','g','LineWidth',2);

                plot(centerPixel(1),centerPixel(2), ...
                    'yo','MarkerFaceColor','y', ...
                    'MarkerSize',6);

                plot(edgePixel(1),edgePixel(2), ...
                    'ro','MarkerFaceColor','r', ...
                    'MarkerSize',6);

                title('Overlay + selected cell boundary');

                %% GRAYSCALE

                nexttile;
                imshow(gray);
                hold on;

                visboundaries(currentMask, ...
                    'Color','g','LineWidth',2);

                title('Outer cell boundary on grayscale');

                %% MEASUREMENTS

                nexttile;
                imshow(gray);
                hold on;

                Mpreview = physicalMeasurements( ...
                    currentMask,px,py);

                drawPhysicalMeasurements( ...
                    Mpreview,px,py);

                title(sprintf( ...
                    'Major %.2f | Minor %.2f | EqD %.2f um', ...
                    Mpreview.majorUm, ...
                    Mpreview.minorUm, ...
                    Mpreview.equivUm));

                %% REVIEW MENU

                choice = menu( ...
                    'Does the proposed boundary fit the cell?', ...
                    'Accept', ...
                    'Adjust ellipse', ...
                    'Trace boundary manually', ...
                    'Retry center and edge', ...
                    'Skip');

                if isgraphics(hReview)
                    close(hReview);
                end

                switch choice

                    case 1
                        accepted = true;
                        reviewDone = true;

                    case 2

                        %% ADJUST ELLIPSE INTERACTIVELY

                        hAdjust = figure( ...
                            'Color','w', ...
                            'Name','Adjust Ellipse', ...
                            'Position',[150 90 950 750]);

                        imshow(RGB);
                        hold on;

                        title({ ...
                            'Drag, resize or rotate the ellipse', ...
                            'Double-click ellipse to finish'});

                        % Use moment-equivalent ellipse of
                        % current mask as editable initial ROI.
                        roiPosition = maskEllipseROI( ...
                            currentMask);

                        roi = drawellipse( ...
                            'Position',roiPosition(1:4), ...
                            'RotationAngle',roiPosition(5), ...
                            'Color','green', ...
                            'LineWidth',2);

                        wait(roi);

                        candidateMask = createMask(roi);

                        if isgraphics(hAdjust)
                            close(hAdjust);
                        end

                        if nnz(candidateMask)>=minimumCellPixels

                            currentMask = candidateMask;

                            measurementMethod = ...
                                "Manually adjusted ellipse";

                        else

                            warning('Ellipse mask too small.');

                        end

                    case 3

                        %% MANUAL FREEHAND OUTER BOUNDARY

                        hManual = figure( ...
                            'Color','w', ...
                            'Position',[100 80 1200 600]);

                        subplot(1,2,1);
                        imshow(RGB);
                        title('Overlay reference');

                        subplot(1,2,2);
                        imshow(gray);
                        hold on;

                        title({ ...
                            'Trace the complete OUTER cell boundary', ...
                            'Double-click to finish'});

                        roi = drawfreehand( ...
                            'Color','cyan', ...
                            'LineWidth',1.5);

                        wait(roi);

                        candidateMask = createMask(roi);

                        if isgraphics(hManual)
                            close(hManual);
                        end

                        if nnz(candidateMask)>=minimumCellPixels

                            candidateMask = ...
                                imfill(candidateMask,'holes');

                            currentMask = ...
                                bwareafilt(candidateMask,1);

                            measurementMethod = ...
                                "Manual freehand boundary";

                        else

                            warning('Freehand mask too small.');

                        end

                    case 4

                        retryPoints = true;
                        reviewDone = true;

                    otherwise

                        skipped = true;
                        reviewDone = true;

                end
            end

            if accepted || skipped
                break;
            end

            if retryPoints
                continue;
            end

        end

        %% =========================================================
        % SAVE ACCEPTED CELL
        %% =========================================================

        if skipped || ~accepted
            fprintf('Skipped Series %d Z=%d\n',s,z);
            continue;
        end

        if nnz(currentMask)<minimumCellPixels
            warning('Accepted mask too small.');
            continue;
        end

        %% CALCULATE PHYSICAL MEASUREMENTS

        M = physicalMeasurements(currentMask,px,py);

        touchesBorder = ...
            any(currentMask(1,:)) || ...
            any(currentMask(end,:)) || ...
            any(currentMask(:,1)) || ...
            any(currentMask(:,end));

        %% STORE RESULTS

        row = table( ...
            s,selectedTime,z, ...
            M.areaUm2, ...
            M.majorUm, ...
            M.minorUm, ...
            M.equivUm, ...
            M.aspectRatio, ...
            M.orientationDeg, ...
            centerPixel(1),centerPixel(2), ...
            edgePixel(1),edgePixel(2), ...
            touchesBorder,measurementMethod, ...
            'VariableNames',{ ...
            'Series','Time','Z', ...
            'Area_um2', ...
            'Major_um', ...
            'Minor_um', ...
            'EquivalentDiameter_um', ...
            'AspectRatio', ...
            'Orientation_deg', ...
            'SelectedCenterX_px','SelectedCenterY_px', ...
            'SelectedEdgeX_px','SelectedEdgeY_px', ...
            'TouchesBorder','Method'});

        seriesResults = [seriesResults;row];

        %% FINAL ANNOTATED IMAGE

        hFinal = figure( ...
            'Color','w', ...
            'Position',[150 100 950 750]);

        imshow(gray);
        hold on;

        % Accepted cell boundary
        visboundaries(currentMask, ...
            'Color','g','LineWidth',1.8);

        % Major axis, minor axis, ellipse and circle
        drawPhysicalMeasurements(M,px,py);

        title({ ...
            sprintf('Series %d | Z=%d',s,z), ...
            sprintf('Major=%.2f | Minor=%.2f | EqD=%.2f um', ...
            M.majorUm,M.minorUm,M.equivUm)});

        prefix = sprintf( ...
            'Series%03d_T%03d_Z%03d', ...
            s,selectedTime,z);

        if savePNG

            exportgraphics(hFinal, ...
                fullfile(seriesDirectory, ...
                [prefix '_Cell_Size.png']), ...
                'Resolution',exportDPI);

        end

        if saveFIG

            savefig(hFinal, ...
                fullfile(seriesDirectory, ...
                [prefix '_Cell_Size.fig']));

        end

        imwrite(uint8(currentMask)*255, ...
            fullfile(seriesDirectory, ...
            [prefix '_Mask.png']));

        save(fullfile(seriesDirectory, ...
            [prefix '_Mask.mat']), ...
            'currentMask','M','px','py', ...
            'centerPixel','edgePixel', ...
            'measurementMethod');

        close(hFinal);

        fprintf('Area = %.3f um^2\n',M.areaUm2);
        fprintf('Major = %.3f um\n',M.majorUm);
        fprintf('Minor = %.3f um\n',M.minorUm);
        fprintf('EqD = %.3f um\n',M.equivUm);

    end

    %% =============================================================
    % SELECT MAXIMUM-AREA ACCEPTED Z SLICE PER SERIES
    %% =============================================================

    if isempty(seriesResults)
        continue;
    end

    writetable(seriesResults, ...
        fullfile(seriesDirectory,'All_Z_Results.csv'));

    allResults = [allResults;seriesResults];

    validIndices = find(~seriesResults.TouchesBorder);

    if isempty(validIndices)

        warning('Series %d has no complete cell masks.',s);
        continue;

    end

    [~,idx] = max( ...
        seriesResults.Area_um2(validIndices));

    best = seriesResults(validIndices(idx),:);

    centralResults = [centralResults;best];

    writetable(best, ...
        fullfile(seriesDirectory, ...
        'Largest_Area_Central_Z.csv'));

    fprintf('\nSERIES %d: MAXIMUM AREA AT Z=%d\n', ...
        s,best.Z);

end

%% =================================================================
% ALL-SERIES SUMMARY
%% =================================================================

if ~isempty(allResults)

    writetable(allResults, ...
        fullfile(outputDirectory, ...
        'All_Series_All_Z.csv'));

end

if ~isempty(centralResults)

    centralResults = sortrows(centralResults,'Series');

    writetable(centralResults, ...
        fullfile(outputDirectory, ...
        'All_Series_Central_Cell_Size.csv'));

    %% GROUPED BAR GRAPH

    sizes = [ ...
        centralResults.Major_um, ...
        centralResults.Minor_um, ...
        centralResults.EquivalentDiameter_um];

    hBarFigure = figure( ...
        'Color','w', ...
        'Position',[150 150 1150 650]);

    hBar = bar(sizes,'grouped');

    hBar(1).FaceColor = [0.96 0.72 0.22];
    hBar(2).FaceColor = [0.18 0.70 0.88];
    hBar(3).FaceColor = [0.85 0.32 0.52];

    xticks(1:height(centralResults));

    xticklabels(compose( ...
        'Series %d',centralResults.Series));

    xlabel('LIF Series');
    ylabel('Cell size (\mum)');

    title('Cell Size: Maximum-Area Accepted Z-Slice');

    legend( ...
        'Major Axis', ...
        'Minor Axis', ...
        'Equivalent Diameter', ...
        'Location','best');

    grid on;
    box on;

    exportgraphics(hBarFigure, ...
        fullfile(outputDirectory, ...
        'Cell_Size_All_Series_BarGraph.png'), ...
        'Resolution',300);

    savefig(hBarFigure, ...
        fullfile(outputDirectory, ...
        'Cell_Size_All_Series_BarGraph.fig'));

    disp(centralResults);

end

fprintf('\nAnalysis complete.\n');
fprintf('Output directory:\n%s\n',outputDirectory);

%% =================================================================
% LOCAL FUNCTIONS
%% =================================================================

function J = normalizeMicroscopy(I)

    I = double(I);
    v = I(isfinite(I));

    if isempty(v)
        J = zeros(size(I));
        return;
    end

    lim = prctile(v,[1 99.5]);

    if lim(2)<=lim(1)
        J = zeros(size(I));
    else
        J = (I-lim(1))/(lim(2)-lim(1));
    end

    J = max(0,min(1,J));
    J(~isfinite(J)) = 0;

end

%% CHLOROPHYLL MAGENTA OVERLAY

function RGB = overlayChlorophyll(G,C,strength)

    alpha = strength*C;

    RGB = zeros([size(G),3]);

    RGB(:,:,1) = G.*(1-alpha)+alpha;
    RGB(:,:,2) = G.*(1-alpha);
    RGB(:,:,3) = G.*(1-alpha)+alpha;

end

%% INITIAL CIRCLE MASK

function B = ellipseMaskFromPosition(position,angle,H,W)

    cx = position(1)+position(3)/2;
    cy = position(2)+position(4)/2;

    a = position(3)/2;
    b = position(4)/2;

    theta = deg2rad(angle);

    t = linspace(0,2*pi,500);

    x = cx+a*cos(t)*cos(theta) ...
          - b*sin(t)*sin(theta);

    y = cy+a*cos(t)*sin(theta) ...
          + b*sin(t)*cos(theta);

    B = poly2mask(x,y,H,W);
    B = imfill(B,'holes');

end

%% EXTRACT AN EDITABLE ELLIPSE INITIALIZATION FROM MASK

function roiPosition = maskEllipseROI(B)

    p = regionprops(B, ...
        'Centroid','MajorAxisLength', ...
        'MinorAxisLength','Orientation');

    if isempty(p)
        error('Cannot initialize ellipse from empty mask.');
    end

    p = p(1);

    cx = p.Centroid(1);
    cy = p.Centroid(2);

    width = p.MajorAxisLength;
    height = p.MinorAxisLength;

    angle = -p.Orientation;

    roiPosition = [ ...
        cx-width/2, ...
        cy-height/2, ...
        width,height,angle];

end

%% PHYSICAL CELL MEASUREMENTS
% Exact physical pixel spacing is used in second moments.

function M = physicalMeasurements(B,px,py)

    [y,x] = find(B);

    if isempty(x)
        error('Empty cell mask.');
    end

    xy = [ ...
        (double(x)-0.5)*px, ...
        (double(y)-0.5)*py];

    center = mean(xy,1);

    X = xy-center;

    covariance = (X'*X)/size(X,1);

    [V,D] = eig(covariance);

    [lambda,idx] = sort(diag(D),'descend');

    V = V(:,idx);
    lambda = max(lambda,0);

    M.majorUm = 4*sqrt(lambda(1));
    M.minorUm = 4*sqrt(lambda(2));

    M.areaUm2 = nnz(B)*px*py;

    M.equivUm = sqrt(4*M.areaUm2/pi);

    M.aspectRatio = ...
        M.majorUm/max(M.minorUm,eps);

    M.orientationDeg = ...
        atan2d(V(2,1),V(1,1));

    M.centroidUm = center;
    M.axisVectors = V;

end

%% DRAW FITTED ELLIPSE, MAJOR AXIS, MINOR AXIS AND CIRCLE

function drawPhysicalMeasurements(M,px,py)

    c = M.centroidUm(:);
    V = M.axisVectors;

    a = M.majorUm/2;
    b = M.minorUm/2;

    t = linspace(0,2*pi,400);

    %% FITTED ELLIPSE

    E = c+V*[a*cos(t);b*sin(t)];

    plot(E(1,:)/px+0.5, ...
         E(2,:)/py+0.5, ...
         'y--','LineWidth',1.5);

    %% MAJOR AXIS

    major = [ ...
        c-a*V(:,1), ...
        c+a*V(:,1)];

    plot(major(1,:)/px+0.5, ...
         major(2,:)/py+0.5, ...
         'c-','LineWidth',2.5);

    %% MINOR AXIS

    minor = [ ...
        c-b*V(:,2), ...
        c+b*V(:,2)];

    plot(minor(1,:)/px+0.5, ...
         minor(2,:)/py+0.5, ...
         '-','Color',[1 0.65 0], ...
         'LineWidth',2.5);

    %% EQUIVALENT DIAMETER CIRCLE

    radius = M.equivUm/2;

    xc = (c(1)+radius*cos(t))/px+0.5;
    yc = (c(2)+radius*sin(t))/py+0.5;

    plot(xc,yc,'r-','LineWidth',1.5);

    %% CENTROID

    plot(c(1)/px+0.5, ...
         c(2)/py+0.5, ...
         'wo','MarkerFaceColor','w', ...
         'MarkerSize',5);

end
