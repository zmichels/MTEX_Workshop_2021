%% plot/save/project string
pname = 'workshop';


%%
%**************************************************************************
%%                   Field data in MTEX
%**************************************************************************

% Fake geographic orientations of foliation plane and lineation
strikeRHR = [110, 200, 330, 20, 33, 45, 35, 36, 29, 27, 40];
dip = [70, 80, 15, 40, 45, 60, 47, 50, 39, 33, 52];
trend = [15, 1, 11, 105, 220, 230, 36, 101, 30, 103, 227];

% Orientations of downward foliation pole and lineation in your EBSD map.
mapPole = -yvector;
mapLin = xvector;

% use the Fabrica function SDT2or to define a fabric orientation
[fabOr, fabRot, strikeV, poleV, linV] = SDT2or(strikeRHR,dip,trend);


% Plot the foliation, lineation 
figure
plot(poleV,'plane','LineColor','k','LineWidth',1,'antipodal','lower')
hold on
% add lineation (should plot in lower hemisphere on the foliation plane)
plot(linV,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower','DisplayName','lineation')
% pole to foliation (should plot in lower hemisphere if working properly
plot(poleV,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower','DisplayName','foliation')

% add a legend to show the DisplayName/s designated above
legend show

% save the figure
saveFigure(sprintf('%s_geo_field_example_Lin_&_Fol.png',pname),'-bestfit')



%% Figures to compare directional analysis with orientation analysis

% compute the mean folitation pole direction and lineation direction
% separately
mFol = mean(poleV);
mLin = mean(linV);

figure,
plot(mFol,'plane','LineColor','k','LineWidth',1,'antipodal','lower')
hold on
% add lineation (should plot in lower hemisphere on the foliation plane)
plot(mLin,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower', 'DisplayName','mean lineation')
% pole to foliation (should plot in lower hemisphere if working properly
plot(mFol,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower','DisplayName','mean foliation')

% compute the the closest vector on foliation plane to lineation
onPlaneV = cross(cross(mFol,mLin),mFol);

% compute te angle between it and mean lineation
angMeanLinFolSep = angle(mLin,onPlaneV,'antipodal')./degree;

% add this point to the plot and show the arc between them
hold on
% add the point
annotate(onPlaneV,'MarkerSize',12,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k')
% draw shortest line/arc between them, give DisplayName the angle
line(mLin,onPlaneV,'LineColor','g','DisplayName', sprintf('~%i°',round(angMeanLinFolSep)))

legend show


% save the figure
saveFigure(sprintf('%s_geo_field_example_meanSep_Lin_&_Fol.png',pname),'-bestfit')

%% Now compare with the fabric **orientations** (computed above)

% NOTE:
% in the function used above, the fabric orientation 'fabOr' is configured
% such that the Miller direction [100] of the orientation matches the
% direction of the lineation, and the foliation pole matches the Miller
% direction [010].
% so...

% lineations
linDir = Miller(1,0,0,'direction',fabOr.CS);
lins = fabOr.*linDir;
% foliations
folDir = Miller(0,1,0,'direction',fabOr.CS);
fols = fabOr*Miller(0,1,0,'direction',fabOr.CS);

% for plotting
hFab = [Miller(0,1,0,'direction',fabOr.CS),Miller(1,0,0,'direction',fabOr.CS)];


% so we could make a 'pole figure' of the foliation and lineation
figure,
plotPDF(fabOr,hFab,'antipodal','lower','smooth','halfwidth',20*degree,'colorrange','equal')
% add a color bar
cb = mtexColorbar('Title','M.U.D.');
% add the field data
hold on
plot(poleV,'plane','LineColor','k','LineWidth',1,'antipodal','lower','add2all')
hold on
% add lineation (should plot in lower hemisphere on the foliation plane)
plot(linV,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower','DisplayName','lineation','add2all')
% pole to foliation (should plot in lower hemisphere if working properly
plot(poleV,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower','DisplayName','foliation','add2all')

% modify the "Miller" labels to reflect the plot
f = gcf;
f.Children(end).Title.String = ['\bf ', 'foliation', ' \rm'];
f.Children(end-1).Title.String = ['\bf ', 'lineation', ' \rm'];

saveFigure(sprintf('%s_field_example_fabPDF_data.png',pname),'-bestfit')

%% lets compute the mean of fabOr 
% compute the mean
mFab = mean(fabOr);
% lets compute the mode of fabOr by first computing an ODF
% ODF
fabODF = calcDensity(fabOr,'halfwidth',10*degree);
% mode
fabMode = calcModes(fabODF,1);

% plot together
figure,
plotPDF(fabOr,hFab,'antipodal','lower','smooth','halfwidth',20*degree,'colorrange','equal')
% add a color bar
cb = mtexColorbar('Title','M.U.D.');
% add the mean
hold on
plotPDF(mFab,hFab,'Marker','s','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','mean','MarkerFaceAlpha',.75)
% add the mode
hold on
plotPDF(fabMode,hFab,'Marker','s','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','mode','MarkerFaceAlpha',.75)


% modify the "Miller" labels to reflect the plot
f = gcf;
f.Children(end).Title.String = ['\bf ', 'foliation', ' \rm'];
f.Children(end-1).Title.String = ['\bf ', 'lineation', ' \rm'];


saveFigure(sprintf('%s_field_example_fabPDF_Mean_Mode.png',pname),'-bestfit')

%% NOTE:
        % The angle between the foliation and lineation is always
        % 90-degrees based on the way the fabric orientation was defined –
        % including means and modes.
        % That is not always the case when analyzing the directions of the
        % fabric elements separately (as demonstrated above)
        
        
%% A comparison: Which to use for "preferred orientation", mode or mean?       
%% Plot each (fabric mean and mode) with the fabric data
%% for the fabric mean
% Plot the foliation, lineation 
figure
plot(fols,'plane','LineColor','k','LineWidth',1,'antipodal','lower')
hold on
% add lineation (should plot in lower hemisphere on the foliation plane)
plot(lins,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower')
% pole to foliation (should plot in lower hemisphere if working properly
plot(fols,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10)

hold on
% add the 'fabric' mean (mFab) of each fabric-element direction
% great circle of foliation plane
plot(mFab*folDir,'plane','LineColor','m','LineWidth',4)
% foliation
plot(mFab*folDir,'Marker','s','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','foliation')
% lineation
plot(mFab*linDir,'Marker','o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','lineation')

% a legend
leg = legend('show');
title(leg,'fabric mean')
saveFigure(sprintf('%s_field_example_fabPDF_Mean__and_data.png',pname),'-bestfit')


%% for the fabric mode
% Plot the foliation, lineation 
figure
plot(fols,'plane','LineColor','k','LineWidth',1,'antipodal','lower')
hold on
% add lineation (should plot in lower hemisphere on the foliation plane)
plot(lins,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower')
% pole to foliation (should plot in lower hemisphere if working properly
plot(fols,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10)

hold on
% add the 'fabric' mode (fabMode) of each fabric-element direction
% great circle of foliation plane
plot(fabMode*folDir,'plane','LineColor','c','LineWidth',4)
% foliation
plot(fabMode*folDir,'Marker','s','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','foliation')
% lineation
plot(fabMode*linDir,'Marker','o','MarkerFaceColor','c','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','lineation')

% a legend
leg = legend('show');
title(leg,'fabric mode')
saveFigure(sprintf('%s_field_example_fabPDF_Mode__and_data.png',pname),'-bestfit')


%%
hold on
% add the 'fabric' mean (mFab) of each fabric-element direction
% great circle of foliation plane
plot(mFab*folDir,'plane','LineColor','m','LineWidth',4)
% foliation
plot(mFab*folDir,'Marker','s','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','foliation')
% lineation
plot(mFab*linDir,'Marker','o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',15,'DisplayName','lineation')

% 
saveFigure(sprintf('%s_field_example_fabPDF_Mean_&_Mode.png',pname),'-bestfit')


%% Directional analysis of individual fabric elements is also possible 
% Example: 

%% First foliation poles separately
% Generically, we can define a set of vectors 'v' to analze...
v1 = fols; % ... in this case, foliation poles.
% setup a spherical grid for analysis
r = plotS2Grid('resolution',1*degree,'antipodal');
% compute the density of axes
kde1 = calcDensity([v1 -v1],r,'antipodal','halfwidth',10*degree);
% find the index of the point on the grid with the maximum density
[~,I1] = max(kde1);
% assign the vector from the "max" point to a new variable
vmax1 = r(I1);
% force to lower hemisphere
vmax1(vmax1.z>0) = -vmax1(vmax1.z>0);

figure,
plot(r,kde1,'antipodal','lower')
hold on
plot(vmax1,'Marker','s','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5,'MarkerSize',15,'DisplayName','foliation')
plot(vmax1,'plane')
saveFigure(sprintf('%s_field_example_folPoles_Maxdensity.png',pname),'-bestfit')


%% Second lineations separately
v2 = lins; % ... in this case, foliation poles.
% compute the density of lineations (use same grid as above)
kde2 = calcDensity([v2 -v2],r,'antipodal','halfwidth',10*degree);
[~,I2] = max(kde2);
vmax2 = r(I2);
vmax2(vmax2.z>0) = -vmax2(vmax2.z>0);

figure,
plot(r,kde2,'antipodal','lower')
hold on
plot(vmax2,'Marker','o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5,'MarkerSize',15,'DisplayName','foliation')
saveFigure(sprintf('%s_field_example_lineations_Maxdensity.png',pname),'-bestfit')

%% check the angle between "max" lineation and foliation pole (is it 90°?)
angleCheck = angle(vmax1,vmax2,'antipodal')./degree
% probably it's not *exactly* 90 degrees as it is supposed to be... :-/


%% Plot the results of separate directional analyses
figure,
% contour the foliation poles in black
plot(fols,'antipodal','lower','contour','LineColor','k','halfwidth',10*degree)
hold on
% add contours of lineation in red
plot(lins,'antipodal','lower','contour','LineColor','r','halfwidth',10*degree)
% add the "max" directions determined from the KDE analysis above
plot(vmax1,'Marker','s','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',15,'LineWidth',3,'DisplayName','foliation')
plot(vmax2,'Marker','o','MarkerFaceColor','w','MarkerEdgeColor','r','MarkerSize',15,'LineWidth',3,'DisplayName','lineation')
plot(vmax1,'antipodal','lower','plane')
hold on
annotation('textbox',[0 .7 0 .3],'String',sprintf('angle = %i°',floor(angleCheck)),'FitBoxToText','on','FontSize',12);
saveFigure(sprintf('%s_field_example_sep_dir_analysis.png',pname),'-bestfit')



%% Close figures
close all








%%
%**************************************************************************
%%             E B S D     D A T A   A N D   A N A L Y S I S
%*************************************************************************

%**************************************************************************
%%              I M P O R T     A N D     R O T A T E
%**************************************************************************

%% Import some data (we are going to use olivine as an example)
mtexdata forsterite
% plot/save/project string
pname = 'workshop';

% phase of interest today
phase = 'forsterite';

% default plotting of the coordinate axes
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

setMTEXpref('bAxisDirection','east');
setMTEXpref('aAxisDirection',''); % undefined

% some setup for plotting 
% orientations of 'Forsterite'
o = ebsd(phase).orientations;
% Crystal Symmerty (needed for defining relevant Miller indices)
cs = o.CS;
% Miller indicies of interest
h = [Miller(1,0,0,'direction',cs),Miller(0,1,0,'direction',cs),Miller(0,0,1,'direction',cs)];


%% Plot the band-contrast map and a pole figure (check 1)
% Confirm a match with with acquisition software
% If not a match, apply coorect rotations to your needs
figure,
plot(ebsd,ebsd.bc,'figsize','large')
mtexColorMap(gray)
saveFigure(sprintf('%s_import_map_check1.png',pname),'-bestfit')

figure,
plotPDF(o,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal','figsize','medium')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
saveFigure(sprintf('%s_import_pf_check1.png',pname),'-bestfit')

%% Hmmm... something looks off compared to acquisition software (example)
close all

%% You determine you need the following rotations to match acquisition...

% NOTE: Also using 'convertEuler2SpatialReferenceFrame' with CTFs
ebsd = rotate(ebsd,rotation('axis',yvector,'angle',180*degree));
ebsd = rotate(ebsd,rotation('axis',zvector,'angle',180*degree),'keepEuler');
o = ebsd(phase).orientations;

%% Check again (check 2)
figure,
plot(ebsd,ebsd.bc,'figsize','large')
mtexColorMap(gray)
saveFigure(sprintf('%s_import_map_check2_corr.png',pname),'-bestfit')

figure,
plotPDF(o,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal','figsize','medium')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
saveFigure(sprintf('%s_import_pf_check2_corr.png',pname),'-bestfit')

%% Now... let's say we also want to rotate it 90° for our fabric reference
close all
% rotate 90°
ebsd = rotate(ebsd,rotation('axis',zvector,'angle',90*degree));
o = ebsd(phase).orientations;

%% Last check (check 3)
figure,
plot(ebsd,ebsd.bc,'figsize','large')
mtexColorMap(gray)
saveFigure(sprintf('%s_import_map_check3_ROT_corr.png',pname),'-bestfit')

figure,
plotPDF(o,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal','figsize','medium')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
saveFigure(sprintf('%s_import_pf_check3_ROT_corr.png',pname),'-bestfit')

%% Hmmm... maybe we want to crop the map
close all

figure,
plot(ebsd,ebsd.bc,'figsize','large')
mtexColorMap(gray)
hold on

% select opposite corners of a rectangle on the plot
% p2 = ginput(2);
% Get XY pairs for all 4 corners
% p = [p2(1,:); diag(flipud(p2))'; p2(2,:); diag(p2)'];

% for wokshop we'll use same points I already clicked...
p = 1.0e+04  *  [0.0081    3.2957
                 1.6599    3.2957
                 1.6599    0.6134
                 0.0081    0.6134];

%% Isolate the EBSD data
ebsd = ebsd(inpolygon(ebsd,p));
o = ebsd(phase).orientations;

%% Plot the working dataset
figure,
plot(ebsd,ebsd.bc,'figsize','large')
mtexColorMap(gray)
saveFigure(sprintf('%s_BC_map.png',pname),'-bestfit')

figure,
plotPDF(o,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal','figsize','medium')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
hold on
plot(yvector,'plane','add2all','LineWidth',3)
saveFigure(sprintf('%s_%s_CPO.png',pname,o.CS.mineral),'-bestfit')



%**************************************************************************
%% C O M P U T E    G R A I N   B O U N D A R I E S
%**************************************************************************

%% Initial grain boundary computation


% in this example, we are going to use only indexed points
% also we will specify a second angle threshold to use for identifying
% intragranular boundaries and their (subgrain boundaries).
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'boundary','tight','angle',[10 1]*degree);

% remove [potentially spurious] small grains (few solutions)
grains = grains(grains.grainSize>9);

% filter the ebsd data by the remaining grainset
ebsd = ebsd(grains);

% recompute grain boundaries
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'boundary','tight','angle',[10 1]*degree);

% smooth grain boundaries for sake of later shape analysis
grains = smooth(grains,5,'moveTriplePoints');

%% Plot the grains
figure,
plot(grains,'figsize','large')
saveFigure(sprintf('%s_grain_reconstruction.png',pname),'-bestfit')


%% Make a histogram of some grain attribute
% let's chose something shape-related like aspect ratio...
figure,
histogram(grains('indexed').aspectRatio)
xlabel('aspect ratio')
ylabel('number of grains')
saveFigure(sprintf('%s_all_grains_aspectRatio.png',pname),'-bestfit')



%% close figs
close all

%**************************************************************************
%%               G R A I N       S .P .O.
%**************************************************************************

% first filter the grains to remove the border grains

% ids of the outer boundary segment
outerBoundary_id = any(grains.boundary.grainId==0,2);

% corresponding grain_id
grain_id = grains.boundary(outerBoundary_id).grainId;

% remove all zeros
grain_id(grain_id==0) = [];


%% fit ellipses to grain hulls
[grains.prop.omega,grains.prop.a,grains.prop.b] = fitEllipse(grains);


% Make vectors of grain-ellipse long axes
grains.prop.spov = rotate(xvector,rotation('axis',zvector,'angle',grains.omega));


% % SPO grains filtered by aspect ratio > 1.4
gSPO = grains(grains.a./grains.b>1.4);
% gSPO = g;


% plot polar histograms (rose diagrams) of grain long axes 
nbins = 360/5;

% plot
figure,
polarhistogram([gSPO(phase).omega; gSPO(phase).omega+pi],nbins,'Normalization','pdf','FaceColor',gSPO(phase).color,'FaceAlpha',.75,'EdgeColor','none');
hold on
prose = polarhistogram([gSPO(phase).omega; gSPO(phase).omega+pi],nbins,'Normalization','pdf');
hold on
prose.DisplayStyle = 'stairs';
prose.EdgeColor = 'k';
drawnow
mineral = grains(phase).mineral;
annotation('textbox',[0 .7 0 .3],'String',sprintf('n = %i',length(gSPO(phase))),'FitBoxToText','on','FontSize',12);

saveFigure(sprintf('%s_%s_longAxisRose.png',pname,mineral),'-bestfit')


%% Invert long and short axes in crystal reference frame and plot them

% inverse orientation
io = inv(gSPO(phase).meanOrientation);

% spo
spo = gSPO(phase).spov;

% long-axes
hl = io.*spo;

% short-axes
hs = io.*rotate(spo,rotation('axis',zvector,'angle',90*degree));

warning off

% plot long-axis in crystal reference
figure,
plot(hl.symmetrise,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal','fundamentalRegion','figSize','medium')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
mtexTitle(sprintf('%s grain long-axes',mineral))
saveFigure(sprintf('%s_%s_longAxisXTAL.png',pname,mineral),'-bestfit')


% plots short-axis in crystal reference
figure,
plot(hs.symmetrise,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal','fundamentalRegion','figSize','medium')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
mtexTitle(sprintf('%s grain short-axes',mineral))
saveFigure(sprintf('%s_%s_shortAxisXTAL.png',pname,mineral),'-bestfit')




%**************************************************************************
%%        S O M E   M I S O R I E N T A T I O N S
%**************************************************************************
%% Low-angle (2–10°) boundary misorientation analysis

% Note, the particular approach used here requires that you specified a
% second angle for inner-boundaries during grain reconstructuion (as above)

% get all the subgrain/inner boundaries from the grainset
subB = grains.innerBoundary;

% boundaries of the phase of interest (forsterite in our example case)
subFo = subB(phase,phase);

% condition for only misorientations with angles of 2–10° 
cond = subFo.misorientation.angle>=1*degree & subFo.misorientation.angle<10*degree;

% the "low-angle" boundaries of interest
labFo = subFo(cond);

%% plot pole figure of misorientation axes in the crystal reference frame
figure,
plotAxisDistribution(labFo.misorientation,'antipodal','lower','smooth','halfwidth',15*degree,'figSize','small')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
saveFigure(sprintf('%s_lowAngleBoundaryAx_XTAL.png',pname),'-bestfit')



%% how to invert to specimen reference frame?

% first get all the crystal orientations along the boundaries
o = ebsd('id',labFo.ebsdId).orientations;

% compute the axis of rotation in specimen reference frame
misAx = axis(o(:,1),o(:,2));
% and angle
misAng = angle(o(:,1),o(:,2));

% plot
figure,
plot(misAx,'antipodal','lower','smooth','halfwidth',10*degree)
cb = mtexColorbar('Title','M.U.D.');
setColorRange([.5 max(cb.Limits)]);
saveFigure(sprintf('%s_lowAngleBoundaryAx_SPEC.png',pname),'-bestfit')


%% Just as an example... lets use LAB max to rotate CPO (maybe kinematic?)
% compute axis density (use same 'r' grid as before)
r = plotS2Grid('resolution',1*degree,'antipodal');
kdeLAB = calcDensity(misAx,r,'antipodal','halfwidth',10*degree);
% find index of point on the grid with the maximum density
[~,labI] = max(kdeLAB);
% axis at that grid point
labMaxSpec = r(labI);


% cross vector
notMax= cross(labMaxSpec,zvector);

% a rotation
rot2labMax = rotation.map(labMaxSpec,zvector,notMax,yvector);

%% a plot
figure,
plotPDF(rot2labMax*o,h,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
saveFigure(sprintf('%s_CP_after_labBased_rotation.png',pname),'-bestfit')







%**************************************************************************
%%      O R G A N I Z I N G    D A T A       ( B Y    S A M P L E )
%**************************************************************************

% Setup a data structure (to contain field and microstrcutrual data)
sample = struct();

% give the sample a name (we'll keep using the string in pname for now)
sammple.name = pname;

% add info about the field data for this sample
% foliation strike and dip
sample.strike = 31;
sample.dip = 39;
% lineation trend
sample.linTrend = 200;

% Since we have a sample and a field fabric, lets specify for our records
% which directions those are in our specimen reference frame (i.e., EBSD
% map).
sample.linSpec = vector3d.X;
sample.folPolSpec = -vector3d.Y;

%% Lets add some orientations and vector versions of eild fabric (similar to above)
% use the Fabrica function SDT2or to define a fabric orientation
[sample.fabOr, sample.fabRot, sample.strikeV, sample.poleV, sample.linV] = SDT2or(sample.strike,sample.dip,sample.linTrend);

%% add some of the microstructural data from the example dataset
% just as an example that you can store anything in a structure to keep it
% all together.
% EBSD
sample.ebsd = ebsd;
% Grains
sample.grains = grains;
% grains used for spo analysis
sample.gSPO = gSPO;

%% This data structure will contain a variety of data for one sample
% first to make a pole figure (as above)
os = sample.ebsd('f').orientations;
cs = o.CS;
hs = [Miller(1,0,0,'direction',cs),Miller(0,1,0,'direction',cs),Miller(0,0,1,'direction',cs)];

figure,
plotPDF(os,hs,'antipodal','lower','smooth','halfwidth',10*degree,'colorrange','equal')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);


%% rotating data from specimen reference frame to "geographic"

% lets define it for this sample and store it in the data structure
% use the rotation 'map' option to specify exactly how two directions in
% the specimen reference frame should "map" to the new reference frame,
% thereby defining a unique rotation between the two. In our example, let's
% say that the sample was cut relative to the lineation measured in the
% field and the pole to the foliation plane. Specifically, the sample was
% cut such that the linetion points towards the positive x-direction in our
% maps, and the downward pole to foliation is aligned with the negative
% y-direction. So we specify the following rotation:
sample.rotSpec2Geo = rotation.map(xvector,sample.linV,-yvector,sample.poleV);



%% Now we plot a "geographic" version of the crystal orientation data

% rotate the rotations using the rotation computed above
og = sample.rotSpec2Geo.*os;

% plot
figure,
plotPDF(og,hs,'antipodal','lower','smooth','halfwidth',15*degree,'colorrange','equal')
cb = mtexColorbar('Title','M.U.D.');
setColorRange([1 max(cb.Limits)]);
hold on
% add structural data from the field for context
plot(fols,'plane','LineColor','k','LineWidth',1,'antipodal','lower','add2all')
% add lineations
plot(lins,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower','add2all')
% poles to foliation
plot(fols,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10,'add2all')
saveFigure(sprintf('%s_Geo_FieldFabs_and_CPO.png',pname),'-bestfit')



%% Now do the same with the specimen low-angle boundaries
sample.labFoMisXTAL = labFo.misorientation;
sample.labFoAxSPEC = misAx;
sample.labFo = labFo;
figure,
plot(sample.rotSpec2Geo.*sample.labFoAxSPEC,'antipodal','lower','smooth','halfwidth',10*degree,'figSize','large')
% mtexColorMap white2black
cb = mtexColorbar('Title','M.U.D.');
setColorRange([.5 max(cb.Limits)]);
hold on
% add structural data from the field for context
plot(fols,'plane','LineColor','k','LineWidth',1,'antipodal','lower')
% lineations
plot(lins,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k',...
    'MarkerSize',10,'antipodal','lower')
% poles to foliation
plot(fols,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerSize',10)
hold on
plot(sample.poleV,'antipodal','lower','plane','LineWidth',6,'LineColor','m')
plot(sample.poleV,'antipodal','lower','plane','LineWidth',3,'LineColor','k')
plot(sample.linV,'antipodal','lower','Marker','o','MarkerSize',20,'MarkerFaceColor','m','MarkerEdgeColor','k')
plot(sample.poleV,'antipodal','lower','Marker','s','MarkerSize',20,'MarkerFaceColor','m','MarkerEdgeColor','k')
saveFigure(sprintf('%s_Geo_FieldFabs_and_lowAngle_boundary_Axes.png',pname),'-bestfit')


