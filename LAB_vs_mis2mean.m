%% 'mis2mean' VS low-angle-boundary misorientations

% Misorientation axes are axes of rotation between any two orientations.
% When computing grains in MTEX, one type of intragranular misorientation
% can be computed – the misorientation of each point to the mean
% orientation of the grain. These misorientations are often used for making
% spatial plots of the ebsd data colored by the misorientation angle
% magnitude to visually assess the intragranular distortion of grains.
% *However*, these 'mis2mean' misorientations are – by definition – *not*
% boundary misorientations (between two *neighboring* orientations). In the
% case of intragranular boundary misorientations, the misorientation axes
% are commonly used to infer the rotation axes of slip systems, by assuming
% that the curvature between neighboring points is associated with slip on
% a discrete plane between the points. Plots of low-angle boundaries
% axes will yield different patterns than mis2mean, and sometimes the
% results are very different. For inferring slip systems, it is therefore
% very important to use a neighbor-neighbor misorientation computation such
% as the one included below.

% load example data
mtexdata forsterite

% phase of interest
phase = 'f';

% compute grains with mis2mean and using low-angle inner-boundary spec
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'boundary','tight','angle',[10 1]*degree);

% remove small grains
grains = grains(grains.grainSize>9);
ebsd = ebsd(grains);

% re-compute grains 
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'boundary','tight','angle',[10 1]*degree);

% smooth grains
grains = smooth(grains,5,'moveTriplePoints');


%% First: mis2mean
% get the misorientations to grain-mean for the phase of interest
mis2mean = ebsd(phase).mis2mean;

% condition for only misorientations with angles of 2–10° 
condM2M = mis2mean.angle>=2*degree & mis2mean.angle<10*degree;

% isolate 2-10° misorientations
mis2meanLow = mis2mean(condM2M);

% plot in crystal reference frame
figure,
plotAxisDistribution(mis2meanLow,'antipodal','lower','smooth','halfwidth',10*degree,'figSize','small')
mtexColorbar('Title','M.U.D.');
mtexTitle('mis2mean axes')


% In specimen/spatial reference frame...
% axis of rotation between each point and the mean orientation of the same grain
ax = axis(ebsd(phase).orientations,grains(ebsd(phase).grainId).meanOrientation);

% plot
figure,
plot(ax,'antipodal','lower','smooth','halfwidth',10*degree)
cb = mtexColorbar('Title','M.U.D.');
mtexTitle('mis2mean axes')

%% Second: Discrete low-angle (2–10°) boundary misorientation analysis

% Note, the particular approach used here requires that you specify a
% second angle for inner-boundaries during grain reconstructuion (as above)

% get all the subgrain/inner boundaries from the grainset
subB = grains.innerBoundary;

% boundaries of the phase of interest (forsterite in our example case)
subFo = subB(phase,phase);

% condition for only misorientations with angles of 2–10° 
condLAB = subFo.misorientation.angle>=2*degree & subFo.misorientation.angle<10*degree;

% the "low-angle" boundaries of interest
labFo = subFo(condLAB);

% plot crystal reference frame 
figure,
plotAxisDistribution(labFo.misorientation,'antipodal','lower','smooth','halfwidth',10*degree,'figSize','small')
mtexColorbar('Title','M.U.D.');
mtexTitle('boundary axes')


% plot specimen/spatial reference frame
% first get all the crystal orientations along the boundaries
oLab = ebsd('id',labFo.ebsdId).orientations;

% compute the axis of rotation in specimen reference frame
misAx = axis(oLab(:,1),oLab(:,2));

% and angle
misAng = angle(oLab(:,1),oLab(:,2));

% plot
figure,
plot(misAx,'antipodal','lower','smooth','halfwidth',10*degree)
mtexColorbar('Title','M.U.D.');
mtexTitle('boundary axes')

%% Plot the mis2mean and low-angle boundaries spatially
figure,
plot(ebsd('indexed'),ebsd('indexed').mis2mean.angle./degree,'figsize','large')
cb = mtexColorbar('Title','\Theta (degrees)');
mtexColorMap(flipud(gray))
hold on
plot(grains.boundary,'LineWidth',1.5,'DisplayName','grain')
plot(labFo,'LineColor','r','LineWidth',3,'DisplayName','subgrain')
