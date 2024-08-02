%% Photonic Material Calculation 
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Modeling parameters                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = (0.4:0.005:0.7)'; % wavelengths for which to calculate scattering in µm

% parameter determining scattering strength of top coat
thetaScatter_topcoat = 3.8; %degrees, average of scatter angle of FWHM in all directions

% parameter determining scattering strength of reflector
thetaScatterX_pp = 0.01; %degrees of scatter angle of FWHM, in horizontal direction
thetaScatterY_pp = 0.01; %degrees of scatter angle of FWHM, in vertical direction

%for ii = 1:size(thetaScatterX_pp,2)
%for jj = 1:size(thetaScatter_topcoat,2)

% parameters to define photonic material scattering characteristics
refrInd = 1.5; % average refractive index of material
lambda0 = 0.675; % wavelength at which material was exposed during holographic manufacture
FWHM0 = 0.005; % FWHM of red filter used in holographic manufacture
%sigma0 = 1/(FWHM0/2.35)/(2*pi); % FWHM to standard deviation conversion
sigma0 = 1;
    lambda_m = lambda0/refrInd; % in µm, wavelength with which photonic material was manufactured 
    kM = 2*pi/lambda_m; % reasonable reference lengthscale for spatial frequencies

% angles defining light incidence direction 
theta_in = 10; 
phi_in = 0;

% observation angle resolution
dthetaObs = 5;
dphiObs = 10;

% resolution of integration variables 
NMesh = 4; 
% NMesh = 1 takes 2 seconds; = 2 takes 18 seconds; = 3 takes 244 seconds
% NMesh = 4 takes 3188 seconds; = 5 takes X seconds
% NMesh = 6 takes X seconds; = 7 takes X seconds

brightness = 50; %adjustment parameter for color conversion 
% Used to adjust brightness in XYZ to RGB conversion for better visibility
% (akin to changing exposure time on camera) 

% Choice of locations to pull out spectra from 
phiChoiceSpec = [0; 0; 0; 0];
thetaChoiceSpec = [10; 20; 30; 50];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Functions needed to model scattering of topcoat and photonic material         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Function definition for the spatial frequency representation of %%%%%%%%%%%%%%%%%%%%%%%
%   a photonic structure produced with Lippmann process on rough shimstock                      
freqRepFunc = @(Kx, Ky, Kz, kM, sigma0, sigmaX_pp, sigmaY_pp) exp(-abs((Kz-kM).^2 + Kx.^2 + Ky.^2 - kM^2).^2/(kM*sigma0)) .* ... % sphere with normal distribution, with Gaussian width determined by FWHM of light exposure
  exp(-1/2 * (Kx./sigmaX_pp).^2) .* ... % Gaussian in X direction based on scattering angle of reflective surface
  exp(-1/2 .* (Ky./sigmaY_pp).^2);  % Gaussian in Y direction based on scattering angle of reflective surface      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for greater accuracy, this should be a sinc function around a sphere


%%% Vector difference between incident and scattered light wave vector %%%%%%%%%%%%%%%%%%%%
matKVecX = @(k, theta_in, phi_in, theta_out, phi_out) ...
            k .* (sind(theta_out).*cosd(phi_out) - sind(theta_in).*cosd(phi_in));

matKVecY = @(k, theta_in, phi_in, theta_out, phi_out) ...
            k .* (sind(theta_out).*sind(phi_out) - sind(theta_in).*sind(phi_in));

matKVecZ = @(k, theta_in, phi_in, theta_out, phi_out) ...
            k .* (cosd(theta_out) + cosd(theta_in));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: written separately for each coordinate to allow using 3D matrices as input for faster calculations


%%% Scattering function rough topcoat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DeltaK = @(k, theta, phi, theta_in, phi_in) k.*sqrt(2*(1-sind(theta).*sind(theta_in).*cosd(phi - phi_in) - cosd(theta).*cosd(theta_in)));
ScatterFuncDeltaK = @(deltaK, sigma) 1/sigma/sqrt(2*pi).*exp(-1/2*deltaK.^2/sigma^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Building inputs: integration variable value matrices, k-Vectors, data arrays       %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Preparing inputs')


% establishing set of integration angles theta1, phi1, theta2, phi2, which
% are uniformly distributed across the relevant angular space (half sphere of
% reflection) 

[V,F] = icosphere(NMesh); % using icosphere script from Wil O.C. Ward 19/03/2015, University of Nottingham, UK

rotAngle = 31.71748; % angle needed to rotate outputs from icosphere to have vertex distribution symetric around normal 

% rotation matrix to align icosphere output with sample surface normal
RotMatrix = [cosd(rotAngle) 0 -sind(rotAngle); sind(rotAngle) 0 cosd(rotAngle); 0 1 0]; 
Vrot = V*RotMatrix; % doing the rotation 

% Picking the upper half sphere of vertices in the next three lines 
VTemplate = Vrot;
VTemplate(find(Vrot(:,3)<-0.0001),:) = NaN; 
VHalf = VTemplate(~isnan(VTemplate(:,1)),:);

% plotting wave vector directions established with icosphere script to
% % check for uniform coverage of integration angles 
% figure()
%     plot3(VHalf(:,1), VHalf(:,2), VHalf(:,3), 'o', 'MarkerSize', 5, 'Color', 'blue')
%     set(gca, 'Fontsize', 18)
%     xlabel('k_x / k')
%     ylabel('k_y / k')
%     zlabel('k_z / k')
%     title('integration variable grid')
%     set(gcf,'Color', 'white')
%     box on 
%     drawnow


% defining integration angles based on processed icosphere output     
theta1 = (acosd(VHalf(:,3)))';
phi1 = (180/pi*atan2(VHalf(:,2),VHalf(:,1)))';
theta2 = theta1;
phi2 = phi1;

% observation angles 
thetaObs = 0:dthetaObs:90;
phiObs = 0:dphiObs:360; 

lambdaMat = repmat(lambda, 1, size(theta1,2));
theta1Mat = repmat(theta1, size(lambda,1), 1);
phi1Mat = repmat(phi1, size(lambda,1), 1);

theta2Mat = theta1Mat;
phi2Mat = phi1Mat;

% k-Vector matrix
kVecsMat = 2*pi./lambdaMat;

    %%% Secondary scattering conversions
    FWHMx_pp = tand(thetaScatterX_pp)*2*kM; 
    FWHMy_pp= tand(thetaScatterY_pp)*2*kM;
    % convert scatter angle of reflector surface to FWHM of Gaussian distribution of intensity in K-space
    
    sigmaX_pp = FWHMx_pp/(2*sqrt(2*log(2))); 
    sigmaY_pp = FWHMy_pp/(2*sqrt(2*log(2)));
    % convert FWHM to standard deviation of Gaussian distribution 
    % FWHM = 2*sigma*sqrt(2*log(2))

    %%% add in scatter component
    FWHM_topcoat = tand(thetaScatter_topcoat) * 2 * max(refrInd*max(kVecsMat,[],'all'));
    % convert scatter angle of topcoat to FWHM of Gaussian distribution of intensity

    sigma_topcoat = FWHM_topcoat/(2*sqrt(2*log(2)));
    % convert FWHM to standard deviation of Gaussian distribution

% temporal data array used in calculations 
tempArray = zeros(size(lambda,1),size(phi2,2));

% array to store spectrum for different observation angles 
SpectrArray = zeros(size(lambda,1),size(thetaObs,2),size(phiObs,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Actual calculation                                     %
% four loops running over the second pair of integration                                  %        
% variables theta2, phi2, and the pairs of observation angles (thetaObs,                  %
% phiObs); theta1, phi1, and lambda are all captured in 3D matrices and                   %
% don't need loops                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting calculations')
totalTime = tic; % counter for total calculation time

% Prepare S_in matrix (only needs to be calculated once, doesn't change
deltaKIncScatter = DeltaK(refrInd*kVecsMat, theta1Mat, phi1Mat, theta_in, phi_in); % deviation of incident scattered light from incident light direction
S_in = ScatterFuncDeltaK(deltaKIncScatter, sigma_topcoat); % intensity associated with incident scattered light 

for aa = 1:1:size(thetaObs,2)
    thetaLoopTime = tic;
    for bb = 1:1:size(phiObs,2) 
        phiLoopTime = tic;
        for ii = 1:1:size(theta2,2)       
                structKVecX = matKVecX(refrInd*kVecsMat, theta2(ii), phi2(ii), theta1Mat, phi1Mat);  % relevant K-vector of photonic material that would enable scattering of this k_in into this k_out          
                structKVecY = matKVecY(refrInd*kVecsMat, theta2(ii), phi2(ii), theta1Mat, phi1Mat);  % relevant K-vector of photonic material that would enable scattering of this k_in into this k_out          
                structKVecZ = matKVecZ(refrInd*kVecsMat, theta2(ii), phi2(ii), theta1Mat, phi1Mat);  % relevant K-vector of photonic material that would enable scattering of this k_in into this k_out                                                     
                temp = S_in.*freqRepFunc(structKVecX,structKVecY,structKVecZ, kM, sigma0, sigmaX_pp, sigmaY_pp); % photonic material scattering behavior 
                tempArray(:,ii) = sum(temp,2); % integration over theta1 and phi1
        end
        deltaKOutScatter = DeltaK(refrInd*kVecsMat, thetaObs(aa), phiObs(bb), theta2Mat, phi2Mat); % deviation of outgoing scattered light from incident light direction
        S_out = ScatterFuncDeltaK(deltaKOutScatter, sigma_topcoat); % intensity associated with outgoing scattered light            
        SpectrArray(:,aa,bb) = sum(S_out.*tempArray, 2);
        %disp (['phiObs = ' num2str(phiObs(bb)) 'º - duration: ' num2str(toc(phiLoopTime), '%.2f') ' seconds'])
    end
    disp (['thetaObs = ' num2str(thetaObs(aa)) 'º - duration: ' num2str(toc(thetaLoopTime), '%.2f') ' seconds'])
end
disp (['Total duration: ' num2str(toc(totalTime), '%.2f') ' seconds'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      Normalization                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SpectrArrayRaw = SpectrArray; % save original data before norming 
TotalIntensityPerWavelength = sum(SpectrArrayRaw, [2 3]); % determining normalization factor
SpectrArray = SpectrArrayRaw/max(TotalIntensityPerWavelength); % Normalization


save_location = '/Users/Andrew/Documents/MATLAB/Examples/R2022a/stats/EvaluateAStandardMultivariateNormalPDFExample/Surface coating project (MSRP Summer 23)/OutputArray Matrices';
%filename = "Parameters_Topcoat"+thetaScatter_topcoat+"_PPX"+thetaScatterX_pp+"_PPY"+thetaScatterY_pp+"_theta"+theta_in+"_phi"+phi_in+"_dtO"+dthetaObs+"_dpO"+dphiObs;
%fig_file = fullfile(save_location , filename);

%save(fig_file,'SpectrArray');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             CIE color coordinates and CIE to RGB conversion                             %
%             after http://www.brucelindbloom.com/index.html?Math.html                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
disp('Starting Spectra to CIE and RGB conversions')
tic

dLambda = lambda(2) - lambda(1); %wavelength step

% Lightsource data b
% Standard illuminant D65 daylight
D65_Data = dlmread('D65.csv');
LightsourceSPDInterpol = interp1(D65_Data(:,1)/1000,D65_Data(:,2),lambda, 'pchip', 0);

% CIE tristimulus curves
CIEXYZ_Data = dlmread('ciexyz31.csv');
CIE_1 = interp1(CIEXYZ_Data(:,1)/1000,CIEXYZ_Data(:,2),lambda, 'pchip', 0);
CIE_2 = interp1(CIEXYZ_Data(:,1)/1000,CIEXYZ_Data(:,3),lambda, 'pchip', 0);
CIE_3 = interp1(CIEXYZ_Data(:,1)/1000,CIEXYZ_Data(:,4),lambda, 'pchip', 0);

% Used to adjust brightness in XYZ to RGB conversion for better visibility
% (akin to changing exposure time on camera) 
%brightness = 100; 

% normalization factor for CIE values
NormFac = squeeze(sum(LightsourceSPDInterpol .* CIE_2 .* dLambda));

% preallocation of space for XYZ and RGB arrays (holding XYZ values and
% RGB values)
XYZ = zeros(size(SpectrArray,2), size(SpectrArray,3), 3);
RGB = zeros(size(SpectrArray,2), size(SpectrArray,3), 3);

for ii = 1:1:size(XYZ,1)
    for jj = 1:1:size(XYZ,2)
        signal = brightness*SpectrArray(:, ii, jj);
        XYZ(ii,jj,1) = 1/NormFac * sum(signal.*LightsourceSPDInterpol.*CIE_1.*dLambda);
        XYZ(ii,jj,2) = 1/NormFac * sum(signal.*LightsourceSPDInterpol.*CIE_2.*dLambda);
        XYZ(ii,jj,3) = 1/NormFac * sum(signal.*LightsourceSPDInterpol.*CIE_3.*dLambda);
    
        RGB(ii,jj,:) = xyz2rgb(XYZ(ii,jj,:),'WhitePoint','d65');
        RGB(ii,jj,:) = min(max(RGB(ii,jj,:),0),1); % bounding rgb data to interval [0,1] 
    end
end  
disp(['Finished CIE and RGB calculations - duration: ' num2str(toc, '%.2f') ' seconds.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                  Creating and saving the output figure                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Creating final outputs')
tic

% location matrices needed for plotting the observed color distributions 
[phiOutMat, thetaOutMat] = meshgrid(phiObs, thetaObs);
xMat = sind(thetaOutMat).*cosd(phiOutMat);
yMat = sind(thetaOutMat).*sind(phiOutMat);
zMat = cosd(thetaOutMat);

% angle markings on plot
phiMarkings = 0:1:360;
thetaMarkings = [20; 40; 60; 90];
xRingMarkings = sind(thetaMarkings)*cosd(phiMarkings);
yRingMarkings = sind(thetaMarkings)*sind(phiMarkings);

phiChoice = [0, 45, 90, 135];
LineMarkings = zeros(size(phiChoice,2), 2, 2);
for ii = 1:1:size(phiChoice,2)
    LineMarkings(ii,:,:) = [cosd(phiChoice(ii)) -sind(phiChoice(ii)); sind(phiChoice(ii)), cosd(phiChoice(ii))] ...
                            * [-1 1; 0 0];
end

% plot of scattering distribution in RGB color
dataFig1 = figure('units', 'points', 'Position', [200, 200, 600, 600]);
    tcolor(xMat,yMat,RGB,'normal') 
    shading flat
    set(gca,'color',[0 0 0])
    set(gca,'xtick',[],'ytick',[])
    hold on
    for ii = 1:1:size(thetaMarkings)-1
        plot(xRingMarkings(ii,:), yRingMarkings(ii,:), 'Color', 'white', 'LineStyle',':', 'LineWidth',1.5)
    end
    plot(xRingMarkings(end,:), yRingMarkings(end,:), 'Color', 'white', 'LineStyle','-', 'LineWidth',1.5)

    for ii = 1:1:size(LineMarkings,1)
        plot(squeeze(LineMarkings(ii,1,:)), squeeze(LineMarkings(ii,2,:)), 'Color', 'white', 'LineStyle',':', 'LineWidth',1.5)
    end
    plot(sind(theta_in)*cosd(phi_in+180), sind(theta_in)*sind(phi_in+180), 'o', 'Markersize', 7, 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'none', 'Linewidth', 2)

    xlim([-1.1 1.1])
    ylim([-1.1 1.1])
    set(gcf,'Color','white')
    axis equal
    drawnow

% Pulling out specific spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will fail, if you choose angle value pairs that have not been simulated;

spectraChoice = zeros(size(SpectrArray,1), size(phiChoiceSpec,1));
for ii = 1:1:size(phiChoiceSpec,1)
    spectraChoice(:,ii) = SpectrArray(:, find(thetaObs == thetaChoiceSpec(ii)), find(phiObs == phiChoiceSpec(ii)));
end


XYZChoice = zeros(size(phiChoiceSpec,1),3);
RGBChoice = zeros(size(phiChoiceSpec,1),3);
XYZChoiceChroma = zeros(size(phiChoiceSpec,1),3);
RGBChoiceChroma = zeros(size(phiChoiceSpec,1),3);

for ii = 1:1:size(phiChoiceSpec,1)
% RGB color of total specular spectrum 
    XYZChoice(ii,1) = 1/NormFac * sum(brightness*spectraChoice(:,ii).*LightsourceSPDInterpol.*CIE_1.*dLambda);
    XYZChoice(ii,2) = 1/NormFac * sum(brightness*spectraChoice(:,ii).*LightsourceSPDInterpol.*CIE_2.*dLambda);
    XYZChoice(ii,3) = 1/NormFac * sum(brightness*spectraChoice(:,ii).*LightsourceSPDInterpol.*CIE_3.*dLambda);
        
    XYZChoiceChroma(ii,:) = XYZChoice(ii,:)./sum(XYZChoice(ii,:));
    
    RGBChoice(ii,:) = xyz2rgb(XYZChoice(ii,:),'WhitePoint','d65');
    RGBChoice(ii,:) = min(max(RGBChoice(ii,:),0),1); % bounding rgb data to interval [0,1]
    
    RGBChoiceChroma(ii,:) = xyz2rgb(XYZChoiceChroma(ii,:),'WhitePoint','d65');
    RGBChoiceChroma(ii,:) = min(max(RGBChoiceChroma(ii,:),0),1); % bounding rgb data to interval [0,1]
end

% generating faces that connect vertices produced by icosphere script to
% plot the distribution of the scattering function S_in 

FHalf = delaunayTriangulation(VHalf(:,1),VHalf(:,2),VHalf(:,3));
FHull = convexHull(FHalf); % faces that connect the vertices as a convex envelope 

index = 1; 
S_in_Slice = (S_in(index,:))';
faceColor = S_in_Slice/max(S_in_Slice); % MK check how this pans out without normalization 
mirrorMatrix = [-1 0 0; 0 1 0; 0 0 1]; % light incidence is opposite of reflected light with respect toi sample normal 


multFreq = 2.5; % range of spatial frequencies for which to evaluate spatial frequency representation of photonic material in temrs of multiples of 2*pi/lambda_m
divFreq = 100; % divider of spatial frequency unit to establish Kx and Ky

%matrices with k-space coordinates 
Kx = -multFreq*kM:kM/divFreq:multFreq*kM;
[KxMat, KyMat, KzMat] = meshgrid(Kx,Kx,Kx);

% frequency representation of photonic material; this is only one of the
% 'lobes', the other one should be in the mirror position, with the mirror
% plane being the xy plane but it is not important for our construction,
% because there is no k_inc, k_sc combination that can reach that for
% elastic scattering 
freqRep = freqRepFunc(KxMat, KyMat, KzMat, kM, sigma0, sigmaX_pp, sigmaY_pp);

% Set the line below kM equal to zero
j = KxMat(1,:,1);
k = find(abs(j-kM)<1e-06);
freqRep(:,:,1:k) = 0;

sliceY = 251;
sliceX = 251;

% plotting a slice of the 3D frequency representation of the photonic
% material 
circX = -kM:kM/100:kM;
circY1 = sqrt(kM^2 - circX.^2)+kM;
circY2 = -sqrt(kM^2 - circX.^2)+kM;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                Plotting everything                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

dataFig1 = figure('units', 'points', 'Position', [100, 100, 1900, 600]);   
    tiledlayout(size(phiChoiceSpec,1),4);
    nexttile(1,[size(phiChoiceSpec,1)/2 1])
        patch('Faces',FHull,'Vertices',VHalf*mirrorMatrix, ...
          'FaceVertexCData', faceColor, 'FaceColor','interp', ...
          'Edgecolor', 'none')
        colormap("gray")
        set(gca, 'Fontsize', 18)
        xlabel('k_x / k')
        ylabel('k_y / k')
        zlabel('k_z / k')
        set(gcf, 'Color', 'white')
        box on

    nexttile(4*size(phiChoiceSpec,1)/2+1,[size(phiChoiceSpec,1)/2 1])
        h = slice(KxMat,KyMat,KzMat,freqRep,[-15:0.5:15],[],[]);
        set(h,'EdgeColor','none',...
        'FaceColor','interp',...
        'FaceAlpha','interp')
        alpha('color')
        xlim([min(KxMat,[],"all"), max(KxMat,[],"all")])
        ylim([min(KyMat,[],"all"), max(KyMat,[],"all")])
        zlim([min(KzMat,[],"all"), max(KzMat,[],"all")])
        c = colorbar();
        ylabel(c,'Amplitude','FontSize',18)    
        set(gca,'FontSize', 18)
        line([0 0], [0 0], [min(KzMat,[],"all") max(KzMat,[],"all")], 'Linestyle', '-', 'Linewidth', 1.5, 'Color', [0 0 0])
        line([min(KxMat,[],"all") max(KxMat,[],"all")], [0 0], [0 0], 'Linestyle', '-', 'Linewidth', 1.5, 'Color', [0 0 0]) 
        line([0 0], [min(KyMat,[],"all") max(KyMat,[],"all")], [0 0], 'Linestyle', '-', 'Linewidth', 1.5, 'Color', [0 0 0]) 
        xlabel('K_x [1/µm]')
        ylabel('K_y [1/µm]')
        zlabel('K_z [1/µm]')
        box on 
    
    nexttile(2,[size(phiChoiceSpec,1) 2]);
        tcolor(xMat,yMat,RGB,'triangles') 
        shading interp   
        % tcolor(xMat,yMat,RGB,'normal') 
        % shading flat
        set(gca,'color',[0 0 0])
        set(gca,'xtick',[],'ytick',[])
        hold on
        for ii = 1:1:size(thetaMarkings)-1
            plot(xRingMarkings(ii,:), yRingMarkings(ii,:), 'Color', 'white', 'LineStyle',':', 'LineWidth',1.5)
        end
        plot(xRingMarkings(end,:), yRingMarkings(end,:), 'Color', 'white', 'LineStyle','-', 'LineWidth',1.5)
    
        for ii = 1:1:size(LineMarkings,1)
            plot(squeeze(LineMarkings(ii,1,:)), squeeze(LineMarkings(ii,2,:)), 'Color', 'white', 'LineStyle',':', 'LineWidth',1.5)
        end
        plot(sind(theta_in)*cosd(phi_in+180), sind(theta_in)*sind(phi_in+180), 'o', 'Markersize', 7, 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'none', 'Linewidth', 2)
        
        for ii = 1:1:size(phiChoiceSpec,1)
            plot(sind(thetaChoiceSpec(ii))*cosd(phiChoiceSpec(ii)), sind(thetaChoiceSpec(ii))*sind(phiChoiceSpec(ii)), 'o', 'Markersize', 7, 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'none', 'Linewidth', 2)
        end

        xlim([-1.1 1.1])
        ylim([-1.1 1.1])
        set(gcf,'Color','white')
        axis equal
    
   for ii = 1:1:size(phiChoiceSpec,1)
       nexttile(4+4*(ii-1))
         plot(lambda, spectraChoice(:,ii), 'Linewidth', 2.5, 'Color', RGBChoice(ii,:))   
         hold on
         %plot(lambda, spectraChoice(:,ii), 'Linewidth', 1.5, 'Color', RGBChoiceChroma(ii,:), 'Linestyle', '--')           
         set(gca, 'Fontsize', 18)
         ylabel('Intensity')
         ylim([0 1.1*max(spectraChoice,[],'all')])
   end
        xlabel('\lambda [nm]')
        box on
   set(gcf,'color', 'white')


finalfilename = ['FullSystem-CenterWL_' num2str(lambda0*1000) '-PPX_' num2str(thetaScatterX_pp) ...
            '-PPY_' num2str(thetaScatterY_pp) '-TCdeg_' num2str(thetaScatter_topcoat) ... 
            '-incTheta_' num2str(theta_in) '-incPhi_' num2str(phi_in) '-ResIntMesh' num2str(NMesh) ... 
            '-brightness' num2str(brightness) '.png'];
finalfig_file = fullfile(save_location , finalfilename);

exportgraphics(dataFig1, finalfig_file, 'Resolution', 600) 

disp(['Data plotting - duration: ' num2str(toc, '%.2f') 's.']);

% end % for looping more than one scatter parameter
% end % for looping more than one mirror parameter


%% OUTPUT SCATTER FUN ONLY
% Prepare S_in matrix (only needs to be calculated once, doesn't change
deltaKIncScatter = DeltaK(refrInd*kVecsMat, theta1Mat, phi1Mat, theta_in, phi_in); % deviation of incident scattered light from incident light direction
S_in = ScatterFuncDeltaK(deltaKIncScatter, sigma); % intensity associated with incident scattered light 


FHalf = delaunayTriangulation(VHalf(:,1),VHalf(:,2),VHalf(:,3));
FHull = convexHull(FHalf); % faces that connect the vertices as a convex envelope 

index = 1; 
S_in_Slice = (S_in(index,:))';
faceColor = S_in_Slice/max(S_in_Slice); % MK check how this pans out without normalization 
mirrorMatrix = [-1 0 0; 0 1 0; 0 0 1]; % light incidence is opposite of reflected light with respect toi sample normal 


multFreq = 2.5; % range of spatial frequencies for which to evaluate spatial frequency representation of photonic material in temrs of multiples of 2*pi/lambda_m
divFreq = 100; % divider of spatial frequency unit to establish Kx and Ky

figure();
        patch('Faces',FHull,'Vertices',VHalf*mirrorMatrix, ...
          'FaceVertexCData', faceColor, 'FaceColor','interp', ...
          'Edgecolor', 'none')
        colormap("gray")
        set(gca, 'Fontsize', 18)
        xlabel('k_x / k')
        ylabel('k_y / k')
        zlabel('k_z / k')
        set(gcf, 'Color', 'white')
        box on