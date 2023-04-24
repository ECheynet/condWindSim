function [u,uCheck,nodes] = condWindSim2D(u0,t,y,z,newY0,newZ0,Cuy,Cuz,varargin)
% [u,uCheck,nodes] = condWindSim(u0,t,y0,newY,Cuy) generate 
% wind speed time series based on preexisting time series in query points, 
% in a 2D plane.
% 
% The function takes the following inputs:
% 
%     u0: [N x Ny] matrix of preexisting time series velocity data in Ny query points
%     t: [1 x N] vector time
%     y: [1 x Ny] vector representing the position of the query points along a horizontal line
%     z: [1 x Ny] vector representing the height of the query points along a vertical line
%     newY0: [1 x Nm] vector of new points where the time series are simulated along a horizontal line
%     newZ0: [1 x Nm] vector of new points where the time series are simulated along a vertical line
%     Cuy: Davenport decay coefficient for the coherence
% 
% The function outputs the following:
% 
%     u: [N x Nm] matrix of conditionally generated time series based on u0
%     uCheck: [N x (Ny + Nm)] matrix of conditionally generated time series based,
% % complemented with the already preexisting query points.
%     nodes: Structure containing information on points where the time series are generated.
% 
% Optional input arguments include:
% 
%     Nwin: the number of window. If set to 1, no moving mean is applied, 
% otherwise, a moving mean is used to compute a time-dependent coherence function.
% 
% The function first combines the target nodes with the reference ones. 
% It then computes a moving mean (optional) to calculate a time-dependent 
% coherence function. It gets the power spectral density (PSD) and phases 
% for each known location. Then it replaces the random phases by the known
% phases from the target locations. 
% Finally, it generates turbulence in the new grid.
% 
% Author: E. Cheynet - UIb - last modified 24-04-2023

%% inputparser

p = inputParser();
p.CaseSensitive = false;
p.addOptional('Nwin',1);
p.parse(varargin{:});
Nwin = p.Results.Nwin ;
%% Combine target nodes with the reference ones

newY = [y(:); newY0(:)];
newZ = [z(:); newZ0(:)];
Nm = numel(newY);
[Ny,N] = size(u0);
%% moving mean (optional)
% The  moving mean is used to compute a time-dependent coherence function

if Nwin==1 % no moving mean
    movMeanU0 = mean(u0,2,'omitnan');
    Twin = [];
elseif isempty(Nwin) % if empty number of window, a 10 min window is used
    Twin = round(600/median(diff(t))); % window duration expressed in number of time step
else
    Twin = round((t(end)/Nwin)/median(diff(t)));% window duration in time step
end

if ~isempty(Twin)
    movMeanU0 = zeros(round(size(u0,2)/Twin),size(u0,1));
    for ii=1:Ny
        movMeanU0(ii,:) = decimate(movingmean(u0(ii,:),Twin,2),Twin,'fir');
    end
    Nwin = size(movMeanU0,1);
end

if isempty(Nwin) || Nwin~=1
    tMovMean = decimate(t,Twin,'fir');
    tMovMean(1)=0;s
elseif Nwin==1
    tMovMean = t;
end

% plot(y,z,'*');hold on;plot(newY,newZ,'o')


%% Get PSD and phases for each known location
dt = median(diff(t));
fs = 1/dt;
clear Su cocoh quadCoh f

if mod(N,2)==0
    Nfreq = round(N/2)+1;
else
    Nfreq = round(N/2);
end

phase = zeros(Nfreq,Ny);
Su = zeros(Ny,Nfreq);

for ii= 1:Ny
    [Su(ii,:),f] = periodogram(detrend(u0(ii,:),'constant'),rectwin(N),N,fs);
    Su(ii,:) = sqrt(Su(ii,:)); % to speed up the computation later on
    Y = fft(detrend(u0(ii,:),'constant'));
    Y = Y(1:Nfreq);
    phase(:,ii) = atan2(imag(Y),real(Y));
end
phase(phase<0)= phase(phase<0) + 2*pi;

%% Replace the random phases by the known phases from the target locations
phi = 2*pi*rand(Nfreq,Nm);% Random phases for unknown locations
indTarget = nan(1,Ny);
for ii=1:Ny
    dummyInd = find(newY(:)==y(ii) & newZ(:)==z(ii) );
    if ~isempty(dummyInd) && ~isnan(dummyInd)
        indTarget(ii) = dummyInd;
        phi(:,indTarget(ii)) =  phase(:,ii);
    end
end

%% generate turbulence in the new grid

dy = abs(newY(:)'-newY(:));
dz = abs(newZ(:)'-newZ(:));


meanU = nan(Nwin,Nm);
u = zeros(numel(t),Nm); % Output
for pp = 1:Nwin % For each wind window (Nwin> 1  if Non-stationary fluctuations)

    % start and ending point for the subsamples
    [~,indStart] = min(abs(tMovMean(pp)-t));
    if pp<Nwin
        [~,indEnd] = min(abs(tMovMean(pp+1)-t));
    else
        indEnd = numel(t);
    end
    % interpolate moving mean at each time window
    F1 = scatteredInterpolant(y(:),z(:),movMeanU0(:,pp));
    F1.Method = 'linear' ;
    F1.ExtrapolationMethod = 'linear'  ;
    meanU(pp,:)  = F1(newY,newZ);
    MeanUCoh = 1/2*(meanU(pp,:)+meanU(pp,:)');

    % Get A
    A = zeros(Nfreq,Nm);
    for ii=2:Nfreq
         [A(ii,:)] = getA(newY,newZ,y,z,Su(:,ii),f(ii),phi(ii,:),indTarget);
    end
    %% Generate the IFFT matrix
    Nu = [A(1:Nfreq,:) ; real(A(Nfreq,:)); conj(flipud(A(2:Nfreq,:)))];    
    % create a dummyU for each window and duration
    dummyU=real(ifft(Nu).*sqrt(Nfreq./(dt)));

    u(indStart:indEnd,:) = dummyU(indStart:indEnd,:) + meanU(pp,:);
    

end
%% Get a copy of u without the reference locations removed
uCheck = u; 
%% remove the target nodes which are not in the newY vector
[ind2remove]= find(ismember([newY(:),newZ(:)],[y(:),z(:)],'row'));
nodes.Zcheck = newZ;
nodes.Ycheck = newY;
u(:,ind2remove)=[];
newZ(ind2remove)= [];
newY(ind2remove)= [];
nodes.U = mean(u,'omitnan');
nodes.Y = newY;
nodes.Z = newZ;
for ii=1:Nm 
    nodes.name{ii} = strcat('N',num2str(ii));
end % names affected at each nodes


%% Nested functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A] = getA2D_all(indTarget,newSu,phi,f,Cy,Cz,meanU,dy,dz)

        if  indTarget>1
            % The frist element of the corrleation matrix governs the
            % phase correlation. So we have to apply the Cholesky
            % decomposition for both row and columns before and after the
            % reference location.

            indB = fliplr(1:indTarget); % backward
            cohU = cohDavenport(f,Cy,Cz,dy(indB,indB),dz(indB,indB),meanU(indB,indB));
            S = (newSu(indB)*newSu(indB)').*cohU;

            G = chol(S,'lower');
            A(indB) = G*exp(1i.*phi(indB)');

            indF = indTarget:Nm; % forward
            cohU = cohDavenport(f,Cy,Cz,dy(indF,indF),dz(indF,indF),meanU(indF,indF));
            S = (newSu(indF)*newSu(indF)').*cohU;
            G = chol(S,'lower');
            A(indF) = G*exp(1i.*phi(indF)');
        else
            % If the first element is already the top left element,
            % no need to change anything
            cohU = cohDavenport(f,Cy,Cz,dy,dz,meanU);
            S = (newSu*newSu').*cohU;
            [L,D]=ldl(S,'lower'); % a LDL decomposition is applied this time
            G = L*sqrt(D);
            A = G*exp(1i.*phi');
        end
    end

    function [coh] = cohDavenport(f,Cy,Cz,dy,dz,U)
        B = (dy.*Cy).^2 + (dz.*Cz).^2;
        coh = exp(-f./U.*sqrt(B));
    end

    function [A] = getA(newY,newZ,y,z,Su,f,phi,indTarget)
        % for each new nodes, put the element with random phases, weighted
        % by the nearest reference locations, expressed as a decaying
        % function following the Davenport model.
        % The first X locations are the reference locations

        % interpolate the target PSD at the new locations
        F0 = scatteredInterpolant(y,z,Su);
        F0.Method =  'linear' ;
        F0.ExtrapolationMethod =  'nearest'  ;
        newSu = F0(newY,newZ);
        % for each target node
        lowTrimat = nan(Nm,Ny); % lower triangular matrix
        Weight = nan(Nm,Ny); % Each weight is defined w.r.t to reference nodes
        for jj=1:Ny
            % Get Cholesky decomposition for all taghet nodes
            [lowTrimat(:,jj)] = getA2D_all(indTarget(jj), newSu,phi,...
                f,Cuy,Cuz,abs(MeanUCoh),dy,dz); % for each target nodes
            % Compute the weight for each target nodes using the Davenport
            % model.
            Weight(:,jj) = cohDavenport(f,Cuy,Cuz,...
                abs(newY(indTarget(jj))'-newY),abs(newZ(indTarget(jj))'-newZ),meanU(pp,:)');
            % No weight at the target location
            Weight(indTarget([1:jj-1,jj+1:end]),jj) = 0;
        end
        A = conj(sum(lowTrimat .* (Weight ./ sum(Weight, 2)), 2)');
    end
%%
end

