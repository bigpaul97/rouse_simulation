function [zeroEigIndex]=MainNUF(traa)
% Non-uniform friction environment

% random number seed
rng(traa)

% In the line below, put in the loop configuration file: 
% Five variables are needed: 
% 1. total num of lattice sites, L
% 2. total num of LEFs, N
% 3. positons of left legs at different times, l_sites_traj (A-by-N matrix)
% 4. position of right legs at different times, r_sites_traj (A-by-N matrix)
% 5. time series corresponding to the positions of LEF legs, ts_traj (A-by-1 vector)
% Examples of loop configuration file are given
load('rosette_24SMC.mat') 

% Load the vector for different friction coefficient as well
% Note that this vector is unitless, and only specifies the relative scale for the unmodified friction coefficient
% The variable name needs to be set to "VMu"
load('VMu_logn(0,log(2)).mat')


%% Adjustable parameters

zeta = 1.08e-6; % unmodified friction coefficient (in Ns/m): 1.08e-6 for mouse; 5.4e-8 for yeast
det = 1; % simulation timestep: 1s for mouse, 0.02s for yeast
taup_limit = 1e1; % total simulation time, in unit of tau_p (approximately)
b = 2e-7; % contour length represented by one bead: 2e-7 nm for mouse; 1e-8 for yeast
bc = 1; % boundary condition: 0, periodic; 1, free; 2, fixed.

%% Physical parameters

% Below are the physical parameters that often don't require a change
kb = 1.38e-23; % Boltzmann constant (J/K)
Temp = 300; % effective temperature (K)
lk = 138e-9; % chromatin Kuhn length (m)
Nk = b/lk; % number of Kuhn length per bead
k = 2*kb*Temp/Nk/lk^2; % spring constant of nearest-neighbor springs (N/m)
k2 = k; % spring constant of non-nearest-neighbor springs, which can be set to other preferred values

%% Initialization

stablept = 2e4; % first 2e4 points in the LEF simulation are discarded to ensure steady state, adjust accordingly
l_sites_traj(1:stablept,:) = [];
r_sites_traj(1:stablept,:) = [];
ts_traj = ts_traj-ts_traj(stablept);
ts_traj(1:stablept) = [];

% find the range of LEF simulation timepoint to be used
realtime_limit = taup_limit*det;
ted = find(ts_traj >= realtime_limit,1,'first'); 

% Produce the time series with regular interval det, with loop configuration points inserted
ts = [ts_traj(1:ted)',ts_traj(1):det:ts_traj(ted)]';
ts(1:ted,2) = 1:ted;
ts = sortrows(ts);
fst = find(ts(:,2)>0,1,'first');
for i = fst+1:size(ts,1)
    if ts(i,2) == 0
        ts(i,2) = ts(i-1,2);
    end
end
ts(1,:) = []; % 2nd column indicates when to update loop info (when the next number appears)
timeM = ts_traj(1):det:ts_traj(ted); % regular timepoints for recording

% Obtain the non-uniform friction matrix
FMatrix = diag(VMu); % Z=zeta*F
ZMatrix = FMatrix.*zeta; % Z matrix

%% Main loop

for loopin = 1:2 % 1: with loop, 2: without loop
    for dmrouse = 1:2 % each dimension is independent (2D here)
        for indts = 1:size(ts,1)
            disp(indts) % print time points
            np = indts+1;

            % delta T
            if indts == 1
                dt = ts(1,1);
            else
                dt = ts(indts,1)-ts(indts-1,1);
            end

            % Obtain the loop-free dynamical matrix
            VKappa = ones(L,1);
            VKappa = VKappa*(2*k);
            VKappaoff = ones(L-1,1)*(-k);
            KKappa = diag(VKappa,0)+diag(VKappaoff,1)+diag(VKappaoff,-1); % Dynamical matrix
            
            % Adjust the dynamical matrix for different boundary conditions
            if bc == 0 % periodic
                KKappa(end,1) = -k;
                KKappa(1,end) = -k;
            elseif bc == 1 % free
                KKappa(1,1) = k;
                KKappa(end,end) = k; 
            elseif bc == 2 % fixed
                % do nothing as KKappa is already for fixed b.c.
            end
    
            % In the looped case, incorporate loops in the dynamical matrix
            if loopin == 1 
                for kapi = 1:N
                    if l_sites_traj(ts(indts,2),kapi) > 0
                        KKappa(l_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi)) = ...
                            KKappa(l_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi))-1*k2;
                        KKappa(r_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi)) = ...
                            KKappa(r_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi))-1*k2;
                        KKappa(l_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi)) = ...
                            KKappa(l_sites_traj(ts(indts,2),kapi),l_sites_traj(ts(indts,2),kapi))+1*k2;
                        KKappa(r_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi)) = ...
                            KKappa(r_sites_traj(ts(indts,2),kapi),r_sites_traj(ts(indts,2),kapi))+1*k2;
                    end                   
                end
            end

            % dimensionless dynamical matrix A
            AMatrix = KKappa./k;
            
            % Diagonalization of -Z^{-1}*K
            [vecs,Lambda] = eig(-KKappa,ZMatrix);

            % Report error if eigenvectors or eigenvalues contain imaginary parts, something is not right
            if sum(sum(imag(vecs)~=0)) > 0 || sum(sum(imag(Lambda)~=0)) > 0
                error('some eigenvectors or eigenvalues are imaginary!!!')
            end
            
            lambda = diag(Lambda); % vector of eigenvalues
            
            % Check for number of zero eigenvalues by considering small numerical residues
            zeroEigIndex = 0;
            numZeroEig = 0;
            for i = 1:L
                if abs(lambda(i)) < 1e-15 % may need to change the threshold in some case
                    lambda(i) = 0;
                    zeroEigIndex = i;
                    numZeroEig = numZeroEig+1;
                end
            end

            % Report error if more than 1 zero eigenvalue was detected
            if numZeroEig > 1
                error('more than 1 zero eigenvalue, something went wrong!!!')
            end

            % check if there are repeated eigenvalues 
            % if so, find the indices of the corresponding eigenvectors and orthogonalize within the sub-eigenspace
            if length(lambda) ~= length(uniquetol(lambda,1e-10)) % tolerance level 1e-10
                [uniqEigenValue,~,ic] = uniquetol(lambda,1e-10); 
                eigenSpace_indexList = [];
                for i = 1:length(uniqEigenValue)
                    if length(find(ic==i)) == 1
                        continue
                    else
                        find_list = find(ic==i); % repeated eigenvalue index in lambda vector
                        for j = 1:length(find_list)
                            eigenSpace_indexList = [eigenSpace_indexList;find_list(j),i];
                        end
                    end
                    
                    % Update eigenvectors so that A^{1/2}*S has orthogonal columns
                    transfmedVecs = sqrtm(AMatrix)*vecs(:,ic==i); % transform to the space of A^{1/2}*{eigenspace}
                    orthVecs = orth(transfmedVecs); % Gram-schmidt orthogonalization
                    newEigVecs = sqrtm(AMatrix)\orthVecs; % transform back to the eigenspace
                    vecs(:,ic==i) = newEigVecs;
                    disp('repeated eigenvalue')
                end
            end

            % Calculate terms in Eq. 36
            for i = 1:L
                SMatrix(:,i) = vecs(:,i); % S-matrix
                expDecay(i) = exp(lambda(i)*dt); % exponential decay in the deterministic term
                expSD(i) = sqrt(1-expDecay(i)^2); % exponential term in the standard deviation (stochastic term)
            end
            
            % Calcualte D_2 matrix
            D2 = SMatrix'*AMatrix*SMatrix;
            
            % Set the zero eigenvalue to a very small positive value for convenience
            % This won't affect result since the eigenmode with zero eigenvalue is set to 0 later
            for i = 1:L
                if D2(i,i) <= 0
                    D2(i,i) = 1e-30;
                end
            end
            
            for i = 1:L
                SSSD(i) = sqrt(kb*Temp./(k*D2(i,i))); % steady-state standard deviation
                TSD(i) = SSSD(i)*expSD(i); % total standard deviation of normal modes 
                % (note: the entry of zero eigenvalue index in TSD is zero by multiplying with expSD)
            end

            % Initial condition for normal coordinates
            if np == 2
                for iniPW = 1:L                       
                    P1(iniPW) = randn(1)*SSSD(iniPW); % initial condition is using mean=0, variance=long-term variance
                end
                
                % set the eigenmode with zero eigenvalue to zero to avoid overall drift
                if zeroEigIndex ~= 0
                    P1(zeroEigIndex) = 0;
                end

                PhyWalk(:,1) = SMatrix*P1'; 
            end
            
            % Update the physical positions of all beads according to the algorithm
            PhyWalk(:,np) = SMatrix*diag(expDecay)*(SMatrix\PhyWalk(:,np-1)) + SMatrix*(TSD.*randn(1,L))';          
        end

        % save the position series for different cases: loop, no loop, physical dimensions
        if loopin == 1 && dmrouse == 1
            dimPhyWalk11 = PhyWalk;
        elseif loopin == 1 && dmrouse == 2
            dimPhyWalk12 = PhyWalk;
        elseif loopin == 2 && dmrouse == 1
            dimPhyWalk21 = PhyWalk;
        elseif loopin == 2 && dmrouse == 2
            dimPhyWalk22 = PhyWalk;
        end
    end
    
end

% Only save positions at the regular time gaps, specified by det
[~,loc] = ismember(timeM,ts(:,1));
x1 = dimPhyWalk11(:,loc);
y1 = dimPhyWalk12(:,loc);
x2 = dimPhyWalk21(:,loc);
y2 = dimPhyWalk22(:,loc);


%% Calculate MSD 

% Pick the beads of interest to calculate MSDs (change accordingly)
% bn1 = 1:25:600;
% bn2 = 13:25:600; % 600-bead simulation (for non-uniform rosette)
bn = 5:10:600; % 600-bead simulation

for jj = 1:length(bn)
    j = bn(jj);
    for i = 1:length(timeM)-1
        MSD1(jj,i) = mean((x1(j,1:end-i)'-x1(j,1+i:end)').^2+(y1(j,1:end-i)'-y1(j,1+i:end)').^2,1); %loops
    end
end

for jj = 1:length(bn)
    j = bn(jj);
    for i = 1:length(timeM)-1
        MSD2(jj,i) = mean((x2(j,1:end-i)'-x2(j,1+i:end)').^2+(y2(j,1:end-i)'-y2(j,1+i:end)').^2,1); %loops
    end
end

%% Tips and bases (the following is the example for calculating MSD for Figure.15)

% for jj=1:length(bn1)
%     j=bn1(jj);
%     for i=1:length(timeM)-1
%         MSD1(jj,i)=mean((x1(j,1:end-i)'-x1(j,1+i:end)').^2+(y1(j,1:end-i)'-y1(j,1+i:end)').^2,1); %loops
%     end
% end
% 
% for jj=1:length(bn2)
%     j=bn2(jj);
%     for i=1:length(timeM)-1
%         MSD2(jj,i)=mean((x1(j,1:end-i)'-x1(j,1+i:end)').^2+(y1(j,1:end-i)'-y1(j,1+i:end)').^2,1); %loops
%     end
% end

% for jj=1:length(bn1)
%     j=bn1(jj);
%     for i=1:length(timeM)-1
%         MSD3(jj,i)=mean((x2(j,1:end-i)'-x2(j,1+i:end)').^2+(y2(j,1:end-i)'-y2(j,1+i:end)').^2,1);%no loop
%     end
% end
% 
% for jj=1:length(bn2)
%     j=bn2(jj);
%     for i=1:length(timeM)-1
%         MSD4(jj,i)=mean((x2(j,1:end-i)'-x2(j,1+i:end)').^2+(y2(j,1:end-i)'-y2(j,1+i:end)').^2,1);%no loop
%     end
% end

%% Delete unused data (undelete as needed)
clear dimPhyWalk11 dimPhyWalk12 dimPhyWalk21 dimPhyWalk22
clear x1 y1 x2 y2
clear vecs Lambda SMatrix expDecayMatrix PhyWalk MMukappa Lcheck

%% Save data

% Specify path and file name
savepath = '';
fname = '';

save([savepath fname]);


end