function [eigenMk] = Main(traa)
% Uniform friction environment

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

%% Commonly adjustable parameters 
% (Note: mouse has bead size of 10kb, yeast of 500 bp, 
% using chromatin compaction fraction of 50bp/nm)

zeta = 1.08e-6; % friction coefficient (in Ns/m): 1.08e-6 for mouse; 5.4e-8 for yeast
det = 1; % simulation timestep: 1s for mouse, 0.02s for yeast
taup_limit = 1e4; % total simulation time, in unit of tau_p (approximately)
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

stablept = 2e4; % first 2e4 points in the LEF simulation are discarded to ensure steady state
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

% Obtain friction matrix
VMu = ones(L,1);
VMu = VMu*(1/zeta);
MMu = diag(VMu,0); % Inverse of the friction matrix


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
                    if l_sites_traj(ts(indts,2),kapi)>0
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
            
            % Diagonalize the Z^{-1}K matrix
            MMukappa = MMu*KKappa;
            [vecs,Lambda] = eig(MMukappa);
            lambda = diag(Lambda); % vector of eigenvalues

            % Make sure the center of mass does not move for cases of periodic and free b.c.
            if (bc == 0) || (bc == 1)
                lambda(1) = 0;
            end
            
            % Calculate terms in Eq.36
            for i = 1:L
                if vecs(1,i)<0
                    SMatrix(:,i) = -vecs(:,i); 
                else
                    SMatrix(:,i) = vecs(:,i); % eigenvectors
                end

                expDecay(i) = exp(-lambda(i)*dt); % exponential decay in the deterministic term
                expSD(i) = sqrt(1-expDecay(i)^2); % exponential term in the standard deviation (stochastic term)
            end
            
            % The eigenvalue matrix scaled by the spring constant
            eigenMk = SMatrix'*KKappa*SMatrix;

            % Set the zero eigenvalue to a very small positive value for convenience
            % This won't affect result since the eigenmode with zero eigenvalue is set to 0 later
            for i = 1:L
                if eigenMk(i,i) <= 0
                    eigenMk(i,i) = 1e-30;
                end
            end
            
            % Calculate terms in Eq.36
            for i=1:L
                SSSD(i) = sqrt(kb*Temp./eigenMk(i,i)); % steady-state standard deviation
                TSD(i) = SSSD(i)*expSD(i); % total standard deviation of normal modes
                % (note: the entry of zero eigenvalue index in TSD is zero by multiplying with expSD)
            end
                        
            % Obtain the initial physical positions of all beads
            if np == 2 
                % initial condition for the first normal coordinate depends on the b.c.
                if (bc == 0) || (bc == 1)
                    P1(1) = 0;
                elseif bc == 2
                    P1(1) = randn(1)*SSSD(1);
                end
                
                % initial condition for the rest of normal coordinates
                for iniPW = 2:L
                    P1(iniPW) = randn(1)*SSSD(iniPW);
                end
                
                % The initial positions of all beads
                PhyWalk(:,1) = SMatrix*P1';
            end
            
            % Update the physical positions of all beads according to the algorithm
            PhyWalk(:,np) = SMatrix*diag(expDecay)*(SMatrix\PhyWalk(:,np-1)) + SMatrix*(TSD.*randn(1,L))';
        end
        
        % save the position series for different cases: loop, no loop, physical dimensions
        if loopin == 1 && dmrouse == 1
            dimPhyWalk11=PhyWalk;
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
% bn2 = 13:25:600;
bn = 5:10:600; % 600-bead simulation (mostly used)

for jj = 1:length(bn)
    j = bn(jj);
    for i = 1:length(timeM)-1
        MSD1(jj,i) = mean((x1(j,1:end-i)'-x1(j,1+i:end)').^2+(y1(j,1:end-i)'-y1(j,1+i:end)').^2,1);
    end
end

for jj = 1:length(bn)
    j = bn(jj);
    for i = 1:length(timeM)-1
        MSD2(jj,i) = mean((x2(j,1:end-i)'-x2(j,1+i:end)').^2+(y2(j,1:end-i)'-y2(j,1+i:end)').^2,1);
    end
end

%% Delete unused data (undelete as needed)
clear dimPhyWalk11 dimPhyWalk12 dimPhyWalk21 dimPhyWalk22
clear NormMod PhyWalk
clear MMukappa vecs Lambda SMatrix SS
clear x1 x2 y1 y2

%% Save data

% Specify path and file name
savepath = '';
fname = ''; 

save([savepath fname]);

end