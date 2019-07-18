% Calculate vector of inverse dynamics joint torques for
% S4PRRR2
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [19x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S4PRRR2_invdynJ_fixb_regmin2vec.m
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR2_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(19,1), zeros(10,1)}
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vr: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1); RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(5) * MDP(5) + RV(7) * MDP(6) + RV(9) * MDP(7) + RV(11) * MDP(8) + RV(14) * MDP(9) + RV(17) * MDP(10); RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8) + RV(15) * MDP(9) + RV(18) * MDP(10); RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10);];
tauJ  = t1;
