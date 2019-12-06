% Calculate vector of inverse dynamics joint torques for
% S5PPPRR2
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [34x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5PPPRR2_invdynJ_fixb_regmin2vec.m
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPPRR2_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(34,1), zeros(13,1)}
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR2_invdynJ_fixb_mdp_slag_vr: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3) + RV(8) * MDP(5) + RV(12) * MDP(6) + RV(25) * MDP(12) + RV(30) * MDP(13); RV(3) * MDP(2) + RV(5) * MDP(3) + RV(9) * MDP(5) + RV(13) * MDP(6) + RV(26) * MDP(12) + RV(31) * MDP(13); RV(6) * MDP(3) + RV(10) * MDP(5) + RV(14) * MDP(6) + RV(27) * MDP(12) + RV(32) * MDP(13); RV(7) * MDP(4) + RV(11) * MDP(5) + RV(15) * MDP(6) + RV(16) * MDP(7) + RV(18) * MDP(8) + RV(20) * MDP(9) + RV(22) * MDP(10) + RV(28) * MDP(12) + RV(33) * MDP(13); RV(17) * MDP(7) + RV(19) * MDP(8) + RV(21) * MDP(9) + RV(23) * MDP(10) + RV(24) * MDP(11) + RV(29) * MDP(12) + RV(34) * MDP(13);];
tauJ = t1;
