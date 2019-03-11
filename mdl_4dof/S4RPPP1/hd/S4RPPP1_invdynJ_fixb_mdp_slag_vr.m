% Calculate vector of inverse dynamics joint torques for
% S4RPPP1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [39x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S4RPPP1_invdynJ_fixb_regmin2vec.m
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPPP1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(39,1), zeros(15,1)}
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPP1_invdynJ_fixb_mdp_slag_vr: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8) + RV(15) * MDP(9) + RV(18) * MDP(10) + RV(21) * MDP(11) + RV(24) * MDP(12) + RV(28) * MDP(13) + RV(32) * MDP(14) + RV(36) * MDP(15); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(22) * MDP(11) + RV(25) * MDP(12) + RV(29) * MDP(13) + RV(33) * MDP(14) + RV(37) * MDP(15); RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(23) * MDP(11) + RV(26) * MDP(12) + RV(30) * MDP(13) + RV(34) * MDP(14) + RV(38) * MDP(15); RV(27) * MDP(12) + RV(31) * MDP(13) + RV(35) * MDP(14) + RV(39) * MDP(15);];
tauJ  = t1;
