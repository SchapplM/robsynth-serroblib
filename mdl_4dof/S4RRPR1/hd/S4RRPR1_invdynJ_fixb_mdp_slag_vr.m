% Calculate vector of inverse dynamics joint torques for
% S4RRPR1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [21x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S4RRPR1_invdynJ_fixb_regmin2vec.m
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(21,1), zeros(10,1)}
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPR1_invdynJ_fixb_mdp_slag_vr: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10); RV(12) * MDP(7); RV(15) * MDP(8) + RV(18) * MDP(9) + RV(21) * MDP(10);];
tauJ  = t1;
