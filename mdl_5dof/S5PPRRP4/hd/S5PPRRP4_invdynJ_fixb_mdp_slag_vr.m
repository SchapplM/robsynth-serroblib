% Calculate vector of inverse dynamics joint torques for
% S5PPRRP4
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [34x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5PPRRP4_invdynJ_fixb_regmin2vec.m
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRP4_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(34,1), zeros(14,1)}
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vr: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(18) * MDP(11) + RV(22) * MDP(12) + RV(30) * MDP(14); RV(3) * MDP(2) + RV(5) * MDP(4) + RV(7) * MDP(5) + RV(19) * MDP(11) + RV(23) * MDP(12) + RV(26) * MDP(13) + RV(31) * MDP(14); RV(4) * MDP(3) + RV(6) * MDP(4) + RV(8) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(13) * MDP(8) + RV(15) * MDP(9) + RV(20) * MDP(11) + RV(24) * MDP(12) + RV(27) * MDP(13) + RV(32) * MDP(14); RV(10) * MDP(6) + RV(12) * MDP(7) + RV(14) * MDP(8) + RV(16) * MDP(9) + RV(17) * MDP(10) + RV(21) * MDP(11) + RV(25) * MDP(12) + RV(28) * MDP(13) + RV(33) * MDP(14); RV(29) * MDP(13) + RV(34) * MDP(14);];
tauJ = t1;
