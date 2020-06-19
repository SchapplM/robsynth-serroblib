% Calculate vector of inverse dynamics joint torques for
% S2RR3
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [9x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S2RR3_invdynJ_fixb_regmin2vec.m
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S2RR3_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(9,1), zeros(6,1)}
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S2RR3_invdynJ_fixb_mdp_slag_vr: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:28
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.02s
% Computational Cost: add. (7->7), mult. (9->9), div. (0->0), fcn. (9->9), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6);];
tauJ = t1;
