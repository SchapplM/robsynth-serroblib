% Calculate vector of inverse dynamics joint torques for
% S1P1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [1x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S1P1_invdynJ_fixb_regmin2vec.m
% MDP [1x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S1P1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S1P1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1), zeros(1,1)}
assert(isreal(MDP) && all(size(MDP) == [1 1]), ...
  'S1P1_invdynJ_fixb_mdp_slag_vr: MDP has to be [1x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:56
% EndTime: 2021-01-14 12:21:56
% DurationCPUTime: 0.02s
% Computational Cost: add. (0->0), mult. (1->1), div. (0->0), fcn. (1->1), ass. (0->1)
t1 = [RV(1) * MDP(1);];
tauJ = t1;
