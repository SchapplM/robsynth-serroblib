% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (82->30), mult. (249->63), div. (0->0), fcn. (164->6), ass. (0->25)
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t64 = (t48 ^ 2 - t50 ^ 2) * MDP(7);
t46 = sin(pkin(7));
t47 = cos(pkin(7));
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t42 = t49 * t46 - t51 * t47;
t38 = t42 * qJD(1);
t43 = t51 * t46 + t49 * t47;
t52 = qJD(4) ^ 2;
t63 = t43 * t52;
t61 = t50 * MDP(6);
t60 = qJD(4) * t48;
t59 = qJD(4) * t50;
t58 = t50 * MDP(12);
t41 = t43 * qJD(3);
t57 = pkin(5) * t52 + qJD(1) * t41;
t34 = -qJD(3) * pkin(3) + t38;
t56 = qJD(4) * (t34 - t38);
t40 = t42 * qJD(3);
t55 = qJD(1) * t40 - t34 * qJD(3);
t54 = -t48 * MDP(11) - t58;
t53 = qJD(3) ^ 2;
t1 = [(t40 * t60 - t50 * t63) * MDP(11) + (t40 * t59 + t48 * t63) * MDP(12) + (-t41 * MDP(4) + t40 * MDP(5) + (-t41 * t50 + t42 * t60) * MDP(11) + (t41 * t48 + t42 * t59) * MDP(12)) * qJD(3); t54 * t52; (-t57 * MDP(11) + MDP(12) * t56 + t52 * MDP(8)) * t50 + (MDP(11) * t56 + t57 * MDP(12) - t52 * MDP(9)) * t48 + ((t50 * MDP(11) - t48 * MDP(12)) * t43 * qJD(1) + (t54 * pkin(3) + 0.2e1 * t48 * t61 - 0.2e1 * t64) * qJD(4)) * qJD(3); t53 * t64 + t55 * t58 + (t55 * MDP(11) - t53 * t61) * t48;];
tauc = t1;
