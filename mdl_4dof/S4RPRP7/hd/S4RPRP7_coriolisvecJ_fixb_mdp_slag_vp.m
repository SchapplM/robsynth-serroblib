% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:18
% EndTime: 2019-12-31 16:47:19
% DurationCPUTime: 0.28s
% Computational Cost: add. (189->67), mult. (384->96), div. (0->0), fcn. (130->2), ass. (0->36)
t60 = cos(qJ(3));
t58 = t60 ^ 2;
t59 = sin(qJ(3));
t87 = (t59 ^ 2 - t58) * MDP(8);
t53 = pkin(3) * t59 - qJ(4) * t60 + qJ(2);
t50 = qJD(1) * t53;
t86 = t50 * MDP(14);
t61 = -pkin(1) - pkin(5);
t85 = -t61 * MDP(17) + MDP(15);
t84 = 0.2e1 * qJD(1);
t83 = 2 * qJD(3);
t71 = qJD(3) * qJ(4);
t56 = qJD(1) * t61 + qJD(2);
t81 = t59 * t56;
t51 = t71 + t81;
t82 = t51 * t60;
t80 = t60 * t56;
t62 = qJD(3) ^ 2;
t79 = t61 * t62;
t76 = qJD(3) * pkin(3);
t75 = t60 * MDP(7);
t66 = pkin(3) * t60 + qJ(4) * t59;
t47 = qJD(3) * t66 - t60 * qJD(4) + qJD(2);
t46 = qJD(1) * t47;
t74 = t50 * MDP(16);
t73 = t50 * qJD(1);
t68 = -qJD(4) + t76;
t49 = -t68 - t80;
t70 = -t49 + t80;
t69 = MDP(6) * qJ(2) + MDP(5);
t65 = qJD(2) * t84 - t79;
t64 = 0.2e1 * t46 - t79;
t63 = qJD(1) ^ 2;
t52 = t66 * qJD(1);
t48 = (qJD(4) + t80) * qJD(3);
t1 = [(t46 * t53 + t50 * t47) * MDP(17) + (0.2e1 * qJD(2) * t69 + t83 * t87) * qJD(1) + (-(t62 * MDP(10)) + t65 * MDP(13) - t64 * MDP(16) + (MDP(12) * qJ(2) * t84 - t51 * t85 + 0.2e1 * t86) * qJD(3)) * t60 + (-(t62 * MDP(9)) + t65 * MDP(12) + t64 * MDP(14) - t85 * t48 + (t74 + (-0.2e1 * MDP(13) * qJ(2) + MDP(16) * t53 - 0.2e1 * t75) * qJD(1) + t85 * t70) * qJD(3)) * t59; -t69 * t63 + ((-MDP(13) + MDP(16)) * t60 - (MDP(12) + MDP(14)) * t59) * (t62 + t63) + (t48 * t59 - t73 + (-t59 * t70 + t82) * qJD(3)) * MDP(17); qJD(4) * MDP(16) * t83 + (t59 * t75 - t87 + (-MDP(12) * t60 + MDP(13) * t59) * qJ(2)) * t63 + ((-t86 + (t51 - t71) * MDP(15) + t52 * MDP(16)) * t60 + (-t52 * MDP(14) + (t49 + t68) * MDP(15) - t74) * t59) * qJD(1) + (t48 * qJ(4) + t51 * qJD(4) - t50 * t52 + (-t82 + (-t49 - t76) * t59) * t56) * MDP(17); t60 * t63 * t59 * MDP(14) + (-t58 * t63 - t62) * MDP(16) + (t60 * t73 + (-t51 + t81) * qJD(3)) * MDP(17);];
tauc = t1;
