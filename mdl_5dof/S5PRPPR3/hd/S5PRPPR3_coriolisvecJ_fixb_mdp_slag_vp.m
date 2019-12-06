% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:57
% DurationCPUTime: 0.24s
% Computational Cost: add. (187->69), mult. (470->105), div. (0->0), fcn. (304->6), ass. (0->46)
t79 = cos(qJ(2));
t92 = qJD(1) * t79;
t67 = qJD(2) * pkin(2) + t92;
t74 = sin(pkin(8));
t75 = cos(pkin(8));
t77 = sin(qJ(2));
t93 = qJD(1) * t77;
t57 = t67 * t74 + t75 * t93;
t54 = qJD(2) * qJ(4) + t57;
t65 = t74 * t79 + t75 * t77;
t61 = t65 * qJD(1);
t104 = -t54 + t61;
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t103 = (t76 ^ 2 - t78 ^ 2) * MDP(10);
t60 = t65 * qJD(2);
t58 = qJD(1) * t60;
t100 = t75 * t79;
t101 = t74 * t77;
t64 = -t100 + t101;
t102 = t58 * t64;
t80 = qJD(5) ^ 2;
t99 = t76 * t80;
t98 = t78 * t80;
t95 = MDP(9) * t78;
t94 = MDP(14) * t78;
t91 = qJD(5) * t76;
t90 = qJD(5) * t78;
t68 = t74 * t93;
t87 = t75 * t92;
t63 = -t68 + t87;
t89 = qJD(4) - t63;
t86 = qJD(2) * t100;
t85 = -pkin(2) * t75 - pkin(3);
t56 = t67 * t75 - t68;
t84 = qJD(2) * t54 - t58;
t71 = pkin(2) * t74 + qJ(4);
t83 = qJD(2) * t71 - t104;
t66 = qJD(2) * t68;
t55 = -t66 + (qJD(4) + t87) * qJD(2);
t82 = t89 * qJD(2) - (-pkin(6) + t85) * t80 + t55;
t81 = qJD(2) ^ 2;
t62 = -qJD(2) * t101 + t86;
t59 = qJD(1) * t86 - t66;
t53 = -qJD(2) * pkin(3) + qJD(4) - t56;
t1 = [(-t56 * t60 + t57 * t62 + t59 * t65 + t102) * MDP(5) + (t53 * t60 + t54 * t62 + t55 * t65 + t102) * MDP(8) + (t60 * t90 - t64 * t99) * MDP(14) + (-t60 * t91 - t64 * t98) * MDP(15) + (-MDP(3) * t77 - MDP(4) * t79) * t81 + (t60 * MDP(6) + t62 * MDP(7) + (t62 * t76 + t65 * t90) * MDP(14) + (t62 * t78 - t65 * t91) * MDP(15)) * qJD(2); (t56 * t61 - t57 * t63 + (-t58 * t75 + t59 * t74) * pkin(2)) * MDP(5) + (-t66 + (0.2e1 * qJD(4) - t63 + t87) * qJD(2)) * MDP(7) + (-t53 * t61 + t89 * t54 + t55 * t71 + t58 * t85) * MDP(8) - MDP(11) * t99 - MDP(12) * t98 + (t82 * t76 + t83 * t90) * MDP(14) + (t82 * t78 - t83 * t91) * MDP(15) + 0.2e1 * (-t95 * t76 + t103) * qJD(2) * qJD(5); (MDP(15) * t76 - t94) * t80; -t81 * MDP(7) + t104 * MDP(8) * qJD(2) + (t76 * MDP(14) + t78 * MDP(15)) * (-t80 - t81); -t81 * t103 - t84 * t94 + (t84 * MDP(15) + t81 * t95) * t76;];
tauc = t1;
