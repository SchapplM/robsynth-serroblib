% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:08
% EndTime: 2019-12-31 16:35:09
% DurationCPUTime: 0.63s
% Computational Cost: add. (470->145), mult. (1060->217), div. (0->0), fcn. (785->10), ass. (0->95)
t60 = cos(qJ(3));
t88 = qJD(2) * qJD(3);
t82 = t60 * t88;
t57 = sin(qJ(3));
t92 = t57 * qJDD(2);
t121 = t82 + t92;
t61 = cos(qJ(2));
t113 = g(3) * t61;
t58 = sin(qJ(2));
t54 = sin(pkin(7));
t55 = cos(pkin(7));
t80 = g(1) * t55 + g(2) * t54;
t74 = t80 * t58;
t120 = t74 - t113;
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t27 = t56 * t60 + t59 * t57;
t50 = qJD(3) + qJD(4);
t118 = t27 * t50;
t119 = qJD(2) * t118;
t91 = t60 * qJDD(2);
t5 = t56 * t92 - t59 * t91 + t119;
t26 = t56 * t57 - t59 * t60;
t71 = t26 * t50;
t100 = qJD(3) * t60;
t117 = -qJD(4) * t60 - t100;
t89 = qJD(1) * qJD(2);
t96 = qJDD(2) * pkin(2);
t23 = -t61 * qJDD(1) + t58 * t89 - t96;
t62 = qJD(3) ^ 2;
t116 = -pkin(5) * t62 + (t80 + t89) * t58 - t113 - t23 + t96;
t115 = pkin(5) + pkin(6);
t114 = g(3) * t58;
t102 = qJD(2) * t60;
t85 = t59 * t102;
t103 = qJD(2) * t57;
t86 = t56 * t103;
t20 = -t85 + t86;
t22 = -t56 * t102 - t59 * t103;
t112 = t22 * t20;
t49 = qJDD(3) + qJDD(4);
t111 = t26 * t49;
t110 = t27 * t49;
t109 = t54 * t61;
t108 = t55 * t61;
t36 = qJD(2) * pkin(5) + t58 * qJD(1);
t81 = pkin(6) * qJD(2) + t36;
t19 = t81 * t60;
t107 = t59 * t19;
t51 = t57 ^ 2;
t106 = -t60 ^ 2 + t51;
t63 = qJD(2) ^ 2;
t105 = t62 + t63;
t104 = qJD(2) * pkin(2);
t101 = qJD(3) * t57;
t99 = qJD(4) * t56;
t98 = qJD(4) * t59;
t95 = t61 * qJD(1);
t94 = qJDD(1) - g(3);
t93 = qJDD(3) * t57;
t90 = t61 * qJDD(2);
t87 = pkin(3) * t101;
t46 = -t60 * pkin(3) - pkin(2);
t84 = qJD(3) * t115;
t83 = t57 * t88;
t79 = g(1) * t54 - g(2) * t55;
t18 = t81 * t57;
t13 = qJD(3) * pkin(3) - t18;
t77 = -t56 * t13 - t107;
t30 = t115 * t57;
t31 = t115 * t60;
t76 = -t59 * t30 - t56 * t31;
t75 = -t56 * t30 + t59 * t31;
t73 = t80 * t61;
t70 = t83 - t91;
t4 = qJD(4) * t85 + t121 * t59 - t50 * t86 + t56 * t91;
t37 = -t95 - t104;
t67 = -pkin(5) * qJDD(3) + (t37 + t95 - t104) * qJD(3);
t24 = qJDD(2) * pkin(5) + t58 * qJDD(1) + t61 * t89;
t66 = -t37 * qJD(2) + t114 - t24 + t73;
t25 = t46 * qJD(2) - t95;
t53 = qJ(3) + qJ(4);
t47 = sin(t53);
t48 = cos(t53);
t7 = qJDD(3) * pkin(3) - pkin(6) * t121 - t36 * t100 - t57 * t24;
t65 = -g(1) * (-t48 * t108 - t54 * t47) - g(2) * (-t48 * t109 + t55 * t47) + t25 * t20 + t19 * t99 + t48 * t114 + (-t19 * t50 - t7) * t56;
t8 = -t70 * pkin(6) - t36 * t101 + t60 * t24;
t64 = -g(1) * (-t47 * t108 + t54 * t48) - g(2) * (-t47 * t109 - t55 * t48) + t77 * qJD(4) + t25 * t22 - t56 * t8 + t59 * t7 + t47 * t114;
t29 = t60 * t84;
t28 = t57 * t84;
t11 = t70 * pkin(3) + t23;
t6 = -t20 ^ 2 + t22 ^ 2;
t2 = -t22 * t50 - t5;
t1 = t20 * t50 + t4;
t3 = [t94, 0, -t63 * t58 + t90, -qJDD(2) * t58 - t63 * t61, 0, 0, 0, 0, 0, (-0.2e1 * t83 + t91) * t61 + (-t105 * t60 - t93) * t58, (-qJDD(3) * t58 - 0.2e1 * t61 * t88) * t60 + (t105 * t58 - t90) * t57, 0, 0, 0, 0, 0, (-t5 - t119) * t61 + ((t56 * t101 + t117 * t59 + t57 * t99) * t50 - t110 + qJD(2) * t20) * t58, (qJD(2) * t71 - t4) * t61 + (-(-t59 * t101 + t117 * t56 - t57 * t98) * t50 + t111 - qJD(2) * t22) * t58; 0, qJDD(2), t94 * t61 + t74, -t94 * t58 + t73, t51 * qJDD(2) + 0.2e1 * t57 * t82, -0.2e1 * t106 * t88 + 0.2e1 * t57 * t91, t62 * t60 + t93, qJDD(3) * t60 - t62 * t57, 0, t116 * t60 + t67 * t57, -t116 * t57 + t67 * t60, t22 * t71 + t4 * t27, t118 * t22 + t20 * t71 - t4 * t26 - t27 * t5, -t50 * t71 + t110, -t118 * t50 - t111, 0, (-t75 * qJD(4) + t56 * t28 - t59 * t29) * t50 + t76 * t49 + t20 * t87 + t46 * t5 + t11 * t26 + t25 * t118 + t120 * t48 + (t118 * t61 - t58 * t20) * qJD(1), -(t76 * qJD(4) - t59 * t28 - t56 * t29) * t50 - t75 * t49 - t22 * t87 + t46 * t4 + t11 * t27 - t25 * t71 - t120 * t47 + (t58 * t22 - t61 * t71) * qJD(1); 0, 0, 0, 0, -t57 * t63 * t60, t106 * t63, t92, t91, qJDD(3), t66 * t57 - t79 * t60, t79 * t57 + t66 * t60, -t112, t6, t1, t2, t49, -(t56 * t18 - t107) * t50 + (-t20 * t103 + t59 * t49 - t50 * t99) * pkin(3) + t64, (-qJD(4) * t13 - t18 * t50 - t8) * t59 + (t22 * t103 - t56 * t49 - t50 * t98) * pkin(3) + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t6, t1, t2, t49, -t77 * t50 + t64, (-t8 + (-qJD(4) + t50) * t13) * t59 + t65;];
tau_reg = t3;
