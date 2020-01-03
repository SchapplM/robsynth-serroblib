% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:49
% EndTime: 2019-12-31 16:36:53
% DurationCPUTime: 1.18s
% Computational Cost: add. (631->189), mult. (1585->305), div. (0->0), fcn. (1277->10), ass. (0->111)
t60 = cos(qJ(3));
t103 = t60 * qJD(2);
t43 = -qJD(4) + t103;
t137 = t43 + qJD(4);
t57 = sin(qJ(3));
t107 = qJD(4) * t57;
t136 = -qJD(2) * t107 + qJDD(3);
t100 = t57 * qJDD(2);
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t8 = (qJD(3) * (qJD(4) + t103) + t100) * t56 - t136 * t59;
t55 = cos(pkin(4));
t115 = qJD(1) * t55;
t53 = sin(pkin(4));
t116 = qJD(1) * t53;
t58 = sin(qJ(2));
t36 = qJD(2) * pkin(6) + t116 * t58;
t20 = t115 * t57 + t60 * t36;
t61 = cos(qJ(2));
t99 = qJD(1) * qJD(2);
t23 = qJDD(2) * pkin(6) + (qJDD(1) * t58 + t61 * t99) * t53;
t101 = qJDD(1) * t55;
t92 = t60 * t101;
t2 = -qJDD(3) * pkin(3) + t20 * qJD(3) + t57 * t23 - t92;
t124 = t53 * t60;
t122 = t55 * t58;
t52 = sin(pkin(8));
t54 = cos(pkin(8));
t25 = t122 * t54 + t52 * t61;
t27 = -t122 * t52 + t54 * t61;
t125 = t53 * t58;
t28 = t125 * t57 - t55 * t60;
t69 = g(1) * (t124 * t52 - t27 * t57) + g(2) * (-t124 * t54 - t25 * t57) - g(3) * t28;
t83 = pkin(3) * t57 - pkin(7) * t60;
t135 = t43 * (pkin(7) * qJD(4) + t83 * qJD(2)) - t2 - t69;
t62 = qJD(3) ^ 2;
t123 = t53 * t61;
t90 = t58 * t99;
t80 = -qJDD(1) * t123 + t53 * t90;
t121 = t55 * t61;
t24 = -t121 * t54 + t52 * t58;
t26 = t121 * t52 + t54 * t58;
t82 = g(1) * t26 + g(2) * t24;
t134 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t62 + t53 * (-g(3) * t61 + t90) - t80 + t82;
t18 = qJD(3) * pkin(7) + t20;
t133 = (pkin(6) * t43 + t18) * qJD(4) + t82;
t104 = t59 * qJD(3);
t98 = qJD(2) * qJD(3);
t89 = t60 * t98;
t7 = qJD(4) * t104 + (t100 + t89) * t59 + t136 * t56;
t131 = t7 * t56;
t113 = qJD(2) * t57;
t31 = t113 * t56 - t104;
t130 = t31 * t43;
t105 = t56 * qJD(3);
t33 = t113 * t59 + t105;
t129 = t33 * t43;
t128 = t33 * t59;
t127 = t43 * t60;
t126 = t53 * t57;
t120 = t57 * t61;
t119 = t60 * t61;
t50 = t57 ^ 2;
t118 = -t60 ^ 2 + t50;
t117 = qJD(2) * pkin(2);
t114 = qJD(2) * t53;
t112 = qJD(3) * t33;
t111 = qJD(3) * t57;
t110 = qJD(3) * t60;
t109 = qJD(4) * t43;
t108 = qJD(4) * t56;
t102 = qJDD(1) - g(3);
t48 = t60 * qJDD(2);
t97 = t58 * t114;
t96 = t61 * t114;
t95 = t43 * t104;
t94 = t61 * t116;
t93 = t57 * t101;
t87 = t57 * t98;
t37 = -t94 - t117;
t85 = -qJD(2) * t37 - t23;
t81 = g(1) * t27 + g(2) * t25;
t38 = -t60 * pkin(3) - t57 * pkin(7) - pkin(2);
t21 = qJD(2) * t38 - t94;
t4 = t59 * t18 + t56 * t21;
t79 = t56 * t18 - t59 * t21;
t30 = qJDD(4) - t48 + t87;
t35 = t83 * qJD(3);
t78 = -t38 * t30 + t35 * t43;
t63 = qJD(2) ^ 2;
t77 = qJDD(2) * t61 - t58 * t63;
t29 = t124 * t58 + t55 * t57;
t15 = -t123 * t59 - t29 * t56;
t75 = t123 * t56 - t29 * t59;
t74 = t119 * t59 + t56 * t58;
t73 = -t119 * t56 + t58 * t59;
t19 = t115 * t60 - t57 * t36;
t72 = -t109 * t59 + t56 * t30;
t71 = t108 * t43 + t59 * t30;
t1 = qJDD(3) * pkin(7) + t19 * qJD(3) + t60 * t23 + t93;
t17 = -qJD(3) * pkin(3) - t19;
t67 = -pkin(6) * t30 + qJD(3) * t17 + qJD(4) * t21 + t1;
t65 = -pkin(7) * t30 + (-t17 - t19) * t43;
t64 = -pkin(6) * qJDD(3) + (t37 + t94 - t117) * qJD(3);
t14 = qJD(3) * t29 + t57 * t96;
t13 = -qJD(3) * t28 + t60 * t96;
t12 = t126 * t52 + t27 * t60;
t10 = -t126 * t54 + t25 * t60;
t6 = qJD(2) * t35 + qJDD(2) * t38 + t80;
t5 = t59 * t6;
t3 = [t102, 0, t77 * t53, (-qJDD(2) * t58 - t61 * t63) * t53, 0, 0, 0, 0, 0, -t14 * qJD(3) - t28 * qJDD(3) + (t60 * t77 - t61 * t87) * t53, -t13 * qJD(3) - t29 * qJDD(3) + (-t57 * t77 - t61 * t89) * t53, 0, 0, 0, 0, 0, -(qJD(4) * t75 - t13 * t56 + t59 * t97) * t43 + t15 * t30 + t14 * t31 + t28 * t8, (qJD(4) * t15 + t13 * t59 + t56 * t97) * t43 + t75 * t30 + t14 * t33 + t28 * t7; 0, qJDD(2), t102 * t123 + t82, -t102 * t125 + t81, t50 * qJDD(2) + 0.2e1 * t60 * t87, -0.2e1 * t118 * t98 + 0.2e1 * t48 * t57, qJDD(3) * t57 + t62 * t60, qJDD(3) * t60 - t62 * t57, 0, t134 * t60 + t64 * t57, -t134 * t57 + t64 * t60, t7 * t59 * t57 + (t104 * t60 - t107 * t56) * t33, (-t31 * t59 - t33 * t56) * t110 + (-t131 - t59 * t8 + (t31 * t56 - t128) * qJD(4)) * t57, (-t7 - t95) * t60 + (t71 + t112) * t57, (t105 * t43 + t8) * t60 + (-qJD(3) * t31 - t72) * t57, -t111 * t43 - t30 * t60, -t79 * t111 - t5 * t60 + (t110 * t31 + t57 * t8) * pkin(6) + (t17 * t107 + t133 * t60 - t78) * t59 + (-(pkin(6) * t111 - qJD(4) * t38) * t43 + t2 * t57 + t67 * t60 - t81) * t56 + (-g(3) * t74 + (-t120 * t31 + t43 * t73) * qJD(1)) * t53, t78 * t56 + (t109 * t38 - t81) * t59 + (pkin(6) * t112 + t67 * t59 + (-t133 + t6) * t56) * t60 + (-t17 * t108 - t4 * qJD(3) + t2 * t59 + (t7 - t95) * pkin(6)) * t57 + (-g(3) * t73 + (-t120 * t33 - t43 * t74) * qJD(1)) * t53; 0, 0, 0, 0, -t57 * t63 * t60, t118 * t63, t100, t48, qJDD(3), t57 * t85 - t69 + t92, g(1) * t12 + g(2) * t10 + g(3) * t29 + t60 * t85 - t93, -t128 * t43 + t131, (t7 + t130) * t59 + (-t8 + t129) * t56, (t127 * t59 - t57 * t33) * qJD(2) + t72, (-t127 * t56 + t57 * t31) * qJD(2) + t71, t43 * t113, -pkin(3) * t8 + t113 * t79 + t135 * t59 - t20 * t31 + t65 * t56, -pkin(3) * t7 + t4 * t113 - t135 * t56 - t20 * t33 + t65 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t31, -t31 ^ 2 + t33 ^ 2, t7 - t130, -t129 - t8, t30, -t56 * t1 + t5 - t17 * t33 - g(1) * (-t12 * t56 + t26 * t59) - g(2) * (-t10 * t56 + t24 * t59) - g(3) * t15 - t137 * t4, -t59 * t1 - t56 * t6 + t17 * t31 - g(1) * (-t12 * t59 - t26 * t56) - g(2) * (-t10 * t59 - t24 * t56) - g(3) * t75 + t137 * t79;];
tau_reg = t3;
