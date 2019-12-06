% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:19
% EndTime: 2019-12-05 15:11:21
% DurationCPUTime: 0.78s
% Computational Cost: add. (688->166), mult. (1592->228), div. (0->0), fcn. (1228->8), ass. (0->113)
t57 = cos(pkin(8));
t121 = qJD(1) * t57;
t60 = sin(qJ(3));
t55 = sin(pkin(8));
t122 = qJD(1) * t55;
t62 = cos(qJ(3));
t97 = t62 * t122;
t35 = qJD(2) * t60 + t97;
t31 = qJD(3) * pkin(6) + t35;
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t20 = -t59 * t121 + t31 * t61;
t148 = qJD(4) * t20;
t77 = t61 * t121 + t31 * t59;
t147 = qJD(5) + t77;
t9 = -qJD(4) * pkin(4) + t147;
t115 = qJD(3) * t61;
t79 = pkin(4) * t61 + qJ(5) * t59 + pkin(3);
t146 = t79 * qJD(3);
t143 = t62 * qJD(2) - t60 * t122;
t145 = qJD(3) * t143;
t144 = t79 * qJDD(3);
t12 = qJD(4) * qJ(5) + t20;
t136 = t55 * t60;
t101 = g(3) * t136;
t142 = -qJD(3) * t35 - t101;
t58 = cos(pkin(7));
t130 = t58 * t62;
t131 = t58 * t60;
t56 = sin(pkin(7));
t132 = t56 * t62;
t133 = t56 * t60;
t86 = -g(2) * (-t57 * t133 - t130) - g(1) * (-t57 * t131 + t132);
t113 = qJDD(4) * pkin(4);
t141 = qJDD(5) - t113;
t137 = t55 * t59;
t135 = t55 * t61;
t134 = t55 * t62;
t64 = qJD(3) ^ 2;
t129 = t62 * t64;
t82 = pkin(4) * t59 - qJ(5) * t61;
t25 = t82 * qJD(4) - t59 * qJD(5);
t128 = t25 - t35;
t53 = t59 ^ 2;
t54 = t61 ^ 2;
t127 = t53 - t54;
t126 = t53 + t54;
t63 = qJD(4) ^ 2;
t125 = t63 + t64;
t124 = qJD(3) * pkin(3);
t123 = pkin(6) * qJDD(4);
t21 = -t143 - t146;
t120 = qJD(3) * t21;
t119 = qJD(3) * t25;
t117 = qJD(3) * t59;
t116 = qJD(3) * t60;
t111 = qJDD(1) - g(3);
t109 = qJDD(1) * t55;
t108 = qJDD(1) * t57;
t107 = qJDD(4) * t59;
t106 = t60 * qJDD(2);
t105 = t61 * qJDD(3);
t104 = t62 * qJDD(3);
t103 = qJD(3) * qJD(4);
t102 = qJDD(4) * qJ(5);
t100 = t61 * t134;
t99 = t59 * t64 * t61;
t98 = qJD(4) * t143 * t59 + t61 * t101 + t35 * t115;
t95 = t55 * t116;
t94 = -g(1) * t56 + g(2) * t58;
t93 = t59 * t103;
t92 = t61 * t103;
t11 = qJDD(3) * pkin(6) + t62 * t109 + t106 + t145;
t91 = t61 * t108 + t59 * t11 + t148;
t30 = -t143 - t124;
t90 = t30 - t124;
t89 = t21 - t146;
t88 = t126 * qJDD(3);
t27 = t57 * t132 - t131;
t29 = t57 * t130 + t133;
t85 = g(1) * t29 + g(2) * t27;
t84 = -t59 * t108 + t61 * t11;
t83 = t12 * t61 + t59 * t9;
t81 = -t60 * t64 + t104;
t80 = -qJDD(3) * t60 - t129;
t78 = qJD(2) * t116 + qJD(3) * t97 - qJDD(2) * t62 + t60 * t109;
t32 = t59 * t134 + t57 * t61;
t75 = t80 * t55;
t74 = -pkin(6) * t63 + t86;
t13 = -t56 * t135 + t27 * t59;
t15 = -t58 * t135 + t29 * t59;
t72 = g(1) * t15 + g(2) * t13 + g(3) * t32 - t91;
t71 = 0.2e1 * qJDD(3) * pkin(3) + t74 - t78;
t4 = t119 + t78 - t144;
t70 = -t4 + t74 + t144;
t14 = t56 * t137 + t27 * t61;
t16 = t58 * t137 + t29 * t61;
t33 = -t57 * t59 + t100;
t69 = g(1) * t16 + g(2) * t14 + g(3) * t33 - t84;
t2 = t102 + (qJD(5) - t77) * qJD(4) + t84;
t3 = t91 + t141;
t68 = t2 * t61 + t3 * t59 + (-t12 * t59 + t61 * t9) * qJD(4);
t67 = t72 + t148;
t17 = -qJD(4) * t100 + (qJD(4) * t57 + t95) * t59;
t66 = qJD(4) * t17 - qJDD(4) * t32 + t93 * t136 + t61 * t75;
t65 = -g(3) * t134 + t68 - t85;
t51 = t59 * qJDD(3);
t36 = t82 * qJD(3);
t18 = -qJD(4) * t32 - t61 * t95;
t6 = (-0.2e1 * t93 + t105) * t62 + (-t125 * t61 - t107) * t60;
t5 = (-qJDD(4) * t60 - 0.2e1 * t62 * t103) * t61 + (t125 * t60 - t104) * t59;
t1 = qJD(4) * t18 + qJDD(4) * t33 + (t80 * t59 - t60 * t92) * t55;
t7 = [t111, -g(3) + (t55 ^ 2 + t57 ^ 2) * qJDD(1), 0, t75, -t81 * t55, 0, 0, 0, 0, 0, t66, -t1, t66, (t32 * t59 + t33 * t61) * qJDD(3) + (-t17 * t59 + t18 * t61 + (t32 * t61 - t33 * t59) * qJD(4)) * qJD(3), t1, t12 * t18 - t17 * t9 + t2 * t33 + t3 * t32 - g(3) + (t62 * t120 + t4 * t60) * t55; 0, qJDD(2) + t94, 0, t81, t80, 0, 0, 0, 0, 0, t6, t5, t6, t126 * t129 + t60 * t88, -t5, (t83 * qJD(3) - t4) * t62 + (t68 + t120) * t60 + t94; 0, 0, qJDD(3), -t78 + t86 - t142, -t111 * t134 - t106 + t85, qJDD(3) * t53 + 0.2e1 * t59 * t92, -0.2e1 * t127 * t103 + 0.2e1 * t59 * t105, t61 * t63 + t107, qJDD(4) * t61 - t59 * t63, 0, (t90 * qJD(4) - t123) * t59 + t71 * t61 + t98, (-t123 + (t143 + t90) * qJD(4)) * t61 + (-t71 + t142) * t59, (t89 * qJD(4) - t123) * t59 + (t70 - t119) * t61 + t98, pkin(6) * t88 - t126 * t145 + t65, (t123 + (-t143 - t89) * qJD(4)) * t61 + (-t128 * qJD(3) + t101 + t70) * t59, t65 * pkin(6) + t128 * t21 - t83 * t143 + (-t4 + t101 + t86) * t79; 0, 0, 0, 0, 0, -t99, t127 * t64, t51, t105, qJDD(4), -t30 * t117 + t67, -t30 * t115 + t69, 0.2e1 * t113 - qJDD(5) + (-t21 * t59 + t36 * t61) * qJD(3) + t67, -t82 * qJDD(3), 0.2e1 * t102 + (t21 * t61 + t36 * t59) * qJD(3) + 0.2e1 * qJD(4) * qJD(5) - t69, t2 * qJ(5) - t3 * pkin(4) - t21 * t36 - t9 * t20 - g(1) * (-pkin(4) * t15 + qJ(5) * t16) - g(2) * (-pkin(4) * t13 + qJ(5) * t14) - g(3) * (-pkin(4) * t32 + qJ(5) * t33) + t147 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t99, t51, -t53 * t64 - t63, -qJD(4) * t12 + t21 * t117 + t141 - t72;];
tau_reg = t7;
