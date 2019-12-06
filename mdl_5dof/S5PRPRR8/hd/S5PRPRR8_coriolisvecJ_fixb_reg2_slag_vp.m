% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:34
% DurationCPUTime: 1.49s
% Computational Cost: add. (1494->228), mult. (3443->341), div. (0->0), fcn. (2325->8), ass. (0->134)
t57 = cos(qJ(4));
t111 = qJD(4) * t57;
t54 = sin(qJ(4));
t115 = qJD(1) * t54;
t59 = -pkin(2) - pkin(7);
t51 = sin(pkin(5));
t116 = qJD(1) * t51;
t58 = cos(qJ(2));
t88 = t58 * t116;
t75 = qJD(3) - t88;
t32 = t59 * qJD(2) + t75;
t52 = cos(pkin(5));
t114 = qJD(2) * t51;
t55 = sin(qJ(2));
t95 = t55 * t114;
t11 = t32 * t111 + (-qJD(4) * t52 + t95) * t115;
t79 = pkin(4) * t57 + pkin(8) * t54;
t33 = t79 * qJD(4) + qJD(3);
t22 = (t33 + t88) * qJD(2);
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t125 = t52 * t57;
t21 = qJD(1) * t125 + t54 * t32;
t15 = qJD(4) * pkin(8) + t21;
t42 = t54 * pkin(4) - t57 * pkin(8) + qJ(3);
t89 = t55 * t116;
t27 = t42 * qJD(2) + t89;
t72 = t53 * t15 - t56 * t27;
t1 = -t72 * qJD(5) + t56 * t11 + t53 * t22;
t103 = t54 * qJD(2);
t46 = qJD(5) + t103;
t151 = t72 * t46 + t1;
t6 = t56 * t15 + t53 * t27;
t2 = -qJD(5) * t6 - t53 * t11 + t56 * t22;
t150 = -t6 * t46 - t2;
t104 = t53 * qJD(4);
t113 = qJD(2) * t57;
t87 = t56 * t113;
t38 = t87 + t104;
t109 = qJD(5) * t38;
t93 = t54 * t104;
t24 = -qJD(2) * t93 + t109;
t102 = t56 * qJD(4);
t106 = qJD(5) * t57;
t91 = t53 * t106;
t66 = t54 * t102 + t91;
t23 = t66 * qJD(2) - qJD(5) * t102;
t81 = qJD(2) * t89;
t41 = t57 * t81;
t12 = t21 * qJD(4) - t41;
t142 = t12 * t57;
t67 = t52 * t115 - t57 * t32;
t62 = -(-t21 * t57 - t54 * t67) * qJD(4) + t11 * t54 - t142;
t124 = t53 * t55;
t123 = t54 * t59;
t29 = t56 * t123 + t53 * t42;
t110 = qJD(4) * t59;
t92 = t57 * t110;
t147 = t29 * qJD(5) - t56 * t33 + t53 * t92 + (-t54 * t124 + t56 * t58) * t116;
t122 = t55 * t56;
t28 = -t53 * t123 + t56 * t42;
t146 = -t28 * qJD(5) - t53 * t33 - t56 * t92 + (t54 * t122 + t53 * t58) * t116;
t127 = t51 * t58;
t30 = t57 * t127 + t52 * t54;
t145 = t12 * t30;
t144 = t12 * t53;
t143 = t12 * t56;
t14 = -qJD(4) * pkin(4) + t67;
t141 = t14 * t53;
t140 = t14 * t56;
t139 = t23 * t53;
t138 = t23 * t54;
t137 = t24 * t54;
t136 = t24 * t56;
t34 = (qJD(3) + t88) * qJD(2);
t135 = t34 * t55;
t36 = t53 * t113 - t102;
t134 = t36 * t46;
t133 = t38 * t36;
t132 = t38 * t46;
t101 = qJD(2) * qJ(3);
t40 = t89 + t101;
t131 = t40 * t58;
t130 = t46 * t53;
t129 = t46 * t54;
t128 = t46 * t56;
t61 = qJD(2) ^ 2;
t126 = t51 * t61;
t121 = t57 * t23;
t120 = t57 * t24;
t49 = t54 ^ 2;
t50 = t57 ^ 2;
t119 = t49 - t50;
t60 = qJD(4) ^ 2;
t118 = -t60 - t61;
t117 = qJD(2) * pkin(2);
t112 = qJD(4) * t54;
t108 = qJD(5) * t53;
t107 = qJD(5) * t56;
t105 = t40 * qJD(2);
t100 = qJD(2) * qJD(4);
t99 = t54 * t127;
t98 = t55 * t126;
t97 = t58 * t126;
t96 = t57 * t61 * t54;
t94 = t58 * t114;
t90 = t56 * t106;
t86 = t57 * t100;
t85 = t46 + t103;
t84 = t57 * t95;
t83 = t54 * t95;
t82 = qJD(5) * t54 + qJD(2);
t80 = t54 * t86;
t78 = -t40 + t89;
t77 = t53 * t6 - t56 * t72;
t76 = -t53 * t72 - t56 * t6;
t70 = qJD(2) * t50 - t129;
t69 = t34 * qJ(3) + t40 * qJD(3);
t31 = -t99 + t125;
t18 = t51 * t122 - t31 * t53;
t19 = t51 * t124 + t31 * t56;
t68 = -pkin(8) * t111 + t14 * t54;
t65 = t78 - t101;
t64 = t75 * qJD(2) - t59 * t60 + t34;
t63 = -t77 * qJD(5) + t1 * t56 - t2 * t53;
t39 = t79 * qJD(2);
t35 = t75 - t117;
t17 = -qJD(4) * t99 + t52 * t111 - t84;
t16 = -t30 * qJD(4) + t83;
t10 = t53 * t39 - t56 * t67;
t9 = t56 * t39 + t53 * t67;
t4 = t18 * qJD(5) + t16 * t56 + t53 * t94;
t3 = -t19 * qJD(5) - t16 * t53 + t56 * t94;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t97, (t135 + (t131 + (t35 - t88) * t55) * qJD(2)) * t51, 0, 0, 0, 0, 0, 0, t54 * t97 + (-t17 + t84) * qJD(4), t57 * t97 + (-t16 - t83) * qJD(4), (-t16 * t54 + t17 * t57 + (-t30 * t54 - t31 * t57) * qJD(4)) * qJD(2), t11 * t31 + t145 + t21 * t16 + t67 * t17 + (t105 * t58 + t135) * t51, 0, 0, 0, 0, 0, 0, t17 * t36 + t18 * t86 + t30 * t24 + t3 * t46, t17 * t38 - t19 * t86 - t30 * t23 - t4 * t46, t18 * t23 - t19 * t24 - t3 * t38 - t4 * t36, t1 * t19 + t14 * t17 + t2 * t18 - t3 * t72 + t6 * t4 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), (-t131 + (-t35 - t117) * t55) * t116 + t69, -0.2e1 * t80, 0.2e1 * t119 * t100, -t60 * t54, 0.2e1 * t80, -t60 * t57, 0, -t111 * t65 + t54 * t64, t112 * t65 + t57 * t64, (t49 + t50) * t81 - t62, (-t131 + (-t21 * t54 + t57 * t67) * t55) * t116 + t62 * t59 + t69, -t121 * t56 - t38 * t66, (t36 * t56 + t38 * t53) * t112 + (t139 - t136 + (t36 * t53 - t38 * t56) * qJD(5)) * t57, -t46 * t91 - t138 + (t38 * t57 + t56 * t70) * qJD(4), t53 * t120 + (t90 - t93) * t36, -t46 * t90 - t137 + (-t36 * t57 - t53 * t70) * qJD(4), t85 * t111, -t147 * t46 + (t2 + (t36 * t59 - t141) * qJD(4)) * t54 + (t36 * t89 + t14 * t107 + t144 - t24 * t59 + (qJD(2) * t28 - t72) * qJD(4)) * t57, t146 * t46 + (-t1 + (t38 * t59 - t140) * qJD(4)) * t54 + (t38 * t89 - t14 * t108 + t143 + t23 * t59 + (-qJD(2) * t29 - t6) * qJD(4)) * t57, t28 * t23 - t29 * t24 + t147 * t38 + t146 * t36 + t77 * t112 + (qJD(5) * t76 - t1 * t53 - t2 * t56) * t57, -t59 * t142 + t1 * t29 + t2 * t28 - t146 * t6 + t147 * t72 + (t110 * t54 + t57 * t89) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t78 * qJD(2), 0, 0, 0, 0, 0, 0, t118 * t54, t118 * t57, 0, t62 - t105, 0, 0, 0, 0, 0, 0, -t120 - t82 * t128 + (-t53 * t57 * t85 + t36 * t54) * qJD(4), t121 + t82 * t130 + (-t57 * t128 + (t38 - t87) * t54) * qJD(4), (-t111 * t36 + t38 * t82 - t137) * t56 + (t111 * t38 + t36 * t82 - t138) * t53, -t77 * qJD(2) + (-qJD(4) * t76 - t12) * t57 + (qJD(4) * t14 + t63) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t119 * t61, 0, -t96, 0, 0, -t105 * t57 + t41, -t78 * t103, 0, 0, t128 * t38 - t139, (-t23 - t134) * t56 + (-t24 - t132) * t53, t46 * t107 + (t54 * t128 + (-t38 + t104) * t57) * qJD(2), t130 * t36 - t136, -t46 * t108 + (-t53 * t129 + (t36 + t102) * t57) * qJD(2), -t46 * t113, -pkin(4) * t24 - t143 - t21 * t36 - t9 * t46 + (-pkin(8) * t128 + t141) * qJD(5) + (t53 * t68 + t57 * t72) * qJD(2), pkin(4) * t23 + t10 * t46 + t144 - t21 * t38 + (pkin(8) * t130 + t140) * qJD(5) + (t56 * t68 + t57 * t6) * qJD(2), t10 * t36 + t9 * t38 + ((-t24 + t109) * pkin(8) + t151) * t56 + ((qJD(5) * t36 - t23) * pkin(8) + t150) * t53, -t12 * pkin(4) + pkin(8) * t63 - t6 * t10 - t14 * t21 + t72 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, -t36 ^ 2 + t38 ^ 2, t134 - t23, -t133, t132 - t24, t86, -t14 * t38 - t150, t14 * t36 - t151, 0, 0;];
tauc_reg = t5;
