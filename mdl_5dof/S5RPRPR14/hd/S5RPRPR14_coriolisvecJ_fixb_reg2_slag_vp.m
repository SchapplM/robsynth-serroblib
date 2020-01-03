% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR14_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:21
% DurationCPUTime: 1.49s
% Computational Cost: add. (2676->256), mult. (5766->357), div. (0->0), fcn. (3765->6), ass. (0->144)
t92 = -pkin(1) - pkin(6);
t68 = qJD(1) * t92 + qJD(2);
t175 = -qJ(4) * qJD(1) + t68;
t140 = cos(pkin(8));
t87 = sin(pkin(8));
t89 = sin(qJ(3));
t91 = cos(qJ(3));
t101 = -t140 * t89 - t87 * t91;
t54 = t101 * qJD(1);
t170 = qJD(5) - t54;
t90 = cos(qJ(5));
t115 = t90 * t170;
t139 = qJD(1) * t89;
t121 = t87 * t139;
t117 = t140 * t91;
t110 = qJD(1) * t117;
t67 = qJD(3) * t110;
t46 = qJD(3) * t121 - t67;
t88 = sin(qJ(5));
t149 = t88 * t46;
t174 = t115 * t170 - t149;
t126 = qJ(4) * qJD(3);
t132 = t89 * qJD(4);
t136 = qJD(3) * t91;
t35 = t68 * t136 + (-t126 * t91 - t132) * qJD(1);
t118 = t140 * t35;
t130 = t91 * qJD(4);
t100 = t126 * t89 - t130;
t137 = qJD(3) * t89;
t95 = qJD(1) * t100 - t137 * t68;
t12 = t87 * t95 + t118;
t59 = t101 * qJD(3);
t47 = qJD(1) * t59;
t125 = qJD(1) * qJD(3);
t120 = t91 * t125;
t83 = qJD(1) * qJD(2);
t64 = pkin(3) * t120 + t83;
t15 = -pkin(4) * t46 - pkin(7) * t47 + t64;
t50 = t175 * t89;
t43 = t140 * t50;
t51 = t175 * t91;
t45 = qJD(3) * pkin(3) + t51;
t21 = t45 * t87 + t43;
t17 = qJD(3) * pkin(7) + t21;
t57 = t110 - t121;
t65 = pkin(3) * t139 + qJD(1) * qJ(2) + qJD(4);
t22 = -pkin(4) * t54 - pkin(7) * t57 + t65;
t6 = t17 * t90 + t22 * t88;
t2 = -qJD(5) * t6 - t88 * t12 + t15 * t90;
t173 = t170 * t6 + t2;
t108 = t17 * t88 - t22 * t90;
t1 = -t108 * qJD(5) + t90 * t12 + t88 * t15;
t172 = t108 * t170 + t1;
t169 = t57 ^ 2;
t168 = 0.2e1 * t83;
t11 = -t140 * t95 + t35 * t87;
t141 = qJ(4) - t92;
t66 = t141 * t89;
t32 = t117 * t141 - t66 * t87;
t167 = t11 * t32;
t150 = t87 * t89;
t62 = t117 - t150;
t166 = t11 * t62;
t165 = t11 * t88;
t131 = t90 * qJD(3);
t135 = qJD(5) * t88;
t18 = -qJD(5) * t131 + t135 * t57 - t47 * t90;
t164 = t18 * t88;
t148 = t88 * t47;
t38 = qJD(3) * t88 + t57 * t90;
t19 = qJD(5) * t38 + t148;
t163 = t19 * t90;
t36 = t57 * t88 - t131;
t162 = t36 * t54;
t161 = t36 * t57;
t160 = t36 * t88;
t159 = t36 * t90;
t158 = t38 * t36;
t157 = t38 * t57;
t156 = t38 * t88;
t155 = t38 * t90;
t154 = t46 * t101;
t153 = t57 * t54;
t152 = t62 * t90;
t151 = t87 * t50;
t14 = t88 * t19;
t41 = t90 * t46;
t93 = qJD(3) ^ 2;
t147 = t93 * t89;
t146 = t93 * t91;
t134 = qJD(5) * t90;
t145 = -t134 * t36 - t14;
t144 = t89 ^ 2 - t91 ^ 2;
t94 = qJD(1) ^ 2;
t143 = -t93 - t94;
t142 = t94 * qJ(2);
t76 = pkin(3) * t89 + qJ(2);
t138 = qJD(1) * t91;
t56 = -qJD(3) * t117 + t137 * t87;
t133 = t56 * qJD(3);
t69 = pkin(3) * t136 + qJD(2);
t128 = qJ(2) * qJD(3);
t124 = t91 * t94 * t89;
t123 = 0.2e1 * qJD(1);
t122 = pkin(3) * t138;
t119 = t141 * t91;
t116 = t88 * t170;
t114 = -qJD(5) * t101 + qJD(1);
t113 = t89 * t120;
t112 = -t108 * t90 + t6 * t88;
t111 = -t108 * t88 - t6 * t90;
t30 = -pkin(4) * t101 - pkin(7) * t62 + t76;
t33 = -t119 * t87 - t140 * t66;
t9 = t30 * t90 - t33 * t88;
t10 = t30 * t88 + t33 * t90;
t107 = t155 + t160;
t106 = t54 * t56 + t154;
t105 = t47 * t62 + t57 * t59;
t104 = -t41 + (t54 * t88 - t135) * t170;
t103 = t134 * t62 + t59 * t88;
t102 = -t135 * t62 + t59 * t90;
t20 = t140 * t45 - t151;
t16 = -qJD(3) * pkin(4) - t20;
t74 = pkin(3) * t87 + pkin(7);
t99 = t16 * t170 + t46 * t74;
t98 = t137 * t141 - t130;
t97 = -t101 * t12 + t20 * t59 - t21 * t56 - t166;
t96 = -qJD(5) * t112 + t1 * t90 - t2 * t88;
t80 = qJ(2) * t168;
t75 = -pkin(3) * t140 - pkin(4);
t53 = t54 ^ 2;
t49 = t59 * qJD(3);
t48 = -qJD(3) * t119 - t132;
t28 = pkin(4) * t57 - pkin(7) * t54 + t122;
t27 = -pkin(4) * t56 - pkin(7) * t59 + t69;
t26 = t140 * t51 - t151;
t25 = t51 * t87 + t43;
t24 = t140 * t48 + t87 * t98;
t23 = -t140 * t98 + t48 * t87;
t8 = t26 * t90 + t28 * t88;
t7 = -t26 * t88 + t28 * t90;
t4 = -qJD(5) * t10 - t88 * t24 + t90 * t27;
t3 = qJD(5) * t9 + t90 * t24 + t88 * t27;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t80, -0.2e1 * t113, 0.2e1 * t144 * t125, -t147, 0.2e1 * t113, -t146, 0, -t92 * t147 + (qJD(2) * t89 + t128 * t91) * t123, -t92 * t146 + (qJD(2) * t91 - t128 * t89) * t123, 0, t80, t105, t101 * t47 + t46 * t62 + t54 * t59 + t56 * t57, t49, t106, t133, 0, -qJD(3) * t23 - t101 * t64 - t46 * t76 - t54 * t69 - t56 * t65, -qJD(3) * t24 + t47 * t76 + t57 * t69 + t59 * t65 + t62 * t64, t23 * t57 + t24 * t54 + t32 * t47 + t33 * t46 - t97, t12 * t33 - t20 * t23 + t21 * t24 + t64 * t76 + t65 * t69 + t167, t102 * t38 - t152 * t18, -(t156 + t159) * t59 + (t164 - t163 + (-t155 + t160) * qJD(5)) * t62, t101 * t18 + t102 * t170 - t38 * t56 - t41 * t62, t103 * t36 + t14 * t62, t101 * t19 - t103 * t170 + t149 * t62 + t36 * t56, -t170 * t56 + t154, -t101 * t2 + t103 * t16 + t108 * t56 + t165 * t62 + t170 * t4 + t19 * t32 + t23 * t36 - t46 * t9, t1 * t101 + t10 * t46 + t102 * t16 + t11 * t152 - t170 * t3 - t18 * t32 + t23 * t38 + t56 * t6, -t10 * t19 + t18 * t9 - t3 * t36 - t38 * t4 - t112 * t59 + (qJD(5) * t111 - t1 * t88 - t2 * t90) * t62, t1 * t10 - t108 * t4 + t16 * t23 + t2 * t9 + t3 * t6 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t142, 0, 0, 0, 0, 0, 0, t143 * t89, t143 * t91, 0, -t142, 0, 0, 0, 0, 0, 0, qJD(1) * t54 + t49, -qJD(1) * t57 + t133, -t105 - t106, -qJD(1) * t65 + t97, 0, 0, 0, 0, 0, 0, -t101 * t149 - t62 * t19 - t59 * t36 + (-t114 * t90 + t56 * t88) * t170, -t101 * t41 + t62 * t18 - t59 * t38 + (t114 * t88 + t56 * t90) * t170, (-t156 + t159) * t56 + t107 * qJD(1) - (qJD(5) * t107 - t163 - t164) * t101, -qJD(1) * t112 - t101 * t96 + t111 * t56 - t16 * t59 - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t144 * t94, 0, -t124, 0, 0, -t91 * t142, t89 * t142, 0, 0, -t153, -t53 + t169, 0, t153, -t67 + (t57 + t121) * qJD(3), 0, qJD(3) * t25 + t122 * t54 - t57 * t65 - t11, -t118 - t65 * t54 + (t150 * t68 + t26) * qJD(3) + (-t91 * pkin(3) * t57 - t100 * t87) * qJD(1), (t21 - t25) * t57 - (-t20 + t26) * t54 + (-t140 * t47 + t46 * t87) * pkin(3), t20 * t25 - t21 * t26 + (-t11 * t140 + t12 * t87 - t138 * t65) * pkin(3), t115 * t38 - t164, (-t18 + t162) * t90 - t170 * t156 + t145, -t157 + t174, t116 * t36 - t163, t104 + t161, -t170 * t57, -t11 * t90 + t19 * t75 - t25 * t36 + t108 * t57 + (-t134 * t74 - t7) * t170 + t99 * t88, t165 - t18 * t75 - t25 * t38 + t57 * t6 + (t135 * t74 + t8) * t170 + t99 * t90, t36 * t8 + t38 * t7 + (-t19 * t74 - t108 * t54 + t1 + (t38 * t74 + t108) * qJD(5)) * t90 + (-t18 * t74 + t54 * t6 - t2 + (t36 * t74 - t6) * qJD(5)) * t88, t108 * t7 + t11 * t75 - t16 * t25 - t6 * t8 + t74 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 + (t57 - t121) * qJD(3), 0.2e1 * t54 * qJD(3), -t53 - t169, t20 * t57 - t21 * t54 + t64, 0, 0, 0, 0, 0, 0, t104 - t161, -t157 - t174, (t18 + t162) * t90 + t38 * t116 + t145, -t16 * t57 + t172 * t88 + t173 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, -t36 ^ 2 + t38 ^ 2, t170 * t36 - t18, -t158, -t148 + (-qJD(5) + t170) * t38, -t46, -t16 * t38 + t173, t16 * t36 - t172, 0, 0;];
tauc_reg = t5;
