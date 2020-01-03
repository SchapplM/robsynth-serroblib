% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:32
% EndTime: 2019-12-31 22:18:41
% DurationCPUTime: 2.60s
% Computational Cost: add. (2843->285), mult. (7841->532), div. (0->0), fcn. (7107->8), ass. (0->144)
t88 = sin(qJ(3));
t187 = -0.4e1 * t88;
t87 = sin(qJ(4));
t150 = qJD(4) * t87;
t92 = cos(qJ(2));
t156 = qJD(2) * t92;
t85 = sin(pkin(5));
t131 = t85 * t156;
t89 = sin(qJ(2));
t167 = t85 * t89;
t86 = cos(pkin(5));
t91 = cos(qJ(3));
t59 = t91 * t167 + t86 * t88;
t32 = t59 * qJD(3) + t88 * t131;
t58 = t88 * t167 - t86 * t91;
t90 = cos(qJ(4));
t101 = t58 * t150 - t90 * t32;
t186 = t101 * pkin(9);
t50 = pkin(7) * t167 + (-pkin(1) * t92 - pkin(2)) * t86;
t25 = t58 * pkin(3) - t59 * pkin(9) + t50;
t166 = t85 * t92;
t142 = pkin(7) * t166;
t178 = pkin(1) * t89;
t51 = t142 + (pkin(8) + t178) * t86;
t52 = (-pkin(2) * t92 - pkin(8) * t89 - pkin(1)) * t85;
t163 = t91 * t51 + t88 * t52;
t27 = -pkin(9) * t166 + t163;
t185 = t87 * t25 + t90 * t27;
t172 = t88 * pkin(9);
t118 = -t91 * pkin(3) - t172;
t70 = -pkin(2) + t118;
t175 = pkin(8) * t91;
t77 = t90 * t175;
t184 = t87 * t70 + t77;
t33 = -qJD(3) * t58 + t91 * t131;
t183 = -qJD(4) * t166 + t33;
t83 = t90 ^ 2;
t160 = t87 ^ 2 - t83;
t123 = t160 * qJD(4);
t174 = pkin(9) * t91;
t117 = pkin(3) * t88 - t174;
t65 = t117 * qJD(3);
t79 = qJD(4) * t90;
t148 = qJD(4) * t91;
t134 = t87 * t148;
t153 = qJD(3) * t90;
t98 = t88 * t153 + t134;
t30 = t98 * pkin(8) - t87 * t65 - t70 * t79;
t113 = pkin(4) * t87 - qJ(5) * t90;
t107 = pkin(8) + t113;
t53 = t107 * t88;
t57 = t113 * qJD(4) - t87 * qJD(5);
t114 = t90 * pkin(4) + t87 * qJ(5);
t66 = -pkin(3) - t114;
t182 = qJD(3) * (-t66 * t91 + t172) - qJD(4) * t53 - t57 * t88;
t154 = qJD(3) * t88;
t130 = t87 * t154;
t31 = pkin(8) * t130 - qJD(4) * t184 + t90 * t65;
t157 = qJD(2) * t89;
t132 = t85 * t157;
t152 = qJD(3) * t91;
t54 = (pkin(2) * t89 - pkin(8) * t92) * t85 * qJD(2);
t55 = -t86 * pkin(1) * t156 + pkin(7) * t132;
t14 = -t52 * t152 + t51 * t154 - t88 * t54 + t91 * t55;
t12 = pkin(9) * t132 - t14;
t56 = (t86 * t178 + t142) * qJD(2);
t18 = t32 * pkin(3) - t33 * pkin(9) + t56;
t4 = -qJD(4) * t185 - t87 * t12 + t90 * t18;
t181 = t114 * qJD(4) - t90 * qJD(5);
t180 = 0.2e1 * t85;
t179 = 0.2e1 * qJD(5);
t177 = pkin(8) * t85;
t176 = pkin(8) * t87;
t173 = t32 * pkin(4);
t19 = -t90 * t132 + t183 * t87 + t59 * t79;
t171 = t19 * t90;
t20 = t87 * t132 - t59 * t150 + t183 * t90;
t170 = t20 * t87;
t34 = t90 * t166 + t59 * t87;
t169 = t34 * t87;
t35 = -t87 * t166 + t59 * t90;
t168 = t35 * t90;
t165 = t88 * t90;
t164 = t90 * t70;
t82 = t88 ^ 2;
t159 = -t91 ^ 2 + t82;
t158 = t32 * qJ(5);
t155 = qJD(3) * t87;
t151 = qJD(3) * t92;
t149 = qJD(4) * t88;
t147 = t58 * qJD(5);
t145 = t91 * qJD(5);
t144 = qJ(5) * qJD(3);
t143 = t87 * t175;
t141 = -0.2e1 * pkin(2) * qJD(3);
t140 = -0.2e1 * pkin(3) * qJD(4);
t139 = pkin(4) * t154;
t138 = pkin(9) * t150;
t137 = pkin(9) * t79;
t80 = t85 ^ 2;
t136 = t80 * t156;
t133 = t90 * t148;
t129 = t87 * t79;
t128 = t88 * t152;
t127 = t90 * t152;
t126 = t88 * t144;
t124 = -t88 * t51 + t91 * t52;
t122 = t159 * qJD(3);
t121 = 0.2e1 * t128;
t120 = t89 * t136;
t119 = t87 * t127;
t26 = pkin(3) * t166 - t124;
t6 = t58 * qJ(5) + t185;
t8 = t90 * t25 - t87 * t27;
t7 = -t58 * pkin(4) - t8;
t116 = -t6 * t87 + t7 * t90;
t111 = -t34 * t90 - t35 * t87;
t36 = -t91 * qJ(5) + t184;
t37 = -t164 + (pkin(4) + t176) * t91;
t110 = -t36 * t87 + t37 * t90;
t15 = -t51 * t152 - t52 * t154 + t91 * t54 + t88 * t55;
t10 = t34 * pkin(4) - t35 * qJ(5) + t26;
t13 = -pkin(3) * t132 - t15;
t5 = t19 * pkin(4) - t20 * qJ(5) - t35 * qJD(5) + t13;
t106 = -t10 * t79 - t5 * t87;
t105 = t10 * t150 - t5 * t90;
t104 = t13 * t87 + t26 * t79;
t103 = -t13 * t90 + t26 * t150;
t102 = t87 * t32 + t58 * t79;
t3 = -t90 * t12 + t27 * t150 - t87 * t18 - t25 * t79;
t100 = t88 * t151 + t91 * t157;
t99 = -t91 * t151 + t88 * t157;
t60 = -t87 * t149 + t127;
t97 = t102 * pkin(9);
t1 = t147 - t3 + t158;
t2 = -t173 - t4;
t94 = t116 * qJD(4) + t1 * t90 + t2 * t87;
t23 = t126 - t30 - t145;
t28 = -t139 - t31;
t93 = t110 * qJD(4) + t23 * t90 + t28 * t87;
t74 = pkin(9) * t133;
t47 = -t143 + t164;
t29 = t107 * t152 + t181 * t88;
t9 = [0, 0, 0, 0.2e1 * t120, 0.2e1 * (-t89 ^ 2 + t92 ^ 2) * t80 * qJD(2), 0.2e1 * t86 * t131, -0.2e1 * t86 * t132, 0, -0.2e1 * t80 * pkin(1) * t157 - 0.2e1 * t56 * t86, -0.2e1 * pkin(1) * t136 + 0.2e1 * t55 * t86, 0.2e1 * t59 * t33, -0.2e1 * t59 * t32 - 0.2e1 * t33 * t58, (t59 * t157 - t33 * t92) * t180, (-t58 * t157 + t32 * t92) * t180, -0.2e1 * t120, 0.2e1 * t50 * t32 + 0.2e1 * t56 * t58 + 0.2e1 * (t124 * t157 - t15 * t92) * t85, 0.2e1 * t50 * t33 + 0.2e1 * t56 * t59 + 0.2e1 * (-t14 * t92 - t157 * t163) * t85, 0.2e1 * t35 * t20, -0.2e1 * t35 * t19 - 0.2e1 * t20 * t34, 0.2e1 * t20 * t58 + 0.2e1 * t35 * t32, -0.2e1 * t19 * t58 - 0.2e1 * t34 * t32, 0.2e1 * t58 * t32, 0.2e1 * t13 * t34 + 0.2e1 * t26 * t19 + 0.2e1 * t8 * t32 + 0.2e1 * t4 * t58, 0.2e1 * t13 * t35 - 0.2e1 * t185 * t32 + 0.2e1 * t26 * t20 + 0.2e1 * t3 * t58, 0.2e1 * t10 * t19 - 0.2e1 * t2 * t58 - 0.2e1 * t7 * t32 + 0.2e1 * t5 * t34, -0.2e1 * t1 * t34 - 0.2e1 * t6 * t19 + 0.2e1 * t2 * t35 + 0.2e1 * t7 * t20, 0.2e1 * t1 * t58 - 0.2e1 * t10 * t20 + 0.2e1 * t6 * t32 - 0.2e1 * t5 * t35, 0.2e1 * t1 * t6 + 0.2e1 * t10 * t5 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, t131, -t132, 0, -t56, t55, t59 * t152 + t33 * t88, -t88 * t32 + t33 * t91 + (-t58 * t91 - t59 * t88) * qJD(3), t99 * t85, t100 * t85, 0, -pkin(2) * t32 + t154 * t50 - t99 * t177 - t56 * t91, -pkin(2) * t33 - t100 * t177 + t152 * t50 + t56 * t88, t165 * t20 + t35 * t60, t111 * t152 + (-t171 - t170 + (-t168 + t169) * qJD(4)) * t88, (t153 * t58 - t20) * t91 + (qJD(3) * t35 - t101) * t88, (-t155 * t58 + t19) * t91 + (-qJD(3) * t34 - t102) * t88, t154 * t58 - t32 * t91, t31 * t58 + t47 * t32 + (-t4 + (pkin(8) * t34 + t26 * t87) * qJD(3)) * t91 + (pkin(8) * t19 + qJD(3) * t8 + t104) * t88, t30 * t58 - t184 * t32 + (-t3 + (pkin(8) * t35 + t26 * t90) * qJD(3)) * t91 + (pkin(8) * t20 - qJD(3) * t185 - t103) * t88, t53 * t19 - t28 * t58 + t29 * t34 - t37 * t32 + (t10 * t155 + t2) * t91 + (-qJD(3) * t7 - t106) * t88, -t36 * t19 + t37 * t20 - t23 * t34 + t28 * t35 + t116 * t152 + (-t1 * t87 + t2 * t90 + (-t6 * t90 - t7 * t87) * qJD(4)) * t88, -t53 * t20 + t23 * t58 - t29 * t35 + t36 * t32 + (-t10 * t153 - t1) * t91 + (qJD(3) * t6 + t105) * t88, t1 * t36 + t10 * t29 + t2 * t37 + t6 * t23 + t7 * t28 + t5 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -0.2e1 * t122, 0, 0, 0, t88 * t141, t91 * t141, 0.2e1 * t128 * t83 - 0.2e1 * t129 * t82, t119 * t187 + 0.2e1 * t82 * t123, 0.2e1 * t134 * t88 + 0.2e1 * t153 * t159, -0.2e1 * t122 * t87 + 0.2e1 * t133 * t88, -0.2e1 * t128, 0.2e1 * t47 * t154 - 0.2e1 * t31 * t91 + 0.2e1 * (t121 * t87 + t79 * t82) * pkin(8), -0.2e1 * t184 * t154 - 0.2e1 * t30 * t91 + 0.2e1 * (t121 * t90 - t150 * t82) * pkin(8), 0.2e1 * (t155 * t53 + t28) * t91 + 0.2e1 * (-qJD(3) * t37 + t29 * t87 + t53 * t79) * t88, 0.2e1 * t110 * t152 + 0.2e1 * (-t23 * t87 + t28 * t90 + (-t36 * t90 - t37 * t87) * qJD(4)) * t88, 0.2e1 * (-t153 * t53 - t23) * t91 + 0.2e1 * (qJD(3) * t36 + t150 * t53 - t29 * t90) * t88, 0.2e1 * t36 * t23 + 0.2e1 * t37 * t28 + 0.2e1 * t53 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t132, t15, t14, t35 * t79 + t170, qJD(4) * t111 - t87 * t19 + t20 * t90, t102, -t101, 0, -pkin(3) * t19 + t103 - t97, -pkin(3) * t20 + t104 + t186, t66 * t19 + t57 * t34 + t105 - t97, (-t171 + t170 + (t168 + t169) * qJD(4)) * pkin(9) + t94, -t66 * t20 - t57 * t35 + t106 - t186, pkin(9) * t94 + t10 * t57 + t5 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t154, 0, -pkin(8) * t152, pkin(8) * t154, -t123 * t88 + t119, t129 * t187 - t152 * t160, t130 - t133, t98, 0, t74 + (-pkin(3) * t90 + t176) * t149 + (t118 * t87 - t77) * qJD(3), (pkin(8) * t165 + t117 * t87) * qJD(4) + (t118 * t90 + t143) * qJD(3), t74 + (t149 * t66 - t29) * t90 - t182 * t87, t93, (-t29 + (t66 * t88 + t174) * qJD(4)) * t87 + t182 * t90, pkin(9) * t93 + t29 * t66 + t53 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t129, -0.2e1 * t123, 0, 0, 0, t87 * t140, t90 * t140, 0.2e1 * t150 * t66 - 0.2e1 * t57 * t90, 0, -0.2e1 * t57 * t87 - 0.2e1 * t66 * t79, 0.2e1 * t66 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t32, t4, t3, t4 + 0.2e1 * t173, -t20 * pkin(4) - t19 * qJ(5) - t34 * qJD(5), 0.2e1 * t147 - t3 + 0.2e1 * t158, -pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t152 * t87 - t79 * t88, t154, t31, t30, t31 + 0.2e1 * t139, (-pkin(4) * t152 - qJ(5) * t149) * t90 + (-t91 * t144 + (pkin(4) * qJD(4) - qJD(5)) * t88) * t87, 0.2e1 * t126 - t30 - 0.2e1 * t145, -t28 * pkin(4) + t23 * qJ(5) + t36 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t150, 0, -t137, t138, -t137, -t181, -t138, -t181 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, qJ(5) * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t20, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t60, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
