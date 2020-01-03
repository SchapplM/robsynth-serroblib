% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:58
% EndTime: 2019-12-31 22:36:05
% DurationCPUTime: 2.19s
% Computational Cost: add. (3380->269), mult. (9223->510), div. (0->0), fcn. (9049->10), ass. (0->151)
t118 = cos(qJ(5));
t111 = t118 ^ 2;
t114 = sin(qJ(5));
t172 = t114 ^ 2 - t111;
t144 = t172 * qJD(5);
t188 = qJD(3) + qJD(4);
t115 = sin(qJ(4));
t116 = sin(qJ(3));
t120 = cos(qJ(3));
t113 = cos(pkin(5));
t112 = sin(pkin(5));
t121 = cos(qJ(2));
t174 = t112 * t121;
t161 = pkin(7) * t174;
t117 = sin(qJ(2));
t184 = pkin(1) * t117;
t75 = t161 + (pkin(8) + t184) * t113;
t76 = (-pkin(2) * t121 - pkin(8) * t117 - pkin(1)) * t112;
t146 = -t116 * t75 + t120 * t76;
t162 = pkin(3) * t174;
t175 = t112 * t117;
t82 = t113 * t116 + t120 * t175;
t39 = -t82 * pkin(9) + t146 - t162;
t119 = cos(qJ(4));
t180 = t116 * t76 + t120 * t75;
t81 = -t113 * t120 + t116 * t175;
t43 = -t81 * pkin(9) + t180;
t40 = t119 * t43;
t132 = t115 * t39 + t40;
t171 = qJD(2) * t117;
t101 = t112 * t171;
t143 = pkin(3) * t101;
t78 = (pkin(2) * t117 - pkin(8) * t121) * t112 * qJD(2);
t170 = qJD(2) * t121;
t79 = -t113 * pkin(1) * t170 + pkin(7) * t101;
t33 = -t180 * qJD(3) + t116 * t79 + t120 * t78;
t150 = t112 * t170;
t64 = -t81 * qJD(3) + t120 * t150;
t25 = -t64 * pkin(9) + t143 + t33;
t168 = qJD(3) * t120;
t169 = qJD(3) * t116;
t32 = -t116 * t78 + t120 * t79 - t76 * t168 + t75 * t169;
t63 = t82 * qJD(3) + t116 * t150;
t28 = -t63 * pkin(9) - t32;
t147 = t115 * t28 - t119 * t25;
t8 = -t132 * qJD(4) - t147;
t187 = 0.2e1 * t112;
t186 = pkin(8) + pkin(9);
t108 = qJD(5) * t118;
t133 = -t115 * t43 + t119 * t39;
t21 = pkin(4) * t174 - t133;
t20 = t21 * t108;
t6 = -pkin(4) * t101 - t8;
t185 = t6 * t114 + t20;
t183 = pkin(8) * t112;
t89 = t115 * t116 - t119 * t120;
t61 = t188 * t89;
t90 = t115 * t120 + t119 * t116;
t182 = t90 * t61;
t152 = qJD(3) * t186;
t141 = t120 * t152;
t142 = t116 * t152;
t96 = t186 * t116;
t97 = t186 * t120;
t66 = -t115 * t96 + t119 * t97;
t47 = t66 * qJD(4) - t115 * t142 + t119 * t141;
t65 = t115 * t97 + t119 * t96;
t60 = t65 * t108;
t181 = t47 * t114 + t60;
t106 = -t119 * pkin(3) - pkin(4);
t166 = qJD(4) * t115;
t157 = pkin(3) * t166;
t179 = t106 * t108 + t114 * t157;
t178 = pkin(3) * qJD(4);
t53 = t115 * t82 + t119 * t81;
t177 = t119 * t53;
t30 = -qJD(4) * t53 - t115 * t63 + t119 * t64;
t54 = -t115 * t81 + t119 * t82;
t48 = t114 * t54 + t118 * t174;
t15 = -qJD(5) * t48 + t101 * t114 + t118 * t30;
t176 = t15 * t114;
t173 = t114 * t118;
t167 = qJD(3) * t121;
t165 = qJD(4) * t119;
t164 = qJD(5) * t114;
t163 = -0.2e1 * pkin(2) * qJD(3);
t160 = pkin(4) * t164;
t159 = pkin(4) * t108;
t158 = pkin(3) * t169;
t156 = pkin(3) * t165;
t155 = t90 * t164;
t154 = t90 * t108;
t153 = t114 * t174;
t19 = t21 * t164;
t59 = t65 * t164;
t107 = -t120 * pkin(3) - pkin(2);
t109 = t112 ^ 2;
t151 = t109 * t170;
t149 = t114 * t108;
t148 = -t6 * t118 + t19;
t145 = -0.4e1 * t90 * t173;
t140 = t117 * t151;
t139 = -t21 * t61 + t6 * t90;
t31 = qJD(4) * t54 + t115 * t64 + t119 * t63;
t138 = t90 * t31 - t61 * t53;
t137 = t47 * t90 - t61 * t65;
t62 = t188 * t90;
t136 = -t61 * t89 + t90 * t62;
t105 = t115 * pkin(3) + pkin(10);
t135 = t105 * t89 - t106 * t90;
t22 = -pkin(10) * t174 + t132;
t74 = pkin(7) * t175 + (-pkin(1) * t121 - pkin(2)) * t113;
t56 = t81 * pkin(3) + t74;
t29 = t53 * pkin(4) - t54 * pkin(10) + t56;
t10 = t114 * t29 + t118 * t22;
t49 = t118 * t54 - t153;
t134 = -t114 * t49 - t118 * t48;
t57 = t89 * pkin(4) - t90 * pkin(10) + t107;
t42 = t114 * t57 + t118 * t66;
t131 = t106 * t164 - t118 * t157;
t18 = t53 * t108 + t114 * t31;
t130 = -t118 * t31 + t53 * t164;
t129 = t118 * t61 + t155;
t128 = -t118 * t62 + t89 * t164;
t7 = -t115 * t25 - t119 * t28 - t39 * t165 + t43 * t166;
t127 = t116 * t167 + t120 * t171;
t126 = t116 * t171 - t120 * t167;
t125 = t62 * pkin(4) + t61 * pkin(10) + t158;
t80 = (t113 * t184 + t161) * qJD(2);
t124 = pkin(10) * t101 - t7;
t52 = t63 * pkin(3) + t80;
t123 = -t105 * t62 - t106 * t61 + (t115 * t90 - t119 * t89) * t178;
t122 = t31 * pkin(4) - t30 * pkin(10) + t52;
t100 = 0.2e1 * t149;
t91 = -0.2e1 * t140;
t88 = -0.2e1 * t144;
t87 = t90 ^ 2;
t51 = t89 * t108 + t114 * t62;
t46 = t115 * t141 + t119 * t142 + t96 * t165 + t97 * t166;
t41 = -t114 * t66 + t118 * t57;
t38 = -t144 * t90 - t61 * t173;
t34 = qJD(5) * t145 + t172 * t61;
t16 = -qJD(5) * t153 - t101 * t118 + t54 * t108 + t114 * t30;
t14 = -qJD(5) * t42 + t114 * t46 + t118 * t125;
t13 = -t57 * t108 - t114 * t125 + t118 * t46 + t66 * t164;
t12 = t49 * t108 + t176;
t9 = -t114 * t22 + t118 * t29;
t3 = qJD(5) * t134 - t114 * t16 + t15 * t118;
t2 = -qJD(5) * t10 - t114 * t124 + t118 * t122;
t1 = -t29 * t108 - t114 * t122 - t118 * t124 + t22 * t164;
t4 = [0, 0, 0, 0.2e1 * t140, 0.2e1 * (-t117 ^ 2 + t121 ^ 2) * t109 * qJD(2), 0.2e1 * t113 * t150, -0.2e1 * t113 * t101, 0, -0.2e1 * t109 * pkin(1) * t171 - 0.2e1 * t80 * t113, -0.2e1 * pkin(1) * t151 + 0.2e1 * t79 * t113, 0.2e1 * t82 * t64, -0.2e1 * t82 * t63 - 0.2e1 * t64 * t81, (-t121 * t64 + t82 * t171) * t187, (t121 * t63 - t81 * t171) * t187, t91, 0.2e1 * t74 * t63 + 0.2e1 * t80 * t81 + 0.2e1 * (-t33 * t121 + t146 * t171) * t112, 0.2e1 * t74 * t64 + 0.2e1 * t80 * t82 + 0.2e1 * (-t32 * t121 - t180 * t171) * t112, 0.2e1 * t54 * t30, -0.2e1 * t30 * t53 - 0.2e1 * t54 * t31, (-t121 * t30 + t54 * t171) * t187, (t121 * t31 - t53 * t171) * t187, t91, 0.2e1 * t56 * t31 + 0.2e1 * t52 * t53 + 0.2e1 * (-t8 * t121 + t133 * t171) * t112, 0.2e1 * t56 * t30 + 0.2e1 * t52 * t54 + 0.2e1 * (-t7 * t121 - t132 * t171) * t112, 0.2e1 * t49 * t15, -0.2e1 * t15 * t48 - 0.2e1 * t49 * t16, 0.2e1 * t15 * t53 + 0.2e1 * t49 * t31, -0.2e1 * t16 * t53 - 0.2e1 * t48 * t31, 0.2e1 * t53 * t31, 0.2e1 * t21 * t16 + 0.2e1 * t2 * t53 + 0.2e1 * t9 * t31 + 0.2e1 * t6 * t48, 0.2e1 * t1 * t53 - 0.2e1 * t10 * t31 + 0.2e1 * t21 * t15 + 0.2e1 * t6 * t49; 0, 0, 0, 0, 0, t150, -t101, 0, -t80, t79, t64 * t116 + t82 * t168, -t116 * t63 + t64 * t120 + (-t116 * t82 - t120 * t81) * qJD(3), t126 * t112, t127 * t112, 0, -pkin(2) * t63 - t80 * t120 - t126 * t183 + t74 * t169, -pkin(2) * t64 + t80 * t116 - t127 * t183 + t74 * t168, t30 * t90 - t54 * t61, -t30 * t89 - t54 * t62 - t138, (t121 * t61 + t90 * t171) * t112, (t121 * t62 - t89 * t171) * t112, 0, t53 * t158 + t107 * t31 + t52 * t89 + t56 * t62 + (t121 * t47 - t65 * t171) * t112, t54 * t158 + t107 * t30 + t52 * t90 - t56 * t61 + (-t121 * t46 - t66 * t171) * t112, -t49 * t155 + (t15 * t90 - t49 * t61) * t118, -t134 * t61 + (-t176 - t118 * t16 + (t114 * t48 - t118 * t49) * qJD(5)) * t90, t118 * t138 + t15 * t89 - t155 * t53 + t49 * t62, -t114 * t138 - t154 * t53 - t16 * t89 - t48 * t62, t31 * t89 + t53 * t62, t114 * t139 + t14 * t53 + t65 * t16 + t2 * t89 + t20 * t90 + t41 * t31 + t47 * t48 + t9 * t62, t1 * t89 - t10 * t62 + t118 * t139 + t13 * t53 + t65 * t15 - t19 * t90 - t42 * t31 + t47 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t116 * t168, 0.2e1 * (-t116 ^ 2 + t120 ^ 2) * qJD(3), 0, 0, 0, t116 * t163, t120 * t163, -0.2e1 * t182, -0.2e1 * t136, 0, 0, 0, 0.2e1 * t107 * t62 + 0.2e1 * t89 * t158, -0.2e1 * t107 * t61 + 0.2e1 * t158 * t90, -0.2e1 * t111 * t182 - 0.2e1 * t149 * t87, 0.2e1 * t87 * t144 - t61 * t145, 0.2e1 * t118 * t136 - 0.2e1 * t155 * t89, -0.2e1 * t114 * t136 - 0.2e1 * t154 * t89, 0.2e1 * t89 * t62, 0.2e1 * t114 * t137 + 0.2e1 * t14 * t89 + 0.2e1 * t41 * t62 + 0.2e1 * t60 * t90, 0.2e1 * t118 * t137 + 0.2e1 * t13 * t89 - 0.2e1 * t42 * t62 - 0.2e1 * t59 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t63, t101, t33, t32, 0, 0, t30, -t31, t101, t119 * t143 + (-t40 + (-t39 + t162) * t115) * qJD(4) - t147, (-t115 * t171 + t121 * t165) * t112 * pkin(3) + t7, t12, t3, t18, -t130, 0, t106 * t16 - t18 * t105 + (-t114 * t177 + t115 * t48) * t178 + t148, t106 * t15 + t130 * t105 + (t115 * t49 - t118 * t177) * t178 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, -t169, 0, -pkin(8) * t168, pkin(8) * t169, 0, 0, -t61, -t62, 0, -t47, t46, t38, t34, t51, -t128, 0, t59 + (-qJD(5) * t135 - t47) * t118 + t123 * t114, t118 * t123 + t135 * t164 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t157, -0.2e1 * t156, t100, t88, 0, 0, 0, 0.2e1 * t131, 0.2e1 * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, t101, t8, t7, t12, t3, t18, -t130, 0, -pkin(4) * t16 - pkin(10) * t18 + t148, -pkin(4) * t15 + pkin(10) * t130 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t62, 0, -t47, t46, t38, t34, t51, -t128, 0, t59 + (pkin(4) * t61 - pkin(10) * t62) * t114 + (-t47 + (-pkin(4) * t90 - pkin(10) * t89) * qJD(5)) * t118, pkin(4) * t129 + pkin(10) * t128 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, -t156, t100, t88, 0, 0, 0, t131 - t160, -t159 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t88, 0, 0, 0, -0.2e1 * t160, -0.2e1 * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, t31, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t114 * t61 - t154, t62, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t164, 0, -t105 * t108 - t114 * t156, t105 * t164 - t118 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t164, 0, -pkin(10) * t108, pkin(10) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
