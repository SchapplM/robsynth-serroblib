% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:15
% EndTime: 2019-03-09 05:07:19
% DurationCPUTime: 1.62s
% Computational Cost: add. (1227->221), mult. (2872->364), div. (0->0), fcn. (2303->8), ass. (0->137)
t92 = sin(qJ(3));
t181 = -0.4e1 * t92;
t93 = cos(qJ(6));
t153 = qJD(6) * t93;
t90 = sin(qJ(6));
t154 = qJD(6) * t90;
t91 = sin(qJ(4));
t157 = qJD(4) * t91;
t94 = cos(qJ(4));
t83 = qJD(4) * t94;
t26 = t91 * t153 - t94 * t154 - t93 * t157 + t90 * t83;
t51 = t90 * t91 + t93 * t94;
t25 = t51 * qJD(6) - t90 * t157 - t93 * t83;
t156 = qJD(4) * t92;
t124 = qJ(5) * t156;
t147 = qJ(5) * qJD(3);
t95 = cos(qJ(3));
t125 = t95 * t147;
t150 = t94 * qJD(5);
t165 = t94 * t125 + t92 * t150;
t149 = t95 * qJD(3);
t133 = t91 * t149;
t39 = t92 * t83 + t133;
t180 = -t39 * pkin(4) - t91 * t124 + t165;
t172 = t92 * pkin(8);
t116 = -t95 * pkin(3) - t172;
t79 = -cos(pkin(10)) * pkin(1) - pkin(2);
t46 = t116 + t79;
t174 = pkin(8) * t95;
t115 = pkin(3) * t92 - t174;
t55 = t115 * qJD(3);
t166 = t46 * t83 + t91 * t55;
t155 = qJD(4) * t95;
t138 = t91 * t155;
t81 = t92 * qJD(3);
t38 = t94 * t81 + t138;
t78 = sin(pkin(10)) * pkin(1) + pkin(7);
t11 = t38 * t78 - t166;
t161 = t91 * qJ(5);
t176 = pkin(4) + pkin(5);
t179 = -t176 * t94 - t161;
t162 = qJ(5) * t94;
t111 = pkin(4) * t91 - t162;
t27 = (t111 + t78) * t92;
t151 = t91 * qJD(5);
t36 = t111 * qJD(4) - t151;
t112 = -t94 * pkin(4) - t161;
t56 = -pkin(3) + t112;
t178 = (-t56 * t95 + t172) * qJD(3) - qJD(4) * t27 - t36 * t92;
t173 = pkin(9) * t92;
t170 = t78 * t91;
t53 = t95 * t170;
t84 = t95 * pkin(4);
t20 = t95 * pkin(5) + t53 + t84 + (-t46 - t173) * t94;
t35 = t91 * t46;
t167 = t94 * t95;
t54 = t78 * t167;
t22 = -t95 * qJ(5) + t35 + t54;
t21 = t91 * t173 + t22;
t110 = t90 * t20 + t93 * t21;
t130 = t78 * t81;
t137 = t94 * t155;
t120 = -t91 * t130 + t78 * t137 + t46 * t157 - t94 * t55;
t139 = t91 * t156;
t5 = pkin(9) * t139 + (-pkin(9) * t167 - t176 * t92) * qJD(3) + t120;
t169 = t92 * t94;
t80 = t92 * t147;
t6 = t80 + (pkin(9) * qJD(4) - qJD(3) * t78) * t169 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t78) * t91) * t95 + t166;
t2 = -t110 * qJD(6) + t93 * t5 - t90 * t6;
t177 = 0.2e1 * qJD(5);
t175 = pkin(8) - pkin(9);
t171 = t56 * t92;
t168 = t93 * t91;
t85 = t91 ^ 2;
t87 = t94 ^ 2;
t164 = t85 - t87;
t86 = t92 ^ 2;
t163 = -t95 ^ 2 + t86;
t160 = qJD(3) * t27;
t159 = qJD(3) * t94;
t158 = qJD(4) * t86;
t152 = qJD(6) * t95;
t148 = t95 * qJD(5);
t146 = -0.2e1 * pkin(3) * qJD(4);
t145 = 0.2e1 * qJD(3) * t79;
t144 = pkin(4) * t81;
t143 = pkin(8) * t157;
t142 = pkin(8) * t83;
t141 = t176 * t91;
t64 = t175 * t94;
t140 = t78 * t158;
t135 = t85 * t149;
t132 = t91 * t83;
t131 = t92 * t149;
t127 = t94 * t149;
t126 = t78 * t149;
t123 = -t94 * t46 + t53;
t122 = t164 * qJD(4);
t121 = t163 * qJD(3);
t119 = t91 * t127;
t118 = t175 * t157;
t117 = qJD(4) * t64;
t23 = t123 + t84;
t109 = t22 * t94 + t23 * t91;
t108 = -t22 * t91 + t23 * t94;
t63 = t175 * t91;
t107 = t90 * t63 + t93 * t64;
t106 = t93 * qJ(5) - t176 * t90;
t104 = -t95 * t26 + t51 * t81;
t1 = -t20 * t153 + t21 * t154 - t90 * t5 - t93 * t6;
t103 = -t141 + t162;
t99 = t112 * qJD(4) + t150;
t7 = t80 - t11 - t148;
t9 = t120 - t144;
t97 = t108 * qJD(4) + t7 * t94 + t9 * t91;
t73 = t87 * t149;
t68 = -0.2e1 * t131;
t67 = pkin(8) * t137;
t61 = t87 * t131;
t52 = -t90 * t94 + t168;
t47 = pkin(3) - t179;
t40 = -t91 * t81 + t137;
t37 = -t127 + t139;
t34 = t51 * t92;
t33 = -t92 * t168 + t90 * t169;
t32 = t90 * qJD(5) + t106 * qJD(6);
t31 = qJ(5) * t154 - t93 * qJD(5) + t153 * t176;
t28 = t103 * qJD(4) + t151;
t24 = (-t78 + t103) * t92;
t18 = -t25 * t95 - t52 * t81;
t17 = t126 - t180;
t16 = t107 * qJD(6) - t93 * t117 - t90 * t118;
t15 = -t90 * t117 + t93 * t118 - t63 * t153 + t64 * t154;
t14 = t51 * t149 + t26 * t92;
t13 = t90 * t127 - t93 * t133 + t25 * t92;
t8 = t179 * t156 + (-t78 - t141) * t149 + t165;
t3 = [0, 0, 0, 0, 0.2e1 * t131, -0.2e1 * t121, 0, 0, 0, t92 * t145, t95 * t145, -0.2e1 * t86 * t132 + 0.2e1 * t61, t119 * t181 + 0.2e1 * t164 * t158, 0.2e1 * t92 * t138 + 0.2e1 * t163 * t159, -0.2e1 * t91 * t121 + 0.2e1 * t92 * t137, t68, 0.2e1 * t94 * t140 + 0.2e1 * t120 * t95 + 0.2e1 * (-t123 + 0.2e1 * t53) * t81, -0.2e1 * t91 * t140 - 0.2e1 * t11 * t95 + 0.2e1 * (-t35 + t54) * t81, 0.2e1 * (t91 * t160 + t9) * t95 + 0.2e1 * (-qJD(3) * t23 + t17 * t91 + t27 * t83) * t92, 0.2e1 * t108 * t149 + 0.2e1 * (-t109 * qJD(4) - t7 * t91 + t9 * t94) * t92, 0.2e1 * (-t27 * t159 - t7) * t95 + 0.2e1 * (qJD(3) * t22 + t27 * t157 - t17 * t94) * t92, 0.2e1 * t27 * t17 + 0.2e1 * t22 * t7 + 0.2e1 * t23 * t9, 0.2e1 * t34 * t14, -0.2e1 * t34 * t13 - 0.2e1 * t14 * t33, 0.2e1 * t14 * t95 - 0.2e1 * t34 * t81, -0.2e1 * t95 * t13 + 0.2e1 * t33 * t81, t68, 0.2e1 * t2 * t95 - 0.2e1 * (t93 * t20 - t90 * t21) * t81 + 0.2e1 * t8 * t33 + 0.2e1 * t24 * t13, 0.2e1 * t1 * t95 + 0.2e1 * t110 * t81 + 0.2e1 * t24 * t14 + 0.2e1 * t8 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t109 * qJD(3) - t17) * t95 + (t97 + t160) * t92, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t61 + 0.2e1 * (t85 - 0.1e1) * t131, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t149, -t81, 0, -t126, t130, -t92 * t122 + t119, t132 * t181 - t135 + t73, -t40, t38, 0, t67 + (-pkin(3) * t94 + t170) * t156 + (t116 * t91 - t54) * qJD(3) (t115 * t91 + t78 * t169) * qJD(4) + (t116 * t94 + t53) * qJD(3), t67 + (t56 * t156 - t17) * t94 - t178 * t91, t97 (-t17 + (t171 + t174) * qJD(4)) * t91 + t178 * t94, t97 * pkin(8) + t17 * t56 + t27 * t36, t14 * t52 - t34 * t25, -t52 * t13 - t14 * t51 + t25 * t33 - t34 * t26, t18, t104, 0, -t16 * t95 - (t93 * t63 - t90 * t64) * t81 + t28 * t33 + t47 * t13 + t8 * t51 + t24 * t26, t107 * t81 + t47 * t14 + t15 * t95 - t24 * t25 + t28 * t34 + t8 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t149, 0, 0, 0, 0, 0, -t38, -t40, -t38, t73 + t135, t40, -t95 * t36 + (t171 + (t85 + t87) * t174) * qJD(3), 0, 0, 0, 0, 0, -t104, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132, -0.2e1 * t122, 0, 0, 0, t91 * t146, t94 * t146, 0.2e1 * t56 * t157 - 0.2e1 * t36 * t94, 0, -0.2e1 * t36 * t91 - 0.2e1 * t56 * t83, 0.2e1 * t56 * t36, -0.2e1 * t52 * t25, 0.2e1 * t25 * t51 - 0.2e1 * t52 * t26, 0, 0, 0, 0.2e1 * t47 * t26 + 0.2e1 * t28 * t51, -0.2e1 * t47 * t25 + 0.2e1 * t28 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39, t81, -t120, t11, -t120 + 0.2e1 * t144 (-pkin(4) * t149 - t124) * t94 + (-t125 + (pkin(4) * qJD(4) - qJD(5)) * t92) * t91, 0.2e1 * t80 - t11 - 0.2e1 * t148, -t9 * pkin(4) + t7 * qJ(5) + t22 * qJD(5), 0, 0, -t14, t13, t81, -t32 * t95 - (-t90 * qJ(5) - t176 * t93) * t81 - t2, t106 * t81 + t31 * t95 - t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t37, -t39, 0, -t37, t180, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t157, 0, -t142, t143, -t142, t99, -t143, t99 * pkin(8), 0, 0, t25, t26, 0, t16, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, qJ(5) * t177, 0, 0, 0, 0, 0, 0.2e1 * t32, -0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t37, 0, t9, 0, 0, 0, 0, 0, -t90 * t152 - t93 * t81, -t93 * t152 + t90 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, t142, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t81, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, -t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
