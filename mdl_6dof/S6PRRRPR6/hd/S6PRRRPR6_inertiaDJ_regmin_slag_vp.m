% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:46
% EndTime: 2019-03-08 23:36:53
% DurationCPUTime: 2.37s
% Computational Cost: add. (1352->257), mult. (3696->447), div. (0->0), fcn. (3255->10), ass. (0->150)
t100 = sin(qJ(4));
t104 = cos(qJ(4));
t105 = cos(qJ(3));
t165 = qJD(4) * t105;
t147 = t104 * t165;
t101 = sin(qJ(3));
t92 = qJD(3) * t101;
t200 = t100 * t92 - t147;
t169 = qJD(3) * t105;
t151 = t100 * t169;
t91 = qJD(4) * t104;
t199 = t101 * t91 + t151;
t198 = -0.4e1 * t101;
t103 = cos(qJ(6));
t164 = qJD(6) * t103;
t167 = qJD(4) * t100;
t99 = sin(qJ(6));
t177 = qJD(6) * t99;
t197 = -t100 * t164 + t103 * t167 + t104 * t177;
t60 = t100 * t99 + t103 * t104;
t196 = qJD(6) * t60 - t103 * t91;
t141 = t104 * t169;
t166 = qJD(4) * t101;
t150 = t100 * t166;
t49 = t141 - t150;
t96 = t104 ^ 2;
t183 = t100 ^ 2 - t96;
t139 = t183 * qJD(4);
t149 = t100 * t165;
t170 = qJD(3) * t104;
t113 = t101 * t170 + t149;
t187 = pkin(9) * t105;
t131 = pkin(3) * t101 - t187;
t62 = t131 * qJD(3);
t186 = t101 * pkin(9);
t132 = -pkin(3) * t105 - t186;
t66 = -pkin(2) + t132;
t185 = t100 * t62 + t66 * t91;
t23 = pkin(8) * t113 - t185;
t176 = t100 * qJ(5);
t190 = pkin(4) + pkin(5);
t194 = -t104 * t190 - t176;
t178 = qJ(5) * t104;
t123 = pkin(4) * t100 - t178;
t118 = pkin(8) + t123;
t41 = t118 * t101;
t162 = t100 * qJD(5);
t46 = qJD(4) * t123 - t162;
t124 = pkin(4) * t104 + t176;
t63 = -pkin(3) - t124;
t193 = qJD(3) * (-t105 * t63 + t186) - qJD(4) * t41 - t101 * t46;
t161 = t104 * qJD(5);
t192 = qJD(4) * t124 - t161;
t174 = t101 * t104;
t89 = qJ(5) * t92;
t12 = t89 + (-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t174 + (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t100) * t105 + t185;
t188 = pkin(8) * t100;
t86 = t105 * t188;
t93 = t105 * pkin(4);
t26 = t105 * pkin(5) + t86 + t93 + (-pkin(10) * t101 - t66) * t104;
t175 = t100 * t101;
t172 = t104 * t105;
t87 = pkin(8) * t172;
t184 = t100 * t66 + t87;
t36 = -qJ(5) * t105 + t184;
t28 = pkin(10) * t175 + t36;
t128 = t103 * t28 + t26 * t99;
t137 = -pkin(8) * t200 - t104 * t62 + t66 * t167;
t9 = pkin(10) * t150 + (-pkin(10) * t172 - t101 * t190) * qJD(3) + t137;
t2 = -qJD(6) * t128 + t103 * t9 - t99 * t12;
t191 = 0.2e1 * qJD(5);
t189 = pkin(9) - pkin(10);
t95 = t101 ^ 2;
t182 = -t105 ^ 2 + t95;
t102 = sin(qJ(2));
t98 = sin(pkin(6));
t181 = t102 * t98;
t106 = cos(qJ(2));
t180 = t106 * t98;
t179 = cos(pkin(6));
t173 = t103 * t100;
t171 = qJD(2) * t102;
t168 = qJD(3) * t106;
t163 = qJD(6) * t105;
t160 = t105 * qJD(5);
t159 = -0.2e1 * pkin(2) * qJD(3);
t158 = -0.2e1 * pkin(3) * qJD(4);
t157 = pkin(4) * t92;
t156 = pkin(9) * t167;
t155 = pkin(9) * t91;
t154 = t100 * t180;
t75 = t189 * t104;
t153 = t98 * t171;
t152 = qJD(2) * t180;
t144 = t100 * t91;
t143 = t101 * t169;
t140 = t104 * t66 - t86;
t138 = t182 * qJD(3);
t136 = 0.2e1 * t143;
t135 = t189 * t167;
t134 = qJD(4) * t75;
t133 = t100 * t141;
t48 = t101 * t179 + t105 * t181;
t33 = t100 * t48 + t104 * t180;
t34 = t104 * t48 - t154;
t127 = t103 * t34 + t33 * t99;
t126 = t103 * t33 - t34 * t99;
t74 = t189 * t100;
t125 = t103 * t75 + t74 * t99;
t122 = -t100 * t34 + t104 * t33;
t37 = -t140 + t93;
t121 = -t100 * t36 + t104 * t37;
t120 = qJ(5) * t103 - t190 * t99;
t1 = -t103 * t12 - t164 * t26 + t177 * t28 - t9 * t99;
t32 = qJD(3) * t48 + t101 * t152;
t47 = t101 * t181 - t105 * t179;
t20 = -t104 * t32 + t167 * t47;
t117 = -t100 * t190 + t178;
t115 = -pkin(8) + t117;
t31 = -qJD(3) * t47 + t105 * t152;
t10 = -qJD(4) * t154 + t100 * t31 - t104 * t153 + t48 * t91;
t112 = t10 * t105 + t32 * t175 + t199 * t47 - t33 * t92;
t11 = -qJD(4) * t33 + t100 * t153 + t31 * t104;
t109 = qJD(4) * t122 + t10 * t100 + t11 * t104;
t18 = -t23 + t89 - t160;
t21 = t137 - t157;
t108 = qJD(4) * t121 + t21 * t100 + t18 * t104;
t80 = -0.2e1 * t143;
t79 = pkin(9) * t147;
t61 = -t104 * t99 + t173;
t57 = pkin(3) - t194;
t45 = t60 * t101;
t44 = -t101 * t173 + t174 * t99;
t43 = qJD(5) * t99 + qJD(6) * t120;
t42 = qJ(5) * t177 - qJD(5) * t103 + t164 * t190;
t40 = qJD(4) * t117 + t162;
t35 = t115 * t101;
t30 = t91 * t99 - t197;
t29 = -t167 * t99 + t196;
t22 = t101 * t192 + t118 * t169;
t19 = t100 * t32 + t47 * t91;
t17 = qJD(6) * t125 - t103 * t134 - t135 * t99;
t16 = t103 * t135 - t134 * t99 - t164 * t74 + t177 * t75;
t15 = t101 * t196 - t103 * t151 + t49 * t99;
t14 = t101 * t197 - t103 * t141 - t199 * t99;
t13 = (qJD(4) * t194 + t161) * t101 + t115 * t169;
t5 = (t170 * t47 + t11) * t105 + (-qJD(3) * t34 - t20) * t101;
t4 = qJD(6) * t126 + t10 * t99 + t11 * t103;
t3 = qJD(6) * t127 - t10 * t103 + t11 * t99;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t33 + 0.2e1 * t11 * t34 + 0.2e1 * t32 * t47, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t153, -t152, 0, 0, 0, 0, 0 (-t101 * t168 - t105 * t171) * t98 (t101 * t171 - t105 * t168) * t98, 0, 0, 0, 0, 0, t112, t5, t112, t122 * t169 + (t10 * t104 - t100 * t11 + (-t100 * t33 - t104 * t34) * qJD(4)) * t101, -t5, t10 * t37 + t11 * t36 + t18 * t34 + t21 * t33 + t22 * t47 + t32 * t41, 0, 0, 0, 0, 0, -t105 * t3 - t126 * t92 - t15 * t47 - t32 * t44, -t105 * t4 + t127 * t92 + t14 * t47 - t32 * t45; 0, 0, 0, 0, t136, -0.2e1 * t138, 0, 0, 0, t101 * t159, t105 * t159, 0.2e1 * t143 * t96 - 0.2e1 * t144 * t95, t133 * t198 + 0.2e1 * t139 * t95, 0.2e1 * t101 * t149 + 0.2e1 * t170 * t182, -0.2e1 * t100 * t138 + 0.2e1 * t101 * t147, t80, 0.2e1 * t137 * t105 + 0.2e1 * t140 * t92 + 0.2e1 * (t100 * t136 + t91 * t95) * pkin(8), -0.2e1 * t23 * t105 - 0.2e1 * t184 * t92 + 0.2e1 * (t104 * t136 - t167 * t95) * pkin(8), 0.2e1 * (qJD(3) * t100 * t41 + t21) * t105 + 0.2e1 * (-qJD(3) * t37 + t22 * t100 + t41 * t91) * t101, 0.2e1 * t121 * t169 + 0.2e1 * (-t100 * t18 + t104 * t21 + (-t100 * t37 - t104 * t36) * qJD(4)) * t101, 0.2e1 * (-t170 * t41 - t18) * t105 + 0.2e1 * (qJD(3) * t36 - t22 * t104 + t167 * t41) * t101, 0.2e1 * t18 * t36 + 0.2e1 * t21 * t37 + 0.2e1 * t22 * t41, -0.2e1 * t45 * t14, 0.2e1 * t14 * t44 - 0.2e1 * t15 * t45, -0.2e1 * t105 * t14 - 0.2e1 * t45 * t92, -0.2e1 * t105 * t15 + 0.2e1 * t44 * t92, t80, 0.2e1 * t2 * t105 - 0.2e1 * (t103 * t26 - t28 * t99) * t92 + 0.2e1 * t13 * t44 + 0.2e1 * t35 * t15, 0.2e1 * t1 * t105 + 0.2e1 * t128 * t92 + 0.2e1 * t13 * t45 - 0.2e1 * t14 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t31, 0, 0, 0, 0, 0, t20, t19, t20, t109, -t19, pkin(9) * t109 + t32 * t63 + t47 * t46, 0, 0, 0, 0, 0, -t30 * t47 - t32 * t60, t29 * t47 - t32 * t61; 0, 0, 0, 0, 0, 0, t169, -t92, 0, -pkin(8) * t169, pkin(8) * t92, -t101 * t139 + t133, t144 * t198 - t169 * t183, t200, t113, 0, t79 + (-pkin(3) * t104 + t188) * t166 + (t100 * t132 - t87) * qJD(3) (pkin(8) * t174 + t100 * t131) * qJD(4) + (t104 * t132 + t86) * qJD(3), t79 + (t166 * t63 - t22) * t104 - t193 * t100, t108 (-t22 + (t101 * t63 + t187) * qJD(4)) * t100 + t193 * t104, pkin(9) * t108 + t22 * t63 + t41 * t46, -t14 * t61 - t29 * t45, t14 * t60 - t15 * t61 + t29 * t44 - t30 * t45, -t105 * t29 - t61 * t92, -t105 * t30 + t60 * t92, 0, -t17 * t105 - (t103 * t74 - t75 * t99) * t92 + t40 * t44 + t57 * t15 + t13 * t60 + t35 * t30, t105 * t16 + t125 * t92 + t13 * t61 - t14 * t57 - t29 * t35 + t40 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t144, -0.2e1 * t139, 0, 0, 0, t100 * t158, t104 * t158, -0.2e1 * t104 * t46 + 0.2e1 * t167 * t63, 0, -0.2e1 * t100 * t46 - 0.2e1 * t63 * t91, 0.2e1 * t63 * t46, -0.2e1 * t61 * t29, 0.2e1 * t29 * t60 - 0.2e1 * t30 * t61, 0, 0, 0, 0.2e1 * t30 * t57 + 0.2e1 * t40 * t60, -0.2e1 * t29 * t57 + 0.2e1 * t40 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -t10, 0, t11, -pkin(4) * t10 + qJ(5) * t11 + qJD(5) * t34, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t199, t92, -t137, t23, -t137 + 0.2e1 * t157 (-pkin(4) * t169 - qJ(5) * t166) * t104 + (-qJ(5) * t169 + (pkin(4) * qJD(4) - qJD(5)) * t101) * t100, -t23 + 0.2e1 * t89 - 0.2e1 * t160, -pkin(4) * t21 + qJ(5) * t18 + qJD(5) * t36, 0, 0, t14, t15, t92, -t43 * t105 - (-qJ(5) * t99 - t103 * t190) * t92 - t2, t42 * t105 + t120 * t92 - t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t167, 0, -t155, t156, -t155, -t192, -t156, -t192 * pkin(9), 0, 0, t29, t30, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, qJ(5) * t191, 0, 0, 0, 0, 0, 0.2e1 * t43, -0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t49, 0, t21, 0, 0, 0, 0, 0, -t103 * t92 - t163 * t99, -t103 * t163 + t92 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, t155, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, -t92, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, -t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
