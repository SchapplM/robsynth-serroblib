% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:32
% EndTime: 2019-03-09 03:42:40
% DurationCPUTime: 3.76s
% Computational Cost: add. (5093->286), mult. (11507->510), div. (0->0), fcn. (10852->10), ass. (0->148)
t105 = sin(pkin(11));
t106 = cos(pkin(11));
t202 = sin(qJ(5));
t204 = cos(qJ(5));
t75 = t204 * t105 + t202 * t106;
t162 = qJD(5) * t204;
t161 = qJD(5) * t202;
t92 = t105 * t161;
t136 = -t106 * t162 + t92;
t109 = cos(qJ(3));
t108 = sin(qJ(3));
t200 = t109 * pkin(3);
t98 = -cos(pkin(10)) * pkin(1) - pkin(2);
t133 = t98 - t200;
t196 = pkin(8) + qJ(4);
t128 = -t196 * t108 + t133;
t97 = sin(pkin(10)) * pkin(1) + pkin(7);
t163 = t105 * t97 + pkin(4);
t119 = t128 * t106 - t163 * t109;
t190 = t109 * t97;
t175 = t106 * t190;
t123 = t105 * t128 + t175;
t23 = t204 * t119 - t202 * t123;
t211 = t75 * t196;
t100 = qJD(3) * t108;
t144 = t202 * t105 - t204 * t106;
t191 = t105 * t162 + t106 * t161;
t159 = t109 * t191;
t210 = t144 * t100 - t159;
t209 = -0.2e1 * t136;
t208 = 0.2e1 * t191;
t24 = t202 * t119 + t204 * t123;
t184 = qJD(3) * t109;
t151 = t202 * t184;
t152 = t204 * t184;
t131 = t105 * t152 + t106 * t151 - t136 * t108;
t192 = t75 * t108;
t206 = t75 * t131 - t192 * t136;
t205 = 0.2e1 * qJD(3);
t203 = cos(qJ(6));
t201 = pkin(3) * t108;
t107 = sin(qJ(6));
t68 = t144 * t108;
t45 = -t107 * t192 - t203 * t68;
t48 = t105 * t151 - t106 * t152 + t191 * t108;
t16 = t45 * qJD(6) - t107 * t48 + t203 * t131;
t145 = t203 * t192;
t44 = -t107 * t68 + t145;
t199 = t44 * t16;
t182 = qJD(6) * t107;
t15 = qJD(6) * t145 + t107 * t131 - t68 * t182 + t203 * t48;
t198 = t45 * t15;
t197 = t68 * t48;
t55 = -t107 * t144 + t203 * t75;
t26 = t55 * qJD(6) - t107 * t136 + t203 * t191;
t132 = t203 * t144;
t54 = t107 * t75 + t132;
t194 = t15 * t54 - t45 * t26;
t193 = t48 * t144 + t68 * t191;
t146 = t204 * t196;
t147 = t196 * t202;
t59 = -t105 * t147 + t106 * t146;
t169 = t105 * t184;
t76 = t97 * t184;
t64 = pkin(4) * t169 + t76;
t188 = t105 * t108;
t83 = t108 * t97;
t70 = pkin(4) * t188 + t83;
t189 = qJ(4) * t109;
t187 = t108 * qJ(4);
t101 = t105 ^ 2;
t102 = t106 ^ 2;
t186 = t101 + t102;
t185 = t108 ^ 2 - t109 ^ 2;
t183 = qJD(4) * t109;
t181 = t108 * qJD(4);
t180 = t98 * t205;
t179 = pkin(5) * t100;
t178 = pkin(5) * t182;
t177 = t105 * t190;
t176 = t106 * t83;
t174 = t191 * pkin(5);
t99 = -t106 * pkin(4) - pkin(3);
t170 = t101 * t184;
t168 = t105 * t181;
t167 = t106 * t100;
t166 = t106 * t184;
t165 = t106 * t181;
t164 = t108 * t184;
t160 = qJD(6) * t203;
t158 = t204 * qJD(4);
t157 = t202 * qJD(4);
t156 = t186 * qJD(4);
t155 = 0.2e1 * t164;
t154 = pkin(5) * t160;
t153 = t105 * t166;
t149 = 0.2e1 * t156;
t148 = -0.2e1 * t185 * qJD(3);
t25 = qJD(6) * t132 + t107 * t191 + t203 * t136 + t75 * t182;
t143 = -t16 * t55 + t25 * t44;
t142 = -t187 - t200;
t141 = -t189 + t201;
t52 = -t165 + (t106 * t141 + t97 * t188) * qJD(3);
t53 = -t168 + (t105 * t141 - t176) * qJD(3);
t140 = -t105 * t52 + t106 * t53;
t130 = t133 - t187;
t56 = t106 * t130 - t177;
t57 = t105 * t130 + t175;
t139 = -t105 * t56 + t106 * t57;
t134 = -t196 * t109 + t201;
t117 = -t165 + (t134 * t106 + t163 * t108) * qJD(3);
t120 = -t168 + (t134 * t105 - t176) * qJD(3);
t11 = -t24 * qJD(5) + t204 * t117 - t202 * t120;
t110 = t48 * pkin(9) + t11 + t179;
t113 = -t109 * pkin(5) + t68 * pkin(9) + t23;
t111 = t203 * t113;
t10 = -t23 * qJD(5) - t202 * t117 - t204 * t120;
t114 = -pkin(9) * t131 - t10;
t20 = -t192 * pkin(9) + t24;
t1 = -qJD(6) * t111 - t107 * t110 - t203 * t114 + t20 * t182;
t135 = t54 * t100 - t109 * t26;
t127 = t192 * t131;
t126 = (qJD(5) * t147 - t158) * t105;
t39 = qJD(5) * t211 + t105 * t157 - t106 * t158;
t125 = -t75 * pkin(9) - t211;
t124 = t107 * t125;
t122 = t203 * t125;
t121 = -t191 * pkin(9) - t39;
t115 = t92 * pkin(9) + (-t157 + (-t204 * pkin(9) - t146) * qJD(5)) * t106 + t126;
t112 = t107 * t113;
t2 = -qJD(6) * t112 - t107 * t114 + t203 * t110 - t20 * t160;
t90 = t102 * t184;
t89 = t105 * t100;
t88 = -0.2e1 * t164;
t63 = t144 * pkin(5) + t99;
t50 = t192 * pkin(5) + t70;
t49 = t75 * t100 + t136 * t109;
t46 = -t144 * pkin(9) + t59;
t40 = (-qJD(5) * t146 - t157) * t106 + t126;
t29 = pkin(5) * t131 + t64;
t22 = t203 * t46 + t124;
t21 = -t107 * t46 + t122;
t18 = t55 * t100 + t109 * t25;
t9 = -qJD(6) * t124 - t107 * t121 + t203 * t115 - t46 * t160;
t8 = -qJD(6) * t122 - t107 * t115 - t203 * t121 + t46 * t182;
t7 = t203 * t20 + t112;
t6 = -t107 * t20 + t111;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t148, 0, t88, 0, 0, t108 * t180, t109 * t180, 0, 0, t102 * t155, -0.4e1 * t108 * t153, t185 * t106 * t205, t101 * t155, t105 * t148, t88, -0.2e1 * t109 * t52 + 0.2e1 * (t56 + 0.2e1 * t177) * t100, 0.2e1 * t109 * t53 + 0.2e1 * (-t57 + 0.2e1 * t175) * t100, 0.2e1 * (-t105 * t53 - t106 * t52) * t108 + 0.2e1 * (-t105 * t57 - t106 * t56) * t184, 0.2e1 * t97 ^ 2 * t164 + 0.2e1 * t56 * t52 + 0.2e1 * t57 * t53, 0.2e1 * t197, 0.2e1 * t68 * t131 + 0.2e1 * t192 * t48, -0.2e1 * t68 * t100 + 0.2e1 * t109 * t48, 0.2e1 * t127, -0.2e1 * t192 * t100 + 0.2e1 * t109 * t131, t88, 0.2e1 * t23 * t100 - 0.2e1 * t11 * t109 + 0.2e1 * t70 * t131 + 0.2e1 * t64 * t192, -0.2e1 * t10 * t109 - 0.2e1 * t24 * t100 - 0.2e1 * t70 * t48 - 0.2e1 * t64 * t68, 0.2e1 * t10 * t192 + 0.2e1 * t11 * t68 - 0.2e1 * t24 * t131 + 0.2e1 * t23 * t48, -0.2e1 * t10 * t24 + 0.2e1 * t11 * t23 + 0.2e1 * t64 * t70, -0.2e1 * t198, 0.2e1 * t44 * t15 - 0.2e1 * t45 * t16, 0.2e1 * t45 * t100 + 0.2e1 * t109 * t15, 0.2e1 * t199, -0.2e1 * t44 * t100 + 0.2e1 * t109 * t16, t88, 0.2e1 * t6 * t100 - 0.2e1 * t109 * t2 + 0.2e1 * t50 * t16 + 0.2e1 * t29 * t44, -0.2e1 * t1 * t109 - 0.2e1 * t7 * t100 - 0.2e1 * t50 * t15 + 0.2e1 * t29 * t45, 0.2e1 * t1 * t44 + 0.2e1 * t15 * t6 - 0.2e1 * t16 * t7 - 0.2e1 * t2 * t45, -0.2e1 * t1 * t7 + 0.2e1 * t2 * t6 + 0.2e1 * t29 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t108 + (t139 * t109 + t185 * t97) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t68 + t70 * t100 - t64 * t109 - t11 * t192 - t23 * t131 - t24 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t45 + t50 * t100 - t109 * t29 - t7 * t15 - t6 * t16 - t2 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t186) * t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t127 - 0.2e1 * t164 + 0.2e1 * t197, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t164 - 0.2e1 * t198 + 0.2e1 * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, 0, -t100, 0, -t76, t97 * t100, 0, 0, t153, t90 - t170, t89, -t153, t167, 0, t105 * t183 + (t142 * t105 - t175) * qJD(3), t106 * t183 + (t142 * t106 + t177) * qJD(3), t140, -pkin(3) * t76 + t140 * qJ(4) + t139 * qJD(4), t68 * t136 - t48 * t75, t193 - t206, t49, t131 * t144 + t192 * t191, -t210, 0, -t100 * t211 - t40 * t109 + t99 * t131 + t64 * t144 + t70 * t191, -t59 * t100 - t39 * t109 - t136 * t70 - t99 * t48 + t64 * t75, t10 * t144 - t11 * t75 - t59 * t131 + t23 * t136 - t24 * t191 + t39 * t192 - t211 * t48 + t40 * t68, -t10 * t59 - t11 * t211 + t23 * t40 - t24 * t39 + t64 * t99, -t15 * t55 - t25 * t45, t143 + t194, t18, t16 * t54 + t26 * t44, -t135, 0, t21 * t100 - t9 * t109 + t63 * t16 + t44 * t174 + t50 * t26 + t29 * t54, -t22 * t100 - t8 * t109 - t63 * t15 + t45 * t174 - t50 * t25 + t29 * t55, t1 * t54 + t15 * t21 - t16 * t22 - t2 * t55 + t25 * t6 - t26 * t7 + t44 * t8 - t45 * t9, -t1 * t22 + t50 * t174 + t2 * t21 + t29 * t63 + t6 * t9 - t7 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t184, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t89, t90 + t170, t108 * t156 + (t186 * t189 - t201) * qJD(3), 0, 0, 0, 0, 0, 0, t210, t49, t193 + t206, t99 * t100 + t131 * t211 - t192 * t40 + t68 * t39 - t48 * t59, 0, 0, 0, 0, 0, 0, t135, t18, -t143 + t194, -pkin(5) * t159 + t63 * t100 - t15 * t22 - t16 * t21 - t44 * t9 - t45 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, qJ(4) * t149, t75 * t209, 0.2e1 * t136 * t144 - 0.2e1 * t75 * t191, 0, t144 * t208, 0, 0, t99 * t208, t99 * t209, -0.2e1 * t136 * t211 + 0.2e1 * t39 * t144 - 0.2e1 * t59 * t191 - 0.2e1 * t40 * t75, -0.2e1 * t211 * t40 - 0.2e1 * t39 * t59, -0.2e1 * t55 * t25, 0.2e1 * t25 * t54 - 0.2e1 * t26 * t55, 0, 0.2e1 * t54 * t26, 0, 0, 0.2e1 * t54 * t174 + 0.2e1 * t63 * t26, 0.2e1 * t55 * t174 - 0.2e1 * t63 * t25, 0.2e1 * t21 * t25 - 0.2e1 * t22 * t26 + 0.2e1 * t54 * t8 - 0.2e1 * t55 * t9, 0.2e1 * t63 * t174 + 0.2e1 * t21 * t9 - 0.2e1 * t22 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t166, 0, t76, 0, 0, 0, 0, 0, 0, t131, -t48, 0, t64, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, -t136, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, -t131, t100, t11, t10, 0, 0, 0, 0, -t15, 0, -t16, t100, t109 * t178 + t203 * t179 + t2 (-t107 * t100 + t109 * t160) * pkin(5) + t1 (t203 * t15 - t107 * t16 + (t107 * t45 - t203 * t44) * qJD(6)) * pkin(5) (t203 * t2 - t1 * t107 + (-t107 * t6 + t203 * t7) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t48, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0 (-t203 * t16 - t107 * t15 + (t107 * t44 + t203 * t45) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, 0, -t191, 0, t40, t39, 0, 0, 0, 0, -t25, 0, -t26, 0, t9, t8 (t203 * t25 - t107 * t26 + (t107 * t55 - t203 * t54) * qJD(6)) * pkin(5) (t203 * t9 - t107 * t8 + (-t107 * t21 + t203 * t22) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t178, -0.2e1 * t154, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t16, t100, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, -t26, 0, t9, t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, -t154, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
