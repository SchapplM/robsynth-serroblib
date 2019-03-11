% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaDJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:30
% EndTime: 2019-03-08 19:01:38
% DurationCPUTime: 3.16s
% Computational Cost: add. (4351->278), mult. (12383->502), div. (0->0), fcn. (13581->14), ass. (0->159)
t101 = sin(qJ(6));
t96 = t101 ^ 2;
t105 = cos(qJ(6));
t97 = t105 ^ 2;
t190 = t96 - t97;
t202 = t190 * qJD(6);
t103 = sin(qJ(4));
t106 = cos(qJ(4));
t179 = cos(pkin(7));
t104 = sin(qJ(3));
t99 = sin(pkin(7));
t185 = t104 * t99;
t66 = -t103 * t185 + t106 * t179;
t201 = qJD(4) + qJD(5);
t100 = sin(pkin(6));
t107 = cos(qJ(3));
t178 = cos(pkin(13));
t127 = t179 * t178;
t180 = cos(pkin(6));
t143 = t99 * t180;
t98 = sin(pkin(13));
t53 = -t107 * t143 + (t104 * t98 - t107 * t127) * t100;
t173 = qJD(3) * t107;
t153 = t99 * t173;
t67 = t103 * t179 + t106 * t185;
t110 = t67 * qJD(4) + t103 * t153;
t169 = t106 * qJD(4);
t170 = t103 * qJD(4);
t59 = -t66 * qJD(4) - t106 * t153;
t200 = t110 * t103 - t59 * t106 - t66 * t169 - t67 * t170;
t199 = t103 ^ 2;
t198 = -pkin(10) - pkin(9);
t197 = cos(qJ(5));
t102 = sin(qJ(5));
t54 = t104 * t143 + (t104 * t127 + t107 * t98) * t100;
t65 = -t100 * t178 * t99 + t180 * t179;
t38 = t103 * t65 + t106 * t54;
t50 = t53 * qJD(3);
t112 = -t38 * qJD(4) + t103 * t50;
t37 = -t103 * t54 + t106 * t65;
t20 = t102 * t37 + t197 * t38;
t25 = t37 * qJD(4) - t106 * t50;
t10 = t20 * qJD(5) + t102 * t25 - t197 * t112;
t19 = t102 * t38 - t197 * t37;
t196 = t10 * t19;
t148 = t197 * t106;
t136 = qJD(4) * t148;
t176 = t102 * t103;
t83 = t198 * t106;
t61 = t198 * t176 - t197 * t83;
t152 = qJD(4) * t198;
t76 = t103 * t152;
t36 = t61 * qJD(5) + t102 * t76 - t198 * t136;
t149 = t197 * t103;
t60 = -t102 * t83 - t198 * t149;
t195 = t36 * t60;
t44 = t102 * t66 + t197 * t67;
t23 = t44 * qJD(5) - t102 * t59 + t197 * t110;
t43 = t102 * t67 - t197 * t66;
t194 = t43 * t23;
t51 = t54 * qJD(3);
t29 = t53 * t51;
t144 = qJD(5) * t197;
t57 = -t106 * t144 + t201 * t176 - t136;
t175 = t102 * t106;
t73 = t149 + t175;
t193 = t73 * t57;
t94 = qJD(6) * t105;
t56 = t60 * t94;
t192 = t36 * t101 + t56;
t172 = qJD(5) * t102;
t163 = pkin(4) * t172;
t92 = -t197 * pkin(4) - pkin(5);
t191 = t101 * t163 + t92 * t94;
t189 = pkin(4) * qJD(5);
t188 = t102 * t19;
t187 = t102 * t43;
t186 = t102 * t60;
t184 = t105 * t57;
t183 = t107 * t99;
t182 = t51 * t107;
t174 = qJD(3) * t104;
t171 = qJD(6) * t101;
t58 = t201 * t73;
t72 = -t148 + t176;
t168 = 0.2e1 * t72 * t58;
t167 = -0.2e1 * pkin(3) * qJD(4);
t166 = pkin(4) * t170;
t165 = pkin(5) * t171;
t164 = pkin(5) * t94;
t160 = t101 * t184;
t159 = t53 * t170;
t158 = t73 * t171;
t157 = t73 * t94;
t55 = t60 * t171;
t154 = t99 * t174;
t93 = -pkin(4) * t106 - pkin(3);
t151 = t101 * t197;
t150 = t105 * t197;
t147 = t107 * t170;
t146 = t101 * t94;
t145 = t103 * t169;
t141 = pkin(4) * t144;
t70 = t73 ^ 2;
t140 = t70 * t146;
t139 = t99 ^ 2 * t104 * t173;
t138 = t53 * t154 - t99 * t182;
t135 = t10 * t43 + t19 * t23;
t134 = t10 * t60 + t19 * t36;
t133 = t10 * t73 - t19 * t57;
t132 = t23 * t60 + t43 * t36;
t131 = t23 * t73 - t43 * t57;
t130 = t36 * t73 - t57 * t60;
t129 = t57 * t72 - t58 * t73;
t91 = pkin(4) * t102 + pkin(11);
t128 = t72 * t91 - t73 * t92;
t15 = -t101 * t20 + t105 * t53;
t16 = t101 * t53 + t105 * t20;
t126 = t101 * t16 + t105 * t15;
t116 = -pkin(5) * t72 + pkin(11) * t73 - t93;
t113 = t105 * t116;
t30 = -t101 * t61 - t113;
t31 = -t101 * t116 + t105 * t61;
t125 = t101 * t31 + t105 * t30;
t121 = t101 * t183 - t105 * t44;
t122 = t101 * t44 + t105 * t183;
t124 = -t101 * t121 - t105 * t122;
t123 = -t105 * t163 + t92 * t171;
t120 = t158 + t184;
t119 = -t105 * t58 + t72 * t171;
t117 = (t96 + t97) * t197;
t115 = (t102 * t73 - t197 * t72) * qJD(5);
t114 = pkin(5) * t58 + pkin(11) * t57 + t166;
t9 = t19 * qJD(5) - t102 * t112 - t197 * t25;
t3 = -t51 * t101 + t105 * t9 + t20 * t171 - t53 * t94;
t4 = -t16 * qJD(6) + t101 * t9 + t105 * t51;
t1 = -t126 * qJD(6) - t4 * t101 - t3 * t105;
t35 = t60 * qJD(5) - t152 * t175 - t197 * t76;
t11 = qJD(6) * t113 - t101 * t114 + t105 * t35 + t61 * t171;
t12 = -t31 * qJD(6) + t101 * t35 + t105 * t114;
t2 = -t125 * qJD(6) - t12 * t101 - t11 * t105;
t22 = t102 * t110 - t66 * t144 + t67 * t172 + t197 * t59;
t13 = t122 * qJD(6) - t101 * t154 + t105 * t22;
t14 = t121 * qJD(6) + t101 * t22 + t105 * t154;
t5 = -t124 * qJD(6) - t14 * t101 - t13 * t105;
t111 = pkin(4) * t115 - t57 * t92 - t58 * t91;
t108 = t25 * t106 - t37 * t169 - t50 * t199;
t86 = -0.2e1 * t146;
t85 = 0.2e1 * t146;
t71 = -0.2e1 * t202;
t64 = t117 * t189;
t42 = t101 * t58 + t72 * t94;
t28 = t73 * t202 + t160;
t26 = -0.4e1 * t73 * t146 + t190 * t57;
t18 = -t105 * t23 + t43 * t171;
t17 = t101 * t23 + t43 * t94;
t7 = -t10 * t105 + t19 * t171;
t6 = t10 * t101 + t19 * t94;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t50 * t54 + 0.2e1 * t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t37 * t112 + 0.2e1 * t38 * t25 + 0.2e1 * t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t20 * t9 + 0.2e1 * t196 + 0.2e1 * t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t15 * t4 - 0.2e1 * t16 * t3 + 0.2e1 * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t104 * t50 - t182 + (t104 * t53 + t107 * t54) * qJD(3)) * t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t110 + t112 * t66 + t25 * t67 - t38 * t59 + t138, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t22 - t9 * t44 + t135 + t138, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t3 - t122 * t4 - t13 * t16 + t14 * t15 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66 * t110 - 0.2e1 * t67 * t59 - 0.2e1 * t139, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t44 * t22 - 0.2e1 * t139 + 0.2e1 * t194, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121 * t13 - 0.2e1 * t122 * t14 + 0.2e1 * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t51 + t159, t103 * t51 + t53 * t169, t108, -t51 * pkin(3) + t108 * pkin(9), 0, 0, 0, 0, 0, 0, t51 * t72 + t53 * t58, t51 * t73 - t53 * t57, -t20 * t58 + t72 * t9 + t133, pkin(4) * t159 - t20 * t35 + t51 * t93 - t61 * t9 + t134, 0, 0, 0, 0, 0, 0, t101 * t133 + t15 * t58 + t157 * t19 + t4 * t72, t105 * t133 - t158 * t19 - t16 * t58 + t3 * t72, t126 * t57 + (t101 * t3 - t105 * t4 + (t101 * t15 - t105 * t16) * qJD(6)) * t73, -t11 * t16 + t12 * t15 - t3 * t31 + t30 * t4 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, -t153, 0, 0, 0, 0, 0, 0, 0, 0 (-t106 * t174 - t147) * t99 (t103 * t174 - t107 * t169) * t99, t200, -pkin(3) * t154 + t200 * pkin(9), 0, 0, 0, 0, 0, 0 (-t107 * t58 + t72 * t174) * t99 (t107 * t57 + t73 * t174) * t99, t22 * t72 - t44 * t58 + t131, -t22 * t61 - t44 * t35 + (-pkin(4) * t147 + t93 * t174) * t99 + t132, 0, 0, 0, 0, 0, 0, t101 * t131 - t122 * t58 + t14 * t72 + t157 * t43, t105 * t131 + t121 * t58 + t13 * t72 - t158 * t43, t124 * t57 + (t101 * t13 - t105 * t14 + (-t101 * t122 + t105 * t121) * qJD(6)) * t73, t11 * t121 - t12 * t122 - t13 * t31 + t14 * t30 + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t145, 0.2e1 * (t106 ^ 2 - t199) * qJD(4), 0, -0.2e1 * t145, 0, 0, t103 * t167, t106 * t167, 0, 0, -0.2e1 * t193, 0.2e1 * t129, 0, t168, 0, 0, 0.2e1 * t72 * t166 + 0.2e1 * t58 * t93, 0.2e1 * t73 * t166 - 0.2e1 * t57 * t93, 0.2e1 * t35 * t72 - 0.2e1 * t58 * t61 + 0.2e1 * t130, 0.2e1 * t166 * t93 - 0.2e1 * t35 * t61 + 0.2e1 * t195, -0.2e1 * t97 * t193 - 0.2e1 * t140, 0.4e1 * t73 * t160 + 0.2e1 * t70 * t202, -0.2e1 * t105 * t129 - 0.2e1 * t158 * t72, -0.2e1 * t96 * t193 + 0.2e1 * t140, 0.2e1 * t101 * t129 - 0.2e1 * t157 * t72, t168, 0.2e1 * t101 * t130 + 0.2e1 * t12 * t72 + 0.2e1 * t30 * t58 + 0.2e1 * t56 * t73, 0.2e1 * t105 * t130 + 0.2e1 * t11 * t72 - 0.2e1 * t31 * t58 - 0.2e1 * t55 * t73, 0.2e1 * t125 * t57 + 0.2e1 * (t101 * t11 - t105 * t12 + (t101 * t30 - t105 * t31) * qJD(6)) * t73, -0.2e1 * t11 * t31 + 0.2e1 * t12 * t30 + 0.2e1 * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0 (-t197 * t10 - t102 * t9 + (t197 * t20 + t188) * qJD(5)) * pkin(4), 0, 0, 0, 0, 0, 0, t7, t6, t1, t10 * t92 + (-t15 * t151 + t150 * t16 + t188) * t189 + t1 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0 (-t197 * t23 - t102 * t22 + (t197 * t44 + t187) * qJD(5)) * pkin(4), 0, 0, 0, 0, 0, 0, t18, t17, t5, t23 * t92 + (-t121 * t150 + t122 * t151 + t187) * t189 + t5 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, 0, -t170, 0, -pkin(9) * t169, pkin(9) * t170, 0, 0, 0, 0, -t57, 0, -t58, 0, -t36, t35 (-t102 * t58 + t197 * t57 + t115) * pkin(4) (-t197 * t36 - t102 * t35 + (t197 * t61 + t186) * qJD(5)) * pkin(4), -t28, t26, t42, t28, -t119, 0, t55 + (-qJD(6) * t128 - t36) * t105 + t111 * t101, t105 * t111 + t128 * t171 + t192, t2, t36 * t92 + (t150 * t31 - t151 * t30 + t186) * t189 + t2 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t163, -0.2e1 * t141, 0, 0, t85, t71, 0, t86, 0, 0, 0.2e1 * t123, 0.2e1 * t191, 0.2e1 * t64, 0.2e1 * (t102 * t92 + t117 * t91) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, t6, t1, -pkin(5) * t10 + pkin(11) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, t5, -pkin(5) * t23 + pkin(11) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, -t58, 0, -t36, t35, 0, 0, -t28, t26, t42, t28, -t119, 0, t55 + (pkin(5) * t57 - pkin(11) * t58) * t101 + (-t36 + (-pkin(5) * t73 - pkin(11) * t72) * qJD(6)) * t105, pkin(5) * t120 + pkin(11) * t119 + t192, t2, -pkin(5) * t36 + pkin(11) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t141, 0, 0, t85, t71, 0, t86, 0, 0, t123 - t165, -t164 + t191, t64 (-pkin(5) * t102 + pkin(11) * t117) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t71, 0, t86, 0, 0, -0.2e1 * t165, -0.2e1 * t164, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, 0, t101 * t57 - t157, t58, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, -t171, 0, -t101 * t141 - t91 * t94, -t105 * t141 + t91 * t171, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, -t171, 0, -pkin(11) * t94, pkin(11) * t171, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
