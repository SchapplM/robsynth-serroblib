% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:38
% EndTime: 2019-03-08 18:43:44
% DurationCPUTime: 2.20s
% Computational Cost: add. (4714->293), mult. (13022->448), div. (0->0), fcn. (11767->14), ass. (0->169)
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t104 = sin(qJ(5));
t107 = cos(qJ(5));
t105 = sin(qJ(3));
t100 = cos(pkin(12));
t101 = cos(pkin(7));
t166 = t100 * t101;
t98 = sin(pkin(6));
t170 = qJD(1) * t98;
t108 = cos(qJ(3));
t96 = sin(pkin(12));
t173 = t108 * t96;
t102 = cos(pkin(6));
t86 = qJD(1) * t102 + qJD(2);
t97 = sin(pkin(7));
t185 = t86 * t97;
t56 = (t105 * t166 + t173) * t170 + t105 * t185;
t99 = cos(pkin(13));
t51 = t99 * t56;
t145 = t108 * t166;
t177 = t105 * t96;
t149 = t98 * t177;
t79 = t108 * t185;
t55 = -qJD(1) * t149 + t145 * t170 + t79;
t53 = qJD(3) * pkin(3) + t55;
t95 = sin(pkin(13));
t33 = t95 * t53 + t51;
t31 = qJD(3) * pkin(9) + t33;
t151 = t100 * t97 * t98;
t68 = -qJD(1) * t151 + t101 * t86 + qJD(4);
t22 = t104 * t68 + t107 * t31;
t20 = qJD(5) * pkin(10) + t22;
t121 = -pkin(5) * t107 - pkin(10) * t104 - pkin(4);
t50 = t95 * t56;
t32 = t99 * t53 - t50;
t25 = qJD(3) * t121 - t32;
t127 = t103 * t20 - t106 * t25;
t199 = (t145 - t177) * t170 + t79;
t48 = t199 * qJD(3);
t49 = t56 * qJD(3);
t28 = t95 * t48 + t99 * t49;
t131 = pkin(5) * t104 - pkin(10) * t107;
t84 = t131 * qJD(5);
t24 = qJD(3) * t84 + t28;
t125 = t104 * t31 - t107 * t68;
t29 = t48 * t99 - t49 * t95;
t7 = -qJD(5) * t125 + t107 * t29;
t1 = -qJD(6) * t127 + t103 * t24 + t106 * t7;
t160 = qJD(3) * t107;
t87 = -qJD(6) + t160;
t203 = -t127 * t87 + t1;
t34 = t55 * t95 + t51;
t202 = t34 - t84;
t6 = t103 * t25 + t106 * t20;
t2 = -qJD(6) * t6 - t103 * t7 + t106 * t24;
t201 = t6 * t87 - t2;
t73 = (t105 * t95 - t108 * t99) * t97;
t200 = qJD(3) * t34 - t28;
t120 = t102 * t97 + t98 * t166;
t60 = t105 * t120 + t98 * t173;
t93 = t104 ^ 2;
t123 = qJD(3) * t93 - t107 * t87;
t156 = qJD(6) * t103;
t142 = t104 * t156;
t154 = t106 * qJD(5);
t198 = -t123 * t154 - t87 * t142;
t155 = qJD(6) * t106;
t141 = t104 * t155;
t152 = qJD(5) * qJD(6);
t157 = qJD(5) * t107;
t67 = qJD(3) * (t103 * t157 + t141) + t103 * t152;
t197 = pkin(3) * t99;
t59 = t108 * t120 - t149;
t39 = t59 * t95 + t60 * t99;
t76 = t101 * t102 - t151;
t26 = t104 * t39 - t76 * t107;
t8 = t22 * qJD(5) + t104 * t29;
t196 = t26 * t8;
t74 = (t105 * t99 + t108 * t95) * t97;
t122 = t107 * t101 - t104 * t74;
t193 = t122 * t8;
t38 = -t59 * t99 + t60 * t95;
t192 = t28 * t38;
t191 = t28 * t73;
t190 = t8 * t103;
t189 = t8 * t106;
t161 = qJD(3) * t104;
t80 = t103 * t161 - t154;
t188 = t80 * t87;
t159 = qJD(5) * t103;
t82 = t106 * t161 + t159;
t187 = t82 * t80;
t186 = t82 * t87;
t158 = qJD(5) * t104;
t89 = pkin(3) * t95 + pkin(9);
t147 = t89 * t158;
t165 = t103 * t107;
t35 = t55 * t99 - t50;
t164 = t106 * t107;
t77 = t121 - t197;
t62 = t103 * t77 + t89 * t164;
t184 = qJD(6) * t62 - t103 * t147 + t202 * t106 - t35 * t165;
t61 = t106 * t77 - t89 * t165;
t183 = -qJD(6) * t61 + t202 * t103 + t106 * t147 + t35 * t164;
t143 = t107 * t154;
t171 = t67 * t106;
t182 = -t104 * t171 - t80 * t143;
t94 = t107 ^ 2;
t181 = t93 - t94;
t19 = -qJD(5) * pkin(5) + t125;
t180 = t103 * t19;
t179 = t103 * t87;
t178 = t104 * t80;
t176 = t106 * t19;
t175 = t106 * t87;
t174 = t107 * t67;
t110 = qJD(3) ^ 2;
t172 = t110 * t97;
t168 = qJD(3) * t35;
t167 = qJD(6) * t80;
t109 = qJD(5) ^ 2;
t163 = t109 * t104;
t162 = t109 * t107;
t153 = qJD(3) * qJD(5);
t146 = t82 * t157;
t144 = t104 * t110 * t107;
t139 = t104 * t153;
t66 = -t106 * t152 + (t142 - t143) * qJD(3);
t138 = t66 * t107 + t82 * t158;
t30 = -qJD(3) * pkin(4) - t32;
t137 = -qJD(3) * t30 - t29;
t136 = -t66 + t167;
t134 = t87 * t141;
t133 = t82 * t141;
t132 = t107 * t139;
t130 = -t103 * t6 + t106 * t127;
t129 = -t103 * t127 - t106 * t6;
t27 = t104 * t76 + t107 * t39;
t14 = t103 * t38 + t106 * t27;
t13 = -t103 * t27 + t106 * t38;
t65 = t101 * t104 + t107 * t74;
t45 = t103 * t73 + t106 * t65;
t44 = -t103 * t65 + t106 * t73;
t126 = -t104 * t125 - t107 * t22;
t119 = -t109 * t89 + t200;
t90 = -pkin(4) - t197;
t118 = qJD(5) * (qJD(3) * t90 + t30 + t35);
t115 = t123 * t103;
t113 = qJD(6) * t130 + t1 * t106 - t2 * t103;
t112 = t8 * t104 + t7 * t107 + (-t104 * t22 + t107 * t125) * qJD(5);
t83 = t131 * qJD(3);
t70 = qJD(3) * t73;
t69 = qJD(3) * t74;
t58 = t60 * qJD(3);
t57 = t59 * qJD(3);
t43 = qJD(5) * t65 - t104 * t70;
t42 = qJD(5) * t122 - t107 * t70;
t37 = t57 * t99 - t58 * t95;
t36 = t57 * t95 + t99 * t58;
t18 = t103 * t83 - t106 * t125;
t17 = t103 * t125 + t106 * t83;
t16 = -qJD(6) * t45 - t103 * t42 + t106 * t69;
t15 = qJD(6) * t44 + t103 * t69 + t106 * t42;
t10 = qJD(5) * t27 + t104 * t37;
t9 = -qJD(5) * t26 + t107 * t37;
t4 = -qJD(6) * t14 - t103 * t9 + t106 * t36;
t3 = qJD(6) * t13 + t103 * t36 + t106 * t9;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58 * qJD(3), -t57 * qJD(3), 0, t48 * t60 - t49 * t59 - t55 * t58 + t56 * t57, 0, 0, 0, 0, 0, 0, -t36 * qJD(3), -t37 * qJD(3), 0, t29 * t39 - t32 * t36 + t33 * t37 + t192, 0, 0, 0, 0, 0, 0, -t10 * qJD(5) + (-t107 * t36 + t158 * t38) * qJD(3), -t9 * qJD(5) + (t104 * t36 + t157 * t38) * qJD(3) (t10 * t104 + t107 * t9 + (-t104 * t27 + t107 * t26) * qJD(5)) * qJD(3), t10 * t125 + t22 * t9 + t27 * t7 + t30 * t36 + t192 + t196, 0, 0, 0, 0, 0, 0, t10 * t80 + t13 * t139 + t26 * t67 - t4 * t87, t10 * t82 - t139 * t14 - t26 * t66 + t3 * t87, t13 * t66 - t14 * t67 - t3 * t80 - t4 * t82, t1 * t14 + t10 * t19 - t127 * t4 + t13 * t2 + t3 * t6 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 * t172, -t108 * t172, 0 (t105 * t48 - t108 * t49 + (-t105 * t55 + t108 * t56) * qJD(3)) * t97, 0, 0, 0, 0, 0, 0, -t69 * qJD(3), t70 * qJD(3), 0, t29 * t74 - t32 * t69 - t33 * t70 + t191, 0, 0, 0, 0, 0, 0, -t43 * qJD(5) + (-t107 * t69 + t158 * t73) * qJD(3), -t42 * qJD(5) + (t104 * t69 + t157 * t73) * qJD(3) (t104 * t43 + t107 * t42 + (-t104 * t65 - t107 * t122) * qJD(5)) * qJD(3), t125 * t43 + t22 * t42 + t30 * t69 + t65 * t7 + t191 - t193, 0, 0, 0, 0, 0, 0, -t122 * t67 + t139 * t44 - t16 * t87 + t43 * t80, t122 * t66 - t139 * t45 + t15 * t87 + t43 * t82, -t15 * t80 - t16 * t82 + t44 * t66 - t45 * t67, t1 * t45 - t127 * t16 + t15 * t6 + t19 * t43 + t2 * t44 - t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t199 + t55) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, t200, -t29 + t168, 0, t32 * t34 - t33 * t35 + (-t28 * t99 + t29 * t95) * pkin(3), 0.2e1 * t132, -0.2e1 * t181 * t153, t162, -0.2e1 * t132, -t163, 0, t104 * t118 + t107 * t119, -t104 * t119 + t107 * t118 (-t93 - t94) * t168 + t112, t112 * t89 + t126 * t35 + t28 * t90 - t30 * t34, t82 * t143 + (-t66 * t106 - t156 * t82) * t104, -t133 + (-t146 + (t66 + t167) * t104) * t103 + t182, t138 - t198, t80 * t141 + (t104 * t67 + t157 * t80) * t103, t134 + t174 + (-t115 - t178) * qJD(5) (-t87 - t160) * t158, t184 * t87 + (-t2 + (t80 * t89 + t180) * qJD(5)) * t107 + (t19 * t155 + t190 - t35 * t80 + t67 * t89 + (qJD(3) * t61 - t127) * qJD(5)) * t104, -t183 * t87 + (t1 + (t82 * t89 + t176) * qJD(5)) * t107 + (-t19 * t156 + t189 - t35 * t82 - t66 * t89 + (-qJD(3) * t62 - t6) * qJD(5)) * t104, t61 * t66 - t62 * t67 + t184 * t82 + t183 * t80 + t130 * t157 + (qJD(6) * t129 - t1 * t103 - t106 * t2) * t104, t19 * t89 * t157 + t1 * t62 + t2 * t61 - t183 * t6 + t184 * t127 + (-t19 * t35 + t8 * t89) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t162, 0, -qJD(5) * t126 + t7 * t104 - t8 * t107, 0, 0, 0, 0, 0, 0, t134 - t174 + (-t115 + t178) * qJD(5), t138 + t198, t133 + (t104 * t136 + t146) * t103 + t182 (-qJD(5) * t129 - t8) * t107 + (qJD(5) * t19 + t113) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t181 * t110, 0, t144, 0, 0, t137 * t104, t137 * t107, 0, 0, -t66 * t103 - t175 * t82 (-t66 + t188) * t106 + (-t67 + t186) * t103, -t87 * t155 + (t87 * t164 + (-t82 + t159) * t104) * qJD(3), -t179 * t80 - t171, t87 * t156 + (-t87 * t165 + (t80 + t154) * t104) * qJD(3), t87 * t161, -pkin(5) * t67 - t189 + t17 * t87 - t22 * t80 + (pkin(10) * t175 + t180) * qJD(6) + (t104 * t127 + (-pkin(10) * t158 - t107 * t19) * t103) * qJD(3), pkin(5) * t66 + t190 - t18 * t87 - t22 * t82 + (-pkin(10) * t179 + t176) * qJD(6) + (-t19 * t164 + (-pkin(10) * t154 + t6) * t104) * qJD(3), t17 * t82 + t18 * t80 + ((qJD(6) * t82 - t67) * pkin(10) + t203) * t106 + (pkin(10) * t136 + t201) * t103, -t8 * pkin(5) + pkin(10) * t113 + t127 * t17 - t6 * t18 - t19 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, -t80 ^ 2 + t82 ^ 2, -t66 - t188, -t187, -t186 - t67, t139, -t19 * t82 - t201, t19 * t80 - t203, 0, 0;];
tauc_reg  = t5;
