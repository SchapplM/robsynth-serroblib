% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:09
% EndTime: 2019-03-09 01:47:16
% DurationCPUTime: 2.38s
% Computational Cost: add. (4668->325), mult. (9425->449), div. (0->0), fcn. (6065->8), ass. (0->185)
t114 = sin(pkin(9));
t120 = cos(qJ(4));
t199 = cos(pkin(10));
t162 = t199 * t120;
t113 = sin(pkin(10));
t118 = sin(qJ(4));
t197 = t113 * t118;
t131 = t162 - t197;
t115 = cos(pkin(9));
t185 = qJD(1) * t115;
t163 = t199 * t118;
t83 = t113 * t120 + t163;
t78 = t83 * qJD(4);
t212 = t114 * t78 + t131 * t185;
t117 = sin(qJ(6));
t119 = cos(qJ(6));
t184 = qJD(1) * t118;
t93 = qJD(1) * t162;
t77 = t113 * t184 - t93;
t191 = qJD(6) - t77;
t159 = t119 * t191;
t70 = qJD(1) * t78;
t236 = -t117 * t70 + t159 * t191;
t181 = qJD(4) * t118;
t109 = t120 * qJD(3);
t176 = qJD(1) * qJD(2);
t164 = t115 * t176;
t200 = qJD(4) * t109 + t120 * t164;
t121 = -pkin(1) - pkin(2);
t91 = t121 * qJD(1) + qJD(2);
t73 = qJ(2) * t185 + t114 * t91;
t66 = -qJD(1) * pkin(7) + t73;
t36 = -t66 * t181 + t200;
t26 = (qJ(5) * t181 - t120 * qJD(5)) * qJD(1) + t36;
t166 = t199 * t26;
t178 = t118 * qJD(3);
t50 = t120 * t66 + t178;
t130 = t50 * qJD(4);
t182 = qJD(2) * t115;
t154 = -qJD(5) + t182;
t134 = t154 * t118;
t168 = qJD(4) * t120 * qJ(5);
t231 = -t130 + (-t134 + t168) * qJD(1);
t10 = t231 * t113 + t166;
t175 = qJD(1) * qJD(4);
t165 = t118 * t175;
t89 = t113 * t165;
t71 = qJD(4) * t93 - t89;
t94 = t114 * t176;
t81 = -pkin(4) * t165 + t94;
t27 = -t70 * pkin(5) + t71 * pkin(8) + t81;
t46 = t178 + (-qJ(5) * qJD(1) + t66) * t120;
t40 = t199 * t46;
t206 = t118 * t66;
t49 = t109 - t206;
t45 = qJ(5) * t184 + t49;
t42 = qJD(4) * pkin(4) + t45;
t16 = t113 * t42 + t40;
t14 = qJD(4) * pkin(8) + t16;
t183 = qJD(1) * t120;
t186 = qJD(1) * t114;
t72 = -qJ(2) * t186 + t115 * t91;
t65 = qJD(1) * pkin(3) - t72;
t54 = pkin(4) * t183 + qJD(5) + t65;
t79 = t83 * qJD(1);
t21 = -t77 * pkin(5) + t79 * pkin(8) + t54;
t6 = t117 * t21 + t119 * t14;
t2 = -qJD(6) * t6 - t117 * t10 + t119 * t27;
t235 = t191 * t6 + t2;
t140 = t117 * t14 - t119 * t21;
t1 = -t140 * qJD(6) + t119 * t10 + t117 * t27;
t234 = t140 * t191 + t1;
t122 = qJD(4) ^ 2;
t123 = qJD(1) ^ 2;
t232 = t114 * (t122 + t123);
t80 = qJD(4) * t162 - t113 * t181;
t217 = t83 * t70;
t145 = -t191 * t80 + t217;
t179 = qJD(6) * t119;
t170 = t83 * t179;
t230 = -t117 * t145 + t170 * t191;
t37 = -t118 * t164 - t130;
t124 = -(t118 * t50 + t120 * t49) * qJD(4) - t37 * t118 + t36 * t120;
t229 = t79 ^ 2;
t228 = 0.2e1 * t94;
t190 = t115 * qJ(2) + t114 * t121;
t85 = -pkin(7) + t190;
t201 = qJ(5) - t85;
t67 = t201 * t120;
t29 = -t113 * t67 - t201 * t163;
t9 = t113 * t26 - t199 * t231;
t227 = t9 * t29;
t68 = t83 * t114;
t226 = t9 * t68;
t225 = t9 * t131;
t224 = t9 * t83;
t177 = t119 * qJD(4);
t55 = -t117 * t79 - t177;
t223 = t55 * t77;
t57 = t117 * qJD(4) - t119 * t79;
t222 = t57 * t55;
t221 = t57 * t79;
t220 = t70 * t131;
t219 = t77 * t79;
t218 = t79 * t55;
t69 = t131 * t114;
t51 = -t119 * t115 - t117 * t69;
t216 = t51 * qJD(6) - t117 * t186 - t212 * t119;
t136 = t117 * t115 - t119 * t69;
t215 = t136 * qJD(6) + t212 * t117 - t119 * t186;
t207 = t117 * t71;
t32 = t57 * qJD(6) - t207;
t214 = -t117 * t32 - t55 * t179;
t213 = -t80 * t114 + t115 * t79;
t211 = t113 * t46;
t210 = t117 * t55;
t209 = t117 * t57;
t205 = t119 * t55;
t204 = t119 * t57;
t180 = qJD(6) * t117;
t31 = -qJD(6) * t177 + t119 * t71 - t79 * t180;
t203 = t31 * t117;
t202 = t32 * t119;
t157 = -t114 * qJ(2) + t115 * t121;
t84 = pkin(3) - t157;
t198 = qJD(1) * t84;
t196 = t115 * t123;
t195 = t122 * t118;
t194 = t122 * t120;
t193 = t78 * qJD(4);
t192 = t80 * qJD(4);
t111 = t118 ^ 2;
t112 = t120 ^ 2;
t189 = t111 - t112;
t188 = t111 + t112;
t174 = pkin(4) * t184;
t173 = 0.2e1 * t176;
t172 = 0.2e1 * t175;
t171 = t83 * t180;
t169 = t118 * t123 * t120;
t161 = qJD(4) * t201;
t160 = t117 * t191;
t158 = -t65 - t182;
t155 = t115 * t172;
t98 = t113 * pkin(4) + pkin(8);
t153 = -qJD(6) * t191 * t98 - t9;
t152 = t120 * t165;
t87 = -pkin(4) * t181 + t114 * qJD(2);
t15 = t199 * t42 - t211;
t13 = -qJD(4) * pkin(5) - t15;
t151 = -t13 * t80 - t224;
t149 = t117 * t6 - t119 * t140;
t148 = t117 * t140 + t119 * t6;
t147 = t131 * t31 + t78 * t57;
t146 = -t131 * t32 + t78 * t55;
t144 = -t80 * t77 - t217;
t143 = t131 * t71 - t78 * t79;
t141 = t72 * t114 - t73 * t115;
t30 = t201 * t197 - t199 * t67;
t75 = t120 * pkin(4) + t84;
t34 = pkin(5) * t131 + t83 * pkin(8) + t75;
t12 = t117 * t34 + t119 * t30;
t11 = -t117 * t30 + t119 * t34;
t137 = t118 * t49 - t120 * t50;
t135 = -t119 * t70 + (t117 * t77 - t180) * t191;
t133 = t191 * t13 + t70 * t98;
t132 = -t122 * t85 + t228;
t129 = qJD(4) * (t158 - t198);
t128 = t120 * t161 - t134;
t127 = t145 * t119 + t171 * t191;
t125 = -t149 * qJD(6) + t1 * t119 - t2 * t117;
t99 = -t199 * pkin(4) - pkin(5);
t76 = t77 ^ 2;
t44 = t118 * t161 + t120 * t154;
t38 = -t79 * pkin(5) - t77 * pkin(8) - t174;
t33 = -t78 * pkin(5) + t80 * pkin(8) + t87;
t20 = t199 * t45 - t211;
t19 = t113 * t128 + t199 * t44;
t18 = t113 * t45 + t40;
t17 = t113 * t44 - t199 * t128;
t8 = t117 * t38 + t119 * t20;
t7 = -t117 * t20 + t119 * t38;
t4 = -t12 * qJD(6) - t117 * t19 + t119 * t33;
t3 = t11 * qJD(6) + t117 * t33 + t119 * t19;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, qJ(2) * t173, 0, 0, 0, 0, 0, 0, t228, 0.2e1 * t164, 0 ((-t114 * t157 + t115 * t190) * qJD(1) - t141) * qJD(2), 0.2e1 * t152, -t189 * t172, -t194, -0.2e1 * t152, t195, 0, t118 * t129 + t120 * t132, -t118 * t132 + t120 * t129, -t188 * t164 - t124, t124 * t85 + (-t137 * t115 + (t65 + t198) * t114) * qJD(2), t71 * t83 + t79 * t80, t143 + t144, -t192, t77 * t78 - t220, t193, 0, -t17 * qJD(4) + t131 * t81 - t54 * t78 - t75 * t70 - t87 * t77, -t19 * qJD(4) - t54 * t80 - t75 * t71 - t87 * t79 - t81 * t83, -t10 * t131 + t15 * t80 + t16 * t78 - t17 * t79 + t19 * t77 - t29 * t71 + t30 * t70 - t224, t10 * t30 - t15 * t17 + t16 * t19 + t54 * t87 + t81 * t75 + t227, t57 * t171 + (t31 * t83 - t57 * t80) * t119 (t205 + t209) * t80 + (-t203 + t202 + (t204 - t210) * qJD(6)) * t83, t127 - t147, -t55 * t170 + (-t32 * t83 - t55 * t80) * t117, t146 + t230, -t191 * t78 - t220, -t11 * t70 + t117 * t151 - t13 * t170 + t131 * t2 + t140 * t78 + t17 * t55 + t191 * t4 + t29 * t32, -t1 * t131 + t119 * t151 + t12 * t70 + t13 * t171 + t17 * t57 - t191 * t3 - t29 * t31 + t6 * t78, t11 * t31 - t12 * t32 - t3 * t55 - t4 * t57 + t149 * t80 + (qJD(6) * t148 + t1 * t117 + t2 * t119) * t83, t1 * t12 + t2 * t11 + t13 * t17 - t140 * t4 + t6 * t3 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t123 * qJ(2), 0, 0, 0, 0, 0, 0, -t114 * t123, -t196, 0, t141 * qJD(1), 0, 0, 0, 0, 0, 0, t118 * t155 - t120 * t232, t118 * t232 + t120 * t155, t188 * t196, t137 * t185 + (qJD(1) * t158 + t124) * t114, 0, 0, 0, 0, 0, 0, t213 * qJD(4) + t115 * t70 + t77 * t186, t212 * qJD(4) + t115 * t71 + t79 * t186, -t212 * t77 + t213 * t79 - t68 * t71 + t69 * t70, t10 * t69 - t81 * t115 + t213 * t15 - t212 * t16 - t54 * t186 + t226, 0, 0, 0, 0, 0, 0, t191 * t215 - t213 * t55 + t68 * t32 - t51 * t70, -t136 * t70 - t191 * t216 - t213 * t57 - t68 * t31, t136 * t32 - t215 * t57 - t216 * t55 + t51 * t31, -t1 * t136 - t213 * t13 - t140 * t215 + t2 * t51 + t216 * t6 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, -t194, 0, -qJD(4) * t137 + t36 * t118 + t37 * t120, 0, 0, 0, 0, 0, 0, -t193, -t192, t143 - t144, t10 * t83 - t15 * t78 + t16 * t80 - t225, 0, 0, 0, 0, 0, 0, t146 - t230, t127 + t147 (-t205 + t209) * t80 + (-t203 - t202 + (t204 + t210) * qJD(6)) * t83, t125 * t83 + t13 * t78 + t148 * t80 - t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t189 * t123, 0, t169, 0, 0 (t65 - t182) * t184, t65 * t183 + (t49 + t206) * qJD(4) - t200, 0, 0, t219, -t76 + t229, t89 + (-t93 - t77) * qJD(4), -t219, 0, 0, t18 * qJD(4) - t174 * t77 + t54 * t79 - t9, -t166 - t54 * t77 + (t113 * t50 + t20) * qJD(4) + (-t113 * t168 + (-pkin(4) * t79 + t113 * t154) * t118) * qJD(1) (-t16 + t18) * t79 + (t15 - t20) * t77 + (t113 * t70 + t199 * t71) * pkin(4), t15 * t18 - t16 * t20 + (t10 * t113 + t54 * t184 - t199 * t9) * pkin(4), t159 * t57 - t203 (-t31 + t223) * t119 - t191 * t209 + t214, t221 + t236, t160 * t55 - t202, t135 - t218, t191 * t79, t117 * t133 + t119 * t153 - t140 * t79 - t18 * t55 - t191 * t7 + t99 * t32, -t117 * t153 + t119 * t133 - t18 * t57 + t191 * t8 - t99 * t31 - t6 * t79, t8 * t55 + t7 * t57 + (-t32 * t98 - t140 * t77 + t1 + (t57 * t98 + t140) * qJD(6)) * t119 + (-t31 * t98 + t6 * t77 - t2 + (t55 * t98 - t6) * qJD(6)) * t117, t125 * t98 - t13 * t18 + t140 * t7 - t6 * t8 + t9 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t79 * qJD(4), t89 + (-t93 + t77) * qJD(4), -t76 - t229, -t15 * t79 - t16 * t77 + t81, 0, 0, 0, 0, 0, 0, t135 + t218, t221 - t236 (t31 + t223) * t119 + t57 * t160 + t214, t234 * t117 + t235 * t119 + t13 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t55 ^ 2 + t57 ^ 2, t191 * t55 - t31, -t222, t207 + (-qJD(6) + t191) * t57, -t70, -t13 * t57 + t235, t13 * t55 - t234, 0, 0;];
tauc_reg  = t5;
