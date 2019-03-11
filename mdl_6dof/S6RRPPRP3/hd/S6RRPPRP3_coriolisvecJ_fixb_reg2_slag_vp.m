% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:33
% EndTime: 2019-03-09 08:35:41
% DurationCPUTime: 2.65s
% Computational Cost: add. (3364->375), mult. (7193->481), div. (0->0), fcn. (3815->4), ass. (0->206)
t145 = sin(qJ(2));
t205 = t145 * qJD(1);
t115 = qJD(5) + t205;
t146 = cos(qJ(5));
t144 = sin(qJ(5));
t206 = t144 * qJD(2);
t147 = cos(qJ(2));
t212 = qJD(1) * t147;
t91 = t146 * t212 + t206;
t228 = t91 * t115;
t201 = qJD(1) * qJD(2);
t188 = t145 * t201;
t64 = qJD(5) * t91 - t144 * t188;
t251 = -t228 - t64;
t203 = t146 * qJD(2);
t90 = t144 * t212 - t203;
t167 = t90 * t115;
t207 = qJD(5) * t147;
t192 = t144 * t207;
t63 = -qJD(5) * t203 + (t145 * t203 + t192) * qJD(1);
t39 = t63 - t167;
t250 = t63 + t167;
t249 = t64 - t228;
t208 = qJD(5) * t146;
t209 = qJD(5) * t144;
t148 = -pkin(2) - pkin(3);
t137 = -pkin(8) + t148;
t157 = pkin(4) * t147 + t137 * t145;
t155 = t157 * qJD(2);
t117 = t147 * t201;
t128 = t145 * qJD(3);
t214 = qJ(3) * t117 + qJD(1) * t128;
t41 = qJD(1) * t155 + t214;
t174 = t145 * pkin(4) + t147 * pkin(8);
t86 = -qJD(1) * pkin(1) - pkin(2) * t212 - qJ(3) * t205;
t68 = pkin(3) * t212 + qJD(4) - t86;
t51 = qJD(1) * t174 + t68;
t119 = qJ(4) * t205;
t199 = pkin(7) * t205;
t173 = qJD(3) + t199;
t165 = -t119 + t173;
t66 = t137 * qJD(2) + t165;
t114 = pkin(7) * t117;
t210 = qJD(2) * t147;
t194 = qJ(4) * t210;
t204 = t145 * qJD(4);
t69 = t114 + (-t194 - t204) * qJD(1);
t184 = -t144 * t41 - t146 * t69 - t51 * t208 + t66 * t209;
t22 = -t144 * t66 + t146 * t51;
t248 = -t22 * t115 - t184;
t23 = t144 * t51 + t146 * t66;
t5 = -qJD(5) * t23 - t144 * t69 + t146 * t41;
t152 = -t63 * qJ(6) + t5;
t177 = pkin(5) * t117;
t1 = t91 * qJD(6) + t152 + t177;
t16 = t91 * qJ(6) + t22;
t13 = t115 * pkin(5) + t16;
t162 = -t64 * qJ(6) + t184;
t2 = t90 * qJD(6) - t162;
t17 = t90 * qJ(6) + t23;
t235 = t17 * t115;
t247 = -(t1 + t235) * t144 - (t115 * t13 - t2) * t146;
t245 = -t23 * t115 - t5;
t138 = qJD(2) * qJ(3);
t124 = pkin(7) * t212;
t95 = -qJ(4) * t212 + t124;
t83 = -t138 - t95;
t195 = qJD(2) * t148;
t110 = qJ(4) * t188;
t136 = qJD(2) * qJD(3);
t211 = qJD(2) * t145;
t159 = pkin(7) * t211 + t147 * qJD(4);
t61 = qJD(1) * t159 - t110 - t136;
t243 = t91 ^ 2;
t242 = t90 * t91;
t241 = pkin(7) - qJ(4);
t143 = qJ(3) + pkin(4);
t240 = t13 - t16;
t221 = t145 * t146;
t163 = pkin(5) * t147 - qJ(6) * t221;
t218 = qJ(6) - t137;
t178 = qJD(5) * t218;
t121 = qJ(3) * t212;
t55 = qJD(1) * t157 + t121;
t34 = -t144 * t95 + t146 * t55;
t239 = -qJD(1) * t163 + t144 * qJD(6) + t146 * t178 - t34;
t202 = t146 * qJD(6);
t223 = t144 * qJ(6);
t35 = t144 * t55 + t146 * t95;
t238 = -t144 * t178 - t205 * t223 + t202 + t35;
t100 = -t147 * pkin(2) - t145 * qJ(3) - pkin(1);
t87 = t147 * pkin(3) - t100;
t65 = t174 + t87;
t103 = t241 * t145;
t92 = t146 * t103;
t44 = t144 * t65 + t92;
t237 = qJD(2) * pkin(2);
t132 = t147 * qJ(4);
t104 = t147 * pkin(7) - t132;
t236 = t104 * t61;
t232 = t61 * t144;
t231 = t61 * t146;
t230 = t63 * t144;
t229 = t64 * t146;
t196 = pkin(5) * t144 - pkin(7);
t227 = -pkin(5) * t209 - t196 * t205 + qJD(3) - t119;
t226 = qJ(6) * t147;
t225 = t115 * t146;
t140 = t147 ^ 2;
t150 = qJD(1) ^ 2;
t224 = t140 * t150;
t222 = t144 * t145;
t220 = t147 * t150;
t149 = qJD(2) ^ 2;
t219 = t149 * t145;
t130 = t149 * t147;
t93 = -t119 + t199;
t217 = qJD(3) + t93;
t216 = -qJD(4) - t68;
t213 = qJ(3) * t210 + t128;
t50 = t155 + t213;
t80 = t241 * t210 - t204;
t200 = t144 * t50 + t146 * t80 + t65 * t208;
t198 = t115 * t222;
t197 = t115 * t221;
t193 = t115 * t209;
t191 = t115 * t208;
t190 = t146 * t207;
t189 = t115 * t212;
t99 = t173 - t237;
t187 = -t99 - t237;
t186 = -t90 * pkin(5) + qJD(6);
t185 = -t144 * t80 + t146 * t50;
t176 = t145 * t195;
t56 = qJD(1) * t176 + t214;
t67 = t176 + t213;
t183 = qJD(1) * t67 + t56;
t182 = qJD(1) * t87 + t68;
t43 = -t144 * t103 + t146 * t65;
t181 = t216 * t145;
t180 = -0.2e1 * pkin(1) * t201;
t179 = qJD(1) * t100 + t86;
t175 = t145 * t117;
t172 = -t13 * t144 + t146 * t17;
t171 = -t144 * t23 - t146 * t22;
t170 = -t144 * t22 + t146 * t23;
t166 = qJD(1) * t140 - t115 * t145;
t70 = pkin(2) * t188 - t214;
t81 = pkin(2) * t211 - t213;
t164 = -pkin(7) * t149 - qJD(1) * t81 - t70;
t31 = -t64 * pkin(5) - t61;
t77 = qJD(2) * pkin(4) - t83;
t49 = t186 + t77;
t161 = -t31 * t144 - t49 * t208;
t160 = -t31 * t146 + t49 * t209;
t158 = -t137 * t210 - t145 * t77;
t79 = -qJ(4) * t211 + t159;
t154 = qJD(5) * t172 + t1 * t146 + t2 * t144;
t153 = qJD(5) * t170 - t144 * t184 + t5 * t146;
t102 = t124 + t138;
t98 = -pkin(7) * t188 + t136;
t151 = t98 * t147 + (t147 * t99 + (-t102 + t124) * t145) * qJD(2);
t139 = t145 ^ 2;
t133 = 0.2e1 * t136;
t127 = t139 * t150;
t113 = t146 * pkin(5) + t143;
t112 = t145 * t220;
t109 = -t127 - t149;
t106 = -0.2e1 * t175;
t105 = 0.2e1 * t175;
t101 = -t127 + t224;
t97 = t218 * t146;
t96 = t218 * t144;
t94 = pkin(2) * t205 - t121;
t89 = (-t139 + t140) * t201;
t88 = t90 ^ 2;
t85 = 0.2e1 * t89;
t82 = -t147 * t196 - t132;
t78 = t148 * t205 + t121;
t74 = t195 + t165;
t71 = (t115 + t205) * t210;
t52 = (t145 * t206 - t190) * pkin(5) - t79;
t42 = -t88 + t243;
t36 = t147 * t223 + t44;
t33 = -t191 + qJD(2) * t90 + (-t147 * t206 - t197) * qJD(1);
t32 = t193 + qJD(2) * t91 + (-t147 * t203 + t198) * qJD(1);
t30 = t145 * pkin(5) + t146 * t226 + t43;
t29 = -t191 + (-t197 + (-t91 - t206) * t147) * qJD(1);
t28 = -t191 + (-t197 + (t91 - t206) * t147) * qJD(1);
t27 = -t193 + (-t198 + (-t90 + t203) * t147) * qJD(1);
t26 = t193 + (t198 + (-t90 - t203) * t147) * qJD(1);
t21 = t144 * t167 - t229;
t20 = t91 * t225 - t230;
t19 = t90 * t190 + (t147 * t64 - t90 * t211) * t144;
t18 = -t91 * t192 + (-t147 * t63 - t91 * t211) * t146;
t15 = t115 * t190 + t64 * t145 + (t144 * t166 + t147 * t90) * qJD(2);
t14 = t115 * t192 + t63 * t145 + (-t146 * t166 - t147 * t91) * qJD(2);
t12 = -t44 * qJD(5) + t185;
t11 = -t103 * t209 + t200;
t10 = qJ(6) * t190 + (-qJ(6) * t211 - qJD(5) * t103 + qJD(6) * t147) * t144 + t200;
t9 = t147 * t202 + t163 * qJD(2) + (-t92 + (-t65 - t226) * t144) * qJD(5) + t185;
t8 = t144 * t39 + t146 * t249;
t7 = t144 * t251 - t250 * t146;
t6 = t144 * t249 - t39 * t146;
t3 = (t144 * t91 + t146 * t90) * t211 + (t230 - t229 + (t144 * t90 - t146 * t91) * qJD(5)) * t147;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t85, t130, t106, -t219, 0, -pkin(7) * t130 + t145 * t180, pkin(7) * t219 + t147 * t180, 0, 0, t105, t130, -0.2e1 * t89, 0, t219, t106, t147 * t164 + t179 * t211, t151, t145 * t164 - t179 * t210, pkin(7) * t151 + t70 * t100 + t86 * t81, t106, t85, -t219, t105, t130, 0, t183 * t145 + (t147 * t182 - t79) * qJD(2), -t183 * t147 + (t145 * t182 + t80) * qJD(2), -t69 * t145 + t61 * t147 + (-t145 * t83 - t147 * t74) * qJD(2) + (-t145 * t80 + t147 * t79 + (-t103 * t147 + t104 * t145) * qJD(2)) * qJD(1), t103 * t69 + t56 * t87 + t67 * t68 + t74 * t80 + t79 * t83 - t236, t18, t3, t14, t19, t15, t71, -t104 * t64 + t12 * t115 + t79 * t90 + (t206 * t77 + t5) * t145 + (-t77 * t208 + t232 + (qJD(1) * t43 + t22) * qJD(2)) * t147, t104 * t63 - t11 * t115 + t79 * t91 + (t203 * t77 + t184) * t145 + (t77 * t209 + t231 + (-qJD(1) * t44 - t23) * qJD(2)) * t147, t11 * t90 + t12 * t91 + t147 * t153 + t171 * t211 - t43 * t63 + t44 * t64, t11 * t23 + t12 * t22 - t184 * t44 + t43 * t5 - t77 * t79 - t236, t18, t3, t14, t19, t15, t71, t9 * t115 - t52 * t90 - t82 * t64 + (t206 * t49 + t1) * t145 + ((qJD(1) * t30 + t13) * qJD(2) + t161) * t147, -t10 * t115 - t52 * t91 + t82 * t63 + (t203 * t49 - t2) * t145 + ((-qJD(1) * t36 - t17) * qJD(2) + t160) * t147, t10 * t90 - t30 * t63 + t36 * t64 + t9 * t91 + (-t13 * t146 - t144 * t17) * t211 + t154 * t147, t1 * t30 + t10 * t17 + t13 * t9 + t2 * t36 + t31 * t82 + t49 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t101, 0, t112, 0, 0, t150 * pkin(1) * t145, pkin(1) * t220, 0, 0, -t112, 0, t101, 0, 0, t112 (-t145 * t86 + t147 * t94) * qJD(1) ((t102 - t138) * t145 + (qJD(3) + t187) * t147) * qJD(1), t133 + (t145 * t94 + t147 * t86) * qJD(1), t98 * qJ(3) + t102 * qJD(3) - t86 * t94 + (t102 * t145 + t147 * t187) * qJD(1) * pkin(7), t112, -t101, 0, -t112, 0, 0, t93 * qJD(2) + t110 + t133 + (t216 * t147 + (-pkin(7) * qJD(2) - t78) * t145) * qJD(1), -t95 * qJD(2) + t114 + ((-qJ(4) * qJD(2) + t78) * t147 + t181) * qJD(1) (-t217 + t74 - t195) * t212, -t61 * qJ(3) + t69 * t148 - t217 * t83 - t68 * t78 - t74 * t95, t20, t7, t28, t21, t26, -t189, -t34 * t115 - t143 * t64 - t231 - t217 * t90 + (-t137 * t225 - t144 * t77) * qJD(5) + (t144 * t158 - t147 * t22) * qJD(1), t35 * t115 + t143 * t63 + t232 - t217 * t91 + (t115 * t137 * t144 - t146 * t77) * qJD(5) + (t146 * t158 + t147 * t23) * qJD(1), -t34 * t91 - t35 * t90 + (t22 * t205 + t137 * t64 + t184 + (-t137 * t91 + t22) * qJD(5)) * t146 + (t23 * t205 + t137 * t63 + t5 + (-t137 * t90 + t23) * qJD(5)) * t144, -t61 * t143 - t22 * t34 - t23 * t35 + t217 * t77 + (qJD(5) * t171 - t5 * t144 - t146 * t184) * t137, t20, t7, t28, t21, t26, -t189, -t113 * t64 - t227 * t90 + t239 * t115 + (-t49 * t222 + (qJD(2) * t96 - t13) * t147) * qJD(1) - t160, t113 * t63 - t227 * t91 + t238 * t115 + (-t49 * t221 + (qJD(2) * t97 + t17) * t147) * qJD(1) + t161, -t238 * t90 + t239 * t91 - t96 * t63 - t97 * t64 - t247, t1 * t96 + t113 * t31 + t239 * t13 - t238 * t17 - t2 * t97 + t227 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, 0, t109, -t102 * qJD(2) + t86 * t205 + t114, 0, 0, 0, 0, 0, 0, t109, t112, 0, t83 * qJD(2) + t114 + (t181 - t194) * qJD(1), 0, 0, 0, 0, 0, 0, t33, t32, t8, -t77 * qJD(2) + t245 * t144 + t146 * t248, 0, 0, 0, 0, 0, 0, t33, t32, t8, -t49 * qJD(2) + t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t117, 0.2e1 * t188, -t127 - t224 (-t147 * t83 + (t74 + t195) * t145) * qJD(1) + t214, 0, 0, 0, 0, 0, 0, t27, t29, t6 (t145 * t170 + t147 * t77) * qJD(1) + t153, 0, 0, 0, 0, 0, 0, t27, t29, t6 (t145 * t172 + t147 * t49) * qJD(1) + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t42, t39, -t242, t249, t117, t77 * t91 - t245, -t77 * t90 - t248, 0, 0, t242, t42, t39, -t242, t249, t117, 0.2e1 * t177 + t235 + (t186 + t49) * t91 + t152, -t243 * pkin(5) + t16 * t115 + (-qJD(6) - t49) * t90 + t162, -pkin(5) * t63 + t240 * t90, t240 * t17 + (t49 * t91 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t250, -t88 - t243, -t13 * t91 - t17 * t90 + t31;];
tauc_reg  = t4;
