% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:42
% EndTime: 2019-03-08 19:37:48
% DurationCPUTime: 2.62s
% Computational Cost: add. (1898->340), mult. (4466->460), div. (0->0), fcn. (3774->12), ass. (0->190)
t133 = sin(qJ(4));
t203 = qJD(2) * qJD(4);
t186 = t133 * t203;
t136 = cos(qJ(4));
t198 = t136 * qJDD(2);
t267 = t186 - t198;
t131 = cos(pkin(6));
t112 = t131 * qJD(1) + qJD(3);
t126 = sin(pkin(11));
t137 = cos(qJ(2));
t128 = sin(pkin(6));
t215 = qJD(1) * t128;
t189 = t137 * t215;
t91 = qJD(2) * pkin(2) + t189;
t129 = cos(pkin(11));
t134 = sin(qJ(2));
t192 = t134 * t215;
t96 = t129 * t192;
t53 = t126 * t91 + t96;
t51 = qJD(2) * pkin(8) + t53;
t240 = t136 * t112 - t133 * t51;
t266 = -qJD(5) + t240;
t214 = qJD(2) * t128;
t187 = qJD(1) * t214;
t202 = qJDD(1) * t128;
t268 = t134 * t202 + t137 * t187;
t26 = -qJD(4) * pkin(4) - t266;
t213 = qJD(2) * t136;
t33 = t133 * t112 + t136 * t51;
t29 = -qJD(4) * qJ(5) - t33;
t135 = cos(qJ(6));
t132 = sin(qJ(6));
t207 = t132 * qJD(4);
t83 = t135 * t213 + t207;
t233 = qJD(4) * t83;
t185 = t136 * t203;
t199 = t133 * qJDD(2);
t159 = t185 + t199;
t80 = qJDD(6) + t159;
t73 = t135 * t80;
t265 = t233 - t73;
t211 = qJD(4) * t133;
t59 = t126 * t189 + t96;
t264 = -pkin(4) * t211 + t59;
t127 = sin(pkin(10));
t130 = cos(pkin(10));
t223 = t131 * t137;
t263 = -t127 * t223 - t130 * t134;
t81 = t134 * t126 - t137 * t129;
t206 = t133 * qJD(2);
t219 = pkin(5) * t206 - t266;
t183 = -t133 * qJ(5) - pkin(3);
t170 = t136 * pkin(4) - t183;
t95 = t126 * t192;
t52 = t129 * t91 - t95;
t34 = -qJD(2) * t170 - t52;
t262 = t34 * t206 + qJDD(5);
t114 = qJD(6) + t206;
t261 = t114 - qJD(6);
t210 = qJD(4) * t136;
t108 = t137 * t202;
t63 = qJDD(2) * pkin(2) - t134 * t187 + t108;
t28 = t126 * t63 + t268 * t129;
t23 = qJDD(2) * pkin(8) + t28;
t260 = t112 * t211 + t133 * t23 + t51 * t210;
t160 = t81 * t131;
t175 = t137 * t126 + t134 * t129;
t36 = -t127 * t175 - t130 * t160;
t39 = t127 * t160 - t130 * t175;
t64 = t81 * t128;
t259 = -g(1) * t39 - g(2) * t36 + g(3) * t64;
t227 = t128 * t136;
t66 = t175 * t131;
t35 = t127 * t81 - t130 * t66;
t17 = t130 * t227 - t133 * t35;
t38 = t127 * t66 + t130 * t81;
t19 = -t127 * t227 - t133 * t38;
t65 = t175 * t128;
t48 = -t131 * t136 + t65 * t133;
t258 = g(1) * t19 + g(2) * t17 + g(3) * t48;
t119 = pkin(5) * t213;
t138 = -pkin(4) - pkin(9);
t21 = t119 - t29;
t257 = t138 * t80 + (t21 - t119 - t33) * t114;
t208 = qJD(6) * t136;
t190 = t135 * t208;
t239 = t132 * t80;
t256 = -t114 * (t133 * t207 - t190) + t136 * t239;
t188 = t132 * t213;
t43 = -qJD(6) * t188 + t132 * qJDD(4) + (qJD(4) * qJD(6) - t267) * t135;
t117 = t126 * pkin(2) + pkin(8);
t139 = qJD(4) ^ 2;
t152 = t117 * t139 - t259;
t205 = t133 * qJD(5);
t166 = -qJ(5) * t210 - t205;
t242 = t166 - t264;
t246 = t129 * pkin(2);
t75 = -t170 - t246;
t27 = -t268 * t126 + t129 * t63;
t172 = pkin(4) * t186 - t27;
t9 = qJD(2) * t166 - qJDD(2) * t170 + t172;
t255 = qJD(2) * t242 + qJDD(2) * t75 + t152 + t9;
t61 = t81 * t214;
t11 = -t65 * t211 + (qJD(4) * t131 - t61) * t136;
t49 = t131 * t133 + t65 * t136;
t60 = qJD(2) * t65;
t254 = qJD(2) * (t133 * t60 + t210 * t64) - t11 * qJD(4) - t49 * qJDD(4) + t199 * t64;
t10 = qJD(4) * t49 - t61 * t133;
t253 = qJD(2) * (-t136 * t60 + t211 * t64) - t10 * qJD(4) - t48 * qJDD(4) - t198 * t64;
t245 = pkin(5) + t117;
t42 = -qJD(6) * t83 + t135 * qJDD(4) + t267 * t132;
t204 = t135 * qJD(4);
t85 = -t188 + t204;
t244 = t42 * t133 + t85 * t210;
t234 = qJ(5) * t136;
t179 = pkin(9) * t133 - t234;
t148 = qJD(4) * t179 - t205;
t243 = -t148 + t264;
t209 = qJD(6) * t132;
t191 = t114 * t209;
t230 = t114 * t135;
t194 = t133 * t230;
t241 = qJD(4) * t194 + t136 * t191;
t238 = t133 * t43;
t237 = t42 * t135;
t236 = t83 * t114;
t235 = t85 * t114;
t232 = qJDD(4) * pkin(4);
t231 = t114 * t132;
t229 = t127 * t134;
t228 = t128 * t133;
t226 = t128 * t137;
t224 = t131 * t134;
t220 = t33 * qJD(4);
t218 = qJDD(1) - g(3);
t124 = t133 ^ 2;
t125 = t136 ^ 2;
t217 = t124 - t125;
t201 = qJDD(4) * qJ(5);
t200 = qJDD(4) * t117;
t109 = t131 * qJDD(1) + qJDD(3);
t195 = t130 * t223;
t140 = qJD(2) ^ 2;
t193 = t133 * t140 * t136;
t77 = t245 * t136;
t182 = -t133 * t109 - t112 * t210 - t136 * t23 + t51 * t211;
t12 = t138 * qJD(4) + t219;
t158 = t138 * t136 + t183;
t30 = qJD(2) * t158 - t52;
t6 = t135 * t12 - t132 * t30;
t7 = t132 * t12 + t135 * t30;
t13 = -t64 * t132 + t48 * t135;
t14 = t48 * t132 + t64 * t135;
t177 = t26 * t133 - t29 * t136;
t173 = -t136 * t109 + t260;
t18 = -t130 * t228 - t136 * t35;
t20 = t127 * t228 - t136 * t38;
t165 = -g(1) * t20 - g(2) * t18 - g(3) * t49;
t164 = g(1) * t38 + g(2) * t35 - g(3) * t65;
t162 = qJDD(5) + t173;
t161 = -g(3) * t131 + (-g(1) * t127 + g(2) * t130) * t128;
t8 = qJD(2) * t148 + qJDD(2) * t158 + t172;
t157 = t259 - t8;
t76 = t245 * t133;
t155 = t76 * t80 + t164;
t151 = qJD(4) * qJD(5) - t182 + t201;
t62 = t129 * t189 - t95;
t150 = t200 + (-qJD(2) * t75 - t34 - t62) * qJD(4);
t118 = -pkin(3) - t246;
t50 = -qJD(2) * pkin(3) - t52;
t149 = -t200 + (qJD(2) * t118 + t50 + t62) * qJD(4);
t147 = t165 - t182;
t146 = -g(1) * t263 - g(3) * t226;
t145 = -t173 + t258;
t5 = t162 - t232;
t144 = t5 * t133 + t151 * t136 + (t133 * t29 + t136 * t26) * qJD(4);
t121 = pkin(4) * t206;
t3 = -t267 * pkin(5) + t151;
t142 = t3 + (qJD(2) * t179 - qJD(6) * t138 + t121) * t114 + t165;
t141 = -qJD(2) * t59 + t152 - t27 + (-pkin(3) + t118) * qJDD(2);
t101 = pkin(2) * t195;
t100 = qJDD(4) * t136 - t139 * t133;
t99 = qJDD(4) * t133 + t139 * t136;
t86 = -qJ(5) * t213 + t121;
t71 = qJD(4) * t77;
t70 = t245 * t211;
t67 = t158 - t246;
t2 = t159 * pkin(5) + t138 * qJDD(4) + t162;
t1 = t135 * t2;
t4 = [t218, 0 (qJDD(2) * t137 - t134 * t140) * t128 (-qJDD(2) * t134 - t137 * t140) * t128, t109 * t131 - t27 * t64 + t28 * t65 - t52 * t60 - t53 * t61 - g(3), 0, 0, 0, 0, 0, t253, t254 (t133 * t48 + t136 * t49) * qJDD(2) + (t10 * t133 + t11 * t136 + (-t133 * t49 + t136 * t48) * qJD(4)) * qJD(2), -t253, -t254, t10 * t26 - t11 * t29 + t151 * t49 + t34 * t60 + t48 * t5 + t64 * t9 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t14 + t10 * t135 - t60 * t132) * t114 + t13 * t80 + t11 * t83 + t49 * t43 -(qJD(6) * t13 + t10 * t132 + t60 * t135) * t114 - t14 * t80 + t11 * t85 + t49 * t42; 0, qJDD(2), t108 - g(2) * (t195 - t229) + t146, -g(1) * (t127 * t224 - t130 * t137) - g(2) * (-t127 * t137 - t130 * t224) - t218 * t134 * t128, -g(2) * t101 + t52 * t59 - t53 * t62 + (g(2) * t229 + t28 * t126 + t27 * t129 + t146) * pkin(2), t124 * qJDD(2) + 0.2e1 * t133 * t185, 0.2e1 * t133 * t198 - 0.2e1 * t203 * t217, t99, t100, 0, t133 * t149 - t136 * t141, t133 * t141 + t136 * t149, t144 + t164 + (-qJD(2) * t62 + qJDD(2) * t117) * (t124 + t125) t150 * t133 + t255 * t136, -t255 * t133 + t150 * t136, t9 * t75 - g(1) * (t263 * pkin(2) - t38 * pkin(8)) - g(2) * (-pkin(2) * t229 - t35 * pkin(8) + t101) - g(3) * (pkin(2) * t226 + t65 * pkin(8)) - t177 * t62 + t242 * t34 + t144 * t117 + t259 * t170, -t85 * t190 + (-t136 * t42 + t211 * t85) * t132 (-t132 * t83 + t135 * t85) * t211 + (t132 * t43 - t237 + (t132 * t85 + t135 * t83) * qJD(6)) * t136, t244 - t256, -t238 + (-t233 - t73) * t136 + t241, t114 * t210 + t80 * t133, -t67 * t239 + t77 * t43 - t70 * t83 + t155 * t135 + (t132 * t157 - t204 * t21 + t1) * t133 + ((-t133 * t62 + t71) * t135 + t243 * t132) * t114 + ((-t132 * t76 - t135 * t67) * t114 - t7 * t133) * qJD(6) + (t6 * qJD(4) + t3 * t135 - t209 * t21 - t62 * t83) * t136, t77 * t42 - t70 * t85 + (-t7 * qJD(4) - t62 * t85) * t136 + (-t21 * t208 - t67 * t80 + (-qJD(6) * t76 + t243) * t114 + (-qJD(6) * t12 + t157) * t133) * t135 + (-(-qJD(6) * t67 + t71) * t114 - t3 * t136 + (t21 * qJD(4) + qJD(6) * t30 + t62 * t114 - t2) * t133 - t155) * t132; 0, 0, 0, 0, t161 + t109, 0, 0, 0, 0, 0, t100, -t99, 0, -t100, t99, qJD(4) * t177 + t133 * t151 - t5 * t136 + t161, 0, 0, 0, 0, 0, t265 * t136 + t238 + t241, t244 + t256; 0, 0, 0, 0, 0, -t193, t217 * t140, t199, t198, qJDD(4), -t206 * t50 + t145 + t220, qJD(4) * t240 - t213 * t50 - t147 (-pkin(4) * t133 + t234) * qJDD(2), -0.2e1 * t232 - t220 + (-qJD(2) * t86 - t109) * t136 - t258 + t260 + t262, 0.2e1 * t201 + (0.2e1 * qJD(5) - t240) * qJD(4) + (t133 * t86 + t136 * t34) * qJD(2) + t147, t151 * qJ(5) - t5 * pkin(4) - t34 * t86 - t26 * t33 - g(1) * (-pkin(4) * t19 + qJ(5) * t20) - g(2) * (-pkin(4) * t17 + qJ(5) * t18) - g(3) * (-pkin(4) * t48 + qJ(5) * t49) + t266 * t29, -t231 * t85 + t237 (-t43 - t235) * t135 + (-t42 + t236) * t132, -t191 + t73 + (-t133 * t231 - t136 * t85) * qJD(2), -qJD(6) * t230 - t239 + (t136 * t83 - t194) * qJD(2), -t114 * t213, qJ(5) * t43 + t142 * t132 + t257 * t135 - t6 * t213 + t219 * t83, qJ(5) * t42 - t257 * t132 + t142 * t135 + t7 * t213 + t219 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, qJDD(4) + t193, -t124 * t140 - t139, t29 * qJD(4) - t145 - t232 + t262, 0, 0, 0, 0, 0, -t114 * t231 - t265, -t114 ^ 2 * t135 - qJD(4) * t85 - t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t83, -t83 ^ 2 + t85 ^ 2, t42 + t236, t235 - t43, t80, -t132 * t8 + t1 - t21 * t85 - g(1) * (t39 * t132 + t19 * t135) - g(2) * (t36 * t132 + t17 * t135) - g(3) * t13 + t261 * t7, -t135 * t8 - t132 * t2 + t21 * t83 - g(1) * (-t19 * t132 + t39 * t135) - g(2) * (-t17 * t132 + t36 * t135) + g(3) * t14 + t261 * t6;];
tau_reg  = t4;
