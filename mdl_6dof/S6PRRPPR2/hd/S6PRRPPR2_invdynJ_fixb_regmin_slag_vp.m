% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:33
% EndTime: 2019-03-08 21:06:43
% DurationCPUTime: 3.72s
% Computational Cost: add. (3033->402), mult. (7223->538), div. (0->0), fcn. (5749->14), ass. (0->216)
t142 = sin(pkin(11));
t145 = cos(pkin(11));
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t106 = t142 * t152 + t145 * t149;
t100 = t106 * qJD(2);
t282 = qJD(6) + t100;
t148 = sin(qJ(6));
t283 = t148 * t282;
t151 = cos(qJ(6));
t217 = qJD(2) * qJD(3);
t203 = t152 * t217;
t204 = t149 * t217;
t54 = qJDD(2) * t106 - t142 * t204 + t145 * t203;
t49 = qJDD(6) + t54;
t46 = t151 * t49;
t177 = -t282 * t283 + t46;
t224 = qJD(2) * t149;
t233 = t145 * t152;
t97 = -qJD(2) * t233 + t142 * t224;
t72 = qJD(3) * t148 - t151 * t97;
t284 = t282 * t72;
t261 = qJ(4) + pkin(8);
t146 = cos(pkin(6));
t144 = sin(pkin(6));
t150 = sin(qJ(2));
t236 = t144 * t150;
t103 = t146 * t152 - t149 * t236;
t153 = cos(qJ(2));
t234 = t144 * t153;
t206 = qJD(1) * t234;
t200 = qJD(3) * t261;
t90 = t152 * qJD(4) - t149 * t200;
t91 = -t149 * qJD(4) - t152 * t200;
t258 = -t106 * t206 + t142 * t90 - t145 * t91;
t277 = -t142 * t149 + t233;
t257 = t142 * t91 + t145 * t90 - t277 * t206;
t92 = t100 ^ 2;
t281 = -t97 ^ 2 - t92;
t225 = qJD(1) * t150;
t209 = t144 * t225;
t194 = t261 * qJD(2) + t209;
t226 = qJD(1) * t146;
t67 = t149 * t226 + t194 * t152;
t59 = t142 * t67;
t66 = -t194 * t149 + t152 * t226;
t32 = t145 * t66 - t59;
t280 = -qJD(5) + t32;
t230 = t148 * t153;
t119 = t144 * t230;
t235 = t144 * t152;
t104 = t146 * t149 + t150 * t235;
t44 = -t145 * t103 + t104 * t142;
t279 = t151 * t44 + t119;
t244 = cos(pkin(10));
t197 = t244 * t153;
t143 = sin(pkin(10));
t238 = t143 * t150;
t93 = -t146 * t197 + t238;
t198 = t244 * t150;
t237 = t143 * t153;
t95 = t146 * t237 + t198;
t191 = g(1) * t95 + g(2) * t93;
t278 = -g(3) * t234 + t191;
t138 = qJ(3) + pkin(11);
t135 = sin(t138);
t136 = cos(t138);
t276 = pkin(4) * t136 + qJ(5) * t135;
t275 = t282 - qJD(6);
t267 = pkin(3) * t145;
t131 = -pkin(4) - t267;
t127 = -pkin(9) + t131;
t251 = t145 * t67;
t63 = qJD(3) * pkin(3) + t66;
t28 = t142 * t63 + t251;
t23 = -qJD(3) * qJ(5) - t28;
t269 = pkin(5) * t97;
t13 = -t23 - t269;
t30 = t142 * t66 + t251;
t274 = t127 * t49 + (t13 - t30 + t269) * t282;
t154 = qJD(3) ^ 2;
t201 = qJDD(1) * t234;
t218 = qJD(1) * qJD(2);
t205 = t150 * t218;
t184 = t144 * t205 - t201;
t273 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t154 + t144 * (-g(3) * t153 + t205) - t184 + t191;
t271 = pkin(4) + pkin(9);
t214 = t152 * qJDD(2);
t215 = t149 * qJDD(2);
t183 = -t142 * t215 + t145 * t214;
t99 = t106 * qJD(3);
t53 = qJD(2) * t99 - t183;
t270 = pkin(4) * t53;
t268 = pkin(3) * t142;
t265 = t100 * pkin(5);
t132 = t152 * pkin(3) + pkin(2);
t179 = -t106 * qJ(5) - t132;
t38 = -t271 * t277 + t179;
t264 = t38 * t49;
t74 = qJD(3) * t151 + t148 * t97;
t263 = t74 * t97;
t262 = t97 * t72;
t216 = t146 * qJDD(1);
t121 = t152 * t216;
t78 = qJDD(2) * pkin(8) + (qJDD(1) * t150 + t153 * t218) * t144;
t164 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t226 + t78;
t181 = t194 * qJD(3);
t19 = qJDD(3) * pkin(3) - t164 * t149 - t152 * t181 + t121;
t20 = (-t181 + t216) * t149 + t164 * t152;
t9 = t142 * t19 + t145 * t20;
t222 = qJD(3) * t149;
t102 = qJD(3) * t233 - t142 * t222;
t260 = pkin(5) * t102 + t258;
t259 = -pkin(5) * t99 + t257;
t94 = t146 * t198 + t237;
t256 = -t93 * t132 + t261 * t94;
t96 = -t146 * t238 + t197;
t255 = -t95 * t132 + t261 * t96;
t254 = qJ(5) * t54;
t253 = qJD(2) * pkin(2);
t252 = t277 * t13;
t250 = t148 * t49;
t248 = t148 * t99;
t246 = t151 * t282;
t220 = qJD(6) * t151;
t221 = qJD(6) * t148;
t25 = -qJD(3) * t221 + t151 * qJDD(3) + t148 * t53 + t97 * t220;
t245 = t25 * t151;
t241 = qJDD(3) * pkin(4);
t239 = t143 * t144;
t231 = t261 * t150;
t229 = t265 - t280;
t228 = qJDD(1) - g(3);
t140 = t149 ^ 2;
t227 = -t152 ^ 2 + t140;
t223 = qJD(2) * t150;
t213 = g(3) * t236;
t212 = qJDD(3) * qJ(5) + t9;
t134 = pkin(3) * t222;
t133 = pkin(3) * t224;
t210 = t151 * t234;
t208 = t144 * t223;
t207 = qJD(2) * t234;
t202 = t153 * t217;
t8 = -t142 * t20 + t145 * t19;
t27 = t145 * t63 - t59;
t199 = t144 * t244;
t196 = t97 * qJ(5) + t133;
t113 = t261 * t149;
t114 = t261 * t152;
t70 = t145 * t113 + t114 * t142;
t193 = t148 * qJDD(3) - t151 * t53;
t192 = (t143 * t235 - t149 * t96) * pkin(3);
t190 = g(1) * t96 + g(2) * t94;
t189 = qJD(5) - t27;
t173 = -t102 * qJ(5) - t106 * qJD(5) + t134;
t36 = pkin(4) * t99 + t173;
t188 = -t36 + t209;
t187 = qJDD(5) - t8;
t12 = -t271 * qJD(3) + t189 + t265;
t87 = -t132 * qJD(2) + qJD(4) - t206;
t161 = -qJ(5) * t100 + t87;
t24 = t271 * t97 + t161;
t5 = t12 * t151 - t148 * t24;
t6 = t12 * t148 + t151 * t24;
t71 = -t113 * t142 + t114 * t145;
t182 = t103 * pkin(3);
t155 = qJD(2) ^ 2;
t180 = qJDD(2) * t153 - t150 * t155;
t4 = -qJD(3) * qJD(5) - t212;
t176 = -g(1) * t143 + t244 * g(2);
t175 = -t148 * t44 + t210;
t56 = -t135 * t199 + t94 * t136;
t58 = t135 * t239 + t136 * t96;
t82 = t135 * t146 + t136 * t236;
t172 = -g(1) * t58 - g(2) * t56 - g(3) * t82;
t171 = -t246 * t282 - t250;
t111 = -t206 - t253;
t169 = -qJD(2) * t111 + t190 - t78;
t168 = (-t94 * t149 - t152 * t199) * pkin(3);
t167 = t278 * t135;
t64 = qJD(3) * t103 + t152 * t207;
t65 = -qJD(3) * t104 - t149 * t207;
t29 = t142 * t64 - t145 * t65;
t31 = t142 * t65 + t145 * t64;
t45 = t103 * t142 + t104 * t145;
t166 = t100 * t29 - t31 * t97 + t44 * t54 - t45 * t53;
t3 = -pkin(5) * t53 - t4;
t42 = pkin(5) * t106 + t70;
t163 = -t13 * t99 + t277 * t3 + t42 * t49 - t190;
t162 = -pkin(8) * qJDD(3) + (t111 + t206 - t253) * qJD(3);
t52 = pkin(3) * t204 - t132 * qJDD(2) + qJDD(4) + t184;
t160 = t3 + (-qJD(6) * t127 + t271 * t100 + t196) * t282 + t172;
t37 = pkin(4) * t97 + t161;
t55 = t94 * t135 + t136 * t199;
t57 = t135 * t96 - t136 * t239;
t81 = t135 * t236 - t136 * t146;
t159 = -g(1) * t57 - g(2) * t55 - g(3) * t81 + t100 * t37 + t187;
t158 = -qJD(5) * t100 - t254 + t52;
t157 = t52 - t278;
t156 = t258 * t100 - t257 * t97 - t53 * t71 + t54 * t70 - t190 - t213;
t129 = qJ(5) + t268;
t109 = t132 * t234;
t88 = qJD(3) * t97;
t50 = -pkin(4) * t277 + t179;
t43 = pkin(5) * t277 + t71;
t39 = pkin(4) * t100 + t196;
t26 = t74 * qJD(6) + t193;
t22 = -qJD(3) * pkin(4) + t189;
t21 = t271 * t99 + t173;
t11 = t158 + t270;
t10 = t271 * t53 + t158;
t7 = t187 - t241;
t2 = t54 * pkin(5) - t271 * qJDD(3) + t187;
t1 = t151 * t2;
t14 = [t228, 0, t180 * t144 (-qJDD(2) * t150 - t153 * t155) * t144, 0, 0, 0, 0, 0, t65 * qJD(3) + t103 * qJDD(3) + (-t149 * t202 + t152 * t180) * t144, -t64 * qJD(3) - t104 * qJDD(3) + (-t149 * t180 - t152 * t202) * t144, t166, -t27 * t29 + t28 * t31 - t44 * t8 + t45 * t9 - g(3) + (-t153 * t52 + t87 * t223) * t144, t166, t29 * qJD(3) + t44 * qJDD(3) + (t153 * t53 - t97 * t223) * t144, t31 * qJD(3) + t45 * qJDD(3) + (-t100 * t223 + t153 * t54) * t144, t22 * t29 - t23 * t31 - t4 * t45 + t44 * t7 - g(3) + (-t11 * t153 + t37 * t223) * t144, 0, 0, 0, 0, 0 (qJD(6) * t175 - t148 * t208 + t151 * t29) * t282 + t279 * t49 + t31 * t72 + t45 * t26 -(qJD(6) * t279 + t148 * t29 + t151 * t208) * t282 + t175 * t49 + t31 * t74 + t45 * t25; 0, qJDD(2), t201 + t278, -t228 * t236 + t190, qJDD(2) * t140 + 0.2e1 * t149 * t203, 0.2e1 * t149 * t214 - 0.2e1 * t227 * t217, qJDD(3) * t149 + t152 * t154, qJDD(3) * t152 - t149 * t154, 0, t162 * t149 + t152 * t273, -t149 * t273 + t162 * t152, -t102 * t27 - t106 * t8 + t277 * t9 - t28 * t99 + t156, t9 * t71 - t8 * t70 - t52 * t132 - g(1) * t255 - g(2) * t256 - g(3) * (t144 * t231 + t109) + (-t209 + t134) * t87 + t257 * t28 - t258 * t27, t102 * t22 + t106 * t7 + t23 * t99 - t277 * t4 + t156, t258 * qJD(3) + qJDD(3) * t70 + t11 * t277 - t136 * t278 + t188 * t97 - t37 * t99 - t50 * t53, t257 * qJD(3) + qJDD(3) * t71 + t188 * t100 - t102 * t37 - t106 * t11 - t50 * t54 + t167, t11 * t50 + t37 * t36 - t4 * t71 + t7 * t70 - g(1) * (-t276 * t95 + t255) - g(2) * (-t276 * t93 + t256) - g(3) * t109 - t257 * t23 + t258 * t22 + (-t37 * t225 - g(3) * (t153 * t276 + t231)) * t144, t74 * t248 - (t25 * t148 + t74 * t220) * t277 (-t148 * t72 + t151 * t74) * t99 - (-t148 * t26 + t245 + (-t148 * t74 - t151 * t72) * qJD(6)) * t277, t282 * t248 + t74 * t102 + t25 * t106 - (t220 * t282 + t250) * t277, t99 * t246 - t72 * t102 - t26 * t106 - (-t221 * t282 + t46) * t277, t102 * t282 + t106 * t49, t1 * t106 + t5 * t102 + t43 * t26 + t259 * t72 + (-t10 * t106 + t135 * t191 - t21 * t282 - t264) * t148 + (t260 * t282 + t163) * t151 + (t225 * t283 - g(3) * (t135 * t230 + t150 * t151)) * t144 + ((-t148 * t42 - t151 * t38) * t282 - t6 * t106 - t148 * t252) * qJD(6), -t6 * t102 + t43 * t25 + t259 * t74 + (-t264 - (qJD(6) * t12 + t10) * t106 - qJD(6) * t252 + (-qJD(6) * t42 + t209 - t21) * t282 + t167) * t151 + (-(-qJD(6) * t24 + t2) * t106 + t213 + (qJD(6) * t38 - t260) * t282 - t163) * t148; 0, 0, 0, 0, -t149 * t155 * t152, t227 * t155, t215, t214, qJDD(3), -g(3) * t103 + t169 * t149 + t176 * t235 + t121, g(3) * t104 + (-t144 * t176 - t216) * t149 + t169 * t152 (-t27 + t32) * t97 + (t28 - t30) * t100 + (-t142 * t53 - t145 * t54) * pkin(3), -g(1) * t192 - g(2) * t168 - g(3) * t182 - t133 * t87 + t8 * t267 + t9 * t268 + t27 * t30 - t28 * t32, -t129 * t53 + t131 * t54 + (-t23 - t30) * t100 + (t22 + t280) * t97, -qJD(3) * t30 + t39 * t97 + (-pkin(4) + t131) * qJDD(3) + t159, qJDD(3) * t129 + t100 * t39 - t37 * t97 + (0.2e1 * qJD(5) - t32) * qJD(3) + t172 + t212, -t4 * t129 + t7 * t131 - t37 * t39 - t22 * t30 - g(1) * (-pkin(4) * t57 + qJ(5) * t58 + t192) - g(2) * (-t55 * pkin(4) + t56 * qJ(5) + t168) - g(3) * (-pkin(4) * t81 + qJ(5) * t82 + t182) + t280 * t23, -t283 * t74 + t245 (-t282 * t74 - t26) * t151 + (-t25 + t284) * t148, t177 + t263, t171 - t262, t282 * t97, t129 * t26 + t160 * t148 + t151 * t274 + t229 * t72 + t5 * t97, t129 * t25 - t148 * t274 + t160 * t151 + t229 * t74 - t6 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, t100 * t27 + t28 * t97 + t157, t281, -0.2e1 * t100 * qJD(3) + t183, -t54 + t88, t270 - t254 - t23 * t97 + (-qJD(5) - t22) * t100 + t157, 0, 0, 0, 0, 0, t171 + t262, -t177 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 + t54, -t100 * t97 + qJDD(3), -t92 - t154, qJD(3) * t23 + t159 - t241, 0, 0, 0, 0, 0, -qJD(3) * t72 + t177, -qJD(3) * t74 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t72 ^ 2 + t74 ^ 2, t25 + t284, t275 * t74 - t193, t49, -t148 * t10 + t1 - t13 * t74 - g(1) * (-t148 * t95 + t151 * t57) - g(2) * (-t148 * t93 + t151 * t55) - g(3) * (t151 * t81 + t119) + t275 * t6, -t151 * t10 - t148 * t2 + t13 * t72 - g(1) * (-t148 * t57 - t151 * t95) - g(2) * (-t148 * t55 - t151 * t93) - g(3) * (-t148 * t81 + t210) + t275 * t5;];
tau_reg  = t14;
