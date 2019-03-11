% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:15
% EndTime: 2019-03-08 21:02:26
% DurationCPUTime: 4.95s
% Computational Cost: add. (4861->457), mult. (11686->656), div. (0->0), fcn. (9474->18), ass. (0->235)
t202 = sin(qJ(3));
t203 = sin(qJ(2));
t197 = sin(pkin(6));
t275 = qJD(1) * t197;
t256 = t203 * t275;
t309 = qJD(3) * pkin(3);
t328 = t202 * t309 - t256;
t195 = sin(pkin(11));
t205 = cos(qJ(3));
t299 = cos(pkin(11));
t156 = t195 * t205 + t299 * t202;
t144 = t156 * qJD(3);
t242 = t299 * t205;
t218 = -t195 * t202 + t242;
t147 = t218 * qJD(3);
t327 = pkin(4) * t144 - qJ(5) * t147 - qJD(5) * t156 + t328;
t315 = qJ(4) + pkin(8);
t247 = qJD(3) * t315;
t135 = t205 * qJD(4) - t202 * t247;
t136 = -t202 * qJD(4) - t205 * t247;
t206 = cos(qJ(2));
t253 = t206 * t275;
t301 = t299 * t135 + t195 * t136 - t218 * t253;
t174 = qJD(2) * t242;
t273 = qJD(2) * t202;
t142 = t195 * t273 - t174;
t134 = qJD(6) + t142;
t326 = t134 - qJD(6);
t145 = t156 * qJD(2);
t194 = sin(pkin(12));
t198 = cos(pkin(12));
t120 = qJD(3) * t194 + t145 * t198;
t121 = t198 * qJD(3) - t145 * t194;
t201 = sin(qJ(6));
t204 = cos(qJ(6));
t59 = t120 * t201 - t121 * t204;
t325 = t134 * t59;
t157 = t194 * t204 + t198 * t201;
t149 = t157 * qJD(6);
t304 = t157 * t142 + t149;
t227 = -t120 * t204 - t121 * t201;
t324 = t134 * t227;
t199 = cos(pkin(6));
t285 = t197 * t203;
t150 = t199 * t205 - t202 * t285;
t312 = -t301 * t194 + t198 * t327;
t311 = t194 * t327 + t301 * t198;
t241 = qJD(2) * t315 + t256;
t274 = qJD(1) * t199;
t112 = t202 * t274 + t205 * t241;
t102 = t195 * t112;
t111 = -t202 * t241 + t205 * t274;
t55 = t299 * t111 - t102;
t261 = pkin(3) * t273;
t82 = pkin(4) * t145 + qJ(5) * t142 + t261;
t28 = t194 * t82 + t198 * t55;
t323 = qJD(5) * t198 - t28;
t27 = -t194 * t55 + t198 * t82;
t322 = -qJD(5) * t194 - t27;
t302 = t135 * t195 - t299 * t136 - t156 * t253;
t300 = cos(pkin(10));
t244 = t300 * t206;
t196 = sin(pkin(10));
t287 = t196 * t203;
t138 = -t199 * t244 + t287;
t245 = t300 * t203;
t286 = t196 * t206;
t140 = t199 * t286 + t245;
t238 = g(1) * t140 + g(2) * t138;
t283 = t197 * t206;
t215 = -g(3) * t283 + t238;
t155 = t194 * t201 - t204 * t198;
t305 = t134 * t155;
t264 = t202 * qJDD(2);
t95 = qJD(2) * t144 - qJDD(2) * t242 + t195 * t264;
t91 = qJDD(6) + t95;
t321 = t305 * t134 - t157 * t91;
t207 = qJD(3) ^ 2;
t249 = qJDD(1) * t283;
t267 = qJD(1) * qJD(2);
t252 = t203 * t267;
t232 = t197 * t252 - t249;
t320 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t207 + t197 * (-g(3) * t206 + t252) - t232 + t238;
t266 = qJD(2) * qJD(3);
t251 = t202 * t266;
t96 = qJD(3) * t174 + qJDD(2) * t156 - t195 * t251;
t78 = -t198 * qJDD(3) + t194 * t96;
t79 = qJDD(3) * t194 + t198 * t96;
t17 = -qJD(6) * t227 + t201 * t79 + t204 * t78;
t137 = t142 ^ 2;
t319 = pkin(3) * t195;
t318 = pkin(9) * t198;
t317 = g(3) * t197;
t179 = qJ(5) + t319;
t316 = pkin(9) + t179;
t265 = t199 * qJDD(1);
t173 = t205 * t265;
t126 = qJDD(2) * pkin(8) + (qJDD(1) * t203 + t206 * t267) * t197;
t210 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t274 + t126;
t225 = t241 * qJD(3);
t42 = qJDD(3) * pkin(3) - t202 * t210 - t205 * t225 + t173;
t43 = (-t225 + t265) * t202 + t210 * t205;
t14 = t195 * t42 + t299 * t43;
t11 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t14;
t183 = pkin(3) * t205 + pkin(2);
t94 = pkin(3) * t251 - t183 * qJDD(2) + qJDD(4) + t232;
t26 = t95 * pkin(4) - t96 * qJ(5) - t145 * qJD(5) + t94;
t7 = t198 * t11 + t194 * t26;
t314 = pkin(5) * t144 - t147 * t318 + t312;
t296 = t147 * t194;
t313 = pkin(9) * t296 - t311;
t106 = t111 + t309;
t243 = t299 * t112;
t51 = t195 * t106 + t243;
t46 = qJD(3) * qJ(5) + t51;
t133 = -t183 * qJD(2) + qJD(4) - t253;
t65 = t142 * pkin(4) - t145 * qJ(5) + t133;
t22 = t194 * t65 + t198 * t46;
t310 = qJD(2) * pkin(2);
t308 = t145 * t59;
t307 = t145 * t227;
t165 = t315 * t202;
t166 = t315 * t205;
t117 = -t195 * t165 + t299 * t166;
t92 = -pkin(4) * t218 - qJ(5) * t156 - t183;
t48 = t198 * t117 + t194 * t92;
t303 = pkin(5) * t296 + t302;
t297 = t142 * t194;
t295 = t156 * t194;
t294 = t156 * t198;
t293 = t179 * t194;
t292 = t179 * t198;
t190 = pkin(12) + qJ(6);
t186 = sin(t190);
t191 = qJ(3) + pkin(11);
t189 = cos(t191);
t291 = t186 * t189;
t188 = cos(t190);
t290 = t188 * t189;
t289 = t189 * t206;
t288 = t196 * t197;
t284 = t197 * t205;
t281 = t315 * t203;
t50 = t299 * t106 - t102;
t45 = -qJD(3) * pkin(4) + qJD(5) - t50;
t280 = -qJD(5) + t45;
t279 = qJDD(1) - g(3);
t139 = t199 * t245 + t286;
t278 = -t138 * t183 + t139 * t315;
t141 = -t199 * t287 + t244;
t277 = -t140 * t183 + t141 * t315;
t192 = t202 ^ 2;
t276 = -t205 ^ 2 + t192;
t272 = qJD(2) * t203;
t269 = qJD(6) * t201;
t268 = qJD(6) * t204;
t263 = t205 * qJDD(2);
t262 = g(3) * t285;
t6 = -t11 * t194 + t198 * t26;
t2 = pkin(5) * t95 - pkin(9) * t79 + t6;
t3 = -pkin(9) * t78 + t7;
t258 = t204 * t2 - t201 * t3;
t257 = t299 * pkin(3);
t255 = t197 * t272;
t254 = qJD(2) * t283;
t250 = t206 * t266;
t21 = -t194 * t46 + t198 * t65;
t246 = t197 * t300;
t47 = -t117 * t194 + t198 * t92;
t53 = t111 * t195 + t243;
t116 = t299 * t165 + t166 * t195;
t240 = -t134 * t304 - t155 * t91;
t239 = (-t141 * t202 + t196 * t284) * pkin(3);
t182 = -t257 - pkin(4);
t237 = g(1) * t141 + g(2) * t139;
t236 = -t7 * t194 - t6 * t198;
t235 = t201 * t2 + t204 * t3;
t13 = -t195 * t43 + t299 * t42;
t15 = pkin(9) * t121 + t22;
t9 = pkin(5) * t142 - pkin(9) * t120 + t21;
t5 = t15 * t204 + t201 * t9;
t234 = t15 * t201 - t204 * t9;
t187 = sin(t191);
t233 = pkin(4) * t189 + qJ(5) * t187;
t31 = -pkin(5) * t218 - pkin(9) * t294 + t47;
t33 = -pkin(9) * t295 + t48;
t231 = -t201 * t33 + t204 * t31;
t230 = t201 * t31 + t204 * t33;
t151 = t199 * t202 + t203 * t284;
t88 = t195 * t150 + t299 * t151;
t66 = -t194 * t88 - t198 * t283;
t67 = -t194 * t283 + t198 * t88;
t229 = -t201 * t67 + t204 * t66;
t228 = t201 * t66 + t204 * t67;
t226 = t150 * pkin(3);
t208 = qJD(2) ^ 2;
t224 = qJDD(2) * t206 - t203 * t208;
t222 = -g(1) * t196 + t300 * g(2);
t152 = t316 * t194;
t220 = pkin(9) * t297 + qJD(6) * t152 - t323;
t153 = t316 * t198;
t219 = pkin(5) * t145 + qJD(6) * t153 + t142 * t318 - t322;
t16 = -t120 * t269 + t121 * t268 - t201 * t78 + t204 * t79;
t100 = t141 * t187 - t189 * t288;
t128 = t187 * t285 - t199 * t189;
t98 = t139 * t187 + t189 * t246;
t217 = g(1) * t100 + g(2) * t98 + g(3) * t128;
t12 = -qJDD(3) * pkin(4) + qJDD(5) - t13;
t216 = -t12 + t217;
t162 = -t253 - t310;
t214 = -qJD(2) * t162 - t126 + t237;
t213 = (-t139 * t202 - t205 * t246) * pkin(3);
t212 = t12 * t156 + t147 * t45 - t237;
t209 = -pkin(8) * qJDD(3) + (t162 + t253 - t310) * qJD(3);
t163 = -t198 * pkin(5) + t182;
t159 = t183 * t283;
t129 = t187 * t199 + t189 * t285;
t110 = -qJD(3) * t151 - t202 * t254;
t109 = qJD(3) * t150 + t205 * t254;
t101 = t141 * t189 + t187 * t288;
t99 = t139 * t189 - t187 * t246;
t87 = -t299 * t150 + t195 * t151;
t86 = t155 * t156;
t85 = t157 * t156;
t77 = pkin(5) * t295 + t116;
t54 = t299 * t109 + t195 * t110;
t52 = t195 * t109 - t299 * t110;
t39 = t194 * t255 + t198 * t54;
t38 = -t194 * t54 + t198 * t255;
t36 = -pkin(5) * t297 + t53;
t35 = t147 * t157 + t268 * t294 - t269 * t295;
t34 = -t147 * t155 - t156 * t149;
t32 = -pkin(5) * t121 + t45;
t8 = t78 * pkin(5) + t12;
t1 = [t279, 0, t224 * t197 (-qJDD(2) * t203 - t206 * t208) * t197, 0, 0, 0, 0, 0, t110 * qJD(3) + t150 * qJDD(3) + (-t202 * t250 + t205 * t224) * t197, -t109 * qJD(3) - t151 * qJDD(3) + (-t202 * t224 - t205 * t250) * t197, -t142 * t54 + t145 * t52 + t87 * t96 - t88 * t95, -t13 * t87 + t14 * t88 - t50 * t52 + t51 * t54 - g(3) + (t133 * t272 - t206 * t94) * t197, -t121 * t52 + t142 * t38 + t66 * t95 + t78 * t87, t120 * t52 - t142 * t39 - t67 * t95 + t79 * t87, -t120 * t38 + t121 * t39 - t66 * t79 - t67 * t78, t12 * t87 + t21 * t38 + t22 * t39 + t45 * t52 + t6 * t66 + t67 * t7 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t228 - t201 * t39 + t204 * t38) * t134 + t229 * t91 + t52 * t59 + t87 * t17 -(qJD(6) * t229 + t201 * t38 + t204 * t39) * t134 - t228 * t91 - t52 * t227 + t87 * t16; 0, qJDD(2), t215 + t249, -t279 * t285 + t237, qJDD(2) * t192 + 0.2e1 * t205 * t251, 0.2e1 * t202 * t263 - 0.2e1 * t276 * t266, qJDD(3) * t202 + t205 * t207, qJDD(3) * t205 - t202 * t207, 0, t209 * t202 + t205 * t320, -t202 * t320 + t209 * t205, t116 * t96 - t117 * t95 - t13 * t156 + t14 * t218 - t301 * t142 - t51 * t144 + t302 * t145 - t50 * t147 - t237 - t262, t14 * t117 - t13 * t116 - t94 * t183 - g(1) * t277 - g(2) * t278 - g(3) * (t197 * t281 + t159) + t301 * t51 - t302 * t50 + t328 * t133, t116 * t78 + t21 * t144 - t6 * t218 + t47 * t95 + t215 * t198 * t189 + (t212 - t262) * t194 + t312 * t142 - t302 * t121, t116 * t79 - t22 * t144 + t7 * t218 - t48 * t95 - t238 * t194 * t189 + t212 * t198 - (-t194 * t289 + t198 * t203) * t317 - t311 * t142 + t302 * t120, -t47 * t79 - t48 * t78 + t236 * t156 + (-t194 * t22 - t198 * t21) * t147 - t312 * t120 + t311 * t121 + t215 * t187, t7 * t48 + t6 * t47 + t12 * t116 - g(1) * (-t140 * t233 + t277) - g(2) * (-t138 * t233 + t278) - g(3) * t159 + t302 * t45 + t311 * t22 + t312 * t21 - (t206 * t233 + t281) * t317, -t16 * t86 - t227 * t34, -t16 * t85 + t17 * t86 + t227 * t35 - t34 * t59, t134 * t34 - t144 * t227 - t16 * t218 - t86 * t91, -t134 * t35 - t144 * t59 + t17 * t218 - t85 * t91, t134 * t144 - t218 * t91, t231 * t91 - t258 * t218 - t234 * t144 + t77 * t17 + t8 * t85 + t32 * t35 - g(1) * (-t140 * t290 + t141 * t186) - g(2) * (-t138 * t290 + t139 * t186) + t303 * t59 - (t186 * t203 + t188 * t289) * t317 + (t313 * t201 + t314 * t204) * t134 + (-t134 * t230 + t218 * t5) * qJD(6), -t230 * t91 + t235 * t218 - t5 * t144 + t77 * t16 - t8 * t86 + t32 * t34 - g(1) * (t140 * t291 + t141 * t188) - g(2) * (t138 * t291 + t139 * t188) - t303 * t227 - (-t186 * t289 + t188 * t203) * t317 + (-t201 * t314 + t204 * t313) * t134 + (-t134 * t231 - t218 * t234) * qJD(6); 0, 0, 0, 0, -t202 * t208 * t205, t276 * t208, t264, t263, qJDD(3), -g(3) * t150 + t202 * t214 + t222 * t284 + t173, g(3) * t151 + (-t197 * t222 - t265) * t202 + t214 * t205 (t51 - t53) * t145 + (-t50 + t55) * t142 + (-t195 * t95 - t299 * t96) * pkin(3), -g(1) * t239 - g(2) * t213 - g(3) * t226 + t13 * t257 - t133 * t261 + t14 * t319 + t50 * t53 - t51 * t55, -t95 * t293 + t53 * t121 - t21 * t145 + t182 * t78 + (t280 * t194 - t27) * t142 + t216 * t198, -t95 * t292 - t53 * t120 + t22 * t145 + t182 * t79 + (t280 * t198 + t28) * t142 - t216 * t194, -g(1) * t101 - g(2) * t99 - g(3) * t129 - t28 * t121 + t27 * t120 + (qJD(5) * t121 - t142 * t21 - t179 * t78 + t7) * t198 + (qJD(5) * t120 - t142 * t22 + t179 * t79 - t6) * t194, t7 * t292 - t6 * t293 + t12 * t182 - t45 * t53 - g(1) * (-pkin(4) * t100 + qJ(5) * t101 + t239) - g(2) * (-t98 * pkin(4) + t99 * qJ(5) + t213) - g(3) * (-pkin(4) * t128 + qJ(5) * t129 + t226) + t323 * t22 + t322 * t21, t157 * t16 + t227 * t305, -t155 * t16 - t157 * t17 + t227 * t304 + t305 * t59, t307 - t321, t240 + t308, -t134 * t145 (-t152 * t204 - t153 * t201) * t91 + t163 * t17 + t8 * t155 + t234 * t145 - t36 * t59 + t304 * t32 + (t201 * t220 - t204 * t219) * t134 + t217 * t188 -(-t152 * t201 + t153 * t204) * t91 + t163 * t16 + t8 * t157 + t5 * t145 + t36 * t227 - t305 * t32 + (t201 * t219 + t204 * t220) * t134 - t217 * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145 ^ 2 - t137, t51 * t142 + t50 * t145 - t215 + t94, t121 * t145 - t137 * t194 + t198 * t95, -t120 * t145 - t137 * t198 - t194 * t95, -t194 * t78 - t198 * t79 + (t120 * t194 + t121 * t198) * t142, -t45 * t145 + (-t194 * t21 + t198 * t22) * t142 - t215 - t236, 0, 0, 0, 0, 0, t240 - t308, t307 + t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120 * t142 + t78, t121 * t142 + t79, -t120 ^ 2 - t121 ^ 2, t120 * t21 - t22 * t121 - t216, 0, 0, 0, 0, 0, t17 - t324, t16 - t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227 * t59, t227 ^ 2 - t59 ^ 2, t16 + t325, -t17 - t324, t91, t32 * t227 - g(1) * (-t101 * t186 + t140 * t188) - g(2) * (t138 * t188 - t186 * t99) - g(3) * (-t129 * t186 - t188 * t283) + t258 + t326 * t5, t32 * t59 - g(1) * (-t101 * t188 - t140 * t186) - g(2) * (-t138 * t186 - t188 * t99) - g(3) * (-t129 * t188 + t186 * t283) - t235 - t326 * t234;];
tau_reg  = t1;
