% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:55
% EndTime: 2019-03-09 15:54:11
% DurationCPUTime: 5.72s
% Computational Cost: add. (8686->473), mult. (23379->649), div. (0->0), fcn. (18331->10), ass. (0->233)
t217 = sin(qJ(3));
t213 = sin(pkin(6));
t221 = cos(qJ(2));
t291 = qJD(1) * t221;
t275 = t213 * t291;
t258 = t217 * t275;
t287 = qJD(3) * t217;
t342 = t258 - t287;
t218 = sin(qJ(2));
t220 = cos(qJ(3));
t293 = qJD(1) * t213;
t276 = t218 * t293;
t215 = cos(pkin(6));
t292 = qJD(1) * t215;
t281 = pkin(1) * t292;
t164 = -pkin(8) * t276 + t221 * t281;
t241 = (pkin(2) * t218 - pkin(9) * t221) * t213;
t165 = qJD(1) * t241;
t266 = -t164 * t217 + t220 * t165;
t323 = -qJ(4) - pkin(9);
t268 = qJD(3) * t323;
t301 = t220 * t221;
t341 = -t217 * qJD(4) + t220 * t268 - (pkin(3) * t218 - qJ(4) * t301) * t293 - t266;
t296 = t220 * t164 + t217 * t165;
t340 = qJ(4) * t258 + t220 * qJD(4) + t217 * t268 - t296;
t212 = sin(pkin(11));
t214 = cos(pkin(11));
t286 = qJD(3) * t220;
t302 = t214 * t220;
t295 = t342 * t212 + t214 * t286 - t275 * t302;
t216 = sin(qJ(6));
t219 = cos(qJ(6));
t200 = qJD(2) + t292;
t259 = t217 * t276;
t147 = -t220 * t200 + t259;
t149 = t200 * t217 + t220 * t276;
t244 = -t147 * t212 + t214 * t149;
t328 = qJD(6) + t244;
t265 = t328 ^ 2;
t282 = qJD(1) * qJD(2);
t269 = t213 * t282;
t257 = t221 * t269;
t122 = -qJD(3) * t259 + t200 * t286 + t220 * t257;
t288 = qJD(2) * t221;
t272 = t217 * t288;
t123 = (t218 * t286 + t272) * t293 + t200 * t287;
t77 = t122 * t214 - t123 * t212;
t231 = -t216 * t77 - t219 * t265;
t108 = t214 * t147 + t149 * t212;
t339 = pkin(5) * t108;
t190 = -qJD(3) + t275;
t82 = -t219 * t108 - t190 * t216;
t338 = t328 * t82;
t311 = t341 * t212 + t340 * t214;
t264 = t328 * t216;
t73 = t219 * t77;
t337 = t328 * t264 - t73;
t336 = t108 * t190;
t177 = t212 * t220 + t214 * t217;
t130 = t177 * t275;
t172 = t177 * qJD(3);
t298 = t130 - t172;
t196 = t218 * t281;
t167 = pkin(8) * t275 + t196;
t335 = -t342 * pkin(3) - t167;
t334 = t244 ^ 2;
t333 = pkin(5) * t244;
t135 = -pkin(2) * t200 - t164;
t106 = pkin(3) * t147 + qJD(4) + t135;
t223 = -qJ(5) * t244 + t106;
t46 = pkin(4) * t108 + t223;
t332 = t244 * t46;
t313 = t340 * t212 - t341 * t214;
t312 = qJ(5) * t276 - t311;
t136 = pkin(9) * t200 + t167;
t162 = (-pkin(2) * t221 - pkin(9) * t218 - pkin(1)) * t213;
t143 = qJD(1) * t162;
t98 = t136 * t220 + t143 * t217;
t75 = -qJ(4) * t147 + t98;
t68 = t212 * t75;
t97 = -t136 * t217 + t220 * t143;
t74 = -qJ(4) * t149 + t97;
t40 = t214 * t74 - t68;
t300 = -qJD(5) + t40;
t303 = t213 * t221;
t161 = pkin(8) * t303 + (pkin(1) * t218 + pkin(9)) * t215;
t297 = t220 * t161 + t217 * t162;
t330 = t295 * qJ(5) + qJD(5) * t177 - t335;
t329 = -qJD(6) + t328;
t304 = t213 * t218;
t174 = -t215 * t220 + t217 * t304;
t273 = t213 * t288;
t128 = -qJD(3) * t174 + t220 * t273;
t175 = t215 * t217 + t220 * t304;
t166 = qJD(2) * t241;
t201 = pkin(8) * t304;
t324 = pkin(1) * t221;
t168 = (t215 * t324 - t201) * qJD(2);
t226 = -t297 * qJD(3) + t220 * t166 - t217 * t168;
t290 = qJD(2) * t218;
t274 = t213 * t290;
t50 = pkin(3) * t274 - qJ(4) * t128 - qJD(4) * t175 + t226;
t127 = qJD(3) * t175 + t213 * t272;
t235 = -t161 * t287 + t162 * t286 + t217 * t166 + t220 * t168;
t55 = -qJ(4) * t127 - qJD(4) * t174 + t235;
t18 = t212 * t50 + t214 * t55;
t15 = -t213 * (qJ(5) * t290 - qJD(5) * t221) - t18;
t207 = -pkin(3) * t214 - pkin(4);
t203 = -pkin(10) + t207;
t318 = t214 * t75;
t66 = -pkin(3) * t190 + t74;
t34 = t212 * t66 + t318;
t29 = qJ(5) * t190 - t34;
t22 = -t29 - t339;
t39 = t212 * t74 + t318;
t327 = t203 * t77 + (t22 - t39 + t339) * t328;
t176 = t212 * t217 - t302;
t277 = -pkin(3) * t220 - pkin(2);
t243 = -qJ(5) * t177 + t277;
t325 = pkin(4) + pkin(10);
t112 = t325 * t176 + t243;
t191 = t323 * t217;
t192 = t323 * t220;
t133 = -t214 * t191 - t192 * t212;
t115 = pkin(5) * t177 + t133;
t260 = t325 * t304;
t193 = t218 * t269;
t156 = qJD(1) * t166;
t157 = qJD(1) * t168;
t227 = -qJD(3) * t98 + t220 * t156 - t217 * t157;
t35 = pkin(3) * t193 - qJ(4) * t122 - qJD(4) * t149 + t227;
t236 = -t136 * t287 + t143 * t286 + t217 * t156 + t220 * t157;
t41 = -qJ(4) * t123 - qJD(4) * t147 + t236;
t13 = t212 * t35 + t214 * t41;
t263 = qJ(5) * t193 - qJD(5) * t190 + t13;
t76 = t122 * t212 + t214 * t123;
t6 = -pkin(5) * t76 + t263;
t326 = t328 * (-t295 * pkin(5) - qJD(1) * t260 + qJD(6) * t112 - t313) - t115 * t77 + t22 * t172 + t6 * t176;
t322 = t298 * pkin(4) + t330;
t267 = -t161 * t217 + t220 * t162;
t80 = -pkin(3) * t303 - qJ(4) * t175 + t267;
t90 = -qJ(4) * t174 + t297;
t54 = t212 * t80 + t214 * t90;
t321 = t298 * pkin(5) - t312;
t320 = t108 * t82;
t319 = t244 * t39;
t284 = qJD(6) * t219;
t285 = qJD(6) * t216;
t42 = t108 * t284 + t190 * t285 + t219 * t193 + t216 * t76;
t316 = t42 * t219;
t84 = t108 * t216 - t190 * t219;
t315 = t84 * t108;
t314 = pkin(4) * t276 + t313;
t310 = t147 * t190;
t309 = t149 * t190;
t308 = t176 * t216;
t307 = t190 * t217;
t306 = t190 * t220;
t209 = t213 ^ 2;
t305 = t209 * qJD(1) ^ 2;
t299 = -t300 + t333;
t158 = pkin(8) * t257 + qJD(2) * t196;
t169 = t215 * pkin(1) * t290 + pkin(8) * t273;
t294 = t218 ^ 2 - t221 ^ 2;
t289 = qJD(2) * t220;
t280 = t218 * t305;
t279 = t219 * t303;
t270 = t209 * t282;
t12 = -t212 * t41 + t214 * t35;
t17 = -t212 * t55 + t214 * t50;
t33 = t214 * t66 - t68;
t53 = -t212 * t90 + t214 * t80;
t262 = t200 + t292;
t261 = 0.2e1 * t270;
t94 = pkin(3) * t123 + t158;
t256 = pkin(3) * t127 + t169;
t49 = pkin(4) * t303 - t53;
t255 = qJD(5) - t33;
t253 = -0.2e1 * pkin(1) * t270;
t229 = -qJ(5) * t77 - qJD(5) * t244 + t94;
t14 = t325 * t76 + t229;
t246 = qJD(2) * t260;
t7 = pkin(5) * t77 - qJD(1) * t246 - t12;
t252 = t219 * t14 + t216 * t7;
t251 = t149 * pkin(3) + qJ(5) * t108;
t134 = t191 * t212 - t192 * t214;
t249 = t133 * t77 - t134 * t76;
t21 = t325 * t190 + t255 + t333;
t26 = t325 * t108 + t223;
t3 = t21 * t219 - t216 * t26;
t4 = t21 * t216 + t219 * t26;
t120 = -t174 * t212 + t175 * t214;
t27 = pkin(5) * t120 + pkin(10) * t303 + t49;
t119 = t214 * t174 + t175 * t212;
t160 = t201 + (-pkin(2) - t324) * t215;
t230 = t174 * pkin(3) + t160;
t225 = -qJ(5) * t120 + t230;
t45 = t325 * t119 + t225;
t248 = -t216 * t45 + t219 * t27;
t247 = t216 * t27 + t219 * t45;
t245 = -t108 ^ 2 - t334;
t48 = qJ(5) * t303 - t54;
t103 = t119 * t219 + t216 * t303;
t239 = t6 + (-qJD(6) * t203 + t244 * t325 + t251) * t328;
t237 = t193 * t216 - t219 * t76;
t118 = t130 * t216 + t219 * t276;
t233 = t172 * t216 + t176 * t284 - t118;
t11 = -pkin(4) * t193 - t12;
t89 = -t127 * t212 + t128 * t214;
t228 = -qJ(5) * t89 - qJD(5) * t120 + t256;
t2 = -qJD(6) * t4 - t216 * t14 + t219 * t7;
t20 = pkin(4) * t76 + t229;
t224 = t22 * qJD(6) * t176 - t112 * t77 + (-qJD(6) * t115 + t298 * t325 + t330) * t328;
t204 = pkin(3) * t212 + qJ(5);
t124 = pkin(4) * t176 + t243;
t117 = -t219 * t130 + t216 * t276;
t116 = -pkin(5) * t176 + t134;
t104 = t119 * t216 - t279;
t88 = t214 * t127 + t128 * t212;
t63 = pkin(4) * t119 + t225;
t60 = pkin(4) * t244 + t251;
t57 = qJD(6) * t103 + t216 * t88 + t219 * t274;
t56 = -qJD(6) * t279 - t219 * t88 + (qJD(6) * t119 + t274) * t216;
t43 = qJD(6) * t84 + t237;
t30 = -pkin(5) * t119 - t48;
t28 = pkin(4) * t190 + t255;
t23 = pkin(4) * t88 + t228;
t19 = t325 * t88 + t228;
t16 = -pkin(4) * t274 - t17;
t9 = -pkin(5) * t88 - t15;
t8 = pkin(5) * t89 - t17 - t246;
t1 = qJD(6) * t3 + t252;
t5 = [0, 0, 0, t218 * t221 * t261, -t294 * t261, t262 * t273, -t262 * t274, 0, -t158 * t215 - t169 * t200 + t218 * t253, -t157 * t215 - t168 * t200 + t221 * t253, t122 * t175 + t128 * t149, -t122 * t174 - t123 * t175 - t127 * t149 - t128 * t147, -t128 * t190 + (-t122 * t221 + (qJD(1) * t175 + t149) * t290) * t213, t127 * t190 + (t123 * t221 + (-qJD(1) * t174 - t147) * t290) * t213 (-t190 * t213 - t209 * t291) * t290, -t226 * t190 + t169 * t147 + t160 * t123 + t158 * t174 + t135 * t127 + (-t227 * t221 + (qJD(1) * t267 + t97) * t290) * t213, t235 * t190 + t169 * t149 + t160 * t122 + t158 * t175 + t135 * t128 + (t236 * t221 + (-qJD(1) * t297 - t98) * t290) * t213, -t108 * t18 - t119 * t13 - t12 * t120 - t17 * t244 - t33 * t89 - t34 * t88 - t53 * t77 - t54 * t76, t106 * t256 + t12 * t53 + t13 * t54 + t33 * t17 + t34 * t18 + t230 * t94, t108 * t15 + t11 * t120 - t119 * t263 + t16 * t244 + t28 * t89 + t29 * t88 + t48 * t76 + t49 * t77, -t108 * t23 - t119 * t20 - t16 * t190 - t46 * t88 - t63 * t76 + (-t11 * t221 + (qJD(1) * t49 + t28) * t290) * t213, -t244 * t23 - t120 * t20 + t15 * t190 - t46 * t89 - t63 * t77 + (-t263 * t221 + (-qJD(1) * t48 - t29) * t290) * t213, t11 * t49 + t15 * t29 + t16 * t28 + t20 * t63 + t23 * t46 - t263 * t48, t104 * t42 + t57 * t84, t103 * t42 - t104 * t43 - t56 * t84 - t57 * t82, t104 * t77 + t120 * t42 + t328 * t57 + t84 * t89, t103 * t77 - t120 * t43 - t328 * t56 - t82 * t89, t120 * t77 + t328 * t89 (-qJD(6) * t247 - t216 * t19 + t219 * t8) * t328 + t248 * t77 + t2 * t120 + t3 * t89 + t9 * t82 + t30 * t43 - t6 * t103 + t22 * t56 -(qJD(6) * t248 + t219 * t19 + t216 * t8) * t328 - t247 * t77 - t1 * t120 - t4 * t89 + t9 * t84 + t30 * t42 + t6 * t104 + t22 * t57; 0, 0, 0, -t221 * t280, t294 * t305 (qJD(2) - t200) * t275, t200 * t276 - t193, 0, pkin(1) * t280 + t167 * t200 - t158, pkin(8) * t193 + t164 * t200 + (-t215 * t282 + t305) * t324, t122 * t217 - t149 * t306 (t122 + t310) * t220 + (-t123 + t309) * t217, -t190 * t286 + (t190 * t301 + (qJD(2) * t217 - t149) * t218) * t293, t190 * t287 + (-t221 * t307 + (t147 + t289) * t218) * t293, t190 * t276, -pkin(2) * t123 - t158 * t220 + t266 * t190 - t167 * t147 + (pkin(9) * t306 + t135 * t217) * qJD(3) + (-t97 * t218 + (-pkin(9) * t290 - t135 * t221) * t217) * t293, -pkin(2) * t122 + t158 * t217 - t296 * t190 - t167 * t149 + (-pkin(9) * t307 + t135 * t220) * qJD(3) + (-t135 * t301 + (-pkin(9) * t289 + t98) * t218) * t293, -t108 * t311 - t12 * t177 - t13 * t176 + t244 * t313 - t295 * t33 + t298 * t34 + t249, t335 * t106 - t12 * t133 + t13 * t134 + t94 * t277 + t311 * t34 - t313 * t33, t108 * t312 + t11 * t177 - t176 * t263 + t244 * t314 + t28 * t295 - t29 * t298 + t249, -t124 * t76 - t176 * t20 + t298 * t46 - t314 * t190 + t322 * t108 + (qJD(2) * t133 - t28) * t276, -t124 * t77 - t177 * t20 - t295 * t46 + t312 * t190 + t322 * t244 + (qJD(2) * t134 + t29) * t276, t11 * t133 + t124 * t20 + t134 * t263 + t314 * t28 + t312 * t29 - t322 * t46, t233 * t84 + t308 * t42, t117 * t84 + t118 * t82 + (-t216 * t82 + t219 * t84) * t172 + (-t216 * t43 + t316 + (-t216 * t84 - t219 * t82) * qJD(6)) * t176, t177 * t42 + t233 * t328 + t295 * t84 + t308 * t77, t176 * t73 - t177 * t43 - t295 * t82 + (t172 * t219 - t176 * t285 + t117) * t328, t77 * t177 + t295 * t328, t116 * t43 - t22 * t117 + t2 * t177 + t224 * t216 - t326 * t219 + t295 * t3 + t321 * t82, -t1 * t177 + t116 * t42 - t22 * t118 + t326 * t216 + t224 * t219 - t295 * t4 + t321 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149 * t147, -t147 ^ 2 + t149 ^ 2, t122 - t310, -t123 - t309, t193, -t135 * t149 - t190 * t98 + t227, t135 * t147 - t190 * t97 - t236, t244 * t34 - t319 + (-t212 * t76 - t214 * t77) * pkin(3) + (t40 - t33) * t108, t33 * t39 - t34 * t40 + (-t106 * t149 + t12 * t214 + t13 * t212) * pkin(3), -t204 * t76 + t207 * t77 - t244 * t29 - t319 + (t28 + t300) * t108, t332 + t108 * t60 + t190 * t39 + (-pkin(4) + t207) * t193 - t12, -t108 * t46 + t190 * t300 + t193 * t204 + t244 * t60 + t263, t11 * t207 + t204 * t263 - t28 * t39 + t29 * t300 - t46 * t60, -t264 * t84 + t316 (-t328 * t84 - t43) * t219 + (-t42 + t338) * t216, t315 - t337, t231 - t320, t328 * t108, t108 * t3 + t204 * t43 + t239 * t216 + t327 * t219 + t299 * t82, -t108 * t4 + t204 * t42 - t327 * t216 + t239 * t219 + t299 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, t108 * t34 + t244 * t33 + t94, t245, t190 * t244 - t76, -t77 - t336, -t108 * t29 - t244 * t28 + t20, 0, 0, 0, 0, 0, t231 + t320, t315 + t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 - t336, -t108 * t244 + t193, -t190 ^ 2 - t334, -t190 * t29 + t11 + t332, 0, 0, 0, 0, 0, t190 * t82 - t216 * t265 + t73, t190 * t84 + t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t82, -t82 ^ 2 + t84 ^ 2, t42 + t338, t329 * t84 - t237, t77, -t22 * t84 + t328 * t4 + t2, t22 * t82 + t329 * t3 - t252;];
tauc_reg  = t5;
