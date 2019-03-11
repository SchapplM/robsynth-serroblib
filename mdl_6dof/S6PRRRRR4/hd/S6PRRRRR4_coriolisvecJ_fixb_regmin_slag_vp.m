% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:01
% EndTime: 2019-03-09 01:00:22
% DurationCPUTime: 8.68s
% Computational Cost: add. (8876->485), mult. (23444->726), div. (0->0), fcn. (19747->14), ass. (0->254)
t212 = sin(pkin(7));
t219 = sin(qJ(3));
t333 = t212 * t219;
t203 = pkin(9) * t333;
t214 = cos(pkin(7));
t224 = cos(qJ(3));
t225 = cos(qJ(2));
t325 = t224 * t225;
t220 = sin(qJ(2));
t329 = t219 * t220;
t242 = -t214 * t329 + t325;
t213 = sin(pkin(6));
t317 = qJD(1) * t213;
t330 = t214 * t224;
t392 = t242 * t317 - (pkin(2) * t330 - t203) * qJD(3);
t258 = pkin(3) * t219 - pkin(10) * t224;
t239 = t258 * qJD(3);
t295 = t220 * t317;
t391 = (t239 - t295) * t212;
t218 = sin(qJ(4));
t223 = cos(qJ(4));
t332 = t212 * t224;
t158 = pkin(9) * t332 + (pkin(2) * t219 + pkin(10)) * t214;
t259 = -pkin(3) * t224 - pkin(10) * t219;
t159 = (-pkin(2) + t259) * t212;
t365 = t223 * t158 + t218 * t159;
t390 = qJD(4) * t365 - t392 * t218 - t391 * t223;
t307 = qJD(4) * t223;
t308 = qJD(4) * t218;
t389 = -t158 * t308 + t159 * t307 + t391 * t218 - t392 * t223;
t171 = -t223 * t214 + t218 * t333;
t309 = qJD(3) * t224;
t290 = t212 * t309;
t134 = -qJD(4) * t171 + t223 * t290;
t311 = qJD(3) * t219;
t291 = t212 * t311;
t388 = -pkin(4) * t291 + pkin(11) * t134 + t390;
t172 = t214 * t218 + t223 * t333;
t288 = t218 * t309;
t135 = qJD(4) * t172 + t212 * t288;
t387 = -pkin(11) * t135 + t389;
t312 = qJD(2) * t224;
t292 = t212 * t312;
t270 = t218 * t292;
t356 = pkin(10) + pkin(11);
t297 = qJD(4) * t356;
t315 = qJD(2) * t212;
t176 = pkin(9) * t315 + t295;
t168 = t219 * t176;
t215 = cos(pkin(6));
t316 = qJD(1) * t215;
t296 = t212 * t316;
t348 = qJD(2) * pkin(2);
t187 = t225 * t317 + t348;
t337 = t187 * t214;
t108 = t224 * (t296 + t337) - t168;
t163 = t258 * t315;
t324 = t223 * t108 + t218 * t163;
t386 = pkin(11) * t270 - t218 * t297 - t324;
t279 = -t108 * t218 + t223 * t163;
t326 = t223 * t224;
t385 = t223 * t297 + (pkin(4) * t219 - pkin(11) * t326) * t315 + t279;
t216 = sin(qJ(6));
t313 = qJD(2) * t214;
t202 = qJD(3) + t313;
t293 = t219 * t315;
t268 = t218 * t293;
t150 = t202 * t223 - t268;
t151 = t202 * t218 + t223 * t293;
t217 = sin(qJ(5));
t222 = cos(qJ(5));
t249 = t150 * t217 + t222 * t151;
t102 = -t222 * t150 + t151 * t217;
t302 = -qJD(6) - t102;
t221 = cos(qJ(6));
t301 = qJD(2) * qJD(3);
t284 = t212 * t301;
t266 = t224 * t284;
t123 = -qJD(4) * t268 + t202 * t307 + t223 * t266;
t124 = (t219 * t307 + t288) * t315 + t202 * t308;
t49 = qJD(5) * t249 + t123 * t217 + t222 * t124;
t44 = t221 * t49;
t195 = -qJD(4) + t292;
t188 = -qJD(5) + t195;
t83 = t221 * t188 + t216 * t249;
t384 = -t302 ^ 2 * t216 + t249 * t83 + t44;
t198 = t219 * t284;
t303 = qJD(6) * t221;
t304 = qJD(6) * t216;
t305 = qJD(5) * t222;
t306 = qJD(5) * t217;
t48 = t222 * t123 - t217 * t124 + t150 * t305 - t151 * t306;
t31 = -t188 * t303 + t216 * t198 + t221 * t48 - t249 * t304;
t29 = t31 * t216;
t380 = t102 * t221;
t85 = -t188 * t216 + t221 * t249;
t383 = t29 + (t303 + t380) * t85;
t352 = t216 * t49 - t302 * t303;
t382 = -t249 * t85 - t302 * t380 + t352;
t201 = t214 * t316;
t114 = t201 + (qJD(2) * t259 - t187) * t212;
t190 = t219 * t296;
t109 = t224 * t176 + t219 * t337 + t190;
t95 = pkin(10) * t202 + t109;
t66 = t114 * t218 + t223 * t95;
t54 = pkin(11) * t150 + t66;
t347 = t217 * t54;
t65 = t223 * t114 - t218 * t95;
t53 = -pkin(11) * t151 + t65;
t41 = -pkin(4) * t195 + t53;
t20 = t222 * t41 - t347;
t16 = pkin(5) * t188 - t20;
t381 = t102 * t16;
t379 = t249 * t102;
t177 = t217 * t218 - t222 * t223;
t364 = qJD(4) + qJD(5);
t132 = t364 * t177;
t137 = t177 * t292;
t322 = -t132 + t137;
t178 = t217 * t223 + t218 * t222;
t321 = (-t292 + t364) * t178;
t327 = t220 * t224;
t328 = t219 * t225;
t244 = t214 * t327 + t328;
t289 = t214 * t311;
t319 = pkin(2) * t289 + pkin(9) * t290 - t244 * t317;
t377 = -t102 ^ 2 + t249 ^ 2;
t64 = pkin(5) * t249 + pkin(12) * t102;
t376 = -t102 * t188 + t48;
t128 = (t239 + t295) * t315;
t231 = t242 * qJD(2);
t269 = t215 * t290;
t75 = (t187 * t330 - t168) * qJD(3) + (t213 * t231 + t269) * qJD(1);
t229 = -qJD(4) * t66 + t223 * t128 - t218 * t75;
t22 = pkin(4) * t198 - pkin(11) * t123 + t229;
t237 = t114 * t307 + t218 * t128 + t223 * t75 - t95 * t308;
t27 = -pkin(11) * t124 + t237;
t282 = -t217 * t22 - t222 * t27 - t41 * t305 + t54 * t306;
t94 = -pkin(3) * t202 - t108;
t74 = -pkin(4) * t150 + t94;
t375 = t102 * t74 + t282;
t254 = t216 * t85 + t221 * t83;
t32 = qJD(6) * t85 - t221 * t198 + t216 * t48;
t374 = -t102 * t254 - t216 * t32 + t31 * t221 - t83 * t303 - t85 * t304;
t277 = -t158 * t218 + t223 * t159;
t81 = -pkin(4) * t332 - pkin(11) * t172 + t277;
t87 = -pkin(11) * t171 + t365;
t371 = -t217 * t388 + t387 * t222 + t81 * t305 - t87 * t306;
t196 = t356 * t218;
t197 = t356 * t223;
t140 = -t196 * t217 + t197 * t222;
t370 = qJD(5) * t140 + t217 * t386 + t222 * t385;
t350 = t217 * t81 + t222 * t87;
t369 = t350 * qJD(5) + t387 * t217 + t222 * t388;
t247 = -t196 * t222 - t197 * t217;
t368 = qJD(5) * t247 - t217 * t385 + t222 * t386;
t265 = -t109 + (-t270 + t308) * pkin(4);
t367 = t302 * t249;
t366 = t219 * t224;
t323 = pkin(4) * t135 + t319;
t346 = t222 * t54;
t21 = t217 * t41 + t346;
t17 = -pkin(12) * t188 + t21;
t35 = pkin(5) * t102 - pkin(12) * t249 + t74;
t255 = t17 * t216 - t221 * t35;
t230 = -qJD(5) * t21 - t217 * t27 + t222 * t22;
t5 = -pkin(5) * t198 - t230;
t363 = t16 * t304 - t5 * t221 + t249 * t255;
t7 = t17 * t221 + t216 * t35;
t362 = t16 * t303 + t5 * t216 + t7 * t249;
t360 = -t249 * t74 + t230;
t358 = -t188 * t249 - t49;
t314 = qJD(2) * t213;
t285 = qJD(1) * t314;
t246 = t214 * t220 * t285;
t76 = qJD(3) * t190 + t176 * t309 + t187 * t289 + t224 * t246 + t285 * t328;
t60 = pkin(4) * t124 + t76;
t13 = pkin(5) * t49 - pkin(12) * t48 + t60;
t4 = pkin(12) * t198 - t282;
t1 = -qJD(6) * t255 + t13 * t216 + t221 * t4;
t354 = -pkin(5) * t291 + t369;
t349 = pkin(5) * t293 + t370;
t341 = t150 * t195;
t340 = t151 * t195;
t339 = t178 * t216;
t338 = t178 * t221;
t336 = t195 * t218;
t335 = t195 * t223;
t209 = t212 ^ 2;
t226 = qJD(2) ^ 2;
t334 = t209 * t226;
t331 = t213 * t226;
t318 = t219 ^ 2 - t224 ^ 2;
t310 = qJD(3) * t223;
t299 = t216 * t332;
t298 = t220 * t331;
t208 = -pkin(4) * t223 - pkin(3);
t294 = t209 * t312;
t206 = pkin(4) * t217 + pkin(12);
t280 = pkin(4) * t151 + qJD(6) * t206 + t64;
t275 = t202 + t313;
t274 = 0.2e1 * t209 * t301;
t273 = t209 * t298;
t271 = t212 * t220 * t314;
t120 = t222 * t171 + t172 * t217;
t121 = -t171 * t217 + t172 * t222;
t157 = t203 + (-pkin(2) * t224 - pkin(3)) * t214;
t122 = pkin(4) * t171 + t157;
t57 = pkin(5) * t120 - pkin(12) * t121 + t122;
t267 = -pkin(12) * t291 - qJD(6) * t57 - t371;
t25 = t217 * t53 + t346;
t264 = pkin(4) * t306 - t25;
t26 = t222 * t53 - t347;
t263 = -pkin(4) * t305 + t26;
t127 = pkin(5) * t177 - pkin(12) * t178 + t208;
t262 = pkin(12) * t293 - qJD(6) * t127 - t368;
t261 = -pkin(5) * t321 + pkin(12) * t322 + qJD(6) * t140 - t265;
t39 = -pkin(12) * t332 + t350;
t62 = -qJD(5) * t120 + t134 * t222 - t135 * t217;
t63 = qJD(5) * t121 + t134 * t217 + t222 * t135;
t260 = -pkin(5) * t63 + pkin(12) * t62 + qJD(6) * t39 - t323;
t256 = -t206 * t49 + t381;
t252 = -t217 * t87 + t222 * t81;
t243 = t214 * t328 + t327;
t131 = t213 * t243 + t215 * t333;
t170 = -t212 * t213 * t225 + t214 * t215;
t96 = -t131 * t218 + t170 * t223;
t97 = t131 * t223 + t170 * t218;
t58 = t217 * t97 - t222 * t96;
t59 = t217 * t96 + t222 * t97;
t245 = t214 * t325 - t329;
t130 = -t213 * t245 - t215 * t332;
t251 = t130 * t221 - t216 * t59;
t250 = t130 * t216 + t221 * t59;
t98 = t121 * t216 + t221 * t332;
t115 = -t137 * t216 - t221 * t293;
t235 = -t216 * t132 + t178 * t303 - t115;
t116 = -t137 * t221 + t216 * t293;
t234 = -t221 * t132 - t178 * t304 - t116;
t146 = -t187 * t212 + t201;
t232 = qJD(3) * (t146 * t212 - t209 * t348);
t2 = -qJD(6) * t7 + t221 * t13 - t216 * t4;
t207 = -pkin(4) * t222 - pkin(5);
t99 = t121 * t221 - t299;
t91 = t269 + (qJD(3) * t245 + t231) * t213;
t90 = t215 * t291 + (qJD(2) * t244 + qJD(3) * t243) * t213;
t47 = qJD(4) * t96 + t218 * t271 + t223 * t91;
t46 = -qJD(4) * t97 - t218 * t91 + t223 * t271;
t38 = pkin(5) * t332 - t252;
t37 = -qJD(6) * t299 + t121 * t303 + t216 * t62 - t221 * t291;
t36 = -qJD(6) * t98 + t216 * t291 + t221 * t62;
t11 = qJD(5) * t59 + t217 * t47 - t222 * t46;
t10 = -qJD(5) * t58 + t217 * t46 + t222 * t47;
t3 = [0, 0, -t298, -t225 * t331, 0, 0, 0, 0, 0, t170 * t198 - t202 * t90 - t224 * t273, t170 * t266 - t202 * t91 + t219 * t273, 0, 0, 0, 0, 0, t124 * t130 - t150 * t90 - t195 * t46 + t198 * t96, t123 * t130 + t151 * t90 + t195 * t47 - t198 * t97, 0, 0, 0, 0, 0, t102 * t90 + t11 * t188 + t130 * t49 - t198 * t58, t10 * t188 + t130 * t48 - t198 * t59 + t249 * t90, 0, 0, 0, 0, 0 -(-qJD(6) * t250 - t10 * t216 + t221 * t90) * t302 + t251 * t49 + t11 * t83 + t58 * t32 (qJD(6) * t251 + t10 * t221 + t216 * t90) * t302 - t250 * t49 + t11 * t85 + t58 * t31; 0, 0, 0, 0, t274 * t366, -t318 * t274, t275 * t290, -t275 * t291, 0, -t319 * t202 - t214 * t76 + t219 * t232, t392 * t202 - t214 * t75 + t224 * t232, t123 * t172 + t134 * t151, -t123 * t171 - t124 * t172 + t134 * t150 - t135 * t151, -t134 * t195 + (-t123 * t224 + (qJD(2) * t172 + t151) * t311) * t212, t135 * t195 + (t124 * t224 + (-qJD(2) * t171 + t150) * t311) * t212 (-t195 * t212 - t294) * t311, t157 * t124 + t94 * t135 + t76 * t171 + t390 * t195 - t319 * t150 + (-t229 * t224 + (qJD(2) * t277 + t65) * t311) * t212, t157 * t123 + t94 * t134 + t76 * t172 + t389 * t195 + t319 * t151 + (t237 * t224 + (-qJD(2) * t365 - t66) * t311) * t212, t121 * t48 + t249 * t62, -t102 * t62 - t120 * t48 - t121 * t49 - t249 * t63, -t188 * t62 + (-t224 * t48 + (qJD(2) * t121 + t249) * t311) * t212, t188 * t63 + (t224 * t49 + (-qJD(2) * t120 - t102) * t311) * t212 (-t188 * t212 - t294) * t311, t60 * t120 + t122 * t49 + t74 * t63 + t369 * t188 + t323 * t102 + (-t230 * t224 + (qJD(2) * t252 + t20) * t311) * t212, t60 * t121 + t122 * t48 + t74 * t62 + t371 * t188 + t323 * t249 + (-t282 * t224 + (-qJD(2) * t350 - t21) * t311) * t212, t31 * t99 + t36 * t85, -t31 * t98 - t32 * t99 - t36 * t83 - t37 * t85, t120 * t31 - t302 * t36 + t49 * t99 + t63 * t85, -t120 * t32 + t302 * t37 - t49 * t98 - t63 * t83, t120 * t49 - t302 * t63 (-t216 * t39 + t221 * t57) * t49 + t2 * t120 - t255 * t63 + t38 * t32 + t5 * t98 + t16 * t37 + t354 * t83 - (t216 * t267 - t221 * t260) * t302 -(t216 * t57 + t221 * t39) * t49 - t1 * t120 - t7 * t63 + t38 * t31 + t5 * t99 + t16 * t36 + t354 * t85 - (t216 * t260 + t221 * t267) * t302; 0, 0, 0, 0, -t334 * t366, t318 * t334 (qJD(3) - t202) * t292, t202 * t293 - t198, 0, t109 * t202 - t146 * t293 - t76, t108 * t202 + (qJD(3) * t176 + t246) * t219 + (-t146 * t315 - qJD(3) * t337 + (-qJD(3) * t212 * t215 - t225 * t314) * qJD(1)) * t224, t123 * t218 - t151 * t335 (t123 - t341) * t223 + (-t124 + t340) * t218, -t195 * t307 + (t195 * t326 + (qJD(3) * t218 - t151) * t219) * t315, t195 * t308 + (-t224 * t336 + (-t150 + t310) * t219) * t315, t195 * t293, -pkin(3) * t124 - t76 * t223 + t279 * t195 + t109 * t150 + (pkin(10) * t335 + t218 * t94) * qJD(4) + (-t219 * t65 + (-pkin(10) * t311 - t224 * t94) * t218) * t315, -pkin(3) * t123 + t76 * t218 - t109 * t151 - t324 * t195 + (-pkin(10) * t336 + t223 * t94) * qJD(4) + (-t94 * t326 + (-pkin(10) * t310 + t66) * t219) * t315, t178 * t48 + t249 * t322, -t102 * t322 - t177 * t48 - t178 * t49 - t249 * t321, -t322 * t188 + (qJD(3) * t178 - t249) * t293, t321 * t188 + (-qJD(3) * t177 + t102) * t293, t188 * t293, t60 * t177 + t208 * t49 + t321 * t74 + t370 * t188 + t265 * t102 + (qJD(3) * t247 - t20) * t293, t60 * t178 + t208 * t48 + t322 * t74 + t368 * t188 + t265 * t249 + (-qJD(3) * t140 + t21) * t293, t234 * t85 + t31 * t338, t115 * t85 + t116 * t83 + t254 * t132 + (-t29 - t221 * t32 + (t216 * t83 - t221 * t85) * qJD(6)) * t178, t177 * t31 - t234 * t302 + t321 * t85 + t338 * t49, -t177 * t32 + t235 * t302 - t321 * t83 - t339 * t49, t177 * t49 - t302 * t321 (t127 * t221 - t140 * t216) * t49 + t2 * t177 - t247 * t32 + t5 * t339 + t349 * t83 - t321 * t255 - (t216 * t262 - t221 * t261) * t302 + t235 * t16 -(t127 * t216 + t140 * t221) * t49 - t1 * t177 - t247 * t31 + t5 * t338 + t349 * t85 - t321 * t7 - (t216 * t261 + t221 * t262) * t302 + t234 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151 * t150, -t150 ^ 2 + t151 ^ 2, t123 + t341, -t124 - t340, t198, -t151 * t94 - t195 * t66 + t229, -t150 * t94 - t195 * t65 - t237, t379, t377, t376, t358, t198, -t188 * t25 + (-t102 * t151 + t188 * t306 + t198 * t222) * pkin(4) + t360, -t188 * t26 + (-t151 * t249 + t188 * t305 - t198 * t217) * pkin(4) + t375, t383, t374, t382, t384, t367, t207 * t32 + t264 * t83 + t256 * t216 - (t216 * t263 - t221 * t280) * t302 + t363, t207 * t31 + t264 * t85 + t256 * t221 - (t216 * t280 + t221 * t263) * t302 + t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t379, t377, t376, t358, t198, -t188 * t21 + t360, -t188 * t20 + t375, t383, t374, t382, t384, t367, -pkin(5) * t32 + (-t20 * t216 + t221 * t64) * t302 - t21 * t83 + t216 * t381 - t352 * pkin(12) + t363, -pkin(5) * t31 - (t20 * t221 + t216 * t64) * t302 - t21 * t85 + t16 * t380 + (-t302 * t304 - t44) * pkin(12) + t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t83, -t83 ^ 2 + t85 ^ 2, -t302 * t83 + t31, -t302 * t85 - t32, t49, -t16 * t85 - t302 * t7 + t2, t16 * t83 + t255 * t302 - t1;];
tauc_reg  = t3;
