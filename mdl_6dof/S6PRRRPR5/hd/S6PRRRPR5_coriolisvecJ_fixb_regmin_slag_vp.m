% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:19
% EndTime: 2019-03-08 23:28:37
% DurationCPUTime: 6.63s
% Computational Cost: add. (7430->431), mult. (20375->642), div. (0->0), fcn. (16940->14), ass. (0->237)
t204 = sin(pkin(7));
t211 = sin(qJ(3));
t316 = t204 * t211;
t195 = pkin(9) * t316;
t207 = cos(pkin(7));
t215 = cos(qJ(3));
t216 = cos(qJ(2));
t308 = t215 * t216;
t212 = sin(qJ(2));
t312 = t211 * t212;
t232 = -t207 * t312 + t308;
t205 = sin(pkin(6));
t301 = qJD(1) * t205;
t313 = t207 * t215;
t361 = t232 * t301 - (pkin(2) * t313 - t195) * qJD(3);
t246 = pkin(3) * t211 - pkin(10) * t215;
t231 = t246 * qJD(3);
t274 = t212 * t301;
t360 = (t231 - t274) * t204;
t210 = sin(qJ(4));
t214 = cos(qJ(4));
t315 = t204 * t215;
t151 = pkin(9) * t315 + (pkin(2) * t211 + pkin(10)) * t207;
t247 = -pkin(3) * t215 - pkin(10) * t211;
t152 = (-pkin(2) + t247) * t204;
t344 = t214 * t151 + t210 * t152;
t359 = qJD(4) * t344 - t361 * t210 - t360 * t214;
t291 = qJD(4) * t214;
t292 = qJD(4) * t210;
t358 = -t151 * t292 + t152 * t291 + t360 * t210 - t361 * t214;
t296 = qJD(2) * t215;
t280 = t204 * t296;
t357 = qJD(4) - t280;
t167 = -t214 * t207 + t210 * t316;
t293 = qJD(3) * t215;
t278 = t204 * t293;
t126 = -qJD(4) * t167 + t214 * t278;
t168 = t207 * t210 + t214 * t316;
t295 = qJD(3) * t211;
t279 = t204 * t295;
t356 = pkin(4) * t279 - t126 * qJ(5) - t168 * qJD(5) - t359;
t276 = t210 * t293;
t127 = qJD(4) * t168 + t204 * t276;
t355 = qJ(5) * t127 + qJD(5) * t167 - t358;
t255 = t210 * t280;
t335 = -qJ(5) - pkin(10);
t270 = qJD(4) * t335;
t299 = qJD(2) * t204;
t171 = pkin(9) * t299 + t274;
t161 = t211 * t171;
t208 = cos(pkin(6));
t300 = qJD(1) * t208;
t332 = qJD(2) * pkin(2);
t182 = t216 * t301 + t332;
t320 = t182 * t207;
t349 = t204 * t300 + t320;
t101 = t215 * t349 - t161;
t156 = t246 * t299;
t306 = t214 * t101 + t210 * t156;
t354 = qJ(5) * t255 + t214 * qJD(5) + t210 * t270 - t306;
t266 = -t210 * t101 + t214 * t156;
t309 = t214 * t215;
t353 = -t210 * qJD(5) + t214 * t270 - (pkin(4) * t211 - qJ(5) * t309) * t299 - t266;
t297 = qJD(2) * t207;
t194 = qJD(3) + t297;
t281 = t211 * t299;
t256 = t210 * t281;
t141 = -t214 * t194 + t256;
t143 = t194 * t210 + t214 * t281;
t203 = sin(pkin(13));
t206 = cos(pkin(13));
t265 = -t206 * t141 - t143 * t203;
t342 = qJD(6) - t265;
t209 = sin(qJ(6));
t213 = cos(qJ(6));
t239 = -t141 * t203 + t206 * t143;
t79 = t209 * t239 - t213 * t357;
t352 = t342 * t79;
t240 = -t209 * t357 - t213 * t239;
t351 = t240 * t342;
t172 = t203 * t210 - t206 * t214;
t129 = t172 * t280;
t166 = t172 * qJD(4);
t350 = t129 - t166;
t310 = t212 * t215;
t311 = t211 * t216;
t234 = t207 * t310 + t311;
t277 = t207 * t295;
t304 = pkin(2) * t277 + pkin(9) * t278 - t234 * t301;
t173 = t203 * t214 + t206 * t210;
t303 = t357 * t173;
t268 = t213 * t342;
t287 = qJD(2) * qJD(3);
t271 = t204 * t287;
t252 = t215 * t271;
t117 = -qJD(4) * t256 + t194 * t291 + t214 * t252;
t118 = (t211 * t291 + t276) * t299 + t194 * t292;
t72 = t117 * t203 + t206 * t118;
t328 = t209 * t72;
t348 = -t268 * t342 - t328;
t334 = t355 * t203 + t356 * t206;
t333 = t356 * t203 - t355 * t206;
t325 = t354 * t203 - t353 * t206;
t324 = t353 * t203 + t354 * t206;
t346 = t211 * t215;
t345 = pkin(4) * t127 + t304;
t102 = t215 * t171 + t211 * t349;
t343 = -t102 + (-t255 + t292) * pkin(4);
t198 = pkin(4) * t203 + pkin(11);
t251 = t211 * t271;
t122 = (t231 + t274) * t299;
t193 = t207 * t300;
t113 = t193 + (qJD(2) * t247 - t182) * t204;
t90 = pkin(10) * t194 + t102;
t62 = t113 * t210 + t214 * t90;
t222 = t232 * qJD(2);
t253 = t208 * t278;
t75 = (t182 * t313 - t161) * qJD(3) + (t205 * t222 + t253) * qJD(1);
t219 = -qJD(4) * t62 + t214 * t122 - t210 * t75;
t14 = pkin(4) * t251 - t117 * qJ(5) - t143 * qJD(5) + t219;
t229 = t113 * t291 + t210 * t122 + t214 * t75 - t292 * t90;
t18 = -qJ(5) * t118 - qJD(5) * t141 + t229;
t5 = t14 * t206 - t18 * t203;
t3 = -pkin(5) * t251 - t5;
t341 = (pkin(4) * t143 + pkin(5) * t239 - pkin(11) * t265 + qJD(6) * t198) * t342 + t3;
t283 = -pkin(4) * t214 - pkin(3);
t121 = pkin(5) * t172 - pkin(11) * t173 + t283;
t190 = t335 * t214;
t273 = t210 * t335;
t131 = -t206 * t190 + t203 * t273;
t340 = (-pkin(5) * t303 + pkin(11) * t350 + qJD(6) * t131 - t343) * t342 - t121 * t72;
t73 = t117 * t206 - t118 * t203;
t34 = -qJD(6) * t240 + t209 * t73 - t213 * t251;
t298 = qJD(2) * t205;
t272 = qJD(1) * t298;
t237 = t207 * t212 * t272;
t254 = t208 * t279;
t76 = qJD(1) * t254 + t171 * t293 + t182 * t277 + t215 * t237 + t272 * t311;
t59 = pkin(4) * t118 + t76;
t20 = pkin(5) * t72 - pkin(11) * t73 + t59;
t61 = t214 * t113 - t210 * t90;
t53 = -qJ(5) * t143 + t61;
t44 = pkin(4) * t357 + t53;
t54 = -qJ(5) * t141 + t62;
t50 = t206 * t54;
t22 = t203 * t44 + t50;
t17 = pkin(11) * t357 + t22;
t89 = -t194 * pkin(3) - t101;
t74 = t141 * pkin(4) + qJD(5) + t89;
t29 = -pkin(5) * t265 - pkin(11) * t239 + t74;
t243 = t17 * t209 - t213 * t29;
t6 = t203 * t14 + t206 * t18;
t4 = pkin(11) * t251 + t6;
t1 = -qJD(6) * t243 + t209 * t20 + t213 * t4;
t338 = -pkin(5) * t279 - t334;
t337 = t79 * t239;
t336 = t240 * t239;
t264 = -t151 * t210 + t214 * t152;
t78 = -pkin(4) * t315 - qJ(5) * t168 + t264;
t85 = -qJ(5) * t167 + t344;
t41 = t203 * t78 + t206 * t85;
t330 = t131 * t72;
t329 = t203 * t54;
t289 = qJD(6) * t213;
t290 = qJD(6) * t209;
t33 = t209 * t251 + t213 * t73 - t239 * t290 + t289 * t357;
t327 = t33 * t209;
t326 = pkin(5) * t281 + t325;
t323 = t141 * t357;
t322 = t143 * t357;
t321 = t173 * t213;
t319 = t357 * t210;
t318 = t357 * t214;
t200 = t204 ^ 2;
t217 = qJD(2) ^ 2;
t317 = t200 * t217;
t314 = t205 * t217;
t302 = t211 ^ 2 - t215 ^ 2;
t294 = qJD(3) * t214;
t288 = qJD(3) - t194;
t285 = t209 * t315;
t284 = t212 * t314;
t262 = t194 + t297;
t261 = 0.2e1 * t200 * t287;
t260 = t200 * t284;
t257 = t204 * t212 * t298;
t114 = t206 * t167 + t168 * t203;
t115 = -t167 * t203 + t168 * t206;
t150 = t195 + (-pkin(2) * t215 - pkin(3)) * t207;
t221 = t167 * pkin(4) + t150;
t55 = t114 * pkin(5) - t115 * pkin(11) + t221;
t249 = -pkin(11) * t279 - qJD(6) * t55 - t333;
t37 = -pkin(11) * t315 + t41;
t83 = t126 * t203 + t206 * t127;
t84 = t126 * t206 - t127 * t203;
t248 = -pkin(5) * t83 + pkin(11) * t84 + qJD(6) * t37 - t345;
t8 = t17 * t213 + t209 * t29;
t21 = t206 * t44 - t329;
t40 = -t203 * t85 + t206 * t78;
t235 = t207 * t308 - t312;
t124 = -t205 * t235 - t208 * t315;
t233 = t207 * t311 + t310;
t125 = t205 * t233 + t208 * t316;
t164 = -t204 * t205 * t216 + t207 * t208;
t91 = -t125 * t210 + t164 * t214;
t92 = t125 * t214 + t164 * t210;
t57 = t203 * t91 + t206 * t92;
t242 = t124 * t213 - t209 * t57;
t241 = t124 * t209 + t213 * t57;
t236 = t213 * t72 + (t209 * t265 - t290) * t342;
t93 = t115 * t209 + t213 * t315;
t109 = -t129 * t209 - t213 * t281;
t227 = -t166 * t209 + t173 * t289 - t109;
t111 = -t129 * t213 + t209 * t281;
t226 = -t166 * t213 - t173 * t290 - t111;
t137 = -t204 * t182 + t193;
t224 = qJD(3) * (t137 * t204 - t200 * t332);
t16 = -pkin(5) * t357 - t21;
t26 = t206 * t53 - t329;
t220 = -t198 * t72 + (t16 + t26) * t342;
t2 = -qJD(6) * t8 + t213 * t20 - t209 * t4;
t218 = -t330 + t3 * t173 + (pkin(11) * t281 - qJD(6) * t121 - t324) * t342;
t199 = -pkin(4) * t206 - pkin(5);
t130 = -t190 * t203 - t206 * t273;
t94 = t115 * t213 - t285;
t87 = t253 + (qJD(3) * t235 + t222) * t205;
t86 = t254 + (qJD(2) * t234 + qJD(3) * t233) * t205;
t56 = t203 * t92 - t206 * t91;
t49 = qJD(4) * t91 + t210 * t257 + t87 * t214;
t48 = -qJD(4) * t92 - t87 * t210 + t214 * t257;
t47 = -qJD(6) * t285 + t115 * t289 + t209 * t84 - t213 * t279;
t46 = -qJD(6) * t93 + t209 * t279 + t213 * t84;
t36 = pkin(5) * t315 - t40;
t25 = t203 * t53 + t50;
t24 = t203 * t48 + t206 * t49;
t23 = t203 * t49 - t206 * t48;
t7 = [0, 0, -t284, -t216 * t314, 0, 0, 0, 0, 0, t164 * t251 - t194 * t86 - t215 * t260, t164 * t252 - t194 * t87 + t211 * t260, 0, 0, 0, 0, 0, t118 * t124 + t141 * t86 + t251 * t91 + t357 * t48, t117 * t124 + t143 * t86 - t251 * t92 - t357 * t49, t23 * t239 + t24 * t265 + t56 * t73 - t57 * t72, t124 * t59 - t21 * t23 + t22 * t24 - t5 * t56 + t57 * t6 + t74 * t86, 0, 0, 0, 0, 0 (-qJD(6) * t241 - t209 * t24 + t86 * t213) * t342 + t242 * t72 + t23 * t79 + t56 * t34 -(qJD(6) * t242 + t86 * t209 + t213 * t24) * t342 - t241 * t72 - t23 * t240 + t56 * t33; 0, 0, 0, 0, t261 * t346, -t302 * t261, t262 * t278, -t262 * t279, 0, -t194 * t304 - t76 * t207 + t211 * t224, t361 * t194 - t75 * t207 + t215 * t224, t117 * t168 + t126 * t143, -t117 * t167 - t118 * t168 - t126 * t141 - t127 * t143, t126 * t357 + (-t117 * t215 + (qJD(2) * t168 + t143) * t295) * t204, -t127 * t357 + (t118 * t215 + (-qJD(2) * t167 - t141) * t295) * t204 (-t200 * t296 + t204 * t357) * t295, t150 * t118 + t89 * t127 + t76 * t167 - t359 * t357 + t304 * t141 + (-t219 * t215 + (qJD(2) * t264 + t61) * t295) * t204, t150 * t117 + t89 * t126 + t76 * t168 - t358 * t357 + t304 * t143 + (t229 * t215 + (-qJD(2) * t344 - t62) * t295) * t204, -t6 * t114 - t5 * t115 - t21 * t84 - t22 * t83 - t239 * t334 + t265 * t333 - t40 * t73 - t41 * t72, t334 * t21 + t333 * t22 + t59 * t221 + t345 * t74 + t5 * t40 + t6 * t41, -t240 * t46 + t33 * t94, t240 * t47 - t33 * t93 - t34 * t94 - t46 * t79, t114 * t33 - t240 * t83 + t342 * t46 + t72 * t94, -t114 * t34 - t342 * t47 - t72 * t93 - t79 * t83, t114 * t72 + t342 * t83 (-t209 * t37 + t213 * t55) * t72 + t2 * t114 - t243 * t83 + t36 * t34 + t3 * t93 + t16 * t47 + (t209 * t249 - t213 * t248) * t342 + t338 * t79 -(t209 * t55 + t213 * t37) * t72 - t1 * t114 - t8 * t83 + t36 * t33 + t3 * t94 + t16 * t46 + (t209 * t248 + t213 * t249) * t342 - t338 * t240; 0, 0, 0, 0, -t317 * t346, t302 * t317, t288 * t280, -t288 * t281, 0, t102 * t194 - t137 * t281 - t76, t101 * t194 + (qJD(3) * t171 + t237) * t211 + (-t137 * t299 - qJD(3) * t320 + (-qJD(3) * t204 * t208 - t216 * t298) * qJD(1)) * t215, t117 * t210 + t143 * t318 (t117 - t323) * t214 + (-t118 - t322) * t210, t357 * t291 + (-t357 * t309 + (qJD(3) * t210 - t143) * t211) * t299, -t357 * t292 + (t215 * t319 + (t141 + t294) * t211) * t299, -t357 * t281, -pkin(3) * t118 - t76 * t214 - t266 * t357 - t102 * t141 + (-pkin(10) * t318 + t210 * t89) * qJD(4) + (-t211 * t61 + (-pkin(10) * t295 - t215 * t89) * t210) * t299, -pkin(3) * t117 + t76 * t210 - t102 * t143 + t306 * t357 + (pkin(10) * t319 + t214 * t89) * qJD(4) + (-t89 * t309 + (-pkin(10) * t294 + t62) * t211) * t299, t130 * t73 - t6 * t172 - t5 * t173 - t21 * t350 - t303 * t22 + t325 * t239 + t324 * t265 - t330, -t5 * t130 + t6 * t131 - t325 * t21 + t324 * t22 + t59 * t283 + t343 * t74, -t226 * t240 + t321 * t33, -t240 * t109 + t111 * t79 - (t209 * t240 - t213 * t79) * t166 + (-t327 - t213 * t34 + (t209 * t79 + t213 * t240) * qJD(6)) * t173, t33 * t172 + t226 * t342 - t240 * t303 + t321 * t72, -t34 * t172 - t173 * t328 - t227 * t342 - t303 * t79, t72 * t172 + t303 * t342, t130 * t34 + t227 * t16 + t2 * t172 + t218 * t209 - t213 * t340 - t243 * t303 + t326 * t79, -t1 * t172 + t130 * t33 + t226 * t16 + t209 * t340 + t218 * t213 - t240 * t326 - t303 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 * t141, -t141 ^ 2 + t143 ^ 2, t117 + t323, -t118 + t322, t251, -t89 * t143 + t357 * t62 + t219, t141 * t89 + t357 * t61 - t229 (-t203 * t72 - t206 * t73) * pkin(4) + (t21 - t26) * t265 + (t22 - t25) * t239, t21 * t25 - t22 * t26 + (-t143 * t74 + t203 * t6 + t206 * t5) * pkin(4), -t240 * t268 + t327 (t33 - t352) * t213 + (-t34 + t351) * t209, t336 - t348, t236 + t337, -t342 * t239, t199 * t34 + t220 * t209 - t213 * t341 + t239 * t243 - t25 * t79, t199 * t33 + t209 * t341 + t220 * t213 + t8 * t239 + t240 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239 ^ 2 - t265 ^ 2, t21 * t239 - t22 * t265 + t59, 0, 0, 0, 0, 0, t236 - t337, t336 + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240 * t79, t240 ^ 2 - t79 ^ 2, t33 + t352, -t34 - t351, t72, t16 * t240 + t342 * t8 + t2, t16 * t79 - t243 * t342 - t1;];
tauc_reg  = t7;
