% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:36
% EndTime: 2019-03-09 12:53:53
% DurationCPUTime: 8.03s
% Computational Cost: add. (7242->582), mult. (14553->708), div. (0->0), fcn. (9291->10), ass. (0->294)
t218 = sin(qJ(2));
t335 = qJD(1) * t218;
t176 = qJD(4) + t335;
t385 = cos(qJ(5));
t296 = t385 * qJD(5);
t411 = t176 * t385 + t296;
t196 = pkin(2) * t335;
t221 = cos(qJ(2));
t360 = qJ(3) * t221;
t273 = pkin(8) * t218 - t360;
t103 = qJD(1) * t273 + t196;
t334 = qJD(1) * t221;
t192 = pkin(7) * t334;
t140 = pkin(3) * t334 + t192;
t220 = cos(qJ(4));
t123 = t220 * t140;
t217 = sin(qJ(4));
t352 = t217 * t218;
t268 = pkin(4) * t221 - pkin(9) * t352;
t327 = qJD(4) * t217;
t224 = -pkin(2) - pkin(8);
t375 = pkin(9) - t224;
t410 = -qJD(1) * t268 + t103 * t217 + t375 * t327 - t123;
t144 = t375 * t220;
t304 = t220 * t335;
t362 = t220 * t103 + t217 * t140;
t409 = pkin(9) * t304 + qJD(4) * t144 + t362;
t319 = qJD(1) * qJD(2);
t294 = t221 * t319;
t316 = t218 * qJDD(1);
t252 = t294 + t316;
t130 = qJDD(4) + t252;
t125 = qJDD(5) + t130;
t216 = sin(qJ(5));
t305 = t385 * t220;
t133 = t216 * t217 - t305;
t159 = qJD(5) + t176;
t323 = qJD(5) * t216;
t326 = qJD(4) * t220;
t370 = -t220 * t323 - t411 * t217 + (-t304 - t326) * t216;
t280 = -t133 * t125 + t159 * t370;
t321 = t217 * qJD(2);
t131 = -t220 * t334 - t321;
t300 = t217 * t334;
t331 = qJD(2) * t220;
t132 = -t300 + t331;
t72 = -t385 * t131 + t132 * t216;
t408 = -qJD(2) * t72 + t280;
t215 = qJ(4) + qJ(5);
t200 = cos(t215);
t219 = sin(qJ(1));
t222 = cos(qJ(1));
t274 = g(1) * t222 + g(2) * t219;
t378 = g(3) * t218;
t392 = t221 * t274 + t378;
t143 = t375 * t217;
t81 = -t143 * t385 - t216 * t144;
t407 = t81 * t125 + t200 * t392;
t386 = pkin(3) + pkin(7);
t258 = t216 * t131 + t132 * t385;
t212 = qJD(2) * qJ(3);
t115 = t212 + t140;
t79 = -pkin(4) * t131 + t115;
t25 = pkin(5) * t72 - qJ(6) * t258 + t79;
t406 = t25 * t72;
t405 = t72 * t79;
t311 = t216 * t352;
t369 = -qJD(1) * t311 - t216 * t327 - t217 * t323 + t220 * t411;
t377 = t258 * t72;
t315 = t221 * qJDD(1);
t266 = qJD(2) * qJD(4) + t315;
t191 = pkin(7) * t335;
t404 = qJD(3) + t191;
t387 = t258 ^ 2;
t403 = -t72 ^ 2 + t387;
t284 = -qJD(4) + t335;
t242 = qJD(2) * t284 - t315;
t318 = qJD(1) * qJD(4);
t271 = t221 * t318 - qJDD(2);
t253 = t271 * t220;
t229 = t242 * t217 - t253;
t325 = qJD(4) * t221;
t302 = t217 * t325;
t303 = t218 * t331;
t248 = t302 + t303;
t308 = t217 * qJDD(2) + t220 * t266;
t233 = qJD(1) * t248 - t308;
t20 = -t131 * t296 + t132 * t323 - t216 * t233 - t385 * t229;
t13 = t159 * t72 - t20;
t37 = pkin(5) * t258 + qJ(6) * t72;
t187 = pkin(4) * t220 + pkin(3);
t363 = pkin(4) * t326 + t187 * t335 + t404;
t257 = t216 * t143 - t144 * t385;
t401 = -qJD(5) * t257 - t216 * t410 + t409 * t385;
t400 = -qJD(5) * t81 + t409 * t216 + t385 * t410;
t109 = t125 * qJ(6);
t148 = t159 * qJD(6);
t399 = t109 + t148;
t134 = t216 * t220 + t217 * t385;
t255 = -t125 * t134 - t159 * t369;
t398 = -qJD(2) * t258 + t255;
t343 = t220 * t222;
t117 = -t217 * t219 + t218 * t343;
t346 = t219 * t220;
t119 = t217 * t222 + t218 * t346;
t396 = -g(1) * t117 - g(2) * t119;
t111 = t125 * pkin(5);
t395 = t111 - qJDD(6);
t199 = sin(t215);
t350 = t218 * t219;
t100 = t199 * t222 + t200 * t350;
t209 = g(3) * t221;
t295 = t218 * t319;
t175 = pkin(2) * t295;
t329 = qJD(3) * t218;
t240 = qJD(2) * t273 - t329;
t201 = t218 * qJ(3);
t291 = -pkin(1) - t201;
t251 = t221 * t224 + t291;
t51 = qJD(1) * t240 + qJDD(1) * t251 + t175;
t174 = pkin(7) * t294;
t188 = pkin(7) * t316;
t293 = qJDD(3) + t174 + t188;
t69 = pkin(3) * t252 + qJDD(2) * t224 + t293;
t320 = pkin(3) * t335 + t404;
t97 = qJD(2) * t224 + t320;
t314 = t217 * t69 + t220 * t51 + t97 * t326;
t89 = t251 * qJD(1);
t11 = pkin(9) * t233 - t327 * t89 + t314;
t53 = -t217 * t89 + t220 * t97;
t43 = -pkin(9) * t132 + t53;
t36 = pkin(4) * t176 + t43;
t54 = t217 * t97 + t220 * t89;
t44 = pkin(9) * t131 + t54;
t241 = -t266 + t295;
t357 = qJD(4) * t89;
t66 = t220 * t69;
t8 = t130 * pkin(4) + t66 + (pkin(9) * t271 - t357) * t220 + (-pkin(9) * t241 - qJD(4) * t97 - t51) * t217;
t299 = t216 * t11 + t44 * t296 + t36 * t323 - t385 * t8;
t348 = t218 * t222;
t98 = t199 * t219 - t200 * t348;
t245 = g(1) * t98 - g(2) * t100 + t200 * t209 - t299;
t232 = t25 * t258 - t245 - t395;
t394 = -t79 * t258 + t245;
t21 = qJD(5) * t258 + t216 * t229 - t385 * t233;
t393 = t159 * t258 - t21;
t391 = t133 * t20 + t258 * t370;
t154 = t386 * t218;
t135 = t217 * t154;
t330 = qJD(2) * t221;
t141 = t386 * t330;
t332 = qJD(2) * t218;
t195 = pkin(2) * t332;
t86 = t195 + t240;
t288 = t220 * t141 - t217 * t86;
t206 = t221 * pkin(2);
t337 = t206 + t201;
t147 = -pkin(1) - t337;
t126 = -pkin(8) * t221 + t147;
t290 = pkin(9) * t221 - t126;
t26 = t268 * qJD(2) + (t220 * t290 - t135) * qJD(4) + t288;
t136 = t220 * t154;
t58 = pkin(4) * t218 + t217 * t290 + t136;
t338 = t220 * t126 + t135;
t344 = t220 * t221;
t64 = -pkin(9) * t344 + t338;
t262 = t216 * t58 + t385 * t64;
t313 = t217 * t141 + t154 * t326 + t220 * t86;
t29 = pkin(9) * t248 - t126 * t327 + t313;
t390 = -qJD(5) * t262 - t216 * t29 + t26 * t385;
t389 = t221 * (qJD(4) + qJD(5));
t292 = t385 * t11 + t216 * t8 + t36 * t296 - t44 * t323;
t1 = t292 + t399;
t366 = t216 * t44;
t16 = t36 * t385 - t366;
t339 = qJD(6) - t16;
t14 = -t159 * pkin(5) + t339;
t312 = t385 * t44;
t17 = t216 * t36 + t312;
t15 = t159 * qJ(6) + t17;
t2 = t299 - t395;
t306 = -g(1) * t348 - g(2) * t350 + t209;
t388 = t1 * t134 + t133 * t2 - t14 * t370 + t15 * t369 + t306;
t383 = g(1) * t219;
t379 = g(2) * t222;
t202 = t217 * pkin(4);
t374 = pkin(5) * t369 - qJ(6) * t370 + qJD(6) * t133 + t363;
t373 = -qJ(6) * t334 - t401;
t372 = pkin(5) * t334 - t400;
t367 = t159 * t17;
t19 = t385 * t43 - t366;
t364 = pkin(4) * t296 + qJD(6) - t19;
t361 = pkin(7) * qJDD(2);
t356 = qJDD(2) * pkin(2);
t354 = t132 * t176;
t353 = t217 * t130;
t351 = t217 * t221;
t349 = t218 * t220;
t226 = qJD(1) ^ 2;
t347 = t218 * t226;
t345 = t219 * t221;
t110 = t220 * t130;
t342 = t221 * t222;
t223 = -pkin(9) - pkin(8);
t341 = t221 * t223;
t340 = t224 * t130;
t183 = qJ(3) + t202;
t155 = t386 * t221;
t213 = t218 ^ 2;
t214 = t221 ^ 2;
t336 = t213 - t214;
t333 = qJD(2) * t131;
t328 = qJD(4) * t131;
t324 = qJD(4) * t224;
t322 = t115 * qJD(4);
t210 = qJDD(2) * qJ(3);
t310 = t217 * t348;
t309 = t221 * t347;
t189 = pkin(7) * t315;
t211 = qJD(2) * qJD(3);
t307 = t189 + t210 + t211;
t177 = pkin(4) * t344;
t114 = t177 + t155;
t301 = t220 * t325;
t298 = g(3) * t337;
t289 = -t217 * t51 + t66;
t287 = -qJD(2) * pkin(2) + qJD(3);
t286 = -qJD(1) * t155 - t115;
t285 = t176 + t335;
t283 = t222 * pkin(1) + pkin(2) * t342 + t219 * pkin(7) + qJ(3) * t348;
t282 = -t188 - t306;
t281 = pkin(3) * t315 + t307;
t139 = t386 * t332;
t18 = t216 * t43 + t312;
t278 = pkin(4) * t323 - t18;
t277 = g(1) * t100 + g(2) * t98;
t101 = -t199 * t350 + t200 * t222;
t99 = t199 * t348 + t200 * t219;
t276 = -g(1) * t101 - g(2) * t99;
t225 = qJD(2) ^ 2;
t275 = pkin(7) * t225 + t379;
t272 = t126 * t176 + t89 * t218;
t142 = t191 + t287;
t153 = -t192 - t212;
t270 = t142 * t221 + t153 * t218;
t269 = t176 ^ 2;
t267 = t291 - t206;
t264 = -t216 * t64 + t385 * t58;
t116 = t267 * qJD(1);
t261 = t116 * t335 + qJDD(3) - t282;
t260 = -0.2e1 * pkin(1) * t319 - t361;
t259 = t216 * t26 + t385 * t29 + t58 * t296 - t323 * t64;
t256 = (-pkin(5) * t200 - qJ(6) * t199) * t221;
t254 = -qJ(3) * t330 - t329;
t250 = 0.2e1 * qJDD(1) * pkin(1) - t275;
t249 = pkin(5) * t199 - qJ(6) * t200 + t202;
t246 = t361 + (-qJD(1) * t147 - t116) * qJD(2);
t244 = g(1) * t99 - g(2) * t101 - t199 * t209 - t292;
t239 = -g(1) * (-t98 * pkin(5) + qJ(6) * t99) - g(2) * (t100 * pkin(5) - qJ(6) * t101);
t108 = t195 + t254;
t67 = qJD(1) * t254 + qJDD(1) * t267 + t175;
t236 = qJD(1) * t108 + qJDD(1) * t147 + t275 + t67;
t235 = t159 * t16 + t244;
t82 = -pkin(4) * t302 + (-pkin(7) - t187) * t332;
t102 = t293 - t356;
t94 = pkin(7) * t295 - t307;
t234 = qJD(2) * t270 + t102 * t218 - t94 * t221;
t231 = t125 * t257 - t199 * t392;
t70 = -qJD(1) * t139 + t281;
t230 = t70 + (-qJ(3) * t318 - t274) * t221 - t378;
t35 = pkin(4) * t308 + qJD(1) * t82 + t281;
t207 = t222 * pkin(7);
t186 = -pkin(4) * t385 - pkin(5);
t182 = pkin(4) * t216 + qJ(6);
t180 = g(1) * t345;
t173 = qJ(3) * t342;
t171 = qJ(3) * t345;
t137 = -qJ(3) * t334 + t196;
t120 = -t217 * t350 + t343;
t118 = t310 + t346;
t105 = t134 * t221;
t104 = t216 * t351 - t221 * t305;
t68 = pkin(5) * t134 + qJ(6) * t133 + t183;
t47 = -pkin(5) * t104 + qJ(6) * t105 + t114;
t46 = -qJD(2) * t311 + t134 * t389 + t385 * t303;
t45 = t133 * t389 + t134 * t332;
t32 = pkin(4) * t132 + t37;
t30 = -t218 * pkin(5) - t264;
t28 = qJ(6) * t218 + t262;
t12 = -pkin(5) * t46 - qJ(6) * t45 + qJD(6) * t105 + t82;
t5 = -pkin(5) * t330 - t390;
t4 = qJ(6) * t330 + qJD(6) * t218 + t259;
t3 = t21 * pkin(5) + t20 * qJ(6) - qJD(6) * t258 + t35;
t6 = [qJDD(1), -t379 + t383, t274, qJDD(1) * t213 + 0.2e1 * t218 * t294, 0.2e1 * t218 * t315 - 0.2e1 * t319 * t336, qJDD(2) * t218 + t221 * t225, qJDD(2) * t221 - t218 * t225, 0, t218 * t260 + t221 * t250 + t180, t260 * t221 + (-t250 - t383) * t218 (t213 + t214) * qJDD(1) * pkin(7) + t234 - t274, t218 * t246 + t221 * t236 - t180, t246 * t221 + (-t236 + t383) * t218, pkin(7) * t234 - g(1) * t207 - g(2) * t283 + t116 * t108 + t67 * t147 - t267 * t383, -t132 * t301 + (t132 * t332 + t221 * t253 - t241 * t351) * t217 (t217 * t131 + t132 * t220) * t332 + (((t132 - t300) * qJD(4) + t308) * t217 + (-t328 + t253 + (t315 + (qJD(4) - 0.2e1 * t335) * qJD(2)) * t217) * t220) * t221, t132 * t330 + (-t176 * t325 - t218 * t271) * t220 + ((-t130 - t316) * t221 + (t176 + t284) * t332) * t217 (t285 * t331 - t308) * t218 + (t285 * t327 - t110 + t333) * t221, t130 * t218 + t176 * t330, t288 * t176 + (-t126 * t217 + t136) * t130 + t289 * t218 + t139 * t131 + t155 * t308 + t70 * t344 - g(1) * t120 - g(2) * t118 + (t53 * t221 + t286 * t349) * qJD(2) + (-t272 * t220 + (-t154 * t176 - t97 * t218 + t221 * t286) * t217) * qJD(4), -t313 * t176 - t338 * t130 - t314 * t218 - t54 * t330 - t139 * t132 + g(1) * t119 - g(2) * t117 + (-t155 * t271 - t221 * t322) * t220 + ((-qJDD(1) * t155 - t70) * t221 + t272 * qJD(4) + (t115 * t218 + t155 * t284) * qJD(2)) * t217, t105 * t20 + t258 * t45, -t104 * t20 + t105 * t21 + t258 * t46 - t45 * t72, -t105 * t125 + t159 * t45 - t20 * t218 + t258 * t330, t104 * t125 + t159 * t46 - t21 * t218 - t330 * t72, t125 * t218 + t159 * t330, -t35 * t104 + t114 * t21 + t264 * t125 + t159 * t390 + t16 * t330 - t299 * t218 - t79 * t46 + t82 * t72 + t276, -t35 * t105 - t114 * t20 - t125 * t262 - t159 * t259 - t17 * t330 - t218 * t292 + t258 * t82 + t79 * t45 + t277, -t104 * t3 + t12 * t72 - t125 * t30 - t14 * t330 - t159 * t5 - t2 * t218 + t21 * t47 - t25 * t46 + t276, -g(2) * t342 + t1 * t104 - t105 * t2 + t14 * t45 + t15 * t46 - t20 * t30 - t21 * t28 + t258 * t5 - t4 * t72 + t180, t1 * t218 + t105 * t3 - t12 * t258 + t125 * t28 + t15 * t330 + t159 * t4 + t20 * t47 - t25 * t45 - t277, t1 * t28 + t15 * t4 + t3 * t47 + t25 * t12 + t2 * t30 + t14 * t5 - g(1) * (pkin(5) * t101 + qJ(6) * t100 + t187 * t222 + t207) - g(2) * (pkin(4) * t310 + pkin(5) * t99 + qJ(6) * t98 - t222 * t341 + t283) + (-g(1) * (-pkin(4) * t352 + t267 + t341) - g(2) * t187) * t219; 0, 0, 0, -t309, t336 * t226, t316, t315, qJDD(2), pkin(1) * t347 + t282, t378 - t189 + (pkin(1) * t226 + t274) * t221 (-pkin(2) * t218 + t360) * qJDD(1) + ((-t153 - t212) * t218 + (-t142 + t287) * t221) * qJD(1), -t137 * t334 + t261 - 0.2e1 * t356, t189 + 0.2e1 * t210 + 0.2e1 * t211 + (qJD(1) * t137 - g(3)) * t218 + (qJD(1) * t116 - t274) * t221, -t94 * qJ(3) - t153 * qJD(3) - t102 * pkin(2) - t116 * t137 - g(1) * (-pkin(2) * t348 + t173) - g(2) * (-pkin(2) * t350 + t171) - t298 - t270 * qJD(1) * pkin(7), -t271 * t220 ^ 2 + (t220 * t242 - t354) * t217 (-t132 * qJD(4) + (-t132 + t331) * t335 - t308) * t220 + (-t328 - t220 * qJDD(2) + t266 * t217 + (0.2e1 * t301 + (-t131 - t321) * t218) * qJD(1)) * t217, -t176 * t327 + t110 + (-t132 * t221 - t176 * t352) * qJD(1), -t176 * t326 - t353 + (-t131 * t221 - t176 * t349) * qJD(1), -t176 * t334, qJ(3) * t308 - t123 * t176 - t53 * t334 - t320 * t131 + (t322 + t340 + (t115 - t212) * t335) * t220 + ((t103 - t324) * t176 + t230) * t217, t362 * t176 + t54 * t334 + t320 * t132 + (qJ(3) * t241 - t115 * t176 - t340) * t217 + (-t176 * t324 + t210 + t230) * t220, t391, t133 * t21 + t134 * t20 - t258 * t369 - t370 * t72, -t258 * t334 + t280, t334 * t72 + t255, -t159 * t334, t35 * t134 + t159 * t400 - t16 * t334 + t183 * t21 + t363 * t72 + t369 * t79 + t231, -t35 * t133 + t159 * t401 + t17 * t334 - t183 * t20 + t363 * t258 + t370 * t79 - t407, t134 * t3 + t14 * t334 - t159 * t372 + t21 * t68 + t25 * t369 + t374 * t72 + t231, t20 * t257 - t21 * t81 + t258 * t372 - t373 * t72 - t388, t133 * t3 - t15 * t334 + t373 * t159 + t20 * t68 - t370 * t25 - t374 * t258 + t407, t1 * t81 + t3 * t68 - t2 * t257 - g(1) * t173 - g(2) * t171 - t298 + t374 * t25 + t373 * t15 + t372 * t14 + (g(3) * t223 - t249 * t274) * t221 + (-g(3) * t249 + t274 * (pkin(2) - t223)) * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, qJDD(2) + t309, -t213 * t226 - t225, qJD(2) * t153 + t174 + t261 - t356, 0, 0, 0, 0, 0, -t217 * t269 + t110 + t333, -qJD(2) * t132 - t220 * t269 - t353, 0, 0, 0, 0, 0, t408, t398, t408, -t134 * t21 - t369 * t72 - t391, -t398, -qJD(2) * t25 + t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132 * t131, -t131 ^ 2 + t132 ^ 2, -t131 * t176 + t229, t233 + t354, t130, g(3) * t344 - t115 * t132 + t289 + (-qJD(4) + t176) * t54 + t396, g(1) * t118 - g(2) * t120 - t115 * t131 + t176 * t53 + (t357 - t209) * t217 - t314, t377, t403, t13, t393, t125, t18 * t159 + (t125 * t385 - t132 * t72 - t159 * t323) * pkin(4) + t394, t19 * t159 + t405 + (-t125 * t216 - t132 * t258 - t159 * t296) * pkin(4) + t244, -t125 * t186 - t159 * t278 - t32 * t72 - t232, -t182 * t21 - t186 * t20 + (t15 + t278) * t258 + (t14 - t364) * t72, t125 * t182 + t159 * t364 + t258 * t32 - t244 + t399 - t406, t1 * t182 + t2 * t186 - t25 * t32 - t14 * t18 - g(3) * (-t177 + t256) + t364 * t15 + (t14 * t323 + t396) * pkin(4) + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t377, t403, t13, t393, t125, t367 + t394, t235 + t405, -t37 * t72 + t111 - t232 + t367, pkin(5) * t20 - qJ(6) * t21 + (t15 - t17) * t258 + (t14 - t339) * t72, t258 * t37 + 0.2e1 * t109 + 0.2e1 * t148 - t235 - t406, -t2 * pkin(5) - g(3) * t256 + t1 * qJ(6) - t14 * t17 + t15 * t339 - t25 * t37 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 + t377, t13, -t159 ^ 2 - t387, -t15 * t159 + t232;];
tau_reg  = t6;
