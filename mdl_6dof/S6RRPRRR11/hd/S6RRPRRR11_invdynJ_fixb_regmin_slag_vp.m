% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:54
% EndTime: 2019-03-09 14:33:11
% DurationCPUTime: 7.28s
% Computational Cost: add. (7406->552), mult. (15508->727), div. (0->0), fcn. (10826->14), ass. (0->291)
t250 = cos(qJ(4));
t245 = sin(qJ(4));
t347 = qJD(2) * t245;
t251 = cos(qJ(2));
t350 = qJD(1) * t251;
t181 = t250 * t350 + t347;
t320 = t245 * t350;
t345 = qJD(2) * t250;
t183 = -t320 + t345;
t244 = sin(qJ(5));
t249 = cos(qJ(5));
t109 = t249 * t181 + t183 * t244;
t243 = sin(qJ(6));
t248 = cos(qJ(6));
t290 = t181 * t244 - t249 * t183;
t292 = t109 * t243 + t248 * t290;
t49 = t248 * t109 - t243 * t290;
t419 = t292 * t49;
t414 = t292 ^ 2 - t49 ^ 2;
t246 = sin(qJ(2));
t351 = qJD(1) * t246;
t214 = qJD(4) + t351;
t207 = qJD(5) + t214;
t199 = qJD(6) + t207;
t332 = qJD(1) * qJD(2);
t319 = t246 * t332;
t330 = t251 * qJDD(1);
t422 = t319 - t330;
t97 = -qJD(4) * t181 + t250 * qJDD(2) + t245 * t422;
t98 = -qJD(4) * t320 + t245 * qJDD(2) + (qJD(2) * qJD(4) - t422) * t250;
t262 = qJD(5) * t290 - t244 * t97 - t249 * t98;
t337 = qJD(5) * t249;
t338 = qJD(5) * t244;
t31 = -t181 * t337 - t183 * t338 - t244 * t98 + t249 * t97;
t335 = qJD(6) * t248;
t336 = qJD(6) * t243;
t8 = -t109 * t335 + t243 * t262 + t248 * t31 + t290 * t336;
t412 = t199 * t49 + t8;
t242 = qJ(4) + qJ(5);
t232 = qJ(6) + t242;
t217 = sin(t232);
t218 = cos(t232);
t247 = sin(qJ(1));
t252 = cos(qJ(1));
t364 = t246 * t252;
t127 = t217 * t364 + t218 * t247;
t366 = t246 * t247;
t129 = -t217 * t366 + t218 * t252;
t253 = -pkin(2) - pkin(8);
t367 = t246 * qJ(3);
t173 = t253 * t251 - pkin(1) - t367;
t132 = t173 * qJD(1);
t223 = pkin(7) * t351;
t416 = qJD(3) + t223;
t333 = pkin(3) * t351 + t416;
t139 = t253 * qJD(2) + t333;
t75 = -t132 * t245 + t250 * t139;
t62 = -pkin(9) * t183 + t75;
t55 = pkin(4) * t214 + t62;
t76 = t132 * t250 + t139 * t245;
t63 = -pkin(9) * t181 + t76;
t59 = t249 * t63;
t25 = t244 * t55 + t59;
t421 = pkin(10) * t109;
t19 = t25 - t421;
t16 = t19 * t336;
t236 = g(3) * t251;
t224 = pkin(7) * t350;
t191 = pkin(3) * t350 + t224;
t239 = qJD(2) * qJ(3);
t162 = t239 + t191;
t117 = pkin(4) * t181 + t162;
t60 = pkin(5) * t109 + t117;
t411 = g(1) * t127 - g(2) * t129 - t217 * t236 + t49 * t60 + t16;
t126 = -t217 * t247 + t218 * t364;
t128 = t217 * t252 + t218 * t366;
t318 = t251 * t332;
t331 = t246 * qJDD(1);
t277 = t318 + t331;
t180 = qJDD(4) + t277;
t172 = qJDD(5) + t180;
t212 = pkin(7) * t318;
t220 = pkin(7) * t331;
t329 = qJDD(3) + t220;
t317 = t212 + t329;
t104 = t277 * pkin(3) + t253 * qJDD(2) + t317;
t213 = pkin(2) * t319;
t381 = qJ(3) * t251;
t295 = pkin(8) * t246 - t381;
t334 = t246 * qJD(3);
t267 = qJD(2) * t295 - t334;
t71 = qJD(1) * t267 + qJDD(1) * t173 + t213;
t308 = t250 * t104 - t245 * t71;
t264 = -qJD(4) * t76 + t308;
t12 = t180 * pkin(4) - t97 * pkin(9) + t264;
t341 = qJD(4) * t250;
t328 = -t245 * t104 - t139 * t341 - t250 * t71;
t342 = qJD(4) * t245;
t17 = -pkin(9) * t98 - t132 * t342 - t328;
t314 = t249 * t12 - t244 * t17;
t265 = -qJD(5) * t25 + t314;
t2 = t172 * pkin(5) - t31 * pkin(10) + t265;
t307 = -t244 * t12 - t249 * t17 - t55 * t337 + t63 * t338;
t3 = pkin(10) * t262 - t307;
t324 = t248 * t2 - t243 * t3;
t425 = -g(1) * t126 - g(2) * t128 + t218 * t236 + t60 * t292 + t324;
t263 = qJD(6) * t292 - t243 * t31 + t248 * t262;
t406 = -t199 * t292 + t263;
t228 = pkin(2) * t351;
t145 = qJD(1) * t295 + t228;
t170 = t250 * t191;
t369 = t245 * t246;
t286 = pkin(4) * t251 - pkin(9) * t369;
t391 = pkin(9) - t253;
t424 = -qJD(1) * t286 + t245 * t145 + t391 * t342 - t170;
t195 = t391 * t250;
t323 = t250 * t351;
t356 = t250 * t145 + t245 * t191;
t423 = pkin(9) * t323 + qJD(4) * t195 + t356;
t184 = t244 * t250 + t245 * t249;
t279 = t246 * t184;
t400 = qJD(4) + qJD(5);
t358 = -qJD(1) * t279 - t400 * t184;
t397 = pkin(3) + pkin(7);
t420 = pkin(10) * t290;
t417 = t290 * t109;
t362 = t249 * t250;
t372 = t244 * t245;
t357 = -t244 * t342 - t245 * t338 + t249 * t323 - t351 * t372 + t362 * t400;
t297 = pkin(2) * t251 + t367;
t285 = pkin(1) + t297;
t163 = t285 * qJD(1);
t415 = qJDD(1) * t285;
t413 = -t109 ^ 2 + t290 ^ 2;
t410 = t109 * t207 + t31;
t57 = t244 * t63;
t24 = t249 * t55 - t57;
t18 = t24 + t420;
t14 = pkin(5) * t207 + t18;
t384 = t248 * t19;
t5 = t243 * t14 + t384;
t409 = -qJD(6) * t5 + t425;
t230 = sin(t242);
t231 = cos(t242);
t141 = t230 * t364 + t231 * t247;
t143 = -t230 * t366 + t231 * t252;
t408 = g(1) * t141 - g(2) * t143 + t109 * t117 - t230 * t236 + t307;
t140 = -t230 * t247 + t231 * t364;
t142 = t230 * t252 + t231 * t366;
t407 = -g(1) * t140 - g(2) * t142 + t117 * t290 + t231 * t236 + t265;
t405 = -t207 * t290 + t262;
t404 = t424 * t249;
t403 = t358 * t248;
t201 = t397 * t246;
t186 = t245 * t201;
t355 = t250 * t173 + t186;
t154 = t250 * t180;
t402 = -t214 * t342 + t154;
t288 = -t362 + t372;
t114 = -t184 * t243 - t248 * t288;
t194 = t391 * t245;
t354 = -t249 * t194 - t244 * t195;
t325 = -pkin(4) * t250 - pkin(3);
t353 = pkin(4) * t341 - t325 * t351 + t416;
t401 = -t194 * t338 + t195 * t337 - t424 * t244 + t423 * t249;
t299 = g(1) * t247 - g(2) * t252;
t300 = g(1) * t252 + g(2) * t247;
t399 = t162 * t214 + t253 * t180;
t398 = t400 * t251;
t392 = g(3) * t246;
t113 = t248 * t184 - t243 * t288;
t390 = -qJD(6) * t113 - t357 * t243 + t403;
t389 = t114 * qJD(6) + t358 * t243 + t357 * t248;
t388 = t249 * t62 - t57;
t187 = t250 * t201;
t316 = pkin(9) * t251 - t173;
t86 = t246 * pkin(4) + t316 * t245 + t187;
t361 = t250 * t251;
t93 = -pkin(9) * t361 + t355;
t386 = t244 * t86 + t249 * t93;
t385 = t357 * pkin(5) + t353;
t383 = t97 * t250;
t382 = pkin(7) * qJDD(2);
t380 = qJDD(2) * pkin(2);
t155 = qJDD(6) + t172;
t379 = t113 * t155;
t378 = t114 * t155;
t377 = t155 * t244;
t376 = t181 * t214;
t375 = t183 * t214;
t219 = pkin(4) * t249 + pkin(5);
t373 = t219 * t155;
t371 = t244 * t248;
t370 = t245 * t180;
t368 = t245 * t251;
t365 = t246 * t250;
t363 = t247 * t250;
t360 = t250 * t252;
t216 = t245 * pkin(4) + qJ(3);
t202 = t397 * t251;
t240 = t246 ^ 2;
t241 = t251 ^ 2;
t352 = t240 - t241;
t349 = qJD(2) * t181;
t348 = qJD(2) * t183;
t346 = qJD(2) * t246;
t344 = qJD(2) * t251;
t343 = qJD(4) * t132;
t340 = qJD(4) * t251;
t339 = qJD(4) * t253;
t255 = qJD(1) ^ 2;
t327 = t246 * t255 * t251;
t221 = pkin(7) * t330;
t237 = qJDD(2) * qJ(3);
t238 = qJD(2) * qJD(3);
t326 = t221 + t237 + t238;
t158 = pkin(4) * t361 + t202;
t321 = t245 * t340;
t315 = qJD(6) * t14 + t3;
t227 = pkin(2) * t346;
t124 = t227 + t267;
t192 = t397 * t344;
t305 = -t245 * t124 + t250 * t192;
t39 = t286 * qJD(2) + (t316 * t250 - t186) * qJD(4) + t305;
t276 = t250 * t124 - t173 * t342 + t245 * t192 + t201 * t341;
t43 = (t246 * t345 + t321) * pkin(9) + t276;
t312 = -t244 * t43 + t249 * t39;
t311 = -t244 * t62 - t59;
t310 = -t244 * t93 + t249 * t86;
t306 = -qJD(2) * pkin(2) + qJD(3);
t304 = t194 * t244 - t249 * t195;
t190 = t397 * t346;
t84 = -pkin(10) * t184 + t354;
t303 = pkin(5) * t350 + t358 * pkin(10) + t354 * qJD(5) + qJD(6) * t84 - t244 * t423 - t404;
t83 = pkin(10) * t288 + t304;
t302 = t357 * pkin(10) - qJD(6) * t83 + t401;
t301 = -t288 * t172 + t358 * t207;
t298 = qJD(6) * t288 - t357;
t296 = pkin(2) * t246 - t381;
t146 = t288 * t251;
t147 = t184 * t251;
t291 = t248 * t146 + t147 * t243;
t82 = t146 * t243 - t147 * t248;
t193 = t223 + t306;
t200 = -t224 - t239;
t289 = t193 * t251 + t200 * t246;
t287 = t214 * t245;
t283 = -0.2e1 * pkin(1) * t332 - t382;
t282 = t244 * t39 + t249 * t43 + t86 * t337 - t93 * t338;
t281 = -t214 * t341 - t370;
t280 = -qJ(3) * t344 - t334;
t278 = -t184 * t172 - t357 * t207;
t274 = pkin(1) * t255 + t300;
t254 = qJD(2) ^ 2;
t273 = pkin(7) * t254 - t299;
t272 = 0.2e1 * t163 * qJD(2) + t382;
t269 = -t251 * t300 - t392;
t268 = 0.2e1 * qJDD(1) * pkin(1) - t273;
t107 = pkin(3) * t330 - qJD(1) * t190 + t326;
t266 = t107 + t269;
t261 = -t163 * t351 - t246 * t300 + t236 + t329;
t118 = -pkin(4) * t321 + (-pkin(7) + t325) * t346;
t54 = pkin(4) * t98 + t107;
t150 = t227 + t280;
t99 = qJD(1) * t280 + t213 - t415;
t259 = qJD(1) * t150 + t273 - t415 + t99;
t135 = pkin(7) * t319 - t326;
t144 = t317 - t380;
t256 = qJD(2) * t289 - t135 * t251 + t144 * t246 - t300;
t188 = -qJ(3) * t350 + t228;
t167 = -t245 * t366 + t360;
t166 = t245 * t252 + t246 * t363;
t165 = t245 * t364 + t363;
t164 = -t245 * t247 + t246 * t360;
t136 = pkin(5) * t184 + t216;
t103 = -pkin(5) * t146 + t158;
t70 = pkin(4) * t183 - pkin(5) * t290;
t65 = t184 * t398 - t288 * t346;
t64 = qJD(2) * t279 + t288 * t398;
t44 = -t65 * pkin(5) + t118;
t34 = pkin(10) * t146 + t386;
t33 = pkin(5) * t246 + pkin(10) * t147 + t310;
t23 = t388 + t420;
t22 = t311 + t421;
t21 = qJD(6) * t82 + t243 * t64 - t248 * t65;
t20 = qJD(6) * t291 + t243 * t65 + t248 * t64;
t13 = -pkin(5) * t262 + t54;
t7 = pkin(10) * t65 + t282;
t6 = pkin(5) * t344 - t64 * pkin(10) - qJD(5) * t386 + t312;
t4 = t14 * t248 - t19 * t243;
t1 = [qJDD(1), t299, t300, qJDD(1) * t240 + 0.2e1 * t246 * t318, 0.2e1 * t246 * t330 - 0.2e1 * t352 * t332, qJDD(2) * t246 + t251 * t254, qJDD(2) * t251 - t246 * t254, 0, t246 * t283 + t251 * t268, -t246 * t268 + t251 * t283 (t240 + t241) * qJDD(1) * pkin(7) + t256, t246 * t272 + t251 * t259, -t246 * t259 + t251 * t272, pkin(7) * t256 - t163 * t150 + (t299 - t99) * t285, -t97 * t368 + (t245 * t346 - t250 * t340) * t183 (-t181 * t245 + t183 * t250) * t346 + (t245 * t98 - t383 + (t181 * t250 + t183 * t245) * qJD(4)) * t251 (t214 * t347 + t97) * t246 + (t281 + t348) * t251 (t214 * t345 - t98) * t246 + (-t349 - t402) * t251, t180 * t246 + t214 * t344, t305 * t214 + (-t173 * t245 + t187) * t180 + t308 * t246 - t190 * t181 + t202 * t98 + t107 * t361 - g(1) * t167 - g(2) * t165 + (-t162 * t365 + t251 * t75) * qJD(2) + (-t162 * t368 - t355 * t214 - t76 * t246) * qJD(4), -t276 * t214 - t355 * t180 - t190 * t183 + t202 * t97 + g(1) * t166 - g(2) * t164 + ((qJD(2) * t162 + t343) * t245 + t328) * t246 + (-qJD(2) * t76 - t107 * t245 - t162 * t341) * t251, -t147 * t31 - t290 * t64, -t109 * t64 + t146 * t31 - t147 * t262 - t290 * t65, -t147 * t172 + t207 * t64 + t246 * t31 - t290 * t344, -t109 * t344 + t146 * t172 + t207 * t65 + t246 * t262, t172 * t246 + t207 * t344, t312 * t207 + t310 * t172 + t314 * t246 + t24 * t344 + t118 * t109 - t158 * t262 - t54 * t146 - t117 * t65 - g(1) * t143 - g(2) * t141 + (-t207 * t386 - t246 * t25) * qJD(5), g(1) * t142 - g(2) * t140 + t117 * t64 - t118 * t290 - t54 * t147 + t158 * t31 - t386 * t172 - t282 * t207 + t307 * t246 - t25 * t344, -t20 * t292 + t8 * t82, -t20 * t49 + t21 * t292 + t263 * t82 + t291 * t8, t155 * t82 + t199 * t20 + t246 * t8 - t292 * t344, t155 * t291 - t199 * t21 + t246 * t263 - t344 * t49, t155 * t246 + t199 * t344 (-t243 * t7 + t248 * t6) * t199 + (-t243 * t34 + t248 * t33) * t155 + t324 * t246 + t4 * t344 + t44 * t49 - t103 * t263 - t13 * t291 + t60 * t21 - g(1) * t129 - g(2) * t127 + ((-t243 * t33 - t248 * t34) * t199 - t5 * t246) * qJD(6), -t5 * t344 + g(1) * t128 - g(2) * t126 + t103 * t8 + t13 * t82 + t16 * t246 + t60 * t20 - t44 * t292 + (-(-qJD(6) * t34 + t6) * t199 - t33 * t155 - t2 * t246) * t243 + (-(qJD(6) * t33 + t7) * t199 - t34 * t155 - t315 * t246) * t248; 0, 0, 0, -t327, t352 * t255, t331, t330, qJDD(2), t246 * t274 - t220 - t236, t251 * t274 - t221 + t392, -t296 * qJDD(1) + ((-t200 - t239) * t246 + (-t193 + t306) * t251) * qJD(1), -t188 * t350 + t261 - 0.2e1 * t380, t221 + 0.2e1 * t237 + 0.2e1 * t238 + (qJD(1) * t188 - g(3)) * t246 + (-qJD(1) * t163 - t300) * t251, -pkin(7) * qJD(1) * t289 - t144 * pkin(2) - g(3) * t297 - t135 * qJ(3) - t200 * qJD(3) + t163 * t188 + t296 * t300, -t183 * t287 + t383 (-t98 - t375) * t250 + (-t97 + t376) * t245 (-t183 * t251 - t214 * t369) * qJD(1) + t402 (t181 * t251 - t214 * t365) * qJD(1) + t281, -t214 * t350, -t75 * t350 + qJ(3) * t98 - t170 * t214 + t333 * t181 + t399 * t250 + ((t145 - t339) * t214 + t266) * t245, qJ(3) * t97 + t356 * t214 + t76 * t350 + t333 * t183 - t399 * t245 + (-t214 * t339 + t266) * t250, -t288 * t31 - t290 * t358, -t358 * t109 - t31 * t184 - t262 * t288 + t290 * t357, t290 * t350 + t301, t109 * t350 + t278, -t207 * t350, t304 * t172 - t216 * t262 + t54 * t184 - t24 * t350 + (t194 * t337 + (qJD(5) * t195 + t423) * t244 + t404) * t207 + t357 * t117 + t353 * t109 + t269 * t230, t358 * t117 - t354 * t172 + t207 * t401 + t216 * t31 + t269 * t231 + t25 * t350 - t288 * t54 - t290 * t353, t8 * t114 - t292 * t390, -t8 * t113 + t114 * t263 + t292 * t389 - t390 * t49, t199 * t390 + t292 * t350 + t378, -t199 * t389 + t350 * t49 - t379, -t199 * t350 (-t243 * t84 + t248 * t83) * t155 - t136 * t263 + t13 * t113 - t4 * t350 + t389 * t60 + t385 * t49 + (t243 * t302 - t248 * t303) * t199 + t269 * t217 -(t243 * t83 + t248 * t84) * t155 + t136 * t8 + t13 * t114 + t5 * t350 + t390 * t60 - t385 * t292 + (t243 * t303 + t248 * t302) * t199 + t269 * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, qJDD(2) + t327, -t240 * t255 - t254, qJD(2) * t200 + t212 + t261 - t380, 0, 0, 0, 0, 0, -t214 * t287 + t154 - t349, -t214 ^ 2 * t250 - t348 - t370, 0, 0, 0, 0, 0, -qJD(2) * t109 + t301, qJD(2) * t290 + t278, 0, 0, 0, 0, 0, t378 - qJD(2) * t49 + (-t184 * t335 + t243 * t298 + t403) * t199, -t379 + qJD(2) * t292 + (t298 * t248 + (qJD(6) * t184 - t358) * t243) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183 * t181, -t181 ^ 2 + t183 ^ 2, t97 + t376, t375 - t98, t180, -g(1) * t164 - g(2) * t166 + g(3) * t361 - t162 * t183 + t214 * t76 + t264, g(1) * t165 - g(2) * t167 + t162 * t181 + t214 * t75 + (t343 - t236) * t245 + t328, -t417, t413, t410, t405, t172, -t311 * t207 + (-t109 * t183 + t172 * t249 - t207 * t338) * pkin(4) + t407, t388 * t207 + (-t172 * t244 + t183 * t290 - t207 * t337) * pkin(4) + t408, -t419, t414, t412, t406, t155, t248 * t373 - (t22 * t248 - t23 * t243) * t199 - t70 * t49 + (-t243 * t377 + (-t243 * t249 - t371) * t199 * qJD(5)) * pkin(4) + ((-pkin(4) * t371 - t219 * t243) * t199 - t5) * qJD(6) + t425, t70 * t292 + (-t373 - t2 + (t22 - (-qJD(5) - qJD(6)) * t244 * pkin(4)) * t199) * t243 + (-pkin(4) * t377 + (-pkin(4) * t337 - qJD(6) * t219 + t23) * t199 - t315) * t248 + t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, t413, t410, t405, t172, t207 * t25 + t407, t207 * t24 + t408, -t419, t414, t412, t406, t155 -(-t18 * t243 - t384) * t199 + (t155 * t248 - t199 * t336 + t290 * t49) * pkin(5) + t409 (-t19 * t199 - t2) * t243 + (t18 * t199 - t315) * t248 + (-t155 * t243 - t199 * t335 - t290 * t292) * pkin(5) + t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t419, t414, t412, t406, t155, t199 * t5 + t409, t199 * t4 - t243 * t2 - t248 * t315 + t411;];
tau_reg  = t1;
