% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPRRP11
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
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:44
% EndTime: 2019-03-09 12:48:59
% DurationCPUTime: 8.16s
% Computational Cost: add. (10959->698), mult. (22703->816), div. (0->0), fcn. (14584->10), ass. (0->342)
t277 = sin(qJ(2));
t386 = qJD(1) * t277;
t253 = pkin(2) * t386;
t280 = cos(qJ(2));
t420 = qJ(3) * t280;
t319 = pkin(8) * t277 - t420;
t152 = qJD(1) * t319 + t253;
t385 = qJD(1) * t280;
t249 = pkin(7) * t385;
t199 = pkin(3) * t385 + t249;
t276 = sin(qJ(4));
t279 = cos(qJ(4));
t108 = -t152 * t276 + t199 * t279;
t407 = t276 * t277;
t315 = pkin(4) * t280 - pkin(9) * t407;
t378 = qJD(4) * t276;
t454 = pkin(2) + pkin(8);
t436 = pkin(9) + t454;
t479 = -qJD(1) * t315 + t378 * t436 - t108;
t109 = t152 * t279 + t199 * t276;
t205 = t436 * t279;
t357 = t279 * t386;
t478 = pkin(9) * t357 + qJD(4) * t205 + t109;
t350 = t276 * t385;
t381 = qJD(2) * t279;
t190 = -t350 + t381;
t275 = sin(qJ(5));
t372 = t276 * qJD(2);
t307 = t279 * t385 + t372;
t452 = cos(qJ(5));
t117 = t190 * t275 + t307 * t452;
t114 = t117 ^ 2;
t311 = t190 * t452 - t275 * t307;
t455 = t311 ^ 2;
t477 = -t114 + t455;
t453 = pkin(3) + pkin(7);
t375 = qJD(5) * t275;
t476 = pkin(4) * t375;
t475 = pkin(4) * t452;
t233 = qJD(4) + t386;
t220 = qJD(5) + t233;
t417 = t117 * t220;
t370 = qJD(1) * qJD(2);
t348 = t277 * t370;
t368 = t280 * qJDD(1);
t110 = qJD(4) * t307 - t279 * qJDD(2) + (-t348 + t368) * t276;
t376 = qJD(4) * t280;
t354 = t276 * t376;
t355 = t277 * t381;
t302 = t354 + t355;
t377 = qJD(4) * t279;
t362 = qJD(2) * t377 + qJDD(2) * t276 + t279 * t368;
t289 = qJD(1) * t302 - t362;
t349 = t452 * qJD(5);
t47 = t110 * t452 + t190 * t375 - t275 * t289 + t307 * t349;
t474 = -t47 + t417;
t473 = t117 * qJ(6);
t472 = t117 * t311;
t191 = t275 * t279 + t276 * t452;
t463 = qJD(4) + qJD(5);
t125 = t463 * t191;
t393 = t191 * t386 + t125;
t409 = t275 * t276;
t392 = -t275 * t378 - t276 * t375 + (qJD(4) * t452 + t349) * t279 + t357 * t452 - t386 * t409;
t411 = t220 * t311;
t48 = qJD(5) * t311 - t110 * t275 - t289 * t452;
t471 = -t48 + t411;
t248 = pkin(7) * t386;
t470 = qJD(3) + t248;
t271 = qJD(2) * qJ(3);
t167 = t271 + t199;
t127 = pkin(4) * t307 + t167;
t274 = qJ(4) + qJ(5);
t257 = sin(t274);
t258 = cos(t274);
t278 = sin(qJ(1));
t281 = cos(qJ(1));
t403 = t277 * t281;
t147 = t257 * t403 + t258 * t278;
t405 = t277 * t278;
t149 = -t257 * t405 + t258 * t281;
t267 = g(3) * t280;
t347 = t280 * t370;
t369 = t277 * qJDD(1);
t305 = t347 + t369;
t188 = qJDD(4) + t305;
t231 = pkin(7) * t347;
t245 = pkin(7) * t369;
t346 = qJDD(3) + t231 + t245;
t113 = pkin(3) * t305 - qJDD(2) * t454 + t346;
t232 = pkin(2) * t348;
t379 = qJD(3) * t277;
t296 = qJD(2) * t319 - t379;
t259 = t277 * qJ(3);
t345 = -pkin(1) - t259;
t304 = -t280 * t454 + t345;
t86 = qJD(1) * t296 + qJDD(1) * t304 + t232;
t138 = t304 * qJD(1);
t371 = pkin(3) * t386 + t470;
t145 = -qJD(2) * t454 + t371;
t89 = t138 * t279 + t145 * t276;
t34 = -qJD(4) * t89 + t113 * t279 - t276 * t86;
t20 = pkin(4) * t188 + pkin(9) * t110 + t34;
t365 = t113 * t276 + t145 * t377 + t279 * t86;
t33 = -t138 * t378 + t365;
t24 = pkin(9) * t289 + t33;
t88 = -t138 * t276 + t145 * t279;
t77 = -pkin(9) * t190 + t88;
t69 = pkin(4) * t233 + t77;
t78 = -pkin(9) * t307 + t89;
t3 = t20 * t275 + t24 * t452 + t349 * t69 - t375 * t78;
t292 = g(1) * t147 - g(2) * t149 - t257 * t267 - t3;
t469 = t117 * t127 + t292;
t261 = t279 * pkin(4);
t244 = t261 + pkin(3);
t204 = t436 * t276;
t432 = t204 * t375 - t205 * t349 + t275 * t479 - t452 * t478;
t129 = -t204 * t452 - t205 * t275;
t431 = -qJD(5) * t129 + t275 * t478 + t452 * t479;
t468 = t233 * t89 + t34;
t467 = qJ(6) * t311;
t382 = qJD(2) * t277;
t198 = t453 * t382;
t246 = pkin(7) * t368;
t269 = qJDD(2) * qJ(3);
t270 = qJD(2) * qJD(3);
t361 = t246 + t269 + t270;
t332 = pkin(3) * t368 + t361;
t115 = -qJD(1) * t198 + t332;
t465 = qJD(4) * t233 * t454 + t115;
t177 = qJDD(5) + t188;
t335 = pkin(4) * t349;
t450 = pkin(4) * t275;
t464 = -t177 * t450 - t220 * t335;
t264 = t280 * pkin(2);
t390 = t264 + t259;
t448 = pkin(8) * t280;
t327 = t390 + t448;
t178 = -pkin(1) - t327;
t215 = t453 * t277;
t193 = t276 * t215;
t124 = t178 * t279 + t193;
t391 = pkin(4) * t377 + t244 * t386 + t470;
t321 = g(1) * t281 + g(2) * t278;
t162 = t177 * pkin(5);
t424 = t47 * qJ(6);
t462 = -qJD(6) * t311 + t162 + t424;
t146 = -t257 * t278 + t258 * t403;
t148 = t257 * t281 + t258 * t405;
t461 = -g(1) * t146 - g(2) * t148 + t258 * t267;
t397 = t279 * t281;
t169 = -t276 * t278 + t277 * t397;
t401 = t278 * t279;
t171 = t276 * t281 + t277 * t401;
t398 = t279 * t280;
t460 = -g(1) * t169 - g(2) * t171 + g(3) * t398;
t15 = t117 * t392 + t191 * t48;
t343 = t20 * t452 - t275 * t24;
t76 = t452 * t78;
t36 = t275 * t69 + t76;
t4 = -qJD(5) * t36 + t343;
t287 = t4 + t461;
t459 = -t127 * t311 + t287;
t358 = t452 * t279;
t192 = t358 - t409;
t14 = -t192 * t47 - t311 * t393;
t342 = -pkin(5) * t117 - qJD(6);
t71 = t127 - t342;
t458 = -t311 * t71 + t461;
t1 = t4 + t462;
t374 = qJD(6) * t117;
t430 = qJ(6) * t48;
t2 = t3 - t374 - t430;
t74 = t275 * t78;
t35 = t452 * t69 - t74;
t25 = t35 - t467;
t22 = pkin(5) * t220 + t25;
t26 = t36 - t473;
t360 = -g(1) * t403 - g(2) * t405 + t267;
t457 = t1 * t192 + t191 * t2 - t22 * t393 + t26 * t392 + t360;
t456 = t191 * t3 + t192 * t4 - t35 * t393 + t36 * t392 + t360;
t282 = -pkin(9) - pkin(8);
t451 = pkin(4) * t190;
t449 = pkin(4) * t276;
t445 = g(1) * t278;
t440 = g(2) * t281;
t439 = g(3) * t277;
t435 = -t25 + t22;
t434 = -qJ(6) * t392 - qJD(6) * t191 + t432;
t433 = -pkin(5) * t385 + qJ(6) * t393 - t192 * qJD(6) + t431;
t44 = t452 * t77 - t74;
t105 = -pkin(9) * t398 + t124;
t194 = t279 * t215;
t344 = pkin(9) * t280 - t178;
t97 = pkin(4) * t277 + t276 * t344 + t194;
t63 = t105 * t452 + t275 * t97;
t427 = t233 * t88;
t425 = t276 * t88;
t423 = pkin(5) * t392 + t391;
t422 = -t117 * t335 - t450 * t48;
t421 = pkin(7) * qJDD(2);
t419 = qJDD(2) * pkin(2);
t414 = t188 * t276;
t413 = t307 * t279;
t412 = t190 * t307;
t268 = -qJ(6) + t282;
t410 = t268 * t280;
t408 = t276 * t307;
t406 = t276 * t280;
t404 = t277 * t279;
t285 = qJD(1) ^ 2;
t402 = t277 * t285;
t400 = t278 * t280;
t399 = t279 * t110;
t161 = t279 * t188;
t396 = t280 * t281;
t395 = t280 * t282;
t394 = t454 * t188;
t202 = pkin(5) * t258 + t261;
t216 = t453 * t280;
t389 = pkin(1) * t281 + pkin(7) * t278;
t272 = t277 ^ 2;
t273 = t280 ^ 2;
t388 = t272 - t273;
t387 = t272 + t273;
t384 = qJD(2) * t307;
t383 = qJD(2) * t190;
t380 = qJD(2) * t280;
t367 = pkin(4) * t406;
t364 = t307 * t404;
t363 = t276 * t403;
t165 = pkin(4) * t398 + t216;
t356 = t277 * t372;
t353 = t279 * t376;
t351 = t220 * t385;
t239 = qJ(3) + t449;
t43 = -t275 * t77 - t76;
t62 = -t105 * t275 + t452 * t97;
t341 = -qJD(2) * pkin(2) + qJD(3);
t340 = qJD(1) * t124 + t89;
t252 = pkin(2) * t382;
t135 = t252 + t296;
t200 = t453 * t380;
t338 = -t135 * t276 + t200 * t279;
t128 = t204 * t275 - t205 * t452;
t337 = -qJD(1) * t216 - t167;
t336 = t233 + t386;
t334 = pkin(2) * t396 + qJ(3) * t403 + t389;
t333 = -t245 - t360;
t331 = t280 * t358;
t330 = t452 * t382;
t329 = t277 * t347;
t237 = g(1) * t400;
t328 = -g(2) * t396 + t237;
t326 = t362 * t276;
t325 = t387 * qJDD(1) * pkin(7);
t284 = qJD(2) ^ 2;
t324 = pkin(7) * t284 + t440;
t323 = g(1) * t148 - g(2) * t146;
t322 = -g(1) * t149 - g(2) * t147;
t320 = t177 * t192 - t220 * t393;
t318 = t190 * t279 - t408;
t203 = t248 + t341;
t214 = -t249 - t271;
t317 = t203 * t280 + t214 * t277;
t316 = t233 * t276;
t314 = t345 - t264;
t312 = -0.2e1 * pkin(1) * t370 - t421;
t168 = t314 * qJD(1);
t310 = t168 * t386 + qJDD(3) - t333;
t309 = -t233 * t377 - t414;
t58 = t315 * qJD(2) + (t279 * t344 - t193) * qJD(4) + t338;
t64 = t135 * t279 - t178 * t378 + t200 * t276 + t215 * t377;
t60 = pkin(9) * t302 + t64;
t10 = -t105 * t375 + t275 * t58 + t349 * t97 + t452 * t60;
t308 = -qJ(3) * t380 - t379;
t306 = -t177 * t191 - t220 * t392;
t303 = 0.2e1 * qJDD(1) * pkin(1) - t324;
t299 = t279 * t348 - t362;
t206 = -pkin(1) - t390;
t298 = t421 + (-qJD(1) * t206 - t168) * qJD(2);
t297 = -t280 * t321 - t439;
t295 = t297 * t257;
t294 = t297 * t258;
t111 = qJD(1) * t308 + qJDD(1) * t314 + t232;
t157 = t252 + t308;
t293 = qJD(1) * t157 + qJDD(1) * t206 + t111 + t324;
t130 = -pkin(4) * t354 + (-pkin(7) - t244) * t382;
t11 = -qJD(5) * t63 - t275 * t60 + t452 * t58;
t141 = pkin(7) * t348 - t361;
t151 = t346 - t419;
t290 = qJD(2) * t317 - t141 * t280 + t151 * t277;
t288 = t292 + t430;
t68 = pkin(4) * t362 + qJD(1) * t130 + t332;
t21 = t48 * pkin(5) + qJDD(6) + t68;
t265 = t281 * pkin(7);
t243 = pkin(5) + t475;
t230 = qJ(3) * t396;
t228 = qJ(3) * t400;
t225 = t280 * t402;
t209 = t388 * t285;
t208 = qJDD(2) * t280 - t277 * t284;
t207 = qJDD(2) * t277 + t280 * t284;
t201 = pkin(5) * t257 + t449;
t196 = pkin(3) + t202;
t195 = -qJ(3) * t385 + t253;
t180 = qJDD(1) * t273 - 0.2e1 * t329;
t179 = qJDD(1) * t272 + 0.2e1 * t329;
t172 = -t276 * t405 + t397;
t170 = t363 + t401;
t154 = t191 * t280;
t153 = t275 * t406 - t331;
t142 = pkin(5) * t191 + t239;
t134 = 0.2e1 * t277 * t368 - 0.2e1 * t370 * t388;
t123 = -t178 * t276 + t194;
t122 = t177 * t277 + t220 * t380;
t112 = -pkin(5) * t153 + t165;
t95 = -qJ(6) * t191 + t129;
t94 = -qJ(6) * t192 + t128;
t85 = pkin(5) * t311 + t451;
t80 = t125 * t280 - t275 * t356 + t279 * t330;
t79 = t276 * t330 + (t406 * t463 + t355) * t275 - t463 * t331;
t65 = -qJD(4) * t124 + t338;
t61 = -pkin(5) * t80 + t130;
t52 = qJ(6) * t153 + t63;
t51 = pkin(5) * t277 + qJ(6) * t154 + t62;
t41 = -qJD(2) * t311 + t306;
t40 = -qJD(2) * t117 + t320;
t38 = -t311 * t385 + t320;
t37 = t117 * t385 + t306;
t28 = t44 - t467;
t27 = t43 + t473;
t17 = -t117 * t80 - t153 * t48;
t16 = t154 * t47 + t311 * t79;
t13 = -t154 * t177 + t220 * t79 - t277 * t47 + t311 * t380;
t12 = -t117 * t380 + t153 * t177 + t220 * t80 - t277 * t48;
t9 = qJ(6) * t80 + qJD(6) * t153 + t10;
t8 = pkin(5) * t380 - t79 * qJ(6) + t154 * qJD(6) + t11;
t7 = -t117 * t79 - t153 * t47 + t154 * t48 + t311 * t80;
t6 = -t14 - t15;
t5 = t117 * t393 + t191 * t47 - t192 * t48 - t311 * t392;
t18 = [0, 0, 0, 0, 0, qJDD(1), -t440 + t445, t321, 0, 0, t179, t134, t207, t180, t208, 0, t277 * t312 + t280 * t303 + t237, t312 * t280 + (-t303 - t445) * t277, -t321 + 0.2e1 * t325, -g(1) * (-pkin(1) * t278 + t265) - g(2) * t389 + (pkin(7) ^ 2 * t387 + pkin(1) ^ 2) * qJDD(1), 0, -t207, -t208, t179, t134, t180, t325 + t290 - t321, t277 * t298 + t280 * t293 - t237, t298 * t280 + (-t293 + t445) * t277, pkin(7) * t290 - g(1) * t265 - g(2) * t334 + t111 * t206 + t168 * t157 - t314 * t445, t110 * t406 + (-t353 + t356) * t190, t318 * t382 + (-t299 * t276 + t399 + (t413 + (t190 - t350) * t276) * qJD(4)) * t280 (t233 * t372 - t110) * t277 + (t309 + t383) * t280, -t289 * t398 - t302 * t307 (t336 * t381 - t362) * t277 + (t336 * t378 - t161 - t384) * t280, t188 * t277 + t233 * t380, t65 * t233 + t123 * t188 - t198 * t307 + t216 * t362 - g(1) * t172 - g(2) * t170 + (t337 * t381 + t34) * t277 + (t88 * qJD(2) + t115 * t279 + t337 * t378) * t280, g(1) * t171 - g(2) * t169 - t110 * t216 - t124 * t188 - t190 * t198 - t233 * t64 + (t167 * t372 - t33) * t277 + (-qJD(2) * t89 - t115 * t276 - t167 * t377) * t280, -t64 * t307 - t124 * t362 - t65 * t190 + t123 * t110 + t237 + (t279 * t340 - t425) * t382 + (-t440 + t34 * t276 - t33 * t279 + (t276 * t340 + t279 * t88) * qJD(4)) * t280, t33 * t124 + t89 * t64 + t34 * t123 + t88 * t65 + t115 * t216 - t167 * t198 - g(1) * (pkin(3) * t281 + t265) - g(2) * (pkin(8) * t396 + t334) + (-g(1) * (t314 - t448) - g(2) * pkin(3)) * t278, t16, t7, t13, t17, t12, t122, t11 * t220 + t117 * t130 - t127 * t80 - t153 * t68 + t165 * t48 + t177 * t62 + t277 * t4 + t35 * t380 + t322, -t10 * t220 + t127 * t79 + t130 * t311 - t154 * t68 - t165 * t47 - t177 * t63 - t277 * t3 - t36 * t380 + t323, -t10 * t117 - t11 * t311 + t153 * t3 + t154 * t4 - t35 * t79 + t36 * t80 + t47 * t62 - t48 * t63 + t328, t3 * t63 + t36 * t10 + t4 * t62 + t35 * t11 + t68 * t165 + t127 * t130 - g(1) * (t244 * t281 + t265) - g(2) * (pkin(4) * t363 - t281 * t395 + t334) + (-g(1) * (-pkin(4) * t407 + t314 + t395) - g(2) * t244) * t278, t16, t7, t13, t17, t12, t122, t1 * t277 + t112 * t48 + t117 * t61 - t153 * t21 + t177 * t51 + t22 * t380 + t220 * t8 - t71 * t80 + t322, -t112 * t47 - t154 * t21 - t177 * t52 - t2 * t277 - t220 * t9 - t26 * t380 + t311 * t61 + t71 * t79 + t323, t1 * t154 - t117 * t9 + t153 * t2 - t22 * t79 + t26 * t80 - t311 * t8 + t47 * t51 - t48 * t52 + t328, t2 * t52 + t26 * t9 + t1 * t51 + t22 * t8 + t21 * t112 + t71 * t61 - g(1) * (t196 * t281 + t265) - g(2) * (t201 * t403 - t268 * t396 + t334) + (-g(1) * (-t201 * t277 + t314 + t410) - g(2) * t196) * t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t209, t369, t225, t368, qJDD(2), pkin(1) * t402 + t333, t439 - t246 + (pkin(1) * t285 + t321) * t280, 0, 0, qJDD(2), -t369, -t368, -t225, t209, t225 (-pkin(2) * t277 + t420) * qJDD(1) + ((-t214 - t271) * t277 + (-t203 + t341) * t280) * qJD(1), -t195 * t385 + t310 - 0.2e1 * t419, t246 + 0.2e1 * t269 + 0.2e1 * t270 + (qJD(1) * t195 - g(3)) * t277 + (qJD(1) * t168 - t321) * t280, -t141 * qJ(3) - t214 * qJD(3) - t151 * pkin(2) - t168 * t195 - g(1) * (-pkin(2) * t403 + t230) - g(2) * (-pkin(2) * t405 + t228) - g(3) * t390 - t317 * qJD(1) * pkin(7), -t190 * t316 - t399, -t362 * t279 + t110 * t276 - t318 * qJD(4) + (t276 * t353 + (t408 + (-t190 + t381) * t279) * t277) * qJD(1), -t233 * t378 + t161 + (-t190 * t280 - t233 * t407) * qJD(1), t326 + t307 * t377 + (-t276 * t302 + t364) * qJD(1) (-t233 * t404 + t280 * t307) * qJD(1) + t309, -t233 * t385, qJ(3) * t362 - t108 * t233 - t88 * t385 + t371 * t307 + (qJD(4) * t167 - t394 + (t167 - t271) * t386) * t279 + (-t439 + (-qJ(3) * qJD(1) * qJD(4) - t321) * t280 + t465) * t276, t89 * t385 - qJ(3) * t110 + t109 * t233 + t371 * t190 + (-t167 * t233 + t394) * t276 + (t297 + t465) * t279, t108 * t190 + t109 * t307 + (-t89 * t386 - t454 * t110 - t34 + (t307 * t454 - t89) * qJD(4)) * t279 + (-t454 * t299 - t33 + t88 * t386 + (t88 - (t190 + t350) * t454) * qJD(4)) * t276 - t360, -g(1) * t230 - g(2) * t228 - g(3) * t327 + t115 * qJ(3) - t88 * t108 - t89 * t109 + t371 * t167 + (-t33 * t276 - t34 * t279 - (t279 * t89 - t425) * qJD(4) + t321 * t277) * t454, t14, t5, t38, t15, t37, -t351, t117 * t391 + t127 * t392 + t128 * t177 + t191 * t68 + t220 * t431 + t239 * t48 - t35 * t385 + t295, -t127 * t393 - t129 * t177 + t192 * t68 - t220 * t432 - t239 * t47 + t311 * t391 + t36 * t385 + t294, -t117 * t432 + t128 * t47 - t129 * t48 - t311 * t431 - t456, t3 * t129 + t4 * t128 + t68 * t239 - g(1) * (t281 * t367 + t230) - g(2) * (t278 * t367 + t228) - g(3) * (t390 - t395) + t432 * t36 + t431 * t35 + (-g(3) * t449 + t321 * (pkin(2) - t282)) * t277 + t391 * t127, t14, t5, t38, t15, t37, -t351, t117 * t423 + t142 * t48 + t177 * t94 + t191 * t21 - t22 * t385 + t220 * t433 + t392 * t71 + t295, -t142 * t47 - t177 * t95 + t192 * t21 - t220 * t434 + t26 * t385 + t311 * t423 - t393 * t71 + t294, -t117 * t434 - t311 * t433 + t47 * t94 - t48 * t95 - t457, t2 * t95 + t1 * t94 + t21 * t142 - g(1) * (t201 * t396 + t230) - g(2) * (t201 * t400 + t228) - g(3) * (t390 - t410) + t423 * t71 + (-g(3) * t201 + t321 * (pkin(2) - t268)) * t277 + t434 * t26 + t433 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, qJDD(2) + t225, -t272 * t285 - t284, qJD(2) * t214 + t231 + t310 - t419, 0, 0, 0, 0, 0, 0, -t233 * t316 + t161 - t384, -t233 ^ 2 * t279 - t383 - t414, -t326 + t399 + (t190 * t276 - t413) * qJD(4) + (-t364 + (t354 + (t190 + t381) * t277) * t276) * qJD(1), -qJD(2) * t167 + t468 * t279 + (t33 - t427) * t276 + t360, 0, 0, 0, 0, 0, 0, t40, t41, t6, -qJD(2) * t127 + t456, 0, 0, 0, 0, 0, 0, t40, t41, t6, -qJD(2) * t71 + t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t412, t190 ^ 2 - t307 ^ 2, t233 * t307 - t110, -t412, t190 * t233 + t289, t188, -t167 * t190 + t460 + t468, g(1) * t170 - g(2) * t172 + t167 * t307 + t427 + (qJD(4) * t138 - t267) * t276 - t365, 0, 0, t472, t477, t474, -t472, t471, t177, -t43 * t220 + (-t117 * t190 + t177 * t452 - t220 * t375) * pkin(4) + t459, t220 * t44 - t311 * t451 + t464 + t469, t47 * t475 + t422 + (t36 + t43 + t476) * t311 + (t44 - t35) * t117, -t35 * t43 - t36 * t44 + (t3 * t275 + t4 * t452 - t127 * t190 + (-t275 * t35 + t36 * t452) * qJD(5) + t460) * pkin(4), t472, t477, t474, -t472, t471, t177, -t85 * t117 + t243 * t177 - t27 * t220 + (-t76 + (-pkin(4) * t220 - t69) * t275) * qJD(5) + t343 + t458 + t462, t117 * t71 + t220 * t28 - t311 * t85 + t288 + t374 + t464, t243 * t47 + t422 + (t26 + t27 + t476) * t311 + (t28 - t22) * t117, t1 * t243 - t26 * t28 - t22 * t27 - t71 * t85 - g(1) * (-t201 * t278 + t202 * t403) - g(2) * (t201 * t281 + t202 * t405) + t202 * t267 + (t2 * t275 + (-t22 * t275 + t26 * t452) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t472, t477, t474, -t472, t471, t177, t36 * t220 + t459, t220 * t35 + t469, 0, 0, t472, t477, t474, -t472, t471, t177, t424 + t26 * t220 + 0.2e1 * t162 + (t342 - t71) * t311 + t287, -pkin(5) * t455 + t220 * t25 + (qJD(6) + t71) * t117 + t288, pkin(5) * t47 - t117 * t435, t435 * t26 + (t1 + t458) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 + t411, -t47 - t417, -t114 - t455, t26 * t117 + t22 * t311 + t21 + t297;];
tau_reg  = t18;
