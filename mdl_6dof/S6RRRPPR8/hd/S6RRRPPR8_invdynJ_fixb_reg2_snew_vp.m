% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRPPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:16:34
% EndTime: 2019-05-07 06:17:14
% DurationCPUTime: 12.36s
% Computational Cost: add. (31661->551), mult. (68923->684), div. (0->0), fcn. (52754->10), ass. (0->329)
t292 = cos(qJ(2));
t288 = sin(qJ(3));
t291 = cos(qJ(3));
t286 = cos(pkin(6));
t335 = qJD(1) * t286 + qJD(2);
t285 = sin(pkin(6));
t289 = sin(qJ(2));
t360 = qJD(1) * t289;
t345 = t285 * t360;
t246 = t288 * t335 + t291 * t345;
t361 = qJD(1) * t285;
t342 = qJD(2) * t361;
t354 = qJDD(1) * t285;
t257 = t289 * t354 + t292 * t342;
t280 = qJDD(1) * t286 + qJDD(2);
t337 = t288 * t257 - t291 * t280;
t191 = qJD(3) * t246 + t337;
t344 = t292 * t361;
t273 = -qJD(3) + t344;
t226 = t246 * t273;
t442 = t191 - t226;
t271 = t273 ^ 2;
t244 = t288 * t345 - t291 * t335;
t422 = t244 ^ 2;
t439 = -t271 - t422;
t202 = t246 * t244;
t258 = -t289 * t342 + t292 * t354;
t252 = -qJDD(3) + t258;
t441 = t202 + t252;
t455 = t441 * t288;
t475 = t439 * t291 + t455;
t494 = t289 * t475;
t527 = pkin(1) * (t292 * t442 - t494);
t436 = t422 - t271;
t176 = t252 - t202;
t478 = t176 * t291;
t131 = -t288 * t436 + t478;
t310 = (-qJD(3) - t273) * t246 - t337;
t479 = t176 * t288;
t513 = (t289 * (t291 * t436 + t479) - t292 * t310) * t285 - t286 * t131;
t161 = (qJD(3) - t273) * t246 + t337;
t347 = -t244 * qJD(3) + t291 * t257 + t288 * t280;
t377 = t273 * t244;
t154 = -t377 - t347;
t398 = t154 * t288;
t421 = t246 ^ 2;
t435 = t422 - t421;
t99 = t154 * t291 + t161 * t288;
t526 = (t289 * (t161 * t291 - t398) - t292 * t435) * t285 + t286 * t99;
t440 = -t271 - t421;
t120 = t291 * t440 + t479;
t123 = -t288 * t440 + t478;
t514 = -pkin(1) * t120 + pkin(8) * (t123 * t292 - t154 * t289);
t493 = t292 * t475;
t454 = t441 * t291;
t474 = t439 * t288 - t454;
t498 = pkin(1) * t474;
t525 = -pkin(8) * (t289 * t442 + t493) + t498;
t524 = pkin(1) * (t123 * t289 + t154 * t292);
t155 = -t377 + t347;
t397 = t155 * t288;
t103 = -t291 * t310 - t397;
t434 = t422 + t421;
t458 = t292 * t434;
t523 = pkin(1) * (t103 * t289 - t458);
t461 = t291 * t155;
t100 = t288 * t310 - t461;
t522 = pkin(2) * t100;
t510 = pkin(2) * t120;
t520 = pkin(9) * t100;
t509 = pkin(9) * t120;
t519 = pkin(9) * t123;
t403 = qJ(4) * t439;
t197 = pkin(3) * t244 - qJ(4) * t246;
t293 = qJD(1) ^ 2;
t414 = sin(qJ(1));
t415 = cos(qJ(1));
t315 = g(1) * t415 + g(2) * t414;
t253 = -t293 * pkin(1) + pkin(8) * t354 - t315;
t332 = -pkin(2) * t292 - pkin(9) * t289;
t256 = t332 * t361;
t314 = g(1) * t414 - g(2) * t415;
t371 = t285 * t293;
t302 = qJDD(1) * pkin(1) + pkin(8) * t371 + t314;
t300 = t286 * t302;
t373 = t285 * t289;
t296 = -g(3) * t373 + t289 * t300;
t334 = t335 ^ 2;
t139 = t280 * pkin(9) - t334 * pkin(2) + (t256 * t361 + t253) * t292 + t296;
t325 = t335 * qJD(1);
t319 = t289 * t325;
t320 = t292 * t325;
t409 = t286 * g(3);
t140 = -t258 * pkin(2) - t257 * pkin(9) - t409 + (pkin(2) * t319 - pkin(9) * t320 - t302) * t285;
t82 = t288 * t139 - t291 * t140;
t313 = t252 * pkin(3) - t271 * qJ(4) + qJDD(4) + t82;
t68 = t246 * t197 + t313;
t299 = -pkin(3) * t441 + t403 - t68;
t497 = pkin(2) * t474;
t517 = -t299 - t497;
t464 = t289 * t434;
t512 = pkin(1) * t100 + pkin(8) * (t103 * t292 + t464);
t496 = pkin(9) * t474;
t495 = pkin(9) * t475;
t471 = pkin(2) * t434;
t507 = pkin(9) * t103 - t471;
t506 = pkin(2) * t154 + t519;
t499 = 2 * qJD(4);
t338 = t289 * t253 - t292 * t300;
t138 = t285 * (g(3) * t292 + t256 * t360) - t280 * pkin(2) - t334 * pkin(9) + t338;
t485 = qJ(4) * t154;
t297 = t191 * pkin(3) + t138 + t485;
t295 = t246 * t499 - t297;
t58 = t295 - (t442 - t226) * pkin(3);
t437 = -t421 + t271;
t477 = t291 * t437 - t455;
t490 = t286 * t477 + (t289 * (-t437 * t288 - t454) - t292 * t155) * t285;
t484 = qJ(4) * t434;
t483 = qJ(5) * t155;
t83 = t291 * t139 + t288 * t140;
t346 = t271 * pkin(3) - t83;
t317 = t244 * t197 + t273 * t499 + t346;
t473 = qJ(4) * (-t176 - t252);
t476 = -pkin(3) * t440 - t317 + t473;
t57 = pkin(3) * t226 + t295 - t485;
t408 = pkin(5) + qJ(4);
t198 = pkin(5) * t246 - pkin(10) * t244;
t211 = pkin(4) * t273 - qJ(5) * t246;
t236 = t422 * pkin(4);
t386 = t191 * qJ(5);
t420 = -2 * qJD(4);
t305 = t386 - t236 + (t420 - t211) * t273 - t346;
t418 = 2 * qJD(5);
t356 = t418 - t197;
t378 = t252 * qJ(4);
t41 = -t378 - t252 * pkin(5) - t271 * pkin(10) + (t198 + t356) * t244 + t305;
t470 = t408 * t41;
t287 = sin(qJ(6));
t290 = cos(qJ(6));
t207 = t244 * t287 - t290 * t273;
t209 = t244 * t290 + t273 * t287;
t164 = t209 * t207;
t188 = qJDD(6) + t347;
t443 = -t164 + t188;
t467 = t287 * t443;
t462 = t290 * t443;
t316 = (t244 * t288 + t246 * t291) * t273;
t372 = t285 * t292;
t453 = t252 * t372 + t286 * t316;
t206 = t207 ^ 2;
t423 = t209 ^ 2;
t119 = -t206 - t423;
t379 = t422 * qJ(5);
t416 = -pkin(4) - pkin(10);
t433 = t246 * t211 + qJDD(5);
t38 = -t379 + t295 + t416 * t191 + (t244 * pkin(5) + (pkin(3) + pkin(10)) * t246) * t273 + t347 * pkin(5) + t433;
t431 = t252 * pkin(4) - t483;
t306 = t313 + t431;
t340 = -pkin(4) * t244 - t197;
t419 = -2 * qJD(5);
t333 = t419 - t340;
t40 = -t271 * pkin(5) + t252 * pkin(10) + (-t198 + t333) * t246 + t306;
t21 = t287 * t40 - t290 * t38;
t22 = t287 * t38 + t290 * t40;
t12 = t287 * t21 + t290 * t22;
t451 = -t119 * t408 + t12;
t405 = t290 * t41;
t237 = qJD(6) + t246;
t339 = t287 * t191 - t290 * t252;
t90 = (qJD(6) + t237) * t209 + t339;
t450 = -t408 * t90 - t405;
t118 = -t207 * qJD(6) + t290 * t191 + t287 * t252;
t382 = t237 * t207;
t323 = t118 - t382;
t406 = t287 * t41;
t449 = -t323 * t408 + t406;
t404 = qJ(4) * t310;
t282 = t285 ^ 2;
t448 = t282 * (-t286 * t293 + t325);
t375 = t273 * t291;
t330 = -t246 * t375 + t288 * t347;
t446 = t286 * t330;
t195 = g(3) * t372 + t338;
t196 = t292 * t253 + t296;
t438 = t289 * t195 + t292 * t196;
t430 = t191 * pkin(4) - t433;
t417 = pkin(3) + pkin(4);
t45 = t246 * t333 + t306;
t298 = t244 * t356 + t305;
t46 = t298 - t378;
t429 = qJ(4) * t46 - t417 * t45;
t427 = t155 * t417 - t404;
t376 = t273 * t288;
t329 = -t291 * t191 - t244 * t376;
t348 = t292 * t202;
t351 = t244 * t375;
t426 = t286 * t329 + (t289 * (t191 * t288 - t351) + t348) * t285;
t212 = t246 * t376;
t425 = (-t212 + t351) * t373 + t453;
t424 = -t440 * t417 + t473;
t234 = t237 ^ 2;
t412 = pkin(3) * t291;
t411 = pkin(8) * t285;
t402 = qJ(4) * t291;
t108 = t164 + t188;
t401 = t108 * t290;
t400 = t138 * t288;
t399 = t138 * t291;
t381 = t237 * t287;
t380 = t237 * t290;
t374 = t282 * t293;
t370 = t287 * t108;
t272 = t292 * t289 * t374;
t255 = t272 + t280;
t368 = t289 * t255;
t254 = -t272 + t280;
t365 = t292 * t254;
t362 = t291 * t347 + t212;
t355 = pkin(3) - t416;
t353 = t288 * t164;
t352 = t291 * t164;
t283 = t289 ^ 2;
t350 = t283 * t374;
t284 = t292 ^ 2;
t349 = t284 * t374;
t343 = -t234 - t423;
t341 = -qJ(4) * t288 - pkin(2);
t54 = t288 * t82 + t291 * t83;
t336 = t292 * t355;
t66 = -t317 - t378;
t331 = -pkin(3) * t68 + qJ(4) * t66;
t328 = -pkin(3) * t155 + t404;
t11 = -t290 * t21 + t287 * t22;
t327 = t288 * t83 - t291 * t82;
t326 = -pkin(1) + t332;
t318 = -t285 * t348 + t362 * t373 + t446;
t311 = (-qJD(6) + t237) * t209 - t339;
t304 = t273 * t211 + t244 * t419 + t236 + t317;
t69 = (-pkin(3) * t273 + t420) * t246 + t297;
t262 = t285 * t320;
t261 = t285 * t319;
t260 = (t283 - t284) * t374;
t259 = -t334 - t349;
t235 = -t350 - t334;
t231 = t246 * t418;
t227 = t285 * t302 + t409;
t222 = t258 - t261;
t221 = t258 + t261;
t220 = -t262 + t257;
t174 = t234 - t423;
t173 = t206 - t234;
t163 = -t206 + t423;
t150 = -t191 - t226;
t132 = -t234 - t206;
t117 = -qJD(6) * t209 - t339;
t110 = (t207 * t290 - t209 * t287) * t237;
t109 = (-t207 * t287 - t209 * t290) * t237;
t105 = qJ(4) * t442 + qJ(5) * t441;
t102 = t150 * t291 + t397;
t94 = -t118 - t382;
t88 = -t118 * t290 + t209 * t381;
t87 = t118 * t287 + t209 * t380;
t86 = t117 * t287 - t207 * t380;
t85 = -t117 * t290 - t207 * t381;
t84 = t110 * t291 + t188 * t288;
t80 = -t173 * t290 + t370;
t79 = t174 * t287 - t462;
t78 = t173 * t287 + t401;
t77 = t174 * t290 + t467;
t76 = -t287 * t343 - t401;
t75 = t290 * t343 - t370;
t74 = qJ(5) * t176 - t154 * t417;
t73 = t290 * t132 - t467;
t72 = t287 * t132 + t462;
t71 = t291 * t88 + t353;
t70 = t291 * t86 - t353;
t67 = t400 + t506;
t65 = -pkin(2) * t161 - t399 + t495;
t64 = t68 + t484;
t63 = -t287 * t94 + t290 * t311;
t62 = t287 * t323 + t290 * t90;
t61 = t287 * t311 + t290 * t94;
t60 = -t287 * t90 + t290 * t323;
t59 = pkin(3) * t434 + t66;
t56 = t288 * t311 + t291 * t80;
t55 = -t288 * t94 + t291 * t79;
t52 = t288 * t76 + t291 * t323;
t51 = t288 * t323 - t291 * t76;
t50 = t288 * t73 + t291 * t90;
t49 = t288 * t90 - t291 * t73;
t48 = t379 + t69 + t430;
t47 = t163 * t288 + t291 * t62;
t44 = t119 * t291 + t288 * t63;
t43 = t119 * t288 - t291 * t63;
t42 = -pkin(2) * t138 + pkin(9) * t54;
t39 = (-t440 - t422) * qJ(5) - t430 + t57;
t37 = t288 * t68 + t291 * t66;
t36 = t288 * t66 - t291 * t68;
t35 = t246 * t340 + t231 - t306 + t483 - t484;
t34 = pkin(9) * t102 + t471 + t54;
t33 = (t439 + t422) * qJ(5) + (t442 + t191) * pkin(4) - t433 - t58;
t32 = t378 - t417 * t434 + (t310 - t191) * qJ(5) + t304;
t31 = t291 * t58 + t341 * t442 + t495;
t30 = -t519 + t288 * t57 - (pkin(2) + t412) * t154;
t29 = -qJ(4) * t48 - qJ(5) * t45;
t28 = t288 * t45 + t291 * t46;
t27 = t288 * t46 - t291 * t45;
t26 = t288 * t64 + t291 * t59 - t507;
t25 = pkin(2) * t442 + t105 * t288 + t291 * t33 - t495;
t24 = t288 * t39 + t291 * t74 - t506;
t23 = -qJ(5) * t63 + t408 * t61;
t19 = -qJ(5) * t46 - t417 * t48;
t18 = -qJ(5) * t323 + t355 * t75 - t405;
t17 = -qJ(5) * t90 + t355 * t72 - t406;
t16 = pkin(9) * t37 + (t341 - t412) * t69;
t15 = t288 * t35 + t291 * t32 + t507;
t14 = -qJ(5) * t76 + t408 * t75 - t22;
t13 = -qJ(5) * t73 + t408 * t72 - t21;
t10 = t12 * t288 + t291 * t41;
t9 = -t12 * t291 + t288 * t41;
t8 = -pkin(2) * t48 + pkin(9) * t28 + t19 * t291 + t288 * t29;
t7 = -qJ(5) * t119 + t355 * t61 + t11;
t6 = pkin(2) * t75 + pkin(9) * t52 + t14 * t288 + t18 * t291;
t5 = pkin(2) * t72 + pkin(9) * t50 + t13 * t288 + t17 * t291;
t4 = pkin(2) * t61 + pkin(9) * t44 + t23 * t288 + t291 * t7;
t3 = -qJ(5) * t12 + t11 * t408;
t2 = -qJ(5) * t41 + t11 * t355;
t1 = pkin(2) * t11 + pkin(9) * t10 + t2 * t291 + t288 * t3;
t20 = [0, 0, 0, 0, 0, qJDD(1), t314, t315, 0, 0, (t285 * t257 + t292 * t448) * t289, t286 * t260 + (t289 * t222 + t292 * (t262 + t257)) * t285, t286 * t220 + (t368 + t292 * (t334 - t350)) * t285, (t285 * t258 - t289 * t448) * t292, t286 * t221 + (t289 * (-t334 + t349) + t365) * t285, t286 * t280, (-t195 + pkin(1) * (t255 * t292 + t259 * t289)) * t286 + (t292 * t227 + pkin(1) * t222 + pkin(8) * (t259 * t292 - t368)) * t285, -t227 * t373 - t286 * t196 + pkin(1) * (-t285 * (t335 * t344 + t257) + (t235 * t292 - t254 * t289) * t286) + (-t235 * t289 - t365) * t411, pkin(1) * ((-t220 * t292 + t221 * t289) * t286 - (-t283 - t284) * t282 * t371) + (t220 * t289 + t221 * t292) * t411 + t438 * t285, pkin(1) * (t285 * t227 + (-t195 * t292 + t196 * t289) * t286) + t438 * t411, t318, -t526, t490, t426, t513, t425, (t65 + pkin(1) * (-t161 * t292 + t494)) * t286 + (t289 * (t400 - t496) + t292 * (t82 - t497) - t498 + pkin(8) * (t161 * t289 + t493)) * t285, (t67 + t524) * t286 + (t289 * (t399 - t509) + t292 * (t83 - t510) + t514) * t285, (t34 + pkin(1) * (t102 * t289 + t458)) * t286 + (-t289 * t327 + pkin(8) * (t102 * t292 - t464) + t326 * (t150 * t288 - t461)) * t285, (t42 + pkin(1) * (-t138 * t292 + t289 * t54)) * t286 + (pkin(8) * (t138 * t289 + t292 * t54) + t326 * t327) * t285, t318, t490, t526, t425, -t513, t426, (t31 - t527) * t286 + (t289 * (-t288 * t58 - t402 * t442 - t496) + t292 * t517 - t525) * t285, (t26 - t523) * t286 + (t289 * (-t288 * t59 + t291 * t64 - t520) + t292 * (-t328 - t522) - t512) * t285, (t30 - t524) * t286 + (t289 * (pkin(3) * t398 + t291 * t57 + t509) + t292 * (-t476 + t510) - t514) * t285, (t16 + pkin(1) * (t289 * t37 - t292 * t69)) * t286 + (t289 * (-pkin(9) * t36 + (pkin(3) * t288 - t402) * t69) + t292 * (-pkin(2) * t36 - t331) - pkin(1) * t36 + pkin(8) * (t289 * t69 + t292 * t37)) * t285, t426, -t526, t513, t446 + (t289 * t362 - t348) * t285, t490, (t244 * t291 - t246 * t288) * t273 * t373 + t453, (t24 - t524) * t286 + (t289 * (-t288 * t74 + t291 * t39 + t509) + t292 * (t304 - t386 - t424 + t510) - t514) * t285, (t25 + t527) * t286 + (t289 * (t105 * t291 - t288 * t33 + t496) + t292 * (t231 + (-t441 - t202) * pkin(4) - t431 - t517) + t525) * t285, (t15 + t523) * t286 + (t289 * (-t288 * t32 + t291 * t35 + t520) + t292 * (-t427 + t522) + t512) * t285, (t8 + pkin(1) * (t28 * t289 - t292 * t48)) * t286 + (t289 * (-pkin(9) * t27 - t19 * t288 + t29 * t291) + t292 * (-pkin(2) * t27 - t429) - pkin(1) * t27 + pkin(8) * (t28 * t292 + t289 * t48)) * t285, t286 * t71 + (t289 * (-t288 * t88 + t352) + t292 * t87) * t285, t286 * t47 + (t289 * (t163 * t291 - t288 * t62) + t292 * t60) * t285, t286 * t55 + (t289 * (-t288 * t79 - t291 * t94) + t292 * t77) * t285, t286 * t70 + (t289 * (-t288 * t86 - t352) - t292 * t85) * t285, t286 * t56 + (t289 * (-t288 * t80 + t291 * t311) + t292 * t78) * t285, t286 * t84 + (t289 * (-t110 * t288 + t188 * t291) + t292 * t109) * t285, (t5 + pkin(1) * (t289 * t50 + t292 * t72)) * t286 + (t289 * (-pkin(9) * t49 + t13 * t291 - t17 * t288) + t292 * (-pkin(2) * t49 + t450) - pkin(1) * t49 + pkin(8) * (-t289 * t72 + t292 * t50) + t73 * t336) * t285, (t6 + pkin(1) * (t289 * t52 + t292 * t75)) * t286 + (t289 * (-pkin(9) * t51 + t14 * t291 - t18 * t288) + t292 * (-pkin(2) * t51 + t449) - pkin(1) * t51 + pkin(8) * (-t289 * t75 + t292 * t52) + t76 * t336) * t285, (t4 + pkin(1) * (t289 * t44 + t292 * t61)) * t286 + (t289 * (-pkin(9) * t43 + t23 * t291 - t288 * t7) + t292 * (-pkin(2) * t43 + t451) - pkin(1) * t43 + pkin(8) * (-t289 * t61 + t292 * t44) + t63 * t336) * t285, (t1 + pkin(1) * (t10 * t289 + t11 * t292)) * t286 + (t289 * (-pkin(9) * t9 - t2 * t288 + t291 * t3) + t292 * (-pkin(2) * t9 - t470) - pkin(1) * t9 + pkin(8) * (t10 * t292 - t11 * t289) + t12 * t336) * t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, t260, t220, t272, t221, t280, -t195, -t196, 0, 0, t330, -t99, t477, t329, -t131, t316, t65, t67, t34, t42, t330, t477, t99, t316, t131, t329, t31, t26, t30, t16, t329, -t99, -t131, t330, t477, t316, t24, t25, t15, t8, t71, t47, t55, t70, t56, t84, t5, t6, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t435, t155, -t202, t310, -t252, -t82, -t83, 0, 0, t202, t155, t435, -t252, -t310, -t202, t299, t328, t476, t331, -t202, -t435, t310, t202, t155, -t252, t298 + t424, t417 * t441 - t403 + t45, t427, t429, -t87, -t60, -t77, t85, -t78, -t109, -t355 * t73 - t450, -t355 * t76 - t449, -t355 * t63 - t451, -t12 * t355 + t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t441, t155, t440, t68, 0, 0, 0, 0, 0, 0, t440, -t441, -t155, t45, 0, 0, 0, 0, 0, 0, t73, t76, t63, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t442, -t434, -t48, 0, 0, 0, 0, 0, 0, t72, t75, t61, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t163, -t94, -t164, t311, t188, -t21, -t22, 0, 0;];
tauJ_reg  = t20;
