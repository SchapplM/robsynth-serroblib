% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 18:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRRPP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:28:09
% EndTime: 2019-05-07 18:28:38
% DurationCPUTime: 11.51s
% Computational Cost: add. (28221->506), mult. (56710->573), div. (0->0), fcn. (40330->8), ass. (0->312)
t312 = sin(qJ(3));
t316 = cos(qJ(3));
t313 = sin(qJ(2));
t375 = qJD(1) * t313;
t281 = -t316 * qJD(2) + t312 * t375;
t302 = t313 * qJDD(1);
t317 = cos(qJ(2));
t369 = qJD(1) * qJD(2);
t365 = t317 * t369;
t345 = t302 + t365;
t330 = -t312 * qJDD(2) - t316 * t345;
t251 = -t281 * qJD(3) - t330;
t282 = qJD(2) * t312 + t316 * t375;
t311 = sin(qJ(4));
t315 = cos(qJ(4));
t257 = t315 * t281 + t282 * t311;
t329 = t316 * qJDD(2) - t312 * t345;
t325 = -t282 * qJD(3) + t329;
t321 = -t257 * qJD(4) + t315 * t251 + t311 * t325;
t299 = qJD(1) * t317 - qJD(3);
t295 = -qJD(4) + t299;
t385 = t295 * t257;
t135 = -t385 - t321;
t259 = -t281 * t311 + t282 * t315;
t256 = t259 ^ 2;
t427 = t295 ^ 2;
t197 = -t427 - t256;
t301 = t313 * t369;
t368 = t317 * qJDD(1);
t285 = -t301 + t368;
t280 = -qJDD(3) + t285;
t277 = -qJDD(4) + t280;
t389 = t259 * t257;
t442 = -t389 + t277;
t458 = t442 * t311;
t136 = t197 * t315 + t458;
t457 = t442 * t315;
t496 = t197 * t311 - t457;
t518 = t136 * t312 + t316 * t496;
t85 = t136 * t316 - t312 * t496;
t536 = pkin(7) * (t135 * t313 + t317 * t518) + pkin(1) * t85;
t530 = pkin(2) * t85;
t528 = pkin(8) * t85;
t429 = t257 ^ 2;
t444 = t256 - t429;
t361 = t311 * t251 - t315 * t325;
t164 = qJD(4) * t259 + t361;
t388 = t259 * t295;
t448 = t164 - t388;
t473 = t135 * t315 + t448 * t311;
t406 = t135 * t311;
t456 = t448 * t315;
t82 = -t456 + t406;
t535 = t317 * t444 - t313 * (t312 * t473 + t316 * t82);
t443 = -t385 + t321;
t116 = t315 * t443;
t123 = (qJD(4) + t295) * t259 + t361;
t72 = t123 * t311 + t116;
t407 = t443 * t311;
t78 = t123 * t315 - t407;
t45 = t312 * t78 + t316 * t72;
t439 = -t429 - t256;
t464 = t313 * t439;
t50 = t312 * t72 - t316 * t78;
t534 = pkin(7) * (t317 * t50 + t464) + pkin(1) * t45;
t531 = pkin(2) * t45;
t529 = pkin(8) * t45;
t520 = -t312 * t82 + t316 * t473;
t470 = pkin(2) * t439;
t526 = pkin(8) * t50 - t470;
t525 = -pkin(2) * t135 + pkin(8) * t518;
t438 = t429 - t427;
t146 = -t311 * t438 + t457;
t447 = t164 + t388;
t477 = t315 * t438 + t458;
t488 = t313 * (t146 * t312 + t316 * t477) + t317 * t447;
t231 = t256 - t427;
t441 = t389 + t277;
t459 = t441 * t311;
t143 = t231 * t315 + t459;
t461 = t315 * t441;
t149 = -t231 * t311 + t461;
t489 = t313 * (t143 * t312 - t149 * t316) - t317 * t443;
t511 = pkin(3) * t136;
t510 = pkin(9) * t136;
t517 = t143 * t316 + t149 * t312;
t515 = t146 * t316 - t312 * t477;
t422 = pkin(3) * t72;
t513 = pkin(9) * t72;
t440 = -t427 - t429;
t474 = t311 * t440 - t461;
t476 = t315 * t440 + t459;
t494 = t312 * t476 + t316 * t474;
t512 = pkin(1) * t494;
t501 = pkin(2) * t494;
t499 = pkin(8) * t494;
t509 = pkin(9) * t496;
t485 = pkin(9) * t476;
t504 = -pkin(3) * t448 + t485;
t469 = pkin(3) * t439;
t503 = -pkin(9) * t78 - t469;
t493 = -t312 * t474 + t316 * t476;
t495 = -pkin(2) * t448 + pkin(8) * t493;
t502 = pkin(3) * t135 - t509;
t500 = pkin(3) * t474;
t484 = pkin(9) * t474;
t409 = t448 * t313;
t490 = pkin(7) * (t317 * t493 + t409) - t512;
t483 = qJ(5) * t135;
t482 = qJ(5) * t439;
t481 = qJ(5) * t442;
t319 = qJD(1) ^ 2;
t314 = sin(qJ(1));
t318 = cos(qJ(1));
t351 = g(1) * t318 + g(2) * t314;
t412 = qJDD(1) * pkin(7);
t276 = -pkin(1) * t319 - t351 + t412;
t352 = -pkin(2) * t317 - pkin(8) * t313;
t359 = t319 * t352 + t276;
t417 = t317 * g(3);
t425 = qJD(2) ^ 2;
t224 = -qJDD(2) * pkin(2) - t425 * pkin(8) + t313 * t359 + t417;
t265 = -pkin(3) * t299 - pkin(9) * t282;
t428 = t281 ^ 2;
t142 = -t325 * pkin(3) - t428 * pkin(9) + t282 * t265 + t224;
t480 = t164 * pkin(4) + t142 + t483;
t424 = pkin(4) + pkin(5);
t467 = qJ(6) * t443;
t453 = -pkin(4) * t441 + qJ(5) * t440;
t386 = t282 * t281;
t333 = -t280 - t386;
t452 = t312 * t333;
t336 = (t257 * t311 + t259 * t315) * t295;
t451 = t312 * t336;
t450 = t316 * t333;
t449 = t316 * t336;
t374 = qJD(5) * t295;
t287 = -0.2e1 * t374;
t373 = qJD(6) * t257;
t446 = 0.2e1 * t373 + t287;
t288 = 0.2e1 * t374;
t445 = -0.2e1 * t373 + t288;
t268 = t281 * t299;
t215 = -t268 + t251;
t226 = pkin(5) * t295 - qJ(6) * t259;
t437 = t259 * t226 + qJDD(6);
t434 = -t164 * pkin(5) + t437;
t211 = (qJD(3) + t299) * t282 - t329;
t384 = t295 * t311;
t227 = t259 * t384;
t383 = t295 * t315;
t367 = t257 * t383;
t348 = -t227 + t367;
t433 = t312 * t348 + t449;
t339 = t164 * t311 - t367;
t349 = -t315 * t164 - t257 * t384;
t432 = t312 * t339 + t316 * t349;
t350 = -t259 * t383 + t311 * t321;
t376 = t315 * t321 + t227;
t71 = t312 * t376 + t316 * t350;
t269 = t317 * t277;
t431 = t313 * (t316 * t348 - t451) + t269;
t202 = t317 * t389;
t430 = t313 * (-t312 * t349 + t316 * t339) + t202;
t413 = t313 * (-t312 * t350 + t316 * t376) - t202;
t279 = t282 ^ 2;
t297 = t299 ^ 2;
t426 = 0.2e1 * t259;
t362 = t314 * g(1) - t318 * g(2);
t275 = qJDD(1) * pkin(1) + t319 * pkin(7) + t362;
t284 = t302 + 0.2e1 * t365;
t346 = -t285 + t301;
t210 = pkin(2) * t346 - pkin(8) * t284 - t275;
t418 = t313 * g(3);
t225 = -pkin(2) * t425 + qJDD(2) * pkin(8) + t317 * t359 - t418;
t166 = -t316 * t210 + t312 * t225;
t97 = t333 * pkin(3) - pkin(9) * t215 - t166;
t167 = t312 * t210 + t316 * t225;
t99 = -pkin(3) * t428 + pkin(9) * t325 + t299 * t265 + t167;
t61 = t311 * t99 - t315 * t97;
t62 = t311 * t97 + t315 * t99;
t36 = t311 * t62 - t315 * t61;
t423 = pkin(3) * t36;
t420 = pkin(4) * t315;
t201 = pkin(4) * t257 - qJ(5) * t259;
t347 = -pkin(4) * t427 - t277 * qJ(5) - t257 * t201 + t62;
t54 = t287 + t347;
t337 = t277 * pkin(4) - qJ(5) * t427 + qJDD(5) + t61;
t390 = t259 * t201;
t56 = t337 + t390;
t416 = -pkin(4) * t56 + qJ(5) * t54;
t415 = t312 * t36;
t414 = t316 * t36;
t404 = t142 * t311;
t403 = t142 * t315;
t395 = t224 * t312;
t394 = t224 * t316;
t242 = t280 - t386;
t392 = t242 * t312;
t391 = t242 * t316;
t382 = t299 * t312;
t381 = t299 * t316;
t298 = t317 * t319 * t313;
t380 = t313 * (qJDD(2) + t298);
t378 = t317 * (-t298 + qJDD(2));
t377 = -pkin(4) * t443 - qJ(5) * t447;
t372 = qJD(3) - t299;
t366 = t317 * t386;
t364 = -qJ(5) * t311 - pkin(3);
t363 = -pkin(5) * t257 - t201;
t37 = t311 * t61 + t315 * t62;
t94 = t166 * t312 + t316 * t167;
t263 = t313 * t276 + t417;
t264 = t276 * t317 - t418;
t360 = t313 * t263 + t317 * t264;
t29 = t311 * t54 - t315 * t56;
t358 = pkin(3) * t29 + t416;
t331 = t277 * pkin(5) + t337 - t467;
t341 = (-0.2e1 * qJD(6) - t363) * t259;
t39 = t341 + t331;
t334 = pkin(5) * t429 - t164 * qJ(6) + t295 * t226 - t347;
t41 = -t334 + t446;
t357 = qJ(5) * t41 - t39 * t424;
t356 = -t62 + t511;
t75 = -t311 * t447 - t116;
t355 = pkin(3) * t75 + t377;
t354 = qJ(5) * t123 + t424 * t443;
t344 = t166 * t316 - t167 * t312;
t343 = -pkin(1) + t352;
t342 = -t61 + t500;
t18 = t311 * t41 - t315 * t39;
t340 = pkin(3) * t18 + t357;
t338 = t354 + t422;
t335 = -pkin(4) * t197 + t347 - t481;
t332 = t335 - t511;
t328 = t453 - t56;
t327 = -t331 + t453;
t326 = -t197 * t424 - t334 - t481;
t324 = t328 + t500;
t323 = t326 - t511;
t322 = qJD(5) * t426 - t480;
t249 = qJD(6) * t426;
t320 = -t390 + t249 + (-t441 - t389) * pkin(5) + t327;
t57 = (-pkin(4) * t295 - 0.2e1 * qJD(5)) * t259 + t480;
t52 = (-t448 + t388) * pkin(4) + t322;
t51 = pkin(4) * t388 + t322 - t483;
t308 = t317 ^ 2;
t307 = t313 ^ 2;
t305 = t308 * t319;
t303 = t307 * t319;
t286 = -0.2e1 * t301 + t368;
t267 = -t279 + t297;
t266 = -t297 + t428;
t261 = -t279 + t428;
t260 = -t279 - t297;
t252 = -t297 - t428;
t241 = t279 + t428;
t216 = t281 * t372 + t330;
t214 = t268 + t251;
t212 = -t282 * t372 + t329;
t199 = -t260 * t312 + t391;
t198 = t260 * t316 + t392;
t191 = t252 * t316 - t452;
t190 = t252 * t312 + t450;
t175 = (t257 * t315 - t259 * t311) * t295;
t156 = -t211 * t316 + t215 * t312;
t92 = -qJ(5) * t448 - qJ(6) * t441;
t91 = t403 - t510;
t84 = t404 - t484;
t81 = -t315 * t447 + t407;
t63 = qJ(6) * t442 - t135 * t424;
t59 = t404 + t502;
t58 = -t403 + t504;
t49 = -t312 * t75 + t316 * t81;
t46 = t312 * t81 + t316 * t75;
t44 = t56 - t482;
t43 = -pkin(4) * t439 + t54;
t42 = qJ(6) * t429 - t434 + t57;
t35 = t51 + (-t197 - t429) * qJ(6) + t434;
t34 = -qJ(5) * t456 - t311 * t52 - t484;
t33 = -pkin(3) * t142 + pkin(9) * t37;
t32 = pkin(4) * t406 + t315 * t51 + t510;
t31 = t259 * t363 + t249 - t331 + t467 + t482;
t30 = t311 * t56 + t315 * t54;
t28 = t315 * t52 + t364 * t448 + t485;
t27 = (-t440 - t429) * qJ(6) + t52 + (-t448 - t164) * pkin(5) + t437;
t26 = t509 + t311 * t51 - (pkin(3) + t420) * t135;
t25 = -t36 + t513;
t24 = -qJ(6) * t123 + t424 * t439 + t334 + t445;
t23 = t37 + t503;
t22 = -t311 * t63 + t315 * t35 + t510;
t21 = -t27 * t311 + t315 * t92 - t484;
t20 = -qJ(5) * t42 - qJ(6) * t39;
t19 = t311 * t39 + t315 * t41;
t17 = -pkin(9) * t75 - t311 * t43 + t315 * t44;
t16 = t311 * t35 + t315 * t63 - t502;
t15 = t316 * t37 - t415;
t14 = t312 * t37 + t414;
t13 = t27 * t315 + t311 * t92 + t504;
t12 = pkin(9) * t81 + t311 * t44 + t315 * t43 - t469;
t11 = -pkin(9) * t29 + (pkin(4) * t311 - qJ(5) * t315) * t57;
t10 = -t29 * t312 + t30 * t316;
t9 = t29 * t316 + t30 * t312;
t8 = -qJ(6) * t41 - t42 * t424;
t7 = -t24 * t311 + t31 * t315 - t513;
t6 = pkin(9) * t30 + (t364 - t420) * t57;
t5 = t24 * t315 + t31 * t311 - t503;
t4 = -t18 * t312 + t19 * t316;
t3 = t18 * t316 + t19 * t312;
t2 = -pkin(9) * t18 + t20 * t315 - t311 * t8;
t1 = -pkin(3) * t42 + pkin(9) * t19 + t20 * t311 + t315 * t8;
t38 = [0, 0, 0, 0, 0, qJDD(1), t362, t351, 0, 0, t284 * t313, t284 * t317 + t286 * t313, t380 + t317 * (-t303 + t425), -t346 * t317, t313 * (t305 - t425) + t378, 0, t317 * t275 + pkin(1) * t286 + pkin(7) * (t317 * (-t305 - t425) - t380), -t313 * t275 - pkin(1) * t284 + pkin(7) * (-t378 - t313 * (-t303 - t425)), pkin(1) * (t303 + t305) + (t307 + t308) * t412 + t360, pkin(1) * t275 + pkin(7) * t360, t313 * (t251 * t316 + t282 * t382) - t366, t313 * (t212 * t316 - t214 * t312) + t317 * t261, t313 * (-t267 * t312 + t450) - t317 * t215, t313 * (-t281 * t381 - t312 * t325) + t366, t313 * (t266 * t316 + t392) + t317 * t211, t317 * t280 + t313 * (t281 * t316 - t282 * t312) * t299, t313 * (-pkin(8) * t190 + t395) + t317 * (-pkin(2) * t190 + t166) - pkin(1) * t190 + pkin(7) * (t191 * t317 - t212 * t313), t313 * (-pkin(8) * t198 + t394) + t317 * (-pkin(2) * t198 + t167) - pkin(1) * t198 + pkin(7) * (t199 * t317 - t216 * t313), t313 * t344 + pkin(7) * (t156 * t317 - t241 * t313) + t343 * (-t211 * t312 - t215 * t316), pkin(7) * (t224 * t313 + t317 * t94) - t343 * t344, t413, -t535, t489, t430, t488, t431, t313 * (-t312 * t58 + t316 * t84 - t499) + t317 * (-t342 - t501) + t490, t313 * (-t312 * t59 + t316 * t91 - t528) + t317 * (-t356 - t530) - t536, t313 * (-t23 * t312 + t25 * t316 + t529) + t317 * (t422 + t531) + t534, t313 * (-pkin(8) * t14 - pkin(9) * t414 - t312 * t33) + t317 * (-pkin(2) * t14 - t423) - pkin(1) * t14 + pkin(7) * (t142 * t313 + t15 * t317), t413, t489, t535, t431, -t488, t430, t313 * (-t28 * t312 + t316 * t34 - t499) + t317 * (-t324 - t501) + t490, t313 * (-pkin(8) * t46 - t12 * t312 + t17 * t316) + t317 * (-pkin(2) * t46 - t355) - pkin(1) * t46 + pkin(7) * (t317 * t49 + t464), t313 * (-t26 * t312 + t316 * t32 + t528) + t317 * (t288 - t332 + t530) + t536, t313 * (-pkin(8) * t9 + t11 * t316 - t312 * t6) + t317 * (-pkin(2) * t9 - t358) - pkin(1) * t9 + pkin(7) * (t10 * t317 + t313 * t57), t413, t535, -t489, t430, t488, t313 * (t175 * t316 - t451) + t269, t313 * (-t13 * t312 + t21 * t316 - t499) - t512 + pkin(7) * t409 + (pkin(5) * t441 + pkin(7) * t493 - t327 + t341 - t500 - t501) * t317, t313 * (-t16 * t312 + t22 * t316 + t528) + t317 * (-t323 + t445 + t530) + t536, t313 * (-t312 * t5 + t316 * t7 - t529) + t317 * (-t338 - t531) - t534, t313 * (-pkin(8) * t3 - t1 * t312 + t2 * t316) + t317 * (-pkin(2) * t3 - t340) - pkin(1) * t3 + pkin(7) * (t313 * t42 + t317 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t298, t303 - t305, t302, t298, t368, qJDD(2), -t263, -t264, 0, 0, t251 * t312 - t282 * t381, t212 * t312 + t214 * t316, t267 * t316 + t452, -t281 * t382 + t316 * t325, t266 * t312 - t391, (t281 * t312 + t282 * t316) * t299, pkin(2) * t212 + pkin(8) * t191 - t394, pkin(2) * t216 + pkin(8) * t199 + t395, pkin(2) * t241 + pkin(8) * t156 + t94, -pkin(2) * t224 + pkin(8) * t94, t71, -t520, -t517, t432, -t515, t433, t312 * t84 + t316 * t58 + t495, t312 * t91 + t316 * t59 - t525, t23 * t316 + t25 * t312 + t526, -pkin(2) * t142 + pkin(8) * t15 - pkin(9) * t415 + t316 * t33, t71, -t517, t520, t433, t515, t432, t28 * t316 + t312 * t34 + t495, pkin(8) * t49 + t12 * t316 + t17 * t312 - t470, t26 * t316 + t312 * t32 + t525, -pkin(2) * t57 + pkin(8) * t10 + t11 * t312 + t316 * t6, t71, t520, t517, t432, -t515, t175 * t312 + t449, t13 * t316 + t21 * t312 + t495, t16 * t316 + t22 * t312 + t525, t312 * t7 + t316 * t5 - t526, -pkin(2) * t42 + pkin(8) * t4 + t1 * t316 + t2 * t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, -t261, t215, -t386, -t211, -t280, -t166, -t167, 0, 0, t389, t444, t443, -t389, -t447, -t277, t342, t356, -t422, t423, t389, t443, -t444, -t277, t447, -t389, t324, t355, t287 + t332, t358, t389, -t444, -t443, -t389, -t447, -t277, t320 + t500, t323 + t446, t338, t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t389, t444, t443, -t389, -t447, -t277, -t61, -t62, 0, 0, t389, t443, -t444, -t277, t447, -t389, t328, t377, t287 + t335, t416, t389, -t444, -t443, -t389, -t447, -t277, t320, t326 + t446, t354, t357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t441, t443, t197, t56, 0, 0, 0, 0, 0, 0, t441, t197, -t443, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448, -t135, t439, -t42;];
tauJ_reg  = t38;
