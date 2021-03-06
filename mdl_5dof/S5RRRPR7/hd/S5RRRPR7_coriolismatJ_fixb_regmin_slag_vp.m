% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:19
% EndTime: 2019-12-31 21:17:33
% DurationCPUTime: 6.19s
% Computational Cost: add. (6684->443), mult. (14026->603), div. (0->0), fcn. (15504->8), ass. (0->367)
t409 = qJD(2) + qJD(3);
t314 = cos(pkin(9));
t316 = sin(qJ(2));
t513 = -pkin(7) - pkin(6);
t290 = t513 * t316;
t501 = sin(qJ(3));
t277 = t501 * t290;
t318 = cos(qJ(2));
t291 = t513 * t318;
t502 = cos(qJ(3));
t400 = t502 * t291;
t368 = -t400 + t277;
t528 = t368 * t314;
t512 = -t528 / 0.2e1;
t274 = t501 * t316 - t502 * t318;
t304 = t501 * t318;
t305 = t502 * t316;
t276 = -t304 - t305;
t525 = t409 * t276;
t544 = t274 * t525;
t313 = sin(pkin(9));
t367 = -t502 * t290 - t501 * t291;
t451 = t313 * t367;
t201 = t314 * t367;
t268 = -t400 / 0.2e1;
t543 = t268 + t277 / 0.2e1;
t317 = cos(qJ(5));
t446 = t317 * t313;
t315 = sin(qJ(5));
t447 = t315 * t314;
t273 = t446 + t447;
t526 = t409 * t273;
t448 = t315 * t313;
t247 = t276 * t448;
t445 = t317 * t314;
t81 = -t247 / 0.2e1 + (-t448 / 0.2e1 + t445) * t276;
t542 = -t81 * qJD(1) + t526;
t308 = -t318 * pkin(2) - pkin(1);
t498 = t274 * pkin(3);
t366 = t276 * qJ(4) + t498;
t187 = t308 + t366;
t529 = t368 * t313;
t103 = t314 * t187 - t529;
t104 = t313 * t187 + t528;
t535 = t104 * t314 / 0.2e1 - t103 * t313 / 0.2e1;
t523 = -t535 + t543;
t541 = qJD(1) * t523;
t540 = qJD(4) * t523;
t522 = t535 + t543;
t539 = t522 * qJD(4);
t538 = (-t103 - t529) * t276;
t537 = (t104 - t528) * t276;
t271 = -t445 + t448;
t168 = -t276 * t445 + t247;
t340 = t273 * t276;
t509 = t273 / 0.2e1;
t510 = -t271 / 0.2e1;
t59 = t168 * t510 + t340 * t509;
t53 = -t59 * qJD(1) + t271 * t526;
t536 = t409 * t367;
t499 = t368 * pkin(3);
t345 = t447 / 0.2e1 + t446 / 0.2e1;
t111 = (t509 - t345) * t274;
t534 = t111 * t409;
t112 = (t509 + t345) * t274;
t533 = t112 * t409;
t428 = qJD(1) * t274;
t398 = t340 * t428;
t185 = t313 * t274;
t356 = -pkin(4) * t185 + t368;
t532 = t356 * t271;
t531 = t356 * t273;
t530 = t367 * t368;
t497 = t314 * pkin(4);
t303 = -pkin(3) - t497;
t503 = t303 / 0.2e1;
t407 = t502 * pkin(2);
t307 = -t407 - pkin(3);
t286 = t307 - t497;
t504 = t286 / 0.2e1;
t389 = t503 + t504;
t527 = t389 * t273;
t311 = t313 ^ 2;
t312 = t314 ^ 2;
t293 = t311 + t312;
t524 = t409 * t293;
t406 = t501 * pkin(2);
t302 = t406 + qJ(4);
t458 = t276 * t302;
t459 = t274 * t307;
t506 = -t276 / 0.2e1;
t508 = -t274 / 0.2e1;
t521 = t498 / 0.2e1 + t458 / 0.2e1 - t459 / 0.2e1 + (t501 * t506 + t502 * t508) * pkin(2);
t418 = t340 * qJD(5);
t429 = qJD(1) * t168;
t519 = t340 * t429 + t409 * t59;
t55 = -t273 * t168 - t271 * t340;
t57 = -t168 ^ 2 + t340 ^ 2;
t518 = t57 * qJD(1) + t409 * t55;
t156 = t271 ^ 2 - t273 ^ 2;
t30 = t55 * qJD(1) + t409 * t156;
t457 = t276 * t313;
t150 = -pkin(4) * t457 + t367;
t167 = t271 * t274;
t327 = t274 * pkin(4) - t529 + (pkin(8) * t276 + t187) * t314;
t76 = pkin(8) * t457 + t104;
t39 = t315 * t327 + t317 * t76;
t517 = t150 * t167 + t356 * t168 + t39 * t276;
t516 = -t340 * qJD(1) + t271 * t409;
t163 = t273 * t274;
t38 = t315 * t76 - t317 * t327;
t515 = -t150 * t163 + t38 * t276 - t356 * t340;
t514 = t276 ^ 2;
t511 = t368 / 0.2e1;
t507 = t274 / 0.2e1;
t500 = pkin(2) * t274;
t310 = t314 * pkin(8);
t496 = t316 * pkin(2);
t200 = -t276 * pkin(3) + t274 * qJ(4);
t188 = t200 + t496;
t117 = t314 * t188 + t451;
t369 = -t276 * pkin(4) + t274 * t310;
t68 = t117 + t369;
t488 = t317 * t68;
t118 = t313 * t188 - t201;
t408 = pkin(8) * t185;
t79 = t408 + t118;
t490 = t315 * t79;
t1 = (t488 - t490) * t274 + t515;
t495 = t1 * qJD(1);
t486 = t317 * t79;
t492 = t315 * t68;
t2 = -(t486 + t492) * t274 + t517;
t494 = t2 * qJD(1);
t119 = t314 * t200 + t451;
t69 = t119 + t369;
t487 = t317 * t69;
t120 = t313 * t200 - t201;
t80 = t408 + t120;
t489 = t315 * t80;
t3 = (t487 - t489) * t274 + t515;
t493 = t3 * qJD(1);
t491 = t315 * t69;
t485 = t317 * t80;
t4 = -(t485 + t491) * t274 + t517;
t484 = t4 * qJD(1);
t481 = t117 * t313;
t480 = t118 * t314;
t479 = t119 * t313;
t478 = t120 * t314;
t475 = t150 * t271;
t474 = t150 * t273;
t362 = t103 * t314 + t104 * t313;
t341 = t362 * t274;
t16 = (-t117 * t314 - t118 * t313) * t276 - t341;
t471 = t16 * qJD(1);
t262 = t314 * t302 + t310;
t387 = (-pkin(8) - t302) * t313;
t172 = t315 * t262 - t317 * t387;
t470 = t172 * t276;
t173 = t317 * t262 + t315 * t387;
t469 = t173 * t276;
t287 = (-pkin(8) - qJ(4)) * t313;
t449 = t314 * qJ(4);
t288 = t310 + t449;
t209 = -t317 * t287 + t315 * t288;
t468 = t209 * t276;
t21 = (-t119 * t314 - t120 * t313) * t276 - t341;
t467 = t21 * qJD(1);
t210 = t315 * t287 + t317 * t288;
t466 = t210 * t276;
t22 = t103 * t117 + t104 * t118 + t530;
t463 = t22 * qJD(1);
t24 = t103 * t119 + t104 * t120 + t530;
t462 = t24 * qJD(1);
t26 = t150 * t340 + t38 * t274;
t461 = t26 * qJD(1);
t27 = t150 * t168 - t39 * t274;
t460 = t27 * qJD(1);
t456 = t286 * t163;
t455 = t286 * t167;
t454 = t303 * t163;
t453 = t303 * t167;
t452 = t313 * qJ(4);
t32 = t537 + (-t118 - t201) * t274;
t444 = t32 * qJD(1);
t33 = t538 + (t117 - t451) * t274;
t443 = t33 * qJD(1);
t34 = t537 + (-t120 - t201) * t274;
t442 = t34 * qJD(1);
t35 = t538 + (t119 - t451) * t274;
t441 = t35 * qJD(1);
t375 = -t406 / 0.2e1;
t322 = (t302 / 0.2e1 + t375 - qJ(4) / 0.2e1) * t276 + (-t407 / 0.2e1 - t307 / 0.2e1 - pkin(3) / 0.2e1) * t274;
t41 = t512 + t528 / 0.2e1 + t322 * t313;
t440 = t41 * qJD(1);
t51 = t362 * t276;
t439 = t51 * qJD(1);
t56 = t168 * t163 + t167 * t340;
t437 = t56 * qJD(1);
t61 = t163 * t274 - t276 * t340;
t434 = t61 * qJD(1);
t62 = t167 * t274 - t168 * t276;
t433 = t62 * qJD(1);
t339 = t293 * t502;
t261 = t339 * pkin(2);
t289 = t293 * qJD(4);
t430 = t261 * qJD(3) + t289;
t427 = qJD(1) * t276;
t426 = qJD(1) * t308;
t425 = qJD(1) * t318;
t424 = qJD(2) * t286;
t423 = qJD(3) * t303;
t422 = qJD(3) * t308;
t421 = qJD(4) * t274;
t102 = t111 * qJD(1);
t95 = t112 * qJD(1);
t390 = -t445 / 0.2e1;
t344 = t390 + t448 / 0.2e1;
t113 = (t271 / 0.2e1 + t344) * t274;
t96 = t113 * qJD(1);
t114 = (t510 + t344) * t274;
t99 = t114 * qJD(1);
t121 = t293 * t514;
t420 = t121 * qJD(1);
t159 = t274 ^ 2 - t514;
t419 = t159 * qJD(1);
t181 = t274 * t496 - t308 * t276;
t417 = t181 * qJD(1);
t182 = -t308 * t274 - t276 * t496;
t416 = t182 * qJD(1);
t269 = t400 / 0.2e1;
t213 = t269 + t268;
t415 = t213 * qJD(1);
t264 = t305 / 0.2e1 + t304 / 0.2e1;
t414 = t264 * qJD(1);
t259 = t271 * qJD(5);
t413 = t273 * qJD(5);
t294 = -t316 ^ 2 + t318 ^ 2;
t412 = t294 * qJD(1);
t411 = t316 * qJD(2);
t410 = t318 * qJD(2);
t405 = pkin(1) * t316 * qJD(1);
t404 = pkin(1) * t425;
t402 = t315 * t502;
t401 = t317 * t502;
t162 = t271 * t276;
t399 = t162 * t428;
t211 = t274 * t427;
t396 = t274 * t426;
t395 = t276 * t426;
t253 = t314 * t421;
t394 = t316 * t425;
t393 = t475 / 0.2e1;
t392 = -t474 / 0.2e1;
t391 = t449 / 0.2e1;
t386 = t502 * qJD(2);
t385 = t502 * qJD(3);
t384 = t501 * qJD(2);
t383 = t501 * qJD(3);
t378 = -qJD(5) - t428;
t377 = pkin(2) * t383;
t376 = pkin(2) * t384;
t374 = t406 / 0.2e1;
t373 = t313 * t402;
t372 = t314 * t402;
t371 = t314 * t211;
t370 = -t401 / 0.2e1;
t174 = (t339 * t302 + t501 * t307) * pkin(2);
t350 = t478 / 0.2e1 - t479 / 0.2e1;
t319 = t350 * t302 + (t367 * t501 / 0.2e1 + t535 * t502) * pkin(2) + t307 * t511;
t351 = -t480 / 0.2e1 + t481 / 0.2e1;
t6 = t499 / 0.2e1 + t351 * qJ(4) + t319;
t363 = t6 * qJD(1) + t174 * qJD(2);
t361 = t480 - t481;
t360 = t478 - t479;
t359 = t458 - t459;
t232 = t293 * t302;
t358 = t232 * qJD(2) - t541;
t28 = (t120 / 0.2e1 - t118 / 0.2e1) * t314 + (-t119 / 0.2e1 + t117 / 0.2e1) * t313;
t357 = -t28 * qJD(1) - t261 * qJD(2);
t355 = -t492 / 0.2e1 - t486 / 0.2e1;
t354 = -t491 / 0.2e1 - t485 / 0.2e1;
t353 = -t490 / 0.2e1 + t488 / 0.2e1;
t352 = -t489 / 0.2e1 + t487 / 0.2e1;
t349 = t172 * t507 + t340 * t504;
t348 = t168 * t504 + t173 * t508;
t347 = t209 * t507 + t340 * t503;
t346 = t168 * t503 + t210 * t508;
t12 = t392 - t348 + t353;
t343 = t12 * qJD(1) - t273 * t424;
t13 = t393 - t349 + t355;
t342 = t13 * qJD(1) + t271 * t424;
t145 = -t264 * qJD(5) + t211;
t338 = t344 * t274;
t321 = t532 / 0.2e1 + t470 / 0.2e1 - t456 / 0.2e1 + (-t313 * t401 - t372) * t500 / 0.2e1 + t340 * t375;
t329 = -t532 / 0.2e1 - t468 / 0.2e1 + t454 / 0.2e1;
t9 = t321 + t329;
t336 = -t9 * qJD(1) - t271 * t376;
t326 = (qJ(4) + t302) * (t311 / 0.2e1 + t312 / 0.2e1);
t171 = t375 + t326;
t285 = t293 * qJ(4);
t335 = -t171 * qJD(2) - t285 * qJD(3) + t541;
t320 = t531 / 0.2e1 + t469 / 0.2e1 + t455 / 0.2e1 - (t314 * t401 - t373) * t500 / 0.2e1 + t168 * t374;
t328 = -t531 / 0.2e1 - t466 / 0.2e1 - t453 / 0.2e1;
t11 = t320 + t328;
t334 = -t11 * qJD(1) - t273 * t376;
t43 = (t511 - t368 / 0.2e1) * t313 + t322 * t314;
t333 = -t43 * qJD(1) - t313 * t376;
t332 = (-t384 - t383) * pkin(2);
t325 = (-t372 / 0.2e1 + t313 * t370) * pkin(2);
t107 = t325 - t527;
t17 = t392 - t346 + t352;
t331 = t17 * qJD(1) + t107 * qJD(2) - t273 * t423;
t324 = (t314 * t370 + t373 / 0.2e1) * pkin(2);
t108 = t389 * t271 + t324;
t18 = t393 - t347 + t354;
t330 = t18 * qJD(1) + t108 * qJD(2) + t271 * t423;
t292 = t313 * t377;
t260 = t273 * qJD(4);
t258 = t271 * qJD(4);
t246 = t273 * t377;
t245 = t271 * t377;
t204 = t271 * t413;
t198 = t409 * t274;
t192 = t409 * t264;
t177 = t185 * qJD(4);
t170 = t374 + t326;
t161 = 0.2e1 * t269 - t277;
t148 = t156 * qJD(5);
t125 = t474 / 0.2e1;
t124 = -t475 / 0.2e1;
t116 = t167 / 0.2e1 + t338;
t115 = -t167 / 0.2e1 + t338;
t110 = t325 + t527;
t109 = t324 + (t286 + t303) * t510;
t101 = t111 * qJD(4);
t100 = t112 * qJD(4);
t98 = t114 * qJD(4);
t97 = t116 * qJD(4);
t84 = t247 / 0.2e1 + (t390 - t344) * t276;
t83 = t273 * t506 + t345 * t276;
t75 = -t413 - t95;
t74 = -t259 - t96;
t73 = -t259 + t99;
t72 = t413 + t102;
t58 = t59 * qJD(5);
t54 = t55 * qJD(5);
t46 = -t113 * qJD(5) - t433;
t45 = -t112 * qJD(5) - t434;
t44 = t276 * t391 + t521 * t314 + t529;
t42 = 0.2e1 * t512 + t276 * t452 / 0.2e1 + t521 * t313;
t40 = -t167 * t429 + t58;
t37 = t115 * qJD(5) - t273 * t525 + t433;
t36 = -t111 * qJD(5) + t271 * t525 + t434;
t29 = t350 - t351;
t25 = t58 + (t526 + t429) * t167;
t23 = t54 - t437;
t20 = t125 + t346 + t352;
t19 = t124 + t347 + t354;
t15 = t125 + t348 + t353;
t14 = t124 + t349 + t355;
t10 = t320 - t328;
t8 = t321 - t329;
t7 = t437 + t54 + t409 * (t273 * t163 - t167 * t271);
t5 = t118 * t391 - t117 * t452 / 0.2e1 - t499 / 0.2e1 + t319;
t31 = [0, 0, 0, t316 * t410, t294 * qJD(2), 0, 0, 0, -pkin(1) * t411, -pkin(1) * t410, t544, t409 * t159, 0, 0, 0, t181 * qJD(2) - t276 * t422, t182 * qJD(2) - t274 * t422, t33 * qJD(2) + t35 * qJD(3) + t253 * t276, t32 * qJD(2) + t34 * qJD(3) - t421 * t457, -t16 * qJD(2) - t21 * qJD(3) + t121 * qJD(4), t22 * qJD(2) + t24 * qJD(3) + t51 * qJD(4), (t167 * t409 + t418) * t168, t57 * qJD(5) + t409 * t56, t274 * t418 + t409 * t62, -t168 * t274 * qJD(5) + t409 * t61, -t544, t1 * qJD(2) + t3 * qJD(3) + t27 * qJD(5) - t162 * t421, t2 * qJD(2) + t4 * qJD(3) + t26 * qJD(5) - t340 * t421; 0, 0, 0, t394, t412, t410, -t411, 0, -pkin(6) * t410 - t405, pkin(6) * t411 - t404, t211, t419, -t198, t525, 0, -qJD(2) * t368 + t161 * qJD(3) + t417, t416 + t536, t443 + (t313 * t359 - t528) * qJD(2) + t42 * qJD(3) - t177, t444 + (t314 * t359 + t529) * qJD(2) + t44 * qJD(3) - t253, qJD(2) * t361 + t29 * qJD(3) - t471, t463 + (t302 * t361 + t307 * t368) * qJD(2) + t5 * qJD(3) + t539, t25, t7, t37, t36, -t145, t495 + (-t456 + t470 + t532) * qJD(2) + t8 * qJD(3) - t100 + t15 * qJD(5), t494 + (t455 + t469 + t531) * qJD(2) + t10 * qJD(3) + t97 + t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t419, -t198, t525, 0, t161 * qJD(2) - qJD(3) * t368 - t395, -t396 + t536, t441 + t42 * qJD(2) + (t313 * t366 - t528) * qJD(3) - t177, t442 + t44 * qJD(2) + (t314 * t366 + t529) * qJD(3) - t253, t29 * qJD(2) + qJD(3) * t360 - t467, t462 + t5 * qJD(2) + (qJ(4) * t360 - t499) * qJD(3) + t539, t25, t7, t37, t36, -t145, t493 + t8 * qJD(2) + (-t454 + t468 + t532) * qJD(3) - t100 + t20 * qJD(5), t484 + t10 * qJD(2) + (t453 + t466 + t531) * qJD(3) + t97 + t19 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 * t409 + t371, (-t313 * t427 - t314 * t409) * t274, t420, t409 * t522 + t439, 0, 0, 0, 0, 0, t84 * qJD(5) - t399 - t533, t83 * qJD(5) + t116 * t409 - t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t519, t518, t115 * t409 - t340 * t378, t168 * t378 - t534, t192, t15 * qJD(2) + t20 * qJD(3) + t84 * qJD(4) - t39 * qJD(5) + t460, t14 * qJD(2) + t19 * qJD(3) + t83 * qJD(4) + t38 * qJD(5) + t461; 0, 0, 0, -t394, -t412, 0, 0, 0, t405, t404, -t211, -t419, 0, 0, 0, t213 * qJD(3) - t417, -t416, t41 * qJD(3) - t443, t43 * qJD(3) - t444, t28 * qJD(3) + t471, t6 * qJD(3) - t463 - t540, t40, t23, t46, t45, t145, t9 * qJD(3) - t12 * qJD(5) - t101 - t495, t11 * qJD(3) - t13 * qJD(5) - t494 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377, -pkin(2) * t385, -t314 * t377, t292, t430, t174 * qJD(3) + t232 * qJD(4), -t204, t148, 0, 0, 0, t286 * t413 + t245, -t259 * t286 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332 + t415, (-t386 - t385) * pkin(2), t314 * t332 + t440, t292 - t333, -t357 + t430, (-pkin(3) * t501 + qJ(4) * t339) * pkin(2) * qJD(3) + t170 * qJD(4) + t363, -t204, t148, 0, 0, 0, t110 * qJD(5) + t245 - t336, t109 * qJD(5) + t246 - t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t524, t170 * qJD(3) + t358, 0, 0, 0, 0, 0, -t102, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t30, t74, t75, -t414, t110 * qJD(3) - t173 * qJD(5) - t343, t109 * qJD(3) + t172 * qJD(5) - t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, -t419, 0, 0, 0, -t213 * qJD(2) + t395, t396, -t41 * qJD(2) - t441, -t43 * qJD(2) - t442, -t28 * qJD(2) + t467, -t6 * qJD(2) - t462 - t540, t40, t23, t46, t45, t145, -t9 * qJD(2) - t17 * qJD(5) - t101 - t493, -t11 * qJD(2) - t18 * qJD(5) - t484 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t376 - t415, pkin(2) * t386, t314 * t376 - t440, t333, t289 + t357, t171 * qJD(4) - t363, -t204, t148, 0, 0, 0, -t107 * qJD(5) + t336, -t108 * qJD(5) + t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, t285 * qJD(4), -t204, t148, 0, 0, 0, t303 * t413, -t303 * t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t524, -t335, 0, 0, 0, 0, 0, -t102, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t30, t74, t75, -t414, -t210 * qJD(5) - t331, t209 * qJD(5) - t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t371, t313 * t211, -t420, t409 * t523 - t439, 0, 0, 0, 0, 0, -t81 * qJD(5) + t399 + t534, t114 * t409 + t398 + t418; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t524, -t171 * qJD(3) - t358, 0, 0, 0, 0, 0, t72, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t524, t335, 0, 0, 0, 0, 0, t72, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t542, -t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t519, -t518, t113 * t409 - t398, t168 * t428 + t533, t192, t12 * qJD(2) + t17 * qJD(3) + t81 * qJD(4) - t460, t13 * qJD(2) + t18 * qJD(3) - qJD(4) * t340 - t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t30, t96, t95, t414, t107 * qJD(3) - t260 + t343, t108 * qJD(3) + t258 + t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t30, t96, t95, t414, -t260 + t331, t258 + t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t542, t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t31;
