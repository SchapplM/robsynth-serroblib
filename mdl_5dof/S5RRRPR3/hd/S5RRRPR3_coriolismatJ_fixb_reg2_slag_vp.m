% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRPR3
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
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:34
% EndTime: 2020-01-03 12:09:44
% DurationCPUTime: 7.24s
% Computational Cost: add. (11112->460), mult. (21044->535), div. (0->0), fcn. (22951->8), ass. (0->356)
t554 = qJD(3) + qJD(5);
t438 = qJD(1) + qJD(2);
t348 = sin(qJ(3));
t497 = cos(pkin(9));
t412 = t497 * t348;
t346 = sin(pkin(9));
t350 = cos(qJ(3));
t473 = t346 * t350;
t315 = t412 + t473;
t347 = sin(qJ(5));
t411 = t497 * t350;
t474 = t346 * t348;
t370 = t411 - t474;
t509 = cos(qJ(5));
t249 = t347 * t315 - t370 * t509;
t487 = t249 ^ 2;
t299 = t509 * t315;
t464 = t347 * t370;
t526 = t299 + t464;
t488 = t526 ^ 2;
t551 = t487 - t488;
t553 = t438 * t551;
t547 = t438 * t526;
t552 = t249 * t547;
t351 = cos(qJ(2));
t504 = t351 * pkin(1);
t277 = t315 * t504;
t426 = t509 * t277;
t278 = t370 * t504;
t465 = t347 * t278;
t366 = -t465 / 0.2e1 - t426 / 0.2e1;
t338 = -t350 * pkin(3) - pkin(2);
t276 = -pkin(4) * t370 + t338;
t269 = t276 - t504;
t415 = t276 / 0.2e1 + t269 / 0.2e1;
t539 = t526 * t415;
t66 = t366 - t539;
t528 = t438 * t315;
t549 = t370 * t528;
t548 = t438 * t249;
t546 = t554 * t249;
t545 = t249 * qJD(4);
t425 = t509 * t278;
t466 = t347 * t277;
t367 = t466 / 0.2e1 - t425 / 0.2e1;
t67 = t249 * t415 + t367;
t544 = t346 / 0.2e1;
t543 = -t497 / 0.2e1;
t89 = t487 + t488;
t542 = t438 * t89;
t395 = t497 * pkin(3) + pkin(4);
t508 = pkin(3) * t346;
t296 = t347 * t508 - t509 * t395;
t297 = t347 * t395 + t509 * t508;
t512 = -t297 / 0.2e1;
t373 = -t526 * t296 / 0.2e1 - t249 * t512;
t506 = t348 * pkin(3);
t339 = t506 / 0.2e1;
t507 = t315 * pkin(4);
t397 = t339 + t507 / 0.2e1;
t110 = t373 + t397;
t538 = t438 * t110;
t396 = t299 / 0.2e1;
t191 = 0.2e1 * t396 + t464;
t535 = t438 * t191;
t364 = t315 * t543 + t370 * t544;
t510 = -t348 / 0.2e1;
t213 = (t510 + t364) * pkin(3);
t534 = t438 * t213;
t312 = t370 ^ 2;
t524 = t315 ^ 2;
t223 = t312 - t524;
t533 = t438 * t223;
t247 = t396 - t299 / 0.2e1;
t532 = t438 * t247;
t254 = t312 + t524;
t530 = t438 * t254;
t529 = t438 * t370;
t344 = t348 ^ 2;
t345 = t350 ^ 2;
t333 = t345 - t344;
t527 = t438 * t333;
t337 = -pkin(2) - t504;
t433 = -t504 / 0.2e1;
t525 = t433 + t337 / 0.2e1;
t349 = sin(qJ(2));
t505 = t349 * pkin(1);
t336 = pkin(7) + t505;
t461 = qJ(4) + t336;
t392 = t461 * t497;
t293 = t350 * t392;
t409 = t346 * t461;
t291 = t348 * t409;
t309 = t370 * pkin(8);
t459 = -t309 + t291;
t192 = -t293 + t459;
t432 = t509 * t192;
t377 = t348 * t392;
t241 = -t350 * t409 - t377;
t310 = t315 * pkin(8);
t193 = -t310 + t241;
t471 = t347 * t193;
t115 = t432 - t471;
t523 = -t115 / 0.2e1;
t503 = -qJ(4) - pkin(7);
t413 = t346 * t503;
t317 = t348 * t413;
t393 = t503 * t497;
t257 = t350 * t393 - t317;
t224 = -t309 + t257;
t429 = t509 * t224;
t379 = t348 * t393;
t259 = t350 * t413 + t379;
t225 = -t310 + t259;
t468 = t347 * t225;
t128 = t429 - t468;
t522 = -t128 / 0.2e1;
t408 = t461 * t350;
t240 = -t346 * t408 - t377;
t357 = t240 - t310;
t173 = t509 * t357;
t521 = -t173 / 0.2e1;
t258 = t503 * t473 + t379;
t358 = t258 - t310;
t206 = t509 * t358;
t520 = -t206 / 0.2e1;
t239 = t291 - t293;
t519 = -t239 / 0.2e1;
t518 = -t241 / 0.2e1;
t516 = -t257 / 0.2e1;
t515 = -t259 / 0.2e1;
t502 = pkin(1) * qJD(1);
t501 = pkin(1) * qJD(2);
t500 = pkin(2) * qJD(2);
t499 = qJD(3) * pkin(3);
t498 = (-t249 * t296 - t297 * t526) * qJD(3);
t355 = t347 * t357;
t292 = t497 * t408;
t194 = t292 - t459;
t430 = t509 * t194;
t118 = t430 + t355;
t493 = t118 * t249;
t470 = t347 * t194;
t116 = -t173 + t470;
t494 = t116 * t526;
t28 = -t493 + t494;
t496 = qJD(1) * t28;
t195 = -t426 - t465;
t196 = t425 - t466;
t65 = -t195 * t526 - t196 * t249;
t495 = qJD(1) * t65;
t431 = t509 * t193;
t472 = t347 * t192;
t117 = t431 + t472;
t279 = t506 + t507;
t17 = -t115 * t116 + t117 * t118 + t269 * t279;
t492 = t17 * qJD(1);
t491 = t240 * t315;
t242 = t292 - t291;
t490 = t242 * t370;
t260 = -t503 * t411 + t317;
t226 = t309 + t260;
t467 = t347 * t226;
t129 = -t206 + t467;
t489 = t526 * t129;
t25 = -t116 * t195 + t118 * t196 + t269 * t505;
t486 = t25 * qJD(1);
t356 = t347 * t358;
t427 = t509 * t226;
t131 = t427 + t356;
t485 = t249 * t131;
t484 = t258 * t315;
t483 = t260 * t370;
t482 = t269 * t526;
t481 = t269 * t249;
t480 = t276 * t526;
t479 = t276 * t249;
t327 = t338 - t504;
t478 = t327 * t370;
t477 = t327 * t315;
t476 = t338 * t370;
t475 = t338 * t315;
t469 = t347 * t224;
t71 = t239 * t240 + t241 * t242 + t327 * t506;
t463 = t71 * qJD(1);
t91 = -t240 * t277 + t242 * t278 + t327 * t505;
t462 = t91 * qJD(1);
t146 = t277 * t315 + t278 * t370;
t253 = t254 * qJD(4);
t460 = t146 * qJD(2) + t253;
t174 = t279 * t249;
t113 = t174 + t482;
t458 = qJD(1) * t113;
t175 = t279 * t526;
t114 = t175 - t481;
t457 = qJD(1) * t114;
t127 = t490 - t491;
t456 = qJD(1) * t127;
t455 = qJD(1) * t146;
t294 = t370 * t506;
t221 = -t294 + t477;
t454 = qJD(1) * t221;
t295 = t315 * t506;
t222 = t295 + t478;
t453 = qJD(1) * t222;
t452 = qJD(1) * t269;
t451 = qJD(1) * t337;
t450 = qJD(2) * t276;
t449 = qJD(3) * t526;
t446 = qJD(5) * t526;
t445 = qJD(5) * t269;
t444 = qJD(5) * t276;
t407 = (t344 + t345) * t351;
t245 = (t336 * t407 + t337 * t349) * pkin(1);
t443 = t245 * qJD(1);
t442 = t247 * qJD(5);
t311 = pkin(1) * t407;
t441 = t311 * qJD(1);
t440 = t315 * qJD(3);
t439 = t348 * qJD(3);
t343 = t350 * qJD(3);
t428 = t509 * t225;
t130 = t428 + t469;
t418 = -t131 / 0.2e1 - t118 / 0.2e1;
t419 = t129 / 0.2e1 + t116 / 0.2e1;
t4 = (t522 + t523 + t418) * t526 + (-t419 - t130 / 0.2e1 - t117 / 0.2e1) * t249;
t88 = t89 * qJD(4);
t437 = t4 * qJD(3) + t88;
t436 = t349 * t501;
t435 = t349 * t502;
t434 = t505 / 0.2e1;
t424 = t249 * t452;
t423 = t526 * t452;
t422 = t348 * t451;
t421 = t350 * t451;
t420 = t370 * t440;
t417 = t258 / 0.2e1 + t240 / 0.2e1;
t416 = t260 / 0.2e1 + t242 / 0.2e1;
t414 = t338 / 0.2e1 + t327 / 0.2e1;
t410 = pkin(1) * t438;
t402 = t249 * t435;
t401 = t526 * t435;
t400 = t370 * t435;
t399 = t315 * t435;
t398 = t348 * t435;
t394 = t349 * t410;
t10 = (-t115 - t118) * t526 + (-t116 - t117) * t249;
t391 = -qJD(1) * t10 - qJD(2) * t4;
t16 = (-t128 - t131) * t526 + (-t129 - t130) * t249;
t390 = -qJD(1) * t4 - qJD(2) * t16;
t21 = -t128 * t129 + t130 * t131 + t279 * t276;
t354 = t415 * t279 + t129 * t523 + t116 * t522 + t117 * t131 / 0.2e1 + t118 * t130 / 0.2e1;
t374 = t195 * t296 / 0.2e1 + t196 * t512;
t6 = t354 + t374;
t389 = t6 * qJD(1) + t21 * qJD(2);
t11 = -t249 * t418 - t419 * t526 + t434;
t33 = -t485 + t489;
t388 = qJD(1) * t11 - qJD(2) * t33;
t23 = (t516 + t519 - t416) * t315 - (t515 + t518 + t417) * t370;
t49 = (-t239 - t242) * t315 - (t240 - t241) * t370;
t387 = -qJD(1) * t49 - qJD(2) * t23;
t70 = (-t257 - t260) * t315 - (t258 - t259) * t370;
t386 = -qJD(1) * t23 - qJD(2) * t70;
t112 = t257 * t258 + t259 * t260 + t338 * t506;
t359 = t240 * t516 + t242 * t515 + t258 * t519 + t260 * t518;
t365 = t277 * t543 + t278 * t544;
t26 = (-t348 * t414 + t365) * pkin(3) + t359;
t385 = -t26 * qJD(1) + t112 * qJD(2);
t119 = t174 + t480;
t39 = -t174 + t66;
t384 = qJD(1) * t39 - qJD(2) * t119;
t120 = t175 - t479;
t40 = -t175 + t67;
t383 = qJD(1) * t40 - qJD(2) * t120;
t143 = t483 - t484;
t58 = t315 * t417 - t370 * t416 + t434;
t382 = -qJD(1) * t58 + qJD(2) * t143;
t361 = (-t473 / 0.2e1 - t412 / 0.2e1) * t504;
t139 = -t315 * t414 + t294 + t361;
t230 = -t294 + t475;
t381 = qJD(1) * t139 - qJD(2) * t230;
t360 = (-t411 / 0.2e1 + t474 / 0.2e1) * t504;
t140 = -t370 * t414 - t295 + t360;
t231 = t295 + t476;
t380 = qJD(1) * t140 - qJD(2) * t231;
t134 = qJD(5) * t191 + t449;
t378 = t433 + pkin(2) / 0.2e1 - t337 / 0.2e1;
t270 = t378 * t348;
t376 = qJD(1) * t270 + t348 * t500;
t271 = t378 * t350;
t375 = qJD(1) * t271 + t350 * t500;
t372 = qJD(1) * t67 + t249 * t450;
t371 = qJD(1) * t66 - t450 * t526;
t369 = -t471 / 0.2e1 + t432 / 0.2e1;
t368 = -t468 / 0.2e1 + t429 / 0.2e1;
t36 = t521 + t431 / 0.2e1 + (t194 / 0.2e1 + t192 / 0.2e1) * t347;
t61 = t520 + t428 / 0.2e1 + (t226 / 0.2e1 + t224 / 0.2e1) * t347;
t363 = -qJD(1) * t36 - qJD(2) * t61 - qJD(3) * t296;
t352 = t355 / 0.2e1 + t430 / 0.2e1;
t35 = t352 + t369;
t353 = t356 / 0.2e1 + t427 / 0.2e1;
t60 = t353 + t368;
t362 = qJD(1) * t35 + qJD(2) * t60 + qJD(3) * t297;
t334 = t348 * t343;
t332 = t348 * t436;
t329 = t333 * qJD(3);
t307 = t315 * qJD(4);
t305 = t370 * qJD(3);
t304 = t370 * qJD(4);
t303 = t438 * t350 * t348;
t300 = t311 * qJD(2);
t287 = t315 * t436;
t286 = t370 * t436;
t283 = t297 * qJD(5);
t282 = t296 * qJD(5);
t273 = (-pkin(2) / 0.2e1 + t525) * t350;
t272 = pkin(2) * t510 + t348 * t525;
t243 = t526 * qJD(4);
t237 = t247 * qJD(3);
t236 = t247 * qJD(4);
t229 = (-t315 * t346 - t370 * t497) * t499;
t228 = t526 * t436;
t227 = t249 * t436;
t212 = pkin(3) * t364 + t339;
t211 = t223 * qJD(3);
t210 = t213 * qJD(3);
t209 = t213 * qJD(4);
t208 = t212 * qJD(3);
t207 = t212 * qJD(4);
t177 = t191 * qJD(4);
t142 = t295 + t476 / 0.2e1 + t478 / 0.2e1 + t360;
t141 = -t294 + t475 / 0.2e1 + t477 / 0.2e1 + t361;
t133 = -qJD(3) * t191 - t446;
t111 = -t373 + t397;
t99 = t110 * qJD(4);
t98 = t111 * qJD(3);
t97 = t111 * qJD(4);
t96 = t110 * qJD(3);
t86 = t546 * t526;
t82 = (t446 + t449) * t249;
t69 = t366 + t539;
t68 = t367 - (t269 + t276) * t249 / 0.2e1;
t64 = t65 * qJD(2);
t63 = -t353 + t368;
t62 = t467 / 0.2e1 + t520 - t428 / 0.2e1 - t469 / 0.2e1;
t59 = t483 / 0.2e1 + t490 / 0.2e1 - t484 / 0.2e1 - t491 / 0.2e1 + t434;
t42 = t175 - t479 / 0.2e1 - t481 / 0.2e1 + t367;
t41 = t174 + t480 / 0.2e1 + t482 / 0.2e1 + t366;
t38 = -t352 + t369;
t37 = t470 / 0.2e1 + t521 - t431 / 0.2e1 - t472 / 0.2e1;
t34 = t554 * t551;
t27 = pkin(3) * t365 - t359 + (t327 + t338) * t339;
t22 = t23 * qJD(3);
t12 = -t485 / 0.2e1 - t493 / 0.2e1 + t489 / 0.2e1 + t494 / 0.2e1 + t434;
t5 = t354 - t374;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t436, -t351 * t501, 0, 0, t334, t329, 0, -t334, 0, 0, t337 * t439 - t350 * t436, t337 * t343 + t332, t300, t245 * qJD(2), t420, t211, 0, -t420, 0, 0, qJD(3) * t221 - t286, qJD(3) * t222 + t287, qJD(3) * t49 + t460, qJD(2) * t91 + qJD(3) * t71 + qJD(4) * t127, -t86, t34, 0, t82, 0, 0, qJD(3) * t113 + t445 * t526 + t227, qJD(3) * t114 - t249 * t445 + t228, qJD(3) * t10 + t64 + t88, qJD(2) * t25 + qJD(3) * t17 + qJD(4) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t394, -t351 * t410, 0, 0, t334, t329, 0, -t334, 0, 0, t272 * qJD(3) - t350 * t394, qJD(3) * t273 + t332 + t398, t300 + t441, t443 + (-pkin(2) * t349 + pkin(7) * t407) * t501, t420, t211, 0, -t420, 0, 0, qJD(3) * t141 - t286 - t400, qJD(3) * t142 + t287 + t399, t22 + t455 + t460, t462 + (-t258 * t277 + t260 * t278 + t338 * t505) * qJD(2) + t27 * qJD(3) + t59 * qJD(4), -t86, t34, 0, t82, 0, 0, qJD(3) * t41 + qJD(5) * t69 + t227 + t402, qJD(3) * t42 + qJD(5) * t68 + t228 + t401, t437 + t64 + t495, t486 + (-t129 * t195 + t131 * t196 + t276 * t505) * qJD(2) + t5 * qJD(3) + t12 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, t527, t343, -t303, -t439, 0, qJD(2) * t272 - t336 * t343 + t422, qJD(2) * t273 + t336 * t439 + t421, 0, 0, t549, t533, t305, -t549, -t440, 0, qJD(2) * t141 + qJD(3) * t239 + t454, qJD(2) * t142 - qJD(3) * t241 + t453, t229 - t387, t463 + t27 * qJD(2) + (t239 * t497 + t241 * t346) * t499 + t207, -t552, t553, -t546, t552, -t134, 0, qJD(2) * t41 + qJD(3) * t115 + qJD(5) * t38 + t458, qJD(2) * t42 - qJD(3) * t117 + qJD(5) * t37 + t457, -t391 + t498, t492 + t5 * qJD(2) + (-t115 * t296 + t117 * t297) * qJD(3) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530, qJD(2) * t59 + t208 + t456, 0, 0, 0, 0, 0, 0, t442, 0, t542, qJD(2) * t12 + t496 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t552, t553, -t546, t552, t133, 0, qJD(2) * t69 + qJD(3) * t38 - qJD(5) * t118 + t236 + t423, qJD(2) * t68 + qJD(3) * t37 + qJD(5) * t116 - t424, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t435, t351 * t502, 0, 0, t334, t329, 0, -t334, 0, 0, -qJD(3) * t270 + t350 * t435, -qJD(3) * t271 - t398, -t441, -t443, t420, t211, 0, -t420, 0, 0, -qJD(3) * t139 + t400, -qJD(3) * t140 - t399, t22 + t253 - t455, -qJD(3) * t26 - qJD(4) * t58 - t462, -t86, t34, 0, t82, 0, 0, -qJD(3) * t39 - qJD(5) * t66 - t402, -qJD(3) * t40 - qJD(5) * t67 - t401, t437 - t495, qJD(3) * t6 - qJD(4) * t11 - t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t329, 0, -t334, 0, 0, -pkin(2) * t439, -pkin(2) * t343, 0, 0, t420, t211, 0, -t420, 0, 0, t230 * qJD(3), t231 * qJD(3), qJD(3) * t70 + t253, qJD(3) * t112 + qJD(4) * t143, -t86, t34, 0, t82, 0, 0, qJD(3) * t119 + t444 * t526, qJD(3) * t120 - t249 * t444, qJD(3) * t16 + t88, qJD(3) * t21 + qJD(4) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, t527, t343, -t303, -t439, 0, -pkin(7) * t343 - t376, pkin(7) * t439 - t375, 0, 0, t549, t533, t305, -t549, -t440, 0, qJD(3) * t257 - t381, -qJD(3) * t259 - t380, t229 - t386, (t257 * t497 + t259 * t346) * t499 + t207 + t385, -t552, t553, -t546, t552, -t134, 0, qJD(3) * t128 + qJD(5) * t63 - t384, -qJD(3) * t130 + qJD(5) * t62 - t383, -t390 + t498, (-t128 * t296 + t130 * t297) * qJD(3) + t97 + t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530, t208 + t382, 0, 0, 0, 0, 0, 0, t442, 0, t542, -t388 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t552, t553, -t546, t552, t133, 0, qJD(3) * t63 - qJD(5) * t131 + t236 - t371, qJD(3) * t62 + qJD(5) * t129 - t372, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303, -t527, 0, t303, 0, 0, qJD(2) * t270 - t422, qJD(2) * t271 - t421, 0, 0, -t549, -t533, 0, t549, 0, 0, qJD(2) * t139 - t307 - t454, qJD(2) * t140 - t304 - t453, t387, qJD(2) * t26 + t209 - t463, t552, -t553, 0, -t552, -t442, 0, qJD(2) * t39 - qJD(5) * t35 - t243 - t458, qJD(2) * t40 + qJD(5) * t36 - t457 + t545, t391, -qJD(2) * t6 - t492 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303, -t527, 0, t303, 0, 0, t376, t375, 0, 0, -t549, -t533, 0, t549, 0, 0, -t307 + t381, -t304 + t380, t386, t209 - t385, t552, -t553, 0, -t552, -t442, 0, -qJD(5) * t60 - t243 + t384, qJD(5) * t61 + t383 + t545, t390, -t389 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t283, t282, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t528, -t529, 0, t534, 0, 0, 0, 0, 0, 0, -t547, t548, 0, -t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t532, 0, -t283 - t362, t282 - t363, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t440, t305, -t530, qJD(2) * t58 - t210 - t456, 0, 0, 0, 0, 0, 0, t134, -t546, -t542, qJD(2) * t11 - t496 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t440, t305, -t530, -t210 - t382, 0, 0, 0, 0, 0, 0, t134, -t546, -t542, t388 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t528, t529, 0, -t534, 0, 0, 0, 0, 0, 0, t547, -t548, 0, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t535, -t548, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t552, -t553, 0, -t552, t237, 0, qJD(2) * t66 + qJD(3) * t35 - t177 - t423, qJD(2) * t67 - qJD(3) * t36 + t424 + t545, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t552, -t553, 0, -t552, t237, 0, qJD(3) * t60 - t177 + t371, -qJD(3) * t61 + t372 + t545, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t532, 0, t362, t363, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t535, t548, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;