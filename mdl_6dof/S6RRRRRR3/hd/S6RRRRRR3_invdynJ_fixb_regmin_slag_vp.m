% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x38]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:44
% EndTime: 2019-03-10 03:43:09
% DurationCPUTime: 11.64s
% Computational Cost: add. (16304->614), mult. (35545->810), div. (0->0), fcn. (27894->18), ass. (0->340)
t355 = cos(qJ(2));
t508 = pkin(7) + pkin(8);
t303 = t508 * t355;
t285 = qJD(1) * t303;
t349 = sin(qJ(3));
t260 = t349 * t285;
t350 = sin(qJ(2));
t301 = t508 * t350;
t283 = qJD(1) * t301;
t506 = cos(qJ(3));
t198 = -t283 * t506 - t260;
t428 = qJD(3) * t506;
t553 = -pkin(2) * t428 + t198;
t509 = qJD(4) + qJD(5);
t429 = qJD(1) * t506;
t451 = qJD(1) * t350;
t542 = -t349 * t451 + t355 * t429;
t552 = t509 - t542;
t341 = qJD(2) + qJD(3);
t551 = t341 * t542;
t345 = qJ(2) + qJ(3);
t337 = cos(t345);
t502 = g(3) * t337;
t335 = sin(t345);
t351 = sin(qJ(1));
t356 = cos(qJ(1));
t402 = g(1) * t356 + g(2) * t351;
t549 = t402 * t335;
t531 = -t549 + t502;
t468 = t349 * t355;
t256 = -qJD(1) * t468 - t350 * t429;
t186 = -pkin(3) * t256 - pkin(9) * t542;
t172 = pkin(2) * t451 + t186;
t348 = sin(qJ(4));
t354 = cos(qJ(4));
t550 = -t354 * t172 + t348 * t553;
t213 = -t256 * t348 - t354 * t341;
t347 = sin(qJ(5));
t353 = cos(qJ(5));
t390 = t256 * t354 - t341 * t348;
t141 = t353 * t213 - t347 * t390;
t346 = sin(qJ(6));
t352 = cos(qJ(6));
t391 = t213 * t347 + t353 * t390;
t393 = t141 * t346 + t352 * t391;
t71 = t352 * t141 - t346 * t391;
t539 = t393 * t71;
t473 = t347 * t354;
t277 = t348 * t353 + t473;
t533 = t552 * t277;
t275 = t347 * t348 - t353 * t354;
t460 = t552 * t275;
t548 = t348 * t172 + t354 * t553;
t528 = t393 ^ 2 - t71 ^ 2;
t418 = qJDD(1) * t506;
t441 = t355 * qJDD(1);
t152 = t349 * t441 + t350 * t418 + t551;
t340 = qJDD(2) + qJDD(3);
t448 = qJD(4) * t354;
t449 = qJD(4) * t348;
t119 = t354 * t152 + t256 * t449 + t348 * t340 + t341 * t448;
t120 = -qJD(4) * t390 + t348 * t152 - t354 * t340;
t367 = qJD(5) * t391 - t347 * t119 - t353 * t120;
t444 = qJD(6) * t352;
t445 = qJD(6) * t346;
t446 = qJD(5) * t353;
t447 = qJD(5) * t347;
t50 = t353 * t119 - t347 * t120 - t213 * t446 + t390 * t447;
t11 = -t141 * t444 + t346 * t367 + t352 * t50 + t391 * t445;
t248 = qJD(4) - t542;
t236 = qJD(5) + t248;
t232 = qJD(6) + t236;
t526 = t232 * t71 + t11;
t344 = qJ(4) + qJ(5);
t338 = qJ(6) + t344;
t325 = sin(t338);
t326 = cos(t338);
t476 = t337 * t351;
t222 = t325 * t356 - t326 * t476;
t475 = t337 * t356;
t224 = t325 * t351 + t326 * t475;
t331 = -pkin(2) * t355 - pkin(1);
t299 = t331 * qJD(1);
t169 = -pkin(3) * t542 + t256 * pkin(9) + t299;
t263 = t506 * t285;
t495 = qJD(2) * pkin(2);
t264 = -t283 + t495;
t191 = t349 * t264 + t263;
t176 = pkin(9) * t341 + t191;
t122 = t354 * t169 - t176 * t348;
t87 = pkin(10) * t390 + t122;
t67 = pkin(4) * t248 + t87;
t123 = t169 * t348 + t176 * t354;
t88 = -pkin(10) * t213 + t123;
t84 = t353 * t88;
t41 = t347 * t67 + t84;
t541 = pkin(11) * t141;
t30 = t41 - t541;
t28 = t30 * t445;
t322 = g(3) * t335;
t190 = t264 * t506 - t260;
t175 = -t341 * pkin(3) - t190;
t134 = t213 * pkin(4) + t175;
t68 = t141 * pkin(5) + t134;
t525 = g(1) * t224 - g(2) * t222 + t326 * t322 + t68 * t71 + t28;
t221 = t325 * t476 + t326 * t356;
t223 = -t325 * t475 + t326 * t351;
t278 = t350 * t506 + t468;
t203 = t341 * t278;
t442 = t350 * qJDD(1);
t396 = t349 * t442 - t355 * t418;
t153 = qJD(1) * t203 + t396;
t151 = qJDD(4) + t153;
t149 = qJDD(5) + t151;
t443 = qJD(1) * qJD(2);
t427 = t350 * t443;
t250 = pkin(2) * t427 + qJDD(1) * t331;
t81 = t153 * pkin(3) - t152 * pkin(9) + t250;
t79 = t354 * t81;
t426 = t355 * t443;
t207 = qJDD(2) * pkin(2) + t508 * (-t426 - t442);
t212 = t508 * (-t427 + t441);
t450 = qJD(3) * t349;
t365 = t349 * t207 + t212 * t506 + t264 * t428 - t285 * t450;
t96 = t340 * pkin(9) + t365;
t17 = t151 * pkin(4) - t119 * pkin(10) - qJD(4) * t123 - t348 * t96 + t79;
t380 = t169 * t448 - t176 * t449 + t348 * t81 + t354 * t96;
t23 = -pkin(10) * t120 + t380;
t423 = t353 * t17 - t347 * t23;
t370 = -qJD(5) * t41 + t423;
t3 = t149 * pkin(5) - t50 * pkin(11) + t370;
t417 = -t347 * t17 - t353 * t23 - t67 * t446 + t88 * t447;
t4 = pkin(11) * t367 - t417;
t437 = t352 * t3 - t346 * t4;
t547 = -g(1) * t223 + g(2) * t221 + t325 * t322 + t68 * t393 + t437;
t369 = qJD(6) * t393 - t346 * t50 + t352 * t367;
t520 = -t232 * t393 + t369;
t339 = t354 * pkin(10);
t407 = -t256 * pkin(4) - t339 * t542;
t327 = pkin(2) * t349 + pkin(9);
t497 = -pkin(10) - t327;
t419 = qJD(4) * t497;
t546 = t354 * t419 - t407 + t550;
t180 = t354 * t186;
t507 = -pkin(9) - pkin(10);
t435 = qJD(4) * t507;
t545 = t190 * t348 + t354 * t435 - t180 - t407;
t481 = t542 * t348;
t440 = pkin(10) * t481;
t544 = -t348 * t419 - t440 + t548;
t458 = t348 * t186 + t354 * t190;
t543 = -t348 * t435 - t440 + t458;
t540 = pkin(11) * t391;
t194 = t352 * t275 + t277 * t346;
t491 = qJD(6) * t194 + t346 * t533 + t352 * t460;
t408 = -t506 * t207 + t349 * t212 + t264 * t450 + t285 * t428;
t97 = -pkin(3) * t340 + t408;
t538 = t97 + t502;
t537 = t391 * t141;
t195 = -t275 * t346 + t277 * t352;
t536 = qJD(6) * t195 - t346 * t460 + t352 * t533;
t197 = -t349 * t283 + t263;
t406 = pkin(2) * t450 - t197;
t385 = -t349 * t350 + t355 * t506;
t202 = t341 * t385;
t471 = t348 * t202;
t530 = t278 * t448 + t471;
t527 = -t141 ^ 2 + t391 ^ 2;
t524 = t141 * t236 + t50;
t82 = t347 * t88;
t40 = t353 * t67 - t82;
t29 = t40 + t540;
t26 = pkin(5) * t236 + t29;
t494 = t352 * t30;
t10 = t346 * t26 + t494;
t523 = -qJD(6) * t10 + t547;
t334 = sin(t344);
t336 = cos(t344);
t226 = t334 * t356 - t336 * t476;
t228 = t334 * t351 + t336 * t475;
t522 = g(1) * t228 - g(2) * t226 + t134 * t141 + t336 * t322 + t417;
t225 = t334 * t476 + t336 * t356;
t227 = -t334 * t475 + t336 * t351;
t521 = -g(1) * t227 + g(2) * t225 + t134 * t391 + t334 * t322 + t370;
t519 = -t236 * t391 + t367;
t518 = t546 * t353;
t517 = t545 * t353;
t181 = t277 * t278;
t189 = -pkin(3) * t385 - pkin(9) * t278 + t331;
t184 = t354 * t189;
t217 = -t349 * t301 + t303 * t506;
t478 = t278 * t354;
t107 = -pkin(4) * t385 - pkin(10) * t478 - t217 * t348 + t184;
t208 = t354 * t217;
t456 = t348 * t189 + t208;
t479 = t278 * t348;
t125 = -pkin(10) * t479 + t456;
t463 = t347 * t107 + t353 * t125;
t332 = pkin(4) * t449;
t516 = pkin(5) * t533 + t332;
t515 = t533 * pkin(11);
t268 = t497 * t348;
t269 = t327 * t354 + t339;
t455 = t347 * t268 + t353 * t269;
t300 = t507 * t348;
t302 = pkin(9) * t354 + t339;
t454 = t347 * t300 + t353 * t302;
t233 = pkin(4) * t481;
t514 = -t233 + t406;
t513 = -t300 * t446 + t302 * t447 - t347 * t545 + t353 * t543;
t512 = -t268 * t446 + t269 * t447 - t347 * t546 + t544 * t353;
t511 = -t256 * pkin(5) - pkin(11) * t460;
t510 = -t506 * t301 - t349 * t303;
t505 = pkin(11) * t277;
t499 = t275 * pkin(5);
t498 = t354 * pkin(4);
t496 = t353 * t87 - t82;
t490 = t119 * t348;
t147 = qJDD(6) + t149;
t489 = t147 * t347;
t488 = t175 * t542;
t487 = t202 * t354;
t486 = t213 * t248;
t485 = t390 * t248;
t484 = t232 * t256;
t483 = t236 * t256;
t482 = t248 * t256;
t480 = t256 * t542;
t328 = pkin(4) * t353 + pkin(5);
t477 = t328 * t147;
t474 = t347 * t352;
t472 = t348 * t151;
t470 = t348 * t351;
t469 = t348 * t356;
t467 = t351 * t354;
t466 = t354 * t151;
t465 = t354 * t356;
t462 = t514 + t516;
t150 = t233 + t191;
t459 = -t150 + t516;
t453 = t332 + t514;
t342 = t350 ^ 2;
t452 = -t355 ^ 2 + t342;
t439 = t350 * t495;
t438 = qJD(4) * pkin(9) * t248;
t330 = -pkin(3) - t498;
t436 = qJD(2) * t508;
t432 = t278 * t449;
t163 = t175 * t448;
t424 = qJD(6) * t26 + t4;
t133 = pkin(3) * t203 - pkin(9) * t202 + t439;
t130 = t354 * t133;
t284 = t350 * t436;
t286 = t355 * t436;
t145 = qJD(3) * t510 - t506 * t284 - t349 * t286;
t35 = -pkin(10) * t487 + t203 * pkin(4) - t348 * t145 + t130 + (-t208 + (pkin(10) * t278 - t189) * t348) * qJD(4);
t377 = t348 * t133 + t354 * t145 + t189 * t448 - t217 * t449;
t44 = -pkin(10) * t530 + t377;
t422 = -t347 * t44 + t353 * t35;
t420 = -t347 * t87 - t84;
t416 = -qJD(4) * t169 - t96;
t415 = t353 * t107 - t125 * t347;
t412 = t353 * t268 - t269 * t347;
t411 = t353 * t300 - t302 * t347;
t410 = t248 * t354;
t329 = -pkin(2) * t506 - pkin(3);
t405 = -t150 + t332;
t155 = t412 - t505;
t404 = -qJD(6) * t155 + t512 + t515;
t267 = t275 * pkin(11);
t156 = -t267 + t455;
t403 = t455 * qJD(5) + qJD(6) * t156 - t347 * t544 + t511 - t518;
t401 = g(1) * t351 - g(2) * t356;
t400 = -t176 * t448 + t79;
t173 = t411 - t505;
t399 = -qJD(6) * t173 + t513 + t515;
t174 = -t267 + t454;
t398 = t454 * qJD(5) + qJD(6) * t174 - t347 * t543 + t511 - t517;
t397 = -pkin(9) * t151 - t488;
t392 = -t151 * t327 - t488;
t182 = t275 * t278;
t126 = t352 * t181 - t182 * t346;
t127 = -t181 * t346 - t182 * t352;
t389 = -t123 * t256 + t348 * t538 + t163;
t388 = t122 * t256 + t175 * t449 + t354 * t549;
t170 = pkin(4) * t479 - t510;
t386 = -0.2e1 * pkin(1) * t443 - pkin(7) * qJDD(2);
t383 = -t432 + t487;
t382 = t107 * t446 - t125 * t447 + t347 * t35 + t353 * t44;
t298 = t329 - t498;
t146 = -t349 * t284 + t286 * t506 - t301 * t450 + t303 * t428;
t58 = pkin(4) * t120 + t97;
t93 = pkin(4) * t530 + t146;
t357 = qJD(2) ^ 2;
t372 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t357 + t401;
t358 = qJD(1) ^ 2;
t371 = pkin(1) * t358 - pkin(7) * qJDD(1) + t402;
t368 = t299 * t256 - t408 - t531;
t21 = -pkin(5) * t367 + t58;
t9 = t352 * t26 - t30 * t346;
t366 = t21 * t194 + t9 * t256 - t326 * t531 + t536 * t68;
t364 = t134 * t533 + t40 * t256 + t58 * t275 - t336 * t531;
t363 = -t10 * t256 + t21 * t195 + t325 * t531 - t491 * t68;
t362 = -t134 * t460 - t41 * t256 + t58 * t277 + t334 * t531;
t359 = g(1) * t475 + g(2) * t476 - t299 * t542 + t322 - t365;
t245 = t337 * t465 + t470;
t244 = -t337 * t469 + t467;
t243 = -t337 * t467 + t469;
t242 = t337 * t470 + t465;
t229 = t330 + t499;
t220 = t298 + t499;
t157 = t256 ^ 2 - t542 ^ 2;
t132 = -t396 + (-qJD(1) * t278 - t256) * t341;
t131 = t152 - t551;
t128 = t181 * pkin(5) + t170;
t121 = -pkin(4) * t390 - pkin(5) * t391;
t64 = t202 * t473 - t347 * t432 - t447 * t479 + (t478 * t509 + t471) * t353;
t63 = -t181 * t509 - t275 * t202;
t62 = t248 * t410 - t256 * t390 + t472;
t61 = -t248 ^ 2 * t348 - t213 * t256 + t466;
t59 = -t390 * t410 + t490;
t55 = -pkin(11) * t181 + t463;
t52 = -pkin(5) * t385 + pkin(11) * t182 + t415;
t46 = -t141 * t256 - t275 * t149 - t236 * t533;
t45 = t277 * t149 - t236 * t460 - t256 * t391;
t39 = t64 * pkin(5) + t93;
t32 = t496 + t540;
t31 = t420 + t541;
t27 = (t119 - t486) * t354 + (-t120 + t485) * t348;
t25 = qJD(6) * t127 + t346 * t63 + t352 * t64;
t24 = -qJD(6) * t126 - t346 * t64 + t352 * t63;
t20 = t50 * t277 + t391 * t460;
t19 = -t194 * t147 - t232 * t536 - t71 * t256;
t18 = t195 * t147 - t232 * t491 - t256 * t393;
t8 = -pkin(11) * t64 + t382;
t7 = t141 * t460 - t50 * t275 + t277 * t367 + t391 * t533;
t6 = t203 * pkin(5) - t63 * pkin(11) - qJD(5) * t463 + t422;
t5 = t11 * t195 + t393 * t491;
t1 = -t11 * t194 + t195 * t369 + t393 * t536 + t491 * t71;
t2 = [qJDD(1), t401, t402, qJDD(1) * t342 + 0.2e1 * t350 * t426, 0.2e1 * t350 * t441 - 0.2e1 * t443 * t452, qJDD(2) * t350 + t355 * t357, qJDD(2) * t355 - t350 * t357, 0, t350 * t386 + t355 * t372, -t350 * t372 + t355 * t386, t152 * t278 - t202 * t256, t152 * t385 - t153 * t278 + t202 * t542 + t203 * t256, t202 * t341 + t278 * t340, -t203 * t341 + t340 * t385, 0, -t146 * t341 + t331 * t153 + t299 * t203 - t250 * t385 + t337 * t401 + t340 * t510 - t439 * t542, -t145 * t341 + t331 * t152 + t299 * t202 - t217 * t340 + t250 * t278 - t256 * t439 - t335 * t401, t119 * t478 - t383 * t390 (-t213 * t354 + t348 * t390) * t202 + (-t490 - t120 * t354 + (t213 * t348 + t354 * t390) * qJD(4)) * t278, -t119 * t385 - t203 * t390 + t248 * t383 + t278 * t466, t120 * t385 - t213 * t203 - t248 * t530 - t278 * t472, -t151 * t385 + t203 * t248 (-t217 * t448 + t130) * t248 + t184 * t151 - t400 * t385 + t122 * t203 + t146 * t213 - t510 * t120 + t278 * t163 - g(1) * t243 - g(2) * t245 + ((-qJD(4) * t189 - t145) * t248 - t217 * t151 - t416 * t385 + t97 * t278 + t175 * t202) * t348, -g(1) * t242 - g(2) * t244 - t119 * t510 - t123 * t203 - t146 * t390 - t151 * t456 + t175 * t383 - t248 * t377 + t380 * t385 + t478 * t97, -t182 * t50 - t391 * t63, -t141 * t63 - t181 * t50 - t182 * t367 + t391 * t64, -t149 * t182 - t203 * t391 + t236 * t63 - t385 * t50, -t141 * t203 - t149 * t181 - t236 * t64 - t367 * t385, -t149 * t385 + t203 * t236, t422 * t236 + t415 * t149 - t423 * t385 + t40 * t203 + t93 * t141 - t170 * t367 + t58 * t181 + t134 * t64 - g(1) * t226 - g(2) * t228 + (-t236 * t463 + t385 * t41) * qJD(5), -g(1) * t225 - g(2) * t227 + t134 * t63 - t149 * t463 + t170 * t50 - t58 * t182 - t41 * t203 - t236 * t382 - t385 * t417 - t391 * t93, t11 * t127 - t24 * t393, -t11 * t126 + t127 * t369 - t24 * t71 + t25 * t393, -t11 * t385 + t127 * t147 - t203 * t393 + t232 * t24, -t126 * t147 - t203 * t71 - t232 * t25 - t369 * t385, -t147 * t385 + t203 * t232 (-t346 * t8 + t352 * t6) * t232 + (-t346 * t55 + t352 * t52) * t147 - t437 * t385 + t9 * t203 + t39 * t71 - t128 * t369 + t21 * t126 + t68 * t25 - g(1) * t222 - g(2) * t224 + ((-t346 * t52 - t352 * t55) * t232 + t10 * t385) * qJD(6), -g(1) * t221 - g(2) * t223 - t10 * t203 + t128 * t11 + t21 * t127 + t68 * t24 - t28 * t385 - t39 * t393 + (-(-qJD(6) * t55 + t6) * t232 - t52 * t147 + t3 * t385) * t346 + (-(qJD(6) * t52 + t8) * t232 - t55 * t147 + t424 * t385) * t352; 0, 0, 0, -t350 * t358 * t355, t452 * t358, t442, t441, qJDD(2), -g(3) * t355 + t350 * t371, g(3) * t350 + t355 * t371, t480, t157, t131, t132, t340, t197 * t341 + (t340 * t506 - t341 * t450 + t451 * t542) * pkin(2) + t368, t198 * t341 + (t256 * t451 - t340 * t349 - t341 * t428) * pkin(2) + t359, t59, t27, t62, t61, t482, t329 * t120 - t538 * t354 + t392 * t348 + t406 * t213 + (-t327 * t448 + t550) * t248 + t388, t329 * t119 + t392 * t354 - t348 * t549 - t406 * t390 + (t327 * t449 + t548) * t248 + t389, t20, t7, t45, t46, t483, t412 * t149 - t298 * t367 + (-t269 * t446 + (-qJD(5) * t268 + t544) * t347 + t518) * t236 + t453 * t141 + t364, -t455 * t149 + t236 * t512 + t298 * t50 - t391 * t453 + t362, t5, t1, t18, t19, t484 (t155 * t352 - t156 * t346) * t147 - t220 * t369 + t462 * t71 + (t346 * t404 - t352 * t403) * t232 + t366 -(t155 * t346 + t156 * t352) * t147 + t220 * t11 - t462 * t393 + (t346 * t403 + t352 * t404) * t232 + t363; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t480, t157, t131, t132, t340, t191 * t341 + t368, t190 * t341 + t359, t59, t27, t62, t61, t482, -pkin(3) * t120 - t180 * t248 - t191 * t213 + (t190 * t248 + t397) * t348 + (-t538 - t438) * t354 + t388, -pkin(3) * t119 + t458 * t248 + t191 * t390 + t397 * t354 + (-t549 + t438) * t348 + t389, t20, t7, t45, t46, t483, t411 * t149 - t330 * t367 + (-t302 * t446 + (-qJD(5) * t300 + t543) * t347 + t517) * t236 + t405 * t141 + t364, -t454 * t149 + t236 * t513 + t330 * t50 - t391 * t405 + t362, t5, t1, t18, t19, t484 (t173 * t352 - t174 * t346) * t147 - t229 * t369 + t459 * t71 + (t346 * t399 - t352 * t398) * t232 + t366 -(t173 * t346 + t174 * t352) * t147 + t229 * t11 - t459 * t393 + (t346 * t398 + t352 * t399) * t232 + t363; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t390 * t213, -t213 ^ 2 + t390 ^ 2, t119 + t486, -t120 - t485, t151, -g(1) * t244 + g(2) * t242 + t123 * t248 + t175 * t390 + (t416 + t322) * t348 + t400, g(1) * t245 - g(2) * t243 + t122 * t248 + t175 * t213 + t322 * t354 - t380, -t537, t527, t524, t519, t149, -t420 * t236 + (t141 * t390 + t149 * t353 - t236 * t447) * pkin(4) + t521, t496 * t236 + (-t149 * t347 - t236 * t446 - t390 * t391) * pkin(4) + t522, -t539, t528, t526, t520, t147, t352 * t477 - (t31 * t352 - t32 * t346) * t232 - t121 * t71 + (-t346 * t489 + (-t346 * t353 - t474) * t232 * qJD(5)) * pkin(4) + ((-pkin(4) * t474 - t328 * t346) * t232 - t10) * qJD(6) + t547, t121 * t393 + (-t477 - t3 + (t31 - (-qJD(5) - qJD(6)) * t347 * pkin(4)) * t232) * t346 + (-pkin(4) * t489 + (-pkin(4) * t446 - qJD(6) * t328 + t32) * t232 - t424) * t352 + t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t537, t527, t524, t519, t149, t41 * t236 + t521, t236 * t40 + t522, -t539, t528, t526, t520, t147 -(-t29 * t346 - t494) * t232 + (t147 * t352 - t232 * t445 + t391 * t71) * pkin(5) + t523 (-t232 * t30 - t3) * t346 + (t232 * t29 - t424) * t352 + (-t147 * t346 - t232 * t444 - t391 * t393) * pkin(5) + t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t539, t528, t526, t520, t147, t10 * t232 + t523, t9 * t232 - t346 * t3 - t352 * t424 + t525;];
tau_reg  = t2;
