% Calculate inertial parameters regressor of coriolis matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:34
% EndTime: 2019-12-05 16:56:44
% DurationCPUTime: 5.57s
% Computational Cost: add. (3742->403), mult. (9895->598), div. (0->0), fcn. (9764->8), ass. (0->344)
t347 = cos(qJ(3));
t346 = cos(qJ(4));
t344 = sin(qJ(3));
t381 = -pkin(3) * t347 - pkin(8) * t344;
t370 = -pkin(2) + t381;
t273 = t346 * t370;
t481 = t346 * qJ(5);
t380 = -t344 * t481 + t273;
t343 = sin(qJ(4));
t519 = pkin(7) * t343;
t173 = (-pkin(4) - t519) * t347 + t380;
t485 = t343 * t347;
t435 = pkin(7) * t485;
t183 = -t380 + t435;
t547 = t173 + t183;
t478 = t346 * t347;
t433 = pkin(7) * t478;
t231 = t343 * t370 + t433;
t513 = t347 * pkin(8);
t517 = t344 * pkin(3);
t304 = -t513 + t517;
t288 = t343 * t304;
t484 = t344 * t346;
t434 = pkin(7) * t484;
t243 = t288 - t434;
t103 = t243 * t347 + (-t231 + 0.2e1 * t433) * t344;
t342 = sin(pkin(5));
t348 = cos(qJ(2));
t488 = t342 * t348;
t402 = t488 / 0.2e1;
t511 = cos(pkin(5));
t393 = t511 * t344;
t345 = sin(qJ(2));
t489 = t342 * t345;
t422 = t347 * t489;
t263 = t393 + t422;
t499 = t263 * t346;
t421 = t343 * t488;
t182 = -t421 + t499;
t531 = t182 / 0.2e1;
t57 = (t343 * t402 + t531 - t499 / 0.2e1) * t344;
t472 = t57 * qJD(1);
t546 = t103 * qJD(2) - t472;
t230 = -t273 + t435;
t289 = t346 * t304;
t486 = t343 * t344;
t328 = pkin(7) * t486;
t242 = t328 + t289;
t102 = t230 * t344 + (t242 - 0.2e1 * t328) * t347;
t403 = -t488 / 0.2e1;
t383 = t346 * t403;
t500 = t263 * t343;
t419 = t346 * t488;
t181 = t419 + t500;
t533 = t181 / 0.2e1;
t58 = (t383 + t533 - t500 / 0.2e1) * t344;
t471 = t58 * qJD(1);
t545 = -t102 * qJD(2) - t471;
t339 = t344 ^ 2;
t479 = t346 * t339;
t501 = t231 * t347;
t175 = -pkin(7) * t479 - t501;
t476 = t347 * t348;
t396 = -t476 / 0.2e1;
t482 = t345 * t346;
t352 = (t343 * t396 + t482 / 0.2e1) * t342;
t262 = t344 * t489 - t511 * t347;
t400 = -t484 / 0.2e1;
t386 = t262 * t400;
t520 = -t347 / 0.2e1;
t405 = t182 * t520;
t357 = t405 + t386;
t78 = t352 + t357;
t467 = t78 * qJD(1);
t544 = qJD(2) * t175 + t467;
t518 = t343 * pkin(4);
t418 = pkin(7) + t518;
t290 = t418 * t344;
t340 = t346 ^ 2;
t104 = pkin(4) * t339 * t340 - t183 * t347 - t290 * t486;
t483 = t345 * t343;
t351 = (t346 * t396 - t483 / 0.2e1) * t342;
t401 = -t486 / 0.2e1;
t505 = t181 * t347;
t358 = -t505 / 0.2e1 + t262 * t401;
t79 = t351 - t358;
t466 = t79 * qJD(1);
t543 = qJD(2) * t104 - t466;
t502 = t230 * t347;
t174 = -t339 * t519 - t502;
t542 = -qJD(2) * t174 + t466;
t423 = t343 * t479;
t424 = t290 * t484;
t184 = -qJ(5) * t486 + t231;
t503 = t184 * t347;
t97 = -pkin(4) * t423 - t424 - t503;
t541 = qJD(2) * t97 + t467;
t477 = t347 * qJ(5);
t186 = -t343 * t477 + t243;
t291 = t418 * t347;
t497 = t291 * t346;
t498 = t290 * t346;
t63 = (t186 + t498) * t347 + (-t184 + t497) * t344;
t540 = t63 * qJD(2) - t472;
t516 = t344 * pkin(4);
t176 = -t346 * t477 + t242 + t516;
t61 = -t173 * t344 + t176 * t347 - t290 * t485 - t291 * t486;
t539 = -t61 * qJD(2) - t471;
t229 = (t346 * t476 + t483) * t342;
t228 = (-t343 * t476 + t482) * t342;
t530 = t228 / 0.2e1;
t47 = (t505 / 0.2e1 - t229 / 0.2e1) * t346 + (t405 + t530) * t343;
t473 = t47 * qJD(1);
t62 = (t242 * t344 - t502) * t346 + (t243 * t344 + t501) * t343;
t538 = t62 * qJD(2) - t473;
t31 = (t173 * t347 + t176 * t344) * t346 + (t186 * t344 + t503) * t343;
t537 = -t31 * qJD(2) + t473;
t338 = t343 ^ 2;
t394 = t338 / 0.2e1 - t340 / 0.2e1;
t320 = t340 - t338;
t443 = t344 * qJD(2);
t417 = t346 * t443;
t536 = qJD(3) * t320 - 0.2e1 * t343 * t417;
t535 = t173 / 0.2e1;
t534 = t176 / 0.2e1;
t532 = -t182 / 0.2e1;
t528 = t263 / 0.2e1;
t512 = pkin(8) + qJ(5);
t300 = t512 * t346;
t527 = -t300 / 0.2e1;
t524 = -t343 / 0.2e1;
t523 = -t344 / 0.2e1;
t522 = -t346 / 0.2e1;
t521 = t346 / 0.2e1;
t515 = t344 * pkin(7);
t514 = t347 * pkin(7);
t395 = t183 / 0.2e1 + t535;
t427 = pkin(4) * t520;
t24 = (t427 + t395) * t346;
t510 = qJD(2) * t24;
t34 = t547 * t486;
t509 = qJD(2) * t34;
t507 = t173 * t346;
t506 = t181 * t343;
t504 = t182 * t346;
t130 = t262 * t343;
t299 = t512 * t343;
t496 = t299 * t344;
t495 = t299 * t347;
t494 = t300 * t344;
t493 = t300 * t347;
t335 = -pkin(4) * t346 - pkin(3);
t491 = t335 * t343;
t490 = t335 * t346;
t487 = t343 * t230;
t480 = t346 * t231;
t38 = (t263 - t504 - t506) * t262;
t475 = t38 * qJD(1);
t420 = t344 * t488;
t41 = -t181 * t228 + t182 * t229 + t262 * t420;
t474 = t41 * qJD(1);
t96 = (t262 * t344 + t263 * t347 - t489) * t488;
t465 = t96 * qJD(1);
t319 = t340 + t338;
t341 = t347 ^ 2;
t321 = t341 - t339;
t286 = t321 * t343;
t461 = qJD(2) * t286;
t287 = t341 * t346 - t479;
t460 = qJD(2) * t287;
t459 = qJD(2) * t342;
t458 = qJD(3) * t343;
t457 = qJD(3) * t346;
t456 = qJD(4) * t182;
t455 = qJD(4) * t184;
t454 = qJD(4) * t300;
t453 = qJD(4) * t343;
t337 = qJD(4) * t346;
t452 = qJD(4) * t347;
t451 = qJD(5) * t343;
t448 = t263 * qJD(3);
t267 = t394 * t344;
t447 = t267 * qJD(4);
t284 = t319 * t339;
t446 = t284 * qJD(2);
t445 = t319 * qJD(3);
t444 = t321 * qJD(2);
t442 = t344 * qJD(3);
t441 = t346 * qJD(5);
t440 = t347 * qJD(2);
t439 = t347 * qJD(3);
t438 = t347 * qJD(5);
t271 = t289 / 0.2e1;
t437 = t271 + t328 / 0.2e1;
t436 = t346 * t518;
t432 = pkin(2) * t443;
t431 = pkin(2) * t440;
t430 = pkin(4) * t453;
t429 = pkin(4) * t337;
t428 = t518 / 0.2e1;
t426 = -t514 / 0.2e1;
t425 = t514 / 0.2e1;
t416 = t343 * t457;
t415 = t346 * t442;
t414 = t343 * t452;
t413 = t346 * t452;
t412 = t344 * t441;
t411 = t348 * t459;
t324 = t343 * t337;
t410 = t343 * t440;
t409 = t343 * t438;
t408 = t344 * t439;
t407 = t344 * t440;
t406 = t346 * t438;
t404 = t290 * t524;
t399 = t484 / 0.2e1;
t398 = -t478 / 0.2e1;
t397 = t477 / 0.2e1;
t392 = t319 * t262;
t391 = pkin(4) * t399;
t390 = -qJD(4) + t440;
t389 = t343 * t415;
t388 = t339 * t324;
t385 = t344 * t403;
t384 = t344 * t402;
t382 = -t344 * t436 + t495 / 0.2e1;
t379 = -t494 / 0.2e1 - t173 / 0.2e1;
t349 = (t184 * t522 + t343 * t535 + t291 / 0.2e1) * t262 - t181 * t176 / 0.2e1 + t186 * t531 + t290 * t528;
t350 = t229 * t527 + t299 * t530 + t335 * t385;
t3 = t349 + t350;
t30 = t173 * t176 + t184 * t186 + t290 * t291;
t378 = t3 * qJD(1) + t30 * qJD(2);
t35 = pkin(4) * t424 - t547 * t184;
t6 = t395 * t182 + (t530 + t386) * pkin(4);
t377 = -qJD(1) * t6 + qJD(2) * t35;
t361 = t228 * t524 + t229 * t521;
t355 = t361 * pkin(8);
t364 = t242 * t533 + t243 * t532;
t4 = pkin(3) * t385 - t263 * t515 / 0.2e1 + t355 + (t480 / 0.2e1 + t487 / 0.2e1 + t426) * t262 + t364;
t64 = pkin(7) ^ 2 * t344 * t347 - t230 * t242 + t231 * t243;
t376 = -t4 * qJD(1) + t64 * qJD(2);
t362 = t181 * t521 + t182 * t524;
t70 = (t402 - t362) * t344;
t85 = (t184 * t343 + t507) * t344;
t375 = -qJD(1) * t70 - qJD(2) * t85;
t374 = -t242 * t343 + t243 * t346;
t196 = t299 * t343 + t300 * t346;
t117 = pkin(4) * t491;
t356 = t335 * t400 + t404;
t12 = t395 * t300 + (t534 + t356) * pkin(4);
t373 = -qJD(2) * t12 + qJD(3) * t117;
t256 = t436 - t491;
t65 = t404 + (-t481 / 0.2e1 + t527) * t347 + (-t490 / 0.2e1 + (0.1e1 - t394) * pkin(4)) * t344 + t437;
t372 = -qJD(2) * t65 - qJD(3) * t256;
t268 = t338 * pkin(4) + t490;
t270 = -t288 / 0.2e1;
t72 = t270 + (t515 / 0.2e1 - t290 / 0.2e1) * t346 + (t397 + t335 * t344 / 0.2e1) * t343 + t382;
t371 = -qJD(2) * t72 + qJD(3) * t268;
t369 = t390 * t344;
t368 = t513 / 0.2e1 - t517 / 0.2e1;
t360 = t368 * t343;
t179 = t288 / 0.2e1 - t360;
t367 = pkin(3) * t457 - qJD(2) * t179;
t359 = t368 * t346;
t180 = -t289 / 0.2e1 + t359;
t366 = pkin(3) * t458 - qJD(2) * t180;
t365 = (t184 / 0.2e1 + t496 / 0.2e1) * t346;
t363 = t506 / 0.2e1 + t504 / 0.2e1;
t247 = t346 * t369;
t207 = -qJD(2) * t267 + t416;
t282 = t417 + t458;
t280 = t343 * t443 - t457;
t189 = qJD(2) * t423 + qJD(3) * t267;
t285 = t320 * t339;
t203 = qJD(2) * t285 + 0.2e1 * t389;
t56 = t426 + t365 + (t427 + t379) * t343;
t353 = t393 / 0.2e1 + t422 / 0.2e1;
t68 = t353 - t363;
t354 = -qJD(1) * t68 + qJD(2) * t56 + qJD(3) * t196;
t333 = -t443 / 0.2e1;
t332 = t443 / 0.2e1;
t331 = t442 / 0.2e1;
t325 = t346 * t440;
t303 = t320 * qJD(4);
t298 = t339 * pkin(7) * t488;
t283 = -t325 + t337;
t281 = t390 * t343;
t272 = (t440 - qJD(4) / 0.2e1) * t344;
t264 = t282 * pkin(4);
t246 = t282 * t347;
t245 = t280 * t347;
t244 = t343 * t369;
t233 = t340 * t408 - t388;
t232 = t338 * t408 + t388;
t222 = t343 * t442 - t460;
t221 = -t413 + t460;
t220 = t414 - t461;
t219 = t415 + t461;
t205 = -qJD(3) * t287 + t344 * t414;
t204 = qJD(3) * t286 + t344 * t413;
t199 = 0.2e1 * t343 * t247;
t188 = -t340 * t407 - t447;
t187 = -t338 * t407 + t447;
t185 = -qJD(4) * t285 - 0.2e1 * t347 * t389;
t154 = t328 + t271 + t359;
t153 = t270 - t360 + t434;
t146 = -t447 + (t340 * t443 + t416) * t347;
t145 = t447 + (t338 * t443 - t416) * t347;
t142 = -0.2e1 * t344 * t324 + t347 * t536;
t132 = t262 * t346;
t81 = t352 - t357;
t80 = t351 + t358;
t73 = t335 * t401 + t498 / 0.2e1 + pkin(7) * t399 + t270 + t343 * t397 - t382;
t71 = t181 * t399 + t182 * t401 + t384;
t69 = t353 + t363;
t66 = t493 / 0.2e1 + qJ(5) * t398 - t356 + t437 + (0.1e1 + t394) * t516;
t60 = t182 * t523 + t263 * t399 + t343 * t384;
t59 = t181 * t523 + t344 * t383 + t486 * t528;
t55 = t343 * t379 + t347 * t428 + t365 + t425;
t46 = t347 * t362 + t361;
t45 = qJD(2) * t79;
t44 = qJD(2) * t78;
t29 = qJD(2) * t81 + qJD(3) * t130 - t456;
t28 = qJD(2) * t80 + qJD(3) * t132 + qJD(4) * t181;
t27 = qJD(2) * t58;
t26 = qJD(2) * t57;
t25 = t183 * t522 - t507 / 0.2e1 + pkin(4) * t398;
t23 = qJD(2) * t60 + qJD(4) * t132 + t343 * t448;
t22 = qJD(2) * t59 + qJD(4) * t130 - t346 * t448;
t21 = qJD(2) * t47;
t20 = -qJD(3) * t57 - qJD(4) * t79;
t19 = -qJD(3) * t58 - qJD(4) * t78;
t18 = qJD(3) * t47;
t17 = pkin(4) * t130;
t13 = pkin(4) * t534 + t290 * t428 + t335 * t391 + t547 * t527;
t11 = t46 * qJD(2) - qJD(3) * t392;
t10 = (t229 * t347 + t339 * t419) * qJD(2) + t60 * qJD(3) + t80 * qJD(4);
t9 = (-t228 * t347 + t339 * t421) * qJD(2) + t59 * qJD(3) + t81 * qJD(4);
t8 = t46 * qJD(3) + (-t228 * t346 - t229 * t343) * t443;
t7 = pkin(4) * t530 + t262 * t391 + t547 * t532;
t5 = t262 * t425 + (pkin(3) * t403 + pkin(7) * t528) * t344 + t355 - t364 - (t480 + t487) * t262 / 0.2e1;
t2 = t349 - t350;
t1 = qJD(2) * t41 + qJD(3) * t38;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345 * t459, -t411, 0, 0, 0, 0, 0, 0, 0, 0, (-t345 * t440 - t348 * t442) * t342, (t345 * t443 - t348 * t439) * t342, (t339 + t341) * t411, t465 + (t298 + (pkin(7) * t341 * t348 - pkin(2) * t345) * t342) * qJD(2), 0, 0, 0, 0, 0, 0, t9, t10, t8, t474 + (-t228 * t230 + t229 * t231 + t298) * qJD(2) + t5 * qJD(3), 0, 0, 0, 0, 0, 0, t9, t10, t8, t474 + (t228 * t173 + t229 * t184 + t290 * t420) * qJD(2) + t2 * qJD(3) + t7 * qJD(4) + t71 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t344 * t411 - t448, t262 * qJD(3) - t347 * t411, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, t11, t475 + t5 * qJD(2) + (-pkin(3) * t263 - pkin(8) * t392) * qJD(3), 0, 0, 0, 0, 0, 0, t22, t23, t11, t475 + t2 * qJD(2) + (-t196 * t262 + t263 * t335) * qJD(3) + t17 * qJD(4) + t69 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t28, 0, 0, 0, 0, 0, 0, 0, 0, t29, t28, 0, -pkin(4) * t456 + qJD(2) * t7 + qJD(3) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t71 + qJD(3) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t465, 0, 0, 0, 0, 0, 0, t19, t20, t18, -qJD(3) * t4 - t474, 0, 0, 0, 0, 0, 0, t19, t20, t18, qJD(3) * t3 - qJD(4) * t6 - qJD(5) * t70 - t474; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t408, t321 * qJD(3), 0, -t408, 0, 0, -pkin(2) * t442, -pkin(2) * t439, 0, 0, t233, t185, t205, t232, t204, -t408, -qJD(3) * t102 - qJD(4) * t175, qJD(3) * t103 + qJD(4) * t174, -qJD(3) * t62, qJD(3) * t64, t233, t185, t205, t232, t204, -t408, -qJD(3) * t61 - qJD(4) * t97 + t344 * t406, qJD(3) * t63 + qJD(4) * t104 - t344 * t409, -qJD(3) * t31 + qJD(4) * t34 + qJD(5) * t284, qJD(3) * t30 + qJD(4) * t35 - qJD(5) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t407, t444, t439, -t407, -t442, 0, -pkin(7) * t439 - t432, pkin(7) * t442 - t431, 0, 0, t146, t142, t222, t145, t219, -t272, (t343 * t381 - t433) * qJD(3) + t154 * qJD(4) + t545, (t346 * t381 + t435) * qJD(3) + t153 * qJD(4) + t546, qJD(3) * t374 - t538, (-pkin(3) * t514 + pkin(8) * t374) * qJD(3) + t376, t146, t142, t222, t145, t219, -t272, (t335 * t485 - t496 - t497) * qJD(3) + t66 * qJD(4) + t409 + t539, (t291 * t343 + t335 * t478 - t494) * qJD(3) + t73 * qJD(4) + t406 + t540, ((t186 + t495) * t346 + (-t176 - t493) * t343) * qJD(3) + t25 * qJD(4) + t537, (-t176 * t299 + t186 * t300 + t291 * t335) * qJD(3) + t13 * qJD(4) + t55 * qJD(5) + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t203, t244, t189, t247, t331, qJD(3) * t154 - qJD(4) * t231 - t544, qJD(3) * t153 + qJD(4) * t230 - t542, 0, 0, -t189, -t203, t244, t189, t247, t331, qJD(3) * t66 - t455 - t541, qJD(3) * t73 + qJD(4) * t183 + t543, qJD(3) * t25 + t344 * t430 + t509, -pkin(4) * t455 + qJD(3) * t13 + t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t245, t446, qJD(3) * t55 + t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, -t21, qJD(2) * t4 - t475, 0, 0, 0, 0, 0, 0, t27, t26, -t21, -qJD(2) * t3 - qJD(5) * t68 - t475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t407, -t444, 0, t407, 0, 0, t432, t431, 0, 0, t188, t199, t221, t187, t220, t272, qJD(4) * t180 - t545, qJD(4) * t179 - t546, t538, -t376, t188, t199, t221, t187, t220, t272, -qJD(4) * t65 - t539, -qJD(4) * t72 - t540, -qJD(4) * t24 - t537, -qJD(4) * t12 + qJD(5) * t56 - t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t303, 0, -t324, 0, 0, -pkin(3) * t453, -pkin(3) * t337, 0, 0, t324, t303, 0, -t324, 0, 0, -t256 * qJD(4), t268 * qJD(4), qJD(5) * t319, qJD(4) * t117 + qJD(5) * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t536, t283, -t207, t281, t333, -pkin(8) * t337 - t366, pkin(8) * t453 - t367, 0, 0, t207, t536, t283, -t207, t281, t333, t372 - t454, qJD(4) * t299 + t371, -t429 - t510, -pkin(4) * t454 + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t445, t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, qJD(2) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t203, -t245, -t189, -t246, t331, -qJD(3) * t180 + t544, -qJD(3) * t179 + t542, 0, 0, t189, t203, -t245, -t189, -t246, t331, qJD(3) * t65 - t412 + t541, qJD(3) * t72 + t344 * t451 - t543, qJD(3) * t24 - t509, -pkin(4) * t412 + qJD(3) * t12 - t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, -t536, t325, t207, -t410, t332, t366, t367, 0, 0, -t207, -t536, t325, t207, -t410, t332, -t372 - t451, -t371 - t441, t510, -pkin(4) * t451 - t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, t280, 0, -t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t70 + qJD(3) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, t244, -t446, -qJD(3) * t56 + t344 * t429 - t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t453, t337, -t445, -t354 + t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, -t280, 0, t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t14;
