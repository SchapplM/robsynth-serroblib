% Calculate inertial parameters regressor of coriolis matrix for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:42
% EndTime: 2019-12-05 17:00:56
% DurationCPUTime: 6.45s
% Computational Cost: add. (3639->460), mult. (9637->625), div. (0->0), fcn. (9484->8), ass. (0->350)
t323 = sin(qJ(4));
t318 = t323 ^ 2;
t326 = cos(qJ(4));
t320 = t326 ^ 2;
t561 = -t318 - t320;
t322 = sin(pkin(5));
t325 = sin(qJ(2));
t477 = t322 * t325;
t285 = t326 * t477;
t327 = cos(qJ(3));
t328 = cos(qJ(2));
t476 = t322 * t328;
t425 = t323 * t476;
t210 = t327 * t425 - t285;
t466 = t326 * t328;
t423 = t327 * t466;
t471 = t325 * t323;
t211 = (t423 + t471) * t322;
t519 = t326 / 0.2e1;
t522 = t323 / 0.2e1;
t560 = t210 * t522 + t211 * t519;
t448 = qJD(4) * t327;
t420 = t323 * t448;
t321 = t327 ^ 2;
t324 = sin(qJ(3));
t319 = t324 ^ 2;
t475 = t323 * t319;
t268 = t321 * t323 - t475;
t442 = t268 * qJD(2);
t559 = t442 - t420;
t438 = t324 * qJD(3);
t558 = t326 * t438 + t442;
t515 = pkin(8) * t324;
t385 = -pkin(3) * t327 - t515;
t363 = -pkin(2) + t385;
t258 = t326 * t363;
t473 = t323 * t327;
t433 = pkin(7) * t473;
t212 = -t258 + t433;
t488 = t212 * t327;
t141 = -pkin(7) * t475 - t488;
t507 = cos(pkin(5));
t246 = t324 * t477 - t507 * t327;
t108 = t246 * t323;
t411 = -t108 / 0.2e1;
t424 = t322 * t466;
t247 = t507 * t324 + t327 * t477;
t483 = t247 * t323;
t153 = t424 + t483;
t495 = t153 * t327;
t351 = -t495 / 0.2e1 + t324 * t411;
t51 = (t423 / 0.2e1 + t471 / 0.2e1) * t322 + t351;
t463 = t51 * qJD(1);
t557 = -qJD(2) * t141 - t463;
t482 = t247 * t326;
t154 = -t425 + t482;
t534 = -t154 / 0.2e1;
t382 = t534 + t482 / 0.2e1;
t408 = -t476 / 0.2e1;
t390 = t323 * t408;
t35 = (t390 + t382) * t324;
t465 = t35 * qJD(1);
t467 = t326 * t327;
t432 = pkin(7) * t467;
t213 = t323 * t363 + t432;
t509 = t327 * pkin(8);
t512 = t324 * pkin(3);
t284 = -t509 + t512;
t270 = t323 * t284;
t472 = t324 * t326;
t228 = -pkin(7) * t472 + t270;
t73 = t228 * t327 + (-t213 + 0.2e1 * t432) * t324;
t556 = t73 * qJD(2) + t465;
t493 = t154 * t327;
t529 = t246 / 0.2e1;
t350 = t493 / 0.2e1 + t472 * t529;
t407 = t476 / 0.2e1;
t389 = t323 * t407;
t53 = t327 * t389 - t285 / 0.2e1 + t350;
t462 = t53 * qJD(1);
t380 = t326 * pkin(4) + t323 * qJ(5);
t249 = t380 * t324;
t481 = t249 * t323;
t474 = t323 * t324;
t303 = pkin(4) * t474;
t470 = t326 * qJ(5);
t232 = t303 + (pkin(7) - t470) * t324;
t485 = t232 * t326;
t487 = t213 * t327;
t62 = t487 + (t481 + t485) * t324;
t555 = -qJD(2) * t62 - t462;
t511 = t326 * pkin(7);
t513 = t323 * pkin(3);
t178 = t323 * (-pkin(2) - t515) + (-qJ(5) + t511 - t513) * t327;
t491 = t178 * t327;
t79 = t232 * t472 + t491;
t554 = qJD(2) * t79 + t462;
t63 = -t232 * t474 + t249 * t472 - t488;
t553 = qJD(2) * t63 + t463;
t315 = t324 * qJ(5);
t184 = t228 + t315;
t517 = pkin(4) * t323;
t283 = -t470 + t517;
t233 = (pkin(7) + t283) * t327;
t368 = t232 * t327 + t233 * t324;
t41 = -t178 * t324 + t184 * t327 + t326 * t368;
t552 = -t41 * qJD(2) - t465;
t458 = t327 * t390 + t285 / 0.2e1;
t339 = t350 + t458;
t551 = qJD(2) * t339 + t108 * qJD(3) - t154 * qJD(4);
t531 = -t210 / 0.2e1;
t330 = (t495 / 0.2e1 - t211 / 0.2e1) * t326 + (-t493 / 0.2e1 + t531) * t323;
t550 = qJD(1) * t330;
t388 = t326 * t408;
t535 = t153 / 0.2e1;
t334 = (t388 - t483 / 0.2e1 + t535) * t324;
t549 = qJD(1) * t334;
t548 = qJD(2) * t330;
t547 = qJD(2) * t334;
t546 = qJD(3) * t330;
t545 = qJD(3) * t334;
t435 = t327 * qJD(3);
t426 = t324 * t476;
t337 = t153 * t210 + t154 * t211 + t246 * t426;
t544 = t337 * qJD(1);
t346 = (-t153 * t323 - t154 * t326 + t247) * t246;
t543 = t346 * qJD(1);
t295 = t320 - t318;
t436 = t327 * qJD(2);
t542 = 0.2e1 * t323 * (qJD(4) + t436) * t472 - t295 * t435;
t523 = -t323 / 0.2e1;
t332 = (t153 * t519 + t154 * t523) * t327 + t560;
t445 = t246 * qJD(3);
t541 = qJD(2) * t332 + t561 * t445;
t404 = t473 / 0.2e1;
t405 = -t473 / 0.2e1;
t406 = t474 / 0.2e1;
t521 = -t324 / 0.2e1;
t333 = t153 * t521 + t247 * t406 + t324 * t388 + (t404 + t405) * t246;
t444 = t247 * qJD(3);
t540 = qJD(2) * t333 + t108 * qJD(4) - t326 * t444;
t439 = t324 * qJD(2);
t539 = qJD(3) * t332 + (t210 * t326 - t211 * t323) * t439;
t538 = qJD(2) * t337 + qJD(3) * t346;
t537 = qJD(3) * t333 + qJD(4) * t339 + (t210 * t327 + t319 * t425) * qJD(2);
t536 = -t153 / 0.2e1;
t533 = t178 / 0.2e1;
t532 = -t184 / 0.2e1;
t530 = t213 / 0.2e1;
t528 = -t249 / 0.2e1;
t527 = -t283 / 0.2e1;
t526 = t283 / 0.2e1;
t525 = -t284 / 0.2e1;
t305 = pkin(7) * t474;
t524 = -t305 / 0.2e1;
t520 = t324 / 0.2e1;
t518 = -t327 / 0.2e1;
t514 = t247 * pkin(7);
t510 = t327 * pkin(7);
t494 = t154 * t212;
t86 = -t494 / 0.2e1;
t10 = t86 + t494 / 0.2e1;
t506 = qJD(1) * t10;
t505 = qJD(2) * t10;
t504 = qJD(2) * t35;
t503 = qJD(2) * t51;
t502 = qJD(2) * t53;
t182 = -t258 + (pkin(7) * t323 + pkin(4)) * t327;
t399 = t212 / 0.2e1 - t182 / 0.2e1;
t401 = qJ(5) * t518;
t19 = t213 * t523 + (t401 + t533) * t323 + (pkin(4) * t518 + t399) * t326;
t489 = t19 * qJD(2);
t486 = t232 * t323;
t484 = t246 * t249;
t110 = t246 * t326;
t273 = -pkin(3) - t380;
t480 = t273 * t323;
t479 = t273 * t324;
t28 = ((t178 - t213) * t326 + (t182 - t212) * t323) * t324;
t478 = t28 * qJD(2);
t469 = t326 * t284;
t468 = t326 * t319;
t71 = (t246 * t324 + t247 * t327 - t477) * t476;
t461 = t71 * qJD(1);
t459 = t561 * pkin(8) * t246;
t269 = t321 * t326 - t468;
t456 = qJD(2) * t269;
t455 = qJD(2) * t322;
t453 = qJD(3) * t323;
t452 = qJD(3) * t326;
t451 = qJD(4) * t212;
t450 = qJD(4) * t323;
t449 = qJD(4) * t326;
t447 = qJD(5) * t323;
t446 = qJD(5) * t327;
t397 = t318 / 0.2e1 - t320 / 0.2e1;
t252 = t397 * t324;
t443 = t252 * qJD(4);
t441 = t295 * qJD(4);
t296 = t321 - t319;
t440 = t296 * qJD(2);
t437 = t326 * qJD(5);
t316 = t324 * pkin(4);
t434 = -t316 + t524;
t431 = pkin(2) * t439;
t430 = pkin(2) * t436;
t429 = pkin(8) * t450;
t428 = pkin(8) * t449;
t427 = t509 / 0.2e1;
t422 = t326 * t439;
t419 = t326 * t448;
t418 = t323 * t437;
t417 = t328 * t455;
t300 = t323 * t449;
t416 = t324 * t447;
t299 = t323 * t452;
t415 = t323 * t436;
t414 = t324 * t435;
t413 = t324 * t436;
t410 = -t110 / 0.2e1;
t409 = t479 / 0.2e1;
t403 = -t472 / 0.2e1;
t402 = -t467 / 0.2e1;
t400 = t533 - t213 / 0.2e1;
t398 = t524 - t316 / 0.2e1;
t227 = t305 + t469;
t240 = (-0.1e1 / 0.2e1 + t397) * t324;
t396 = qJD(2) * t240 - t299;
t187 = qJD(2) * t252 - t299;
t280 = t323 * qJD(2) * t468;
t164 = qJD(3) * t252 + t280;
t395 = pkin(3) * t408;
t302 = -qJD(4) + t436;
t393 = t323 * t422;
t392 = t324 * t299;
t391 = t319 * t300;
t383 = 0.2e1 * t392;
t381 = t283 * t520 + t232 / 0.2e1;
t379 = -t273 * t327 + t515;
t348 = t560 * pkin(8);
t331 = t407 * t479 + t348;
t186 = -t227 - t316;
t338 = t186 * t536 + t154 * t532 - t247 * t232 / 0.2e1;
t1 = (t178 * t519 - t233 / 0.2e1 + t182 * t522) * t246 + t331 + t338;
t11 = t178 * t184 + t182 * t186 + t232 * t233;
t378 = -t1 * qJD(1) + t11 * qJD(2);
t16 = -t178 * t212 + t182 * t213 + t232 * t249;
t359 = pkin(4) * t531 + t211 * qJ(5) / 0.2e1;
t3 = -t484 / 0.2e1 + t399 * t154 + t400 * t153 + t359;
t377 = -t3 * qJD(1) + t16 * qJD(2);
t44 = pkin(7) ^ 2 * t324 * t327 - t212 * t227 + t213 * t228;
t353 = t227 * t535 + t228 * t534;
t5 = t324 * t395 + t514 * t521 + t348 + (t213 * t519 + t212 * t522 - t510 / 0.2e1) * t246 + t353;
t376 = -t5 * qJD(1) + t44 * qJD(2);
t21 = -t182 * t467 - t186 * t472 + (t184 * t324 + t491) * t323;
t375 = -t21 * qJD(2) + t550;
t43 = (t227 * t324 - t488) * t326 + (t228 * t324 + t487) * t323;
t374 = -t43 * qJD(2) + t550;
t72 = t212 * t324 + (t227 - 0.2e1 * t305) * t327;
t373 = -t72 * qJD(2) - t549;
t42 = -t182 * t324 + t186 * t327 + t323 * t368;
t372 = t42 * qJD(2) - t549;
t371 = qJD(3) * t35 + qJD(4) * t51;
t370 = t184 * t326 + t186 * t323;
t369 = -t227 * t323 + t228 * t326;
t142 = -pkin(7) * t468 - t487;
t50 = -t350 + t458;
t366 = qJD(1) * t50 + qJD(2) * t142;
t158 = t273 * t326 + t283 * t323;
t254 = -t270 / 0.2e1;
t287 = pkin(8) * t404;
t38 = t287 - t481 / 0.2e1 - t485 / 0.2e1 - t315 + t254 + (t480 / 0.2e1 + (t527 + pkin(7) / 0.2e1) * t326) * t324;
t365 = -qJD(2) * t38 + qJD(3) * t158;
t159 = t283 * t326 - t480;
t358 = t409 + t427;
t347 = t528 + t358;
t354 = t381 * t323;
t40 = t354 + (t525 + t347) * t326 + t434;
t364 = -qJD(2) * t40 + qJD(3) * t159;
t362 = t302 * t324;
t361 = t427 - t512 / 0.2e1;
t360 = qJ(5) * t532 + t186 * pkin(4) / 0.2e1;
t357 = -t470 / 0.2e1 + t517 / 0.2e1;
t253 = t270 / 0.2e1;
t288 = pkin(8) * t405;
t151 = pkin(3) * t406 + t253 + t288;
t356 = pkin(3) * t452 - qJD(2) * t151;
t152 = (t525 + t361) * t326;
t355 = pkin(3) * t453 - qJD(2) * t152;
t58 = t486 / 0.2e1 + (t525 + t358) * t326 + t398;
t352 = qJD(2) * t58 + t273 * t453;
t237 = t326 * t362;
t312 = t320 * t319;
t266 = t318 * t319 - t312;
t180 = -qJD(2) * t266 + t383;
t223 = -qJD(3) * t295 + 0.2e1 * t393;
t349 = qJD(3) * t268 + t324 * t419;
t34 = (t382 + t389) * t324;
t52 = t327 * t388 - t322 * t471 / 0.2e1 + t351;
t344 = qJD(2) * (t211 * t327 + t319 * t424) + qJD(3) * t34 + qJD(4) * t52;
t343 = qJD(2) * t52 + qJD(3) * t110 + qJD(4) * t153;
t342 = -qJD(4) * t266 + t327 * t383;
t12 = (t527 + t357) * t246;
t329 = (-t400 * t323 - t399 * t326) * pkin(8) + t232 * t526 + t249 * t273 / 0.2e1;
t8 = t329 + t360;
t341 = t273 * t283 * qJD(3) - t12 * qJD(1) + t8 * qJD(2);
t340 = qJD(2) * t34 + qJD(4) * t110 + t323 * t444;
t336 = -qJD(4) * t380 + t437;
t279 = t321 + t312;
t335 = qJD(2) * t279 + t392 - t448;
t310 = -t439 / 0.2e1;
t309 = t439 / 0.2e1;
t308 = t438 / 0.2e1;
t301 = t326 * t436;
t282 = t319 * pkin(7) * t476;
t281 = t326 * t416;
t278 = t302 * qJ(5);
t265 = -t301 + t449;
t264 = t302 * t323;
t257 = (t436 - qJD(4) / 0.2e1) * t324;
t255 = t469 / 0.2e1;
t248 = qJD(3) * t318 + t393;
t241 = t318 * t521 + t320 * t520 + t521;
t236 = (t422 + t453) * t327;
t235 = (-t323 * t439 + t452) * t327;
t234 = t323 * t362;
t215 = t320 * t414 - t391;
t214 = t318 * t414 + t391;
t204 = t323 * t438 - t456;
t203 = -t419 + t456;
t202 = t213 * qJD(4);
t193 = t323 * t237;
t181 = -qJD(3) * t269 + t324 * t420;
t163 = -t320 * t413 - t443;
t162 = -t318 * t413 + t443;
t125 = t326 * t361 + t255 + t305;
t124 = t288 + t254 + (t513 / 0.2e1 + t511) * t324;
t115 = -t443 + (t320 * t439 + t299) * t327;
t114 = t443 + (t318 * t439 - t299) * t327;
t59 = pkin(8) * t402 + t273 * t403 - t486 / 0.2e1 - t469 / 0.2e1 + t398;
t39 = t347 * t326 + t255 + t354 - t434;
t37 = t287 + t315 + pkin(7) * t403 + t253 - t381 * t326 + (t409 + t528) * t323;
t20 = -t212 * t326 / 0.2e1 + t178 * t523 + t182 * t519 + pkin(4) * t402 + (t530 + t401) * t323;
t13 = (t357 + t526) * t246;
t9 = t10 * qJD(4);
t7 = t329 - t360;
t6 = t213 * t410 + t212 * t411 + t510 * t529 + (t514 / 0.2e1 + t395) * t324 + t348 - t353;
t4 = t178 * t536 + t86 + t484 / 0.2e1 + t154 * t182 / 0.2e1 + t153 * t530 + t359;
t2 = t178 * t410 + t182 * t411 + t233 * t529 + t331 - t338;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, t538, 0, 0, 0, 0, 0, 0, 0, 0, 0, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t325 * t455, -t417, 0, 0, 0, 0, 0, 0, 0, 0, (-t325 * t436 - t328 * t438) * t322, (t325 * t439 - t328 * t435) * t322, (t319 + t321) * t417, t461 + (t282 + (pkin(7) * t321 * t328 - pkin(2) * t325) * t322) * qJD(2), 0, 0, 0, 0, 0, 0, t537, t344, t539, t544 + (t210 * t212 + t211 * t213 + t282) * qJD(2) + t6 * qJD(3) + t9, 0, 0, 0, 0, 0, 0, t537, t539, -t344, t544 + (t211 * t178 + t210 * t182 + t232 * t426) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) - t339 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324 * t417 - t444, -t327 * t417 + t445, 0, 0, 0, 0, 0, 0, 0, 0, t540, t340, t541, t543 + t6 * qJD(2) + (-pkin(3) * t247 + t459) * qJD(3), 0, 0, 0, 0, 0, 0, t540, t541, -t340, t543 + t2 * qJD(2) + (t247 * t273 + t459) * qJD(3) + t13 * qJD(4) - t108 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t551, t343, 0, t505, 0, 0, 0, 0, 0, 0, t551, 0, -t343, t4 * qJD(2) + t13 * qJD(3) + (-pkin(4) * t154 - qJ(5) * t153) * qJD(4) + t154 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t551; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t461, 0, 0, 0, 0, 0, 0, -qJD(4) * t50 - t545, t371, t546, -qJD(3) * t5 - t544 + t9, 0, 0, 0, 0, 0, 0, qJD(4) * t53 - t545, t546, -t371, -qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t53 - t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t414, t296 * qJD(3), 0, -t414, 0, 0, -pkin(2) * t438, -pkin(2) * t435, 0, 0, t215, -t342, t181, t214, t349, -t414, -qJD(3) * t72 - qJD(4) * t142, qJD(3) * t73 + qJD(4) * t141, -qJD(3) * t43, qJD(3) * t44, t215, t181, t342, -t414, -t349, t214, qJD(3) * t42 + qJD(4) * t62 - t319 * t418, -qJD(3) * t21 - qJD(4) * t28 + t327 * t416, -qJD(3) * t41 - qJD(4) * t63 + qJD(5) * t279, qJD(3) * t11 + qJD(4) * t16 - qJD(5) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t413, t440, t435, -t413, -t438, 0, -pkin(7) * t435 - t431, pkin(7) * t438 - t430, 0, 0, t115, -t542, t204, t114, t558, -t257, (t323 * t385 - t432) * qJD(3) + t125 * qJD(4) + t373, (t326 * t385 + t433) * qJD(3) + t124 * qJD(4) + t556, qJD(3) * t369 + t374, (-pkin(3) * t510 + pkin(8) * t369) * qJD(3) + t376, t115, t204, t542, -t257, -t558, t114, (-t233 * t326 - t323 * t379) * qJD(3) + t39 * qJD(4) + t241 * qJD(5) + t372, qJD(3) * t370 + t20 * qJD(4) + t375, (-t233 * t323 + t326 * t379) * qJD(3) + t37 * qJD(4) + t281 + t552, (pkin(8) * t370 + t233 * t273) * qJD(3) + t7 * qJD(4) + t59 * qJD(5) + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, -t180, t234, t164, t237, t308, qJD(3) * t125 - t202 - t366, qJD(3) * t124 + t451 - t557, 0, t506, -t164, t234, t180, t308, -t237, t164, qJD(3) * t39 - t202 - t555, -t478 + t20 * qJD(3) + (-t324 * t470 + t303) * qJD(4) - t416, qJD(3) * t37 - t446 - t451 - t553, t7 * qJD(3) + (-pkin(4) * t213 - qJ(5) * t212) * qJD(4) + t178 * qJD(5) + t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t241 - t280, t234, t335, qJD(3) * t59 + qJD(4) * t178 - t554; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547, -t504, -t548, qJD(2) * t5 - t543, 0, 0, 0, 0, 0, 0, t547, -t548, t504, qJD(2) * t1 - qJD(4) * t12 - t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, -t440, 0, t413, 0, 0, t431, t430, 0, 0, t163, 0.2e1 * t193, t203, t162, -t559, t257, qJD(4) * t152 - t373, qJD(4) * t151 - t556, -t374, -t376, t163, t203, -0.2e1 * t193, t257, t559, t162, qJD(4) * t40 - qJD(5) * t240 - t372, -qJD(4) * t19 - t327 * t437 - t375, qJD(4) * t38 + t281 - t552, qJD(4) * t8 - qJD(5) * t58 - t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, t441, 0, -t300, 0, 0, -pkin(3) * t450, -pkin(3) * t449, 0, 0, t300, 0, -t441, 0, 0, -t300, -qJD(4) * t159 + t418, 0, -qJD(4) * t158 + qJD(5) * t318, (qJD(4) * t283 - t447) * t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, -t223, t265, t187, t264, t310, -t355 - t428, -t356 + t429, 0, 0, -t187, t265, t223, t310, -t264, t187, -t364 - t428, t336 - t489, -t365 - t429, pkin(8) * t336 + t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t396, t265, t248, -t352 + t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t50, -t503, 0, -t505, 0, 0, 0, 0, 0, 0, -t502, 0, t503, qJD(2) * t3 + qJD(3) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t180, t235, -t164, -t236, t308, -qJD(3) * t152 + t366, -qJD(3) * t151 + t557, 0, -t506, t164, t235, -t180, t308, t236, -t164, -qJD(3) * t40 + t555, qJD(3) * t19 + t478, -qJD(3) * t38 - t446 + t553, -qJ(5) * t446 - qJD(3) * t8 - t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t223, t301, -t187, -t415, t309, t355, t356, 0, 0, t187, t301, -t223, t309, t415, -t187, t364, t489, t365, -t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t302, -t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t502; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t240 + t280, t235, -t335, qJ(5) * t448 + qJD(3) * t58 + t554; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, t301, -t248, t352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t14;
