% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:34
% EndTime: 2019-03-08 23:03:41
% DurationCPUTime: 4.05s
% Computational Cost: add. (5200->370), mult. (13500->528), div. (0->0), fcn. (10639->12), ass. (0->198)
t449 = sin(qJ(4));
t450 = sin(qJ(3));
t453 = cos(qJ(4));
t454 = cos(qJ(3));
t412 = t449 * t454 + t450 * t453;
t563 = pkin(8) + pkin(9);
t519 = qJD(3) * t563;
t413 = t450 * t519;
t414 = t454 * t519;
t418 = t563 * t450;
t419 = t563 * t454;
t478 = t418 * t449 - t419 * t453;
t455 = cos(qJ(2));
t446 = sin(pkin(6));
t531 = qJD(1) * t446;
t517 = t455 * t531;
t577 = qJD(4) * t478 + t412 * t517 + t413 * t449 - t453 * t414;
t411 = t449 * t450 - t453 * t454;
t525 = qJD(4) * t449;
t547 = t418 * t453;
t576 = qJD(4) * t547 - t411 * t517 + t453 * t413 + t449 * t414 + t419 * t525;
t442 = qJD(3) + qJD(4);
t452 = cos(qJ(6));
t448 = sin(qJ(6));
t528 = qJD(2) * t450;
t513 = t449 * t528;
t526 = qJD(2) * t454;
t404 = -t453 * t526 + t513;
t512 = t453 * t528;
t406 = -t449 * t526 - t512;
t445 = sin(pkin(12));
t558 = cos(pkin(12));
t468 = -t445 * t404 - t406 * t558;
t549 = t468 * t448;
t355 = -t452 * t442 + t549;
t500 = -t558 * t404 + t406 * t445;
t565 = qJD(6) - t500;
t575 = t355 * t565;
t357 = t442 * t448 + t452 * t468;
t574 = t357 * t565;
t383 = t442 * t412;
t573 = -qJ(5) * t383 - qJD(5) * t411 - t576;
t382 = t442 * t411;
t572 = qJ(5) * t382 - qJD(5) * t412 + t577;
t521 = qJD(2) * qJD(3);
t509 = t454 * t521;
t571 = qJD(4) * t526 + t509;
t451 = sin(qJ(2));
t518 = t451 * t531;
t559 = qJD(3) * pkin(3);
t466 = t450 * t559 - t518;
t499 = t565 * t452;
t370 = -t442 * t513 + t453 * t571;
t510 = t450 * t521;
t494 = qJD(4) * t512 + t449 * t571 + t453 * t510;
t345 = t370 * t445 + t494 * t558;
t541 = t448 * t345;
t570 = -t565 * t499 - t541;
t569 = MDP(6) * (t450 ^ 2 - t454 ^ 2);
t504 = qJD(2) * t563 + t518;
t447 = cos(pkin(6));
t530 = qJD(1) * t447;
t386 = -t504 * t450 + t454 * t530;
t529 = qJD(2) * t446;
t511 = qJD(1) * t529;
t491 = t455 * t511;
t358 = qJD(3) * t386 + t454 * t491;
t380 = t386 + t559;
t567 = (qJD(4) * t380 + t358) * t453;
t398 = t406 * qJ(5);
t387 = t450 * t530 + t454 * t504;
t377 = t449 * t387;
t502 = t453 * t380 - t377;
t339 = t398 + t502;
t566 = pkin(4) * t383 + t466;
t359 = -qJD(3) * t387 - t450 * t491;
t503 = t449 * t359 - t387 * t525;
t301 = -qJ(5) * t494 - t404 * qJD(5) + t503 + t567;
t379 = t453 * t387;
t480 = -t380 * t449 - t379;
t481 = -t449 * t358 + t453 * t359;
t462 = qJD(4) * t480 + t481;
t458 = -qJ(5) * t370 + qJD(5) * t406 + t462;
t294 = t301 * t558 + t445 * t458;
t330 = pkin(4) * t442 + t339;
t557 = qJ(5) * t404;
t340 = -t480 - t557;
t546 = t445 * t340;
t312 = t330 * t558 - t546;
t310 = -t442 * pkin(5) - t312;
t439 = -pkin(3) * t454 - pkin(2);
t399 = qJD(2) * t439 - t517;
t371 = pkin(4) * t404 + qJD(5) + t399;
t321 = -pkin(5) * t500 - pkin(10) * t468 + t371;
t375 = t411 * t558 + t412 * t445;
t376 = -t445 * t411 + t412 * t558;
t477 = pkin(4) * t411 + t439;
t331 = pkin(5) * t375 - pkin(10) * t376 + t477;
t348 = -t382 * t558 - t445 * t383;
t293 = t301 * t445 - t458 * t558;
t363 = -qJ(5) * t411 - t478;
t465 = -qJ(5) * t412 - t419 * t449 - t547;
t335 = t363 * t558 + t445 * t465;
t486 = t293 * t376 - t335 * t345;
t537 = t445 * t572 + t558 * t573;
t564 = -(qJD(6) * t321 + t294) * t375 + t310 * t348 + (-qJD(6) * t331 - t537) * t565 + t486;
t562 = pkin(4) * t406;
t561 = pkin(3) * qJD(4);
t560 = qJD(2) * pkin(2);
t556 = t310 * t500;
t555 = t310 * t376;
t524 = qJD(6) * t448;
t346 = t370 * t558 - t445 * t494;
t523 = qJD(6) * t452;
t536 = t452 * t346 + t442 * t523;
t319 = -t468 * t524 + t536;
t554 = t319 * t448;
t553 = t331 * t345;
t552 = t346 * t448;
t551 = t355 * t468;
t550 = t357 * t468;
t548 = t399 * t406;
t545 = t445 * t449;
t544 = t446 * t451;
t543 = t446 * t455;
t457 = qJD(2) ^ 2;
t542 = t446 * t457;
t456 = qJD(3) ^ 2;
t540 = t450 * t456;
t341 = t452 * t345;
t539 = t454 * t456;
t538 = t445 * t573 - t558 * t572;
t332 = t558 * t340;
t313 = t445 * t330 + t332;
t535 = t453 * t386 - t377;
t344 = t398 + t535;
t501 = -t386 * t449 - t379;
t473 = t501 + t557;
t505 = t558 * t449;
t534 = -t344 * t445 + t558 * t473 + (t445 * t453 + t505) * t561;
t533 = -t558 * t344 - t445 * t473 + (t453 * t558 - t545) * t561;
t400 = pkin(3) * t510 + t451 * t511;
t438 = pkin(3) * t453 + pkin(4);
t397 = pkin(3) * t505 + t445 * t438;
t527 = qJD(2) * t451;
t440 = pkin(3) * t528;
t520 = t451 * t542;
t515 = t446 * t527;
t514 = t455 * t529;
t507 = -pkin(3) * t442 - t380;
t328 = pkin(5) * t468 - pkin(10) * t500 - t562;
t393 = pkin(10) + t397;
t496 = qJD(6) * t393 + t328 + t440;
t435 = pkin(4) * t445 + pkin(10);
t495 = qJD(6) * t435 + t328;
t493 = t450 * t514;
t492 = t454 * t514;
t311 = pkin(10) * t442 + t313;
t296 = t311 * t452 + t321 * t448;
t490 = t293 * t448 + t296 * t468 + t310 * t523;
t347 = -t382 * t445 + t383 * t558;
t489 = pkin(5) * t347 - pkin(10) * t348 + t566;
t485 = -t345 * t393 - t556;
t484 = -t345 * t435 - t556;
t483 = t311 * t448 - t321 * t452;
t482 = t312 * t500 + t313 * t468;
t402 = t447 * t454 - t450 * t544;
t403 = t447 * t450 + t454 * t544;
t479 = t402 * t453 - t403 * t449;
t365 = t402 * t449 + t403 * t453;
t476 = t341 + (t448 * t500 - t524) * t565;
t337 = t365 * t558 + t445 * t479;
t475 = -t337 * t448 - t452 * t543;
t474 = -t337 * t452 + t448 * t543;
t472 = -t293 * t452 + t310 * t524 + t468 * t483;
t471 = t399 * t404 - t503;
t470 = t348 * t452 - t376 * t524;
t469 = qJD(2) * t560;
t396 = -pkin(3) * t545 + t438 * t558;
t353 = pkin(4) * t494 + t400;
t464 = -0.2e1 * qJD(3) * t560;
t320 = qJD(6) * t357 + t552;
t461 = -t406 * t404 * MDP(12) - t565 * t468 * MDP(25) + ((t319 - t575) * t452 + (-t320 - t574) * t448) * MDP(22) + (t476 + t551) * MDP(24) + (-t550 - t570) * MDP(23) + (t357 * t499 + t554) * MDP(21) + (t404 * t442 + t370) * MDP(14) + (-t406 * t442 - t494) * MDP(15) + (-t404 ^ 2 + t406 ^ 2) * MDP(13);
t436 = -pkin(4) * t558 - pkin(5);
t392 = -pkin(5) - t396;
t385 = -qJD(3) * t403 - t493;
t384 = qJD(3) * t402 + t492;
t336 = t365 * t445 - t479 * t558;
t334 = t363 * t445 - t465 * t558;
t325 = -qJD(4) * t365 - t384 * t449 + t385 * t453;
t324 = qJD(4) * t479 + t384 * t453 + t385 * t449;
t315 = t339 * t558 - t546;
t314 = t339 * t445 + t332;
t307 = t345 * pkin(5) - t346 * pkin(10) + t353;
t304 = t452 * t307;
t303 = t324 * t558 + t445 * t325;
t302 = t324 * t445 - t325 * t558;
t1 = [-MDP(3) * t520 - t455 * MDP(4) * t542 + (-t454 * t520 + (t385 - t493) * qJD(3)) * MDP(10) + (t450 * t520 + (-t384 - t492) * qJD(3)) * MDP(11) + (t325 * t442 + (t404 * t527 - t455 * t494) * t446) * MDP(17) + (-t324 * t442 + (-t370 * t455 - t406 * t527) * t446) * MDP(18) + (t302 * t468 + t303 * t500 + t336 * t346 - t337 * t345) * MDP(19) + (t293 * t336 + t294 * t337 - t302 * t312 + t303 * t313 + (-t353 * t455 + t371 * t527) * t446) * MDP(20) + ((qJD(6) * t474 - t303 * t448 + t452 * t515) * t565 + t475 * t345 + t302 * t355 + t336 * t320) * MDP(26) + (-(qJD(6) * t475 + t303 * t452 + t448 * t515) * t565 + t474 * t345 + t302 * t357 + t336 * t319) * MDP(27); 0.2e1 * t450 * MDP(5) * t509 - 0.2e1 * t521 * t569 + MDP(7) * t539 - MDP(8) * t540 + (-pkin(8) * t539 + t450 * t464) * MDP(10) + (pkin(8) * t540 + t454 * t464) * MDP(11) + (t370 * t412 + t382 * t406) * MDP(12) + (-t370 * t411 + t382 * t404 + t406 * t383 - t412 * t494) * MDP(13) + (t399 * t383 + t400 * t411 + t466 * t404 + t439 * t494) * MDP(17) + (t439 * t370 - t399 * t382 + t400 * t412 - t466 * t406) * MDP(18) + (-t294 * t375 - t312 * t348 - t313 * t347 + t334 * t346 + t468 * t538 + t500 * t537 + t486) * MDP(19) + (t293 * t334 + t294 * t335 - t538 * t312 + t537 * t313 + t353 * t477 + t371 * t566) * MDP(20) + (t319 * t376 * t452 + t357 * t470) * MDP(21) + ((-t355 * t452 - t357 * t448) * t348 + (-t554 - t320 * t452 + (t355 * t448 - t357 * t452) * qJD(6)) * t376) * MDP(22) + (t319 * t375 + t341 * t376 + t347 * t357 + t470 * t565) * MDP(23) + (-t376 * t541 - t320 * t375 - t347 * t355 + (-t348 * t448 - t376 * t523) * t565) * MDP(24) + (t345 * t375 + t347 * t565) * MDP(25) + (-t483 * t347 + t304 * t375 + t334 * t320 + t538 * t355 + (t553 + t489 * t565 + (-t311 * t375 - t335 * t565 + t555) * qJD(6)) * t452 + t564 * t448) * MDP(26) + (-t296 * t347 + t334 * t319 + t538 * t357 + (-t553 - (-qJD(6) * t311 + t307) * t375 - qJD(6) * t555 + (qJD(6) * t335 - t489) * t565) * t448 + t564 * t452) * MDP(27) + (-t382 * MDP(14) - t383 * MDP(15) + MDP(17) * t577 + t576 * MDP(18)) * t442; (-t345 * t397 - t346 * t396 + t468 * t534 + t500 * t533 + t482) * MDP(19) + (t294 * t397 - t293 * t396 - t371 * (t440 - t562) + t533 * t313 - t534 * t312) * MDP(20) + (t392 * t320 + t485 * t448 + t534 * t355 + (-t448 * t533 - t452 * t496) * t565 + t472) * MDP(26) + t454 * t469 * MDP(11) + (t392 * t319 + t485 * t452 + t534 * t357 + (t448 * t496 - t452 * t533) * t565 + t490) * MDP(27) + t457 * t569 + (t535 * t442 + t406 * t440 + (qJD(4) * t507 - t358) * t453 + t471) * MDP(18) + (-t501 * t442 - t404 * t440 + t548 + (t449 * t507 - t379) * qJD(4) + t481) * MDP(17) + t461 + (-t457 * t454 * MDP(5) + MDP(10) * t469) * t450; (-t314 * t468 - t315 * t500 + (-t345 * t445 - t346 * t558) * pkin(4) + t482) * MDP(19) + (-t442 * t480 + t462 + t548) * MDP(17) + (t442 * t502 + t471 - t567) * MDP(18) + (t312 * t314 - t313 * t315 + (-t293 * t558 + t294 * t445 + t371 * t406) * pkin(4)) * MDP(20) + (-t314 * t355 + t436 * t320 + t484 * t448 + (t315 * t448 - t452 * t495) * t565 + t472) * MDP(26) + (-t314 * t357 + t436 * t319 + t484 * t452 + (t315 * t452 + t448 * t495) * t565 + t490) * MDP(27) + t461; (-t468 ^ 2 - t500 ^ 2) * MDP(19) + (t312 * t468 - t313 * t500 + t353) * MDP(20) + (t476 - t551) * MDP(26) + (-t550 + t570) * MDP(27); t357 * t355 * MDP(21) + (-t355 ^ 2 + t357 ^ 2) * MDP(22) + (t536 + t575) * MDP(23) + (-t552 + t574) * MDP(24) + t345 * MDP(25) + (-t294 * t448 + t296 * t565 - t310 * t357 + t304) * MDP(26) + (-t294 * t452 - t307 * t448 + t310 * t355 - t483 * t565) * MDP(27) + (-MDP(23) * t549 - MDP(24) * t357 - MDP(26) * t296 + MDP(27) * t483) * qJD(6);];
tauc  = t1;
