% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:25
% EndTime: 2019-03-08 22:00:36
% DurationCPUTime: 5.52s
% Computational Cost: add. (4006->421), mult. (10350->590), div. (0->0), fcn. (8216->12), ass. (0->204)
t438 = sin(pkin(12));
t444 = sin(qJ(3));
t505 = qJD(2) * t444;
t440 = cos(pkin(12));
t448 = cos(qJ(3));
t523 = t440 * t448;
t407 = qJD(2) * t523 - t438 * t505;
t549 = qJD(5) + qJD(6);
t564 = t407 - t549;
t445 = sin(qJ(2));
t439 = sin(pkin(6));
t508 = qJD(1) * t439;
t492 = t445 * t508;
t504 = qJD(3) * t444;
t563 = pkin(3) * t504 - t492;
t546 = -qJ(4) - pkin(8);
t483 = qJD(3) * t546;
t404 = qJD(4) * t448 + t444 * t483;
t405 = -qJD(4) * t444 + t448 * t483;
t351 = t404 * t440 + t405 * t438;
t415 = t438 * t444 - t523;
t449 = cos(qJ(2));
t491 = t449 * t508;
t384 = t415 * t491;
t514 = t351 + t384;
t416 = t438 * t448 + t440 * t444;
t408 = t416 * qJD(3);
t411 = t415 * qJD(3);
t562 = pkin(4) * t408 + pkin(9) * t411 + t563;
t420 = qJD(2) * pkin(8) + t492;
t477 = qJ(4) * qJD(2) + t420;
t441 = cos(pkin(6));
t507 = qJD(1) * t441;
t490 = t444 * t507;
t380 = t448 * t477 + t490;
t367 = t438 * t380;
t428 = t448 * t507;
t379 = -t444 * t477 + t428;
t372 = qJD(3) * pkin(3) + t379;
t315 = t372 * t440 - t367;
t312 = -qJD(3) * pkin(4) - t315;
t409 = t416 * qJD(2);
t443 = sin(qJ(5));
t447 = cos(qJ(5));
t498 = t447 * qJD(3);
t387 = t409 * t443 - t498;
t308 = pkin(5) * t387 + t312;
t446 = cos(qJ(6));
t389 = qJD(3) * t443 + t409 * t447;
t442 = sin(qJ(6));
t535 = t389 * t442;
t326 = t387 * t446 + t535;
t561 = t308 * t326;
t401 = qJD(5) - t407;
t397 = qJD(6) + t401;
t560 = t326 * t397;
t465 = t387 * t442 - t389 * t446;
t559 = t397 * t465;
t522 = t442 * t447;
t419 = t443 * t446 + t522;
t511 = t564 * t419;
t503 = qJD(5) * t443;
t531 = t407 * t443;
t558 = t503 - t531;
t496 = qJD(2) * qJD(3);
t484 = t448 * t496;
t485 = t444 * t496;
t399 = -t438 * t485 + t440 * t484;
t339 = qJD(5) * t498 + t399 * t447 - t409 * t503;
t398 = qJD(2) * t408;
t524 = t440 * t380;
t316 = t372 * t438 + t524;
t313 = qJD(3) * pkin(9) + t316;
t493 = -pkin(3) * t448 - pkin(2);
t400 = qJD(2) * t493 + qJD(4) - t491;
t332 = -pkin(4) * t407 - pkin(9) * t409 + t400;
t301 = t313 * t447 + t332 * t443;
t470 = qJD(4) + t491;
t338 = (-t420 * t444 + t428) * qJD(3) + (-qJ(4) * t504 + t448 * t470) * qJD(2);
t551 = (-t420 * t448 - t490) * qJD(3) + (-qJ(4) * qJD(3) * t448 - t444 * t470) * qJD(2);
t306 = t440 * t338 + t438 * t551;
t506 = qJD(2) * t439;
t489 = t445 * t506;
t406 = pkin(3) * t485 + qJD(1) * t489;
t324 = pkin(4) * t398 - pkin(9) * t399 + t406;
t322 = t447 * t324;
t453 = -qJD(5) * t301 - t306 * t443 + t322;
t283 = pkin(5) * t398 - pkin(10) * t339 + t453;
t340 = qJD(5) * t389 + t399 * t443;
t502 = qJD(5) * t447;
t458 = t306 * t447 - t313 * t503 + t324 * t443 + t332 * t502;
t284 = -pkin(10) * t340 + t458;
t481 = t283 * t446 - t442 * t284;
t557 = t308 * t465 + t481;
t556 = t398 * MDP(25) + (-t326 ^ 2 + t465 ^ 2) * MDP(22) - t326 * MDP(21) * t465;
t555 = MDP(5) * t444;
t554 = MDP(6) * (t444 ^ 2 - t448 ^ 2);
t359 = t419 * t416;
t553 = -t384 * t443 + t447 * t562;
t515 = t404 * t438 - t405 * t440 - t416 * t491;
t365 = pkin(4) * t415 - pkin(9) * t416 + t493;
t423 = t546 * t444;
t424 = t546 * t448;
t386 = t423 * t438 - t424 * t440;
t552 = -t365 * t502 + t386 * t503 - t443 * t562 - t447 * t514;
t550 = MDP(10) * t444 + MDP(11) * t448;
t418 = t442 * t443 - t446 * t447;
t512 = t564 * t418;
t548 = -t397 * t512 - t398 * t419;
t479 = t339 * t442 + t340 * t446;
t291 = -qJD(6) * t465 + t479;
t431 = pkin(3) * t438 + pkin(9);
t547 = pkin(10) + t431;
t545 = qJD(2) * pkin(2);
t300 = -t313 * t443 + t332 * t447;
t294 = -pkin(10) * t389 + t300;
t289 = pkin(5) * t401 + t294;
t544 = t289 * t446;
t295 = -pkin(10) * t387 + t301;
t543 = t295 * t446;
t542 = t326 * t409;
t541 = t465 * t409;
t540 = t339 * t443;
t539 = t387 * t401;
t538 = t387 * t409;
t537 = t389 * t401;
t536 = t389 * t409;
t533 = t398 * t443;
t530 = t411 * t443;
t529 = t411 * t447;
t528 = t416 * t443;
t527 = t416 * t447;
t526 = t439 * t445;
t525 = t439 * t449;
t450 = qJD(3) ^ 2;
t521 = t444 * t450;
t373 = t447 * t386;
t390 = t447 * t398;
t520 = t448 * t450;
t519 = pkin(10) * t529 + pkin(5) * t408 - t351 * t443 + (-t373 + (pkin(10) * t416 - t365) * t443) * qJD(5) + t553;
t486 = t416 * t502;
t461 = t486 - t530;
t518 = pkin(10) * t461 + t552;
t517 = pkin(5) * t461 + t515;
t320 = t379 * t440 - t367;
t352 = pkin(3) * t505 + pkin(4) * t409 - pkin(9) * t407;
t516 = t320 * t447 + t352 * t443;
t513 = t365 * t443 + t373;
t501 = qJD(6) * t442;
t500 = qJD(6) * t446;
t494 = t339 * t446 - t340 * t442 - t387 * t500;
t432 = -pkin(3) * t440 - pkin(4);
t488 = t449 * t506;
t487 = t416 * t503;
t482 = qJD(5) * t547;
t292 = t295 * t501;
t480 = t442 * t283 - t292;
t305 = t338 * t438 - t440 * t551;
t318 = t379 * t438 + t524;
t385 = -t423 * t440 - t424 * t438;
t478 = t401 * t447;
t476 = qJD(6) * t289 + t284;
t475 = pkin(5) * t558 - t318;
t474 = t397 * t511 - t418 * t398;
t345 = t447 * t352;
t414 = t547 * t447;
t472 = pkin(5) * t409 + qJD(6) * t414 - t320 * t443 + t345 + (-pkin(10) * t407 + t482) * t447;
t413 = t547 * t443;
t471 = -pkin(10) * t531 + qJD(6) * t413 + t443 * t482 + t516;
t286 = t289 * t442 + t543;
t469 = t305 * t416 - t386 * t398;
t362 = t447 * t365;
t307 = pkin(5) * t415 - pkin(10) * t527 - t386 * t443 + t362;
t309 = -pkin(10) * t528 + t513;
t468 = t307 * t442 + t309 * t446;
t412 = t441 * t444 + t448 * t526;
t462 = t441 * t448 - t444 * t526;
t358 = t412 * t440 + t438 * t462;
t336 = -t358 * t443 - t447 * t525;
t463 = -t358 * t447 + t443 * t525;
t467 = t336 * t446 + t442 * t463;
t466 = t336 * t442 - t446 * t463;
t464 = -t401 * t558 + t390;
t460 = -t487 - t529;
t290 = -t389 * t501 + t494;
t456 = t312 * t401 - t398 * t431;
t454 = -0.2e1 * qJD(3) * t545;
t451 = qJD(2) ^ 2;
t422 = -pkin(5) * t447 + t432;
t378 = -qJD(3) * t412 - t444 * t488;
t377 = qJD(3) * t462 + t448 * t488;
t366 = t398 * t415;
t360 = t418 * t416;
t357 = t412 * t438 - t440 * t462;
t349 = pkin(5) * t528 + t385;
t319 = t377 * t440 + t378 * t438;
t317 = t377 * t438 - t378 * t440;
t303 = -t411 * t522 - t442 * t487 - t501 * t528 + (t527 * t549 - t530) * t446;
t302 = -t359 * t549 + t418 * t411;
t298 = qJD(5) * t463 - t319 * t443 + t447 * t489;
t297 = qJD(5) * t336 + t319 * t447 + t443 * t489;
t293 = pkin(5) * t340 + t305;
t285 = -t295 * t442 + t544;
t1 = [(t317 * t409 + t319 * t407 + t357 * t399 - t358 * t398) * MDP(12) + (t305 * t357 + t306 * t358 - t315 * t317 + t316 * t319) * MDP(13) + (t298 * t401 + t317 * t387 + t336 * t398 + t340 * t357) * MDP(19) + (-t297 * t401 + t317 * t389 + t339 * t357 + t398 * t463) * MDP(20) + ((-qJD(6) * t466 - t297 * t442 + t298 * t446) * t397 + t467 * t398 + t317 * t326 + t357 * t291) * MDP(26) + (-(qJD(6) * t467 + t297 * t446 + t298 * t442) * t397 - t466 * t398 - t317 * t465 + t357 * t290) * MDP(27) + (MDP(10) * t378 - MDP(11) * t377) * qJD(3) + (-t406 * t449 * MDP(13) + (t400 * t445 * MDP(13) - qJD(3) * t449 * t550) * qJD(2) + (-t449 * MDP(4) + (-MDP(10) * t448 + MDP(11) * t444 - MDP(3)) * t445) * t451) * t439; 0.2e1 * t484 * t555 - 0.2e1 * t496 * t554 + MDP(7) * t520 - MDP(8) * t521 + (-pkin(8) * t520 + t444 * t454) * MDP(10) + (pkin(8) * t521 + t448 * t454) * MDP(11) + (-t306 * t415 + t315 * t411 - t316 * t408 + t385 * t399 + t407 * t514 + t409 * t515 + t469) * MDP(12) + (t305 * t385 + t306 * t386 - t515 * t315 + t514 * t316 + t400 * t563 + t406 * t493) * MDP(13) + (t339 * t527 + t389 * t460) * MDP(14) + (-(-t387 * t447 - t389 * t443) * t411 + (-t540 - t340 * t447 + (t387 * t443 - t389 * t447) * qJD(5)) * t416) * MDP(15) + (t339 * t415 + t389 * t408 + t390 * t416 + t401 * t460) * MDP(16) + (-t340 * t415 - t387 * t408 - t398 * t528 - t401 * t461) * MDP(17) + (t401 * t408 + t366) * MDP(18) + (t362 * t398 + (-t313 * t502 + t322) * t415 + t300 * t408 + t385 * t340 + t312 * t486 + (-t386 * t502 + t553) * t401 + t515 * t387 + ((-qJD(5) * t365 - t351) * t401 + (-qJD(5) * t332 - t306) * t415 - t312 * t411 + t469) * t443) * MDP(19) + (-t301 * t408 + t305 * t527 + t460 * t312 + t385 * t339 + t515 * t389 - t513 * t398 + t401 * t552 - t458 * t415) * MDP(20) + (-t290 * t360 - t302 * t465) * MDP(21) + (-t290 * t359 + t291 * t360 - t302 * t326 + t303 * t465) * MDP(22) + (t290 * t415 + t302 * t397 - t360 * t398 - t408 * t465) * MDP(23) + (-t291 * t415 - t303 * t397 - t326 * t408 - t359 * t398) * MDP(24) + (t397 * t408 + t366) * MDP(25) + ((t307 * t446 - t309 * t442) * t398 + t481 * t415 + t285 * t408 + t349 * t291 + t293 * t359 + t308 * t303 + (t442 * t518 + t446 * t519) * t397 + t517 * t326 + (-t286 * t415 - t397 * t468) * qJD(6)) * MDP(26) + (-t468 * t398 - (t476 * t446 + t480) * t415 - t286 * t408 + t349 * t290 - t293 * t360 + t308 * t302 + ((-qJD(6) * t307 + t518) * t446 + (qJD(6) * t309 - t519) * t442) * t397 - t517 * t465) * MDP(27); ((t315 - t320) * t407 + (-t398 * t438 - t399 * t440) * pkin(3)) * MDP(12) + (t315 * t318 - t316 * t320 + (-t305 * t440 + t306 * t438 - t400 * t505) * pkin(3)) * MDP(13) + (t389 * t478 + t540) * MDP(14) + ((t339 - t539) * t447 + (-t340 - t537) * t443) * MDP(15) + (t401 * t478 + t533 - t536) * MDP(16) + (t464 + t538) * MDP(17) + (-t305 * t447 - t318 * t387 + t432 * t340 + (-t431 * t502 - t345) * t401 + (t320 * t401 + t456) * t443) * MDP(19) + (t305 * t443 - t318 * t389 + t432 * t339 + (t431 * t503 + t516) * t401 + t456 * t447) * MDP(20) + (t290 * t419 - t465 * t512) * MDP(21) + (-t290 * t418 - t291 * t419 - t326 * t512 - t465 * t511) * MDP(22) + (t541 - t548) * MDP(23) + (t474 + t542) * MDP(24) + ((-t413 * t446 - t414 * t442) * t398 + t422 * t291 + t293 * t418 + (t442 * t471 - t446 * t472) * t397 + t475 * t326 - t511 * t308) * MDP(26) + (-(-t413 * t442 + t414 * t446) * t398 + t422 * t290 + t293 * t419 + (t442 * t472 + t446 * t471) * t397 - t475 * t465 + t512 * t308) * MDP(27) + t550 * qJD(2) * t545 + (-t448 * t555 + t554) * t451 - ((-t316 + t318) * MDP(12) + t401 * MDP(18) + t300 * MDP(19) - t301 * MDP(20) + t397 * MDP(25) + t285 * MDP(26) - t286 * MDP(27)) * t409; (-t407 ^ 2 - t409 ^ 2) * MDP(12) + (t315 * t409 - t316 * t407 + t406) * MDP(13) + (t464 - t538) * MDP(19) + (-t401 ^ 2 * t447 - t533 - t536) * MDP(20) + (t474 - t542) * MDP(26) + (t541 + t548) * MDP(27); t389 * t387 * MDP(14) + (-t387 ^ 2 + t389 ^ 2) * MDP(15) + (t339 + t539) * MDP(16) + (-t340 + t537) * MDP(17) + t398 * MDP(18) + (t301 * t401 - t312 * t389 + t453) * MDP(19) + (t300 * t401 + t312 * t387 - t458) * MDP(20) + (t290 + t560) * MDP(23) + (-t291 - t559) * MDP(24) + (-(-t294 * t442 - t543) * t397 - t286 * qJD(6) + (-t326 * t389 - t397 * t501 + t398 * t446) * pkin(5) + t557) * MDP(26) + (t561 + t292 + (-t295 * t397 - t283) * t442 + (t294 * t397 - t476) * t446 + (t389 * t465 - t397 * t500 - t398 * t442) * pkin(5)) * MDP(27) + t556; (t494 + t560) * MDP(23) + (-t479 - t559) * MDP(24) + (t286 * t397 + t557) * MDP(26) + (-t446 * t284 + t285 * t397 - t480 + t561) * MDP(27) + (-MDP(23) * t535 + MDP(24) * t465 - MDP(26) * t286 - MDP(27) * t544) * qJD(6) + t556;];
tauc  = t1;
