% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:55
% EndTime: 2019-03-09 08:24:03
% DurationCPUTime: 5.53s
% Computational Cost: add. (2299->496), mult. (5703->640), div. (0->0), fcn. (3602->6), ass. (0->225)
t566 = MDP(15) - MDP(20);
t565 = MDP(17) + MDP(19);
t501 = (qJD(1) * qJD(2));
t564 = -2 * t501;
t440 = cos(qJ(2));
t434 = t440 ^ 2;
t438 = sin(qJ(2));
t563 = MDP(5) * (t438 ^ 2 - t434);
t435 = sin(pkin(9));
t516 = qJD(1) * t438;
t493 = t435 * t516;
t436 = cos(pkin(9));
t503 = t436 * qJD(2);
t376 = t493 - t503;
t549 = pkin(3) + qJ(5);
t562 = t376 * t549;
t561 = t436 * t549;
t484 = t440 * t501;
t402 = t435 * t484;
t559 = qJ(5) * t402 + t376 * qJD(5);
t507 = qJD(4) * t435;
t515 = qJD(1) * t440;
t421 = pkin(7) * t515;
t492 = t435 * t515;
t520 = pkin(3) * t492 + t421;
t558 = -qJD(5) * t436 - t507 - t520;
t532 = t436 * t438;
t557 = pkin(4) * t532 + t440 * qJ(5);
t506 = qJD(4) * t440;
t513 = qJD(2) * t438;
t556 = -qJ(4) * t513 + t506;
t553 = -MDP(16) + MDP(21);
t512 = qJD(2) * t440;
t423 = pkin(7) * t512;
t489 = t435 * t512;
t519 = pkin(3) * t489 + t423;
t552 = -t438 * (qJD(4) * t436 - qJD(5) * t435) + t519;
t372 = t376 ^ 2;
t491 = t436 * t516;
t514 = qJD(2) * t435;
t378 = t491 + t514;
t551 = t378 ^ 2;
t550 = pkin(4) + pkin(8);
t548 = -pkin(5) - qJ(4);
t547 = qJD(2) * pkin(2);
t546 = qJ(3) * t438;
t545 = qJ(5) * t435;
t392 = -pkin(2) * t440 - pkin(1) - t546;
t367 = t392 * qJD(1);
t400 = qJD(2) * qJ(3) + t421;
t330 = t367 * t436 - t435 * t400;
t320 = pkin(3) * t515 + qJD(4) - t330;
t544 = t320 * t438;
t331 = t435 * t367 + t436 * t400;
t322 = qJ(4) * t515 - t331;
t543 = t322 * t438;
t439 = cos(qJ(6));
t437 = sin(qJ(6));
t538 = t376 * t437;
t326 = -t439 * t378 + t538;
t414 = -qJD(6) + t515;
t542 = t326 * t414;
t328 = t376 * t439 + t378 * t437;
t541 = t328 * t414;
t465 = pkin(2) * t438 - qJ(3) * t440;
t360 = qJD(2) * t465 - qJD(3) * t438;
t540 = t360 * t436;
t539 = t376 * t378;
t384 = t465 * qJD(1);
t537 = t384 * t436;
t536 = t414 * t440;
t442 = qJD(1) ^ 2;
t535 = t434 * t442;
t534 = t435 * t438;
t533 = t435 * t440;
t426 = t436 * qJ(3);
t531 = t436 * t439;
t530 = t436 * t440;
t441 = qJD(2) ^ 2;
t529 = t438 * t441;
t528 = t440 * t441;
t527 = t440 * t442;
t464 = qJ(4) * t436 - t545;
t451 = t440 * t464;
t526 = qJD(1) * t451 + t558;
t380 = t435 * t439 + t436 * t437;
t362 = t380 * qJD(6);
t450 = t440 * t380;
t525 = -qJD(1) * t450 + t362;
t381 = t435 * t437 - t531;
t363 = t381 * qJD(6);
t490 = t436 * t515;
t524 = -t437 * t492 + t439 * t490 + t363;
t365 = t435 * t384;
t417 = qJ(4) * t516;
t523 = t365 + t417;
t352 = t360 * qJD(1);
t420 = pkin(7) * t516;
t388 = (qJD(3) - t420) * qJD(2);
t319 = t435 * t352 + t436 * t388;
t455 = t436 * t548 + t545;
t445 = t440 * t455;
t522 = -qJD(1) * t445 + t558;
t413 = pkin(7) * t484;
t521 = pkin(3) * t402 + t413;
t518 = pkin(3) * t534 + t438 * pkin(7);
t412 = pkin(7) * t530;
t349 = t435 * t392 + t412;
t396 = (pkin(4) + qJ(3)) * t435;
t397 = t436 * pkin(4) + t426;
t511 = qJD(3) * t376;
t510 = qJD(3) * t378;
t509 = qJD(3) * t436;
t508 = qJD(4) * t378;
t505 = qJD(6) * t437;
t504 = qJD(6) * t439;
t502 = t438 * MDP(27);
t427 = t440 * qJD(5);
t500 = pkin(8) * t530;
t410 = pkin(7) * t533;
t499 = pkin(3) * t513;
t498 = pkin(7) * t513;
t497 = t435 * t550;
t415 = t438 * t501;
t496 = qJ(4) * t415 + t319;
t472 = t436 * t484;
t495 = t378 * t504 + t439 * t402 + t437 * t472;
t494 = -pkin(7) * t435 - pkin(3);
t488 = t440 * t503;
t487 = t414 * t505;
t486 = t414 * t504;
t485 = t548 * t440;
t483 = qJ(4) * t435 + pkin(2);
t480 = -qJD(3) + t547;
t468 = -t420 + t480;
t482 = t468 - t547;
t481 = pkin(1) * t564;
t318 = t352 * t436 - t435 * t388;
t456 = pkin(4) * t472 + qJD(1) * t427 - t318;
t285 = (-t438 * t549 + t500) * t501 + t456;
t473 = t440 * t497;
t286 = (-t506 + (pkin(5) * t438 - t473) * qJD(2)) * qJD(1) + t496;
t479 = -t285 * t437 + t439 * t286;
t477 = pkin(4) * t490 - t537;
t348 = t392 * t436 - t410;
t447 = qJ(4) * t378 + t468;
t300 = t447 - t562;
t364 = t483 + t561;
t476 = -qJD(2) * t364 + t300;
t475 = t376 + t503;
t347 = -qJ(4) * t490 + t520;
t474 = t347 + t507;
t471 = qJD(1) * t485;
t470 = t494 * t438;
t469 = -pkin(4) * t514 - qJD(4);
t342 = qJ(4) * t440 - t349;
t467 = pkin(4) * t488 + t427 - t540;
t432 = t440 * pkin(3);
t343 = -t348 + t432;
t353 = t435 * t360;
t466 = t353 - t556;
t463 = t285 * t439 + t286 * t437;
t449 = qJ(5) * t515 + t320;
t288 = t378 * t550 + t449;
t289 = -t376 * t550 + qJD(5) + t331 + t471;
t283 = t288 * t439 + t289 * t437;
t462 = t288 * t437 - t289 * t439;
t296 = pkin(4) * t378 + t449;
t301 = -pkin(4) * t376 + qJD(5) - t322;
t461 = -t296 * t436 + t301 * t435;
t310 = t410 + t432 + (pkin(8) * t438 - t392) * t436 + t557;
t311 = -t438 * t497 + t349 + t485;
t460 = t310 * t439 + t311 * t437;
t459 = -t310 * t437 + t311 * t439;
t345 = -pkin(7) * t491 + t365;
t337 = -t436 * t498 + t353;
t457 = (-qJ(5) + t494) * t438;
t371 = pkin(8) * t436 + t397;
t444 = t457 + t500;
t454 = -qJD(1) * t444 + qJD(3) * t435 + qJD(6) * t371 - t477;
t370 = pkin(8) * t435 + t396;
t443 = (-pkin(7) * t436 + pkin(5)) * t438 - t473;
t453 = -qJD(1) * t443 - qJD(6) * t370 + t509 - t523;
t452 = -pkin(4) * t533 - pkin(7) * t532;
t297 = -t376 * t505 + t495;
t448 = t402 * t437 - t439 * t472;
t446 = -qJ(4) * t472 + t521;
t312 = t446 - t508;
t298 = qJD(6) * t328 + t448;
t401 = qJD(3) * t492;
t389 = -pkin(3) * t436 - t483;
t358 = t380 * t438;
t357 = t437 * t534 - t438 * t531;
t355 = -qJ(4) * t532 + t518;
t354 = t435 * t548 - pkin(2) - t561;
t344 = pkin(7) * t493 + t537;
t341 = t438 * t464 - t518;
t336 = t435 * t498 + t540;
t335 = qJD(1) * t470 - t537;
t334 = -t345 - t417;
t333 = (-qJ(4) * t512 - qJD(4) * t438) * t436 + t519;
t329 = -pkin(4) * t534 - t342;
t325 = t438 * t455 + t518;
t324 = qJD(2) * t470 - t540;
t323 = t343 + t557;
t321 = qJD(1) * t452 + t523;
t317 = -t337 + t556;
t316 = pkin(3) * t376 - t447;
t314 = qJD(2) * t450 - t363 * t438;
t313 = t362 * t438 + t437 * t489 - t439 * t488;
t309 = qJD(1) * t457 + t477;
t308 = -pkin(3) * t415 - t318;
t307 = qJD(2) * t451 - t552;
t306 = qJD(2) * t452 + t466;
t305 = qJD(1) * t506 - t496;
t302 = qJD(2) * t457 + t467;
t299 = qJD(2) * t445 + t552;
t295 = qJD(2) * t443 + t466;
t294 = qJD(2) * t444 + t467;
t293 = t469 * t515 + t496;
t292 = -t312 - t559;
t291 = -t415 * t549 + t456;
t290 = t378 * t548 - t468 + t562;
t287 = t471 * t503 - t508 + t521 + t559;
t1 = [(-t312 * t532 - t333 * t378 + (qJD(1) * t317 + t305) * t440 + (-t316 * t530 - t543 + (-t342 * t438 - t355 * t530) * qJD(1)) * qJD(2)) * MDP(17) + (-t312 * t534 - t333 * t376 + (-qJD(1) * t324 - t308) * t440 + (-t316 * t533 + t544 + (t343 * t438 - t355 * t533) * qJD(1)) * qJD(2)) * MDP(16) + (t292 * t532 + t307 * t378 + (-qJD(1) * t306 - t293) * t440 + (t300 * t530 + t301 * t438 + (t329 * t438 + t341 * t530) * qJD(1)) * qJD(2)) * MDP(19) + (-t292 * t534 - t307 * t376 + (qJD(1) * t302 + t291) * t440 + (-t300 * t533 - t296 * t438 + (-t323 * t438 - t341 * t533) * qJD(1)) * qJD(2)) * MDP(21) - MDP(7) * t529 + (pkin(7) * t529 + t440 * t481) * MDP(10) + (-pkin(7) * t528 + t438 * t481) * MDP(9) + (-t336 * t378 - t337 * t376 + (-t318 * t436 - t319 * t435) * t438 + (-t330 * t436 - t331 * t435 + (-t348 * t436 - t349 * t435) * qJD(1)) * t512) * MDP(13) + (t317 * t376 + t324 * t378 + (t305 * t435 + t308 * t436) * t438 + (t320 * t436 + t322 * t435 + (t342 * t435 + t343 * t436) * qJD(1)) * t512) * MDP(15) + (-t302 * t378 + t306 * t376 + (-t291 * t436 + t293 * t435) * t438 + ((-t323 * t436 + t329 * t435) * qJD(1) + t461) * t512) * MDP(20) + (-t297 * t440 - t314 * t414 + (qJD(1) * t358 + t328) * t513) * MDP(25) + (t298 * t440 + t313 * t414 + (-qJD(1) * t357 - t326) * t513) * MDP(26) + ((t294 * t439 + t295 * t437) * t414 + t463 * t440 + t299 * t328 + t325 * t297 + t287 * t358 + t290 * t314 + (t414 * t459 - t440 * t462) * qJD(6) + (-qJD(1) * t460 - t283) * t513) * MDP(29) + (t291 * t323 + t292 * t341 + t293 * t329 + t296 * t302 + t300 * t307 + t301 * t306) * MDP(22) + (-t414 - t515) * qJD(2) * t502 + 0.2e1 * t440 * MDP(4) * t415 + (t305 * t342 + t308 * t343 + t312 * t355 + t316 * t333 + t317 * t322 + t320 * t324) * MDP(18) + (t297 * t358 + t314 * t328) * MDP(23) + (-t297 * t357 - t298 * t358 - t313 * t328 - t314 * t326) * MDP(24) + (-(-t294 * t437 + t295 * t439) * t414 - t479 * t440 + t299 * t326 + t325 * t298 + t287 * t357 + t290 * t313 + (t283 * t440 + t414 * t460) * qJD(6) + (qJD(1) * t459 - t462) * t513) * MDP(28) + ((-qJD(1) * t336 - t318) * t440 + ((pkin(7) * t376 - t435 * t468) * t440 + (t330 + (t348 + 0.2e1 * t410) * qJD(1)) * t438) * qJD(2)) * MDP(11) + (t318 * t348 + t319 * t349 + t330 * t336 + t331 * t337 + (-t468 + t420) * t423) * MDP(14) + ((qJD(1) * t337 + t319) * t440 + ((pkin(7) * t378 - t436 * t468) * t440 + (-t331 + (-t349 + 0.2e1 * t412) * qJD(1)) * t438) * qJD(2)) * MDP(12) + t563 * t564 + MDP(6) * t528; -t438 * MDP(4) * t527 + t442 * t563 + t401 * MDP(11) + (t344 * t378 + t345 * t376 + (t330 * t515 + t319 - t511) * t436 + (t331 * t515 - t318 + t510) * t435) * MDP(13) + (-t330 * t344 - t331 * t345 + (-t330 * t435 + t331 * t436) * qJD(3) + (-t318 * t435 + t319 * t436) * qJ(3) + t482 * t421) * MDP(14) + (-t334 * t376 - t335 * t378 + (-t320 * t515 - t305 - t511) * t436 + (-t322 * t515 + t308 + t510) * t435) * MDP(15) + (t312 * t436 + t474 * t376 - t401) * MDP(16) + (-t312 * t435 + t474 * t378) * MDP(17) + (-t305 * t426 + t312 * t389 - t316 * t347 - t320 * t335 + (-t334 - t509) * t322 + (qJ(3) * t308 + qJD(3) * t320 - qJD(4) * t316) * t435) * MDP(18) + (t292 * t435 - t526 * t378) * MDP(19) + (-t291 * t435 - t293 * t436 + t309 * t378 - t321 * t376 + (t376 * t436 - t378 * t435) * qJD(3) + ((-t396 * t436 + t397 * t435) * qJD(2) - t461) * t515) * MDP(20) + (t292 * t436 + t526 * t376 + t401) * MDP(21) + (t291 * t396 + t292 * t364 + t293 * t397 - t296 * t309 - t301 * t321 - t526 * t300 + (t296 * t435 + t301 * t436) * qJD(3)) * MDP(22) + (t297 * t381 + t328 * t525) * MDP(23) + (t297 * t380 - t298 * t381 - t326 * t525 - t328 * t524) * MDP(24) + (-t525 * t414 + (qJD(2) * t381 - t328) * t516) * MDP(25) + (t524 * t414 + (qJD(2) * t380 + t326) * t516) * MDP(26) + (-t287 * t380 + t354 * t298 + (t437 * t454 - t439 * t453) * t414 + t522 * t326 + t524 * t290 + ((-t370 * t437 + t371 * t439) * qJD(2) + t462) * t516) * MDP(28) + (t287 * t381 + t354 * t297 + (t437 * t453 + t439 * t454) * t414 + t522 * t328 + t525 * t290 + (-(t370 * t439 + t371 * t437) * qJD(2) + t283) * t516) * MDP(29) + (((-qJ(3) * t514 - t330) * t438 + (-pkin(7) * t475 + t435 * t482 + t344) * t440) * MDP(11) + ((-qJ(3) * t503 + t331) * t438 + (-t345 + (-t378 + t514) * pkin(7) + (t468 - t480) * t436) * t440) * MDP(12) + (-t544 + t335 * t440 + (t316 * t440 + (-t389 * t440 + t546) * qJD(2)) * t435) * MDP(16) + (t543 - t334 * t440 + (qJ(3) * t513 + (-qJD(2) * t389 - qJD(3) + t316) * t440) * t436) * MDP(17) + ((qJD(2) * t397 - t301) * t438 + (t321 + (-qJD(3) - t476) * t436) * t440) * MDP(19) + ((-qJD(2) * t396 + t296) * t438 + (t435 * t476 - t309) * t440) * MDP(21) + t414 * t502) * qJD(1) + (MDP(9) * t438 * t442 + MDP(10) * t527) * pkin(1); (t330 * t378 + t331 * t376 + t413) * MDP(14) + (-t322 * t376 + (-qJD(4) - t320) * t378 + t446) * MDP(18) + (t301 * t376 + (-qJD(4) - t296) * t378 + t446 + t559) * MDP(22) + (t298 - t541) * MDP(28) + (t297 + t542) * MDP(29) + (MDP(11) + t553) * (-t378 * t515 + t402) + (MDP(12) - t565) * t475 * t515 + (-MDP(13) - t566) * (t551 + t372); (t316 * t378 - t318) * MDP(18) + (-t300 * t378 + t456) * MDP(22) + (t326 * t378 + t486) * MDP(28) + (t328 * t378 - t487) * MDP(29) + t566 * (-t376 + t503) * t515 + ((-t322 * t440 - t499) * MDP(18) + (-qJ(5) * t513 + t301 * t440 - t499) * MDP(22) + (-t437 * t513 - t439 * t536) * MDP(28) + (t437 * t536 - t439 * t513) * MDP(29)) * qJD(1) + t553 * (-t415 + t539) + t565 * (-t551 - t535); (t415 + t539) * MDP(19) + t402 * MDP(20) + (-t372 - t535) * MDP(21) + (t300 * t376 + t496) * MDP(22) + (-t326 * t376 + t487) * MDP(28) + (-t328 * t376 + t486) * MDP(29) + ((t439 * MDP(28) - t437 * MDP(29)) * t513 + (t378 * MDP(20) + (-t296 + t469) * MDP(22) + (-t437 * MDP(28) - t439 * MDP(29)) * t414) * t440) * qJD(1); t328 * t326 * MDP(23) + (-t326 ^ 2 + t328 ^ 2) * MDP(24) + (t495 - t542) * MDP(25) + (-t448 - t541) * MDP(26) + MDP(27) * t415 + (-t283 * t414 - t290 * t328 + t479) * MDP(28) + (t290 * t326 + t414 * t462 - t463) * MDP(29) + (-MDP(25) * t538 - MDP(26) * t328 - MDP(28) * t283 + MDP(29) * t462) * qJD(6);];
tauc  = t1;
