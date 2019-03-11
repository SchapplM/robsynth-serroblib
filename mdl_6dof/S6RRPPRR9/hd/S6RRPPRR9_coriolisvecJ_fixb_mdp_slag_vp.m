% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:20
% EndTime: 2019-03-09 09:32:32
% DurationCPUTime: 6.96s
% Computational Cost: add. (3463->528), mult. (9044->695), div. (0->0), fcn. (6438->8), ass. (0->221)
t445 = cos(pkin(6));
t528 = qJD(1) * t445;
t429 = qJD(2) + t528;
t449 = sin(qJ(5));
t452 = cos(qJ(5));
t450 = sin(qJ(2));
t444 = sin(pkin(6));
t529 = qJD(1) * t444;
t506 = t450 * t529;
t380 = t429 * t452 + t449 * t506;
t453 = cos(qJ(2));
t511 = qJD(1) * qJD(2);
t500 = t444 * t511;
t411 = t453 * t500;
t348 = qJD(5) * t380 - t452 * t411;
t527 = qJD(1) * t453;
t505 = t444 * t527;
t408 = qJD(5) + t505;
t551 = t380 * t408;
t571 = -t348 - t551;
t401 = t452 * t506;
t378 = t429 * t449 - t401;
t377 = qJD(6) + t378;
t442 = t450 ^ 2;
t443 = t453 ^ 2;
t570 = (t442 - t443) * MDP(5);
t448 = sin(qJ(6));
t451 = cos(qJ(6));
t341 = t380 * t448 - t451 * t408;
t569 = t341 * t408;
t558 = pkin(1) * t450;
t434 = t445 * t558;
t547 = t444 * t453;
t565 = pkin(8) * t547 + t434;
t382 = -t445 * qJ(3) - t565;
t366 = pkin(3) * t547 - t382;
t340 = pkin(4) * t547 - pkin(9) * t445 + t366;
t446 = qJ(3) - pkin(9);
t447 = pkin(2) + qJ(4);
t484 = -t447 * t453 - pkin(1);
t352 = (-t446 * t450 + t484) * t444;
t568 = t449 * t340 + t452 * t352;
t522 = qJD(5) * t449;
t347 = qJD(5) * t401 + t449 * t411 - t429 * t522;
t552 = t378 * t408;
t567 = -t347 + t552;
t419 = t429 * pkin(2);
t566 = -t429 * qJ(4) - t419;
t510 = pkin(1) * t528;
t386 = pkin(8) * t506 - t453 * t510;
t514 = qJD(3) + t386;
t387 = pkin(8) * t505 + t450 * t510;
t370 = pkin(3) * t505 + t387;
t512 = qJD(4) + t370;
t564 = -MDP(13) - MDP(16);
t414 = t429 * qJ(3);
t513 = pkin(4) * t505 + t512;
t318 = -pkin(9) * t429 + t414 + t513;
t331 = qJD(1) * t352;
t300 = t318 * t449 + t331 * t452;
t523 = qJD(3) * t450;
t455 = (-t523 + (-qJD(2) * t446 - qJD(4)) * t453) * t444;
t486 = t450 * t500;
t402 = pkin(2) * t486;
t535 = qJ(4) * t486 + t402;
t313 = qJD(1) * t455 + t535;
t410 = t429 * qJD(3);
t559 = pkin(1) * t445;
t509 = qJD(2) * t559;
t489 = qJD(1) * t509;
t534 = pkin(8) * t486 - t453 * t489;
t359 = -t410 + t534;
t548 = t444 * t450;
t490 = (-pkin(3) - pkin(4)) * t548;
t480 = qJD(1) * t490;
t322 = qJD(2) * t480 - t359;
t456 = -qJD(5) * t300 - t313 * t449 + t452 * t322;
t290 = pkin(5) * t486 - t456;
t563 = t377 * (pkin(5) * t380 + pkin(10) * t377) + t290;
t350 = t480 - t386;
t516 = qJD(3) - t350;
t317 = -t516 - t566;
t526 = qJD(2) * t450;
t562 = (-t317 * t453 + t446 * t526) * t529 - t317 * qJD(5);
t504 = t444 * t526;
t422 = pkin(2) * t504;
t533 = qJ(4) * t504 + t422;
t321 = t455 + t533;
t560 = pkin(3) + pkin(8);
t497 = t444 * (-pkin(4) - t560);
t426 = t453 * t509;
t437 = t445 * qJD(3);
t532 = t426 + t437;
t333 = t497 * t526 + t532;
t561 = -qJD(5) * t568 - t321 * t449 + t333 * t452;
t428 = t429 ^ 2;
t557 = qJ(3) * t450;
t461 = -qJD(6) * t380 - t486;
t518 = qJD(6) * t451;
t538 = t451 * t347 + t408 * t518;
t303 = t448 * t461 + t538;
t556 = t303 * t448;
t555 = t341 * t377;
t554 = t341 * t380;
t343 = t380 * t451 + t408 * t448;
t553 = t343 * t380;
t550 = t408 * t452;
t441 = t444 ^ 2;
t549 = t441 * qJD(1) ^ 2;
t546 = t448 * t348;
t545 = t448 * t377;
t544 = t449 * t453;
t543 = t450 * t453;
t542 = t451 * t348;
t496 = t451 * t377;
t541 = t451 * t453;
t540 = t452 * t303;
t483 = pkin(5) * t452 + pkin(10) * t449;
t539 = (-pkin(4) - t483) * t505 - qJD(5) * t483 - t512;
t415 = qJ(4) * t506;
t427 = pkin(2) * t506;
t349 = -t446 * t505 + t415 + t427;
t537 = t452 * t349 + t449 * t350;
t381 = pkin(8) * t411 + t450 * t489;
t530 = qJ(3) * qJD(2);
t525 = qJD(2) * t452;
t524 = qJD(2) * t453;
t521 = qJD(5) * t452;
t298 = pkin(10) * t408 + t300;
t520 = qJD(6) * t298;
t519 = qJD(6) * t448;
t369 = -pkin(3) * t506 - t386;
t515 = qJD(3) - t369;
t508 = t408 * t544;
t507 = t453 * t549;
t384 = -t445 * pkin(2) + pkin(8) * t548 - t453 * t559;
t503 = t444 * t524;
t502 = t408 * t522;
t501 = t441 * t511;
t463 = t452 * t313 + t318 * t521 + t449 * t322 - t331 * t522;
t289 = -pkin(10) * t486 + t463;
t302 = pkin(5) * t378 - pkin(10) * t380 + t317;
t494 = -qJD(6) * t302 - t289;
t398 = pkin(5) * t449 - pkin(10) * t452 + t447;
t493 = -pkin(10) * t506 - qJD(6) * t398 + t537;
t491 = -qJD(4) - t530;
t488 = t450 * t507;
t336 = pkin(3) * t411 - t429 * qJD(4) + t381;
t487 = t445 * qJ(4) - t384;
t485 = -0.2e1 * pkin(1) * t501;
t482 = t387 * t429 - t381;
t375 = (-t448 * t450 + t449 * t541) * t529;
t481 = -t451 * t522 - t375;
t292 = t298 * t451 + t302 * t448;
t291 = -t298 * t448 + t302 * t451;
t306 = pkin(10) * t547 + t568;
t339 = t490 + t487;
t390 = t445 * t452 + t449 * t548;
t470 = -t445 * t449 + t452 * t548;
t311 = -pkin(5) * t470 - pkin(10) * t390 + t339;
t479 = t306 * t451 + t311 * t448;
t478 = -t306 * t448 + t311 * t451;
t299 = t318 * t452 - t331 * t449;
t476 = t340 * t452 - t352 * t449;
t388 = t565 * qJD(2);
t474 = t381 * t445 + t388 * t429;
t473 = -pkin(8) * t504 + t426;
t385 = -qJ(3) * t505 + t427;
t372 = -t429 * t505 + t411;
t472 = (-qJD(2) + t429) * t506;
t471 = -t390 * t448 + t444 * t541;
t363 = t390 * t451 + t448 * t547;
t468 = -t377 * t518 - t546;
t467 = -t377 * t519 + t542;
t464 = t408 * t343;
t462 = t452 * t321 + t449 * t333 + t340 * t521 - t352 * t522;
t383 = (-pkin(2) * t453 - pkin(1) - t557) * t444;
t460 = (-qJ(3) * t524 - t523) * t444;
t335 = -pkin(3) * t486 - t359;
t365 = (t484 - t557) * t444;
t297 = -pkin(5) * t408 - t299;
t459 = -pkin(10) * t348 + (t297 + t299) * t377;
t323 = -pkin(4) * t411 - t336;
t458 = (t453 * t491 - t523) * t444;
t436 = t445 * qJD(4);
t334 = t436 + (t453 * t497 - t434) * qJD(2);
t412 = qJD(3) * t505;
t394 = t449 * t486;
t376 = -t437 - t473;
t374 = (t448 * t544 + t450 * t451) * t529;
t373 = qJD(1) * t383;
t371 = t422 + t460;
t368 = -t414 - t387;
t367 = t385 + t415;
t364 = -t419 + t514;
t361 = qJD(5) * t390 - t452 * t503;
t360 = qJD(5) * t470 + t449 * t503;
t358 = qJD(1) * t460 + t402;
t357 = -t436 + (t547 * t560 + t434) * qJD(2);
t356 = -t504 * t560 + t532;
t355 = pkin(3) * t548 - t487;
t354 = t373 * t506;
t353 = qJD(1) * t365;
t346 = t414 + t512;
t332 = t353 * t505;
t329 = t458 + t533;
t326 = t515 + t566;
t324 = qJD(1) * t458 + t535;
t310 = qJD(6) * t471 + t360 * t451 - t448 * t504;
t309 = qJD(6) * t363 + t360 * t448 + t451 * t504;
t307 = pkin(5) * t506 + t349 * t449 - t350 * t452;
t305 = -pkin(5) * t547 - t476;
t304 = qJD(6) * t343 + t347 * t448 + t451 * t486;
t301 = pkin(5) * t361 - pkin(10) * t360 + t334;
t296 = pkin(5) * t348 - pkin(10) * t347 + t323;
t295 = t451 * t296;
t294 = pkin(5) * t504 - t561;
t293 = -pkin(10) * t504 + t462;
t288 = -qJD(6) * t292 - t289 * t448 + t295;
t287 = qJD(6) * t291 + t289 * t451 + t296 * t448;
t1 = [(-t429 * t473 + t445 * t534 + t453 * t485) * MDP(10) + ((-t373 * t526 + t358 * t453 + (t371 * t453 - t383 * t526) * qJD(1)) * t444 + t474) * MDP(12) + (-t336 * t445 - t357 * t429 + (t353 * t526 - t324 * t453 + (-t329 * t453 + t365 * t526) * qJD(1)) * t444) * MDP(17) + (-t359 * t445 - t376 * t429 + (-t373 * t524 - t358 * t450 + (-t371 * t450 - t383 * t524) * qJD(1)) * t444) * MDP(13) + (t335 * t445 + t356 * t429 + (-t353 * t524 - t324 * t450 + (-t329 * t450 - t365 * t524) * qJD(1)) * t444) * MDP(16) + (t360 * t408 + (t347 * t453 + (-qJD(1) * t390 - t380) * t526) * t444) * MDP(21) + (t450 * t485 - t474) * MDP(9) + (MDP(6) * t503 - MDP(7) * t504) * (t429 + t528) + (t347 * t390 + t360 * t380) * MDP(19) + (t358 * t383 + t359 * t382 + t364 * t388 + t368 * t376 + t371 * t373 + t381 * t384) * MDP(14) + (t303 * t363 + t310 * t343) * MDP(26) + (t324 * t365 + t326 * t357 + t329 * t353 + t335 * t366 + t336 * t355 + t346 * t356) * MDP(18) + (-t359 * t453 + t381 * t450 + (t364 * t453 + t368 * t450) * qJD(2) + (-t376 * t453 + t388 * t450 + (t382 * t450 + t384 * t453) * qJD(2)) * qJD(1)) * t444 * MDP(11) + (t335 * t453 + t336 * t450 + (t326 * t453 - t346 * t450) * qJD(2) + (t356 * t453 + t357 * t450 + (t355 * t453 - t366 * t450) * qJD(2)) * qJD(1)) * t444 * MDP(15) + (t303 * t471 - t304 * t363 - t309 * t343 - t310 * t341) * MDP(27) + (t347 * t470 - t348 * t390 - t360 * t378 - t361 * t380) * MDP(20) + (-t348 * t470 + t361 * t377) * MDP(30) + (-t303 * t470 + t310 * t377 + t343 * t361 + t348 * t363) * MDP(28) + (-t361 * t408 + (-t348 * t453 + (-qJD(1) * t470 + t378) * t526) * t444) * MDP(22) + (-(qJD(6) * t478 + t293 * t451 + t301 * t448) * t377 - t479 * t348 + t287 * t470 - t292 * t361 + t294 * t343 + t305 * t303 + t290 * t363 + t297 * t310) * MDP(32) + ((-qJD(6) * t479 - t293 * t448 + t301 * t451) * t377 + t478 * t348 - t288 * t470 + t291 * t361 + t294 * t341 + t305 * t304 - t290 * t471 + t297 * t309) * MDP(31) + (t304 * t470 - t309 * t377 - t341 * t361 + t348 * t471) * MDP(29) + (t561 * t408 + t334 * t378 + t339 * t348 - t323 * t470 + t317 * t361 + (t456 * t453 + (-qJD(1) * t476 - t299) * t526) * t444) * MDP(24) + (-t462 * t408 + t334 * t380 + t339 * t347 + t323 * t390 + t317 * t360 + (-t463 * t453 + (qJD(1) * t568 + t300) * t526) * t444) * MDP(25) + (-t408 * t444 - t441 * t527) * MDP(23) * t526 + 0.2e1 * (MDP(4) * t543 - t570) * t501; (t549 * t558 + t482) * MDP(9) + (t341 * t375 + t343 * t374 + (t341 * t451 + t343 * t448) * t522 + (-t556 - t304 * t451 + (t341 * t448 - t343 * t451) * qJD(6)) * t452) * MDP(27) + (-t408 * t521 + t394 + (-t378 * t450 - t453 * t550) * t529) * MDP(22) + (t348 * t449 + t377 * t550) * MDP(30) + (t398 * t542 - t297 * t374 - t307 * t341 + (t493 * t448 - t539 * t451) * t377 + ((-qJD(3) * t448 - t446 * t518) * t377 - t446 * t546 + t288 + (-t297 * t448 + t446 * t341) * qJD(5)) * t449 + (t291 * t505 + t297 * t518 - qJD(3) * t341 + t290 * t448 - t446 * t304 + (-t446 * t545 + t291) * qJD(5)) * t452) * MDP(31) + (-t398 * t546 - t297 * t375 - t307 * t343 + (t539 * t448 + t493 * t451) * t377 + (-(qJD(3) * t451 - t446 * t519) * t377 - t446 * t542 - t287 + (-t297 * t451 + t446 * t343) * qJD(5)) * t449 + (-t292 * t505 - t297 * t519 - qJD(3) * t343 + t290 * t451 - t446 * t303 + (-t446 * t496 - t292) * qJD(5)) * t452) * MDP(32) + (t451 * t540 + (-t452 * t519 + t481) * t343) * MDP(26) + (-t369 * t429 + t332 + 0.2e1 * t410 + (-pkin(3) * qJD(2) + t367) * t506 - t534) * MDP(16) + (pkin(1) * t507 - t386 * t429 + t534) * MDP(10) + (t412 + ((-qJD(2) * t447 - t326 - t369) * t453 + (t346 - t370 + t491) * t450) * t529) * MDP(15) + (t514 * t429 + (t373 * t453 + t385 * t450) * t529 - t359) * MDP(13) + (t512 * t429 + (-t353 * t450 + t367 * t453) * t529 - t336) * MDP(17) + (-t502 + (-t508 + (t380 - t525) * t450) * t529) * MDP(21) + (t412 + ((-pkin(2) * qJD(2) - t364 + t386) * t453 + (-t368 - t387 - t530) * t450) * t529) * MDP(11) + (-pkin(2) * t381 - qJ(3) * t359 - t364 * t387 - t368 * t514 - t373 * t385) * MDP(14) + (qJ(3) * t335 - t326 * t512 - t336 * t447 + t346 * t515 - t353 * t367) * MDP(18) + (-t385 * t505 + t354 - t482) * MDP(12) - MDP(4) * t488 + (t303 * t449 + t481 * t377 + (t464 + t467) * t452) * MDP(28) + MDP(7) * t472 + t372 * MDP(6) + (t567 * t449 + t452 * t571) * MDP(20) + (t347 * t452 - t449 * t551) * MDP(19) + t549 * t570 + t408 * MDP(23) * t506 + (-t300 * t506 + t323 * t452 + t447 * t347 + (-t446 * t521 + t537) * t408 + t513 * t380 + (-qJD(3) * t408 + t562) * t449) * MDP(25) + (t299 * t506 + t447 * t348 + t513 * t378 + (t323 + (-qJD(5) * t446 + t349) * t408) * t449 + (t408 * t516 - t562) * t452) * MDP(24) + (-t304 * t449 + (t448 * t522 + t374) * t377 + (t468 - t569) * t452) * MDP(29); (t368 * t429 + t354 + t381) * MDP(14) + (-t346 * t429 + t353 * t506 + t336) * MDP(18) + t571 * MDP(24) + t567 * MDP(25) + (-t542 + t554) * MDP(31) + (t546 + t553) * MDP(32) + ((MDP(12) - MDP(17)) * t543 + t564 * t442) * t549 + (MDP(31) * t545 + MDP(32) * t496) * t377 + t564 * t428 + (MDP(11) + MDP(15)) * t372; MDP(15) * t472 - MDP(16) * t488 + (-t443 * t549 - t428) * MDP(17) + (t326 * t429 + t332 + t335) * MDP(18) + (-t502 - t378 * t429 + (-t450 * t525 - t508) * t529) * MDP(24) + (-t380 * t429 - t408 * t550 + t394) * MDP(25) + (-t452 * t304 + (-t451 * t429 - t448 * t550) * t377 + (t468 + t569) * t449) * MDP(31) + (-t540 + (t448 * t429 - t451 * t550) * t377 + (t464 - t467) * t449) * MDP(32); -t378 ^ 2 * MDP(20) + (t347 + t552) * MDP(21) + (-t348 + t551) * MDP(22) - MDP(23) * t486 + (t300 * t408 + t456) * MDP(24) + (t299 * t408 + t317 * t378 - t463) * MDP(25) + (t343 * t496 + t556) * MDP(26) + ((t303 - t555) * t451 + (-t343 * t377 - t304) * t448) * MDP(27) + (t377 * t496 + t546 - t553) * MDP(28) + (-t377 ^ 2 * t448 + t542 + t554) * MDP(29) + (-pkin(5) * t304 - t300 * t341 + t459 * t448 - t451 * t563) * MDP(31) + (-pkin(5) * t303 - t300 * t343 + t448 * t563 + t459 * t451) * MDP(32) + (MDP(19) * t378 + t380 * MDP(20) - t317 * MDP(24) - t377 * MDP(30) - t291 * MDP(31) + t292 * MDP(32)) * t380; -t341 ^ 2 * MDP(27) + (t538 + t555) * MDP(28) + t348 * MDP(30) + (t292 * t377 + t295) * MDP(31) + (t291 * t377 + t297 * t341) * MDP(32) + (MDP(26) * t341 + t343 * MDP(27) + t377 * MDP(29) - t297 * MDP(31)) * t343 + (MDP(29) * t461 - MDP(31) * t520 + MDP(32) * t494) * t451 + (t461 * MDP(28) + (-qJD(6) * t408 - t347) * MDP(29) + t494 * MDP(31) + (-t296 + t520) * MDP(32)) * t448;];
tauc  = t1;
