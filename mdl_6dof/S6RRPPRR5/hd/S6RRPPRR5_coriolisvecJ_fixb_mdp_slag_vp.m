% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR5
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
%   see S6RRPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:11:15
% EndTime: 2019-03-09 09:11:24
% DurationCPUTime: 6.37s
% Computational Cost: add. (3437->535), mult. (8958->713), div. (0->0), fcn. (6386->8), ass. (0->208)
t448 = cos(pkin(6));
t518 = t448 * qJD(1);
t428 = qJD(2) + t518;
t454 = cos(qJ(5));
t451 = sin(qJ(5));
t452 = sin(qJ(2));
t447 = sin(pkin(6));
t530 = qJD(1) * t447;
t508 = t452 * t530;
t489 = t451 * t508;
t567 = t454 * t428 + t489;
t370 = qJD(6) + t567;
t570 = t370 ^ 2;
t374 = -t428 * t451 + t454 * t508;
t455 = cos(qJ(2));
t529 = qJD(1) * t455;
t507 = t447 * t529;
t402 = qJD(5) + t507;
t554 = t374 * t402;
t555 = t567 * t402;
t547 = t447 * t455;
t562 = pkin(1) * t452;
t470 = pkin(8) * t547 + t448 * t562;
t566 = t470 * qJD(2);
t377 = qJD(1) * t566;
t548 = t447 * t452;
t423 = qJD(4) * t548;
t513 = qJD(1) * qJD(2);
t502 = t447 * t513;
t487 = t455 * t502;
t537 = qJ(4) * t487 + qJD(1) * t423;
t335 = t377 - t537;
t526 = qJD(2) * t455;
t505 = t447 * t526;
t347 = -qJ(4) * t505 - t423 + t566;
t445 = t452 ^ 2;
t446 = t455 ^ 2;
t569 = (t445 - t446) * MDP(5);
t379 = -t447 * pkin(1) - pkin(2) * t547 - qJ(3) * t548;
t360 = pkin(3) * t547 - t379;
t474 = (pkin(4) * t455 - pkin(9) * t452) * t447;
t333 = t474 + t360;
t378 = t448 * qJ(3) + t470;
t359 = -qJ(4) * t547 + t378;
t349 = -pkin(9) * t448 + t359;
t568 = t451 * t333 + t454 * t349;
t512 = pkin(1) * t518;
t421 = t455 * t512;
t382 = -pkin(8) * t508 + t421;
t516 = qJD(3) - t382;
t367 = -pkin(1) * t530 - pkin(2) * t507 - qJ(3) * t508;
t348 = pkin(3) * t507 + qJD(4) - t367;
t322 = qJD(1) * t474 + t348;
t383 = pkin(8) * t507 + t452 * t512;
t364 = -qJ(4) * t507 + t383;
t406 = t428 * qJ(3);
t344 = t406 + t364;
t330 = -pkin(9) * t428 + t344;
t302 = t322 * t451 + t330 * t454;
t456 = -pkin(2) - pkin(3);
t444 = pkin(4) - t456;
t464 = t447 * (-pkin(9) * t455 - t444 * t452);
t460 = qJD(2) * t464;
t424 = qJD(3) * t548;
t536 = qJ(3) * t487 + qJD(1) * t424;
t316 = qJD(1) * t460 + t536;
t488 = t452 * t502;
t397 = qJ(4) * t488;
t401 = qJD(2) * t421;
t403 = t428 * qJD(3);
t525 = qJD(4) * t455;
t528 = qJD(2) * t452;
t323 = t397 + t401 + t403 + (-pkin(8) * t528 - t525) * t530;
t459 = -qJD(5) * t302 + t454 * t316 - t451 * t323;
t292 = pkin(5) * t488 - t459;
t565 = t370 * (pkin(5) * t374 + pkin(10) * t370) + t292;
t533 = qJ(3) * t505 + t424;
t321 = t460 + t533;
t563 = pkin(1) * t448;
t511 = qJD(2) * t563;
t422 = t455 * t511;
t436 = t448 * qJD(3);
t336 = t422 + t436 + (-t525 + (-pkin(8) + qJ(4)) * t528) * t447;
t564 = -qJD(5) * t568 + t321 * t454 - t336 * t451;
t561 = pkin(2) * qJD(2);
t560 = pkin(8) * qJD(2);
t450 = sin(qJ(6));
t465 = -qJD(6) * t374 - t488;
t523 = qJD(5) * t454;
t538 = t428 * t523 - t454 * t487;
t342 = -qJD(5) * t489 - t538;
t453 = cos(qJ(6));
t520 = qJD(6) * t453;
t541 = t453 * t342 + t402 * t520;
t305 = t450 * t465 + t541;
t559 = t305 * t450;
t338 = t374 * t450 - t453 * t402;
t558 = t338 * t370;
t557 = t348 * t452;
t556 = t367 * t452;
t497 = t370 * t453;
t553 = t382 * t428;
t552 = t383 * t428;
t551 = t402 * t451;
t550 = t402 * t454;
t443 = t447 ^ 2;
t549 = t443 * qJD(1) ^ 2;
t466 = t451 * t526 + t452 * t523;
t462 = t466 * t447;
t524 = qJD(5) * t451;
t504 = t428 * t524;
t343 = qJD(1) * t462 - t504;
t546 = t450 * t343;
t545 = t451 * t452;
t544 = t453 * t343;
t543 = t454 * t455;
t542 = -t364 + t402 * (pkin(5) * t451 - pkin(10) * t454);
t410 = qJ(3) * t507;
t331 = qJD(1) * t464 + t410;
t407 = qJ(4) * t508;
t362 = t407 + t382;
t539 = t451 * t331 + t454 * t362;
t535 = 0.2e1 * t403 + t401;
t531 = qJ(3) * qJD(2);
t527 = qJD(2) * t454;
t300 = pkin(10) * t402 + t302;
t522 = qJD(6) * t300;
t521 = qJD(6) * t450;
t519 = t402 * MDP(24);
t517 = qJD(3) - t362;
t514 = MDP(25) * qJD(5);
t510 = pkin(8) * t526;
t509 = t455 * t549;
t380 = -t448 * pkin(2) + pkin(8) * t548 - t455 * t563;
t506 = t447 * t528;
t503 = t443 * t513;
t491 = t456 * t548;
t363 = qJD(1) * t491 + t410;
t501 = -t363 - t560;
t381 = pkin(2) * t508 - t410;
t500 = t381 - t560;
t468 = t451 * t316 + t322 * t523 + t454 * t323 - t330 * t524;
t291 = -pkin(10) * t488 + t468;
t358 = -t428 * pkin(2) + t516;
t328 = -t428 * pkin(3) + t358 - t407;
t319 = pkin(4) * t428 - t328;
t303 = pkin(5) * t567 - pkin(10) * t374 + t319;
t495 = -qJD(6) * t303 - t291;
t393 = pkin(5) * t454 + pkin(10) * t451 + t444;
t494 = -pkin(10) * t508 - qJD(6) * t393 + t539;
t493 = t452 * t511;
t490 = t452 * t509;
t486 = -0.2e1 * pkin(1) * t503;
t368 = (t450 * t543 + t452 * t453) * t530;
t484 = t450 * t523 + t368;
t369 = (-t450 * t452 + t453 * t543) * t530;
t483 = -t453 * t523 - t369;
t350 = -t448 * pkin(3) - qJ(4) * t548 + t380;
t482 = qJD(2) * t491;
t294 = t300 * t453 + t303 * t450;
t293 = -t300 * t450 + t303 * t453;
t308 = pkin(10) * t547 + t568;
t341 = t448 * pkin(4) - t350;
t386 = t447 * t545 + t448 * t454;
t387 = -t448 * t451 + t454 * t548;
t313 = pkin(5) * t386 - pkin(10) * t387 + t341;
t481 = t308 * t453 + t313 * t450;
t480 = -t308 * t450 + t313 * t453;
t301 = t322 * t454 - t330 * t451;
t478 = t333 * t454 - t349 * t451;
t340 = t374 * t453 + t402 * t450;
t476 = -t377 * t448 - t428 * t566;
t475 = -pkin(8) * t506 + t422;
t473 = -t387 * t450 + t453 * t547;
t357 = t387 * t453 + t450 * t547;
t472 = t370 * t520 + t546;
t471 = t370 * t521 - t544;
t469 = -pkin(8) * t488 + t401;
t467 = t451 * t321 + t333 * t523 + t454 * t336 - t349 * t524;
t299 = -pkin(5) * t402 - t301;
t461 = -pkin(10) * t343 + (t299 + t301) * t370;
t449 = qJ(3) - pkin(9);
t458 = (-t319 * t455 + t449 * t528) * t530 - t319 * qJD(5);
t394 = t451 * t488;
t388 = t428 * t507;
t371 = t436 + t475;
t366 = -t388 + t487;
t365 = pkin(2) * t506 - t533;
t361 = t406 + t383;
t355 = -qJD(5) * t386 + t454 * t505;
t354 = -t448 * t524 + t462;
t353 = t403 + t469;
t351 = pkin(2) * t488 - t536;
t346 = t482 + t533;
t334 = qJD(1) * t482 + t536;
t312 = qJD(6) * t473 + t355 * t453 - t450 * t506;
t311 = qJD(6) * t357 + t355 * t450 + t453 * t506;
t309 = pkin(5) * t508 - t331 * t454 + t362 * t451;
t307 = -pkin(5) * t547 - t478;
t306 = qJD(6) * t340 + t342 * t450 + t453 * t488;
t304 = pkin(5) * t354 - pkin(10) * t355 - t347;
t298 = pkin(5) * t343 - pkin(10) * t342 - t335;
t297 = t453 * t298;
t296 = pkin(5) * t506 - t564;
t295 = -pkin(10) * t506 + t467;
t290 = -qJD(6) * t294 - t291 * t450 + t297;
t289 = qJD(6) * t293 + t291 * t453 + t298 * t450;
t1 = [(-t428 * t475 - t448 * t469 + t455 * t486) * MDP(10) + (t452 * t486 + t476) * MDP(9) + (-(qJD(6) * t480 + t295 * t453 + t304 * t450) * t370 - t481 * t343 - t289 * t386 - t294 * t354 + t296 * t340 + t307 * t305 + t292 * t357 + t299 * t312) * MDP(32) + (t305 * t386 + t312 * t370 + t340 * t354 + t343 * t357) * MDP(28) + (t343 * t386 + t354 * t370) * MDP(30) + (t342 * t387 + t355 * t374) * MDP(19) + (t305 * t357 + t312 * t340) * MDP(26) + (t323 * t359 + t328 * t347 + t334 * t360 + t335 * t350 + t336 * t344 + t346 * t348) * MDP(18) + (-t323 * t455 - t335 * t452 + (-t328 * t455 + t344 * t452) * qJD(2) + (-t336 * t455 - t347 * t452 + (-t350 * t455 + t359 * t452) * qJD(2)) * qJD(1)) * t447 * MDP(17) + (t353 * t455 + t377 * t452 + (t358 * t455 - t361 * t452) * qJD(2) + (t371 * t455 + t566 * t452 + (-t378 * t452 + t380 * t455) * qJD(2)) * qJD(1)) * t447 * MDP(12) + (-t402 * t447 - t443 * t529) * MDP(23) * t528 + (MDP(6) * t505 - MDP(7) * t506) * (t428 + t518) + 0.2e1 * (t452 * t455 * MDP(4) - t569) * t503 + ((-qJD(6) * t481 - t295 * t450 + t304 * t453) * t370 + t480 * t343 + t290 * t386 + t293 * t354 + t296 * t338 + t307 * t306 - t292 * t473 + t299 * t311) * MDP(31) + (-t306 * t386 - t311 * t370 - t338 * t354 + t343 * t473) * MDP(29) + (t305 * t473 - t306 * t357 - t311 * t340 - t312 * t338) * MDP(27) + (-t342 * t386 - t343 * t387 - t354 * t374 - t355 * t567) * MDP(20) + (t564 * t402 - t347 * t567 + t341 * t343 - t335 * t386 + t319 * t354 + (t459 * t455 + (-qJD(1) * t478 - t301) * t528) * t447) * MDP(24) + (-t354 * t402 + (-t343 * t455 + (qJD(1) * t386 + t567) * t528) * t447) * MDP(22) + (t351 * t379 + t353 * t378 + t358 * t566 + t361 * t371 + t365 * t367 + t377 * t380) * MDP(14) + (-t467 * t402 - t347 * t374 + t341 * t342 - t335 * t387 + t319 * t355 + (-t468 * t455 + (qJD(1) * t568 + t302) * t528) * t447) * MDP(25) + (t355 * t402 + (t342 * t455 + (-qJD(1) * t387 - t374) * t528) * t447) * MDP(21) + (-t335 * t448 - t347 * t428 + (-t348 * t528 + t334 * t455 + (t346 * t455 - t360 * t528) * qJD(1)) * t447) * MDP(15) + ((t367 * t528 - t351 * t455 + (-t365 * t455 + t379 * t528) * qJD(1)) * t447 + t476) * MDP(11) + (t353 * t448 + t371 * t428 + (-t367 * t526 - t351 * t452 + (-t365 * t452 - t379 * t526) * qJD(1)) * t447) * MDP(13) + (t323 * t448 + t336 * t428 + (t348 * t526 + t334 * t452 + (t346 * t452 + t360 * t526) * qJD(1)) * t447) * MDP(16); (t305 * t454 + t483 * t370 + (-t340 * t402 + t471) * t451) * MDP(28) + (-t306 * t454 + t484 * t370 + (t338 * t402 + t472) * t451) * MDP(29) + t366 * MDP(6) + (t335 * t451 + t444 * t342 + t364 * t374 + (t449 * t524 + t539) * t402 + (-qJD(3) * t402 + t458) * t454) * MDP(25) + (-t302 * MDP(25) + t301 * MDP(24) + (-qJD(2) + t428) * MDP(7) + t402 * MDP(23)) * t508 + (t552 + (-t493 + (t455 * t500 - t556) * t447) * qJD(1)) * MDP(11) + (t364 * t428 + (-t493 + (t455 * t501 + t557) * t447) * qJD(1) + t537) * MDP(15) + (-pkin(2) * t377 + qJ(3) * t353 - t358 * t383 + t361 * t516 - t367 * t381) * MDP(14) + (qJ(3) * t323 - t328 * t364 + t335 * t456 + t344 * t517 - t348 * t363) * MDP(18) + (t393 * t544 - t299 * t368 - t309 * t338 + (t494 * t450 - t542 * t453) * t370 + ((-qJD(3) * t450 - t449 * t520) * t370 - t449 * t546 + t290 + (-t299 * t450 + t449 * t338) * qJD(5)) * t454 + (-t293 * t507 - t299 * t520 + qJD(3) * t338 - t292 * t450 + t449 * t306 + (t449 * t450 * t370 - t293) * qJD(5)) * t451) * MDP(31) + (-t342 * t451 - t374 * t550) * MDP(19) + (t343 * t454 - t370 * t551) * MDP(30) + (pkin(1) * t509 - t469 + t553) * MDP(10) + ((-t342 + t555) * t454 + (t343 + t554) * t451) * MDP(20) + (-t393 * t546 - t299 * t369 - t309 * t340 + (t542 * t450 + t494 * t453) * t370 + (-(qJD(3) * t453 - t449 * t521) * t370 - t449 * t544 - t289 + (-t299 * t453 + t449 * t340) * qJD(5)) * t454 + (t294 * t507 + t299 * t521 + qJD(3) * t340 - t292 * t453 + t449 * t305 + (t449 * t497 + t294) * qJD(5)) * t451) * MDP(32) + (-t305 * t451 * t453 + (t451 * t521 + t483) * t340) * MDP(26) - MDP(4) * t490 + (-pkin(8) * t487 + t552 + (-t448 * t513 + t549) * t562) * MDP(9) + (t338 * t369 + t340 * t368 + (t338 * t453 + t340 * t450) * t523 + (t559 + t306 * t453 + (-t338 * t450 + t340 * t453) * qJD(6)) * t451) * MDP(27) + t549 * t569 + (t444 * t343 + t364 * t567 + (-t335 + (-qJD(5) * t449 - t331) * t402) * t454 + (-t402 * t517 + t458) * t451) * MDP(24) + (((-qJD(4) - t348) * t455 + t501 * t452) * MDP(16) + (t374 * t452 - t402 * t543) * MDP(21) + ((t361 - t383 - t531) * t452 + (-t358 + t516 - t561) * t455) * MDP(12) + ((-t344 + t364 + t531) * t452 + (-qJD(2) * t456 + t328 - t517) * t455) * MDP(17) + (t455 * t551 + (-t567 + t527) * t452) * MDP(22) + (t367 * t455 + t452 * t500) * MDP(13)) * t530 + t402 * t524 * MDP(22) + (-t362 * t428 + t397 + t535) * MDP(16) + (-t402 * t523 + t394) * MDP(21) + (-t553 + t535) * MDP(13); -t361 * t428 * MDP(14) + (-t344 * t428 - t537) * MDP(18) + (t504 - t554) * MDP(24) + (t538 + t555) * MDP(25) + (t338 * t374 - t544) * MDP(31) + (t340 * t374 + t546) * MDP(32) + (-MDP(11) - MDP(15)) * t490 + (MDP(13) + MDP(16)) * (-t428 ^ 2 - t445 * t549) + (MDP(12) - MDP(17)) * t366 + ((MDP(14) + MDP(18)) * t493 + ((t510 + t556) * MDP(14) + (t510 - t557) * MDP(18) - t466 * MDP(24) + t514 * t545) * t447) * qJD(1) + (MDP(31) * t450 + MDP(32) * t453) * t570; t388 * MDP(16) + t536 * MDP(18) + t394 * MDP(25) + (-t306 * MDP(31) - t305 * MDP(32) - t402 * t514) * t454 + (-t445 - t446) * MDP(17) * t549 + (-MDP(31) * t484 + MDP(32) * t483) * t370 + (-qJD(5) * t519 + (qJD(5) * t338 - t472) * MDP(31) + (qJD(5) * t340 + t471) * MDP(32)) * t451 + ((-MDP(25) * t550 + qJD(2) * MDP(16) + t344 * MDP(18) + (t338 * MDP(31) + t340 * MDP(32) - t519) * t451) * t455 + ((-qJD(2) - t428) * MDP(15) + (-pkin(3) * qJD(2) + t328 - t561) * MDP(18) + (-t567 - t527) * MDP(24) - t374 * MDP(25)) * t452) * t530; -t567 ^ 2 * MDP(20) + (t342 + t555) * MDP(21) + (-t466 * t530 + t504 + t554) * MDP(22) - MDP(23) * t488 + (t302 * t402 + t459) * MDP(24) + (t301 * t402 + t319 * t567 - t468) * MDP(25) + (t340 * t497 + t559) * MDP(26) + ((t305 - t558) * t453 + (-t340 * t370 - t306) * t450) * MDP(27) + (t370 * t497 + t546) * MDP(28) + (-t450 * t570 + t544) * MDP(29) + (-pkin(5) * t306 - t302 * t338 + t461 * t450 - t453 * t565) * MDP(31) + (-pkin(5) * t305 - t302 * t340 + t450 * t565 + t461 * t453) * MDP(32) + (MDP(19) * t567 + t374 * MDP(20) - t319 * MDP(24) - t340 * MDP(28) + t338 * MDP(29) - t370 * MDP(30) - t293 * MDP(31) + t294 * MDP(32)) * t374; -t338 ^ 2 * MDP(27) + (t541 + t558) * MDP(28) + t343 * MDP(30) + (t294 * t370 + t297) * MDP(31) + (t293 * t370 + t299 * t338) * MDP(32) + (MDP(26) * t338 + t340 * MDP(27) + t370 * MDP(29) - t299 * MDP(31)) * t340 + (MDP(29) * t465 - MDP(31) * t522 + MDP(32) * t495) * t453 + (t465 * MDP(28) + (-qJD(6) * t402 - t342) * MDP(29) + t495 * MDP(31) + (-t298 + t522) * MDP(32)) * t450;];
tauc  = t1;
