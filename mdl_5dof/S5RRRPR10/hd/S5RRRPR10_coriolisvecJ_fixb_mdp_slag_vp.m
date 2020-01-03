% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:30:03
% EndTime: 2019-12-31 21:30:15
% DurationCPUTime: 6.46s
% Computational Cost: add. (5263->419), mult. (14283->597), div. (0->0), fcn. (11121->10), ass. (0->193)
t457 = sin(pkin(5));
t465 = cos(qJ(2));
t529 = qJD(1) * t465;
t511 = t457 * t529;
t569 = qJD(3) - t511;
t461 = sin(qJ(3));
t462 = sin(qJ(2));
t464 = cos(qJ(3));
t531 = qJD(1) * t457;
t512 = t462 * t531;
t459 = cos(pkin(5));
t530 = qJD(1) * t459;
t518 = pkin(1) * t530;
t413 = -pkin(7) * t512 + t465 * t518;
t477 = (pkin(2) * t462 - pkin(8) * t465) * t457;
t414 = qJD(1) * t477;
t496 = -t413 * t461 + t464 * t414;
t554 = -qJ(4) - pkin(8);
t503 = qJD(3) * t554;
t539 = t464 * t465;
t568 = -(pkin(3) * t462 - qJ(4) * t539) * t531 - t496 - qJD(4) * t461 + t464 * t503;
t490 = t461 * t511;
t534 = t464 * t413 + t461 * t414;
t567 = -qJ(4) * t490 - qJD(4) * t464 - t461 * t503 + t534;
t463 = cos(qJ(5));
t460 = sin(qJ(5));
t447 = qJD(2) + t530;
t491 = t461 * t512;
t396 = -t464 * t447 + t491;
t398 = t447 * t461 + t464 * t512;
t456 = sin(pkin(10));
t458 = cos(pkin(10));
t479 = -t396 * t456 + t458 * t398;
t549 = t479 * t460;
t342 = -t463 * t569 + t549;
t498 = -t458 * t396 - t398 * t456;
t558 = qJD(5) - t498;
t566 = t342 * t558;
t344 = t460 * t569 + t463 * t479;
t565 = t344 * t558;
t425 = t456 * t461 - t458 * t464;
t564 = t569 * t425;
t426 = t456 * t464 + t458 * t461;
t533 = t569 * t426;
t495 = t558 * t463;
t519 = qJD(1) * qJD(2);
t504 = t457 * t519;
t488 = t465 * t504;
t524 = qJD(3) * t464;
t372 = -qJD(3) * t491 + t447 * t524 + t464 * t488;
t526 = qJD(2) * t465;
t508 = t461 * t526;
t525 = qJD(3) * t461;
t373 = (t462 * t524 + t508) * t531 + t447 * t525;
t338 = t372 * t456 + t458 * t373;
t541 = t460 * t338;
t563 = -t558 * t495 - t541;
t453 = t457 ^ 2;
t562 = -0.2e1 * t453 * t519;
t560 = (t462 ^ 2 - t465 ^ 2) * MDP(5);
t538 = t567 * t456 + t568 * t458;
t536 = t568 * t456 - t567 * t458;
t542 = t457 * t465;
t556 = pkin(1) * t462;
t410 = pkin(7) * t542 + (pkin(8) + t556) * t459;
t411 = (-pkin(2) * t465 - pkin(8) * t462 - pkin(1)) * t457;
t535 = t464 * t410 + t461 * t411;
t444 = t462 * t518;
t416 = pkin(7) * t511 + t444;
t559 = -t416 + (-t490 + t525) * pkin(3);
t385 = pkin(8) * t447 + t416;
t392 = qJD(1) * t411;
t354 = t385 * t464 + t392 * t461;
t415 = qJD(2) * t477;
t405 = qJD(1) * t415;
t543 = t457 * t462;
t448 = pkin(7) * t543;
t555 = pkin(1) * t465;
t417 = (t459 * t555 - t448) * qJD(2);
t406 = qJD(1) * t417;
t467 = -t354 * qJD(3) + t464 * t405 - t461 * t406;
t489 = t462 * t504;
t304 = pkin(3) * t489 - qJ(4) * t372 - qJD(4) * t398 + t467;
t476 = -t385 * t525 + t392 * t524 + t461 * t405 + t464 * t406;
t308 = -qJ(4) * t373 - qJD(4) * t396 + t476;
t288 = t304 * t458 - t308 * t456;
t286 = -pkin(4) * t489 - t288;
t451 = pkin(3) * t456 + pkin(9);
t557 = (pkin(3) * t398 + pkin(4) * t479 - pkin(9) * t498 + qJD(5) * t451) * t558 + t286;
t339 = t372 * t458 - t373 * t456;
t522 = qJD(5) * t463;
t515 = t463 * t339 + t460 * t489 + t522 * t569;
t523 = qJD(5) * t460;
t309 = -t479 * t523 + t515;
t553 = t309 * t460;
t337 = -qJ(4) * t396 + t354;
t552 = t337 * t456;
t551 = t342 * t479;
t550 = t344 * t479;
t548 = t396 * t569;
t547 = t398 * t569;
t546 = t426 * t463;
t545 = t569 * t464;
t466 = qJD(1) ^ 2;
t544 = t453 * t466;
t332 = t458 * t337;
t540 = t461 * t569;
t334 = t463 * t338;
t289 = t456 * t304 + t458 * t308;
t422 = -t459 * t464 + t461 * t543;
t509 = t457 * t526;
t377 = -t422 * qJD(3) + t464 * t509;
t423 = t459 * t461 + t464 * t543;
t468 = -t535 * qJD(3) + t464 * t415 - t417 * t461;
t528 = qJD(2) * t462;
t510 = t457 * t528;
t315 = pkin(3) * t510 - qJ(4) * t377 - qJD(4) * t423 + t468;
t376 = t423 * qJD(3) + t457 * t508;
t475 = -t410 * t525 + t411 * t524 + t461 * t415 + t464 * t417;
t319 = -qJ(4) * t376 - qJD(4) * t422 + t475;
t295 = t456 * t315 + t458 * t319;
t353 = -t385 * t461 + t464 * t392;
t336 = -qJ(4) * t398 + t353;
t330 = pkin(3) * t569 + t336;
t303 = t456 * t330 + t332;
t497 = -t410 * t461 + t464 * t411;
t341 = -pkin(3) * t542 - qJ(4) * t423 + t497;
t349 = -qJ(4) * t422 + t535;
t318 = t456 * t341 + t458 * t349;
t537 = pkin(4) * t512 - t538;
t407 = pkin(7) * t488 + qJD(2) * t444;
t418 = t459 * pkin(1) * t528 + pkin(7) * t509;
t527 = qJD(2) * t464;
t521 = qJD(2) - t447;
t516 = t460 * t542;
t513 = -pkin(3) * t464 - pkin(2);
t506 = t554 * t461;
t287 = pkin(9) * t489 + t289;
t351 = pkin(3) * t373 + t407;
t297 = pkin(4) * t338 - pkin(9) * t339 + t351;
t502 = -t287 * t460 + t463 * t297;
t501 = t339 * t460 - t463 * t489;
t500 = t564 * t460 - t463 * t512;
t499 = t460 * t512 + t564 * t463;
t492 = t453 * t462 * t465 * MDP(4);
t487 = pkin(3) * t376 + t418;
t486 = pkin(1) * t562;
t374 = pkin(4) * t425 - pkin(9) * t426 + t513;
t485 = pkin(9) * t512 - qJD(5) * t374 - t536;
t441 = t554 * t464;
t383 = -t458 * t441 + t456 * t506;
t484 = -t533 * pkin(4) - t564 * pkin(9) + qJD(5) * t383 - t559;
t483 = t287 * t463 + t297 * t460;
t299 = pkin(9) * t569 + t303;
t384 = -pkin(2) * t447 - t413;
t359 = pkin(3) * t396 + qJD(4) + t384;
t311 = -pkin(4) * t498 - pkin(9) * t479 + t359;
t291 = t299 * t463 + t311 * t460;
t482 = t299 * t460 - t311 * t463;
t314 = -pkin(9) * t542 + t318;
t369 = t458 * t422 + t423 * t456;
t370 = -t422 * t456 + t423 * t458;
t409 = t448 + (-pkin(2) - t555) * t459;
t470 = pkin(3) * t422 + t409;
t327 = pkin(4) * t369 - pkin(9) * t370 + t470;
t481 = t314 * t463 + t327 * t460;
t480 = -t314 * t460 + t327 * t463;
t294 = t315 * t458 - t319 * t456;
t302 = t330 * t458 - t552;
t317 = t341 * t458 - t349 * t456;
t478 = t334 + (t460 * t498 - t523) * t558;
t356 = t370 * t460 + t463 * t542;
t474 = t426 * t522 - t500;
t473 = -t426 * t523 - t499;
t298 = -pkin(4) * t569 - t302;
t307 = t336 * t458 - t552;
t469 = -t451 * t338 + (t298 + t307) * t558;
t452 = -pkin(3) * t458 - pkin(4);
t382 = -t441 * t456 - t458 * t506;
t357 = t370 * t463 - t516;
t348 = -t376 * t456 + t377 * t458;
t347 = t458 * t376 + t377 * t456;
t321 = -qJD(5) * t516 + t348 * t460 + t370 * t522 - t463 * t510;
t320 = -t356 * qJD(5) + t348 * t463 + t460 * t510;
t313 = pkin(4) * t542 - t317;
t310 = t344 * qJD(5) + t501;
t306 = t336 * t456 + t332;
t300 = pkin(4) * t347 - pkin(9) * t348 + t487;
t293 = pkin(9) * t510 + t295;
t292 = -pkin(4) * t510 - t294;
t285 = -t291 * qJD(5) + t502;
t284 = -t482 * qJD(5) + t483;
t1 = [0.2e1 * t492 * t519 + t560 * t562 + (-t407 * t459 - t418 * t447 + t462 * t486) * MDP(9) + (-t406 * t459 - t417 * t447 + t465 * t486) * MDP(10) + (t372 * t423 + t377 * t398) * MDP(11) + (-t372 * t422 - t373 * t423 - t376 * t398 - t377 * t396) * MDP(12) + (t377 * t569 + (-t372 * t465 + (qJD(1) * t423 + t398) * t528) * t457) * MDP(13) + (-t376 * t569 + (t373 * t465 + (-qJD(1) * t422 - t396) * t528) * t457) * MDP(14) + (-t453 * t529 + t457 * t569) * MDP(15) * t528 + (t468 * t569 + t418 * t396 + t409 * t373 + t407 * t422 + t384 * t376 + (-t467 * t465 + (t497 * qJD(1) + t353) * t528) * t457) * MDP(16) + (-t475 * t569 + t418 * t398 + t409 * t372 + t407 * t423 + t384 * t377 + (t476 * t465 + (-t535 * qJD(1) - t354) * t528) * t457) * MDP(17) + (-t288 * t370 - t289 * t369 - t294 * t479 + t295 * t498 - t302 * t348 - t303 * t347 - t317 * t339 - t318 * t338) * MDP(18) + (t288 * t317 + t289 * t318 + t302 * t294 + t303 * t295 + t351 * t470 + t359 * t487) * MDP(19) + (t309 * t357 + t320 * t344) * MDP(20) + (-t309 * t356 - t310 * t357 - t320 * t342 - t321 * t344) * MDP(21) + (t309 * t369 + t320 * t558 + t338 * t357 + t344 * t347) * MDP(22) + (-t310 * t369 - t321 * t558 - t338 * t356 - t342 * t347) * MDP(23) + (t338 * t369 + t347 * t558) * MDP(24) + ((-t481 * qJD(5) - t293 * t460 + t300 * t463) * t558 + t480 * t338 + t285 * t369 - t482 * t347 + t292 * t342 + t313 * t310 + t286 * t356 + t298 * t321) * MDP(25) + (-(t480 * qJD(5) + t293 * t463 + t300 * t460) * t558 - t481 * t338 - t284 * t369 - t291 * t347 + t292 * t344 + t313 * t309 + t286 * t357 + t298 * t320) * MDP(26) + (MDP(6) * t509 - MDP(7) * t510) * (t447 + t530); -t466 * t492 + t544 * t560 + t521 * MDP(6) * t511 + (t416 * t447 + t544 * t556 - t407) * MDP(9) + (pkin(7) * t489 + t413 * t447 + (-t459 * t519 + t544) * t555) * MDP(10) + (t372 * t461 + t398 * t545) * MDP(11) + ((t372 - t548) * t464 + (-t373 - t547) * t461) * MDP(12) + (t569 * t524 + (-t569 * t539 + (qJD(2) * t461 - t398) * t462) * t531) * MDP(13) + (-t569 * t525 + (t465 * t540 + (t396 + t527) * t462) * t531) * MDP(14) + (-pkin(2) * t373 - t407 * t464 - t496 * t569 - t416 * t396 + (-pkin(8) * t545 + t384 * t461) * qJD(3) + (-t353 * t462 + (-pkin(8) * t528 - t384 * t465) * t461) * t531) * MDP(16) + (-pkin(2) * t372 + t407 * t461 + t534 * t569 - t416 * t398 + (pkin(8) * t540 + t384 * t464) * qJD(3) + (-t384 * t539 + (-pkin(8) * t527 + t354) * t462) * t531) * MDP(17) + (-t288 * t426 - t289 * t425 + t564 * t302 - t533 * t303 - t383 * t338 + t339 * t382 - t538 * t479 + t536 * t498) * MDP(18) + (-t288 * t382 + t289 * t383 + t538 * t302 + t536 * t303 + t351 * t513 + t559 * t359) * MDP(19) + (t309 * t546 + t473 * t344) * MDP(20) + (t500 * t344 + t499 * t342 + (-t553 - t310 * t463 + (t342 * t460 - t344 * t463) * qJD(5)) * t426) * MDP(21) + (t309 * t425 + t426 * t334 + t533 * t344 + t473 * t558) * MDP(22) + (-t310 * t425 - t533 * t342 - t426 * t541 - t474 * t558) * MDP(23) + (t338 * t425 + t533 * t558) * MDP(24) + ((t374 * t463 - t383 * t460) * t338 + t285 * t425 + t382 * t310 + t286 * t460 * t426 + (t485 * t460 - t484 * t463) * t558 + t537 * t342 - t533 * t482 + t474 * t298) * MDP(25) + (-(t374 * t460 + t383 * t463) * t338 - t284 * t425 + t382 * t309 + t286 * t546 + (t484 * t460 + t485 * t463) * t558 + t537 * t344 - t533 * t291 + t473 * t298) * MDP(26) + (-MDP(15) * t569 - t521 * MDP(7)) * t512; t398 * t396 * MDP(11) + (-t396 ^ 2 + t398 ^ 2) * MDP(12) + (t372 + t548) * MDP(13) + (-t373 + t547) * MDP(14) + MDP(15) * t489 + (t354 * t569 - t384 * t398 + t467) * MDP(16) + (t353 * t569 + t384 * t396 - t476) * MDP(17) + ((-t338 * t456 - t339 * t458) * pkin(3) + (t302 - t307) * t498 + (t303 - t306) * t479) * MDP(18) + (t302 * t306 - t303 * t307 + (t288 * t458 + t289 * t456 - t359 * t398) * pkin(3)) * MDP(19) + (t344 * t495 + t553) * MDP(20) + ((t309 - t566) * t463 + (-t310 - t565) * t460) * MDP(21) + (-t550 - t563) * MDP(22) + (t478 + t551) * MDP(23) - t558 * t479 * MDP(24) + (-t306 * t342 + t452 * t310 + t469 * t460 - t557 * t463 + t479 * t482) * MDP(25) + (t291 * t479 - t306 * t344 + t452 * t309 + t557 * t460 + t469 * t463) * MDP(26); (-t479 ^ 2 - t498 ^ 2) * MDP(18) + (t302 * t479 - t303 * t498 + t351) * MDP(19) + (t478 - t551) * MDP(25) + (-t550 + t563) * MDP(26); t344 * t342 * MDP(20) + (-t342 ^ 2 + t344 ^ 2) * MDP(21) + (t515 + t566) * MDP(22) + (-t501 + t565) * MDP(23) + t338 * MDP(24) + (t291 * t558 - t298 * t344 + t502) * MDP(25) + (t298 * t342 - t482 * t558 - t483) * MDP(26) + (-MDP(22) * t549 - MDP(23) * t344 - MDP(25) * t291 + MDP(26) * t482) * qJD(5);];
tauc = t1;
