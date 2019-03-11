% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:47:03
% EndTime: 2019-03-08 22:47:13
% DurationCPUTime: 5.56s
% Computational Cost: add. (4838->447), mult. (12103->608), div. (0->0), fcn. (8707->10), ass. (0->200)
t465 = sin(qJ(3));
t468 = cos(qJ(3));
t484 = pkin(3) * t465 - pkin(9) * t468;
t429 = t484 * qJD(3);
t435 = -pkin(3) * t468 - pkin(9) * t465 - pkin(2);
t464 = sin(qJ(4));
t466 = sin(qJ(2));
t467 = cos(qJ(4));
t516 = qJD(4) * t467;
t461 = sin(pkin(6));
t526 = qJD(1) * t461;
t469 = cos(qJ(2));
t546 = t468 * t469;
t581 = -(t464 * t466 + t467 * t546) * t526 + t464 * t429 + t435 * t516;
t522 = qJD(2) * t468;
t499 = t464 * t522;
t518 = qJD(4) * t464;
t580 = -t499 + t518;
t520 = qJD(3) * t465;
t568 = pkin(8) * t464;
t579 = t467 * t429 + t520 * t568 - (-t464 * t546 + t466 * t467) * t526;
t547 = t467 * t468;
t450 = pkin(8) * t547;
t478 = pkin(4) * t465 - qJ(5) * t547;
t515 = qJD(5) * t467;
t578 = -t465 * t515 + t478 * qJD(3) + (-t450 + (qJ(5) * t465 - t435) * t464) * qJD(4) + t579;
t549 = t465 * t467;
t577 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t549 + (-qJD(5) * t465 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t468) * t464 + t581;
t513 = t467 * qJD(3);
t523 = qJD(2) * t465;
t425 = t464 * t523 - t513;
t521 = qJD(3) * t464;
t427 = t467 * t523 + t521;
t460 = sin(pkin(11));
t462 = cos(pkin(11));
t369 = t462 * t425 + t427 * t460;
t447 = -qJD(4) + t522;
t576 = t369 * t447;
t421 = t460 * t467 + t462 * t464;
t408 = t421 * qJD(4);
t534 = t421 * t522 - t408;
t554 = t462 * t467;
t533 = t460 * t580 - t462 * t516 + t522 * t554;
t519 = qJD(3) * t468;
t494 = t464 * t519;
t495 = t465 * t516;
t575 = t494 + t495;
t480 = -t425 * t460 + t462 * t427;
t574 = t480 ^ 2;
t573 = MDP(5) * t465;
t458 = t465 ^ 2;
t572 = MDP(6) * (-t468 ^ 2 + t458);
t504 = t466 * t526;
t431 = qJD(2) * pkin(8) + t504;
t463 = cos(pkin(6));
t525 = qJD(1) * t468;
t397 = -t465 * t431 + t463 * t525;
t389 = -qJD(3) * pkin(3) - t397;
t362 = pkin(4) * t425 + qJD(5) + t389;
t321 = pkin(5) * t369 - qJ(6) * t480 + t362;
t571 = t321 * t480;
t541 = -t577 * t460 + t578 * t462;
t540 = t578 * t460 + t577 * t462;
t428 = t484 * qJD(2);
t487 = -t397 * t464 + t467 * t428;
t350 = qJD(2) * t478 + t487;
t535 = t467 * t397 + t464 * t428;
t357 = -qJ(5) * t499 + t535;
t567 = -qJ(5) - pkin(9);
t489 = qJD(4) * t567;
t405 = t464 * t489 + t515;
t474 = -qJD(5) * t464 + t467 * t489;
t536 = (-t350 + t474) * t462 + (t357 - t405) * t460;
t553 = t463 * t465;
t445 = qJD(1) * t553;
t398 = t468 * t431 + t445;
t570 = pkin(4) * t580 - t398;
t569 = MDP(19) + MDP(22);
t566 = qJD(2) * pkin(2);
t390 = qJD(3) * pkin(9) + t398;
t503 = t469 * t526;
t399 = qJD(2) * t435 - t503;
t552 = t464 * t399;
t354 = t390 * t467 + t552;
t340 = -qJ(5) * t425 + t354;
t334 = t462 * t340;
t352 = -t390 * t464 + t467 * t399;
t339 = -qJ(5) * t427 + t352;
t313 = t339 * t460 + t334;
t565 = t313 * t480;
t564 = t340 * t460;
t524 = qJD(2) * t461;
t500 = t469 * t524;
t485 = t465 * t500;
t365 = qJD(1) * t485 + qJD(3) * t445 + t431 * t519;
t563 = t365 * t464;
t562 = t365 * t467;
t561 = t389 * t464;
t517 = qJD(4) * t465;
t496 = t464 * t517;
t511 = qJD(2) * qJD(3);
t491 = t468 * t511;
t510 = qJD(3) * qJD(4);
t530 = (t491 + t510) * t467;
t392 = -qJD(2) * t496 + t530;
t560 = t392 * t464;
t559 = t425 * t447;
t558 = t427 * t447;
t557 = t447 * t467;
t556 = t461 * t466;
t555 = t461 * t469;
t551 = t464 * t465;
t550 = t464 * t468;
t470 = qJD(3) ^ 2;
t548 = t465 * t470;
t545 = t468 * t470;
t364 = -t431 * t520 + (qJD(3) * t463 + t500) * t525;
t396 = (t429 + t504) * qJD(2);
t488 = t464 * t364 - t467 * t396;
t472 = -qJD(4) * t354 - t488;
t492 = t465 * t511;
t299 = pkin(4) * t492 - qJ(5) * t392 - qJD(5) * t427 + t472;
t490 = t464 * t510;
t393 = t575 * qJD(2) + t490;
t476 = t467 * t364 - t390 * t518 + t464 * t396 + t399 * t516;
t302 = -qJ(5) * t393 - qJD(5) * t425 + t476;
t544 = -t462 * t299 + t460 * t302;
t294 = t460 * t299 + t462 * t302;
t543 = qJ(6) * t520 - qJD(6) * t468 + t540;
t542 = -pkin(5) * t520 - t541;
t320 = t460 * t350 + t462 * t357;
t317 = qJ(6) * t523 + t320;
t361 = t462 * t405 + t460 * t474;
t539 = t317 - t361;
t538 = pkin(5) * t523 - t536;
t537 = t534 * pkin(5) - t533 * qJ(6) + qJD(6) * t421 - t570;
t326 = -pkin(4) * t447 + t339;
t308 = t460 * t326 + t334;
t423 = t467 * t435;
t374 = -qJ(5) * t549 + t423 + (-pkin(4) - t568) * t468;
t529 = t464 * t435 + t450;
t380 = -qJ(5) * t551 + t529;
t342 = t460 * t374 + t462 * t380;
t528 = pkin(4) * t551 + t465 * pkin(8);
t514 = qJD(6) * t447;
t314 = t339 * t462 - t564;
t512 = qJD(6) - t314;
t508 = qJ(6) * t492 + t294;
t506 = t575 * pkin(4) + pkin(8) * t519;
t505 = -pkin(4) * t467 - pkin(3);
t501 = t466 * t524;
t498 = t468 * t513;
t497 = t447 * t518;
t493 = t567 * t464;
t348 = t392 * t460 + t462 * t393;
t486 = t465 * t503;
t349 = t392 * t462 - t393 * t460;
t436 = t567 * t467;
t386 = -t436 * t460 - t462 * t493;
t387 = -t462 * t436 + t460 * t493;
t483 = -t387 * t348 + t349 * t386 - t361 * t369;
t432 = -t503 - t566;
t482 = -t432 - t503;
t307 = t326 * t462 - t564;
t341 = t374 * t462 - t380 * t460;
t420 = t460 * t464 - t554;
t479 = qJD(2) * t458 - t447 * t468;
t345 = pkin(4) * t393 + t365;
t412 = t468 * t556 + t553;
t378 = -t412 * t464 - t467 * t555;
t477 = -t412 * t467 + t464 * t555;
t411 = -t463 * t468 + t465 * t556;
t292 = -pkin(5) * t492 + t544;
t473 = qJD(3) * (-t482 - t566);
t295 = pkin(5) * t348 - qJ(6) * t349 - qJD(6) * t480 + t345;
t471 = qJD(2) ^ 2;
t453 = -pkin(4) * t462 - pkin(5);
t451 = pkin(4) * t460 + qJ(6);
t404 = -t460 * t551 + t462 * t549;
t403 = t421 * t465;
t377 = qJD(3) * t412 + t485;
t376 = -qJD(3) * t411 + t468 * t500;
t366 = pkin(5) * t420 - qJ(6) * t421 + t505;
t359 = t408 * t465 + t460 * t494 - t462 * t498;
t358 = t420 * t517 - t421 * t519;
t356 = pkin(5) * t403 - qJ(6) * t404 + t528;
t344 = t378 * t460 - t462 * t477;
t343 = -t462 * t378 - t460 * t477;
t337 = pkin(5) * t468 - t341;
t336 = -qJ(6) * t468 + t342;
t333 = qJD(4) * t378 + t376 * t467 + t464 * t501;
t332 = qJD(4) * t477 - t376 * t464 + t467 * t501;
t322 = pkin(4) * t427 + pkin(5) * t480 + qJ(6) * t369;
t315 = -pkin(5) * t358 + qJ(6) * t359 - qJD(6) * t404 + t506;
t312 = t332 * t460 + t333 * t462;
t310 = -t462 * t332 + t333 * t460;
t305 = -qJ(6) * t447 + t308;
t303 = pkin(5) * t447 + qJD(6) - t307;
t291 = t508 - t514;
t1 = [(-t332 * t447 + t377 * t425 + t393 * t411) * MDP(17) + (t333 * t447 + t377 * t427 + t392 * t411) * MDP(18) + (t294 * t344 - t307 * t310 + t308 * t312 + t343 * t544 + t345 * t411 + t362 * t377) * MDP(20) + (t310 * t447 + t348 * t411 + t369 * t377) * MDP(21) + (-t312 * t447 - t349 * t411 - t377 * t480) * MDP(23) + (t291 * t344 + t292 * t343 + t295 * t411 + t303 * t310 + t305 * t312 + t321 * t377) * MDP(24) + (-t377 * MDP(10) - t376 * MDP(11) + (MDP(17) * t378 + MDP(18) * t477 - MDP(21) * t343 + MDP(23) * t344) * t523) * qJD(3) + ((-MDP(10) * t465 - MDP(11) * t468) * t469 * t511 + (-t469 * MDP(4) + (-MDP(10) * t468 + MDP(11) * t465 - MDP(3)) * t466) * t471) * t461 + t569 * (t310 * t480 - t312 * t369 + t343 * t349 - t344 * t348); 0.2e1 * t491 * t573 - 0.2e1 * t511 * t572 + MDP(7) * t545 - MDP(8) * t548 + (-pkin(8) * t545 + t465 * t473) * MDP(10) + (pkin(8) * t548 + t468 * t473) * MDP(11) + (t392 * t549 + (-t496 + t498) * t427) * MDP(12) + ((-t425 * t467 - t427 * t464) * t519 + (-t560 - t393 * t467 + (t425 * t464 - t427 * t467) * qJD(4)) * t465) * MDP(13) + (t447 * t496 - t392 * t468 + (t427 * t465 + t467 * t479) * qJD(3)) * MDP(14) + (t447 * t495 + t393 * t468 + (-t425 * t465 - t464 * t479) * qJD(3)) * MDP(15) + (-t447 - t522) * MDP(16) * t520 + ((t435 * t518 - t579) * t447 + ((pkin(8) * t425 + t561) * qJD(3) + (t552 + (pkin(8) * t447 + t390) * t467) * qJD(4) + t488) * t468 + (-t425 * t503 + t389 * t516 + pkin(8) * t393 + t563 + ((-pkin(8) * t550 + t423) * qJD(2) + t352) * qJD(3)) * t465) * MDP(17) + (t581 * t447 + (t389 * t513 + (qJD(3) * t427 - t497) * pkin(8) + t476) * t468 + (-t427 * t503 - t389 * t518 + pkin(8) * t392 + t562 + (-pkin(8) * t557 - qJD(2) * t529 - t354) * qJD(3)) * t465) * MDP(18) + (-t294 * t403 + t307 * t359 + t308 * t358 - t341 * t349 - t342 * t348 - t369 * t540 + t404 * t544 - t480 * t541) * MDP(19) + (t294 * t342 - t544 * t341 + t345 * t528 + (-t486 + t506) * t362 + t540 * t308 + t541 * t307) * MDP(20) + (t292 * t468 + t295 * t403 + t315 * t369 - t321 * t358 + t348 * t356 + t542 * t447 + (-t369 * t503 + (-qJD(2) * t337 - t303) * qJD(3)) * t465) * MDP(21) + (-t291 * t403 + t292 * t404 - t303 * t359 + t305 * t358 - t336 * t348 + t337 * t349 - t369 * t543 + t480 * t542) * MDP(22) + (-t291 * t468 - t295 * t404 - t315 * t480 + t321 * t359 - t349 * t356 - t543 * t447 + (t480 * t503 + (qJD(2) * t336 + t305) * qJD(3)) * t465) * MDP(23) + (t291 * t336 + t292 * t337 + t295 * t356 + (t315 - t486) * t321 + t543 * t305 + t542 * t303) * MDP(24); (qJD(3) * t398 - t432 * t523 - t365) * MDP(10) + t482 * t522 * MDP(11) + (-t427 * t557 + t560) * MDP(12) + ((t392 + t559) * t467 + (-t393 + t558) * t464) * MDP(13) + (-t447 * t516 + (t447 * t547 + (-t427 + t521) * t465) * qJD(2)) * MDP(14) + (t497 + (-t447 * t550 + (t425 + t513) * t465) * qJD(2)) * MDP(15) + t447 * MDP(16) * t523 + (-pkin(3) * t393 - t562 + t487 * t447 - t398 * t425 + (pkin(9) * t557 + t561) * qJD(4) + (-t352 * t465 + (-pkin(9) * t520 - t389 * t468) * t464) * qJD(2)) * MDP(17) + (-pkin(3) * t392 + t563 - t535 * t447 - t398 * t427 + (-pkin(9) * t447 * t464 + t389 * t467) * qJD(4) + (-t389 * t547 + (-pkin(9) * t513 + t354) * t465) * qJD(2)) * MDP(18) + (-t294 * t420 + t307 * t533 + t308 * t534 + t320 * t369 + t421 * t544 - t480 * t536 + t483) * MDP(19) + (t294 * t387 + t544 * t386 + t345 * t505 + t570 * t362 + (t361 - t320) * t308 + t536 * t307) * MDP(20) + (t295 * t420 + t348 * t366 + t538 * t447 - t537 * t369 - t534 * t321 + (-qJD(3) * t386 + t303) * t523) * MDP(21) + (-t291 * t420 + t292 * t421 - t303 * t533 + t305 * t534 + t317 * t369 + t480 * t538 + t483) * MDP(22) + (-t295 * t421 - t349 * t366 + t539 * t447 + t537 * t480 + t533 * t321 + (qJD(3) * t387 - t305) * t523) * MDP(23) + (t291 * t387 + t292 * t386 + t295 * t366 + t303 * t538 - t305 * t539 - t321 * t537) * MDP(24) + (-t468 * t573 + t572) * t471; t427 * t425 * MDP(12) + (-t425 ^ 2 + t427 ^ 2) * MDP(13) + (t530 - t559) * MDP(14) + (-t490 - t558) * MDP(15) + (-t354 * t447 - t389 * t427 + t472) * MDP(17) + (-t352 * t447 + t389 * t425 - t476) * MDP(18) + (t308 * t480 - t565) * MDP(19) + (t307 * t313 - t308 * t314) * MDP(20) + (-t313 * t447 - t544 - t571) * MDP(21) + (t305 * t480 - t348 * t451 + t349 * t453 - t565) * MDP(22) + (t314 * t447 + t322 * t480 + t508 - 0.2e1 * t514) * MDP(23) + (t291 * t451 + t292 * t453 - t303 * t313 + t305 * t512 - t321 * t322) * MDP(24) + (-MDP(15) * t494 + ((-t464 * MDP(14) - t467 * MDP(15)) * qJD(4) + (MDP(16) + (pkin(5) - t453) * MDP(21) + t451 * MDP(23)) * qJD(3)) * t465) * qJD(2) + ((-t348 * t460 - t349 * t462) * MDP(19) + (t294 * t460 - t362 * t427 - t462 * t544) * MDP(20)) * pkin(4) + ((-t307 + t314) * MDP(19) - t322 * MDP(21) + (t303 - t512) * MDP(22) - t321 * MDP(23)) * t369; (t307 * t480 + t308 * t369 + t345) * MDP(20) + (-t447 * t480 + t348) * MDP(21) + (-t349 - t576) * MDP(23) + (-t303 * t480 + t305 * t369 + t295) * MDP(24) + t569 * (-t369 ^ 2 - t574); (t369 * t480 - t492) * MDP(21) + (t349 - t576) * MDP(22) + (-t447 ^ 2 - t574) * MDP(23) + (t305 * t447 + t292 + t571) * MDP(24);];
tauc  = t1;
