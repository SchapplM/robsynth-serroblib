% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:25:04
% EndTime: 2019-03-09 05:25:15
% DurationCPUTime: 6.82s
% Computational Cost: add. (4654->448), mult. (10414->632), div. (0->0), fcn. (7152->8), ass. (0->201)
t490 = cos(qJ(6));
t488 = sin(qJ(4));
t491 = cos(qJ(4));
t542 = t491 * qJD(3);
t492 = cos(qJ(3));
t552 = qJD(1) * t492;
t458 = t488 * t552 - t542;
t533 = t491 * t552;
t551 = qJD(3) * t488;
t460 = t533 + t551;
t485 = sin(pkin(10));
t486 = cos(pkin(10));
t508 = -t458 * t486 - t460 * t485;
t569 = t490 * t508;
t402 = t458 * t485 - t460 * t486;
t487 = sin(qJ(6));
t585 = t402 * t487;
t349 = t569 + t585;
t489 = sin(qJ(3));
t553 = qJD(1) * t489;
t478 = qJD(4) + t553;
t473 = qJD(6) + t478;
t586 = t349 * t473;
t546 = qJD(4) * t492;
t529 = t488 * t546;
t531 = t489 * t542;
t497 = -t529 - t531;
t413 = qJD(1) * t497 + qJD(4) * t542;
t493 = -pkin(1) - pkin(7);
t475 = qJD(1) * t493 + qJD(2);
t464 = pkin(3) * t489 - pkin(8) * t492 + qJ(2);
t443 = t464 * qJD(1);
t463 = t489 * t475;
t446 = qJD(3) * pkin(8) + t463;
t396 = t443 * t488 + t446 * t491;
t516 = pkin(3) * t492 + pkin(8) * t489;
t452 = qJD(3) * t516 + qJD(2);
t432 = t452 * qJD(1);
t418 = t491 * t432;
t496 = -t396 * qJD(4) + t418;
t549 = qJD(3) * t492;
t325 = -qJ(5) * t413 - qJD(5) * t460 + (pkin(4) * qJD(1) - t475 * t488) * t549 + t496;
t527 = t488 * t553;
t414 = -qJD(3) * t527 + t460 * qJD(4);
t530 = t492 * t542;
t547 = qJD(4) * t491;
t548 = qJD(4) * t488;
t499 = t488 * t432 + t443 * t547 - t446 * t548 + t475 * t530;
t330 = -qJ(5) * t414 - qJD(5) * t458 + t499;
t308 = t486 * t325 - t330 * t485;
t364 = t413 * t486 - t414 * t485;
t541 = qJD(1) * qJD(3);
t526 = t492 * t541;
t306 = pkin(5) * t526 - pkin(9) * t364 + t308;
t309 = t485 * t325 + t486 * t330;
t363 = -t413 * t485 - t414 * t486;
t307 = pkin(9) * t363 + t309;
t395 = t491 * t443 - t446 * t488;
t369 = -qJ(5) * t460 + t395;
t357 = pkin(4) * t478 + t369;
t370 = -qJ(5) * t458 + t396;
t576 = t486 * t370;
t328 = t485 * t357 + t576;
t594 = pkin(9) * t508;
t318 = t328 + t594;
t544 = qJD(6) * t487;
t317 = t318 * t544;
t580 = t475 * t492;
t447 = -qJD(3) * pkin(3) - t580;
t408 = pkin(4) * t458 + qJD(5) + t447;
t354 = -pkin(5) * t508 + t408;
t604 = -t487 * t306 - t490 * t307 - t354 * t349 + t317;
t595 = -t490 * t402 + t487 * t508;
t603 = MDP(27) * t526 + (-t349 ^ 2 + t595 ^ 2) * MDP(24) - t349 * t595 * MDP(23);
t587 = t595 * t473;
t462 = t516 * qJD(1);
t445 = t491 * t462;
t589 = -qJ(5) - pkin(8);
t523 = qJD(4) * t589;
t574 = t488 * t492;
t538 = t475 * t574;
t572 = t489 * t491;
t601 = t538 - t445 - (pkin(4) * t492 + qJ(5) * t572) * qJD(1) - qJD(5) * t488 + t491 * t523;
t545 = qJD(5) * t491;
t568 = t491 * t492;
t557 = t488 * t462 + t475 * t568;
t600 = qJ(5) * t527 - t488 * t523 - t545 + t557;
t522 = t490 * t306 - t487 * t307;
t599 = -t354 * t595 + t522;
t598 = pkin(9) * t402;
t454 = t485 * t491 + t486 * t488;
t439 = t454 * qJD(4);
t440 = t454 * qJD(1);
t597 = t489 * t440 + t439;
t575 = t486 * t491;
t507 = t485 * t488 - t575;
t591 = qJD(4) * t507;
t596 = -t485 * t527 + t553 * t575 - t591;
t540 = 0.2e1 * qJD(1);
t484 = t492 ^ 2;
t593 = MDP(8) * (t489 ^ 2 - t484);
t560 = t485 * t600 + t486 * t601;
t559 = t485 * t601 - t486 * t600;
t571 = t489 * t493;
t556 = t488 * t464 + t491 * t571;
t521 = -t490 * t363 + t364 * t487;
t313 = qJD(6) * t595 + t521;
t590 = pkin(4) * t485;
t495 = qJD(1) ^ 2;
t588 = qJ(2) * t495;
t584 = t413 * t488;
t583 = t447 * t488;
t582 = t458 * t478;
t581 = t460 * t478;
t579 = t478 * t488;
t578 = t478 * t489;
t577 = t478 * t491;
t365 = t485 * t370;
t573 = t488 * t493;
t327 = t486 * t357 - t365;
t316 = pkin(5) * t478 + t327 + t598;
t570 = t490 * t316;
t567 = t492 * t493;
t494 = qJD(3) ^ 2;
t565 = t493 * t494;
t437 = t491 * t452;
t524 = pkin(4) - t573;
t341 = qJ(5) * t531 + t437 - t556 * qJD(4) + (qJ(5) * t548 + qJD(3) * t524 - t545) * t492;
t528 = t491 * t546;
t535 = t488 * t452 + t464 * t547 + t493 * t530;
t343 = -qJ(5) * t528 + (-qJD(5) * t492 + (qJ(5) * qJD(3) - qJD(4) * t493) * t489) * t488 + t535;
t315 = t485 * t341 + t486 * t343;
t509 = -t454 * t487 - t490 * t507;
t564 = qJD(6) * t509 - t597 * t487 + t596 * t490;
t401 = t454 * t490 - t487 * t507;
t563 = qJD(6) * t401 + t596 * t487 + t597 * t490;
t333 = t486 * t369 - t365;
t428 = t454 * t492;
t562 = t507 * qJD(1) - qJD(3) * t428 + t489 * t591;
t430 = t507 * t492;
t561 = -qJD(3) * t430 - t439 * t489 - t440;
t451 = t491 * t464;
t399 = -qJ(5) * t568 + t489 * t524 + t451;
t409 = -qJ(5) * t574 + t556;
t352 = t485 * t399 + t486 * t409;
t558 = -t463 + t597 * pkin(5) + (t527 + t548) * pkin(4);
t468 = t589 * t488;
t469 = t589 * t491;
t411 = t485 * t468 - t486 * t469;
t550 = qJD(3) * t489;
t537 = t488 * t571;
t536 = qJD(6) * t569 + t487 * t363 + t490 * t364;
t534 = -pkin(4) * t491 - pkin(3);
t532 = t488 * t550;
t525 = MDP(18) * t552;
t390 = pkin(4) * t414 + t475 * t550;
t314 = t486 * t341 - t343 * t485;
t332 = -t369 * t485 - t576;
t351 = t486 * t399 - t409 * t485;
t410 = t486 * t468 + t469 * t485;
t520 = pkin(4) * t574 - t567;
t519 = t458 + t542;
t518 = -t460 + t551;
t517 = qJD(4) * t489 + qJD(1);
t471 = t489 * t526;
t389 = -pkin(9) * t507 + t411;
t515 = pkin(5) * t552 + t596 * pkin(9) + qJD(6) * t389 - t560;
t388 = -pkin(9) * t454 + t410;
t514 = t597 * pkin(9) - qJD(6) * t388 - t559;
t429 = t507 * t489;
t513 = -qJD(6) * t429 - t562;
t427 = t454 * t489;
t512 = qJD(6) * t427 - t561;
t304 = t487 * t316 + t490 * t318;
t335 = pkin(5) * t489 + pkin(9) * t430 + t351;
t336 = -pkin(9) * t428 + t352;
t511 = t335 * t487 + t336 * t490;
t510 = -t490 * t428 + t430 * t487;
t378 = -t428 * t487 - t430 * t490;
t506 = qJD(1) * t484 - t578;
t480 = pkin(4) * t486 + pkin(5);
t505 = t480 * t487 + t490 * t590;
t504 = t480 * t490 - t487 * t590;
t503 = -pkin(8) * t549 + t447 * t489;
t500 = t493 * t550 + (t528 - t532) * pkin(4);
t312 = t402 * t544 + t536;
t421 = pkin(5) * t507 + t534;
t407 = pkin(5) * t428 + t520;
t385 = t454 * t546 - t485 * t532 + t486 * t531;
t383 = t454 * t550 + t492 * t591;
t371 = pkin(4) * t460 - pkin(5) * t402;
t353 = -pkin(5) * t383 + t500;
t334 = -pkin(5) * t363 + t390;
t322 = qJD(6) * t378 - t490 * t383 - t385 * t487;
t321 = qJD(6) * t510 + t383 * t487 - t385 * t490;
t320 = t333 + t598;
t319 = t332 - t594;
t311 = pkin(9) * t383 + t315;
t310 = pkin(5) * t549 + pkin(9) * t385 + t314;
t303 = -t318 * t487 + t570;
t1 = [0.2e1 * t541 * t593 + (-t489 * t565 + (qJ(2) * t549 + qJD(2) * t489) * t540) * MDP(12) + (-t492 * t565 + (-qJ(2) * t550 + qJD(2) * t492) * t540) * MDP(13) + (t413 * t568 + t460 * t497) * MDP(14) + ((t458 * t491 + t460 * t488) * t550 + (-t584 - t414 * t491 + (t458 * t488 - t460 * t491) * qJD(4)) * t492) * MDP(15) + (-t478 * t529 + t413 * t489 + (t460 * t492 + t491 * t506) * qJD(3)) * MDP(16) + (-t478 * t528 - t414 * t489 + (-t458 * t492 - t488 * t506) * qJD(3)) * MDP(17) + (t478 * t549 + t471) * MDP(18) + (-t414 * t567 + t418 * t489 + t437 * t478 + (-t396 * t489 + t447 * t568 - t478 * t556) * qJD(4) + ((t458 * t493 - t583) * t489 + (-t478 * t573 + (t451 - t537) * qJD(1) + t395) * t492) * qJD(3)) * MDP(19) + (-(-qJD(4) * t537 + t535) * t478 - t499 * t489 + (-t493 * t413 - t447 * t548) * t492 + ((-qJD(1) * t556 - t396) * t492 + (t493 * t460 + (-t447 + t580) * t491) * t489) * qJD(3)) * MDP(20) + (t308 * t430 - t309 * t428 + t314 * t402 + t315 * t508 + t327 * t385 + t328 * t383 - t351 * t364 + t352 * t363) * MDP(21) + (t308 * t351 + t309 * t352 + t327 * t314 + t328 * t315 + t390 * t520 + t408 * t500) * MDP(22) + (t312 * t378 + t321 * t595) * MDP(23) + (t312 * t510 - t313 * t378 + t321 * t349 - t322 * t595) * MDP(24) + (t312 * t489 + t321 * t473 + (qJD(1) * t378 + t595) * t549) * MDP(25) + (-t313 * t489 - t322 * t473 + (qJD(1) * t510 + t349) * t549) * MDP(26) + (t473 * t549 + t471) * MDP(27) + ((t310 * t490 - t311 * t487) * t473 + t522 * t489 - t353 * t349 + t407 * t313 - t334 * t510 + t354 * t322 + (-t304 * t489 - t473 * t511) * qJD(6) + ((t335 * t490 - t336 * t487) * qJD(1) + t303) * t549) * MDP(28) + (t407 * t312 + t317 * t489 + t354 * t321 + t334 * t378 + t353 * t595 + (-(-qJD(6) * t336 + t310) * t473 - t306 * t489) * t487 + (-(qJD(6) * t335 + t311) * t473 - (qJD(6) * t316 + t307) * t489) * t490 + (-qJD(1) * t511 - t304) * t549) * MDP(29) - 0.2e1 * MDP(7) * t471 + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t540 + (-MDP(10) * t492 - MDP(9) * t489) * t494; -t495 * MDP(5) - MDP(6) * t588 + (-t414 * t492 - t517 * t577 + (t458 * t489 + (-t478 - t553) * t574) * qJD(3)) * MDP(19) + (-t413 * t492 + t517 * t579 + (-t478 * t568 + (t460 - t533) * t489) * qJD(3)) * MDP(20) + (-t363 * t429 + t364 * t427 + t402 * t562 + t508 * t561) * MDP(21) + (-t308 * t427 - t309 * t429 + t327 * t562 + t328 * t561 - t390 * t492 + t408 * t550) * MDP(22) + (-t492 * t313 + (t487 * t512 - t490 * t513) * t473 + ((-t427 * t490 + t429 * t487) * t552 - t489 * t349) * qJD(3)) * MDP(28) + (-t492 * t312 + (t487 * t513 + t490 * t512) * t473 + (-(-t427 * t487 - t429 * t490) * t552 + t489 * t595) * qJD(3)) * MDP(29) + (t489 * MDP(12) + t492 * MDP(13)) * (-t494 - t495); t489 * MDP(13) * t588 + (t460 * t577 + t584) * MDP(14) + ((t413 - t582) * t491 + (-t414 - t581) * t488) * MDP(15) + (t478 * t547 + (t478 * t572 + t492 * t518) * qJD(1)) * MDP(16) + (-t478 * t548 + (-t488 * t578 + t492 * t519) * qJD(1)) * MDP(17) - t478 * t525 + (-pkin(3) * t414 - t445 * t478 + (t478 * t574 - t489 * t519) * t475 + (-pkin(8) * t577 + t583) * qJD(4) + (-t395 * t492 + t488 * t503) * qJD(1)) * MDP(19) + (-pkin(3) * t413 + t557 * t478 + t518 * t463 + (pkin(8) * t579 + t447 * t491) * qJD(4) + (t396 * t492 + t491 * t503) * qJD(1)) * MDP(20) + (-t308 * t454 - t309 * t507 - t596 * t327 - t597 * t328 + t363 * t411 - t364 * t410 + t560 * t402 + t559 * t508) * MDP(21) + (t309 * t411 + t308 * t410 + t390 * t534 + (pkin(4) * t579 - t463) * t408 + t559 * t328 + t560 * t327) * MDP(22) + (t312 * t401 + t564 * t595) * MDP(23) + (t312 * t509 - t313 * t401 + t349 * t564 - t563 * t595) * MDP(24) + (t564 * t473 + (qJD(3) * t401 - t595) * t552) * MDP(25) + (-t563 * t473 + (qJD(3) * t509 - t349) * t552) * MDP(26) - t473 * MDP(27) * t552 + (t421 * t313 - t334 * t509 + (t487 * t514 - t490 * t515) * t473 + t563 * t354 - t558 * t349 + ((t388 * t490 - t389 * t487) * qJD(3) - t303) * t552) * MDP(28) + (t421 * t312 + t334 * t401 + (t487 * t515 + t490 * t514) * t473 + t564 * t354 + t558 * t595 + (-(t388 * t487 + t389 * t490) * qJD(3) + t304) * t552) * MDP(29) + (-t593 + (-qJ(2) * MDP(12) + t489 * MDP(7)) * t492) * t495; t460 * t458 * MDP(14) + (-t458 ^ 2 + t460 ^ 2) * MDP(15) + (t413 + t582) * MDP(16) + (-t414 + t581) * MDP(17) + qJD(3) * t525 + (-qJD(3) * t538 + t396 * t478 - t447 * t460 + t496) * MDP(19) + (t395 * t478 + t447 * t458 - t499) * MDP(20) + ((t363 * t485 - t364 * t486) * pkin(4) + (t327 - t333) * t508 + (-t328 - t332) * t402) * MDP(21) + (-t327 * t332 - t328 * t333 + (t308 * t486 + t309 * t485 - t408 * t460) * pkin(4)) * MDP(22) + (t312 - t586) * MDP(25) + (-t313 + t587) * MDP(26) + (t504 * t526 - (t319 * t490 - t320 * t487) * t473 + t371 * t349 + (-t473 * t505 - t304) * qJD(6) + t599) * MDP(28) + (-t505 * t526 + (t319 * t487 + t320 * t490) * t473 - t371 * t595 + (-t473 * t504 - t570) * qJD(6) + t604) * MDP(29) + t603; (-t402 ^ 2 - t508 ^ 2) * MDP(21) + (-t327 * t402 - t328 * t508 + t390) * MDP(22) + (t313 + t587) * MDP(28) + (t312 + t586) * MDP(29); (t536 - t586) * MDP(25) + (-t521 + t587) * MDP(26) + (t304 * t473 + t599) * MDP(28) + (t303 * t473 + t604) * MDP(29) + (MDP(25) * t585 - MDP(26) * t595 - MDP(28) * t304 - MDP(29) * t570) * qJD(6) + t603;];
tauc  = t1;
