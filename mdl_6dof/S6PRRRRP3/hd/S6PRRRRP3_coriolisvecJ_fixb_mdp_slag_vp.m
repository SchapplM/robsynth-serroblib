% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:18
% EndTime: 2019-03-09 00:11:29
% DurationCPUTime: 6.10s
% Computational Cost: add. (4492->456), mult. (11091->623), div. (0->0), fcn. (8161->10), ass. (0->201)
t485 = sin(qJ(3));
t489 = cos(qJ(3));
t505 = pkin(3) * t485 - pkin(9) * t489;
t444 = t505 * qJD(3);
t449 = -pkin(3) * t489 - pkin(9) * t485 - pkin(2);
t484 = sin(qJ(4));
t486 = sin(qJ(2));
t488 = cos(qJ(4));
t541 = qJD(4) * t488;
t481 = sin(pkin(6));
t552 = qJD(1) * t481;
t490 = cos(qJ(2));
t570 = t489 * t490;
t613 = -(t484 * t486 + t488 * t570) * t552 + t484 * t444 + t449 * t541;
t546 = qJD(3) * t485;
t590 = pkin(8) * t484;
t612 = t488 * t444 + t546 * t590 - (-t484 * t570 + t486 * t488) * t552;
t571 = t488 * t489;
t471 = pkin(8) * t571;
t500 = pkin(4) * t485 - pkin(10) * t571;
t611 = t500 * qJD(3) + (-t471 + (pkin(10) * t485 - t449) * t484) * qJD(4) + t612;
t543 = qJD(4) * t484;
t545 = qJD(3) * t488;
t544 = qJD(3) * t489;
t521 = t484 * t544;
t522 = t485 * t541;
t604 = t521 + t522;
t610 = -t604 * pkin(10) + (-t485 * t545 - t489 * t543) * pkin(8) + t613;
t530 = t486 * t552;
t447 = qJD(2) * pkin(8) + t530;
t482 = cos(pkin(6));
t551 = qJD(1) * t489;
t407 = -t485 * t447 + t482 * t551;
t441 = t505 * qJD(2);
t510 = -t407 * t484 + t488 * t441;
t591 = pkin(9) + pkin(10);
t531 = qJD(4) * t591;
t609 = qJD(2) * t500 + t488 * t531 + t510;
t548 = qJD(2) * t489;
t525 = t484 * t548;
t558 = t488 * t407 + t484 * t441;
t606 = -pkin(10) * t525 + t484 * t531 + t558;
t549 = qJD(2) * t485;
t434 = -t484 * t549 + t545;
t487 = cos(qJ(5));
t547 = qJD(3) * t484;
t435 = t488 * t549 + t547;
t483 = sin(qJ(5));
t582 = t435 * t483;
t378 = -t487 * t434 + t582;
t376 = t378 ^ 2;
t537 = qJD(2) * qJD(3);
t519 = t485 * t537;
t502 = t434 * t483 + t487 * t435;
t592 = t502 ^ 2;
t608 = MDP(23) * t519 + (-t376 + t592) * MDP(20);
t607 = qJ(6) * t378;
t577 = t483 * t484;
t436 = -t487 * t488 + t577;
t535 = qJD(4) + qJD(5);
t539 = qJD(5) * t487;
t560 = -t436 * t548 - t487 * t541 - t488 * t539 + t535 * t577;
t437 = t483 * t488 + t484 * t487;
t387 = t535 * t437;
t559 = -t437 * t548 + t387;
t518 = t489 * t537;
t605 = qJD(3) * qJD(4) + t518;
t603 = MDP(5) * t485;
t479 = t485 ^ 2;
t602 = MDP(6) * (-t489 ^ 2 + t479);
t601 = qJ(6) * t502;
t400 = -qJD(3) * pkin(3) - t407;
t368 = -pkin(4) * t434 + t400;
t339 = pkin(5) * t378 + qJD(6) + t368;
t600 = t339 * t502;
t599 = t368 * t502;
t469 = -qJD(4) + t548;
t459 = -qJD(5) + t469;
t598 = t459 * t502;
t597 = t609 * t487;
t596 = t610 * t483 - t611 * t487;
t433 = t488 * t449;
t573 = t485 * t488;
t383 = -pkin(10) * t573 + t433 + (-pkin(4) - t590) * t489;
t554 = t484 * t449 + t471;
t575 = t484 * t485;
t392 = -pkin(10) * t575 + t554;
t540 = qJD(5) * t483;
t595 = t383 * t539 - t392 * t540 + t611 * t483 + t610 * t487;
t578 = t482 * t485;
t464 = qJD(1) * t578;
t408 = t489 * t447 + t464;
t506 = -t408 + (-t525 + t543) * pkin(4);
t561 = t483 * t383 + t487 * t392;
t453 = t591 * t484;
t454 = t591 * t488;
t555 = -t483 * t453 + t487 * t454;
t594 = -t453 * t539 - t454 * t540 - t609 * t483 - t606 * t487;
t520 = t488 * t544;
t542 = qJD(4) * t485;
t523 = t484 * t542;
t593 = t520 - t523;
t517 = qJD(2) * t542;
t403 = -t484 * t517 + t605 * t488;
t532 = t605 * t484 + t488 * t517;
t511 = t403 * t483 + t487 * t532;
t333 = qJD(5) * t502 + t511;
t589 = qJD(2) * pkin(2);
t550 = qJD(2) * t481;
t526 = t490 * t550;
t507 = t485 * t526;
t372 = qJD(1) * t507 + qJD(3) * t464 + t447 * t544;
t588 = t372 * t484;
t587 = t372 * t488;
t586 = t400 * t484;
t585 = t403 * t484;
t584 = t434 * t469;
t583 = t435 * t469;
t581 = t469 * t488;
t580 = t481 * t486;
t579 = t481 * t490;
t401 = qJD(3) * pkin(9) + t408;
t529 = t490 * t552;
t410 = qJD(2) * t449 - t529;
t576 = t484 * t410;
t362 = t401 * t488 + t576;
t349 = pkin(10) * t434 + t362;
t343 = t483 * t349;
t574 = t484 * t489;
t491 = qJD(3) ^ 2;
t572 = t485 * t491;
t345 = t487 * t349;
t569 = t489 * t491;
t356 = t387 * t485 + t483 * t521 - t487 * t520;
t414 = t436 * t485;
t568 = pkin(5) * t546 + qJ(6) * t356 - qJD(5) * t561 + qJD(6) * t414 - t596;
t357 = -t540 * t575 + (t535 * t573 + t521) * t487 + t593 * t483;
t413 = t437 * t485;
t567 = -qJ(6) * t357 - qJD(6) * t413 + t595;
t361 = -t401 * t484 + t488 * t410;
t348 = -pkin(10) * t435 + t361;
t338 = -pkin(4) * t469 + t348;
t323 = t487 * t338 - t343;
t316 = t323 - t601;
t315 = -pkin(5) * t459 + t316;
t566 = t315 - t316;
t565 = -t559 * qJ(6) - qJD(6) * t436 + t594;
t564 = -pkin(5) * t549 + t560 * qJ(6) - t555 * qJD(5) - qJD(6) * t437 + t606 * t483 - t597;
t563 = t487 * t348 - t343;
t445 = pkin(4) * t575 + t485 * pkin(8);
t538 = t502 * MDP(19);
t533 = t487 * t403 + t434 * t539 - t483 * t532;
t409 = t604 * pkin(4) + pkin(8) * t544;
t476 = -pkin(4) * t488 - pkin(3);
t527 = t486 * t550;
t524 = t469 * t543;
t371 = -t447 * t546 + (qJD(3) * t482 + t526) * t551;
t406 = (t444 + t530) * qJD(2);
t513 = t484 * t371 - t488 * t406;
t493 = -qJD(4) * t362 - t513;
t322 = pkin(4) * t519 - pkin(10) * t403 + t493;
t498 = t488 * t371 - t401 * t543 + t484 * t406 + t410 * t541;
t327 = -pkin(10) * t532 + t498;
t516 = t487 * t322 - t483 * t327;
t515 = -t348 * t483 - t345;
t512 = t487 * t383 - t392 * t483;
t509 = -t487 * t453 - t454 * t483;
t508 = -t483 * t322 - t487 * t327 - t338 * t539 + t349 * t540;
t448 = -t529 - t589;
t504 = -t448 - t529;
t324 = t338 * t483 + t345;
t420 = t489 * t580 + t578;
t390 = -t420 * t484 - t488 * t579;
t499 = -t420 * t488 + t484 * t579;
t352 = t390 * t487 + t483 * t499;
t353 = t390 * t483 - t487 * t499;
t501 = qJD(2) * t479 - t469 * t489;
t419 = -t482 * t489 + t485 * t580;
t332 = t435 * t540 - t533;
t355 = pkin(4) * t532 + t372;
t495 = qJD(3) * (-t504 - t589);
t325 = t333 * pkin(5) + t355;
t494 = -t324 * qJD(5) + t516;
t492 = qJD(2) ^ 2;
t475 = pkin(4) * t487 + pkin(5);
t389 = qJD(3) * t420 + t507;
t388 = -qJD(3) * t419 + t489 * t526;
t370 = -qJ(6) * t436 + t555;
t369 = -qJ(6) * t437 + t509;
t342 = qJD(4) * t390 + t388 * t488 + t484 * t527;
t341 = qJD(4) * t499 - t388 * t484 + t488 * t527;
t335 = -qJ(6) * t413 + t561;
t334 = -pkin(5) * t489 + qJ(6) * t414 + t512;
t319 = t563 - t601;
t318 = t515 + t607;
t317 = t324 - t607;
t314 = -qJD(5) * t353 + t341 * t487 - t342 * t483;
t313 = qJD(5) * t352 + t341 * t483 + t342 * t487;
t310 = -qJ(6) * t333 - qJD(6) * t378 - t508;
t309 = pkin(5) * t519 + qJ(6) * t332 - qJD(6) * t502 + t494;
t1 = [(-t341 * t469 - t389 * t434 + t419 * t532) * MDP(17) + (t342 * t469 + t389 * t435 + t403 * t419) * MDP(18) + (-t314 * t459 + t333 * t419 + t378 * t389) * MDP(24) + (t313 * t459 - t332 * t419 + t389 * t502) * MDP(25) + (-t313 * t378 - t314 * t502 + t332 * t352 - t333 * t353) * MDP(26) + (t309 * t352 + t310 * t353 + t313 * t317 + t314 * t315 + t325 * t419 + t339 * t389) * MDP(27) + (-t389 * MDP(10) - t388 * MDP(11) + (MDP(17) * t390 + MDP(18) * t499 + MDP(24) * t352 - MDP(25) * t353) * t549) * qJD(3) + ((-MDP(10) * t485 - MDP(11) * t489) * t490 * t537 + (-t490 * MDP(4) + (-MDP(10) * t489 + MDP(11) * t485 - MDP(3)) * t486) * t492) * t481; 0.2e1 * t518 * t603 - 0.2e1 * t537 * t602 + MDP(7) * t569 - MDP(8) * t572 + (-pkin(8) * t569 + t485 * t495) * MDP(10) + (pkin(8) * t572 + t489 * t495) * MDP(11) + (t403 * t573 + t593 * t435) * MDP(12) + ((t434 * t488 - t435 * t484) * t544 + (-t488 * t532 - t585 + (-t484 * t434 - t435 * t488) * qJD(4)) * t485) * MDP(13) + (t469 * t523 - t403 * t489 + (t435 * t485 + t488 * t501) * qJD(3)) * MDP(14) + (t469 * t522 + t532 * t489 + (t434 * t485 - t484 * t501) * qJD(3)) * MDP(15) + ((t449 * t543 - t612) * t469 + ((-pkin(8) * t434 + t586) * qJD(3) + (t576 + (pkin(8) * t469 + t401) * t488) * qJD(4) + t513) * t489 + (pkin(8) * t532 + t588 + t400 * t541 + t434 * t529 + ((-pkin(8) * t574 + t433) * qJD(2) + t361) * qJD(3)) * t485) * MDP(17) + (t613 * t469 + (t400 * t545 + (qJD(3) * t435 - t524) * pkin(8) + t498) * t489 + (-t435 * t529 - t400 * t543 + pkin(8) * t403 + t587 + (-pkin(8) * t581 - qJD(2) * t554 - t362) * qJD(3)) * t485) * MDP(18) + (t332 * t414 - t356 * t502) * MDP(19) + (t332 * t413 + t333 * t414 + t356 * t378 - t357 * t502) * MDP(20) + (t332 * t489 + t356 * t459) * MDP(21) + (t333 * t489 + t357 * t459) * MDP(22) + (-t516 * t489 + t409 * t378 + t445 * t333 + t355 * t413 + t368 * t357 + t596 * t459 + (t324 * t489 + t459 * t561) * qJD(5) + (-t378 * t529 + (qJD(2) * t512 + t323) * qJD(3)) * t485) * MDP(24) + (-t508 * t489 + t409 * t502 - t445 * t332 - t355 * t414 - t368 * t356 + t595 * t459 + (-t502 * t529 + (-qJD(2) * t561 - t324) * qJD(3)) * t485) * MDP(25) + (t309 * t414 - t310 * t413 + t315 * t356 - t317 * t357 + t332 * t334 - t333 * t335 - t378 * t567 - t502 * t568) * MDP(26) + (t310 * t335 + t309 * t334 + t325 * (pkin(5) * t413 + t445) + (pkin(5) * t357 - t485 * t529 + t409) * t339 + t567 * t317 + t568 * t315) * MDP(27) + ((-t469 - t548) * MDP(16) + (-qJD(2) * t414 + t502) * MDP(21) + (-qJD(2) * t413 - t378) * MDP(22) + (-t459 - t548) * MDP(23)) * t546; (qJD(3) * t408 - t372) * MDP(10) + t504 * t548 * MDP(11) + (-t435 * t581 + t585) * MDP(12) + ((t403 - t584) * t488 + (-t532 + t583) * t484) * MDP(13) + (-t469 * t541 + (t469 * t571 + (-t435 + t547) * t485) * qJD(2)) * MDP(14) + (t524 + (-t469 * t574 + (-t434 + t545) * t485) * qJD(2)) * MDP(15) + (-pkin(3) * t532 - t587 + t510 * t469 + t408 * t434 + (pkin(9) * t581 + t586) * qJD(4) + (-t361 * t485 + (-pkin(9) * t546 - t400 * t489) * t484) * qJD(2)) * MDP(17) + (-pkin(3) * t403 + t588 - t558 * t469 - t408 * t435 + (-pkin(9) * t469 * t484 + t400 * t488) * qJD(4) + (-t400 * t571 + (-pkin(9) * t545 + t362) * t485) * qJD(2)) * MDP(18) + (-t332 * t437 - t502 * t560) * MDP(19) + (t332 * t436 - t333 * t437 + t378 * t560 - t502 * t559) * MDP(20) + (t476 * t333 + t355 * t436 + t559 * t368 + t506 * t378) * MDP(24) + (-t476 * t332 + t355 * t437 - t560 * t368 + t502 * t506) * MDP(25) + (-t309 * t437 - t310 * t436 + t315 * t560 - t317 * t559 + t332 * t369 - t333 * t370 - t378 * t565 - t502 * t564) * MDP(26) + (t310 * t370 + t309 * t369 + t325 * (pkin(5) * t436 + t476) + (pkin(5) * t559 + t506) * t339 + t565 * t317 + t564 * t315) * MDP(27) + (t560 * MDP(21) + t559 * MDP(22) + (t454 * t539 + (-qJD(5) * t453 - t606) * t483 + t597) * MDP(24) + t594 * MDP(25)) * t459 + (-t448 * MDP(10) + t469 * MDP(16) + (qJD(3) * t437 - t502) * MDP(21) + (-qJD(3) * t436 + t378) * MDP(22) + t459 * MDP(23) + (qJD(3) * t509 - t323) * MDP(24) + (-qJD(3) * t555 + t324) * MDP(25)) * t549 + (-t489 * t603 + t602) * t492; -t435 * t434 * MDP(12) + (-t434 ^ 2 + t435 ^ 2) * MDP(13) + (t403 + t584) * MDP(14) + (-t532 - t583) * MDP(15) + MDP(16) * t519 + (-t362 * t469 - t400 * t435 + t493) * MDP(17) + (-t361 * t469 - t400 * t434 - t498) * MDP(18) + t378 * t538 + (-t459 * t378 - t332) * MDP(21) + (-t333 - t598) * MDP(22) + (t515 * t459 - t599 + (-t378 * t435 + t459 * t540 + t487 * t519) * pkin(4) + t494) * MDP(24) + (-t563 * t459 + t368 * t378 + (-t435 * t502 + t459 * t539 - t483 * t519) * pkin(4) + t508) * MDP(25) + (-t315 * t378 + t317 * t502 + t318 * t502 + t319 * t378 + t332 * t475 + (-t333 * t483 + (-t378 * t487 + t483 * t502) * qJD(5)) * pkin(4)) * MDP(26) + (-pkin(5) * t600 + t309 * t475 - t315 * t318 - t317 * t319 + (t310 * t483 - t339 * t435 + (-t315 * t483 + t317 * t487) * qJD(5)) * pkin(4)) * MDP(27) + t608; t533 * MDP(21) + (-t511 - t598) * MDP(22) + (-t324 * t459 + t516 - t599) * MDP(24) + (-t323 * t459 + t508) * MDP(25) + t566 * MDP(27) * t317 + (t332 * MDP(26) + (t309 - t600) * MDP(27)) * pkin(5) + (-t459 * MDP(21) + t368 * MDP(25) - MDP(26) * t566 + t538) * t378 + (-MDP(21) * t582 - t502 * MDP(22) - MDP(24) * t324) * qJD(5) + t608; (-t376 - t592) * MDP(26) + (t315 * t502 + t317 * t378 + t325) * MDP(27);];
tauc  = t1;
