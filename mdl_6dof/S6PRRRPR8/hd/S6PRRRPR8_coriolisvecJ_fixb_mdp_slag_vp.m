% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:35
% EndTime: 2019-03-08 23:53:49
% DurationCPUTime: 8.09s
% Computational Cost: add. (4925->534), mult. (13435->735), div. (0->0), fcn. (10782->12), ass. (0->223)
t488 = sin(pkin(7));
t494 = sin(qJ(3));
t621 = t488 * t494;
t481 = pkin(9) * t621;
t490 = cos(pkin(7));
t498 = cos(qJ(3));
t499 = cos(qJ(2));
t612 = t498 * t499;
t495 = sin(qJ(2));
t616 = t494 * t495;
t518 = -t490 * t616 + t612;
t489 = sin(pkin(6));
t600 = qJD(1) * t489;
t619 = t490 * t498;
t605 = t518 * t600 - (pkin(2) * t619 - t481) * qJD(3);
t544 = pkin(3) * t494 - pkin(10) * t498;
t517 = t544 * qJD(3);
t575 = t495 * t600;
t650 = (t517 - t575) * t488;
t493 = sin(qJ(4));
t497 = cos(qJ(4));
t596 = qJD(2) * t490;
t558 = qJD(3) + t596;
t598 = qJD(2) * t488;
t574 = t494 * t598;
t434 = t493 * t558 + t497 * t574;
t427 = qJD(6) + t434;
t614 = t495 * t498;
t615 = t494 * t499;
t520 = t490 * t614 + t615;
t592 = qJD(3) * t498;
t571 = t488 * t592;
t594 = qJD(3) * t494;
t603 = t490 * pkin(2) * t594 + pkin(9) * t571 - t520 * t600;
t595 = qJD(2) * t498;
t573 = t488 * t595;
t474 = -qJD(4) + t573;
t458 = pkin(9) * t598 + t575;
t491 = cos(pkin(6));
t599 = qJD(1) * t491;
t576 = t488 * t599;
t470 = t494 * t576;
t638 = qJD(2) * pkin(2);
t468 = t499 * t600 + t638;
t624 = t468 * t490;
t390 = t498 * t458 + t494 * t624 + t470;
t379 = pkin(10) * t558 + t390;
t479 = t490 * t599;
t545 = -pkin(3) * t498 - pkin(10) * t494;
t397 = t479 + (qJD(2) * t545 - t468) * t488;
t350 = t379 * t493 - t497 * t397;
t586 = -qJD(5) - t350;
t584 = qJD(2) * qJD(3);
t566 = t488 * t584;
t547 = t494 * t566;
t536 = pkin(4) * t547;
t450 = t494 * t458;
t508 = t518 * qJD(2);
t549 = t491 * t571;
t365 = (t468 * t619 - t450) * qJD(3) + (t489 * t508 + t549) * qJD(1);
t410 = (t517 + t575) * t598;
t589 = qJD(4) * t497;
t590 = qJD(4) * t493;
t554 = t493 * t365 + t379 * t589 + t397 * t590 - t497 * t410;
t334 = -t536 + t554;
t351 = t497 * t379 + t493 * t397;
t348 = qJ(5) * t474 - t351;
t649 = -t348 * t474 + t334;
t620 = t488 * t498;
t441 = pkin(9) * t620 + (pkin(2) * t494 + pkin(10)) * t490;
t442 = (-pkin(2) + t545) * t488;
t648 = t441 * t589 + t442 * t590 - t493 * t605 - t650 * t497;
t647 = -t441 * t590 + t442 * t589 + t650 * t493 - t497 * t605;
t578 = t493 * t621;
t548 = qJD(4) * t578;
t413 = -t490 * t589 - t497 * t571 + t548;
t455 = t490 * t493 + t497 * t621;
t646 = qJ(5) * t413 - qJD(5) * t455 + t603;
t546 = t498 * t566;
t591 = qJD(4) * t434;
t406 = t493 * t546 + t591;
t550 = t493 * t573;
t645 = qJD(5) * t493 + t390 + (t550 - t590) * pkin(4);
t585 = pkin(5) * t434 - t586;
t526 = t497 * t558;
t432 = t493 * t574 - t526;
t639 = pkin(5) * t432;
t340 = -t348 - t639;
t405 = qJD(2) * t548 - qJD(4) * t526 - t497 * t546;
t641 = pkin(4) + pkin(11);
t644 = t641 * t405 + (t340 - t351 + t639) * t427;
t389 = t498 * (t576 + t624) - t450;
t485 = t488 ^ 2;
t643 = (-t494 * t498 * MDP(5) + (t494 ^ 2 - t498 ^ 2) * MDP(6)) * t485;
t610 = -t488 * (qJ(5) * t594 - qJD(5) * t498) - t647;
t642 = t434 ^ 2;
t501 = qJD(2) ^ 2;
t640 = pkin(5) + pkin(10);
t636 = qJ(5) * t432;
t635 = qJ(5) * t497;
t378 = -pkin(3) * t558 - t389;
t502 = -t434 * qJ(5) + t378;
t352 = t432 * pkin(4) + t502;
t633 = t352 * t434;
t492 = sin(qJ(6));
t496 = cos(qJ(6));
t587 = qJD(6) * t496;
t577 = t492 * t406 + t432 * t587 + t496 * t547;
t588 = qJD(6) * t492;
t356 = t474 * t588 + t577;
t632 = t356 * t496;
t623 = t474 * t492;
t402 = -t496 * t432 - t623;
t631 = t402 * t427;
t630 = t402 * t474;
t404 = t432 * t492 - t474 * t496;
t629 = t404 * t427;
t628 = t404 * t474;
t627 = t432 * t434;
t626 = t432 * t474;
t625 = t434 * t474;
t622 = t474 * t497;
t618 = t492 * t405;
t617 = t493 * t498;
t399 = t496 * t405;
t613 = t497 * t498;
t414 = qJD(4) * t455 + t493 * t571;
t611 = -pkin(5) * t414 - t610;
t572 = t488 * t594;
t609 = -pkin(4) * t572 + t648;
t608 = pkin(4) * t414 + t646;
t607 = qJ(5) * t589 - t573 * t635 + t645;
t445 = t544 * t598;
t606 = t497 * t389 + t493 * t445;
t604 = t497 * t441 + t493 * t442;
t602 = -t640 * t590 - (-pkin(5) * t617 + qJ(5) * t494) * t598 - t606;
t597 = qJD(2) * t489;
t593 = qJD(3) * t497;
t583 = pkin(10) * t474 * t493;
t582 = pkin(10) * t622;
t581 = pkin(10) * t594;
t580 = pkin(10) * t593;
t476 = t640 * t497;
t426 = -t468 * t488 + t479;
t570 = t426 * t598;
t568 = qJD(3) * t624;
t567 = qJD(1) * t597;
t565 = -qJ(5) * t493 - pkin(3);
t534 = t641 * t572;
t332 = -pkin(5) * t405 - qJD(2) * t534 + t554;
t528 = t490 * t495 * t567;
t366 = qJD(3) * t470 + t458 * t592 + t494 * t568 + t498 * t528 + t567 * t615;
t506 = qJ(5) * t405 - qJD(5) * t434 + t366;
t335 = t406 * t641 + t506;
t564 = t496 * t332 - t335 * t492;
t563 = -t493 * t389 + t445 * t497;
t562 = -t493 * t441 + t442 * t497;
t561 = t427 * t492;
t560 = t427 * t496;
t555 = -t497 * t365 + t379 * t590 - t397 * t589 - t493 * t410;
t551 = t488 * t495 * t597;
t384 = pkin(4) * t620 - t562;
t424 = t492 * t574 - t496 * t550;
t543 = t496 * t590 + t424;
t425 = (t492 * t617 + t494 * t496) * t598;
t542 = t492 * t590 - t425;
t454 = -t497 * t490 + t578;
t440 = t481 + (-pkin(2) * t498 - pkin(3)) * t490;
t507 = -qJ(5) * t455 + t440;
t369 = t454 * t641 + t507;
t541 = pkin(5) * t413 + qJD(6) * t369 + t534 - t648;
t364 = pkin(5) * t455 + pkin(11) * t620 + t384;
t540 = -qJD(6) * t364 - t414 * t641 - t646;
t457 = -t497 * t641 + t565;
t539 = qJD(6) * t457 + (pkin(5) * t613 - t494 * t641) * t598 - t563 - qJD(4) * t476;
t475 = t640 * t493;
t538 = -qJD(6) * t475 + t645 + t474 * (pkin(11) * t493 - t635);
t535 = -qJD(5) * t474 - t555;
t533 = t332 * t492 + t335 * t496;
t337 = t474 * t641 + t585;
t346 = t432 * t641 + t502;
t329 = t337 * t496 - t346 * t492;
t330 = t337 * t492 + t346 * t496;
t519 = t490 * t615 + t614;
t412 = t489 * t519 + t491 * t621;
t453 = -t488 * t489 * t499 + t490 * t491;
t381 = t412 * t493 - t453 * t497;
t521 = t490 * t612 - t616;
t411 = -t489 * t521 - t491 * t620;
t531 = t381 * t496 - t411 * t492;
t530 = t381 * t492 + t411 * t496;
t382 = t412 * t497 + t453 * t493;
t529 = qJ(5) * t547;
t383 = qJ(5) * t620 - t604;
t523 = t426 * t488 - t485 * t638;
t522 = -t454 * t492 + t496 * t620;
t415 = t454 * t496 + t492 * t620;
t513 = -t351 * t474 - t554;
t333 = -t529 - t535;
t331 = -pkin(5) * t406 - t333;
t511 = t331 + (t427 * t641 + t636) * t427;
t510 = -t496 * t406 + t492 * t547;
t505 = -t405 - t626;
t472 = -pkin(4) * t497 + t565;
t401 = t405 * t493;
t388 = pkin(4) * t434 + t636;
t385 = t405 * t455;
t380 = pkin(4) * t454 + t507;
t377 = t549 + (qJD(3) * t521 + t508) * t489;
t376 = t491 * t572 + (qJD(2) * t520 + qJD(3) * t519) * t489;
t370 = -pkin(5) * t454 - t383;
t368 = qJD(6) * t522 + t414 * t496 - t492 * t572;
t367 = qJD(6) * t415 + t414 * t492 + t496 * t572;
t363 = -pkin(4) * t574 - t563;
t359 = -qJ(5) * t574 - t606;
t357 = qJD(6) * t404 + t510;
t347 = pkin(4) * t474 - t586;
t342 = t493 * t551 - t412 * t590 + (qJD(4) * t453 + t377) * t497;
t341 = qJD(4) * t382 + t377 * t493 - t497 * t551;
t336 = pkin(4) * t406 + t506;
t327 = -qJD(6) * t330 + t564;
t326 = qJD(6) * t329 + t533;
t1 = [(-t376 * t558 + t453 * t547) * MDP(10) + (-t377 * t558 + t453 * t546) * MDP(11) + (t341 * t434 - t342 * t432 - t381 * t405 - t382 * t406) * MDP(19) + (-t333 * t382 + t334 * t381 + t336 * t411 + t341 * t347 - t342 * t348 + t352 * t376) * MDP(22) + ((-qJD(6) * t530 + t341 * t496 - t376 * t492) * t427 - t531 * t405 + t342 * t402 + t382 * t357) * MDP(28) + (-(qJD(6) * t531 + t341 * t492 + t376 * t496) * t427 + t530 * t405 + t342 * t404 + t382 * t356) * MDP(29) + (-MDP(4) * t499 + (-MDP(3) + (-MDP(10) * t498 + MDP(11) * t494) * t485) * t495) * t501 * t489 + (-MDP(17) + MDP(20)) * (-t341 * t474 - t376 * t432 + t381 * t547 - t406 * t411) + (-MDP(18) + MDP(21)) * (-t342 * t474 - t376 * t434 + t382 * t547 + t405 * t411); ((-qJD(2) * t603 - t366) * t490 + (t494 * t523 - t603) * qJD(3)) * MDP(10) + ((qJD(2) * t605 - t365) * t490 + (t498 * t523 + t605) * qJD(3)) * MDP(11) + (-t413 * t434 - t385) * MDP(12) + (t405 * t454 - t406 * t455 + t413 * t432 - t414 * t434) * MDP(13) + (t413 * t474 + (t405 * t498 + (qJD(2) * t455 + t434) * t594) * t488) * MDP(14) + (t414 * t474 + (t406 * t498 + (-qJD(2) * t454 - t432) * t594) * t488) * MDP(15) + (-t474 * t488 - t485 * t595) * MDP(16) * t594 + (t366 * t454 + t378 * t414 + t440 * t406 + t648 * t474 + t603 * t432 + (t554 * t498 + (qJD(2) * t562 - t350) * t594) * t488) * MDP(17) + (t366 * t455 - t378 * t413 - t440 * t405 + t647 * t474 + t603 * t434 + (-t555 * t498 + (-qJD(2) * t604 - t351) * t594) * t488) * MDP(18) + (t333 * t454 + t334 * t455 - t347 * t413 + t348 * t414 + t383 * t406 - t384 * t405 + t432 * t610 + t434 * t609) * MDP(19) + (-t336 * t454 - t352 * t414 - t380 * t406 - t609 * t474 - t608 * t432 + (-t334 * t498 + (qJD(2) * t384 + t347) * t594) * t488) * MDP(20) + (-t336 * t455 + t352 * t413 + t380 * t405 + t610 * t474 - t608 * t434 + (t333 * t498 + (-qJD(2) * t383 - t348) * t594) * t488) * MDP(21) + (t333 * t383 + t334 * t384 + t336 * t380 + t347 * t609 + t348 * t610 + t352 * t608) * MDP(22) + (-t356 * t522 + t367 * t404) * MDP(23) + (t356 * t415 + t357 * t522 - t367 * t402 + t368 * t404) * MDP(24) + (t356 * t455 + t367 * t427 - t404 * t413 + t405 * t522) * MDP(25) + (-t357 * t455 + t368 * t427 + t402 * t413 - t405 * t415) * MDP(26) + (-t413 * t427 - t385) * MDP(27) + (-(t364 * t496 - t369 * t492) * t405 + t327 * t455 - t329 * t413 + t370 * t357 - t331 * t415 - t340 * t368 + (t492 * t540 - t496 * t541) * t427 + t611 * t402) * MDP(28) + ((t364 * t492 + t369 * t496) * t405 - t326 * t455 + t330 * t413 + t370 * t356 - t331 * t522 + t340 * t367 + (t492 * t541 + t496 * t540) * t427 + t611 * t404) * MDP(29) - 0.2e1 * t643 * t584 + (MDP(7) * t571 - MDP(8) * t572) * (qJD(3) + 0.2e1 * t596); (-(-t457 * t492 + t475 * t496) * t405 + t327 * t493 + t476 * t357 + (t492 * t538 - t496 * t539) * t427 + t602 * t402 - t543 * t340 + (-t329 * t474 + t331 * t496 - t340 * t588) * t497) * MDP(28) + ((t457 * t496 + t475 * t492) * t405 - t326 * t493 + t476 * t356 + (t492 * t539 + t496 * t538) * t427 + t602 * t404 + t542 * t340 + (t330 * t474 - t331 * t492 - t340 * t587) * t497) * MDP(29) + (t390 * t558 - t494 * t570 - t366) * MDP(10) + (t389 * t558 + (qJD(3) * t458 + t528) * t494 + (-t570 - t568 + (-qJD(3) * t488 * t491 - t499 * t597) * qJD(1)) * t498) * MDP(11) + (-t434 * t622 - t401) * MDP(12) + ((-t405 + t626) * t497 + (-t406 + t625) * t493) * MDP(13) + (-t474 * t589 + (t474 * t613 + (qJD(3) * t493 - t434) * t494) * t598) * MDP(14) + (t474 * t590 + (-t474 * t617 + (t432 + t593) * t494) * t598) * MDP(15) + (-pkin(3) * t406 - t366 * t497 + t563 * t474 - t390 * t432 + (t378 * t493 + t582) * qJD(4) + (t350 * t494 + (-t378 * t498 - t581) * t493) * t598) * MDP(17) + (pkin(3) * t405 + t366 * t493 - t606 * t474 - t390 * t434 + (t378 * t497 - t583) * qJD(4) + (-t378 * t613 + (t351 - t580) * t494) * t598) * MDP(18) + (-t359 * t432 - t363 * t434 + (-t333 - t474 * t347 + (-t406 + t591) * pkin(10)) * t497 + ((qJD(4) * t432 - t405) * pkin(10) + t649) * t493) * MDP(19) + (t336 * t497 + t363 * t474 - t406 * t472 + t607 * t432 + (-t352 * t493 - t582) * qJD(4) + (-t347 * t494 + (t352 * t498 + t581) * t493) * t598) * MDP(20) + (-t336 * t493 - t359 * t474 + t405 * t472 + t607 * t434 + (-t352 * t497 + t583) * qJD(4) + (t352 * t613 + (t348 + t580) * t494) * t598) * MDP(21) + (t336 * t472 - t347 * t363 - t348 * t359 - t607 * t352 + (-t333 * t497 + t334 * t493 + (t347 * t497 + t348 * t493) * qJD(4)) * pkin(10)) * MDP(22) + (-t356 * t492 * t497 + (-t497 * t587 + t542) * t404) * MDP(23) + (t402 * t425 + t404 * t424 + (-t402 * t492 + t404 * t496) * t590 + (-t632 + t357 * t492 + (t402 * t496 + t404 * t492) * qJD(6)) * t497) * MDP(24) + (t356 * t493 + t542 * t427 + (-t427 * t587 + t618 - t628) * t497) * MDP(25) + (-t357 * t493 + t543 * t427 + (t427 * t588 + t399 + t630) * t497) * MDP(26) + (-t427 * t622 - t401) * MDP(27) + t474 * MDP(16) * t574 + ((-MDP(7) * t498 + MDP(8) * t494) * t488 * t490 + t643) * t501; MDP(12) * t627 + (-t432 ^ 2 + t642) * MDP(13) + t505 * MDP(14) + (-t625 - t406) * MDP(15) + MDP(16) * t547 + (-t378 * t434 + t513) * MDP(17) + (t350 * t474 + t378 * t432 + t555) * MDP(18) + (pkin(4) * t405 - qJ(5) * t406 + (-t348 - t351) * t434 + (t347 + t586) * t432) * MDP(19) + (t388 * t432 - t513 - 0.2e1 * t536 + t633) * MDP(20) + (-t352 * t432 + t388 * t434 + t474 * t586 + 0.2e1 * t529 + t535) * MDP(21) + (-pkin(4) * t334 - qJ(5) * t333 - t347 * t351 + t348 * t586 - t352 * t388) * MDP(22) + (-t404 * t561 + t632) * MDP(23) + ((-t357 - t629) * t496 + (-t356 + t631) * t492) * MDP(24) + (t404 * t432 - t427 * t561 - t399) * MDP(25) + (-t402 * t432 - t427 * t560 + t618) * MDP(26) + t427 * t432 * MDP(27) + (qJ(5) * t357 + t329 * t432 + t585 * t402 + t511 * t492 + t496 * t644) * MDP(28) + (qJ(5) * t356 - t330 * t432 + t585 * t404 - t492 * t644 + t511 * t496) * MDP(29); t505 * MDP(19) + (t547 - t627) * MDP(20) + (-t474 ^ 2 - t642) * MDP(21) + (t633 + t649) * MDP(22) + (-t399 + t630) * MDP(28) + (t618 + t628) * MDP(29) + (-MDP(28) * t561 - MDP(29) * t560) * t427; t404 * t402 * MDP(23) + (-t402 ^ 2 + t404 ^ 2) * MDP(24) + (t577 + t631) * MDP(25) + (-t510 + t629) * MDP(26) - t405 * MDP(27) + (t330 * t427 - t340 * t404 + t564) * MDP(28) + (t329 * t427 + t340 * t402 - t533) * MDP(29) + (MDP(25) * t623 - MDP(26) * t404 - MDP(28) * t330 - MDP(29) * t329) * qJD(6);];
tauc  = t1;
