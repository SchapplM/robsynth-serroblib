% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:57
% EndTime: 2019-03-09 08:59:14
% DurationCPUTime: 11.88s
% Computational Cost: add. (9591->515), mult. (29159->709), div. (0->0), fcn. (24316->12), ass. (0->234)
t541 = sin(pkin(6));
t549 = cos(qJ(2));
t641 = cos(pkin(11));
t594 = t641 * t549;
t580 = t541 * t594;
t525 = qJD(1) * t580;
t540 = sin(pkin(11));
t546 = sin(qJ(2));
t612 = qJD(2) * t541;
t598 = t546 * t612;
t581 = qJD(1) * t598;
t503 = qJD(2) * t525 - t540 * t581;
t539 = sin(pkin(12));
t542 = cos(pkin(12));
t545 = sin(qJ(5));
t548 = cos(qJ(5));
t522 = t539 * t548 + t542 * t545;
t655 = t503 * t522;
t562 = t540 * t549 + t546 * t641;
t614 = qJD(1) * t541;
t512 = t562 * t614;
t543 = cos(pkin(6));
t613 = qJD(1) * t543;
t531 = qJD(2) + t613;
t472 = -t512 * t542 - t531 * t539;
t626 = t539 * t512;
t661 = t531 * t542 - t626;
t662 = -t472 * t548 + t545 * t661;
t388 = qJD(5) * t662 + t655;
t547 = cos(qJ(6));
t386 = t547 * t388;
t544 = sin(qJ(6));
t460 = t548 * t661;
t423 = t472 * t545 + t460;
t416 = qJD(6) - t423;
t664 = t416 ^ 2;
t665 = -t544 * t664 + t386;
t599 = t546 * t614;
t508 = t540 * t599 - t525;
t504 = qJD(5) + t508;
t663 = t423 * t504;
t621 = t548 * t542;
t521 = t539 * t545 - t621;
t617 = t504 * t521;
t616 = t504 * t522;
t658 = pkin(2) * t598;
t657 = MDP(4) * t546;
t656 = MDP(5) * (t546 ^ 2 - t549 ^ 2);
t643 = pkin(8) + qJ(3);
t597 = t643 * t546;
t583 = t541 * t597;
t645 = pkin(1) * t549;
t491 = (pkin(2) + t645) * t543 - t583;
t646 = pkin(1) * t546;
t606 = t543 * t646;
t624 = t541 * t549;
t506 = t624 * t643 + t606;
t449 = t540 * t491 + t641 * t506;
t440 = qJ(4) * t543 + t449;
t625 = t541 * t546;
t514 = t540 * t625 - t580;
t515 = t562 * t541;
t579 = (-pkin(2) * t549 - pkin(1)) * t541;
t455 = pkin(3) * t514 - qJ(4) * t515 + t579;
t394 = -t440 * t539 + t542 * t455;
t486 = t515 * t542 + t539 * t543;
t374 = pkin(4) * t514 - pkin(9) * t486 + t394;
t395 = t542 * t440 + t539 * t455;
t485 = t515 * t539 - t543 * t542;
t380 = -pkin(9) * t485 + t395;
t654 = t545 * t374 + t548 * t380;
t493 = t506 * qJD(1);
t482 = t540 * t493;
t605 = t543 * t645;
t529 = qJD(1) * t605;
t492 = -qJD(1) * t583 + t529;
t443 = t492 * t641 - t482;
t452 = pkin(2) * t599 + pkin(3) * t512 + qJ(4) * t508;
t392 = -t443 * t539 + t542 * t452;
t644 = pkin(9) * t542;
t375 = pkin(4) * t512 + t508 * t644 + t392;
t393 = t542 * t443 + t539 * t452;
t631 = t508 * t539;
t382 = pkin(9) * t631 + t393;
t533 = pkin(2) * t540 + qJ(4);
t642 = pkin(9) + t533;
t519 = t642 * t539;
t520 = t642 * t542;
t564 = -t519 * t548 - t520 * t545;
t653 = qJD(4) * t521 - qJD(5) * t564 + t545 * t375 + t548 * t382;
t468 = -t519 * t545 + t520 * t548;
t652 = -qJD(4) * t522 - qJD(5) * t468 - t375 * t548 + t382 * t545;
t510 = qJD(2) * t515;
t502 = qJD(1) * t510;
t476 = pkin(2) * t531 + t492;
t595 = t641 * t493;
t432 = t540 * t476 + t595;
t427 = qJ(4) * t531 + t432;
t563 = qJD(1) * t579;
t516 = qJD(3) + t563;
t439 = pkin(3) * t508 - qJ(4) * t512 + t516;
t383 = -t427 * t539 + t542 * t439;
t364 = pkin(4) * t508 + pkin(9) * t472 + t383;
t384 = t542 * t427 + t539 * t439;
t370 = pkin(9) * t661 + t384;
t343 = t364 * t545 + t370 * t548;
t527 = qJD(2) * t529;
t553 = (-qJD(2) * t597 + qJD(3) * t549) * t541;
t463 = qJD(1) * t553 + t527;
t478 = -qJD(2) * t506 - qJD(3) * t625;
t464 = t478 * qJD(1);
t411 = t641 * t463 + t540 * t464;
t403 = qJD(4) * t531 + t411;
t526 = pkin(2) * t581;
t412 = pkin(3) * t502 - qJ(4) * t503 - qJD(4) * t512 + t526;
t368 = -t403 * t539 + t542 * t412;
t623 = t542 * t503;
t357 = pkin(4) * t502 - pkin(9) * t623 + t368;
t369 = t542 * t403 + t539 * t412;
t632 = t503 * t539;
t360 = -pkin(9) * t632 + t369;
t591 = -t548 * t357 + t360 * t545;
t647 = -qJD(5) * t343 - t591;
t333 = -pkin(5) * t502 - t647;
t651 = t416 * (pkin(5) * t662 + pkin(10) * t416) + t333;
t400 = t504 * t544 + t547 * t662;
t387 = qJD(5) * t460 + t503 * t621 + (qJD(5) * t472 - t632) * t545;
t590 = t387 * t544 - t547 * t502;
t351 = t400 * qJD(6) + t590;
t634 = t662 * t544;
t398 = -t547 * t504 + t634;
t650 = -t351 * t521 - t398 * t616;
t588 = t512 * t544 + t547 * t617;
t609 = qJD(6) * t544;
t556 = t522 * t609 + t588;
t649 = t522 * t386 - t416 * t556;
t648 = -t502 * t522 + t504 * t617;
t505 = t508 ^ 2;
t608 = qJD(6) * t547;
t600 = t547 * t387 + t544 * t502 + t504 * t608;
t350 = -t609 * t662 + t600;
t640 = t350 * t544;
t638 = t398 * t416;
t637 = t400 * t416;
t636 = t423 * t512;
t635 = t662 * t512;
t511 = (-t540 * t546 + t594) * t612;
t630 = t511 * t539;
t629 = t522 * t547;
t536 = t541 ^ 2;
t550 = qJD(1) ^ 2;
t627 = t536 * t550;
t622 = t544 * t388;
t530 = qJD(2) * t605;
t477 = t530 + t553;
t430 = t641 * t477 + t540 * t478;
t418 = qJD(4) * t543 + t430;
t428 = pkin(3) * t510 - qJ(4) * t511 - qJD(4) * t515 + t658;
t377 = t542 * t418 + t539 * t428;
t618 = pkin(5) * t512 - t652;
t611 = qJD(5) * t545;
t610 = qJD(5) * t548;
t604 = t536 * t646;
t602 = t549 * t627;
t601 = MDP(6) * t624;
t596 = qJD(1) * qJD(2) * t536;
t559 = t545 * t357 + t548 * t360 + t364 * t610 - t370 * t611;
t332 = pkin(10) * t502 + t559;
t410 = t463 * t540 - t641 * t464;
t397 = pkin(4) * t632 + t410;
t341 = pkin(5) * t388 - pkin(10) * t387 + t397;
t592 = -t332 * t544 + t547 * t341;
t589 = -t547 * t512 + t544 * t617;
t376 = -t418 * t539 + t542 * t428;
t429 = t477 * t540 - t641 * t478;
t442 = t492 * t540 + t595;
t585 = t416 * t547;
t584 = t531 + t613;
t582 = t549 * t596;
t535 = -t641 * pkin(2) - pkin(3);
t578 = t350 * t521 + t400 * t616;
t577 = -t521 * t502 - t504 * t616;
t401 = pkin(4) * t630 + t429;
t409 = -pkin(4) * t631 + t442;
t524 = -t542 * pkin(4) + t535;
t461 = t521 * pkin(5) - t522 * pkin(10) + t524;
t576 = pkin(10) * t512 - qJD(6) * t461 + t653;
t575 = -pkin(5) * t616 - pkin(10) * t617 + qJD(6) * t468 + t409;
t431 = t476 * t641 - t482;
t448 = t491 * t641 - t540 * t506;
t574 = t332 * t547 + t341 * t544;
t339 = pkin(10) * t504 + t343;
t426 = -t531 * pkin(3) + qJD(4) - t431;
t396 = -pkin(4) * t661 + t426;
t352 = -pkin(5) * t423 - pkin(10) * t662 + t396;
t337 = t339 * t547 + t352 * t544;
t573 = t339 * t544 - t352 * t547;
t345 = pkin(10) * t514 + t654;
t441 = -t543 * pkin(3) - t448;
t406 = t485 * pkin(4) + t441;
t433 = t548 * t485 + t486 * t545;
t434 = -t485 * t545 + t486 * t548;
t358 = t433 * pkin(5) - t434 * pkin(10) + t406;
t572 = t345 * t547 + t358 * t544;
t571 = -t345 * t544 + t358 * t547;
t363 = pkin(4) * t510 - t511 * t644 + t376;
t367 = -pkin(9) * t630 + t377;
t570 = t363 * t548 - t367 * t545;
t342 = t364 * t548 - t370 * t545;
t569 = t374 * t548 - t380 * t545;
t566 = -t383 * t539 + t384 * t542;
t565 = t426 * t511 + t441 * t503;
t405 = t434 * t547 + t514 * t544;
t404 = t434 * t544 - t514 * t547;
t561 = -pkin(8) * t624 - t606;
t560 = -pkin(8) * t581 + t527;
t558 = t545 * t363 + t548 * t367 + t374 * t610 - t380 * t611;
t557 = t522 * t608 - t589;
t555 = t561 * t531;
t338 = -pkin(5) * t504 - t342;
t554 = -pkin(10) * t388 + (t338 + t342) * t416;
t552 = -t533 * t502 + t535 * t503 + (-qJD(4) + t426) * t508;
t551 = -t416 * t557 - t522 * t622;
t391 = qJD(5) * t434 + t511 * t522;
t390 = -qJD(5) * t433 - t511 * t521;
t354 = qJD(6) * t405 + t390 * t544 - t510 * t547;
t353 = -qJD(6) * t404 + t390 * t547 + t510 * t544;
t346 = pkin(5) * t391 - pkin(10) * t390 + t401;
t344 = -pkin(5) * t514 - t569;
t335 = -pkin(5) * t510 + qJD(5) * t654 - t570;
t334 = pkin(10) * t510 + t558;
t331 = -qJD(6) * t337 + t592;
t330 = -qJD(6) * t573 + t574;
t1 = [((t555 + (t543 * t561 - 0.2e1 * t604) * qJD(1)) * MDP(9) + t584 * t601) * qJD(2) + (-(qJD(6) * t571 + t334 * t547 + t346 * t544) * t416 - t572 * t388 - t330 * t433 - t337 * t391 + t335 * t400 + t344 * t350 + t333 * t405 + t338 * t353) * MDP(30) + (-t410 * t448 + t411 * t449 - t431 * t429 + t432 * t430 + (t516 + t563) * t658) * MDP(12) - t584 * MDP(7) * t598 + (-t387 * t433 - t388 * t434 + t390 * t423 - t391 * t662) * MDP(18) + (-t388 * t514 - t391 * t504 + t423 * t510 - t433 * t502) * MDP(20) + (t570 * t504 + t569 * t502 - t591 * t514 + t342 * t510 - t401 * t423 + t406 * t388 + t397 * t433 + t396 * t391 + (-t343 * t514 - t504 * t654) * qJD(5)) * MDP(22) + (t387 * t514 + t390 * t504 + t434 * t502 + t510 * t662) * MDP(19) + (t387 * t434 + t390 * t662) * MDP(17) + (-t343 * t510 + t406 * t387 + t396 * t390 + t397 * t434 + t401 * t662 - t502 * t654 - t504 * t558 - t514 * t559) * MDP(23) - 0.2e1 * t596 * t656 + 0.2e1 * t582 * t657 + (t502 * t514 + t504 * t510) * MDP(21) + (t410 * t515 - t411 * t514 + t429 * t512 - t430 * t508 - t431 * t511 - t432 * t510 - t448 * t503 - t449 * t502) * MDP(11) + (t368 * t394 + t369 * t395 + t376 * t383 + t377 * t384 + t410 * t441 + t426 * t429) * MDP(16) + (t350 * t433 + t353 * t416 + t388 * t405 + t391 * t400) * MDP(26) + (-t351 * t433 - t354 * t416 - t388 * t404 - t391 * t398) * MDP(27) + (t388 * t433 + t391 * t416) * MDP(28) + (t350 * t405 + t353 * t400) * MDP(24) + (-t350 * t404 - t351 * t405 - t353 * t398 - t354 * t400) * MDP(25) + ((-qJD(6) * t572 - t334 * t544 + t346 * t547) * t416 + t571 * t388 + t331 * t433 - t573 * t391 + t335 * t398 + t344 * t351 + t333 * t404 + t338 * t354) * MDP(29) + (-t369 * t514 - t377 * t508 - t384 * t510 - t395 * t502 + t410 * t486 - t429 * t472 + t542 * t565) * MDP(14) + (-0.2e1 * pkin(1) * t582 - (-pkin(8) * t598 + t530) * t531 - t560 * t543) * MDP(10) + (t368 * t514 + t376 * t508 + t383 * t510 + t394 * t502 + t410 * t485 - t429 * t661 + t539 * t565) * MDP(13) + (t377 * t661 - t369 * t485 + t376 * t472 - t368 * t486 + (-t383 * t542 - t384 * t539) * t511 + (-t394 * t542 - t395 * t539) * t503) * MDP(15); ((t432 - t442) * t512 + (-t431 + t443) * t508 + (-t502 * t540 - t503 * t641) * pkin(2)) * MDP(11) + (-t635 - t648) * MDP(19) - t504 * t512 * MDP(21) + (-t383 * t392 - t384 * t393 + t410 * t535 - t426 * t442 + (-t368 * t539 + t369 * t542) * t533 + t566 * qJD(4)) * MDP(16) + (-MDP(7) * t599 + qJD(1) * t601) * (qJD(2) - t531) + (-t387 * t521 - t388 * t522 - t423 * t617 - t616 * t662) * MDP(18) + (-t342 * t512 + t524 * t388 + t616 * t396 + t397 * t521 + t409 * t423 + t502 * t564 + t504 * t652) * MDP(22) + (t343 * t512 + t524 * t387 - t617 * t396 + t397 * t522 - t409 * t662 - t468 * t502 + t504 * t653) * MDP(23) + (t387 * t522 - t617 * t662) * MDP(17) + (t578 + t649) * MDP(26) + (t551 + t650) * MDP(27) + t627 * t656 - t602 * t657 + (-(t461 * t544 + t468 * t547) * t388 - t330 * t521 - t564 * t350 + t333 * t629 + (t544 * t575 + t547 * t576) * t416 + t618 * t400 - t616 * t337 - t556 * t338) * MDP(30) + ((t461 * t547 - t468 * t544) * t388 + t331 * t521 - t564 * t351 + t333 * t544 * t522 + (t544 * t576 - t547 * t575) * t416 + t618 * t398 - t616 * t573 + t557 * t338) * MDP(29) + (t384 * t512 + t393 * t508 + t410 * t539 + t442 * t472 + t542 * t552) * MDP(14) + (pkin(1) * t602 + (-pkin(8) * t599 + t529) * t531 - t560) * MDP(10) + (t550 * t604 + (qJD(2) * t561 - t555) * qJD(1)) * MDP(9) + (t388 * t521 + t416 * t616) * MDP(28) + (t350 * t629 - t400 * t556) * MDP(24) + (t577 - t636) * MDP(20) + (t589 * t400 + t588 * t398 + (-t640 - t351 * t547 + (t398 * t544 - t400 * t547) * qJD(6)) * t522) * MDP(25) + (t431 * t442 - t432 * t443 + (-t410 * t641 + t411 * t540 - t516 * t599) * pkin(2)) * MDP(12) + (-t383 * t512 - t392 * t508 - t410 * t542 + t442 * t661 + t552 * t539) * MDP(13) + (-t392 * t472 + t393 * t626 + (qJD(4) * t661 - t383 * t508 - t393 * t531 + t369) * t542 + (-qJD(4) * t472 - t384 * t508 - t368) * t539) * MDP(15); (-t512 ^ 2 - t505) * MDP(11) + (t431 * t512 + t432 * t508 + t526) * MDP(12) + (t542 * t502 - t505 * t539 + t512 * t661) * MDP(13) + (t472 * t512 - t502 * t539 - t505 * t542) * MDP(14) + ((-t472 * t539 + t542 * t661) * t508 + (-t539 ^ 2 - t542 ^ 2) * t503) * MDP(15) + (t368 * t542 + t369 * t539 - t426 * t512 + t508 * t566) * MDP(16) + (t577 + t636) * MDP(22) + (-t635 + t648) * MDP(23) + (t551 - t650) * MDP(29) + (t578 - t649) * MDP(30); (-t472 * t508 + t632) * MDP(13) + (t508 * t661 + t623) * MDP(14) + (-t472 ^ 2 - t661 ^ 2) * MDP(15) + (-t383 * t472 - t384 * t661 + t410) * MDP(16) + (t504 * t662 + t388) * MDP(22) + (t387 + t663) * MDP(23) + (-t398 * t662 + t665) * MDP(29) + (-t400 * t662 - t547 * t664 - t622) * MDP(30); -t423 ^ 2 * MDP(18) + (t387 - t663) * MDP(19) - MDP(20) * t655 + t502 * MDP(21) + (t343 * t504 + t647) * MDP(22) + (t342 * t504 - t396 * t423 - t559) * MDP(23) + (t400 * t585 + t640) * MDP(24) + ((t350 - t638) * t547 + (-t351 - t637) * t544) * MDP(25) + (t416 * t585 + t622) * MDP(26) + t665 * MDP(27) + (-pkin(5) * t351 - t343 * t398 + t554 * t544 - t547 * t651) * MDP(29) + (-pkin(5) * t350 - t343 * t400 + t544 * t651 + t554 * t547) * MDP(30) + (-t423 * MDP(17) + (-qJD(5) + t504) * MDP(20) - t396 * MDP(22) - t400 * MDP(26) + t398 * MDP(27) - t416 * MDP(28) + t573 * MDP(29) + t337 * MDP(30) + MDP(18) * t662) * t662; t400 * t398 * MDP(24) + (-t398 ^ 2 + t400 ^ 2) * MDP(25) + (t600 + t638) * MDP(26) + (-t590 + t637) * MDP(27) + t388 * MDP(28) + (t337 * t416 - t338 * t400 + t592) * MDP(29) + (t338 * t398 - t416 * t573 - t574) * MDP(30) + (-MDP(26) * t634 - MDP(27) * t400 - MDP(29) * t337 + MDP(30) * t573) * qJD(6);];
tauc  = t1;
