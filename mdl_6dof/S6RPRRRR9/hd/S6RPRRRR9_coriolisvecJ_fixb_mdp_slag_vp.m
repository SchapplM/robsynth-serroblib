% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:35
% EndTime: 2019-03-09 07:25:52
% DurationCPUTime: 10.56s
% Computational Cost: add. (5648->511), mult. (12617->709), div. (0->0), fcn. (8913->8), ass. (0->214)
t527 = sin(qJ(3));
t531 = cos(qJ(3));
t497 = pkin(3) * t527 - pkin(8) * t531 + qJ(2);
t469 = t497 * qJD(1);
t532 = -pkin(1) - pkin(7);
t511 = qJD(1) * t532 + qJD(2);
t496 = t527 * t511;
t475 = qJD(3) * pkin(8) + t496;
t526 = sin(qJ(4));
t530 = cos(qJ(4));
t416 = t530 * t469 - t475 * t526;
t604 = qJD(1) * t531;
t581 = t530 * t604;
t603 = qJD(3) * t526;
t486 = t581 + t603;
t398 = -pkin(9) * t486 + t416;
t605 = qJD(1) * t527;
t515 = qJD(4) + t605;
t393 = pkin(4) * t515 + t398;
t417 = t469 * t526 + t475 * t530;
t593 = t530 * qJD(3);
t484 = t526 * t604 - t593;
t399 = -pkin(9) * t484 + t417;
t529 = cos(qJ(5));
t397 = t529 * t399;
t525 = sin(qJ(5));
t353 = t393 * t525 + t397;
t427 = t529 * t484 + t486 * t525;
t665 = pkin(10) * t427;
t348 = t353 - t665;
t524 = sin(qJ(6));
t595 = qJD(6) * t524;
t346 = t348 * t595;
t528 = cos(qJ(6));
t551 = t484 * t525 - t529 * t486;
t644 = t551 * t524;
t381 = -t528 * t427 + t644;
t639 = t511 * t531;
t476 = -qJD(3) * pkin(3) - t639;
t438 = pkin(4) * t484 + t476;
t386 = pkin(5) * t427 + t438;
t673 = -t386 * t381 + t346;
t552 = t427 * t524 + t528 * t551;
t592 = qJD(1) * qJD(3);
t575 = t531 * t592;
t672 = MDP(32) * t575 + (-t381 ^ 2 + t552 ^ 2) * MDP(29) + t381 * MDP(28) * t552;
t598 = qJD(4) * t531;
t577 = t526 * t598;
t541 = -t527 * t593 - t577;
t442 = qJD(1) * t541 + qJD(4) * t593;
t602 = qJD(3) * t527;
t580 = t526 * t602;
t443 = -qJD(1) * t580 + t486 * qJD(4);
t596 = qJD(5) * t529;
t597 = qJD(5) * t525;
t367 = t529 * t442 - t525 * t443 - t484 * t596 - t486 * t597;
t559 = pkin(3) * t531 + pkin(8) * t527;
t482 = qJD(3) * t559 + qJD(2);
t455 = t482 * qJD(1);
t447 = t530 * t455;
t539 = -t417 * qJD(4) + t447;
t601 = qJD(3) * t531;
t358 = -pkin(9) * t442 + (pkin(4) * qJD(1) - t511 * t526) * t601 + t539;
t578 = t531 * t593;
t599 = qJD(4) * t530;
t600 = qJD(4) * t526;
t544 = t526 * t455 + t469 * t599 - t475 * t600 + t511 * t578;
t360 = -pkin(9) * t443 + t544;
t572 = t529 * t358 - t525 * t360;
t538 = -t353 * qJD(5) + t572;
t335 = pkin(5) * t575 - pkin(10) * t367 + t538;
t536 = qJD(5) * t551 - t442 * t525 - t529 * t443;
t561 = -t525 * t358 - t529 * t360 - t393 * t596 + t399 * t597;
t336 = pkin(10) * t536 - t561;
t671 = -t524 * t335 - t528 * t336 + t673;
t594 = qJD(6) * t528;
t584 = t528 * t367 - t427 * t594 + t524 * t536;
t341 = t551 * t595 + t584;
t509 = qJD(5) + t515;
t571 = t367 * t524 - t528 * t536;
t537 = qJD(6) * t552 - t571;
t503 = qJD(6) + t509;
t661 = t503 * t552;
t662 = t381 * t503;
t670 = MDP(25) * t575 + (-t427 ^ 2 + t551 ^ 2) * MDP(22) + (t427 * t509 + t367) * MDP(23) + (-t509 * t551 + t536) * MDP(24) - t427 * MDP(21) * t551 + (t537 - t661) * MDP(31) + (t341 - t662) * MDP(30) + t672;
t573 = t528 * t335 - t524 * t336;
t656 = t386 * t552 + t573;
t493 = t559 * qJD(1);
t474 = t530 * t493;
t646 = pkin(8) + pkin(9);
t583 = qJD(4) * t646;
t633 = t526 * t531;
t586 = t511 * t633;
t631 = t527 * t530;
t590 = pkin(9) * t631;
t667 = -t586 + t474 + (pkin(4) * t531 + t590) * qJD(1) + t530 * t583;
t582 = t526 * t605;
t627 = t530 * t531;
t613 = t526 * t493 + t511 * t627;
t666 = pkin(9) * t582 + t526 * t583 + t613;
t664 = pkin(10) * t551;
t489 = t525 * t530 + t526 * t529;
t647 = qJD(4) + qJD(5);
t659 = t647 * t489;
t488 = t525 * t526 - t529 * t530;
t650 = t527 * t488;
t615 = -qJD(1) * t650 - t647 * t488;
t467 = t489 * qJD(1);
t614 = t527 * t467 + t659;
t576 = t530 * t598;
t657 = t576 - t580;
t655 = t427 * t438 + t561;
t654 = t438 * t551 + t538;
t591 = 0.2e1 * qJD(1);
t523 = t531 ^ 2;
t651 = MDP(8) * (t527 ^ 2 - t523);
t481 = t530 * t497;
t632 = t526 * t532;
t574 = pkin(4) - t632;
t424 = -pkin(9) * t627 + t527 * t574 + t481;
t630 = t527 * t532;
t508 = t530 * t630;
t611 = t526 * t497 + t508;
t437 = -pkin(9) * t633 + t611;
t616 = t525 * t424 + t529 * t437;
t560 = -t496 + (t582 + t600) * pkin(4);
t649 = t667 * t529;
t504 = t646 * t526;
t505 = t646 * t530;
t612 = -t525 * t504 + t529 * t505;
t648 = t504 * t596 + t505 * t597 + t525 * t667 + t666 * t529;
t534 = qJD(1) ^ 2;
t645 = qJ(2) * t534;
t643 = t442 * t526;
t642 = t476 * t526;
t641 = t484 * t515;
t640 = t486 * t515;
t638 = t515 * t526;
t637 = t515 * t527;
t636 = t515 * t530;
t635 = t524 * t525;
t395 = t525 * t399;
t634 = t525 * t528;
t352 = t529 * t393 - t395;
t347 = t352 + t664;
t345 = pkin(5) * t509 + t347;
t629 = t528 * t345;
t628 = t528 * t348;
t626 = t531 * t532;
t533 = qJD(3) ^ 2;
t624 = t532 * t533;
t433 = t528 * t488 + t489 * t524;
t623 = -qJD(6) * t433 - t524 * t614 + t528 * t615;
t434 = -t488 * t524 + t489 * t528;
t622 = qJD(6) * t434 + t524 * t615 + t528 * t614;
t459 = t488 * t531;
t621 = qJD(3) * t459 + t527 * t659 + t467;
t620 = t488 * qJD(1) - t489 * t601 + t647 * t650;
t619 = t529 * t398 - t395;
t617 = pkin(5) * t614 + t560;
t588 = pkin(4) * qJD(5) * t503;
t585 = t526 * t630;
t520 = -pkin(4) * t530 - pkin(3);
t412 = pkin(4) * t443 + t511 * t602;
t462 = t530 * t482;
t375 = t462 + (-t508 + (pkin(9) * t531 - t497) * t526) * qJD(4) + (t531 * t574 + t590) * qJD(3);
t540 = -qJD(4) * t585 + t526 * t482 + t497 * t599 + t532 * t578;
t384 = -pkin(9) * t657 + t540;
t570 = t529 * t375 - t384 * t525;
t569 = -t398 * t525 - t397;
t568 = t529 * t424 - t437 * t525;
t566 = -t529 * t504 - t505 * t525;
t483 = pkin(4) * t633 - t626;
t565 = t484 + t593;
t564 = -t486 + t603;
t563 = qJD(6) * t345 + t336;
t562 = qJD(4) * t527 + qJD(1);
t507 = t527 * t575;
t411 = -pkin(10) * t488 + t612;
t558 = pkin(5) * t604 + pkin(10) * t615 + qJD(5) * t612 + qJD(6) * t411 - t525 * t666 + t649;
t410 = -pkin(10) * t489 + t566;
t557 = pkin(10) * t614 - qJD(6) * t410 + t648;
t456 = t489 * t527;
t556 = qJD(6) * t456 + t621;
t555 = -qJD(6) * t650 - t620;
t338 = t524 * t345 + t628;
t363 = pkin(5) * t527 + pkin(10) * t459 + t568;
t457 = t489 * t531;
t371 = -pkin(10) * t457 + t616;
t554 = t363 * t524 + t371 * t528;
t406 = t528 * t457 - t459 * t524;
t407 = -t457 * t524 - t459 * t528;
t550 = qJD(1) * t523 - t637;
t519 = pkin(4) * t529 + pkin(5);
t549 = pkin(4) * t634 + t519 * t524;
t548 = -pkin(4) * t635 + t519 * t528;
t547 = -pkin(8) * t601 + t476 * t527;
t444 = pkin(4) * t657 + t532 * t602;
t543 = t525 * t375 + t529 * t384 + t424 * t596 - t437 * t597;
t452 = pkin(5) * t488 + t520;
t432 = pkin(5) * t457 + t483;
t400 = pkin(4) * t486 - pkin(5) * t551;
t392 = -t597 * t633 + (t627 * t647 - t580) * t529 + t541 * t525;
t390 = qJD(3) * t650 - t531 * t659;
t372 = pkin(5) * t392 + t444;
t351 = -pkin(5) * t536 + t412;
t350 = t619 + t664;
t349 = t569 + t665;
t344 = qJD(6) * t407 + t390 * t524 + t528 * t392;
t343 = -qJD(6) * t406 + t390 * t528 - t392 * t524;
t340 = -pkin(10) * t392 + t543;
t339 = pkin(5) * t601 - pkin(10) * t390 - qJD(5) * t616 + t570;
t337 = -t348 * t524 + t629;
t1 = [(-MDP(10) * t531 - MDP(9) * t527) * t533 + (-t515 * t576 - t443 * t527 + (-t484 * t531 - t526 * t550) * qJD(3)) * MDP(17) + (-t515 * t577 + t442 * t527 + (t486 * t531 + t530 * t550) * qJD(3)) * MDP(16) + (t442 * t627 + t486 * t541) * MDP(14) + ((t484 * t530 + t486 * t526) * t602 + (-t643 - t443 * t530 + (t484 * t526 - t486 * t530) * qJD(4)) * t531) * MDP(15) + 0.2e1 * t592 * t651 + ((t339 * t528 - t340 * t524) * t503 + t573 * t527 - t372 * t381 - t432 * t537 + t351 * t406 + t386 * t344 + (-t338 * t527 - t503 * t554) * qJD(6) + ((t363 * t528 - t371 * t524) * qJD(1) + t337) * t601) * MDP(33) + (t537 * t527 - t344 * t503 + (-qJD(1) * t406 + t381) * t601) * MDP(31) + (-t341 * t406 + t343 * t381 + t344 * t552 + t407 * t537) * MDP(29) + (-t540 * t515 - t544 * t527 + (-t532 * t442 - t476 * t600) * t531 + ((-qJD(1) * t611 - t417) * t531 + (t532 * t486 + (-t476 + t639) * t530) * t527) * qJD(3)) * MDP(20) + (-t443 * t626 + t447 * t527 + t462 * t515 + (-t417 * t527 + t476 * t627 - t515 * t611) * qJD(4) + ((t484 * t532 - t642) * t527 + (-t515 * t632 + (t481 - t585) * qJD(1) + t416) * t531) * qJD(3)) * MDP(19) + (t536 * t527 - t392 * t509 + (-qJD(1) * t457 - t427) * t601) * MDP(24) + (-t367 * t457 - t390 * t427 + t392 * t551 - t459 * t536) * MDP(22) + (t570 * t509 + t572 * t527 + t444 * t427 - t483 * t536 + t412 * t457 + t438 * t392 + (-t353 * t527 - t509 * t616) * qJD(5) + (qJD(1) * t568 + t352) * t601) * MDP(26) + (-t531 * t624 + (-qJ(2) * t602 + qJD(2) * t531) * t591) * MDP(13) + (-t527 * t624 + (qJ(2) * t601 + qJD(2) * t527) * t591) * MDP(12) + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t591 + (t515 * t601 + t507) * MDP(18) + (t509 * t601 + t507) * MDP(25) + (t503 * t601 + t507) * MDP(32) + (t341 * t407 - t343 * t552) * MDP(28) + (-t543 * t509 + t561 * t527 - t444 * t551 + t483 * t367 - t412 * t459 + t438 * t390 + (-qJD(1) * t616 - t353) * t601) * MDP(27) + (t367 * t527 + t390 * t509 + (-qJD(1) * t459 - t551) * t601) * MDP(23) + (-t367 * t459 - t390 * t551) * MDP(21) + (t341 * t527 + t343 * t503 + (qJD(1) * t407 - t552) * t601) * MDP(30) + (t432 * t341 + t386 * t343 + t346 * t527 + t351 * t407 - t372 * t552 + (-(-qJD(6) * t371 + t339) * t503 - t335 * t527) * t524 + (-(qJD(6) * t363 + t340) * t503 - t563 * t527) * t528 + (-qJD(1) * t554 - t338) * t601) * MDP(34) - 0.2e1 * MDP(7) * t507; -t534 * MDP(5) - MDP(6) * t645 + (-t443 * t531 - t562 * t636 + (t484 * t527 + (-t515 - t605) * t633) * qJD(3)) * MDP(19) + (-t442 * t531 + t562 * t638 + (-t515 * t627 + (t486 - t581) * t527) * qJD(3)) * MDP(20) + (t536 * t531 + t620 * t509 + (t427 * t527 - t456 * t604) * qJD(3)) * MDP(26) + (-t367 * t531 + t621 * t509 + (-t527 * t551 + t604 * t650) * qJD(3)) * MDP(27) + (t531 * t537 + (t524 * t556 - t528 * t555) * t503 + ((-t456 * t528 + t524 * t650) * t604 - t527 * t381) * qJD(3)) * MDP(33) + (-t531 * t341 + (t524 * t555 + t528 * t556) * t503 + (-(-t456 * t524 - t528 * t650) * t604 - t527 * t552) * qJD(3)) * MDP(34) + (t527 * MDP(12) + t531 * MDP(13)) * (-t533 - t534); (t412 * t488 + t560 * t427 + t614 * t438 - t520 * t536) * MDP(26) + (t367 * t489 - t551 * t615) * MDP(21) + (-t367 * t488 - t427 * t615 + t489 * t536 + t551 * t614) * MDP(22) + (t486 * t636 + t643) * MDP(14) + (t515 * t599 + (t515 * t631 + t531 * t564) * qJD(1)) * MDP(16) + (-t515 * t600 + (-t526 * t637 + t531 * t565) * qJD(1)) * MDP(17) + (-pkin(3) * t442 + t613 * t515 + t564 * t496 + (pkin(8) * t638 + t476 * t530) * qJD(4) + (t417 * t531 + t530 * t547) * qJD(1)) * MDP(20) + ((t442 - t641) * t530 + (-t443 - t640) * t526) * MDP(15) + (-pkin(3) * t443 - t474 * t515 + (t515 * t633 - t527 * t565) * t511 + (-pkin(8) * t636 + t642) * qJD(4) + (-t416 * t531 + t526 * t547) * qJD(1)) * MDP(19) + (t520 * t367 + t412 * t489 + t615 * t438 - t551 * t560) * MDP(27) + (t351 * t433 - t381 * t617 + t622 * t386 - t452 * t537) * MDP(33) + (t452 * t341 + t351 * t434 + t623 * t386 - t552 * t617) * MDP(34) + (t341 * t434 - t552 * t623) * MDP(28) + (-t341 * t433 + t381 * t623 + t434 * t537 + t552 * t622) * MDP(29) + t527 * MDP(13) * t645 + ((-t505 * t596 + (qJD(5) * t504 + t666) * t525 - t649) * MDP(26) - t614 * MDP(24) + t615 * MDP(23) + t648 * MDP(27)) * t509 + ((t524 * t557 - t528 * t558) * MDP(33) - t622 * MDP(31) + (t524 * t558 + t528 * t557) * MDP(34) + t623 * MDP(30)) * t503 + ((qJD(3) * t566 - t352) * MDP(26) + (-qJD(3) * t488 + t427) * MDP(24) + (qJD(3) * t489 + t551) * MDP(23) + (-qJD(3) * t612 + t353) * MDP(27) + ((t410 * t528 - t411 * t524) * qJD(3) - t337) * MDP(33) + (-qJD(3) * t433 - t381) * MDP(31) + (-(t410 * t524 + t411 * t528) * qJD(3) + t338) * MDP(34) + (qJD(3) * t434 + t552) * MDP(30) - t515 * MDP(18) - t509 * MDP(25) - t503 * MDP(32)) * t604 + (-t651 + (-qJ(2) * MDP(12) + t527 * MDP(7)) * t531) * t534; (-t443 + t640) * MDP(17) + (t416 * t515 + t476 * t484 - t544) * MDP(20) + (t442 + t641) * MDP(16) + (-qJD(3) * t586 + t417 * t515 - t476 * t486 + t539) * MDP(19) + (t548 * t575 - (t349 * t528 - t350 * t524) * t503 + t400 * t381 + (-t524 * t529 - t634) * t588 + (-t503 * t549 - t338) * qJD(6) + t656) * MDP(33) + (-t484 ^ 2 + t486 ^ 2) * MDP(15) + (-t569 * t509 + (-t427 * t486 - t509 * t597 + t529 * t575) * pkin(4) + t654) * MDP(26) + (t619 * t509 + (t486 * t551 - t509 * t596 - t525 * t575) * pkin(4) + t655) * MDP(27) + (-t549 * t575 + (t349 * t524 + t350 * t528) * t503 + t400 * t552 - (t528 * t529 - t635) * t588 + (-t503 * t548 - t629) * qJD(6) + t671) * MDP(34) + t486 * t484 * MDP(14) + MDP(18) * t575 + t670; (t353 * t509 + t654) * MDP(26) + (t352 * t509 + t655) * MDP(27) + (-(-t347 * t524 - t628) * t503 - t338 * qJD(6) + (-t381 * t551 - t503 * t595 + t528 * t575) * pkin(5) + t656) * MDP(33) + ((-t348 * t503 - t335) * t524 + (t347 * t503 - t563) * t528 + (-t503 * t594 - t524 * t575 - t551 * t552) * pkin(5) + t673) * MDP(34) + t670; (t584 - t662) * MDP(30) + (-t571 - t661) * MDP(31) + (t338 * t503 + t656) * MDP(33) + (t337 * t503 + t671) * MDP(34) + (MDP(30) * t644 + MDP(31) * t552 - MDP(33) * t338 - MDP(34) * t629) * qJD(6) + t672;];
tauc  = t1;
