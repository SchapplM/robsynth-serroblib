% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:44
% EndTime: 2019-03-09 07:02:59
% DurationCPUTime: 9.93s
% Computational Cost: add. (5806->496), mult. (13697->667), div. (0->0), fcn. (9797->10), ass. (0->220)
t540 = sin(pkin(11)) * pkin(1) + pkin(7);
t526 = t540 * qJD(1);
t556 = sin(qJ(3));
t547 = t556 * qJD(2);
t560 = cos(qJ(3));
t482 = t560 * t526 + t547;
t472 = qJD(3) * pkin(8) + t482;
t541 = -cos(pkin(11)) * pkin(1) - pkin(2);
t503 = -pkin(3) * t560 - pkin(8) * t556 + t541;
t476 = t503 * qJD(1);
t555 = sin(qJ(4));
t559 = cos(qJ(4));
t420 = -t472 * t555 + t559 * t476;
t631 = qJD(3) * t555;
t634 = qJD(1) * t556;
t511 = t559 * t634 + t631;
t412 = -pkin(9) * t511 + t420;
t633 = qJD(1) * t560;
t538 = -qJD(4) + t633;
t405 = -pkin(4) * t538 + t412;
t662 = t555 * t476;
t421 = t472 * t559 + t662;
t622 = t559 * qJD(3);
t509 = t555 * t634 - t622;
t413 = -pkin(9) * t509 + t421;
t558 = cos(qJ(5));
t411 = t558 * t413;
t554 = sin(qJ(5));
t369 = t405 * t554 + t411;
t450 = t558 * t509 + t511 * t554;
t699 = pkin(10) * t450;
t363 = t369 - t699;
t553 = sin(qJ(6));
t624 = qJD(6) * t553;
t361 = t363 * t624;
t557 = cos(qJ(6));
t580 = t509 * t554 - t558 * t511;
t675 = t580 * t553;
t403 = -t557 * t450 + t675;
t684 = qJD(2) * t560 - t556 * t526;
t471 = -qJD(3) * pkin(3) - t684;
t447 = pkin(4) * t509 + t471;
t406 = pkin(5) * t450 + t447;
t706 = -t406 * t403 + t361;
t581 = t450 * t553 + t557 * t580;
t621 = qJD(1) * qJD(3);
t609 = t556 * t621;
t705 = MDP(30) * t609 + (-t403 ^ 2 + t581 ^ 2) * MDP(27) + t403 * MDP(26) * t581;
t608 = t560 * t621;
t628 = qJD(4) * t555;
t611 = t556 * t628;
t620 = qJD(3) * qJD(4);
t465 = -qJD(1) * t611 + (t608 + t620) * t559;
t627 = qJD(4) * t559;
t610 = t556 * t627;
t629 = qJD(3) * t560;
t613 = t555 * t629;
t567 = t610 + t613;
t466 = qJD(1) * t567 + t555 * t620;
t625 = qJD(5) * t558;
t626 = qJD(5) * t554;
t385 = t558 * t465 - t554 * t466 - t509 * t625 - t511 * t626;
t473 = t684 * qJD(3);
t586 = pkin(3) * t556 - pkin(8) * t560;
t521 = t586 * qJD(3);
t502 = qJD(1) * t521;
t594 = t555 * t473 - t559 * t502;
t565 = -qJD(4) * t421 - t594;
t374 = pkin(4) * t609 - pkin(9) * t465 + t565;
t616 = t559 * t473 + t476 * t627 + t555 * t502;
t571 = -t472 * t628 + t616;
t377 = -pkin(9) * t466 + t571;
t603 = t558 * t374 - t554 * t377;
t566 = -qJD(5) * t369 + t603;
t351 = pkin(5) * t609 - pkin(10) * t385 + t566;
t563 = qJD(5) * t580 - t465 * t554 - t558 * t466;
t590 = -t554 * t374 - t558 * t377 - t405 * t625 + t413 * t626;
t352 = pkin(10) * t563 - t590;
t704 = -t553 * t351 - t557 * t352 + t706;
t623 = qJD(6) * t557;
t617 = t557 * t385 - t450 * t623 + t553 * t563;
t357 = t580 * t624 + t617;
t532 = -qJD(5) + t538;
t602 = t385 * t553 - t557 * t563;
t564 = qJD(6) * t581 - t602;
t525 = -qJD(6) + t532;
t695 = t525 * t581;
t696 = t403 * t525;
t703 = MDP(23) * t609 + (-t450 ^ 2 + t580 ^ 2) * MDP(20) + (-t450 * t532 + t385) * MDP(21) + (t532 * t580 + t563) * MDP(22) - t450 * MDP(19) * t580 + (t564 + t695) * MDP(29) + (t357 + t696) * MDP(28) + t705;
t605 = t557 * t351 - t553 * t352;
t692 = t406 * t581 + t605;
t655 = t559 * t560;
t576 = pkin(4) * t556 - pkin(9) * t655;
t518 = t586 * qJD(1);
t593 = t559 * t518 - t555 * t684;
t680 = pkin(8) + pkin(9);
t615 = qJD(4) * t680;
t700 = qJD(1) * t576 + t559 * t615 + t593;
t614 = t555 * t633;
t643 = t555 * t518 + t559 * t684;
t693 = pkin(9) * t614 - t555 * t615 - t643;
t698 = pkin(10) * t580;
t512 = t554 * t555 - t558 * t559;
t572 = t512 * t560;
t681 = qJD(4) + qJD(5);
t645 = qJD(1) * t572 - t681 * t512;
t513 = t554 * t559 + t555 * t558;
t644 = (-t633 + t681) * t513;
t691 = t447 * t450 + t590;
t690 = t447 * t580 + t566;
t687 = MDP(5) * t556;
t549 = t556 ^ 2;
t686 = MDP(6) * (-t560 ^ 2 + t549);
t485 = t513 * t556;
t685 = t700 * t558;
t488 = t559 * t503;
t659 = t556 * t559;
t665 = t540 * t555;
t435 = -pkin(9) * t659 + t488 + (-pkin(4) - t665) * t560;
t517 = t540 * t655;
t639 = t555 * t503 + t517;
t661 = t555 * t556;
t442 = -pkin(9) * t661 + t639;
t647 = t554 * t435 + t558 * t442;
t587 = -t482 + (-t614 + t628) * pkin(4);
t529 = t680 * t555;
t530 = t680 * t559;
t640 = -t554 * t529 + t558 * t530;
t683 = -t529 * t625 - t530 * t626 - t554 * t700 + t693 * t558;
t682 = t560 * t622 - t611;
t679 = t564 * t560;
t415 = -qJD(3) * t572 - t485 * t681;
t416 = -t626 * t661 + (t659 * t681 + t613) * t558 + t682 * t554;
t486 = t512 * t556;
t433 = t557 * t485 - t486 * t553;
t364 = -qJD(6) * t433 + t415 * t557 - t416 * t553;
t678 = t364 * t525;
t677 = t563 * t560;
t676 = t415 * t532;
t674 = t465 * t555;
t673 = t466 * t560;
t672 = t471 * t555;
t671 = t471 * t559;
t474 = qJD(3) * t547 + t526 * t629;
t670 = t474 * t555;
t669 = t474 * t559;
t668 = t509 * t538;
t667 = t511 * t538;
t666 = t538 * t559;
t664 = t553 * t554;
t409 = t554 * t413;
t663 = t554 * t557;
t660 = t555 * t560;
t561 = qJD(3) ^ 2;
t658 = t556 * t561;
t368 = t558 * t405 - t409;
t362 = t368 + t698;
t359 = -pkin(5) * t532 + t362;
t657 = t557 * t359;
t656 = t557 * t363;
t654 = t560 * t561;
t434 = -t485 * t553 - t486 * t557;
t365 = qJD(6) * t434 + t415 * t553 + t557 * t416;
t653 = t365 * t525 - t433 * t609;
t454 = t557 * t512 + t513 * t553;
t652 = -qJD(6) * t454 - t553 * t644 + t557 * t645;
t455 = -t512 * t553 + t513 * t557;
t651 = qJD(6) * t455 + t553 * t645 + t557 * t644;
t650 = t416 * t532 - t485 * t609;
t649 = t558 * t412 - t409;
t646 = pkin(5) * t644 + t587;
t642 = t503 * t627 + t555 * t521;
t630 = qJD(3) * t556;
t641 = t559 * t521 + t630 * t665;
t491 = pkin(4) * t661 + t556 * t540;
t527 = qJD(1) * t541;
t618 = pkin(4) * qJD(5) * t525;
t458 = pkin(4) * t567 + t540 * t629;
t546 = -pkin(4) * t559 - pkin(3);
t606 = MDP(16) * t630;
t604 = -t357 * t560 - t581 * t630;
t601 = -t385 * t560 - t580 * t630;
t392 = t576 * qJD(3) + (-t517 + (pkin(9) * t556 - t503) * t555) * qJD(4) + t641;
t394 = (-t556 * t622 - t560 * t628) * t540 - t567 * pkin(9) + t642;
t600 = t558 * t392 - t394 * t554;
t599 = -t412 * t554 - t411;
t598 = t558 * t435 - t442 * t554;
t596 = -t465 * t560 + t511 * t630;
t595 = t538 * t540 + t472;
t592 = -t558 * t529 - t530 * t554;
t591 = qJD(6) * t359 + t352;
t589 = t538 * t611;
t588 = t538 * t610;
t428 = pkin(4) * t466 + t474;
t438 = -pkin(10) * t512 + t640;
t585 = pkin(5) * t634 + pkin(10) * t645 + qJD(5) * t640 + qJD(6) * t438 + t554 * t693 + t685;
t437 = -pkin(10) * t513 + t592;
t584 = -pkin(10) * t644 + qJD(6) * t437 + t683;
t354 = t553 * t359 + t656;
t378 = -pkin(5) * t560 + pkin(10) * t486 + t598;
t379 = -pkin(10) * t485 + t647;
t583 = t378 * t553 + t379 * t557;
t579 = qJD(1) * t549 - t538 * t560;
t577 = 0.2e1 * qJD(3) * t527;
t545 = pkin(4) * t558 + pkin(5);
t575 = pkin(4) * t663 + t545 * t553;
t574 = -pkin(4) * t664 + t545 * t557;
t570 = t554 * t392 + t558 * t394 + t435 * t625 - t442 * t626;
t568 = t579 * t555;
t479 = pkin(5) * t512 + t546;
t443 = pkin(5) * t485 + t491;
t422 = pkin(4) * t511 - pkin(5) * t580;
t389 = pkin(5) * t416 + t458;
t370 = -pkin(5) * t563 + t428;
t367 = t649 + t698;
t366 = t599 + t699;
t356 = -pkin(10) * t416 + t570;
t355 = pkin(5) * t630 - pkin(10) * t415 - qJD(5) * t647 + t600;
t353 = -t363 * t553 + t657;
t1 = [(t443 * t357 - t361 * t560 + t406 * t364 + t370 * t434 - t389 * t581 + ((-qJD(6) * t379 + t355) * t525 + t351 * t560) * t553 + ((qJD(6) * t378 + t356) * t525 + t591 * t560) * t557) * MDP(32) + (t357 * t434 - t364 * t581) * MDP(26) + (t491 * t385 + t447 * t415 - t428 * t486 - t458 * t580 + t570 * t532 - t590 * t560) * MDP(25) + (-t385 * t486 - t415 * t580) * MDP(19) + (t650 - t677) * MDP(22) + (t653 - t679) * MDP(29) + (-t486 * t609 + t601 - t676) * MDP(21) + (t588 + t673 + (-t509 * t556 - t568) * qJD(3)) * MDP(15) + (t465 * t659 + t511 * t682) * MDP(12) + (t434 * t609 + t604 - t678) * MDP(28) + (-t357 * t433 + t364 * t403 + t365 * t581 + t434 * t564) * MDP(27) + (-(t355 * t557 - t356 * t553) * t525 - t605 * t560 - t389 * t403 - t443 * t564 + t370 * t433 + t406 * t365 + (t354 * t560 + t525 * t583) * qJD(6)) * MDP(31) + ((-t532 - t633) * MDP(23) + (-t525 - t633) * MDP(30) + t403 * MDP(29) + (-qJD(1) * t647 - t369) * MDP(25) + ((t378 * t557 - t379 * t553) * qJD(1) + t353) * MDP(31) + (qJD(1) * t598 + t368) * MDP(24) + (-qJD(1) * t583 - t354) * MDP(32) - t450 * MDP(22)) * t630 - MDP(8) * t658 + (t540 * t658 + t560 * t577) * MDP(11) + (-t540 * t654 + t556 * t577) * MDP(10) + (-t538 - t633) * t606 + (t579 * t622 + t589 + t596) * MDP(14) + (-t385 * t485 - t415 * t450 + t416 * t580 - t486 * t563) * MDP(20) + (-t600 * t532 - t603 * t560 + t458 * t450 - t491 * t563 + t428 * t485 + t447 * t416 + (t369 * t560 + t532 * t647) * qJD(5)) * MDP(24) - 0.2e1 * t621 * t686 + 0.2e1 * t608 * t687 + MDP(7) * t654 + (t642 * t538 + (-t595 * t628 + (t511 * t540 + t671) * qJD(3) + t616) * t560 + (-t471 * t628 + t540 * t465 + t669 + (-qJD(1) * t639 - t540 * t666 - t421) * qJD(3)) * t556) * MDP(18) + (-(-t503 * t628 + t641) * t538 + ((t509 * t540 + t672) * qJD(3) + (t559 * t595 + t662) * qJD(4) + t594) * t560 + (t471 * t627 + t540 * t466 + t670 + ((-t540 * t660 + t488) * qJD(1) + t420) * qJD(3)) * t556) * MDP(17) + ((-t509 * t559 - t511 * t555) * t629 + (-t674 - t466 * t559 + (t509 * t555 - t511 * t559) * qJD(4)) * t556) * MDP(13); (t588 - t673) * MDP(17) + (-t589 + t596) * MDP(18) + (t650 + t677) * MDP(24) + (t601 + t676) * MDP(25) + (t653 + t679) * MDP(31) + (t604 + t678) * MDP(32) + (-MDP(10) * t556 - MDP(11) * t560) * t561 + (-t579 * MDP(18) * t559 - MDP(17) * t568 + (t509 * MDP(17) + t450 * MDP(24) - t403 * MDP(31) + (MDP(25) * t486 - MDP(32) * t434) * qJD(1)) * t556) * qJD(3); (-t538 * t627 + (t538 * t655 + (-t511 + t631) * t556) * qJD(1)) * MDP(14) + (-t385 * t512 - t450 * t645 + t513 * t563 + t580 * t644) * MDP(20) + (t546 * t385 + t428 * t513 + t645 * t447 - t580 * t587) * MDP(25) + (t370 * t454 - t403 * t646 + t651 * t406 - t479 * t564) * MDP(31) + (t357 * t455 - t581 * t652) * MDP(26) + (t479 * t357 + t370 * t455 + t652 * t406 - t581 * t646) * MDP(32) + (-t357 * t454 + t403 * t652 + t455 * t564 + t581 * t651) * MDP(27) + (qJD(3) * t482 - t474) * MDP(10) + (t428 * t512 + t644 * t447 + t587 * t450 - t546 * t563) * MDP(24) + (t385 * t513 - t580 * t645) * MDP(19) - t527 * t633 * MDP(11) + (t538 * t628 + (-t538 * t660 + (t509 + t622) * t556) * qJD(1)) * MDP(15) + ((t465 + t668) * t559 + (-t466 + t667) * t555) * MDP(13) + (-pkin(3) * t465 + t670 - t643 * t538 - t482 * t511 + (-pkin(8) * t538 * t555 + t671) * qJD(4) + (-t471 * t655 + (-pkin(8) * t622 + t421) * t556) * qJD(1)) * MDP(18) + (-pkin(3) * t466 - t669 + t593 * t538 - t482 * t509 + (pkin(8) * t666 + t672) * qJD(4) + (-t420 * t556 + (-pkin(8) * t630 - t471 * t560) * t555) * qJD(1)) * MDP(17) + (-t511 * t666 + t674) * MDP(12) + (-t560 * t687 + t686) * qJD(1) ^ 2 + (-t645 * MDP(21) + t683 * MDP(25) + t644 * MDP(22) + (t530 * t625 + (-qJD(5) * t529 + t693) * t554 + t685) * MDP(24)) * t532 + ((t553 * t584 + t557 * t585) * MDP(31) + t651 * MDP(29) + (-t553 * t585 + t557 * t584) * MDP(32) - t652 * MDP(28)) * t525 + ((qJD(3) * t513 + t580) * MDP(21) + (-qJD(3) * t640 + t369) * MDP(25) + ((t437 * t557 - t438 * t553) * qJD(3) - t353) * MDP(31) + (-qJD(3) * t454 - t403) * MDP(29) + (-(t437 * t553 + t438 * t557) * qJD(3) + t354) * MDP(32) + (qJD(3) * t455 + t581) * MDP(28) - t527 * MDP(10) + (-qJD(3) * t512 + t450) * MDP(22) + (qJD(3) * t592 - t368) * MDP(24) + t538 * MDP(16) + t532 * MDP(23) + t525 * MDP(30)) * t634; (-t466 - t667) * MDP(15) + t511 * t509 * MDP(12) + (-t575 * t609 - (t366 * t553 + t367 * t557) * t525 + t422 * t581 + (t557 * t558 - t664) * t618 + (t525 * t574 - t657) * qJD(6) + t704) * MDP(32) + (t599 * t532 + (-t450 * t511 + t532 * t626 + t558 * t609) * pkin(4) + t690) * MDP(24) + (-t649 * t532 + (t511 * t580 + t532 * t625 - t554 * t609) * pkin(4) + t691) * MDP(25) + (-t420 * t538 + t471 * t509 - t571) * MDP(18) + (-t421 * t538 - t471 * t511 + t565) * MDP(17) + qJD(1) * t606 + (-t509 ^ 2 + t511 ^ 2) * MDP(13) + (t574 * t609 + (t366 * t557 - t367 * t553) * t525 + t422 * t403 - (-t553 * t558 - t663) * t618 + (t525 * t575 - t354) * qJD(6) + t692) * MDP(31) + (t465 - t668) * MDP(14) + t703; (-t369 * t532 + t690) * MDP(24) + (-t368 * t532 + t691) * MDP(25) + ((-t362 * t553 - t656) * t525 - t354 * qJD(6) + (-t403 * t580 + t525 * t624 + t557 * t609) * pkin(5) + t692) * MDP(31) + ((t363 * t525 - t351) * t553 + (-t362 * t525 - t591) * t557 + (t525 * t623 - t553 * t609 - t580 * t581) * pkin(5) + t706) * MDP(32) + t703; (t617 + t696) * MDP(28) + (-t602 + t695) * MDP(29) + (-t354 * t525 + t692) * MDP(31) + (-t353 * t525 + t704) * MDP(32) + (MDP(28) * t675 + MDP(29) * t581 - MDP(31) * t354 - MDP(32) * t657) * qJD(6) + t705;];
tauc  = t1;
