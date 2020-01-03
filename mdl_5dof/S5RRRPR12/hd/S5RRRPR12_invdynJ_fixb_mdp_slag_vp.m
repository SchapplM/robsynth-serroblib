% Calculate vector of inverse dynamics joint torques for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:57
% EndTime: 2019-12-31 21:41:20
% DurationCPUTime: 14.08s
% Computational Cost: add. (7301->631), mult. (18308->880), div. (0->0), fcn. (14549->14), ass. (0->252)
t618 = cos(qJ(3));
t734 = cos(pkin(5));
t670 = t734 * qJD(1);
t645 = t670 + qJD(2);
t614 = sin(qJ(3));
t615 = sin(qJ(2));
t611 = sin(pkin(5));
t701 = qJD(1) * t611;
t682 = t615 * t701;
t662 = t614 * t682;
t527 = -t618 * t645 + t662;
t520 = qJD(5) + t527;
t529 = t614 * t645 + t618 * t682;
t619 = cos(qJ(2));
t700 = qJD(1) * t619;
t681 = t611 * t700;
t581 = -qJD(3) + t681;
t610 = sin(pkin(10));
t612 = cos(pkin(10));
t485 = t529 * t610 + t612 * t581;
t617 = cos(qJ(5));
t487 = t529 * t612 - t581 * t610;
t613 = sin(qJ(5));
t732 = t487 * t613;
t751 = -t617 * t485 - t732;
t753 = t751 * t520;
t661 = pkin(1) * t670;
t547 = pkin(7) * t681 + t615 * t661;
t752 = qJD(4) * t614 + t547 + t581 * (pkin(3) * t614 - qJ(4) * t618);
t646 = t485 * t613 - t487 * t617;
t750 = t520 * t646;
t690 = t615 * qJDD(1);
t674 = t611 * t690;
t692 = qJD(1) * qJD(2);
t675 = t611 * t692;
t748 = -t619 * t675 - t674;
t616 = sin(qJ(1));
t741 = cos(qJ(1));
t657 = t734 * t741;
t556 = t615 * t657 + t616 * t619;
t683 = t611 * t741;
t503 = t556 * t618 - t614 * t683;
t555 = t615 * t616 - t619 * t657;
t607 = pkin(10) + qJ(5);
t604 = sin(t607);
t605 = cos(t607);
t747 = t503 * t604 - t555 * t605;
t746 = t503 * t605 + t555 * t604;
t666 = t734 * qJDD(1);
t595 = t666 + qJDD(2);
t630 = qJD(3) * t645;
t698 = qJD(2) * t619;
t677 = t614 * t698;
t696 = qJD(3) * t618;
t451 = -t618 * t595 + t611 * (qJD(1) * (t615 * t696 + t677) + t614 * t690) + t614 * t630;
t544 = -pkin(7) * t682 + t619 * t661;
t655 = pkin(2) * t615 - pkin(8) * t619;
t545 = t655 * t701;
t707 = t618 * t544 + t614 * t545;
t462 = qJ(4) * t682 + t707;
t745 = t612 * t462 + t752 * t610;
t697 = qJD(3) * t614;
t688 = pkin(8) * t697;
t710 = -t752 * t612 + (t462 + t688) * t610;
t719 = t611 * t615;
t596 = pkin(7) * t719;
t671 = t619 * t734;
t744 = pkin(1) * t671 - t596;
t691 = qJDD(1) * t619;
t594 = t611 * t691;
t660 = t615 * t675;
t541 = qJDD(3) - t594 + t660;
t743 = -pkin(3) * t541 + qJDD(4);
t620 = qJD(1) ^ 2;
t739 = pkin(8) * t615;
t557 = t615 * t741 + t616 * t671;
t738 = g(1) * t557;
t737 = g(2) * t555;
t736 = pkin(9) + qJ(4);
t735 = pkin(8) * qJD(3);
t733 = qJ(4) * t451;
t731 = t527 * t581;
t730 = t527 * t610;
t729 = t529 * t581;
t726 = t581 * t614;
t725 = t595 * MDP(8);
t724 = t604 * t618;
t723 = t605 * t618;
t606 = t611 ^ 2;
t722 = t606 * t620;
t721 = t610 * t613;
t720 = t610 * t618;
t718 = t611 * t616;
t717 = t611 * t618;
t716 = t611 * t619;
t715 = t612 * t614;
t714 = t612 * t618;
t713 = t618 * t619;
t642 = qJD(2) * t661;
t658 = pkin(1) * t666;
t685 = -pkin(7) * t594 - t615 * t658 - t619 * t642;
t625 = -pkin(7) * t660 - t685;
t471 = pkin(8) * t595 + t625;
t639 = t655 * qJD(2);
t643 = -pkin(2) * t619 - pkin(1) - t739;
t474 = (qJD(1) * t639 + qJDD(1) * t643) * t611;
t513 = pkin(8) * t645 + t547;
t540 = t643 * t611;
t519 = qJD(1) * t540;
t635 = t618 * t471 + t614 * t474 - t513 * t697 + t519 * t696;
t396 = qJ(4) * t541 - qJD(4) * t581 + t635;
t450 = -qJD(3) * t662 + t614 * t595 + (t630 - t748) * t618;
t664 = pkin(7) * t748 - t615 * t642 + t619 * t658;
t472 = -pkin(2) * t595 - t664;
t398 = pkin(3) * t451 - qJ(4) * t450 - qJD(4) * t529 + t472;
t384 = t612 * t396 + t610 * t398;
t672 = t615 * t734;
t703 = pkin(1) * t672 + pkin(7) * t716;
t539 = pkin(8) * t734 + t703;
t546 = t611 * t639;
t548 = t744 * qJD(2);
t634 = -t539 * t697 + t540 * t696 + t614 * t546 + t618 * t548;
t699 = qJD(2) * t615;
t417 = (qJ(4) * t699 - qJD(4) * t619) * t611 + t634;
t554 = t614 * t734 + t615 * t717;
t499 = qJD(3) * t554 + t611 * t677;
t553 = t614 * t719 - t618 * t734;
t678 = t611 * t698;
t500 = -qJD(3) * t553 + t618 * t678;
t549 = t703 * qJD(2);
t423 = t499 * pkin(3) - t500 * qJ(4) - t554 * qJD(4) + t549;
t392 = t612 * t417 + t610 * t423;
t512 = -pkin(2) * t645 - t544;
t437 = t527 * pkin(3) - t529 * qJ(4) + t512;
t448 = t618 * t513 + t614 * t519;
t440 = -qJ(4) * t581 + t448;
t405 = t610 * t437 + t612 * t440;
t447 = -t614 * t513 + t519 * t618;
t473 = pkin(3) * t529 + qJ(4) * t527;
t414 = t612 * t447 + t610 * t473;
t538 = -pkin(2) * t734 - t744;
t459 = t553 * pkin(3) - t554 * qJ(4) + t538;
t708 = t618 * t539 + t614 * t540;
t460 = -qJ(4) * t716 + t708;
t412 = t610 * t459 + t612 * t460;
t514 = -t612 * t682 + t681 * t720;
t515 = (t610 * t615 + t612 * t713) * t701;
t563 = -t617 * t612 + t721;
t564 = t610 * t617 + t612 * t613;
t695 = qJD(5) * t614;
t712 = t514 * t613 - t515 * t617 - t563 * t696 - t564 * t695;
t694 = qJD(5) * t617;
t711 = -t617 * t514 - t515 * t613 + t564 * t696 + t694 * t715 - t695 * t721;
t709 = -t612 * t688 - t745;
t706 = t520 * t563;
t705 = t520 * t564;
t530 = t614 * t544;
t463 = -pkin(3) * t682 - t545 * t618 + t530;
t684 = pkin(4) * t610 + pkin(8);
t704 = -pkin(4) * t514 + t684 * t696 - t463;
t641 = pkin(3) * t618 + qJ(4) * t614 + pkin(2);
t526 = pkin(8) * t714 - t610 * t641;
t608 = t615 ^ 2;
t702 = -t619 ^ 2 + t608;
t439 = pkin(3) * t581 + qJD(4) - t447;
t693 = -qJD(4) + t439;
t689 = 0.2e1 * t606;
t687 = t619 * t722;
t431 = t450 * t610 - t612 * t541;
t432 = t450 * t612 + t541 * t610;
t686 = -t613 * t431 + t617 * t432 - t485 * t694;
t680 = t614 * t700;
t679 = t611 * t699;
t676 = t619 * t692;
t383 = -t396 * t610 + t612 * t398;
t379 = pkin(4) * t451 - pkin(9) * t432 + t383;
t380 = -pkin(9) * t431 + t384;
t669 = t617 * t379 - t613 * t380;
t391 = -t417 * t610 + t612 * t423;
t668 = t617 * t431 + t613 * t432;
t404 = t612 * t437 - t440 * t610;
t411 = t612 * t459 - t460 * t610;
t413 = -t447 * t610 + t612 * t473;
t667 = -t614 * t539 + t540 * t618;
t665 = t614 * t471 - t618 * t474 + t513 * t696 + t519 * t697;
t656 = t611 * t620 * t734;
t502 = t556 * t614 + t618 * t683;
t558 = -t616 * t672 + t619 * t741;
t506 = t558 * t614 - t616 * t717;
t654 = -g(1) * t502 + g(2) * t506;
t653 = g(1) * t558 + g(2) * t556;
t461 = pkin(3) * t716 - t667;
t501 = -pkin(9) * t610 * t614 + t526;
t652 = pkin(4) * t611 * t680 - pkin(9) * t515 + qJD(5) * t501 - (pkin(4) * t614 - pkin(9) * t714) * qJD(3) - t710;
t561 = t612 * t641;
t490 = -pkin(9) * t715 - t561 + (-pkin(8) * t610 - pkin(4)) * t618;
t651 = -pkin(9) * t514 - qJD(5) * t490 - (-pkin(8) * t715 - pkin(9) * t720) * qJD(3) + t745;
t649 = t613 * t379 + t617 * t380;
t390 = pkin(4) * t527 - pkin(9) * t487 + t404;
t397 = -pkin(9) * t485 + t405;
t381 = t390 * t617 - t397 * t613;
t382 = t390 * t613 + t397 * t617;
t498 = t554 * t612 - t610 * t716;
t400 = pkin(4) * t553 - pkin(9) * t498 + t411;
t497 = t554 * t610 + t612 * t716;
t406 = -pkin(9) * t497 + t412;
t648 = t400 * t617 - t406 * t613;
t647 = t400 * t613 + t406 * t617;
t441 = t617 * t497 + t498 * t613;
t442 = -t497 * t613 + t498 * t617;
t644 = 0.2e1 * t670 + qJD(2);
t640 = -t539 * t696 - t540 * t697 + t546 * t618 - t614 * t548;
t577 = t736 * t612;
t637 = pkin(9) * t527 * t612 + pkin(4) * t529 + qJD(4) * t610 + qJD(5) * t577 + t413;
t576 = t736 * t610;
t636 = pkin(9) * t730 - qJD(4) * t612 + qJD(5) * t576 + t414;
t387 = -qJD(5) * t732 + t686;
t633 = t611 * (t666 + t595);
t632 = g(1) * t506 + g(2) * t502 + g(3) * t553;
t507 = t558 * t618 + t614 * t718;
t631 = -g(1) * t507 - g(2) * t503 - g(3) * t554;
t399 = t665 + t743;
t628 = g(3) * t716 - t737 - t738;
t627 = -g(3) * t719 - t653;
t626 = -t399 + t632;
t624 = -t472 - t628;
t623 = -pkin(8) * t541 - t512 * t581;
t422 = -pkin(3) * t679 - t640;
t621 = t632 - t665;
t388 = -qJD(5) * t646 + t668;
t603 = -pkin(4) * t612 - pkin(3);
t565 = t684 * t614;
t543 = t563 * t614;
t542 = t564 * t614;
t525 = -pkin(8) * t720 - t561;
t479 = t500 * t612 + t610 * t679;
t478 = t500 * t610 - t612 * t679;
t457 = t507 * t605 + t557 * t604;
t456 = -t507 * t604 + t557 * t605;
t444 = qJDD(5) + t451;
t433 = -pkin(4) * t730 + t448;
t430 = pkin(4) * t497 + t461;
t416 = pkin(4) * t485 + t439;
t407 = pkin(4) * t478 + t422;
t402 = qJD(5) * t442 + t617 * t478 + t479 * t613;
t401 = -qJD(5) * t441 - t478 * t613 + t479 * t617;
t389 = -pkin(9) * t478 + t392;
t386 = pkin(4) * t431 + t399;
t385 = pkin(4) * t499 - pkin(9) * t479 + t391;
t377 = -t382 * qJD(5) + t669;
t376 = t381 * qJD(5) + t649;
t1 = [(-t388 * t553 - t402 * t520 - t441 * t444 + t499 * t751) * MDP(25) + (-t387 * t441 - t388 * t442 + t401 * t751 + t402 * t646) * MDP(23) + ((-qJD(5) * t647 + t385 * t617 - t389 * t613) * t520 + t648 * t444 + t377 * t553 + t381 * t499 - t407 * t751 + t430 * t388 + t386 * t441 + t416 * t402 + g(1) * t746 - g(2) * t457) * MDP(27) + (t444 * t553 + t499 * t520) * MDP(26) + (-t450 * t553 - t451 * t554 - t499 * t529 - t500 * t527) * MDP(12) + (t450 * t554 + t500 * t529) * MDP(11) + (-t640 * t581 + t667 * t541 + t549 * t527 + t538 * t451 + t472 * t553 + t512 * t499 + g(1) * t503 - g(2) * t507 + (t447 * t699 + t619 * t665) * t611) * MDP(16) + (t384 * t412 + t405 * t392 + t383 * t411 + t404 * t391 + t399 * t461 + t439 * t422 - g(1) * (-t616 * pkin(1) - t556 * pkin(2) - pkin(3) * t503 + pkin(7) * t683 - t555 * pkin(8) - qJ(4) * t502) - g(2) * (pkin(1) * t741 + t558 * pkin(2) + t507 * pkin(3) + pkin(7) * t718 + t557 * pkin(8) + t506 * qJ(4))) * MDP(21) + (t387 * t553 + t401 * t520 + t442 * t444 - t499 * t646) * MDP(24) + (t387 * t442 - t401 * t646) * MDP(22) + (-t392 * t527 - t412 * t451 - t384 * t553 - t405 * t499 + t422 * t487 + t461 * t432 + t399 * t498 + t439 * t479 - g(1) * (t503 * t610 - t555 * t612) - g(2) * (-t507 * t610 + t557 * t612)) * MDP(19) + (t391 * t527 + t411 * t451 + t383 * t553 + t404 * t499 + t422 * t485 + t461 * t431 + t399 * t497 + t439 * t478 - g(1) * (-t503 * t612 - t555 * t610) - g(2) * (t507 * t612 + t557 * t610)) * MDP(18) + (-t407 * t646 + t430 * t387 + t386 * t442 + t416 * t401 - (qJD(5) * t648 + t385 * t613 + t389 * t617) * t520 - t647 * t444 - t376 * t553 - t382 * t499 - g(1) * t747 - g(2) * t456) * MDP(28) + (-t383 * t498 - t384 * t497 - t391 * t487 - t392 * t485 - t404 * t479 - t405 * t478 - t411 * t432 - t412 * t431 - t654) * MDP(20) + (t499 * t581 - t541 * t553 + (t451 * t619 - t527 * t699) * t611) * MDP(14) + (-t500 * t581 + t541 * t554 + (-t450 * t619 + t529 * t699) * t611) * MDP(13) + (g(1) * t616 - g(2) * t741) * MDP(2) + (g(1) * t741 + g(2) * t616) * MDP(3) + t734 * t725 + (-t548 * t645 - t703 * t595 - t625 * t734 - g(1) * t555 + g(2) * t557 + (-t676 - t690) * pkin(1) * t689) * MDP(10) + (-t549 * t645 - t596 * t595 + t664 * t734 + g(1) * t556 - g(2) * t558 + (t595 * t671 + (-t615 * t692 + t691) * t689) * pkin(1)) * MDP(9) + (t634 * t581 - t708 * t541 + t549 * t529 + t538 * t450 + t472 * t554 + t512 * t500 + (-t448 * t699 + t619 * t635) * t611 + t654) * MDP(17) + (t619 * t633 - t644 * t679) * MDP(7) + (t615 * t633 + t644 * t678) * MDP(6) + ((qJDD(1) * t608 + 0.2e1 * t615 * t676) * MDP(4) + 0.2e1 * (t619 * t690 - t692 * t702) * MDP(5)) * t606 + qJDD(1) * MDP(1) + (-t541 * t619 - t581 * t699) * t611 * MDP(15); (-t387 * t542 + t388 * t543 + t646 * t711 + t712 * t751) * MDP(23) + (t388 * t618 - t444 * t542 - t520 * t711 - t726 * t751) * MDP(25) + ((t490 * t617 - t501 * t613) * t444 - t377 * t618 + t381 * t697 + t565 * t388 + t386 * t542 - g(1) * (-t557 * t723 + t558 * t604) - g(2) * (-t555 * t723 + t556 * t604) + (t613 * t651 - t617 * t652) * t520 - t704 * t751 + t711 * t416 + (-t381 * t680 - g(3) * (t604 * t615 + t605 * t713)) * t611) * MDP(27) + t725 + t702 * MDP(5) * t722 + (-t444 * t618 - t520 * t726) * MDP(26) + (-t387 * t618 - t444 * t543 + t520 * t712 + t646 * t726) * MDP(24) + (t565 * t387 - t386 * t543 - (t490 * t613 + t501 * t617) * t444 + t376 * t618 - t382 * t697 - g(1) * (t557 * t724 + t558 * t605) - g(2) * (t555 * t724 + t556 * t605) + (t613 * t652 + t617 * t651) * t520 - t704 * t646 + t712 * t416 + (t382 * t680 - g(3) * (-t604 * t713 + t605 * t615)) * t611) * MDP(28) + (-t387 * t543 - t646 * t712) * MDP(22) + (t450 * t614 - t618 * t729) * MDP(11) + (t581 * t697 + t541 * t618 + (t527 * t615 - t619 * t726) * t701) * MDP(14) + (t615 * t656 + t594) * MDP(7) + (-t581 * t696 + t541 * t614 + (-t529 * t615 + t581 * t713) * t701) * MDP(13) + (-t619 * t656 + t674) * MDP(6) + (pkin(1) * t615 * t722 + t547 * t645 - t628 + t664) * MDP(9) + (-pkin(2) * t450 - t707 * t581 + t448 * t682 - t547 * t529 + t623 * t618 + (-t581 * t735 - t624) * t614) * MDP(17) + (-t447 * t682 - pkin(2) * t451 - t547 * t527 - t530 * t581 + t623 * t614 + ((t545 + t735) * t581 + t624) * t618) * MDP(16) + (-t439 * t515 - t526 * t451 - t463 * t487 - t709 * t527 + t627 * t612 + (t384 + (pkin(8) * t487 + t439 * t612) * qJD(3) + t628 * t610) * t618 + (pkin(8) * t432 + t399 * t612 + t405 * t581) * t614) * MDP(19) + (-t439 * t514 + t525 * t451 - t463 * t485 + t710 * t527 + t627 * t610 + (-t383 + (pkin(8) * t485 + t439 * t610) * qJD(3) - t628 * t612) * t618 + (pkin(8) * t431 + t399 * t610 - t404 * t581) * t614) * MDP(18) + (t404 * t515 + t405 * t514 - t431 * t526 - t432 * t525 - t710 * t487 - t709 * t485 + (-t404 * t612 - t405 * t610) * t696 + (-t383 * t612 - t384 * t610 - t628) * t614) * MDP(20) + (pkin(1) * t687 + t544 * t645 + (pkin(7) * t692 + g(3)) * t719 + t653 + t685) * MDP(10) + (t383 * t525 + t384 * t526 - t439 * t463 + t709 * t405 + t710 * t404 + t641 * t738 + t641 * t737 + (t399 * t614 + t439 * t696 - t653) * pkin(8) - g(3) * (t619 * t641 + t739) * t611) * MDP(21) + ((t450 + t731) * t618 + (-t451 + t729) * t614) * MDP(12) + t581 * MDP(15) * t682 - t615 * MDP(4) * t687; -t527 ^ 2 * MDP(12) + (t450 - t731) * MDP(13) + (-t451 - t729) * MDP(14) + t541 * MDP(15) + (-t448 * t581 + t621) * MDP(16) + (-t447 * t581 + t512 * t527 - t631 - t635) * MDP(17) + (-t610 * t733 - pkin(3) * t431 - t448 * t485 + (t610 * t693 - t413) * t527 + t626 * t612) * MDP(18) + (-t612 * t733 - pkin(3) * t432 - t448 * t487 + (t612 * t693 + t414) * t527 - t626 * t610) * MDP(19) + (t413 * t487 + t414 * t485 + (-qJ(4) * t431 - qJD(4) * t485 - t404 * t527 + t384) * t612 + (qJ(4) * t432 + qJD(4) * t487 - t405 * t527 - t383) * t610 + t631) * MDP(20) + (-t404 * t413 - t405 * t414 - t439 * t448 + (-t404 * t610 + t405 * t612) * qJD(4) + t626 * pkin(3) + (-t383 * t610 + t384 * t612 + t631) * qJ(4)) * MDP(21) + (t387 * t564 + t646 * t706) * MDP(22) + (-t387 * t563 - t388 * t564 + t646 * t705 - t706 * t751) * MDP(23) + (t444 * t564 - t520 * t706) * MDP(24) + (-t444 * t563 - t520 * t705) * MDP(25) + ((-t576 * t617 - t577 * t613) * t444 + t603 * t388 + t386 * t563 + t433 * t751 + (t613 * t636 - t617 * t637) * t520 + t705 * t416 + t632 * t605) * MDP(27) + (t603 * t387 + t386 * t564 - (-t576 * t613 + t577 * t617) * t444 + t433 * t646 + (t613 * t637 + t617 * t636) * t520 - t706 * t416 - t632 * t604) * MDP(28) + (MDP(11) * t527 + t529 * MDP(12) - t512 * MDP(16) - t404 * MDP(18) + t405 * MDP(19) + MDP(24) * t646 - MDP(25) * t751 - t520 * MDP(26) - t381 * MDP(27) + t382 * MDP(28)) * t529; (t487 * t527 + t431) * MDP(18) + (-t485 * t527 + t432) * MDP(19) + (-t485 ^ 2 - t487 ^ 2) * MDP(20) + (t404 * t487 + t405 * t485 - t621 + t743) * MDP(21) + (t388 - t750) * MDP(27) + (t387 + t753) * MDP(28); t646 * t751 * MDP(22) + (t646 ^ 2 - t751 ^ 2) * MDP(23) + (t686 - t753) * MDP(24) + (-t668 - t750) * MDP(25) + t444 * MDP(26) + (t382 * t520 + t416 * t646 - g(1) * t456 + g(2) * t747 - g(3) * (-t554 * t604 - t605 * t716) + t669) * MDP(27) + (-t416 * t751 + t381 * t520 + g(1) * t457 + g(2) * t746 - g(3) * (-t554 * t605 + t604 * t716) - t649) * MDP(28) + (-MDP(24) * t732 + MDP(25) * t646 - MDP(27) * t382 - t381 * MDP(28)) * qJD(5);];
tau = t1;
