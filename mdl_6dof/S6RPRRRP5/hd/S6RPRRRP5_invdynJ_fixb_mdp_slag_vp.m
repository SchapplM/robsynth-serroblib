% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:43
% EndTime: 2019-03-09 06:12:56
% DurationCPUTime: 10.38s
% Computational Cost: add. (11710->581), mult. (28344->693), div. (0->0), fcn. (22473->14), ass. (0->252)
t594 = cos(qJ(5));
t687 = qJD(5) * t594;
t589 = cos(pkin(10));
t595 = cos(qJ(3));
t703 = t595 * t589;
t588 = sin(pkin(10));
t592 = sin(qJ(3));
t709 = t588 * t592;
t530 = -t703 + t709;
t521 = t530 * qJD(1);
t531 = t588 * t595 + t589 * t592;
t522 = t531 * qJD(1);
t591 = sin(qJ(4));
t744 = cos(qJ(4));
t620 = -t744 * t521 - t591 * t522;
t776 = t620 * t594;
t792 = t687 - t776;
t678 = qJD(1) * qJD(3);
t659 = t592 * t678;
t675 = qJDD(1) * t595;
t773 = qJDD(1) * t592 + t595 * t678;
t665 = t588 * t675 + t589 * t773;
t502 = -t588 * t659 + t665;
t666 = -t588 * t773 - t589 * t659;
t613 = t589 * t675 + t666;
t653 = t591 * t502 - t744 * t613;
t619 = -t591 * t521 + t522 * t744;
t754 = t619 * qJD(4);
t448 = t653 + t754;
t446 = qJDD(5) + t448;
t439 = t594 * t446;
t681 = qJD(5) - t620;
t590 = sin(qJ(5));
t786 = t681 * t590;
t791 = -t681 * t786 + t439;
t673 = qJDD(3) + qJDD(4);
t688 = qJD(5) * t590;
t605 = t502 * t744 + t591 * t613;
t782 = t620 * qJD(4);
t600 = t605 + t782;
t674 = qJD(3) + qJD(4);
t774 = qJD(5) * t674 + t600;
t432 = -t590 * t673 - t594 * t774 + t619 * t688;
t430 = t432 * t590;
t431 = t432 * t594;
t489 = t590 * t674 + t594 * t619;
t438 = t590 * t446;
t626 = t681 * t792 + t438;
t433 = t590 * t774 - t594 * t673 + t619 * t687;
t487 = t590 * t619 - t594 * t674;
t627 = -t590 * t433 - t487 * t792;
t684 = t620 * qJD(3);
t685 = t619 * qJD(3);
t720 = t489 * t619;
t790 = (t626 - t720) * MDP(24) + t673 * MDP(19) + (t685 - t653) * MDP(18) - t620 ^ 2 * MDP(16) + (-MDP(15) * t620 + MDP(16) * t619 - MDP(26) * t681) * t619 + (-t684 + t605) * MDP(17) + (t489 * t792 - t430) * MDP(22) + (-t489 * t786 - t431 + t627) * MDP(23);
t789 = t681 ^ 2;
t734 = pkin(7) + qJ(2);
t552 = t734 * t588;
t532 = qJD(1) * t552;
t553 = t734 * t589;
t533 = qJD(1) * t553;
t631 = t532 * t592 - t533 * t595;
t486 = -pkin(8) * t521 - t631;
t483 = t744 * t486;
t760 = -t595 * t532 - t533 * t592;
t485 = -pkin(8) * t522 + t760;
t484 = qJD(3) * pkin(3) + t485;
t456 = t591 * t484 + t483;
t450 = pkin(9) * t674 + t456;
t664 = pkin(2) * t589 + pkin(1);
t540 = -qJD(1) * t664 + qJD(2);
t505 = pkin(3) * t521 + t540;
t457 = -pkin(4) * t620 - pkin(9) * t619 + t505;
t419 = t450 * t594 + t457 * t590;
t412 = qJ(6) * t681 + t419;
t788 = t412 * t681;
t482 = t591 * t486;
t455 = t484 * t744 - t482;
t449 = -pkin(4) * t674 - t455;
t426 = t487 * pkin(5) - t489 * qJ(6) + t449;
t787 = t426 * t681;
t593 = sin(qJ(1));
t596 = cos(qJ(1));
t639 = g(1) * t596 + g(2) * t593;
t778 = t426 * t620;
t777 = t449 * t620;
t587 = pkin(10) + qJ(3);
t582 = qJ(4) + t587;
t570 = sin(t582);
t712 = t570 * t596;
t713 = t570 * t593;
t775 = g(1) * t712 + g(2) * t713;
t730 = qJDD(1) * pkin(1);
t579 = qJDD(2) - t730;
t758 = g(1) * t593 - g(2) * t596;
t628 = -t579 + t758;
t637 = t594 * pkin(5) + t590 * qJ(6);
t571 = cos(t582);
t772 = t639 * t571;
t471 = pkin(4) * t619 - pkin(9) * t620;
t679 = qJD(1) * qJD(2);
t746 = qJDD(1) * t734 + t679;
t509 = t746 * t588;
t510 = t746 * t589;
t652 = -t595 * t509 - t592 * t510;
t444 = qJDD(3) * pkin(3) - pkin(8) * t502 + qJD(3) * t631 + t652;
t632 = -t592 * t509 + t595 * t510;
t447 = t613 * pkin(8) + qJD(3) * t760 + t632;
t660 = qJD(4) * t744;
t689 = qJD(4) * t591;
t607 = t591 * t444 + t447 * t744 + t484 * t660 - t486 * t689;
t564 = g(3) * t570;
t667 = t564 + t772;
t771 = -t505 * t620 - t607 + t667;
t636 = pkin(5) * t590 - qJ(6) * t594;
t770 = pkin(5) * t688 - qJ(6) * t687 - qJD(6) * t590 - t620 * t636;
t765 = pkin(5) * t619;
t762 = qJ(6) * t619;
t723 = t487 * t619;
t405 = pkin(9) * t673 + t607;
t677 = qJDD(1) * t589;
t409 = -pkin(2) * t677 - pkin(3) * t613 + t448 * pkin(4) - pkin(9) * t600 + t579;
t615 = t594 * t405 + t590 * t409 - t450 * t688 + t457 * t687;
t731 = qJ(6) * t446;
t395 = qJD(6) * t681 + t615 + t731;
t646 = t590 * t405 - t594 * t409 + t450 * t687 + t457 * t688;
t742 = pkin(5) * t446;
t397 = qJDD(6) + t646 - t742;
t761 = t395 * t594 + t397 * t590;
t458 = t591 * t485 + t483;
t642 = pkin(3) * t689 - t458;
t528 = t595 * t552;
t651 = -t553 * t592 - t528;
t493 = -pkin(8) * t531 + t651;
t693 = -t592 * t552 + t595 * t553;
t494 = -pkin(8) * t530 + t693;
t467 = t591 * t493 + t494 * t744;
t504 = -t591 * t530 + t531 * t744;
t508 = pkin(3) * t530 - t664;
t618 = -t530 * t744 - t591 * t531;
t468 = -pkin(4) * t618 - pkin(9) * t504 + t508;
t696 = t594 * t467 + t590 * t468;
t759 = t571 * pkin(4) + t570 * pkin(9);
t757 = t681 * t688 - t439;
t418 = -t450 * t590 + t457 * t594;
t680 = qJD(6) - t418;
t411 = -pkin(5) * t681 + t680;
t755 = t411 * t590 + t412 * t594;
t753 = qJ(2) * qJDD(1);
t692 = t775 * t594;
t752 = t411 * t619 + t426 * t688 + t692;
t736 = g(3) * t590;
t668 = -t571 * t736 + t590 * t775;
t645 = -t744 * t444 + t591 * t447 + t484 * t689 + t486 * t660;
t406 = -pkin(4) * t673 + t645;
t398 = t433 * pkin(5) + t432 * qJ(6) - t489 * qJD(6) + t406;
t729 = t398 * t590;
t751 = -t412 * t619 + t668 - t729;
t750 = -t418 * t619 + t449 * t688 + t692;
t749 = t406 * t590 + t419 * t619 + t449 * t687 - t668;
t737 = g(3) * t571;
t748 = -t505 * t619 - t645 - t737 + t775;
t745 = t489 ^ 2;
t524 = t531 * qJD(3);
t743 = pkin(3) * t524;
t741 = pkin(9) * t446;
t733 = pkin(9) * qJD(5);
t732 = MDP(5) * t588;
t726 = t419 * t681;
t725 = t433 * t594;
t576 = pkin(3) * t591 + pkin(9);
t724 = t446 * t576;
t722 = t487 * t590;
t721 = t489 * t487;
t719 = t489 * t590;
t718 = t489 * t594;
t715 = t504 * t594;
t707 = t590 * t593;
t705 = t593 * t594;
t704 = t594 * t596;
t702 = t596 * t590;
t701 = -t642 - t770;
t698 = t594 * t455 + t590 * t471;
t459 = t485 * t744 - t482;
t462 = pkin(3) * t522 + t471;
t697 = t594 * t459 + t590 * t462;
t694 = -t456 + t770;
t691 = t588 ^ 2 + t589 ^ 2;
t690 = qJD(3) * t522;
t672 = t744 * pkin(3);
t663 = qJD(1) * t709;
t661 = t576 * t688;
t656 = -t398 - t737;
t655 = -t406 - t737;
t654 = t691 * qJD(1) ^ 2;
t647 = pkin(3) * t660;
t644 = t571 * t637 + t759;
t643 = 0.2e1 * t691;
t514 = t571 * t707 + t704;
t516 = t571 * t702 - t705;
t641 = -g(1) * t514 + g(2) * t516;
t515 = t571 * t705 - t702;
t517 = t571 * t704 + t707;
t640 = g(1) * t515 - g(2) * t517;
t634 = t411 * t594 - t412 * t590;
t633 = -t724 - t777;
t630 = pkin(4) + t637;
t581 = cos(t587);
t568 = pkin(3) * t581;
t629 = t568 + t664 + t759;
t624 = t620 * t786 - t757;
t623 = t411 * t687 - t667 + t761;
t580 = sin(t587);
t622 = t639 * t580;
t621 = t493 * t744 - t591 * t494;
t523 = t530 * qJD(3);
t472 = qJD(4) * t618 - t523 * t744 - t591 * t524;
t617 = t472 * t590 + t504 * t687;
t616 = -t472 * t594 + t504 * t688;
t610 = -qJD(3) * t528 + qJD(2) * t703 + (-qJD(2) * t588 - qJD(3) * t553) * t592;
t474 = -pkin(8) * t524 + t610;
t602 = -t531 * qJD(2) - qJD(3) * t693;
t475 = pkin(8) * t523 + t602;
t421 = qJD(4) * t621 + t474 * t744 + t591 * t475;
t473 = qJD(4) * t504 - t591 * t523 + t524 * t744;
t435 = pkin(4) * t473 - pkin(9) * t472 + t743;
t614 = t594 * t421 + t590 * t435 - t467 * t688 + t468 * t687;
t609 = g(1) * t516 + g(2) * t514 + t570 * t736 - t646;
t608 = t643 * t679 - t639;
t606 = qJD(5) * t634 + t761;
t604 = t426 * t489 + qJDD(6) - t609;
t603 = -g(1) * t517 - g(2) * t515 - t564 * t594 + t615;
t422 = qJD(4) * t467 + t591 * t474 - t475 * t744;
t599 = t639 * t570 * t630 - pkin(9) * t772;
t584 = -pkin(8) - t734;
t577 = -t672 - pkin(4);
t539 = -qJDD(1) * t664 + qJDD(2);
t526 = -t672 - t630;
t491 = qJDD(2) - t666 * pkin(3) + (-pkin(1) + (-pkin(3) * t595 - pkin(2)) * t589) * qJDD(1);
t460 = pkin(5) * t489 + qJ(6) * t487;
t436 = t504 * t636 - t621;
t424 = pkin(5) * t618 + t467 * t590 - t468 * t594;
t423 = -qJ(6) * t618 + t696;
t417 = t455 * t590 - t471 * t594 - t765;
t416 = t698 + t762;
t415 = t459 * t590 - t462 * t594 - t765;
t414 = t697 + t762;
t413 = t487 * t681 - t432;
t401 = t636 * t472 + (qJD(5) * t637 - qJD(6) * t594) * t504 + t422;
t400 = -pkin(5) * t473 + qJD(5) * t696 + t421 * t590 - t435 * t594;
t399 = qJ(6) * t473 - qJD(6) * t618 + t614;
t1 = [(qJD(3) * t602 + qJDD(3) * t651 + t540 * t524 + t539 * t530 + t581 * t758 + t613 * t664) * MDP(13) + (-qJD(3) * t610 - qJDD(3) * t693 - t502 * t664 - t540 * t523 + t539 * t531 - t580 * t758) * MDP(14) + (-t399 * t487 + t400 * t489 - t423 * t433 - t424 * t432 + t758 * t570 + t634 * t472 + (-qJD(5) * t755 - t395 * t590 + t397 * t594) * t504) * MDP(30) + t758 * MDP(2) + (pkin(1) * t628 + (t691 * t753 + t608) * qJ(2)) * MDP(7) + (t643 * t753 + t608) * MDP(6) + (t395 * t423 + t412 * t399 + t398 * t436 + t426 * t401 + t397 * t424 + t411 * t400 - g(1) * (-pkin(5) * t515 - qJ(6) * t514) - g(2) * (pkin(5) * t517 + qJ(6) * t516) + (g(1) * t584 - g(2) * t629) * t596 + (g(1) * t629 + g(2) * t584) * t593) * MDP(32) + t639 * MDP(3) + (-qJD(3) * t523 + qJDD(3) * t531) * MDP(10) + (t502 * t531 - t522 * t523) * MDP(8) + (t472 * t619 + t504 * t600) * MDP(15) + (-g(1) * t713 + g(2) * t712 - t421 * t674 - t467 * t673 + t505 * t472 + t491 * t504 + t508 * t600 + t619 * t743) * MDP(21) + (MDP(4) * t589 - t732) * (t628 + t730) + (t646 * t618 + t418 * t473 + t422 * t487 - t621 * t433 + ((-qJD(5) * t467 + t435) * t681 + t468 * t446 + t449 * qJD(5) * t504) * t594 + ((-qJD(5) * t468 - t421) * t681 - t467 * t446 + t406 * t504 + t449 * t472) * t590 + t640) * MDP(27) + (t406 * t715 - t419 * t473 + t422 * t489 + t432 * t621 - t446 * t696 - t449 * t616 - t614 * t681 + t615 * t618 + t641) * MDP(28) + (-t446 * t618 + t473 * t681) * MDP(26) + (t432 * t618 + t439 * t504 + t473 * t489 - t616 * t681) * MDP(24) + (t433 * t618 - t438 * t504 - t473 * t487 - t617 * t681) * MDP(25) + (-t395 * t618 - t398 * t715 + t399 * t681 - t401 * t489 + t412 * t473 + t423 * t446 + t426 * t616 + t432 * t436 - t641) * MDP(31) + (t397 * t618 - t400 * t681 + t401 * t487 - t411 * t473 - t424 * t446 + t426 * t617 + t433 * t436 + t504 * t729 + t640) * MDP(29) + (-t473 * t674 + t618 * t673) * MDP(18) + qJDD(1) * MDP(1) + (-qJD(3) * t524 - qJDD(3) * t530) * MDP(11) + (-t502 * t530 + t523 * t521 - t522 * t524 + t531 * t613) * MDP(9) + (t472 * t674 + t504 * t673) * MDP(17) + (-t432 * t715 - t489 * t616) * MDP(22) + ((-t487 * t594 - t719) * t472 + (t430 - t725 + (-t718 + t722) * qJD(5)) * t504) * MDP(23) + (-t422 * t674 + t508 * t448 + t505 * t473 - t491 * t618 + t571 * t758 - t620 * t743 + t621 * t673) * MDP(20) + (-t504 * t448 + t472 * t620 - t473 * t619 + t600 * t618) * MDP(16); -MDP(4) * t677 + qJDD(1) * t732 - MDP(6) * t654 + (-qJ(2) * t654 - t628) * MDP(7) + (-t613 + t690) * MDP(13) + ((-t521 - t663) * qJD(3) + t665) * MDP(14) + (t685 + t653 + 0.2e1 * t754) * MDP(20) + (t605 + t684 + 0.2e1 * t782) * MDP(21) + (t624 - t723) * MDP(27) + (-t594 * t789 - t438 - t720) * MDP(28) + (-t723 + t791) * MDP(29) + (t681 * t719 + t431 + t627) * MDP(30) + (t626 + t720) * MDP(31) + (-t426 * t619 + (-t397 + t788) * t594 + (t411 * t681 + t395) * t590 - t758) * MDP(32); (t459 * t674 + (-t522 * t619 - t591 * t673 - t660 * t674) * pkin(3) + t771) * MDP(21) + (t398 * t526 - t412 * t414 - t411 * t415 - g(3) * (t568 + t644) - t701 * t426 + (t660 * t755 + t622) * pkin(3) + t606 * t576 + t599) * MDP(32) + (-g(3) * t581 - t540 * t522 + t622 + t652) * MDP(13) + ((t521 - t663) * qJD(3) + t665) * MDP(10) + (g(3) * t580 + t540 * t521 + t581 * t639 - t632) * MDP(14) + qJDD(3) * MDP(12) + (-t521 ^ 2 + t522 ^ 2) * MDP(9) + (t613 + t690) * MDP(11) + (t624 + t723) * MDP(25) + (t414 * t487 - t415 * t489 + (-t487 * t647 - t411 * t620 + (qJD(5) * t489 - t433) * t576) * t594 + (t489 * t647 + t412 * t620 - t432 * t576 + (t487 * t576 - t412) * qJD(5)) * t590 + t623) * MDP(30) + (t458 * t674 + (t522 * t620 + t673 * t744 - t674 * t689) * pkin(3) + t748) * MDP(20) + (t526 * t433 + t656 * t594 + (-t724 - t778) * t590 - t701 * t487 + (-t576 * t687 - t590 * t647 + t415) * t681 + t752) * MDP(29) + (t577 * t433 + t655 * t594 + t633 * t590 + t642 * t487 + ((-qJD(5) * t576 - t462) * t594 + (-t647 + t459) * t590) * t681 + t750) * MDP(27) + (-t577 * t432 + t633 * t594 + t642 * t489 + (-t594 * t647 + t661 + t697) * t681 + t749) * MDP(28) + (t526 * t432 + (-t414 - t661) * t681 + t701 * t489 + (t647 * t681 + t724 - t787) * t594 + t751) * MDP(31) + t522 * t521 * MDP(8) + t790; (t456 * t674 + t748) * MDP(20) + (t455 * t674 + t771) * MDP(21) + (t723 + t791) * MDP(25) + (-pkin(4) * t433 - t456 * t487 + (t455 * t681 - t741 - t777) * t590 + ((-t471 - t733) * t681 + t655) * t594 + t750) * MDP(27) + (pkin(4) * t432 + pkin(9) * t757 - t449 * t776 - t456 * t489 + t681 * t698 + t749) * MDP(28) + (t417 * t681 - t433 * t630 + (-t741 - t778) * t590 + t694 * t487 + (-t681 * t733 + t656) * t594 + t752) * MDP(29) + (-t411 * t776 + t416 * t487 - t417 * t489 - t412 * t786 + (-t430 - t725 + (t718 + t722) * qJD(5)) * pkin(9) + t623) * MDP(30) + (-t432 * t630 + (-pkin(9) * t688 - t416) * t681 - t694 * t489 + (t741 - t787) * t594 + t751) * MDP(31) + (pkin(9) * t606 - g(3) * t644 - t398 * t630 - t411 * t417 - t412 * t416 + t426 * t694 + t599) * MDP(32) + t790; MDP(22) * t721 + (-t487 ^ 2 + t745) * MDP(23) + t413 * MDP(24) + (t489 * t681 - t433) * MDP(25) + t446 * MDP(26) + (-t449 * t489 + t609 + t726) * MDP(27) + (t418 * t681 + t449 * t487 - t603) * MDP(28) + (-t460 * t487 - t604 + t726 + 0.2e1 * t742) * MDP(29) + (pkin(5) * t432 - qJ(6) * t433 + (t412 - t419) * t489 + (t411 - t680) * t487) * MDP(30) + (0.2e1 * t731 - t426 * t487 + t460 * t489 + (0.2e1 * qJD(6) - t418) * t681 + t603) * MDP(31) + (t395 * qJ(6) - t397 * pkin(5) - t426 * t460 - t411 * t419 - g(1) * (-pkin(5) * t516 + qJ(6) * t517) - g(2) * (-pkin(5) * t514 + qJ(6) * t515) + t636 * t564 + t680 * t412) * MDP(32); (t721 - t446) * MDP(29) + t413 * MDP(30) + (-t745 - t789) * MDP(31) + (t604 - t742 - t788) * MDP(32);];
tau  = t1;
