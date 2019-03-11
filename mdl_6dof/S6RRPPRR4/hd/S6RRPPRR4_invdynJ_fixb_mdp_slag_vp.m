% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:21
% EndTime: 2019-03-09 09:05:35
% DurationCPUTime: 11.10s
% Computational Cost: add. (7143->652), mult. (20271->845), div. (0->0), fcn. (16702->12), ass. (0->294)
t602 = sin(pkin(11));
t603 = sin(pkin(6));
t608 = sin(qJ(2));
t711 = qJD(1) * t608;
t684 = t603 * t711;
t604 = cos(pkin(11));
t612 = cos(qJ(2));
t721 = t612 * t604;
t689 = t603 * t721;
t545 = -qJD(1) * t689 + t602 * t684;
t605 = cos(pkin(6));
t712 = qJD(1) * t605;
t586 = qJD(2) + t712;
t607 = sin(qJ(5));
t611 = cos(qJ(5));
t514 = -t611 * t545 + t586 * t607;
t643 = t602 * t612 + t604 * t608;
t552 = t643 * t603;
t549 = qJD(1) * t552;
t785 = qJD(5) + t549;
t786 = t514 * t785;
t513 = qJD(6) + t514;
t610 = cos(qJ(6));
t704 = qJD(6) * t610;
t516 = t545 * t607 + t586 * t611;
t696 = qJDD(1) * t612;
t678 = t603 * t696;
t564 = t604 * t678;
t633 = t643 * qJD(2);
t697 = qJDD(1) * t608;
t499 = -t564 + (qJD(1) * t633 + t602 * t697) * t603;
t698 = qJDD(1) * t605;
t585 = qJDD(2) + t698;
t670 = -t611 * t499 + t585 * t607;
t439 = qJD(5) * t516 + t670;
t436 = qJDD(6) + t439;
t606 = sin(qJ(6));
t729 = t606 * t436;
t640 = -t513 * t704 - t729;
t744 = t516 * t606;
t469 = -t610 * t785 + t744;
t782 = t469 * t785;
t784 = t640 + t782;
t596 = pkin(2) * t612 + pkin(1);
t657 = t596 * t603;
t783 = qJDD(1) * t657;
t667 = t611 * t785;
t730 = t605 * t612;
t590 = pkin(1) * t730;
t757 = pkin(8) + qJ(3);
t680 = t757 * t608;
t659 = t603 * t680;
t527 = pkin(2) * t605 + t590 - t659;
t731 = t605 * t608;
t589 = pkin(1) * t731;
t734 = t603 * t612;
t714 = pkin(8) * t734 + t589;
t540 = qJ(3) * t734 + t714;
t480 = t527 * t604 - t602 * t540;
t768 = pkin(3) + pkin(9);
t447 = pkin(4) * t552 - t605 * t768 - t480;
t736 = t603 * t608;
t551 = t602 * t736 - t689;
t754 = qJ(4) * t552;
t628 = -t657 - t754;
t463 = t551 * t768 + t628;
t781 = t607 * t447 + t611 * t463;
t742 = t549 * t607;
t487 = -t545 * t606 + t610 * t742;
t707 = qJD(5) * t607;
t780 = t610 * t707 + t487;
t613 = cos(qJ(1));
t733 = t603 * t613;
t778 = g(2) * t733 - g(3) * t605;
t529 = (t734 * t757 + t589) * qJD(1);
t521 = t602 * t529;
t583 = qJD(1) * t590;
t528 = -qJD(1) * t659 + t583;
t477 = t528 * t604 - t521;
t701 = qJD(4) - t477;
t777 = -MDP(12) - MDP(16);
t517 = pkin(2) * t586 + t528;
t732 = t604 * t529;
t466 = t602 * t517 + t732;
t455 = -t586 * qJ(4) - t466;
t763 = pkin(4) * t545;
t434 = -t455 - t763;
t710 = qJD(2) * t608;
t683 = t603 * t710;
t658 = qJD(1) * t683;
t563 = t602 * t658;
t699 = qJD(1) * qJD(2);
t679 = t612 * t699;
t500 = -t563 + (qJDD(1) * t643 + t604 * t679) * t603;
t498 = qJDD(5) + t500;
t595 = -pkin(2) * t604 - pkin(3);
t592 = -pkin(9) + t595;
t776 = -t434 * t785 - t592 * t498;
t465 = t517 * t604 - t521;
t654 = qJD(4) - t465;
t762 = pkin(4) * t549;
t429 = -t586 * t768 + t654 + t762;
t554 = -qJD(1) * t657 + qJD(3);
t619 = -qJ(4) * t549 + t554;
t446 = t545 * t768 + t619;
t411 = t429 * t607 + t446 * t611;
t692 = pkin(1) * t696;
t580 = t605 * t692;
t693 = pkin(1) * qJD(2) * t605;
t661 = qJD(1) * t693;
t677 = qJD(2) * t757;
t709 = qJD(3) * t608;
t462 = -t608 * t661 + pkin(2) * t585 + t580 + (-qJDD(1) * t680 + (-t612 * t677 - t709) * qJD(1)) * t603;
t627 = qJD(3) * t612 - t608 * t677;
t685 = pkin(8) * t678 + qJDD(1) * t589 + t612 * t661;
t472 = (qJ(3) * t696 + qJD(1) * t627) * t603 + t685;
t421 = t462 * t604 - t602 * t472;
t653 = qJDD(4) - t421;
t408 = pkin(4) * t500 - t585 * t768 + t653;
t695 = pkin(2) * t658 + qJDD(3);
t629 = -qJ(4) * t500 - qJD(4) * t549 + t695;
t616 = t629 - t783;
t415 = t499 * t768 + t616;
t673 = -t611 * t408 + t415 * t607;
t771 = t411 * qJD(5) + t673;
t393 = -pkin(5) * t498 + t771;
t609 = sin(qJ(1));
t737 = t602 * t608;
t561 = -t721 + t737;
t636 = t561 * t605;
t507 = t609 * t636 - t613 * t643;
t735 = t603 * t609;
t488 = -t507 * t611 - t607 * t735;
t504 = -t609 * t643 - t613 * t636;
t490 = -t504 * t611 + t607 * t733;
t645 = t551 * t611 - t605 * t607;
t632 = g(1) * t488 + g(2) * t490 + g(3) * t645;
t775 = t513 * (pkin(5) * t516 + pkin(10) * t513) + t393 + t632;
t593 = pkin(2) * t602 + qJ(4);
t559 = pkin(5) * t607 - pkin(10) * t611 + t593;
t656 = pkin(5) * t611 + pkin(10) * t607;
t774 = t513 * ((-pkin(4) - t656) * t549 - qJD(5) * t656 - t701) - t559 * t436;
t715 = -t602 * t730 - t604 * t731;
t503 = t609 * t561 + t613 * t715;
t642 = t504 * t607 + t611 * t733;
t773 = -t503 * t610 + t606 * t642;
t508 = -t613 * t561 + t609 * t715;
t772 = t503 * t606 + t610 * t642;
t539 = t549 ^ 2;
t599 = t603 ^ 2;
t769 = 0.2e1 * t599;
t767 = pkin(1) * t599;
t766 = pkin(3) * t499;
t765 = pkin(3) * t551;
t764 = pkin(3) * t585;
t724 = t609 * t612;
t726 = t608 * t613;
t557 = -t605 * t724 - t726;
t761 = g(1) * t557;
t760 = g(1) * t609;
t758 = g(3) * t612;
t756 = MDP(6) * t603;
t755 = MDP(7) * t603;
t706 = qJD(5) * t611;
t438 = t607 * t499 + t545 * t706 + t611 * t585 - t586 * t707;
t686 = t610 * t438 + t606 * t498 + t704 * t785;
t705 = qJD(6) * t606;
t403 = -t516 * t705 + t686;
t753 = t403 * t606;
t471 = t516 * t610 + t606 * t785;
t672 = t438 * t606 - t610 * t498;
t404 = qJD(6) * t471 + t672;
t752 = t404 * t607;
t751 = t469 * t513;
t750 = t471 * t513;
t747 = t513 * t606;
t665 = t513 * t610;
t746 = t514 * t545;
t745 = t516 * t545;
t743 = t545 * t586;
t740 = t585 * MDP(8);
t614 = qJD(1) ^ 2;
t738 = t599 * t614;
t728 = t607 * t498;
t727 = t608 * t609;
t723 = t610 * t436;
t722 = t611 * t403;
t720 = t612 * t613;
t476 = t528 * t602 + t732;
t448 = t476 - t763;
t675 = pkin(2) * t684 + qJ(4) * t545;
t453 = t549 * t768 + t675;
t717 = t607 * t448 + t611 * t453;
t422 = t602 * t462 + t604 * t472;
t584 = t612 * t693;
t518 = t603 * t627 + t584;
t681 = t757 * t603;
t519 = -t603 * t709 + (-t612 * t681 - t589) * qJD(2);
t460 = t604 * t518 + t602 * t519;
t481 = t602 * t527 + t604 * t540;
t600 = t608 ^ 2;
t713 = -t612 ^ 2 + t600;
t708 = qJD(5) * t592;
t703 = qJD(2) - t586;
t702 = t762 + t701;
t694 = t608 * t767;
t691 = pkin(8) * t697;
t690 = t612 * t738;
t688 = t605 * t720;
t687 = t403 * t607 + (t549 * t611 + t706) * t471;
t452 = -t605 * qJD(4) - t460;
t474 = -t605 * qJ(4) - t481;
t635 = t607 * t408 + t611 * t415 + t429 * t706 - t446 * t707;
t392 = pkin(10) * t498 + t635;
t660 = t585 * qJ(4) + t586 * qJD(4) + t422;
t409 = -pkin(4) * t499 + t660;
t395 = pkin(5) * t439 - pkin(10) * t438 + t409;
t674 = -t606 * t392 + t610 * t395;
t486 = t610 * t545 + t606 * t742;
t671 = t486 * t513 + t707 * t747;
t459 = t518 * t602 - t604 * t519;
t553 = pkin(2) * t731 - t681;
t669 = -t553 * t609 + t613 * t596;
t668 = t785 ^ 2;
t666 = t785 * t516;
t664 = t586 + t712;
t663 = qJD(1) * t703;
t662 = t585 + t698;
t655 = g(1) * t613 + g(2) * t609;
t652 = t610 * t392 + t606 * t395;
t401 = pkin(10) * t785 + t411;
t418 = pkin(5) * t514 - pkin(10) * t516 + t434;
t397 = t401 * t610 + t418 * t606;
t651 = t401 * t606 - t418 * t610;
t417 = pkin(10) * t552 + t781;
t451 = -pkin(4) * t551 - t474;
t525 = t551 * t607 + t605 * t611;
t424 = -pkin(5) * t645 - pkin(10) * t525 + t451;
t650 = t417 * t610 + t424 * t606;
t649 = -t417 * t606 + t424 * t610;
t410 = t429 * t611 - t446 * t607;
t547 = t603 * t633;
t548 = qJD(2) * t689 - t602 * t683;
t638 = pkin(2) * t683 - qJ(4) * t548 - qJD(4) * t552;
t435 = t547 * t768 + t638;
t437 = pkin(4) * t548 + t459;
t648 = -t435 * t607 + t437 * t611;
t647 = t447 * t611 - t463 * t607;
t484 = t525 * t610 + t552 * t606;
t483 = t525 * t606 - t552 * t610;
t644 = -t553 * t613 - t596 * t609;
t430 = -pkin(4) * t547 - t452;
t641 = t611 * t498 - t607 * t668;
t639 = t513 * t705 - t723;
t634 = t611 * t435 + t607 * t437 + t447 * t706 - t463 * t707;
t631 = -g(1) * t508 + g(2) * t503 - g(3) * t552;
t630 = -g(1) * t507 - g(2) * t504 + g(3) * t551;
t626 = -t667 * t785 - t728;
t624 = t409 + t631;
t623 = t714 * t586;
t622 = -t611 * t705 - t780;
t400 = -pkin(5) * t785 - t410;
t621 = -pkin(10) * t436 + (t400 + t410) * t513;
t620 = t459 * t549 - t603 * t655;
t618 = qJD(6) * t513 * t592 - t631;
t617 = (-pkin(10) * t545 - qJD(6) * t559 + t717) * t513 + t630;
t473 = pkin(3) * t545 + t619;
t615 = t473 * t549 - t630 + t653;
t570 = pkin(2) * t688;
t558 = -t605 * t727 + t720;
t556 = -t605 * t726 - t724;
t555 = -t688 + t727;
t489 = -t507 * t607 + t611 * t735;
t485 = t628 + t765;
t482 = pkin(3) * t549 + t675;
t479 = qJD(5) * t525 - t547 * t611;
t478 = qJD(5) * t645 + t547 * t607;
t475 = -pkin(3) * t605 - t480;
t458 = pkin(3) * t547 + t638;
t454 = -pkin(3) * t586 + t654;
t441 = t489 * t610 + t508 * t606;
t440 = -t489 * t606 + t508 * t610;
t426 = -qJD(6) * t483 + t478 * t610 + t548 * t606;
t425 = qJD(6) * t484 + t478 * t606 - t548 * t610;
t423 = t616 + t766;
t420 = t653 - t764;
t416 = -pkin(5) * t552 - t647;
t413 = pkin(5) * t545 - t448 * t611 + t453 * t607;
t405 = pkin(5) * t479 - pkin(10) * t478 + t430;
t399 = -pkin(5) * t548 + qJD(5) * t781 - t648;
t398 = pkin(10) * t548 + t634;
t391 = -qJD(6) * t397 + t674;
t390 = -qJD(6) * t651 + t652;
t1 = [(t423 * t485 + t473 * t458 - t660 * t474 + t455 * t452 + t420 * t475 + t454 * t459 - g(1) * (pkin(3) * t503 + qJ(4) * t504 + t644) - g(2) * (pkin(3) * t508 - qJ(4) * t507 + t669)) * MDP(16) + (g(1) * t503 + g(2) * t508 + t420 * t605 - t423 * t551 - t458 * t545 + t459 * t586 - t473 * t547 + t475 * t585 - t485 * t499) * MDP(14) + t655 * MDP(3) + (-g(2) * t613 + t760) * MDP(2) + (t420 * t552 + t452 * t545 + t454 * t548 + t455 * t547 + t474 * t499 + t475 * t500 - t551 * t660 + t620) * MDP(13) + (-g(1) * t504 + g(2) * t507 - t423 * t552 - t452 * t586 - t458 * t549 - t473 * t548 - t474 * t585 - t485 * t500 + t605 * t660) * MDP(15) + ((-qJD(6) * t650 - t398 * t606 + t405 * t610) * t513 + t649 * t436 - t391 * t645 - t651 * t479 + t399 * t469 + t416 * t404 + t393 * t483 + t400 * t425 - g(1) * t772 - g(2) * t441) * MDP(29) + (-(qJD(6) * t649 + t398 * t610 + t405 * t606) * t513 - t650 * t436 + t390 * t645 - t397 * t479 + t399 * t471 + t416 * t403 + t393 * t484 + t400 * t426 + g(1) * t773 - g(2) * t440) * MDP(30) + (t438 * t645 - t439 * t525 - t478 * t514 - t479 * t516) * MDP(18) + (-t403 * t645 + t426 * t513 + t436 * t484 + t471 * t479) * MDP(26) + (-t436 * t645 + t479 * t513) * MDP(28) + (t404 * t645 - t425 * t513 - t436 * t483 - t469 * t479) * MDP(27) + (t422 * t481 + t466 * t460 + t421 * t480 - t465 * t459 - g(1) * t644 - g(2) * t669 + (pkin(2) * t554 * t710 + (-t695 + t783) * t596) * t603) * MDP(12) + (qJDD(1) * t600 + 0.2e1 * t608 * t679) * t599 * MDP(4) + (-t421 * t552 - t422 * t551 - t460 * t545 - t465 * t548 - t466 * t547 - t480 * t500 - t481 * t499 + t620) * MDP(11) + (t692 * t769 + (-pkin(8) * t736 + t590) * t585 + (-t603 * t691 + t580) * t605 - g(1) * t556 - g(2) * t558 + (-t623 + (-t605 * t714 - 0.2e1 * t694) * qJD(1)) * qJD(2)) * MDP(9) + (t612 * t662 - t664 * t710) * t755 + (qJD(2) * t612 * t664 + t608 * t662) * t756 + t605 * t740 + (-(-pkin(8) * t683 + t584) * t586 - t714 * t585 - (-pkin(8) * t658 + t685) * t605 - g(1) * t555 - g(2) * t557 + 0.2e1 * (-t679 - t697) * t767) * MDP(10) + (t608 * t696 - t699 * t713) * MDP(5) * t769 + qJDD(1) * MDP(1) + (t438 * t525 + t478 * t516) * MDP(17) + (-t403 * t483 - t404 * t484 - t425 * t471 - t426 * t469) * MDP(25) + (t403 * t484 + t426 * t471) * MDP(24) + (-t439 * t552 - t479 * t785 + t498 * t645 - t514 * t548) * MDP(20) + (t648 * t785 + t647 * t498 - t673 * t552 + t410 * t548 + t430 * t514 + t451 * t439 - t409 * t645 + t434 * t479 - g(1) * t642 - g(2) * t489 + (-t411 * t552 - t781 * t785) * qJD(5)) * MDP(22) + (g(1) * t490 - g(2) * t488 + t409 * t525 - t411 * t548 + t430 * t516 + t434 * t478 + t451 * t438 - t498 * t781 - t552 * t635 - t634 * t785) * MDP(23) + (t498 * t552 + t548 * t785) * MDP(21) + (t438 * t552 + t478 * t785 + t498 * t525 + t516 * t548) * MDP(19); (t660 * t593 + t420 * t595 - t473 * t482 - t454 * t476 - g(1) * (pkin(2) * t557 + pkin(3) * t507 + qJ(4) * t508) - g(2) * (-pkin(2) * t727 + pkin(3) * t504 - qJ(4) * t503 + t570) - g(3) * (pkin(2) * t734 + t754 - t765) - t701 * t455) * MDP(16) + (-g(2) * t570 + t465 * t476 - t466 * t477 + (t422 * t602 + t421 * t604 - t761 + g(2) * t727 + (-t554 * t711 - t758) * t603) * pkin(2)) * MDP(12) + (t469 * t487 + t471 * t486 + (t469 * t610 + t471 * t606) * t707 + (-t753 - t404 * t610 + (t469 * t606 - t471 * t610) * qJD(6)) * t611) * MDP(25) + (t641 + t745) * MDP(19) + (t626 - t746) * MDP(20) + (pkin(1) * t690 + g(1) * t558 - g(2) * t556 + t583 * t586 + (pkin(8) * t663 + g(3)) * t736 - t685) * MDP(10) + (-t400 * t486 - t413 * t469 - t774 * t610 + t617 * t606 + (-t592 * t729 + t391 + (-t400 * t606 + t469 * t592) * qJD(5) - t618 * t610) * t607 + (t400 * t704 + t393 * t606 - t651 * t549 - t592 * t404 + (-t592 * t747 - t651) * qJD(5)) * t611) * MDP(29) - t608 * MDP(4) * t690 + (t438 * t611 - t607 * t666) * MDP(17) + (-t400 * t487 - t413 * t471 + t774 * t606 + t617 * t610 + (-t592 * t723 - t390 + (-t400 * t610 + t471 * t592) * qJD(5) + t618 * t606) * t607 + (-t400 * t705 + t393 * t610 - t397 * t549 - t592 * t403 + (-t592 * t665 - t397) * qJD(5)) * t611) * MDP(30) + (-t752 + (t640 - t782) * t611 + t671) * MDP(27) + (t436 * t607 + t513 * t667) * MDP(28) + (t513 * t622 + t611 * t723 + t687) * MDP(26) + (t471 * t622 + t610 * t722) * MDP(24) + (-t499 * t593 + t500 * t595 + (-t455 - t476) * t549 + (t454 - t701) * t545) * MDP(13) + (-t473 * t545 + t482 * t549 + t585 * t593 + t586 * t701 + t631 + t660) * MDP(15) + ((t466 - t476) * t549 + (-t465 + t477) * t545 + (-t499 * t602 - t500 * t604) * pkin(2)) * MDP(11) + (-t476 * t586 + t482 * t545 + (-pkin(3) + t595) * t585 + t615) * MDP(14) + (t614 * t694 - t761 + g(2) * t555 + t580 + (-t691 - t758) * t603 + (-qJD(2) * t714 + t623) * qJD(1)) * MDP(9) + (-t703 * t711 + t696) * t755 + (t612 * t663 + t697) * t756 + t740 + ((-t439 - t666) * t611 + (-t438 + t786) * t607) * MDP(18) + t785 * t545 * MDP(21) + (t410 * t545 + t593 * t439 + t702 * t514 + (-t448 * t785 - t776) * t611 + ((t453 - t708) * t785 + t624) * t607) * MDP(22) + (t593 * t438 + t717 * t785 - t411 * t545 + t702 * t516 + t776 * t607 + (-t708 * t785 + t624) * t611) * MDP(23) + t713 * MDP(5) * t738; (t465 * t549 + t466 * t545 + t695 + t778) * MDP(12) + (-t549 * t586 + t564) * MDP(14) + (t563 + t743) * MDP(15) + (-t454 * t549 - t455 * t545 + t629 + t766 + t778) * MDP(16) + (-t728 + t746) * MDP(22) + (t745 + (t707 + t742) * t785) * MDP(23) + (t671 + t752) * MDP(29) + (t513 * t780 + t687) * MDP(30) + (-MDP(22) * t668 - t498 * MDP(23) + t784 * MDP(29) + t639 * MDP(30)) * t611 + (t777 * t760 + (-MDP(14) * t643 - MDP(15) * t721) * t699 + (-MDP(14) * t737 - MDP(15) * t643 + t596 * t777) * qJDD(1)) * t603 + (MDP(11) + MDP(13)) * (-t545 ^ 2 - t539); (t500 + t743) * MDP(13) + (-t545 * t549 + t585) * MDP(14) + (-t586 ^ 2 - t539) * MDP(15) + (t455 * t586 + t615 - t764) * MDP(16) + (-t514 * t586 + t641) * MDP(22) + (-t516 * t586 + t626) * MDP(23) + (-t611 * t404 + (-t610 * t586 - t606 * t667) * t513 + t784 * t607) * MDP(29) + (-t722 + (t606 * t586 - t610 * t667) * t513 + (t471 * t785 + t639) * t607) * MDP(30); -t514 ^ 2 * MDP(18) + (t438 + t786) * MDP(19) - t670 * MDP(20) + t498 * MDP(21) + (t411 * t785 - t632 - t771) * MDP(22) + (g(1) * t489 - g(2) * t642 + g(3) * t525 + t410 * t785 + t434 * t514 - t635) * MDP(23) + (t471 * t665 + t753) * MDP(24) + ((t403 - t751) * t610 + (-t404 - t750) * t606) * MDP(25) + (t513 * t665 + t729) * MDP(26) + (-t513 ^ 2 * t606 + t723) * MDP(27) + (-pkin(5) * t404 - t411 * t469 + t621 * t606 - t610 * t775) * MDP(29) + (-pkin(5) * t403 - t411 * t471 + t606 * t775 + t621 * t610) * MDP(30) + (t514 * MDP(17) + (-qJD(5) + t785) * MDP(20) - t434 * MDP(22) - t471 * MDP(26) + t469 * MDP(27) - t513 * MDP(28) + t651 * MDP(29) + t397 * MDP(30) + t516 * MDP(18)) * t516; t471 * t469 * MDP(24) + (-t469 ^ 2 + t471 ^ 2) * MDP(25) + (t686 + t751) * MDP(26) + (-t672 + t750) * MDP(27) + t436 * MDP(28) + (-g(1) * t440 - g(2) * t773 + g(3) * t483 + t397 * t513 - t400 * t471 + t674) * MDP(29) + (g(1) * t441 - g(2) * t772 + g(3) * t484 + t400 * t469 - t513 * t651 - t652) * MDP(30) + (-MDP(26) * t744 - MDP(27) * t471 - MDP(29) * t397 + MDP(30) * t651) * qJD(6);];
tau  = t1;
