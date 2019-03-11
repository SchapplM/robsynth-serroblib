% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP12_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:53
% EndTime: 2019-03-09 12:54:07
% DurationCPUTime: 10.40s
% Computational Cost: add. (7335->679), mult. (14658->814), div. (0->0), fcn. (9291->10), ass. (0->297)
t632 = cos(qJ(4));
t629 = sin(qJ(4));
t732 = t629 * qJD(2);
t633 = cos(qJ(2));
t748 = qJD(1) * t633;
t543 = -t632 * t748 - t732;
t710 = t629 * t748;
t744 = qJD(2) * t632;
t544 = -t710 + t744;
t628 = sin(qJ(5));
t797 = cos(qJ(5));
t483 = -t797 * t543 + t544 * t628;
t671 = t628 * t543 + t544 * t797;
t825 = t483 * t671;
t630 = sin(qJ(2));
t749 = qJD(1) * t630;
t608 = pkin(2) * t749;
t786 = qJ(3) * t633;
t684 = pkin(8) * t630 - t786;
t514 = qJD(1) * t684 + t608;
t604 = pkin(7) * t748;
t552 = pkin(3) * t748 + t604;
t535 = t632 * t552;
t777 = t629 * t630;
t679 = pkin(4) * t633 - pkin(9) * t777;
t739 = qJD(4) * t629;
t636 = -pkin(2) - pkin(8);
t788 = pkin(9) - t636;
t824 = qJD(1) * t679 - t514 * t629 - t788 * t739 + t535;
t556 = t788 * t632;
t714 = t632 * t749;
t756 = t632 * t514 + t629 * t552;
t823 = pkin(9) * t714 + qJD(4) * t556 + t756;
t696 = -qJD(4) + t749;
t725 = qJDD(1) * t633;
t656 = qJD(2) * t696 - t725;
t728 = qJD(1) * qJD(4);
t683 = t633 * t728 - qJDD(2);
t666 = t683 * t632;
t641 = t656 * t629 - t666;
t737 = qJD(4) * t633;
t712 = t629 * t737;
t713 = t630 * t744;
t660 = t712 + t713;
t677 = qJD(2) * qJD(4) + t725;
t719 = t629 * qJDD(2) + t632 * t677;
t645 = qJD(1) * t660 - t719;
t707 = t797 * qJD(5);
t735 = qJD(5) * t628;
t431 = -t543 * t707 + t544 * t735 - t628 * t645 - t797 * t641;
t588 = qJD(4) + t749;
t571 = qJD(5) + t588;
t424 = t483 * t571 - t431;
t432 = qJD(5) * t671 + t628 * t641 - t797 * t645;
t729 = qJD(1) * qJD(2);
t705 = t633 * t729;
t726 = qJDD(1) * t630;
t665 = t705 + t726;
t542 = qJDD(4) + t665;
t537 = qJDD(5) + t542;
t799 = t671 ^ 2;
t822 = t537 * MDP(26) + (t571 * t671 - t432) * MDP(25) + MDP(22) * t825 + (-t483 ^ 2 + t799) * MDP(23) + t424 * MDP(24);
t555 = t788 * t629;
t492 = -t555 * t797 - t628 * t556;
t627 = qJ(4) + qJ(5);
t612 = cos(t627);
t631 = sin(qJ(1));
t634 = cos(qJ(1));
t685 = g(1) * t634 + g(2) * t631;
t790 = g(3) * t630;
t805 = t633 * t685 + t790;
t821 = t492 * t537 + t612 * t805;
t798 = pkin(3) + pkin(7);
t624 = qJD(2) * qJ(3);
t527 = t624 + t552;
t490 = -pkin(4) * t543 + t527;
t436 = pkin(5) * t483 - qJ(6) * t671 + t490;
t820 = t436 * t483;
t819 = t483 * t490;
t738 = qJD(4) * t632;
t810 = t797 * qJD(4) + t707;
t488 = -t628 * t738 - t629 * t810 - t632 * t735;
t690 = t797 * t749;
t502 = -t628 * t714 - t629 * t690;
t758 = t488 + t502;
t724 = t628 * t777;
t757 = -qJD(1) * t724 - t628 * t739 - t629 * t735 + (t690 + t810) * t632;
t603 = pkin(7) * t749;
t817 = qJD(3) + t603;
t448 = pkin(5) * t671 + qJ(6) * t483;
t670 = t628 * t555 - t556 * t797;
t814 = -qJD(5) * t670 + t628 * t824 + t797 * t823;
t813 = -qJD(5) * t492 + t628 * t823 - t797 * t824;
t566 = t798 * t630;
t548 = t632 * t566;
t613 = t630 * qJ(3);
t618 = t633 * pkin(2);
t752 = t618 + t613;
t559 = -pkin(1) - t752;
t538 = -pkin(8) * t633 + t559;
t702 = pkin(9) * t633 - t538;
t469 = pkin(4) * t630 + t629 * t702 + t548;
t547 = t629 * t566;
t755 = t632 * t538 + t547;
t769 = t632 * t633;
t475 = -pkin(9) * t769 + t755;
t812 = t628 * t469 + t797 * t475;
t521 = t537 * qJ(6);
t560 = t571 * qJD(6);
t811 = t521 + t560;
t599 = pkin(4) * t632 + pkin(3);
t753 = pkin(4) * t738 + t599 * t749 + t817;
t768 = t632 * t634;
t529 = -t629 * t631 + t630 * t768;
t771 = t631 * t632;
t531 = t629 * t634 + t630 * t771;
t809 = -g(1) * t529 - g(2) * t531;
t808 = -MDP(28) + MDP(31);
t523 = t537 * pkin(5);
t807 = t523 - qJDD(6);
t611 = sin(t627);
t773 = t630 * t634;
t509 = t611 * t631 - t612 * t773;
t775 = t630 * t631;
t511 = t611 * t634 + t612 * t775;
t621 = g(3) * t633;
t706 = t630 * t729;
t587 = pkin(2) * t706;
t742 = qJD(3) * t630;
t652 = qJD(2) * t684 - t742;
t703 = -pkin(1) - t613;
t663 = t633 * t636 + t703;
t462 = qJD(1) * t652 + qJDD(1) * t663 + t587;
t586 = pkin(7) * t705;
t600 = pkin(7) * t726;
t704 = qJDD(3) + t586 + t600;
t480 = pkin(3) * t665 + qJDD(2) * t636 + t704;
t477 = t632 * t480;
t731 = pkin(3) * t749 + t817;
t508 = qJD(2) * t636 + t731;
t655 = -t677 + t706;
t500 = t663 * qJD(1);
t741 = qJD(4) * t500;
t419 = t542 * pkin(4) + t477 + (pkin(9) * t683 - t741) * t632 + (-pkin(9) * t655 - qJD(4) * t508 - t462) * t629;
t721 = t632 * t462 + t629 * t480 + t508 * t738;
t422 = pkin(9) * t645 - t500 * t739 + t721;
t464 = -t500 * t629 + t632 * t508;
t454 = -pkin(9) * t544 + t464;
t447 = pkin(4) * t588 + t454;
t465 = t500 * t632 + t508 * t629;
t455 = pkin(9) * t543 + t465;
t694 = -t797 * t419 + t628 * t422 + t447 * t735 + t455 * t707;
t653 = g(1) * t509 - g(2) * t511 + t612 * t621 - t694;
t642 = t436 * t671 - t653 - t807;
t806 = -t490 * t671 + t653;
t715 = t797 * t632;
t545 = t628 * t629 - t715;
t804 = t431 * t545 + t671 * t758;
t745 = qJD(2) * t630;
t607 = pkin(2) * t745;
t497 = t607 + t652;
t743 = qJD(2) * t633;
t553 = t798 * t743;
t699 = -t497 * t629 + t632 * t553;
t437 = t679 * qJD(2) + (t632 * t702 - t547) * qJD(4) + t699;
t720 = t632 * t497 + t629 * t553 + t566 * t738;
t440 = pkin(9) * t660 - t538 * t739 + t720;
t803 = -qJD(5) * t812 + t437 * t797 - t628 * t440;
t802 = t633 * (qJD(4) + qJD(5));
t695 = t628 * t419 + t797 * t422 + t447 * t707 - t455 * t735;
t412 = t695 + t811;
t413 = t694 - t807;
t779 = t628 * t455;
t427 = t447 * t797 - t779;
t730 = qJD(6) - t427;
t425 = -t571 * pkin(5) + t730;
t716 = t797 * t455;
t428 = t628 * t447 + t716;
t426 = t571 * qJ(6) + t428;
t546 = t628 * t632 + t629 * t797;
t717 = -g(1) * t773 - g(2) * t775 + t621;
t801 = t412 * t546 + t413 * t545 - t425 * t758 + t426 * t757 + t717;
t795 = g(1) * t631;
t791 = g(2) * t634;
t614 = t629 * pkin(4);
t787 = pkin(7) * qJDD(2);
t785 = qJDD(2) * pkin(2);
t784 = t428 * t571;
t781 = t537 * t546;
t780 = t544 * t588;
t494 = t545 * t537;
t778 = t629 * t542;
t776 = t629 * t633;
t774 = t630 * t632;
t638 = qJD(1) ^ 2;
t772 = t630 * t638;
t770 = t631 * t633;
t522 = t632 * t542;
t767 = t633 * t634;
t635 = -pkin(9) - pkin(8);
t766 = t633 * t635;
t765 = t636 * t542;
t595 = qJ(3) + t614;
t764 = pkin(5) * t757 - qJ(6) * t758 + qJD(6) * t545 + t753;
t763 = -qJ(6) * t748 - t814;
t762 = pkin(5) * t748 - t813;
t759 = t571 * t488 - t494;
t430 = t454 * t797 - t779;
t754 = pkin(4) * t707 + qJD(6) - t430;
t567 = t798 * t633;
t625 = t630 ^ 2;
t626 = t633 ^ 2;
t751 = t625 - t626;
t747 = qJD(2) * t483;
t746 = qJD(2) * t543;
t740 = qJD(4) * t543;
t736 = qJD(4) * t636;
t734 = t527 * qJD(4);
t733 = t544 * qJD(2);
t622 = qJDD(2) * qJ(3);
t723 = t629 * t773;
t722 = t633 * t772;
t601 = pkin(7) * t725;
t623 = qJD(2) * qJD(3);
t718 = t601 + t622 + t623;
t589 = pkin(4) * t769;
t526 = t589 + t567;
t711 = t632 * t737;
t709 = g(3) * t752;
t701 = -qJD(2) * pkin(2) + qJD(3);
t700 = -t629 * t462 + t477;
t698 = -qJD(1) * t567 - t527;
t697 = t588 + t749;
t693 = t634 * pkin(1) + pkin(2) * t767 + t631 * pkin(7) + qJ(3) * t773;
t692 = -t600 - t717;
t691 = pkin(3) * t725 + t718;
t551 = t798 * t745;
t429 = t628 * t454 + t716;
t689 = pkin(4) * t735 - t429;
t637 = qJD(2) ^ 2;
t688 = pkin(7) * t637 + t791;
t687 = g(1) * t511 + g(2) * t509;
t510 = t611 * t773 + t612 * t631;
t512 = -t611 * t775 + t612 * t634;
t686 = -g(1) * t512 - g(2) * t510;
t682 = t500 * t630 + t538 * t588;
t554 = t603 + t701;
t565 = -t604 - t624;
t681 = t554 * t633 + t565 * t630;
t678 = t703 - t618;
t676 = -0.2e1 * pkin(1) * t729 - t787;
t674 = t469 * t797 - t628 * t475;
t528 = t678 * qJD(1);
t669 = t528 * t749 + qJDD(3) - t692;
t668 = (-pkin(5) * t612 - qJ(6) * t611) * t633;
t667 = -qJ(3) * t743 - t742;
t664 = t628 * t437 + t797 * t440 + t469 * t707 - t475 * t735;
t662 = 0.2e1 * qJDD(1) * pkin(1) - t688;
t661 = pkin(5) * t611 - qJ(6) * t612 + t614;
t658 = t787 + (-qJD(1) * t559 - t528) * qJD(2);
t654 = g(1) * t510 - g(2) * t512 - t611 * t621 - t695;
t650 = -g(1) * (-t509 * pkin(5) + qJ(6) * t510) - g(2) * (t511 * pkin(5) - qJ(6) * t512);
t478 = qJD(1) * t667 + qJDD(1) * t678 + t587;
t519 = t607 + t667;
t649 = qJD(1) * t519 + qJDD(1) * t559 + t478 + t688;
t493 = -pkin(4) * t712 + (-pkin(7) - t599) * t745;
t647 = t427 * t571 + t654;
t505 = pkin(7) * t706 - t718;
t513 = t704 - t785;
t646 = qJD(2) * t681 - t505 * t633 + t513 * t630;
t644 = t537 * t670 - t611 * t805;
t481 = -qJD(1) * t551 + t691;
t643 = t481 + (-qJ(3) * t728 - t685) * t633 - t790;
t446 = pkin(4) * t719 + qJD(1) * t493 + t691;
t619 = t634 * pkin(7);
t598 = -pkin(4) * t797 - pkin(5);
t594 = pkin(4) * t628 + qJ(6);
t592 = g(1) * t770;
t585 = qJ(3) * t767;
t583 = qJ(3) * t770;
t549 = -qJ(3) * t748 + t608;
t532 = -t629 * t775 + t768;
t530 = t723 + t771;
t516 = t546 * t633;
t515 = t628 * t776 - t633 * t715;
t479 = pkin(5) * t546 + qJ(6) * t545 + t595;
t458 = -pkin(5) * t515 + qJ(6) * t516 + t526;
t457 = -qJD(2) * t724 + t546 * t802 + t797 * t713;
t456 = t545 * t802 + t546 * t745;
t443 = pkin(4) * t544 + t448;
t441 = -t630 * pkin(5) - t674;
t439 = qJ(6) * t630 + t812;
t423 = -pkin(5) * t457 - qJ(6) * t456 + qJD(6) * t516 + t493;
t416 = -pkin(5) * t743 - t803;
t415 = qJ(6) * t743 + qJD(6) * t630 + t664;
t414 = t432 * pkin(5) + t431 * qJ(6) - qJD(6) * t671 + t446;
t1 = [((t543 * t629 + t544 * t632) * t745 + (((t544 - t710) * qJD(4) + t719) * t629 + (-t740 + t666 + (t725 + (qJD(4) - 0.2e1 * t749) * qJD(2)) * t629) * t632) * t633) * MDP(16) + (-t428 * t743 - t526 * t431 - t446 * t516 + t490 * t456 + t493 * t671 - t537 * t812 - t571 * t664 - t630 * t695 + t687) * MDP(28) + (t630 * t676 + t633 * t662 + t592) * MDP(9) + (t630 * t658 + t633 * t649 - t592) * MDP(12) + (-t544 * t711 + (t630 * t733 + t633 * t666 - t655 * t776) * t629) * MDP(15) + (-g(2) * t767 + t412 * t515 - t413 * t516 - t415 * t483 + t416 * t671 + t425 * t456 + t426 * t457 - t431 * t441 - t432 * t439 + t592) * MDP(30) + (-t431 * t630 + t456 * t571 - t516 * t537 + t671 * t743) * MDP(24) + (t412 * t630 + t414 * t516 + t415 * t571 - t423 * t671 + t426 * t743 + t431 * t458 - t436 * t456 + t439 * t537 - t687) * MDP(31) + (t431 * t516 + t456 * t671) * MDP(22) + (-t431 * t515 + t432 * t516 - t456 * t483 + t457 * t671) * MDP(23) + (-t720 * t588 - t755 * t542 - t721 * t630 - t465 * t743 - t551 * t544 + g(1) * t531 - g(2) * t529 + (-t567 * t683 - t633 * t734) * t632 + ((-qJDD(1) * t567 - t481) * t633 + t682 * qJD(4) + (t527 * t630 + t567 * t696) * qJD(2)) * t629) * MDP(21) + (t699 * t588 + (-t538 * t629 + t548) * t542 + t700 * t630 + t551 * t543 + t567 * t719 + t481 * t769 - g(1) * t532 - g(2) * t530 + (t464 * t633 + t698 * t774) * qJD(2) + (-t682 * t632 + (-t508 * t630 - t566 * t588 + t633 * t698) * t629) * qJD(4)) * MDP(20) + ((t625 + t626) * qJDD(1) * pkin(7) + t646 - t685) * MDP(11) + t685 * MDP(3) + ((t697 * t744 - t719) * t630 + (t697 * t739 - t522 + t746) * t633) * MDP(18) + qJDD(1) * MDP(1) + (t412 * t439 + t426 * t415 + t414 * t458 + t436 * t423 + t413 * t441 + t425 * t416 - g(1) * (pkin(5) * t512 + qJ(6) * t511 + t599 * t634 + t619) - g(2) * (pkin(4) * t723 + pkin(5) * t510 + qJ(6) * t509 - t634 * t766 + t693) + (-g(1) * (-pkin(4) * t777 + t678 + t766) - g(2) * t599) * t631) * MDP(32) + (t427 * t743 + t526 * t432 - t446 * t515 - t490 * t457 + t493 * t483 + t674 * t537 + t571 * t803 - t694 * t630 + t686) * MDP(27) + (pkin(7) * t646 - g(1) * t619 - g(2) * t693 + t478 * t559 + t528 * t519 - t678 * t795) * MDP(14) + (t676 * t633 + (-t662 - t795) * t630) * MDP(10) + (t658 * t633 + (-t649 + t795) * t630) * MDP(13) + (-t791 + t795) * MDP(2) + 0.2e1 * (t630 * t725 - t729 * t751) * MDP(5) + (-t432 * t630 + t457 * t571 - t483 * t743 + t515 * t537) * MDP(25) + (-t413 * t630 - t414 * t515 - t416 * t571 + t423 * t483 - t425 * t743 + t432 * t458 - t436 * t457 - t441 * t537 + t686) * MDP(29) + (t542 * t630 + t588 * t743) * MDP(19) + (t537 * t630 + t571 * t743) * MDP(26) + (t633 * t733 + (-t588 * t737 - t630 * t683) * t632 + ((-t542 - t726) * t633 + (t588 + t696) * t745) * t629) * MDP(17) + (qJDD(1) * t625 + 0.2e1 * t630 * t705) * MDP(4) + (qJDD(2) * t630 + t633 * t637) * MDP(6) + (qJDD(2) * t633 - t630 * t637) * MDP(7); ((-pkin(2) * t630 + t786) * qJDD(1) + ((-t565 - t624) * t630 + (-t554 + t701) * t633) * qJD(1)) * MDP(11) + (t414 * t546 + t432 * t479 + t436 * t757 + t483 * t764 - t571 * t762 + t644) * MDP(29) + (t595 * t432 + t446 * t546 + t753 * t483 + t757 * t490 + t571 * t813 + t644) * MDP(27) + (t431 * t670 - t432 * t492 - t483 * t763 + t671 * t762 - t801) * MDP(30) + (t412 * t492 + t414 * t479 - t413 * t670 - g(1) * t585 - g(2) * t583 - t709 + t764 * t436 + t763 * t426 + t762 * t425 + (g(3) * t635 - t661 * t685) * t633 + (-g(3) * t661 + t685 * (pkin(2) - t635)) * t630) * MDP(32) + (t431 * t546 + t432 * t545 - t483 * t758 - t671 * t757) * MDP(23) + (-t549 * MDP(12) - t588 * MDP(19) - t464 * MDP(20) + t465 * MDP(21) - MDP(24) * t671 + t483 * MDP(25) - t571 * MDP(26) - t427 * MDP(27) + t428 * MDP(28) + t425 * MDP(29) - t426 * MDP(31)) * t748 + (t601 + 0.2e1 * t622 + 0.2e1 * t623 + (qJD(1) * t549 - g(3)) * t630 + (qJD(1) * t528 - t685) * t633) * MDP(13) + (-t595 * t431 - t446 * t545 + t758 * t490 + t571 * t814 + t753 * t671 - t821) * MDP(28) + (t414 * t545 + t431 * t479 - t758 * t436 + t763 * t571 - t764 * t671 + t821) * MDP(31) + (pkin(1) * t772 + t692) * MDP(9) + qJDD(2) * MDP(8) + (-t505 * qJ(3) - t565 * qJD(3) - t513 * pkin(2) - t528 * t549 - g(1) * (-pkin(2) * t773 + t585) - g(2) * (-pkin(2) * t775 + t583) - t709 - t681 * qJD(1) * pkin(7)) * MDP(14) + (-t588 * t739 + t522 + (-t544 * t633 - t588 * t777) * qJD(1)) * MDP(17) + (-t588 * t738 - t778 + (-t543 * t633 - t588 * t774) * qJD(1)) * MDP(18) + t804 * MDP(22) + (t669 - 0.2e1 * t785) * MDP(12) + MDP(7) * t725 + MDP(6) * t726 - MDP(4) * t722 + (t790 - t601 + (pkin(1) * t638 + t685) * t633) * MDP(10) + (t756 * t588 + t731 * t544 + (qJ(3) * t655 - t527 * t588 - t765) * t629 + (-t588 * t736 + t622 + t643) * t632) * MDP(21) + (qJ(3) * t719 - t535 * t588 - t731 * t543 + (t734 + t765 + (t527 - t624) * t749) * t632 + ((t514 - t736) * t588 + t643) * t629) * MDP(20) + ((-t544 * qJD(4) + (-t544 + t744) * t749 - t719) * t632 + (-t740 - t632 * qJDD(2) + t677 * t629 + (0.2e1 * t711 + (-t543 - t732) * t630) * qJD(1)) * t629) * MDP(16) + (-t683 * t632 ^ 2 + (t632 * t656 - t780) * t629) * MDP(15) + (t502 * t571 + t759) * MDP(24) + (-t571 * t757 - t781) * MDP(25) + t751 * MDP(5) * t638; MDP(11) * t726 + (qJDD(2) + t722) * MDP(12) + (-t625 * t638 - t637) * MDP(13) + (qJD(2) * t565 + t586 + t669 - t785) * MDP(14) + (t522 + t746) * MDP(20) + (-t733 - t778) * MDP(21) + (-t747 + t759) * MDP(27) + (-t747 - t494) * MDP(29) + (-t432 * t546 - t483 * t757 - t804) * MDP(30) + (-qJD(2) * t436 + t801) * MDP(32) + (t502 * MDP(27) + MDP(29) * t758 + t757 * t808) * t571 + t808 * (qJD(2) * t671 + t781) + (-t629 * MDP(20) - t632 * MDP(21)) * t588 ^ 2; -t544 * t543 * MDP(15) + (-t543 ^ 2 + t544 ^ 2) * MDP(16) + (-t543 * t588 + t641) * MDP(17) + (t645 + t780) * MDP(18) + t542 * MDP(19) + (g(3) * t769 - t527 * t544 + t700 + (-qJD(4) + t588) * t465 + t809) * MDP(20) + (g(1) * t530 - g(2) * t532 + t464 * t588 - t527 * t543 + (t741 - t621) * t629 - t721) * MDP(21) + (t429 * t571 + (-t483 * t544 + t537 * t797 - t571 * t735) * pkin(4) + t806) * MDP(27) + (t430 * t571 + t819 + (-t628 * t537 - t544 * t671 - t571 * t707) * pkin(4) + t654) * MDP(28) + (-t443 * t483 - t537 * t598 - t571 * t689 - t642) * MDP(29) + (-t431 * t598 - t432 * t594 + (t426 + t689) * t671 + (t425 - t754) * t483) * MDP(30) + (t443 * t671 + t537 * t594 + t571 * t754 - t654 + t811 - t820) * MDP(31) + (t412 * t594 + t413 * t598 - t436 * t443 - t425 * t429 - g(3) * (-t589 + t668) + t754 * t426 + (t425 * t735 + t809) * pkin(4) + t650) * MDP(32) + t822; (t784 + t806) * MDP(27) + (t647 + t819) * MDP(28) + (-t448 * t483 + t523 - t642 + t784) * MDP(29) + (pkin(5) * t431 - qJ(6) * t432 + (t426 - t428) * t671 + (t425 - t730) * t483) * MDP(30) + (t448 * t671 + 0.2e1 * t521 + 0.2e1 * t560 - t647 - t820) * MDP(31) + (-t413 * pkin(5) - g(3) * t668 + t412 * qJ(6) - t425 * t428 + t426 * t730 - t436 * t448 + t650) * MDP(32) + t822; (-t537 + t825) * MDP(29) + t424 * MDP(30) + (-t571 ^ 2 - t799) * MDP(31) + (-t426 * t571 + t642) * MDP(32);];
tau  = t1;
