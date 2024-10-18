% Calculate vector of inverse dynamics joint torques for
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR14_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:42
% EndTime: 2024-09-27 18:43:49
% DurationCPUTime: 4.72s
% Computational Cost: add. (6336->447), mult. (11303->585), div. (0->0), fcn. (8397->26), ass. (0->263)
t631 = sin(pkin(5));
t639 = cos(qJ(4));
t640 = cos(qJ(3));
t748 = t639 * t640;
t714 = t631 * t748;
t634 = sin(qJ(4));
t635 = sin(qJ(3));
t755 = t634 * t635;
t801 = -t631 * t755 + t714;
t675 = t634 * t640 + t635 * t639;
t800 = t675 * t631;
t636 = sin(qJ(2));
t641 = cos(qJ(2));
t632 = cos(pkin(5));
t735 = qJD(3) * t640;
t709 = t632 * t735;
t759 = t632 * t635;
t773 = pkin(1) * qJD(1);
t799 = (-t636 * t759 + t640 * t641) * t773 - pkin(2) * t709;
t624 = qJDD(1) + qJDD(2);
t626 = qJD(1) + qJD(2);
t786 = qJD(3) + qJD(4);
t646 = t786 * t675;
t451 = -t624 * t714 + t631 * (t624 * t755 + t646 * t626);
t721 = t631 * (-pkin(8) - pkin(9));
t692 = t635 * t721;
t798 = -qJD(3) * t692 + t799;
t758 = t632 * t640;
t661 = pkin(1) * (-t635 * t641 - t636 * t758);
t541 = qJD(1) * t661;
t602 = pkin(2) * t759;
t797 = (t640 * t721 - t602) * qJD(3) - t541;
t796 = t801 * t626;
t638 = cos(qJ(5));
t531 = t800 * t626;
t633 = sin(qJ(5));
t770 = t531 * t633;
t470 = t638 * t796 - t770;
t763 = t626 * t632;
t589 = qJD(3) + t763;
t586 = qJD(4) + t589;
t580 = qJD(5) + t586;
t794 = t470 * t580;
t678 = -t638 * t531 - t633 * t796;
t793 = t580 * t678;
t781 = pkin(2) * t626;
t585 = t641 * t773 + t781;
t555 = t585 * t758;
t723 = t636 * t773;
t764 = t626 * t631;
t566 = pkin(8) * t764 + t723;
t687 = pkin(9) * t764 + t566;
t483 = -t687 * t635 + t555;
t472 = pkin(3) * t589 + t483;
t720 = t585 * t759;
t484 = t687 * t640 + t720;
t480 = t639 * t484;
t679 = -t472 * t634 - t480;
t775 = pkin(10) * t796;
t430 = -t679 + t775;
t732 = qJD(5) * t633;
t428 = t430 * t732;
t762 = t626 * t640;
t528 = (-pkin(3) * t762 - t585) * t631;
t475 = -pkin(4) * t796 + t528;
t629 = qJ(3) + qJ(4);
t708 = pkin(5) - t629;
t688 = -qJ(5) + t708;
t670 = sin(t688);
t615 = pkin(5) + t629;
t608 = qJ(5) + t615;
t725 = sin(t608) / 0.2e1;
t558 = t725 - t670 / 0.2e1;
t622 = qJ(5) + t629;
t610 = cos(t622);
t630 = qJ(1) + qJ(2);
t617 = sin(t630);
t619 = cos(t630);
t505 = t558 * t617 - t610 * t619;
t671 = cos(t688);
t677 = -t558 * t619 - t610 * t617;
t726 = cos(t608) / 0.2e1;
t652 = -t475 * t470 - g(1) * t505 - g(2) * t677 - g(3) * (t726 - t671 / 0.2e1) + t428;
t609 = sin(t622);
t654 = t671 / 0.2e1 + t726;
t500 = t609 * t617 - t619 * t654;
t503 = -t619 * t609 - t617 * t654;
t450 = t624 * t800 + t786 * t796;
t767 = t624 * t632;
t587 = qJDD(3) + t767;
t584 = qJDD(4) + t587;
t783 = pkin(1) * t641;
t613 = qJDD(1) * t783;
t784 = pkin(1) * t636;
t722 = qJD(2) * t784;
t782 = pkin(2) * t624;
t550 = -qJD(1) * t722 + t613 + t782;
t535 = t550 * t758;
t728 = qJDD(1) * t636;
t737 = qJD(2) * t641;
t768 = t624 * t631;
t539 = pkin(8) * t768 + (qJD(1) * t737 + t728) * pkin(1);
t438 = pkin(3) * t587 + t535 + (-pkin(9) * t768 - t539) * t635 - t484 * qJD(3);
t712 = t640 * t539 + t550 * t759 + t585 * t709;
t736 = qJD(3) * t635;
t446 = -t566 * t736 + t712;
t711 = t626 * t736;
t766 = t624 * t640;
t777 = pkin(9) * t631;
t442 = (-t711 + t766) * t777 + t446;
t650 = t679 * qJD(4) + t639 * t438 - t634 * t442;
t404 = pkin(4) * t584 - pkin(10) * t450 + t650;
t733 = qJD(4) * t639;
t734 = qJD(4) * t634;
t689 = -t634 * t438 - t639 * t442 - t472 * t733 + t484 * t734;
t407 = -pkin(10) * t451 - t689;
t704 = t638 * t404 - t633 * t407;
t648 = t475 * t678 - g(1) * t503 + g(2) * t500 - g(3) * (t725 + t670 / 0.2e1) + t704;
t571 = qJDD(5) + t584;
t561 = t571 * MDP(25);
t790 = t561 + t470 * MDP(21) * t678 + (-t470 ^ 2 + t678 ^ 2) * MDP(22);
t789 = t631 * (-pkin(3) * t736 + t723);
t603 = pkin(2) * t758;
t621 = t632 * pkin(3);
t524 = t603 + t621 + t692;
t760 = t631 * t640;
t600 = pkin(9) * t760;
t742 = pkin(8) * t760 + t602;
t537 = t600 + t742;
t745 = t634 * t524 + t639 * t537;
t788 = -qJD(4) * t745 + t634 * t798 + t797 * t639;
t787 = -t524 * t733 + t537 * t734 - t797 * t634 + t639 * t798;
t612 = pkin(2) + t783;
t583 = t612 * t758;
t588 = pkin(8) * t631 + t784;
t706 = -t588 - t777;
t498 = t706 * t635 + t583 + t621;
t582 = t612 * t759;
t743 = t640 * t588 + t582;
t506 = t600 + t743;
t746 = t634 * t498 + t639 * t506;
t523 = t531 * pkin(10);
t478 = t634 * t484;
t702 = t639 * t472 - t478;
t429 = -t523 + t702;
t729 = qJD(3) - t589;
t785 = g(3) * t631 + t729 * t566;
t703 = t450 * t633 + t638 * t451;
t417 = -t678 * qJD(5) + t703;
t780 = pkin(3) * t634;
t779 = pkin(3) * t640;
t778 = pkin(4) * t801;
t490 = t786 * (t748 - t755) * t631;
t776 = pkin(10) * t490;
t772 = MDP(9) * t631;
t427 = pkin(4) * t586 + t429;
t771 = t427 * t638;
t625 = t631 ^ 2;
t765 = t625 * t640;
t761 = t631 * t635;
t757 = t633 * t404;
t756 = t633 * t571;
t754 = t634 * t638;
t751 = t635 * t640;
t750 = t638 * t430;
t749 = t638 * t571;
t747 = t639 * t483 - t478;
t744 = t640 * pkin(1) * t737 + t612 * t709;
t741 = g(1) * t619 + g(2) * t617;
t627 = t635 ^ 2;
t740 = -t640 ^ 2 + t627;
t739 = MDP(10) * t631;
t731 = qJD(5) * t638;
t570 = t584 * MDP(18);
t730 = t587 * MDP(11);
t599 = cos(t615) / 0.2e1;
t604 = cos(t708);
t727 = t604 / 0.2e1 + t599;
t724 = sin(t615) / 0.2e1;
t718 = t626 * t761;
t713 = t638 * t450 - t633 * t451 + t731 * t796;
t710 = t631 * t736;
t707 = -t585 - t781;
t491 = t646 * t631;
t593 = pkin(3) * t710;
t476 = pkin(4) * t491 + t593;
t705 = t632 * pkin(4) - pkin(10) * t800;
t701 = -t483 * t634 - t480;
t700 = t639 * t498 - t506 * t634;
t699 = t639 * t524 - t537 * t634;
t698 = t587 + t767;
t697 = t589 + t763;
t696 = qJD(5) * t427 + t407;
t695 = qJD(1) * (-qJD(2) + t626);
t694 = qJD(2) * (-qJD(1) - t626);
t693 = t632 * t722;
t686 = g(1) * t617 - g(2) * t619 + t613;
t685 = sin(t708);
t576 = (-pkin(2) - t779) * t631;
t452 = t699 + t705;
t489 = t491 * pkin(10);
t684 = -qJD(5) * t452 + t489 + t787;
t543 = t801 * pkin(10);
t454 = t543 + t745;
t683 = qJD(5) * t454 + t776 - t788;
t556 = (-t612 - t779) * t631;
t682 = -t633 * t427 - t750;
t443 = t700 + t705;
t444 = t543 + t746;
t681 = t443 * t638 - t444 * t633;
t680 = t443 * t633 + t444 * t638;
t487 = t633 * t800 - t638 * t801;
t488 = t633 * t801 + t638 * t800;
t564 = t724 - t685 / 0.2e1;
t618 = cos(t629);
t676 = -t564 * t619 - t617 * t618;
t515 = t564 * t617 - t618 * t619;
t673 = qJD(3) * (-t612 * t626 - t585);
t672 = -t631 * t723 + t476;
t546 = -t617 * t640 - t619 * t759;
t548 = -t617 * t759 + t619 * t640;
t669 = -g(1) * t546 - g(2) * t548 + (-t635 * t539 + t535 + (-t566 * t640 - t720) * qJD(3)) * t632 + t550 * t765;
t668 = t626 * t723 + t782;
t473 = (t706 * qJD(3) - t693) * t635 + t744;
t657 = qJD(2) * t661;
t474 = t657 + (t706 * t640 - t582) * qJD(3);
t667 = t639 * t473 + t634 * t474 + t498 * t733 - t506 * t734;
t416 = -t531 * t732 + t713;
t665 = t682 * qJD(5);
t663 = t612 * t624 - t626 * t722;
t545 = t617 * t635 - t619 * t758;
t547 = -t617 * t758 - t619 * t635;
t662 = -g(1) * t545 - g(2) * t547 - t446 * t632;
t485 = t626 * t593 + (-pkin(3) * t766 - t550) * t631;
t432 = t488 * qJD(5) + t490 * t633 + t638 * t491;
t433 = pkin(4) * t451 + t485;
t660 = -g(1) * t677 + g(2) * t505 + (t665 + t704) * t632 + t475 * t432 + t433 * t487;
t659 = -g(1) * t676 + g(2) * t515 - t485 * t801 + t528 * t491 + t650 * t632;
t431 = -t487 * qJD(5) + t490 * t638 - t491 * t633;
t656 = -g(1) * t500 - g(2) * t503 - (t696 * t638 - t428 + t757) * t632 + t475 * t431 + t433 * t488;
t616 = sin(t629);
t510 = t616 * t617 - t619 * t727;
t513 = -t619 * t616 - t617 * t727;
t655 = -g(1) * t510 - g(2) * t513 + t485 * t800 + t528 * t490 + t632 * t689;
t651 = -qJD(4) * t746 - t473 * t634 + t639 * t474;
t647 = -g(1) * t515 - g(2) * t676 - g(3) * (t599 - t604 / 0.2e1) - t528 * t796 + t689;
t645 = -t531 * t796 * MDP(14) + (t416 - t794) * MDP(23) + (-t417 - t793) * MDP(24) + (-t586 * t796 + t450) * MDP(16) + (t531 * t586 - t451) * MDP(17) + (t531 ^ 2 - t796 ^ 2) * MDP(15) + t570 + t790;
t644 = (-t416 * t487 - t417 * t488 + t431 * t470 + t432 * t678) * MDP(22) + (t416 * t488 - t431 * t678) * MDP(21) + (t431 * t580 + t488 * t571) * MDP(23) + (-t432 * t580 - t487 * t571) * MDP(24) + (t450 * t801 - t451 * t800 + t490 * t796 - t491 * t531) * MDP(15) + (t490 * t586 + t584 * t800) * MDP(16) + (-t491 * t586 + t584 * t801) * MDP(17) + (t450 * t800 + t490 * t531) * MDP(14) + (t698 * t640 - t697 * t736) * t739 + (t698 * t635 + t697 * t735) * t772 + t624 * MDP(4) + (0.2e1 * (-t740 * t626 * qJD(3) + t624 * t751) * MDP(8) + (t624 * t627 + 0.2e1 * t640 * t711) * MDP(7)) * t625 + (t450 * MDP(16) - t451 * MDP(17) + t416 * MDP(23) - t417 * MDP(24) + t561 + t570 + t730) * t632;
t643 = -g(1) * t513 + g(2) * t510 - g(3) * (t724 + t685 / 0.2e1) - t528 * t531 + t650;
t642 = cos(qJ(1));
t637 = sin(qJ(1));
t611 = pkin(3) * t639 + pkin(4);
t594 = t631 * t722;
t554 = t594 + t593;
t507 = t576 - t778;
t497 = t556 - t778;
t496 = pkin(3) * t718 + pkin(4) * t531;
t462 = t476 + t594;
t435 = -t523 + t747;
t434 = t701 - t775;
t414 = t651 - t776;
t413 = -t489 + t667;
t1 = [(t556 * t450 + t554 * t531 - t746 * t584 - t667 * t586 + t655) * MDP(20) + t644 + qJDD(1) * MDP(1) + (g(1) * t637 - g(2) * t642) * MDP(2) + (g(1) * t642 + g(2) * t637) * MDP(3) + ((t624 * t641 + t636 * t694) * pkin(1) + t686) * MDP(5) + (((-qJDD(1) - t624) * t636 + t641 * t694) * pkin(1) + t741) * MDP(6) + ((-qJD(3) * t743 + t657) * t589 + (-t588 * t635 + t583) * t587 + (t635 * t673 + t663 * t640) * t625 + t669) * MDP(12) + (-((-qJD(3) * t588 - t693) * t635 + t744) * t589 - t743 * t587 + (t640 * t673 + (-t550 - t663) * t635) * t625 + t662) * MDP(13) + ((-t680 * qJD(5) - t413 * t633 + t414 * t638) * t580 + t681 * t571 - t462 * t470 + t497 * t417 + t660) * MDP(26) + (-(t681 * qJD(5) + t413 * t638 + t414 * t633) * t580 - t680 * t571 - t462 * t678 + t497 * t416 + t656) * MDP(27) + (t556 * t451 - t554 * t796 + t700 * t584 + t651 * t586 + t659) * MDP(19); ((t452 * t638 - t454 * t633) * t571 + t507 * t417 + (t684 * t633 - t683 * t638) * t580 - t672 * t470 + t660) * MDP(26) + (-(t452 * t633 + t454 * t638) * t571 + t507 * t416 + (t683 * t633 + t684 * t638) * t580 - t672 * t678 + t656) * MDP(27) + (t576 * t451 + t699 * t584 + t586 * t788 + t789 * t796 + t659) * MDP(19) + (-t742 * t587 + (pkin(8) * t710 + t799) * t589 + (t707 * t735 + (-t550 - t668) * t635) * t625 + t662) * MDP(13) + (t576 * t450 - t531 * t789 - t745 * t584 + t586 * t787 + t655) * MDP(20) + ((-pkin(8) * t761 + t603) * t587 - t541 * t589 + t668 * t765 + (t707 * t635 * t625 - t742 * t589) * qJD(3) + t669) * MDP(12) + t644 + ((t641 * t695 - t728) * pkin(1) + t741) * MDP(6) + (t695 * t784 + t686) * MDP(5); (t611 * t749 - (t434 * t638 - t435 * t633) * t580 + t496 * t470 + (-t634 * t756 + (-t633 * t639 - t754) * t580 * qJD(4)) * pkin(3) + ((-pkin(3) * t754 - t611 * t633) * t580 + t682) * qJD(5) + t648) * MDP(26) + t730 + (-g(1) * t547 + g(2) * t545 + t535 - t785 * t640 + (-t539 + (t625 * t626 - t729 * t632) * t585) * t635) * MDP(12) + (t625 * t585 * t762 + g(1) * t548 - g(2) * t546 + t555 * t589 + t635 * t785 - t712) * MDP(13) + (-t701 * t586 + (t584 * t639 - t586 * t734 + t718 * t796) * pkin(3) + t643) * MDP(19) + t645 + (t496 * t678 + (-t611 * t571 - t404 + (t434 - (-qJD(4) - qJD(5)) * t780) * t580) * t633 + (-t571 * t780 + (-pkin(3) * t733 - qJD(5) * t611 + t435) * t580 - t696) * t638 + t652) * MDP(27) + (-t729 * t635 * t626 + t766) * t739 + (t624 * t635 + t729 * t762) * t772 + (t747 * t586 + (-t531 * t718 - t634 * t584 - t586 * t733) * pkin(3) + t647) * MDP(20) + (-MDP(7) * t751 + MDP(8) * t740) * t626 ^ 2 * t625; (t702 * t586 + t647) * MDP(20) + (-t679 * t586 + t643) * MDP(19) + (-(-t429 * t633 - t750) * t580 + t665 + (t470 * t531 - t580 * t732 + t749) * pkin(4) + t648) * MDP(26) + t645 + ((-t430 * t580 - t404) * t633 + (t429 * t580 - t696) * t638 + (t531 * t678 - t580 * t731 - t756) * pkin(4) + t652) * MDP(27); (t713 - t794) * MDP(23) + (-t703 - t793) * MDP(24) + (-t682 * t580 + t648) * MDP(26) + (-t638 * t407 - t757 + (-t430 * t633 + t771) * t580 + t652) * MDP(27) + (-MDP(23) * t770 + t678 * MDP(24) + t682 * MDP(26) - MDP(27) * t771) * qJD(5) + t790;];
tau = t1;
