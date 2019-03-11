% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:43
% EndTime: 2019-03-09 13:26:04
% DurationCPUTime: 15.17s
% Computational Cost: add. (11607->633), mult. (27195->825), div. (0->0), fcn. (21487->18), ass. (0->282)
t687 = sin(pkin(11));
t688 = cos(pkin(11));
t693 = sin(qJ(2));
t698 = cos(qJ(2));
t648 = t687 * t698 + t688 * t693;
t637 = t648 * qJD(1);
t692 = sin(qJ(4));
t697 = cos(qJ(4));
t761 = t697 * qJD(2);
t597 = t637 * t692 - t761;
t599 = qJD(2) * t692 + t637 * t697;
t691 = sin(qJ(5));
t696 = cos(qJ(5));
t531 = t696 * t597 + t599 * t691;
t695 = cos(qJ(6));
t690 = sin(qJ(6));
t723 = t597 * t691 - t696 * t599;
t813 = t723 * t690;
t484 = -t695 * t531 + t813;
t636 = t648 * qJD(2);
t758 = qJDD(1) * t698;
t759 = qJDD(1) * t693;
t726 = -t687 * t759 + t688 * t758;
t576 = qJD(1) * t636 + qJDD(4) - t726;
t571 = qJDD(5) + t576;
t562 = qJDD(6) + t571;
t724 = t531 * t690 + t695 * t723;
t854 = t562 * MDP(31) + (-t484 ^ 2 + t724 ^ 2) * MDP(28) + t484 * MDP(27) * t724;
t768 = qJD(1) * t693;
t794 = t688 * t698;
t635 = qJD(1) * t794 - t687 * t768;
t829 = qJD(4) + qJD(5);
t853 = t635 - t829;
t760 = qJD(1) * qJD(2);
t750 = t698 * t760;
t751 = t693 * t760;
t582 = qJDD(1) * t648 - t687 * t751 + t688 * t750;
t767 = qJD(4) * t692;
t523 = qJD(4) * t761 + t692 * qJDD(2) + t697 * t582 - t637 * t767;
t524 = qJD(4) * t599 - t697 * qJDD(2) + t582 * t692;
t764 = qJD(5) * t696;
t765 = qJD(5) * t691;
t460 = t696 * t523 - t691 * t524 - t597 * t764 - t599 * t765;
t703 = qJD(5) * t723 - t523 * t691 - t696 * t524;
t762 = qJD(6) * t695;
t755 = t695 * t460 - t531 * t762 + t690 * t703;
t763 = qJD(6) * t690;
t431 = t723 * t763 + t755;
t625 = qJD(4) - t635;
t619 = qJD(5) + t625;
t744 = t460 * t690 - t695 * t703;
t704 = qJD(6) * t724 - t744;
t615 = qJD(6) + t619;
t842 = t615 * t724;
t843 = t484 * t615;
t852 = t571 * MDP(24) + (-t531 ^ 2 + t723 ^ 2) * MDP(21) + (t531 * t619 + t460) * MDP(22) + (-t619 * t723 + t703) * MDP(23) - t531 * MDP(20) * t723 + (t704 - t842) * MDP(30) + (t431 - t843) * MDP(29) + t854;
t674 = pkin(2) * t698 + pkin(1);
t658 = -qJD(1) * t674 + qJD(3);
t543 = -pkin(3) * t635 - pkin(8) * t637 + t658;
t821 = qJ(3) + pkin(7);
t660 = t821 * t693;
t652 = qJD(1) * t660;
t820 = qJD(2) * pkin(2);
t644 = -t652 + t820;
t661 = t821 * t698;
t653 = qJD(1) * t661;
t795 = t688 * t653;
t580 = t687 * t644 + t795;
t568 = qJD(2) * pkin(8) + t580;
t504 = t697 * t543 - t568 * t692;
t487 = -pkin(9) * t599 + t504;
t467 = pkin(4) * t625 + t487;
t505 = t543 * t692 + t568 * t697;
t488 = -pkin(9) * t597 + t505;
t478 = t696 * t488;
t450 = t467 * t691 + t478;
t845 = pkin(10) * t531;
t444 = t450 - t845;
t442 = t444 * t763;
t640 = t687 * t653;
t579 = t644 * t688 - t640;
t567 = -qJD(2) * pkin(3) - t579;
t525 = pkin(4) * t597 + t567;
t471 = pkin(5) * t531 + t525;
t686 = qJ(4) + qJ(5);
t682 = qJ(6) + t686;
t671 = sin(t682);
t672 = cos(t682);
t699 = cos(qJ(1));
t683 = qJ(2) + pkin(11);
t676 = cos(t683);
t694 = sin(qJ(1));
t801 = t676 * t694;
t602 = t671 * t699 - t672 * t801;
t800 = t676 * t699;
t604 = t671 * t694 + t672 * t800;
t675 = sin(t683);
t825 = g(3) * t675;
t838 = g(1) * t604 - g(2) * t602 - t471 * t484 + t672 * t825 + t442;
t650 = t691 * t692 - t696 * t697;
t774 = t853 * t650;
t791 = t691 * t697;
t651 = t692 * t696 + t791;
t773 = t853 * t651;
t601 = t671 * t801 + t672 * t699;
t603 = -t671 * t800 + t672 * t694;
t581 = -qJD(2) * t637 + t726;
t712 = pkin(2) * t751 - qJDD(1) * t674 + qJDD(3);
t506 = -pkin(3) * t581 - pkin(8) * t582 + t712;
t503 = t697 * t506;
t749 = qJD(2) * t821;
t634 = -qJD(3) * t693 - t698 * t749;
t575 = qJDD(2) * pkin(2) + qJD(1) * t634 - qJDD(1) * t660;
t633 = qJD(3) * t698 - t693 * t749;
t584 = qJD(1) * t633 + qJDD(1) * t661;
t518 = t687 * t575 + t688 * t584;
t516 = qJDD(2) * pkin(8) + t518;
t435 = pkin(4) * t576 - pkin(9) * t523 - qJD(4) * t505 - t516 * t692 + t503;
t766 = qJD(4) * t697;
t717 = t692 * t506 + t697 * t516 + t543 * t766 - t568 * t767;
t440 = -pkin(9) * t524 + t717;
t746 = t696 * t435 - t691 * t440;
t705 = -qJD(5) * t450 + t746;
t425 = pkin(5) * t571 - pkin(10) * t460 + t705;
t735 = -t691 * t435 - t696 * t440 - t467 * t764 + t488 * t765;
t426 = pkin(10) * t703 - t735;
t747 = t695 * t425 - t690 * t426;
t837 = -g(1) * t603 + g(2) * t601 + t471 * t724 + t671 * t825 + t747;
t555 = pkin(2) * t768 + pkin(3) * t637 - pkin(8) * t635;
t547 = t697 * t555;
t587 = -t652 * t688 - t640;
t667 = pkin(2) * t687 + pkin(8);
t822 = pkin(9) + t667;
t748 = qJD(4) * t822;
t849 = pkin(4) * t637 - t587 * t692 + t547 + (-pkin(9) * t635 + t748) * t697;
t776 = t692 * t555 + t697 * t587;
t806 = t635 * t692;
t848 = -pkin(9) * t806 + t692 * t748 + t776;
t517 = t575 * t688 - t687 * t584;
t515 = -qJDD(2) * pkin(3) - t517;
t733 = g(1) * t699 + g(2) * t694;
t711 = -g(3) * t676 + t675 * t733;
t847 = -qJD(4) * t667 * t625 - t515 + t711;
t846 = t767 - t806;
t844 = pkin(10) * t723;
t588 = -t650 * t690 + t651 * t695;
t780 = qJD(6) * t588 + t690 * t774 - t695 * t773;
t753 = t648 * t766;
t647 = t687 * t693 - t794;
t639 = t647 * qJD(2);
t805 = t639 * t692;
t839 = t753 - t805;
t681 = cos(t686);
t797 = t681 * t694;
t680 = sin(t686);
t798 = t680 * t699;
t606 = -t676 * t797 + t798;
t796 = t681 * t699;
t799 = t680 * t694;
t608 = t676 * t796 + t799;
t836 = g(1) * t608 - g(2) * t606 + t525 * t531 + t681 * t825 + t735;
t605 = t676 * t799 + t796;
t607 = -t676 * t798 + t797;
t835 = -g(1) * t607 + g(2) * t605 + t525 * t723 + t680 * t825 + t705;
t563 = t651 * t648;
t578 = pkin(3) * t647 - pkin(8) * t648 - t674;
t566 = t697 * t578;
t596 = -t660 * t687 + t661 * t688;
t802 = t648 * t697;
t497 = pkin(4) * t647 - pkin(9) * t802 - t596 * t692 + t566;
t589 = t697 * t596;
t775 = t692 * t578 + t589;
t803 = t648 * t692;
t508 = -pkin(9) * t803 + t775;
t779 = t691 * t497 + t696 * t508;
t585 = -t652 * t687 + t795;
t734 = pkin(4) * t846 - t585;
t831 = t849 * t696;
t645 = t822 * t692;
t646 = t822 * t697;
t772 = -t691 * t645 + t696 * t646;
t830 = t645 * t764 + t646 * t765 + t691 * t849 + t848 * t696;
t586 = t695 * t650 + t651 * t690;
t781 = -qJD(6) * t586 + t690 * t773 + t695 * t774;
t828 = -t562 * t588 - t615 * t781;
t827 = -t571 * t651 - t619 * t774;
t826 = pkin(4) * t691;
t823 = g(3) * t698;
t476 = t691 * t488;
t449 = t696 * t467 - t476;
t443 = t449 + t844;
t441 = pkin(5) * t619 + t443;
t819 = t441 * t695;
t818 = t484 * t637;
t817 = t724 * t637;
t816 = t523 * t692;
t815 = t531 * t637;
t814 = t723 * t637;
t810 = t597 * t625;
t809 = t597 * t637;
t808 = t599 * t625;
t807 = t599 * t637;
t804 = t639 * t697;
t793 = t690 * t562;
t792 = t691 * t695;
t790 = t692 * t576;
t789 = t692 * t694;
t788 = t692 * t699;
t787 = t694 * t697;
t786 = t695 * t444;
t785 = t695 * t562;
t561 = t697 * t576;
t784 = t697 * t699;
t783 = t696 * t487 - t476;
t777 = -pkin(5) * t773 + t734;
t684 = t693 ^ 2;
t771 = -t698 ^ 2 + t684;
t757 = t693 * t820;
t668 = -pkin(2) * t688 - pkin(3);
t754 = t648 * t767;
t556 = pkin(3) * t636 + pkin(8) * t639 + t757;
t548 = t697 * t556;
t554 = t633 * t688 + t634 * t687;
t455 = pkin(9) * t804 + pkin(4) * t636 - t554 * t692 + t548 + (-t589 + (pkin(9) * t648 - t578) * t692) * qJD(4);
t716 = t697 * t554 + t692 * t556 + t578 * t766 - t596 * t767;
t463 = -pkin(9) * t839 + t716;
t745 = t696 * t455 - t463 * t691;
t743 = -t487 * t691 - t478;
t742 = t696 * t497 - t508 * t691;
t553 = t633 * t687 - t688 * t634;
t739 = -t696 * t645 - t646 * t691;
t595 = t688 * t660 + t661 * t687;
t738 = t625 * t697;
t737 = -qJD(4) * t543 - t516;
t736 = qJD(6) * t441 + t426;
t732 = g(1) * t694 - g(2) * t699;
t731 = -t586 * t562 - t615 * t780;
t730 = -t650 * t571 + t619 * t773;
t551 = pkin(4) * t803 + t595;
t729 = -t568 * t766 + t503;
t539 = -pkin(10) * t650 + t772;
t728 = pkin(5) * t637 + pkin(10) * t774 + qJD(5) * t772 + qJD(6) * t539 - t691 * t848 + t831;
t538 = -pkin(10) * t651 + t739;
t727 = -pkin(10) * t773 - qJD(6) * t538 + t830;
t430 = t690 * t441 + t786;
t564 = t650 * t648;
t509 = t695 * t563 - t564 * t690;
t510 = -t563 * t690 - t564 * t695;
t659 = -pkin(4) * t697 + t668;
t514 = pkin(4) * t839 + t553;
t722 = -t625 * t846 + t561;
t721 = -0.2e1 * pkin(1) * t760 - pkin(7) * qJDD(2);
t719 = -t754 - t804;
t715 = t691 * t455 + t696 * t463 + t497 * t764 - t508 * t765;
t713 = t567 * t625 - t667 * t576;
t464 = pkin(4) * t524 + t515;
t700 = qJD(2) ^ 2;
t708 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t700 + t732;
t701 = qJD(1) ^ 2;
t707 = pkin(1) * t701 - pkin(7) * qJDD(1) + t733;
t673 = pkin(4) * t696 + pkin(5);
t629 = t676 * t784 + t789;
t628 = -t676 * t788 + t787;
t627 = -t676 * t787 + t788;
t626 = t676 * t789 + t784;
t600 = pkin(5) * t650 + t659;
t512 = pkin(5) * t563 + t551;
t511 = pkin(4) * t599 - pkin(5) * t723;
t470 = -t639 * t791 - t691 * t754 - t765 * t803 + (t802 * t829 - t805) * t696;
t469 = -t563 * t829 + t650 * t639;
t459 = pkin(5) * t470 + t514;
t452 = -pkin(10) * t563 + t779;
t451 = pkin(5) * t647 + pkin(10) * t564 + t742;
t446 = t783 + t844;
t445 = t743 + t845;
t438 = qJD(6) * t510 + t469 * t690 + t695 * t470;
t437 = -qJD(6) * t509 + t469 * t695 - t470 * t690;
t436 = -pkin(5) * t703 + t464;
t429 = -t444 * t690 + t819;
t428 = -pkin(10) * t470 + t715;
t427 = pkin(5) * t636 - pkin(10) * t469 - qJD(5) * t779 + t745;
t1 = [(-t460 * t564 - t469 * t723) * MDP(20) + (t460 * t647 + t469 * t619 - t564 * t571 - t636 * t723) * MDP(22) + (-g(1) * t605 - g(2) * t607 - t450 * t636 + t551 * t460 - t464 * t564 + t525 * t469 - t514 * t723 - t571 * t779 - t619 * t715 + t647 * t735) * MDP(26) + (t518 * t596 + t580 * t554 - t517 * t595 - t579 * t553 - t712 * t674 + t658 * t757 - g(1) * (-t674 * t694 + t699 * t821) - g(2) * (t674 * t699 + t694 * t821)) * MDP(12) + (-t460 * t563 - t469 * t531 + t470 * t723 - t564 * t703) * MDP(21) + (-t470 * t619 - t531 * t636 - t563 * t571 + t647 * t703) * MDP(23) + (t514 * t531 - t551 * t703 + t464 * t563 + t525 * t470 + t745 * t619 + t742 * t571 + t746 * t647 + t449 * t636 - g(1) * t606 - g(2) * t608 + (-t450 * t647 - t619 * t779) * qJD(5)) * MDP(25) + (t431 * t510 - t437 * t724) * MDP(27) + (t431 * t647 + t437 * t615 + t510 * t562 - t636 * t724) * MDP(29) + (-g(1) * t601 - g(2) * t603 - t430 * t636 + t512 * t431 + t436 * t510 + t471 * t437 + t442 * t647 - t459 * t724 + (-(-qJD(6) * t452 + t427) * t615 - t451 * t562 - t425 * t647) * t690 + (-(qJD(6) * t451 + t428) * t615 - t452 * t562 - t736 * t647) * t695) * MDP(33) + (-t438 * t615 + t484 * t636 - t509 * t562 + t647 * t704) * MDP(30) + ((t427 * t695 - t428 * t690) * t615 + (t451 * t695 - t452 * t690) * t562 + t747 * t647 + t429 * t636 - t459 * t484 - t512 * t704 + t436 * t509 + t471 * t438 - g(1) * t602 - g(2) * t604 + ((-t451 * t690 - t452 * t695) * t615 - t430 * t647) * qJD(6)) * MDP(32) + (-t431 * t509 + t437 * t484 + t438 * t724 + t510 * t704) * MDP(28) + (-t524 * t647 - t597 * t636 - t625 * t839 - t648 * t790) * MDP(16) + (t576 * t647 + t625 * t636) * MDP(17) + (t571 * t647 + t619 * t636) * MDP(24) + (t562 * t647 + t615 * t636) * MDP(31) + (t693 * t721 + t698 * t708) * MDP(9) + (-t693 * t708 + t698 * t721) * MDP(10) + t732 * MDP(2) + t733 * MDP(3) + qJDD(1) * MDP(1) + (qJDD(1) * t684 + 0.2e1 * t693 * t750) * MDP(4) + (-t517 * t648 - t518 * t647 + t553 * t637 + t554 * t635 + t579 * t639 - t580 * t636 + t581 * t596 + t582 * t595 - t733) * MDP(11) + (-(-t597 * t697 - t599 * t692) * t639 + (-t816 - t524 * t697 + (t597 * t692 - t599 * t697) * qJD(4)) * t648) * MDP(14) + ((-t596 * t766 + t548) * t625 + t566 * t576 + t729 * t647 + t504 * t636 + t553 * t597 + t595 * t524 + t567 * t753 - g(1) * t627 - g(2) * t629 + ((-qJD(4) * t578 - t554) * t625 - t596 * t576 + t737 * t647 + t515 * t648 - t567 * t639) * t692) * MDP(18) + 0.2e1 * (t693 * t758 - t760 * t771) * MDP(5) + (t523 * t647 + t561 * t648 + t599 * t636 + t625 * t719) * MDP(15) + (t523 * t802 + t599 * t719) * MDP(13) + (-g(1) * t626 - g(2) * t628 - t505 * t636 + t515 * t802 + t595 * t523 + t553 * t599 + t567 * t719 - t576 * t775 - t625 * t716 - t647 * t717) * MDP(19) + (qJDD(2) * t693 + t698 * t700) * MDP(6) + (qJDD(2) * t698 - t693 * t700) * MDP(7); ((t579 - t587) * t635 + (t581 * t687 - t582 * t688) * pkin(2)) * MDP(11) - ((-t580 + t585) * MDP(11) + t615 * MDP(31) + t619 * MDP(24) + t625 * MDP(17) + t504 * MDP(18) + t449 * MDP(25) - t505 * MDP(19) + t429 * MDP(32) - t430 * MDP(33) - t450 * MDP(26)) * t637 + (t460 * t651 - t723 * t774) * MDP(20) + (g(3) * t693 + t698 * t707) * MDP(10) + (t659 * t460 + t464 * t651 + t774 * t525 - t772 * t571 + t619 * t830 - t680 * t711 - t723 * t734) * MDP(26) + (-(t538 * t690 + t539 * t695) * t562 + t600 * t431 + t436 * t588 + (t690 * t728 + t695 * t727) * t615 - t777 * t724 + t781 * t471 - t711 * t671) * MDP(33) + (-MDP(4) * t693 * t698 + MDP(5) * t771) * t701 + (t431 * t588 - t724 * t781) * MDP(27) + ((t538 * t695 - t539 * t690) * t562 - t600 * t704 + t436 * t586 + (t690 * t727 - t695 * t728) * t615 - t777 * t484 + t780 * t471 + t711 * t672) * MDP(32) + (-t431 * t586 + t484 * t781 + t588 * t704 + t724 * t780) * MDP(28) + (t722 + t809) * MDP(16) + (t730 + t815) * MDP(23) + MDP(7) * t758 + MDP(6) * t759 + (-t460 * t650 - t531 * t774 + t651 * t703 - t723 * t773) * MDP(21) + (-t659 * t703 + t464 * t650 + t739 * t571 + (-t646 * t764 + (qJD(5) * t645 + t848) * t691 - t831) * t619 + t734 * t531 - t773 * t525 + t711 * t681) * MDP(25) + (t693 * t707 - t823) * MDP(9) + (t579 * t585 - t580 * t587 + (-t823 + t517 * t688 + t518 * t687 + (-qJD(1) * t658 + t733) * t693) * pkin(2)) * MDP(12) + (t599 * t738 + t816) * MDP(13) + (t731 - t818) * MDP(30) + ((t523 - t810) * t697 + (-t524 - t808) * t692) * MDP(14) + (t814 - t827) * MDP(22) + (t817 - t828) * MDP(29) + (t668 * t524 - t547 * t625 - t585 * t597 + (t587 * t625 + t713) * t692 + t847 * t697) * MDP(18) + (t668 * t523 - t585 * t599 + t776 * t625 - t692 * t847 + t713 * t697) * MDP(19) + qJDD(2) * MDP(8) + (t625 * t738 + t790 - t807) * MDP(15); (-t635 ^ 2 - t637 ^ 2) * MDP(11) + (t579 * t637 - t580 * t635 + t712 - t732) * MDP(12) + (t722 - t809) * MDP(18) + (-t625 ^ 2 * t697 - t790 - t807) * MDP(19) + (t730 - t815) * MDP(25) + (t814 + t827) * MDP(26) + (t731 + t818) * MDP(32) + (t817 + t828) * MDP(33); (t783 * t619 + (-t571 * t691 + t599 * t723 - t619 * t764) * pkin(4) + t836) * MDP(26) + (t511 * t724 + (-t673 * t562 - t425 + (t445 - (-qJD(5) - qJD(6)) * t826) * t615) * t690 + (-t562 * t826 + (-pkin(4) * t764 - qJD(6) * t673 + t446) * t615 - t736) * t695 + t838) * MDP(33) + (-t743 * t619 + (-t531 * t599 + t571 * t696 - t619 * t765) * pkin(4) + t835) * MDP(25) + t576 * MDP(17) + (-t597 ^ 2 + t599 ^ 2) * MDP(14) + (g(1) * t629 - g(2) * t627 + t504 * t625 + t567 * t597 + t697 * t825 - t717) * MDP(19) + (-g(1) * t628 + g(2) * t626 + t505 * t625 - t567 * t599 + (t737 + t825) * t692 + t729) * MDP(18) + (t523 + t810) * MDP(15) + (-t524 + t808) * MDP(16) + t599 * t597 * MDP(13) + (t673 * t785 - (t445 * t695 - t446 * t690) * t615 + t511 * t484 + (-t691 * t793 + (-t690 * t696 - t792) * t615 * qJD(5)) * pkin(4) + ((-pkin(4) * t792 - t673 * t690) * t615 - t430) * qJD(6) + t837) * MDP(32) + t852; (t450 * t619 + t835) * MDP(25) + (t449 * t619 + t836) * MDP(26) + (-(-t443 * t690 - t786) * t615 - t430 * qJD(6) + (-t484 * t723 - t615 * t763 + t785) * pkin(5) + t837) * MDP(32) + ((-t444 * t615 - t425) * t690 + (t443 * t615 - t736) * t695 + (-t615 * t762 - t723 * t724 - t793) * pkin(5) + t838) * MDP(33) + t852; (t755 - t843) * MDP(29) + (-t744 - t842) * MDP(30) + (t430 * t615 + t837) * MDP(32) + (-t690 * t425 - t695 * t426 + t429 * t615 + t838) * MDP(33) + (MDP(29) * t813 + MDP(30) * t724 - MDP(32) * t430 - MDP(33) * t819) * qJD(6) + t854;];
tau  = t1;
