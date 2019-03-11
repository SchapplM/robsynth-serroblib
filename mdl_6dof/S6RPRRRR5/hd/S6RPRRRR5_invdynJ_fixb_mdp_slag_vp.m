% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:10:12
% EndTime: 2019-03-09 07:10:30
% DurationCPUTime: 13.39s
% Computational Cost: add. (11062->594), mult. (26777->746), div. (0->0), fcn. (22233->18), ass. (0->278)
t702 = cos(pkin(11));
t710 = cos(qJ(3));
t811 = t710 * t702;
t701 = sin(pkin(11));
t706 = sin(qJ(3));
t817 = t701 * t706;
t641 = -t811 + t817;
t629 = t641 * qJD(1);
t642 = t701 * t710 + t702 * t706;
t630 = t642 * qJD(1);
t705 = sin(qJ(4));
t853 = cos(qJ(4));
t596 = t853 * t629 + t630 * t705;
t858 = qJD(5) + qJD(6);
t895 = t596 + t858;
t703 = sin(qJ(6));
t704 = sin(qJ(5));
t708 = cos(qJ(6));
t709 = cos(qJ(5));
t647 = t703 * t704 - t708 * t709;
t891 = t895 * t647;
t786 = qJD(1) * qJD(3);
t770 = t706 * t786;
t784 = qJDD(1) * t710;
t869 = qJDD(1) * t706 + t710 * t786;
t775 = t701 * t784 + t702 * t869;
t603 = -t701 * t770 + t775;
t776 = -t701 * t869 - t702 * t770;
t726 = t702 * t784 + t776;
t733 = -t705 * t629 + t630 * t853;
t527 = qJD(4) * t733 + t705 * t603 - t853 * t726;
t524 = qJDD(5) + t527;
t523 = qJDD(6) + t524;
t867 = -qJD(5) - t596;
t591 = qJD(6) - t867;
t816 = t703 * t709;
t648 = t704 * t708 + t816;
t882 = t648 * t523 - t591 * t891;
t894 = t895 * t648;
t793 = qJD(5) * t709;
t875 = t596 * t709;
t893 = t793 + t875;
t794 = qJD(5) * t704;
t876 = t596 * t704;
t892 = t794 + t876;
t698 = pkin(11) + qJ(3);
t691 = qJ(4) + t698;
t680 = sin(t691);
t707 = sin(qJ(1));
t711 = cos(qJ(1));
t750 = g(1) * t711 + g(2) * t707;
t890 = t750 * t680;
t889 = pkin(10) * t876;
t845 = pkin(7) + qJ(2);
t662 = t845 * t701;
t643 = qJD(1) * t662;
t663 = t845 * t702;
t644 = qJD(1) * t663;
t737 = t643 * t706 - t644 * t710;
t583 = -pkin(8) * t629 - t737;
t577 = t705 * t583;
t861 = -t710 * t643 - t644 * t706;
t582 = -pkin(8) * t630 + t861;
t538 = t582 * t853 - t577;
t771 = qJD(4) * t853;
t870 = -pkin(3) * t771 + t538;
t888 = t892 * pkin(5);
t887 = pkin(5) * t733 + pkin(10) * t875;
t695 = qJDD(3) + qJDD(4);
t787 = qJD(1) * qJD(2);
t857 = qJDD(1) * t845 + t787;
t618 = t857 * t701;
t619 = t857 * t702;
t761 = -t710 * t618 - t706 * t619;
t521 = qJDD(3) * pkin(3) - pkin(8) * t603 + qJD(3) * t737 + t761;
t738 = -t706 * t618 + t710 * t619;
t525 = t726 * pkin(8) + qJD(3) * t861 + t738;
t579 = qJD(3) * pkin(3) + t582;
t795 = qJD(4) * t705;
t754 = -t853 * t521 + t705 * t525 + t579 * t795 + t583 * t771;
t470 = -pkin(4) * t695 + t754;
t699 = qJD(3) + qJD(4);
t586 = t699 * t704 + t709 * t733;
t526 = t853 * t603 - t629 * t771 - t630 * t795 + t705 * t726;
t763 = t526 * t704 - t709 * t695;
t501 = qJD(5) * t586 + t763;
t461 = pkin(5) * t501 + t470;
t534 = t579 * t853 - t577;
t528 = -t699 * pkin(4) - t534;
t584 = -t709 * t699 + t704 * t733;
t505 = t584 * pkin(5) + t528;
t681 = cos(t691);
t700 = qJ(5) + qJ(6);
t692 = sin(t700);
t848 = g(3) * t692;
t886 = t461 * t648 - t505 * t891 + t681 * t848 - t692 * t890;
t693 = cos(t700);
t847 = g(3) * t693;
t885 = t461 * t647 + t505 * t894 - t681 * t847 + t693 * t890;
t561 = pkin(4) * t733 + pkin(9) * t596;
t546 = pkin(3) * t630 + t561;
t884 = -t709 * t546 + t704 * t870;
t500 = t709 * t526 + t704 * t695 + t699 * t793 - t733 * t794;
t883 = t500 * t709 - t704 * t501 - t584 * t893;
t748 = -t647 * t523 - t591 * t894;
t516 = t704 * t524;
t881 = -t867 * t893 + t516;
t791 = qJD(6) * t708;
t779 = t708 * t500 - t703 * t501 - t584 * t791;
t792 = qJD(6) * t703;
t465 = -t586 * t792 + t779;
t739 = t584 * t703 - t708 * t586;
t764 = t500 * t703 + t708 * t501;
t466 = -qJD(6) * t739 + t764;
t498 = t500 * t704;
t832 = t586 * t703;
t540 = t708 * t584 + t832;
t880 = t695 * MDP(19) + (t586 * t893 + t498) * MDP(22) + t465 * t648 * MDP(29) + (-t465 * t647 - t648 * t466 + t540 * t891) * MDP(30) + (t891 * MDP(29) + MDP(30) * t894) * t739;
t879 = t528 * t596;
t878 = t540 * t591;
t877 = t591 * t739;
t849 = g(3) * t681;
t872 = t470 + t849;
t868 = MDP(4) * t702 - MDP(5) * t701;
t578 = t853 * t583;
t535 = t705 * t579 + t578;
t529 = pkin(9) * t699 + t535;
t682 = -pkin(2) * t702 - pkin(1);
t656 = qJD(1) * t682 + qJD(2);
t608 = pkin(3) * t629 + t656;
t536 = pkin(4) * t596 - pkin(9) * t733 + t608;
t487 = t529 * t709 + t536 * t704;
t481 = -pkin(10) * t584 + t487;
t476 = t481 * t792;
t819 = t693 * t707;
t820 = t692 * t711;
t613 = -t681 * t819 + t820;
t818 = t693 * t711;
t821 = t692 * t707;
t615 = t681 * t818 + t821;
t866 = g(1) * t615 - g(2) * t613 + t505 * t540 + t680 * t847 + t476;
t612 = t681 * t821 + t818;
t614 = -t681 * t820 + t819;
t716 = t705 * t521 + t525 * t853 + t579 * t771 - t583 * t795;
t469 = t695 * pkin(9) + t716;
t587 = qJDD(2) - t776 * pkin(3) + (-pkin(1) + (-pkin(3) * t710 - pkin(2)) * t702) * qJDD(1);
t479 = t527 * pkin(4) - t526 * pkin(9) + t587;
t478 = t709 * t479;
t453 = pkin(5) * t524 - pkin(10) * t500 - qJD(5) * t487 - t469 * t704 + t478;
t728 = t709 * t469 + t704 * t479 - t529 * t794 + t536 * t793;
t454 = -pkin(10) * t501 + t728;
t765 = t708 * t453 - t703 * t454;
t865 = -g(1) * t614 + g(2) * t612 + t505 * t739 + t680 * t848 + t765;
t864 = t523 * MDP(33) + (-t540 ^ 2 + t739 ^ 2) * MDP(30) - t540 * MDP(29) * t739;
t605 = -t705 * t641 + t642 * t853;
t562 = t648 * t605;
t537 = t705 * t582 + t578;
t752 = pkin(3) * t795 - t537;
t802 = -t706 * t662 + t710 * t663;
t860 = t704 * t546 + t709 * t870;
t859 = qJ(2) * qJDD(1);
t856 = -t733 * MDP(15) + MDP(16) * t596 - t608 * MDP(21);
t486 = -t529 * t704 + t709 * t536;
t480 = -pkin(10) * t586 + t486;
t472 = -pkin(5) * t867 + t480;
t840 = t472 * t708;
t457 = -t481 * t703 + t840;
t839 = t481 * t708;
t458 = t472 * t703 + t839;
t855 = MDP(16) * t733 - MDP(20) * t608 + MDP(26) * t867 - MDP(27) * t486 + MDP(28) * t487 - MDP(33) * t591 - MDP(34) * t457 + MDP(35) * t458;
t854 = -pkin(9) - pkin(10);
t632 = t642 * qJD(3);
t852 = pkin(3) * t632;
t675 = g(3) * t680;
t846 = t709 * pkin(5);
t685 = pkin(3) * t705 + pkin(9);
t844 = -pkin(10) - t685;
t841 = qJDD(1) * pkin(1);
t838 = t540 * t733;
t837 = t739 * t733;
t631 = t641 * qJD(3);
t732 = -t641 * t853 - t705 * t642;
t564 = qJD(4) * t732 - t631 * t853 - t705 * t632;
t836 = t564 * t704;
t835 = t564 * t709;
t834 = t584 * t733;
t833 = t586 * t733;
t831 = t586 * t704;
t830 = t733 * t699;
t827 = t596 * t699;
t824 = t605 * t704;
t823 = t605 * t709;
t815 = t704 * t707;
t814 = t704 * t711;
t813 = t707 * t709;
t517 = t709 * t524;
t639 = t710 * t662;
t760 = -t663 * t706 - t639;
t589 = -pkin(8) * t642 + t760;
t590 = -pkin(8) * t641 + t802;
t552 = t705 * t589 + t590 * t853;
t548 = t709 * t552;
t812 = t709 * t711;
t805 = t709 * t534 + t704 * t561;
t616 = pkin(3) * t641 + t682;
t553 = -pkin(4) * t732 - pkin(9) * t605 + t616;
t803 = t704 * t553 + t548;
t801 = t752 + t888;
t800 = t701 ^ 2 + t702 ^ 2;
t796 = qJD(3) * t630;
t781 = qJD(5) * pkin(9) * t867;
t519 = t528 * t793;
t780 = t704 * t872 + t519;
t777 = t528 * t794 + t709 * t890;
t774 = qJD(5) * t854;
t773 = qJD(1) * t817;
t772 = t605 * t794;
t767 = qJD(5) * t844;
t766 = t800 * qJD(1) ^ 2;
t759 = t867 * t704;
t758 = -qJD(5) * t536 - t469;
t757 = qJD(6) * t472 + t454;
t753 = 0.2e1 * t800;
t686 = -pkin(3) * t853 - pkin(4);
t751 = -t535 + t888;
t749 = g(1) * t707 - g(2) * t711;
t747 = -t529 * t793 + t478;
t694 = t709 * pkin(10);
t637 = t685 * t709 + t694;
t746 = qJD(6) * t637 - t709 * t767 - t884 + t887;
t559 = t709 * t561;
t666 = pkin(9) * t709 + t694;
t745 = qJD(6) * t666 - t534 * t704 - t709 * t774 + t559 + t887;
t636 = t844 * t704;
t744 = -qJD(6) * t636 - t704 * t767 + t860 + t889;
t665 = t854 * t704;
t743 = -qJD(6) * t665 - t704 * t774 + t805 + t889;
t742 = -pkin(9) * t524 + t879;
t740 = -t524 * t685 + t879;
t736 = t867 * t892 + t517;
t734 = t589 * t853 - t705 * t590;
t731 = t605 * t793 + t836;
t730 = -t772 + t835;
t720 = -qJD(3) * t639 + qJD(2) * t811 + (-qJD(2) * t701 - qJD(3) * t663) * t706;
t568 = -pkin(8) * t632 + t720;
t715 = -t642 * qJD(2) - qJD(3) * t802;
t569 = pkin(8) * t631 + t715;
t490 = qJD(4) * t734 + t568 * t853 + t705 * t569;
t565 = qJD(4) * t605 - t705 * t631 + t632 * t853;
t504 = pkin(4) * t565 - pkin(9) * t564 + t852;
t727 = t709 * t490 + t704 * t504 - t552 * t794 + t553 * t793;
t725 = -t749 - t841;
t722 = -t754 - t849 + t890;
t718 = t753 * t787 - t750;
t491 = qJD(4) * t552 + t705 * t568 - t569 * t853;
t714 = t681 * t750 + t675 - t716;
t690 = cos(t698);
t689 = sin(t698);
t688 = qJDD(2) - t841;
t687 = -pkin(4) - t846;
t664 = t686 - t846;
t655 = qJDD(1) * t682 + qJDD(2);
t626 = t681 * t812 + t815;
t625 = -t681 * t814 + t813;
t624 = -t681 * t813 + t814;
t623 = t681 * t815 + t812;
t563 = t647 * t605;
t550 = t709 * t553;
t513 = pkin(5) * t824 - t734;
t503 = t709 * t504;
t489 = -pkin(10) * t824 + t803;
t483 = -pkin(5) * t732 - pkin(10) * t823 - t552 * t704 + t550;
t474 = t564 * t816 - t703 * t772 - t792 * t824 + (t823 * t858 + t836) * t708;
t473 = -t562 * t858 - t647 * t564;
t471 = pkin(5) * t731 + t491;
t456 = -pkin(10) * t731 + t727;
t455 = -pkin(10) * t835 + pkin(5) * t565 - t490 * t704 + t503 + (-t548 + (pkin(10) * t605 - t553) * t704) * qJD(5);
t1 = [(-t603 * t641 + t631 * t629 - t630 * t632 + t642 * t726) * MDP(9) + (t526 * t732 - t527 * t605 - t564 * t596 - t565 * t733) * MDP(16) + ((t455 * t708 - t456 * t703) * t591 + (t483 * t708 - t489 * t703) * t523 - t765 * t732 + t457 * t565 + t471 * t540 + t513 * t466 + t461 * t562 + t505 * t474 - g(1) * t613 - g(2) * t615 + ((-t483 * t703 - t489 * t708) * t591 + t458 * t732) * qJD(6)) * MDP(34) + (-t565 * t699 + t695 * t732) * MDP(18) + (-t523 * t732 + t565 * t591) * MDP(33) + (t466 * t732 - t474 * t591 - t523 * t562 - t540 * t565) * MDP(32) + (-t465 * t563 - t473 * t739) * MDP(29) + (-t465 * t562 + t466 * t563 - t473 * t540 + t474 * t739) * MDP(30) + (-g(1) * t612 - g(2) * t614 - t458 * t565 - t461 * t563 + t513 * t465 - t471 * t739 + t505 * t473 - t476 * t732 + (-(-qJD(6) * t489 + t455) * t591 - t483 * t523 + t453 * t732) * t703 + (-(qJD(6) * t483 + t456) * t591 - t489 * t523 + t757 * t732) * t708) * MDP(35) + (-t465 * t732 + t473 * t591 - t523 * t563 - t565 * t739) * MDP(31) + (-t491 * t699 + t527 * t616 + t565 * t608 - t587 * t732 + t596 * t852 + t681 * t749 + t695 * t734) * MDP(20) + t868 * (-t688 - t725) + ((-t584 * t709 - t831) * t564 + (-t498 - t501 * t709 + (t584 * t704 - t586 * t709) * qJD(5)) * t605) * MDP(23) + (-t524 * t732 - t565 * t867) * MDP(26) + (-g(1) * t623 - g(2) * t625 + t470 * t823 - t487 * t565 + t491 * t586 - t500 * t734 - t524 * t803 + t528 * t730 + t727 * t867 + t728 * t732) * MDP(28) + (-(-t552 * t793 + t503) * t867 + t550 * t524 - t747 * t732 + t486 * t565 + t491 * t584 - t734 * t501 + t605 * t519 - g(1) * t624 - g(2) * t626 + (-(-qJD(5) * t553 - t490) * t867 - t552 * t524 - t758 * t732 + t470 * t605 + t528 * t564) * t704) * MDP(27) + (-t500 * t732 + t517 * t605 + t565 * t586 - t730 * t867) * MDP(24) + (t501 * t732 - t516 * t605 - t565 * t584 + t731 * t867) * MDP(25) + (-qJD(3) * t720 - qJDD(3) * t802 + t682 * t603 - t656 * t631 + t655 * t642 - t689 * t749) * MDP(14) + (-qJD(3) * t631 + qJDD(3) * t642) * MDP(10) + (t603 * t642 - t630 * t631) * MDP(8) + (-t490 * t699 + t526 * t616 - t552 * t695 + t564 * t608 + t587 * t605 - t680 * t749 + t733 * t852) * MDP(21) + (t526 * t605 + t564 * t733) * MDP(15) + (t500 * t823 + t586 * t730) * MDP(22) + t749 * MDP(2) + t750 * MDP(3) + (qJD(3) * t715 + qJDD(3) * t760 + t656 * t632 + t655 * t641 - t682 * t726 + t690 * t749) * MDP(13) + ((-t688 + t749) * pkin(1) + (t800 * t859 + t718) * qJ(2)) * MDP(7) + (t753 * t859 + t718) * MDP(6) + qJDD(1) * MDP(1) + (t564 * t699 + t605 * t695) * MDP(17) + (-qJD(3) * t632 - qJDD(3) * t641) * MDP(11); -MDP(6) * t766 + (-qJ(2) * t766 + qJDD(2) + t725) * MDP(7) + (-t726 + t796) * MDP(13) + ((-t629 - t773) * qJD(3) + t775) * MDP(14) + (t527 + t830) * MDP(20) + (t526 - t827) * MDP(21) + (t736 - t834) * MDP(27) + (-t709 * t867 ^ 2 - t516 - t833) * MDP(28) + (t748 - t838) * MDP(34) + (t837 - t882) * MDP(35) - t868 * qJDD(1); t855 * t733 + (t748 + t838) * MDP(32) + t880 + (t736 + t834) * MDP(25) + ((t636 * t708 - t637 * t703) * t523 + t664 * t466 + (t703 * t744 - t708 * t746) * t591 + t801 * t540 + t885) * MDP(34) + (-(t636 * t703 + t637 * t708) * t523 + t664 * t465 + (t703 * t746 + t708 * t744) * t591 - t801 * t739 + t886) * MDP(35) + (t837 + t882) * MDP(31) + (-t527 + t830) * MDP(18) + (t686 * t500 + t740 * t709 - t704 * t890 + t752 * t586 - (t685 * t794 + t860) * t867 + t780) * MDP(28) + (t831 * t867 + t883) * MDP(23) + (-t833 + t881) * MDP(24) + (t686 * t501 - t872 * t709 + t740 * t704 + t752 * t584 - (-t685 * t793 + t884) * t867 + t777) * MDP(27) + (t538 * t699 + (-t630 * t733 - t695 * t705 - t699 * t771) * pkin(3) + t714) * MDP(21) + (t537 * t699 + (-t596 * t630 + t695 * t853 - t699 * t795) * pkin(3) + t722) * MDP(20) + (t726 + t796) * MDP(11) + ((t629 - t773) * qJD(3) + t775) * MDP(10) + (g(3) * t689 + t656 * t629 + t690 * t750 - t738) * MDP(14) + (-g(3) * t690 - t656 * t630 + t689 * t750 + t761) * MDP(13) + (t526 + t827) * MDP(17) + qJDD(3) * MDP(12) + t630 * t629 * MDP(8) - t856 * t596 + (-t629 ^ 2 + t630 ^ 2) * MDP(9); t526 * MDP(17) - t527 * MDP(18) + (t535 * t699 + t722) * MDP(20) + (t534 * t699 + t714) * MDP(21) + (t586 * t759 + t883) * MDP(23) + t881 * MDP(24) + (-t759 * t867 + t517) * MDP(25) + (-pkin(4) * t501 - t535 * t584 + t559 * t867 + (-t534 * t867 + t742) * t704 + (-t872 + t781) * t709 + t777) * MDP(27) + (-pkin(4) * t500 - t805 * t867 - t535 * t586 + t742 * t709 + (-t890 - t781) * t704 + t780) * MDP(28) + t882 * MDP(31) + t748 * MDP(32) + ((t665 * t708 - t666 * t703) * t523 + t687 * t466 + (t703 * t743 - t708 * t745) * t591 + t751 * t540 + t885) * MDP(34) + (-(t665 * t703 + t666 * t708) * t523 + t687 * t465 + (t703 * t745 + t708 * t743) * t591 - t751 * t739 + t886) * MDP(35) + (t699 * MDP(17) - t856) * t596 + (MDP(18) * t699 - MDP(24) * t586 + MDP(25) * t584 + MDP(31) * t739 + MDP(32) * t540 + t855) * t733 + t880; t586 * t584 * MDP(22) + (-t584 ^ 2 + t586 ^ 2) * MDP(23) + (-t584 * t867 + t500) * MDP(24) + (-t763 + (-qJD(5) - t867) * t586) * MDP(25) + t524 * MDP(26) + (-g(1) * t625 + g(2) * t623 - t487 * t867 - t528 * t586 + (t758 + t675) * t704 + t747) * MDP(27) + (g(1) * t626 - g(2) * t624 - t486 * t867 + t528 * t584 + t675 * t709 - t728) * MDP(28) + (t465 + t878) * MDP(31) + (-t466 - t877) * MDP(32) + (-(-t480 * t703 - t839) * t591 - t458 * qJD(6) + (t523 * t708 - t540 * t586 - t591 * t792) * pkin(5) + t865) * MDP(34) + ((-t481 * t591 - t453) * t703 + (t480 * t591 - t757) * t708 + (-t523 * t703 + t586 * t739 - t591 * t791) * pkin(5) + t866) * MDP(35) + t864; (t779 + t878) * MDP(31) + (-t764 - t877) * MDP(32) + (t458 * t591 + t865) * MDP(34) + (-t703 * t453 - t708 * t454 + t457 * t591 + t866) * MDP(35) + (-MDP(31) * t832 + MDP(32) * t739 - MDP(34) * t458 - MDP(35) * t840) * qJD(6) + t864;];
tau  = t1;
