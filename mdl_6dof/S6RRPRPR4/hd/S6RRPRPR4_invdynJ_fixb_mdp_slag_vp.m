% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:39
% EndTime: 2019-03-09 10:26:56
% DurationCPUTime: 13.04s
% Computational Cost: add. (13269->661), mult. (37664->902), div. (0->0), fcn. (31290->16), ass. (0->318)
t706 = sin(pkin(11));
t707 = sin(pkin(6));
t714 = sin(qJ(2));
t815 = qJD(1) * t714;
t789 = t707 * t815;
t709 = cos(pkin(11));
t718 = cos(qJ(2));
t828 = t718 * t709;
t793 = t707 * t828;
t644 = qJD(1) * t793 - t706 * t789;
t638 = qJD(4) - t644;
t710 = cos(pkin(6));
t834 = t710 * t714;
t687 = pkin(1) * t834;
t837 = t707 * t718;
t863 = pkin(8) + qJ(3);
t630 = (t837 * t863 + t687) * qJD(1);
t622 = t706 * t630;
t870 = pkin(1) * t710;
t688 = t718 * t870;
t682 = qJD(1) * t688;
t786 = t863 * t714;
t761 = t707 * t786;
t629 = -qJD(1) * t761 + t682;
t563 = t629 * t709 - t622;
t745 = t706 * t718 + t709 * t714;
t735 = qJD(1) * t745;
t647 = t707 * t735;
t583 = pkin(2) * t789 + pkin(3) * t647 - pkin(9) * t644;
t717 = cos(qJ(4));
t567 = t717 * t583;
t713 = sin(qJ(4));
t694 = pkin(2) * t706 + pkin(9);
t826 = qJ(5) + t694;
t769 = qJD(4) * t826;
t890 = -pkin(4) * t647 - t567 + (qJ(5) * t644 - t769) * t717 + (-qJD(5) + t563) * t713;
t824 = t717 * t563 + t713 * t583;
t847 = t644 * t713;
t889 = -qJ(5) * t847 - qJD(5) * t717 + t713 * t769 + t824;
t812 = qJD(4) * t713;
t888 = t812 - t847;
t716 = cos(qJ(6));
t712 = sin(qJ(6));
t816 = qJD(1) * t710;
t685 = qJD(2) + t816;
t612 = t647 * t713 - t717 * t685;
t705 = sin(pkin(12));
t708 = cos(pkin(12));
t746 = -t647 * t717 - t685 * t713;
t747 = -t612 * t705 - t708 * t746;
t855 = t747 * t712;
t524 = -t716 * t638 + t855;
t771 = -t708 * t612 + t705 * t746;
t875 = qJD(6) - t771;
t887 = t524 * t875;
t526 = t638 * t712 + t716 * t747;
t886 = t526 * t875;
t664 = t705 * t717 + t708 * t713;
t823 = t638 * t664;
t662 = t705 * t713 - t708 * t717;
t878 = t638 * t662;
t767 = t875 * t716;
t814 = qJD(2) * t714;
t788 = t707 * t814;
t760 = qJD(1) * t788;
t806 = qJD(1) * qJD(2);
t785 = t718 * t806;
t597 = -t706 * t760 + (qJDD(1) * t745 + t709 * t785) * t707;
t805 = qJDD(1) * t710;
t684 = qJDD(2) + t805;
t811 = qJD(4) * t717;
t530 = t717 * t597 - t647 * t812 + t713 * t684 + t685 * t811;
t531 = -qJD(4) * t746 + t597 * t713 - t717 * t684;
t493 = -t530 * t705 - t708 * t531;
t492 = qJDD(6) - t493;
t833 = t712 * t492;
t885 = -t875 * t767 - t833;
t494 = t530 * t708 - t531 * t705;
t803 = qJDD(1) * t718;
t783 = t707 * t803;
t668 = t709 * t783;
t734 = t745 * qJD(2);
t804 = qJDD(1) * t714;
t784 = t706 * t804;
t595 = qJDD(4) - t668 + (qJD(1) * t734 + t784) * t707;
t778 = t494 * t712 - t716 * t595;
t472 = qJD(6) * t526 + t778;
t884 = t472 * t662 + t524 * t823;
t651 = t745 * t710;
t663 = t706 * t714 - t828;
t715 = sin(qJ(1));
t719 = cos(qJ(1));
t604 = t651 * t719 - t663 * t715;
t702 = qJ(4) + pkin(12);
t699 = sin(t702);
t700 = cos(t702);
t836 = t707 * t719;
t578 = -t604 * t700 + t699 * t836;
t738 = t710 * t663;
t603 = -t715 * t745 - t719 * t738;
t883 = t578 * t712 - t603 * t716;
t882 = t578 * t716 + t603 * t712;
t698 = pkin(2) * t718 + pkin(1);
t758 = t698 * t707;
t802 = pkin(2) * t760 + qJDD(3);
t881 = qJDD(1) * t758 - t802;
t628 = pkin(2) * t710 + t688 - t761;
t818 = pkin(8) * t837 + t687;
t639 = qJ(3) * t837 + t818;
t573 = t706 * t628 + t709 * t639;
t561 = pkin(9) * t710 + t573;
t839 = t707 * t714;
t649 = t706 * t839 - t793;
t650 = t745 * t707;
t585 = pkin(3) * t649 - pkin(9) * t650 - t758;
t825 = t717 * t561 + t713 * t585;
t821 = t889 * t705 + t708 * t890;
t819 = t705 * t890 - t889 * t708;
t625 = t650 * t713 - t710 * t717;
t626 = t650 * t717 + t710 * t713;
t554 = -t625 * t705 + t626 * t708;
t642 = t649 * t716;
t877 = -t554 * t712 + t642;
t835 = t709 * t630;
t562 = t629 * t706 + t835;
t876 = pkin(4) * t888 - t562;
t605 = t715 * t651 + t663 * t719;
t838 = t707 * t715;
t586 = t605 * t713 + t717 * t838;
t740 = t604 * t713 + t717 * t836;
t874 = -g(1) * t586 + g(2) * t740 + g(3) * t625;
t615 = pkin(2) * t685 + t629;
t549 = t706 * t615 + t835;
t544 = pkin(9) * t685 + t549;
t654 = -qJD(1) * t758 + qJD(3);
t559 = -pkin(3) * t644 - pkin(9) * t647 + t654;
t516 = t544 * t717 + t559 * t713;
t798 = pkin(1) * t803;
t681 = t710 * t798;
t800 = qJD(2) * t870;
t762 = qJD(1) * t800;
t781 = qJD(2) * t863;
t813 = qJD(3) * t714;
t547 = -t714 * t762 + pkin(2) * t684 + t681 + (-qJDD(1) * t786 + (-t718 * t781 - t813) * qJD(1)) * t707;
t728 = qJD(3) * t718 - t714 * t781;
t790 = pkin(8) * t783 + qJDD(1) * t687 + t718 * t762;
t556 = (qJ(3) * t803 + qJD(1) * t728) * t707 + t790;
t514 = t706 * t547 + t709 * t556;
t508 = pkin(9) * t684 + t514;
t596 = t668 + (-qJD(2) * t735 - t784) * t707;
t522 = -pkin(3) * t596 - pkin(9) * t597 - t881;
t777 = -t713 * t508 + t717 * t522;
t722 = -qJD(4) * t516 + t777;
t458 = pkin(4) * t595 - qJ(5) * t530 + qJD(5) * t746 + t722;
t737 = -t717 * t508 - t713 * t522 + t544 * t812 - t559 * t811;
t460 = -qJ(5) * t531 - qJD(5) * t612 - t737;
t447 = t458 * t708 - t460 * t705;
t445 = -pkin(5) * t595 - t447;
t693 = pkin(4) * t705 + pkin(10);
t873 = t875 * (-pkin(4) * t746 + pkin(5) * t747 - pkin(10) * t771 + qJD(6) * t693) + g(1) * (t605 * t699 + t700 * t838) + g(2) * (-t604 * t699 - t700 * t836) + g(3) * (-t650 * t699 + t700 * t710) + t445;
t701 = t707 ^ 2;
t872 = 0.2e1 * t701;
t871 = pkin(1) * t701;
t869 = pkin(2) * t709;
t829 = t715 * t718;
t830 = t714 * t719;
t659 = -t710 * t829 - t830;
t867 = g(1) * t659;
t866 = g(1) * t715;
t864 = g(3) * t718;
t862 = MDP(6) * t707;
t861 = MDP(7) * t707;
t809 = qJD(6) * t716;
t791 = t716 * t494 + t712 * t595 + t638 * t809;
t810 = qJD(6) * t712;
t471 = -t747 * t810 + t791;
t860 = t471 * t712;
t500 = -qJ(5) * t612 + t516;
t858 = t500 * t705;
t857 = t524 * t747;
t856 = t526 * t747;
t851 = t612 * t638;
t850 = t612 * t647;
t849 = t746 * t638;
t848 = t746 * t647;
t846 = t649 * t712;
t845 = t664 * t716;
t844 = t684 * MDP(8);
t843 = t698 * t715;
t842 = t700 * t712;
t841 = t700 * t716;
t720 = qJD(1) ^ 2;
t840 = t701 * t720;
t497 = t708 * t500;
t832 = t713 * t595;
t831 = t714 * t715;
t490 = t716 * t492;
t827 = t718 * t719;
t448 = t705 * t458 + t708 * t460;
t646 = t663 * t707 * qJD(2);
t571 = -qJD(4) * t625 - t646 * t717;
t645 = t707 * t734;
t683 = t718 * t800;
t616 = t707 * t728 + t683;
t787 = t863 * t707;
t617 = -t707 * t813 + (-t718 * t787 - t687) * qJD(2);
t546 = t616 * t709 + t617 * t706;
t584 = pkin(2) * t788 + pkin(3) * t645 + pkin(9) * t646;
t774 = -t546 * t713 + t717 * t584;
t470 = pkin(4) * t645 - qJ(5) * t571 - qJD(4) * t825 - qJD(5) * t626 + t774;
t570 = qJD(4) * t626 - t646 * t713;
t736 = t717 * t546 - t561 * t812 + t713 * t584 + t585 * t811;
t478 = -qJ(5) * t570 - qJD(5) * t625 + t736;
t452 = t705 * t470 + t708 * t478;
t515 = -t544 * t713 + t717 * t559;
t499 = qJ(5) * t746 + t515;
t489 = pkin(4) * t638 + t499;
t464 = t705 * t489 + t497;
t773 = -t561 * t713 + t717 * t585;
t503 = pkin(4) * t649 - qJ(5) * t626 + t773;
t511 = -qJ(5) * t625 + t825;
t480 = t705 * t503 + t708 * t511;
t820 = pkin(5) * t647 - t821;
t703 = t714 ^ 2;
t817 = -t718 ^ 2 + t703;
t808 = qJD(2) - t685;
t801 = t714 * t871;
t797 = pkin(8) * t804;
t796 = t664 * t833;
t795 = t664 * t490;
t794 = t718 * t840;
t792 = t710 * t827;
t697 = pkin(4) * t717 + pkin(3);
t782 = g(2) * t836 - g(3) * t710;
t780 = t826 * t713;
t446 = pkin(10) * t595 + t448;
t513 = t547 * t709 - t706 * t556;
t507 = -pkin(3) * t684 - t513;
t484 = pkin(4) * t531 + qJDD(5) + t507;
t456 = -pkin(5) * t493 - pkin(10) * t494 + t484;
t779 = -t712 * t446 + t716 * t456;
t776 = -t716 * t647 + t712 * t878;
t775 = t647 * t712 + t716 * t878;
t770 = t604 * t717 - t713 * t836;
t545 = t616 * t706 - t709 * t617;
t548 = t615 * t709 - t622;
t572 = t628 * t709 - t706 * t639;
t768 = t638 * t717;
t766 = t685 + t816;
t764 = qJD(1) * t808;
t763 = t684 + t805;
t653 = pkin(2) * t834 - t787;
t759 = pkin(4) * t707 * t713 - t653;
t757 = g(1) * t719 + g(2) * t715;
t756 = t471 * t662 + t526 * t823;
t744 = -t697 - t869;
t601 = pkin(5) * t662 - pkin(10) * t664 + t744;
t754 = pkin(10) * t647 - qJD(6) * t601 - t819;
t661 = t826 * t717;
t611 = t708 * t661 - t705 * t780;
t753 = -pkin(5) * t823 - pkin(10) * t878 + qJD(6) * t611 - t876;
t752 = t716 * t446 + t712 * t456;
t462 = pkin(10) * t638 + t464;
t543 = -pkin(3) * t685 - t548;
t523 = pkin(4) * t612 + qJD(5) + t543;
t483 = -pkin(5) * t771 - pkin(10) * t747 + t523;
t454 = t462 * t716 + t483 * t712;
t751 = t462 * t712 - t483 * t716;
t451 = t470 * t708 - t478 * t705;
t475 = pkin(10) * t649 + t480;
t553 = t708 * t625 + t626 * t705;
t560 = -pkin(3) * t710 - t572;
t724 = pkin(4) * t625 + t560;
t487 = pkin(5) * t553 - pkin(10) * t554 + t724;
t750 = t475 * t716 + t487 * t712;
t749 = -t475 * t712 + t487 * t716;
t463 = t489 * t708 - t858;
t479 = t503 * t708 - t511 * t705;
t533 = t554 * t716 + t846;
t743 = t717 * t595 - t638 * t888;
t742 = t490 + (t712 * t771 - t810) * t875;
t741 = pkin(4) * t570 + t545;
t733 = t543 * t638 - t694 * t595;
t606 = t715 * t738 - t719 * t745;
t732 = g(1) * t606 + g(2) * t603 - g(3) * t649;
t731 = t664 * t809 - t776;
t730 = t664 * t810 + t775;
t725 = t818 * t685;
t461 = -pkin(5) * t638 - t463;
t468 = t499 * t708 - t858;
t723 = -t693 * t492 + (t461 + t468) * t875;
t721 = qJD(4) * t638 * t694 + t507 + t732;
t711 = -qJ(5) - pkin(9);
t696 = -pkin(3) - t869;
t695 = -pkin(4) * t708 - pkin(5);
t677 = t719 * t698;
t675 = pkin(2) * t792;
t660 = -t710 * t831 + t827;
t658 = -t710 * t830 - t829;
t657 = -t792 + t831;
t619 = t650 * t700 + t699 * t710;
t610 = t661 * t705 + t708 * t780;
t587 = -t605 * t717 + t713 * t838;
t580 = -t605 * t700 + t699 * t838;
t528 = t580 * t716 - t606 * t712;
t527 = -t580 * t712 - t606 * t716;
t521 = -t570 * t705 + t571 * t708;
t520 = t708 * t570 + t571 * t705;
t486 = qJD(6) * t533 + t521 * t712 - t645 * t716;
t485 = qJD(6) * t877 + t521 * t716 + t645 * t712;
t474 = -pkin(5) * t649 - t479;
t467 = t499 * t705 + t497;
t466 = pkin(5) * t520 - pkin(10) * t521 + t741;
t450 = pkin(10) * t645 + t452;
t449 = -pkin(5) * t645 - t451;
t444 = -qJD(6) * t454 + t779;
t443 = -qJD(6) * t751 + t752;
t1 = [(-t513 * t650 - t514 * t649 + t545 * t647 + t546 * t644 + t548 * t646 - t549 * t645 - t572 * t597 + t573 * t596 - t707 * t757) * MDP(11) + (t798 * t872 + (-pkin(8) * t839 + t688) * t684 + (-t707 * t797 + t681) * t710 - g(1) * t658 - g(2) * t660 + (-t725 + (-t710 * t818 - 0.2e1 * t801) * qJD(1)) * qJD(2)) * MDP(9) + (t492 * t553 + t520 * t875) * MDP(26) + (t471 * t553 + t485 * t875 + t492 * t533 + t520 * t526) * MDP(24) + (t448 * t480 + t464 * t452 + t447 * t479 + t463 * t451 + t484 * t724 + t523 * t741 - g(1) * (-t603 * t711 - t604 * t697 + t719 * t759 - t843) - g(2) * (-t605 * t697 + t606 * t711 + t715 * t759 + t677)) * MDP(21) + (-t530 * t625 - t531 * t626 + t570 * t746 - t571 * t612) * MDP(14) + (t530 * t626 - t571 * t746) * MDP(13) + (-g(1) * t740 - g(2) * t586 + t507 * t626 - t516 * t645 + t560 * t530 + t543 * t571 - t545 * t746 - t595 * t825 - t638 * t736 + t649 * t737) * MDP(19) + (t530 * t649 + t571 * t638 + t595 * t626 - t645 * t746) * MDP(15) + (-(-pkin(8) * t788 + t683) * t685 - t818 * t684 - (-pkin(8) * t760 + t790) * t710 - g(1) * t657 - g(2) * t659 + 0.2e1 * (-t785 - t804) * t871) * MDP(10) + (-g(1) * t603 + g(2) * t606 - t447 * t554 - t448 * t553 - t451 * t747 + t452 * t771 - t463 * t521 - t464 * t520 - t479 * t494 + t480 * t493) * MDP(20) + (-t472 * t553 - t486 * t875 + t492 * t877 - t520 * t524) * MDP(25) + (t471 * t877 - t472 * t533 - t485 * t524 - t486 * t526) * MDP(23) + ((-qJD(6) * t750 - t450 * t712 + t466 * t716) * t875 + t749 * t492 + t444 * t553 - t751 * t520 + t449 * t524 + t474 * t472 - t445 * t877 + t461 * t486 - g(1) * t882 - g(2) * t528) * MDP(27) + (t774 * t638 + t773 * t595 + t777 * t649 + t515 * t645 + t545 * t612 + t560 * t531 + t507 * t625 + t543 * t570 + g(1) * t770 - g(2) * t587 + (-t516 * t649 - t638 * t825) * qJD(4)) * MDP(18) + (-(qJD(6) * t749 + t450 * t716 + t466 * t712) * t875 - t750 * t492 - t443 * t553 - t454 * t520 + t449 * t526 + t474 * t471 + t445 * t533 + t461 * t485 + g(1) * t883 - g(2) * t527) * MDP(28) + t757 * MDP(3) + (-g(2) * t719 + t866) * MDP(2) + (t471 * t533 + t485 * t526) * MDP(22) + (t514 * t573 + t549 * t546 + t513 * t572 - t548 * t545 - g(1) * (-t653 * t719 - t843) - g(2) * (-t653 * t715 + t677) + (pkin(2) * t654 * t814 + t698 * t881) * t707) * MDP(12) + (qJD(2) * t718 * t766 + t714 * t763) * t862 + (t718 * t763 - t766 * t814) * t861 + t710 * t844 + (qJDD(1) * t703 + 0.2e1 * t714 * t785) * t701 * MDP(4) + (t714 * t803 - t806 * t817) * MDP(5) * t872 + qJDD(1) * MDP(1) + (-t531 * t649 - t570 * t638 - t595 * t625 - t612 * t645) * MDP(16) + (t595 * t649 + t638 * t645) * MDP(17); (t516 * t647 + t696 * t530 + t562 * t746 + t638 * t824 + t713 * t721 + t717 * t733) * MDP(19) + (-t515 * t647 + t696 * t531 - t562 * t612 - t567 * t638 + (t563 * t638 + t733) * t713 - t721 * t717) * MDP(18) + ((t549 - t562) * t647 + (t548 - t563) * t644 + (t596 * t706 - t597 * t709) * pkin(2)) * MDP(11) + (t720 * t801 - t867 + g(2) * t657 + t681 + (-t797 - t864) * t707 + (-qJD(2) * t818 + t725) * qJD(1)) * MDP(9) + (-t730 * t875 + t756 + t795) * MDP(24) + (t448 * t611 - t447 * t610 + t484 * t744 - g(1) * (pkin(2) * t659 + t605 * t711 + t606 * t697) - g(2) * (-pkin(2) * t831 + t603 * t697 - t604 * t711 + t675) - g(3) * (pkin(2) * t837 - t649 * t697 - t650 * t711) + t876 * t523 + t819 * t464 + t821 * t463) * MDP(21) + (t530 * t713 - t746 * t768) * MDP(13) + ((t530 - t851) * t717 + (-t531 + t849) * t713) * MDP(14) + t844 + (t743 + t850) * MDP(16) + (t492 * t662 + t823 * t875) * MDP(26) + ((t601 * t716 - t611 * t712) * t492 + t444 * t662 + t610 * t472 + t445 * t712 * t664 - g(1) * (-t605 * t712 + t606 * t841) - g(2) * (t603 * t841 + t604 * t712) - g(3) * (-t649 * t841 + t650 * t712) + (t712 * t754 - t716 * t753) * t875 + t820 * t524 - t823 * t751 + t731 * t461) * MDP(27) + (-(t601 * t712 + t611 * t716) * t492 - t443 * t662 + t610 * t471 + t445 * t845 - g(1) * (-t605 * t716 - t606 * t842) - g(2) * (-t603 * t842 + t604 * t716) - g(3) * (t649 * t842 + t650 * t716) + (t712 * t753 + t716 * t754) * t875 + t820 * t526 - t823 * t454 - t730 * t461) * MDP(28) + (g(1) * t605 - g(2) * t604 - g(3) * t650 - t447 * t664 - t448 * t662 + t463 * t878 - t464 * t823 + t493 * t611 + t494 * t610 - t747 * t821 + t771 * t819) * MDP(20) + (t776 * t526 + t775 * t524 + (-t860 - t472 * t716 + (t524 * t712 - t526 * t716) * qJD(6)) * t664) * MDP(23) - t714 * MDP(4) * t794 - t638 * t647 * MDP(17) + (pkin(1) * t794 + g(1) * t660 - g(2) * t658 + t682 * t685 + (pkin(8) * t764 + g(3)) * t839 - t790) * MDP(10) + (-t731 * t875 - t796 - t884) * MDP(25) + (t638 * t768 + t832 + t848) * MDP(15) + (-g(2) * t675 + t548 * t562 - t549 * t563 + (t514 * t706 + t513 * t709 - t867 + g(2) * t831 + (-t654 * t815 - t864) * t707) * pkin(2)) * MDP(12) + (t471 * t845 - t526 * t730) * MDP(22) + (-t808 * t815 + t803) * t861 + (t718 * t764 + t804) * t862 + t817 * MDP(5) * t840; (-t644 ^ 2 - t647 ^ 2) * MDP(11) + (t548 * t647 - t549 * t644 + t782 + t802) * MDP(12) + (t743 - t850) * MDP(18) + (-t638 ^ 2 * t717 - t832 + t848) * MDP(19) + (t493 * t664 + t494 * t662 + t747 * t823 - t771 * t878) * MDP(20) + (-t447 * t662 + t448 * t664 - t463 * t823 - t464 * t878 - t523 * t647 + t782) * MDP(21) + (-t796 + t884) * MDP(27) + (t756 - t795) * MDP(28) + ((-pkin(1) * qJDD(1) - pkin(2) * t803 - t866) * MDP(12) - MDP(21) * t866) * t707 + (-MDP(27) * t731 + MDP(28) * t730) * t875; -t746 * t612 * MDP(13) + (-t612 ^ 2 + t746 ^ 2) * MDP(14) + (t530 + t851) * MDP(15) + (-t531 - t849) * MDP(16) + t595 * MDP(17) + (t516 * t638 + t543 * t746 + t722 + t874) * MDP(18) + (g(1) * t587 + g(2) * t770 + g(3) * t626 + t515 * t638 + t543 * t612 + t737) * MDP(19) + ((t493 * t705 - t494 * t708) * pkin(4) + (t463 - t468) * t771 + (t464 - t467) * t747) * MDP(20) + (t463 * t467 - t464 * t468 + (t447 * t708 + t448 * t705 + t523 * t746 + t874) * pkin(4)) * MDP(21) + (t526 * t767 + t860) * MDP(22) + ((t471 - t887) * t716 + (-t472 - t886) * t712) * MDP(23) + (-t856 - t885) * MDP(24) + (t742 + t857) * MDP(25) - t875 * t747 * MDP(26) + (-t467 * t524 + t695 * t472 + t723 * t712 - t716 * t873 + t747 * t751) * MDP(27) + (t454 * t747 - t467 * t526 + t695 * t471 + t712 * t873 + t723 * t716) * MDP(28); (-t747 ^ 2 - t771 ^ 2) * MDP(20) + (t463 * t747 - t464 * t771 + t484 + t732) * MDP(21) + (t742 - t857) * MDP(27) + (-t856 + t885) * MDP(28); t526 * t524 * MDP(22) + (-t524 ^ 2 + t526 ^ 2) * MDP(23) + (t791 + t887) * MDP(24) + (-t778 + t886) * MDP(25) + t492 * MDP(26) + (t454 * t875 - t461 * t526 - g(1) * t527 - g(2) * t883 - g(3) * (-t619 * t712 + t642) + t779) * MDP(27) + (-t751 * t875 + t461 * t524 + g(1) * t528 - g(2) * t882 - g(3) * (-t619 * t716 - t846) - t752) * MDP(28) + (-MDP(24) * t855 - MDP(25) * t526 - MDP(27) * t454 + MDP(28) * t751) * qJD(6);];
tau  = t1;
