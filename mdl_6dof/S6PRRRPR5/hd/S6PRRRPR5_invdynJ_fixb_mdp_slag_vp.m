% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:46
% EndTime: 2019-03-08 23:29:09
% DurationCPUTime: 17.12s
% Computational Cost: add. (9595->666), mult. (24442->948), div. (0->0), fcn. (20926->18), ass. (0->306)
t722 = sin(pkin(7));
t731 = sin(qJ(3));
t735 = cos(qJ(3));
t764 = t722 * (pkin(3) * t731 - pkin(10) * t735);
t723 = sin(pkin(6));
t732 = sin(qJ(2));
t869 = t723 * t732;
t816 = qJD(1) * t869;
t924 = qJD(3) * t764 - t722 * t816;
t726 = cos(pkin(7));
t736 = cos(qJ(2));
t858 = t735 * t736;
t863 = t731 * t732;
t761 = -t726 * t863 + t858;
t637 = t761 * t723;
t872 = t722 * t731;
t707 = pkin(9) * t872;
t864 = t726 * t735;
t906 = pkin(2) * t864 - t707;
t923 = qJD(1) * t637 - t906 * qJD(3);
t730 = sin(qJ(4));
t734 = cos(qJ(4));
t865 = t726 * t731;
t871 = t722 * t735;
t846 = pkin(2) * t865 + pkin(9) * t871;
t646 = pkin(10) * t726 + t846;
t778 = -pkin(3) * t735 - pkin(10) * t731;
t647 = (-pkin(2) + t778) * t722;
t907 = t734 * t646 + t730 * t647;
t922 = qJD(4) * t907 - t923 * t730 - t734 * t924;
t836 = qJD(4) * t734;
t837 = qJD(4) * t730;
t921 = -t646 * t837 + t647 * t836 + t730 * t924 - t923 * t734;
t840 = qJD(2) * t735;
t814 = t722 * t840;
t920 = qJD(4) - t814;
t668 = -t734 * t726 + t730 * t872;
t838 = qJD(3) * t735;
t811 = t722 * t838;
t783 = t734 * t811;
t607 = -qJD(4) * t668 + t783;
t821 = t734 * t872;
t669 = t726 * t730 + t821;
t839 = qJD(3) * t731;
t812 = t722 * t839;
t919 = pkin(4) * t812 - qJ(5) * t607 - qJD(5) * t669 - t922;
t809 = t730 * t838;
t608 = qJD(4) * t669 + t722 * t809;
t918 = -qJ(5) * t608 - qJD(5) * t668 + t921;
t842 = qJD(2) * t722;
t671 = pkin(9) * t842 + t816;
t727 = cos(pkin(6));
t844 = qJD(1) * t727;
t817 = t722 * t844;
t843 = qJD(1) * t736;
t683 = qJD(2) * pkin(2) + t723 * t843;
t878 = t683 * t726;
t579 = -t731 * t671 + t735 * (t817 + t878);
t653 = qJD(2) * t764;
t640 = t734 * t653;
t728 = -qJ(5) - pkin(10);
t803 = qJD(4) * t728;
t859 = t734 * t735;
t917 = -t640 - (pkin(4) * t731 - qJ(5) * t859) * t842 + t734 * t803 + (t579 - qJD(5)) * t730;
t784 = t730 * t814;
t853 = t734 * t579 + t730 * t653;
t916 = -qJ(5) * t784 - qJD(5) * t734 - t730 * t803 + t853;
t733 = cos(qJ(6));
t729 = sin(qJ(6));
t841 = qJD(2) * t726;
t706 = qJD(3) + t841;
t815 = t731 * t842;
t785 = t730 * t815;
t633 = -t734 * t706 + t785;
t635 = t706 * t730 + t734 * t815;
t721 = sin(pkin(13));
t724 = cos(pkin(13));
t768 = -t633 * t721 + t724 * t635;
t884 = t768 * t729;
t548 = -t733 * t920 + t884;
t795 = -t724 * t633 - t635 * t721;
t902 = qJD(6) - t795;
t915 = t548 * t902;
t550 = t729 * t920 + t733 * t768;
t914 = t550 * t902;
t672 = t721 * t730 - t724 * t734;
t913 = t920 * t672;
t912 = pkin(9) * qJDD(2) * t722 + (qJD(2) * t843 + qJDD(1) * t732) * t723 + qJD(3) * t817;
t861 = t732 * t735;
t862 = t731 * t736;
t763 = t726 * t861 + t862;
t636 = t763 * t723;
t810 = t726 * t839;
t848 = pkin(2) * t810 + pkin(9) * t811 - qJD(1) * t636;
t673 = t721 * t734 + t724 * t730;
t847 = t920 * t673;
t829 = qJDD(2) * t726;
t705 = qJDD(3) + t829;
t561 = qJD(2) * t783 - qJD(4) * t785 + qJDD(2) * t821 + t730 * t705 + t706 * t836;
t828 = qJDD(2) * t731;
t562 = -t734 * t705 + t706 * t837 + t722 * (qJD(2) * (t731 * t836 + t809) + t730 * t828);
t521 = -t561 * t721 - t724 * t562;
t520 = qJDD(6) - t521;
t792 = t902 * t733;
t911 = -t520 * t729 - t902 * t792;
t855 = -t918 * t721 + t724 * t919;
t854 = t721 * t919 + t918 * t724;
t852 = t721 * t916 + t724 * t917;
t850 = t721 * t917 - t724 * t916;
t762 = t726 * t862 + t861;
t604 = t723 * t762 + t727 * t872;
t868 = t723 * t736;
t661 = -t722 * t868 + t726 * t727;
t569 = -t604 * t730 + t661 * t734;
t570 = t604 * t734 + t661 * t730;
t525 = t569 * t721 + t570 * t724;
t903 = t726 * t858 - t863;
t603 = -t723 * t903 - t727 * t871;
t600 = t603 * t733;
t909 = -t525 * t729 + t600;
t908 = pkin(4) * t608 + t848;
t892 = sin(pkin(12));
t801 = t892 * t732;
t725 = cos(pkin(12));
t866 = t725 * t736;
t662 = t727 * t866 - t801;
t800 = t892 * t736;
t867 = t725 * t732;
t663 = t727 * t867 + t800;
t870 = t723 * t725;
t823 = t722 * t870;
t566 = t662 * t865 + t663 * t735 - t731 * t823;
t664 = -t727 * t800 - t867;
t665 = -t727 * t801 + t866;
t802 = t723 * t892;
t779 = t722 * t802;
t568 = t665 * t735 + (t664 * t726 + t779) * t731;
t605 = -t662 * t722 - t726 * t870;
t606 = -t664 * t722 + t726 * t802;
t905 = -g(1) * (-t568 * t730 + t606 * t734) - g(2) * (-t566 * t730 + t605 * t734) - g(3) * t569;
t777 = g(1) * t665 + g(2) * t663;
t747 = -g(3) * t869 - t777;
t580 = t735 * t671 + t683 * t865 + t731 * t817;
t904 = -t580 + (-t784 + t837) * pkin(4);
t901 = pkin(4) * t562 + qJDD(5);
t827 = qJDD(2) * t735;
t703 = t722 * t827;
t831 = qJD(2) * qJD(3);
t807 = t731 * t831;
t649 = t722 * t807 + qJDD(4) - t703;
t704 = qJDD(1) * t868;
t813 = qJD(2) * t869;
t782 = qJD(1) * t813;
t643 = qJDD(2) * pkin(2) + t704 - t782;
t830 = qJDD(1) * t727;
t805 = t722 * t830;
t744 = t643 * t865 - t671 * t839 + t731 * t805 + t735 * t912 + t838 * t878;
t518 = pkin(10) * t705 + t744;
t564 = pkin(10) * t706 + t580;
t702 = t726 * t844;
t596 = t702 + (qJD(2) * t778 - t683) * t722;
t529 = t564 * t734 + t596 * t730;
t701 = t726 * t830;
t758 = t807 - t827;
t806 = t735 * t831;
t759 = t806 + t828;
t547 = t701 + (pkin(3) * t758 - pkin(10) * t759 - t643) * t722;
t543 = t734 * t547;
t741 = -qJD(4) * t529 - t730 * t518 + t543;
t471 = pkin(4) * t649 - qJ(5) * t561 - qJD(5) * t635 + t741;
t757 = -t734 * t518 - t730 * t547 + t564 * t837 - t596 * t836;
t473 = -qJ(5) * t562 - qJD(5) * t633 - t757;
t464 = t471 * t724 - t473 * t721;
t462 = -pkin(5) * t649 - t464;
t712 = pkin(4) * t721 + pkin(11);
t718 = qJ(4) + pkin(13);
t715 = sin(t718);
t716 = cos(t718);
t900 = t902 * (pkin(4) * t635 + pkin(5) * t768 - pkin(11) * t795 + qJD(6) * t712) + g(1) * (-t568 * t715 + t606 * t716) + g(2) * (-t566 * t715 + t605 * t716) + g(3) * (-t604 * t715 + t661 * t716) + t462;
t897 = pkin(3) * t705;
t891 = MDP(7) * t722;
t522 = t561 * t724 - t562 * t721;
t834 = qJD(6) * t733;
t819 = t733 * t522 + t729 * t649 + t834 * t920;
t835 = qJD(6) * t729;
t486 = -t768 * t835 + t819;
t890 = t486 * t729;
t516 = -qJ(5) * t633 + t529;
t889 = t516 * t721;
t886 = t548 * t768;
t885 = t550 * t768;
t883 = t603 * t729;
t882 = t633 * t920;
t881 = t635 * t920;
t880 = t673 * t729;
t879 = t673 * t733;
t877 = t705 * MDP(9);
t876 = t715 * t722;
t875 = t716 * t729;
t874 = t716 * t733;
t510 = t724 * t516;
t737 = qJD(2) ^ 2;
t860 = t732 * t737;
t514 = t733 * t520;
t857 = qJDD(1) - g(3);
t465 = t721 * t471 + t724 * t473;
t856 = -pkin(5) * t812 - t855;
t528 = -t564 * t730 + t734 * t596;
t515 = -qJ(5) * t635 + t528;
t504 = pkin(4) * t920 + t515;
t481 = t721 * t504 + t510;
t794 = -t646 * t730 + t734 * t647;
t546 = -pkin(4) * t871 - qJ(5) * t669 + t794;
t554 = -qJ(5) * t668 + t907;
t501 = t721 * t546 + t724 * t554;
t851 = pkin(5) * t815 - t852;
t719 = t731 ^ 2;
t845 = -t735 ^ 2 + t719;
t833 = qJD(3) - t706;
t822 = t729 * t871;
t714 = pkin(4) * t734 + pkin(3);
t808 = t728 * t730;
t463 = pkin(11) * t649 + t465;
t765 = t643 * t864 - t671 * t838 - t683 * t810 - t731 * t912 + t735 * t805;
t519 = -t765 - t897;
t490 = t519 + t901;
t469 = -pkin(5) * t521 - pkin(11) * t522 + t490;
t799 = -t729 * t463 + t733 * t469;
t798 = t522 * t729 - t733 * t649;
t797 = t729 * t913 - t733 * t815;
t796 = t729 * t815 + t733 * t913;
t791 = t706 + t841;
t789 = t705 + t829;
t786 = t722 * t813;
t597 = t724 * t668 + t669 * t721;
t598 = -t668 * t721 + t669 * t724;
t645 = t707 + (-pkin(2) * t735 - pkin(3)) * t726;
t745 = pkin(4) * t668 + t645;
t523 = pkin(5) * t597 - pkin(11) * t598 + t745;
t775 = -pkin(11) * t812 - qJD(6) * t523 - t854;
t599 = pkin(5) * t672 - pkin(11) * t673 - t714;
t774 = pkin(11) * t815 - qJD(6) * t599 - t850;
t497 = -pkin(11) * t871 + t501;
t552 = t607 * t721 + t724 * t608;
t553 = t607 * t724 - t608 * t721;
t773 = -pkin(5) * t552 + pkin(11) * t553 + qJD(6) * t497 - t908;
t696 = t728 * t734;
t616 = -t724 * t696 + t721 * t808;
t772 = -pkin(5) * t847 - pkin(11) * t913 + qJD(6) * t616 - t904;
t771 = t733 * t463 + t729 * t469;
t479 = pkin(11) * t920 + t481;
t563 = -pkin(3) * t706 - t579;
t540 = pkin(4) * t633 + qJD(5) + t563;
t491 = -pkin(5) * t795 - pkin(11) * t768 + t540;
t467 = t479 * t733 + t491 * t729;
t770 = t479 * t729 - t491 * t733;
t480 = t504 * t724 - t889;
t769 = t525 * t733 + t883;
t500 = t546 * t724 - t554 * t721;
t766 = t514 + (t729 * t795 - t835) * t902;
t571 = t598 * t729 + t733 * t871;
t565 = -t662 * t864 + t663 * t731 + t735 * t823;
t567 = -t664 * t864 + t665 * t731 - t735 * t779;
t755 = g(1) * t567 + g(2) * t565 + g(3) * t603;
t754 = -g(1) * t568 - g(2) * t566 - g(3) * t604;
t584 = t662 * t731 + t663 * t864;
t586 = t664 * t731 + t665 * t864;
t753 = g(1) * t586 + g(2) * t584 + g(3) * t636;
t585 = t662 * t735 - t663 * t865;
t587 = t664 * t735 - t665 * t865;
t752 = g(1) * t587 + g(2) * t585 + g(3) * t637;
t751 = t673 * t834 - t797;
t750 = -t673 * t835 - t796;
t743 = -pkin(10) * t649 + t563 * t920;
t478 = -pkin(5) * t920 - t480;
t485 = t515 * t724 - t889;
t742 = -t712 * t520 + (t478 + t485) * t902;
t740 = -pkin(10) * qJD(4) * t920 - t519 + t755;
t739 = t755 + t765;
t717 = t722 ^ 2;
t713 = -pkin(4) * t724 - pkin(5);
t623 = -t683 * t722 + t702;
t615 = -t696 * t721 - t724 * t808;
t601 = -t643 * t722 + t701;
t590 = t637 * t716 + t869 * t876;
t572 = t598 * t733 - t822;
t559 = t727 * t811 + (t761 * qJD(2) + qJD(3) * t903) * t723;
t558 = t727 * t812 + (qJD(2) * t763 + qJD(3) * t762) * t723;
t557 = t604 * t716 + t661 * t715;
t545 = t587 * t716 + t665 * t876;
t544 = t585 * t716 + t663 * t876;
t533 = t568 * t716 + t606 * t715;
t531 = t566 * t716 + t605 * t715;
t524 = -t724 * t569 + t570 * t721;
t509 = qJD(4) * t569 + t559 * t734 + t730 * t786;
t508 = -qJD(4) * t570 - t559 * t730 + t734 * t786;
t507 = -qJD(6) * t822 + t553 * t729 + t598 * t834 - t733 * t812;
t506 = -qJD(6) * t571 + t553 * t733 + t729 * t812;
t496 = pkin(5) * t871 - t500;
t487 = qJD(6) * t550 + t798;
t484 = t515 * t721 + t510;
t483 = t508 * t721 + t509 * t724;
t482 = -t724 * t508 + t509 * t721;
t461 = -t467 * qJD(6) + t799;
t460 = -qJD(6) * t770 + t771;
t1 = [t857 * MDP(1) + (-t558 * t706 - t603 * t705) * MDP(10) + (-t559 * t706 - t604 * t705) * MDP(11) + (t508 * t920 + t558 * t633 + t562 * t603 + t569 * t649) * MDP(17) + (-t509 * t920 + t558 * t635 + t561 * t603 - t570 * t649) * MDP(18) + (t482 * t768 + t483 * t795 + t521 * t525 + t522 * t524) * MDP(19) + (-t464 * t524 + t465 * t525 - t480 * t482 + t481 * t483 + t490 * t603 + t540 * t558 - g(3)) * MDP(20) + ((-qJD(6) * t769 - t483 * t729 + t558 * t733) * t902 + t909 * t520 + t482 * t548 + t524 * t487) * MDP(26) + (-(qJD(6) * t909 + t483 * t733 + t558 * t729) * t902 - t769 * t520 + t482 * t550 + t524 * t486) * MDP(27) + (MDP(10) * t758 + MDP(11) * t759) * t722 * t661 + ((qJDD(2) * t736 - t860) * MDP(3) + (-qJDD(2) * t732 - t736 * t737) * MDP(4) + (-MDP(10) * t735 + MDP(11) * t731) * t717 * t860) * t723; (t519 * t668 + t645 * t562 + t563 * t608 + t848 * t633 + t794 * t649 - t752 * t734 - t920 * t922) * MDP(17) + (-t464 * t598 - t465 * t597 - t480 * t553 - t481 * t552 - t500 * t522 + t501 * t521 - t768 * t855 + t795 * t854 - t753) * MDP(19) + ((-t497 * t729 + t523 * t733) * t520 + t461 * t597 - t770 * t552 + t496 * t487 + t462 * t571 + t478 * t507 - g(1) * (t545 * t733 + t586 * t729) - g(2) * (t544 * t733 + t584 * t729) - g(3) * (t590 * t733 + t636 * t729) + (t729 * t775 - t733 * t773) * t902 + t856 * t548) * MDP(26) + (t731 * t789 + t791 * t838) * t891 + t726 * t877 + (-(t497 * t733 + t523 * t729) * t520 - t460 * t597 - t467 * t552 + t496 * t486 + t462 * t572 + t478 * t506 - g(1) * (-t545 * t729 + t586 * t733) - g(2) * (-t544 * t729 + t584 * t733) - g(3) * (-t590 * t729 + t636 * t733) + (t729 * t773 + t733 * t775) * t902 + t856 * t550) * MDP(27) + (t519 * t669 + t645 * t561 + t563 * t607 + t848 * t635 - t907 * t649 + t752 * t730 - t920 * t921) * MDP(18) + (-t608 * t920 - t649 * t668) * MDP(15) + (t607 * t920 + t649 * t669) * MDP(14) + (-t846 * t705 + t923 * t706 - t744 * t726 + t753) * MDP(11) + (t906 * t705 - t848 * t706 + t765 * t726 - t752) * MDP(10) + (t465 * t501 + t464 * t500 + t490 * t745 - g(1) * (pkin(2) * t664 - t586 * t728 + t587 * t714) - g(2) * (pkin(2) * t662 - t584 * t728 + t585 * t714) - g(3) * (pkin(2) * t868 - t636 * t728 + t637 * t714) + t908 * t540 + t854 * t481 + t855 * t480) * MDP(20) + (-g(1) * t664 - g(2) * t662 - g(3) * t868 + t704) * MDP(3) + (-t857 * t869 + t777) * MDP(4) + qJDD(2) * MDP(2) + (-t561 * t668 - t562 * t669 - t607 * t633 - t608 * t635) * MDP(13) + (t561 * t669 + t607 * t635) * MDP(12) + (-t487 * t597 - t507 * t902 - t520 * t571 - t548 * t552) * MDP(24) + (t486 * t597 + t506 * t902 + t520 * t572 + t550 * t552) * MDP(23) + (t520 * t597 + t552 * t902) * MDP(25) + (-t486 * t571 - t487 * t572 - t506 * t548 - t507 * t550) * MDP(22) + (t486 * t572 + t506 * t550) * MDP(21) + ((-(-t564 * t836 + t543) * t735 + t528 * t839 + (-(-qJD(4) * t596 - t518) * t735 + t747) * t730) * MDP(17) + (-t529 * t839 + t734 * t747 - t735 * t757) * MDP(18) + (t562 * t735 - t633 * t839) * MDP(15) + (-t561 * t735 + t635 * t839) * MDP(14) + (t601 * t731 + t623 * t838) * MDP(11) + (-t601 * t735 + t623 * t839) * MDP(10) + t747 * (pkin(4) * t730 + pkin(9)) * MDP(20) + (t735 * t789 - t791 * t839) * MDP(8) + (-t649 * t735 + t839 * t920) * MDP(16)) * t722 + ((-pkin(2) * t758 + t735 * t782) * MDP(10) + (-pkin(2) * t759 - t731 * t782) * MDP(11) + 0.2e1 * (t731 * t827 - t831 * t845) * MDP(6) + (qJDD(2) * t719 + 0.2e1 * t731 * t806) * MDP(5)) * t717; (t833 * t840 + t828) * t891 + (-t815 * t833 + t703) * MDP(8) + t877 + (t580 * t706 - t623 * t815 + t739) * MDP(10) + (t579 * t706 - t623 * t814 - t744 - t754) * MDP(11) + (t561 * t730 + t734 * t881) * MDP(12) + ((t561 - t882) * t734 + (-t562 - t881) * t730) * MDP(13) + (t920 * t836 + t649 * t730 + (-t635 * t731 - t859 * t920) * t842) * MDP(14) + (-t920 * t837 + t649 * t734 + (t730 * t735 * t920 + t633 * t731) * t842) * MDP(15) - t920 * MDP(16) * t815 + (-t528 * t815 - pkin(3) * t562 - t580 * t633 - t640 * t920 + (t579 * t920 + t743) * t730 + t740 * t734) * MDP(17) + (-pkin(3) * t561 + t529 * t815 - t580 * t635 - t730 * t740 + t734 * t743 + t853 * t920) * MDP(18) + (-t464 * t673 - t465 * t672 + t480 * t913 - t847 * t481 + t521 * t616 + t522 * t615 - t852 * t768 + t850 * t795 + t754) * MDP(19) + (t465 * t616 - t464 * t615 - t490 * t714 - g(1) * (-t567 * t714 - t568 * t728) - g(2) * (-t565 * t714 - t566 * t728) - g(3) * (-t603 * t714 - t604 * t728) + t904 * t540 + t850 * t481 + t852 * t480) * MDP(20) + (t486 * t879 + t550 * t750) * MDP(21) + (t797 * t550 + t796 * t548 + (-t890 - t487 * t733 + (t548 * t729 - t550 * t733) * qJD(6)) * t673) * MDP(22) + (t486 * t672 + t514 * t673 + t550 * t847 + t750 * t902) * MDP(23) + (-t487 * t672 - t520 * t880 - t548 * t847 - t751 * t902) * MDP(24) + (t520 * t672 + t847 * t902) * MDP(25) + ((t599 * t733 - t616 * t729) * t520 + t461 * t672 + t615 * t487 + t462 * t880 - g(1) * (-t567 * t874 + t568 * t729) - g(2) * (-t565 * t874 + t566 * t729) - g(3) * (-t603 * t874 + t604 * t729) + (t729 * t774 - t733 * t772) * t902 + t851 * t548 - t847 * t770 + t751 * t478) * MDP(26) + (-(t599 * t729 + t616 * t733) * t520 - t460 * t672 + t615 * t486 + t462 * t879 - g(1) * (t567 * t875 + t568 * t733) - g(2) * (t565 * t875 + t566 * t733) - g(3) * (t603 * t875 + t604 * t733) + (t729 * t772 + t733 * t774) * t902 + t851 * t550 - t847 * t467 + t750 * t478) * MDP(27) + (-MDP(5) * t731 * t735 + MDP(6) * t845) * t717 * t737; t635 * t633 * MDP(12) + (-t633 ^ 2 + t635 ^ 2) * MDP(13) + (t561 + t882) * MDP(14) + (-t562 + t881) * MDP(15) + t649 * MDP(16) + (t529 * t920 - t563 * t635 + t741 + t905) * MDP(17) + (t563 * t633 + t528 * t920 - g(1) * (-t568 * t734 - t606 * t730) - g(2) * (-t566 * t734 - t605 * t730) + g(3) * t570 + t757) * MDP(18) + ((t521 * t721 - t522 * t724) * pkin(4) + (t480 - t485) * t795 + (t481 - t484) * t768) * MDP(19) + (t480 * t484 - t481 * t485 + (t464 * t724 + t465 * t721 - t540 * t635 + t905) * pkin(4)) * MDP(20) + (t550 * t792 + t890) * MDP(21) + ((t486 - t915) * t733 + (-t487 - t914) * t729) * MDP(22) + (-t885 - t911) * MDP(23) + (t766 + t886) * MDP(24) - t902 * t768 * MDP(25) + (-t484 * t548 + t713 * t487 + t742 * t729 - t733 * t900 + t768 * t770) * MDP(26) + (t467 * t768 - t484 * t550 + t713 * t486 + t729 * t900 + t742 * t733) * MDP(27); (-t768 ^ 2 - t795 ^ 2) * MDP(19) + (t480 * t768 - t481 * t795 - t739 - t897 + t901) * MDP(20) + (t766 - t886) * MDP(26) + (-t885 + t911) * MDP(27); t550 * t548 * MDP(21) + (-t548 ^ 2 + t550 ^ 2) * MDP(22) + (t819 + t915) * MDP(23) + (-t798 + t914) * MDP(24) + t520 * MDP(25) + (t467 * t902 - t478 * t550 - g(1) * (-t533 * t729 + t567 * t733) - g(2) * (-t531 * t729 + t565 * t733) - g(3) * (-t557 * t729 + t600) + t799) * MDP(26) + (-t770 * t902 + t478 * t548 - g(1) * (-t533 * t733 - t567 * t729) - g(2) * (-t531 * t733 - t565 * t729) - g(3) * (-t557 * t733 - t883) - t771) * MDP(27) + (-MDP(23) * t884 - MDP(24) * t550 - MDP(26) * t467 + MDP(27) * t770) * qJD(6);];
tau  = t1;
