% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:40
% EndTime: 2019-03-09 22:00:03
% DurationCPUTime: 16.27s
% Computational Cost: add. (15789->651), mult. (36939->831), div. (0->0), fcn. (28256->18), ass. (0->304)
t792 = cos(qJ(3));
t787 = sin(qJ(3));
t788 = sin(qJ(2));
t892 = qJD(1) * t788;
t870 = t787 * t892;
t793 = cos(qJ(2));
t891 = qJD(1) * t793;
t686 = t792 * t891 - t870;
t687 = -t787 * t891 - t792 * t892;
t786 = sin(qJ(4));
t791 = cos(qJ(4));
t651 = t791 * t686 + t687 * t786;
t647 = qJD(6) - t651;
t783 = sin(pkin(11));
t784 = cos(pkin(11));
t785 = sin(qJ(6));
t790 = cos(qJ(6));
t700 = t783 * t785 - t790 * t784;
t973 = t647 * t700;
t878 = qJD(3) + qJD(4);
t770 = qJD(2) + t878;
t838 = t686 * t786 - t791 * t687;
t635 = -t784 * t770 + t783 * t838;
t637 = t770 * t783 + t784 * t838;
t926 = t637 * t785;
t971 = -t790 * t635 - t926;
t974 = t647 * t971;
t701 = t783 * t790 + t784 * t785;
t684 = t701 * qJD(6);
t972 = -t701 * t651 + t684;
t967 = t651 * t783;
t642 = pkin(5) * t967;
t877 = pkin(10) * t967;
t704 = t787 * t788 - t792 * t793;
t775 = t793 * pkin(2);
t931 = pkin(1) + t775;
t674 = pkin(3) * t704 - t931;
t727 = qJD(1) * t931;
t970 = qJDD(1) * t931;
t680 = t687 * pkin(9);
t945 = pkin(7) + pkin(8);
t729 = t945 * t793;
t713 = qJD(1) * t729;
t688 = t787 * t713;
t728 = t945 * t788;
t711 = qJD(1) * t728;
t930 = qJD(2) * pkin(2);
t696 = -t711 + t930;
t859 = t792 * t696 - t688;
t632 = t680 + t859;
t778 = qJD(2) + qJD(3);
t622 = pkin(3) * t778 + t632;
t692 = t792 * t713;
t837 = -t696 * t787 - t692;
t937 = pkin(9) * t686;
t633 = -t837 + t937;
t628 = t791 * t633;
t577 = t786 * t622 + t628;
t574 = qJ(5) * t770 + t577;
t666 = -pkin(3) * t686 - t727;
t589 = -pkin(4) * t651 - qJ(5) * t838 + t666;
t534 = -t574 * t783 + t784 * t589;
t515 = -pkin(5) * t651 - pkin(10) * t637 + t534;
t535 = t784 * t574 + t783 * t589;
t520 = -pkin(10) * t635 + t535;
t494 = t515 * t785 + t520 * t790;
t879 = qJDD(1) * t793;
t882 = qJD(1) * qJD(2);
t868 = t793 * t882;
t880 = qJDD(1) * t788;
t825 = -t868 - t880;
t881 = qJD(1) * qJD(3);
t962 = -t793 * t881 + t825;
t638 = -t778 * t870 + t787 * t879 - t792 * t962;
t776 = qJDD(2) + qJDD(3);
t662 = qJDD(2) * pkin(2) + t825 * t945;
t869 = t788 * t882;
t824 = -t869 + t879;
t665 = t945 * t824;
t808 = qJD(3) * t837 + t792 * t662 - t787 * t665;
t559 = pkin(3) * t776 - pkin(9) * t638 + t808;
t890 = qJD(3) * t787;
t682 = t713 * t890;
t855 = -qJD(3) * t696 - t665;
t563 = -t682 + (pkin(9) * t962 + t662) * t787 + ((-t788 * t881 + t824) * pkin(9) - t855) * t792;
t887 = qJD(4) * t791;
t888 = qJD(4) * t786;
t851 = -t791 * t559 + t786 * t563 + t622 * t888 + t633 * t887;
t767 = qJDD(4) + t776;
t953 = -pkin(4) * t767 + qJDD(5);
t510 = t851 + t953;
t705 = t787 * t793 + t788 * t792;
t661 = t778 * t705;
t800 = qJD(1) * t661;
t797 = -t704 * qJDD(1) - t800;
t565 = t791 * t638 + t686 * t887 + t687 * t888 + t786 * t797;
t555 = t565 * t783 - t784 * t767;
t497 = pkin(5) * t555 + t510;
t626 = t786 * t633;
t576 = t622 * t791 - t626;
t573 = -pkin(4) * t770 + qJD(5) - t576;
t549 = pkin(5) * t635 + t573;
t782 = qJ(2) + qJ(3);
t774 = qJ(4) + t782;
t761 = cos(t774);
t777 = pkin(11) + qJ(6);
t768 = sin(t777);
t760 = sin(t774);
t789 = sin(qJ(1));
t794 = cos(qJ(1));
t846 = g(1) * t794 + g(2) * t789;
t834 = t760 * t846;
t933 = g(3) * t768;
t801 = t494 * t838 + t497 * t701 - t549 * t973 + t761 * t933 - t768 * t834;
t493 = t515 * t790 - t520 * t785;
t769 = cos(t777);
t932 = g(3) * t769;
t913 = t760 * t794;
t914 = t760 * t789;
t963 = g(1) * t913 + g(2) * t914;
t804 = -t493 * t838 + t497 * t700 + t549 * t972 - t761 * t932 + t769 * t963;
t556 = t565 * t784 + t767 * t783;
t885 = qJD(6) * t790;
t875 = -t785 * t555 + t790 * t556 - t635 * t885;
t886 = qJD(6) * t785;
t508 = -t637 * t886 + t875;
t840 = t635 * t785 - t637 * t790;
t862 = t790 * t555 + t785 * t556;
t509 = -qJD(6) * t840 + t862;
t566 = qJD(4) * t838 + t638 * t786 - t791 * t797;
t564 = qJDD(6) + t566;
t969 = t767 * MDP(22) - t566 * MDP(21) - t651 ^ 2 * MDP(19) + (-t651 * t770 + t565) * MDP(20) + (-MDP(18) * t651 + MDP(19) * t838 + MDP(21) * t770) * t838 + (-t700 * t564 - t838 * t971) * MDP(32) + (-t508 * t700 - t971 * t973) * MDP(30) + (MDP(29) * t508 - MDP(30) * t509 + MDP(31) * t564) * t701 + (MDP(29) * t973 + MDP(30) * t972 + t838 * MDP(31)) * t840 + (-MDP(31) * t973 - MDP(32) * t972 - MDP(33) * t838) * t647;
t968 = t647 * t840;
t772 = cos(t782);
t896 = t761 * pkin(4) + t760 * qJ(5);
t873 = pkin(3) * t772 + t896;
t961 = -g(3) * t761 + t963;
t612 = pkin(4) * t838 - qJ(5) * t651;
t910 = t761 * t794;
t911 = t761 * t789;
t874 = g(1) * t910 + g(2) * t911 + g(3) * t760;
t947 = -(qJD(4) * t622 + t563) * t791 - t786 * t559 + t633 * t888;
t803 = -t666 * t651 + t874 + t947;
t957 = pkin(5) * t838;
t845 = g(1) * t789 - g(2) * t794;
t956 = t845 * t760;
t585 = t632 * t791 - t626;
t944 = pkin(3) * t687;
t597 = t612 - t944;
t537 = -t585 * t783 + t784 * t597;
t740 = pkin(3) * t887 + qJD(5);
t955 = -t740 * t783 - t537;
t538 = t784 * t585 + t783 * t597;
t864 = t740 * t784 - t538;
t858 = t711 * t787 - t692;
t640 = t858 - t937;
t899 = -t792 * t711 - t688;
t641 = t680 + t899;
t591 = t640 * t786 + t641 * t791;
t765 = pkin(2) * t892;
t592 = t597 + t765;
t539 = -t591 * t783 + t784 * t592;
t751 = t786 * t787 * pkin(2);
t763 = pkin(2) * t792 + pkin(3);
t889 = qJD(3) * t792;
t948 = t791 * pkin(2) * t889 - t751 * t878 + t763 * t887;
t654 = qJD(5) + t948;
t954 = -t654 * t783 - t539;
t540 = t784 * t591 + t783 * t592;
t863 = t654 * t784 - t540;
t583 = t632 * t786 + t628;
t848 = pkin(3) * t888 - t583;
t905 = t787 * t791;
t900 = t640 * t791 - t641 * t786 + t763 * t888 + (t787 * t887 + (t786 * t792 + t905) * qJD(3)) * pkin(2);
t898 = -t787 * t728 + t792 * t729;
t771 = sin(t782);
t847 = -pkin(3) * t771 - pkin(4) * t760;
t812 = -t534 * t838 + (-t510 + t961) * t784;
t912 = t761 * t783;
t806 = t535 * t838 + g(3) * t912 + (t510 - t834) * t783;
t817 = -t851 + t961;
t810 = -t666 * t838 + t817;
t941 = pkin(3) * t791;
t938 = pkin(5) * t784;
t773 = t784 * pkin(10);
t929 = qJ(5) * t784;
t506 = qJ(5) * t767 + qJD(5) * t770 - t947;
t750 = pkin(2) * t869;
t617 = pkin(3) * t800 + qJDD(1) * t674 + t750;
t513 = t566 * pkin(4) - t565 * qJ(5) - qJD(5) * t838 + t617;
t489 = t784 * t506 + t783 * t513;
t487 = t489 * t784;
t928 = t573 * t651;
t658 = t791 * t704 + t705 * t786;
t660 = t778 * t704;
t593 = -qJD(4) * t658 - t660 * t791 - t661 * t786;
t927 = t593 * t783;
t923 = t651 * t784;
t659 = -t704 * t786 + t705 * t791;
t920 = t659 * t783;
t919 = t659 * t784;
t895 = pkin(2) * t905 + t786 * t763;
t678 = qJ(5) + t895;
t918 = t678 * t784;
t754 = pkin(3) * t786 + qJ(5);
t915 = t754 * t784;
t909 = t768 * t789;
t908 = t768 * t794;
t907 = t769 * t789;
t906 = t769 * t794;
t594 = qJD(4) * t659 - t660 * t786 + t791 * t661;
t766 = t788 * t930;
t653 = pkin(3) * t661 + t766;
t525 = pkin(4) * t594 - qJ(5) * t593 - qJD(5) * t659 + t653;
t871 = qJD(2) * t945;
t712 = t788 * t871;
t714 = t793 * t871;
t823 = -t792 * t712 - t787 * t714 - t728 * t889 - t729 * t890;
t600 = -pkin(9) * t661 + t823;
t807 = -qJD(3) * t898 + t712 * t787 - t792 * t714;
t601 = pkin(9) * t660 + t807;
t857 = -t792 * t728 - t729 * t787;
t644 = -pkin(9) * t705 + t857;
t645 = -pkin(9) * t704 + t898;
t839 = t644 * t791 - t645 * t786;
t528 = qJD(4) * t839 + t600 * t791 + t601 * t786;
t499 = t783 * t525 + t784 * t528;
t542 = t784 * t576 + t783 * t612;
t607 = pkin(4) * t658 - qJ(5) * t659 + t674;
t609 = t644 * t786 + t645 * t791;
t548 = t783 * t607 + t784 * t609;
t901 = -t642 + t900;
t780 = t788 ^ 2;
t894 = -t793 ^ 2 + t780;
t883 = -qJD(5) + t573;
t755 = -pkin(4) - t938;
t866 = -pkin(2) * t788 + t847;
t488 = -t506 * t783 + t784 * t513;
t485 = pkin(5) * t566 - pkin(10) * t556 + t488;
t486 = -pkin(10) * t555 + t489;
t865 = t790 * t485 - t486 * t785;
t498 = t784 * t525 - t528 * t783;
t541 = -t576 * t783 + t784 * t612;
t547 = t784 * t607 - t609 * t783;
t856 = t763 * t791 - t751;
t853 = t487 - t874;
t850 = -pkin(10) * t923 + t957;
t679 = -pkin(4) - t856;
t849 = -t642 + t848;
t844 = t485 * t785 + t486 * t790;
t843 = -t488 * t783 + t487;
t532 = pkin(5) * t658 - pkin(10) * t919 + t547;
t536 = -pkin(10) * t920 + t548;
t842 = t532 * t790 - t536 * t785;
t841 = t532 * t785 + t536 * t790;
t836 = t931 + t873;
t835 = t534 * t923 + t535 * t967 + t853;
t833 = t845 * t761;
t832 = -0.2e1 * pkin(1) * t882 - pkin(7) * qJDD(2);
t695 = t773 + t915;
t831 = qJD(6) * t695 + t850 - t955;
t664 = t773 + t918;
t830 = qJD(6) * t664 + t850 - t954;
t694 = (-pkin(10) - t754) * t783;
t829 = -qJD(6) * t694 - t864 - t877;
t663 = (-pkin(10) - t678) * t783;
t828 = -qJD(6) * t663 - t863 - t877;
t722 = t773 + t929;
t827 = qJD(5) * t783 + qJD(6) * t722 - t651 * t773 + t541 + t957;
t721 = (-pkin(10) - qJ(5)) * t783;
t826 = -qJD(5) * t784 - qJD(6) * t721 + t542 - t877;
t819 = -t566 * t678 + t651 * t654 - t928;
t818 = -t566 * t754 + t651 * t740 - t928;
t795 = qJD(2) ^ 2;
t816 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t795 + t845;
t796 = qJD(1) ^ 2;
t815 = pkin(1) * t796 - pkin(7) * qJDD(1) + t846;
t813 = t510 * t659 + t573 * t593 - t846;
t529 = qJD(4) * t609 + t600 * t786 - t601 * t791;
t802 = g(3) * t771 - t787 * t662 + t727 * t686 + t772 * t846 + t792 * t855 + t682;
t799 = -g(3) * t772 - t727 * t687 + t771 * t846 + t808;
t798 = t687 * t686 * MDP(11) + (-t686 * t778 + t638) * MDP(13) + (-t687 * t778 + t797) * MDP(14) + (-t686 ^ 2 + t687 ^ 2) * MDP(12) + t776 * MDP(15) + t969;
t779 = -pkin(9) - t945;
t762 = -pkin(4) - t941;
t719 = qJ(5) * t910;
t718 = qJ(5) * t911;
t717 = t755 - t941;
t681 = t750 - t970;
t672 = t679 - t938;
t671 = t761 * t906 + t909;
t670 = -t761 * t908 + t907;
t669 = -t761 * t907 + t908;
t668 = t761 * t909 + t906;
t667 = t765 - t944;
t614 = t700 * t659;
t613 = t701 * t659;
t572 = pkin(5) * t920 - t839;
t558 = t577 + t642;
t517 = t593 * t701 + t885 * t919 - t886 * t920;
t516 = -t593 * t700 - t659 * t684;
t514 = pkin(5) * t927 + t529;
t492 = -pkin(10) * t927 + t499;
t491 = pkin(5) * t594 - t593 * t773 + t498;
t1 = [(t788 * t832 + t793 * t816) * MDP(9) + (-t788 * t816 + t793 * t832) * MDP(10) + (t565 * t659 + t593 * t838) * MDP(18) + t845 * MDP(2) + t846 * MDP(3) + (-t528 * t770 + t565 * t674 + t593 * t666 - t609 * t767 + t617 * t659 + t653 * t838 - t956) * MDP(24) + (-t498 * t637 - t499 * t635 - t547 * t556 - t548 * t555 + t956 + (-t488 * t784 - t489 * t783) * t659 + (-t534 * t784 - t535 * t783) * t593) * MDP(27) + (t638 * t705 + t660 * t687) * MDP(11) + (-t638 * t704 - t660 * t686 + t687 * t661 + t705 * t797) * MDP(12) + (-t660 * t778 + t705 * t776) * MDP(13) + (t508 * t658 + t516 * t647 - t564 * t614 - t594 * t840) * MDP(31) + (-(t491 * t785 + t492 * t790) * t647 - t841 * t564 - t844 * t658 - t494 * t594 - t514 * t840 + t572 * t508 - t497 * t614 + t549 * t516 - g(1) * t668 - g(2) * t670 + (-t493 * t658 - t647 * t842) * qJD(6)) * MDP(35) + (-t508 * t614 - t516 * t840) * MDP(29) + (t488 * t547 + t489 * t548 + t534 * t498 + t535 * t499 - t510 * t839 + t573 * t529 + (g(1) * t779 - g(2) * t836) * t794 + (g(1) * t836 + g(2) * t779) * t789) * MDP(28) + (-t565 * t658 - t566 * t659 + t593 * t651 - t594 * t838) * MDP(19) + (t488 * t658 - t498 * t651 + t529 * t635 + t534 * t594 + t547 * t566 - t555 * t839 + t783 * t813 + t784 * t833) * MDP(25) + (-t529 * t770 + t566 * t674 + t594 * t666 + t617 * t658 - t651 * t653 + t767 * t839 + t833) * MDP(23) + (-t489 * t658 + t499 * t651 + t529 * t637 - t535 * t594 - t548 * t566 - t556 * t839 + t784 * t813 - t845 * t912) * MDP(26) + (-t638 * t931 + t727 * t660 + t681 * t705 - t687 * t766 - t771 * t845 - t776 * t898 - t778 * t823) * MDP(17) + (-t594 * t770 - t658 * t767) * MDP(21) + (t593 * t770 + t659 * t767) * MDP(20) + (qJDD(1) * t780 + 0.2e1 * t788 * t868) * MDP(4) + (qJDD(2) * t788 + t793 * t795) * MDP(6) + (qJDD(2) * t793 - t788 * t795) * MDP(7) + (t564 * t658 + t594 * t647) * MDP(33) + qJDD(1) * MDP(1) + (-t686 * t766 + t772 * t845 + t776 * t857 + t778 * t807 + (t681 - t970) * t704 - 0.2e1 * t727 * t661) * MDP(16) + 0.2e1 * (t788 * t879 - t882 * t894) * MDP(5) + (-t661 * t778 - t704 * t776) * MDP(14) + ((t491 * t790 - t492 * t785) * t647 + t842 * t564 + t865 * t658 + t493 * t594 - t514 * t971 + t572 * t509 + t497 * t613 + t549 * t517 - g(1) * t669 - g(2) * t671 + (-t494 * t658 - t647 * t841) * qJD(6)) * MDP(34) + (-t508 * t613 + t509 * t614 + t516 * t971 + t517 * t840) * MDP(30) + (-t509 * t658 - t517 * t647 - t564 * t613 + t594 * t971) * MDP(32); (t539 * t651 + t555 * t679 + t635 * t900 + t783 * t819 + t812) * MDP(25) + (-t555 * t918 + t539 * t637 - t863 * t635 + (t556 * t678 + t637 * t654 - t488) * t783 + t835) * MDP(27) + (-t895 * t767 - t667 * t838 + (t591 - t948) * t770 + t803) * MDP(24) + (-g(3) * t793 + t788 * t815) * MDP(9) + (t510 * t679 - g(1) * (t794 * t866 + t719) - g(2) * (t789 * t866 + t718) - g(3) * (t775 + t873) + t843 * t678 + t900 * t573 + t863 * t535 + t954 * t534) * MDP(28) + MDP(7) * t879 + MDP(6) * t880 + (t651 * t667 + t767 * t856 - t770 * t900 + t810) * MDP(23) + qJDD(2) * MDP(8) + (-t540 * t651 + t556 * t679 + t637 * t900 + t784 * t819 + t806) * MDP(26) + (g(3) * t788 + t793 * t815) * MDP(10) + ((t663 * t790 - t664 * t785) * t564 + t672 * t509 + (t785 * t828 - t790 * t830) * t647 - t901 * t971 + t804) * MDP(34) + (t899 * t778 + (t687 * t892 - t776 * t787 - t778 * t889) * pkin(2) + t802) * MDP(17) + (-t858 * t778 + (t686 * t892 + t776 * t792 - t778 * t890) * pkin(2) + t799) * MDP(16) + (-(t663 * t785 + t664 * t790) * t564 + t672 * t508 + (t785 * t830 + t790 * t828) * t647 - t901 * t840 + t801) * MDP(35) + t798 + (-MDP(4) * t788 * t793 + MDP(5) * t894) * t796; (t537 * t651 + t555 * t762 + t635 * t848 + t783 * t818 + t812) * MDP(25) + (t583 * t770 + (-t651 * t687 + t767 * t791 - t770 * t888) * pkin(3) + t810) * MDP(23) + (-t555 * t915 + t537 * t637 - t864 * t635 + (t556 * t754 + t637 * t740 - t488) * t783 + t835) * MDP(27) + ((t694 * t790 - t695 * t785) * t564 + t717 * t509 + (t785 * t829 - t790 * t831) * t647 - t849 * t971 + t804) * MDP(34) + (t585 * t770 + (t687 * t838 - t767 * t786 - t770 * t887) * pkin(3) + t803) * MDP(24) + (t510 * t762 - g(1) * (t794 * t847 + t719) - g(2) * (t789 * t847 + t718) - g(3) * t873 + t843 * t754 + t848 * t573 + t864 * t535 + t955 * t534) * MDP(28) + (-t538 * t651 + t556 * t762 + t637 * t848 + t784 * t818 + t806) * MDP(26) + (-(t694 * t785 + t695 * t790) * t564 + t717 * t508 + (t785 * t831 + t790 * t829) * t647 - t849 * t840 + t801) * MDP(35) + (-t778 * t837 + t799) * MDP(16) + (t778 * t859 + t802) * MDP(17) + t798; (t577 * t770 + t810) * MDP(23) + (t576 * t770 + t803) * MDP(24) + (-qJ(5) * t566 * t783 - pkin(4) * t555 - t577 * t635 - (t783 * t883 - t541) * t651 + t812) * MDP(25) + (-t566 * t929 - pkin(4) * t556 - t577 * t637 - (t784 * t883 + t542) * t651 + t806) * MDP(26) + (t541 * t637 + t542 * t635 + (-qJ(5) * t555 - qJD(5) * t635 + t534 * t651) * t784 + (qJ(5) * t556 + qJD(5) * t637 + t535 * t651 - t488) * t783 + t853) * MDP(27) + (-t510 * pkin(4) - t535 * t542 - t534 * t541 - t573 * t577 - g(1) * (-pkin(4) * t913 + t719) - g(2) * (-pkin(4) * t914 + t718) - g(3) * t896 + (-t534 * t783 + t535 * t784) * qJD(5) + t843 * qJ(5)) * MDP(28) + ((t721 * t790 - t722 * t785) * t564 + t755 * t509 + t558 * t971 + (t785 * t826 - t790 * t827) * t647 + t804) * MDP(34) + (-(t721 * t785 + t722 * t790) * t564 + t755 * t508 + t558 * t840 + (t785 * t827 + t790 * t826) * t647 + t801) * MDP(35) + t969; (-t637 * t651 + t555) * MDP(25) + (t635 * t651 + t556) * MDP(26) + (-t635 ^ 2 - t637 ^ 2) * MDP(27) + (t534 * t637 + t535 * t635 - t817 + t953) * MDP(28) + (t509 - t968) * MDP(34) + (t508 + t974) * MDP(35); t840 * t971 * MDP(29) + (t840 ^ 2 - t971 ^ 2) * MDP(30) + (t875 - t974) * MDP(31) + (-t862 - t968) * MDP(32) + t564 * MDP(33) + (-g(1) * t670 + g(2) * t668 + t494 * t647 + t549 * t840 + t760 * t933 + t865) * MDP(34) + (g(1) * t671 - g(2) * t669 + t493 * t647 - t549 * t971 + t760 * t932 - t844) * MDP(35) + (-MDP(31) * t926 + MDP(32) * t840 - MDP(34) * t494 - MDP(35) * t493) * qJD(6);];
tau  = t1;
