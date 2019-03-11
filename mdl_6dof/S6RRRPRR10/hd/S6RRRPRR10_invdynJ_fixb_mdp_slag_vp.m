% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:22:21
% EndTime: 2019-03-09 19:22:48
% DurationCPUTime: 19.78s
% Computational Cost: add. (7766->730), mult. (16698->929), div. (0->0), fcn. (11983->12), ass. (0->289)
t773 = cos(qJ(2));
t751 = t773 * qJDD(1);
t768 = sin(qJ(2));
t864 = qJD(1) * qJD(2);
t959 = -t768 * t864 + t751;
t687 = qJDD(3) - t959;
t675 = -qJDD(5) + t687;
t662 = -qJDD(6) + t675;
t767 = sin(qJ(3));
t772 = cos(qJ(3));
t867 = t772 * qJD(2);
t881 = qJD(1) * t768;
t689 = t767 * t881 - t867;
t857 = t772 * t881;
t878 = qJD(2) * t767;
t691 = t857 + t878;
t766 = sin(qJ(5));
t771 = cos(qJ(5));
t617 = t689 * t766 + t691 * t771;
t765 = sin(qJ(6));
t770 = cos(qJ(6));
t809 = -t771 * t689 + t691 * t766;
t972 = -t770 * t617 + t765 * t809;
t924 = t617 * t765;
t987 = t770 * t809 + t924;
t999 = -t662 * MDP(33) - (-t972 ^ 2 + t987 ^ 2) * MDP(30) - t987 * MDP(29) * t972;
t847 = t773 * t864;
t863 = qJDD(1) * t768;
t874 = qJD(3) * t768;
t976 = qJD(1) * t874 - qJDD(2);
t605 = -qJD(3) * t867 + (-t847 - t863) * t772 + t976 * t767;
t880 = qJD(1) * t773;
t606 = t767 * (qJD(2) * (qJD(3) + t880) + t863) + t976 * t772;
t871 = qJD(5) * t771;
t872 = qJD(5) * t766;
t528 = -t771 * t605 + t766 * t606 + t689 * t871 - t691 * t872;
t838 = -t605 * t766 - t771 * t606;
t529 = qJD(5) * t617 + t838;
t869 = qJD(6) * t770;
t861 = t770 * t528 - t765 * t529 - t809 * t869;
t870 = qJD(6) * t765;
t501 = -t617 * t870 + t861;
t841 = t528 * t765 + t770 * t529;
t502 = -qJD(6) * t972 + t841;
t734 = -qJD(3) + t880;
t865 = -qJD(5) - t734;
t977 = qJD(6) - t865;
t980 = t977 * t972;
t981 = t987 * t977;
t998 = (t501 + t981) * MDP(31) - (t502 + t980) * MDP(32) + (t617 ^ 2 - t809 ^ 2) * MDP(23) + t809 * MDP(22) * t617 - (t809 * t865 - t528) * MDP(24) - t675 * MDP(26) + t999;
t873 = qJD(3) * t772;
t993 = t772 * t880 - t873;
t858 = t767 * t880;
t875 = qJD(3) * t767;
t992 = t858 - t875;
t803 = pkin(2) * t773 + pkin(8) * t768 + pkin(1);
t676 = t803 * qJD(1);
t748 = pkin(7) * t880;
t717 = qJD(2) * pkin(8) + t748;
t623 = -t772 * t676 - t767 * t717;
t866 = qJD(4) - t623;
t716 = -qJD(2) * pkin(2) + pkin(7) * t881;
t600 = t689 * pkin(3) - t691 * qJ(4) + t716;
t567 = -pkin(4) * t689 - t600;
t540 = pkin(5) * t809 + t567;
t774 = cos(qJ(1));
t905 = t774 * t767;
t769 = sin(qJ(1));
t909 = t769 * t772;
t667 = t773 * t905 - t909;
t906 = t772 * t774;
t668 = t767 * t769 + t773 * t906;
t762 = qJ(5) + qJ(6);
t754 = sin(t762);
t755 = cos(t762);
t586 = t667 * t754 + t668 * t755;
t806 = t754 * t767 + t755 * t772;
t911 = t767 * t773;
t665 = t769 * t911 + t906;
t907 = t772 * t773;
t666 = t769 * t907 - t905;
t811 = t665 * t754 + t666 * t755;
t942 = pkin(3) + pkin(4);
t986 = pkin(9) * t691 - t866;
t556 = t734 * t942 - t986;
t624 = -t767 * t676 + t772 * t717;
t577 = pkin(9) * t689 + t624;
t722 = t734 * qJ(4);
t563 = t577 - t722;
t522 = t556 * t766 + t563 * t771;
t826 = pkin(2) * t768 - pkin(8) * t773;
t704 = t826 * qJD(2);
t631 = qJD(1) * t704 - qJDD(1) * t803;
t660 = pkin(7) * t959 + qJDD(2) * pkin(8);
t829 = -t772 * t631 + t767 * t660 - t676 * t875 + t717 * t873;
t804 = qJDD(4) + t829;
t515 = pkin(9) * t605 - t687 * t942 + t804;
t674 = t687 * qJ(4);
t720 = t734 * qJD(4);
t792 = t767 * t631 + t772 * t660 - t676 * t873 - t717 * t875;
t531 = t674 - t720 + t792;
t517 = pkin(9) * t606 + t531;
t842 = t771 * t515 - t766 * t517;
t780 = -t522 * qJD(5) + t842;
t496 = -pkin(5) * t675 - pkin(10) * t528 + t780;
t791 = -t766 * t515 - t771 * t517 - t556 * t871 + t563 * t872;
t498 = -pkin(10) * t529 - t791;
t521 = t771 * t556 - t563 * t766;
t982 = pkin(10) * t617;
t511 = t521 - t982;
t509 = -pkin(5) * t865 + t511;
t967 = pkin(10) * t809;
t512 = t522 - t967;
t830 = t765 * t496 + t770 * t498 + t509 * t869 - t512 * t870;
t931 = g(3) * t768;
t974 = g(1) * t586 + g(2) * t811 + t540 * t987 + t806 * t931 - t830;
t941 = pkin(8) - pkin(9);
t719 = t941 * t772;
t859 = -pkin(7) * t767 - pkin(3);
t785 = -pkin(9) * t907 + (-pkin(4) + t859) * t768;
t701 = t826 * qJD(1);
t916 = t701 * t772;
t989 = qJD(1) * t785 - qJD(3) * t719 - t916;
t670 = t767 * t701;
t892 = qJ(4) * t881 + t670;
t910 = t768 * t772;
t988 = (-pkin(7) * t910 + pkin(9) * t911) * qJD(1) + t892 + t941 * t875;
t602 = t667 * t771 - t668 * t766;
t912 = t767 * t771;
t656 = t766 * t910 - t768 * t912;
t956 = t665 * t771 - t666 * t766;
t985 = -g(1) * t602 - g(2) * t956 + g(3) * t656 - t567 * t617 + t780;
t585 = t667 * t755 - t668 * t754;
t843 = t770 * t496 - t765 * t498;
t913 = t767 * t768;
t955 = t665 * t755 - t666 * t754;
t983 = g(2) * t955 - t540 * t972 - g(3) * (t754 * t910 - t755 * t913) - t843 + g(1) * t585;
t693 = t766 * t767 + t771 * t772;
t794 = t773 * t693;
t953 = qJD(3) - qJD(5);
t896 = -qJD(1) * t794 + t953 * t693;
t895 = -t766 * t993 + t767 * t871 + t771 * t992 - t772 * t872;
t978 = qJ(4) * t993 - t767 * qJD(4) - t748;
t930 = g(3) * t773;
t823 = g(1) * t774 + g(2) * t769;
t951 = t768 * t823;
t945 = t930 - t951;
t736 = pkin(7) * t911;
t759 = t773 * pkin(3);
t610 = pkin(4) * t773 + t736 + t759 + (-pkin(9) * t768 + t803) * t772;
t737 = pkin(7) * t907;
t886 = -t767 * t803 + t737;
t634 = -qJ(4) * t773 + t886;
t622 = pkin(9) * t913 + t634;
t897 = t766 * t610 + t771 * t622;
t862 = t942 * t767;
t894 = -qJD(3) * t862 + t858 * t942 - t978;
t844 = -qJ(4) * t766 - t771 * t942;
t964 = -qJD(5) * t844 + t766 * t577 + t771 * t986;
t705 = qJ(4) * t771 - t766 * t942;
t963 = -qJD(5) * t705 - t771 * t577 + t766 * t986;
t962 = t989 * t771;
t718 = t941 * t767;
t890 = t766 * t718 + t771 * t719;
t961 = -t718 * t871 + t719 * t872 + t766 * t989 + t771 * t988;
t876 = qJD(2) * t773;
t958 = t767 * t876 + t768 * t873;
t850 = t767 * t874;
t852 = t773 * t867;
t957 = -t850 + t852;
t756 = t767 * qJ(4);
t952 = t772 * pkin(3) + pkin(2) + t756;
t939 = pkin(8) * t687;
t949 = t600 * t734 + t939;
t603 = t667 * t766 + t668 * t771;
t657 = t693 * t768;
t810 = t665 * t766 + t666 * t771;
t944 = -g(1) * t603 - g(2) * t810 - g(3) * t657 - t567 * t809 - t791;
t943 = t691 ^ 2;
t940 = pkin(3) * t687;
t929 = pkin(8) * qJD(3);
t927 = t512 * t765;
t590 = -t722 + t624;
t926 = t590 * t734;
t925 = t605 * t767;
t923 = t624 * t734;
t920 = t689 * t691;
t919 = t689 * t734;
t918 = t691 * t734;
t917 = t691 * t772;
t915 = t765 * t766;
t914 = t766 * t770;
t908 = t770 * t512;
t904 = t963 + t967;
t903 = t964 + t982;
t805 = t766 * t772 - t912;
t620 = t770 * t693 - t765 * t805;
t902 = -qJD(6) * t620 - t765 * t895 + t770 * t896;
t621 = -t693 * t765 - t770 * t805;
t901 = qJD(6) * t621 + t765 * t896 + t770 * t895;
t900 = pkin(5) * t895 + t894;
t893 = -pkin(3) * t992 + t978;
t891 = t767 * t704 - t803 * t873;
t630 = t691 * pkin(3) + t689 * qJ(4);
t889 = (g(1) * t906 + g(2) * t909) * t768;
t887 = qJ(4) * t852 + qJD(4) * t910;
t760 = t768 ^ 2;
t884 = -t773 ^ 2 + t760;
t879 = qJD(2) * t691;
t877 = qJD(2) * t768;
t860 = qJ(4) * t877 + t891;
t745 = pkin(7) * t863;
t661 = -qJDD(2) * pkin(2) + pkin(7) * t847 + t745;
t855 = t734 * t878;
t854 = t734 * t867;
t851 = t734 * t875;
t821 = -qJD(3) * t737 + t704 * t772 + t803 * t875;
t545 = pkin(9) * t850 + qJD(2) * t785 - t821;
t546 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t910 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t767) * t773 + t860;
t840 = t771 * t545 - t546 * t766;
t837 = t771 * t610 - t622 * t766;
t834 = t771 * t718 - t719 * t766;
t833 = -t772 * t803 - t736;
t832 = t865 ^ 2;
t682 = t772 * pkin(4) + t952;
t828 = -pkin(7) - t862;
t579 = -pkin(4) * t691 - t630;
t827 = t768 * t859;
t825 = -g(1) * t665 + g(2) * t667;
t824 = g(1) * t666 - g(2) * t668;
t822 = g(1) * t769 - g(2) * t774;
t595 = -pkin(10) * t693 + t890;
t820 = -pkin(5) * t881 + pkin(10) * t896 + qJD(5) * t890 + qJD(6) * t595 - t766 * t988 + t962;
t594 = pkin(10) * t805 + t834;
t819 = pkin(10) * t895 - qJD(6) * t594 + t961;
t818 = t771 * t977;
t817 = qJD(3) * t716 - t939;
t500 = t765 * t509 + t908;
t536 = pkin(5) * t773 - pkin(10) * t657 + t837;
t537 = -pkin(10) * t656 + t897;
t816 = t536 * t770 - t537 * t765;
t815 = t536 * t765 + t537 * t770;
t588 = pkin(3) * t734 + t866;
t814 = t588 * t772 - t590 * t767;
t581 = t770 * t656 + t657 * t765;
t582 = -t656 * t765 + t657 * t770;
t700 = -pkin(5) + t844;
t808 = t700 * t765 + t705 * t770;
t801 = -t734 * t929 + t930;
t532 = t606 * pkin(3) + t605 * qJ(4) - t691 * qJD(4) + t661;
t798 = -0.2e1 * pkin(1) * t864 - pkin(7) * qJDD(2);
t797 = t687 * t767 - t734 * t873;
t796 = t687 * t772 + t851;
t793 = -t532 - t801;
t790 = t766 * t545 + t771 * t546 + t610 * t871 - t622 * t872;
t777 = qJD(1) ^ 2;
t788 = pkin(1) * t777 + t823;
t732 = qJ(4) * t910;
t633 = t768 * t828 + t732;
t518 = -pkin(4) * t606 - t532;
t787 = -t773 * t823 - t931;
t776 = qJD(2) ^ 2;
t784 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t776 + t822;
t783 = g(1) * t667 + g(2) * t665 + g(3) * t913 - t829;
t779 = t600 * t691 + qJDD(4) - t783;
t778 = g(1) * t668 + g(2) * t666 + g(3) * t910 - t623 * t734 - t792;
t558 = (-t772 * t942 - t756) * t874 + t828 * t876 + t887;
t649 = -t732 + (pkin(3) * t767 + pkin(7)) * t768;
t635 = t759 - t833;
t632 = pkin(5) * t693 + t682;
t629 = qJD(1) * t827 - t916;
t628 = -pkin(7) * t857 + t892;
t578 = pkin(5) * t656 + t633;
t573 = pkin(3) * t958 + pkin(7) * t876 + qJ(4) * t850 - t887;
t566 = qJD(2) * t827 - t821;
t562 = -t605 - t919;
t561 = -qJD(4) * t773 + (-t768 * t867 - t773 * t875) * pkin(7) + t860;
t560 = t768 * t805 * t953 + qJD(2) * t794;
t559 = qJD(5) * t657 + t766 * t957 - t771 * t958;
t541 = -pkin(5) * t617 + t579;
t533 = t804 - t940;
t523 = pkin(5) * t559 + t558;
t508 = qJD(6) * t582 + t770 * t559 + t560 * t765;
t507 = -qJD(6) * t581 - t559 * t765 + t560 * t770;
t505 = pkin(5) * t529 + t518;
t504 = -pkin(10) * t559 + t790;
t503 = -pkin(5) * t877 - pkin(10) * t560 - qJD(5) * t897 + t840;
t499 = t509 * t770 - t927;
t1 = [(-t528 * t656 - t529 * t657 - t559 * t617 - t560 * t809) * MDP(23) + (-t529 * t773 + t559 * t865 + t656 * t675 + t809 * t877) * MDP(25) + (-t675 * t773 + t865 * t877) * MDP(26) + (t528 * t773 - t560 * t865 - t617 * t877 - t657 * t675) * MDP(24) + ((t503 * t770 - t504 * t765) * t977 - t816 * t662 + t843 * t773 - t499 * t877 + t523 * t987 + t578 * t502 + t505 * t581 + t540 * t508 + g(1) * t811 - g(2) * t586 + (-t500 * t773 - t815 * t977) * qJD(6)) * MDP(34) + (-t502 * t773 - t508 * t977 + t581 * t662 + t877 * t987) * MDP(32) + (-t501 * t581 - t502 * t582 - t507 * t987 + t508 * t972) * MDP(30) + (-t561 * t734 - t573 * t691 + t605 * t649 + t634 * t687 + (-t600 * t867 - t531) * t773 + (qJD(2) * t590 - t532 * t772 + t600 * t875) * t768 - t825) * MDP(20) + (-t662 * t773 - t877 * t977) * MDP(33) + (t501 * t582 - t507 * t972) * MDP(29) + (t501 * t773 + t507 * t977 - t582 * t662 + t877 * t972) * MDP(31) + (-t840 * t865 - t837 * t675 + t842 * t773 - t521 * t877 + t558 * t809 + t633 * t529 + t518 * t656 + t567 * t559 + g(1) * t810 - g(2) * t603 + (-t522 * t773 + t865 * t897) * qJD(5)) * MDP(27) + (g(1) * t956 - g(2) * t602 + t518 * t657 + t522 * t877 + t633 * t528 + t558 * t617 + t567 * t560 + t675 * t897 + t773 * t791 + t790 * t865) * MDP(28) + (-(qJD(6) * t816 + t503 * t765 + t504 * t770) * t977 + t815 * t662 - t830 * t773 + t500 * t877 - t523 * t972 + t578 * t501 + t505 * t582 + t540 * t507 + g(1) * t955 - g(2) * t585) * MDP(35) + (t531 * t634 + t590 * t561 + t532 * t649 + t600 * t573 + t533 * t635 + t588 * t566 - g(1) * (-pkin(3) * t666 - qJ(4) * t665) - g(2) * (pkin(3) * t668 + qJ(4) * t667) + (-g(1) * pkin(7) - g(2) * t803) * t774 + (-g(2) * pkin(7) + g(1) * t803) * t769) * MDP(21) + 0.2e1 * (t751 * t768 - t864 * t884) * MDP(5) + (t891 * t734 - t886 * t687 + (t716 * t867 + (-t851 + t879) * pkin(7) + t792) * t773 + (-t716 * t875 - t624 * qJD(2) + t661 * t772 + (-t605 - t854) * pkin(7)) * t768 + t825) * MDP(17) + ((t606 + t855) * t773 + (-qJD(2) * t689 - t797) * t768) * MDP(14) + (-t561 * t689 + t566 * t691 - t605 * t635 - t606 * t634 + t814 * t876 + (-t531 * t767 + t533 * t772 + (-t588 * t767 - t590 * t772) * qJD(3) + t822) * t768) * MDP(19) + (-t687 * t773 - t734 * t877) * MDP(15) + (t566 * t734 + t573 * t689 + t606 * t649 - t635 * t687 + (t600 * t878 + t533) * t773 + (-qJD(2) * t588 + t532 * t767 + t600 * t873) * t768 + t824) * MDP(18) + ((t605 - t854) * t773 + (t796 + t879) * t768) * MDP(13) + t822 * MDP(2) + t823 * MDP(3) + ((-t689 * t772 - t691 * t767) * t876 + (t925 - t606 * t772 + (t689 * t767 - t917) * qJD(3)) * t768) * MDP(12) + (qJDD(1) * t760 + 0.2e1 * t768 * t847) * MDP(4) + (-t821 * t734 + t833 * t687 + ((pkin(7) * t689 + t716 * t767) * qJD(2) + t829) * t773 + (t716 * t873 + t623 * qJD(2) + t661 * t767 + (t606 - t855) * pkin(7)) * t768 + t824) * MDP(16) + (t768 * t798 + t773 * t784) * MDP(9) + (-t768 * t784 + t773 * t798) * MDP(10) + qJDD(1) * MDP(1) + (-t605 * t910 + t691 * t957) * MDP(11) + (qJDD(2) * t768 + t773 * t776) * MDP(6) + (qJDD(2) * t773 - t768 * t776) * MDP(7) + (t528 * t657 + t560 * t617) * MDP(22); (-MDP(4) * t768 * t773 + MDP(5) * t884) * t777 + (-t528 * t693 + t529 * t805 - t617 * t895 - t809 * t896) * MDP(23) + (-t528 * t805 + t617 * t896) * MDP(22) + (t675 * t805 - t865 * t896) * MDP(24) + (t675 * t693 + t865 * t895) * MDP(25) + ((t689 * t768 - t734 * t911) * qJD(1) + t796) * MDP(14) + (-t501 * t620 - t502 * t621 + t901 * t972 - t902 * t987) * MDP(30) + (t734 * MDP(15) + t588 * MDP(18) - t590 * MDP(20) + t617 * MDP(24) - MDP(25) * t809 - MDP(26) * t865 + t521 * MDP(27) - t522 * MDP(28) - MDP(31) * t972 - MDP(32) * t987 + MDP(33) * t977 + t499 * MDP(34) - t500 * MDP(35)) * t881 + (-(t594 * t770 - t595 * t765) * t662 + t632 * t502 + t505 * t620 + (t765 * t819 - t770 * t820) * t977 + t900 * t987 + t901 * t540 - t945 * t806) * MDP(34) + (-t588 * t629 - t590 * t628 + t893 * t600 + (qJD(3) * t814 + t531 * t772 + t533 * t767 + t787) * pkin(8) + (-t532 - t945) * t952) * MDP(21) + (-t834 * t675 + t682 * t529 - g(3) * t794 - (-t719 * t871 + (-qJD(5) * t718 + t988) * t766 - t962) * t865 + t894 * t809 + t895 * t567 + (t518 + t951) * t693) * MDP(27) + (-t621 * t662 + t902 * t977) * MDP(31) + (t620 * t662 - t901 * t977) * MDP(32) + (t501 * t621 - t902 * t972) * MDP(29) + ((t594 * t765 + t595 * t770) * t662 + t632 * t501 + t505 * t621 + (t765 * t820 + t770 * t819) * t977 - t900 * t972 + t902 * t540 + t945 * (t754 * t772 - t755 * t767)) * MDP(35) + (-t605 * t952 + t628 * t734 - t893 * t691 + t949 * t772 + (t793 + t951) * t767) * MDP(20) + (-t606 * t952 - t629 * t734 + t893 * t689 - t767 * t949 + t793 * t772 + t889) * MDP(18) + (t682 * t528 + t896 * t567 + t894 * t617 + t890 * t675 - t865 * t961 + (-t518 + t945) * t805) * MDP(28) + ((-t605 + t919) * t772 + (-t606 + t918) * t767) * MDP(12) + ((-t691 * t768 + t734 * t907) * qJD(1) + t797) * MDP(13) + (-pkin(2) * t606 + t817 * t767 + (-t930 - t661 + (t701 + t929) * t734) * t772 + (-t716 * t911 - t623 * t768 + (-t689 * t773 + t734 * t913) * pkin(7)) * qJD(1) + t889) * MDP(16) + (t768 * t788 - t745 - t930) * MDP(9) + (t931 + (-pkin(7) * qJDD(1) + t788) * t773) * MDP(10) + MDP(7) * t751 + MDP(6) * t863 + (-t734 * t917 - t925) * MDP(11) + (t628 * t689 - t629 * t691 + (t531 - t734 * t588 + (qJD(3) * t691 - t606) * pkin(8)) * t772 + (t533 + t926 + (qJD(3) * t689 - t605) * pkin(8)) * t767 + t787) * MDP(19) + qJDD(2) * MDP(8) + (pkin(2) * t605 - t670 * t734 + t817 * t772 + (-t716 * t907 + t624 * t768 + (-t691 * t773 + t734 * t910) * pkin(7)) * qJD(1) + (t661 + t801 - t951) * t767) * MDP(17); (t617 * t865 + t529) * MDP(25) + (-t579 * t617 + t705 * t675 - t865 * t964 + t944) * MDP(28) + (t531 * qJ(4) - t533 * pkin(3) - t600 * t630 - t588 * t624 - g(1) * (-pkin(3) * t667 + qJ(4) * t668) - g(2) * (-pkin(3) * t665 + qJ(4) * t666) - g(3) * (-pkin(3) * t913 + t732) + t866 * t590) * MDP(21) + (-t630 * t689 - t779 - t923 + 0.2e1 * t940) * MDP(18) + (t808 * t662 + t541 * t972 + ((-qJD(6) * t700 + t903) * t770 + (qJD(6) * t705 - t904) * t765) * t977 - t974) * MDP(35) + (-t579 * t809 - t844 * t675 - t865 * t963 - t985) * MDP(27) + (-t689 ^ 2 + t943) * MDP(12) + (-t691 * t716 + t783 - t923) * MDP(16) + MDP(11) * t920 + (-t606 - t918) * MDP(14) + (-t600 * t689 + t630 * t691 + 0.2e1 * t674 - 0.2e1 * t720 - t778) * MDP(20) + (t689 * t716 + t778) * MDP(17) + (pkin(3) * t605 - qJ(4) * t606 + (t590 - t624) * t691 + (t588 - t866) * t689) * MDP(19) + (-(t700 * t770 - t705 * t765) * t662 - t541 * t987 + (t765 * t903 + t770 * t904) * t977 + (-t808 * t977 + t500) * qJD(6) + t983) * MDP(34) + t687 * MDP(15) + t562 * MDP(13) - t998; (-t687 + t920) * MDP(18) + t562 * MDP(19) + (-t734 ^ 2 - t943) * MDP(20) + (t779 + t926 - t940) * MDP(21) + (-t675 * t771 - t691 * t809 - t766 * t832) * MDP(27) + (-t617 * t691 + t766 * t675 - t771 * t832) * MDP(28) + (-(t770 * t771 - t915) * t662 - t691 * t987 + (-t765 * t818 - t914 * t977) * t977) * MDP(34) + ((t765 * t771 + t914) * t662 + t691 * t972 + (-t770 * t818 + t915 * t977) * t977) * MDP(35); (-t838 + (-qJD(5) - t865) * t617) * MDP(25) + (-t522 * t865 + t985) * MDP(27) + (-t521 * t865 - t944) * MDP(28) + (-(-t511 * t765 - t908) * t977 - t500 * qJD(6) + (-t617 * t987 - t662 * t770 - t870 * t977) * pkin(5) - t983) * MDP(34) + ((t511 * t770 - t927) * t977 + (t617 * t972 + t662 * t765 - t869 * t977) * pkin(5) + t974) * MDP(35) + t998; (t861 + t981) * MDP(31) + (-t841 - t980) * MDP(32) + (t500 * t977 - t983) * MDP(34) + (t499 * t977 + t974) * MDP(35) + (-MDP(31) * t924 + MDP(32) * t972 - MDP(34) * t500) * qJD(6) + t999;];
tau  = t1;
