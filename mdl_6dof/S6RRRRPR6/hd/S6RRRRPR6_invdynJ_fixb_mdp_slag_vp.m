% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR6
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
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:41
% EndTime: 2019-03-09 22:24:05
% DurationCPUTime: 18.03s
% Computational Cost: add. (15616->707), mult. (35685->921), div. (0->0), fcn. (26570->16), ass. (0->307)
t831 = cos(qJ(3));
t827 = sin(qJ(2));
t832 = cos(qJ(2));
t945 = t831 * t832;
t859 = pkin(3) * t827 - pkin(9) * t945;
t974 = pkin(8) + pkin(9);
t902 = qJD(3) * t974;
t875 = pkin(2) * t827 - pkin(8) * t832;
t764 = t875 * qJD(1);
t826 = sin(qJ(3));
t920 = qJD(1) * t827;
t901 = t826 * t920;
t925 = pkin(7) * t901 + t831 * t764;
t1019 = qJD(1) * t859 + t831 * t902 + t925;
t741 = t826 * t764;
t949 = t827 * t831;
t950 = t826 * t832;
t1014 = t741 + (-pkin(7) * t949 - pkin(9) * t950) * qJD(1) + t826 * t902;
t829 = cos(qJ(6));
t907 = t831 * qJD(2);
t755 = t901 - t907;
t916 = qJD(2) * t826;
t757 = t831 * t920 + t916;
t825 = sin(qJ(4));
t830 = cos(qJ(4));
t691 = t830 * t755 + t757 * t825;
t822 = sin(pkin(11));
t823 = cos(pkin(11));
t863 = t755 * t825 - t830 * t757;
t985 = -t691 * t823 + t822 * t863;
t946 = t829 * t985;
t824 = sin(qJ(6));
t984 = t691 * t822 + t823 * t863;
t959 = t984 * t824;
t1006 = t946 + t959;
t812 = t832 * qJDD(1);
t906 = qJD(1) * qJD(2);
t979 = -t827 * t906 + t812;
t752 = qJDD(3) - t979;
t747 = qJDD(4) + t752;
t733 = qJDD(6) + t747;
t998 = t824 * t985 - t829 * t984;
t1018 = t733 * MDP(31) + (-t1006 ^ 2 + t998 ^ 2) * MDP(28) - t1006 * MDP(27) * t998;
t773 = -pkin(2) * t832 - pkin(8) * t827 - pkin(1);
t748 = t773 * qJD(1);
t919 = qJD(1) * t832;
t807 = pkin(7) * t919;
t779 = qJD(2) * pkin(8) + t807;
t699 = t831 * t748 - t779 * t826;
t666 = -pkin(9) * t757 + t699;
t795 = -qJD(3) + t919;
t657 = -pkin(3) * t795 + t666;
t700 = t748 * t826 + t779 * t831;
t667 = -pkin(9) * t755 + t700;
t662 = t825 * t667;
t610 = t830 * t657 - t662;
t993 = qJ(5) * t863;
t584 = t610 + t993;
t785 = -qJD(4) + t795;
t577 = -pkin(4) * t785 + t584;
t664 = t830 * t667;
t611 = t657 * t825 + t664;
t994 = qJ(5) * t691;
t585 = t611 - t994;
t953 = t823 * t585;
t537 = t822 * t577 + t953;
t982 = pkin(10) * t985;
t529 = t537 + t982;
t908 = qJD(6) * t824;
t528 = t529 * t908;
t778 = -qJD(2) * pkin(2) + pkin(7) * t920;
t710 = pkin(3) * t755 + t778;
t652 = pkin(4) * t691 + qJD(5) + t710;
t590 = -pkin(5) * t985 + t652;
t821 = qJ(3) + qJ(4);
t801 = pkin(11) + qJ(6) + t821;
t792 = sin(t801);
t793 = cos(t801);
t833 = cos(qJ(1));
t828 = sin(qJ(1));
t948 = t828 * t832;
t707 = t792 * t833 - t793 * t948;
t944 = t832 * t833;
t709 = t792 * t828 + t793 * t944;
t966 = g(3) * t827;
t1017 = g(1) * t709 - g(2) * t707 - t590 * t1006 + t793 * t966 + t528;
t892 = t832 * t906;
t905 = qJDD(1) * t827;
t912 = qJD(3) * t827;
t989 = -qJD(1) * t912 + qJDD(2);
t681 = qJD(3) * t907 + (t892 + t905) * t831 + t989 * t826;
t682 = t826 * (qJD(2) * (qJD(3) + t919) + t905) - t989 * t831;
t909 = qJD(4) * t830;
t910 = qJD(4) * t825;
t601 = t830 * t681 - t825 * t682 - t755 * t909 - t757 * t910;
t838 = qJD(4) * t863 - t681 * t825 - t830 * t682;
t554 = -t601 * t822 + t823 * t838;
t555 = t601 * t823 + t822 * t838;
t903 = qJD(6) * t946 + t824 * t554 + t829 * t555;
t517 = t908 * t984 + t903;
t889 = -t829 * t554 + t555 * t824;
t518 = qJD(6) * t998 + t889;
t776 = -qJD(6) + t785;
t960 = t1006 * t776;
t961 = t998 * t776;
t1016 = (t517 + t960) * MDP(29) + t747 * MDP(22) + (-t691 ^ 2 + t863 ^ 2) * MDP(19) + (-t691 * t785 + t601) * MDP(20) + (t785 * t863 + t838) * MDP(21) - t691 * MDP(18) * t863 + (-t518 - t961) * MDP(30) + t1018;
t758 = t825 * t826 - t830 * t831;
t850 = t832 * t758;
t977 = qJD(3) + qJD(4);
t1008 = -qJD(1) * t850 + t977 * t758;
t759 = t825 * t831 + t826 * t830;
t932 = (-t919 + t977) * t759;
t767 = t875 * qJD(2);
t704 = qJD(1) * t767 + qJDD(1) * t773;
t696 = t831 * t704;
t731 = pkin(7) * t979 + qJDD(2) * pkin(8);
t586 = pkin(3) * t752 - pkin(9) * t681 - qJD(3) * t700 - t731 * t826 + t696;
t911 = qJD(3) * t831;
t913 = qJD(3) * t826;
t849 = t826 * t704 + t831 * t731 + t748 * t911 - t779 * t913;
t597 = -pkin(9) * t682 + t849;
t888 = t830 * t586 - t825 * t597;
t839 = -qJD(4) * t611 + t888;
t523 = pkin(4) * t747 - qJ(5) * t601 + qJD(5) * t863 + t839;
t879 = -t825 * t586 - t830 * t597 - t657 * t909 + t667 * t910;
t525 = qJ(5) * t838 - qJD(5) * t691 - t879;
t511 = t823 * t523 - t525 * t822;
t509 = pkin(5) * t747 - pkin(10) * t555 + t511;
t512 = t822 * t523 + t823 * t525;
t510 = pkin(10) * t554 + t512;
t1012 = -t824 * t509 - t829 * t510 + t1017;
t1009 = t1019 * t830;
t780 = t974 * t826;
t781 = t974 * t831;
t1007 = -t1014 * t830 - t1019 * t825 - t780 * t909 - t781 * t910;
t706 = t792 * t948 + t793 * t833;
t708 = -t792 * t944 + t793 * t828;
t890 = t829 * t509 - t824 * t510;
t999 = -g(1) * t708 + g(2) * t706 - t590 * t998 + t792 * t966 + t890;
t1004 = pkin(5) * t984;
t995 = pkin(10) * t984;
t926 = -t825 * t780 + t830 * t781;
t1003 = -pkin(4) * t920 + t1008 * qJ(5) - qJD(4) * t926 - qJD(5) * t759 + t1014 * t825 - t1009;
t1002 = -t932 * qJ(5) - qJD(5) * t758 + t1007;
t805 = pkin(7) * t905;
t732 = -qJDD(2) * pkin(2) + pkin(7) * t892 + t805;
t873 = g(1) * t833 + g(2) * t828;
t965 = g(3) * t832;
t845 = t827 * t873 - t965;
t1001 = qJD(3) * pkin(8) * t795 - t732 + t845;
t973 = pkin(3) * t826;
t876 = pkin(3) * t913 - t919 * t973 - t807;
t1000 = pkin(4) * t932 + t876;
t997 = pkin(4) * t863;
t992 = t1008 * t822 - t823 * t932;
t991 = -t1008 * t823 - t822 * t932;
t914 = qJD(2) * t832;
t899 = t826 * t914;
t990 = t827 * t911 + t899;
t814 = sin(t821);
t815 = cos(t821);
t717 = t814 * t833 - t815 * t948;
t719 = t814 * t828 + t815 * t944;
t987 = g(1) * t719 - g(2) * t717 + t691 * t710 + t815 * t966 + t879;
t716 = t814 * t948 + t815 * t833;
t718 = -t814 * t944 + t815 * t828;
t976 = -g(1) * t718 + g(2) * t716 + t814 * t966;
t986 = t710 * t863 + t839 + t976;
t725 = t759 * t827;
t939 = -t1002 * t822 + t1003 * t823;
t938 = t1002 * t823 + t1003 * t822;
t754 = t831 * t773;
t971 = pkin(7) * t826;
t698 = -pkin(9) * t949 + t754 + (-pkin(3) - t971) * t832;
t797 = pkin(7) * t945;
t924 = t826 * t773 + t797;
t951 = t826 * t827;
t705 = -pkin(9) * t951 + t924;
t934 = t825 * t698 + t830 * t705;
t886 = -t666 * t825 - t664;
t591 = t886 + t994;
t936 = t830 * t666 - t662;
t592 = t936 + t993;
t952 = t823 * t825;
t962 = pkin(3) * qJD(4);
t930 = -t823 * t591 + t592 * t822 + (-t822 * t830 - t952) * t962;
t954 = t822 * t825;
t929 = -t822 * t591 - t823 * t592 + (t823 * t830 - t954) * t962;
t978 = -t826 * t912 + t832 * t907;
t972 = pkin(4) * t822;
t817 = t831 * pkin(3);
t964 = pkin(2) + t817;
t771 = pkin(4) * t814 + t973;
t963 = pkin(7) + t771;
t958 = t681 * t826;
t957 = t755 * t795;
t956 = t757 * t795;
t955 = t757 * t831;
t578 = t822 * t585;
t536 = t823 * t577 - t578;
t527 = -pkin(5) * t785 + t536 + t995;
t947 = t829 * t527;
t943 = -t930 - t982;
t942 = t929 - t995;
t650 = -qJD(2) * t850 - t725 * t977;
t726 = t758 * t827;
t915 = qJD(2) * t827;
t927 = t831 * t767 + t915 * t971;
t645 = t859 * qJD(2) + (-t797 + (pkin(9) * t827 - t773) * t826) * qJD(3) + t927;
t928 = t826 * t767 + t773 * t911;
t649 = -t990 * pkin(9) + (-t827 * t907 - t832 * t913) * pkin(7) + t928;
t887 = t830 * t645 - t649 * t825;
t540 = pkin(4) * t915 - qJ(5) * t650 - qJD(4) * t934 + qJD(5) * t726 + t887;
t651 = -t910 * t951 + (t949 * t977 + t899) * t830 + t978 * t825;
t848 = t825 * t645 + t830 * t649 + t698 * t909 - t705 * t910;
t544 = -qJ(5) * t651 - qJD(5) * t725 + t848;
t520 = t822 * t540 + t823 * t544;
t688 = -t758 * t823 - t759 * t822;
t689 = -t758 * t822 + t759 * t823;
t866 = t829 * t688 - t689 * t824;
t941 = qJD(6) * t866 + t824 * t992 + t829 * t991;
t641 = t688 * t824 + t689 * t829;
t940 = qJD(6) * t641 + t824 * t991 - t829 * t992;
t543 = t823 * t584 - t578;
t937 = -pkin(5) * t992 + t1000;
t884 = t830 * t698 - t705 * t825;
t623 = -pkin(4) * t832 + qJ(5) * t726 + t884;
t628 = -qJ(5) * t725 + t934;
t566 = t822 * t623 + t823 * t628;
t883 = -t830 * t780 - t781 * t825;
t675 = -qJ(5) * t759 + t883;
t676 = -qJ(5) * t758 + t926;
t626 = t822 * t675 + t823 * t676;
t772 = pkin(4) * t815 + t817;
t768 = pkin(3) * t951 + t827 * pkin(7);
t819 = t827 ^ 2;
t923 = -t832 ^ 2 + t819;
t918 = qJD(2) * t755;
t917 = qJD(2) * t757;
t711 = pkin(3) * t990 + pkin(7) * t914;
t900 = t795 * t907;
t897 = t795 * t913;
t896 = t795 * t911;
t519 = t823 * t540 - t544 * t822;
t542 = -t584 * t822 - t953;
t565 = t823 * t623 - t628 * t822;
t625 = t823 * t675 - t676 * t822;
t882 = -qJD(3) * t748 - t731;
t881 = qJD(6) * t527 + t510;
t878 = pkin(4) * t758 - t964;
t803 = pkin(3) * t830 + pkin(4);
t728 = -pkin(3) * t954 + t823 * t803;
t877 = pkin(4) * t725 + t768;
t874 = pkin(3) * t757 - t997;
t872 = g(1) * t828 - g(2) * t833;
t871 = t779 * t911 - t696;
t600 = pkin(10) * t688 + t626;
t870 = pkin(5) * t920 + pkin(10) * t991 + qJD(6) * t600 - t939;
t599 = -pkin(10) * t689 + t625;
t869 = pkin(10) * t992 + qJD(6) * t599 + t938;
t868 = -pkin(8) * t752 + qJD(3) * t778;
t514 = t824 * t527 + t829 * t529;
t668 = -t725 * t823 + t726 * t822;
t669 = -t725 * t822 - t726 * t823;
t867 = t829 * t668 - t669 * t824;
t620 = t668 * t824 + t669 * t829;
t763 = pkin(2) + t772;
t818 = -qJ(5) - t974;
t862 = t763 * t832 - t818 * t827;
t860 = pkin(4) * t651 + t711;
t798 = pkin(4) * t823 + pkin(5);
t858 = t798 * t824 + t829 * t972;
t857 = t798 * t829 - t824 * t972;
t855 = pkin(1) + t862;
t854 = -0.2e1 * pkin(1) * t906 - pkin(7) * qJDD(2);
t853 = t752 * t826 - t896;
t852 = t752 * t831 + t897;
t835 = qJD(1) ^ 2;
t846 = pkin(1) * t835 + t873;
t655 = pkin(3) * t682 + t732;
t834 = qJD(2) ^ 2;
t842 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t834 + t872;
t567 = -pkin(4) * t838 + qJDD(5) + t655;
t739 = t826 * t828 + t831 * t944;
t738 = -t826 * t944 + t828 * t831;
t737 = t826 * t833 - t828 * t945;
t736 = t826 * t948 + t831 * t833;
t729 = pkin(3) * t952 + t803 * t822;
t721 = pkin(5) + t728;
t656 = -pkin(5) * t688 + t878;
t633 = -pkin(5) * t668 + t877;
t603 = -t997 - t1004;
t598 = t874 - t1004;
t595 = t650 * t823 - t651 * t822;
t594 = -t650 * t822 - t651 * t823;
t562 = -pkin(5) * t594 + t860;
t557 = pkin(10) * t668 + t566;
t556 = -pkin(5) * t832 - pkin(10) * t669 + t565;
t533 = qJD(6) * t620 - t829 * t594 + t595 * t824;
t532 = qJD(6) * t867 + t594 * t824 + t595 * t829;
t531 = t543 + t995;
t530 = t542 - t982;
t526 = -pkin(5) * t554 + t567;
t516 = pkin(10) * t594 + t520;
t515 = pkin(5) * t915 - pkin(10) * t595 + t519;
t513 = -t529 * t824 + t947;
t1 = [(-g(1) * t716 - g(2) * t718 + t768 * t601 - t611 * t915 + t710 * t650 - t655 * t726 - t711 * t863 - t747 * t934 + t785 * t848 - t832 * t879) * MDP(24) + (-t601 * t832 - t650 * t785 - t726 * t747 - t863 * t915) * MDP(20) + (-t601 * t726 - t650 * t863) * MDP(18) + (-(-t773 * t913 + t927) * t795 + t754 * t752 - g(1) * t737 - g(2) * t739 + ((t896 + t918) * pkin(7) + (-pkin(7) * t752 + qJD(2) * t778 - t882) * t826 + t871) * t832 + (pkin(7) * t682 + qJD(2) * t699 + t732 * t826 + t778 * t911) * t827) * MDP(16) + (t928 * t795 - t924 * t752 - g(1) * t736 - g(2) * t738 + (t778 * t907 + (-t897 + t917) * pkin(7) + t849) * t832 + (-t778 * t913 - t700 * qJD(2) + t732 * t831 + (t681 - t900) * pkin(7)) * t827) * MDP(17) + (t512 * t566 + t537 * t520 + t511 * t565 + t536 * t519 + t567 * t877 + t652 * t860 + (-g(1) * t963 - g(2) * t855) * t833 + (g(1) * t855 - g(2) * t963) * t828) * MDP(26) + ((-t755 * t831 - t757 * t826) * t914 + (-t958 - t682 * t831 + (t755 * t826 - t955) * qJD(3)) * t827) * MDP(12) + (qJDD(1) * t819 + 0.2e1 * t827 * t892) * MDP(4) + (t827 * t854 + t832 * t842) * MDP(9) + (-t827 * t842 + t832 * t854) * MDP(10) + 0.2e1 * (t812 * t827 - t906 * t923) * MDP(5) + ((-t681 - t900) * t832 + (t852 + t917) * t827) * MDP(13) + ((t795 * t916 + t682) * t832 + (-t853 - t918) * t827) * MDP(14) + (-t514 * t915 - g(1) * t706 - g(2) * t708 + t633 * t517 + t526 * t620 - t528 * t832 + t590 * t532 + t562 * t998 + ((-qJD(6) * t557 + t515) * t776 - t556 * t733 + t509 * t832) * t824 + ((qJD(6) * t556 + t516) * t776 - t557 * t733 + t881 * t832) * t829) * MDP(33) + (-t517 * t832 - t532 * t776 + t620 * t733 + t915 * t998) * MDP(29) + (t517 * t620 + t532 * t998) * MDP(27) + (t681 * t949 + t757 * t978) * MDP(11) + (t1006 * t915 + t518 * t832 + t533 * t776 + t733 * t867) * MDP(30) + (-(t515 * t829 - t516 * t824) * t776 + (t556 * t829 - t557 * t824) * t733 - t890 * t832 + t513 * t915 - t562 * t1006 + t633 * t518 - t526 * t867 + t590 * t533 - g(1) * t707 - g(2) * t709 + (-(-t556 * t824 - t557 * t829) * t776 + t514 * t832) * qJD(6)) * MDP(32) + (t1006 * t532 + t517 * t867 - t518 * t620 - t533 * t998) * MDP(28) + (-t752 * t832 - t795 * t915) * MDP(15) + (-t747 * t832 - t785 * t915) * MDP(22) + (-t733 * t832 - t776 * t915) * MDP(31) + t872 * MDP(2) + t873 * MDP(3) + (qJDD(2) * t827 + t832 * t834) * MDP(6) + (qJDD(2) * t832 - t827 * t834) * MDP(7) + qJDD(1) * MDP(1) + (-t601 * t725 - t650 * t691 + t651 * t863 - t726 * t838) * MDP(19) + (t651 * t785 - t691 * t915 - t725 * t747 - t832 * t838) * MDP(21) + (-t887 * t785 + t884 * t747 - t888 * t832 + t610 * t915 + t711 * t691 - t768 * t838 + t655 * t725 + t710 * t651 - g(1) * t717 - g(2) * t719 + (t611 * t832 + t785 * t934) * qJD(4)) * MDP(23) + (-t511 * t669 + t512 * t668 + t519 * t984 + t520 * t985 - t536 * t595 + t537 * t594 + t554 * t566 - t555 * t565 + t827 * t872) * MDP(25); (t733 * t866 + t776 * t940) * MDP(30) + (t827 * t846 - t805 - t965) * MDP(9) + (t966 + (-pkin(7) * qJDD(1) + t846) * t832) * MDP(10) + (-t795 * t955 + t958) * MDP(11) + ((t681 + t957) * t831 + (-t682 + t956) * t826) * MDP(12) + ((t755 * t827 - t795 * t950) * qJD(1) + t852) * MDP(14) + ((-t757 * t827 + t795 * t945) * qJD(1) + t853) * MDP(13) + (t512 * t626 + t511 * t625 + t567 * t878 - g(3) * t862 + t1000 * t652 + t938 * t537 + t939 * t536 + t873 * (t763 * t827 + t818 * t832)) * MDP(26) + (t517 * t641 + t941 * t998) * MDP(27) + (-(t599 * t824 + t600 * t829) * t733 + t656 * t517 + t526 * t641 + (-t824 * t870 + t829 * t869) * t776 + t941 * t590 + t937 * t998 - t845 * t792) * MDP(33) + (-t747 * t758 + t785 * t932) * MDP(21) + ((t599 * t829 - t600 * t824) * t733 + t656 * t518 - t526 * t866 + (t824 * t869 + t829 * t870) * t776 + t940 * t590 - t937 * t1006 + t845 * t793) * MDP(32) + (t795 * MDP(15) + MDP(20) * t863 + t691 * MDP(21) + t785 * MDP(22) - t610 * MDP(23) + t611 * MDP(24) - MDP(29) * t998 - MDP(30) * t1006 + t776 * MDP(31) - t513 * MDP(32) + t514 * MDP(33)) * t920 + (t1006 * t941 + t517 * t866 - t518 * t641 - t940 * t998) * MDP(28) + (t883 * t747 + t964 * t838 + t655 * t758 + (t781 * t909 + (-qJD(4) * t780 - t1014) * t825 + t1009) * t785 + t932 * t710 + t876 * t691 + t845 * t815) * MDP(23) + (-MDP(4) * t827 * t832 + MDP(5) * t923) * t835 + (-pkin(2) * t682 + t925 * t795 + t868 * t826 + (-t699 * t827 + (-pkin(7) * t755 - t778 * t826) * t832) * qJD(1) + t1001 * t831) * MDP(16) + (-pkin(2) * t681 - t741 * t795 + t868 * t831 + (-t778 * t945 + t700 * t827 + (-t757 * t832 + t795 * t949) * pkin(7)) * qJD(1) - t1001 * t826) * MDP(17) + (t641 * t733 - t776 * t941) * MDP(29) + MDP(6) * t905 + MDP(7) * t812 + (t1008 * t863 + t601 * t759) * MDP(18) + (t1007 * t785 - t1008 * t710 - t964 * t601 + t655 * t759 - t926 * t747 - t845 * t814 - t876 * t863) * MDP(24) + (t1008 * t785 + t747 * t759) * MDP(20) + (t1008 * t691 - t601 * t758 + t759 * t838 + t863 * t932) * MDP(19) + qJDD(2) * MDP(8) + (-t511 * t689 + t512 * t688 - t536 * t991 + t537 * t992 + t554 * t626 - t555 * t625 - t873 * t832 + t938 * t985 + t939 * t984 - t966) * MDP(25); (t554 * t729 - t555 * t728 + (t536 + t929) * t985 + (-t537 + t930) * t984) * MDP(25) + (t886 * t785 + (-t691 * t757 + t747 * t830 + t785 * t910) * pkin(3) + t986) * MDP(23) + (-t936 * t785 + (-t747 * t825 + t757 * t863 + t785 * t909) * pkin(3) + t987) * MDP(24) + (t512 * t729 + t511 * t728 - t652 * t874 - g(1) * (-t771 * t944 + t772 * t828) - g(2) * (-t771 * t948 - t772 * t833) + t771 * t966 + t929 * t537 + t930 * t536) * MDP(26) + (-g(1) * t738 + g(2) * t736 - t700 * t795 - t757 * t778 + (t882 + t966) * t826 - t871) * MDP(16) + (t681 - t957) * MDP(13) + (g(1) * t739 - g(2) * t737 + g(3) * t949 - t699 * t795 + t755 * t778 - t849) * MDP(17) + (-t598 * t998 + (-t721 * t733 - t509 + (-qJD(6) * t729 - t943) * t776) * t824 + (-t729 * t733 + (qJD(6) * t721 + t942) * t776 - t881) * t829 + t1017) * MDP(33) + ((t721 * t829 - t729 * t824) * t733 + t598 * t1006 + (t942 * t824 + t943 * t829) * t776 + (-(-t721 * t824 - t729 * t829) * t776 - t514) * qJD(6) + t999) * MDP(32) + (-t682 - t956) * MDP(14) + t757 * t755 * MDP(11) + t752 * MDP(15) + (-t755 ^ 2 + t757 ^ 2) * MDP(12) + t1016; (-t611 * t785 + t986) * MDP(23) + (-t610 * t785 + t987) * MDP(24) + ((t554 * t822 - t555 * t823) * pkin(4) + (t536 - t543) * t985 + (-t537 - t542) * t984) * MDP(25) + (-t536 * t542 - t537 * t543 + (t511 * t823 + t512 * t822 + t652 * t863 + t976) * pkin(4)) * MDP(26) + (t857 * t733 + (t530 * t829 - t531 * t824) * t776 + t603 * t1006 + (t776 * t858 - t514) * qJD(6) + t999) * MDP(32) + (-t858 * t733 - (t530 * t824 + t531 * t829) * t776 - t603 * t998 + (t776 * t857 - t947) * qJD(6) + t1012) * MDP(33) + t1016; (-t984 ^ 2 - t985 ^ 2) * MDP(25) + (-t536 * t984 - t537 * t985 + t567 - t845) * MDP(26) + (t518 - t961) * MDP(32) + (t517 - t960) * MDP(33); (t903 + t960) * MDP(29) + (-t889 - t961) * MDP(30) + (-t514 * t776 + t999) * MDP(32) + (-t513 * t776 + t1012) * MDP(33) + (MDP(29) * t959 - MDP(30) * t998 - MDP(32) * t514 - MDP(33) * t947) * qJD(6) + t1018;];
tau  = t1;
