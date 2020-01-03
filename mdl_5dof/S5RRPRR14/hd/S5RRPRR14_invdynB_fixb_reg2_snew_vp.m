% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRR14_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:39:01
% EndTime: 2019-12-31 20:39:32
% DurationCPUTime: 27.86s
% Computational Cost: add. (192002->816), mult. (441786->1317), div. (0->0), fcn. (353012->12), ass. (0->595)
t875 = cos(pkin(5));
t864 = t875 * qJD(1) + qJD(2);
t872 = sin(pkin(10));
t874 = cos(pkin(10));
t873 = sin(pkin(5));
t878 = sin(qJ(2));
t977 = qJD(1) * t878;
t959 = t873 * t977;
t830 = t874 * t864 - t872 * t959;
t831 = t872 * t864 + t874 * t959;
t877 = sin(qJ(4));
t881 = cos(qJ(4));
t791 = -t881 * t830 + t877 * t831;
t882 = cos(qJ(2));
t976 = qJD(1) * t882;
t958 = t873 * t976;
t856 = -qJD(4) + t958;
t780 = t791 * t856;
t973 = qJDD(1) * t878;
t838 = (qJD(2) * t976 + t973) * t873;
t953 = t875 * qJDD(1) + qJDD(2);
t804 = t874 * t838 + t872 * t953;
t954 = t872 * t838 - t874 * t953;
t896 = t791 * qJD(4) - t881 * t804 + t877 * t954;
t1031 = t896 - t780;
t793 = t877 * t830 + t881 * t831;
t876 = sin(qJ(5));
t880 = cos(qJ(5));
t764 = t876 * t793 + t880 * t856;
t766 = t880 * t793 - t876 * t856;
t708 = t766 * t764;
t956 = t877 * t804 + t881 * t954;
t719 = -t793 * qJD(4) - t956;
t897 = qJDD(5) - t719;
t1024 = -t708 + t897;
t1030 = t1024 * t876;
t1029 = t1024 * t880;
t869 = t873 ^ 2;
t884 = qJD(1) ^ 2;
t1005 = t869 * t884;
t1010 = t830 * t831;
t972 = t873 * qJDD(1);
t839 = -qJD(2) * t959 + t882 * t972;
t898 = -t839 + t1010;
t1028 = t872 * t898;
t1027 = t874 * t898;
t1011 = t793 * t791;
t833 = -qJDD(4) + t839;
t894 = -t833 - t1011;
t1026 = t877 * t894;
t1025 = t881 * t894;
t815 = t830 * t958;
t776 = -t815 - t804;
t816 = t831 * t958;
t772 = t816 + t954;
t846 = t864 * t959;
t810 = t839 - t846;
t847 = t864 * t958;
t807 = t847 + t838;
t788 = qJD(5) + t791;
t728 = t788 * t764;
t962 = t764 * qJD(5) + t876 * t833 + t880 * t896;
t631 = t962 + t728;
t957 = t880 * t833 - t876 * t896;
t626 = (qJD(5) - t788) * t766 + t957;
t689 = (qJD(4) + t856) * t793 + t956;
t762 = t764 ^ 2;
t763 = t766 ^ 2;
t786 = t788 ^ 2;
t789 = t791 ^ 2;
t790 = t793 ^ 2;
t1022 = t830 ^ 2;
t827 = t831 ^ 2;
t1021 = t856 ^ 2;
t1020 = t864 ^ 2;
t1019 = 2 * qJD(3);
t1018 = pkin(2) * t878;
t1017 = pkin(2) * t882;
t1016 = pkin(4) * t877;
t1015 = pkin(7) * t873;
t1014 = t875 * g(3);
t1013 = t788 * t876;
t1012 = t788 * t880;
t1009 = t838 * qJ(3);
t1008 = t856 * t877;
t1007 = t856 * t881;
t1006 = t873 * t1005;
t879 = sin(qJ(1));
t883 = cos(qJ(1));
t859 = t883 * g(1) + t879 * g(2);
t834 = -t884 * pkin(1) + pkin(7) * t972 - t859;
t951 = -qJ(3) * t878 - t1017;
t978 = qJD(1) * t873;
t835 = t951 * t978;
t1001 = t873 * t878;
t858 = t879 * g(1) - t883 * g(2);
t893 = qJDD(1) * pkin(1) + t884 * t1015 + t858;
t891 = t875 * t893;
t890 = -g(3) * t1001 + t878 * t891;
t889 = qJ(3) * t953 + t890;
t887 = -t1020 * pkin(2) + (t835 * t978 + t834) * t882 + t889;
t950 = -qJ(3) * t882 + t1018;
t979 = qJD(1) * t864;
t888 = -t839 * pkin(2) - t1014 - t1009 + (t950 * t979 - t893) * t873;
t664 = t830 * t1019 + t872 * t888 + t874 * t887;
t805 = -pkin(3) * t958 - t831 * pkin(8);
t637 = -t1022 * pkin(3) - pkin(8) * t954 + t805 * t958 + t664;
t971 = t831 * t1019;
t663 = t872 * t887 - t874 * t888 + t971;
t886 = t898 * pkin(3) + t776 * pkin(8) - t663;
t560 = t877 * t637 - t881 * t886;
t561 = t881 * t637 + t877 * t886;
t505 = -t881 * t560 + t877 * t561;
t1004 = t872 * t505;
t955 = t878 * t834 - t882 * t891;
t743 = qJDD(3) - t953 * pkin(2) - t1020 * qJ(3) + (g(3) * t882 + t835 * t977) * t873 + t955;
t1003 = t872 * t743;
t777 = t839 + t1010;
t1002 = t872 * t777;
t1000 = t873 * t882;
t999 = t874 * t505;
t998 = t874 * t743;
t997 = t874 * t777;
t733 = t791 * pkin(4) - t793 * pkin(9);
t536 = t833 * pkin(4) - t1021 * pkin(9) + t793 * t733 + t560;
t996 = t876 * t536;
t654 = t708 + t897;
t995 = t876 * t654;
t669 = pkin(3) * t954 - t1022 * pkin(8) + t831 * t805 + t743;
t994 = t877 * t669;
t725 = t833 - t1011;
t993 = t877 * t725;
t817 = t873 * t893 + t1014;
t992 = t878 * t817;
t988 = t878 * t882;
t967 = t869 * t988;
t855 = t884 * t967;
t836 = -t855 + t953;
t991 = t878 * t836;
t837 = t855 + t953;
t990 = t878 * t837;
t989 = t878 * t839;
t987 = t880 * t536;
t986 = t880 * t654;
t985 = t881 * t669;
t984 = t881 * t725;
t983 = t882 * t817;
t982 = t882 * t834;
t981 = t882 * t836;
t980 = t882 * t837;
t537 = -t1021 * pkin(4) - t833 * pkin(9) - t791 * t733 + t561;
t569 = t1031 * pkin(9) + (-t856 * t793 - t719) * pkin(4) + t669;
t500 = t880 * t537 + t876 * t569;
t970 = t882 * t1011;
t969 = t882 * t1010;
t870 = t878 ^ 2;
t968 = t870 * t1005;
t966 = t877 * t708;
t965 = t878 * t1011;
t964 = t878 * t1010;
t963 = t881 * t708;
t961 = -pkin(4) * t881 - pkin(3);
t960 = t864 * t978;
t499 = t876 * t537 - t880 * t569;
t506 = t877 * t560 + t881 * t561;
t593 = t872 * t663 + t874 * t664;
t819 = -t879 * t858 - t883 * t859;
t854 = t883 * qJDD(1) - t879 * t884;
t952 = -pkin(6) * t854 - t879 * g(3);
t826 = -t968 - t1020;
t794 = -t878 * t826 - t981;
t949 = pkin(7) * t794 - t992;
t871 = t882 ^ 2;
t863 = t871 * t1005;
t843 = -t863 - t1020;
t798 = t882 * t843 - t990;
t948 = pkin(7) * t798 + t983;
t455 = t876 * t499 + t880 * t500;
t442 = t877 * t455 - t881 * t536;
t443 = t881 * t455 + t877 * t536;
t412 = -t872 * t442 + t874 * t443;
t454 = -t880 * t499 + t876 * t500;
t947 = t412 * t878 - t454 * t882;
t457 = t874 * t506 - t1004;
t946 = t457 * t878 - t669 * t882;
t630 = -t728 + t962;
t557 = -t626 * t880 - t876 * t630;
t671 = t762 + t763;
t532 = t877 * t557 + t881 * t671;
t533 = t881 * t557 - t877 * t671;
t475 = -t872 * t532 + t874 * t533;
t555 = -t626 * t876 + t880 * t630;
t945 = t475 * t878 - t555 * t882;
t627 = (-qJD(5) - t788) * t766 - t957;
t556 = t880 * t627 + t631 * t876;
t707 = -t763 + t762;
t538 = t877 * t556 + t881 * t707;
t539 = t881 * t556 - t877 * t707;
t482 = -t872 * t538 + t874 * t539;
t554 = -t876 * t627 + t631 * t880;
t944 = t482 * t878 + t554 * t882;
t687 = -t786 - t762;
t588 = t880 * t687 - t1030;
t540 = t877 * t588 + t881 * t627;
t541 = t881 * t588 - t877 * t627;
t488 = -t872 * t540 + t874 * t541;
t587 = t876 * t687 + t1029;
t943 = t488 * t878 - t587 * t882;
t697 = -t763 - t786;
t595 = -t876 * t697 - t986;
t542 = t877 * t595 + t881 * t631;
t543 = t881 * t595 - t877 * t631;
t491 = -t872 * t542 + t874 * t543;
t594 = t880 * t697 - t995;
t942 = t491 * t878 - t594 * t882;
t724 = -t763 + t786;
t605 = -t876 * t724 + t1029;
t550 = t877 * t605 + t881 * t630;
t552 = t881 * t605 - t877 * t630;
t497 = -t872 * t550 + t874 * t552;
t603 = -t880 * t724 - t1030;
t941 = t497 * t878 + t603 * t882;
t723 = t762 - t786;
t606 = t880 * t723 - t995;
t551 = t877 * t606 + t881 * t626;
t553 = t881 * t606 - t877 * t626;
t498 = -t872 * t551 + t874 * t553;
t604 = -t876 * t723 - t986;
t940 = t498 * t878 + t604 * t882;
t665 = -t766 * qJD(5) - t957;
t620 = t764 * t1012 - t876 * t665;
t583 = t877 * t620 + t963;
t585 = t881 * t620 - t966;
t525 = -t872 * t583 + t874 * t585;
t619 = -t764 * t1013 - t880 * t665;
t939 = t525 * t878 + t619 * t882;
t622 = -t766 * t1013 - t880 * t962;
t584 = t877 * t622 - t963;
t586 = t881 * t622 + t966;
t526 = -t872 * t584 + t874 * t586;
t621 = -t766 * t1012 + t876 * t962;
t938 = t526 * t878 + t621 * t882;
t657 = (-t764 * t880 + t766 * t876) * t788;
t611 = t877 * t657 - t881 * t897;
t612 = t881 * t657 + t877 * t897;
t545 = -t872 * t611 + t874 * t612;
t656 = (t764 * t876 + t766 * t880) * t788;
t937 = t545 * t878 + t656 * t882;
t688 = (qJD(4) - t856) * t793 + t956;
t613 = -t1031 * t881 - t877 * t688;
t615 = t1031 * t877 - t881 * t688;
t548 = -t872 * t613 + t874 * t615;
t734 = -t790 + t789;
t936 = t548 * t878 + t734 * t882;
t693 = t780 + t896;
t614 = -t689 * t877 + t881 * t693;
t616 = -t689 * t881 - t877 * t693;
t549 = -t872 * t614 + t874 * t616;
t702 = -t789 - t790;
t935 = t549 * t878 - t702 * t882;
t731 = -t1021 - t789;
t658 = t877 * t731 + t1025;
t659 = t881 * t731 - t1026;
t582 = -t872 * t658 + t874 * t659;
t934 = t582 * t878 - t688 * t882;
t933 = t593 * t878 - t743 * t882;
t754 = -t790 - t1021;
t674 = t881 * t754 + t993;
t675 = -t877 * t754 + t984;
t598 = -t872 * t674 + t874 * t675;
t932 = t1031 * t882 + t598 * t878;
t771 = -t790 + t1021;
t683 = t881 * t771 + t1026;
t685 = -t877 * t771 + t1025;
t609 = -t872 * t683 + t874 * t685;
t931 = t609 * t878 + t693 * t882;
t770 = t789 - t1021;
t684 = t877 * t770 - t984;
t686 = t881 * t770 + t993;
t610 = -t872 * t684 + t874 * t686;
t930 = t610 * t878 + t689 * t882;
t716 = (t791 * t877 + t793 * t881) * t856;
t717 = (t791 * t881 - t793 * t877) * t856;
t652 = -t872 * t716 + t874 * t717;
t929 = t652 * t878 + t833 * t882;
t592 = -t874 * t663 + t872 * t664;
t773 = t816 - t954;
t775 = -t815 + t804;
t711 = t874 * t773 - t872 * t775;
t795 = -t827 + t1022;
t928 = t711 * t878 + t795 * t882;
t712 = -t772 * t874 - t872 * t776;
t767 = t827 + t1022;
t927 = t712 * t878 + t767 * t882;
t787 = -t863 - t1022;
t722 = t874 * t787 - t1028;
t926 = t722 * t878 + t773 * t882;
t813 = -t827 - t863;
t738 = -t872 * t813 + t997;
t925 = t738 * t878 - t775 * t882;
t812 = -t827 + t863;
t739 = -t872 * t812 + t1027;
t924 = t739 * t878 + t776 * t882;
t811 = -t863 + t1022;
t740 = t874 * t811 + t1002;
t923 = t740 * t878 + t772 * t882;
t784 = g(3) * t1000 + t955;
t785 = t890 + t982;
t922 = -t882 * t784 + t878 * t785;
t732 = t878 * t784 + t882 * t785;
t921 = t807 * t882 + t810 * t878;
t808 = -t847 + t838;
t809 = t839 + t846;
t920 = -t808 * t882 + t809 * t878;
t919 = t826 * t882 - t991;
t842 = t863 - t1020;
t918 = t842 * t878 + t981;
t841 = -t968 + t1020;
t917 = t841 * t882 + t990;
t916 = t843 * t878 + t980;
t818 = t883 * t858 - t879 * t859;
t915 = -t875 * t884 + t979;
t914 = t873 * t953;
t678 = -t791 * t1008 + t881 * t719;
t679 = -t791 * t1007 - t877 * t719;
t601 = -t872 * t678 + t874 * t679;
t913 = t601 * t878 + t970;
t680 = -t793 * t1007 - t877 * t896;
t681 = t793 * t1008 - t881 * t896;
t602 = -t872 * t680 + t874 * t681;
t912 = t602 * t878 - t970;
t759 = t815 * t874 + t872 * t954;
t911 = t759 * t878 - t969;
t761 = t874 * t804 + t816 * t872;
t910 = t761 * t878 + t969;
t401 = pkin(8) * t443 + (-pkin(9) * t877 + t961) * t454;
t407 = -pkin(8) * t442 + (-pkin(9) * t881 + t1016) * t454;
t411 = t874 * t442 + t872 * t443;
t387 = -qJ(3) * t411 - t872 * t401 + t874 * t407;
t397 = -pkin(2) * t411 - pkin(3) * t442 + pkin(4) * t536 - pkin(9) * t455;
t403 = t882 * t412 + t878 * t454;
t909 = pkin(7) * t403 + t387 * t878 + t397 * t882;
t448 = -pkin(9) * t555 - t454;
t418 = pkin(8) * t533 + t877 * t448 + t961 * t555;
t428 = -pkin(8) * t532 + t555 * t1016 + t881 * t448;
t474 = t874 * t532 + t872 * t533;
t399 = -qJ(3) * t474 - t872 * t418 + t874 * t428;
t416 = -pkin(2) * t474 - pkin(3) * t532 - pkin(4) * t671 - pkin(9) * t557 - t455;
t460 = t882 * t475 + t878 * t555;
t908 = pkin(7) * t460 + t399 * t878 + t416 * t882;
t467 = -pkin(4) * t587 + t499;
t509 = -pkin(9) * t587 + t996;
t431 = -pkin(3) * t587 + pkin(8) * t541 + t881 * t467 + t877 * t509;
t438 = -pkin(8) * t540 - t877 * t467 + t881 * t509;
t487 = t874 * t540 + t872 * t541;
t405 = -qJ(3) * t487 - t872 * t431 + t874 * t438;
t437 = -pkin(2) * t487 - pkin(3) * t540 - pkin(4) * t627 - pkin(9) * t588 + t987;
t464 = t882 * t488 + t878 * t587;
t907 = pkin(7) * t464 + t405 * t878 + t437 * t882;
t468 = -pkin(4) * t594 + t500;
t510 = -pkin(9) * t594 + t987;
t432 = -pkin(3) * t594 + pkin(8) * t543 + t881 * t468 + t877 * t510;
t441 = -pkin(8) * t542 - t877 * t468 + t881 * t510;
t490 = t874 * t542 + t872 * t543;
t406 = -qJ(3) * t490 - t872 * t432 + t874 * t441;
t440 = -pkin(2) * t490 - pkin(3) * t542 - pkin(4) * t631 - pkin(9) * t595 - t996;
t465 = t882 * t491 + t878 * t594;
t906 = pkin(7) * t465 + t406 * t878 + t440 * t882;
t456 = t872 * t506 + t999;
t494 = -pkin(3) * t669 + pkin(8) * t506;
t424 = -pkin(8) * t999 - qJ(3) * t456 - t872 * t494;
t436 = -pkin(2) * t456 - pkin(3) * t505;
t453 = t882 * t457 + t878 * t669;
t905 = pkin(7) * t453 + t424 * t878 + t436 * t882;
t469 = -pkin(3) * t702 + pkin(8) * t616 + t506;
t478 = -pkin(8) * t614 - t505;
t547 = t874 * t614 + t872 * t616;
t433 = -qJ(3) * t547 - t872 * t469 + t874 * t478;
t516 = -pkin(2) * t547 - pkin(3) * t614;
t529 = t882 * t549 + t878 * t702;
t904 = pkin(7) * t529 + t433 * t878 + t516 * t882;
t563 = -pkin(3) * t688 + pkin(8) * t659 - t985;
t581 = t874 * t658 + t872 * t659;
t589 = -pkin(8) * t658 + t994;
t489 = -qJ(3) * t581 - t872 * t563 + t874 * t589;
t503 = -pkin(2) * t581 - pkin(3) * t658 + t560;
t558 = t882 * t582 + t878 * t688;
t903 = pkin(7) * t558 + t489 * t878 + t503 * t882;
t564 = pkin(3) * t1031 + pkin(8) * t675 + t994;
t596 = -pkin(8) * t674 + t985;
t597 = t874 * t674 + t872 * t675;
t501 = -qJ(3) * t597 - t872 * t564 + t874 * t596;
t511 = -pkin(2) * t597 - pkin(3) * t674 + t561;
t565 = -t1031 * t878 + t882 * t598;
t902 = pkin(7) * t565 + t501 * t878 + t511 * t882;
t721 = t872 * t787 + t1027;
t625 = t872 * (t835 * t958 + t889 + t982) - t874 * (-qJ(3) * t847 - t1009 - t817) + t971 + (-t872 * t1020 + t874 * t810 - t721) * pkin(2);
t662 = -qJ(3) * t721 + t1003;
t682 = t882 * t722 - t878 * t773;
t901 = pkin(7) * t682 + t625 * t882 + t662 * t878;
t735 = t874 * t813 + t1002;
t633 = -pkin(2) * t735 + t664;
t672 = -qJ(3) * t735 + t998;
t698 = t882 * t738 + t878 * t775;
t900 = pkin(7) * t698 + t633 * t882 + t672 * t878;
t755 = t878 * t808 + t882 * t809;
t899 = pkin(7) * t755 + t732;
t710 = -t772 * t872 + t874 * t776;
t568 = -qJ(3) * t710 - t592;
t670 = t882 * t712 - t878 * t767;
t895 = pkin(7) * t670 - t710 * t1017 + t568 * t878;
t575 = t882 * t593 + t878 * t743;
t892 = pkin(7) * t575 + t592 * t951;
t853 = t879 * qJDD(1) + t883 * t884;
t845 = -t863 - t968;
t844 = -t863 + t968;
t840 = -pkin(6) * t853 + t883 * g(3);
t824 = t875 * t882 * t839;
t823 = t839 * t1000;
t814 = (t870 + t871) * t960;
t806 = (t973 + (qJD(2) + t864) * t976) * t873;
t800 = t882 * t838 - t870 * t960;
t799 = -t871 * t960 - t989;
t797 = t882 * t842 - t991;
t796 = -t878 * t841 + t980;
t783 = (t882 * t1006 + t807 * t875) * t878;
t782 = t824 + (-t875 * t960 - t1006) * t988;
t769 = (-t830 * t874 - t831 * t872) * t958;
t768 = (-t830 * t872 + t831 * t874) * t958;
t760 = t872 * t804 - t816 * t874;
t758 = t815 * t872 - t874 * t954;
t756 = -t878 * t807 + t882 * t810;
t753 = t873 * t810 + t875 * t916;
t752 = -t873 * t809 + t875 * t918;
t751 = -t873 * t808 + t875 * t917;
t750 = -t875 * t810 + t873 * t916;
t749 = -t873 * t806 + t875 * t919;
t748 = t875 * t806 + t873 * t919;
t747 = t882 * t769 - t989;
t746 = -t873 * t844 + t875 * t921;
t745 = -t873 * t845 + t875 * t920;
t744 = t875 * t845 + t873 * t920;
t737 = t872 * t811 - t997;
t736 = t874 * t812 + t1028;
t730 = t882 * t761 - t964;
t729 = t882 * t759 + t964;
t714 = -t879 * t753 + t883 * t798;
t713 = t883 * t753 + t879 * t798;
t709 = t872 * t773 + t874 * t775;
t706 = -t879 * t749 + t883 * t794;
t705 = t883 * t749 + t879 * t794;
t704 = t873 * t817 + t875 * t922;
t703 = -t875 * t817 + t873 * t922;
t701 = t875 * t878 * t769 - t873 * t768 + t824;
t700 = t882 * t740 - t878 * t772;
t699 = t882 * t739 - t878 * t776;
t696 = -t879 * t745 + t883 * t755;
t695 = t883 * t745 + t879 * t755;
t694 = t882 * t711 - t878 * t795;
t677 = -t873 * t760 + t875 * t910;
t676 = -t873 * t758 + t875 * t911;
t673 = -t992 + (-t750 * t873 - t753 * t875) * pkin(7);
t668 = -t983 + (-t748 * t873 - t749 * t875) * pkin(7);
t667 = -pkin(1) * t750 + t873 * t784 + t875 * t948;
t660 = -pkin(1) * t748 + t873 * t785 + t875 * t949;
t651 = t874 * t716 + t872 * t717;
t650 = -pkin(2) * t775 + qJ(3) * t738 + t1003;
t649 = t875 * pkin(7) * t732 - pkin(1) * t703;
t648 = -t879 * t704 + t883 * t732;
t647 = t883 * t704 + t879 * t732;
t646 = -t873 * t737 + t875 * t923;
t645 = -t873 * t736 + t875 * t924;
t644 = -t873 * t735 + t875 * t925;
t643 = t875 * t735 + t873 * t925;
t642 = pkin(2) * t773 + qJ(3) * t722 - t998;
t641 = t882 * t652 - t878 * t833;
t640 = -pkin(1) * t744 + t875 * t899;
t639 = (-t703 * t873 - t704 * t875) * pkin(7);
t638 = (-t744 * t873 - t745 * t875) * pkin(7) - t922;
t636 = -t873 * t721 + t875 * t926;
t635 = t875 * t721 + t873 * t926;
t632 = -t873 * t709 + t875 * t928;
t624 = -t873 * t710 + t875 * t927;
t623 = t875 * t710 + t873 * t927;
t608 = t874 * t684 + t872 * t686;
t607 = t874 * t683 + t872 * t685;
t600 = t874 * t680 + t872 * t681;
t599 = t874 * t678 + t872 * t679;
t591 = -t879 * t644 + t883 * t698;
t590 = t883 * t644 + t879 * t698;
t580 = t882 * t602 + t965;
t579 = t882 * t601 - t965;
t578 = -t879 * t636 + t883 * t682;
t577 = t883 * t636 + t879 * t682;
t576 = -pkin(2) * t743 + qJ(3) * t593;
t574 = -t879 * t624 + t883 * t670;
t573 = t883 * t624 + t879 * t670;
t572 = -t873 * t651 + t875 * t929;
t571 = t882 * t610 - t878 * t689;
t570 = t882 * t609 - t878 * t693;
t562 = pkin(2) * t767 + qJ(3) * t712 + t593;
t546 = t874 * t613 + t872 * t615;
t544 = t874 * t611 + t872 * t612;
t534 = t882 * t548 - t878 * t734;
t531 = -t873 * t600 + t875 * t912;
t530 = -t873 * t599 + t875 * t913;
t528 = -t873 * t608 + t875 * t930;
t527 = -t873 * t607 + t875 * t931;
t524 = t874 * t584 + t872 * t586;
t523 = t874 * t583 + t872 * t585;
t522 = -t873 * t597 + t875 * t932;
t521 = t875 * t597 + t873 * t932;
t520 = t882 * t545 - t878 * t656;
t519 = -t873 * t592 + t875 * t933;
t518 = t875 * t592 + t873 * t933;
t517 = -t878 * t633 + t882 * t672 + (-t643 * t873 - t644 * t875) * pkin(7);
t515 = -t873 * t581 + t875 * t934;
t514 = t875 * t581 + t873 * t934;
t513 = -t878 * t625 + t882 * t662 + (-t635 * t873 - t636 * t875) * pkin(7);
t512 = -pkin(1) * t643 - t873 * t650 + t875 * t900;
t508 = t882 * t526 - t878 * t621;
t507 = t882 * t525 - t878 * t619;
t504 = -pkin(1) * t635 - t873 * t642 + t875 * t901;
t502 = t710 * t1018 + t882 * t568 + (-t623 * t873 - t624 * t875) * pkin(7);
t496 = t874 * t551 + t872 * t553;
t495 = t874 * t550 + t872 * t552;
t493 = -t879 * t519 + t883 * t575;
t492 = t883 * t519 + t879 * t575;
t486 = pkin(2) * t1031 + qJ(3) * t598 + t874 * t564 + t872 * t596;
t485 = -t879 * t522 + t883 * t565;
t484 = t883 * t522 + t879 * t565;
t483 = -t873 * t546 + t875 * t936;
t481 = t874 * t538 + t872 * t539;
t480 = -t873 * t547 + t875 * t935;
t479 = t875 * t547 + t873 * t935;
t477 = -pkin(2) * t688 + qJ(3) * t582 + t874 * t563 + t872 * t589;
t476 = -t873 * t544 + t875 * t937;
t473 = -t879 * t515 + t883 * t558;
t472 = t883 * t515 + t879 * t558;
t471 = t882 * t498 - t878 * t604;
t470 = t882 * t497 - t878 * t603;
t466 = -pkin(1) * t623 - t873 * t562 + t875 * t895;
t463 = -t873 * t524 + t875 * t938;
t462 = -t873 * t523 + t875 * t939;
t461 = t882 * t482 - t878 * t554;
t459 = -t879 * t480 + t883 * t529;
t458 = t883 * t480 + t879 * t529;
t452 = t950 * t592 + (-t518 * t873 - t519 * t875) * pkin(7);
t451 = -pkin(1) * t518 - t873 * t576 + t875 * t892;
t450 = -t873 * t496 + t875 * t940;
t449 = -t873 * t495 + t875 * t941;
t447 = -t873 * t490 + t875 * t942;
t446 = t875 * t490 + t873 * t942;
t445 = -t873 * t487 + t875 * t943;
t444 = t875 * t487 + t873 * t943;
t439 = -t873 * t481 + t875 * t944;
t435 = -t873 * t474 + t875 * t945;
t434 = t875 * t474 + t873 * t945;
t430 = -pkin(2) * t702 + qJ(3) * t549 + t874 * t469 + t872 * t478;
t429 = t882 * t501 - t878 * t511 + (-t521 * t873 - t522 * t875) * pkin(7);
t427 = t882 * t489 - t878 * t503 + (-t514 * t873 - t515 * t875) * pkin(7);
t426 = -t879 * t447 + t883 * t465;
t425 = t883 * t447 + t879 * t465;
t423 = -t873 * t456 + t875 * t946;
t422 = t875 * t456 + t873 * t946;
t421 = -t879 * t445 + t883 * t464;
t420 = t883 * t445 + t879 * t464;
t419 = -pkin(2) * t669 - pkin(8) * t1004 + qJ(3) * t457 + t874 * t494;
t417 = -pkin(1) * t521 - t873 * t486 + t875 * t902;
t415 = -pkin(1) * t514 - t873 * t477 + t875 * t903;
t414 = -t879 * t435 + t883 * t460;
t413 = t883 * t435 + t879 * t460;
t410 = -t879 * t423 + t883 * t453;
t409 = t883 * t423 + t879 * t453;
t408 = t882 * t433 - t878 * t516 + (-t479 * t873 - t480 * t875) * pkin(7);
t404 = -pkin(2) * t594 + qJ(3) * t491 + t874 * t432 + t872 * t441;
t402 = -pkin(2) * t587 + qJ(3) * t488 + t874 * t431 + t872 * t438;
t400 = -pkin(1) * t479 - t873 * t430 + t875 * t904;
t398 = -pkin(2) * t555 + qJ(3) * t475 + t874 * t418 + t872 * t428;
t396 = -t873 * t411 + t875 * t947;
t395 = t875 * t411 + t873 * t947;
t394 = t882 * t424 - t878 * t436 + (-t422 * t873 - t423 * t875) * pkin(7);
t393 = t882 * t406 - t878 * t440 + (-t446 * t873 - t447 * t875) * pkin(7);
t392 = t882 * t405 - t878 * t437 + (-t444 * t873 - t445 * t875) * pkin(7);
t391 = -pkin(1) * t422 - t873 * t419 + t875 * t905;
t390 = t882 * t399 - t878 * t416 + (-t434 * t873 - t435 * t875) * pkin(7);
t389 = -t879 * t396 + t883 * t403;
t388 = t883 * t396 + t879 * t403;
t386 = -pkin(1) * t446 - t873 * t404 + t875 * t906;
t385 = -pkin(1) * t444 - t873 * t402 + t875 * t907;
t384 = -pkin(2) * t454 + qJ(3) * t412 + t874 * t401 + t872 * t407;
t383 = -pkin(1) * t434 - t873 * t398 + t875 * t908;
t382 = t882 * t387 - t878 * t397 + (-t395 * t873 - t396 * t875) * pkin(7);
t381 = -pkin(1) * t395 - t873 * t384 + t875 * t909;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t853, -t854, 0, t819, 0, 0, 0, 0, 0, 0, t714, t706, t696, t648, 0, 0, 0, 0, 0, 0, t578, t591, t574, t493, 0, 0, 0, 0, 0, 0, t473, t485, t459, t410, 0, 0, 0, 0, 0, 0, t421, t426, t414, t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t854, -t853, 0, t818, 0, 0, 0, 0, 0, 0, t713, t705, t695, t647, 0, 0, 0, 0, 0, 0, t577, t590, t573, t492, 0, 0, 0, 0, 0, 0, t472, t484, t458, t409, 0, 0, 0, 0, 0, 0, t420, t425, t413, t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t750, t748, t744, t703, 0, 0, 0, 0, 0, 0, t635, t643, t623, t518, 0, 0, 0, 0, 0, 0, t514, t521, t479, t422, 0, 0, 0, 0, 0, 0, t444, t446, t434, t395; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t854, 0, -t853, 0, t952, -t840, -t818, -pkin(6) * t818, -t879 * t783 + t883 * t800, -t879 * t746 + t883 * t756, -t879 * t751 + t883 * t796, -t879 * t782 + t883 * t799, -t879 * t752 + t883 * t797, t883 * t814 + t879 * t914, -pkin(6) * t713 - t879 * t667 + t883 * t673, -pkin(6) * t705 - t879 * t660 + t883 * t668, -pkin(6) * t695 + t883 * t638 - t879 * t640, -pkin(6) * t647 + t883 * t639 - t879 * t649, -t879 * t677 + t883 * t730, -t879 * t632 + t883 * t694, -t879 * t645 + t883 * t699, -t879 * t676 + t883 * t729, -t879 * t646 + t883 * t700, -t879 * t701 + t883 * t747, -pkin(6) * t577 - t879 * t504 + t883 * t513, -pkin(6) * t590 - t879 * t512 + t883 * t517, -pkin(6) * t573 - t879 * t466 + t883 * t502, -pkin(6) * t492 - t879 * t451 + t883 * t452, -t879 * t531 + t883 * t580, -t879 * t483 + t883 * t534, -t879 * t527 + t883 * t570, -t879 * t530 + t883 * t579, -t879 * t528 + t883 * t571, -t879 * t572 + t883 * t641, -pkin(6) * t472 - t879 * t415 + t883 * t427, -pkin(6) * t484 - t879 * t417 + t883 * t429, -pkin(6) * t458 - t879 * t400 + t883 * t408, -pkin(6) * t409 - t879 * t391 + t883 * t394, -t879 * t463 + t883 * t508, -t879 * t439 + t883 * t461, -t879 * t449 + t883 * t470, -t879 * t462 + t883 * t507, -t879 * t450 + t883 * t471, -t879 * t476 + t883 * t520, -pkin(6) * t420 - t879 * t385 + t883 * t392, -pkin(6) * t425 - t879 * t386 + t883 * t393, -pkin(6) * t413 - t879 * t383 + t883 * t390, -pkin(6) * t388 - t879 * t381 + t883 * t382; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t853, 0, t854, 0, t840, t952, t819, pkin(6) * t819, t883 * t783 + t879 * t800, t883 * t746 + t879 * t756, t883 * t751 + t879 * t796, t883 * t782 + t879 * t799, t883 * t752 + t879 * t797, t879 * t814 - t883 * t914, pkin(6) * t714 + t883 * t667 + t879 * t673, pkin(6) * t706 + t883 * t660 + t879 * t668, pkin(6) * t696 + t879 * t638 + t883 * t640, pkin(6) * t648 + t879 * t639 + t883 * t649, t883 * t677 + t879 * t730, t883 * t632 + t879 * t694, t883 * t645 + t879 * t699, t883 * t676 + t879 * t729, t883 * t646 + t879 * t700, t883 * t701 + t879 * t747, pkin(6) * t578 + t883 * t504 + t879 * t513, pkin(6) * t591 + t883 * t512 + t879 * t517, pkin(6) * t574 + t883 * t466 + t879 * t502, pkin(6) * t493 + t883 * t451 + t879 * t452, t883 * t531 + t879 * t580, t883 * t483 + t879 * t534, t883 * t527 + t879 * t570, t883 * t530 + t879 * t579, t883 * t528 + t879 * t571, t883 * t572 + t879 * t641, pkin(6) * t473 + t883 * t415 + t879 * t427, pkin(6) * t485 + t883 * t417 + t879 * t429, pkin(6) * t459 + t883 * t400 + t879 * t408, pkin(6) * t410 + t883 * t391 + t879 * t394, t883 * t463 + t879 * t508, t883 * t439 + t879 * t461, t883 * t449 + t879 * t470, t883 * t462 + t879 * t507, t883 * t450 + t879 * t471, t883 * t476 + t879 * t520, pkin(6) * t421 + t883 * t385 + t879 * t392, pkin(6) * t426 + t883 * t386 + t879 * t393, pkin(6) * t414 + t883 * t383 + t879 * t390, pkin(6) * t389 + t883 * t381 + t879 * t382; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t858, t859, 0, 0, (t869 * t882 * t915 + t838 * t873) * t878, t875 * t844 + t873 * t921, t875 * t808 + t873 * t917, -t915 * t967 + t823, t875 * t809 + t873 * t918, t875 * t953, pkin(1) * t753 - t875 * t784 + t873 * t948, pkin(1) * t749 - t875 * t785 + t873 * t949, pkin(1) * t745 + t873 * t899, pkin(1) * t704 + t732 * t1015, t875 * t760 + t873 * t910, t875 * t709 + t873 * t928, t875 * t736 + t873 * t924, t875 * t758 + t873 * t911, t875 * t737 + t873 * t923, t1001 * t769 + t875 * t768 + t823, pkin(1) * t636 + t875 * t642 + t873 * t901, pkin(1) * t644 + t875 * t650 + t873 * t900, pkin(1) * t624 + t875 * t562 + t873 * t895, pkin(1) * t519 + t875 * t576 + t873 * t892, t875 * t600 + t873 * t912, t875 * t546 + t873 * t936, t875 * t607 + t873 * t931, t875 * t599 + t873 * t913, t875 * t608 + t873 * t930, t875 * t651 + t873 * t929, pkin(1) * t515 + t875 * t477 + t873 * t903, pkin(1) * t522 + t875 * t486 + t873 * t902, pkin(1) * t480 + t875 * t430 + t873 * t904, pkin(1) * t423 + t875 * t419 + t873 * t905, t875 * t524 + t873 * t938, t875 * t481 + t873 * t944, t875 * t495 + t873 * t941, t875 * t523 + t873 * t939, t875 * t496 + t873 * t940, t875 * t544 + t873 * t937, pkin(1) * t445 + t875 * t402 + t873 * t907, pkin(1) * t447 + t875 * t404 + t873 * t906, pkin(1) * t435 + t875 * t398 + t873 * t908, pkin(1) * t396 + t875 * t384 + t873 * t909;];
tauB_reg = t1;
