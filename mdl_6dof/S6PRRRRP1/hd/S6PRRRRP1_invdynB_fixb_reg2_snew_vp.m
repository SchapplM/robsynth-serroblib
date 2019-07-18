% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6PRRRRP1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:31:45
% EndTime: 2019-05-05 09:32:15
% DurationCPUTime: 26.76s
% Computational Cost: add. (171238->853), mult. (339481->1322), div. (0->0), fcn. (252168->12), ass. (0->636)
t877 = sin(qJ(5));
t878 = sin(qJ(4));
t879 = sin(qJ(3));
t882 = cos(qJ(4));
t883 = cos(qJ(3));
t834 = (t878 * t883 + t879 * t882) * qJD(2);
t881 = cos(qJ(5));
t977 = qJD(3) + qJD(4);
t800 = t834 * t877 - t881 * t977;
t802 = t881 * t834 + t877 * t977;
t753 = t802 * t800;
t983 = qJD(2) * qJD(3);
t967 = t883 * t983;
t981 = qJDD(2) * t879;
t841 = t967 + t981;
t968 = t879 * t983;
t979 = qJDD(2) * t883;
t904 = t968 - t979;
t957 = t878 * t841 + t882 * t904;
t756 = -t834 * qJD(4) - t957;
t754 = qJDD(5) - t756;
t959 = -t754 + t753;
t1012 = t959 * t877;
t1011 = t959 * t881;
t1042 = t959 * pkin(5);
t986 = qJD(2) * t879;
t832 = -t882 * t883 * qJD(2) + t878 * t986;
t757 = -t832 * qJD(4) + t882 * t841 - t878 * t904;
t822 = t977 * t832;
t1046 = -t822 + t757;
t873 = sin(pkin(6));
t875 = cos(pkin(6));
t872 = sin(pkin(11));
t874 = cos(pkin(11));
t961 = g(1) * t872 - t874 * g(2);
t995 = g(3) - qJDD(1);
t1045 = -t873 * t995 + t875 * t961;
t788 = t834 * t832;
t976 = qJDD(3) + qJDD(4);
t1038 = -t788 + t976;
t1044 = t1038 * t878;
t1043 = t1038 * t882;
t975 = t977 ^ 2;
t1041 = t872 * t995;
t1040 = t874 * t995;
t710 = -t800 * qJD(5) + t881 * t757 + t877 * t976;
t824 = qJD(5) + t832;
t766 = t824 * t800;
t1039 = -t766 + t710;
t886 = qJD(2) ^ 2;
t857 = t879 * t886 * t883;
t849 = qJDD(3) + t857;
t958 = t757 * t877 - t881 * t976;
t675 = (qJD(5) - t824) * t802 + t958;
t848 = g(1) * t874 + g(2) * t872;
t797 = -t874 * t848 - t872 * t961;
t796 = -t872 * t848 + t874 * t961;
t798 = t800 ^ 2;
t799 = t802 ^ 2;
t823 = t824 ^ 2;
t830 = t832 ^ 2;
t831 = t834 ^ 2;
t730 = -t823 - t798;
t637 = t730 * t877 - t1011;
t1037 = pkin(4) * t637;
t694 = t753 + t754;
t1014 = t694 * t877;
t751 = -t799 - t823;
t641 = t751 * t881 - t1014;
t1036 = pkin(4) * t641;
t1035 = pkin(4) * t878;
t1034 = pkin(7) * t873;
t1033 = pkin(7) * t875;
t679 = -t766 - t710;
t607 = -t675 * t881 - t679 * t877;
t713 = -t798 - t799;
t574 = t607 * t878 - t713 * t882;
t575 = t607 * t882 + t713 * t878;
t509 = t574 * t883 + t575 * t879;
t1032 = pkin(8) * t509;
t638 = t730 * t881 + t1012;
t674 = (qJD(5) + t824) * t802 + t958;
t585 = t638 * t878 - t674 * t882;
t586 = t638 * t882 + t674 * t878;
t522 = t585 * t883 + t586 * t879;
t1031 = pkin(8) * t522;
t1013 = t694 * t881;
t642 = -t751 * t877 - t1013;
t588 = -t1039 * t882 + t642 * t878;
t589 = t1039 * t878 + t642 * t882;
t528 = t588 * t883 + t589 * t879;
t1030 = pkin(8) * t528;
t1029 = pkin(9) * t574;
t1028 = pkin(9) * t585;
t1027 = pkin(9) * t588;
t605 = -t675 * t877 + t679 * t881;
t1026 = pkin(10) * t605;
t1025 = pkin(10) * t637;
t1024 = pkin(10) * t641;
t510 = -t574 * t879 + t575 * t883;
t880 = sin(qJ(2));
t884 = cos(qJ(2));
t936 = t510 * t880 - t605 * t884;
t439 = -t509 * t873 + t875 * t936;
t483 = t510 * t884 + t605 * t880;
t401 = t439 * t874 + t483 * t872;
t1023 = qJ(1) * t401;
t523 = -t585 * t879 + t586 * t883;
t934 = t523 * t880 - t637 * t884;
t455 = -t522 * t873 + t875 * t934;
t499 = t523 * t884 + t637 * t880;
t413 = t455 * t874 + t499 * t872;
t1022 = qJ(1) * t413;
t529 = -t588 * t879 + t589 * t883;
t933 = t529 * t880 - t641 * t884;
t462 = -t528 * t873 + t875 * t933;
t502 = t529 * t884 + t641 * t880;
t417 = t462 * t874 + t502 * t872;
t1021 = qJ(1) * t417;
t818 = t873 * t961 + t875 * t995;
t771 = t1045 * t880 - t884 * t848;
t888 = -t886 * pkin(2) + qJDD(2) * pkin(8) + t771;
t729 = -t879 * t818 + t883 * t888;
t852 = qJD(3) * pkin(3) - pkin(9) * t986;
t870 = t883 ^ 2;
t867 = t870 * t886;
t690 = -pkin(3) * t867 - pkin(9) * t904 - qJD(3) * t852 + t729;
t728 = t883 * t818 + t879 * t888;
t887 = -t728 + (-t841 + t967) * pkin(9) + t849 * pkin(3);
t620 = t882 * t690 + t878 * t887;
t784 = pkin(4) * t832 - pkin(10) * t834;
t603 = -pkin(4) * t975 + pkin(10) * t976 - t832 * t784 + t620;
t770 = -t1045 * t884 - t880 * t848;
t761 = -qJDD(2) * pkin(2) - t886 * pkin(8) + t770;
t716 = t904 * pkin(3) - pkin(9) * t867 + t852 * t986 + t761;
t956 = t834 * t977;
t618 = -t1046 * pkin(10) + (-t756 + t956) * pkin(4) + t716;
t988 = -t877 * t603 + t881 * t618;
t970 = t710 * qJ(6) - t988;
t906 = -qJ(6) * t766 - t970;
t985 = qJD(6) * t802;
t503 = t906 - 0.2e1 * t985 - t1042;
t1020 = t503 * t877;
t1019 = t503 * t881;
t619 = t690 * t878 - t882 * t887;
t552 = -t619 * t882 + t620 * t878;
t1018 = t552 * t879;
t1017 = t552 * t883;
t602 = -t976 * pkin(4) - t975 * pkin(10) + t784 * t834 + t619;
t1016 = t602 * t877;
t1015 = t602 * t881;
t1010 = t716 * t878;
t1009 = t716 * t882;
t1008 = t761 * t879;
t1007 = t761 * t883;
t778 = t788 + t976;
t1006 = t778 * t878;
t1005 = t778 * t882;
t1004 = t818 * t880;
t1003 = t818 * t884;
t1002 = t824 * t877;
t1001 = t824 * t881;
t842 = -0.2e1 * t968 + t979;
t803 = t842 * t883;
t1000 = t849 * t879;
t850 = qJDD(3) - t857;
t999 = t850 * t879;
t998 = t850 * t883;
t869 = t879 ^ 2;
t997 = t869 * t886;
t994 = pkin(1) * t439 + t483 * t1034;
t993 = pkin(1) * t455 + t499 * t1034;
t992 = pkin(1) * t462 + t502 * t1034;
t991 = -pkin(2) * t605 + pkin(8) * t510;
t990 = -pkin(2) * t637 + pkin(8) * t523;
t989 = -pkin(2) * t641 + pkin(8) * t529;
t547 = t881 * t603 + t877 * t618;
t987 = t869 + t870;
t982 = qJDD(2) * t873;
t980 = qJDD(2) * t880;
t978 = qJDD(2) * t884;
t974 = t878 * t753;
t973 = t882 * t753;
t972 = t880 * t788;
t971 = t884 * t788;
t969 = -pkin(4) * t882 - pkin(3);
t438 = t509 * t875 + t873 * t936;
t966 = -pkin(1) * t438 + t483 * t1033;
t454 = t522 * t875 + t873 * t934;
t965 = -pkin(1) * t454 + t499 * t1033;
t461 = t528 * t875 + t873 * t933;
t964 = -pkin(1) * t461 + t502 * t1033;
t963 = -pkin(3) * t637 + pkin(9) * t586;
t962 = -pkin(3) * t641 + pkin(9) * t589;
t553 = t619 * t878 + t882 * t620;
t656 = t728 * t879 + t883 * t729;
t955 = t880 * t857;
t954 = t884 * t857;
t952 = t878 * t956;
t951 = t878 * t822;
t950 = t882 * t956;
t949 = t882 * t822;
t655 = t728 * t883 - t729 * t879;
t843 = t987 * qJDD(2);
t846 = t867 + t997;
t794 = t843 * t884 - t846 * t880;
t948 = pkin(7) * t794 + t655 * t880;
t844 = -t880 * t886 + t978;
t947 = -pkin(7) * t844 - t1004;
t915 = t884 * t886 + t980;
t946 = -pkin(7) * t915 + t1003;
t709 = -qJD(5) * t802 - t958;
t760 = pkin(5) * t824 - qJ(6) * t802;
t905 = t709 * qJ(6) - 0.2e1 * qJD(6) * t800 - t760 * t824 + t547;
t493 = -qJ(6) * t675 + (-t713 - t798) * pkin(5) + t905;
t792 = 0.2e1 * t985;
t494 = t792 + (-t679 + t766) * qJ(6) + t1042 + t970;
t426 = -t493 * t877 + t494 * t881 - t1026;
t562 = -pkin(4) * t605 - pkin(5) * t679;
t573 = pkin(9) * t575;
t396 = -pkin(3) * t605 + t426 * t878 + t562 * t882 + t573;
t397 = t426 * t882 - t562 * t878 - t1029;
t372 = -t396 * t879 + t397 * t883 - t1032;
t891 = -pkin(2) * t509 - pkin(3) * t574 + pkin(4) * t713 - pkin(10) * t607;
t394 = -t493 * t881 - t494 * t877 + t891;
t945 = t372 * t880 + t394 * t884;
t476 = -t1037 + t792 - t906 + 0.2e1 * t1042;
t548 = -t709 * pkin(5) - t798 * qJ(6) + t760 * t802 + qJDD(6) + t602;
t514 = -pkin(5) * t674 + qJ(6) * t730 - t548;
t485 = qJ(6) * t1011 - t514 * t877 - t1025;
t407 = t476 * t882 + t485 * t878 + t963;
t410 = -t476 * t878 + t485 * t882 - t1028;
t376 = -t407 * t879 + t410 * t883 - t1031;
t890 = -pkin(2) * t522 - pkin(3) * t585 + pkin(4) * t674 - pkin(10) * t638;
t427 = -qJ(6) * t1012 - t514 * t881 + t890;
t944 = t376 * t880 + t427 * t884;
t543 = -qJ(6) * t751 + t548;
t611 = -pkin(5) * t1039 - qJ(6) * t694;
t490 = t543 * t881 - t611 * t877 - t1024;
t491 = -t1036 + (-t751 - t798) * pkin(5) + t905;
t409 = t490 * t878 + t491 * t882 + t962;
t419 = t490 * t882 - t491 * t878 - t1027;
t379 = -t409 * t879 + t419 * t883 - t1030;
t889 = -pkin(2) * t528 - pkin(3) * t588 + pkin(4) * t1039 - pkin(10) * t642;
t428 = -t543 * t877 - t611 * t881 + t889;
t943 = t379 * t880 + t428 * t884;
t477 = t547 * t877 + t881 * t988;
t470 = -t477 - t1026;
t420 = t470 * t878 + t605 * t969 + t573;
t432 = t1035 * t605 + t470 * t882 - t1029;
t382 = -t420 * t879 + t432 * t883 - t1032;
t478 = t547 * t881 - t877 * t988;
t408 = -t478 + t891;
t942 = t382 * t880 + t408 * t884;
t511 = -pkin(5) * t798 + t905;
t457 = t511 * t881 - t1020;
t429 = t457 * t878 - t548 * t882;
t430 = t457 * t882 + t548 * t878;
t387 = -t429 * t879 + t430 * t883;
t456 = t511 * t877 + t1019;
t941 = t387 * t880 - t456 * t884;
t512 = -t988 - t1037;
t549 = t1016 - t1025;
t441 = t512 * t882 + t549 * t878 + t963;
t458 = -t512 * t878 + t549 * t882 - t1028;
t389 = -t441 * t879 + t458 * t883 - t1031;
t446 = t890 + t1015;
t940 = t389 * t880 + t446 * t884;
t513 = t547 - t1036;
t551 = t1015 - t1024;
t444 = t513 * t882 + t551 * t878 + t962;
t459 = -t513 * t878 + t551 * t882 - t1027;
t390 = -t444 * t879 + t459 * t883 - t1030;
t447 = t889 - t1016;
t939 = t390 * t880 + t447 * t884;
t468 = t478 * t878 - t602 * t882;
t469 = t478 * t882 + t602 * t878;
t404 = -t468 * t879 + t469 * t883;
t938 = t404 * t880 - t477 * t884;
t496 = t553 * t883 - t1018;
t937 = t496 * t880 - t716 * t884;
t608 = -t1039 * t877 - t674 * t881;
t752 = -t799 + t798;
t582 = t608 * t878 + t752 * t882;
t583 = t608 * t882 - t752 * t878;
t516 = -t582 * t879 + t583 * t883;
t606 = -t1039 * t881 + t674 * t877;
t935 = t516 * t880 + t606 * t884;
t765 = -t799 + t823;
t652 = -t765 * t877 - t1011;
t592 = t652 * t878 + t679 * t882;
t594 = t652 * t882 - t679 * t878;
t532 = -t592 * t879 + t594 * t883;
t650 = -t765 * t881 + t1012;
t932 = t532 * t880 + t650 * t884;
t764 = t798 - t823;
t653 = t764 * t881 - t1014;
t593 = t653 * t878 + t675 * t882;
t595 = t653 * t882 - t675 * t878;
t533 = -t593 * t879 + t595 * t883;
t651 = -t764 * t877 - t1013;
t931 = t533 * t880 + t651 * t884;
t667 = t1001 * t800 - t709 * t877;
t630 = t667 * t878 + t973;
t632 = t667 * t882 - t974;
t565 = -t630 * t879 + t632 * t883;
t666 = -t1002 * t800 - t709 * t881;
t930 = t565 * t880 + t666 * t884;
t669 = -t1002 * t802 + t710 * t881;
t631 = t669 * t878 - t973;
t633 = t669 * t882 + t974;
t566 = -t631 * t879 + t633 * t883;
t668 = -t1001 * t802 - t710 * t877;
t929 = t566 * t880 + t668 * t884;
t701 = (-t800 * t881 + t802 * t877) * t824;
t658 = t701 * t878 - t754 * t882;
t659 = t701 * t882 + t754 * t878;
t591 = -t658 * t879 + t659 * t883;
t700 = (t800 * t877 + t802 * t881) * t824;
t928 = t591 * t880 + t700 * t884;
t737 = (0.2e1 * qJD(4) + qJD(3)) * t834 + t957;
t660 = t1046 * t882 - t737 * t878;
t662 = -t1046 * t878 - t737 * t882;
t598 = -t660 * t879 + t662 * t883;
t787 = -t831 + t830;
t927 = t598 * t880 + t787 * t884;
t739 = qJD(3) * t834 - t957;
t742 = -t822 - t757;
t661 = t739 * t878 + t742 * t882;
t663 = t739 * t882 - t742 * t878;
t599 = -t661 * t879 + t663 * t883;
t759 = -t830 - t831;
t926 = t599 * t880 - t759 * t884;
t772 = -t975 - t830;
t714 = t772 * t878 + t1043;
t715 = t772 * t882 - t1044;
t644 = -t714 * t879 + t715 * t883;
t925 = t644 * t880 - t737 * t884;
t924 = t656 * t880 - t761 * t884;
t813 = -t831 - t975;
t743 = t813 * t882 - t1006;
t744 = -t813 * t878 - t1005;
t665 = -t743 * t879 + t744 * t883;
t923 = -t1046 * t884 + t665 * t880;
t820 = -t831 + t975;
t747 = t820 * t882 + t1044;
t749 = -t820 * t878 + t1043;
t672 = -t747 * t879 + t749 * t883;
t922 = t672 * t880 + t742 * t884;
t819 = t830 - t975;
t748 = t819 * t878 + t1005;
t750 = t819 * t882 - t1006;
t673 = -t748 * t879 + t750 * t883;
t921 = t673 * t880 - t739 * t884;
t711 = t770 * t884 - t771 * t880;
t712 = t770 * t880 + t771 * t884;
t840 = 0.2e1 * t967 + t981;
t790 = -t840 * t879 + t803;
t847 = t867 - t997;
t920 = t790 * t880 + t847 * t884;
t885 = qJD(3) ^ 2;
t856 = -t867 - t885;
t810 = t856 * t883 - t1000;
t919 = t810 * t880 + t842 * t884;
t854 = -t885 - t997;
t812 = -t854 * t879 - t998;
t918 = t812 * t880 - t840 * t884;
t827 = t915 * t875;
t917 = t827 * t874 + t844 * t872;
t782 = t827 * t872 - t844 * t874;
t916 = t843 * t880 + t846 * t884;
t838 = t987 * t983;
t914 = -qJDD(3) * t884 + t838 * t880;
t724 = t882 * t756 + t951;
t725 = -t878 * t756 + t949;
t648 = -t724 * t879 + t725 * t883;
t913 = t648 * t880 + t971;
t726 = t878 * t757 + t950;
t727 = t882 * t757 - t952;
t649 = -t726 * t879 + t727 * t883;
t912 = t649 * t880 - t971;
t911 = (-t438 * t873 - t439 * t875) * pkin(7);
t910 = (-t454 * t873 - t455 * t875) * pkin(7);
t909 = (-t461 * t873 - t462 * t875) * pkin(7);
t855 = t867 - t885;
t809 = t855 * t883 - t999;
t908 = t809 * t880 - t883 * t978;
t839 = t883 * t849;
t853 = t885 - t997;
t811 = -t853 * t879 + t839;
t907 = t811 * t880 - t879 * t978;
t815 = -t870 * t983 + t879 * t904;
t903 = t815 * t880 - t954;
t816 = t841 * t883 - t869 * t983;
t902 = t816 * t880 + t954;
t762 = -t951 - t950;
t763 = -t949 + t952;
t699 = -t762 * t879 + t763 * t883;
t901 = t880 * t699 - t884 * t976;
t471 = -pkin(5) * t548 + qJ(6) * t511;
t392 = -pkin(10) * t456 - qJ(6) * t1019 - t471 * t877;
t415 = -pkin(4) * t456 - pkin(5) * t503;
t365 = -pkin(3) * t456 + pkin(9) * t430 + t392 * t878 + t415 * t882;
t369 = -pkin(9) * t429 + t392 * t882 - t415 * t878;
t386 = t429 * t883 + t430 * t879;
t345 = -pkin(8) * t386 - t365 * t879 + t369 * t883;
t362 = -pkin(2) * t386 - pkin(3) * t429 + pkin(4) * t548 - pkin(10) * t457 + qJ(6) * t1020 - t471 * t881;
t377 = t387 * t884 + t456 * t880;
t900 = pkin(7) * t377 + t345 * t880 + t362 * t884;
t384 = pkin(9) * t469 + (-pkin(10) * t878 + t969) * t477;
t393 = -pkin(9) * t468 + (-pkin(10) * t882 + t1035) * t477;
t403 = t468 * t883 + t469 * t879;
t358 = -pkin(8) * t403 - t384 * t879 + t393 * t883;
t380 = -pkin(2) * t403 - pkin(3) * t468 + pkin(4) * t602 - pkin(10) * t478;
t391 = t404 * t884 + t477 * t880;
t899 = pkin(7) * t391 + t358 * t880 + t380 * t884;
t495 = t553 * t879 + t1017;
t544 = -pkin(3) * t716 + pkin(9) * t553;
t436 = -pkin(8) * t495 - pkin(9) * t1017 - t544 * t879;
t467 = -pkin(2) * t495 - pkin(3) * t552;
t484 = t496 * t884 + t716 * t880;
t898 = pkin(7) * t484 + t436 * t880 + t467 * t884;
t526 = -pkin(3) * t759 + pkin(9) * t663 + t553;
t539 = -pkin(9) * t661 - t552;
t597 = t661 * t883 + t663 * t879;
t448 = -pkin(8) * t597 - t526 * t879 + t539 * t883;
t554 = -pkin(2) * t597 - pkin(3) * t661;
t578 = t599 * t884 + t759 * t880;
t897 = pkin(7) * t578 + t448 * t880 + t554 * t884;
t613 = -pkin(3) * t737 + pkin(9) * t715 - t1009;
t643 = t714 * t883 + t715 * t879;
t645 = -pkin(9) * t714 + t1010;
t540 = -pkin(8) * t643 - t613 * t879 + t645 * t883;
t550 = -pkin(2) * t643 - pkin(3) * t714 + t619;
t610 = t644 * t884 + t737 * t880;
t896 = pkin(7) * t610 + t540 * t880 + t550 * t884;
t616 = -pkin(3) * t1046 + pkin(9) * t744 + t1010;
t657 = -pkin(9) * t743 + t1009;
t664 = t743 * t883 + t744 * t879;
t545 = -pkin(8) * t664 - t616 * t879 + t657 * t883;
t555 = -pkin(2) * t664 - pkin(3) * t743 + t620;
t621 = t1046 * t880 + t665 * t884;
t895 = pkin(7) * t621 + t545 * t880 + t555 * t884;
t806 = t856 * t879 + t839;
t696 = -pkin(2) * t806 + t728;
t722 = -pkin(8) * t806 + t1008;
t768 = t810 * t884 - t842 * t880;
t894 = pkin(7) * t768 + t696 * t884 + t722 * t880;
t808 = t854 * t883 - t999;
t697 = -pkin(2) * t808 + t729;
t723 = -pkin(8) * t808 + t1007;
t769 = t812 * t884 + t840 * t880;
t893 = pkin(7) * t769 + t697 * t884 + t723 * t880;
t628 = t656 * t884 + t761 * t880;
t892 = pkin(7) * t628 - (-pkin(2) * t884 - pkin(8) * t880) * t655;
t828 = t844 * t875;
t826 = t844 * t873;
t825 = t915 * t873;
t817 = qJDD(3) * t880 + t838 * t884;
t807 = t853 * t883 + t1000;
t805 = t855 * t879 + t998;
t804 = (t841 + t967) * t879;
t795 = t914 * t875;
t789 = t840 * t883 + t842 * t879;
t786 = t916 * t875;
t785 = t916 * t873;
t783 = -t828 * t872 - t874 * t915;
t781 = t828 * t874 - t872 * t915;
t776 = t816 * t884 - t955;
t775 = t815 * t884 + t955;
t774 = t811 * t884 + t879 * t980;
t773 = t809 * t884 + t880 * t979;
t758 = t790 * t884 - t847 * t880;
t746 = -t1003 + (t825 * t873 + t827 * t875) * pkin(7);
t745 = -t1004 + (-t826 * t873 - t828 * t875) * pkin(7);
t736 = -t786 * t872 + t794 * t874;
t735 = t786 * t874 + t794 * t872;
t734 = -t804 * t873 + t875 * t902;
t733 = -t803 * t873 + t875 * t903;
t732 = -t807 * t873 + t875 * t907;
t731 = -t805 * t873 + t875 * t908;
t720 = -t808 * t873 + t875 * t918;
t719 = -t806 * t873 + t875 * t919;
t718 = t808 * t875 + t873 * t918;
t717 = t806 * t875 + t873 * t919;
t708 = -t789 * t873 + t875 * t920;
t707 = pkin(2) * t842 + pkin(8) * t810 - t1007;
t706 = -pkin(2) * t840 + pkin(8) * t812 + t1008;
t702 = t712 * t875;
t698 = t762 * t883 + t763 * t879;
t692 = -pkin(1) * t826 + t770 * t873 + t875 * t946;
t691 = pkin(1) * t825 + t771 * t873 + t875 * t947;
t686 = t884 * t699 + t880 * t976;
t685 = -t711 * t875 + t818 * t873;
t684 = -t711 * t873 - t818 * t875;
t683 = -t720 * t872 + t769 * t874;
t682 = -t719 * t872 + t768 * t874;
t681 = t720 * t874 + t769 * t872;
t680 = t719 * t874 + t768 * t872;
t671 = t748 * t883 + t750 * t879;
t670 = t747 * t883 + t749 * t879;
t647 = t726 * t883 + t727 * t879;
t646 = t724 * t883 + t725 * t879;
t639 = pkin(2) * t846 + pkin(8) * t843 + t656;
t635 = t649 * t884 + t972;
t634 = t648 * t884 - t972;
t629 = -pkin(2) * t761 + pkin(8) * t656;
t627 = -pkin(1) * t684 + t1033 * t712;
t626 = t673 * t884 + t739 * t880;
t625 = t672 * t884 - t742 * t880;
t624 = -t685 * t872 + t712 * t874;
t623 = t685 * t874 + t712 * t872;
t622 = -t873 * t698 + t875 * t901;
t612 = t655 * t884 + (-t785 * t873 - t786 * t875) * pkin(7);
t609 = (-t684 * t873 - t685 * t875) * pkin(7);
t596 = t660 * t883 + t662 * t879;
t590 = t658 * t883 + t659 * t879;
t581 = t598 * t884 - t787 * t880;
t580 = -t697 * t880 + t723 * t884 + (-t718 * t873 - t720 * t875) * pkin(7);
t579 = -t696 * t880 + t722 * t884 + (-t717 * t873 - t719 * t875) * pkin(7);
t577 = -t671 * t873 + t875 * t921;
t576 = -t670 * t873 + t875 * t922;
t572 = -t647 * t873 + t875 * t912;
t571 = -t646 * t873 + t875 * t913;
t570 = -t664 * t873 + t875 * t923;
t569 = t664 * t875 + t873 * t923;
t568 = t655 * t873 + t875 * t924;
t567 = -t655 * t875 + t873 * t924;
t564 = t631 * t883 + t633 * t879;
t563 = t630 * t883 + t632 * t879;
t561 = -pkin(1) * t718 - t706 * t873 + t875 * t893;
t560 = -pkin(1) * t717 - t707 * t873 + t875 * t894;
t559 = t591 * t884 - t700 * t880;
t558 = -pkin(1) * t785 - t639 * t873 + t875 * t948;
t557 = -t643 * t873 + t875 * t925;
t556 = t643 * t875 + t873 * t925;
t542 = t566 * t884 - t668 * t880;
t541 = t565 * t884 - t666 * t880;
t538 = -t568 * t872 + t628 * t874;
t537 = t568 * t874 + t628 * t872;
t536 = -pkin(2) * t1046 + pkin(8) * t665 + t616 * t883 + t657 * t879;
t535 = -t570 * t872 + t621 * t874;
t534 = t570 * t874 + t621 * t872;
t531 = t593 * t883 + t595 * t879;
t530 = t592 * t883 + t594 * t879;
t525 = -t596 * t873 + t875 * t927;
t524 = -pkin(2) * t737 + pkin(8) * t644 + t613 * t883 + t645 * t879;
t520 = -t597 * t873 + t875 * t926;
t519 = t597 * t875 + t873 * t926;
t518 = -t557 * t872 + t610 * t874;
t517 = t557 * t874 + t610 * t872;
t515 = t582 * t883 + t583 * t879;
t507 = -t590 * t873 + t875 * t928;
t506 = t590 * t875 + t873 * t928;
t505 = t533 * t884 - t651 * t880;
t504 = t532 * t884 - t650 * t880;
t492 = t516 * t884 - t606 * t880;
t489 = -t564 * t873 + t875 * t929;
t488 = -t563 * t873 + t875 * t930;
t487 = t564 * t875 + t873 * t929;
t486 = t563 * t875 + t873 * t930;
t480 = -t520 * t872 + t578 * t874;
t479 = t520 * t874 + t578 * t872;
t475 = -(pkin(2) * t880 - pkin(8) * t884) * t655 + (-t567 * t873 - t568 * t875) * pkin(7);
t474 = -pkin(1) * t567 - t629 * t873 + t875 * t892;
t473 = -t507 * t872 + t559 * t874;
t472 = t507 * t874 + t559 * t872;
t466 = -t531 * t873 + t875 * t931;
t465 = -t530 * t873 + t875 * t932;
t464 = t531 * t875 + t873 * t931;
t463 = t530 * t875 + t873 * t932;
t452 = -t489 * t872 + t542 * t874;
t451 = -t488 * t872 + t541 * t874;
t450 = t489 * t874 + t542 * t872;
t449 = t488 * t874 + t541 * t872;
t445 = -pkin(2) * t759 + pkin(8) * t599 + t526 * t883 + t539 * t879;
t443 = -t515 * t873 + t875 * t935;
t442 = t515 * t875 + t873 * t935;
t440 = t545 * t884 - t555 * t880 + (-t569 * t873 - t570 * t875) * pkin(7);
t435 = -t495 * t873 + t875 * t937;
t434 = t495 * t875 + t873 * t937;
t433 = t540 * t884 - t550 * t880 + (-t556 * t873 - t557 * t875) * pkin(7);
t431 = -pkin(2) * t716 + pkin(8) * t496 - pkin(9) * t1018 + t544 * t883;
t425 = -t466 * t872 + t505 * t874;
t424 = -t465 * t872 + t504 * t874;
t423 = t466 * t874 + t505 * t872;
t422 = t465 * t874 + t504 * t872;
t421 = -pkin(1) * t569 - t536 * t873 + t875 * t895;
t418 = -t462 * t872 + t502 * t874;
t416 = qJ(1) * t418;
t414 = -t455 * t872 + t499 * t874;
t412 = qJ(1) * t414;
t411 = -pkin(1) * t556 - t524 * t873 + t875 * t896;
t406 = -t443 * t872 + t492 * t874;
t405 = t443 * t874 + t492 * t872;
t402 = -t439 * t872 + t483 * t874;
t400 = qJ(1) * t402;
t399 = -t435 * t872 + t484 * t874;
t398 = t435 * t874 + t484 * t872;
t395 = t448 * t884 - t554 * t880 + (-t519 * t873 - t520 * t875) * pkin(7);
t388 = t444 * t883 + t459 * t879 + t989;
t385 = t441 * t883 + t458 * t879 + t990;
t383 = -pkin(1) * t519 - t445 * t873 + t875 * t897;
t381 = t420 * t883 + t432 * t879 + t991;
t378 = t409 * t883 + t419 * t879 + t989;
t375 = t407 * t883 + t410 * t879 + t990;
t374 = -t403 * t873 + t875 * t938;
t373 = t403 * t875 + t873 * t938;
t371 = t436 * t884 - t467 * t880 + (-t434 * t873 - t435 * t875) * pkin(7);
t370 = t396 * t883 + t397 * t879 + t991;
t368 = -pkin(1) * t434 - t431 * t873 + t875 * t898;
t367 = t390 * t884 - t447 * t880 + t909;
t366 = t389 * t884 - t446 * t880 + t910;
t364 = -t386 * t873 + t875 * t941;
t363 = t386 * t875 + t873 * t941;
t361 = -t374 * t872 + t391 * t874;
t360 = t374 * t874 + t391 * t872;
t359 = t379 * t884 - t428 * t880 + t909;
t357 = t376 * t884 - t427 * t880 + t910;
t356 = t382 * t884 - t408 * t880 + t911;
t355 = -pkin(2) * t477 + pkin(8) * t404 + t384 * t883 + t393 * t879;
t354 = -t388 * t873 + t875 * t939 + t964;
t353 = -t385 * t873 + t875 * t940 + t965;
t352 = t372 * t884 - t394 * t880 + t911;
t351 = -t364 * t872 + t377 * t874;
t350 = t364 * t874 + t377 * t872;
t349 = -t378 * t873 + t875 * t943 + t964;
t348 = -t381 * t873 + t875 * t942 + t966;
t347 = -t375 * t873 + t875 * t944 + t965;
t346 = -t370 * t873 + t875 * t945 + t966;
t344 = -pkin(2) * t456 + pkin(8) * t387 + t365 * t883 + t369 * t879;
t343 = t358 * t884 - t380 * t880 + (-t373 * t873 - t374 * t875) * pkin(7);
t342 = -pkin(1) * t373 - t355 * t873 + t875 * t899;
t341 = t345 * t884 - t362 * t880 + (-t363 * t873 - t364 * t875) * pkin(7);
t340 = -pkin(1) * t363 - t344 * t873 + t875 * t900;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t797, 0, 0, 0, 0, 0, 0, t783, t782, 0, t624, 0, 0, 0, 0, 0, 0, t682, t683, t736, t538, 0, 0, 0, 0, 0, 0, t518, t535, t480, t399, 0, 0, 0, 0, 0, 0, t414, t418, t402, t361, 0, 0, 0, 0, 0, 0, t414, t418, t402, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t796, 0, 0, 0, 0, 0, 0, t781, -t917, 0, t623, 0, 0, 0, 0, 0, 0, t680, t681, t735, t537, 0, 0, 0, 0, 0, 0, t517, t534, t479, t398, 0, 0, 0, 0, 0, 0, t413, t417, t401, t360, 0, 0, 0, 0, 0, 0, t413, t417, t401, t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t995, 0, 0, 0, 0, 0, 0, t826, -t825, 0, t684, 0, 0, 0, 0, 0, 0, t717, t718, t785, t567, 0, 0, 0, 0, 0, 0, t556, t569, t519, t434, 0, 0, 0, 0, 0, 0, t454, t461, t438, t373, 0, 0, 0, 0, 0, 0, t454, t461, t438, t363; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1041, -t1040, -t796, -qJ(1) * t796, 0, 0, -t782, 0, t783, t872 * t982, -qJ(1) * t781 - t692 * t872 + t745 * t874, qJ(1) * t917 - t691 * t872 + t746 * t874, -t702 * t872 + t711 * t874, -qJ(1) * t623 + t609 * t874 - t627 * t872, -t734 * t872 + t776 * t874, -t708 * t872 + t758 * t874, -t732 * t872 + t774 * t874, -t733 * t872 + t775 * t874, -t731 * t872 + t773 * t874, -t795 * t872 + t817 * t874, -qJ(1) * t680 - t560 * t872 + t579 * t874, -qJ(1) * t681 - t561 * t872 + t580 * t874, -qJ(1) * t735 - t558 * t872 + t612 * t874, -qJ(1) * t537 - t474 * t872 + t475 * t874, -t572 * t872 + t635 * t874, -t525 * t872 + t581 * t874, -t576 * t872 + t625 * t874, -t571 * t872 + t634 * t874, -t577 * t872 + t626 * t874, -t622 * t872 + t686 * t874, -qJ(1) * t517 - t411 * t872 + t433 * t874, -qJ(1) * t534 - t421 * t872 + t440 * t874, -qJ(1) * t479 - t383 * t872 + t395 * t874, -qJ(1) * t398 - t368 * t872 + t371 * t874, t452, t406, t424, t451, t425, t473, -t353 * t872 + t366 * t874 - t1022, -t354 * t872 + t367 * t874 - t1021, -t348 * t872 + t356 * t874 - t1023, -qJ(1) * t360 - t342 * t872 + t343 * t874, t452, t406, t424, t451, t425, t473, -t347 * t872 + t357 * t874 - t1022, -t349 * t872 + t359 * t874 - t1021, -t346 * t872 + t352 * t874 - t1023, -qJ(1) * t350 - t340 * t872 + t341 * t874; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1040, -t1041, t797, qJ(1) * t797, 0, 0, t917, 0, t781, -t874 * t982, qJ(1) * t783 + t692 * t874 + t745 * t872, qJ(1) * t782 + t691 * t874 + t746 * t872, t702 * t874 + t711 * t872, qJ(1) * t624 + t609 * t872 + t627 * t874, t734 * t874 + t776 * t872, t708 * t874 + t758 * t872, t732 * t874 + t774 * t872, t733 * t874 + t775 * t872, t731 * t874 + t773 * t872, t795 * t874 + t817 * t872, qJ(1) * t682 + t560 * t874 + t579 * t872, qJ(1) * t683 + t561 * t874 + t580 * t872, qJ(1) * t736 + t558 * t874 + t612 * t872, qJ(1) * t538 + t474 * t874 + t475 * t872, t572 * t874 + t635 * t872, t525 * t874 + t581 * t872, t576 * t874 + t625 * t872, t571 * t874 + t634 * t872, t577 * t874 + t626 * t872, t622 * t874 + t686 * t872, qJ(1) * t518 + t411 * t874 + t433 * t872, qJ(1) * t535 + t421 * t874 + t440 * t872, qJ(1) * t480 + t383 * t874 + t395 * t872, qJ(1) * t399 + t368 * t874 + t371 * t872, t450, t405, t422, t449, t423, t472, t353 * t874 + t366 * t872 + t412, t354 * t874 + t367 * t872 + t416, t348 * t874 + t356 * t872 + t400, qJ(1) * t361 + t342 * t874 + t343 * t872, t450, t405, t422, t449, t423, t472, t347 * t874 + t357 * t872 + t412, t349 * t874 + t359 * t872 + t416, t346 * t874 + t352 * t872 + t400, qJ(1) * t351 + t340 * t874 + t341 * t872; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t961, t848, 0, 0, 0, 0, t825, 0, t826, t875 * qJDD(2), pkin(1) * t828 - t770 * t875 + t873 * t946, -pkin(1) * t827 - t771 * t875 + t873 * t947, t712 * t873, pkin(1) * t685 + t1034 * t712, t804 * t875 + t873 * t902, t789 * t875 + t873 * t920, t807 * t875 + t873 * t907, t803 * t875 + t873 * t903, t805 * t875 + t873 * t908, t914 * t873, pkin(1) * t719 + t707 * t875 + t873 * t894, pkin(1) * t720 + t706 * t875 + t873 * t893, pkin(1) * t786 + t639 * t875 + t873 * t948, pkin(1) * t568 + t629 * t875 + t873 * t892, t647 * t875 + t873 * t912, t596 * t875 + t873 * t927, t670 * t875 + t873 * t922, t646 * t875 + t873 * t913, t671 * t875 + t873 * t921, t875 * t698 + t873 * t901, pkin(1) * t557 + t524 * t875 + t873 * t896, pkin(1) * t570 + t536 * t875 + t873 * t895, pkin(1) * t520 + t445 * t875 + t873 * t897, pkin(1) * t435 + t431 * t875 + t873 * t898, t487, t442, t463, t486, t464, t506, t385 * t875 + t873 * t940 + t993, t388 * t875 + t873 * t939 + t992, t381 * t875 + t873 * t942 + t994, pkin(1) * t374 + t355 * t875 + t873 * t899, t487, t442, t463, t486, t464, t506, t375 * t875 + t873 * t944 + t993, t378 * t875 + t873 * t943 + t992, t370 * t875 + t873 * t945 + t994, pkin(1) * t364 + t344 * t875 + t873 * t900;];
tauB_reg  = t1;