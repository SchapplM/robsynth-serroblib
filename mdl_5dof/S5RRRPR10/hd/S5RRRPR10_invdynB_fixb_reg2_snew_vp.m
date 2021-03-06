% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPR10_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:30:16
% EndTime: 2019-12-31 21:30:46
% DurationCPUTime: 28.63s
% Computational Cost: add. (202963->814), mult. (452634->1322), div. (0->0), fcn. (361324->12), ass. (0->599)
t902 = cos(qJ(2));
t893 = sin(pkin(5));
t898 = sin(qJ(2));
t996 = qJD(1) * t898;
t980 = t893 * t996;
t991 = t893 * qJDD(1);
t858 = -qJD(2) * t980 + t902 * t991;
t895 = cos(pkin(5));
t884 = qJD(1) * t895 + qJD(2);
t865 = t884 * t980;
t829 = t858 - t865;
t1054 = t895 * t829;
t897 = sin(qJ(3));
t901 = cos(qJ(3));
t847 = t901 * t884 - t897 * t980;
t848 = t897 * t884 + t901 * t980;
t892 = sin(pkin(10));
t894 = cos(pkin(10));
t813 = t892 * t847 + t894 * t848;
t995 = qJD(1) * t902;
t979 = t893 * t995;
t877 = -qJD(3) + t979;
t896 = sin(qJ(5));
t900 = cos(qJ(5));
t785 = t813 * t896 + t900 * t877;
t787 = t900 * t813 - t896 * t877;
t727 = t787 * t785;
t992 = qJDD(1) * t898;
t857 = (qJD(2) * t995 + t992) * t893;
t972 = qJDD(1) * t895 + qJDD(2);
t973 = t857 * t897 - t901 * t972;
t798 = -t848 * qJD(3) - t973;
t799 = t847 * qJD(3) + t901 * t857 + t897 * t972;
t975 = -t894 * t798 + t799 * t892;
t966 = qJDD(5) + t975;
t1047 = -t727 + t966;
t1053 = t1047 * t896;
t1052 = t1047 * t900;
t1027 = t813 * t877;
t709 = t975 - t1027;
t811 = -t894 * t847 + t848 * t892;
t1028 = t813 * t811;
t852 = -qJDD(3) + t858;
t915 = -t852 - t1028;
t1051 = t892 * t915;
t1050 = t894 * t915;
t1026 = t847 * t848;
t914 = -t852 + t1026;
t1049 = t897 * t914;
t1048 = t901 * t914;
t743 = t798 * t892 + t799 * t894;
t684 = -t785 * qJD(5) + t900 * t743 - t896 * t852;
t808 = qJD(5) + t811;
t741 = t808 * t785;
t644 = -t741 + t684;
t794 = t811 * t877;
t712 = t743 + t794;
t830 = t847 * t877;
t772 = -t830 - t799;
t770 = -t830 + t799;
t866 = t884 * t979;
t826 = t866 + t857;
t976 = t743 * t896 + t900 * t852;
t641 = (qJD(5) - t808) * t787 + t976;
t767 = (qJD(3) + t877) * t848 + t973;
t783 = t785 ^ 2;
t784 = t787 ^ 2;
t807 = t808 ^ 2;
t809 = t811 ^ 2;
t810 = t813 ^ 2;
t1046 = t847 ^ 2;
t846 = t848 ^ 2;
t1045 = t877 ^ 2;
t1044 = t884 ^ 2;
t1043 = 2 * qJD(4);
t1042 = pkin(2) * t898;
t1041 = pkin(2) * t902;
t1040 = pkin(4) * t892;
t1039 = pkin(7) * t893;
t1038 = t857 * pkin(8);
t1037 = t895 * g(3);
t748 = pkin(4) * t811 - pkin(9) * t813;
t899 = sin(qJ(1));
t903 = cos(qJ(1));
t880 = t903 * g(1) + t899 * g(2);
t904 = qJD(1) ^ 2;
t853 = -t904 * pkin(1) + pkin(7) * t991 - t880;
t971 = -pkin(8) * t898 - t1041;
t997 = qJD(1) * t893;
t856 = t971 * t997;
t1019 = t893 * t898;
t879 = t899 * g(1) - t903 * g(2);
t912 = qJDD(1) * pkin(1) + t904 * t1039 + t879;
t911 = t895 * t912;
t910 = -g(3) * t1019 + t898 * t911;
t909 = pkin(8) * t972 + t910;
t907 = -t1044 * pkin(2) + (t856 * t997 + t853) * t902 + t909;
t970 = -pkin(8) * t902 + t1042;
t998 = qJD(1) * t884;
t908 = -t858 * pkin(2) - t1038 - t1037 + (t970 * t998 - t912) * t893;
t693 = t897 * t908 + t901 * t907;
t822 = -pkin(3) * t877 - qJ(4) * t848;
t650 = -pkin(3) * t1046 + qJ(4) * t798 + t822 * t877 + t693;
t692 = t897 * t907 - t901 * t908;
t906 = pkin(3) * t914 + qJ(4) * t772 - t692;
t978 = t650 * t892 - t894 * t906;
t541 = t852 * pkin(4) - t1045 * pkin(9) + (t1043 + t748) * t813 + t978;
t1036 = t541 * t896;
t974 = t898 * t853 - t902 * t911;
t759 = -t972 * pkin(2) - t1044 * pkin(8) + (g(3) * t902 + t856 * t996) * t893 + t974;
t676 = -t798 * pkin(3) - t1046 * qJ(4) + t848 * t822 + qJDD(4) + t759;
t1035 = t676 * t892;
t1034 = t676 * t894;
t737 = t852 - t1028;
t1033 = t737 * t892;
t1032 = t737 * t894;
t1031 = t759 * t897;
t1030 = t808 * t896;
t1029 = t808 * t900;
t1025 = t877 * t892;
t1024 = t877 * t894;
t1023 = t877 * t897;
t1022 = t877 * t901;
t889 = t893 ^ 2;
t1021 = t889 * t898;
t1020 = t889 * t904;
t1018 = t893 * t902;
t1017 = t895 * t898;
t667 = t727 + t966;
t1016 = t896 * t667;
t564 = t813 * t1043 + t978;
t565 = -0.2e1 * qJD(4) * t811 + t894 * t650 + t892 * t906;
t509 = -t564 * t894 + t565 * t892;
t1015 = t897 * t509;
t789 = t852 + t1026;
t1014 = t897 * t789;
t832 = t893 * t912 + t1037;
t1013 = t898 * t832;
t1012 = t898 * t852;
t1000 = t902 * t904;
t876 = t1000 * t1021;
t854 = -t876 + t972;
t1011 = t898 * t854;
t855 = t876 + t972;
t1010 = t898 * t855;
t1009 = t900 * t541;
t1008 = t900 * t667;
t1007 = t901 * t509;
t1006 = t901 * t759;
t1005 = t901 * t789;
t1004 = t902 * t832;
t1003 = t902 * t853;
t1002 = t902 * t854;
t1001 = t902 * t855;
t542 = -pkin(4) * t1045 - pkin(9) * t852 - t748 * t811 + t565;
t582 = t709 * pkin(4) - pkin(9) * t712 + t676;
t508 = t900 * t542 + t896 * t582;
t890 = t898 ^ 2;
t891 = t902 ^ 2;
t999 = t890 + t891;
t990 = t892 * t727;
t989 = t894 * t727;
t988 = t902 * t1028;
t987 = t902 * t1026;
t986 = t898 * t1028;
t985 = t898 * t1026;
t984 = t890 * t1020;
t983 = t891 * t1020;
t982 = -pkin(4) * t894 - pkin(3);
t981 = t884 * t997;
t507 = t542 * t896 - t900 * t582;
t510 = t564 * t892 + t894 * t565;
t614 = t692 * t897 + t901 * t693;
t834 = -t899 * t879 - t903 * t880;
t874 = t903 * qJDD(1) - t899 * t904;
t969 = -pkin(6) * t874 - t899 * g(3);
t842 = -t984 - t1044;
t814 = -t898 * t842 - t1002;
t968 = pkin(7) * t814 - t1013;
t862 = -t983 - t1044;
t819 = t902 * t862 - t1010;
t967 = pkin(7) * t819 + t1004;
t466 = t896 * t507 + t900 * t508;
t454 = t466 * t892 - t541 * t894;
t455 = t466 * t894 + t541 * t892;
t425 = -t897 * t454 + t901 * t455;
t465 = -t900 * t507 + t896 * t508;
t965 = t425 * t898 - t465 * t902;
t468 = t901 * t510 - t1015;
t964 = t468 * t898 - t676 * t902;
t645 = -t741 - t684;
t576 = -t641 * t900 - t896 * t645;
t688 = t783 + t784;
t546 = t576 * t892 + t688 * t894;
t547 = t576 * t894 - t688 * t892;
t490 = -t897 * t546 + t901 * t547;
t574 = -t641 * t896 + t900 * t645;
t963 = t490 * t898 - t574 * t902;
t642 = (-qJD(5) - t808) * t787 - t976;
t575 = t900 * t642 - t896 * t644;
t726 = -t784 + t783;
t551 = t575 * t892 + t726 * t894;
t552 = t575 * t894 - t726 * t892;
t498 = -t897 * t551 + t901 * t552;
t573 = -t896 * t642 - t900 * t644;
t962 = t498 * t898 + t573 * t902;
t700 = -t807 - t783;
t604 = t900 * t700 - t1053;
t553 = t604 * t892 + t642 * t894;
t554 = t604 * t894 - t642 * t892;
t504 = -t897 * t553 + t901 * t554;
t603 = t896 * t700 + t1052;
t961 = t504 * t898 - t603 * t902;
t718 = -t784 - t807;
t609 = -t896 * t718 - t1008;
t555 = t609 * t892 - t644 * t894;
t556 = t609 * t894 + t644 * t892;
t506 = -t897 * t555 + t901 * t556;
t608 = t900 * t718 - t1016;
t960 = t506 * t898 - t608 * t902;
t736 = -t784 + t807;
t621 = -t896 * t736 + t1052;
t566 = t621 * t892 + t645 * t894;
t568 = t621 * t894 - t645 * t892;
t515 = -t897 * t566 + t901 * t568;
t619 = -t900 * t736 - t1053;
t959 = t515 * t898 + t619 * t902;
t735 = t783 - t807;
t622 = t900 * t735 - t1016;
t567 = t622 * t892 + t641 * t894;
t569 = t622 * t894 - t641 * t892;
t516 = -t897 * t567 + t901 * t569;
t620 = -t896 * t735 - t1008;
t958 = t516 * t898 + t620 * t902;
t683 = -qJD(5) * t787 - t976;
t638 = t785 * t1029 - t896 * t683;
t598 = t638 * t892 + t989;
t600 = t638 * t894 - t990;
t536 = -t897 * t598 + t901 * t600;
t637 = -t785 * t1030 - t900 * t683;
t957 = t536 * t898 + t637 * t902;
t640 = -t787 * t1030 + t900 * t684;
t599 = t640 * t892 - t989;
t601 = t640 * t894 + t990;
t537 = -t897 * t599 + t901 * t601;
t639 = -t787 * t1029 - t896 * t684;
t956 = t537 * t898 + t639 * t902;
t673 = (-t785 * t900 + t787 * t896) * t808;
t627 = t673 * t892 - t894 * t966;
t628 = t673 * t894 + t892 * t966;
t558 = -t897 * t627 + t901 * t628;
t672 = (t785 * t896 + t787 * t900) * t808;
t955 = t558 * t898 + t672 * t902;
t629 = -t709 * t892 + t712 * t894;
t631 = -t709 * t894 - t712 * t892;
t562 = -t897 * t629 + t901 * t631;
t749 = -t810 + t809;
t954 = t562 * t898 + t749 * t902;
t710 = t975 + t1027;
t714 = -t743 + t794;
t630 = -t710 * t892 + t714 * t894;
t632 = -t710 * t894 - t714 * t892;
t563 = -t897 * t630 + t901 * t632;
t719 = -t809 - t810;
t953 = t563 * t898 - t719 * t902;
t745 = -t1045 - t809;
t674 = t745 * t892 + t1050;
t675 = t745 * t894 - t1051;
t595 = -t897 * t674 + t901 * t675;
t952 = t595 * t898 - t709 * t902;
t777 = -t810 - t1045;
t690 = t777 * t894 + t1033;
t691 = -t777 * t892 + t1032;
t612 = -t897 * t690 + t901 * t691;
t951 = t612 * t898 - t712 * t902;
t950 = t614 * t898 - t759 * t902;
t793 = -t810 + t1045;
t696 = t793 * t894 + t1051;
t698 = -t793 * t892 + t1050;
t617 = -t897 * t696 + t901 * t698;
t949 = t617 * t898 + t714 * t902;
t792 = t809 - t1045;
t697 = t792 * t892 - t1032;
t699 = t792 * t894 + t1033;
t618 = -t897 * t697 + t901 * t699;
t948 = t618 * t898 + t710 * t902;
t613 = -t901 * t692 + t897 * t693;
t768 = (-qJD(3) + t877) * t848 - t973;
t707 = t901 * t768 - t897 * t770;
t815 = -t846 + t1046;
t947 = t707 * t898 + t815 * t902;
t708 = -t767 * t901 - t897 * t772;
t788 = t846 + t1046;
t946 = t708 * t898 + t788 * t902;
t802 = -t1045 - t1046;
t733 = t901 * t802 - t1049;
t945 = t733 * t898 + t768 * t902;
t816 = -t846 - t1045;
t747 = -t897 * t816 + t1005;
t944 = t747 * t898 - t770 * t902;
t824 = -t846 + t1045;
t752 = -t897 * t824 + t1048;
t943 = t752 * t898 + t772 * t902;
t823 = -t1045 + t1046;
t753 = t901 * t823 + t1014;
t942 = t753 * t898 + t767 * t902;
t803 = g(3) * t1018 + t974;
t804 = t910 + t1003;
t941 = -t902 * t803 + t898 * t804;
t744 = t898 * t803 + t902 * t804;
t940 = t826 * t902 + t829 * t898;
t827 = -t866 + t857;
t828 = t858 + t865;
t939 = -t827 * t902 + t828 * t898;
t938 = t842 * t902 - t1011;
t861 = t983 - t1044;
t937 = t861 * t898 + t1002;
t860 = -t984 + t1044;
t936 = t860 * t902 + t1010;
t935 = t862 * t898 + t1001;
t833 = t903 * t879 - t899 * t880;
t934 = -t895 * t904 + t998;
t933 = t893 * t972;
t701 = -t811 * t1025 - t894 * t975;
t702 = -t811 * t1024 + t892 * t975;
t625 = -t897 * t701 + t901 * t702;
t931 = t625 * t898 + t988;
t703 = -t813 * t1024 + t743 * t892;
t704 = t813 * t1025 + t743 * t894;
t626 = -t897 * t703 + t901 * t704;
t930 = t626 * t898 - t988;
t762 = t847 * t1022 - t897 * t798;
t929 = t762 * t898 - t987;
t764 = t848 * t1023 + t901 * t799;
t928 = t764 * t898 + t987;
t412 = qJ(4) * t455 + (-pkin(9) * t892 + t982) * t465;
t420 = -qJ(4) * t454 + (-pkin(9) * t894 + t1040) * t465;
t424 = t901 * t454 + t897 * t455;
t398 = -pkin(8) * t424 - t897 * t412 + t901 * t420;
t410 = -pkin(2) * t424 - pkin(3) * t454 + pkin(4) * t541 - pkin(9) * t466;
t415 = t902 * t425 + t898 * t465;
t927 = pkin(7) * t415 + t398 * t898 + t410 * t902;
t457 = -pkin(9) * t574 - t465;
t435 = qJ(4) * t547 + t457 * t892 + t574 * t982;
t441 = -qJ(4) * t546 + t574 * t1040 + t457 * t894;
t489 = t901 * t546 + t897 * t547;
t413 = -pkin(8) * t489 - t897 * t435 + t901 * t441;
t429 = -pkin(2) * t489 - pkin(3) * t546 - pkin(4) * t688 - pkin(9) * t576 - t466;
t473 = t902 * t490 + t898 * t574;
t926 = pkin(7) * t473 + t413 * t898 + t429 * t902;
t481 = -pkin(4) * t603 + t507;
t521 = -pkin(9) * t603 + t1036;
t444 = -pkin(3) * t603 + qJ(4) * t554 + t481 * t894 + t521 * t892;
t450 = -qJ(4) * t553 - t481 * t892 + t521 * t894;
t503 = t901 * t553 + t897 * t554;
t418 = -pkin(8) * t503 - t897 * t444 + t901 * t450;
t452 = -pkin(2) * t503 - pkin(3) * t553 - pkin(4) * t642 - pkin(9) * t604 + t1009;
t478 = t902 * t504 + t898 * t603;
t925 = pkin(7) * t478 + t418 * t898 + t452 * t902;
t482 = -pkin(4) * t608 + t508;
t522 = -pkin(9) * t608 + t1009;
t446 = -pkin(3) * t608 + qJ(4) * t556 + t482 * t894 + t522 * t892;
t451 = -qJ(4) * t555 - t482 * t892 + t522 * t894;
t505 = t901 * t555 + t897 * t556;
t419 = -pkin(8) * t505 - t897 * t446 + t901 * t451;
t453 = -pkin(2) * t505 - pkin(3) * t555 + pkin(4) * t644 - pkin(9) * t609 - t1036;
t480 = t902 * t506 + t898 * t608;
t924 = pkin(7) * t480 + t419 * t898 + t453 * t902;
t467 = t897 * t510 + t1007;
t492 = -pkin(3) * t676 + qJ(4) * t510;
t433 = -pkin(8) * t467 - qJ(4) * t1007 - t897 * t492;
t447 = -pkin(2) * t467 - pkin(3) * t509;
t464 = t902 * t468 + t898 * t676;
t923 = pkin(7) * t464 + t433 * t898 + t447 * t902;
t479 = -pkin(3) * t719 + qJ(4) * t632 + t510;
t485 = -qJ(4) * t630 - t509;
t561 = t901 * t630 + t897 * t632;
t445 = -pkin(8) * t561 - t897 * t479 + t901 * t485;
t530 = -pkin(2) * t561 - pkin(3) * t630;
t545 = t902 * t563 + t898 * t719;
t922 = pkin(7) * t545 + t445 * t898 + t530 * t902;
t577 = -pkin(3) * t709 + qJ(4) * t675 - t1034;
t594 = t901 * t674 + t897 * t675;
t602 = -qJ(4) * t674 + t1035;
t499 = -pkin(8) * t594 - t897 * t577 + t901 * t602;
t512 = -pkin(2) * t594 - pkin(3) * t674 + t564;
t572 = t902 * t595 + t898 * t709;
t921 = pkin(7) * t572 + t499 * t898 + t512 * t902;
t579 = -pkin(3) * t712 + qJ(4) * t691 + t1035;
t607 = -qJ(4) * t690 + t1034;
t611 = t901 * t690 + t897 * t691;
t511 = -pkin(8) * t611 - t897 * t579 + t901 * t607;
t523 = -pkin(2) * t611 - pkin(3) * t690 + t565;
t583 = t902 * t612 + t712 * t898;
t920 = pkin(7) * t583 + t511 * t898 + t523 * t902;
t732 = t897 * t802 + t1048;
t651 = t897 * (t856 * t979 + t1003 + t909) - t901 * (-pkin(8) * t866 - t1038 - t832) + (-t897 * t1044 + t829 * t901 - t732) * pkin(2);
t679 = -pkin(8) * t732 + t1031;
t682 = t902 * t733 - t898 * t768;
t919 = pkin(7) * t682 + t651 * t902 + t679 * t898;
t746 = t901 * t816 + t1014;
t654 = -pkin(2) * t746 + t693;
t685 = -pkin(8) * t746 + t1006;
t686 = t902 * t747 + t770 * t898;
t918 = pkin(7) * t686 + t654 * t902 + t685 * t898;
t778 = t898 * t827 + t902 * t828;
t917 = pkin(7) * t778 + t744;
t706 = -t767 * t897 + t901 * t772;
t584 = -pkin(8) * t706 - t613;
t669 = t902 * t708 - t898 * t788;
t916 = pkin(7) * t669 - t706 * t1041 + t584 * t898;
t592 = t902 * t614 + t898 * t759;
t913 = pkin(7) * t592 + t613 * t971;
t888 = t893 * t889;
t873 = t899 * qJDD(1) + t903 * t904;
t864 = t999 * t1020;
t863 = (t890 - t891) * t1020;
t859 = -pkin(6) * t873 + t903 * g(3);
t836 = t895 * t902 * t852;
t835 = t852 * t1018;
t831 = t999 * t981;
t825 = (t992 + (qJD(2) + t884) * t995) * t893;
t821 = t902 * t857 - t890 * t981;
t820 = -t898 * t858 - t891 * t981;
t818 = t902 * t861 - t1011;
t817 = -t898 * t860 + t1001;
t801 = (t888 * t1000 + t826 * t895) * t898;
t800 = (-t888 * t898 * t904 + t1054) * t902;
t782 = (-t847 * t901 - t848 * t897) * t877;
t781 = (-t847 * t897 + t848 * t901) * t877;
t779 = -t898 * t826 + t902 * t829;
t776 = t893 * t829 + t895 * t935;
t775 = -t893 * t828 + t895 * t937;
t774 = -t893 * t827 + t895 * t936;
t773 = t893 * t935 - t1054;
t766 = -t893 * t825 + t895 * t938;
t765 = t895 * t825 + t893 * t938;
t763 = -t848 * t1022 + t897 * t799;
t761 = t847 * t1023 + t901 * t798;
t760 = t902 * t782 - t1012;
t758 = -t893 * t863 + t895 * t940;
t757 = t893 * t864 + t895 * t939;
t756 = -t895 * t864 + t893 * t939;
t751 = t897 * t823 - t1005;
t750 = t901 * t824 + t1049;
t731 = (t811 * t894 - t813 * t892) * t877;
t730 = (t811 * t892 + t813 * t894) * t877;
t729 = -t899 * t776 + t903 * t819;
t728 = t903 * t776 + t899 * t819;
t725 = t902 * t764 - t985;
t724 = t902 * t762 + t985;
t723 = -t899 * t766 + t903 * t814;
t722 = t903 * t766 + t899 * t814;
t721 = t893 * t832 + t895 * t941;
t720 = -t895 * t832 + t893 * t941;
t717 = t782 * t1017 - t893 * t781 + t836;
t716 = -t899 * t757 + t903 * t778;
t715 = t903 * t757 + t899 * t778;
t705 = t897 * t768 + t901 * t770;
t695 = t902 * t753 - t898 * t767;
t694 = t902 * t752 - t898 * t772;
t687 = -t1013 + (-t773 * t893 - t776 * t895) * pkin(7);
t681 = -t1004 + (-t765 * t893 - t766 * t895) * pkin(7);
t680 = -pkin(1) * t773 + t893 * t803 + t895 * t967;
t678 = t902 * t707 - t898 * t815;
t677 = -pkin(1) * t765 + t893 * t804 + t895 * t968;
t671 = -t893 * t763 + t895 * t928;
t670 = -t893 * t761 + t895 * t929;
t665 = -t897 * t730 + t901 * t731;
t664 = t901 * t730 + t897 * t731;
t663 = pkin(7) * t744 * t895 - pkin(1) * t720;
t662 = -t899 * t721 + t903 * t744;
t661 = t903 * t721 + t899 * t744;
t660 = t902 * t665 - t1012;
t659 = -pkin(1) * t756 + t895 * t917;
t658 = -pkin(2) * t770 + pkin(8) * t747 + t1031;
t657 = -t893 * t751 + t895 * t942;
t656 = -t893 * t750 + t895 * t943;
t655 = (-t720 * t893 - t721 * t895) * pkin(7);
t653 = pkin(2) * t768 + pkin(8) * t733 - t1006;
t652 = (-t756 * t893 - t757 * t895) * pkin(7) - t941;
t649 = -t893 * t746 + t895 * t944;
t648 = t895 * t746 + t893 * t944;
t636 = -t893 * t732 + t895 * t945;
t635 = t895 * t732 + t893 * t945;
t624 = t901 * t703 + t897 * t704;
t623 = t901 * t701 + t897 * t702;
t616 = t901 * t697 + t897 * t699;
t615 = t901 * t696 + t897 * t698;
t610 = -t893 * t705 + t895 * t947;
t606 = -t893 * t706 + t895 * t946;
t605 = t895 * t706 + t893 * t946;
t597 = t902 * t626 + t986;
t596 = t902 * t625 - t986;
t593 = -pkin(2) * t759 + pkin(8) * t614;
t591 = -t899 * t649 + t903 * t686;
t590 = t903 * t649 + t899 * t686;
t589 = t665 * t1017 - t893 * t664 + t836;
t588 = t902 * t618 - t898 * t710;
t587 = t902 * t617 - t898 * t714;
t586 = -t899 * t636 + t903 * t682;
t585 = t903 * t636 + t899 * t682;
t578 = pkin(2) * t788 + pkin(8) * t708 + t614;
t571 = -t899 * t606 + t903 * t669;
t570 = t903 * t606 + t899 * t669;
t560 = t901 * t629 + t897 * t631;
t557 = t901 * t627 + t897 * t628;
t550 = t902 * t562 - t898 * t749;
t549 = -t893 * t624 + t895 * t930;
t548 = -t893 * t623 + t895 * t931;
t544 = -t893 * t613 + t895 * t950;
t543 = t895 * t613 + t893 * t950;
t539 = -t893 * t616 + t895 * t948;
t538 = -t893 * t615 + t895 * t949;
t535 = t901 * t599 + t897 * t601;
t534 = t901 * t598 + t897 * t600;
t533 = t902 * t558 - t898 * t672;
t532 = -t893 * t611 + t895 * t951;
t531 = t895 * t611 + t893 * t951;
t529 = -t893 * t594 + t895 * t952;
t528 = t895 * t594 + t893 * t952;
t527 = -t898 * t654 + t902 * t685 + (-t648 * t893 - t649 * t895) * pkin(7);
t526 = -t898 * t651 + t902 * t679 + (-t635 * t893 - t636 * t895) * pkin(7);
t525 = t902 * t537 - t898 * t639;
t524 = t902 * t536 - t898 * t637;
t520 = -pkin(1) * t648 - t893 * t658 + t895 * t918;
t519 = -t899 * t544 + t903 * t592;
t518 = t903 * t544 + t899 * t592;
t517 = -pkin(1) * t635 - t893 * t653 + t895 * t919;
t514 = t901 * t567 + t897 * t569;
t513 = t901 * t566 + t897 * t568;
t502 = -t893 * t560 + t895 * t954;
t501 = -t899 * t532 + t903 * t583;
t500 = t903 * t532 + t899 * t583;
t497 = t901 * t551 + t897 * t552;
t496 = -t893 * t561 + t895 * t953;
t495 = t895 * t561 + t893 * t953;
t494 = -pkin(2) * t712 + pkin(8) * t612 + t901 * t579 + t897 * t607;
t493 = t706 * t1042 + t902 * t584 + (-t605 * t893 - t606 * t895) * pkin(7);
t491 = -t893 * t557 + t895 * t955;
t488 = -pkin(2) * t709 + pkin(8) * t595 + t901 * t577 + t897 * t602;
t487 = t902 * t516 - t898 * t620;
t486 = t902 * t515 - t898 * t619;
t484 = -t899 * t529 + t903 * t572;
t483 = t903 * t529 + t899 * t572;
t477 = -pkin(1) * t605 - t893 * t578 + t895 * t916;
t476 = -t893 * t535 + t895 * t956;
t475 = -t893 * t534 + t895 * t957;
t474 = t902 * t498 - t898 * t573;
t472 = -t899 * t496 + t903 * t545;
t471 = t903 * t496 + t899 * t545;
t470 = t970 * t613 + (-t543 * t893 - t544 * t895) * pkin(7);
t469 = -pkin(1) * t543 - t893 * t593 + t895 * t913;
t463 = -t893 * t514 + t895 * t958;
t462 = -t893 * t513 + t895 * t959;
t461 = -t893 * t505 + t895 * t960;
t460 = t895 * t505 + t893 * t960;
t459 = -t893 * t503 + t895 * t961;
t458 = t895 * t503 + t893 * t961;
t456 = -t893 * t497 + t895 * t962;
t449 = -t893 * t489 + t895 * t963;
t448 = t895 * t489 + t893 * t963;
t443 = -pkin(2) * t719 + pkin(8) * t563 + t901 * t479 + t897 * t485;
t442 = t902 * t511 - t898 * t523 + (-t531 * t893 - t532 * t895) * pkin(7);
t440 = t902 * t499 - t898 * t512 + (-t528 * t893 - t529 * t895) * pkin(7);
t439 = -t899 * t461 + t903 * t480;
t438 = t903 * t461 + t899 * t480;
t437 = -t899 * t459 + t903 * t478;
t436 = t903 * t459 + t899 * t478;
t434 = -pkin(1) * t531 - t893 * t494 + t895 * t920;
t432 = -t893 * t467 + t895 * t964;
t431 = t895 * t467 + t893 * t964;
t430 = -pkin(2) * t676 + pkin(8) * t468 - qJ(4) * t1015 + t901 * t492;
t428 = -t899 * t449 + t903 * t473;
t427 = t903 * t449 + t899 * t473;
t426 = -pkin(1) * t528 - t893 * t488 + t895 * t921;
t423 = t902 * t445 - t898 * t530 + (-t495 * t893 - t496 * t895) * pkin(7);
t422 = -t899 * t432 + t903 * t464;
t421 = t903 * t432 + t899 * t464;
t417 = -pkin(2) * t608 + pkin(8) * t506 + t901 * t446 + t897 * t451;
t416 = -pkin(2) * t603 + pkin(8) * t504 + t901 * t444 + t897 * t450;
t414 = -pkin(1) * t495 - t893 * t443 + t895 * t922;
t411 = -pkin(2) * t574 + pkin(8) * t490 + t901 * t435 + t897 * t441;
t409 = -t893 * t424 + t895 * t965;
t408 = t895 * t424 + t893 * t965;
t407 = t902 * t419 - t898 * t453 + (-t460 * t893 - t461 * t895) * pkin(7);
t406 = t902 * t418 - t898 * t452 + (-t458 * t893 - t459 * t895) * pkin(7);
t405 = t902 * t433 - t898 * t447 + (-t431 * t893 - t432 * t895) * pkin(7);
t404 = -pkin(1) * t431 - t893 * t430 + t895 * t923;
t403 = t902 * t413 - t898 * t429 + (-t448 * t893 - t449 * t895) * pkin(7);
t402 = -t899 * t409 + t903 * t415;
t401 = t903 * t409 + t899 * t415;
t400 = -pkin(1) * t460 - t893 * t417 + t895 * t924;
t399 = -pkin(1) * t458 - t893 * t416 + t895 * t925;
t397 = -pkin(2) * t465 + pkin(8) * t425 + t901 * t412 + t897 * t420;
t396 = -pkin(1) * t448 - t893 * t411 + t895 * t926;
t395 = t902 * t398 - t898 * t410 + (-t408 * t893 - t409 * t895) * pkin(7);
t394 = -pkin(1) * t408 - t893 * t397 + t895 * t927;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t873, -t874, 0, t834, 0, 0, 0, 0, 0, 0, t729, t723, t716, t662, 0, 0, 0, 0, 0, 0, t586, t591, t571, t519, 0, 0, 0, 0, 0, 0, t484, t501, t472, t422, 0, 0, 0, 0, 0, 0, t437, t439, t428, t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t874, -t873, 0, t833, 0, 0, 0, 0, 0, 0, t728, t722, t715, t661, 0, 0, 0, 0, 0, 0, t585, t590, t570, t518, 0, 0, 0, 0, 0, 0, t483, t500, t471, t421, 0, 0, 0, 0, 0, 0, t436, t438, t427, t401; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t773, t765, t756, t720, 0, 0, 0, 0, 0, 0, t635, t648, t605, t543, 0, 0, 0, 0, 0, 0, t528, t531, t495, t431, 0, 0, 0, 0, 0, 0, t458, t460, t448, t408; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t874, 0, -t873, 0, t969, -t859, -t833, -pkin(6) * t833, -t899 * t801 + t903 * t821, -t899 * t758 + t903 * t779, -t899 * t774 + t903 * t817, -t899 * t800 + t903 * t820, -t899 * t775 + t903 * t818, t903 * t831 + t899 * t933, -pkin(6) * t728 - t899 * t680 + t903 * t687, -pkin(6) * t722 - t899 * t677 + t903 * t681, -pkin(6) * t715 + t903 * t652 - t899 * t659, -pkin(6) * t661 + t903 * t655 - t899 * t663, -t899 * t671 + t903 * t725, -t899 * t610 + t903 * t678, -t899 * t656 + t903 * t694, -t899 * t670 + t903 * t724, -t899 * t657 + t903 * t695, -t899 * t717 + t903 * t760, -pkin(6) * t585 - t899 * t517 + t903 * t526, -pkin(6) * t590 - t899 * t520 + t903 * t527, -pkin(6) * t570 - t899 * t477 + t903 * t493, -pkin(6) * t518 - t899 * t469 + t903 * t470, -t899 * t549 + t903 * t597, -t899 * t502 + t903 * t550, -t899 * t538 + t903 * t587, -t899 * t548 + t903 * t596, -t899 * t539 + t903 * t588, -t899 * t589 + t903 * t660, -pkin(6) * t483 - t899 * t426 + t903 * t440, -pkin(6) * t500 - t899 * t434 + t903 * t442, -pkin(6) * t471 - t899 * t414 + t903 * t423, -pkin(6) * t421 - t899 * t404 + t903 * t405, -t899 * t476 + t903 * t525, -t899 * t456 + t903 * t474, -t899 * t462 + t903 * t486, -t899 * t475 + t903 * t524, -t899 * t463 + t903 * t487, -t899 * t491 + t903 * t533, -pkin(6) * t436 - t899 * t399 + t903 * t406, -pkin(6) * t438 - t899 * t400 + t903 * t407, -pkin(6) * t427 - t899 * t396 + t903 * t403, -pkin(6) * t401 - t899 * t394 + t903 * t395; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t873, 0, t874, 0, t859, t969, t834, pkin(6) * t834, t903 * t801 + t899 * t821, t903 * t758 + t899 * t779, t903 * t774 + t899 * t817, t903 * t800 + t899 * t820, t903 * t775 + t899 * t818, t899 * t831 - t903 * t933, pkin(6) * t729 + t903 * t680 + t899 * t687, pkin(6) * t723 + t903 * t677 + t899 * t681, pkin(6) * t716 + t899 * t652 + t903 * t659, pkin(6) * t662 + t899 * t655 + t903 * t663, t903 * t671 + t899 * t725, t903 * t610 + t899 * t678, t903 * t656 + t899 * t694, t903 * t670 + t899 * t724, t903 * t657 + t899 * t695, t903 * t717 + t899 * t760, pkin(6) * t586 + t903 * t517 + t899 * t526, pkin(6) * t591 + t903 * t520 + t899 * t527, pkin(6) * t571 + t903 * t477 + t899 * t493, pkin(6) * t519 + t903 * t469 + t899 * t470, t903 * t549 + t899 * t597, t903 * t502 + t899 * t550, t903 * t538 + t899 * t587, t903 * t548 + t899 * t596, t903 * t539 + t899 * t588, t903 * t589 + t899 * t660, pkin(6) * t484 + t903 * t426 + t899 * t440, pkin(6) * t501 + t903 * t434 + t899 * t442, pkin(6) * t472 + t903 * t414 + t899 * t423, pkin(6) * t422 + t903 * t404 + t899 * t405, t903 * t476 + t899 * t525, t903 * t456 + t899 * t474, t903 * t462 + t899 * t486, t903 * t475 + t899 * t524, t903 * t463 + t899 * t487, t903 * t491 + t899 * t533, pkin(6) * t437 + t903 * t399 + t899 * t406, pkin(6) * t439 + t903 * t400 + t899 * t407, pkin(6) * t428 + t903 * t396 + t899 * t403, pkin(6) * t402 + t903 * t394 + t899 * t395; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t879, t880, 0, 0, (t889 * t902 * t934 + t857 * t893) * t898, t895 * t863 + t893 * t940, t895 * t827 + t893 * t936, (-t1021 * t934 + t858 * t893) * t902, t895 * t828 + t893 * t937, t895 * t972, pkin(1) * t776 - t895 * t803 + t893 * t967, pkin(1) * t766 - t895 * t804 + t893 * t968, pkin(1) * t757 + t893 * t917, pkin(1) * t721 + t744 * t1039, t895 * t763 + t893 * t928, t895 * t705 + t893 * t947, t895 * t750 + t893 * t943, t895 * t761 + t893 * t929, t895 * t751 + t893 * t942, t1019 * t782 + t895 * t781 + t835, pkin(1) * t636 + t895 * t653 + t893 * t919, pkin(1) * t649 + t895 * t658 + t893 * t918, pkin(1) * t606 + t895 * t578 + t893 * t916, pkin(1) * t544 + t895 * t593 + t893 * t913, t895 * t624 + t893 * t930, t895 * t560 + t893 * t954, t895 * t615 + t893 * t949, t895 * t623 + t893 * t931, t895 * t616 + t893 * t948, t1019 * t665 + t895 * t664 + t835, pkin(1) * t529 + t895 * t488 + t893 * t921, pkin(1) * t532 + t895 * t494 + t893 * t920, pkin(1) * t496 + t895 * t443 + t893 * t922, pkin(1) * t432 + t895 * t430 + t893 * t923, t895 * t535 + t893 * t956, t895 * t497 + t893 * t962, t895 * t513 + t893 * t959, t895 * t534 + t893 * t957, t895 * t514 + t893 * t958, t895 * t557 + t893 * t955, pkin(1) * t459 + t895 * t416 + t893 * t925, pkin(1) * t461 + t895 * t417 + t893 * t924, pkin(1) * t449 + t895 * t411 + t893 * t926, pkin(1) * t409 + t895 * t397 + t893 * t927;];
tauB_reg = t1;
