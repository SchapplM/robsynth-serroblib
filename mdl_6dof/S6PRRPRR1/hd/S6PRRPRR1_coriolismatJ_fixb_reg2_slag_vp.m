% Calculate inertial parameters regressor of coriolis matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:14
% EndTime: 2019-03-08 21:54:52
% DurationCPUTime: 33.07s
% Computational Cost: add. (24409->755), mult. (53211->1051), div. (0->0), fcn. (64791->12), ass. (0->603)
t801 = qJD(3) + qJD(5);
t599 = sin(qJ(3));
t601 = cos(qJ(3));
t596 = sin(pkin(6));
t955 = sin(qJ(2));
t796 = t596 * t955;
t941 = cos(pkin(6));
t643 = t599 * t941 + t601 * t796;
t939 = sin(pkin(12));
t622 = t939 * t643;
t644 = t599 * t796 - t601 * t941;
t940 = cos(pkin(12));
t624 = t940 * t644;
t1016 = -t622 - t624;
t473 = t643 * t940 - t939 * t644;
t598 = sin(qJ(5));
t956 = cos(qJ(5));
t1043 = t598 * t1016 + t956 * t473;
t1031 = t956 * t1016;
t1036 = t598 * t473;
t329 = t1031 - t1036;
t600 = cos(qJ(6));
t1058 = t1043 * t600;
t597 = sin(qJ(6));
t602 = cos(qJ(2));
t898 = t596 * t602;
t276 = -t597 * t898 + t1058;
t925 = t276 * t600;
t894 = t597 * t1043;
t275 = t600 * t898 + t894;
t927 = t275 * t597;
t1068 = (t1043 - t925 - t927) * t329;
t1069 = t1068 * qJD(1);
t745 = t939 * t601;
t748 = t940 * t599;
t566 = t748 + t745;
t561 = t956 * t566;
t746 = t939 * t599;
t747 = t940 * t601;
t565 = -t746 + t747;
t883 = t598 * t565;
t1004 = t561 + t883;
t766 = -t898 / 0.2e1;
t1012 = t1004 * t766;
t629 = t566 * t898;
t528 = t956 * t629;
t534 = t565 * t898;
t884 = t598 * t534;
t869 = -t884 / 0.2e1 - t528 / 0.2e1;
t1029 = t1012 + t869;
t1067 = qJD(2) * t1029 - t1043 * t801;
t984 = -t1058 / 0.2e1;
t1065 = t1058 / 0.2e1;
t594 = t600 ^ 2;
t592 = t597 ^ 2;
t962 = t592 / 0.2e1;
t750 = t594 / 0.2e1 + t962;
t1064 = t750 * t329;
t868 = t592 + t594;
t1063 = t868 * t329;
t517 = -t956 * t565 + t598 * t566;
t1019 = t600 * t517;
t1037 = t1019 / 0.2e1;
t1038 = -t1019 / 0.2e1;
t1042 = t1038 + t1037;
t1055 = qJD(6) * t1042;
t1025 = t1004 ^ 2;
t903 = t517 ^ 2;
t1026 = t1025 - t903;
t1052 = t1026 * t600;
t1057 = qJD(2) * t1052;
t1062 = t1055 + t1057;
t1061 = 0.2e1 * t597;
t944 = -qJ(4) - pkin(8);
t723 = t944 * t939;
t572 = t599 * t723;
t724 = t944 * t940;
t530 = t601 * t724 - t572;
t947 = t565 * pkin(9);
t632 = t530 - t947;
t441 = t956 * t632;
t677 = t599 * t724;
t532 = t601 * t723 + t677;
t945 = t566 * pkin(9);
t462 = t532 - t945;
t889 = t598 * t462;
t319 = -t441 + t889;
t1060 = -t319 / 0.2e1;
t790 = t940 * pkin(3);
t727 = t790 + pkin(4);
t789 = t939 * pkin(3);
t558 = t598 * t789 - t956 * t727;
t552 = -pkin(5) + t558;
t971 = -t552 / 0.2e1;
t1059 = t1043 * pkin(5);
t1053 = t1026 * t597;
t1056 = qJD(2) * t1053;
t982 = -t1031 / 0.2e1;
t533 = -t944 * t747 + t572;
t1054 = t530 + t533;
t958 = -t600 / 0.2e1;
t959 = t597 / 0.2e1;
t1018 = (t275 * t958 + t276 * t959) * t517;
t626 = t598 * t629;
t791 = t956 * t534;
t405 = t791 - t626;
t892 = t597 * t405;
t379 = t600 * t796 - t892;
t739 = t597 * t796;
t878 = t600 * t405;
t380 = t739 + t878;
t957 = t600 / 0.2e1;
t960 = -t597 / 0.2e1;
t872 = t379 * t960 + t380 * t957;
t1028 = t1018 - t872;
t1051 = qJD(1) * t1028;
t1030 = t1012 - t869;
t1050 = qJD(1) * t1030;
t1027 = t1018 + t872;
t1049 = qJD(2) * t1027;
t1048 = qJD(2) * t1028;
t1046 = qJD(2) * t1030;
t1045 = t1026 * qJD(2);
t728 = t561 / 0.2e1;
t1002 = t728 + t883 / 0.2e1;
t1008 = t1004 * qJD(2);
t1033 = t517 * t1008;
t1044 = qJD(6) * t1002 + t1033;
t404 = t528 + t884;
t1041 = t872 * pkin(10) - t404 * pkin(5) / 0.2e1;
t1040 = pkin(10) * t1064 - t1059 / 0.2e1;
t559 = t598 * t727 + t789 * t956;
t553 = pkin(10) + t559;
t1039 = t404 * t971 - t872 * t553;
t1011 = t597 * t1004;
t763 = t1011 / 0.2e1;
t980 = t1016 / 0.2e1;
t1020 = t594 * t517;
t495 = -t1020 / 0.2e1;
t496 = t1020 / 0.2e1;
t1021 = t1004 * pkin(5);
t949 = t517 * pkin(10);
t377 = t1021 + t949;
t1024 = -t517 / 0.2e1;
t975 = t517 / 0.2e1;
t1023 = -t1004 / 0.2e1;
t1022 = t1004 / 0.2e1;
t954 = pkin(10) * t1004;
t586 = -pkin(3) * t601 - pkin(2);
t544 = -pkin(4) * t565 + t586;
t725 = pkin(5) * t517 - t954;
t631 = t544 + t725;
t531 = t745 * t944 + t677;
t618 = t531 - t945;
t616 = t598 * t618;
t463 = t533 + t947;
t794 = t956 * t463;
t322 = t794 + t616;
t896 = t597 * t322;
t146 = -t600 * t631 + t896;
t881 = t600 * t322;
t147 = t597 * t631 + t881;
t1017 = (t146 * t600 - t147 * t597) * t517;
t1015 = t801 * t517;
t1014 = qJD(2) * t517;
t1013 = qJD(4) * t517;
t1009 = t1002 * qJD(2);
t1007 = 0.2e1 * t1004;
t581 = t594 - t592;
t1006 = t581 * t801;
t579 = t955 * t602 * t596 ^ 2;
t361 = 0.2e1 * t1038;
t1005 = t552 + t558;
t593 = t599 ^ 2;
t595 = t601 ^ 2;
t1003 = t593 + t595;
t627 = t598 * t632;
t795 = t956 * t462;
t321 = t795 + t627;
t882 = t600 * t321;
t591 = t599 * pkin(3);
t946 = t566 * pkin(4);
t547 = t591 + t946;
t323 = t547 + t377;
t895 = t597 * t323;
t149 = t882 + t895;
t933 = t149 * t600;
t880 = t600 * t323;
t897 = t597 * t321;
t148 = t880 - t897;
t936 = t148 * t597;
t1000 = t933 / 0.2e1 - t936 / 0.2e1;
t966 = t558 / 0.2e1;
t751 = t966 + t971;
t965 = t559 / 0.2e1;
t969 = -t553 / 0.2e1;
t628 = (t969 + t965) * t1004 + t751 * t517;
t999 = -t954 / 0.2e1 + t628;
t890 = t597 * t600;
t336 = (t1024 + t975) * t890;
t344 = (t962 - t594 / 0.2e1) * t1004;
t830 = t344 * qJD(6);
t998 = -qJD(5) * t336 + t830;
t335 = t361 * t597;
t997 = -qJD(5) * t335 + t830;
t996 = qJD(3) * t336 + t830;
t995 = -qJD(3) * t335 + t830;
t590 = qJD(6) * t600;
t583 = t597 * t590;
t831 = t336 * qJD(2);
t994 = t831 - t583;
t993 = t831 + t583;
t786 = qJD(2) * t890;
t130 = t1025 * t786 + t344 * t801;
t991 = t566 ^ 2;
t990 = t148 / 0.2e1;
t989 = -t149 / 0.2e1;
t988 = -t275 / 0.2e1;
t987 = t276 / 0.2e1;
t442 = t956 * t618;
t888 = t598 * t463;
t320 = -t442 + t888;
t986 = t320 / 0.2e1;
t985 = t322 / 0.2e1;
t983 = t1043 / 0.2e1;
t981 = -t473 / 0.2e1;
t979 = t473 / 0.2e1;
t970 = t552 / 0.2e1;
t968 = t553 / 0.2e1;
t967 = -t558 / 0.2e1;
t964 = t565 / 0.2e1;
t963 = -t566 / 0.2e1;
t953 = t319 * pkin(5);
t633 = t1022 * t1043;
t52 = t1023 * t1043 + t633;
t943 = t52 * qJD(3);
t942 = qJD(3) * pkin(3);
t938 = t146 * t597;
t937 = t147 * t600;
t935 = t148 * t600;
t934 = t149 * t597;
t879 = t600 * t377;
t920 = t320 * t597;
t186 = t879 + t920;
t932 = t186 * t597;
t931 = t186 * t600;
t893 = t597 * t377;
t919 = t320 * t600;
t187 = t893 - t919;
t930 = t187 * t597;
t929 = t187 * t600;
t928 = t275 * t517;
t926 = t276 * t517;
t924 = t319 * t320;
t923 = t319 * t597;
t922 = t319 * t600;
t921 = t320 * t1004;
t916 = t329 * t404;
t913 = t379 * t600;
t911 = t380 * t597;
t720 = t750 * t517;
t634 = t1004 * t970 - t553 * t720;
t672 = -t935 / 0.2e1 - t934 / 0.2e1;
t40 = t634 + t672;
t909 = t40 * qJD(2);
t908 = t404 * t320;
t906 = t404 * t597;
t905 = t404 * t600;
t902 = t52 * qJD(2);
t53 = -t275 * t379 + t276 * t380 - t916;
t901 = t53 * qJD(1);
t640 = -pkin(10) * t720 - t1021 / 0.2e1;
t671 = t931 / 0.2e1 + t930 / 0.2e1;
t54 = t640 - t671;
t900 = t54 * qJD(2);
t349 = t597 * t517;
t876 = t600 * t1004;
t96 = t1043 * t405 - t579 - t916;
t873 = t96 * qJD(1);
t752 = 0.2e1 * t1022;
t356 = t752 * t600;
t849 = qJD(3) * t600;
t871 = t356 * qJD(5) + t1004 * t849;
t870 = t626 / 0.2e1 - t791 / 0.2e1;
t508 = t592 * t517;
t485 = -t508 / 0.2e1;
t768 = t517 * t962;
t714 = t496 + t768;
t190 = t495 + t485 + t714;
t867 = qJD(2) * t190;
t207 = t903 + t1025;
t195 = t207 * t597;
t865 = qJD(2) * t195;
t197 = t207 * t600;
t863 = qJD(2) * t197;
t860 = qJD(2) * t344;
t370 = t581 * t1025;
t859 = qJD(2) * t370;
t857 = qJD(2) * t544;
t856 = qJD(2) * t596;
t855 = qJD(2) * t601;
t852 = qJD(3) * t1004;
t850 = qJD(3) * t597;
t848 = qJD(4) * t1004;
t844 = qJD(5) * t1004;
t843 = qJD(5) * t544;
t842 = qJD(5) * t597;
t841 = qJD(5) * t600;
t840 = qJD(6) * t597;
t754 = t979 + t981;
t117 = (t980 + t622 / 0.2e1 + t624 / 0.2e1) * t565 - t754 * t566;
t839 = t117 * qJD(2);
t486 = t508 / 0.2e1;
t193 = t496 + t486 + t714;
t837 = t193 * qJD(2);
t198 = -t1016 * t629 + t473 * t534 - t579;
t836 = t198 * qJD(1);
t835 = t207 * qJD(2);
t661 = t1004 * t967 + t559 * t975;
t587 = t591 / 0.2e1;
t740 = t587 + t946 / 0.2e1;
t215 = t661 + t740;
t833 = t215 * qJD(2);
t829 = t1011 * qJD(2);
t753 = t1022 + t1023;
t347 = t753 * t597;
t828 = t347 * qJD(2);
t348 = t752 * t597;
t827 = t348 * qJD(2);
t826 = t349 * qJD(2);
t351 = t1024 * t1061;
t341 = t351 * qJD(2);
t355 = t753 * t600;
t825 = t355 * qJD(2);
t824 = t356 * qJD(2);
t823 = t1019 * qJD(2);
t360 = 0.2e1 * t1037;
t822 = t360 * qJD(2);
t821 = t361 * qJD(2);
t369 = t508 + t1020;
t820 = t369 * qJD(2);
t392 = (t1003 - 0.1e1) * t579;
t819 = t392 * qJD(1);
t424 = 0.2e1 * t728 + t883;
t817 = t424 * qJD(2);
t630 = (t939 * t964 + t940 * t963) * pkin(3);
t456 = -t591 / 0.2e1 + t630;
t816 = t456 * qJD(2);
t564 = t565 ^ 2;
t460 = t564 - t991;
t815 = t460 * qJD(2);
t514 = t728 - t561 / 0.2e1;
t812 = t514 * qJD(2);
t811 = t514 * qJD(5);
t527 = t564 + t991;
t808 = t527 * qJD(2);
t549 = t559 * qJD(5);
t807 = t565 * qJD(2);
t806 = t566 * qJD(2);
t805 = t566 * qJD(3);
t582 = t595 - t593;
t804 = t582 * qJD(2);
t803 = t599 * qJD(3);
t802 = t601 * qJD(3);
t799 = pkin(2) * t599 * qJD(2);
t798 = pkin(2) * t855;
t788 = t592 * t1008;
t787 = t594 * t1008;
t785 = t597 * t849;
t784 = t597 * t841;
t783 = qJD(6) * t517 * t1004;
t780 = t1004 * t1014;
t779 = t565 * t806;
t778 = t565 * t805;
t777 = t602 * t856;
t776 = t599 * t802;
t775 = t600 * t1008;
t772 = t329 * t1022;
t317 = t329 * t959;
t318 = t329 * t957;
t771 = -t906 / 0.2e1;
t770 = t905 / 0.2e1;
t765 = t894 / 0.2e1;
t764 = -t892 / 0.2e1;
t762 = t517 * t960;
t761 = -t1011 / 0.2e1;
t759 = -t878 / 0.2e1;
t757 = t986 + t321 / 0.2e1;
t756 = t985 + t1060;
t749 = qJD(2) * t955;
t459 = t868 * t558;
t742 = t1003 * t602;
t741 = t801 * t600;
t737 = t1025 * t583;
t736 = t1004 * t786;
t734 = t898 * t975;
t733 = t517 * t766;
t731 = t495 + t768;
t730 = t596 * t749;
t729 = t796 / 0.2e1;
t726 = pkin(5) / 0.2e1 + t751;
t719 = -0.2e1 * t736;
t718 = 0.2e1 * t736;
t717 = t441 / 0.2e1 - t889 / 0.2e1;
t716 = -t442 / 0.2e1 + t888 / 0.2e1;
t712 = t597 * t741;
t711 = t801 * t890;
t11 = -t146 * t148 + t147 * t149 + t924;
t666 = t1043 * t986 + t1060 * t329;
t673 = t937 / 0.2e1 + t938 / 0.2e1;
t609 = t148 * t988 + t149 * t987 + t329 * t673 + t666;
t2 = t609 + t1039;
t710 = t2 * qJD(1) + t11 * qJD(2);
t13 = -t146 * t186 + t147 * t187 + t320 * t322;
t612 = -(t985 - t673) * t329 + t186 * t988 + t187 * t987 + t320 * t983;
t4 = t612 - t1041;
t709 = t4 * qJD(1) + t13 * qJD(2);
t708 = t517 * t736;
t12 = (t934 + t935) * t1004 + t1017;
t707 = -t12 * qJD(2) + t1051;
t614 = t321 * t983 + t329 * t985 + t547 * t766 + t666;
t662 = t404 * t966 + t405 * t965;
t14 = -t614 + t662;
t67 = t321 * t322 + t544 * t547 + t924;
t706 = -t14 * qJD(1) + t67 * qJD(2);
t16 = (t930 + t931) * t1004 + t1017;
t705 = -t16 * qJD(2) + t1051;
t692 = t1004 * t319 - t320 * t517;
t22 = -t1004 * t146 + t148 * t517 + t597 * t692;
t620 = t1023 * t275 + t597 * t633;
t31 = t770 + t620;
t704 = t31 * qJD(1) + t22 * qJD(2);
t23 = -t1004 * t147 - t149 * t517 + t600 * t692;
t619 = t1023 * t276 + t600 * t633;
t36 = t771 + t619;
t703 = t36 * qJD(1) + t23 * qJD(2);
t24 = (-t146 + t896) * t1004 + (t186 - t920) * t517;
t675 = (t988 + t765) * t1004;
t29 = t770 + t675;
t702 = t29 * qJD(1) + t24 * qJD(2);
t25 = (-t147 + t881) * t1004 + (-t187 - t919) * t517;
t674 = (-t276 / 0.2e1 + t1065) * t1004;
t34 = t771 + t674;
t701 = t34 * qJD(1) + t25 * qJD(2);
t28 = t921 - (t937 + t938) * t517;
t668 = t925 / 0.2e1 + t927 / 0.2e1;
t625 = -t517 * t668 - t772;
t664 = -t913 / 0.2e1 - t911 / 0.2e1;
t42 = t625 + t664;
t700 = qJD(1) * t42 + qJD(2) * t28;
t44 = -t1004 * t322 - t321 * t517 + t692;
t699 = t52 * qJD(1) + t44 * qJD(2);
t82 = -t1011 * t320 + t146 * t517;
t652 = -t1023 * t329 + t729;
t93 = t759 - t928 / 0.2e1 - t652 * t597;
t698 = qJD(1) * t93 - qJD(2) * t82;
t83 = -t147 * t517 + t320 * t876;
t92 = t764 + t926 / 0.2e1 + t652 * t600;
t697 = qJD(1) * t92 - qJD(2) * t83;
t695 = t933 - t936;
t694 = t929 - t932;
t690 = -t1004 * t553 - t552 * t517;
t160 = t984 + t1065;
t617 = t954 / 0.2e1 + pkin(5) * t1024 + t628;
t47 = t597 * t617 - t600 * t756;
t689 = qJD(1) * t160 + qJD(2) * t47;
t246 = t530 * t531 + t532 * t533 + t586 * t591;
t603 = t1054 * t980 + t531 * t981 + t532 * t979 + t766 * t591;
t621 = -t629 / 0.2e1;
t611 = t534 * t789 / 0.2e1 + t621 * t790;
t90 = -t603 + t611;
t688 = -t90 * qJD(1) + t246 * qJD(2);
t687 = t1004 * (-qJD(6) - t1014);
t665 = t1024 * t1043 - t772;
t111 = t729 - t665;
t113 = -t322 * t517 + t921;
t686 = qJD(1) * t111 - qJD(2) * t113;
t188 = -t1054 * t566 + (-t531 + t532) * t565;
t685 = -t117 * qJD(1) - t188 * qJD(2);
t605 = t1016 * t963 + t473 * t964;
t247 = t729 - t605;
t334 = -t531 * t566 + t533 * t565;
t684 = -qJD(1) * t247 + qJD(2) * t334;
t251 = t1004 * t544 + t517 * t547;
t683 = qJD(2) * t251 + t1050;
t252 = t1004 * t547 - t517 * t544;
t613 = t791 / 0.2e1 + t598 * t621;
t281 = -t733 + t613;
t682 = qJD(1) * t281 + qJD(2) * t252;
t635 = -t745 / 0.2e1 - t748 / 0.2e1;
t437 = (t566 / 0.2e1 + t635) * t898;
t505 = -t565 * t591 + t566 * t586;
t681 = -qJD(1) * t437 + qJD(2) * t505;
t636 = -t747 / 0.2e1 + t746 / 0.2e1;
t438 = (t964 + t636) * t898;
t506 = t565 * t586 + t566 * t591;
t680 = -qJD(1) * t438 + qJD(2) * t506;
t679 = qJD(5) * t424 + t852;
t676 = t949 / 0.2e1 + t1021 / 0.2e1;
t670 = t929 / 0.2e1 - t932 / 0.2e1;
t667 = t320 * t965 + t322 * t970;
t660 = t1004 * t971 + t517 * t968;
t279 = t734 + t613;
t659 = -qJD(1) * t279 + t517 * t857;
t658 = -t1004 * t857 - t1050;
t655 = (t844 + t852) * t517;
t654 = t1015 * t1004;
t653 = t377 / 0.2e1 + t676;
t610 = -(-t553 * t750 + t965) * t329 - t668 * t558 + t1043 * t970;
t10 = t610 - t1040;
t249 = -t459 * t553 + t552 * t559;
t6 = t953 / 0.2e1 + (pkin(10) * t989 + t147 * t967 + t187 * t968) * t600 + (pkin(10) * t990 + t146 * t967 + t186 * t969) * t597 + t667;
t651 = t10 * qJD(1) + t6 * qJD(2) + t249 * qJD(3);
t648 = t323 / 0.2e1 + t660;
t19 = (t187 / 0.2e1 + t989) * t600 + (-t186 / 0.2e1 + t990) * t597;
t647 = -qJD(2) * t19 + qJD(3) * t459;
t608 = t616 / 0.2e1 + t794 / 0.2e1;
t150 = t608 + t717;
t646 = qJD(2) * t150 + qJD(3) * t559;
t615 = -t627 / 0.2e1 - t795 / 0.2e1;
t151 = -t615 + t716;
t179 = t982 + t1031 / 0.2e1 + t754 * t598;
t645 = -qJD(1) * t179 - qJD(2) * t151 - qJD(3) * t558;
t425 = t726 * t597;
t88 = t653 * t600;
t642 = pkin(5) * t842 + qJD(2) * t88 + qJD(3) * t425;
t426 = t726 * t600;
t86 = t653 * t597;
t641 = pkin(5) * t841 - qJD(2) * t86 + qJD(3) * t426;
t71 = t597 * t757 - t600 * t648;
t639 = -qJD(2) * t71 - t552 * t850;
t157 = (t983 - t1043 / 0.2e1) * t597;
t50 = t597 * t756 + t600 * t617;
t638 = -qJD(1) * t157 - qJD(2) * t50 - t559 * t850;
t69 = t597 * t648 + t600 * t757;
t637 = -qJD(2) * t69 - t552 * t849;
t584 = t599 * t855;
t578 = t581 * qJD(6);
t563 = t565 * qJD(3);
t548 = t558 * qJD(5);
t539 = t597 * t549;
t461 = -0.2e1 * t1004 * t583;
t455 = t587 + t630;
t452 = t459 * qJD(5);
t440 = t566 * t766 + t635 * t898;
t439 = t565 * t766 + t636 * t898;
t428 = pkin(5) * t958 + t1005 * t957;
t427 = pkin(5) * t960 + t1005 * t959;
t391 = t719 + t1006;
t390 = t718 - t1006;
t389 = -t905 / 0.2e1;
t388 = t906 / 0.2e1;
t368 = t801 * t1002;
t354 = t763 + t761;
t353 = 0.2e1 * t763;
t352 = -t349 / 0.2e1 - t762;
t343 = t355 * qJD(5);
t340 = t352 * qJD(6);
t339 = t351 * qJD(6);
t331 = t341 - t840;
t299 = t711 - t860;
t298 = -t712 + t860;
t282 = -t733 + t870;
t280 = t734 + t870;
t248 = t729 + t605;
t216 = -t661 + t740;
t192 = t496 + t485 + t731;
t191 = t495 + t486 + t731;
t180 = t1036 + 0.2e1 * t982;
t165 = t329 * t958 - t318;
t164 = -0.2e1 * t318;
t163 = t329 * t960 - t317;
t162 = -0.2e1 * t317;
t161 = 0.2e1 * t984;
t156 = t1043 * t959 + t765;
t153 = -t608 + t717;
t152 = t615 + t716;
t116 = t117 * qJD(3);
t112 = t729 + t665;
t95 = -t926 / 0.2e1 - t1004 * t318 + t764 + t600 * t729;
t94 = t928 / 0.2e1 - t329 * t761 + t759 - t739 / 0.2e1;
t91 = t603 + t611;
t89 = t920 + t879 / 0.2e1 - t676 * t600;
t87 = t919 - t893 / 0.2e1 + t676 * t597;
t80 = 0.2e1 * t1064;
t72 = t920 / 0.2e1 - t897 / 0.2e1 + t880 / 0.2e1 - t660 * t600;
t70 = t919 / 0.2e1 - t882 / 0.2e1 - t895 / 0.2e1 + t660 * t597;
t55 = t640 + t671;
t49 = t896 / 0.2e1 + pkin(5) * t1037 + t923 / 0.2e1 + t999 * t600;
t48 = -t881 / 0.2e1 - pkin(5) * t762 - t922 / 0.2e1 + t999 * t597;
t41 = t625 - t664;
t39 = t634 - t672;
t35 = t388 + t619;
t33 = t388 + t674;
t32 = t389 + t620;
t30 = t389 + t675;
t18 = t670 + t1000;
t15 = t614 + t662;
t9 = t610 + t1040;
t5 = -t953 / 0.2e1 - t673 * t558 + t670 * t553 + t667 + t1000 * pkin(10);
t3 = t612 + t1041;
t1 = t609 - t1039;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t198, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t53 - t1068 * t801; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t730, -t777, 0, 0, 0, 0, 0, 0, 0, 0 (-t601 * t749 - t602 * t803) * t596 (t599 * t749 - t602 * t802) * t596, t742 * t856, t819 + (-pkin(2) * t955 + pkin(8) * t742) * t856, 0, 0, 0, 0, 0, 0, t440 * qJD(3) - t565 * t730, t439 * qJD(3) + t566 * t730 (t534 * t565 + t566 * t629) * qJD(2) + t116, t836 + (-t531 * t629 + t534 * t533 + t586 * t796) * qJD(2) + t91 * qJD(3) + t248 * qJD(4), 0, 0, 0, 0, 0, 0, t1029 * t801 + t517 * t730, t282 * qJD(3) + t280 * qJD(5) + t1004 * t730 (t1004 * t404 - t405 * t517) * qJD(2) + t943, t873 + (t405 * t322 + t544 * t796 + t908) * qJD(2) + t15 * qJD(3) + t112 * qJD(4), 0, 0, 0, 0, 0, 0 (t1011 * t404 + t379 * t517) * qJD(2) + t32 * qJD(3) + t30 * qJD(5) + t95 * qJD(6) (-t380 * t517 + t404 * t876) * qJD(2) + t35 * qJD(3) + t33 * qJD(5) + t94 * qJD(6) (-t911 - t913) * t1008 + t801 * t1027, t901 + (-t146 * t379 + t147 * t380 + t908) * qJD(2) + t1 * qJD(3) + t41 * qJD(4) + t3 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t643 - t599 * t777, qJD(3) * t644 - t601 * t777, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t440 - qJD(3) * t473, qJD(2) * t439 - qJD(3) * t1016, t839, t91 * qJD(2) + (t1016 * t939 - t473 * t940) * t942, 0, 0, 0, 0, 0, 0, t1067, qJD(2) * t282 - qJD(3) * t329 + qJD(5) * t180, t902, t15 * qJD(2) + (t1043 * t558 + t329 * t559) * qJD(3), 0, 0, 0, 0, 0, 0, qJD(2) * t32 + qJD(5) * t161 + qJD(6) * t163 - t1043 * t849, qJD(2) * t35 + qJD(5) * t156 + qJD(6) * t165 + t1043 * t850, qJD(3) * t1063 + t80 * qJD(5) + t1049, -t1069 + t1 * qJD(2) + (t1043 * t552 + t1063 * t553) * qJD(3) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, qJD(2) * t280 + qJD(3) * t180 - qJD(5) * t329, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t30 + qJD(3) * t161 + qJD(6) * t162 - t1043 * t841, qJD(2) * t33 + qJD(3) * t156 + qJD(6) * t164 + t1043 * t842, t80 * qJD(3) + qJD(5) * t1063 + t1049, -t1069 + t3 * qJD(2) + t9 * qJD(3) + (pkin(10) * t1063 - t1059) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t95 + qJD(3) * t163 + qJD(5) * t162 - qJD(6) * t276, qJD(2) * t94 + qJD(3) * t165 + qJD(5) * t164 + qJD(6) * t275, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t819, 0, 0, 0, 0, 0, 0, -t437 * qJD(3), -t438 * qJD(3), t116, -qJD(3) * t90 - qJD(4) * t247 - t836, 0, 0, 0, 0, 0, 0, t801 * t1030, qJD(3) * t281 + qJD(5) * t279, t943, -qJD(3) * t14 - qJD(4) * t111 - t873, 0, 0, 0, 0, 0, 0, qJD(3) * t31 + qJD(5) * t29 - qJD(6) * t92, qJD(3) * t36 + qJD(5) * t34 - qJD(6) * t93, t801 * t1028, qJD(3) * t2 + qJD(4) * t42 + qJD(5) * t4 - t901; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t776, t582 * qJD(3), 0, -t776, 0, 0, -pkin(2) * t803, -pkin(2) * t802, 0, 0, t778, t460 * qJD(3), 0, -t778, 0, 0, t505 * qJD(3), t506 * qJD(3), qJD(3) * t188 + qJD(4) * t527, qJD(3) * t246 + qJD(4) * t334, -t654, -t801 * t1026, 0, t655, 0, 0, qJD(3) * t251 + t1004 * t843, qJD(3) * t252 - t517 * t843, qJD(3) * t44 + qJD(4) * t207, qJD(3) * t67 + qJD(4) * t113, -t594 * t654 - t737, t1015 * t1061 * t876 - qJD(6) * t370, t1052 * t801 - t597 * t783, -t592 * t654 + t737, -t1053 * t801 - t600 * t783, t655, qJD(3) * t22 + qJD(4) * t195 + qJD(5) * t24 + qJD(6) * t83, qJD(3) * t23 + qJD(4) * t197 + qJD(5) * t25 + qJD(6) * t82, -qJD(3) * t12 - qJD(5) * t16, qJD(3) * t11 + qJD(4) * t28 + qJD(5) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t584, t804, t802, -t584, -t803, 0, -pkin(8) * t802 - t799, pkin(8) * t803 - t798, 0, 0, t779, t815, t563, -t779, -t805, 0, qJD(3) * t530 + t681, -qJD(3) * t532 + t680 (-t565 * t940 - t566 * t939) * t942 - t685 (t530 * t940 + t532 * t939) * t942 + t455 * qJD(4) + t688, -t780, -t1045, -t1015, t1033, -t679, 0, -qJD(3) * t319 + qJD(5) * t153 + t683, -qJD(3) * t321 + qJD(5) * t152 + t682 (-t1004 * t559 - t558 * t517) * qJD(3) + t699 (t319 * t558 + t321 * t559) * qJD(3) + t216 * qJD(4) + t706 -(t785 + t787) * t517 - t997, 0.2e1 * t708 + (t508 - t1020) * qJD(3) + t191 * qJD(5) + t461, qJD(5) * t353 + t1004 * t850 + t1062 -(-t785 + t788) * t517 + t997, t340 - t1056 + t871, t1044 (t597 * t690 - t922) * qJD(3) + t48 * qJD(5) + t72 * qJD(6) + t704 (t600 * t690 + t923) * qJD(3) + t49 * qJD(5) + t70 * qJD(6) + t703, qJD(3) * t695 + t18 * qJD(5) + t707 (t319 * t552 + t553 * t695) * qJD(3) + t39 * qJD(4) + t5 * qJD(5) + t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t808, qJD(3) * t455 + t684, 0, 0, 0, 0, 0, 0, t811, 0, t835, qJD(3) * t216 - t686, 0, 0, 0, 0, 0, 0, t340 - t343 + t865, qJD(5) * t354 + t1055 + t863, t192 * qJD(5), qJD(3) * t39 + qJD(5) * t55 + t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1033, -t1045, -t1015, t1033, -qJD(3) * t424 - t844, 0, qJD(3) * t153 + qJD(4) * t514 - qJD(5) * t322 - t658, qJD(3) * t152 + qJD(5) * t320 - t659, 0, 0 (-t784 - t787) * t517 - t995, t191 * qJD(3) + t461 + (-qJD(5) * t581 + t718) * t517, qJD(3) * t353 + t1004 * t842 + t1062 (t784 - t788) * t517 + t995, qJD(3) * t356 + t1004 * t841 - t1056, t1044, t48 * qJD(3) - t355 * qJD(4) + (t597 * t725 - t881) * qJD(5) + t89 * qJD(6) + t702, t49 * qJD(3) + t354 * qJD(4) + (t600 * t725 + t896) * qJD(5) + t87 * qJD(6) + t701, t18 * qJD(3) + t192 * qJD(4) + qJD(5) * t694 + t705, t5 * qJD(3) + t55 * qJD(4) + (-pkin(5) * t322 + pkin(10) * t694) * qJD(5) + t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t1007 * t712 - t859, t1042 * t801 + t597 * t687, t130, qJD(3) * t352 + t600 * t687, t368, qJD(3) * t72 + qJD(4) * t352 + qJD(5) * t89 - qJD(6) * t147 - t697, qJD(3) * t70 + qJD(4) * t1042 + qJD(5) * t87 + qJD(6) * t146 - t698, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t437 * qJD(2), t438 * qJD(2), -t839, qJD(2) * t90, 0, 0, 0, 0, 0, 0, -t1046, -qJD(2) * t281 + qJD(5) * t179, -t902, qJD(2) * t14, 0, 0, 0, 0, 0, 0, -qJD(2) * t31 + qJD(5) * t160, -qJD(2) * t36 + qJD(5) * t157, -t1048, -qJD(2) * t2 + qJD(5) * t10 + t1069; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t584, -t804, 0, t584, 0, 0, t799, t798, 0, 0, -t779, -t815, 0, t779, 0, 0, -qJD(4) * t566 - t681, -qJD(4) * t565 - t680, t685, qJD(4) * t456 - t688, t780, t1045, 0, -t1033, -t811, 0, -qJD(5) * t150 - t683 - t848, qJD(5) * t151 + t1013 - t682, -t699, -qJD(4) * t215 - t706, t594 * t780 - t998, qJD(5) * t190 + t461 - 0.2e1 * t708, -qJD(5) * t347 + qJD(6) * t360 - t1057, t592 * t780 + t998, t339 - t343 + t1056, -t1044, qJD(5) * t47 + qJD(6) * t71 - t600 * t848 - t704, qJD(4) * t1011 + qJD(5) * t50 + qJD(6) * t69 - t703, -qJD(4) * t369 + qJD(5) * t19 - t707, qJD(4) * t40 + qJD(5) * t6 - t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t549, t548, 0, 0, t583, t578, 0, -t583, 0, 0, -t549 * t600 + t552 * t840, t552 * t590 + t539, -t452, t249 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t806, -t807, 0, t816, 0, 0, 0, 0, 0, 0, -t1008, t1014, 0, -t833, 0, 0, 0, 0, 0, 0, -t775, t829, -t820, t909; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t812, 0, -t549 - t646, t548 - t645, 0, 0, t993, t578 + t867, -t828, -t993, -t825, 0, qJD(6) * t427 - t559 * t741 + t689, qJD(6) * t428 + t539 - t638, -t452 - t647 (-pkin(5) * t559 - pkin(10) * t459) * qJD(5) + t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t391, t590 + t822, t298, t331, -t1009, qJD(5) * t427 - t553 * t590 - t639, qJD(5) * t428 + t553 * t840 - t637, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t111 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t805, t563, -t808, -qJD(3) * t456 - t684, 0, 0, 0, 0, 0, 0, t679, -t1015, -t835, qJD(3) * t215 + t686, 0, 0, 0, 0, 0, 0, t339 - t865 + t871, -qJD(3) * t1011 - qJD(5) * t348 + qJD(6) * t361 - t863, qJD(3) * t369 + qJD(5) * t193, -qJD(3) * t40 - qJD(5) * t54 - t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t806, t807, 0, -t816, 0, 0, 0, 0, 0, 0, t1008, -t1014, 0, t833, 0, 0, 0, 0, 0, 0, t775, -t829, t820, -t909; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t817, -t1014, 0, 0, 0, 0, 0, 0, 0, 0, t824, -t827, t837, -t900; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, -t590 + t821, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1046, -qJD(2) * t279 - qJD(3) * t179, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t29 - qJD(3) * t160, -qJD(2) * t34 - qJD(3) * t157, -t1048, -qJD(2) * t4 - qJD(3) * t10 + t1069; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1033, t1045, 0, -t1033, t514 * qJD(3), 0, qJD(3) * t150 - qJD(4) * t424 + t658, -qJD(3) * t151 + t1013 + t659, 0, 0, t1033 * t594 - t996, -qJD(3) * t190 + t517 * t719 + t461, qJD(3) * t347 + qJD(6) * t1019 - t1057, t1033 * t592 + t996, qJD(3) * t355 - qJD(6) * t349 + t1056, -t1044, -qJD(3) * t47 - qJD(4) * t356 - qJD(6) * t88 - t702, -qJD(3) * t50 + qJD(4) * t348 + qJD(6) * t86 - t701, -qJD(3) * t19 - qJD(4) * t193 - t705, -qJD(3) * t6 + qJD(4) * t54 - t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t812, 0, t646, t645, 0, 0, -t994, t578 - t867, t828, t994, t825, 0, -qJD(6) * t425 + t559 * t849 - t689, -qJD(6) * t426 + t638, t647, -t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t817, t1014, 0, 0, 0, 0, 0, 0, 0, 0, -t824, t827, -t837, t900; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583, t578, 0, -t583, 0, 0, -pkin(5) * t840, -pkin(5) * t590, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t391, t590 + t823, t298, -t826 - t840, -t1009, -pkin(10) * t590 - t642, pkin(10) * t840 - t641, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t92, qJD(2) * t93, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t1007 * t711 + t859, -qJD(3) * t360 - qJD(5) * t1019 + t1033 * t597, -t130, -qJD(3) * t351 + qJD(5) * t349 + t1033 * t600, t368, -qJD(3) * t71 - qJD(4) * t351 + qJD(5) * t88 + t697, -qJD(3) * t69 - qJD(4) * t361 - qJD(5) * t86 + t698, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, t390, -t822, t299, -t341, t1009, qJD(5) * t425 + t639, qJD(5) * t426 + t637, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341, -t821, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, t390, -t823, t299, t826, t1009, t642, t641, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
