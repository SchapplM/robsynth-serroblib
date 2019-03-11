% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRP9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:34:02
% EndTime: 2019-03-09 12:34:53
% DurationCPUTime: 30.58s
% Computational Cost: add. (29336->949), mult. (69472->1315), div. (0->0), fcn. (78822->10), ass. (0->684)
t616 = sin(pkin(11));
t618 = cos(pkin(11));
t619 = cos(pkin(6));
t617 = sin(pkin(6));
t622 = sin(qJ(2));
t915 = t617 * t622;
t561 = -t616 * t915 + t618 * t619;
t562 = t616 * t619 + t618 * t915;
t621 = sin(qJ(4));
t970 = cos(qJ(4));
t448 = -t970 * t561 + t621 * t562;
t623 = cos(qJ(5));
t1021 = t623 * t448;
t430 = -t1021 / 0.2e1;
t772 = t1021 / 0.2e1;
t1031 = t772 + t430;
t1032 = qJD(5) * t1031;
t620 = sin(qJ(5));
t1022 = t620 * t448;
t780 = t1022 / 0.2e1;
t808 = t970 * t618;
t902 = t621 * t616;
t678 = t808 - t902;
t624 = cos(qJ(2));
t914 = t617 * t624;
t783 = -t914 / 0.2e1;
t742 = t621 * t783;
t747 = t970 * t914;
t975 = t618 / 0.2e1;
t414 = t616 * t742 - t678 * t783 + t747 * t975;
t870 = qJD(1) * t624;
t805 = t617 * t870;
t745 = t448 * t805;
t1030 = t414 * qJD(2) - t745;
t809 = t970 * t616;
t901 = t621 * t618;
t581 = t809 + t901;
t976 = -t616 / 0.2e1;
t413 = t581 * t783 + t618 * t742 + t747 * t976;
t537 = t970 * t562;
t903 = t621 * t561;
t1014 = t537 + t903;
t746 = t1014 * t805;
t1029 = t413 * qJD(2) - t746;
t868 = qJD(2) * t581;
t873 = qJD(1) * t1014;
t1028 = t868 + t873;
t830 = t1014 * qJD(4);
t1024 = t448 / 0.2e1;
t609 = -pkin(3) * t618 - pkin(2);
t738 = -pkin(4) * t678 - pkin(10) * t581;
t668 = t609 + t738;
t960 = pkin(9) + qJ(3);
t758 = t960 * t616;
t740 = t621 * t758;
t589 = t960 * t618;
t810 = t970 * t589;
t502 = t810 - t740;
t907 = t620 * t502;
t322 = -t623 * t668 + t907;
t921 = t581 * t623;
t287 = -qJ(6) * t921 - t322;
t963 = t678 * pkin(5);
t257 = t287 - t963;
t605 = pkin(8) * t914;
t606 = t619 * t622 * pkin(1);
t569 = t605 + t606;
t534 = qJ(3) * t619 + t569;
t731 = -pkin(2) * t624 - qJ(3) * t622;
t535 = (-pkin(1) + t731) * t617;
t406 = t618 * t534 + t616 * t535;
t331 = pkin(9) * t561 + t406;
t812 = t970 * t331;
t405 = -t616 * t534 + t618 * t535;
t310 = -pkin(3) * t914 - t562 * pkin(9) + t405;
t904 = t621 * t310;
t207 = t812 + t904;
t198 = -pkin(10) * t914 + t207;
t969 = pkin(1) * t624;
t568 = pkin(8) * t915 - t619 * t969;
t538 = -pkin(2) * t619 + t568;
t459 = -pkin(3) * t561 + t538;
t739 = pkin(4) * t448 - pkin(10) * t1014;
t633 = t459 + t739;
t110 = t198 * t620 - t623 * t633;
t898 = t623 * t1014;
t383 = -t620 * t914 + t898;
t93 = -qJ(6) * t383 - t110;
t965 = t448 * pkin(5);
t75 = t93 + t965;
t986 = -t678 / 0.2e1;
t1027 = t1024 * t257 + t75 * t986;
t206 = -t970 * t310 + t331 * t621;
t197 = pkin(4) * t914 + t206;
t1026 = t197 / 0.2e1;
t1025 = -t448 / 0.2e1;
t1023 = t1014 / 0.2e1;
t615 = t623 ^ 2;
t977 = -t615 / 0.2e1;
t614 = t620 ^ 2;
t979 = -t614 / 0.2e1;
t737 = (t977 + t979) * t448;
t869 = qJD(2) * t678;
t874 = qJD(1) * t448;
t1020 = -t869 + t874;
t909 = t620 * t1014;
t381 = t623 * t914 + t909;
t1019 = t1014 * t381;
t1018 = t1014 * t383;
t927 = t448 * t581;
t251 = t1014 * t986 + t927 / 0.2e1;
t1017 = t251 * qJD(4);
t831 = t448 * qJD(4);
t804 = t678 * t868;
t676 = qJD(1) * t251 - t804;
t677 = qJD(2) * t251 + t448 * t873;
t765 = t287 / 0.2e1 - t257 / 0.2e1;
t111 = t623 * t198 + t620 * t633;
t94 = -t381 * qJ(6) + t111;
t1016 = t765 * t94;
t1015 = -t257 + t287;
t367 = t383 * t620;
t900 = t623 * t381;
t243 = t367 - t900;
t528 = t537 / 0.2e1;
t443 = t528 + t903 / 0.2e1;
t575 = t678 ^ 2;
t576 = t581 ^ 2;
t1013 = -t576 - t575;
t823 = t576 - t575;
t1012 = t383 ^ 2;
t1011 = t448 ^ 2;
t1010 = -pkin(4) / 0.2e1;
t1009 = t75 / 0.2e1;
t293 = pkin(4) * t1014 + pkin(10) * t448;
t290 = t623 * t293;
t913 = t620 * t206;
t756 = t290 + t913;
t97 = pkin(5) * t1014 + qJ(6) * t1021 + t756;
t1008 = t97 / 0.2e1;
t815 = t616 * t914;
t515 = pkin(3) * t815 + t569;
t517 = t581 * t914;
t518 = t678 * t914;
t308 = pkin(4) * t517 - pkin(10) * t518 + t515;
t303 = t623 * t308;
t566 = (pkin(2) * t622 - qJ(3) * t624) * t617;
t457 = t618 * t566 + t616 * t568;
t376 = (-pkin(9) * t618 * t624 + pkin(3) * t622) * t617 + t457;
t345 = t621 * t376;
t458 = t616 * t566 - t618 * t568;
t401 = -pkin(9) * t815 + t458;
t392 = t970 * t401;
t890 = t392 + t345;
t236 = pkin(10) * t915 + t890;
t912 = t620 * t236;
t154 = t303 - t912;
t895 = t623 * t518;
t461 = t620 * t915 + t895;
t124 = pkin(5) * t517 - qJ(6) * t461 + t154;
t1007 = t124 / 0.2e1;
t233 = t623 * t236;
t302 = t620 * t308;
t155 = t233 + t302;
t594 = t623 * t915;
t905 = t620 * t518;
t460 = -t594 + t905;
t139 = -qJ(6) * t460 + t155;
t1006 = -t139 / 0.2e1;
t473 = t623 * t678;
t961 = t581 * pkin(4);
t962 = t678 * pkin(10);
t491 = t961 - t962;
t484 = t623 * t491;
t501 = t589 * t621 + t970 * t758;
t924 = t501 * t620;
t755 = t484 + t924;
t264 = pkin(5) * t581 - qJ(6) * t473 + t755;
t1004 = t264 / 0.2e1;
t897 = t623 * t502;
t323 = t620 * t668 + t897;
t470 = t620 * t581;
t288 = -qJ(6) * t470 + t323;
t1002 = t288 / 0.2e1;
t791 = -t367 / 0.2e1;
t1001 = -t381 / 0.2e1;
t1000 = t381 / 0.2e1;
t999 = -t383 / 0.2e1;
t998 = t383 / 0.2e1;
t991 = -t484 / 0.2e1;
t990 = t501 / 0.2e1;
t989 = -t502 / 0.2e1;
t988 = -t537 / 0.2e1;
t987 = t678 / 0.2e1;
t985 = -t581 / 0.2e1;
t984 = t581 / 0.2e1;
t959 = -qJ(6) - pkin(10);
t591 = t959 * t620;
t983 = t591 / 0.2e1;
t592 = t959 * t623;
t982 = t592 / 0.2e1;
t981 = -t592 / 0.2e1;
t610 = -pkin(5) * t623 - pkin(4);
t980 = t610 / 0.2e1;
t978 = t614 / 0.2e1;
t974 = -t620 / 0.2e1;
t973 = t620 / 0.2e1;
t972 = -t623 / 0.2e1;
t971 = t623 / 0.2e1;
t968 = pkin(5) * t383;
t967 = pkin(5) * t461;
t966 = pkin(5) * t620;
t964 = t460 * pkin(5);
t950 = t94 * t623;
t954 = t75 * t620;
t958 = -t954 / 0.2e1 + t950 / 0.2e1;
t957 = t75 - t93;
t956 = pkin(5) * qJD(5);
t955 = pkin(5) * qJD(6);
t953 = t75 * t623;
t952 = t94 * t678;
t951 = t94 * t620;
t153 = pkin(5) * t381 + t197;
t16 = t153 * t968 - t94 * t957;
t949 = qJD(1) * t16;
t17 = t957 * t381;
t948 = qJD(1) * t17;
t29 = t153 * t1014 - (t950 - t954) * t448;
t947 = qJD(1) * t29;
t40 = -t381 * t94 - t383 * t75;
t946 = qJD(1) * t40;
t45 = t110 * t448 - t197 * t381;
t945 = qJD(1) * t45;
t46 = -t111 * t448 + t197 * t383;
t944 = qJD(1) * t46;
t943 = t124 * t623;
t203 = t623 * t206;
t289 = t620 * t293;
t892 = t203 - t289;
t109 = qJ(6) * t1022 - t892;
t158 = -pkin(5) * t1022 + t207;
t13 = t109 * t94 + t153 * t158 + t75 * t97;
t942 = t13 * qJD(1);
t941 = t139 * t620;
t391 = t621 * t401;
t811 = t970 * t376;
t733 = -t391 + t811;
t235 = -pkin(4) * t915 - t733;
t194 = t235 + t964;
t18 = t124 * t75 + t139 * t94 + t153 * t194;
t940 = t18 * qJD(1);
t19 = -t109 * t381 - t97 * t383 + (t951 + t953) * t448;
t939 = t19 * qJD(1);
t938 = t197 * t623;
t937 = t257 * t620;
t26 = -t124 * t383 - t139 * t381 - t460 * t94 - t461 * t75;
t936 = t26 * qJD(1);
t935 = t264 * t623;
t934 = t288 * t623;
t468 = t620 * t678;
t483 = t620 * t491;
t490 = t501 * t623;
t889 = t490 - t483;
t294 = -qJ(6) * t468 - t889;
t933 = t294 * t620;
t32 = -t110 * t1014 + t207 * t381 + (t290 + (-t197 + t206) * t620) * t448;
t932 = t32 * qJD(1);
t33 = -t111 * t1014 + t207 * t383 + (t892 - t938) * t448;
t931 = t33 * qJD(1);
t34 = -t110 * t517 + t154 * t448 + t197 * t460 + t235 * t381;
t930 = t34 * qJD(1);
t35 = -t111 * t517 - t155 * t448 + t197 * t461 + t235 * t383;
t929 = t35 * qJD(1);
t928 = t383 * t623;
t926 = t461 * t620;
t925 = t461 * t623;
t923 = t561 * t618;
t922 = t562 * t616;
t920 = t591 * t620;
t919 = t591 * t623;
t918 = t592 * t620;
t917 = t592 * t623;
t612 = t617 ^ 2;
t916 = t612 * t622;
t911 = t620 * t381;
t908 = t620 * t460;
t906 = t620 * t517;
t896 = t623 * t517;
t66 = t206 * t915 - t515 * t448 - t459 * t517 + t733 * t914;
t894 = t66 * qJD(1);
t67 = t515 * t1014 + t459 * t518 + (-t207 * t622 + t624 * t890) * t617;
t893 = t67 * qJD(1);
t891 = -t937 / 0.2e1 + t934 / 0.2e1;
t888 = t809 / 0.2e1 + t901 / 0.2e1;
t599 = t616 ^ 2 + t618 ^ 2;
t602 = t614 + t615;
t603 = t615 - t614;
t141 = -t206 * t914 - t459 * t448;
t887 = qJD(1) * t141;
t142 = -t1014 * t459 - t207 * t914;
t886 = qJD(1) * t142;
t157 = t243 * t448;
t885 = qJD(1) * t157;
t160 = -t1011 * t620 + t1019;
t884 = qJD(1) * t160;
t161 = -t1022 * t448 - t1019;
t883 = qJD(1) * t161;
t162 = -t1011 * t623 + t1018;
t882 = qJD(1) * t162;
t163 = t1021 * t448 + t1018;
t881 = qJD(1) * t163;
t224 = -t405 * t562 + t406 * t561;
t880 = qJD(1) * t224;
t879 = qJD(1) * t1022;
t763 = 0.2e1 * t1024;
t275 = t763 * t620;
t878 = qJD(1) * t275;
t305 = t448 * t915 - t517 * t914;
t877 = qJD(1) * t305;
t306 = -t1014 * t915 + t518 * t914;
t876 = qJD(1) * t306;
t875 = qJD(1) * t383;
t514 = t888 * t914;
t872 = qJD(1) * t514;
t871 = qJD(1) * t619;
t867 = qJD(2) * t609;
t866 = qJD(2) * t616;
t865 = qJD(2) * t618;
t864 = qJD(2) * t623;
t863 = qJD(2) * t624;
t862 = qJD(3) * t623;
t861 = qJD(3) * t624;
t860 = qJD(4) * t581;
t859 = qJD(4) * t620;
t858 = qJD(4) * t623;
t857 = qJD(5) * t381;
t856 = qJD(5) * t448;
t855 = qJD(5) * t620;
t854 = qJD(5) * t623;
t641 = t1025 * t609 + t459 * t987 + t501 * t783;
t764 = -t345 / 0.2e1 - t392 / 0.2e1;
t133 = t641 - t764;
t853 = t133 * qJD(1);
t640 = t459 * t984 + t609 * t1023 + t502 * t914 / 0.2e1;
t693 = t391 / 0.2e1 - t811 / 0.2e1;
t135 = t640 + t693;
t852 = t135 * qJD(1);
t144 = -t457 * t562 + t458 * t561 + (-t405 * t618 - t406 * t616) * t914;
t851 = t144 * qJD(1);
t773 = t900 / 0.2e1;
t666 = (t773 + t791) * t678;
t778 = t908 / 0.2e1;
t789 = t925 / 0.2e1;
t147 = t778 + t789 - t666;
t850 = t147 * qJD(1);
t152 = t405 * t457 + t406 * t458 + t538 * t569;
t849 = t152 * qJD(1);
t156 = (t900 + t367) * t448;
t848 = t156 * qJD(1);
t159 = -t381 * t461 - t383 * t460;
t847 = t159 * qJD(1);
t769 = t473 / 0.2e1;
t786 = t581 * t977;
t662 = t383 * t769 + t448 * t786;
t173 = -t926 / 0.2e1 + t662;
t846 = t173 * qJD(1);
t785 = -t915 / 0.2e1;
t792 = t381 * t987;
t180 = t792 - t895 / 0.2e1 + (-t927 / 0.2e1 + t785) * t620;
t845 = t180 * qJD(1);
t671 = t383 * t987 + t430 * t581;
t732 = t594 / 0.2e1 - t905 / 0.2e1;
t181 = t671 - t732;
t844 = t181 * qJD(1);
t186 = -t405 * t915 + t457 * t914 - t538 * t815 + t569 * t561;
t843 = t186 * qJD(1);
t187 = t569 * t562 + (-t406 * t622 + (t538 * t618 + t458) * t624) * t617;
t842 = t187 * qJD(1);
t192 = -t381 * t517 - t448 * t460;
t841 = t192 * qJD(1);
t193 = t383 * t517 + t448 * t461;
t840 = t193 * qJD(1);
t223 = -t1014 * t517 - t448 * t518;
t839 = t223 * qJD(1);
t281 = 0.2e1 * t430;
t838 = t281 * qJD(1);
t412 = (t985 + t888) * t914;
t837 = t412 * qJD(1);
t836 = t412 * qJD(2);
t835 = t413 * qJD(1);
t833 = t414 * qJD(1);
t415 = (t987 - t808 / 0.2e1 + t902 / 0.2e1) * t914;
t832 = t415 * qJD(1);
t829 = t468 * qJD(2);
t510 = pkin(1) * t916 + t569 * t619;
t828 = t510 * qJD(1);
t511 = -t568 * t619 + t612 * t969;
t827 = t511 * qJD(1);
t826 = t514 * qJD(5);
t567 = t678 * qJD(4);
t580 = (-t622 ^ 2 + t624 ^ 2) * t612;
t825 = t580 * qJD(1);
t824 = t617 * qJD(4);
t822 = pkin(5) * t921;
t821 = pkin(5) * t855;
t820 = t968 / 0.2e1;
t819 = t966 / 0.2e1;
t818 = -t963 / 0.2e1;
t817 = -t93 / 0.2e1 + t1009;
t816 = t623 * t470;
t813 = -t950 / 0.2e1;
t807 = t383 * t874;
t803 = t615 * t868;
t802 = t581 * t864;
t801 = t617 * t861;
t800 = t624 * t824;
t799 = t620 * t858;
t798 = t581 * t855;
t797 = t581 * t854;
t796 = t581 * t567;
t795 = t612 * t870;
t794 = qJD(2) * t915;
t793 = t620 * t854;
t790 = -t925 / 0.2e1;
t788 = -t921 / 0.2e1;
t787 = t921 / 0.2e1;
t784 = t915 / 0.2e1;
t782 = -t911 / 0.2e1;
t781 = -t1022 / 0.2e1;
t779 = -t909 / 0.2e1;
t777 = -t906 / 0.2e1;
t776 = t906 / 0.2e1;
t775 = -t470 / 0.2e1;
t774 = t470 / 0.2e1;
t771 = t898 / 0.2e1;
t770 = -t896 / 0.2e1;
t768 = t1026 - t206 / 0.2e1;
t767 = -t203 / 0.2e1 + t289 / 0.2e1;
t766 = -t233 / 0.2e1 - t302 / 0.2e1;
t762 = 0.2e1 * t1023;
t761 = t483 / 0.2e1 - t490 / 0.2e1;
t760 = -t605 / 0.2e1 - t606 / 0.2e1;
t757 = 0.2e1 * t816;
t754 = t829 - t855;
t753 = pkin(5) * t787;
t752 = t620 * t818;
t751 = qJD(2) + t871;
t750 = -qJD(5) - t874;
t749 = qJD(3) + t867;
t748 = -qJD(5) + t869;
t744 = t863 * t916;
t743 = t622 * t795;
t741 = t581 * t979 + t786;
t736 = qJD(2) * t757;
t735 = qJD(4) * t757;
t734 = t303 / 0.2e1 - t912 / 0.2e1;
t730 = -qJD(4) + t805;
t435 = pkin(5) * t470 + t501;
t436 = pkin(5) * t468 + t502;
t626 = t109 * t1002 + t153 * t436 / 0.2e1 + t158 * t435 / 0.2e1 + t75 * t1004 + t94 * t294 / 0.2e1 + t257 * t1008;
t649 = t124 * t983 + t139 * t981 + t194 * t980;
t1 = -t626 + t649;
t53 = t257 * t264 + t288 * t294 + t435 * t436;
t729 = -t1 * qJD(1) + t53 * qJD(2);
t652 = t288 * t1024 + t109 * t985 - t952 / 0.2e1;
t657 = t97 * t985 + t1027;
t681 = t460 * t981 + t461 * t983;
t685 = t1001 * t294 + t264 * t999;
t4 = (t1006 + t657) * t623 + (t1007 + t652) * t620 + t681 + t685;
t723 = t257 * t623 + t288 * t620;
t57 = (t933 + t935) * t581 + t723 * t678;
t728 = t4 * qJD(1) - t57 * qJD(2);
t5 = -t1016 + t817 * t288 + (t153 * t788 + t435 * t999 + t1007) * pkin(5);
t59 = t1015 * t288 + t435 * t822;
t727 = -qJD(1) * t5 + qJD(2) * t59;
t58 = t1015 * t470;
t638 = -t381 * t765 + t470 * t817;
t8 = t967 / 0.2e1 + t638;
t726 = -qJD(1) * t8 + qJD(2) * t58;
t14 = (t965 / 0.2e1 + t817) * t623;
t689 = t818 - t765;
t54 = t689 * t623;
t725 = qJD(1) * t14 + qJD(2) * t54;
t20 = (pkin(5) * t1024 + t817) * t620;
t64 = t689 * t620;
t724 = qJD(1) * t20 + qJD(2) * t64;
t722 = -t457 * t616 + t458 * t618;
t721 = -t918 + t919;
t106 = t435 * t581 - (-t934 + t937) * t678;
t687 = t1023 * t435 + t153 * t984;
t12 = (t288 * t1025 + t952 / 0.2e1 - t124 / 0.2e1) * t623 + (t1006 + t1027) * t620 + t687;
t720 = -qJD(1) * t12 - qJD(2) * t106;
t680 = -t917 / 0.2e1 - t920 / 0.2e1;
t636 = t581 * t980 + t678 * t680;
t684 = t935 / 0.2e1 + t933 / 0.2e1;
t125 = t636 - t684;
t635 = t1014 * t980 - t448 * t680;
t688 = t109 * t973 + t97 * t971;
t41 = t635 - t688;
t719 = -qJD(1) * t41 - qJD(2) * t125;
t131 = -t484 * t678 + (-t322 + t907) * t581;
t650 = t1023 * t322 + t110 * t984 + t381 * t989;
t690 = -pkin(10) * t517 / 0.2e1 + t207 * t985;
t692 = t1010 * t460 + t235 * t972;
t22 = t448 * t991 + t290 * t987 + (-t678 * t768 + t690) * t620 + t650 + t692;
t718 = -t22 * qJD(1) + t131 * qJD(2);
t132 = (-t323 + t897) * t581 - (t889 - t490) * t678;
t625 = t1023 * t323 + t1025 * t889 + t111 * t984 + t383 * t989 + t892 * t987;
t691 = t1010 * t461 + t235 * t973;
t23 = (t197 * t986 + t448 * t990 + t690) * t623 + t625 + t691;
t717 = -t23 * qJD(1) + t132 * qJD(2);
t145 = t723 * t581;
t634 = t964 / 0.2e1 + pkin(4) * t785 + t693;
t686 = t1000 * t288 + t257 * t998;
t27 = (t951 / 0.2e1 + t953 / 0.2e1) * t581 + t634 + t686;
t716 = -qJD(1) * t27 - qJD(2) * t145;
t225 = -t322 * t678 - t470 * t501;
t631 = t1025 * t322 + t110 * t987 + t197 * t774 + t381 * t990;
t37 = t631 + t766;
t715 = -qJD(1) * t37 + qJD(2) * t225;
t226 = t323 * t678 + t501 * t921;
t630 = t1025 * t323 + t111 * t987 + t197 * t787 + t383 * t990;
t36 = -t630 + t734;
t714 = qJD(1) * t36 - qJD(2) * t226;
t713 = t748 * t623;
t644 = (t779 + t1001) * t581 - t448 * t468;
t117 = t770 + t644;
t407 = t823 * t620;
t712 = -qJD(1) * t117 + qJD(2) * t407;
t661 = t1014 * t984 - t678 * t763;
t629 = t381 * t984 + t620 * t661;
t119 = t770 + t629;
t408 = t1013 * t620;
t711 = qJD(1) * t119 - qJD(2) * t408;
t643 = (t771 + t998) * t581 + t448 * t473;
t121 = t777 + t643;
t409 = t823 * t623;
t710 = -qJD(1) * t121 - qJD(2) * t409;
t628 = t383 * t984 + t623 * t661;
t123 = t776 + t628;
t482 = t1013 * t623;
t709 = qJD(1) * t123 - qJD(2) * t482;
t165 = -0.2e1 * t737;
t564 = t614 * t678;
t565 = t615 * t678;
t479 = -t564 - t565;
t708 = qJD(1) * t165 + qJD(2) * t479;
t175 = -t1014 * t581 - t448 * t678;
t228 = -t1014 ^ 2 + t1011;
t707 = qJD(1) * t228 + qJD(2) * t175;
t706 = qJD(1) * t175 - qJD(2) * t823;
t627 = (t923 / 0.2e1 + t922 / 0.2e1) * qJ(3) + t405 * t976 + t406 * t975;
t178 = t627 + t760;
t586 = t599 * qJ(3);
t705 = qJD(1) * t178 + qJD(2) * t586;
t683 = -t928 / 0.2e1 + t782;
t196 = -t903 / 0.2e1 + t988 + t683;
t396 = -t888 + t741;
t704 = qJD(1) * t196 + qJD(2) * t396;
t215 = (-t911 - t928) * t581;
t480 = t602 * t576;
t703 = qJD(1) * t215 - qJD(2) * t480;
t278 = t762 * t620;
t702 = -qJD(1) * t278 - qJD(2) * t470;
t701 = -qJD(1) * t1021 + qJD(2) * t473;
t698 = qJD(1) * t443 + qJD(2) * t888;
t446 = t922 + t923;
t474 = t561 ^ 2 + t562 ^ 2;
t697 = qJD(1) * t474 + qJD(2) * t446;
t696 = qJD(1) * t446 + qJD(2) * t599;
t695 = qJD(1) * t243 + qJD(4) * t602;
t694 = -t962 / 0.2e1 + t961 / 0.2e1;
t682 = t381 * t981 + t383 * t983;
t679 = t773 + t367 / 0.2e1;
t283 = t762 * t623;
t675 = qJD(1) * t283 + t802;
t674 = t802 + t875;
t673 = pkin(4) * t999 + pkin(10) * t430;
t672 = t694 * t623;
t670 = t904 / 0.2e1 + t812 / 0.2e1;
t394 = t610 * t966;
t43 = t765 * t592 + (t435 * t974 + t610 * t788 + t1004) * pkin(5);
t9 = -t817 * t592 + (t153 * t974 + t610 * t999 + t1008) * pkin(5);
t667 = -qJD(1) * t9 - qJD(2) * t43 + qJD(4) * t394;
t632 = -t810 / 0.2e1 + t752 + t740 / 0.2e1;
t642 = (t918 / 0.2e1 - t919 / 0.2e1) * t581 + t891;
t100 = t632 + t642;
t30 = t813 + (-t965 / 0.2e1 + t1009) * t620 + t670 + t682;
t513 = -t917 - t920;
t665 = -qJD(1) * t30 + qJD(2) * t100 + qJD(4) * t513;
t151 = (-t911 + t928) * t581;
t167 = -t900 + 0.2e1 * t791;
t380 = t381 ^ 2;
t185 = t380 - t1012;
t664 = qJD(1) * t185 - qJD(2) * t151 + qJD(4) * t167;
t252 = t380 + t1012;
t663 = qJD(1) * t252 - qJD(2) * t215 + qJD(4) * t243;
t651 = t448 * t816 - t678 * t679;
t108 = t778 + t790 + t651;
t660 = t108 * qJD(1) - t678 * t736;
t639 = t694 * t620 + t490 / 0.2e1;
t217 = t639 + t761;
t648 = pkin(4) * t1000 + t938 / 0.2e1 + pkin(10) * t780;
t47 = t648 + t767;
t659 = pkin(4) * t858 - qJD(1) * t47 - qJD(2) * t217;
t219 = t991 - t672;
t49 = -t290 / 0.2e1 + t768 * t620 + t673;
t658 = pkin(4) * t859 - qJD(1) * t49 - qJD(2) * t219;
t214 = t679 * t581;
t240 = t782 + t928 / 0.2e1;
t656 = qJD(2) * t214 - qJD(4) * t240 + t381 * t875;
t467 = (t978 + t977) * t581;
t655 = qJD(1) * t240 - qJD(2) * t467 + t799;
t654 = qJD(5) * t443 + t677;
t653 = qJD(5) * t888 + t676;
t481 = t603 * t576;
t647 = qJD(1) * t151 + qJD(2) * t481 + t735;
t646 = -qJD(1) * t167 - qJD(4) * t603 + t736;
t645 = t576 * t620 * t864 + qJD(1) * t214 + qJD(4) * t467;
t637 = (qJD(2) * t731 + t861) * t617;
t593 = qJD(2) * t784;
t560 = t888 * qJD(4);
t558 = t581 * t858;
t539 = (t795 - t824 / 0.2e1) * t622;
t524 = -0.2e1 * t581 * t793;
t498 = t896 / 0.2e1;
t463 = t468 * qJD(5);
t462 = t467 * qJD(5);
t444 = t528 + t988;
t438 = -t908 / 0.2e1;
t437 = t446 * qJD(3);
t404 = t415 * qJD(2);
t395 = t741 + t888;
t313 = qJD(2) * t514 + qJD(4) * t443;
t301 = (-t674 - t859) * pkin(5);
t284 = t1014 * t972 + t771;
t277 = t1014 * t973 + t779;
t276 = t781 + t780;
t239 = t243 * qJD(6);
t237 = t240 * qJD(5);
t220 = t924 + t484 / 0.2e1 - t672;
t218 = t639 - t761;
t213 = t215 * qJD(6);
t208 = t214 * qJD(5);
t195 = t683 + t443;
t182 = t671 + t732;
t179 = t448 * t775 + t792 + t895 / 0.2e1 + t620 * t784;
t177 = t627 - t760;
t174 = t175 * qJD(4);
t172 = t926 / 0.2e1 + t662;
t166 = t167 * qJD(5);
t164 = t1025 * t615 - t448 * t978 - t737;
t150 = t151 * qJD(5);
t146 = t438 + t790 - t666;
t136 = t640 - t693;
t134 = t641 + t764;
t126 = t636 + t684;
t122 = t777 + t628;
t120 = t776 + t643;
t118 = t498 + t629;
t116 = t498 + t644;
t107 = t438 + t789 + t651;
t99 = -t632 + t642;
t65 = t287 * t973 - t934 / 0.2e1 + t752 + t891;
t55 = (t765 + t818) * t623;
t50 = t197 * t973 + t913 / 0.2e1 + t290 / 0.2e1 + t673;
t48 = t648 - t767;
t44 = pkin(5) * t1004 + t257 * t982 + t287 * t981 + t435 * t819 + t610 * t753;
t42 = t635 + t688;
t39 = t630 + t734;
t38 = -t631 + t766;
t31 = pkin(5) * t781 + t670 - t682 + t958;
t28 = t75 * t788 + t775 * t94 + t634 - t686;
t25 = pkin(10) * t770 + t197 * t769 + t207 * t787 + t430 * t501 - t625 + t691;
t24 = pkin(10) * t777 + t1024 * t755 + t1026 * t468 + t207 * t774 + t501 * t781 + t756 * t986 - t650 + t692;
t21 = pkin(5) * t780 + t93 * t973 + t813 + t958;
t15 = t93 * t971 - t953 / 0.2e1 + pkin(5) * t772;
t11 = t941 / 0.2e1 + t943 / 0.2e1 - (t813 + t954 / 0.2e1) * t678 - t891 * t448 + t687;
t10 = pkin(5) * t1008 + t153 * t819 + t610 * t820 + t75 * t982 + t93 * t981;
t7 = -t967 / 0.2e1 + t638;
t6 = t93 * t1002 - t75 * t288 / 0.2e1 + t435 * t820 + t153 * t753 + pkin(5) * t1007 + t1016;
t3 = t124 * t974 + t139 * t971 + t620 * t652 + t623 * t657 - t681 + t685;
t2 = t626 + t649;
t51 = [0, 0, 0, t744, t580 * qJD(2), t617 * t619 * t863, -t619 * t794, 0, -t510 * qJD(2), -t511 * qJD(2), -t186 * qJD(2) + t562 * t801, t187 * qJD(2) + t561 * t801, qJD(2) * t144 + qJD(3) * t474, qJD(2) * t152 + qJD(3) * t224 (qJD(2) * t518 - t831) * t1014, qJD(2) * t223 + qJD(4) * t228, -t306 * qJD(2) + t448 * t800, -t305 * qJD(2) + t1014 * t800, -t744, -t66 * qJD(2) - t142 * qJD(4) + t1014 * t801, t67 * qJD(2) + t141 * qJD(4) - t448 * t801 (qJD(2) * t461 - t623 * t831 - t857) * t383, qJD(2) * t159 + qJD(4) * t156 + qJD(5) * t185, qJD(2) * t193 + qJD(4) * t162 - t381 * t856, qJD(2) * t192 - qJD(4) * t160 - t383 * t856 (qJD(2) * t517 + t830) * t448, qJD(2) * t34 - qJD(3) * t161 + qJD(4) * t32 + qJD(5) * t46, qJD(2) * t35 + qJD(3) * t163 + qJD(4) * t33 + qJD(5) * t45, qJD(2) * t26 - qJD(3) * t157 + qJD(4) * t19 + qJD(5) * t17 + qJD(6) * t252, qJD(2) * t18 + qJD(3) * t29 + qJD(4) * t13 + qJD(5) * t16 + qJD(6) * t40; 0, 0, 0, t743, t825, t751 * t914, -t751 * t915, 0, -qJD(2) * t569 - t828, qJD(2) * t568 - t827, -t569 * t865 + t616 * t637 - t843, t569 * t866 + t618 * t637 + t842, qJD(2) * t722 + t437 + t851, t849 + (-t569 * pkin(2) + qJ(3) * t722) * qJD(2) + t177 * qJD(3), t1028 * t518 - t1017, t839 + (-t517 * t581 + t518 * t678) * qJD(2) + t174, -qJD(4) * t415 + t581 * t794 - t876, -qJD(4) * t412 + t678 * t794 - t877, -t539, -t894 + (-t501 * t915 - t515 * t678 + t517 * t609) * qJD(2) - t413 * qJD(3) + t136 * qJD(4), t893 + (-t502 * t915 + t515 * t581 + t518 * t609) * qJD(2) + t414 * qJD(3) + t134 * qJD(4), t172 * qJD(4) + t461 * t674 - t208, t847 + t107 * qJD(4) - t150 + (-t460 * t623 - t926) * t868, t840 + (-t461 * t678 + t581 * t896) * qJD(2) + t120 * qJD(4) + t179 * qJD(5), t841 + (t460 * t678 - t470 * t517) * qJD(2) + t116 * qJD(4) + t182 * qJD(5), t1020 * t517 + t1017 + t826, t930 + (-t154 * t678 + t235 * t470 - t322 * t517 + t460 * t501) * qJD(2) + t118 * qJD(3) + t24 * qJD(4) + t39 * qJD(5), t929 + (t155 * t678 + t235 * t921 - t323 * t517 + t461 * t501) * qJD(2) + t122 * qJD(3) + t25 * qJD(4) + t38 * qJD(5), t936 + (-t257 * t461 - t288 * t460 + (-t941 - t943) * t581) * qJD(2) + t146 * qJD(3) + t3 * qJD(4) + t7 * qJD(5) - t213, t940 + (t124 * t257 + t139 * t288 + t194 * t435) * qJD(2) + t11 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) + t28 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(1) * t562 + t866) * t914 (qJD(1) * t561 + t865) * t914, t697, qJD(2) * t177 + t880, 0, 0, 0, 0, 0, t444 * qJD(4) - t1029, t1030, 0, 0, 0, 0, 0, qJD(2) * t118 + qJD(4) * t284 + qJD(5) * t276 - t883, qJD(2) * t122 + qJD(4) * t277 + t1032 + t881, qJD(2) * t146 + qJD(4) * t164 - t885, qJD(2) * t11 + qJD(4) * t42 + qJD(5) * t21 + qJD(6) * t195 + t947; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t677, t707, t448 * t730 - t404, t1014 * t730 - t836, t593, qJD(2) * t136 + qJD(3) * t444 - qJD(4) * t207 - t886, qJD(2) * t134 + qJD(4) * t206 + t887, t172 * qJD(2) + t237 + (-t859 - t875) * t1021, t107 * qJD(2) - t603 * t831 + t166 + t848, qJD(2) * t120 + t620 * t830 + t1032 + t882, qJD(2) * t116 + t623 * t830 - t884, t654, t932 + t24 * qJD(2) + t284 * qJD(3) + (-t207 * t623 + t620 * t739) * qJD(4) + t50 * qJD(5), t931 + t25 * qJD(2) + t277 * qJD(3) + (t207 * t620 + t623 * t739) * qJD(4) + t48 * qJD(5), t939 + t3 * qJD(2) + t164 * qJD(3) + (t109 * t623 + t448 * t721 - t97 * t620) * qJD(4) + t15 * qJD(5) + t239, t942 + t2 * qJD(2) + t42 * qJD(3) + (-t109 * t592 + t158 * t610 + t591 * t97) * qJD(4) + t10 * qJD(5) + t31 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t656, t664, t179 * qJD(2) + qJD(4) * t1031 + t381 * t750, t182 * qJD(2) + t383 * t750, t313, qJD(2) * t39 + qJD(3) * t276 + qJD(4) * t50 - qJD(5) * t111 + t944, qJD(2) * t38 + qJD(3) * t1031 + qJD(4) * t48 + qJD(5) * t110 + t945, pkin(5) * t857 + qJD(2) * t7 + qJD(4) * t15 + t948, qJD(2) * t6 + qJD(3) * t21 + qJD(4) * t10 - t94 * t956 + t949; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t663, qJD(2) * t28 + qJD(3) * t195 + qJD(4) * t31 + t946; 0, 0, 0, -t743, -t825, -t619 * t805, t871 * t915, 0, t828, t827, t843, -t842, t437 - t851, qJD(3) * t178 - t849, -t518 * t873 - t1017, t174 - t839, -qJD(4) * t414 + t876, -qJD(4) * t413 + t877, t539, -qJD(3) * t412 + qJD(4) * t135 + t894, qJD(3) * t415 + qJD(4) * t133 - t893, qJD(4) * t173 - t461 * t875 - t208, qJD(4) * t108 - t150 - t847, qJD(4) * t121 + qJD(5) * t180 - t840, qJD(4) * t117 + qJD(5) * t181 - t841, -t517 * t874 + t1017 - t826, qJD(3) * t119 - qJD(4) * t22 - qJD(5) * t36 - t930, qJD(3) * t123 - qJD(4) * t23 - qJD(5) * t37 - t929, qJD(3) * t147 + qJD(4) * t4 + qJD(5) * t8 - t213 - t936, qJD(3) * t12 - qJD(4) * t1 - qJD(5) * t5 - qJD(6) * t27 - t940; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599 * qJD(3), t586 * qJD(3), t796, -t823 * qJD(4), 0, 0, 0, t609 * t860, t609 * t567, -t576 * t793 + t615 * t796, -qJD(5) * t481 - t678 * t735, qJD(4) * t409 + t678 * t798, -qJD(4) * t407 + t678 * t797, -t796, -qJD(3) * t408 + qJD(4) * t131 + qJD(5) * t226, -qJD(3) * t482 + qJD(4) * t132 + qJD(5) * t225, -qJD(4) * t57 - qJD(5) * t58 + qJD(6) * t480, qJD(3) * t106 + qJD(4) * t53 + qJD(5) * t59 - qJD(6) * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t696, t705, 0, 0, 0, 0, 0, -t837, t832, 0, 0, 0, 0, 0, t711, t709, t850, qJD(4) * t126 + qJD(5) * t65 + qJD(6) * t395 - t720; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t676, t706, t567 - t833, -t835 - t860, qJD(1) * t785, -qJD(4) * t502 + t581 * t867 + t852, qJD(4) * t501 + t678 * t867 + t853, t846 - t462 - (-t799 - t803) * t678 (-t564 + t565) * qJD(4) + t524 + t660, t581 * t859 - t710, t558 - t712, t653 (t620 * t738 - t897) * qJD(4) + t220 * qJD(5) + t718 (t623 * t738 + t907) * qJD(4) + t218 * qJD(5) + t717 (-t264 * t620 + t294 * t623 - t678 * t721) * qJD(4) + t55 * qJD(5) + t728, t126 * qJD(3) + (t264 * t591 - t294 * t592 + t436 * t610) * qJD(4) + t44 * qJD(5) + t99 * qJD(6) + t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t645, -t647, t470 * t748 + t845, t581 * t713 + t844, t560 - t872, qJD(4) * t220 - qJD(5) * t323 - t714, qJD(4) * t218 + qJD(5) * t322 + t715, pkin(5) * t798 + qJD(4) * t55 - t726, qJD(3) * t65 + qJD(4) * t44 - t288 * t956 + t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t703, qJD(3) * t395 + qJD(4) * t99 + t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t562 * t805, -t561 * t805, -t697, -qJD(2) * t178 - t880, 0, 0, 0, 0, 0, -t746 + t836 + t830, -t404 + t745 - t831, 0, 0, 0, 0, 0, -qJD(2) * t119 + qJD(4) * t283 - qJD(5) * t275 + t883, -qJD(2) * t123 - qJD(4) * t278 + qJD(5) * t281 - t881, -qJD(2) * t147 + qJD(4) * t165 + t885, -qJD(2) * t12 - qJD(4) * t41 - qJD(5) * t20 + qJD(6) * t196 - t947; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t696, -t705, 0, 0, 0, 0, 0, t837 + t860, t567 - t832, 0, 0, 0, 0, 0, t463 + t558 - t711, -qJD(4) * t470 + t678 * t854 - t709, qJD(4) * t479 - t850, -qJD(4) * t125 - qJD(5) * t64 + qJD(6) * t396 + t720; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1028, -t1020, 0, 0, 0, 0, 0, t675, t702, t708, t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t754 - t878, t713 + t838, 0, -t724 - t821; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t677, -t707, t1030, t1029, t593, -qJD(2) * t135 - qJD(3) * t1014 + t886, -qJD(2) * t133 + qJD(3) * t448 - t887, -qJD(2) * t173 + t623 * t807 + t237, -qJD(2) * t108 + t166 - t848, -qJD(2) * t121 + qJD(5) * t1021 - t882, -qJD(2) * t117 - qJD(5) * t1022 + t884, -t654, qJD(2) * t22 - qJD(3) * t283 + qJD(5) * t49 - t932, qJD(2) * t23 + qJD(3) * t278 + qJD(5) * t47 - t931, -qJD(2) * t4 - qJD(3) * t165 - qJD(5) * t14 + t239 - t939, qJD(2) * t1 + qJD(3) * t41 - qJD(5) * t9 - qJD(6) * t30 - t942; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t676, -t706, t833, t835, qJD(1) * t784, -t581 * t749 - t852, -t678 * t749 - t853, -t678 * t803 - t462 - t846, t524 - t660, -qJD(5) * t473 + t710, t463 + t712, -t653, qJD(5) * t219 - t581 * t862 - t718, qJD(3) * t470 + qJD(5) * t217 - t717, -qJD(3) * t479 - qJD(5) * t54 - t728, qJD(3) * t125 - qJD(5) * t43 + qJD(6) * t100 - t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1028, t1020, 0, 0, 0, 0, 0, -t675, -t702, -t708, -t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t793, t603 * qJD(5), 0, 0, 0, -pkin(4) * t855, -pkin(4) * t854, qJD(6) * t602, qJD(5) * t394 + qJD(6) * t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t655, -t646, -t701 + t854, t754 - t879, -t698, -pkin(10) * t854 - t658, pkin(10) * t855 - t659, -pkin(5) * t854 - t725, t592 * t956 + t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t695, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t656, -t664, -qJD(2) * t180 - qJD(4) * t1021 + t381 * t874, -qJD(2) * t181 + qJD(4) * t1022 + t807, t313, qJD(2) * t36 + qJD(3) * t275 - qJD(4) * t49 - t944, qJD(2) * t37 - qJD(3) * t281 - qJD(4) * t47 - t945, -qJD(2) * t8 + qJD(4) * t14 - t948, qJD(2) * t5 + qJD(3) * t20 + qJD(4) * t9 - t383 * t955 - t949; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t645, t647, qJD(4) * t473 - t620 * t804 - t845, -qJD(4) * t468 - t678 * t802 - t844, t560 + t872, -qJD(3) * t468 - qJD(4) * t219 + t714, -qJD(4) * t217 - t678 * t862 - t715, qJD(4) * t54 + t726, qJD(3) * t64 + qJD(4) * t43 - qJD(6) * t822 - t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t829 + t878, -t678 * t864 - t838, 0, t724; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t655, t646, t701, -t829 + t879, t698, t658, t659, t725, -t620 * t955 - t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t663, qJD(2) * t27 - qJD(3) * t196 + qJD(4) * t30 + t383 * t956 - t946; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t703, pkin(5) * t797 - qJD(3) * t396 - qJD(4) * t100 - t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t695, -t665 + t821; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t51;
