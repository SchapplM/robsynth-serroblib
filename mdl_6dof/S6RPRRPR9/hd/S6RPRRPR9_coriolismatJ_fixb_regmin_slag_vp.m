% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:38
% EndTime: 2019-03-09 05:32:30
% DurationCPUTime: 31.63s
% Computational Cost: add. (43629->948), mult. (119722->1378), div. (0->0), fcn. (141058->14), ass. (0->698)
t638 = sin(pkin(6));
t639 = cos(pkin(12));
t642 = sin(qJ(3));
t637 = sin(pkin(7));
t998 = cos(pkin(6));
t788 = t998 * t637;
t774 = t642 * t788;
t997 = cos(pkin(7));
t789 = t642 * t997;
t1008 = cos(qJ(3));
t636 = sin(pkin(12));
t833 = t1008 * t636;
t536 = t774 + (t639 * t789 + t833) * t638;
t950 = t638 * t639;
t586 = t637 * t950 - t997 * t998;
t641 = sin(qJ(4));
t644 = cos(qJ(4));
t454 = t536 * t641 + t586 * t644;
t456 = t536 * t644 - t586 * t641;
t635 = sin(pkin(13));
t996 = cos(pkin(13));
t356 = t454 * t996 + t456 * t635;
t643 = cos(qJ(6));
t1055 = t643 * t356;
t1060 = t1055 / 0.2e1;
t1061 = -t1055 / 0.2e1;
t1062 = t1060 + t1061;
t1063 = qJD(6) * t1062;
t951 = t637 * t642;
t588 = t641 * t997 + t644 * t951;
t708 = t641 * t951 - t644 * t997;
t653 = t588 * t996 - t635 * t708;
t1033 = -t653 / 0.2e1;
t1037 = t653 / 0.2e1;
t1059 = t1033 + t1037;
t1005 = -qJ(5) - pkin(10);
t614 = t1005 * t644;
t787 = t996 * t641;
t550 = -t1005 * t787 - t614 * t635;
t1056 = t550 / 0.2e1;
t1057 = -t550 / 0.2e1;
t792 = t1056 + t1057;
t1043 = t356 / 0.2e1;
t483 = t588 * t635 + t708 * t996;
t1058 = t483 / 0.2e1;
t442 = t996 * t456;
t955 = t635 * t454;
t1053 = t442 - t955;
t1045 = t1053 / 0.2e1;
t640 = sin(qJ(6));
t190 = t640 * t1053;
t954 = t635 * t641;
t688 = t1005 * t954 - t614 * t996;
t1054 = -t688 / 0.2e1;
t967 = t688 * t640;
t966 = t688 * t643;
t786 = t996 * t644;
t599 = -t786 + t954;
t795 = 0.2e1 * t1043;
t765 = t795 * t599;
t591 = t599 ^ 2;
t953 = t635 * t644;
t601 = t787 + t953;
t592 = t601 ^ 2;
t1052 = -t592 - t591;
t631 = t640 ^ 2;
t633 = t643 ^ 2;
t622 = t633 - t631;
t1051 = qJD(4) * t622;
t882 = qJD(3) * t643;
t827 = t640 * t882;
t780 = t601 * t827;
t1050 = 0.2e1 * t780 - t1051;
t772 = t997 * t1008;
t952 = t636 * t638;
t534 = -t1008 * t788 + t642 * t952 - t772 * t950;
t1049 = t534 ^ 2;
t836 = pkin(1) * t998;
t621 = t639 * t836;
t537 = t998 * pkin(2) + t621 + (-pkin(9) * t997 - qJ(2)) * t952;
t568 = (-pkin(9) * t636 * t637 - pkin(2) * t639 - pkin(1)) * t638;
t434 = -t537 * t637 + t568 * t997;
t766 = pkin(3) * t534 - pkin(10) * t536;
t328 = t434 + t766;
t497 = t537 * t789;
t914 = qJ(2) * t950 + t636 * t836;
t526 = (t950 * t997 + t788) * pkin(9) + t914;
t499 = t1008 * t526;
t542 = t568 * t951;
t371 = t499 + t497 + t542;
t345 = -pkin(10) * t586 + t371;
t181 = -t328 * t644 + t345 * t641;
t161 = -qJ(5) * t456 - t181;
t152 = pkin(4) * t534 + t161;
t182 = t328 * t641 + t345 * t644;
t162 = -qJ(5) * t454 + t182;
t956 = t635 * t162;
t79 = t152 * t996 - t956;
t1048 = -t79 / 0.2e1;
t159 = t996 * t162;
t80 = t152 * t635 + t159;
t1047 = t80 / 0.2e1;
t100 = t161 * t996 - t956;
t1046 = t100 / 0.2e1;
t1044 = -t356 / 0.2e1;
t628 = -pkin(4) * t644 - pkin(3);
t707 = pkin(5) * t599 - pkin(11) * t601 + t628;
t418 = t640 * t707 + t966;
t1042 = t418 / 0.2e1;
t1041 = -t456 / 0.2e1;
t1040 = t456 / 0.2e1;
t834 = t637 * t1008;
t463 = t640 * t653 + t643 * t834;
t1039 = -t463 / 0.2e1;
t464 = -t640 * t834 + t643 * t653;
t1038 = -t464 / 0.2e1;
t1036 = -t483 / 0.2e1;
t1032 = -t536 / 0.2e1;
t771 = t1008 * t996;
t835 = t635 * t1008;
t562 = (t641 * t771 + t644 * t835) * t637;
t1025 = t562 / 0.2e1;
t1024 = -t599 / 0.2e1;
t1023 = t599 / 0.2e1;
t1022 = -t601 / 0.2e1;
t1021 = t601 / 0.2e1;
t626 = pkin(4) * t635 + pkin(11);
t1020 = t626 / 0.2e1;
t627 = -pkin(4) * t996 - pkin(5);
t1019 = -t627 / 0.2e1;
t1018 = t627 / 0.2e1;
t1017 = -t628 / 0.2e1;
t1016 = -t635 / 0.2e1;
t1015 = t635 / 0.2e1;
t1014 = -t640 / 0.2e1;
t1013 = t640 / 0.2e1;
t1012 = -t641 / 0.2e1;
t1011 = t641 / 0.2e1;
t1010 = -t643 / 0.2e1;
t1009 = t643 / 0.2e1;
t1007 = t456 * pkin(4);
t1006 = t641 * pkin(4);
t1004 = pkin(4) * qJD(4);
t269 = -t534 * t643 + t190;
t370 = t526 * t642 - t537 * t772 - t568 * t834;
t344 = pkin(3) * t586 + t370;
t226 = pkin(4) * t454 + t344;
t652 = pkin(5) * t356 - pkin(11) * t1053 + t226;
t68 = pkin(11) * t534 + t80;
t39 = t640 * t68 - t643 * t652;
t67 = -pkin(5) * t534 - t79;
t169 = pkin(5) * t1053 + pkin(11) * t356 + t1007;
t931 = t643 * t169;
t948 = t640 * t100;
t757 = t931 - t948;
t943 = t640 * t356;
t99 = t161 * t635 + t159;
t5 = -t1053 * t39 + t99 * t269 + t356 * t757 - t67 * t943;
t1003 = t5 * qJD(1);
t271 = t1053 * t643 + t534 * t640;
t40 = t640 * t652 + t643 * t68;
t933 = t643 * t100;
t946 = t640 * t169;
t758 = t933 + t946;
t6 = -t1053 * t40 - t1055 * t67 + t99 * t271 - t356 * t758;
t1002 = t6 * qJD(1);
t412 = t644 * t534;
t427 = pkin(3) * t536 + pkin(10) * t534;
t416 = t644 * t427;
t935 = t641 * t370;
t188 = pkin(4) * t536 + qJ(5) * t412 + t416 + t935;
t409 = t641 * t534;
t369 = t644 * t370;
t415 = t641 * t427;
t918 = t369 - t415;
t204 = qJ(5) * t409 - t918;
t133 = t188 * t996 - t204 * t635;
t117 = -pkin(5) * t536 - t133;
t524 = t536 * t643;
t408 = t409 * t635 - t534 * t786;
t940 = t640 * t408;
t363 = -t524 + t940;
t407 = t601 * t534;
t514 = pkin(4) * t409;
t298 = -t514 + t371;
t170 = -pkin(5) * t407 - pkin(11) * t408 + t298;
t930 = t643 * t170;
t134 = t188 * t635 + t204 * t996;
t118 = pkin(11) * t536 + t134;
t947 = t640 * t118;
t49 = t930 - t947;
t9 = t117 * t269 + t356 * t49 + t363 * t67 + t39 * t407;
t1001 = t9 * qJD(1);
t1000 = t99 * t640;
t999 = t99 * t643;
t25 = -t269 * t67 + t356 * t39;
t995 = qJD(1) * t25;
t26 = t271 * t67 - t356 * t40;
t994 = qJD(1) * t26;
t34 = -t1053 * t79 - t356 * t80;
t993 = qJD(1) * t34;
t650 = t1037 * t269 + t1039 * t1053;
t843 = t637 * t952;
t781 = t641 * t843;
t773 = t636 * t789;
t832 = t1008 * t639;
t570 = (-t773 + t832) * t638;
t919 = t644 * t570;
t493 = t781 + t919;
t842 = t644 * t952;
t934 = t641 * t570;
t711 = t637 * t842 - t934;
t390 = t493 * t635 - t711 * t996;
t977 = t390 * t643;
t51 = t977 / 0.2e1 + t650;
t992 = qJD(1) * t51;
t649 = t1037 * t271 + t1038 * t1053;
t978 = t390 * t640;
t54 = -t978 / 0.2e1 + t649;
t991 = qJD(1) * t54;
t924 = t643 * t408;
t968 = t536 * t640;
t364 = t924 + t968;
t932 = t643 * t118;
t945 = t640 * t170;
t50 = t932 + t945;
t10 = t117 * t271 - t356 * t50 + t364 * t67 + t40 * t407;
t990 = t10 * qJD(1);
t22 = t100 * t80 + t1007 * t226 - t79 * t99;
t989 = t22 * qJD(1);
t23 = -t1053 * t133 - t134 * t356 + t407 * t80 - t408 * t79;
t988 = t23 * qJD(1);
t24 = t133 * t79 + t134 * t80 + t226 * t298;
t987 = t24 * qJD(1);
t986 = t271 * t1053;
t985 = t271 * t640;
t984 = t271 * t643;
t391 = t493 * t996 + t635 * t711;
t734 = t636 * t772;
t949 = t639 * t642;
t569 = (t734 + t949) * t638;
t33 = t226 * t569 - t390 * t79 + t391 * t80;
t983 = t33 * qJD(1);
t982 = t344 * t644;
t981 = t1053 * t269;
t980 = t356 * t601;
t979 = t364 * t640;
t563 = (-t641 * t835 + t644 * t771) * t637;
t937 = t640 * t563;
t495 = t643 * t951 - t937;
t662 = t1025 * t269 - t1039 * t407 + t1043 * t495 + t1058 * t363;
t920 = t643 * t569;
t942 = t640 * t391;
t754 = t920 - t942;
t520 = t640 * t601;
t803 = -t520 / 0.2e1;
t666 = t1024 * t754 + t390 * t803;
t44 = t662 + t666;
t976 = t44 * qJD(1);
t921 = t643 * t563;
t496 = t640 * t951 + t921;
t661 = t1025 * t271 - t1038 * t407 + t1044 * t496 + t1058 * t364;
t926 = t643 * t391;
t936 = t640 * t569;
t755 = t926 + t936;
t958 = t601 * t643;
t815 = t958 / 0.2e1;
t667 = t1024 * t755 + t390 * t815;
t46 = t661 - t667;
t975 = t46 * qJD(1);
t974 = t463 * t599;
t973 = t464 * t599;
t47 = -t181 * t536 + t371 * t454 + (t416 + (-t344 + t370) * t641) * t534;
t972 = t47 * qJD(1);
t48 = -t182 * t536 + t371 * t456 + (t918 - t982) * t534;
t971 = t48 * qJD(1);
t970 = t483 * t601;
t969 = t653 * t599;
t929 = t643 * t269;
t718 = t929 / 0.2e1 + t985 / 0.2e1;
t680 = t1055 * t520 + t599 * t718;
t727 = t1009 * t364 + t1014 * t363;
t56 = t680 - t727;
t965 = t56 * qJD(1);
t964 = t562 * t640;
t963 = t562 * t643;
t962 = t569 * t644;
t961 = t588 * t534;
t960 = t599 * t601;
t959 = t601 * t269;
t659 = t1025 * t1053 - t1033 * t407 + t1044 * t563 + t1058 * t408;
t726 = t1021 * t390 + t1024 * t391;
t62 = t659 - t726;
t957 = t62 * qJD(1);
t944 = t640 * t269;
t941 = t640 * t407;
t525 = pkin(5) * t601 + pkin(11) * t599 + t1006;
t939 = t640 * t525;
t938 = t550 * t640;
t925 = t643 * t407;
t923 = t643 * t525;
t922 = t550 * t643;
t713 = t638 * t734;
t916 = t713 / 0.2e1 + t638 * t949 / 0.2e1;
t809 = -t950 / 0.2e1;
t915 = -t713 / 0.2e1 + t642 * t809;
t623 = -t641 ^ 2 + t644 ^ 2;
t845 = t356 * t943;
t107 = -t845 - t981;
t913 = qJD(1) * t107;
t108 = -t845 + t981;
t912 = qJD(1) * t108;
t844 = t356 * t1055;
t109 = -t844 + t986;
t911 = qJD(1) * t109;
t110 = t844 + t986;
t910 = qJD(1) * t110;
t720 = -t942 / 0.2e1 + t920 / 0.2e1;
t724 = t1038 * t356 + t1058 * t271;
t119 = t720 - t724;
t909 = qJD(1) * t119;
t719 = -t936 / 0.2e1 - t926 / 0.2e1;
t725 = t1036 * t269 + t1043 * t463;
t120 = t719 - t725;
t908 = qJD(1) * t120;
t126 = t269 * t390 + t356 * t754;
t907 = qJD(1) * t126;
t127 = t271 * t390 - t356 * t755;
t906 = qJD(1) * t127;
t131 = t181 * t534 - t344 * t454;
t905 = qJD(1) * t131;
t132 = -t182 * t534 + t344 * t456;
t904 = qJD(1) * t132;
t144 = t1053 * t390 - t356 * t391;
t903 = qJD(1) * t144;
t728 = -t1033 * t356 + t1036 * t1053;
t153 = t728 + t916;
t902 = qJD(1) * t153;
t202 = -t370 * t586 - t434 * t534;
t901 = qJD(1) * t202;
t203 = t371 * t586 + t434 * t536;
t900 = qJD(1) * t203;
t899 = qJD(1) * t271;
t275 = -t1049 * t644 + t456 * t536;
t898 = qJD(1) * t275;
t282 = t454 * t569 + t534 * t711;
t897 = qJD(1) * t282;
t283 = t456 * t569 - t493 * t534;
t896 = qJD(1) * t283;
t895 = qJD(1) * t356;
t399 = t534 * t843 + t569 * t586;
t894 = qJD(1) * t399;
t400 = t536 * t843 + t570 * t586;
t893 = qJD(1) * t400;
t589 = t787 / 0.2e1 + t953 / 0.2e1;
t406 = t589 * t534;
t892 = qJD(1) * t406;
t891 = qJD(1) * t456;
t890 = qJD(1) * t534;
t889 = qJD(1) * t586;
t888 = qJD(3) * t534;
t887 = qJD(3) * t586;
t886 = qJD(3) * t599;
t885 = qJD(3) * t601;
t884 = qJD(3) * t641;
t883 = qJD(3) * t642;
t881 = qJD(3) * t644;
t880 = qJD(4) * t534;
t879 = qJD(4) * t640;
t878 = qJD(4) * t641;
t877 = qJD(4) * t643;
t876 = qJD(4) * t644;
t875 = qJD(5) * t643;
t874 = qJD(6) * t356;
t873 = qJD(6) * t640;
t872 = qJD(6) * t643;
t124 = -t929 - t985;
t101 = t124 * t356;
t871 = t101 * qJD(1);
t102 = -t269 * t364 - t271 * t363;
t870 = t102 * qJD(1);
t135 = t269 * t407 - t356 * t363;
t869 = t135 * qJD(1);
t136 = -t271 * t407 + t356 * t364;
t868 = t136 * qJD(1);
t819 = t271 * t1024;
t691 = -t1021 * t356 * t633 + t643 * t819;
t138 = -t979 / 0.2e1 + t691;
t867 = t138 * qJD(1);
t820 = t269 * t1024;
t141 = t820 - t924 / 0.2e1 + (-t980 / 0.2e1 + t1032) * t640;
t866 = t141 * qJD(1);
t710 = t1061 * t601 + t819;
t761 = t524 / 0.2e1 - t940 / 0.2e1;
t142 = t710 - t761;
t865 = t142 * qJD(1);
t195 = 0.2e1 * t1061;
t864 = t195 * qJD(1);
t810 = t951 / 0.2e1;
t651 = t1032 * t708 + t454 * t810;
t198 = t962 / 0.2e1 + t651;
t863 = t198 * qJD(1);
t753 = t454 * t644 + t456 * t641;
t224 = t753 * t534;
t862 = t224 * qJD(1);
t274 = -t1049 * t641 + t454 * t536;
t861 = t274 * qJD(1);
t775 = t1008 * t1040;
t779 = t842 / 0.2e1;
t801 = -t934 / 0.2e1;
t276 = t801 + t961 / 0.2e1 + (t779 + t775) * t637;
t860 = t276 * qJD(1);
t777 = -t834 / 0.2e1;
t648 = -t708 * t534 / 0.2e1 + t454 * t777;
t692 = -t919 / 0.2e1 - t781 / 0.2e1;
t277 = t648 + t692;
t859 = t277 * qJD(1);
t811 = -t951 / 0.2e1;
t709 = t588 * t536 / 0.2e1 + t456 * t811;
t816 = t569 * t1011;
t301 = t816 + t709;
t858 = t301 * qJD(1);
t376 = -t536 ^ 2 + t1049;
t857 = t376 * qJD(1);
t839 = -t997 / 0.2e1;
t689 = t536 * t839 + t586 * t811;
t392 = t689 + t915;
t856 = t392 * qJD(1);
t776 = t834 / 0.2e1;
t678 = t534 * t839 + t586 * t776;
t394 = (t832 / 0.2e1 - t773 / 0.2e1) * t638 + t678;
t855 = t394 * qJD(1);
t854 = t406 * qJD(6);
t853 = t409 * qJD(1);
t852 = t412 * qJD(1);
t472 = (t914 * t639 + (qJ(2) * t952 - t621) * t636) * t638;
t851 = t472 * qJD(1);
t768 = t789 / 0.2e1;
t529 = t774 / 0.2e1 + (t833 / 0.2e1 + t639 * t768) * t638;
t850 = t529 * qJD(1);
t590 = (t636 ^ 2 + t639 ^ 2) * t638 ^ 2;
t849 = t590 * qJD(1);
t848 = t592 - t591;
t847 = t1007 / 0.2e1;
t846 = t1006 / 0.2e1;
t841 = t67 * t1013;
t840 = t67 * t1009;
t838 = -t996 / 0.2e1;
t837 = t67 / 0.2e1 + t1046;
t831 = t456 * t890;
t830 = t599 * t885;
t829 = t633 * t885;
t828 = t601 * t882;
t826 = t640 * t877;
t825 = t599 * t872;
t824 = t536 * t890;
t823 = t536 * t888;
t822 = qJD(4) * t960;
t821 = t640 * t872;
t585 = t601 * t877;
t818 = t271 * t1021;
t817 = t1053 * t1021;
t814 = t356 * t1020;
t813 = t269 * t1019;
t812 = t271 * t1018;
t808 = -t947 / 0.2e1;
t807 = t943 / 0.2e1;
t806 = t941 / 0.2e1;
t805 = -t941 / 0.2e1;
t804 = -t937 / 0.2e1;
t802 = t520 / 0.2e1;
t800 = -t932 / 0.2e1;
t797 = t925 / 0.2e1;
t796 = -t921 / 0.2e1;
t794 = -t369 / 0.2e1 + t415 / 0.2e1;
t793 = t1058 + t1036;
t791 = qJD(4) * t1008;
t790 = t638 * t998;
t785 = -qJD(3) + t889;
t784 = -qJD(4) - t890;
t783 = -qJD(6) - t895;
t782 = -qJD(6) - t886;
t778 = qJD(3) * t834;
t770 = qJD(1) * t790;
t769 = qJD(2) * t790;
t767 = t170 / 0.2e1 + t67 * t1022;
t764 = t793 * t599;
t763 = 0.2e1 * t640 * t585;
t760 = t814 + t169 / 0.2e1;
t13 = (-t100 + t79) * t356 + (-t80 + t99) * t1053;
t98 = t1053 * t1059;
t759 = t13 * qJD(1) + t98 * qJD(2);
t756 = -t1053 * t626 - t356 * t627;
t752 = -t599 * t627 - t601 * t626;
t751 = t456 * t777;
t665 = t100 * t1033 + t1036 * t99 + t1037 * t79 + t1058 * t80;
t704 = t1015 * t391 + t390 * t838;
t18 = (t637 * t775 + t704) * pkin(4) + t665;
t750 = t18 * qJD(1);
t647 = (-t298 * t1008 / 0.2e1 + t226 * t642 / 0.2e1) * t637 + t133 * t1036 + t134 * t1037 + t562 * t1048 + t563 * t1047;
t679 = t1017 * t569 + t1054 * t391 + t1057 * t390;
t21 = t647 + t679;
t299 = -t1008 * t637 ^ 2 * t642 + t483 * t562 + t563 * t653;
t749 = qJD(1) * t21 + qJD(2) * t299;
t200 = t1059 * t601 + t764;
t748 = qJD(1) * t98 + qJD(3) * t200;
t469 = t848 * t640;
t731 = t640 * t765;
t672 = t1053 * t803 + t731 - t959 / 0.2e1;
t84 = t797 + t672;
t747 = -qJD(1) * t84 + qJD(3) * t469;
t470 = t1052 * t640;
t673 = t1053 * t802 + t731 + t959 / 0.2e1;
t86 = t797 + t673;
t746 = qJD(1) * t86 - qJD(3) * t470;
t471 = t848 * t643;
t657 = (t817 - t765) * t643 + t818;
t88 = t806 + t657;
t745 = -qJD(1) * t88 - qJD(3) * t471;
t528 = t1052 * t643;
t656 = (t817 + t765) * t643 + t818;
t90 = t805 + t656;
t744 = qJD(1) * t90 - qJD(3) * t528;
t743 = t782 * t643;
t104 = 0.2e1 * t1045 * t601 + t765;
t125 = t1053 ^ 2 + t356 ^ 2;
t742 = qJD(1) * t125 + qJD(3) * t104;
t741 = qJD(1) * t104 - qJD(3) * t1052;
t705 = -t1015 * t356 + t1053 * t838;
t163 = (t1041 + t705) * pkin(4);
t702 = t1016 * t599 + t601 * t838;
t478 = (t1012 + t702) * pkin(4);
t740 = qJD(1) * t163 + qJD(3) * t478;
t739 = qJD(1) * t190 + qJD(3) * t520;
t191 = t795 * t640;
t517 = t640 * t599;
t172 = qJD(1) * t191 + qJD(3) * t517;
t194 = 0.2e1 * t1060;
t522 = t643 * t599;
t738 = -qJD(1) * t194 - qJD(3) * t522;
t244 = t454 ^ 2 - t456 ^ 2;
t737 = qJD(1) * t244 - qJD(3) * t753;
t736 = -qJD(1) * t753 + qJD(3) * t623;
t349 = -t955 / 0.2e1 + t442 / 0.2e1;
t735 = qJD(1) * t349 + qJD(3) * t589;
t733 = t636 * t770;
t732 = t639 * t770;
t677 = pkin(3) * t454 / 0.2e1 + t982 / 0.2e1 + pkin(10) * t409 / 0.2e1;
t145 = t677 + t794;
t730 = pkin(3) * t881 - qJD(1) * t145;
t712 = pkin(3) * t1041 - pkin(10) * t412 / 0.2e1;
t147 = -t416 / 0.2e1 + (t344 / 0.2e1 - t370 / 0.2e1) * t641 + t712;
t729 = pkin(3) * t884 - qJD(1) * t147;
t722 = t1019 * t601 + t1020 * t599;
t721 = t810 - t970 / 0.2e1;
t717 = t497 / 0.2e1 + t499 / 0.2e1 - t514 / 0.2e1 + t542 / 0.2e1;
t350 = t1012 * t454 + t1040 * t644;
t716 = -qJD(3) * t350 + t454 * t891;
t715 = qJD(1) * t350 + t641 * t881;
t714 = qJD(4) * t529 + t824;
t706 = t134 * t1015 + t133 * t996 / 0.2e1;
t703 = t1015 * t563 + t562 * t838;
t701 = t525 / 0.2e1 + t722;
t417 = -t643 * t707 + t967;
t646 = t417 * t1045 + t269 * t1054 + t757 * t1024 + (t923 + t938) * t1044 + t550 * t807 + t599 * t841;
t670 = t1010 * t117 + t1018 * t363 + t626 * t806;
t1 = (t39 / 0.2e1 - t1000 / 0.2e1) * t601 + t646 + t670;
t690 = t1021 * t653 + t764;
t655 = t1022 * t463 + t640 * t690;
t177 = t963 / 0.2e1 + t655;
t183 = (-t417 + t967) * t601 + t923 * t599;
t700 = -t1 * qJD(1) + t177 * qJD(2) + t183 * qJD(3);
t660 = t792 * t653;
t155 = (t641 * t776 + t703) * pkin(4) + t660;
t351 = t1006 * t628;
t664 = t100 * t1054 + t79 * t688 / 0.2e1 + t80 * t1056 + t99 * t1057;
t7 = (t1012 * t226 + t1017 * t456 + t706) * pkin(4) + t664;
t699 = -t7 * qJD(1) - t155 * qJD(2) + t351 * qJD(3);
t654 = t1022 * t464 + t643 * t690;
t180 = -t964 / 0.2e1 + t654;
t184 = (-t418 + t966) * t601 - t939 * t599;
t645 = t1053 * t1042 + t271 * t1054 + t758 * t1023 + (-t922 + t939) * t1043 + t550 * t1060 + t599 * t840;
t671 = t1013 * t117 + t1018 * t364 + t626 * t797;
t2 = (t40 / 0.2e1 - t999 / 0.2e1) * t601 + t645 + t671;
t698 = -t2 * qJD(1) + t180 * qJD(2) + t184 * qJD(3);
t658 = t1045 * t688 + t1053 * t1054;
t687 = (-t1016 * t407 + t408 * t838) * pkin(4);
t11 = (t1047 - t99 / 0.2e1) * t601 + (t1046 + t1048) * t599 + t687 + t658;
t697 = t11 * qJD(1) - t200 * qJD(2);
t681 = t1023 * t40 + t1042 * t356 + t1057 * t271;
t14 = t643 * t767 + t681 + t808;
t252 = t804 + t973 / 0.2e1 + t721 * t643;
t332 = -t418 * t599 + t550 * t958;
t696 = qJD(1) * t14 + qJD(2) * t252 - qJD(3) * t332;
t682 = t1024 * t39 + t1044 * t417 + t1056 * t269;
t15 = -t640 * t767 + t682 + t800;
t253 = t796 - t974 / 0.2e1 - t721 * t640;
t331 = t417 * t599 - t520 * t550;
t695 = -qJD(1) * t15 - qJD(2) * t253 + qJD(3) * t331;
t663 = t1021 * t79 + t1023 * t80 + t1053 * t1057 - t1054 * t356;
t31 = t663 + t717;
t361 = t969 / 0.2e1 + t721;
t396 = t550 * t601 - t599 * t688;
t694 = -qJD(1) * t31 - qJD(2) * t361 + qJD(3) * t396;
t106 = (-t944 + t984) * t601;
t130 = t269 ^ 2 - t271 ^ 2;
t693 = qJD(1) * t130 - qJD(3) * t106 + qJD(4) * t124;
t151 = t718 * t601;
t158 = -t944 / 0.2e1 + t984 / 0.2e1;
t686 = qJD(3) * t151 - qJD(4) * t158 + t269 * t899;
t516 = (t631 / 0.2e1 - t633 / 0.2e1) * t601;
t685 = qJD(1) * t158 - qJD(3) * t516 + t826;
t176 = t1053 * t1023 + t980 / 0.2e1;
t684 = qJD(3) * t176 + qJD(6) * t349 + t1053 * t895;
t683 = qJD(1) * t176 + qJD(6) * t589 + t830;
t527 = t622 * t592;
t676 = qJD(1) * t106 + qJD(3) * t527 + t763;
t675 = -qJD(1) * t124 + t1050;
t674 = qJD(1) * t151 + qJD(4) * t516 + t592 * t827;
t27 = t640 * t760 + t643 * t837 + t813;
t291 = t640 * t701 + t643 * t792;
t384 = t793 * t643;
t669 = -qJD(1) * t27 + qJD(2) * t384 - qJD(3) * t291 - t627 * t877;
t29 = t640 * t837 - t643 * t760 + t812;
t293 = t640 * t792 - t643 * t701;
t383 = t793 * t640;
t668 = -qJD(1) * t29 + qJD(2) * t383 - qJD(3) * t293 - t627 * t879;
t587 = t589 * qJD(4);
t571 = -0.2e1 * t601 * t821;
t515 = t529 * qJD(3);
t501 = t517 * qJD(6);
t500 = t516 * qJD(6);
t477 = pkin(4) * t702 + t846;
t404 = -t925 / 0.2e1;
t395 = t1008 * t809 + t768 * t952 + t678;
t393 = -t689 + t915;
t386 = (t1009 - t1010) * t483;
t385 = (t1013 - t1014) * t483;
t362 = -t969 / 0.2e1 + t970 / 0.2e1 + t810;
t348 = t350 * qJD(4);
t343 = t1053 * t877;
t302 = t816 - t709;
t297 = (qJD(1) * t1053 + t885) * t643;
t294 = t550 * t1013 + t938 / 0.2e1 + t923 / 0.2e1 - t722 * t643;
t292 = t550 * t1009 + t922 / 0.2e1 - t939 / 0.2e1 + t722 * t640;
t279 = -t961 / 0.2e1 + t751 + t801 + t637 * t779;
t278 = -t648 + t692;
t255 = -t973 / 0.2e1 + t483 * t815 + t804 + t643 * t810;
t254 = t974 / 0.2e1 + t483 * t803 + t796 + t640 * t811;
t225 = t753 * qJD(4);
t205 = -qJD(3) * t406 + qJD(4) * t349;
t199 = -t962 / 0.2e1 + t651;
t197 = t200 * qJD(4);
t192 = t1014 * t356 + t807;
t187 = t192 * qJD(6);
t186 = t191 * qJD(6);
t179 = t964 / 0.2e1 + t654;
t178 = -t963 / 0.2e1 + t655;
t175 = t176 * qJD(4);
t171 = -t172 - t873;
t164 = pkin(4) * t705 + t847;
t157 = t158 * qJD(6);
t156 = pkin(4) * t703 + t1006 * t777 - t660;
t154 = -t728 + t916;
t149 = t151 * qJD(6);
t148 = t344 * t1011 + t935 / 0.2e1 + t416 / 0.2e1 + t712;
t146 = t677 - t794;
t143 = t710 + t761;
t140 = t356 * t803 + t820 + t924 / 0.2e1 + t968 / 0.2e1;
t137 = t979 / 0.2e1 + t691;
t123 = t124 * qJD(6);
t122 = t720 + t724;
t121 = t719 + t725;
t105 = t106 * qJD(6);
t103 = t104 * qJD(5);
t97 = t98 * qJD(4);
t89 = t806 + t656;
t87 = t805 + t657;
t85 = t404 + t673;
t83 = t404 + t672;
t61 = t659 + t726;
t55 = t680 + t727;
t53 = t978 / 0.2e1 + t649;
t52 = -t977 / 0.2e1 + t650;
t45 = t661 + t667;
t43 = t662 - t666;
t32 = -t663 + t717;
t30 = t626 * t1061 + t812 + t841 - t948 / 0.2e1 + t931 / 0.2e1;
t28 = t640 * t814 + t813 + t840 - t933 / 0.2e1 - t946 / 0.2e1;
t20 = t647 - t679;
t19 = -t665 + (t704 + t751) * pkin(4);
t17 = t67 * t815 + t808 + t930 / 0.2e1 - t681;
t16 = t67 * t803 + t800 - t945 / 0.2e1 - t682;
t12 = t100 * t1024 + t1021 * t99 + t1022 * t80 + t1023 * t79 - t658 + t687;
t8 = pkin(4) * t706 + t226 * t846 + t628 * t847 - t664;
t4 = t1022 * t40 + t815 * t99 - t645 + t671;
t3 = t1022 * t39 + t802 * t99 - t646 + t670;
t35 = [0, 0, 0, -t636 * t769, -t639 * t769, t590 * qJD(2), t472 * qJD(2), -t823, t376 * qJD(3), t534 * t887, t536 * t887, 0, qJD(2) * t399 + qJD(3) * t203, qJD(2) * t400 + qJD(3) * t202 (-qJD(4) * t454 - t534 * t881) * t456, qJD(3) * t224 + qJD(4) * t244, qJD(3) * t275 - t454 * t880, -qJD(3) * t274 - t456 * t880, t823, qJD(2) * t282 + qJD(3) * t47 + qJD(4) * t132, qJD(2) * t283 + qJD(3) * t48 + qJD(4) * t131, qJD(2) * t144 + qJD(3) * t23 + qJD(4) * t13 + qJD(5) * t125, qJD(2) * t33 + qJD(3) * t24 + qJD(4) * t22 + qJD(5) * t34 (qJD(3) * t364 - qJD(6) * t269 - t356 * t877) * t271, qJD(3) * t102 - qJD(4) * t101 + qJD(6) * t130, qJD(3) * t136 + qJD(4) * t109 - t269 * t874, qJD(3) * t135 - qJD(4) * t108 - t271 * t874 (-qJD(3) * t407 + qJD(4) * t1053) * t356, qJD(2) * t126 + qJD(3) * t9 + qJD(4) * t5 - qJD(5) * t107 + qJD(6) * t26, qJD(2) * t127 + qJD(3) * t10 + qJD(4) * t6 + qJD(5) * t110 + qJD(6) * t25; 0, 0, 0, -t733, -t732, t849, t851, 0, 0, 0, 0, 0, qJD(3) * t393 + t894, qJD(3) * t395 + t893, 0, 0, 0, 0, 0, qJD(3) * t199 + qJD(4) * t279 + t897, qJD(3) * t302 + qJD(4) * t278 + t896, qJD(3) * t61 + t903 + t97, t983 + (t390 * t483 + t391 * t653 - t569 * t834) * qJD(2) + t20 * qJD(3) + t19 * qJD(4) + t154 * qJD(5), 0, 0, 0, 0, 0, qJD(3) * t43 + qJD(4) * t52 + qJD(6) * t122 + t907, qJD(3) * t45 + qJD(4) * t53 + qJD(6) * t121 + t906; 0, 0, 0, 0, 0, 0, 0, -t824, t857, t785 * t534, t785 * t536, 0, qJD(2) * t393 - qJD(3) * t371 + t900, qJD(2) * t395 + qJD(3) * t370 + t901, t348 + (-t884 - t891) * t412, -t623 * t888 - t225 + t862, t536 * t884 + t898, t536 * t881 - t861, t714, t972 + t199 * qJD(2) + (-t371 * t644 + t641 * t766) * qJD(3) + t148 * qJD(4), t971 + t302 * qJD(2) + (t371 * t641 + t644 * t766) * qJD(3) + t146 * qJD(4), t988 + t61 * qJD(2) + (-t133 * t601 - t134 * t599 + t407 * t688 + t408 * t550) * qJD(3) + t12 * qJD(4) + t103, t987 + t20 * qJD(2) + (-t133 * t550 + t134 * t688 + t298 * t628) * qJD(3) + t8 * qJD(4) + t32 * qJD(5), t137 * qJD(4) - t149 + (t828 + t899) * t364, t870 + t55 * qJD(4) - t105 + (-t363 * t643 - t979) * t885, t868 + (t364 * t599 - t601 * t925) * qJD(3) + t87 * qJD(4) + t140 * qJD(6), t869 + (-t363 * t599 + t407 * t520) * qJD(3) + t83 * qJD(4) + t143 * qJD(6), -t854 + t175 - (t886 + t895) * t407, t1001 + t43 * qJD(2) + (t117 * t520 + t363 * t550 + t407 * t417 + t49 * t599) * qJD(3) + t3 * qJD(4) + t85 * qJD(5) + t17 * qJD(6), t990 + t45 * qJD(2) + (t117 * t958 + t364 * t550 + t407 * t418 - t50 * t599) * qJD(3) + t4 * qJD(4) + t89 * qJD(5) + t16 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t716, t737, t784 * t454, t784 * t456, t515, qJD(2) * t279 + qJD(3) * t148 - qJD(4) * t182 + t904, qJD(2) * t278 + qJD(3) * t146 + qJD(4) * t181 + t905, t12 * qJD(3) + (-t1053 * t635 + t356 * t996) * t1004 + t759, t989 + t19 * qJD(2) + t8 * qJD(3) + (t100 * t635 - t99 * t996) * t1004 + t164 * qJD(5), t137 * qJD(3) + t157 - (t879 + t899) * t1055, t55 * qJD(3) - t1051 * t356 + t123 - t871, qJD(3) * t87 + t1053 * t879 + t1063 + t911, qJD(3) * t83 + t187 + t343 - t912, t684, t1003 + t52 * qJD(2) + t3 * qJD(3) + (t640 * t756 - t999) * qJD(4) + t30 * qJD(6), t1002 + t53 * qJD(2) + t4 * qJD(3) + (t643 * t756 + t1000) * qJD(4) + t28 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t742, qJD(2) * t154 + qJD(3) * t32 + qJD(4) * t164 + t993, 0, 0, 0, 0, 0, qJD(3) * t85 + t187 - t913, qJD(3) * t89 + t1063 + t910; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t686, t693, t140 * qJD(3) + qJD(4) * t1062 + t269 * t783, t143 * qJD(3) + t192 * qJD(4) + t271 * t783, t205, qJD(2) * t122 + qJD(3) * t17 + qJD(4) * t30 + qJD(5) * t192 - qJD(6) * t40 + t994, qJD(2) * t121 + qJD(3) * t16 + qJD(4) * t28 + qJD(5) * t1062 + qJD(6) * t39 + t995; 0, 0, 0, t733, t732, -t849, -t851, 0, 0, 0, 0, 0, -qJD(3) * t392 - t894, qJD(3) * t394 - t893, 0, 0, 0, 0, 0, qJD(3) * t198 - qJD(4) * t276 - t897, -qJD(3) * t301 - qJD(4) * t277 - t896, qJD(3) * t62 - t903 + t97, qJD(3) * t21 - qJD(4) * t18 - qJD(5) * t153 - t983, 0, 0, 0, 0, 0, qJD(3) * t44 + qJD(4) * t51 - qJD(6) * t119 - t907, qJD(3) * t46 + qJD(4) * t54 - qJD(6) * t120 - t906; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t299, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t637 * t883 - t856, -t778 + t855, 0, 0, 0, 0, 0, t863 + (-t641 * t791 - t642 * t881) * t637, -t858 + (t641 * t883 - t644 * t791) * t637, t957 + (t562 * t601 - t563 * t599) * qJD(3) + t197 (t550 * t562 + t563 * t688 + t628 * t951) * qJD(3) + t156 * qJD(4) + t362 * qJD(5) + t749, 0, 0, 0, 0, 0, t976 + (t495 * t599 + t520 * t562) * qJD(3) + t178 * qJD(4) + t255 * qJD(6), t975 + (-t496 * t599 + t562 * t958) * qJD(3) + t179 * qJD(4) + t254 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t588 - t641 * t778 - t860, qJD(4) * t708 - t644 * t778 - t859, t748, t156 * qJD(3) + (-t483 * t635 - t653 * t996) * t1004 - t750, 0, 0, 0, 0, 0, qJD(3) * t178 + qJD(6) * t385 - t653 * t877 + t992, qJD(3) * t179 + qJD(6) * t386 + t653 * t879 + t991; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t362 - t902, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t255 + qJD(4) * t385 - qJD(6) * t464 - t909, qJD(3) * t254 + qJD(4) * t386 + qJD(6) * t463 - t908; 0, 0, 0, 0, 0, 0, 0, t824, -t857, -t534 * t889, -t536 * t889, 0, qJD(2) * t392 - t900, -qJD(2) * t394 - t901, t644 * t831 + t348, -t225 - t862, qJD(4) * t412 - t898, -qJD(4) * t409 + t861, -t714, -qJD(2) * t198 + qJD(4) * t147 - t972, qJD(2) * t301 + qJD(4) * t145 - t971, -qJD(2) * t62 - qJD(4) * t11 + t103 - t988, -qJD(2) * t21 - qJD(4) * t7 - qJD(5) * t31 - t987, qJD(4) * t138 - t364 * t899 - t149, qJD(4) * t56 - t105 - t870, qJD(4) * t88 + qJD(6) * t141 - t868, qJD(4) * t84 + qJD(6) * t142 - t869, t407 * t895 + t175 + t854, -qJD(2) * t44 - qJD(4) * t1 + qJD(5) * t86 - qJD(6) * t14 - t1001, -qJD(2) * t46 - qJD(4) * t2 + qJD(5) * t90 - qJD(6) * t15 - t990; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, -t855, 0, 0, 0, 0, 0, -t863, t858, t197 - t957, -qJD(4) * t155 - qJD(5) * t361 - t749, 0, 0, 0, 0, 0, qJD(4) * t177 - qJD(6) * t252 - t976, qJD(4) * t180 - qJD(6) * t253 - t975; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t641 * t876, t623 * qJD(4), 0, 0, 0, -pkin(3) * t878, -pkin(3) * t876, -qJD(5) * t1052, qJD(4) * t351 + qJD(5) * t396, -t592 * t821 - t633 * t822, -qJD(6) * t527 + t599 * t763, qJD(4) * t471 - t873 * t960, -qJD(4) * t469 - t601 * t825, t822, qJD(4) * t183 - qJD(5) * t470 + qJD(6) * t332, qJD(4) * t184 - qJD(5) * t528 + qJD(6) * t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t715, t736, t852 + t876, -t853 - t878, -t850, -pkin(10) * t876 - t729, pkin(10) * t878 - t730 (t599 * t996 - t601 * t635) * t1004 - t697 (-t550 * t635 - t688 * t996) * t1004 + t477 * qJD(5) + t699, t867 - t500 + (-t826 - t829) * t599, t1050 * t599 + t571 + t965, t601 * t879 - t745, t585 - t747, t683 (t640 * t752 - t966) * qJD(4) + t294 * qJD(6) + t700 (t643 * t752 + t967) * qJD(4) + t292 * qJD(6) + t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t741, qJD(4) * t477 + t694, 0, 0, 0, 0, 0, t746, t744; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t674, -t676, t520 * t782 + t866, t601 * t743 + t865, t587 + t892, qJD(4) * t294 - qJD(6) * t418 - t696, qJD(4) * t292 + qJD(6) * t417 + t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, -t737, -qJD(3) * t412 + t454 * t890, qJD(3) * t409 + t831, t515, qJD(2) * t276 - qJD(3) * t147 - t904, qJD(2) * t277 - qJD(3) * t145 - t905, qJD(3) * t11 - t759, qJD(2) * t18 + qJD(3) * t7 + qJD(5) * t163 - t989, -qJD(3) * t138 + t1055 * t899 + t157, -qJD(3) * t56 + t123 + t871, -qJD(3) * t88 + qJD(6) * t194 - t911, -qJD(3) * t84 - t186 + t912, -t684, -qJD(2) * t51 + qJD(3) * t1 + qJD(6) * t29 - t1053 * t875 - t1003, -qJD(2) * t54 + qJD(3) * t2 + qJD(5) * t190 + qJD(6) * t27 - t1002; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t860, t859, -t748, qJD(3) * t155 + t750, 0, 0, 0, 0, 0, -qJD(3) * t177 - qJD(6) * t383 - t992, -qJD(3) * t180 - qJD(6) * t384 - t991; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t715, -t736, -t852, t853, t850, t729, t730, t697, qJD(5) * t478 - t699, t599 * t829 - t500 - t867, -0.2e1 * t599 * t780 + t571 - t965, qJD(6) * t522 + t745, -t501 + t747, -t683, qJD(6) * t293 - t601 * t875 - t700, qJD(5) * t520 + qJD(6) * t291 - t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t821, t622 * qJD(6), 0, 0, 0, t627 * t873, t627 * t872; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t740, 0, 0, 0, 0, 0, -t297, t739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t685, -t675, -t738 + t872, t171, -t735, -t626 * t872 - t668, t626 * t873 - t669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t742, qJD(2) * t153 + qJD(3) * t31 - qJD(4) * t163 - t993, 0, 0, 0, 0, 0, -qJD(3) * t86 - t186 + t343 + t913, -qJD(3) * t90 - qJD(4) * t190 + qJD(6) * t195 - t910; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t361 + t902, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t741, -qJD(4) * t478 - t694, 0, 0, 0, 0, 0, -t501 + t585 - t746, -qJD(4) * t520 - t744 - t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t740, 0, 0, 0, 0, 0, t297, -t739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t743 + t864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t686, -t693, -qJD(3) * t141 - qJD(4) * t194 + t269 * t895, -qJD(3) * t142 + qJD(4) * t191 + t271 * t895, t205, qJD(2) * t119 + qJD(3) * t14 - qJD(4) * t29 + qJD(5) * t191 - t994, qJD(2) * t120 + qJD(3) * t15 - qJD(4) * t27 - qJD(5) * t195 - t995; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t252 + qJD(4) * t383 + t909, qJD(3) * t253 + qJD(4) * t384 + t908; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t674, t676, -qJD(4) * t522 + t640 * t830 - t866, qJD(4) * t517 + t599 * t828 - t865, t587 - t892, -qJD(4) * t293 + qJD(5) * t517 + t696, -qJD(4) * t291 + t599 * t875 - t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t685, t675, t738, t172, t735, t668, t669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t599 * t882 - t864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t35;
