% Calculate minimal parameter regressor of coriolis matrix for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRRPRP11_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:07
% EndTime: 2019-03-09 17:49:56
% DurationCPUTime: 28.76s
% Computational Cost: add. (17635->990), mult. (40746->1336), div. (0->0), fcn. (42533->8), ass. (0->706)
t1003 = pkin(3) + pkin(10);
t596 = cos(qJ(5));
t524 = (-qJ(6) - t1003) * t596;
t980 = t524 / 0.2e1;
t594 = sin(qJ(3));
t597 = cos(qJ(3));
t592 = sin(pkin(6));
t595 = sin(qJ(2));
t894 = t592 * t595;
t950 = cos(pkin(6));
t485 = t594 * t950 + t597 * t894;
t984 = t485 / 0.2e1;
t598 = cos(qJ(2));
t893 = t592 * t598;
t795 = t596 * t893;
t517 = t594 * t795;
t593 = sin(qJ(5));
t438 = t593 * t894 - t517;
t877 = t595 * t596;
t878 = t594 * t598;
t439 = (t593 * t878 + t877) * t592;
t566 = pkin(3) * t894;
t791 = pkin(1) * t950;
t493 = pkin(8) * t894 - t598 * t791;
t471 = t594 * t493;
t494 = (pkin(2) * t595 - pkin(9) * t598) * t592;
t871 = t597 * t494;
t733 = t471 + t871;
t304 = -t566 - t733;
t870 = t597 * t598;
t235 = (pkin(4) * t870 - pkin(10) * t595) * t592 + t304;
t229 = t596 * t235;
t949 = qJ(4) * t597;
t702 = pkin(10) * t594 - t949;
t495 = pkin(8) * t893 + t595 * t791;
t797 = t592 * t878;
t793 = pkin(3) * t797 + t495;
t303 = t702 * t893 + t793;
t891 = t593 * t303;
t148 = t229 - t891;
t796 = t592 * t870;
t133 = pkin(5) * t796 - qJ(6) * t439 + t148;
t1001 = t133 / 0.2e1;
t228 = t593 * t235;
t287 = t596 * t303;
t149 = t287 + t228;
t143 = -qJ(6) * t438 + t149;
t972 = t593 / 0.2e1;
t867 = t1001 * t596 + t143 * t972;
t884 = t593 * t1003;
t523 = -qJ(6) * t593 - t884;
t982 = -t523 / 0.2e1;
t1027 = -t438 * t982 + t439 * t980 + t867;
t906 = t485 * t597;
t761 = t906 / 0.2e1;
t572 = t950 * t597;
t483 = t594 * t894 - t572;
t880 = t594 * t483;
t329 = -t880 / 0.2e1 + t761;
t850 = qJD(1) * t485;
t357 = t483 * t850;
t1018 = -qJD(2) * t329 + t357;
t806 = t597 * qJD(2);
t561 = t594 * t806;
t1014 = qJD(1) * t329 + t561;
t1026 = t1003 * t483;
t883 = t594 * qJ(4);
t956 = t597 * pkin(3);
t703 = -t883 - t956;
t526 = -pkin(2) + t703;
t453 = -pkin(2) * t950 + t493;
t907 = t485 * qJ(4);
t613 = t453 - t907;
t960 = t483 * pkin(3);
t256 = t613 + t960;
t929 = t256 * t594;
t1017 = t526 * t984 + t929 / 0.2e1;
t1025 = t1003 * t597 + t883;
t1024 = pkin(4) + pkin(9);
t199 = t613 + t1026;
t963 = pkin(9) * t595;
t707 = -pkin(2) * t598 - t963;
t455 = (-pkin(1) + t707) * t592;
t428 = t597 * t455;
t454 = pkin(9) * t950 + t495;
t882 = t594 * t454;
t297 = -t428 + t882;
t230 = -t485 * pkin(4) - t297;
t567 = pkin(3) * t893;
t615 = pkin(10) * t893 - t230 + t567;
t122 = t199 * t596 + t593 * t615;
t274 = t593 * t893;
t908 = t483 * t596;
t386 = t274 + t908;
t101 = qJ(6) * t386 + t122;
t941 = t101 * t596;
t121 = t199 * t593 - t596 * t615;
t909 = t483 * t593;
t387 = t795 - t909;
t100 = qJ(6) * t387 - t121;
t82 = t485 * pkin(5) + t100;
t952 = t82 * t593;
t1023 = t485 * (t941 - t952);
t497 = -pkin(2) - t1025;
t539 = t1024 * t594;
t513 = t596 * t539;
t378 = t497 * t593 - t513;
t885 = t593 * t597;
t335 = qJ(6) * t885 - t378;
t958 = t594 * pkin(5);
t317 = t335 + t958;
t923 = t317 * t593;
t768 = -t923 / 0.2e1;
t887 = t593 * t539;
t379 = t497 * t596 + t887;
t873 = t596 * t597;
t336 = -qJ(6) * t873 + t379;
t920 = t336 * t596;
t1022 = (t920 / 0.2e1 + t768) * t485;
t163 = (-t920 + t923) * t594;
t999 = -t317 / 0.2e1;
t741 = t335 / 0.2e1 + t999;
t1021 = t741 * t101;
t312 = t329 * qJD(3);
t844 = qJD(2) * t598;
t778 = t592 * t844;
t531 = t597 * t778;
t1020 = t531 * t594 + t312;
t1019 = -t317 + t335;
t774 = t941 / 0.2e1;
t1016 = t774 - t941 / 0.2e1;
t757 = t894 / 0.2e1;
t476 = t594 * t757 - t572 / 0.2e1;
t1015 = -qJD(5) * t476 - t1018;
t805 = t597 * qJD(5);
t1013 = t805 / 0.2e1 + t1014;
t638 = (t958 / 0.2e1 - t741) * t593;
t1012 = qJD(2) * t638;
t587 = t593 ^ 2;
t589 = t596 ^ 2;
t556 = t587 - t589;
t731 = qJD(3) * t556;
t1011 = qJD(5) * t638;
t1004 = t485 ^ 2;
t1010 = t1004 * t593;
t1009 = t1004 * t596;
t872 = t597 * t454;
t881 = t594 * t455;
t298 = t872 + t881;
t959 = t483 * pkin(4);
t231 = t298 - t959;
t554 = qJ(4) * t893;
t203 = t231 - t554;
t755 = t893 / 0.2e1;
t583 = t597 * pkin(9);
t584 = t597 * pkin(4);
t540 = t583 + t584;
t977 = -t540 / 0.2e1;
t759 = t485 * t977;
t971 = -t594 / 0.2e1;
t1008 = -t203 * t971 + t597 * (-t1003 * t755 - t230 / 0.2e1) - t759;
t834 = qJD(4) * t597;
t1007 = qJD(3) * t1025 - t834;
t259 = t297 + t567;
t924 = t259 * t597;
t258 = t298 - t554;
t925 = t258 * t594;
t968 = -t597 / 0.2e1;
t970 = t594 / 0.2e1;
t1006 = t297 * t968 + t298 * t970 + t924 / 0.2e1 - t925 / 0.2e1;
t1005 = t387 ^ 2;
t218 = t596 * t231;
t910 = t483 * qJ(4);
t260 = t1003 * t485 + t910;
t115 = -t483 * pkin(5) + t218 + (-qJ(6) * t485 - t260) * t593;
t1002 = t115 / 0.2e1;
t1000 = t229 / 0.2e1;
t585 = pkin(3) * t594;
t501 = t585 + t702;
t515 = t540 * t596;
t955 = t597 * pkin(5);
t325 = t955 + t515 + (-qJ(6) * t594 - t501) * t593;
t998 = t325 / 0.2e1;
t996 = -t336 / 0.2e1;
t995 = t336 / 0.2e1;
t994 = -t386 / 0.2e1;
t993 = t386 / 0.2e1;
t992 = -t387 / 0.2e1;
t991 = t387 / 0.2e1;
t990 = t438 / 0.2e1;
t989 = t439 / 0.2e1;
t988 = -t471 / 0.2e1;
t482 = pkin(5) * t873 + t540;
t987 = t482 / 0.2e1;
t986 = -t483 / 0.2e1;
t985 = -t485 / 0.2e1;
t518 = t597 * t274;
t983 = -t518 / 0.2e1;
t981 = t523 / 0.2e1;
t979 = -t526 / 0.2e1;
t978 = -t539 / 0.2e1;
t976 = -t566 / 0.2e1;
t573 = pkin(5) * t593 + qJ(4);
t975 = t573 / 0.2e1;
t974 = t589 / 0.2e1;
t973 = -t593 / 0.2e1;
t969 = -t596 / 0.2e1;
t967 = t597 / 0.2e1;
t966 = t1003 / 0.2e1;
t965 = pkin(5) * t387;
t964 = pkin(5) * t439;
t962 = t386 * pkin(5);
t961 = t438 * pkin(5);
t957 = t596 * pkin(5);
t954 = pkin(5) * qJD(5);
t953 = pkin(5) * qJD(6);
t951 = t100 - t82;
t156 = t203 - t962;
t13 = t101 * t951 - t156 * t965;
t948 = qJD(1) * t13;
t15 = t951 * t386;
t947 = qJD(1) * t15;
t31 = -t156 * t893 - t1023;
t946 = qJD(1) * t31;
t38 = t101 * t386 + t387 * t82;
t945 = qJD(1) * t38;
t51 = t121 * t485 + t203 * t386;
t944 = qJD(1) * t51;
t52 = -t122 * t485 - t203 * t387;
t943 = qJD(1) * t52;
t942 = t100 * t593;
t940 = t115 * t596;
t217 = t593 * t231;
t252 = t596 * t260;
t865 = t252 + t217;
t876 = t596 * qJ(6);
t134 = t485 * t876 + t865;
t938 = t134 * t593;
t792 = -pkin(4) - t957;
t177 = t485 * t792 - t297;
t14 = t101 * t134 + t115 * t82 + t156 * t177;
t937 = t14 * qJD(1);
t553 = qJ(4) * t894;
t472 = t594 * t494;
t473 = t597 * t493;
t864 = t473 - t472;
t302 = -t553 + t864;
t257 = -pkin(4) * t797 - t302;
t189 = t257 + t961;
t16 = t101 * t143 + t133 * t82 + t156 * t189;
t935 = t16 * qJD(1);
t934 = t203 * t593;
t933 = t203 * t596;
t21 = t115 * t387 + t134 * t386 + t1023;
t932 = t21 * qJD(1);
t22 = -t101 * t438 + t133 * t387 + t143 * t386 - t439 * t82;
t931 = t22 * qJD(1);
t930 = t256 * t485;
t928 = t256 * t597;
t927 = t257 * t593;
t926 = t257 * t596;
t892 = t593 * t260;
t734 = t218 - t892;
t32 = t121 * t483 - t230 * t386 + (t734 - t933) * t485;
t922 = t32 * qJD(1);
t33 = t122 * t483 - t230 * t387 + (-t865 + t934) * t485;
t921 = t33 * qJD(1);
t36 = -t121 * t796 + t148 * t485 + t203 * t438 - t257 * t386;
t919 = t36 * qJD(1);
t37 = -t122 * t796 - t149 * t485 + t203 * t439 - t257 * t387;
t918 = t37 * qJD(1);
t917 = t386 * t594;
t916 = t387 * t593;
t915 = t387 * t594;
t914 = t439 * t596;
t349 = pkin(3) * t485 + t910;
t45 = t256 * t349 - t258 * t297 + t259 * t298;
t913 = t45 * qJD(1);
t912 = t453 * t597;
t352 = -t554 * t597 + t793;
t46 = t256 * t352 - t258 * t302 + t259 * t304;
t911 = t46 * qJD(1);
t328 = t485 * t593;
t742 = t297 / 0.2e1 - t259 / 0.2e1;
t743 = t258 / 0.2e1 - t298 / 0.2e1;
t756 = -t893 / 0.2e1;
t49 = (pkin(3) * t756 + t742) * t597 + (qJ(4) * t756 + t743) * t594;
t905 = t49 * qJD(1);
t904 = t523 * t594;
t903 = t523 * t596;
t902 = t524 * t593;
t901 = t524 * t594;
t53 = (-t258 + t298) * t485 + (-t259 + t297) * t483;
t899 = t53 * qJD(1);
t538 = t585 - t949;
t898 = t538 * t483;
t586 = t592 ^ 2;
t591 = t598 ^ 2;
t897 = t586 * t591;
t896 = t586 * t595;
t59 = t302 * t483 + t304 * t485 + (t924 - t925) * t893;
t895 = t59 * qJD(1);
t890 = t593 * t386;
t889 = t593 * t438;
t888 = t593 * t501;
t886 = t593 * t594;
t879 = t594 * t596;
t875 = t596 * t386;
t874 = t596 * t387;
t65 = -t352 * t485 + (t258 * t595 + (t302 - t928) * t598) * t592;
t869 = t65 * qJD(1);
t66 = -t352 * t483 + (t259 * t595 + (-t304 - t929) * t598) * t592;
t868 = t66 * qJD(1);
t477 = t596 * t501;
t514 = t540 * t593;
t863 = t477 + t514;
t555 = t587 + t589;
t588 = t594 ^ 2;
t590 = t597 ^ 2;
t557 = t590 - t588;
t798 = t298 * t893;
t108 = -t349 * t483 - t798 - t930;
t862 = qJD(1) * t108;
t799 = t297 * t893;
t109 = t256 * t483 - t349 * t485 + t799;
t861 = qJD(1) * t109;
t144 = -t258 * t893 - t930;
t860 = qJD(1) * t144;
t161 = (-t875 + t916) * t485;
t859 = qJD(1) * t161;
t173 = -t453 * t483 - t799;
t858 = qJD(1) * t173;
t174 = -t453 * t485 - t798;
t857 = qJD(1) * t174;
t206 = -t386 * t893 - t1010;
t856 = qJD(1) * t206;
t225 = t387 * t483 + t1010;
t855 = qJD(1) * t225;
t226 = -t386 * t483 + t1009;
t854 = qJD(1) * t226;
t371 = -t483 * t894 + t594 * t897;
t853 = qJD(1) * t371;
t372 = -t485 * t894 + t597 * t897;
t852 = qJD(1) * t372;
t851 = qJD(1) * t387;
t849 = qJD(1) * t598;
t848 = qJD(2) * t588;
t847 = qJD(2) * t592;
t846 = qJD(2) * t594;
t845 = qJD(2) * t596;
t843 = qJD(3) * t297;
t842 = qJD(3) * t298;
t841 = qJD(3) * t593;
t840 = qJD(3) * t594;
t839 = qJD(3) * t596;
t838 = qJD(3) * t597;
t837 = qJD(4) * t593;
t836 = qJD(4) * t594;
t835 = qJD(4) * t596;
t833 = qJD(5) * t386;
t832 = qJD(5) * t485;
t831 = qJD(5) * t593;
t830 = qJD(5) * t594;
t829 = qJD(5) * t596;
t828 = qJD(5) * t1003;
t123 = t297 * t894 - t453 * t797 - t483 * t495 + t733 * t893;
t827 = t123 * qJD(1);
t124 = t495 * t485 + (-t298 * t595 + (-t864 + t912) * t598) * t592;
t826 = t124 * qJD(1);
t406 = t873 * t328;
t767 = t917 / 0.2e1;
t132 = -t406 + (-t915 / 0.2e1 + t990) * t596 + (t767 + t989) * t593;
t825 = t132 * qJD(1);
t154 = t386 * t439 + t387 * t438;
t824 = t154 * qJD(1);
t158 = -t874 + t890;
t162 = t158 * t485;
t823 = t162 * qJD(1);
t753 = -t886 / 0.2e1;
t762 = -t906 / 0.2e1;
t630 = t387 * t753 + t587 * t762;
t764 = -t914 / 0.2e1;
t179 = t764 + t630;
t822 = t179 * qJD(1);
t648 = t757 + t761;
t765 = t915 / 0.2e1;
t185 = t765 - t517 / 0.2e1 + t648 * t593;
t821 = t185 * qJD(1);
t752 = t886 / 0.2e1;
t629 = (t877 / 0.2e1 + t598 * t752) * t592;
t748 = -t873 / 0.2e1;
t639 = t485 * t748 + t767;
t187 = t629 - t639;
t820 = t187 * qJD(1);
t190 = t386 * t796 - t438 * t485;
t819 = t190 * qJD(1);
t191 = -t387 * t796 + t439 * t485;
t818 = t191 * qJD(1);
t430 = t485 * t886;
t301 = t430 / 0.2e1 + t485 * t752;
t251 = t597 * t795 - t301;
t817 = t251 * qJD(1);
t255 = t387 * t893 + t1009;
t816 = t255 * qJD(1);
t237 = -t483 * t597 - t485 * t594;
t269 = t237 * t893;
t815 = t269 * qJD(1);
t431 = t485 * t879;
t356 = -t518 - t431;
t812 = t356 * qJD(1);
t395 = pkin(1) * t896 + t495 * t950;
t811 = t395 * qJD(1);
t396 = pkin(1) * t586 * t598 - t493 * t950;
t810 = t396 * qJD(1);
t474 = t483 * qJD(3);
t809 = t483 * qJD(4);
t502 = (-t595 ^ 2 + t591) * t586;
t808 = t502 * qJD(1);
t807 = t592 * qJD(3);
t804 = pkin(5) * t885;
t803 = pkin(5) * t831;
t802 = pkin(9) * t840;
t801 = -t965 / 0.2e1;
t800 = t957 / 0.2e1;
t75 = -t952 / 0.2e1;
t794 = -t100 / 0.2e1 + t82 / 0.2e1;
t790 = t593 * t850;
t789 = t592 * t849;
t788 = t593 * t806;
t787 = t596 * t806;
t786 = t598 * t807;
t785 = t593 * t839;
t784 = t593 * t838;
t783 = qJD(4) * t893;
t782 = t593 * t805;
t781 = t596 * t805;
t780 = t586 * t849;
t779 = t595 * t847;
t562 = t594 * t838;
t777 = t485 * t836;
t564 = t596 * t838;
t776 = t593 * t829;
t563 = t594 * t845;
t766 = t916 / 0.2e1;
t763 = t328 / 0.2e1;
t758 = -t894 / 0.2e1;
t754 = t890 / 0.2e1;
t751 = -t885 / 0.2e1;
t750 = t885 / 0.2e1;
t749 = -t879 / 0.2e1;
t747 = t873 / 0.2e1;
t746 = t485 * t966;
t745 = -t217 / 0.2e1 - t252 / 0.2e1;
t744 = -t228 / 0.2e1 - t287 / 0.2e1;
t740 = t472 / 0.2e1 - t473 / 0.2e1;
t739 = t583 / 0.2e1 + t584 / 0.2e1;
t738 = t974 + t587 / 0.2e1;
t737 = t952 / 0.2e1 - t942 / 0.2e1 + t1016;
t736 = t75 + t942 / 0.2e1 + t1016;
t735 = t950 * qJD(1);
t732 = t515 - t888;
t441 = t880 / 0.2e1;
t292 = t441 - t648;
t730 = qJD(1) * t292 - t561;
t560 = t593 * t846;
t729 = qJD(1) * t328 + t560;
t727 = pkin(5) * t751;
t726 = -pkin(5) * t328 / 0.2e1;
t725 = pkin(5) * t763;
t724 = pkin(9) * t755;
t723 = qJD(5) + t850;
t722 = qJD(5) + t846;
t721 = t483 * t789;
t720 = t593 * t564;
t719 = t844 * t896;
t718 = t595 * t780;
t530 = t597 * t789;
t717 = t593 * t787;
t716 = t594 * t756;
t715 = t597 * t756;
t714 = t596 * t755;
t713 = -pkin(4) / 0.2e1 - t957 / 0.2e1;
t712 = t874 / 0.2e1 + t754;
t711 = t553 + t740;
t709 = t592 * t735;
t708 = t950 * t847;
t706 = 0.2e1 * t717;
t705 = t471 / 0.2e1 + t871 / 0.2e1;
t704 = t596 * t850 + t563;
t534 = -qJD(3) + t789;
t338 = t594 * t876 + t863;
t481 = (-pkin(9) + t792) * t594;
t604 = t101 * t338 / 0.2e1 + t317 * t1002 + t134 * t995 + t156 * t481 / 0.2e1 + t177 * t987 + t82 * t998;
t622 = t133 * t980 + t143 * t981 + t189 * t975;
t1 = -t604 + t622;
t92 = t317 * t325 + t336 * t338 + t481 * t482;
t701 = -qJD(1) * t1 + qJD(2) * t92;
t699 = t597 * t714;
t698 = t907 - t1026;
t103 = -t325 * t885 + t338 * t873 + t163;
t601 = t82 * t753 + t594 * t774 + (t115 * t972 + t134 * t969) * t597 + t1022 + t325 * t991 + t338 * t993;
t4 = t601 + t1027;
t696 = t4 * qJD(1) - t103 * qJD(2);
t104 = t1019 * t336 - t482 * t804;
t5 = t794 * t336 - t1021 + (t156 * t750 + t387 * t987 + t1001) * pkin(5);
t695 = -qJD(1) * t5 + qJD(2) * t104;
t11 = t726 + t736;
t694 = qJD(1) * t11 - t1012;
t693 = -t302 * t597 + t304 * t594;
t691 = t735 + qJD(2);
t610 = t386 * t741 + t794 * t873;
t10 = t964 / 0.2e1 + t610;
t105 = t1019 * t873;
t690 = -qJD(1) * t10 + qJD(2) * t105;
t19 = t725 + t737;
t689 = qJD(1) * t19 + t1012;
t159 = (-t378 - t513) * t597 + (t732 - t515) * t594;
t602 = t121 * t967 + t378 * t986 + t386 * t978 + t732 * t985 + t734 * t971;
t657 = qJ(4) * t990 + t927 / 0.2e1;
t24 = t1008 * t596 + t602 + t657;
t688 = -t24 * qJD(1) + t159 * qJD(2);
t160 = t863 * t594 - t540 * t886 + (t379 - t887) * t597;
t603 = t122 * t967 + t379 * t986 + t387 * t978 + t863 * t984 + t865 * t970;
t656 = qJ(4) * t989 + t926 / 0.2e1;
t23 = -t1008 * t593 + t603 + t656;
t687 = -t23 * qJD(1) - t160 * qJD(2);
t605 = t101 * t749 + t482 * t756 + t752 * t82 - t1022;
t18 = t605 - t867;
t686 = -qJD(1) * t18 - qJD(2) * t163;
t164 = -t317 * t885 + t336 * t873;
t608 = t553 / 0.2e1 + t961 / 0.2e1 + pkin(4) * t716 + t740;
t652 = t317 * t992 + t336 * t994;
t27 = (t774 + t75) * t597 + t608 + t652;
t685 = -qJD(1) * t27 - qJD(2) * t164;
t270 = -t378 * t594 + t540 * t873;
t607 = t121 * t971 + t203 * t747 + t378 * t985 + t386 * t977;
t40 = t607 + t744;
t684 = -qJD(1) * t40 - qJD(2) * t270;
t271 = -t379 * t594 - t540 * t885;
t623 = t122 * t970 + t379 * t984 + t540 * t991;
t39 = t1000 + (-t303 / 0.2e1 + t203 * t967) * t593 + t623;
t683 = -qJD(1) * t39 + qJD(2) * t271;
t393 = t526 * t597 + t538 * t594;
t621 = t928 / 0.2e1 + t483 * t979 + t538 * t984;
t77 = (t349 / 0.2e1 + pkin(9) * t756) * t594 + t621 + t711;
t682 = qJD(1) * t77 + qJD(2) * t393;
t394 = -t526 * t594 + t538 * t597;
t680 = t724 - t494 / 0.2e1;
t79 = -t566 + t988 + t898 / 0.2e1 + (-t349 / 0.2e1 + t680) * t597 + t1017;
t681 = qJD(1) * t79 - qJD(2) * t394;
t669 = t738 * t906;
t113 = (-t917 / 0.2e1 + t989) * t596 + (t765 + t990) * t593 + t669;
t475 = t555 * t597 * t594;
t679 = qJD(1) * t113 + qJD(2) * t475;
t169 = t983 - t431 + (-t908 / 0.2e1 + t994) * t597;
t512 = t557 * t596;
t678 = -qJD(1) * t169 - qJD(2) * t512;
t170 = -t430 + (t714 - t909 / 0.2e1 + t991) * t597;
t509 = t557 * t593;
t677 = -qJD(1) * t170 - qJD(2) * t509;
t353 = t387 * t885;
t207 = t386 * t873 - t353;
t510 = t555 * t590;
t676 = qJD(1) * t207 - qJD(2) * t510;
t272 = t483 ^ 2 - t1004;
t675 = qJD(1) * t272 + qJD(2) * t237;
t674 = qJD(1) * t237 + qJD(2) * t557;
t195 = t755 + t712;
t533 = -0.1e1 / 0.2e1 - t738;
t673 = qJD(1) * t195 + qJD(3) * t533;
t227 = t874 + t890;
t672 = qJD(1) * t227 - qJD(3) * t555;
t671 = qJD(1) * t274 - t841;
t670 = -t830 - t848;
t667 = -t474 + t721;
t666 = t531 - t721;
t665 = pkin(2) * t985 + t453 * t970;
t664 = -t302 * qJ(4) / 0.2e1 - t304 * pkin(3) / 0.2e1;
t663 = t476 * qJD(1) - t806 / 0.2e1;
t662 = t592 * t707;
t612 = pkin(2) * t483 / 0.2e1 + t912 / 0.2e1 + pkin(9) * t716;
t165 = t612 + t740;
t661 = pkin(2) * t806 - qJD(1) * t165;
t626 = t597 * t680 + t988;
t167 = t626 + t665;
t660 = pkin(2) * t846 - qJD(1) * t167;
t397 = t485 * t530;
t659 = qJD(5) * t715 - t397;
t658 = qJ(4) * t992 + t933 / 0.2e1;
t655 = -t940 / 0.2e1 - t938 / 0.2e1;
t653 = -t256 * t538 / 0.2e1 + t349 * t979;
t649 = -t903 / 0.2e1 + t902 / 0.2e1;
t647 = t754 - t874 / 0.2e1;
t646 = t594 * t966 - t949 / 0.2e1;
t43 = (t594 * t743 + t597 * t742) * pkin(9) + t653 + t664;
t645 = -t526 * t538 * qJD(2) + t43 * qJD(1);
t644 = t592 * (-t526 * t598 + t963);
t135 = t976 + t626 + t1017;
t643 = qJD(1) * t135 + t526 * t846;
t642 = -t788 + t839;
t641 = (t839 - t851) * t485;
t640 = -t596 * t723 - t563;
t637 = pkin(9) * t715 - t1017 - t705;
t253 = t573 * t957;
t47 = -t741 * t523 + (t482 * t969 + t573 * t750 + t998) * pkin(5);
t7 = t794 * t523 + (t156 * t969 + t387 * t975 + t1002) * pkin(5);
t636 = -qJD(1) * t7 - qJD(2) * t47 + qJD(3) * t253;
t635 = t646 * t596;
t606 = (t524 * t967 + t996) * t593 + (t523 * t968 + t999) * t596;
t137 = (pkin(9) / 0.2e1 - t713) * t594 + t606;
t609 = t101 * t973 + t386 * t981 + t387 * t980 + t82 * t969;
t611 = t713 * t485 + t428 / 0.2e1 - t882 / 0.2e1;
t29 = -t609 + t611;
t384 = t523 * t593 + t524 * t596;
t634 = -qJD(1) * t29 + qJD(2) * t137 - qJD(3) * t384;
t142 = (t901 / 0.2e1 - t338 / 0.2e1) * t593 + (-t904 / 0.2e1 + t955 / 0.2e1 - t325 / 0.2e1) * t596 + t739;
t600 = t649 * t485 - t554 / 0.2e1 - t962 / 0.2e1 - t959 / 0.2e1 + t881 / 0.2e1 + t872 / 0.2e1 + t573 * t756;
t35 = t600 + t655;
t633 = qJD(1) * t35 + qJD(2) * t142 + qJD(3) * t573;
t153 = t353 / 0.2e1 + (t766 + t875) * t597;
t383 = t386 ^ 2;
t175 = t383 - t1005;
t632 = qJD(1) * t175 - qJD(2) * t153 - qJD(3) * t158;
t238 = t383 + t1005;
t631 = qJD(1) * t238 - qJD(2) * t207 - qJD(3) * t227;
t350 = (t501 / 0.2e1 + t646) * t593;
t54 = -t218 / 0.2e1 + (t746 + t260 / 0.2e1) * t593 + t658;
t628 = -qJ(4) * t839 - qJD(1) * t54 - qJD(2) * t350;
t351 = t477 / 0.2e1 + t635;
t614 = qJ(4) * t993 - t934 / 0.2e1 + t596 * t746;
t56 = t614 - t745;
t627 = qJ(4) * t841 - qJD(1) * t56 - qJD(2) * t351;
t202 = t647 * t597;
t221 = t875 / 0.2e1 + t766;
t625 = qJD(2) * t202 - qJD(3) * t221 + t386 * t851;
t499 = (t974 - t587 / 0.2e1) * t597;
t624 = -qJD(1) * t221 + qJD(2) * t499 + t785;
t620 = qJD(3) * t703 + t834;
t511 = t556 * t590;
t619 = -qJD(1) * t153 - qJD(2) * t511 + 0.2e1 * t720;
t618 = -qJD(1) * t158 + t706 + t731;
t617 = t590 * t593 * t845 - qJD(1) * t202 - qJD(3) * t499;
t415 = t1004 + t897;
t616 = qJD(1) * t415 + t485 * t846 - t786;
t579 = pkin(9) * t838;
t574 = t838 / 0.2e1;
t541 = qJD(2) * t757;
t532 = 0.1e1 / 0.2e1 - t738;
t525 = 0.2e1 * t597 * t776;
t503 = t534 * qJ(4);
t488 = -t530 + t838;
t487 = t534 * t596;
t480 = t499 * qJD(5);
t456 = (t780 - t807 / 0.2e1) * t595;
t398 = -t594 * t850 - t848;
t377 = t531 / 0.2e1 - t476 * qJD(3);
t299 = (-t642 + t851) * pkin(5);
t291 = t441 + t762 + t757;
t289 = -t514 - t477 / 0.2e1 + t635;
t288 = t515 - t888 / 0.2e1 + t646 * t593;
t236 = (t531 - t474) * t485;
t234 = t237 * qJD(3);
t219 = t227 * qJD(6);
t215 = t221 * qJD(5);
t201 = t207 * qJD(6);
t196 = t756 + t712;
t194 = t202 * qJD(5);
t188 = t629 + t639;
t186 = t485 * t750 + t765 + t593 * t758 + t517 / 0.2e1;
t178 = t914 / 0.2e1 + t630;
t172 = t386 * t967 + t483 * t747 + t431 + t983;
t171 = t387 * t968 + t483 * t750 + t430 + t699;
t168 = t597 * t724 + t665 + t705;
t166 = t612 - t740;
t157 = t158 * qJD(5);
t151 = t153 * qJD(5);
t141 = pkin(5) * t747 + t338 * t972 + t594 * t649 + t596 * t998 + t739;
t138 = pkin(5) * t749 + t1024 * t971 + t606;
t136 = t976 + t637;
t131 = t438 * t969 + t439 * t973 + t594 * t647 - t406;
t112 = -t889 / 0.2e1 + t764 + (-t875 / 0.2e1 + t766) * t594 + t669;
t107 = pkin(5) * t752 + t335 * t972 + t768;
t80 = -t898 / 0.2e1 + t349 * t967 - t566 + t637;
t78 = t349 * t971 + t594 * t724 - t621 + t711;
t69 = pkin(5) * t753 - t593 * t741;
t57 = t614 + t745;
t55 = t593 * t746 - t892 / 0.2e1 + t218 / 0.2e1 + t658;
t50 = (-t883 / 0.2e1 - t956 / 0.2e1) * t893 + t1006;
t48 = pkin(5) * t998 + t317 * t982 + t335 * t981 + t482 * t800 + t573 * t727;
t44 = pkin(9) * t1006 - t653 + t664;
t42 = t203 * t751 - t891 / 0.2e1 + t1000 - t623;
t41 = -t607 + t744;
t34 = t600 - t655;
t30 = t609 + t611;
t28 = t101 * t748 + t750 * t82 + t608 - t652;
t26 = t203 * t752 + t230 * t751 + t540 * t763 - t715 * t884 - t603 + t656;
t25 = -t1003 * t699 + t203 * t749 + t230 * t747 + t596 * t759 - t602 + t657;
t20 = t725 + t736;
t17 = t605 + t867;
t12 = t726 + t737;
t9 = -t964 / 0.2e1 + t610;
t8 = pkin(5) * t1002 + t100 * t981 + t156 * t800 + t573 * t801 + t82 * t982;
t6 = pkin(5) * t1001 + t100 * t995 + t156 * t727 + t482 * t801 + t82 * t996 + t1021;
t3 = t601 - t1027;
t2 = t604 + t622;
t58 = [0, 0, 0, t719, t502 * qJD(2), t598 * t708, -t595 * t708, 0, -t395 * qJD(2), -t396 * qJD(2), t236, qJD(2) * t269 + qJD(3) * t272, -qJD(2) * t372 + t483 * t786, qJD(2) * t371 + t485 * t786, -t719, -qJD(2) * t123 - qJD(3) * t174, qJD(2) * t124 + qJD(3) * t173, qJD(2) * t59 + qJD(3) * t53 + t483 * t783, qJD(2) * t66 + qJD(3) * t108 + t485 * t809, qJD(2) * t65 + qJD(3) * t109 + qJD(4) * t415, qJD(2) * t46 + qJD(3) * t45 + qJD(4) * t144 (-qJD(2) * t439 - t485 * t841 - t833) * t387, qJD(2) * t154 + qJD(3) * t162 + qJD(5) * t175, qJD(2) * t191 + qJD(3) * t225 + t386 * t832, qJD(2) * t190 + qJD(3) * t226 + t387 * t832, t236, qJD(2) * t36 + qJD(3) * t32 - qJD(4) * t206 + qJD(5) * t52, qJD(2) * t37 + qJD(3) * t33 + qJD(4) * t255 + qJD(5) * t51, qJD(2) * t22 + qJD(3) * t21 + qJD(4) * t161 + qJD(5) * t15 + qJD(6) * t238, qJD(2) * t16 + qJD(3) * t14 + qJD(4) * t31 + qJD(5) * t13 + qJD(6) * t38; 0, 0, 0, t718, t808, t691 * t893, -t691 * t894, 0, -qJD(2) * t495 - t811, qJD(2) * t493 - t810, t397 + t1020, t557 * t778 + t234 + t815, t594 * t779 - t852, t597 * t779 + t853, -t456, -t827 + (-t495 * t597 + t594 * t662) * qJD(2) + t168 * qJD(3), t826 + (t495 * t594 + t597 * t662) * qJD(2) + t166 * qJD(3), qJD(2) * t693 + qJD(3) * t50 + t895, t868 + (t352 * t597 + t594 * t644) * qJD(2) + t80 * qJD(3) + t291 * qJD(4), t869 + (-t352 * t594 + t597 * t644) * qJD(2) + t78 * qJD(3) + t777, t911 + (pkin(9) * t693 + t352 * t526) * qJD(2) + t44 * qJD(3) + t136 * qJD(4), t178 * qJD(3) - t194 + (-t788 - t851) * t439, t824 + t131 * qJD(3) - t151 + (t889 - t914) * t806, t818 + (-t274 * t590 + t439 * t594) * qJD(2) + t171 * qJD(3) + t188 * qJD(5), t819 + (-t438 * t594 - t590 * t795) * qJD(2) + t172 * qJD(3) + t186 * qJD(5), -t659 + t1020, t919 + (t148 * t594 + t540 * t438 + (-t378 * t893 + t926) * t597) * qJD(2) + t25 * qJD(3) + t301 * qJD(4) + t42 * qJD(5), t918 + (-t149 * t594 + t540 * t439 + (-t379 * t893 - t927) * t597) * qJD(2) + t26 * qJD(3) + t596 * t777 + t41 * qJD(5), t931 + (-t317 * t439 - t336 * t438 + (t133 * t593 - t143 * t596) * t597) * qJD(2) + t3 * qJD(3) + t112 * qJD(4) + t9 * qJD(5) - t201, t935 + (t133 * t317 + t143 * t336 + t189 * t482) * qJD(2) + t2 * qJD(3) + t17 * qJD(4) + t6 * qJD(5) + t28 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1018, t675, t667, t534 * t485, t541, qJD(2) * t168 - t842 - t857, qJD(2) * t166 + t843 + t858, t899 + t50 * qJD(2) + (-t907 + t960) * qJD(3) - t809, qJD(2) * t80 + t842 + t862, qJD(2) * t78 - t783 - t843 + t861, t913 + t44 * qJD(2) + (-pkin(3) * t298 - qJ(4) * t297) * qJD(3) + t258 * qJD(4), t178 * qJD(2) + t593 * t641 + t215, t131 * qJD(2) - t485 * t731 - t157 + t823, qJD(2) * t171 - t474 * t596 + t855, qJD(2) * t172 + t474 * t593 + t854, t1015, t922 + t25 * qJD(2) + (t230 * t593 - t596 * t698) * qJD(3) - t386 * qJD(4) + t55 * qJD(5), t921 + t26 * qJD(2) + (t230 * t596 + t593 * t698) * qJD(3) - qJD(4) * t387 + t57 * qJD(5), t932 + t3 * qJD(2) + (-t940 - t938 + (-t902 + t903) * t485) * qJD(3) + t12 * qJD(5) - t219, t937 + t2 * qJD(2) + (t115 * t524 + t134 * t523 + t177 * t573) * qJD(3) + t34 * qJD(4) + t8 * qJD(5) + t30 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t667, qJD(2) * t291 + t357, t616, qJD(2) * t136 + qJD(3) * t258 + t860, 0, 0, 0, 0, 0, qJD(2) * t301 - qJD(3) * t386 - t856, -qJD(3) * t387 + t485 * t563 + t816, qJD(2) * t112 + t859, qJD(2) * t17 + qJD(3) * t34 + qJD(5) * t20 + qJD(6) * t196 + t946; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t625, t632, t188 * qJD(2) + t386 * t723, t186 * qJD(2) + t387 * t723, t377, qJD(2) * t42 + qJD(3) * t55 - qJD(5) * t122 + t943, qJD(2) * t41 + qJD(3) * t57 + qJD(5) * t121 + t944, -pkin(5) * t833 + qJD(2) * t9 + qJD(3) * t12 + t947, qJD(2) * t6 + qJD(3) * t8 + qJD(4) * t20 - t101 * t954 + t948; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t631, qJD(2) * t28 + qJD(3) * t30 + qJD(4) * t196 + t945; 0, 0, 0, -t718, -t808, -t598 * t709, t595 * t709, 0, t811, t810, -t397 + t312, t234 - t815, -t597 * t786 + t852, t594 * t786 - t853, t456, qJD(3) * t167 + t827, qJD(3) * t165 - t826, -qJD(3) * t49 - t597 * t783 - t895, -qJD(3) * t79 + qJD(4) * t292 - t868, -qJD(3) * t77 + t777 - t869, -qJD(3) * t43 - qJD(4) * t135 - t911, qJD(3) * t179 + t439 * t851 - t194, qJD(3) * t132 - t151 - t824, -qJD(3) * t170 - qJD(5) * t187 - t818, -qJD(3) * t169 + qJD(5) * t185 - t819, t312 + t659, -qJD(3) * t24 - qJD(4) * t251 - qJD(5) * t39 - t919, -qJD(3) * t23 - qJD(4) * t356 - qJD(5) * t40 - t918, qJD(3) * t4 + qJD(4) * t113 + qJD(5) * t10 - t201 - t931, -qJD(3) * t1 + qJD(4) * t18 - qJD(5) * t5 - qJD(6) * t27 - t935; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, t557 * qJD(3), 0, 0, 0, -pkin(2) * t840, -pkin(2) * t838, 0, qJD(3) * t394 - t594 * t834, -qJD(3) * t393 + qJD(4) * t588 (qJD(3) * t538 - t836) * t526, -t562 * t587 + t590 * t776, -qJD(5) * t511 - 0.2e1 * t594 * t720, -qJD(3) * t509 - t594 * t781, -qJD(3) * t512 + t594 * t782, t562, qJD(3) * t159 + qJD(5) * t271 + t588 * t837, -qJD(3) * t160 - qJD(5) * t270 + t588 * t835, -qJD(3) * t103 + qJD(4) * t475 - qJD(5) * t105 + qJD(6) * t510, qJD(3) * t92 + qJD(4) * t163 + qJD(5) * t104 - qJD(6) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1014, t674, t488, t534 * t594, qJD(1) * t758, -t579 - t660, -t661 + t802, t620 - t905, t579 - t681, -t682 - t802, pkin(9) * t620 - t645, t822 - t480 + (-t587 * t806 + t785) * t594, t825 + t525 + (-0.2e1 * t717 - t731) * t594, t564 + t677, t678 - t784, t1013, t288 * qJD(5) - t1007 * t596 - t539 * t841 + t688, t289 * qJD(5) + t1007 * t593 - t539 * t839 + t687 ((-t325 + t904) * t596 + (-t338 - t901) * t593) * qJD(3) + t69 * qJD(5) + t696 (t325 * t524 + t338 * t523 + t481 * t573) * qJD(3) + t141 * qJD(4) + t48 * qJD(5) + t138 * qJD(6) + t701; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488, t730, -t398, t579 - t643, 0, 0, 0, 0, 0, t593 * t848 + t564 - t817, t588 * t845 - t784 - t812, t679, qJD(3) * t141 + qJD(5) * t107 - t686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t617, t619, -t722 * t873 - t820, t722 * t885 + t821, -t530 / 0.2e1 + t574, qJD(3) * t288 - qJD(5) * t379 + t683, qJD(3) * t289 + qJD(5) * t378 + t684, pkin(5) * t781 + qJD(3) * t69 - t690, qJD(3) * t48 + qJD(4) * t107 - t336 * t954 + t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t676, qJD(3) * t138 + t685; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1018, -t675, t666 (-t846 - t850) * t893, t541, -qJD(2) * t167 + t857, -qJD(2) * t165 - t858, qJD(2) * t49 - t899, qJD(2) * t79 - t862, qJD(2) * t77 - t783 - t861, -qJ(4) * t783 + qJD(2) * t43 - t913, -qJD(2) * t179 + t387 * t790 + t215, -qJD(2) * t132 - t157 - t823, qJD(2) * t170 - t485 * t831 - t855, qJD(2) * t169 - t485 * t829 - t854, -t1015, qJD(2) * t24 - qJD(4) * t274 + qJD(5) * t54 - t922, qJD(2) * t23 + qJD(5) * t56 - t596 * t783 - t921, -qJD(2) * t4 - qJD(5) * t11 - t219 - t932, qJD(2) * t1 + qJD(4) * t35 - qJD(5) * t7 - qJD(6) * t29 - t937; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1014, -t674, t530, -t594 * t789, qJD(1) * t757, t660, t661, t905, t681, t682, t645, t561 * t587 - t480 - t822, t594 * t706 + t525 - t825, -t593 * t830 - t677, -t594 * t829 - t678, -t1013, qJD(5) * t350 - t688, qJD(5) * t351 - t687, -t696 + t1011, qJD(4) * t142 - qJD(5) * t47 + qJD(6) * t137 - t701; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t776, t556 * qJD(5), 0, 0, 0, qJ(4) * t829 + t837, -qJ(4) * t831 + t835, qJD(6) * t555, qJD(4) * t573 + qJD(5) * t253 - qJD(6) * t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t534, -t503, 0, 0, 0, 0, 0, -t671, -t487, 0, qJD(6) * t532 + t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t624, t618, -t593 * t723 - t560, t640, t663, t593 * t828 - t628, t596 * t828 - t627, -t694 + t803, -t523 * t954 + t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t672, qJD(4) * t532 + t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666, -qJD(2) * t292 - t357, -t616, qJ(4) * t786 + qJD(2) * t135 - t860, 0, 0, 0, 0, 0, qJD(2) * t251 + qJD(3) * t274 - qJD(5) * t328 + t856, -t816 + t356 * qJD(2) + (t786 - t832) * t596, -qJD(2) * t113 - t859, -qJD(2) * t18 - qJD(3) * t35 - qJD(5) * t19 + qJD(6) * t195 - t946; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530, -t730, t398, t643, 0, 0, 0, 0, 0, t593 * t670 + t817, t596 * t670 + t812, -t679, -qJD(3) * t142 - t1011 + t686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t534, t503, 0, 0, 0, 0, 0, t671, t487, 0, qJD(6) * t533 - t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t729 - t831, t640, 0, -t689 - t803; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t625, -t632, t187 * qJD(2) + (-qJD(1) * t386 + t841) * t485, -t185 * qJD(2) + t641, t377, qJD(2) * t39 - qJD(3) * t54 + qJD(4) * t328 - t943, qJD(2) * t40 - qJD(3) * t56 + t485 * t835 - t944, -qJD(2) * t10 + qJD(3) * t11 - t947, qJD(2) * t5 + qJD(3) * t7 + qJD(4) * t19 + t387 * t953 - t948; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t617, -t619, t820 + (t787 + t841) * t594, t594 * t642 - t821, t530 / 0.2e1 + t574, -qJD(3) * t350 + t593 * t836 - t683, -qJD(3) * t351 + t594 * t835 - t684, -qJD(3) * t638 + t690, qJD(3) * t47 + qJD(4) * t638 + qJD(6) * t804 - t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t624, -t618, t560 + t790, t704, -t663, t628, t627, t694, -t596 * t953 - t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t729, t704, 0, t689; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t631, qJD(2) * t27 + qJD(3) * t29 - qJD(4) * t195 - t387 * t954 - t945; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t676, -pkin(5) * t782 - qJD(3) * t137 - t685; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t672, pkin(5) * t829 - qJD(4) * t533 - t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t58;