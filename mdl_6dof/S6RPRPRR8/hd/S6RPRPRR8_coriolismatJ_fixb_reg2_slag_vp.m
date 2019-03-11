% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRR8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:39
% EndTime: 2019-03-09 04:00:03
% DurationCPUTime: 19.96s
% Computational Cost: add. (24443->699), mult. (44123->927), div. (0->0), fcn. (50714->8), ass. (0->530)
t628 = sin(qJ(3));
t630 = cos(qJ(3));
t933 = sin(pkin(10));
t934 = cos(pkin(10));
t594 = t628 * t934 + t630 * t933;
t962 = -t594 / 0.2e1;
t629 = cos(qJ(5));
t613 = pkin(3) * t933 + pkin(8);
t946 = pkin(9) + t613;
t591 = t946 * t629;
t626 = sin(qJ(6));
t627 = sin(qJ(5));
t716 = t946 * t627;
t952 = cos(qJ(6));
t711 = t626 * t591 + t952 * t716;
t967 = t711 / 0.2e1;
t759 = qJD(5) + qJD(6);
t467 = t629 * t594;
t740 = t952 * t627;
t434 = -t467 * t626 - t594 * t740;
t462 = t627 * t594;
t540 = t626 * t462;
t616 = t952 * t629;
t702 = t594 * t616;
t438 = t702 - t540;
t458 = t591 * t952 - t626 * t716;
t965 = t458 / 0.2e1;
t968 = -t458 / 0.2e1;
t718 = t965 + t968;
t966 = -t711 / 0.2e1;
t719 = t966 + t967;
t1008 = -t718 * t434 - t719 * t438;
t612 = t933 * t628;
t713 = t934 * t630;
t596 = t713 - t612;
t465 = t627 * t596;
t975 = -pkin(1) - pkin(7);
t758 = qJ(4) - t975;
t679 = t758 * t934;
t590 = t630 * t679;
t678 = t758 * t933;
t494 = -t628 * t678 + t590;
t399 = pkin(5) * t465 + t494;
t541 = t626 * t465;
t440 = t596 * t616 - t541;
t854 = t626 * t629;
t601 = t740 + t854;
t959 = -t601 / 0.2e1;
t614 = -pkin(3) * t934 - pkin(4);
t602 = -t629 * pkin(5) + t614;
t997 = -t602 / 0.2e1;
t1007 = -t399 * t959 - t440 * t997 + t458 * t962;
t864 = t601 * t596;
t512 = -t864 / 0.2e1;
t855 = t626 * t627;
t599 = -t616 + t855;
t960 = -t599 / 0.2e1;
t999 = t399 * t960 + t512 * t602 + t594 * t967;
t1006 = t438 / 0.2e1;
t1005 = t438 * t601;
t866 = t601 * t864;
t893 = t440 * t599;
t227 = -t866 / 0.2e1 - t893 / 0.2e1;
t1002 = t759 * t227;
t665 = t630 * t678;
t636 = -t628 * t679 - t665;
t615 = t628 * pkin(3) + qJ(2);
t650 = t594 * pkin(4) - t596 * pkin(8) + t615;
t297 = t627 * t636 - t629 * t650;
t469 = t629 * t596;
t249 = -pkin(9) * t469 - t297;
t949 = t594 * pkin(5);
t633 = t249 + t949;
t632 = t626 * t633;
t298 = t627 * t650 + t629 * t636;
t250 = -pkin(9) * t465 + t298;
t744 = t952 * t250;
t140 = t744 + t632;
t858 = t626 * t249;
t159 = -t744 - t858;
t1001 = t159 + t140;
t220 = t952 * t633;
t857 = t626 * t250;
t139 = -t220 + t857;
t745 = t952 * t249;
t160 = t745 - t857;
t1000 = t160 + t139;
t947 = t630 * pkin(3);
t474 = t596 * pkin(4) + t594 * pkin(8) + t947;
t443 = t627 * t474;
t703 = t758 * t628;
t496 = t703 * t933 - t590;
t479 = t629 * t496;
t300 = t479 + t443;
t253 = pkin(9) * t462 + t300;
t743 = t952 * t253;
t700 = -t743 / 0.2e1;
t444 = t629 * t474;
t853 = t627 * t496;
t299 = t444 - t853;
t951 = pkin(5) * t596;
t229 = pkin(9) * t467 + t299 + t951;
t859 = t626 * t229;
t658 = -t859 / 0.2e1 + t700;
t50 = t658 - t999;
t746 = t952 * t229;
t856 = t626 * t253;
t657 = -t856 / 0.2e1 + t746 / 0.2e1;
t52 = t657 + t1007;
t49 = t657 - t1007;
t998 = 0.2e1 * t627;
t973 = t434 / 0.2e1;
t958 = t601 / 0.2e1;
t768 = t596 * qJD(1);
t506 = t594 * t768;
t585 = t713 / 0.2e1 - t612 / 0.2e1;
t985 = t585 * qJD(5) + t506;
t755 = -t952 / 0.2e1;
t865 = t601 * t594;
t956 = -t626 / 0.2e1;
t646 = (-t438 * t755 - t865 * t956) * pkin(5);
t989 = t139 / 0.2e1;
t721 = t160 / 0.2e1 + t989;
t974 = t140 / 0.2e1;
t722 = t974 + t159 / 0.2e1;
t3 = t440 * t718 + t599 * t721 + t601 * t722 + t719 * t864 + t646;
t598 = t599 ^ 2;
t586 = t598 / 0.2e1;
t988 = -t598 / 0.2e1;
t266 = t586 + t988;
t786 = t266 * qJD(4);
t994 = t3 * qJD(1) - t786;
t797 = qJD(3) * t601;
t993 = -qJD(1) * t227 + t599 * t797;
t804 = qJD(1) * t440;
t992 = qJD(3) * t227 - t864 * t804;
t699 = -t616 / 0.2e1;
t963 = t541 / 0.2e1;
t320 = t963 + (t699 + t960) * t596;
t784 = t320 * qJD(2);
t991 = t711 * t759 - t784;
t990 = t759 * t458;
t712 = t469 * t998;
t954 = -t627 / 0.2e1;
t941 = qJD(3) * pkin(3);
t987 = (-t594 * t934 + t596 * t933) * t941;
t707 = t759 * t599;
t592 = t594 ^ 2;
t593 = t596 ^ 2;
t984 = -t593 - t592;
t757 = t593 - t592;
t760 = t593 / 0.2e1;
t698 = t592 / 0.2e1 + t760;
t983 = qJD(6) * t585 + t985;
t656 = t854 / 0.2e1 + t740 / 0.2e1;
t317 = (t958 - t656) * t596;
t785 = t317 * qJD(2);
t796 = qJD(3) * t602;
t982 = qJD(1) * t49 - t601 * t796 + t785;
t754 = t952 / 0.2e1;
t953 = -t629 / 0.2e1;
t33 = (t864 * t954 + (t599 * t953 + t754) * t596) * pkin(5) + t49;
t950 = pkin(5) * t627;
t429 = t599 * t950 + t601 * t602;
t981 = qJD(1) * t33 - qJD(3) * t429 + t785;
t971 = -t438 / 0.2e1;
t672 = t1006 * t434 - t865 * t971;
t980 = qJD(3) * t672;
t147 = t746 - t856;
t148 = t743 + t859;
t14 = -t139 * t438 + t140 * t865 - t147 * t440 - t148 * t864;
t978 = t14 * qJD(1) + t672 * qJD(2);
t976 = t601 ^ 2;
t972 = t865 / 0.2e1;
t964 = t494 / 0.2e1;
t961 = -t596 / 0.2e1;
t625 = t629 ^ 2;
t957 = -t625 / 0.2e1;
t955 = t626 / 0.2e1;
t948 = t626 * pkin(5);
t943 = pkin(5) * qJD(5);
t741 = t952 * t599;
t863 = t601 * t626;
t659 = t863 / 0.2e1 - t741 / 0.2e1;
t9 = -t721 * t438 - t722 * t434 + (t629 * t760 + t659) * pkin(5);
t939 = t9 * qJD(1);
t495 = -t934 * t703 - t665;
t908 = t300 * t629;
t911 = t299 * t627;
t912 = t298 * t629;
t913 = t297 * t627;
t37 = (t912 / 0.2e1 + t913 / 0.2e1 - t495 / 0.2e1) * t596 + (t908 / 0.2e1 - t911 / 0.2e1 + t964) * t594;
t935 = t37 * qJD(3);
t750 = -t949 / 0.2e1;
t631 = -t632 / 0.2e1 + t626 * t750;
t43 = t858 / 0.2e1 + t631;
t931 = qJD(1) * t43;
t756 = pkin(5) * t469;
t905 = t399 * t440;
t55 = t159 * t594 + t756 * t864 + t905;
t930 = qJD(1) * t55;
t906 = t399 * t864;
t56 = -t160 * t594 + t440 * t756 - t906;
t929 = qJD(1) * t56;
t64 = t139 * t594 - t906;
t928 = qJD(1) * t64;
t65 = -t140 * t594 + t905;
t927 = qJD(1) * t65;
t881 = t494 * t596;
t93 = t881 + (-t912 - t913) * t594;
t926 = qJD(1) * t93;
t923 = t139 * t599;
t918 = t140 * t601;
t17 = -t1000 * t864 - t1001 * t440;
t917 = t17 * qJD(1);
t400 = -pkin(5) * t462 + t495;
t20 = -t139 * t147 + t140 * t148 + t399 * t400;
t916 = t20 * qJD(1);
t23 = -t139 * t159 + t140 * t160 + t399 * t756;
t915 = t23 * qJD(1);
t29 = -t139 * t865 - t140 * t438 + t399 * t596;
t914 = t29 * qJD(1);
t910 = t299 * t629;
t909 = t300 * t627;
t38 = t918 + t923;
t907 = t38 * qJD(1);
t901 = t865 * t440;
t900 = t865 * t594;
t899 = t864 * t438;
t897 = t438 * t594;
t896 = t438 * t599;
t894 = t440 * t596;
t892 = t440 * t601;
t891 = t440 * t626;
t662 = -t220 / 0.2e1 + t952 * t750;
t45 = t745 / 0.2e1 + t662;
t890 = t45 * qJD(1);
t888 = t711 * t864;
t884 = t672 * qJD(1);
t883 = t494 * t495;
t882 = t494 * t594;
t880 = t494 * t627;
t879 = t494 * t629;
t878 = t495 * t596;
t877 = t495 * t627;
t876 = t495 * t629;
t875 = t496 * t594;
t689 = t297 * t629 - t298 * t627;
t53 = (t909 + t910) * t596 + t689 * t594;
t874 = t53 * qJD(1);
t873 = t594 * t596;
t872 = t594 * t614;
t871 = t596 * t864;
t339 = t596 * t599;
t870 = t596 * t614;
t869 = t599 * t864;
t868 = t599 * t594;
t867 = t601 * t865;
t860 = t613 * t594;
t73 = (-t297 + t877) * t596 + (t299 - t880) * t594;
t852 = t73 * qJD(1);
t74 = (-t298 + t876) * t596 + (-t300 - t879) * t594;
t851 = t74 * qJD(1);
t624 = t627 ^ 2;
t717 = t957 - t624 / 0.2e1;
t644 = t717 * t860 + t870 / 0.2e1;
t674 = -t910 / 0.2e1 - t909 / 0.2e1;
t95 = t644 + t674;
t850 = t95 * qJD(1);
t150 = -t599 * t718 + t601 * t719;
t841 = t150 * qJD(5);
t175 = t869 - t892;
t840 = t759 * t175;
t671 = t438 * t962 - t894 / 0.2e1;
t696 = t616 / 0.2e1 - t855 / 0.2e1;
t204 = -t671 + t696;
t316 = (t958 + t656) * t594;
t839 = t204 * qJD(2) + t316 * qJD(4);
t673 = t434 * t962 + t871 / 0.2e1;
t205 = -t656 - t673;
t724 = -t868 / 0.2e1;
t829 = -t540 / 0.2e1 + t702 / 0.2e1;
t319 = t724 + t829;
t838 = t205 * qJD(2) + t319 * qJD(4);
t206 = -t656 + t673;
t723 = t868 / 0.2e1;
t322 = t723 + t829;
t837 = t206 * qJD(2) + t322 * qJD(4);
t207 = t671 + t696;
t325 = -t865 / 0.2e1 + t656 * t594;
t836 = t207 * qJD(2) + t325 * qJD(4);
t833 = t759 * t266;
t831 = t759 * t316;
t830 = t759 * t325;
t828 = t540 / 0.2e1 + t594 * t699;
t827 = t624 + t625;
t606 = t625 - t624;
t761 = -t593 / 0.2e1;
t648 = -t1006 * t438 + t434 * t972 + t761;
t107 = -t976 / 0.2e1 + t988 + t648;
t826 = qJD(1) * t107;
t108 = (t973 + t972) * t596;
t825 = qJD(1) * t108;
t109 = (t971 + t1006) * t596;
t824 = qJD(1) * t109;
t130 = -t867 + t896;
t823 = qJD(1) * t130;
t177 = t899 - t901;
t822 = qJD(1) * t177;
t821 = qJD(1) * t689;
t193 = t297 * t594 - t465 * t494;
t820 = qJD(1) * t193;
t194 = -t298 * t594 + t469 * t494;
t819 = qJD(1) * t194;
t818 = qJD(1) * t204;
t817 = qJD(1) * t205;
t221 = -t871 + t900;
t816 = qJD(1) * t221;
t222 = t871 + t900;
t815 = qJD(1) * t222;
t223 = t894 - t897;
t814 = qJD(1) * t223;
t224 = t894 + t897;
t813 = qJD(1) * t224;
t228 = -t866 + t893;
t812 = qJD(1) * t228;
t256 = -t594 * t636 + t881;
t811 = qJD(1) * t256;
t337 = -0.1e1 / 0.2e1 - t698;
t326 = t337 * t627;
t810 = qJD(1) * t326;
t383 = t337 * t629;
t809 = qJD(1) * t383;
t401 = t757 * t627;
t808 = qJD(1) * t401;
t402 = t984 * t627;
t807 = qJD(1) * t402;
t403 = t757 * t629;
t806 = qJD(1) * t403;
t805 = qJD(1) * t864;
t441 = t594 * t947 + t615 * t596;
t803 = qJD(1) * t441;
t442 = -t615 * t594 + t596 * t947;
t802 = qJD(1) * t442;
t478 = t984 * t629;
t801 = qJD(1) * t478;
t800 = qJD(1) * t615;
t799 = qJD(2) * t594;
t798 = qJD(3) * t266;
t795 = qJD(3) * t629;
t794 = qJD(4) * t629;
t793 = qJD(5) * t627;
t792 = qJD(5) * t629;
t791 = qJD(6) * t602;
t635 = t636 * t596;
t166 = t875 / 0.2e1 + t635 / 0.2e1 - t878 / 0.2e1 + t882 / 0.2e1;
t790 = t166 * qJD(1);
t167 = -t635 - t875 + t878 - t882;
t789 = t167 * qJD(1);
t178 = t899 + t901;
t788 = t178 * qJD(1);
t655 = t592 * t717 + t761;
t234 = t655 + t717;
t787 = t234 * qJD(1);
t302 = t316 * qJD(1);
t318 = t723 + t828;
t304 = t318 * qJD(1);
t306 = t319 * qJD(1);
t783 = t337 * qJD(1);
t782 = t339 * qJD(1);
t341 = 0.2e1 * t512;
t781 = t341 * qJD(1);
t645 = t933 * t962 + t934 * t961;
t425 = (-t630 / 0.2e1 + t645) * pkin(3);
t780 = t425 * qJD(1);
t779 = t757 * qJD(1);
t461 = (t624 / 0.2e1 + t957) * t596;
t778 = t461 * qJD(5);
t777 = t462 * qJD(1);
t776 = t465 * qJD(1);
t775 = t467 * qJD(1);
t577 = t624 * t594;
t578 = t625 * t594;
t475 = t577 + t578;
t774 = t475 * qJD(1);
t477 = t827 * t596;
t773 = t477 * qJD(1);
t772 = t984 * qJD(1);
t771 = t585 * qJD(1);
t769 = t594 * qJD(1);
t580 = t594 * qJD(3);
t583 = t596 * qJD(3);
t605 = t628 ^ 2 - t630 ^ 2;
t767 = t605 * qJD(1);
t766 = t628 * qJD(1);
t765 = t628 * qJD(3);
t764 = t630 * qJD(1);
t763 = t630 * qJD(3);
t753 = t951 / 0.2e1;
t752 = -t950 / 0.2e1;
t751 = t950 / 0.2e1;
t225 = t892 / 0.2e1 + t869 / 0.2e1;
t748 = t225 * qJD(3);
t747 = qJD(3) * t975;
t742 = t952 * t864;
t739 = qJ(2) * t766;
t738 = qJ(2) * t764;
t736 = t599 * t769;
t735 = t601 * t769;
t734 = t627 * t769;
t732 = t627 * t795;
t731 = t594 * t793;
t730 = t594 * t792;
t505 = t594 * t583;
t729 = t627 * t792;
t728 = t629 * t769;
t727 = t629 * t768;
t726 = t628 * t763;
t725 = t629 * t761;
t720 = -t443 / 0.2e1 - t479 / 0.2e1;
t715 = t952 * qJD(5);
t714 = t952 * qJD(6);
t709 = t759 * t594;
t708 = t759 * t601;
t706 = t629 * t753;
t704 = -qJD(5) - t769;
t701 = t593 * t729;
t697 = qJD(3) * t712;
t695 = -t465 * qJD(3) - t730;
t694 = -qJD(6) + t704;
t54 = -t297 * t299 + t298 * t300 + t883;
t693 = t54 * qJD(1) + t37 * qJD(2);
t688 = t908 - t911;
t687 = -t596 * t613 - t872;
t27 = -t139 * t596 + t147 * t594 - t399 * t865 + t400 * t864;
t686 = t27 * qJD(1) + t108 * qJD(2);
t28 = -t140 * t596 - t148 * t594 - t399 * t438 + t400 * t440;
t685 = t28 * qJD(1) + t109 * qJD(2);
t643 = -t865 * t967 - t438 * t965 + t596 * t602 / 0.2e1;
t675 = t147 * t599 / 0.2e1 + t148 * t959;
t26 = t643 + t675;
t684 = -qJD(1) * t26 + qJD(2) * t225;
t257 = (-0.1e1 + t827) * t873;
t683 = -t37 * qJD(1) - t257 * qJD(2);
t647 = (-t438 * t955 - t755 * t865) * pkin(5);
t18 = t599 * t722 - t601 * t721 + t647;
t682 = -t18 * qJD(1) + t150 * qJD(3);
t32 = (t440 * t954 + (t601 * t953 + t956) * t596) * pkin(5) + t50;
t430 = -t599 * t602 + t601 * t950;
t681 = qJD(1) * t32 - qJD(3) * t430;
t677 = t704 * t629;
t195 = t496 * t636 + t615 * t947 + t883;
t676 = t195 * qJD(1) + t166 * qJD(2);
t208 = -t440 ^ 2 + t864 ^ 2;
t69 = qJD(1) * t208 + qJD(3) * t175;
t431 = t598 - t976;
t137 = qJD(1) * t175 + qJD(3) * t431;
t670 = t458 * t959 + t711 * t960;
t669 = t860 / 0.2e1 - t870 / 0.2e1;
t667 = qJD(1) * t50 + t599 * t796;
t666 = t596 * t677;
t640 = t669 * t627 + t879 / 0.2e1;
t183 = t640 - t720;
t664 = -qJD(1) * t183 - t614 * t795;
t652 = t669 * t629;
t185 = -t444 / 0.2e1 - t652 + (t496 / 0.2e1 + t964) * t627;
t663 = -qJD(3) * t614 * t627 - qJD(1) * t185;
t405 = -qJD(1) * t461 + t732;
t661 = t147 * t754 + t148 * t955;
t660 = t891 / 0.2e1 - t742 / 0.2e1;
t334 = qJD(1) * t593 * t627 * t629 + qJD(3) * t461;
t476 = t606 * t593;
t654 = qJD(1) * t476 + t697;
t653 = qJD(1) * t712 - qJD(3) * t606;
t161 = -t434 * t864 + t438 * t440 - t873;
t634 = t864 * t989 + t440 * t974 + t147 * t973 + t148 * t1006 + t399 * t594 / 0.2e1 + t400 * t961;
t6 = t634 + t670;
t651 = t6 * qJD(1) + t161 * qJD(2) + t225 * qJD(4);
t639 = t1000 * t968 + t159 * t967 + t711 * t974;
t1 = (t399 * t954 + t469 * t997 + t661) * pkin(5) + t639;
t174 = t602 * t950;
t41 = (t465 / 0.2e1 + t660) * pkin(5) + t1008;
t641 = -t1 * qJD(1) - t41 * qJD(2) + t174 * qJD(3) + t150 * qJD(4);
t623 = qJ(2) * qJD(2);
t622 = qJD(1) * qJ(2);
t611 = t628 * t764;
t573 = t585 * qJD(3);
t569 = t629 * t583;
t451 = t462 * qJD(5);
t424 = t947 / 0.2e1 + t645 * pkin(3);
t423 = -t777 - t793;
t384 = t592 * t953 + t725 + t629 / 0.2e1;
t342 = t864 / 0.2e1 + t512;
t336 = 0.1e1 / 0.2e1 - t698;
t327 = t627 * t698 + t954;
t324 = -t596 * t656 + t512;
t323 = t339 / 0.2e1 + t596 * t699 + t963;
t321 = t724 + t828;
t233 = t655 - t717;
t232 = -t708 - t302;
t231 = t707 - t306;
t230 = -t304 - t707;
t197 = qJD(3) * t316 + t440 * t769;
t196 = qJD(3) * t318 + t769 * t864;
t186 = t880 / 0.2e1 - t853 / 0.2e1 + t444 / 0.2e1 - t652;
t184 = t640 + t720;
t164 = t166 * qJD(3);
t154 = qJD(3) * t320 + t817;
t153 = qJD(3) * t317 + t818;
t134 = qJD(3) * t325 + t440 * t694;
t133 = qJD(3) * t321 + t694 * t864;
t106 = t976 / 0.2e1 + t586 + t648;
t105 = t109 * qJD(3);
t104 = t108 * qJD(3);
t94 = t644 - t674;
t72 = qJD(3) * t324 - t438 * t759 - t818;
t71 = qJD(3) * t323 - t434 * t759 - t817;
t51 = t658 + t999;
t46 = t857 - t745 / 0.2e1 + t662;
t44 = -t744 - t858 / 0.2e1 + t631;
t42 = pkin(5) * t660 + t596 * t752 - t1008;
t35 = t440 * t751 + t601 * t706 + t700 + (-t951 / 0.2e1 - t229 / 0.2e1) * t626 + t999;
t34 = t599 * t706 - t864 * t752 + t952 * t753 + t52;
t25 = t643 - t675;
t19 = t1000 * t958 + t1001 * t960 + t647;
t10 = (t659 + t725) * pkin(5) + t1001 * t973 + t1000 * t1006;
t5 = t634 - t670;
t4 = -t864 * t966 + t160 * t960 + t159 * t959 + t646 - t918 / 0.2e1 - t888 / 0.2e1 - t923 / 0.2e1;
t2 = pkin(5) * t661 + t399 * t751 + t602 * t706 - t639;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t623, -t726, t605 * qJD(3), 0, t726, 0, 0, qJ(2) * t763 + qJD(2) * t628, -qJ(2) * t765 + qJD(2) * t630, 0, t623, -t505, -t757 * qJD(3), 0, t505, 0, 0, qJD(3) * t441 + t799, qJD(2) * t596 + qJD(3) * t442, qJD(3) * t167 - qJD(4) * t984, qJD(2) * t615 + qJD(3) * t195 + qJD(4) * t256, -t505 * t625 - t701, -qJD(5) * t476 + t594 * t697, qJD(3) * t403 - t596 * t731, -t505 * t624 + t701, -qJD(3) * t401 - t596 * t730, t505, qJD(3) * t73 - qJD(4) * t402 + qJD(5) * t194 + t629 * t799, qJD(3) * t74 - qJD(4) * t478 + qJD(5) * t193 - t627 * t799, -qJD(2) * t477 - qJD(3) * t53, -qJD(2) * t689 + qJD(3) * t54 + qJD(4) * t93 (-qJD(3) * t438 - t759 * t864) * t440, qJD(3) * t178 + t208 * t759, qJD(3) * t223 - t709 * t864 (-qJD(3) * t865 + t440 * t759) * t864, qJD(3) * t221 - t440 * t709, t505, qJD(3) * t27 + qJD(4) * t222 + qJD(5) * t55 + qJD(6) * t65 - t599 * t799, qJD(3) * t28 + qJD(4) * t224 + qJD(5) * t56 + qJD(6) * t64 - t601 * t799, qJD(2) * t228 + qJD(3) * t14 + qJD(4) * t177 + qJD(5) * t17, qJD(2) * t38 + qJD(3) * t20 + qJD(4) * t29 + qJD(5) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t622, 0, 0, 0, 0, 0, 0, t766, t764, 0, t622, 0, 0, 0, 0, 0, 0, t769, t768, 0, qJD(4) * t336 + t164 + t800, 0, 0, 0, 0, 0, 0, qJD(5) * t384 + t728, qJD(5) * t327 - t734, -t773, qJD(4) * t233 - t821 + t935, 0, 0, 0, 0, 0, 0, t207 * t759 + t104 - t736, t206 * t759 + t105 - t735, t812 + t980, t907 + (-t434 * t599 + t1005) * qJD(2) + t5 * qJD(3) + t106 * qJD(4) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t611, t767, -t765, t611, -t763, 0, -t628 * t747 + t738, -t630 * t747 - t739, 0, 0, -t506, -t779, -t580, t506, -t583, 0, -qJD(3) * t495 + t803, -qJD(3) * t496 + t802, t789 - t987 (-t495 * t934 + t496 * t933) * t941 + t424 * qJD(4) + t676, -t778 + (-t625 * t768 - t732) * t594 (t577 - t578) * qJD(3) + (-qJD(5) + t769) * t712, t583 * t627 + t806, t778 + (-t624 * t768 + t732) * t594, t569 - t808, t985, t852 + (t627 * t687 - t876) * qJD(3) + t186 * qJD(5), t851 + (t629 * t687 + t877) * qJD(3) + t184 * qJD(5), qJD(3) * t688 - t874 (t495 * t614 + t613 * t688) * qJD(3) + t94 * qJD(4) + t693 -(t797 + t804) * t438 + t1002, t788 + (t867 + t896) * qJD(3) + t840, t321 * t759 + t583 * t601 + t814 -(qJD(3) * t599 + t805) * t865 - t1002, -t583 * t599 + t816 + t830, t983 (t400 * t599 - t596 * t711 - t602 * t865) * qJD(3) + t34 * qJD(5) + t52 * qJD(6) + t686 (t400 * t601 - t438 * t602 - t458 * t596) * qJD(3) + t342 * qJD(4) + t35 * qJD(5) + t51 * qJD(6) + t685 (-t147 * t601 - t148 * t599 - t438 * t711 + t458 * t865) * qJD(3) + t4 * qJD(5) + t978, t916 + t5 * qJD(2) + (-t147 * t711 + t148 * t458 + t400 * t602) * qJD(3) + t25 * qJD(4) + t2 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t772, qJD(2) * t336 + qJD(3) * t424 + t811, 0, 0, 0, 0, 0, 0, -t807, -t801, 0, qJD(2) * t233 + qJD(3) * t94 + t926, 0, 0, 0, 0, 0, 0, t815 + t830, qJD(3) * t342 + t322 * t759 + t813, t822, t914 + t106 * qJD(2) + t25 * qJD(3) + (-t599 * t865 - t1005) * qJD(4) + t19 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t334, -t654, t704 * t465, t334, t666, t573, qJD(2) * t384 + qJD(3) * t186 - qJD(5) * t298 + t819, qJD(2) * t327 + qJD(3) * t184 + qJD(5) * t297 + t820, 0, 0, t992, t69, t133, -t992, t134, t573, qJD(3) * t34 + qJD(5) * t159 + qJD(6) * t44 + t836 + t930, qJD(3) * t35 - qJD(5) * t160 + qJD(6) * t46 + t837 + t929, t917 + t4 * qJD(3) + (t742 - t891) * t943, t915 + t10 * qJD(2) + t2 * qJD(3) + t19 * qJD(4) + (t159 * t952 + t160 * t626) * t943; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t992, t69, t133, -t992, t134, t573, qJD(3) * t52 + qJD(5) * t44 - qJD(6) * t140 + t836 + t927, qJD(3) * t51 + qJD(5) * t46 + qJD(6) * t139 + t837 + t928, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t622, 0, 0, 0, 0, 0, 0, -t766, -t764, 0, -t622, 0, 0, 0, 0, 0, 0, -t769, -t768, 0, qJD(4) * t337 + t164 - t800, 0, 0, 0, 0, 0, 0, qJD(5) * t383 - t728, -qJD(5) * t326 + t734, t773, qJD(4) * t234 + t821 + t935, 0, 0, 0, 0, 0, 0, -t204 * t759 + t104 + t736, -t205 * t759 + t105 + t735, -t812 + t980, qJD(3) * t6 + qJD(4) * t107 - qJD(5) * t9 - t907; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t765, -t763, 0, 0, 0, 0, 0, 0, 0, 0, -t580, -t583, 0, t790 + t987, 0, 0, 0, 0, 0, 0, -qJD(5) * t465 - t580 * t629, -qJD(5) * t469 + t580 * t627, t477 * qJD(3) (t477 * t613 + t872) * qJD(3) - t683, 0, 0, 0, 0, 0, 0, t324 * t759 + t580 * t599 + t825, t323 * t759 + t580 * t601 + t824, -qJD(3) * t228 + t884 (t440 * t458 + t594 * t602 + t888) * qJD(3) + t42 * qJD(5) + t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t783, 0, 0, 0, 0, 0, 0, 0, 0, 0, t787, 0, 0, 0, 0, 0, 0, 0, 0, 0, t748 + t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t695 + t809, -qJD(3) * t469 + t731 - t810, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, 0, -t939 + t42 * qJD(3) + (t434 * t626 - t438 * t952) * t943; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t611, -t767, 0, -t611, 0, 0, -t738, t739, 0, 0, t506, t779, 0, -t506, 0, 0, -qJD(4) * t596 - t803, qJD(4) * t594 - t802, -t789, qJD(4) * t425 - t676, t506 * t625 - t778, t666 * t998, qJD(5) * t467 - t806, t506 * t624 + t778, -t451 + t808, -t985, qJD(5) * t185 - t596 * t794 - t852, qJD(4) * t465 + qJD(5) * t183 - t851, -qJD(4) * t475 + t874, qJD(4) * t95 - t693, t438 * t804 + t1002, -t788 + t840, -t318 * t759 - t814, t805 * t865 - t1002, -t816 - t831, -t983, qJD(4) * t339 - qJD(5) * t33 - qJD(6) * t49 - t686, -qJD(4) * t341 - qJD(5) * t32 - qJD(6) * t50 - t685, qJD(4) * t130 - qJD(5) * t3 - t978, -qJD(2) * t6 + qJD(4) * t26 - qJD(5) * t1 - t916; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t790, 0, 0, 0, 0, 0, 0, 0, 0, 0, t683, 0, 0, 0, 0, 0, 0, -t317 * t759 - t825, -t320 * t759 - t824, -t884, -qJD(5) * t41 - t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t729, t606 * qJD(5), 0, -t729, 0, 0, t614 * t793, t614 * t792, 0, 0, -t599 * t708, t759 * t431, 0, t601 * t707, 0, 0, qJD(5) * t429 + t601 * t791, qJD(5) * t430 - t599 * t791, 0, qJD(5) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t768, t769, 0, t780, 0, 0, 0, 0, 0, 0, -t727, t776, -t774, t850, 0, 0, 0, 0, 0, 0, t782, -t781, t823 + t833, -t684 + t841; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t405, -t653, t775 + t792, -t405, t423, -t771, -t613 * t792 - t663, t613 * t793 - t664, 0, 0, -t993, t137, t230, t993, t232, -t771, -t981 - t990, -t681 + t991 (t741 - t863) * t943 - t994 (-t458 * t952 - t626 * t711) * t943 + t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t993, t137, t230, t993, t232, -t771, -t982 - t990, -t667 + t991, t786, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583, -t580, t772, -qJD(2) * t337 - qJD(3) * t425 - t811, 0, 0, 0, 0, 0, 0, -t451 + t569 + t807, t695 + t801, t475 * qJD(3), -qJD(2) * t234 - qJD(3) * t95 - t926, 0, 0, 0, 0, 0, 0, -qJD(3) * t339 - t815 - t831, qJD(3) * t341 - t319 * t759 - t813, -qJD(3) * t130 - t822, -qJD(2) * t107 - qJD(3) * t26 - qJD(5) * t18 - t914; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t787, 0, 0, 0, 0, 0, 0, 0, 0, 0, t748 - t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t768, -t769, 0, -t780, 0, 0, 0, 0, 0, 0, t727, -t776, t774, -t850, 0, 0, 0, 0, 0, 0, -t782, t781, -t823 + t833, t684 + t841; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t423, t677, 0, 0, 0, 0, 0, 0, 0, 0, t232, t231, t798 (-t599 * t626 - t601 * t952) * t943 + t682; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t231, t798, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t654, -qJD(3) * t467 + t506 * t627, -t334, qJD(3) * t462 + t594 * t727, t573, -qJD(2) * t383 - qJD(3) * t185 + qJD(4) * t462 - t819, qJD(2) * t326 - qJD(3) * t183 + t594 * t794 - t820, 0, 0, -t992, -t69, t196, t992, t197, t573, qJD(3) * t33 + qJD(6) * t43 + t839 - t930, qJD(3) * t32 + qJD(6) * t45 + t838 - t929, qJD(3) * t3 - t917, qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t18 - t915; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t809, t810, 0, 0, 0, 0, 0, 0, 0, 0, t153, t154, 0, qJD(3) * t41 + t939; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405, t653, -t775, t405, t777, t771, t663, t664, 0, 0, t993, -t137, t304, -t993, t302, t771, t981, t784 + t681, t994, -t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t777, t728, 0, 0, 0, 0, 0, 0, 0, 0, t302, t306, -t798, -t682; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t948, -pkin(5) * t714, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t759 * t948 + t931, t890 + (-t715 - t714) * pkin(5), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t992, -t69, t196, t992, t197, t573, qJD(3) * t49 - qJD(5) * t43 + t839 - t927, qJD(3) * t50 - qJD(5) * t45 + t838 - t928, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t154, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t993, -t137, t304, -t993, t302, t771, t982, t784 + t667, -t786, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t306, -t798, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t626 * t943 - t931, pkin(5) * t715 - t890, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
