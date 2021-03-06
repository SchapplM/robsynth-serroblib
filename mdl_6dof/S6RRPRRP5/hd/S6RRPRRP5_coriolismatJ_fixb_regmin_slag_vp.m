% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRP5
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
% cmat_reg [(6*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:29
% EndTime: 2019-03-09 12:06:04
% DurationCPUTime: 18.26s
% Computational Cost: add. (25249->904), mult. (63286->1239), div. (0->0), fcn. (70076->10), ass. (0->658)
t881 = cos(pkin(11));
t546 = -pkin(2) * t881 - pkin(3);
t560 = sin(qJ(4));
t563 = cos(qJ(4));
t658 = -t563 * pkin(4) - t560 * pkin(10);
t504 = t658 + t546;
t562 = cos(qJ(5));
t491 = t562 * t504;
t559 = sin(qJ(5));
t557 = sin(pkin(11));
t545 = pkin(2) * t557 + pkin(9);
t813 = t563 * t545;
t736 = t559 * t813;
t433 = -t491 + t736;
t821 = t560 * t562;
t739 = qJ(6) * t821;
t398 = -t433 - t739;
t925 = t398 / 0.2e1;
t389 = -t739 + t491 + (-t545 * t559 - pkin(5)) * t563;
t927 = -t389 / 0.2e1;
t685 = t925 + t927;
t882 = cos(pkin(6));
t683 = t882 * t560;
t558 = sin(pkin(6));
t561 = sin(qJ(2));
t681 = t881 * t561;
t564 = cos(qJ(2));
t832 = t557 * t564;
t488 = (-t681 - t832) * t558;
t837 = t488 * t563;
t440 = t683 - t837;
t830 = t558 * t564;
t831 = t558 * t561;
t486 = t557 * t831 - t830 * t881;
t328 = t440 * t559 - t486 * t562;
t726 = pkin(1) * t882;
t543 = t564 * t726;
t892 = pkin(8) + qJ(3);
t476 = -t831 * t892 + t543;
t445 = pkin(2) * t882 + t476;
t669 = t561 * t726;
t477 = t830 * t892 + t669;
t682 = t881 * t477;
t319 = t557 * t445 + t682;
t291 = pkin(9) * t882 + t319;
t817 = t563 * t291;
t727 = -pkin(2) * t564 - pkin(1);
t337 = t486 * pkin(3) + t488 * pkin(9) + t558 * t727;
t822 = t560 * t337;
t210 = t817 + t822;
t898 = t486 * pkin(10);
t186 = t210 + t898;
t456 = t557 * t477;
t318 = t445 * t881 - t456;
t290 = -pkin(3) * t882 - t318;
t838 = t488 * t560;
t438 = -t563 * t882 - t838;
t659 = t438 * pkin(4) - t440 * pkin(10);
t566 = t290 + t659;
t98 = t562 * t186 + t559 * t566;
t81 = -t328 * qJ(6) + t98;
t944 = t685 * t81;
t485 = t486 ^ 2;
t943 = t485 * t563;
t552 = t558 ^ 2;
t942 = t552 * t561;
t209 = t291 * t560 - t563 * t337;
t185 = -pkin(4) * t486 + t209;
t941 = (t209 / 0.2e1 - t185 / 0.2e1) * t559;
t940 = -t389 + t398;
t818 = t562 * t563;
t695 = -t818 / 0.2e1;
t839 = t488 * t559;
t939 = t486 * t695 - t839 / 0.2e1;
t371 = t563 * t486;
t692 = t371 / 0.2e1;
t938 = t562 * t692 + t839 / 0.2e1;
t907 = -t560 / 0.2e1;
t706 = t438 * t907;
t903 = t563 / 0.2e1;
t618 = t440 * t903 + t706;
t937 = t618 * qJD(4);
t747 = t563 * qJD(2);
t541 = t560 * t747;
t675 = t618 * qJD(1) + t541;
t754 = t438 * qJD(1);
t610 = -qJD(2) * t618 + t440 * t754;
t842 = t486 * t559;
t845 = t440 * t562;
t330 = t842 + t845;
t936 = t330 ^ 2;
t935 = t438 ^ 2;
t97 = t186 * t559 - t562 * t566;
t80 = -qJ(6) * t330 - t97;
t899 = t438 * pkin(5);
t64 = t80 + t899;
t934 = t64 / 0.2e1;
t339 = t476 * t557 + t682;
t893 = t563 * pkin(10);
t895 = t560 * pkin(4);
t526 = -t893 + t895;
t229 = -t486 * t526 + t339;
t224 = t562 * t229;
t340 = t476 * t881 - t456;
t327 = t563 * t340;
t745 = pkin(2) * t831;
t354 = -pkin(3) * t488 + pkin(9) * t486 + t745;
t347 = t560 * t354;
t809 = t327 + t347;
t198 = -pkin(10) * t488 + t809;
t829 = t559 * t198;
t113 = t224 - t829;
t366 = -t371 * t562 - t839;
t368 = t560 * t486;
t86 = -pkin(5) * t368 - qJ(6) * t366 + t113;
t933 = t86 / 0.2e1;
t258 = t562 * t438;
t289 = pkin(4) * t440 + pkin(10) * t438;
t276 = t562 * t289;
t828 = t559 * t209;
t677 = t276 + t828;
t107 = pkin(5) * t440 + qJ(6) * t258 + t677;
t932 = t107 / 0.2e1;
t132 = pkin(5) * t328 + t185;
t931 = t132 / 0.2e1;
t930 = t276 / 0.2e1;
t311 = t330 * t559;
t709 = -t311 / 0.2e1;
t929 = -t330 / 0.2e1;
t928 = -t366 / 0.2e1;
t926 = t389 / 0.2e1;
t735 = t562 * t813;
t434 = t504 * t559 + t735;
t825 = t559 * t560;
t399 = -qJ(6) * t825 + t434;
t924 = -t399 / 0.2e1;
t923 = t399 / 0.2e1;
t515 = t562 * t526;
t517 = t545 * t825;
t808 = t517 + t515;
t414 = pkin(5) * t560 - qJ(6) * t818 + t808;
t922 = t414 / 0.2e1;
t921 = -t438 / 0.2e1;
t920 = t438 / 0.2e1;
t919 = t440 / 0.2e1;
t918 = -t486 / 0.2e1;
t896 = t559 * pkin(5);
t684 = t545 + t896;
t494 = t684 * t560;
t917 = t494 / 0.2e1;
t514 = t559 * t526;
t916 = t514 / 0.2e1;
t891 = -qJ(6) - pkin(10);
t523 = t891 * t559;
t915 = t523 / 0.2e1;
t524 = t891 * t562;
t914 = t524 / 0.2e1;
t913 = -t524 / 0.2e1;
t912 = -t545 / 0.2e1;
t550 = -pkin(5) * t562 - pkin(4);
t911 = t550 / 0.2e1;
t553 = t559 ^ 2;
t910 = t553 / 0.2e1;
t555 = t562 ^ 2;
t909 = -t555 / 0.2e1;
t908 = -t559 / 0.2e1;
t906 = t560 / 0.2e1;
t905 = t562 / 0.2e1;
t904 = -t563 / 0.2e1;
t902 = pkin(5) * t330;
t901 = pkin(5) * t366;
t480 = t488 * t562;
t365 = -t371 * t559 + t480;
t900 = t365 * pkin(5);
t897 = t488 * pkin(4);
t894 = t563 * pkin(5);
t699 = -t821 / 0.2e1;
t704 = -t825 / 0.2e1;
t890 = t64 * t699 + t81 * t704;
t889 = t64 - t80;
t888 = pkin(2) * qJD(2);
t887 = pkin(5) * qJD(5);
t886 = pkin(5) * qJD(6);
t885 = t64 * t562;
t884 = t81 * t559;
t883 = t81 * t562;
t16 = t132 * t902 - t81 * t889;
t880 = qJD(1) * t16;
t17 = t889 * t328;
t879 = qJD(1) * t17;
t622 = t365 * t926 + t366 * t924;
t193 = t562 * t198;
t223 = t559 * t229;
t114 = t193 + t223;
t92 = -qJ(6) * t365 + t114;
t628 = t86 * t908 + t905 * t92;
t326 = t560 * t340;
t816 = t563 * t354;
t676 = -t326 + t816;
t197 = -t676 + t897;
t150 = t197 + t900;
t712 = t150 * t904;
t30 = t712 + (t486 * t917 + t628) * t560 + t622;
t878 = qJD(1) * t30;
t38 = -t328 * t81 - t330 * t64;
t877 = qJD(1) * t38;
t43 = -t185 * t328 + t438 * t97;
t876 = qJD(1) * t43;
t44 = t185 * t330 - t438 * t98;
t875 = qJD(1) * t44;
t874 = t107 * t559;
t254 = t559 * t438;
t204 = t562 * t209;
t275 = t559 * t289;
t810 = t204 - t275;
t121 = qJ(6) * t254 - t810;
t873 = t121 * t562;
t15 = t132 * t150 + t64 * t86 + t81 * t92;
t872 = t15 * qJD(1);
t163 = -pkin(5) * t254 + t210;
t18 = t107 * t64 + t121 * t81 + t132 * t163;
t871 = t18 * qJD(1);
t870 = t185 * t562;
t19 = -t328 * t92 - t330 * t86 - t365 * t81 - t366 * t64;
t869 = t19 * qJD(1);
t868 = t197 * t559;
t867 = t197 * t562;
t22 = -t107 * t330 - t121 * t328 + (t884 + t885) * t438;
t866 = t22 * qJD(1);
t29 = -t132 * t368 - t365 * t64 + t366 * t81;
t865 = t29 * qJD(1);
t864 = t290 * t560;
t863 = t290 * t563;
t862 = t328 * t366;
t861 = t328 * t563;
t860 = t330 * t365;
t859 = t330 * t562;
t858 = t330 * t563;
t34 = t113 * t438 + t185 * t365 + t197 * t328 + t368 * t97;
t857 = t34 * qJD(1);
t35 = -t114 * t438 + t185 * t366 + t197 * t330 + t368 * t98;
t856 = t35 * qJD(1);
t36 = t210 * t328 - t97 * t440 + (t276 + (-t185 + t209) * t559) * t438;
t855 = t36 * qJD(1);
t854 = t365 * t438;
t853 = t365 * t562;
t852 = t366 * t438;
t851 = t366 * t559;
t850 = t366 * t562;
t37 = t210 * t330 - t98 * t440 + (t810 - t870) * t438;
t849 = t37 * qJD(1);
t848 = t434 * t440;
t847 = t440 * t488;
t846 = t440 * t560;
t45 = t209 * t488 + t339 * t438 + (t676 - t864) * t486;
t844 = t45 * qJD(1);
t46 = t210 * t488 + t339 * t440 + (-t809 - t863) * t486;
t843 = t46 * qJD(1);
t841 = t485 * t560;
t840 = t488 * t438;
t836 = t545 * t328;
t835 = t552 * t564;
t834 = t553 * t563;
t554 = t560 ^ 2;
t833 = t554 * t559;
t827 = t559 * t328;
t826 = t559 * t365;
t824 = t559 * t563;
t823 = t560 * t328;
t820 = t562 * t328;
t819 = t562 * t554;
t815 = t563 * t365;
t814 = t563 * t438;
t88 = -(-t319 + t339) * t488 + (t318 - t340) * t486;
t812 = t88 * qJD(1);
t732 = t545 * t821;
t657 = t514 - t732;
t430 = -qJ(6) * t824 + t657;
t495 = t684 * t563;
t707 = t389 * t908;
t124 = (t399 * t905 + t707 - t495 / 0.2e1) * t563 + (t414 * t908 + t430 * t905 + t917) * t560;
t665 = t927 - t894 / 0.2e1;
t142 = (t925 + t665) * t821;
t811 = t124 * qJD(4) + t142 * qJD(5);
t536 = t555 - t553;
t535 = t555 + t553;
t556 = t563 ^ 2;
t537 = t556 - t554;
t115 = t209 * t486 - t290 * t438;
t807 = qJD(1) * t115;
t116 = -t210 * t486 + t290 * t440;
t806 = qJD(1) * t116;
t127 = t860 - t862;
t805 = qJD(1) * t127;
t738 = t328 * t368;
t159 = -t738 - t854;
t804 = qJD(1) * t159;
t737 = t330 * t368;
t162 = -t737 - t852;
t803 = qJD(1) * t162;
t164 = t318 * t488 - t319 * t486;
t802 = qJD(1) * t164;
t694 = t818 / 0.2e1;
t589 = t330 * t694 + t555 * t706;
t166 = -t851 / 0.2e1 + t589;
t801 = qJD(1) * t166;
t465 = t368 / 0.2e1;
t697 = -t820 / 0.2e1;
t613 = t697 + t311 / 0.2e1;
t172 = t560 * t613 + t465;
t800 = qJD(1) * t172;
t173 = t328 * t440 - t559 * t935;
t799 = qJD(1) * t173;
t174 = t330 * t440 - t562 * t935;
t798 = qJD(1) * t174;
t667 = t438 * t704;
t710 = t861 / 0.2e1;
t603 = t710 + t667;
t178 = t603 + t938;
t797 = qJD(1) * t178;
t696 = -t258 / 0.2e1;
t666 = t560 * t696;
t633 = t480 / 0.2e1 + t666;
t705 = -t842 / 0.2e1;
t183 = (t330 / 0.2e1 + t705) * t563 + t633;
t796 = qJD(1) * t183;
t237 = t840 - t841;
t795 = qJD(1) * t237;
t238 = -t840 - t841;
t794 = qJD(1) * t238;
t239 = -t847 - t943;
t793 = qJD(1) * t239;
t250 = -t847 + t943;
t792 = qJD(1) * t250;
t791 = qJD(1) * t330;
t790 = qJD(1) * t486;
t789 = qJD(1) * t563;
t788 = qJD(2) * t486;
t787 = qJD(3) * t563;
t786 = qJD(4) * t438;
t785 = qJD(4) * t486;
t784 = qJD(4) * t559;
t783 = qJD(4) * t560;
t782 = qJD(4) * t562;
t781 = qJD(4) * t563;
t780 = qJD(5) * t328;
t779 = qJD(5) * t438;
t778 = qJD(5) * t559;
t777 = qJD(5) * t562;
t776 = qJD(5) * t563;
t112 = pkin(2) * t727 * t942 - t318 * t339 + t319 * t340;
t775 = t112 * qJD(1);
t341 = -t826 / 0.2e1;
t624 = (-t861 / 0.2e1 + t928) * t562;
t702 = t824 / 0.2e1;
t120 = t330 * t702 + t341 + t624;
t774 = t120 * qJD(1);
t128 = -t860 - t862;
t773 = t128 * qJD(1);
t449 = t486 * t699;
t693 = -t814 / 0.2e1;
t599 = (t693 - t846 / 0.2e1) * t559;
t567 = t438 * t702 + t599 + t823 / 0.2e1;
t137 = t449 + t567;
t772 = t137 * qJD(1);
t153 = (t820 + t311) * t438;
t771 = t153 * qJD(1);
t157 = (t851 - t853) * t560;
t770 = t157 * qJD(1);
t160 = t738 - t854;
t769 = t160 * qJD(1);
t161 = -t737 + t852;
t768 = t161 * qJD(1);
t703 = t825 / 0.2e1;
t602 = t438 * t703 + t710;
t179 = t602 + t939;
t767 = t179 * qJD(1);
t654 = t929 + t705;
t181 = t563 * t654 + t633;
t766 = t181 * qJD(1);
t461 = t486 * t833;
t196 = t461 - t815;
t765 = t196 * qJD(1);
t212 = (t845 / 0.2e1 + t654) * t560;
t764 = t212 * qJD(1);
t651 = t814 + t846;
t220 = t651 * t486;
t763 = t220 * qJD(1);
t348 = t366 * t563;
t462 = t486 * t819;
t246 = t462 - t348;
t762 = t246 * qJD(1);
t761 = t254 * qJD(1);
t760 = t258 * qJD(1);
t601 = t557 * t918 + t881 * t488 / 0.2e1;
t321 = (-t831 / 0.2e1 + t601) * pkin(2);
t759 = t321 * qJD(1);
t758 = t368 * qJD(1);
t757 = t371 * qJD(1);
t756 = t838 * qJD(1);
t393 = t488 ^ 2 + t485;
t755 = t393 * qJD(1);
t593 = -pkin(8) * t830 - t669;
t446 = pkin(1) * t942 - t593 * t882;
t753 = t446 * qJD(1);
t503 = pkin(8) * t831 - t543;
t447 = pkin(1) * t835 - t503 * t882;
t752 = t447 * qJD(1);
t484 = (t681 / 0.2e1 + t832 / 0.2e1) * t558;
t751 = t484 * qJD(1);
t509 = (-t561 ^ 2 + t564 ^ 2) * t552;
t750 = t509 * qJD(1);
t749 = t560 * qJD(2);
t748 = t560 * qJD(5);
t746 = pkin(5) * t821;
t744 = t902 / 0.2e1;
t743 = -t899 / 0.2e1;
t742 = t899 / 0.2e1;
t741 = t896 / 0.2e1;
t740 = -t80 / 0.2e1 + t934;
t734 = t561 * t835;
t733 = t559 * t821;
t731 = t562 * t814;
t730 = t64 * t908;
t729 = t80 * t905;
t728 = t884 / 0.2e1;
t725 = t330 * t754;
t723 = t562 * t749;
t722 = t559 * t782;
t721 = t560 * t782;
t720 = t559 * t748;
t719 = t559 * t776;
t718 = t562 * t748;
t717 = t562 * t776;
t716 = t486 * t789;
t715 = t488 * t789;
t714 = t559 * t777;
t713 = t560 * t781;
t711 = t185 * t559 / 0.2e1;
t303 = t330 * t906;
t708 = -t858 / 0.2e1;
t701 = -t823 / 0.2e1;
t700 = -t368 / 0.2e1;
t698 = t821 / 0.2e1;
t691 = -t813 / 0.2e1;
t690 = t813 / 0.2e1;
t689 = -t193 / 0.2e1 - t223 / 0.2e1;
t688 = t204 / 0.2e1 - t275 / 0.2e1;
t686 = -t327 / 0.2e1 - t347 / 0.2e1;
t680 = t882 * qJD(1);
t679 = -0.2e1 * t733;
t678 = 0.2e1 * t733;
t674 = pkin(5) * t718;
t673 = pkin(5) * t698;
t672 = -qJD(4) - t790;
t671 = -qJD(5) - t754;
t670 = -qJD(5) + t747;
t668 = qJD(1) * t734;
t450 = t486 * t698;
t664 = t754 + qJD(5) / 0.2e1;
t663 = t898 / 0.2e1 - t210 / 0.2e1;
t662 = t558 * t680;
t661 = qJD(2) * t558 * t882;
t660 = t563 * t679;
t656 = t224 / 0.2e1 - t829 / 0.2e1;
t655 = -t326 / 0.2e1 + t816 / 0.2e1;
t653 = t559 * t700 + t303;
t155 = (t389 * t563 + t414 * t560) * t562 + (t399 * t563 + t430 * t560) * t559;
t579 = t121 * t907 + t399 * t920 + t81 * t904;
t580 = t107 * t907 + t389 * t920 + t64 * t904;
t617 = t365 * t913 + t366 * t915;
t619 = t414 * t929 - t430 * t328 / 0.2e1;
t4 = (-t92 / 0.2e1 + t580) * t562 + (t933 + t579) * t559 + t617 + t619;
t652 = t4 * qJD(1) - t155 * qJD(2);
t650 = -t486 * t546 + t488 * t545;
t455 = -t523 * t559 - t524 * t562;
t649 = t680 + qJD(2);
t569 = -t328 * t685 + t740 * t825;
t10 = t901 / 0.2e1 + t569;
t152 = t940 * t825;
t648 = -qJD(1) * t10 + qJD(2) * t152;
t125 = (t894 / 0.2e1 + t685) * t562;
t13 = (t742 + t740) * t562;
t647 = qJD(1) * t13 - qJD(2) * t125;
t577 = (t729 + t728) * t560 + t890;
t21 = (t708 + t365 / 0.2e1) * pkin(5) + t577;
t646 = -qJD(1) * t21 - qJD(2) * t142;
t597 = t433 * t919 + t808 * t921;
t631 = -pkin(4) * t365 / 0.2e1 - t867 / 0.2e1;
t23 = (t930 - t836 / 0.2e1 + t941) * t563 + (t97 / 0.2e1 + (t545 * t920 + t663) * t559) * t560 + t597 + t631;
t244 = t433 * t560 + (-t517 + t515) * t563;
t645 = -t23 * qJD(1) - t244 * qJD(2);
t236 = (t389 * t562 + t399 * t559) * t560;
t571 = t900 / 0.2e1 + t897 / 0.2e1 - t655;
t620 = t328 * t923 + t330 * t926;
t27 = (t728 + t885 / 0.2e1) * t560 + t571 + t620;
t644 = -qJD(1) * t27 - qJD(2) * t236;
t616 = t330 * t912 - t870 / 0.2e1;
t630 = pkin(4) * t928 + t868 / 0.2e1;
t24 = t438 * t916 + t848 / 0.2e1 + (t616 + t688) * t563 + (t98 / 0.2e1 + t663 * t562) * t560 + t630;
t245 = t514 * t563 + (-t434 + t735) * t560;
t643 = -t24 * qJD(1) + t245 * qJD(2);
t374 = -t433 * t563 - t545 * t833;
t627 = t433 * t921 + t903 * t97;
t40 = (t836 / 0.2e1 + t711) * t560 + t627 + t689;
t642 = qJD(1) * t40 - qJD(2) * t374;
t375 = -t434 * t563 - t545 * t819;
t626 = t434 * t920 + t904 * t98;
t39 = t560 * t616 + t626 + t656;
t641 = qJD(1) * t39 + qJD(2) * t375;
t640 = t670 * t560;
t568 = t559 * t693 + t599 + t701;
t136 = t450 + t568;
t512 = t537 * t559;
t639 = qJD(1) * t136 + qJD(2) * t512;
t139 = -t731 + (-t845 / 0.2e1 + t654) * t560;
t513 = t556 * t562 - t819;
t638 = -qJD(1) * t139 - qJD(2) * t513;
t208 = (-t827 - t859) * t560;
t510 = t535 * t554;
t637 = qJD(1) * t208 - qJD(2) * t510;
t240 = -t440 ^ 2 + t935;
t636 = qJD(1) * t240 - qJD(2) * t651;
t635 = -qJD(1) * t651 + qJD(2) * t537;
t216 = t311 - t820;
t634 = qJD(1) * t216 + qJD(4) * t535;
t632 = t893 / 0.2e1 - t895 / 0.2e1;
t431 = -t837 / 0.2e1 + t683 / 0.2e1;
t629 = t431 * qJD(1) + t749 / 0.2e1;
t625 = t559 * t692 - t480 / 0.2e1 + t666;
t623 = t864 / 0.2e1 + t546 * t919;
t621 = t365 * t915 + t366 * t914;
t615 = t822 / 0.2e1 + t817 / 0.2e1;
t614 = t697 + t709;
t572 = t863 / 0.2e1 + t546 * t921 + t545 * t465;
t108 = t572 - t686;
t612 = -qJD(1) * t108 - t546 * t747;
t110 = t326 / 0.2e1 + (t486 * t912 - t354 / 0.2e1) * t563 + t623;
t611 = -qJD(1) * t110 - t546 * t749;
t609 = -qJD(4) * t484 + t488 * t790;
t608 = -t723 - t784;
t607 = pkin(4) * t929 + pkin(10) * t696;
t606 = -t748 / 0.2e1 + t675;
t605 = t632 * t559;
t604 = t632 * t562;
t353 = t550 * t896;
t7 = -t740 * t524 + (t132 * t908 + t550 * t929 + t932) * pkin(5);
t99 = t685 * t524 + (t494 * t908 + t550 * t699 + t922) * pkin(5);
t600 = -qJD(1) * t7 - qJD(2) * t99 + qJD(4) * t353;
t499 = t524 * t704;
t598 = -t499 + (t523 * t907 + t923) * t562;
t596 = t873 / 0.2e1 - t874 / 0.2e1 + t931;
t565 = t107 * t926 + t121 * t923 + t495 * t931 + t163 * t917 + t64 * t922 + t81 * t430 / 0.2e1;
t584 = t150 * t911 + t86 * t915 + t913 * t92;
t1 = -t565 + t584;
t147 = t389 * t414 + t399 * t430 + t494 * t495;
t595 = -t1 * qJD(1) + t147 * qJD(2) + t124 * qJD(3);
t149 = t399 * t940 + t494 * t746;
t5 = -t944 + t740 * t399 + (t132 * t699 + t494 * t929 + t933) * pkin(5);
t594 = -qJD(1) * t5 + qJD(2) * t149 + qJD(3) * t142;
t585 = (t883 / 0.2e1 + t730 - t163 / 0.2e1) * t563;
t12 = t585 + (t486 * t911 + t596) * t560 + t621;
t474 = (-0.1e1 + t535) * t563 * t560;
t592 = t12 * qJD(1) + t124 * qJD(2) + t474 * qJD(3);
t145 = (-t827 + t859) * t560;
t151 = -t820 + 0.2e1 * t709;
t325 = t328 ^ 2;
t156 = t325 - t936;
t591 = qJD(1) * t156 - qJD(2) * t145 + qJD(4) * t151;
t217 = t325 + t936;
t590 = qJD(1) * t217 - qJD(2) * t208 + qJD(4) * t216;
t342 = t826 / 0.2e1;
t397 = t438 * t733;
t102 = t559 * t708 + t342 + t397 + t624;
t588 = t102 * qJD(1) + qJD(2) * t660;
t441 = t916 - t605;
t576 = pkin(4) * t328 / 0.2e1 + t870 / 0.2e1 + pkin(10) * t254 / 0.2e1;
t47 = t576 - t688;
t587 = pkin(4) * t782 - qJD(1) * t47 - qJD(2) * t441;
t442 = -t515 / 0.2e1 + t604;
t49 = -t276 / 0.2e1 - t941 + t607;
t586 = pkin(4) * t784 - qJD(1) * t49 - qJD(2) * t442;
t207 = t614 * t560;
t215 = -t827 / 0.2e1 + t859 / 0.2e1;
t583 = -qJD(2) * t207 - qJD(4) * t215 + t328 * t791;
t505 = (t910 + t909) * t560;
t582 = qJD(1) * t215 - qJD(2) * t505 + t722;
t581 = qJD(5) * t431 + t610;
t578 = t330 * t915 + t328 * t913 - t883 / 0.2e1;
t511 = t536 * t554;
t575 = qJD(1) * t145 + qJD(2) * t511 + qJD(4) * t678;
t574 = -qJD(1) * t151 + qJD(2) * t678 - qJD(4) * t536;
t573 = qJD(2) * t559 * t819 - qJD(1) * t207 + qJD(4) * t505;
t188 = t559 * t665 + t598 + t691;
t32 = (t743 + t934) * t559 + t578 + t615;
t492 = t553 * t907 + (0.1e1 / 0.2e1 + t909) * t560;
t570 = -qJD(1) * t32 + qJD(2) * t188 - qJD(3) * t492 + qJD(4) * t455;
t548 = t555 * t563;
t547 = t783 / 0.2e1;
t540 = t559 * t783;
t520 = qJD(5) * t679;
t501 = t505 * qJD(5);
t493 = t555 * t906 + (t910 + 0.1e1 / 0.2e1) * t560;
t475 = t484 * qJD(2);
t473 = t488 * t747;
t394 = -pkin(5) * t824 - t524 * t703 - t499;
t391 = t517 + t515 / 0.2e1 + t604;
t390 = t732 - t514 / 0.2e1 - t605;
t369 = t700 + t465;
t363 = t369 * qJD(4);
t362 = t368 * qJD(4);
t346 = t850 / 0.2e1;
t334 = -t758 - t783;
t320 = t745 / 0.2e1 + t601 * pkin(2);
t315 = qJD(2) * t700 + t431 * qJD(4);
t249 = (t608 - t791) * pkin(5);
t225 = t651 * qJD(4);
t214 = t216 * qJD(6);
t213 = t440 * t699 + t653;
t211 = t215 * qJD(5);
t205 = t208 * qJD(6);
t202 = t207 * qJD(5);
t187 = pkin(5) * t702 + t598 + t690 + t707;
t184 = t858 / 0.2e1 + t625;
t182 = t708 + t625;
t180 = t602 + t938;
t177 = t603 + t939;
t171 = (t918 + t613) * t560;
t165 = t851 / 0.2e1 + t589;
t148 = t151 * qJD(5);
t143 = t145 * qJD(5);
t140 = t440 * t698 + t653 + t731;
t138 = t450 + t567;
t135 = t449 + t568;
t126 = pkin(5) * t695 + t562 * t685;
t119 = t563 * t613 + t342 + t346;
t111 = t486 * t691 + t623 + t655;
t109 = t572 + t686;
t101 = t563 * t614 + t341 + t346 + t397;
t100 = pkin(5) * t922 + t389 * t914 + t398 * t913 + t494 * t741 + t550 * t673;
t50 = t711 + t828 / 0.2e1 + t930 + t607;
t48 = t576 + t688;
t42 = t185 * t698 + t303 * t545 - t626 + t656;
t41 = t185 * t704 + t545 * t701 - t627 + t689;
t33 = t559 * t743 - t578 + t615 + t730;
t31 = t494 * t700 + t560 * t628 - t622 + t712;
t28 = t571 - t620 + t890;
t26 = t657 * t921 - t848 / 0.2e1 + t810 * t904 + t98 * t907 + t330 * t690 + t545 * t666 + t210 * t698 + t185 * t694 + pkin(10) * t450 + t630;
t25 = pkin(10) * t465 * t559 + t185 * t702 + t210 * t703 + t328 * t690 + t545 * t667 + t677 * t904 + t907 * t97 - t597 + t631;
t20 = pkin(5) * t708 - t900 / 0.2e1 + t577;
t14 = t729 - t885 / 0.2e1 + t562 * t742;
t11 = t550 * t700 + t560 * t596 + t585 - t621;
t9 = -t901 / 0.2e1 + t569;
t8 = pkin(5) * t932 + t132 * t741 + t550 * t744 + t64 * t914 + t80 * t913;
t6 = pkin(5) * t933 + t132 * t673 + t494 * t744 + t64 * t924 + t80 * t923 + t944;
t3 = t559 * t579 + t562 * t580 - t617 + t619 + t628;
t2 = t565 + t584;
t51 = [0, 0, 0, qJD(2) * t734, t509 * qJD(2), t564 * t661, -t561 * t661, 0, -t446 * qJD(2), -t447 * qJD(2), qJD(2) * t88 + qJD(3) * t393, qJD(2) * t112 + qJD(3) * t164 (-t486 * t747 - t786) * t440, qJD(2) * t220 + qJD(4) * t240, qJD(2) * t239 - t438 * t785, -qJD(2) * t238 - t440 * t785, -t488 * t788, qJD(2) * t45 - qJD(3) * t237 + qJD(4) * t116, qJD(2) * t46 + qJD(3) * t250 + qJD(4) * t115 (qJD(2) * t366 - t438 * t782 - t780) * t330, qJD(2) * t128 + qJD(4) * t153 + qJD(5) * t156, qJD(2) * t161 + qJD(4) * t174 - t328 * t779, qJD(2) * t160 - qJD(4) * t173 - t330 * t779 (qJD(4) * t440 - t486 * t749) * t438, qJD(2) * t34 + qJD(3) * t159 + qJD(4) * t36 + qJD(5) * t44, qJD(2) * t35 + qJD(3) * t162 + qJD(4) * t37 + qJD(5) * t43, qJD(2) * t19 + qJD(3) * t127 + qJD(4) * t22 + qJD(5) * t17 + qJD(6) * t217, qJD(2) * t15 + qJD(3) * t29 + qJD(4) * t18 + qJD(5) * t16 + qJD(6) * t38; 0, 0, 0, t668, t750, t649 * t830, -t649 * t831, 0, qJD(2) * t593 - t753, qJD(2) * t503 - t752, t812 + (t486 * t881 + t488 * t557) * t888, t775 + (-t339 * t881 + t340 * t557) * t888 + t320 * qJD(3), t937 + (-qJD(1) * t440 - t749) * t371, -t537 * t788 - t225 + t763, -t488 * t749 + t793, t363 - t473 - t794, -t609, t844 + (-t339 * t563 + t560 * t650) * qJD(2) + t111 * qJD(4), t843 + (t339 * t560 + t563 * t650) * qJD(2) + t109 * qJD(4), qJD(4) * t165 + t202 + (t723 + t791) * t366, t773 + t101 * qJD(4) - t143 + (-t851 - t853) * t749, t768 + (-t462 - t348) * qJD(2) + t140 * qJD(4) + t177 * qJD(5), t769 + (t461 + t815) * qJD(2) + t135 * qJD(4) + t184 * qJD(5), -t937 + (-t664 + t747) * t368, -t113 * t747 + t857 + t25 * qJD(4) + t42 * qJD(5) + (t365 * t545 + t433 * t486 + t868) * t749, t114 * t747 + t856 + t26 * qJD(4) + t41 * qJD(5) + (t366 * t545 + t434 * t486 + t867) * t749, t869 + (-t365 * t399 - t366 * t389 + (-t559 * t92 - t562 * t86) * t560) * qJD(2) + t3 * qJD(4) + t9 * qJD(5) - t205, t872 + (t150 * t494 + t389 * t86 + t399 * t92) * qJD(2) + t31 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) + t28 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t755, qJD(2) * t320 + t802, 0, 0, 0, 0, 0, t363 - t795, t792, 0, 0, 0, 0, 0, qJD(4) * t138 + qJD(5) * t182 + t804, qJD(4) * t213 + qJD(5) * t180 + t803, qJD(4) * t119 + t805, t865 + t31 * qJD(2) + t11 * qJD(4) + t20 * qJD(5) + t171 * qJD(6) + (t371 + t826 + t850) * qJD(3) * t560; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t610, t636, t672 * t438, t369 * qJD(2) + t440 * t672, t475, qJD(2) * t111 + qJD(3) * t369 - qJD(4) * t210 + t806, qJD(2) * t109 + qJD(4) * t209 + t807, qJD(2) * t165 + t211 + (-t784 - t791) * t258, t101 * qJD(2) - t536 * t786 + t148 + t771, qJD(2) * t140 + t440 * t784 + t798, qJD(2) * t135 + t440 * t782 - t799, t581, t855 + t25 * qJD(2) + t138 * qJD(3) + (-t210 * t562 + t559 * t659) * qJD(4) + t50 * qJD(5), t849 + t26 * qJD(2) + t213 * qJD(3) + (t210 * t559 + t562 * t659) * qJD(4) + t48 * qJD(5), t866 + t3 * qJD(2) + t119 * qJD(3) + (-t874 + t873 + (t523 * t562 - t524 * t559) * t438) * qJD(4) + t14 * qJD(5) + t214, t871 + t2 * qJD(2) + t11 * qJD(3) + (t107 * t523 - t121 * t524 + t163 * t550) * qJD(4) + t8 * qJD(5) + t33 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t583, t591, qJD(2) * t177 + t328 * t671, qJD(2) * t184 + t330 * t671, t315, qJD(2) * t42 + qJD(3) * t182 + qJD(4) * t50 - qJD(5) * t98 + t875, qJD(2) * t41 + qJD(3) * t180 + qJD(4) * t48 + qJD(5) * t97 + t876, pkin(5) * t780 + qJD(2) * t9 + qJD(4) * t14 + t879, qJD(2) * t6 + qJD(3) * t20 + qJD(4) * t8 - t81 * t887 + t880; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t590, qJD(2) * t28 + qJD(3) * t171 + qJD(4) * t33 + t877; 0, 0, 0, -t668, -t750, -t564 * t662, t561 * t662, 0, t753, t752, -t812, qJD(3) * t321 - t775, t440 * t716 + t937, -t225 - t763, qJD(4) * t371 - t793, -t362 + t794, t609, qJD(4) * t110 + t488 * t787 - t844, -qJD(3) * t838 + qJD(4) * t108 - t843, qJD(4) * t166 - t366 * t791 + t202, qJD(4) * t102 - t143 - t773, -qJD(4) * t139 + qJD(5) * t178 - t768, qJD(4) * t136 + qJD(5) * t183 - t769, t368 * t664 - t937, -qJD(3) * t196 - qJD(4) * t23 - qJD(5) * t39 - t857, -qJD(3) * t246 - qJD(4) * t24 - qJD(5) * t40 - t856, -qJD(3) * t157 + qJD(4) * t4 + qJD(5) * t10 - t205 - t869, -qJD(3) * t30 - qJD(4) * t1 - qJD(5) * t5 - qJD(6) * t27 - t872; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t713, t537 * qJD(4), 0, 0, 0, t546 * t783, t546 * t781, -t554 * t714 + t555 * t713, qJD(4) * t660 - qJD(5) * t511, -qJD(4) * t513 + t560 * t719, qJD(4) * t512 + t560 * t717, -t713, -qJD(4) * t244 - qJD(5) * t375, qJD(4) * t245 + qJD(5) * t374, -qJD(4) * t155 - qJD(5) * t152 + qJD(6) * t510, qJD(4) * t147 + qJD(5) * t149 - qJD(6) * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t759, 0, 0, 0, 0, 0, t715, -t756, 0, 0, 0, 0, 0, -t765, -t762, -t770, t811 - t878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t675, t635, t757 + t781, t334, -t751, -t545 * t781 - t611, t545 * t783 - t612, t801 - t501 + (t555 * t749 + t722) * t563 (t548 - t834) * qJD(4) + t520 + t588, t540 + t638, t639 + t721, -t606 (t559 * t658 - t735) * qJD(4) + t391 * qJD(5) + t645 (t562 * t658 + t736) * qJD(4) + t390 * qJD(5) + t643 ((-t523 * t563 + t430) * t562 + (t524 * t563 - t414) * t559) * qJD(4) + t126 * qJD(5) + t652 (t414 * t523 - t430 * t524 + t495 * t550) * qJD(4) + t100 * qJD(5) + t187 * qJD(6) + t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t573, -t575, t559 * t640 + t797, t562 * t640 + t796, qJD(1) * t465 + t547, qJD(4) * t391 - qJD(5) * t434 - t641, qJD(4) * t390 + qJD(5) * t433 - t642, pkin(5) * t720 + qJD(4) * t126 - t648, qJD(4) * t100 - t399 * t887 + t594; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t637, qJD(4) * t187 + t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t755, -qJD(2) * t321 - t802, 0, 0, 0, 0, 0, -t362 - t473 + t795, qJD(2) * t838 - t486 * t781 - t792, 0, 0, 0, 0, 0, qJD(2) * t196 + qJD(4) * t137 + qJD(5) * t181 - t804, qJD(2) * t246 - qJD(4) * t212 + qJD(5) * t179 - t803, qJD(2) * t157 + qJD(4) * t120 - t805, qJD(2) * t30 + qJD(4) * t12 + qJD(5) * t21 + qJD(6) * t172 - t865; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t759, 0, 0, 0, 0, 0, -t715, t756, 0, 0, 0, 0, 0, t765, t762, t770, t811 + t878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t672 * t563, 0, 0, 0, 0, 0, -t719 - t721 + t772, t540 - t717 - t764, t774 + (t548 + t834) * qJD(4) (t455 * t563 + t550 * t560) * qJD(4) + t394 * qJD(5) + t493 * qJD(6) + t592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t559 * t781 - t718 + t766, -t562 * t781 + t720 + t767, 0, qJD(4) * t394 - t646 - t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t493 + t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t610, -t636, -qJD(2) * t371 + t486 * t754, qJD(2) * t368 + t440 * t790, t475, -qJD(2) * t110 + qJD(3) * t368 - t806, -qJD(2) * t108 + t486 * t787 - t807, -qJD(2) * t166 + t562 * t725 + t211, -qJD(2) * t102 + t148 - t771, qJD(2) * t139 + qJD(5) * t258 - t798, -qJD(2) * t136 - qJD(5) * t254 + t799, -t581, qJD(2) * t23 - qJD(3) * t137 + qJD(5) * t49 - t855, qJD(2) * t24 + qJD(3) * t212 + qJD(5) * t47 - t849, -qJD(2) * t4 - qJD(3) * t120 - qJD(5) * t13 + t214 - t866, qJD(2) * t1 - qJD(3) * t12 - qJD(5) * t7 - qJD(6) * t32 - t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t675, -t635, -t757, t758, t751, t611, t612, -t541 * t555 - t501 - t801, t520 - t588, -t638 - t717, -t639 + t719, t606, qJD(5) * t442 - t645, qJD(5) * t441 - t643, qJD(5) * t125 - t652, -qJD(5) * t99 + qJD(6) * t188 - t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t758, t716, 0, 0, 0, 0, 0, -t772, t764, -t774, -qJD(6) * t492 - t592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, t536 * qJD(5), 0, 0, 0, -pkin(4) * t778, -pkin(4) * t777, qJD(6) * t535, qJD(5) * t353 + qJD(6) * t455; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t582, -t574, -t562 * t670 + t760, t559 * t670 - t761, -t629, -pkin(10) * t777 - t586, pkin(10) * t778 - t587, -pkin(5) * t777 - t647, t524 * t887 + t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t634, t570; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583, -t591, -qJD(2) * t178 - qJD(4) * t258 + t328 * t754, -qJD(2) * t183 + qJD(4) * t254 + t725, t315, qJD(2) * t39 - qJD(3) * t181 - qJD(4) * t49 - t875, qJD(2) * t40 - qJD(3) * t179 - qJD(4) * t47 - t876, -qJD(2) * t10 + qJD(4) * t13 - t879, qJD(2) * t5 - qJD(3) * t21 + qJD(4) * t7 - t330 * t886 - t880; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t573, t575, -t797 + (-t559 * t749 + t782) * t563, t563 * t608 - t796, qJD(1) * t700 + t547, -qJD(4) * t442 + t641, -qJD(4) * t441 + t642, -qJD(4) * t125 + t648, qJD(4) * t99 - qJD(6) * t746 - t594; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t766, -t767, 0, t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t582, t574, t562 * t747 - t760, -t559 * t747 + t761, t629, t586, t587, t647, -t559 * t886 - t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t590, qJD(2) * t27 - qJD(3) * t172 + qJD(4) * t32 + t330 * t887 - t877; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t637, -qJD(4) * t188 - t644 + t674; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t492 - t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t634, pkin(5) * t778 - t570; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t51;
