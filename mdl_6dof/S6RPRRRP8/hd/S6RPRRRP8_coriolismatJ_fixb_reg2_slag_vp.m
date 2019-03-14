% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRP8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:47
% EndTime: 2019-03-09 06:25:10
% DurationCPUTime: 19.92s
% Computational Cost: add. (17079->831), mult. (30168->922), div. (0->0), fcn. (31631->6), ass. (0->595)
t607 = cos(qJ(5));
t605 = sin(qJ(4));
t606 = sin(qJ(3));
t590 = t605 * t606;
t608 = cos(qJ(4));
t609 = cos(qJ(3));
t856 = t608 * t609;
t560 = -t590 + t856;
t557 = t560 ^ 2;
t558 = t605 * t609 + t608 * t606;
t976 = t558 ^ 2;
t699 = t976 / 0.2e1 + t557 / 0.2e1;
t1024 = t607 * t699;
t604 = sin(qJ(5));
t1023 = t699 * t604;
t973 = -pkin(1) - pkin(7);
t778 = -pkin(8) + t973;
t561 = t778 * t606;
t713 = t778 * t609;
t714 = t605 * t561 - t608 * t713;
t1005 = t714 * t604;
t941 = t560 * pkin(5);
t1021 = -t1005 / 0.2e1 - t941 / 0.2e1;
t777 = qJD(3) + qJD(4);
t991 = t777 * t607;
t1006 = t604 * t991;
t1015 = t560 * t1006;
t603 = t607 ^ 2;
t535 = t603 * t557;
t414 = t535 + t976;
t805 = qJD(5) * t558;
t1022 = qJD(1) * t414 + t1015 + t805;
t376 = t607 * t560;
t374 = t604 * t560;
t527 = pkin(5) * t374;
t691 = -qJ(6) * t376 + t527;
t249 = t691 + t714;
t901 = t249 * t604;
t234 = t901 / 0.2e1;
t859 = t607 * qJ(6);
t947 = pkin(5) * t604;
t565 = -t859 + t947;
t871 = t565 * t604;
t734 = t871 / 0.2e1;
t1020 = t234 + (t734 - pkin(5) / 0.2e1) * t560 + t1021;
t554 = t608 * t561;
t689 = t605 * t713;
t989 = t554 + t689;
t1002 = t989 * t607;
t965 = -t1002 / 0.2e1;
t1004 = t714 * t607;
t1019 = -t1004 / 0.2e1;
t380 = t1005 / 0.2e1;
t551 = t560 * qJ(6);
t942 = t560 * pkin(4);
t943 = t558 * pkin(9);
t429 = t942 + t943;
t948 = pkin(3) * t609;
t392 = t429 + t948;
t339 = t604 * t392;
t845 = t339 / 0.2e1 + t1019;
t1017 = -t551 - t845;
t862 = t604 * qJ(6);
t938 = t607 * pkin(5);
t682 = t862 + t938;
t564 = -pkin(4) - t682;
t937 = t608 * pkin(3);
t550 = t564 - t937;
t512 = t550 * t604;
t556 = t564 * t604;
t837 = t512 / 0.2e1 + t556 / 0.2e1;
t870 = t565 * t607;
t1016 = t870 - t837;
t803 = qJD(5) * t607;
t641 = (0.1e1 / 0.2e1 + t699) * t607;
t984 = qJD(1) * t641;
t998 = -t374 * t777 - t558 * t803 - t984;
t1003 = t989 * t604;
t588 = pkin(3) * t606 + qJ(2);
t940 = t560 * pkin(9);
t945 = t558 * pkin(4);
t692 = -t940 + t945;
t633 = t588 + t692;
t224 = -t607 * t633 + t1003;
t225 = t604 * t633 + t1002;
t597 = t604 / 0.2e1;
t950 = t607 / 0.2e1;
t649 = t224 * t597 + t225 * t950;
t962 = -t989 / 0.2e1;
t1014 = (t962 + t649) * t560;
t882 = t558 * qJ(6);
t179 = t225 + t882;
t944 = t558 * pkin(5);
t180 = t224 - t944;
t652 = t179 * t950 + t180 * t597;
t370 = t604 * t558;
t372 = t607 * t558;
t990 = -pkin(5) * t370 + qJ(6) * t372 + t989;
t969 = -t990 / 0.2e1;
t1013 = (t969 + t652) * t560;
t1001 = t990 * t604;
t1012 = (-t180 + t1001) * t560;
t1011 = (-t224 + t1003) * t560;
t1010 = (-t225 + t1002) * t560;
t1000 = t990 * t607;
t1009 = t249 * t372 + (t179 - t1000) * t560;
t946 = t989 * pkin(4);
t1008 = t249 * t990;
t888 = t714 * t989;
t412 = t777 * t558;
t899 = t990 * t564;
t602 = t604 ^ 2;
t829 = t602 + t603;
t858 = t607 * t392;
t230 = t858 + t1005;
t188 = -t230 - t941;
t916 = t188 * t604;
t231 = -t1004 + t339;
t186 = t551 + t231;
t917 = t186 * t607;
t651 = t917 / 0.2e1 + t916 / 0.2e1;
t997 = -t899 / 0.2e1 - t651 * pkin(9);
t774 = t557 - t976;
t300 = t774 * t604;
t798 = t300 * qJD(1);
t126 = t560 * t991 - t798;
t996 = t777 * t714;
t766 = t937 / 0.2e1;
t960 = t550 / 0.2e1;
t696 = t766 + t960;
t958 = -t564 / 0.2e1;
t994 = t558 * (t958 + t696);
t939 = t605 * pkin(3);
t591 = pkin(9) + t939;
t764 = t591 / 0.2e1 - pkin(9) / 0.2e1;
t993 = t560 * t764;
t874 = t560 * t605;
t877 = t558 * t608;
t934 = qJD(3) * pkin(3);
t992 = (-t874 + t877) * t934;
t580 = t603 - t602;
t795 = t370 * qJD(1);
t804 = qJD(5) * t604;
t318 = t795 + t804;
t182 = -qJD(5) * t370 + t798;
t732 = t859 / 0.2e1;
t630 = (t732 - t947 / 0.2e1) * t937;
t957 = t564 / 0.2e1;
t721 = t957 + t960;
t247 = -t565 * t721 + t630;
t361 = t682 * t560;
t726 = -t224 / 0.2e1 + t180 / 0.2e1;
t621 = t726 * t607 + (-t179 / 0.2e1 + t225 / 0.2e1) * t604;
t956 = t565 / 0.2e1;
t740 = t249 * t956;
t613 = pkin(9) * t621 + t361 * t957 + t740;
t397 = t604 * t429;
t244 = -t1004 + t397;
t198 = t551 + t244;
t857 = t607 * t429;
t243 = t857 + t1005;
t199 = -t243 - t941;
t971 = t199 / 0.2e1;
t657 = -t198 * qJ(6) / 0.2e1 + pkin(5) * t971;
t6 = t613 + t657;
t961 = -t527 / 0.2e1;
t268 = t961 + (t732 + t956) * t560;
t800 = t268 * qJD(2);
t872 = t565 * t564;
t988 = -t6 * qJD(1) + t247 * qJD(3) - qJD(4) * t872 + t800;
t612 = t361 * t960 + t591 * t621 + t740;
t972 = -t186 / 0.2e1;
t658 = qJ(6) * t972 + t188 * pkin(5) / 0.2e1;
t4 = t612 + t658;
t884 = t550 * t565;
t987 = -t4 * qJD(1) - qJD(3) * t884 + t800;
t902 = t244 * t607;
t903 = t243 * t604;
t647 = t902 / 0.2e1 - t903 / 0.2e1;
t963 = t714 / 0.2e1;
t20 = t1014 + (t963 + t647) * t558;
t880 = t558 * t560;
t187 = (-0.1e1 + t829) * t880;
t801 = t187 * qJD(2);
t986 = t20 * qJD(1) + t801;
t904 = t231 * t607;
t905 = t230 * t604;
t648 = -t904 / 0.2e1 + t905 / 0.2e1;
t18 = t1014 + (t963 - t648) * t558;
t985 = t18 * qJD(1) + t801;
t817 = qJD(2) * t641;
t791 = t558 * qJD(1);
t982 = -t791 - qJD(5);
t864 = t602 * t560;
t525 = pkin(9) * t864;
t863 = t603 * t560;
t528 = pkin(9) * t863;
t981 = -t525 / 0.2e1 - t528 / 0.2e1;
t748 = t560 * t803;
t980 = t300 * t777 + t558 * t748;
t953 = t602 / 0.2e1;
t368 = (t953 - t603 / 0.2e1) * t560;
t752 = qJD(1) * t604 * t607;
t485 = t557 * t752;
t156 = t368 * t777 + t485;
t790 = t560 * qJD(1);
t745 = t558 * t790;
t427 = t607 * t745;
t979 = t370 * t777 + t427;
t750 = t558 * t804;
t978 = t376 * t777 - t750;
t704 = t560 * t752;
t325 = t580 * t777 - 0.2e1 * t704;
t754 = t607 * t791;
t977 = -qJD(5) * t641 - t754;
t114 = -0.2e1 * t604 * t376 * (-qJD(5) + t791) + t580 * t412;
t116 = t982 * t374;
t975 = 0.2e1 * t607 * t116;
t974 = pkin(4) / 0.2e1;
t970 = t249 / 0.2e1;
t968 = -t361 / 0.2e1;
t967 = -t392 / 0.2e1;
t966 = t1004 / 0.2e1;
t964 = -t429 / 0.2e1;
t727 = -t554 / 0.2e1;
t959 = -t560 / 0.2e1;
t592 = -pkin(4) - t937;
t955 = -t592 / 0.2e1;
t954 = t592 / 0.2e1;
t598 = -t604 / 0.2e1;
t952 = t605 / 0.2e1;
t951 = -t607 / 0.2e1;
t914 = t199 * t604;
t915 = t198 * t607;
t650 = t915 / 0.2e1 + t914 / 0.2e1;
t12 = t1013 + (t970 + t650) * t558;
t8 = t1013 + (t970 + t651) * t558;
t949 = t8 * qJD(3) + t12 * qJD(4);
t723 = -t714 / 0.2e1 + t963;
t724 = t989 / 0.2e1 + t962;
t86 = t558 * t723 + t560 * t724;
t936 = t86 * qJD(3);
t935 = pkin(3) * qJD(4);
t933 = t8 * qJD(1);
t656 = t862 / 0.2e1 + t938 / 0.2e1;
t907 = t225 * t604;
t909 = t224 * t607;
t918 = t180 * t607;
t920 = t179 * t604;
t30 = -t909 / 0.2e1 - t920 / 0.2e1 + t907 / 0.2e1 + t918 / 0.2e1 + t656 * t558;
t932 = t30 * qJD(5);
t29 = -t907 / 0.2e1 + (t882 / 0.2e1 + t179 / 0.2e1) * t604 + (t944 / 0.2e1 - t726) * t607;
t931 = -t29 * qJD(5) + t372 * qJD(6);
t37 = (t198 / 0.2e1 + t972) * t607 + (t971 - t188 / 0.2e1) * t604;
t930 = qJD(1) * t37;
t42 = (t244 / 0.2e1 - t231 / 0.2e1) * t607 + (-t243 / 0.2e1 + t230 / 0.2e1) * t604;
t929 = qJD(1) * t42;
t908 = t225 * t558;
t61 = -t908 + (t249 * t607 + t361 * t604) * t560;
t928 = qJD(1) * t61;
t910 = t224 * t558;
t62 = -t910 + (-t361 * t607 + t901) * t560;
t927 = qJD(1) * t62;
t921 = t179 * t558;
t75 = -t249 * t376 + t921;
t926 = qJD(1) * t75;
t76 = -t918 + t920;
t925 = qJD(1) * t76;
t674 = -t907 + t909;
t924 = qJD(1) * t674;
t923 = t12 * qJD(1);
t617 = t361 * t959 + t558 * t621;
t16 = t617 - t656;
t922 = t16 * qJD(1);
t21 = -t179 * t224 + t180 * t225 + t249 * t361;
t912 = t21 * qJD(1);
t145 = t180 * t372;
t22 = t145 - t188 * t376 + (t186 * t560 - t921) * t604;
t911 = t22 * qJD(1);
t23 = t145 - t199 * t376 + (t198 * t560 - t921) * t604;
t906 = t23 * qJD(1);
t25 = ((t179 - t225) * t607 + (t180 - t224) * t604) * t560;
t900 = t25 * qJD(1);
t26 = t29 * qJD(1);
t31 = t186 * t558 + t1009;
t894 = t31 * qJD(1);
t32 = t1012 + (-t188 - t901) * t558;
t893 = t32 * qJD(1);
t33 = t198 * t558 + t1009;
t892 = t33 * qJD(1);
t34 = t1012 + (-t199 - t901) * t558;
t891 = t34 * qJD(1);
t642 = t674 * t558;
t35 = (t230 * t607 + t231 * t604) * t560 + t642;
t890 = t35 * qJD(1);
t39 = (t243 * t607 + t244 * t604) * t560 + t642;
t889 = t39 * qJD(1);
t883 = t550 * t607;
t881 = t558 * t550;
t879 = t558 * t564;
t878 = t558 * t592;
t56 = t1011 + (t230 - t1005) * t558;
t876 = t56 * qJD(1);
t875 = t560 * t591;
t873 = t564 * t607;
t57 = t1010 + (-t231 - t1004) * t558;
t869 = t57 * qJD(1);
t58 = t1011 + (t243 - t1005) * t558;
t868 = t58 * qJD(1);
t59 = t1010 + (-t244 - t1004) * t558;
t867 = t59 * qJD(1);
t866 = t591 * t558;
t865 = t592 * t560;
t855 = t86 * qJD(1);
t553 = t856 / 0.2e1 - t590 / 0.2e1;
t835 = -t864 / 0.2e1 + t863 / 0.2e1;
t296 = -t553 + t835;
t852 = t296 * qJD(6);
t295 = t553 + t835;
t851 = t295 * qJD(6);
t843 = t858 / 0.2e1 + t380;
t841 = t397 / 0.2e1 + t1019;
t840 = -t397 / 0.2e1 + t966;
t379 = t857 / 0.2e1;
t839 = t379 + t380;
t472 = t591 * t864;
t473 = t591 * t863;
t838 = t472 + t473;
t834 = t525 + t528;
t719 = t829 * t608;
t552 = pkin(3) * t719;
t529 = t552 * qJD(4);
t793 = t552 * qJD(3);
t833 = t529 + t793;
t832 = t829 * pkin(9) * t937;
t810 = qJD(4) * t607;
t813 = qJD(3) * t607;
t831 = (t810 + t813) * t604;
t770 = t605 * t935;
t576 = t604 * t770;
t596 = t602 * qJD(6);
t830 = t596 - t576;
t110 = -t374 * t714 + t910;
t828 = qJD(1) * t110;
t111 = t376 * t714 - t908;
t827 = qJD(1) * t111;
t252 = t598 - t1023;
t826 = qJD(1) * t252;
t288 = t597 + t1023;
t824 = qJD(1) * t288;
t301 = t774 * t607;
t823 = qJD(1) * t301;
t336 = t558 * t948 + t588 * t560;
t822 = qJD(1) * t336;
t337 = -t588 * t558 + t560 * t948;
t821 = qJD(1) * t337;
t395 = t557 * t602 - t535;
t820 = qJD(1) * t395;
t411 = t727 + t554 / 0.2e1;
t819 = qJD(1) * t411;
t818 = qJD(1) * t588;
t291 = t951 + t1024;
t816 = qJD(2) * t291;
t815 = qJD(2) * t558;
t814 = qJD(3) * t604;
t812 = qJD(4) * t588;
t811 = qJD(4) * t604;
t809 = qJD(5) * t224;
t807 = qJD(5) * t376;
t806 = qJD(5) * t395;
t802 = qJD(6) * t604;
t797 = t774 * qJD(1);
t796 = t368 * qJD(5);
t346 = t372 * qJD(1);
t396 = t829 * t560;
t794 = t396 * qJD(1);
t792 = t553 * qJD(1);
t547 = t558 * qJD(6);
t579 = t606 ^ 2 - t609 ^ 2;
t789 = t579 * qJD(1);
t788 = t580 * qJD(5);
t787 = t606 * qJD(1);
t786 = t606 * qJD(3);
t785 = t607 * qJD(6);
t784 = t609 * qJD(1);
t783 = t609 * qJD(3);
t775 = t18 * qJD(3) + t20 * qJD(4);
t773 = t605 * t934;
t772 = pkin(9) * t804;
t771 = pkin(9) * t803;
t769 = t941 / 0.2e1;
t768 = pkin(9) * t598;
t767 = -t937 / 0.2e1;
t765 = t974 + t955;
t728 = t376 / 0.2e1;
t737 = t565 * t959;
t269 = qJ(6) * t728 + t737 + t961;
t763 = t269 * qJD(5) + t374 * qJD(6) + t801;
t235 = -t901 / 0.2e1;
t729 = -t376 / 0.2e1;
t730 = t372 / 0.2e1;
t762 = t550 * t729 + t591 * t730 + t235;
t761 = pkin(9) * t730 + t564 * t729 + t235;
t760 = -t268 * qJD(5) - t801;
t759 = qJD(3) * t973;
t758 = qJ(2) * t787;
t757 = qJ(2) * t784;
t756 = t588 * t791;
t755 = t604 * t791;
t753 = t588 * t790;
t751 = t604 * t815;
t747 = t591 * t804;
t746 = t591 * t803;
t586 = t604 * t803;
t744 = t560 * t802;
t583 = t604 * t785;
t743 = t606 * t783;
t736 = -t875 / 0.2e1;
t735 = -t874 / 0.2e1;
t731 = -t372 / 0.2e1;
t725 = t969 + t990 / 0.2e1;
t722 = t472 / 0.2e1 + t473 / 0.2e1;
t720 = pkin(3) * t777;
t639 = (t939 / 0.2e1 - t764) * t560;
t47 = t725 * t607 + (t639 - t994) * t604;
t575 = t607 * t773;
t718 = -qJD(1) * t47 + t575;
t708 = pkin(3) * t735;
t479 = t607 * t708;
t49 = t479 + t725 * t604 + (t993 + t994) * t607;
t574 = t604 * t773;
t717 = -qJD(1) * t49 + t574;
t478 = t728 * t939;
t697 = t767 + t955;
t638 = (-pkin(4) / 0.2e1 + t697) * t558;
t74 = t478 + t724 * t604 + (t638 - t993) * t607;
t716 = -qJD(1) * t74 - t574;
t71 = t965 + t1002 / 0.2e1 + (t639 + t638) * t604;
t715 = -qJD(1) * t71 + t575;
t413 = t777 * t560;
t711 = t777 * t604;
t709 = t607 * t770;
t707 = t560 * t768;
t703 = t557 * t586;
t698 = t968 - t943 / 0.2e1;
t694 = t591 * t719;
t690 = (t603 / 0.2e1 + t953) * t608;
t688 = t968 - t866 / 0.2e1;
t686 = t427 + t748;
t245 = qJD(1) * t295 + t831;
t285 = qJD(1) * t368 - t831;
t685 = t558 * t413;
t274 = t777 * t880;
t681 = t879 + t940;
t13 = t179 * t186 + t180 * t188 + t1008;
t680 = t13 * qJD(1) + t8 * qJD(2);
t14 = t179 * t198 + t180 * t199 + t1008;
t679 = t14 * qJD(1) + t12 * qJD(2);
t38 = -t224 * t230 + t225 * t231 + t888;
t678 = t38 * qJD(1) + t18 * qJD(2);
t40 = -t224 * t243 + t225 * t244 + t888;
t677 = t40 * qJD(1) + t20 * qJD(2);
t676 = t916 + t917;
t675 = t914 + t915;
t673 = t904 - t905;
t672 = t902 - t903;
t671 = t875 + t881;
t670 = -t875 - t878;
t112 = t588 * t948;
t668 = t112 * qJD(1) + t86 * qJD(2);
t407 = t871 + t883;
t653 = (t737 - t249 / 0.2e1) * t607;
t614 = t653 + (t560 * t960 + t688) * t604;
t53 = t614 + t1017;
t667 = -qJD(1) * t53 + qJD(3) * t407;
t408 = -t512 + t870;
t404 = t550 * t728;
t43 = t404 + (t967 + t688) * t607 + t1020;
t666 = -qJD(1) * t43 + qJD(3) * t408;
t665 = pkin(3) * t690;
t663 = qJD(5) * t565 - t802;
t662 = t583 - t709;
t64 = t769 + t762 + t843;
t646 = -qJD(1) * t64 + t550 * t814;
t645 = -t374 * qJD(5) - t558 * t991;
t618 = (t866 / 0.2e1 - t865 / 0.2e1) * t604 + t966;
t100 = t618 + t845;
t644 = -qJD(1) * t100 - t592 * t813;
t459 = t591 * t731;
t102 = t459 + (t865 / 0.2e1 + t967) * t607 + t723 * t604;
t643 = -qJD(1) * t102 - t592 * t814;
t305 = qJD(5) * t553 + t745;
t640 = t560 * t734 + t234 + 0.2e1 * t769;
t637 = -0.2e1 * t1015;
t636 = 0.2e1 * t1015;
t611 = t650 * t591 + (t249 * t952 + t608 * t652) * pkin(3) + t990 * t960;
t2 = t611 + t997;
t302 = (t550 * t605 + t694) * pkin(3);
t622 = t708 + t722 + t981;
t68 = (t960 + t958 + t665) * t558 + t622;
t635 = t2 * qJD(1) + t68 * qJD(2) + t302 * qJD(3);
t610 = t647 * t591 + (t608 * t649 + t714 * t952) * pkin(3) + t989 * t954;
t10 = t946 / 0.2e1 + t648 * pkin(9) + t610;
t359 = (t592 * t605 + t694) * pkin(3);
t83 = (t954 + t974 + t665) * t558 + t622;
t634 = t10 * qJD(1) + t83 * qJD(2) + t359 * qJD(3);
t571 = t604 * t767;
t277 = t571 + t1016;
t446 = -t556 + t870;
t437 = t564 * t728;
t50 = t437 + (t964 + t698) * t607 + t1020;
t632 = -qJD(1) * t50 + qJD(3) * t277 + qJD(4) * t446;
t572 = t607 * t766;
t278 = t607 * t721 + t572 + t871;
t445 = t871 + t873;
t615 = t653 + (t560 * t957 + t698) * t604;
t55 = -t551 + t615 + t840;
t631 = -qJD(1) * t55 + qJD(3) * t278 + qJD(4) * t445;
t620 = (t943 / 0.2e1 + t942 / 0.2e1) * t604 + t966;
t104 = t620 + t841;
t573 = t607 * t767;
t468 = t607 * t765 + t573;
t628 = pkin(4) * t810 - qJD(1) * t104 + qJD(3) * t468;
t501 = pkin(9) * t731;
t106 = t501 + (-t942 / 0.2e1 + t964) * t607;
t467 = t604 * t765 + t571;
t627 = pkin(4) * t811 - qJD(1) * t106 + qJD(3) * t467;
t570 = t604 * t766;
t321 = t570 + t837;
t66 = t769 + t761 + t839;
t625 = -qJD(1) * t66 + qJD(3) * t321 + t564 * t811;
t624 = -qJD(5) * t682 + t785;
t619 = (-t877 / 0.2e1 + t874 / 0.2e1) * pkin(3) + t736;
t616 = (t558 * t690 + t735) * pkin(3) + t722 - t981;
t601 = qJ(2) * qJD(2);
t600 = qJD(1) * qJ(2);
t587 = t606 * t784;
t486 = t560 * t583;
t484 = t982 * qJ(6);
t470 = pkin(4) * t951 + t592 * t950 + t573;
t469 = pkin(4) * t598 + t592 * t597 + t571;
t409 = t602 * t777 + t704;
t393 = t777 * t553;
t363 = t396 * qJD(2);
t322 = t570 - t837;
t320 = t346 + t803;
t304 = 0.2e1 * t727 - t689;
t289 = t597 - t1023;
t280 = t571 - t1016;
t279 = -t871 - t883 / 0.2e1 - t873 / 0.2e1 + t572;
t255 = t950 - t1024;
t254 = t598 + t1023;
t248 = t872 / 0.2e1 + t884 / 0.2e1 + t630;
t242 = t603 * t745 - t796;
t241 = t602 * t745 + t796;
t226 = t777 * t396;
t223 = t225 * qJD(5);
t213 = t636 - t820;
t212 = t637 + t820;
t197 = -t603 * t685 - t703;
t196 = -t602 * t685 + t703;
t184 = qJD(5) * t372 - t823;
t138 = -t372 * t777 + t604 * t745;
t137 = -t796 + (-t603 * t790 - t1006) * t558;
t136 = t796 + (-t602 * t790 + t1006) * t558;
t127 = t560 * t711 + t823;
t113 = t301 * t777 - t560 * t750;
t109 = -t412 * t604 + t807;
t107 = pkin(4) * t729 + t379 + 0.2e1 * t380 + t501;
t105 = t620 + t840;
t103 = t592 * t728 + t380 + t459 + t843;
t101 = t618 - t845;
t84 = t777 * t187;
t82 = t878 / 0.2e1 - t945 / 0.2e1 + t616;
t73 = t478 + t1003 + pkin(9) * t729 + pkin(4) * t730 + (t558 * t697 + t736) * t607;
t72 = 0.2e1 * t965 + t707 + t370 * t974 + (-t878 / 0.2e1 + t619) * t604;
t70 = t558 * t711 - t807;
t67 = t881 / 0.2e1 + t879 / 0.2e1 + t616;
t65 = -t857 / 0.2e1 + t761 + t1021;
t63 = -t858 / 0.2e1 + t762 + t1021;
t54 = t551 + t615 + t841;
t52 = t614 - t1017;
t51 = t607 * t698 + t437 + t640 + t839;
t48 = t479 - t1001 + pkin(9) * t728 + t564 * t730 + (t875 / 0.2e1 + t696 * t558) * t607;
t46 = -t1000 + t370 * t958 + t707 + (-t881 / 0.2e1 + t619) * t604;
t44 = t607 * t688 + t404 + t640 + t843;
t41 = t647 - t648;
t36 = t650 + t651;
t24 = t624 - t26;
t15 = t617 + t656;
t9 = pkin(9) * t904 / 0.2e1 + t230 * t768 - t946 / 0.2e1 + t610;
t5 = t613 - t657;
t3 = t612 - t658;
t1 = t611 - t997;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t601, -t743, t579 * qJD(3), 0, t743, 0, 0, qJ(2) * t783 + qJD(2) * t606, -qJ(2) * t786 + qJD(2) * t609, 0, t601, -t685, -t777 * t774, 0, t274, 0, 0, qJD(3) * t336 + t560 * t812 + t815, qJD(2) * t560 + qJD(3) * t337 - t558 * t812, 0, qJD(2) * t588 + qJD(3) * t112, t197, t558 * t636 + t806, t113, t196, -t980, t274, qJD(3) * t56 + qJD(4) * t58 + qJD(5) * t111 + t607 * t815, qJD(3) * t57 + qJD(4) * t59 + qJD(5) * t110 - t751, -qJD(3) * t35 - qJD(4) * t39 - t363, -qJD(2) * t674 + qJD(3) * t38 + qJD(4) * t40, t197, t113, t558 * t637 - t806, t274, t980, t196, qJD(3) * t32 + qJD(4) * t34 + qJD(5) * t61 + (-t557 * t802 + t815) * t607, -qJD(3) * t22 - qJD(4) * t23 - qJD(5) * t25 - t558 * t744 - t363, qJD(3) * t31 + qJD(4) * t33 + qJD(5) * t62 + qJD(6) * t414 + t751, qJD(2) * t76 + qJD(3) * t13 + qJD(4) * t14 + qJD(5) * t21 + qJD(6) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t600, 0, 0, 0, 0, 0, 0, t787, t784, 0, t600, 0, 0, 0, 0, 0, 0, t791, t790, 0, t818 + t936, 0, 0, 0, 0, 0, 0, -qJD(5) * t291 + t754, qJD(5) * t254 - t755, -t794, t775 - t924, 0, 0, 0, 0, 0, 0, qJD(5) * t255 + t754, -t794, qJD(5) * t289 + t755, qJD(5) * t15 + qJD(6) * t291 + t925 + t949; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t587, t789, -t786, t587, -t783, 0, -t606 * t759 + t757, -t609 * t759 - t758, 0, 0, -t745, -t797, -t412, t745, -t413, 0, -qJD(3) * t989 + qJD(4) * t304 + t822, t821 + t996, t992 (-t605 * t714 - t608 * t989) * t934 + t668, t137, -t114, t127, t136, t126, t305, t876 + (t604 * t670 - t1002) * qJD(3) + t72 * qJD(4) + t103 * qJD(5), t869 + (t607 * t670 + t1003) * qJD(3) + t73 * qJD(4) + t101 * qJD(5), qJD(3) * t673 + t41 * qJD(4) - t890 (t591 * t673 + t592 * t989) * qJD(3) + t9 * qJD(4) + t678, t137, t127, t114, t305, -t126, t136, t893 + (-t604 * t671 - t1000) * qJD(3) + t46 * qJD(4) + t44 * qJD(5) + t852, qJD(3) * t676 + t36 * qJD(4) - t911 + t932, t894 + (t607 * t671 - t1001) * qJD(3) + t48 * qJD(4) + t52 * qJD(5) + t486 (t550 * t990 + t591 * t676) * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + t63 * qJD(6) + t680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t745, -t797, -t412, t745, -t413, 0, qJD(3) * t304 - qJD(4) * t989 + t753, -t756 + t996, 0, 0, t137, -t114, t127, t136, t126, t305, t868 + t72 * qJD(3) + (t604 * t692 - t1002) * qJD(4) + t107 * qJD(5), t867 + t73 * qJD(3) + (t607 * t692 + t1003) * qJD(4) + t105 * qJD(5), t41 * qJD(3) + qJD(4) * t672 - t889, t9 * qJD(3) + (pkin(9) * t672 - t946) * qJD(4) + t677, t137, t127, t114, t305, -t126, t136, t891 + t46 * qJD(3) + (-t604 * t681 - t1000) * qJD(4) + t51 * qJD(5) + t852, t36 * qJD(3) + qJD(4) * t675 - t906 + t932, t892 + t48 * qJD(3) + (t607 * t681 - t1001) * qJD(4) + t54 * qJD(5) + t486, t1 * qJD(3) + (pkin(9) * t675 + t899) * qJD(4) + t5 * qJD(5) + t65 * qJD(6) + t679; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, t212, t116, t156, -t686, t393, qJD(3) * t103 + qJD(4) * t107 - t223 - t816 + t827, qJD(2) * t254 + qJD(3) * t101 + qJD(4) * t105 + t809 + t828, 0, 0, -t156, t116, t213, t393, t686, t156, qJD(2) * t255 + qJD(3) * t44 + qJD(4) * t51 - t223 + t928, qJD(5) * t691 + t30 * t777 - t744 - t900, qJD(2) * t289 + qJD(3) * t52 + qJD(4) * t54 + t547 - t809 + t927, t912 + t15 * qJD(2) + t3 * qJD(3) + t5 * qJD(4) + (-pkin(5) * t225 - qJ(6) * t224) * qJD(5) + t179 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296 * t777 - t485, t116, t1022, qJD(3) * t63 + qJD(4) * t65 + qJD(5) * t179 + t816 + t926; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t600, 0, 0, 0, 0, 0, 0, -t787, -t784, 0, -t600, 0, 0, 0, 0, 0, 0, -t791, -t790, 0, -t818 + t936, 0, 0, 0, 0, 0, 0, t977, -qJD(5) * t252 + t755, t794, t775 + t924, 0, 0, 0, 0, 0, 0, t977, t794, -qJD(5) * t288 - t755, qJD(5) * t16 + qJD(6) * t641 - t925 + t949; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t786, -t783, 0, 0, 0, 0, 0, 0, 0, 0, -t412, -t413, 0, t855 - t992, 0, 0, 0, 0, 0, 0, t645, t70, t226 (t838 + t878) * qJD(3) + t82 * qJD(4) + t985, 0, 0, 0, 0, 0, 0, t645, t226, t109, t933 + (t838 + t881) * qJD(3) + t67 * qJD(4) + t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t412, -t413, 0, 0, 0, 0, 0, 0, 0, 0, t645, t70, t226, t82 * qJD(3) + (t834 - t945) * qJD(4) + t986, 0, 0, 0, 0, 0, 0, t645, t226, t109, t923 + t67 * qJD(3) + (t834 + t879) * qJD(4) + t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t998, -t826 - t978, 0, 0, 0, 0, 0, 0, 0, 0, t998, 0, -t824 + t978, t269 * t777 + t558 * t624 + t922; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t587, -t789, 0, -t587, 0, 0, -t757, t758, 0, 0, t745, t797, 0, -t745, 0, 0, qJD(4) * t411 - t822, -t821, 0, -t668, t242, t975, t184, t241, t182, -t305, qJD(4) * t71 + qJD(5) * t102 - t876, qJD(4) * t74 + qJD(5) * t100 - t869, qJD(4) * t42 + t890, qJD(4) * t10 - t678, t242, t184, -t975, -t305, -t182, t241, qJD(4) * t47 + qJD(5) * t43 + t851 - t893, qJD(4) * t37 + t911 + t931, qJD(4) * t49 + qJD(5) * t53 + t486 - t894, qJD(4) * t2 + qJD(5) * t4 + qJD(6) * t64 - t680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t855, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t83 - t985, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t68 + t760 - t933; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t770, -t608 * t935, 0, 0, t586, t788, 0, -t586, 0, 0, t592 * t804 - t709, t592 * t803 + t576, t529, t359 * qJD(4), t586, 0, -t788, 0, 0, -t586, -qJD(5) * t408 + t662, t529, -qJD(5) * t407 + t830, qJD(4) * t302 + t550 * t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t605 * t720 + t819, -t608 * t720, 0, 0, t586, t788, 0, -t586, 0, 0, qJD(5) * t469 - t709 - t715, qJD(5) * t470 + t576 - t716, t833 + t929 (-pkin(4) * t939 + t832) * qJD(4) + t634, t586, 0, -t788, 0, 0, -t586, qJD(5) * t280 + t662 - t718, t833 + t930, qJD(5) * t279 - t717 + t830 (t564 * t939 + t832) * qJD(4) + t248 * qJD(5) + t322 * qJD(6) + t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t285, t325, t320, t285, -t318, -t792, qJD(4) * t469 - t643 - t746, qJD(4) * t470 - t644 + t747, 0, 0, -t285, t320, -t325, -t792, t318, t285, qJD(4) * t280 - t666 - t746, t24, qJD(4) * t279 - t667 - t747, t248 * qJD(4) + t591 * t624 - t987; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, t320, t409, qJD(4) * t322 - t646 + t746; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t745, t797, 0, -t745, 0, 0, -qJD(3) * t411 - t753, t756, 0, 0, t242, t975, t184, t241, t182, -t305, -qJD(3) * t71 + qJD(5) * t106 - t868, -qJD(3) * t74 + qJD(5) * t104 - t867, -qJD(3) * t42 + t889, -qJD(3) * t10 - t677, t242, t184, -t975, -t305, -t182, t241, -qJD(3) * t47 + qJD(5) * t50 + t851 - t891, -qJD(3) * t37 + t906 + t931, -qJD(3) * t49 + qJD(5) * t55 + t486 - t892, -qJD(3) * t2 + qJD(5) * t6 + qJD(6) * t66 - t679; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t83 - t986, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t68 + t760 - t923; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773 - t819, t608 * t934, 0, 0, t586, t788, 0, -t586, 0, 0, -qJD(5) * t467 + t715, -qJD(5) * t468 + t716, -t793 - t929, -t634, t586, 0, -t788, 0, 0, -t586, -qJD(5) * t277 + t583 + t718, -t793 - t930, -qJD(5) * t278 + t596 + t717, -qJD(5) * t247 - qJD(6) * t321 - t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t586, t788, 0, -t586, 0, 0, -pkin(4) * t804, -pkin(4) * t803, 0, 0, t586, 0, -t788, 0, 0, -t586, -qJD(5) * t446 + t583, 0, -qJD(5) * t445 + t596, t663 * t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t285, t325, t320, t285, -t318, -t792, -t627 - t771, -t628 + t772, 0, 0, -t285, t320, -t325, -t792, t318, t285, -t632 - t771, t24, -t631 - t772, pkin(9) * t624 - t988; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, t320, t409, -t625 + t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t213, t138, -t156, t979, t393, -qJD(3) * t102 - qJD(4) * t106 + t817 - t827, qJD(2) * t252 - qJD(3) * t100 - qJD(4) * t104 - t828, 0, 0, t156, t138, t212, t393, -t979, -t156, -qJD(3) * t43 - qJD(4) * t50 + t817 - t928, t29 * t777 + t900, qJD(2) * t288 - qJD(3) * t53 - qJD(4) * t55 + t547 - t927, qJ(6) * t547 - qJD(2) * t16 - qJD(3) * t4 - qJD(4) * t6 - t912; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t984, t826, 0, 0, 0, 0, 0, 0, 0, 0, t984, 0, t824, t268 * t777 - t922; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t325, -t346, -t285, t795, t792, qJD(4) * t467 + t643, qJD(4) * t468 + t644, 0, 0, t285, -t346, t325, t792, -t795, -t285, qJD(4) * t277 + t666, t26, qJD(4) * t278 + t667, qJD(4) * t247 + t987; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t325, -t346, -t285, t795, t792, t627, t628, 0, 0, t285, -t346, t325, t792, -t795, -t285, t632, t26, t631, t988; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t982, -t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295 * t777 + t485, t138, -t1022, -qJ(6) * t805 - qJD(3) * t64 - qJD(4) * t66 - t817 - t926; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t984; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, -t346, -t409, qJD(4) * t321 + t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, -t346, -t409, t625; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t982, t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;