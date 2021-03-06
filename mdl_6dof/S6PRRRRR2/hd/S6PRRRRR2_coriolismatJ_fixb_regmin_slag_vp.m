% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x32]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRRRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:27
% EndTime: 2019-03-09 00:46:03
% DurationCPUTime: 20.70s
% Computational Cost: add. (12385->717), mult. (28710->972), div. (0->0), fcn. (33476->12), ass. (0->606)
t781 = qJD(3) + qJD(4);
t566 = sin(qJ(3));
t921 = sin(qJ(4));
t761 = t921 * t566;
t570 = cos(qJ(3));
t922 = cos(qJ(4));
t762 = t922 * t570;
t521 = t761 - t762;
t931 = -t521 / 0.2e1;
t569 = cos(qJ(5));
t913 = t569 * pkin(5);
t558 = -pkin(4) - t913;
t1019 = t558 / 0.2e1;
t563 = sin(pkin(6));
t567 = sin(qJ(2));
t863 = t563 * t567;
t907 = cos(pkin(6));
t640 = t566 * t863 - t570 * t907;
t466 = t922 * t640;
t639 = t566 * t907 + t570 * t863;
t974 = t921 * t639;
t998 = t974 + t466;
t1014 = t998 / 0.2e1;
t565 = sin(qJ(5));
t568 = cos(qJ(6));
t838 = t568 * t565;
t564 = sin(qJ(6));
t852 = t564 * t569;
t672 = t838 + t852;
t364 = t672 * t521;
t1018 = t1014 * t364;
t551 = t921 * t570;
t554 = t922 * t566;
t504 = t554 / 0.2e1 + t551 / 0.2e1;
t523 = -t554 - t551;
t571 = cos(qJ(2));
t862 = t563 * t571;
t737 = t862 / 0.2e1;
t310 = -t504 * t862 + t523 * t737;
t999 = t922 * t639 - t921 * t640;
t803 = qJD(4) * t999;
t807 = qJD(3) * t999;
t1017 = t310 * qJD(2) - t803 - t807;
t837 = t568 * t569;
t853 = t564 * t565;
t958 = t837 - t853;
t363 = t958 * t523;
t1016 = -t363 / 0.2e1;
t930 = t521 / 0.2e1;
t775 = t922 * pkin(3);
t557 = -t775 - pkin(4);
t535 = t557 - t913;
t1015 = t535 / 0.2e1;
t1013 = t999 / 0.2e1;
t365 = t672 * t523;
t939 = pkin(10) + pkin(11);
t537 = t939 * t565;
t539 = t939 * t569;
t673 = -t568 * t537 - t564 * t539;
t1008 = t365 * t1019 + t673 * t931;
t632 = -t762 / 0.2e1 + t761 / 0.2e1;
t309 = t521 * t737 + t632 * t862;
t1012 = t309 * qJD(2) + t781 * t998;
t940 = pkin(8) + pkin(9);
t538 = t940 * t566;
t540 = t940 * t570;
t698 = t922 * t538 + t540 * t921;
t1000 = t698 * t569;
t749 = t1000 / 0.2e1;
t1011 = t365 * t1013;
t1010 = t363 * t1013;
t1001 = t698 * t565;
t390 = t569 * t521;
t501 = t523 * pkin(5);
t1009 = pkin(11) * t390 + t1001 - t501;
t780 = qJD(5) + qJD(6);
t387 = t565 * t521;
t776 = pkin(11) * t387;
t914 = t523 * pkin(4);
t915 = t521 * pkin(10);
t413 = -t914 + t915;
t397 = t565 * t413;
t824 = t1000 - t397;
t212 = t776 - t824;
t841 = t568 * t212;
t732 = -t841 / 0.2e1;
t398 = t569 * t413;
t171 = t398 + t1009;
t860 = t564 * t171;
t660 = -t860 / 0.2e1 + t732;
t1007 = t660 - t1008;
t431 = -t564 * t537 + t568 * t539;
t995 = t1016 * t558 + t431 * t931;
t774 = t921 * pkin(3);
t556 = t774 + pkin(10);
t912 = pkin(11) + t556;
t502 = t912 * t565;
t503 = t912 * t569;
t719 = t568 * t502 + t564 * t503;
t997 = t1015 * t365 + t719 * t930;
t383 = -t564 * t502 + t568 * t503;
t996 = t1016 * t535 + t383 * t931;
t929 = -t523 / 0.2e1;
t981 = t998 * t569;
t1004 = t981 / 0.2e1;
t975 = t781 * t523;
t1003 = t521 * t975;
t811 = qJD(2) * t523;
t755 = t521 * t811;
t312 = -t504 * qJD(5) + t755;
t932 = t672 / 0.2e1;
t933 = -t958 / 0.2e1;
t156 = t363 * t933 + t365 * t932;
t717 = t781 * t672;
t89 = qJD(2) * t156 + t717 * t958;
t919 = pkin(3) * t566;
t393 = t413 + t919;
t377 = t565 * t393;
t825 = -t377 + t1000;
t199 = t776 - t825;
t842 = t568 * t199;
t733 = -t842 / 0.2e1;
t378 = t569 * t393;
t170 = t378 + t1009;
t861 = t564 * t170;
t661 = -t861 / 0.2e1 + t733;
t994 = t661 - t997;
t857 = t564 * t199;
t736 = -t857 / 0.2e1;
t846 = t568 * t170;
t659 = t736 + t846 / 0.2e1;
t993 = t659 - t996;
t856 = t564 * t212;
t735 = -t856 / 0.2e1;
t845 = t568 * t171;
t658 = t735 + t845 / 0.2e1;
t992 = t658 - t995;
t559 = -pkin(3) * t570 - pkin(2);
t917 = t521 * pkin(4);
t700 = t523 * pkin(10) + t917;
t638 = t559 + t700;
t699 = -t921 * t538 + t540 * t922;
t978 = t699 * t565;
t215 = -t569 * t638 + t978;
t991 = (t215 - t978) * t523;
t977 = t699 * t569;
t216 = t565 * t638 + t977;
t990 = (t216 - t977) * t523;
t308 = (t931 + t632) * t862;
t789 = t308 * qJD(1);
t989 = t698 * t781 - t789;
t988 = t781 * t699;
t987 = t780 * t383;
t986 = t780 * t431;
t985 = t780 * t673;
t984 = 0.2e1 * t523;
t982 = t521 * t780;
t665 = -pkin(5) * t387 + t699;
t980 = t665 * t958;
t979 = t665 * t672;
t411 = t781 * t521;
t972 = 0.2e1 * t999;
t513 = t523 ^ 2;
t779 = -t521 ^ 2 + t513;
t877 = t523 * t565;
t328 = -pkin(5) * t877 + t698;
t895 = t328 * t672;
t750 = -t895 / 0.2e1;
t34 = t750 + t993;
t599 = t998 * t672;
t593 = -t599 / 0.2e1;
t655 = t852 / 0.2e1 + t838 / 0.2e1;
t636 = t655 * t998;
t116 = t593 + t636;
t793 = t116 * qJD(1);
t806 = qJD(3) * t535;
t971 = qJD(2) * t34 - t672 * t806 + t793;
t878 = t523 * t556;
t879 = t521 * t557;
t966 = t917 / 0.2e1 + t878 / 0.2e1 - t879 / 0.2e1 + (t921 * t929 + t922 * t931) * pkin(3);
t230 = qJD(6) * t504 - t312;
t307 = (t929 - t504) * t862;
t790 = t307 * qJD(1);
t810 = qJD(2) * t559;
t965 = t523 * t810 + t790;
t375 = t521 * t919 - t523 * t559;
t964 = -qJD(2) * t375 + t790;
t763 = t569 * t922;
t705 = -t763 / 0.2e1;
t764 = t565 * t922;
t709 = t564 * t764;
t610 = (t568 * t705 + t709 / 0.2e1) * pkin(3);
t725 = t1019 + t1015;
t227 = -t725 * t958 + t610;
t896 = t328 * t958;
t751 = -t896 / 0.2e1;
t41 = t751 + t1007;
t606 = t998 * t565;
t597 = t564 * t606;
t575 = -t597 / 0.2e1 + t568 * t1004;
t600 = t998 * t958;
t594 = -t600 / 0.2e1;
t117 = t594 + t575;
t792 = t117 * qJD(1);
t802 = qJD(4) * t558;
t963 = qJD(2) * t41 + qJD(3) * t227 - t802 * t958 + t792;
t918 = pkin(5) * t565;
t496 = t672 * t918;
t201 = -t496 + t227;
t924 = -t569 / 0.2e1;
t925 = t565 / 0.2e1;
t591 = (t363 * t925 + (t564 / 0.2e1 - t672 * t924) * t523) * pkin(5) + t751;
t21 = t591 + t1007;
t865 = t558 * t958;
t374 = t496 + t865;
t962 = t21 * qJD(2) + t201 * qJD(3) - t374 * qJD(4) + t792;
t495 = t958 * t918;
t870 = t535 * t672;
t338 = -t495 + t870;
t923 = t569 / 0.2e1;
t590 = (t365 * t925 + (-t568 / 0.2e1 - t958 * t923) * t523) * pkin(5) + t750;
t8 = t590 + t993;
t961 = t8 * qJD(2) - t338 * qJD(3) + t793;
t706 = -t764 / 0.2e1;
t611 = (t564 * t705 + t568 * t706) * pkin(3);
t226 = -t672 * t725 + t611;
t40 = t750 + t992;
t596 = t606 / 0.2e1;
t574 = t1004 * t564 + t568 * t596;
t115 = t593 + t574;
t794 = t115 * qJD(1);
t960 = qJD(2) * t40 + qJD(3) * t226 - t672 * t802 + t794;
t200 = t495 + t226;
t22 = t590 + t992;
t864 = t558 * t672;
t373 = -t495 + t864;
t959 = t22 * qJD(2) + t200 * qJD(3) - t373 * qJD(4) + t794;
t767 = t565 * t863;
t437 = t521 * t862;
t834 = t569 * t437;
t372 = t767 - t834;
t839 = t568 * t372;
t848 = t565 * t437;
t370 = t569 * t863 + t848;
t855 = t564 * t370;
t657 = -t855 / 0.2e1 - t839 / 0.2e1;
t851 = t565 * t999;
t271 = t569 * t862 + t851;
t836 = t569 * t999;
t272 = -t565 * t862 + t836;
t87 = t568 * t271 + t564 * t272;
t954 = -t998 * t365 / 0.2e1 + t87 * t931;
t576 = t657 + t954;
t957 = qJD(1) * t576;
t840 = t568 * t370;
t854 = t564 * t372;
t656 = -t854 / 0.2e1 + t840 / 0.2e1;
t88 = t564 * t271 - t568 * t272;
t952 = t1014 * t363 + t88 * t931;
t577 = t656 + t952;
t956 = qJD(1) * t577;
t367 = t958 * t521;
t955 = -t1014 * t367 + t88 * t929;
t953 = t87 * t929 + t1018;
t142 = -t363 ^ 2 + t365 ^ 2;
t97 = t363 * t672 + t365 * t958;
t11 = qJD(2) * t142 + t781 * t97;
t817 = qJD(2) * t363;
t46 = t156 * t781 - t365 * t817;
t348 = -t672 ^ 2 + t958 ^ 2;
t71 = qJD(2) * t97 + t348 * t781;
t876 = t523 * t569;
t195 = pkin(11) * t876 - t215;
t916 = t521 * pkin(5);
t165 = t195 + t916;
t196 = pkin(11) * t877 + t216;
t843 = t568 * t196;
t78 = t564 * t165 + t843;
t951 = -t328 * t367 - t665 * t363 + t78 * t523;
t164 = t568 * t165;
t858 = t564 * t196;
t77 = -t164 + t858;
t950 = -t328 * t364 - t665 * t365 + t77 * t523;
t562 = t569 ^ 2;
t941 = t565 ^ 2;
t386 = (-t941 / 0.2e1 + t562 / 0.2e1) * t523;
t847 = t565 * t569;
t757 = qJD(2) * t847;
t949 = t386 * t781 + t513 * t757;
t778 = -t562 + t941;
t350 = -0.2e1 * t523 * t757 + t778 * t781;
t654 = -t837 / 0.2e1 + t853 / 0.2e1;
t635 = t654 * t998;
t118 = t594 - t635;
t948 = qJD(2) * t576 + t118 * qJD(3) + t117 * qJD(4);
t947 = qJD(2) * t577 + t116 * qJD(3) + t115 * qJD(4);
t264 = t599 / 0.2e1;
t121 = t264 + t574;
t122 = t264 + t636;
t588 = t656 - t952;
t945 = qJD(2) * t588 + t122 * qJD(3) + t121 * qJD(4) + t780 * t88;
t263 = t600 / 0.2e1;
t119 = t263 + t575;
t120 = t263 - t635;
t589 = t657 - t954;
t944 = qJD(2) * t589 + t120 * qJD(3) + t119 * qJD(4) + t780 * t87;
t791 = t118 * qJD(1);
t943 = t719 * t780 - t791;
t938 = -t164 / 0.2e1;
t937 = -t165 / 0.2e1;
t928 = t523 / 0.2e1;
t927 = t556 / 0.2e1;
t926 = -t565 / 0.2e1;
t920 = pkin(3) * t521;
t911 = t780 * t97;
t910 = pkin(5) * qJD(5);
t909 = pkin(5) * qJD(6);
t769 = -t916 / 0.2e1;
t703 = t769 + t195 / 0.2e1;
t26 = (t937 + t703) * t564;
t906 = qJD(2) * t26;
t28 = t568 * t703 + t938;
t905 = qJD(2) * t28;
t98 = -t363 * t364 - t365 * t367;
t903 = qJD(2) * t98;
t902 = t271 * t521;
t901 = t272 * t521;
t898 = t328 * t363;
t897 = t328 * t365;
t891 = t719 * t523;
t889 = t383 * t523;
t887 = t673 * t523;
t882 = t431 * t523;
t436 = t523 * t862;
t881 = t436 * t565;
t880 = t436 * t569;
t874 = t535 * t364;
t872 = t535 * t367;
t871 = t535 * t958;
t868 = t558 * t364;
t866 = t558 * t367;
t859 = t564 * t195;
t844 = t568 * t195;
t826 = t780 * t156;
t157 = t364 * t521 - t365 * t523;
t821 = qJD(2) * t157;
t158 = t363 * t523 - t367 * t521;
t820 = qJD(2) * t158;
t305 = t779 * t565;
t819 = qJD(2) * t305;
t306 = t779 * t569;
t818 = qJD(2) * t306;
t376 = -t521 * t559 - t523 * t919;
t815 = qJD(2) * t376;
t814 = qJD(2) * t386;
t396 = t778 * t513;
t813 = qJD(2) * t396;
t812 = qJD(2) * t521;
t809 = qJD(2) * t567;
t808 = qJD(2) * t570;
t805 = qJD(3) * t565;
t804 = qJD(3) * t569;
t801 = qJD(4) * t559;
t800 = qJD(4) * t565;
t799 = qJD(4) * t569;
t798 = qJD(5) * t565;
t797 = qJD(5) * t569;
t796 = qJD(6) * t535;
t795 = qJD(6) * t558;
t222 = (t932 + t655) * t521;
t213 = t222 * qJD(2);
t223 = (t933 + t654) * t521;
t214 = t223 * qJD(2);
t788 = t779 * qJD(2);
t787 = t387 * qJD(2);
t381 = t390 * qJD(2);
t786 = t504 * qJD(2);
t546 = -t566 ^ 2 + t570 ^ 2;
t784 = t546 * qJD(2);
t783 = t566 * qJD(3);
t782 = t570 * qJD(3);
t777 = pkin(5) * t876;
t773 = pkin(2) * t566 * qJD(2);
t772 = pkin(2) * t808;
t771 = -t918 / 0.2e1;
t768 = t501 / 0.2e1;
t236 = t896 / 0.2e1;
t743 = t876 / 0.2e1;
t766 = -pkin(5) * t672 * t743 + t363 * t771 + t236;
t237 = t895 / 0.2e1;
t765 = t365 * t771 + t237 + t958 * t777 / 0.2e1;
t759 = t521 * t810;
t756 = qJD(5) * t521 * t523;
t754 = t563 * t809;
t753 = qJD(2) * t862;
t549 = t565 * t797;
t752 = t566 * t808;
t748 = t436 * t933;
t747 = t436 * t932;
t746 = t881 / 0.2e1;
t745 = -t880 / 0.2e1;
t744 = t958 * t931;
t738 = t863 / 0.2e1;
t734 = t848 / 0.2e1;
t731 = -t836 / 0.2e1;
t730 = t834 / 0.2e1;
t728 = -t377 / 0.2e1 + t749;
t727 = t397 / 0.2e1 - t1000 / 0.2e1;
t724 = t922 * qJD(3);
t723 = t922 * qJD(4);
t722 = t921 * qJD(3);
t721 = t921 * qJD(4);
t720 = pkin(5) * t780;
t715 = t781 * t569;
t714 = t780 * t672;
t713 = qJD(5) + t812;
t712 = pkin(3) * t721;
t711 = pkin(3) * t722;
t710 = -t774 / 0.2e1;
t702 = t768 - t170 / 0.2e1;
t701 = t768 - t171 / 0.2e1;
t697 = t565 * t715;
t695 = t781 * t847;
t572 = t1011 + (t564 * t981 + t568 * t606) * t931;
t583 = -t748 + t953;
t13 = t572 + t583;
t4 = (t845 - t856) * t521 + t950;
t694 = -t13 * qJD(1) + t4 * qJD(2);
t607 = t1011 - t1018;
t14 = t583 + t607;
t2 = (t846 - t857) * t521 + t950;
t693 = -t14 * qJD(1) + t2 * qJD(2);
t573 = t1010 + (-t568 * t981 + t597) * t930;
t582 = -t747 - t955;
t15 = t573 + t582;
t5 = -(t841 + t860) * t521 + t951;
t692 = -t15 * qJD(1) + t5 * qJD(2);
t608 = t744 * t998 + t1010;
t16 = t582 + t608;
t3 = -(t842 + t861) * t521 + t951;
t691 = -t16 * qJD(1) + t3 * qJD(2);
t688 = qJD(6) + t713;
t686 = -t748 - t953;
t685 = -t747 + t955;
t85 = -t843 - t859;
t38 = t365 * t777 + t521 * t85 - t898;
t684 = -qJD(2) * t38 + t956;
t86 = t844 - t858;
t39 = t363 * t777 - t521 * t86 + t897;
t683 = -qJD(2) * t39 + t957;
t58 = -t521 * t78 - t898;
t682 = -qJD(2) * t58 + t956;
t57 = t521 * t77 + t897;
t681 = -qJD(2) * t57 + t957;
t664 = (t271 / 0.2e1 - t851 / 0.2e1) * t523;
t61 = t745 + t664;
t75 = t398 * t521 + t991;
t680 = t61 * qJD(1) + t75 * qJD(2);
t584 = t999 * t929;
t581 = t271 * t928 + t565 * t584;
t63 = t745 + t581;
t69 = t378 * t521 + t991;
t679 = t63 * qJD(1) + t69 * qJD(2);
t663 = (t272 / 0.2e1 + t731) * t523;
t66 = t746 + t663;
t76 = t990 + (t824 - t1000) * t521;
t678 = t66 * qJD(1) + t76 * qJD(2);
t580 = t272 * t928 + t569 * t584;
t68 = t746 + t580;
t70 = t990 + (t825 - t1000) * t521;
t677 = t68 * qJD(1) + t70 * qJD(2);
t339 = t496 + t871;
t7 = t591 + t994;
t675 = -t7 * qJD(2) + t339 * qJD(3);
t674 = t878 - t879;
t138 = t215 * t521 + t698 * t877;
t586 = t928 * t998 + t738;
t92 = t730 - t902 / 0.2e1 - t586 * t565;
t671 = qJD(1) * t92 - qJD(2) * t138;
t139 = -t216 * t521 - t698 * t876;
t91 = t734 + t901 / 0.2e1 + t586 * t569;
t670 = qJD(1) * t91 - qJD(2) * t139;
t667 = -t557 / 0.2e1 - t775 / 0.2e1;
t585 = (t927 + t710 - pkin(10) / 0.2e1) * t523 + (-pkin(4) / 0.2e1 + t667) * t521;
t81 = t565 * t585;
t669 = t81 * qJD(2);
t668 = t523 * t713;
t666 = t915 / 0.2e1 - t914 / 0.2e1;
t662 = t521 * t927 + t557 * t928;
t35 = t751 + t994;
t652 = qJD(2) * t35 - t806 * t958;
t651 = t569 * t668;
t643 = pkin(4) / 0.2e1 + t667;
t642 = t666 * t569;
t641 = t697 * t984;
t637 = t662 * t569;
t613 = t565 * t666 + t749;
t127 = t613 + t727;
t447 = t643 * t569;
t630 = pkin(4) * t799 - qJD(2) * t127 + qJD(3) * t447;
t129 = -t398 / 0.2e1 - t642;
t446 = t643 * t565;
t629 = pkin(4) * t800 - qJD(2) * t129 + qJD(3) * t446;
t628 = (-t722 - t721) * pkin(3);
t612 = t565 * t662 + t749;
t123 = t612 - t728;
t592 = t1014 - t974 / 0.2e1 - t466 / 0.2e1;
t185 = t592 * t569;
t625 = qJD(1) * t185 - qJD(2) * t123 - t557 * t804;
t125 = -t378 / 0.2e1 - t637;
t181 = t592 * t565;
t624 = qJD(1) * t181 - qJD(2) * t125 - t557 * t805;
t623 = t980 / 0.2e1 + t887 / 0.2e1 + t868 / 0.2e1;
t622 = -t979 / 0.2e1 - t882 / 0.2e1 + t866 / 0.2e1;
t579 = -t980 / 0.2e1 + t891 / 0.2e1 - t874 / 0.2e1 + (-t564 * t763 - t568 * t764) * t920 / 0.2e1 + t365 * t710;
t31 = t579 + t623;
t616 = -t31 * qJD(2) + t711 * t958;
t578 = t979 / 0.2e1 + t889 / 0.2e1 - t872 / 0.2e1 - (t568 * t763 - t709) * t920 / 0.2e1 + t363 * t710;
t33 = t578 + t622;
t615 = -t33 * qJD(2) - t672 * t711;
t84 = t569 * t585;
t614 = -t84 * qJD(2) - t565 * t711;
t228 = t865 / 0.2e1 + t871 / 0.2e1 + t610;
t229 = t864 / 0.2e1 + t870 / 0.2e1 + t611;
t321 = t998 * t925;
t544 = t565 * t712;
t536 = t778 * qJD(5);
t485 = t672 * t712;
t484 = t958 * t712;
t449 = pkin(3) * t705 + pkin(4) * t924 + t557 * t923;
t448 = pkin(3) * t706 + pkin(4) * t926 + t557 * t925;
t419 = t880 / 0.2e1;
t418 = -t881 / 0.2e1;
t394 = t781 * t504;
t379 = t386 * qJD(5);
t333 = t381 + t797;
t332 = -t787 - t798;
t287 = t308 * qJD(2);
t285 = t307 * qJD(2);
t274 = t695 - t814;
t273 = -t697 + t814;
t258 = t958 * t714;
t256 = 0.2e1 * t565 * t651;
t225 = t521 * t655 - t672 * t930;
t224 = t521 * t654 - t744;
t221 = -t562 * t755 - t379;
t203 = t496 + t228;
t202 = -t495 + t229;
t194 = qJD(5) * t390 - t818;
t193 = -qJD(5) * t387 + t819;
t192 = t780 * t348;
t191 = -t924 * t998 + t1004;
t190 = 0.2e1 * t1004;
t189 = -t926 * t998 + t321;
t188 = t321 + t596;
t187 = t924 * t999 + t731;
t182 = t851 / 0.2e1 + t999 * t925;
t167 = -t714 - t213;
t166 = t780 * t958 - t214;
t149 = -t379 + (t562 * t811 - t697) * t521;
t148 = t972 * t932;
t147 = t972 * t933;
t141 = -t565 * t975 + t818;
t140 = -t523 * t715 - t819;
t137 = (qJD(5) - t812) * t847 * t984 + t778 * t411;
t130 = t1001 + t398 / 0.2e1 - t642;
t128 = t613 - t727;
t126 = t1001 + t378 / 0.2e1 - t637;
t124 = t612 + t728;
t94 = -t901 / 0.2e1 + t981 * t929 + t734 + t569 * t738;
t93 = t902 / 0.2e1 + t523 * t321 + t730 - t767 / 0.2e1;
t83 = pkin(10) * t743 + t966 * t569 + t978;
t82 = -t977 + pkin(10) * t877 / 0.2e1 + t966 * t565;
t80 = t223 * t781 - t365 * t812;
t79 = t222 * t781 - t363 * t812;
t74 = -t223 * t780 - t820;
t73 = -t222 * t780 - t821;
t67 = t418 + t580;
t65 = t418 + t663;
t64 = t419 + t581;
t62 = t419 + t664;
t60 = t224 * t781 + t365 * t688;
t59 = t225 * t781 + t363 * t688;
t48 = -t367 * t817 + t826;
t45 = t224 * t780 - t672 * t975 + t820;
t44 = t225 * t780 - t958 * t975 + t821;
t43 = t237 + t658 + t995;
t42 = t236 + t660 + t1008;
t37 = t237 + t659 + t996;
t36 = t236 + t661 + t997;
t32 = t578 - t622;
t30 = t579 - t623;
t29 = t568 * t769 + t858 + t938 - t844 / 0.2e1;
t27 = -t843 - t859 / 0.2e1 + (t769 + t937) * t564;
t25 = -(t717 - t817) * t367 + t826;
t24 = t701 * t564 + t1008 + t732 + t766;
t23 = -t701 * t568 + t735 + t765 + t995;
t20 = -t608 + t685;
t19 = -t573 + t685;
t18 = -t607 + t686;
t17 = -t572 + t686;
t10 = t702 * t564 + t733 + t766 + t997;
t9 = -t702 * t568 + t736 + t765 + t996;
t6 = -t903 + t911;
t1 = t903 + t781 * (t364 * t672 - t367 * t958) + t911;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t754, -t753, 0, 0, 0, 0, 0 (-t567 * t808 - t571 * t783) * t563 (t566 * t809 - t571 * t782) * t563, 0, 0, 0, 0, 0, t310 * t781 + t521 * t754, t309 * t781 - t523 * t754, 0, 0, 0, 0, 0 (t370 * t521 + t436 * t877) * qJD(2) + t64 * qJD(3) + t62 * qJD(4) + t94 * qJD(5) (-t372 * t521 + t436 * t876) * qJD(2) + t67 * qJD(3) + t65 * qJD(4) + t93 * qJD(5), 0, 0, 0, 0, 0 ((t840 - t854) * t521 + t436 * t365) * qJD(2) + t18 * qJD(3) + t17 * qJD(4) + t780 * t588 (-(t839 + t855) * t521 + t436 * t363) * qJD(2) + t20 * qJD(3) + t19 * qJD(4) + t780 * t589; 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t639 - t566 * t753, qJD(3) * t640 - t570 * t753, 0, 0, 0, 0, 0, t1017, t1012, 0, 0, 0, 0, 0, qJD(2) * t64 + qJD(4) * t187 + qJD(5) * t189 - t804 * t999, qJD(2) * t67 + qJD(4) * t182 + qJD(5) * t191 + t805 * t999, 0, 0, 0, 0, 0, qJD(2) * t18 + qJD(4) * t147 + t122 * t780 - t807 * t958, qJD(2) * t20 + qJD(4) * t148 + t120 * t780 + t672 * t807; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1017, t1012, 0, 0, 0, 0, 0, qJD(2) * t62 + qJD(3) * t187 + qJD(5) * t188 - t799 * t999, qJD(2) * t65 + qJD(3) * t182 + qJD(5) * t190 + t800 * t999, 0, 0, 0, 0, 0, qJD(2) * t17 + qJD(3) * t147 + t121 * t780 - t803 * t958, qJD(2) * t19 + qJD(3) * t148 + t119 * t780 + t672 * t803; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t94 + qJD(3) * t189 + qJD(4) * t188 - qJD(5) * t272, qJD(2) * t93 + qJD(3) * t191 + qJD(4) * t190 + qJD(5) * t271, 0, 0, 0, 0, 0, t945, t944; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t945, t944; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t781 * t307, -t781 * t308, 0, 0, 0, 0, 0, qJD(3) * t63 + qJD(4) * t61 - qJD(5) * t91, qJD(3) * t68 + qJD(4) * t66 - qJD(5) * t92, 0, 0, 0, 0, 0, -qJD(3) * t14 - qJD(4) * t13 - t577 * t780, -qJD(3) * t16 - qJD(4) * t15 - t576 * t780; 0, 0, 0, 0, t566 * t782, t546 * qJD(3), 0, 0, 0, -pkin(2) * t783, -pkin(2) * t782, t1003, -t781 * t779, 0, 0, 0, qJD(3) * t375 - t523 * t801, qJD(3) * t376 - t521 * t801, t1003 * t562 - t513 * t549, qJD(5) * t396 - t521 * t641, t306 * t781 + t565 * t756, -t305 * t781 + t569 * t756, -t1003, qJD(3) * t69 + qJD(4) * t75 + qJD(5) * t139, qJD(3) * t70 + qJD(4) * t76 + qJD(5) * t138 (-t365 * t780 + t367 * t781) * t363, t142 * t780 + t781 * t98, t158 * t781 + t365 * t982, t157 * t781 + t363 * t982, -t1003, qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t38 + qJD(6) * t58, qJD(3) * t3 + qJD(4) * t5 + qJD(5) * t39 + qJD(6) * t57; 0, 0, 0, 0, t752, t784, t782, -t783, 0, -pkin(8) * t782 - t773, pkin(8) * t783 - t772, t755, -t788, -t411, t975, 0, -t964 - t988, t815 + t989, t149, t137, t141, t140, -t312 (t565 * t674 - t977) * qJD(3) + t82 * qJD(4) + t126 * qJD(5) + t679 (t569 * t674 + t978) * qJD(3) + t83 * qJD(4) + t124 * qJD(5) + t677, t25, t1, t45, t44, t230 (-t874 + t891 - t980) * qJD(3) + t30 * qJD(4) + t9 * qJD(5) + t37 * qJD(6) + t693 (-t872 + t889 + t979) * qJD(3) + t32 * qJD(4) + t10 * qJD(5) + t36 * qJD(6) + t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t755, -t788, -t411, t975, 0, -t965 - t988, -t759 + t989, t149, t137, t141, t140, -t312, t82 * qJD(3) + (t565 * t700 - t977) * qJD(4) + t130 * qJD(5) + t680, t83 * qJD(3) + (t569 * t700 + t978) * qJD(4) + t128 * qJD(5) + t678, t25, t1, t45, t44, t230, t30 * qJD(3) + (-t868 - t887 - t980) * qJD(4) + t23 * qJD(5) + t43 * qJD(6) + t694, t32 * qJD(3) + (-t866 + t882 + t979) * qJD(4) + t24 * qJD(5) + t42 * qJD(6) + t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t949, t695 * t984 + t813, t565 * t668, t651, t394, qJD(3) * t126 + qJD(4) * t130 - qJD(5) * t216 - t670, qJD(3) * t124 + qJD(4) * t128 + qJD(5) * t215 - t671, t46, t11, t60, t59, t394, qJD(3) * t9 + qJD(4) * t23 + qJD(5) * t85 + qJD(6) * t27 - t684, qJD(3) * t10 + qJD(4) * t24 - qJD(5) * t86 + qJD(6) * t29 - t683; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t11, t60, t59, t394, qJD(3) * t37 + qJD(4) * t43 + qJD(5) * t27 - qJD(6) * t78 - t682, qJD(3) * t36 + qJD(4) * t42 + qJD(5) * t29 + qJD(6) * t77 - t681; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, t287, 0, 0, 0, 0, 0, -qJD(2) * t63 - qJD(5) * t181, -qJD(2) * t68 - qJD(5) * t185, 0, 0, 0, 0, 0, qJD(2) * t14 - t116 * t780, qJD(2) * t16 - t118 * t780; 0, 0, 0, 0, -t752, -t784, 0, 0, 0, t773, t772, -t755, t788, 0, 0, 0, t964, t789 - t815, t221, t256, t194, t193, t312, qJD(4) * t81 + qJD(5) * t125 - t679, qJD(4) * t84 + qJD(5) * t123 - t677, t48, t6, t74, t73, -t230, qJD(4) * t31 - qJD(5) * t8 - qJD(6) * t34 - t693, qJD(4) * t33 - qJD(5) * t7 - qJD(6) * t35 - t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t712, -pkin(3) * t723, t549, -t536, 0, 0, 0, t557 * t798 - t569 * t712, t557 * t797 + t544, t258, t192, 0, 0, 0, qJD(5) * t338 + t672 * t796 - t484, qJD(5) * t339 + t796 * t958 + t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t628 (-t724 - t723) * pkin(3), t549, -t536, 0, 0, 0, t448 * qJD(5) + t569 * t628 + t669, t449 * qJD(5) + t544 - t614, t258, t192, 0, 0, 0, t202 * qJD(5) + t229 * qJD(6) - t484 - t616, t203 * qJD(5) + t228 * qJD(6) + t485 - t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, -t350, t333, t332, -t786, qJD(4) * t448 - t556 * t797 - t624, qJD(4) * t449 + t556 * t798 - t625, t89, t71, t166, t167, -t786, t202 * qJD(4) - t961 - t987, t203 * qJD(4) + t675 + t943; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t71, t166, t167, -t786, qJD(4) * t229 - t971 - t987, qJD(4) * t228 - t652 + t943; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, t287, 0, 0, 0, 0, 0, -qJD(2) * t61, -qJD(2) * t66, 0, 0, 0, 0, 0, qJD(2) * t13 - t115 * t780, qJD(2) * t15 - t117 * t780; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t755, t788, 0, 0, 0, t965, t789 + t759, t221, t256, t194, t193, t312, -qJD(3) * t81 + qJD(5) * t129 - t680, -qJD(3) * t84 + qJD(5) * t127 - t678, t48, t6, t74, t73, -t230, -qJD(3) * t31 - qJD(5) * t22 - qJD(6) * t40 - t694, -qJD(3) * t33 - qJD(5) * t21 - qJD(6) * t41 - t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t711, pkin(3) * t724, t549, -t536, 0, 0, 0, -t446 * qJD(5) + t569 * t711 - t669, -t447 * qJD(5) + t614, t258, t192, 0, 0, 0, -t200 * qJD(5) - t226 * qJD(6) + t616, -t201 * qJD(5) - t227 * qJD(6) + t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t549, -t536, 0, 0, 0, -pkin(4) * t798, -pkin(4) * t797, t258, t192, 0, 0, 0, qJD(5) * t373 + t672 * t795, qJD(5) * t374 + t795 * t958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, -t350, t333, t332, -t786, -pkin(10) * t797 - t629, pkin(10) * t798 - t630, t89, t71, t166, t167, -t786, -t959 - t986, -t962 - t985; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t71, t166, t167, -t786, -t960 - t986, -t963 - t985; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t91 + qJD(3) * t181, qJD(2) * t92 + qJD(3) * t185, 0, 0, 0, 0, 0, t947, t948; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t949, -t641 - t813, -t390 * t781 - t565 * t755, t387 * t781 - t569 * t755, t394, -qJD(3) * t125 - qJD(4) * t129 + t670, -qJD(3) * t123 - qJD(4) * t127 + t671, -t46, -t11, t80, t79, t394, qJD(3) * t8 + qJD(4) * t22 + qJD(6) * t26 + t684, qJD(3) * t7 + qJD(4) * t21 + qJD(6) * t28 + t683; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, t350, -t381, t787, t786, qJD(4) * t446 + t624, qJD(4) * t447 + t625, -t89, -t71, t214, t213, t786, qJD(4) * t200 + t961, qJD(4) * t201 - t675 + t791; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, t350, -t381, t787, t786, t629, t630, -t89, -t71, t214, t213, t786, t959, t962; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t564 * t909, -t568 * t909; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t564 * t720 + t906, -t568 * t720 + t905; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t947, t948; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t11, t80, t79, t394, qJD(3) * t34 + qJD(4) * t40 - qJD(5) * t26 + t682, qJD(3) * t35 + qJD(4) * t41 - qJD(5) * t28 + t681; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t71, t214, t213, t786, qJD(4) * t226 + t971, qJD(4) * t227 + t652 + t791; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t71, t214, t213, t786, t960, t963; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t564 * t910 - t906, t568 * t910 - t905; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t12;
