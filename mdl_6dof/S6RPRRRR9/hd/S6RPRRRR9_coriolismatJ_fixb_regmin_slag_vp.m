% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x34]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:22
% EndTime: 2019-03-09 07:26:14
% DurationCPUTime: 34.34s
% Computational Cost: add. (12422->766), mult. (26707->987), div. (0->0), fcn. (29009->8), ass. (0->616)
t999 = qJD(5) + qJD(6);
t1008 = qJD(4) + t999;
t571 = cos(qJ(4));
t926 = pkin(8) + pkin(9);
t509 = t926 * t571;
t903 = cos(qJ(5));
t493 = t903 * t509;
t568 = sin(qJ(4));
t507 = t926 * t568;
t567 = sin(qJ(5));
t841 = t567 * t507;
t394 = t493 - t841;
t554 = t903 * t571;
t840 = t567 * t568;
t485 = -t554 + t840;
t312 = -pkin(10) * t485 + t394;
t566 = sin(qJ(6));
t570 = cos(qJ(6));
t552 = t567 * t571;
t553 = t903 * t568;
t623 = t553 + t552;
t690 = t507 * t903 + t509 * t567;
t936 = -pkin(10) * t623 - t690;
t171 = t570 * t312 + t566 * t936;
t1005 = t1008 * t171;
t1002 = t566 * t312 - t570 * t936;
t1004 = t1008 * t1002;
t569 = sin(qJ(3));
t907 = -t569 / 0.2e1;
t1003 = t171 * t907;
t906 = t569 / 0.2e1;
t1000 = t906 * t1002;
t750 = qJD(4) + qJD(5);
t967 = qJD(6) + t750;
t572 = cos(qJ(3));
t512 = t572 * t554;
t837 = t568 * t572;
t444 = t567 * t837 - t512;
t431 = t570 * t444;
t419 = t431 / 0.2e1;
t446 = t623 * t572;
t848 = t566 * t446;
t681 = t419 + t848 / 0.2e1;
t470 = t570 * t485;
t845 = t566 * t623;
t366 = -t470 - t845;
t947 = t572 * t366;
t965 = t947 / 0.2e1 + t681;
t471 = t570 * t623;
t460 = -t471 / 0.2e1;
t846 = t566 * t485;
t680 = t460 + t846 / 0.2e1;
t445 = t623 * t569;
t432 = t570 * t445;
t838 = t568 * t569;
t511 = t567 * t838;
t448 = t554 * t569 - t511;
t847 = t566 * t448;
t650 = t432 + t847;
t613 = t650 * t569;
t433 = t570 * t446;
t850 = t566 * t444;
t302 = -t433 + t850;
t948 = t572 * t302;
t972 = t613 / 0.2e1 - t948 / 0.2e1;
t963 = t680 - t972;
t977 = t963 * qJD(1);
t45 = qJD(3) * t965 + t977;
t966 = -t947 / 0.2e1 + t681;
t991 = qJD(3) * t966 - t977;
t410 = t850 / 0.2e1;
t748 = t410 - t433 / 0.2e1;
t691 = -t471 + t846;
t971 = t572 * t691;
t973 = t971 / 0.2e1 + t748;
t988 = qJD(3) * t973;
t974 = -t971 / 0.2e1 + t748;
t987 = qJD(3) * t974;
t408 = t847 / 0.2e1;
t749 = t408 + t432 / 0.2e1;
t970 = t691 * t569;
t975 = t970 / 0.2e1 + t749;
t986 = qJD(3) * t975;
t773 = t974 * qJD(2);
t976 = -t970 / 0.2e1 + t749;
t985 = t976 * qJD(1);
t755 = t569 * qJD(1);
t938 = -t431 - t848;
t68 = qJD(3) * t976 + t755 * t938;
t960 = t366 ^ 2 - t691 ^ 2;
t982 = qJD(3) * t960;
t849 = t566 * t445;
t434 = t570 * t448;
t917 = -t434 / 0.2e1;
t682 = t917 + t849 / 0.2e1;
t952 = t366 * t569;
t961 = t952 / 0.2e1 + t682;
t981 = qJD(3) * t961;
t962 = -t952 / 0.2e1 + t682;
t978 = t962 * qJD(1);
t106 = t963 * qJD(2);
t964 = t680 + t972;
t108 = t964 * qJD(2);
t771 = t965 * qJD(2);
t69 = qJD(3) * t962 - t302 * t755;
t901 = pkin(4) * t571;
t556 = -pkin(3) - t901;
t437 = pkin(5) * t485 + t556;
t573 = -pkin(1) - pkin(7);
t902 = pkin(4) * t568;
t695 = -t573 + t902;
t665 = t695 * t572;
t897 = t446 * pkin(5);
t602 = t665 + t897;
t918 = -t691 / 0.2e1;
t924 = t938 / 0.2e1;
t940 = t437 * t924 + t602 * t918;
t958 = t366 / 0.2e1;
t969 = t938 * t958;
t945 = t302 ^ 2 - t938 ^ 2;
t968 = qJD(1) * t945;
t563 = t569 ^ 2;
t957 = t563 / 0.2e1;
t797 = qJD(1) * t938;
t953 = t302 * t797;
t790 = qJD(3) * t691;
t951 = t366 * t790;
t210 = t602 * t302;
t589 = t602 * t938;
t944 = 0.2e1 * t568;
t942 = t569 * t750;
t588 = t602 * t366;
t873 = t437 * t302;
t939 = t588 / 0.2e1 + t873 / 0.2e1;
t937 = -t434 + t849;
t672 = pkin(3) * t569 - pkin(8) * t572;
t498 = qJ(2) + t672;
t482 = t571 * t498;
t836 = t568 * t573;
t696 = pkin(4) - t836;
t824 = t571 * t572;
t743 = pkin(9) * t824;
t358 = t569 * t696 + t482 - t743;
t335 = t903 * t358;
t823 = t571 * t573;
t733 = t569 * t823;
t402 = t498 * t568 + t733;
t383 = -pkin(9) * t837 + t402;
t843 = t567 * t383;
t213 = -t335 + t843;
t898 = t444 * pkin(10);
t637 = -t213 + t898;
t608 = t566 * t637;
t593 = t608 / 0.2e1;
t399 = t566 * pkin(5);
t737 = -t399 / 0.2e1;
t361 = t903 * t383;
t214 = t567 * t358 + t361;
t896 = t446 * pkin(10);
t168 = t214 - t896;
t834 = t570 * t168;
t707 = -t834 / 0.2e1;
t152 = t569 * pkin(5) + t637;
t861 = t566 * t152;
t813 = -t861 / 0.2e1 + t707;
t639 = t569 * t737 + t813;
t23 = t593 + t834 / 0.2e1 + t639;
t916 = t434 / 0.2e1;
t283 = t917 + t916;
t765 = t283 * qJD(2);
t935 = t23 * qJD(1) - t399 * qJD(4) + t765;
t557 = t569 * t573;
t401 = t557 * t568 - t482;
t382 = -t401 - t743;
t731 = t903 * t382;
t223 = t731 - t843;
t173 = t223 + t898;
t844 = t567 * t382;
t222 = -t361 - t844;
t172 = t222 + t896;
t703 = t172 / 0.2e1 + t168 / 0.2e1;
t670 = t703 * t570;
t741 = t903 * pkin(4);
t679 = t741 + pkin(5);
t641 = t566 * t679;
t839 = t567 * t570;
t466 = pkin(4) * t839 + t641;
t868 = t466 * t569;
t19 = t868 / 0.2e1 + t670 + (-t173 / 0.2e1 + t152 / 0.2e1) * t566;
t934 = qJD(1) * t19 + qJD(4) * t466 - t765;
t606 = t552 / 0.2e1 + t553 / 0.2e1;
t912 = t623 / 0.2e1;
t340 = (t912 - t606) * t572;
t761 = t340 * qJD(2);
t787 = qJD(3) * t556;
t892 = t572 * pkin(3);
t893 = t569 * pkin(8);
t508 = t892 + t893;
t494 = t571 * t508;
t363 = t571 * t569 * pkin(9) + t572 * t696 + t494;
t347 = t903 * t363;
t491 = t568 * t508;
t817 = t572 * t573;
t514 = t571 * t817;
t386 = pkin(9) * t838 + t491 + t514;
t842 = t567 * t386;
t667 = t347 / 0.2e1 - t842 / 0.2e1;
t863 = t556 * t444;
t601 = t863 / 0.2e1 + t667;
t635 = t665 / 0.2e1;
t875 = t394 * t569;
t91 = t875 / 0.2e1 - t623 * t635 + t601;
t933 = qJD(1) * t91 - t623 * t787 + t761;
t336 = t485 * t902 + t556 * t623;
t740 = t902 / 0.2e1;
t584 = -t394 * t907 - t446 * t740;
t683 = t741 / 0.2e1;
t911 = -t623 / 0.2e1;
t914 = -t485 / 0.2e1;
t57 = (t695 * t911 + t901 * t914 + t683) * t572 + t584 + t601;
t932 = qJD(1) * t57 - qJD(3) * t336 + t761;
t257 = t437 * t366;
t899 = pkin(5) * t623;
t130 = t691 * t899 - t257;
t614 = t1002 * t569;
t346 = t567 * t363;
t375 = t903 * t386;
t803 = t375 + t346;
t169 = pkin(10) * t445 + t803;
t833 = t570 * t169;
t706 = -t833 / 0.2e1;
t692 = t347 - t842;
t891 = t572 * pkin(5);
t157 = pkin(10) * t448 + t692 + t891;
t860 = t566 * t157;
t629 = -t860 / 0.2e1 + t706;
t636 = -t665 / 0.2e1;
t714 = -t873 / 0.2e1;
t905 = -t572 / 0.2e1;
t7 = -t614 / 0.2e1 + t714 + t366 * t636 + (t444 * t918 - t446 * t958 + t566 * t905 + t911 * t938) * pkin(5) + t629;
t931 = qJD(1) * t7 + qJD(3) * t130 + t771;
t438 = t899 + t902;
t124 = -t438 * t691 + t257;
t744 = pkin(4) * t824;
t900 = pkin(5) * t444;
t642 = t744 - t900;
t915 = -t438 / 0.2e1;
t919 = t691 / 0.2e1;
t575 = t642 * t919 + t915 * t938 - t1000;
t586 = t466 * t905 + t629;
t2 = t714 - t588 / 0.2e1 + t575 + t586;
t930 = qJD(1) * t2 - qJD(3) * t124 + t771;
t256 = t437 * t691;
t129 = t366 * t899 + t256;
t858 = t566 * t169;
t710 = -t858 / 0.2e1;
t835 = t570 * t157;
t628 = t710 + t835 / 0.2e1;
t574 = t628 - t940;
t738 = t899 / 0.2e1;
t739 = -t900 / 0.2e1;
t579 = -t171 * t906 - t302 * t738 - t366 * t739;
t734 = t891 / 0.2e1;
t8 = t570 * t734 + t574 - t579;
t929 = qJD(1) * t8 + qJD(3) * t129 + t773;
t576 = -t302 * t915 + t642 * t958 - t1003;
t510 = t570 * t679;
t894 = t567 * pkin(4);
t545 = t566 * t894;
t465 = t545 - t510;
t713 = t465 * t905;
t1 = t713 + t574 + t576;
t123 = -t366 * t438 - t256;
t928 = qJD(1) * t1 - qJD(3) * t123 + t773;
t562 = t568 ^ 2;
t564 = t571 ^ 2;
t530 = t564 - t562;
t693 = t824 * t944;
t603 = qJD(1) * t693 - qJD(3) * t530;
t708 = t840 / 0.2e1;
t909 = -t512 / 0.2e1;
t338 = t909 + (t708 + t914) * t572;
t762 = t338 * qJD(2);
t927 = t690 * t750 - t762;
t925 = t302 / 0.2e1;
t678 = -t361 / 0.2e1;
t913 = t485 / 0.2e1;
t677 = -t493 / 0.2e1;
t910 = t511 / 0.2e1;
t904 = t572 / 0.2e1;
t890 = pkin(5) * qJD(5);
t889 = pkin(5) * qJD(6);
t64 = -t608 - t834;
t33 = t302 * t900 + t64 * t569 + t589;
t888 = qJD(1) * t33;
t607 = t570 * t637;
t859 = t566 * t168;
t65 = t607 - t859;
t34 = t569 * t65 + t900 * t938 - t210;
t887 = qJD(1) * t34;
t832 = t570 * t172;
t856 = t566 * t173;
t66 = t832 - t856;
t35 = -t302 * t642 + t66 * t569 + t589;
t886 = qJD(1) * t35;
t831 = t570 * t173;
t857 = t566 * t172;
t67 = t831 + t857;
t36 = t569 * t67 - t642 * t938 - t210;
t885 = qJD(1) * t36;
t136 = t570 * t152;
t60 = -t136 + t859;
t42 = -t569 * t60 - t210;
t884 = qJD(1) * t42;
t61 = t834 + t861;
t43 = -t61 * t569 + t589;
t883 = qJD(1) * t43;
t881 = t937 * t569;
t876 = t690 * t569;
t872 = t445 * t569;
t871 = t446 * t572;
t870 = t448 * t569;
t869 = t465 * t569;
t730 = t903 * t566;
t472 = (t730 + t839) * pkin(4);
t867 = t472 * t569;
t473 = t570 * t741 - t545;
t866 = t473 * t569;
t484 = -pkin(4) * t838 + t557;
t52 = t692 * t569 + t484 * t446 + (-t445 * t695 - t213) * t572;
t865 = t52 * qJD(1);
t53 = t214 * t572 + t444 * t484 + t448 * t665 + t569 * t803;
t864 = t53 * qJD(1);
t862 = t556 * t446;
t565 = t572 ^ 2;
t555 = t565 * t571;
t825 = t571 * t563;
t821 = t572 * t938;
t818 = t572 * t444;
t63 = (-t432 / 0.2e1 - t847 / 0.2e1 + t650 / 0.2e1) * t572;
t816 = t63 * qJD(1);
t83 = t302 * t937 + t650 * t938;
t815 = t83 * qJD(1);
t684 = -t741 / 0.2e1;
t615 = -t335 / 0.2e1 + t569 * t684;
t97 = t731 / 0.2e1 + t615;
t814 = t97 * qJD(1);
t711 = t859 / 0.2e1;
t812 = -t136 / 0.2e1 + t711;
t142 = t444 * t623 + t446 * t485;
t811 = t750 * t142;
t810 = -t608 / 0.2e1 + t707;
t809 = t711 - t607 / 0.2e1;
t234 = t444 * t913 - t446 * t912;
t808 = t750 * t234;
t529 = t563 - t565;
t802 = qJD(1) * qJ(2);
t117 = t222 * t569 + (-t444 * t695 + t446 * t901) * t572;
t801 = qJD(1) * t117;
t373 = t446 * t665;
t118 = t223 * t569 + t444 * t744 + t373;
t800 = qJD(1) * t118;
t121 = -t213 * t569 + t373;
t799 = qJD(1) * t121;
t122 = -t214 * t569 - t444 * t665;
t798 = qJD(1) * t122;
t348 = -t401 * t569 - t565 * t836;
t796 = qJD(1) * t348;
t349 = -t402 * t569 - t565 * t823;
t795 = qJD(1) * t349;
t794 = qJD(1) * t444;
t489 = t529 * t568;
t793 = qJD(1) * t489;
t490 = -t555 + t825;
t792 = qJD(1) * t490;
t791 = qJD(2) * t569;
t789 = qJD(3) * t437;
t788 = qJD(3) * t623;
t786 = qJD(3) * t571;
t785 = qJD(4) * t568;
t784 = qJD(4) * t571;
t783 = qJD(5) * t556;
t587 = t938 * t904 - t881 / 0.2e1;
t456 = -t845 / 0.2e1;
t747 = t456 - t470 / 0.2e1;
t109 = t587 + t747;
t782 = t109 * qJD(1);
t583 = t821 / 0.2e1 + t937 * t907;
t110 = t583 + t747;
t781 = t110 * qJD(1);
t143 = (-t849 / 0.2e1 + t916 + t937 / 0.2e1) * t572;
t778 = t143 * qJD(1);
t153 = t948 + t613;
t777 = t153 * qJD(1);
t154 = t821 + t881;
t776 = t154 * qJD(1);
t212 = -t444 * t445 + t446 * t448;
t770 = t212 * qJD(1);
t732 = t568 * t817;
t237 = -t401 * t572 + (t494 + t732) * t569;
t769 = t237 * qJD(1);
t238 = t402 * t572 + (-t514 + t491) * t569;
t768 = t238 * qJD(1);
t630 = -t870 / 0.2e1 + t818 / 0.2e1;
t666 = t554 / 0.2e1 - t840 / 0.2e1;
t249 = -t630 + t666;
t767 = t249 * qJD(1);
t631 = -t872 / 0.2e1 - t871 / 0.2e1;
t250 = t631 - t606;
t766 = t250 * qJD(1);
t305 = -t871 + t872;
t764 = t305 * qJD(1);
t306 = -t818 - t870;
t763 = t306 * qJD(1);
t339 = (t912 + t606) * t569;
t318 = t339 * qJD(1);
t676 = -t554 / 0.2e1;
t341 = t910 + (t676 + t913) * t569;
t320 = t341 * qJD(1);
t760 = t366 * qJD(6);
t759 = t691 * qJD(6);
t742 = 0.1e1 / 0.2e1 + t957;
t449 = (-t565 / 0.2e1 - t742) * t568;
t758 = t449 * qJD(1);
t450 = t555 / 0.2e1 + t742 * t571;
t757 = t450 * qJD(1);
t756 = t529 * qJD(1);
t754 = t569 * qJD(3);
t753 = t572 * qJD(1);
t752 = t572 * qJD(3);
t751 = t572 * qJD(4);
t39 = t302 * t958 + t366 * t925 + t691 * t924 + t919 * t938;
t49 = t302 * t366 + t691 * t938;
t746 = t49 * qJD(6) + t39 * t750;
t89 = t302 * t918 + t969;
t90 = -t302 * t919 + t969;
t745 = t90 * qJD(6) + t750 * t89;
t736 = -t570 * pkin(5) / 0.2e1;
t729 = qJ(2) * t755;
t728 = qJ(2) * t753;
t727 = t568 * t786;
t726 = t568 * t752;
t725 = t571 * t752;
t724 = t569 * t785;
t723 = t568 * t751;
t722 = t569 * t784;
t721 = t571 * t751;
t720 = t366 * t755;
t719 = t691 * t755;
t718 = t485 * t755;
t717 = t623 * t755;
t716 = t568 * t784;
t544 = t569 * t752;
t715 = t569 * t753;
t712 = t485 * t904;
t709 = -t856 / 0.2e1;
t705 = -t831 / 0.2e1;
t700 = -t346 / 0.2e1 - t375 / 0.2e1;
t698 = t903 * qJD(4);
t697 = t903 * qJD(5);
t694 = pkin(5) * t999;
t688 = t750 * t366;
t687 = t750 * t623;
t686 = -qJD(4) - t755;
t685 = -qJD(6) - t755;
t675 = pkin(4) * t907 - t358 / 0.2e1;
t674 = t755 + qJD(4) / 0.2e1;
t673 = t734 + t157 / 0.2e1;
t669 = qJD(3) * t693;
t362 = -pkin(5) * t445 + t484;
t5 = (t835 - t858) * t569 - t362 * t302 - t650 * t897 + (-t650 * t695 - t60) * t572;
t664 = qJD(1) * t5 + qJD(2) * t63;
t663 = qJD(5) - t686;
t660 = t485 * t636 - t862 / 0.2e1 + t700;
t6 = (t833 + t860) * t569 + t61 * t572 - t362 * t938 - t602 * t937;
t659 = -qJD(1) * t6 - qJD(2) * t143;
t11 = qJD(3) * t39 + t968;
t13 = qJD(1) * t39 + t982;
t648 = qJD(3) * t49 + t968;
t647 = qJD(1) * t49 + t982;
t337 = -t485 * t556 + t623 * t902;
t585 = t444 * t740 + t690 * t907;
t622 = t862 / 0.2e1 + t700;
t56 = (-t894 / 0.2e1 + t901 * t911 + t695 * t913) * t572 + t585 + t622;
t645 = -qJD(1) * t56 + qJD(3) * t337;
t390 = t677 + t493 / 0.2e1;
t95 = t678 + t361 / 0.2e1 + (t382 / 0.2e1 + t675) * t567;
t644 = t95 * qJD(1) + t390 * qJD(3);
t643 = t686 * t572;
t253 = -t444 ^ 2 + t446 ^ 2;
t87 = qJD(1) * t253 + qJD(3) * t142;
t307 = t485 ^ 2 - t623 ^ 2;
t99 = qJD(1) * t142 + qJD(3) * t307;
t640 = t893 / 0.2e1 + t892 / 0.2e1;
t638 = t569 * t736 + t812;
t634 = qJD(5) / 0.2e1 + t674;
t609 = t640 * t568;
t384 = t491 / 0.2e1 + t609;
t633 = pkin(3) * t786 - qJD(1) * t384;
t610 = t640 * t571;
t385 = -t494 / 0.2e1 - t610;
t632 = pkin(3) * qJD(3) * t568 - qJD(1) * t385;
t627 = -t857 / 0.2e1 + t705;
t626 = t709 + t832 / 0.2e1;
t621 = -qJD(3) * t90 - t953;
t46 = qJD(3) * t89 + t953;
t620 = -qJD(1) * t90 + t951;
t50 = qJD(1) * t89 - t951;
t92 = -t876 / 0.2e1 + t485 * t635 + t622;
t619 = qJD(1) * t92 + t485 * t787;
t617 = t571 * t643;
t126 = -qJD(3) * t234 - t446 * t794;
t150 = qJD(1) * t234 - t485 * t788;
t475 = (t562 / 0.2e1 - t564 / 0.2e1) * t572;
t616 = -qJD(1) * t475 + t727;
t612 = -t623 * t636 - t863 / 0.2e1 + t667;
t611 = t967 * t302;
t605 = qJD(1) * t555 * t568 + qJD(3) * t475;
t488 = t530 * t565;
t604 = qJD(1) * t488 + t669;
t594 = -t566 * t703 + t705;
t20 = -t869 / 0.2e1 + t136 / 0.2e1 + t594;
t599 = qJD(1) * t20 - qJD(4) * t465;
t592 = t607 / 0.2e1;
t25 = t592 - t859 / 0.2e1 + t638;
t400 = t510 / 0.2e1 + (t684 + pkin(5) / 0.2e1) * t570;
t597 = t25 * qJD(1) - t400 * qJD(4);
t27 = t709 + t867 / 0.2e1 + t593 + t670;
t596 = qJD(1) * t27 + qJD(4) * t472;
t28 = t866 / 0.2e1 + t592 + t594;
t595 = qJD(1) * t28 + qJD(4) * t473;
t577 = t940 + t1003;
t15 = -t577 + t628;
t591 = qJD(1) * t15 + t691 * t789 + t773;
t578 = t437 * t925 + t602 * t958 + t1000;
t16 = -t578 + t629;
t590 = qJD(1) * t16 - t366 * t789 + t771;
t551 = -t753 / 0.2e1;
t550 = t753 / 0.2e1;
t549 = t752 / 0.2e1;
t543 = t571 * t755;
t542 = t568 * t754;
t541 = t568 * t755;
t483 = t674 * t572;
t469 = t475 * qJD(4);
t468 = t473 * qJD(5);
t467 = t472 * qJD(5);
t452 = -t825 / 0.2e1 - t555 / 0.2e1 + t571 / 0.2e1;
t451 = (-0.1e1 / 0.2e1 + t957 + t565 / 0.2e1) * t568;
t443 = t466 * qJD(6);
t442 = t465 * qJD(6);
t441 = t634 * t572;
t396 = (qJD(6) / 0.2e1 + t634) * t572;
t381 = t736 + t545 - t510 / 0.2e1 + t570 * t684;
t380 = t737 - t641 / 0.2e1 + (-t839 - t730 / 0.2e1) * pkin(4);
t345 = t569 * t606 - t623 * t906;
t344 = t572 * t708 + t712 + t909;
t343 = t485 * t907 + t569 * t676 + t910;
t342 = -t572 * t606 - t623 * t904;
t324 = -t732 + t494 / 0.2e1 - t610;
t323 = -t514 - t491 / 0.2e1 + t609;
t314 = 0.2e1 * t677 + t841;
t259 = 0.2e1 * t460 + t846;
t258 = 0.2e1 * t456 - t470;
t252 = t630 + t666;
t251 = -t631 - t606;
t246 = t252 * qJD(2);
t245 = t251 * qJD(2);
t244 = t250 * qJD(2);
t243 = t249 * qJD(2);
t242 = qJD(3) * t341 + t446 * t755;
t241 = qJD(3) * t339 - t444 * t755;
t229 = 0.2e1 * t419 + t848;
t228 = 0.2e1 * t917 + t849;
t227 = 0.2e1 * t410 - t433;
t226 = 0.2e1 * t408 + t432;
t217 = -t687 - t318;
t216 = -t485 * t750 - t320;
t139 = t143 * qJD(3);
t135 = qJD(3) * t338 + t766;
t134 = qJD(3) * t340 + t767;
t128 = t343 * qJD(3) - t446 * t663;
t127 = t345 * qJD(3) + t444 * t663;
t116 = -t583 + t747;
t113 = -t587 + t747;
t107 = t113 * qJD(2);
t105 = t109 * qJD(2);
t102 = t342 * qJD(3) - t448 * t750 - t767;
t101 = t344 * qJD(3) + t445 * t750 - t766;
t98 = t843 - t731 / 0.2e1 + t615;
t96 = 0.2e1 * t678 - t844 / 0.2e1 + t675 * t567;
t94 = -t875 / 0.2e1 + t612;
t93 = t876 / 0.2e1 + t660;
t62 = t63 * qJD(3);
t59 = -t585 + t660 + (-t571 * t623 + t567) * pkin(4) * t905;
t58 = t572 * t683 + t712 * t901 - t584 + t612;
t55 = t258 * qJD(6) + t688 - t978;
t54 = t259 * qJD(6) + t691 * t750 - t985;
t44 = qJD(6) * t283 + t782 + t987;
t41 = t227 * qJD(6) + t302 * t663 + t981;
t40 = t229 * qJD(6) - t663 * t938 + t986;
t32 = t226 * qJD(6) + t650 * t750 + t991;
t31 = t228 * qJD(6) + t750 * t937 - t782 + t988;
t30 = -t866 / 0.2e1 + t627 + t809;
t29 = -t867 / 0.2e1 + t626 + t810;
t26 = t638 + t809;
t24 = t639 + t810;
t22 = -t868 / 0.2e1 + t626 + t813;
t21 = t869 / 0.2e1 + t627 + t812;
t18 = t577 + t628;
t17 = t578 + t629;
t10 = t614 / 0.2e1 + t938 * t738 - t691 * t739 + t706 - t673 * t566 + t939;
t9 = t570 * t673 + t579 + t710 + t940;
t4 = -t575 + t586 + t939;
t3 = t713 - t576 + t628 + t940;
t12 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t544, t529 * qJD(3), 0, 0, 0, qJ(2) * t752 + t791, -qJ(2) * t754 + qJD(2) * t572, -t544 * t564 - t565 * t716, -qJD(4) * t488 + t569 * t669, -qJD(3) * t490 - t569 * t723, qJD(3) * t489 - t569 * t721, t544, qJD(3) * t237 + qJD(4) * t349 + t571 * t791, -qJD(3) * t238 - qJD(4) * t348 - t568 * t791 (qJD(3) * t448 + t446 * t750) * t444, t212 * qJD(3) + t253 * t750, t306 * qJD(3) - t446 * t942, t305 * qJD(3) + t444 * t942, t544, qJD(3) * t52 + qJD(4) * t117 + qJD(5) * t122 - t485 * t791, -qJD(3) * t53 - qJD(4) * t118 - qJD(5) * t121 - t623 * t791 (qJD(3) * t937 + t611) * t938, qJD(3) * t83 + t945 * t967, t154 * qJD(3) + t569 * t611, -t569 * t938 * t967 + t153 * qJD(3), t544, qJD(3) * t5 + qJD(4) * t35 + qJD(5) * t33 + qJD(6) * t43 + t366 * t791, -qJD(3) * t6 - qJD(4) * t36 - qJD(5) * t34 - qJD(6) * t42 + t691 * t791; 0, 0, 0, 0, qJD(1), t802, 0, 0, 0, 0, 0, t755, t753, 0, 0, 0, 0, 0, qJD(4) * t452 + t543, qJD(4) * t451 - t541, 0, 0, 0, 0, 0, t252 * t750 - t718, t251 * t750 - t717, 0, 0, 0, 0, 0, t116 * qJD(6) + t113 * t750 + t62 + t720, t964 * t967 - t139 + t719; 0, 0, 0, 0, 0, 0, -t715, t756, -t754, -t752, 0, -t573 * t754 + t728, -t573 * t752 - t729, -t469 + (-t564 * t753 - t727) * t569, t569 * t603 - 0.2e1 * t572 * t716, t726 - t792, t725 + t793, t483, t769 + (t568 * t672 - t733) * qJD(3) + t324 * qJD(4), -t768 + (-pkin(8) * t824 + (pkin(3) * t571 + t836) * t569) * qJD(3) + t323 * qJD(4) (-t788 + t794) * t448 + t808, t770 + (t445 * t623 + t448 * t485) * qJD(3) + t811, t343 * t750 + t623 * t752 + t763, t345 * t750 - t485 * t752 + t764, t441, t865 + (-t445 * t556 + t484 * t485 - t572 * t690) * qJD(3) + t58 * qJD(4) + t94 * qJD(5), -t864 + (-t394 * t572 - t448 * t556 + t484 * t623) * qJD(3) + t59 * qJD(4) + t93 * qJD(5) (-t790 + t797) * t937 + t745, t815 + (t366 * t937 - t650 * t691) * qJD(3) + t746, -t691 * t752 + t961 * t967 + t776, t366 * t752 + t967 * t975 + t777, t396 (-t1002 * t572 - t362 * t366 - t437 * t650) * qJD(3) + t3 * qJD(4) + t9 * qJD(5) + t18 * qJD(6) + t664 (-t171 * t572 - t362 * t691 + t437 * t937) * qJD(3) + t4 * qJD(4) + t10 * qJD(5) + t17 * qJD(6) + t659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t605, -t604, t568 * t643, t617, t549, qJD(2) * t452 + qJD(3) * t324 - qJD(4) * t402 + t795, qJD(2) * t451 + qJD(3) * t323 + qJD(4) * t401 - t796, -t126, t87, t128, t127, t549, qJD(3) * t58 + qJD(4) * t222 + qJD(5) * t96 + t246 + t801, qJD(3) * t59 - qJD(4) * t223 + qJD(5) * t98 + t245 - t800, t46, t11, t41, t40, t549, qJD(3) * t3 + qJD(4) * t66 + qJD(5) * t29 + qJD(6) * t22 + t107 + t886, qJD(3) * t4 - qJD(4) * t67 + qJD(5) * t30 + qJD(6) * t21 + t108 - t885; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t87, t128, t127, t549, qJD(3) * t94 + qJD(4) * t96 - qJD(5) * t214 + t246 + t798, qJD(3) * t93 + qJD(4) * t98 + qJD(5) * t213 + t245 - t799, t46, t11, t41, t40, t549, qJD(3) * t9 + qJD(4) * t29 + qJD(5) * t64 + qJD(6) * t24 + t107 + t888, qJD(3) * t10 + qJD(4) * t30 - qJD(5) * t65 + qJD(6) * t26 + t108 - t887; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t621, t648, t227 * t750 - t302 * t685 + t981, t229 * t750 + t685 * t938 + t986, t549, qJD(2) * t116 + qJD(3) * t18 + qJD(4) * t22 + qJD(5) * t24 - qJD(6) * t61 + t883, qJD(3) * t17 + qJD(4) * t21 + qJD(5) * t26 + qJD(6) * t60 + t108 - t884; 0, 0, 0, 0, -qJD(1), -t802, 0, 0, 0, 0, 0, -t755, -t753, 0, 0, 0, 0, 0, -qJD(4) * t450 - t543, -qJD(4) * t449 + t541, 0, 0, 0, 0, 0, -t249 * t750 + t718, -t250 * t750 + t717, 0, 0, 0, 0, 0, -t110 * qJD(6) - t109 * t750 + t62 - t720, -t963 * t967 - t139 - t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t754, -t752, 0, 0, 0, 0, 0, -t571 * t754 - t723, t542 - t721, 0, 0, 0, 0, 0, t342 * t750 + t485 * t754, t344 * t750 + t623 * t754, 0, 0, 0, 0, 0, -t366 * t754 + t967 * t973 + t816, -t691 * t754 + t966 * t967 - t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t722 - t726 - t757, t724 - t725 - t758, 0, 0, 0, 0, 0, t102, t101, 0, 0, 0, 0, 0, t31, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t101, 0, 0, 0, 0, 0, t31, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6) * t937 + t228 * t750 - t781 + t988, qJD(6) * t650 + t226 * t750 + t991; 0, 0, 0, 0, 0, 0, t715, -t756, 0, 0, 0, -t728, t729, t564 * t715 - t469, t617 * t944, t722 + t792, -t724 - t793, -t483, qJD(4) * t385 - t769, qJD(4) * t384 + t768, -t448 * t794 + t808, -t770 + t811, -t341 * t750 - t763, -t339 * t750 - t764, -t441, -qJD(4) * t57 - qJD(5) * t91 - t865, -qJD(4) * t56 - qJD(5) * t92 + t864, -t797 * t937 + t745, t746 - t815, -t962 * t967 - t776, -t967 * t976 - t777, -t396, -qJD(4) * t1 - qJD(5) * t8 - qJD(6) * t15 - t664, -qJD(4) * t2 - qJD(5) * t7 - qJD(6) * t16 - t659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t750 * t340, -t750 * t338, 0, 0, 0, 0, 0, -t967 * t974 - t816, -t965 * t967 + t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, t530 * qJD(4), 0, 0, 0, -pkin(3) * t785, -pkin(3) * t784, -t485 * t687, t750 * t307, 0, 0, 0, qJD(4) * t336 + t623 * t783, qJD(4) * t337 - t485 * t783 -(t688 + t760) * t691, t967 * t960, 0, 0, 0, qJD(4) * t123 - qJD(5) * t129 - t437 * t759, qJD(4) * t124 - qJD(5) * t130 + t437 * t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t616, -t603, t543 + t784, -t541 - t785, t551, -pkin(8) * t784 - t632, pkin(8) * t785 - t633, t150, t99, t216, t217, t551, -qJD(4) * t394 + qJD(5) * t314 - t932, t645 + t927, t50, t13, t55, t54, t551, -t928 - t1005, -t930 + t1004; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t99, t216, t217, t551, qJD(4) * t314 - qJD(5) * t394 - t933, -t619 + t927, t50, t13, t55, t54, t551, -t1005 - t929, t1004 - t931; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t620, t647, t258 * t750 + t760 - t978, t259 * t750 + t759 - t985, t551, -t1005 - t591, t1004 - t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, t604 (t568 * t753 - t786) * t569, t571 * t715 + t542, t549, qJD(2) * t450 - qJD(3) * t385 - t795, qJD(2) * t449 - qJD(3) * t384 + t796, t126, -t87, t242, t241, t549, qJD(3) * t57 + qJD(5) * t95 + t243 - t801, qJD(3) * t56 + qJD(5) * t97 + t244 + t800, -t46, -t11, t69, t68, t549, qJD(3) * t1 - qJD(5) * t27 - qJD(6) * t19 + t105 - t886, qJD(3) * t2 - qJD(5) * t28 - qJD(6) * t20 + t106 + t885; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t757, t758, 0, 0, 0, 0, 0, t134, t135, 0, 0, 0, 0, 0, t44, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t616, t603, -t543, t541, t550, t632, t633, -t150, -t99, t320, t318, t550, qJD(5) * t390 + t932, t762 - t645, -t50, -t13, t978, t985, t550, t928, t930; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t894, -pkin(4) * t697, 0, 0, 0, 0, 0, -t467 - t443, -t468 + t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t750 * t894 + t644, t814 + (-t698 - t697) * pkin(4), 0, 0, 0, 0, 0, qJD(6) * t380 - t467 - t596, qJD(6) * t381 - t468 - t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t380 - t443 - t934, qJD(5) * t381 + t442 - t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t87, t242, t241, t549, qJD(3) * t91 - qJD(4) * t95 + t243 - t798, qJD(3) * t92 - qJD(4) * t97 + t244 + t799, -t46, -t11, t69, t68, t549, qJD(3) * t8 + qJD(4) * t27 + qJD(6) * t23 + t105 - t888, qJD(3) * t7 + qJD(4) * t28 + qJD(6) * t25 + t106 + t887; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t135, 0, 0, 0, 0, 0, t44, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t99, t320, t318, t550, -qJD(4) * t390 + t933, t762 + t619, -t50, -t13, t978, t985, t550, t929, t931; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t894 - t644, pkin(4) * t698 - t814, 0, 0, 0, 0, 0, -qJD(6) * t399 + t596, -qJD(6) * t400 + t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t566 * t889, -t570 * t889; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t566 * t694 + t935, -t570 * t694 + t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t621, -t648, t69, t68, t549, qJD(2) * t110 + qJD(3) * t15 + qJD(4) * t19 - qJD(5) * t23 - t883, qJD(3) * t16 + qJD(4) * t20 - qJD(5) * t25 + t106 + t884; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t283 * t750 + t781 + t987, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t620, -t647, t978, t985, t550, t591, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t399 + t934, qJD(5) * t400 + t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t566 * t890 - t935, t570 * t890 - t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t12;
