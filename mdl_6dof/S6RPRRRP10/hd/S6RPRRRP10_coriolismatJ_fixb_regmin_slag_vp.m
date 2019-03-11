% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRRP10
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
% cmat_reg [(6*%NQJ)%x31]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRP10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:36
% EndTime: 2019-03-09 06:33:05
% DurationCPUTime: 16.62s
% Computational Cost: add. (12632->783), mult. (25252->954), div. (0->0), fcn. (25336->6), ass. (0->574)
t558 = cos(qJ(3));
t867 = t558 * pkin(3);
t556 = sin(qJ(3));
t869 = t556 * pkin(8);
t495 = t867 + t869;
t557 = cos(qJ(4));
t480 = t557 * t495;
t555 = sin(qJ(4));
t559 = -pkin(1) - pkin(7);
t791 = t555 * t559;
t668 = pkin(4) - t791;
t325 = t557 * t556 * pkin(9) + t558 * t668 + t480;
t554 = sin(qJ(5));
t306 = t554 * t325;
t478 = t555 * t495;
t782 = t558 * t559;
t500 = t557 * t782;
t793 = t555 * t556;
t365 = pkin(9) * t793 + t478 + t500;
t877 = cos(qJ(5));
t340 = t877 * t365;
t770 = t306 / 0.2e1 + t340 / 0.2e1;
t792 = t555 * t558;
t498 = t554 * t792;
t537 = t877 * t557;
t440 = t537 * t558 - t498;
t807 = t440 * qJ(6);
t697 = t877 * t555;
t794 = t554 * t557;
t475 = t697 + t794;
t783 = t558 * t475;
t874 = t783 * pkin(5);
t640 = -t807 + t874;
t876 = pkin(4) * t555;
t667 = -t559 + t876;
t642 = t667 * t558;
t223 = t642 + t640;
t875 = pkin(4) * t557;
t540 = -pkin(3) - t875;
t803 = t475 * qJ(6);
t795 = t554 * t555;
t471 = -t537 + t795;
t872 = t471 * pkin(5);
t639 = -t803 + t872;
t309 = t540 + t639;
t675 = -t783 / 0.2e1;
t890 = -t471 / 0.2e1;
t912 = t223 * t890 + t309 * t675;
t907 = -pkin(9) - pkin(8);
t496 = t907 * t557;
t665 = -t554 * t496 - t907 * t697;
t937 = t665 * t556;
t947 = t937 / 0.2e1;
t951 = t770 + t912 + t947;
t800 = t540 * t783;
t948 = -t937 / 0.2e1;
t950 = t800 / 0.2e1 - t770 + t948;
t949 = t770 - t912 + t948;
t618 = -t642 / 0.2e1;
t68 = t947 + t471 * t618 - t800 / 0.2e1 - t770;
t717 = qJD(4) + qJD(5);
t649 = t558 * t907 + qJ(2);
t609 = t649 * t557;
t868 = t557 * pkin(3);
t622 = t668 + t868;
t568 = t556 * t622 + t609;
t291 = t877 * t568;
t870 = t556 * pkin(3);
t648 = -pkin(8) * t558 + t870;
t621 = qJ(2) + t648;
t788 = t557 * t559;
t705 = t556 * t788;
t396 = t555 * t621 + t705;
t353 = -pkin(9) * t792 + t396;
t797 = t554 * t353;
t170 = -t291 + t797;
t164 = -pkin(5) * t556 + t170;
t541 = t556 * t559;
t706 = t555 * t541;
t352 = -t706 + (t649 + t870) * t557;
t698 = t877 * t352;
t183 = t698 - t797;
t946 = t164 + t183;
t576 = -t496 * t877 + t795 * t907;
t925 = t576 * t556;
t350 = -t925 / 0.2e1;
t887 = t475 / 0.2e1;
t891 = t440 / 0.2e1;
t604 = t223 * t887 + t309 * t891;
t923 = t350 + t604;
t324 = t877 * t353;
t565 = t554 * t568;
t171 = t324 + t565;
t562 = -t565 / 0.2e1 - t324 / 0.2e1;
t871 = t554 * pkin(4);
t709 = -t871 / 0.2e1;
t561 = t556 * t709 + t562;
t798 = t554 * t352;
t644 = t324 / 0.2e1 + t798 / 0.2e1;
t83 = t561 - t644;
t785 = t558 * t440;
t497 = t554 * t793;
t439 = t537 * t556 - t497;
t808 = t439 * t556;
t598 = t808 / 0.2e1 + t785 / 0.2e1;
t643 = t537 / 0.2e1 - t795 / 0.2e1;
t922 = t643 - t598;
t931 = t922 * qJD(2);
t945 = t83 * qJD(4) - t171 * qJD(5) + t931;
t182 = t324 + t798;
t944 = -t182 * qJD(4) + t83 * qJD(5) + t931;
t699 = t877 * t325;
t796 = t554 * t365;
t768 = t699 / 0.2e1 - t796 / 0.2e1;
t866 = t558 * pkin(5);
t641 = -t866 / 0.2e1 - t768;
t551 = t556 ^ 2;
t943 = t551 / 0.2e1;
t942 = -t665 / 0.2e1;
t941 = t665 / 0.2e1;
t894 = t783 / 0.2e1;
t349 = t925 / 0.2e1;
t82 = t561 + t644;
t572 = t598 + t643;
t915 = t572 * qJD(2);
t940 = t82 * qJD(5) + t915;
t939 = -t82 * qJD(4) + t915;
t938 = t439 * t717;
t518 = t697 / 0.2e1;
t582 = t794 / 0.2e1 + t518;
t296 = (t887 - t582) * t558;
t917 = qJD(1) * t572;
t936 = t296 * qJD(3) + t917;
t910 = t440 ^ 2;
t397 = t551 + t910;
t663 = t717 * t556;
t740 = qJD(3) * t475;
t934 = qJD(1) * t397 + t440 * t740 + t663;
t302 = -t558 * t582 + t675;
t933 = t302 * qJD(3) - t917 - t938;
t730 = t296 * qJD(2);
t924 = t717 * t576;
t932 = -t730 - t924;
t723 = t556 * qJD(1);
t693 = t471 * t723;
t930 = t717 * t922 - t693;
t929 = t717 * t665;
t928 = 0.2e1 * t555;
t927 = -t576 / 0.2e1;
t436 = t475 * t556;
t926 = t436 * t717;
t617 = t642 / 0.2e1;
t67 = t471 * t617 + t950;
t784 = t558 * t471;
t677 = -t784 / 0.2e1;
t883 = -t537 / 0.2e1;
t765 = t498 / 0.2e1 + t558 * t883;
t299 = t677 + t765;
t728 = t299 * qJD(2);
t739 = qJD(3) * t540;
t921 = qJD(1) * t67 + t471 * t739 + t728;
t293 = -t471 * t540 + t475 * t876;
t710 = -t875 / 0.2e1;
t889 = t471 / 0.2e1;
t892 = -t440 / 0.2e1;
t49 = t876 * t892 + (t475 * t710 + t667 * t889 + t709) * t558 + t950;
t920 = qJD(1) * t49 - qJD(3) * t293 + t728;
t354 = pkin(5) * t475 + qJ(6) * t471;
t827 = t309 * t471;
t137 = -t354 * t475 + t827;
t256 = pkin(5) * t440 + qJ(6) * t783;
t888 = -t475 / 0.2e1;
t601 = t256 * t888 + t354 * t892;
t787 = t558 * qJ(6);
t27 = -t601 + t787 + t951;
t676 = t784 / 0.2e1;
t715 = t877 / 0.2e1;
t653 = t558 * t715;
t764 = -t498 / 0.2e1 + t557 * t653;
t298 = t676 + t764;
t729 = t298 * qJD(2);
t919 = qJD(1) * t27 - qJD(3) * t137 + t729;
t318 = t354 + t876;
t125 = -t318 * t475 + t827;
t789 = t557 * t558;
t716 = pkin(4) * t789;
t228 = t256 + t716;
t603 = t228 * t887 + t318 * t891;
t530 = qJ(6) + t871;
t884 = t530 / 0.2e1;
t906 = qJ(6) / 0.2e1;
t647 = (t884 + t906) * t558;
t23 = t647 + t603 + t951;
t918 = qJD(1) * t23 - qJD(3) * t125 + t729;
t713 = t877 * pkin(4);
t660 = t713 / 0.2e1;
t588 = t291 / 0.2e1 + t556 * t660;
t85 = -t698 / 0.2e1 + t588;
t707 = t85 * qJD(1) + qJD(4) * t713;
t914 = t576 * qJD(6);
t812 = t783 * t665;
t819 = t576 * t440;
t913 = -t812 / 0.2e1 - t819 / 0.2e1;
t786 = t558 * t783;
t813 = t436 * t556;
t654 = t813 / 0.2e1 + t786 / 0.2e1;
t911 = -t572 * t717 + t693;
t550 = t555 ^ 2;
t552 = t557 ^ 2;
t512 = t552 - t550;
t666 = t789 * t928;
t579 = qJD(1) * t666 - qJD(3) * t512;
t470 = t475 ^ 2;
t909 = -pkin(5) / 0.2e1;
t908 = pkin(5) / 0.2e1;
t163 = t324 + t554 * t609 + (t554 * t622 + qJ(6)) * t556;
t905 = t163 / 0.2e1;
t904 = -t223 / 0.2e1;
t903 = -t309 / 0.2e1;
t896 = t576 / 0.2e1;
t895 = t436 / 0.2e1;
t893 = t439 / 0.2e1;
t886 = t497 / 0.2e1;
t885 = -t530 / 0.2e1;
t539 = -t713 - pkin(5);
t882 = t539 / 0.2e1;
t881 = t554 / 0.2e1;
t879 = t556 / 0.2e1;
t878 = -t558 / 0.2e1;
t873 = t439 * pkin(5);
t863 = pkin(4) * qJD(5);
t814 = t436 * qJ(6);
t613 = t814 / 0.2e1 + t873 / 0.2e1;
t672 = t170 / 0.2e1 - t164 / 0.2e1;
t674 = t905 - t171 / 0.2e1;
t7 = -t471 * t672 + t475 * t674 + t613;
t862 = t7 * qJD(1);
t831 = t223 * t440;
t838 = t182 * t556;
t43 = t228 * t783 + t831 - t838;
t857 = qJD(1) * t43;
t832 = t223 * t783;
t837 = t183 * t556;
t44 = -t228 * t440 + t832 + t837;
t856 = qJD(1) * t44;
t841 = t170 * t556;
t45 = -t256 * t440 + t832 - t841;
t855 = qJD(1) * t45;
t840 = t171 * t556;
t46 = t256 * t783 + t831 - t840;
t854 = qJD(1) * t46;
t560 = -t562 + (qJ(6) + t530) * t879;
t64 = -t560 + t644;
t853 = qJD(1) * t64;
t74 = t163 * t556 - t831;
t852 = qJD(1) * t74;
t851 = qJD(1) * t82;
t850 = qJD(4) * t85;
t849 = qJD(5) * t85;
t802 = t475 * t530;
t804 = t471 * t539;
t596 = t804 / 0.2e1 + t802 / 0.2e1;
t671 = t183 / 0.2e1 + t164 / 0.2e1;
t673 = t905 - t182 / 0.2e1;
t829 = t228 * t558;
t13 = t829 / 0.2e1 - t671 * t439 + t673 * t436 + t596;
t848 = t13 * qJD(1);
t844 = t164 * t439;
t846 = t163 * t436;
t778 = -t846 / 0.2e1 + t844 / 0.2e1;
t571 = -t170 * t439 / 0.2e1 + t171 * t895 + t256 * t878 + t778;
t611 = t872 / 0.2e1 - t803 / 0.2e1;
t16 = t571 + t611;
t847 = t16 * qJD(1);
t845 = t163 * t475;
t843 = t164 * t471;
t766 = t340 + t306;
t167 = t766 + t787;
t595 = t699 - t796;
t169 = -t595 - t866;
t469 = -pkin(4) * t793 + t541;
t222 = -pkin(5) * t436 + qJ(6) * t439 + t469;
t17 = t163 * t167 + t164 * t169 + t222 * t223;
t842 = t17 * qJD(1);
t18 = -t163 * t170 + t164 * t171 + t223 * t256;
t839 = t18 * qJD(1);
t19 = t163 * t183 + t164 * t182 + t223 * t228;
t836 = t19 * qJD(1);
t20 = -t167 * t783 + t169 * t440 - t844 + t846;
t835 = t20 * qJD(1);
t21 = (-t163 + t171) * t440 + (-t164 + t170) * t783;
t834 = t21 * qJD(1);
t22 = (-t163 + t182) * t440 - t946 * t783;
t833 = t22 * qJD(1);
t826 = t309 * t475;
t31 = t163 * t558 + t167 * t556 - t222 * t440 + t223 * t439;
t825 = t31 * qJD(1);
t32 = -t164 * t558 - t169 * t556 + t222 * t783 - t223 * t436;
t824 = t32 * qJD(1);
t821 = t665 * t558;
t817 = t576 * t558;
t41 = t595 * t556 + t469 * t783 + (-t436 * t667 - t170) * t558;
t816 = t41 * qJD(1);
t42 = t171 * t558 + t439 * t642 - t469 * t440 + t556 * t766;
t815 = t42 * qJD(1);
t811 = t783 * t475;
t810 = t783 * t539;
t809 = t439 * t539;
t806 = t440 * t530;
t805 = t471 * t440;
t801 = t530 * t436;
t799 = t540 * t440;
t553 = t558 ^ 2;
t538 = t553 * t557;
t790 = t557 * t551;
t57 = t843 + t845;
t781 = t57 * qJD(1);
t123 = -t475 * t440 + t471 * t783;
t779 = t717 * t123;
t777 = t439 * qJD(6);
t186 = -t811 / 0.2e1 - t805 / 0.2e1;
t776 = t717 * t186;
t205 = t654 + t582;
t545 = t556 * qJD(6);
t775 = t205 * qJD(2) + t545;
t207 = -t654 + t582;
t774 = t207 * qJD(2) + t545;
t771 = t717 * t296;
t679 = t789 / 0.2e1;
t763 = t558 * t518 + t554 * t679;
t762 = t795 / 0.2e1 + t883;
t761 = -t794 / 0.2e1 - t697 / 0.2e1;
t511 = t551 - t553;
t760 = qJ(6) * qJD(5);
t759 = qJD(1) * qJ(2);
t101 = -t838 + (t440 * t667 + t783 * t875) * t558;
t758 = qJD(1) * t101;
t338 = t783 * t642;
t102 = -t440 * t716 + t338 + t837;
t757 = qJD(1) * t102;
t103 = t338 - t841;
t756 = qJD(1) * t103;
t104 = t440 * t642 - t840;
t755 = qJD(1) * t104;
t187 = t805 - t811;
t754 = qJD(1) * t187;
t202 = t762 - t598;
t753 = qJD(1) * t202;
t750 = qJD(1) * t205;
t206 = -t654 + t761;
t749 = qJD(1) * t206;
t250 = -t786 + t813;
t748 = qJD(1) * t250;
t251 = t785 - t808;
t747 = qJD(1) * t251;
t395 = -t557 * t621 + t706;
t307 = -t395 * t556 - t553 * t791;
t746 = qJD(1) * t307;
t308 = -t396 * t556 - t553 * t788;
t745 = qJD(1) * t308;
t744 = qJD(1) * t440;
t474 = t511 * t555;
t743 = qJD(1) * t474;
t477 = -t538 + t790;
t742 = qJD(1) * t477;
t741 = qJD(2) * t556;
t738 = qJD(3) * t557;
t737 = qJD(4) * t555;
t736 = qJD(4) * t557;
t735 = qJD(5) * t540;
t734 = qJD(6) * t475;
t166 = t436 * t440 + t439 * t783;
t733 = t166 * qJD(1);
t704 = t555 * t782;
t188 = -t395 * t558 + (t480 + t704) * t556;
t732 = t188 * qJD(1);
t189 = t396 * t558 + (-t500 + t478) * t556;
t731 = t189 * qJD(1);
t294 = (t887 + t582) * t556;
t263 = t294 * qJD(1);
t297 = t886 + (t883 + t889) * t556;
t266 = t297 * qJD(1);
t727 = t783 * qJD(6);
t714 = 0.1e1 / 0.2e1 + t943;
t442 = (-t553 / 0.2e1 - t714) * t555;
t726 = t442 * qJD(1);
t443 = t538 / 0.2e1 + t714 * t557;
t725 = t443 * qJD(1);
t459 = t471 * qJD(6);
t724 = t511 * qJD(1);
t722 = t556 * qJD(3);
t721 = t558 * qJD(1);
t720 = t558 * qJD(3);
t719 = t558 * qJD(4);
t544 = qJD(5) * t713;
t718 = t544 + qJD(6);
t712 = t554 * t863;
t711 = t876 / 0.2e1;
t701 = t302 * t717 + t471 * t722;
t696 = qJ(2) * t723;
t695 = qJ(2) * t721;
t694 = t783 * t744;
t692 = t471 * t741;
t691 = t475 * t741;
t690 = t471 * t740;
t689 = t475 * t722;
t688 = t555 * t738;
t687 = t555 * t720;
t686 = t557 * t720;
t685 = t556 * t737;
t684 = t555 * t719;
t683 = t556 * t736;
t682 = t557 * t719;
t681 = t555 * t736;
t528 = t556 * t720;
t680 = t556 * t721;
t678 = t354 * t878;
t670 = t942 + t941;
t669 = t896 + t927;
t664 = t717 * t475;
t662 = pkin(4) * t679;
t661 = -qJD(4) - t723;
t659 = -t845 / 0.2e1 - t843 / 0.2e1 + t913;
t656 = t475 * t617 + t799 / 0.2e1 + t768;
t651 = t768 + t866;
t650 = t723 + qJD(4) / 0.2e1;
t646 = qJD(3) * t666;
t599 = t809 / 0.2e1 - t801 / 0.2e1;
t5 = t440 * t669 + t471 * t671 + t475 * t673 + t670 * t783 - t599;
t638 = t5 * qJD(1);
t501 = qJD(5) - t661;
t636 = -t799 / 0.2e1 + t768;
t563 = t163 * t891 + t164 * t894 + t167 * t893 + t169 * t895 + t222 * t878 + t223 * t879;
t597 = t576 * t888 + t665 * t890;
t10 = t563 + t597;
t155 = t436 * t783 + t439 * t440 - t556 * t558;
t635 = t10 * qJD(1) + t155 * qJD(2);
t587 = t660 + t882 + t908;
t610 = t885 + t871 / 0.2e1 + t906;
t106 = -t471 * t587 + t475 * t610;
t93 = t440 * t610 - t587 * t783;
t634 = qJD(1) * t93 + qJD(3) * t106;
t124 = t318 * t471 + t826;
t570 = t228 * t889 + t318 * t894 + t923;
t25 = (t882 + t909) * t558 + t570 - t768;
t633 = -qJD(1) * t25 - qJD(3) * t124;
t136 = t354 * t471 + t826;
t583 = t349 - t604;
t602 = t256 * t889 + t354 * t894;
t29 = t583 - t602 + t651;
t631 = qJD(1) * t29 - qJD(3) * t136;
t292 = t471 * t876 + t475 * t540;
t586 = -t783 * t711 + t349;
t50 = (t471 * t710 + t667 * t888 + t660) * t558 + t586 + t636;
t629 = qJD(1) * t50 - qJD(3) * t292;
t84 = -t797 + t698 / 0.2e1 + t588;
t627 = -qJD(4) * t84 + qJD(5) * t170;
t626 = -qJD(4) * t183 - qJD(5) * t84;
t625 = t661 * t558;
t214 = t783 ^ 2 - t910;
t62 = qJD(1) * t214 + qJD(3) * t123;
t257 = t471 ^ 2 - t470;
t86 = qJD(1) * t123 + qJD(3) * t257;
t620 = t869 / 0.2e1 + t867 / 0.2e1;
t616 = t167 * t906 + t169 * t909;
t615 = t182 * t908 - t183 * qJ(6) / 0.2e1;
t614 = qJ(6) * t941 + t576 * t908;
t612 = -t874 / 0.2e1 + t807 / 0.2e1;
t584 = t620 * t555;
t361 = t478 / 0.2e1 + t584;
t608 = pkin(3) * t738 - qJD(1) * t361;
t585 = t620 * t557;
t362 = -t480 / 0.2e1 - t585;
t607 = pkin(3) * qJD(3) * t555 - qJD(1) * t362;
t606 = t167 * t884 + t169 * t882;
t605 = t256 * t903 + t354 * t904;
t600 = -t810 / 0.2e1 - t806 / 0.2e1;
t66 = t475 * t618 + t349 + t636;
t593 = qJD(1) * t66 - t475 * t739;
t592 = t557 * t625;
t177 = t878 - t186;
t591 = qJD(1) * t177 + t690;
t107 = qJD(3) * t186 - t694;
t131 = -qJD(1) * t186 + t690;
t460 = (t550 / 0.2e1 - t552 / 0.2e1) * t558;
t589 = -qJD(1) * t460 + t688;
t581 = qJD(1) * t538 * t555 + qJD(3) * t460;
t473 = t512 * t553;
t580 = qJD(1) * t473 + t646;
t564 = t163 * t941 + t182 * t942 + t228 * t903 + t318 * t904 + t927 * t946;
t1 = t564 + t606;
t569 = t318 * t878 - t436 * t669 + t439 * t670;
t36 = t569 + t600;
t58 = t309 * t318;
t578 = -t1 * qJD(1) + t36 * qJD(2) + t58 * qJD(3);
t3 = t576 * t672 + t665 * t674 + t605 + t616;
t38 = t678 - t612;
t59 = t309 * t354;
t577 = -t3 * qJD(1) + t38 * qJD(2) + t59 * qJD(3);
t295 = t675 + t763;
t39 = t641 + t923;
t574 = qJD(1) * t39 + qJD(2) * t295 + t309 * t740;
t567 = (t163 * t715 + t164 * t881) * pkin(4) + t170 * t885 + t171 * t882;
t12 = t567 + t615;
t418 = (t530 * t877 + t539 * t554) * pkin(4);
t566 = (t576 * t715 + t665 * t881) * pkin(4) + t665 * t885 + t576 * t882;
t56 = t566 + t614;
t95 = t436 * t610 + t439 * t587;
t573 = t12 * qJD(1) + t95 * qJD(2) + t56 * qJD(3) + t418 * qJD(4);
t549 = qJ(6) * qJD(6);
t534 = -t721 / 0.2e1;
t533 = t721 / 0.2e1;
t532 = t720 / 0.2e1;
t527 = t557 * t723;
t526 = t555 * t722;
t525 = t555 * t723;
t510 = t530 * qJD(6);
t468 = t501 * qJ(6);
t467 = t650 * t558;
t454 = t460 * qJD(4);
t452 = t475 * t723;
t445 = -t790 / 0.2e1 - t538 / 0.2e1 + t557 / 0.2e1;
t444 = (-0.1e1 / 0.2e1 + t943 + t553 / 0.2e1) * t555;
t433 = (qJD(5) / 0.2e1 + t650) * t558;
t336 = t440 * t734;
t305 = -t436 / 0.2e1 + t582 * t556;
t304 = t886 + (t883 + t890) * t556;
t303 = t894 + t763;
t301 = t676 + t765;
t300 = t677 + t764;
t286 = -t704 + t480 / 0.2e1 - t585;
t285 = -t500 - t478 / 0.2e1 + t584;
t227 = qJD(3) * t470 + t475 * t744;
t209 = t598 + t762;
t208 = t654 + t761;
t199 = t208 * qJD(2);
t197 = t206 * qJD(2);
t193 = qJD(3) * t294 + t440 * t723;
t192 = qJD(3) * t297 + t723 * t783;
t178 = t878 + t186;
t176 = -t664 - t263;
t175 = -t471 * t717 - t266;
t118 = qJD(3) * t299 + t749;
t117 = qJD(3) * t298 + t750;
t114 = qJD(3) * t305 - t440 * t501;
t113 = qJD(3) * t304 - t501 * t783;
t105 = (t475 * t881 + t877 * t890) * pkin(4) - t596 + t611;
t94 = (t436 * t881 + t439 * t715) * pkin(4) + t599 - t613;
t92 = (t440 * t881 - t715 * t783) * pkin(4) + t600 - t612;
t89 = qJD(3) * t301 - t749 + t926;
t88 = qJD(3) * t300 - t750 - t926;
t69 = t350 + t656;
t65 = t560 + t644;
t55 = t566 - t614;
t54 = qJD(4) * t871 - t851;
t53 = -t717 * t871 + t851;
t52 = t440 * t711 + t475 * t662 + t558 * t709 + t68;
t51 = pkin(4) * t653 + t471 * t662 - t586 + t656;
t40 = t583 + t641;
t37 = t678 + t612;
t35 = t569 - t600;
t30 = t602 + t651 + t923;
t28 = t601 + t787 + t949;
t26 = t539 * t878 + t570 - t641;
t24 = t647 - t603 + t949;
t15 = t571 - t611;
t14 = t183 * t893 - t829 / 0.2e1 + t182 * t895 + t596 + t778;
t11 = t567 - t615;
t9 = t563 - t597;
t8 = t170 * t889 + t171 * t887 + t613 + t659 - t913;
t6 = t182 * t887 + t183 * t890 + t576 * t891 - t783 * t942 - t599 + t659;
t4 = t163 * t942 + t164 * t896 + t170 * t927 + t171 * t941 - t605 + t616;
t2 = -t564 + t606;
t33 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t528, t511 * qJD(3), 0, 0, 0, qJ(2) * t720 + t741, -qJ(2) * t722 + qJD(2) * t558, -t528 * t552 - t553 * t681, -qJD(4) * t473 + t556 * t646, -qJD(3) * t477 - t556 * t684, qJD(3) * t474 - t556 * t682, t528, qJD(3) * t188 + qJD(4) * t308 + t557 * t741, -qJD(3) * t189 - qJD(4) * t307 - t555 * t741 (-qJD(3) * t439 - t717 * t783) * t440, qJD(3) * t166 + t214 * t717, qJD(3) * t251 - t663 * t783, qJD(3) * t250 - t440 * t663, t528, qJD(3) * t41 + qJD(4) * t101 + qJD(5) * t104 - t692, -qJD(3) * t42 - qJD(4) * t102 - qJD(5) * t103 - t691, qJD(3) * t32 + qJD(4) * t43 + qJD(5) * t46 - t440 * t727 - t692, qJD(2) * t187 + qJD(3) * t20 + qJD(4) * t22 + qJD(5) * t21 - t545 * t783, qJD(3) * t31 + qJD(4) * t44 + qJD(5) * t45 + qJD(6) * t397 + t691, qJD(2) * t57 + qJD(3) * t17 + qJD(4) * t19 + qJD(5) * t18 + qJD(6) * t74; 0, 0, 0, 0, qJD(1), t759, 0, 0, 0, 0, 0, t723, t721, 0, 0, 0, 0, 0, qJD(4) * t445 + t527, qJD(4) * t444 - t525, 0, 0, 0, 0, 0, t930, t208 * t717 - t452, t930, t754, t207 * t717 + t452, t781 + (t471 * t436 + t475 * t439) * qJD(2) + t9 * qJD(3) + t14 * qJD(4) + t15 * qJD(5) + t209 * qJD(6); 0, 0, 0, 0, 0, 0, -t680, t724, -t722, -t720, 0, -t559 * t722 + t695, -t559 * t720 - t696, -t454 + (-t552 * t721 - t688) * t556, t556 * t579 - 0.2e1 * t558 * t681, t687 - t742, t686 + t743, t467, t732 + (t555 * t648 - t705) * qJD(3) + t286 * qJD(4), -t731 + (-pkin(8) * t789 + (t791 + t868) * t556) * qJD(3) + t285 * qJD(4) (-t740 - t744) * t439 + t776, t733 + (t436 * t475 + t439 * t471) * qJD(3) + t779, t304 * t717 + t475 * t720 + t747, t305 * t717 - t471 * t720 + t748, t433, t816 + (-t436 * t540 + t469 * t471 - t821) * qJD(3) + t51 * qJD(4) + t69 * qJD(5), -t815 + (-t439 * t540 + t469 * t475 - t817) * qJD(3) + t52 * qJD(4) + t68 * qJD(5), t824 + (t222 * t471 - t309 * t436 - t821) * qJD(3) + t26 * qJD(4) + t30 * qJD(5) + t178 * qJD(6), t835 + (-t167 * t471 + t169 * t475 + t436 * t576 - t439 * t665) * qJD(3) + t6 * qJD(4) + t8 * qJD(5) + t304 * qJD(6), t825 + (-t222 * t475 + t309 * t439 + t817) * qJD(3) + t24 * qJD(4) + t28 * qJD(5) + t336, t842 + t9 * qJD(2) + (t167 * t576 + t169 * t665 + t222 * t309) * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t40 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t581, -t580, t555 * t625, t592, t532, qJD(2) * t445 + qJD(3) * t286 - qJD(4) * t396 + t745, qJD(2) * t444 + qJD(3) * t285 + qJD(4) * t395 - t746, t107, t62, t113, t114, t532, qJD(3) * t51 + t758 + t944, qJD(3) * t52 + t199 + t626 - t757, qJD(3) * t26 + t857 + t944, t833 + t6 * qJD(3) + (-t806 - t810) * qJD(4) + t92 * qJD(5) - t727, qJD(3) * t24 - t626 + t774 + t856, t836 + t14 * qJD(2) + t2 * qJD(3) + (t182 * t539 + t183 * t530) * qJD(4) + t11 * qJD(5) + t65 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t62, t113, t114, t532, qJD(3) * t69 + t755 + t945, qJD(3) * t68 + t199 + t627 - t756, qJD(3) * t30 + t854 + t945, t8 * qJD(3) + t92 * qJD(4) + qJD(5) * t640 - t727 + t834, qJD(3) * t28 - t627 + t774 + t855, t839 + t15 * qJD(2) + t4 * qJD(3) + t11 * qJD(4) + (-pkin(5) * t171 - qJ(6) * t170) * qJD(5) + t163 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t178 - t694, t113, t934, qJD(2) * t209 + qJD(3) * t40 + qJD(4) * t65 + qJD(5) * t163 + t852; 0, 0, 0, 0, -qJD(1), -t759, 0, 0, 0, 0, 0, -t723, -t721, 0, 0, 0, 0, 0, -qJD(4) * t443 - t527, -qJD(4) * t442 + t525, 0, 0, 0, 0, 0, t911, -t206 * t717 + t452, t911, -t754, -t205 * t717 - t452, qJD(3) * t10 - qJD(4) * t13 + qJD(5) * t16 - qJD(6) * t202 - t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t722, -t720, 0, 0, 0, 0, 0, -t557 * t722 - t684, t526 - t682, 0, 0, 0, 0, 0, t701, t301 * t717 + t689, t701, -qJD(3) * t187, t300 * t717 - t689 (t309 * t556 + t812 + t819) * qJD(3) + t35 * qJD(4) + t37 * qJD(5) + t303 * qJD(6) + t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t683 - t687 - t725, t685 - t686 - t726, 0, 0, 0, 0, 0, t933, t89, t933, 0, t88, -t848 + t35 * qJD(3) + (-t801 + t809) * qJD(4) + t94 * qJD(5) + t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t933, t89, t933, 0, t88, t847 + t37 * qJD(3) + t94 * qJD(4) + (-t814 - t873) * qJD(5) + t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t303 - t753 + t938; 0, 0, 0, 0, 0, 0, t680, -t724, 0, 0, 0, -t695, t696, t552 * t680 - t454, t592 * t928, t683 + t742, -t685 - t743, -t467, qJD(4) * t362 - t732, qJD(4) * t361 + t731, t439 * t744 + t776, -t733 + t779, -t297 * t717 - t747, -t294 * t717 - t748, -t433, -qJD(4) * t50 - qJD(5) * t66 - t816, -qJD(4) * t49 - qJD(5) * t67 + t815, qJD(4) * t25 - qJD(5) * t29 - qJD(6) * t177 - t824, -qJD(4) * t5 - qJD(5) * t7 - qJD(6) * t297 - t835, -qJD(4) * t23 - qJD(5) * t27 + t336 - t825, -qJD(2) * t10 - qJD(4) * t1 - qJD(5) * t3 - qJD(6) * t39 - t842; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t771, -t717 * t299, -t771, 0, -t717 * t298, qJD(4) * t36 + qJD(5) * t38 - qJD(6) * t295 - t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t681, t512 * qJD(4), 0, 0, 0, -pkin(3) * t737, -pkin(3) * t736, -t471 * t664, t717 * t257, 0, 0, 0, qJD(4) * t292 + t475 * t735, qJD(4) * t293 - t471 * t735, qJD(4) * t124 + qJD(5) * t136 - t459 * t475, 0, qJD(4) * t125 + qJD(5) * t137 + qJD(6) * t470, qJD(4) * t58 + qJD(5) * t59 - t309 * t734; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t589, -t579, t527 + t736, -t525 - t737, t534, -pkin(8) * t736 - t607, pkin(8) * t737 - t608, -t131, t86, t175, t176, t534, -t629 + t932, t929 - t920, -t633 + t932 (-t802 - t804) * qJD(4) + t105 * qJD(5) - t638 - t459, -t929 - t918 (-t530 * t665 + t539 * t576) * qJD(4) + t55 * qJD(5) + t914 + t578; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t86, t175, t176, t534, -t593 + t932, t929 - t921, -t631 + t932, t105 * qJD(4) + qJD(5) * t639 - t459 - t862, -t929 - t919, t55 * qJD(4) + (-pkin(5) * t576 - qJ(6) * t665) * qJD(5) + t914 + t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t591, t175, t227, -t574 + t924; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t581, t580 (t555 * t721 - t738) * t556, t557 * t680 + t526, t532, qJD(2) * t443 - qJD(3) * t362 - t745, qJD(2) * t442 - qJD(3) * t361 + t746, -t107, -t62, t192, t193, t532, qJD(3) * t50 - t758 + t940, qJD(3) * t49 + t197 + t757 - t849, -qJD(3) * t25 - t857 + t940, qJD(3) * t5 + qJD(5) * t93 - t833, qJD(3) * t23 + t775 + t849 - t856, qJD(2) * t13 + qJD(3) * t1 + qJD(5) * t12 - qJD(6) * t64 - t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t725, t726, 0, 0, 0, 0, 0, t936, t118, t936, 0, t117, -qJD(3) * t36 + qJD(5) * t95 + t848; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t589, t579, -t527, t525, t533, t607, t608, t131, -t86, t266, t263, t533, t629 + t730, t920, t633 + t730, qJD(5) * t106 + t638, t918, qJD(5) * t56 - t578; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t712, -t544, -t712, 0, t718, qJD(5) * t418 + t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t544 - t707, t53, t634, t707 + t718 (-pkin(5) * t554 + qJ(6) * t877) * t863 + t510 + t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, t530 * t717 - t853; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t62, t192, t193, t532, qJD(3) * t66 - t755 + t939, qJD(3) * t67 + t197 + t756 + t850, qJD(3) * t29 - t854 + t939, qJD(3) * t7 - qJD(4) * t93 - t834, qJD(3) * t27 + t775 - t850 - t855, qJ(6) * t545 - qJD(2) * t16 + qJD(3) * t3 - qJD(4) * t12 - t839; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t936, t118, t936, 0, t117, -qJD(3) * t38 - qJD(4) * t95 - t847; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t86, t266, t263, t533, t593 + t730, t921, t631 + t730, -qJD(4) * t106 + t862, t919, -qJD(4) * t56 - t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t707, t54, -t634, qJD(6) - t707, t549 - t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t177 + t694, t192, -t934, qJD(2) * t202 + qJD(3) * t39 + qJD(4) * t64 - t556 * t760 - t852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t295 + t753; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t591, t266, -t227, t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501, -qJD(4) * t530 - t760 + t853; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501, -t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t33;
