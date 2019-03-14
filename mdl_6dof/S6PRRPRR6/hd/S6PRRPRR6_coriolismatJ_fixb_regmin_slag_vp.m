% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:41
% EndTime: 2019-03-08 22:28:26
% DurationCPUTime: 25.03s
% Computational Cost: add. (15127->862), mult. (41157->1290), div. (0->0), fcn. (48631->14), ass. (0->641)
t579 = sin(pkin(13));
t582 = cos(pkin(13));
t583 = cos(pkin(7));
t580 = sin(pkin(7));
t586 = sin(qJ(3));
t893 = t580 * t586;
t524 = -t579 * t893 + t582 * t583;
t525 = t579 * t583 + t582 * t893;
t585 = sin(qJ(5));
t929 = cos(qJ(5));
t410 = -t929 * t524 + t585 * t525;
t588 = cos(qJ(6));
t968 = t588 * t410;
t383 = -t968 / 0.2e1;
t974 = t968 / 0.2e1 + t383;
t978 = qJD(6) * t974;
t781 = t929 * t582;
t872 = t585 * t579;
t656 = t781 - t872;
t720 = t781 / 0.2e1;
t589 = cos(qJ(3));
t892 = t580 * t589;
t760 = -t892 / 0.2e1;
t724 = t585 * t760;
t372 = t579 * t724 - t656 * t760 + t720 * t892;
t834 = qJD(2) * t589;
t778 = t580 * t834;
t728 = t410 * t778;
t977 = t372 * qJD(3) - t728;
t782 = t929 * t579;
t871 = t585 * t582;
t538 = t782 + t871;
t371 = t582 * t724 + (t538 + t782) * t760;
t504 = t929 * t525;
t873 = t585 * t524;
t963 = t504 + t873;
t729 = t963 * t778;
t976 = t371 * qJD(3) - t729;
t832 = qJD(3) * t538;
t837 = qJD(2) * t963;
t975 = t832 + t837;
t803 = t963 * qJD(5);
t569 = pkin(9) * t892;
t570 = t583 * t586 * pkin(2);
t530 = t569 + t570;
t501 = qJ(4) * t583 + t530;
t710 = -pkin(3) * t589 - qJ(4) * t586;
t502 = (-pkin(2) + t710) * t580;
t362 = -t579 * t501 + t582 * t502;
t277 = -pkin(4) * t892 - t525 * pkin(10) + t362;
t363 = t582 * t501 + t579 * t502;
t296 = pkin(10) * t524 + t363;
t159 = -t929 * t277 + t585 * t296;
t152 = pkin(5) * t892 + t159;
t973 = t152 / 0.2e1;
t953 = t410 / 0.2e1;
t972 = t579 / 0.2e1;
t937 = -t582 / 0.2e1;
t951 = t963 / 0.2e1;
t584 = sin(qJ(6));
t881 = t584 * t963;
t342 = t588 * t892 + t881;
t971 = t342 * t963;
t859 = t588 * t963;
t344 = -t584 * t892 + t859;
t970 = t344 * t963;
t969 = t584 * t410;
t833 = qJD(3) * t656;
t838 = qJD(2) * t410;
t967 = -t833 + t838;
t912 = t410 * t538;
t940 = -t656 / 0.2e1;
t216 = t963 * t940 + t912 / 0.2e1;
t966 = t216 * qJD(5);
t804 = t410 * qJD(5);
t581 = sin(pkin(6));
t930 = cos(qJ(2));
t787 = t930 * t586;
t587 = sin(qJ(2));
t869 = t587 * t589;
t789 = t583 * t869;
t490 = (t787 + t789) * t581;
t786 = t930 * t589;
t870 = t586 * t587;
t491 = (-t583 * t870 + t786) * t581;
t891 = t581 * t587;
t791 = t580 * t891;
t423 = -t491 * t579 + t582 * t791;
t424 = t491 * t582 + t579 * t791;
t663 = t423 * t972 + t424 * t937;
t965 = -t490 * pkin(3) / 0.2e1 - t663 * qJ(4);
t777 = t656 * t832;
t654 = qJD(2) * t216 - t777;
t655 = qJD(3) * t216 + t410 * t837;
t828 = qJD(3) * t588;
t774 = t584 * t828;
t714 = 0.2e1 * t538 * t774;
t532 = t656 ^ 2;
t533 = t538 ^ 2;
t962 = -t533 - t532;
t798 = t533 - t532;
t577 = t584 ^ 2;
t578 = t588 ^ 2;
t568 = t578 - t577;
t736 = qJD(5) * t568;
t961 = t714 - t736;
t960 = t410 ^ 2;
t922 = cos(pkin(6));
t738 = t922 * t580;
t460 = t586 * t738 + (t583 * t787 + t869) * t581;
t788 = t930 * t581;
t633 = -t580 * t788 + t583 * t922;
t340 = t460 * t582 + t579 * t633;
t604 = -t460 * t579 + t582 * t633;
t203 = t340 * t585 - t604 * t929;
t959 = t203 / 0.2e1;
t958 = -t342 / 0.2e1;
t957 = t342 / 0.2e1;
t956 = -t344 / 0.2e1;
t955 = t344 / 0.2e1;
t954 = -t410 / 0.2e1;
t952 = -t963 / 0.2e1;
t562 = t588 * t893;
t481 = t656 * t892;
t877 = t584 * t481;
t426 = -t562 + t877;
t950 = -t426 / 0.2e1;
t730 = t583 * t788;
t790 = t581 * t870;
t459 = t790 + (-t730 - t738) * t589;
t949 = -t459 / 0.2e1;
t948 = t459 / 0.2e1;
t947 = t460 / 0.2e1;
t924 = t538 * pkin(5);
t925 = t656 * pkin(11);
t461 = t924 - t925;
t946 = t461 / 0.2e1;
t923 = pkin(10) + qJ(4);
t556 = t923 * t582;
t739 = t923 * t579;
t471 = t556 * t929 - t585 * t739;
t945 = -t471 / 0.2e1;
t480 = t538 * t892;
t944 = -t480 / 0.2e1;
t943 = -t524 / 0.2e1;
t942 = t525 / 0.2e1;
t941 = t656 / 0.2e1;
t939 = -t538 / 0.2e1;
t938 = t538 / 0.2e1;
t936 = t583 / 0.2e1;
t935 = t584 / 0.2e1;
t934 = -t586 / 0.2e1;
t933 = -t588 / 0.2e1;
t932 = t588 / 0.2e1;
t931 = t589 / 0.2e1;
t928 = pkin(2) * t589;
t927 = t410 * pkin(11);
t204 = t340 * t929 + t585 * t604;
t888 = t584 * t204;
t148 = -t459 * t588 + t888;
t921 = t148 * t656;
t866 = t588 * t204;
t149 = t459 * t584 + t866;
t920 = t149 * t656;
t784 = t929 * t423;
t874 = t585 * t424;
t246 = -t784 + t874;
t919 = t246 * t584;
t918 = t246 * t588;
t269 = t459 * t538;
t917 = t269 * t584;
t916 = t269 * t588;
t915 = t340 * t582;
t914 = t344 * t584;
t913 = t344 * t588;
t855 = t588 * t481;
t427 = t584 * t893 + t855;
t909 = t427 * t584;
t908 = t460 * t584;
t907 = t460 * t588;
t470 = t556 * t585 + t739 * t929;
t906 = t470 * t584;
t905 = t470 * t588;
t904 = t490 * t656;
t903 = t490 * t538;
t902 = t490 * t579;
t901 = t490 * t582;
t900 = t490 * t584;
t899 = t490 * t588;
t898 = t524 * t582;
t897 = t525 * t579;
t896 = t538 * t588;
t442 = t584 * t538;
t861 = t588 * t342;
t657 = t861 / 0.2e1 + t914 / 0.2e1;
t620 = t442 * t968 - t656 * t657;
t662 = t427 * t932 + t584 * t950;
t54 = t620 - t662;
t895 = t54 * qJD(2);
t575 = t580 ^ 2;
t894 = t575 * t586;
t890 = t584 * t159;
t527 = (pkin(3) * t586 - qJ(4) * t589) * t580;
t529 = pkin(9) * t893 - t583 * t928;
t417 = t582 * t527 + t579 * t529;
t339 = (-pkin(10) * t582 * t589 + pkin(4) * t586) * t580 + t417;
t309 = t585 * t339;
t418 = t579 * t527 - t582 * t529;
t792 = t579 * t892;
t359 = -pkin(10) * t792 + t418;
t354 = t929 * t359;
t850 = t354 + t309;
t202 = pkin(11) * t893 + t850;
t889 = t584 * t202;
t243 = pkin(5) * t963 + t927;
t887 = t584 * t243;
t783 = t929 * t424;
t875 = t585 * t423;
t247 = t783 + t875;
t886 = t584 * t247;
t479 = pkin(4) * t792 + t530;
t267 = pkin(5) * t480 - pkin(11) * t481 + t479;
t885 = t584 * t267;
t270 = t656 * t459;
t884 = t584 * t270;
t883 = t584 * t342;
t880 = t584 * t461;
t879 = t584 * t471;
t878 = t584 * t480;
t440 = t584 * t656;
t876 = t585 * t359;
t868 = t588 * t159;
t867 = t588 * t202;
t865 = t588 * t243;
t864 = t588 * t247;
t863 = t588 * t267;
t862 = t588 * t270;
t858 = t588 * t461;
t857 = t588 * t471;
t856 = t588 * t480;
t445 = t588 * t656;
t596 = t604 * t579;
t81 = (t460 + t596 - t915) * t459;
t854 = t81 * qJD(1);
t82 = t340 * t424 + t423 * t604 + t459 * t490;
t853 = t82 * qJD(1);
t694 = t861 + t914;
t99 = t694 * t410;
t852 = t99 * qJD(2);
t851 = -t152 + t159;
t721 = t787 / 0.2e1;
t750 = t869 / 0.2e1;
t723 = t581 * t750;
t849 = t581 * t721 + t583 * t723;
t722 = -t788 / 0.2e1;
t848 = -t581 * t789 / 0.2e1 + t586 * t722;
t565 = t579 ^ 2 + t582 ^ 2;
t103 = -t584 * t960 + t971;
t847 = qJD(2) * t103;
t104 = -t410 * t969 - t971;
t846 = qJD(2) * t104;
t105 = -t588 * t960 + t970;
t845 = qJD(2) * t105;
t106 = t410 * t968 + t970;
t844 = qJD(2) * t106;
t843 = qJD(2) * t969;
t742 = 0.2e1 * t953;
t227 = t742 * t584;
t842 = qJD(2) * t227;
t265 = t410 * t893 - t480 * t892;
t841 = qJD(2) * t265;
t266 = t481 * t892 - t893 * t963;
t840 = qJD(2) * t266;
t839 = qJD(2) * t344;
t531 = t782 / 0.2e1 + t871 / 0.2e1;
t478 = t531 * t892;
t836 = qJD(2) * t478;
t835 = qJD(2) * t583;
t573 = -pkin(4) * t582 - pkin(3);
t831 = qJD(3) * t573;
t830 = qJD(3) * t579;
t829 = qJD(3) * t582;
t827 = qJD(3) * t589;
t826 = qJD(4) * t588;
t825 = qJD(4) * t589;
t824 = qJD(5) * t538;
t823 = qJD(5) * t584;
t822 = qJD(5) * t588;
t821 = qJD(6) * t410;
t820 = qJD(6) * t584;
t819 = qJD(6) * t588;
t100 = -t342 * t427 - t344 * t426;
t818 = t100 * qJD(2);
t745 = t445 / 0.2e1;
t764 = -t912 / 0.2e1;
t635 = t344 * t745 + t578 * t764;
t122 = -t909 / 0.2e1 + t635;
t817 = t122 * qJD(2);
t762 = -t893 / 0.2e1;
t765 = t342 * t941;
t130 = t765 - t855 / 0.2e1 + (t764 + t762) * t584;
t816 = t130 * qJD(2);
t652 = t344 * t941 + t383 * t538;
t711 = t562 / 0.2e1 - t877 / 0.2e1;
t131 = t652 - t711;
t815 = t131 * qJD(2);
t144 = -t342 * t480 - t410 * t426;
t814 = t144 * qJD(2);
t145 = t344 * t480 + t410 * t427;
t813 = t145 * qJD(2);
t175 = -t410 * t481 - t480 * t963;
t812 = t175 * qJD(2);
t233 = 0.2e1 * t383;
t811 = t233 * qJD(2);
t630 = t939 + t531;
t370 = t630 * t892;
t810 = t370 * qJD(2);
t809 = t370 * qJD(3);
t808 = t371 * qJD(2);
t806 = t372 * qJD(2);
t644 = t720 - t872 / 0.2e1;
t629 = t940 + t644;
t373 = t629 * t892;
t805 = t373 * qJD(2);
t802 = t440 * qJD(3);
t801 = t478 * qJD(6);
t528 = t656 * qJD(5);
t537 = (-t586 ^ 2 + t589 ^ 2) * t575;
t800 = t537 * qJD(2);
t799 = t580 * qJD(5);
t797 = pkin(5) * t957;
t796 = pkin(5) * t956;
t795 = pkin(5) * t950;
t794 = -pkin(5) * t427 / 0.2e1;
t793 = t927 / 0.2e1;
t785 = t929 * t339;
t780 = t344 * t838;
t776 = t578 * t832;
t775 = t538 * t828;
t773 = t580 * t825;
t772 = t589 * t799;
t771 = t584 * t822;
t770 = t656 * t819;
t769 = t538 * t528;
t768 = t575 * t834;
t767 = qJD(3) * t893;
t766 = t584 * t819;
t522 = t538 * t822;
t763 = t896 / 0.2e1;
t761 = t893 / 0.2e1;
t759 = t892 / 0.2e1;
t758 = -t889 / 0.2e1;
t757 = t884 / 0.2e1;
t756 = -t969 / 0.2e1;
t755 = -t881 / 0.2e1;
t754 = -t878 / 0.2e1;
t753 = t878 / 0.2e1;
t752 = -t442 / 0.2e1;
t751 = t870 / 0.2e1;
t749 = -t867 / 0.2e1;
t748 = t862 / 0.2e1;
t747 = t859 / 0.2e1;
t746 = -t856 / 0.2e1;
t744 = t973 - t159 / 0.2e1;
t743 = -t309 / 0.2e1 - t354 / 0.2e1;
t741 = 0.2e1 * t951;
t740 = -t569 / 0.2e1 - t570 / 0.2e1;
t737 = t565 * t459;
t735 = t802 - t820;
t734 = qJD(3) + t835;
t733 = -qJD(6) - t838;
t732 = qJD(4) + t831;
t731 = -qJD(6) + t833;
t727 = t827 * t894;
t726 = t586 * t768;
t719 = t793 + t243 / 0.2e1;
t718 = t738 / 0.2e1;
t717 = pkin(5) * t410 - pkin(11) * t963;
t716 = -pkin(5) * t656 - pkin(11) * t538;
t715 = 0.2e1 * t584 * t522;
t713 = t267 / 0.2e1 + t152 * t939;
t712 = t203 * t939 + t947;
t709 = -qJD(5) + t778;
t213 = t884 + t907;
t606 = t148 * t944 + t213 * t953 - t269 * t957 + t426 * t959;
t695 = -t886 + t899;
t607 = t246 * t752 + t695 * t941;
t14 = t606 + t607;
t658 = t785 - t876;
t201 = -pkin(5) * t893 - t658;
t160 = t585 * t277 + t296 * t929;
t153 = -pkin(11) * t892 + t160;
t505 = -pkin(3) * t583 + t529;
t425 = -pkin(4) * t524 + t505;
t603 = t425 + t717;
t55 = t584 * t153 - t588 * t603;
t93 = t863 - t889;
t7 = t152 * t426 + t201 * t342 + t410 * t93 - t480 * t55;
t708 = t14 * qJD(1) + t7 * qJD(2);
t214 = -t862 + t908;
t605 = t149 * t944 + t214 * t954 - t269 * t955 + t427 * t959;
t696 = t864 + t900;
t608 = t246 * t763 + t696 * t941;
t16 = t605 - t608;
t56 = t588 * t153 + t584 * t603;
t94 = t867 + t885;
t8 = t152 * t427 + t201 * t344 - t410 * t94 - t480 * t56;
t707 = t16 * qJD(1) + t8 * qJD(2);
t667 = t148 * t952 + t204 * t957;
t17 = t918 / 0.2e1 + t667;
t5 = t160 * t342 - t55 * t963 + (t584 * t851 + t865) * t410;
t706 = t17 * qJD(1) + t5 * qJD(2);
t665 = t149 * t952 + t204 * t955;
t20 = -t919 / 0.2e1 + t665;
t6 = t160 * t344 - t56 * t963 + (t588 * t851 - t887) * t410;
t705 = t20 * qJD(1) + t6 * qJD(2);
t664 = t362 * t972 + t363 * t937;
t590 = (t530 / 0.2e1 + t664) * t459 + t340 * t418 / 0.2e1 + t505 * t947 + t604 * t417 / 0.2e1;
t26 = t590 - t965;
t92 = t362 * t417 + t363 * t418 + t505 * t530;
t704 = t26 * qJD(1) + t92 * qJD(2);
t27 = -t152 * t342 + t410 * t55;
t661 = -t900 / 0.2e1 - t864 / 0.2e1;
t668 = t148 * t953 + t203 * t958;
t34 = t661 - t668;
t703 = qJD(1) * t34 - qJD(2) * t27;
t28 = t152 * t344 - t410 * t56;
t660 = t899 / 0.2e1 - t886 / 0.2e1;
t666 = t149 * t954 + t203 * t955;
t33 = t660 - t666;
t702 = qJD(1) * t33 - qJD(2) * t28;
t41 = t159 * t893 - t479 * t410 - t425 * t480 + t658 * t892;
t600 = (t203 * t934 - t269 * t931) * t580 + t480 * t948 + t410 * t947;
t44 = t904 / 0.2e1 + t600;
t701 = t44 * qJD(1) - t41 * qJD(2);
t42 = t479 * t963 + t425 * t481 + (-t160 * t586 + t589 * t850) * t580;
t599 = (t204 * t934 - t270 * t931) * t580 + t481 * t948 + t963 * t947;
t46 = -t903 / 0.2e1 + t599;
t700 = t46 * qJD(1) + t42 * qJD(2);
t597 = t580 * t604;
t659 = t897 / 0.2e1 + t898 / 0.2e1;
t591 = -t659 * t459 + (t597 * t937 - t340 * t579 * t580 / 0.2e1) * t589;
t47 = t591 + t663;
t85 = -t417 * t525 + t418 * t524 + (-t362 * t582 - t363 * t579) * t892;
t699 = t47 * qJD(1) + t85 * qJD(2);
t646 = -t875 / 0.2e1 - t783 / 0.2e1;
t651 = t203 * t760 + t410 * t949;
t57 = -t646 + t651;
t83 = -t159 * t892 - t425 * t410;
t698 = -qJD(1) * t57 - qJD(2) * t83;
t645 = -t874 / 0.2e1 + t784 / 0.2e1;
t650 = t204 * t759 + t948 * t963;
t59 = -t645 + t650;
t84 = -t160 * t892 - t425 * t963;
t697 = -qJD(1) * t59 + qJD(2) * t84;
t693 = -t417 * t579 + t418 * t582;
t140 = -t362 * t893 + t417 * t892 - t505 * t792 + t530 * t524;
t592 = t460 * t943 + t586 * t597 / 0.2e1;
t86 = t901 / 0.2e1 + t592;
t692 = -t86 * qJD(1) + t140 * qJD(2);
t141 = t530 * t525 + (-t363 * t586 + (t505 * t582 + t418) * t589) * t580;
t649 = t340 * t762 + t460 * t942;
t89 = -t902 / 0.2e1 + t649;
t691 = -t89 * qJD(1) - t141 * qJD(2);
t197 = t630 * t459;
t612 = t425 * t938 + t471 * t759 + t573 * t951;
t647 = -t876 / 0.2e1 + t785 / 0.2e1;
t75 = t612 - t647;
t690 = t197 * qJD(1) - t75 * qJD(2);
t198 = t629 * t459;
t613 = t425 * t941 + t470 * t760 + t573 * t954;
t73 = t613 - t743;
t689 = t198 * qJD(1) - t73 * qJD(2);
t365 = t798 * t584;
t615 = (t755 + t958) * t538 - t410 * t440;
t62 = t746 + t615;
t688 = -qJD(2) * t62 + qJD(3) * t365;
t366 = t962 * t584;
t634 = -t656 * t742 + t938 * t963;
t602 = t342 * t938 + t584 * t634;
t64 = t746 + t602;
t687 = qJD(2) * t64 - qJD(3) * t366;
t367 = t798 * t588;
t614 = (t747 + t955) * t538 + t410 * t445;
t66 = t754 + t614;
t686 = -qJD(2) * t66 - qJD(3) * t367;
t453 = t962 * t588;
t601 = t344 * t938 + t588 * t634;
t68 = t753 + t601;
t685 = qJD(2) * t68 - qJD(3) * t453;
t684 = t731 * t588;
t594 = t340 * t943 + t604 * t942;
t123 = t594 + t849;
t176 = -t362 * t525 + t363 * t524;
t683 = -qJD(1) * t123 + qJD(2) * t176;
t616 = t633 * t580;
t609 = -t616 / 0.2e1;
t595 = t460 * t936 + t586 * t609;
t256 = t595 + t848;
t473 = -pkin(2) * t894 - t530 * t583;
t682 = qJD(1) * t256 - qJD(2) * t473;
t257 = (t581 * t751 + t949) * t583 + (t722 + t609) * t589;
t472 = -t529 * t583 + t575 * t928;
t681 = qJD(1) * t257 + qJD(2) * t472;
t126 = -t410 * t656 - t538 * t963;
t188 = -t963 ^ 2 + t960;
t680 = qJD(2) * t188 + qJD(3) * t126;
t679 = qJD(2) * t126 - qJD(3) * t798;
t230 = t741 * t584;
t678 = -qJD(2) * t230 - qJD(3) * t442;
t677 = -qJD(2) * t968 + qJD(3) * t445;
t496 = t504 / 0.2e1;
t401 = t496 + t873 / 0.2e1;
t674 = qJD(2) * t401 + qJD(3) * t531;
t408 = -t897 - t898;
t446 = t524 ^ 2 + t525 ^ 2;
t673 = qJD(2) * t446 - qJD(3) * t408;
t672 = qJD(2) * t408 - qJD(3) * t565;
t671 = -t925 / 0.2e1 + t924 / 0.2e1;
t670 = (-t148 / 0.2e1 + t888 / 0.2e1) * t538;
t669 = (-t149 / 0.2e1 + t866 / 0.2e1) * t538;
t235 = t741 * t588;
t653 = qJD(2) * t235 + t775;
t648 = t946 + t671;
t611 = pkin(11) * t944 + t160 * t939 - t656 * t744;
t641 = t573 + t716;
t288 = -t588 * t641 + t879;
t623 = t288 * t951 + t342 * t945 + t55 * t938;
t639 = t201 / 0.2e1 + t410 * t946 + t243 * t940;
t1 = t584 * t611 - t588 * t639 + t623 + t795;
t21 = -t916 / 0.2e1 + t670;
t71 = -t461 * t445 + (-t288 + t879) * t538;
t643 = t21 * qJD(1) - t1 * qJD(2) + t71 * qJD(3);
t289 = t584 * t641 + t857;
t621 = t289 * t951 + t344 * t945 + t56 * t938;
t2 = t584 * t639 + t588 * t611 + t621 + t794;
t24 = t917 / 0.2e1 + t669;
t72 = t461 * t440 + (-t289 + t857) * t538;
t642 = t24 * qJD(1) - t2 * qJD(2) + t72 * qJD(3);
t180 = t289 * t656 + t470 * t896;
t37 = t757 - t920 / 0.2e1 + t712 * t588;
t622 = t289 * t953 + t470 * t956 + t56 * t940;
t9 = t588 * t713 + t622 + t758;
t640 = qJD(1) * t37 + qJD(2) * t9 - qJD(3) * t180;
t624 = t288 * t954 + t470 * t957 + t55 * t941;
t10 = -t584 * t713 + t624 + t749;
t179 = -t288 * t656 - t442 * t470;
t38 = t748 + t921 / 0.2e1 - t712 * t584;
t638 = -qJD(1) * t38 - qJD(2) * t10 + qJD(3) * t179;
t139 = t342 ^ 2 - t344 ^ 2;
t91 = (-t883 + t913) * t538;
t637 = qJD(2) * t139 - qJD(3) * t91 - qJD(5) * t694;
t598 = qJ(4) * t659 - t664;
t128 = t598 + t740;
t593 = t915 / 0.2e1 - t596 / 0.2e1;
t155 = t723 + (t730 / 0.2e1 + t718) * t586 - t593;
t549 = t565 * qJ(4);
t636 = qJD(1) * t155 - qJD(2) * t128 - qJD(3) * t549;
t171 = t648 * t584;
t29 = t584 * t719 + t588 * t744 + t797;
t632 = pkin(5) * t822 - qJD(2) * t29 - qJD(3) * t171;
t173 = t648 * t588;
t31 = t584 * t744 - t588 * t719 + t796;
t631 = pkin(5) * t823 - qJD(2) * t31 + qJD(3) * t173;
t162 = t657 * t538;
t206 = -t883 / 0.2e1 + t913 / 0.2e1;
t628 = qJD(3) * t162 - qJD(5) * t206 + t342 * t839;
t439 = (t577 / 0.2e1 - t578 / 0.2e1) * t538;
t627 = qJD(2) * t206 - qJD(3) * t439 + t771;
t626 = qJD(6) * t401 + t655;
t625 = qJD(6) * t531 + t654;
t452 = t568 * t533;
t619 = qJD(2) * t91 + qJD(3) * t452 + t715;
t618 = qJD(2) * t694 + t961;
t617 = qJD(2) * t162 + qJD(5) * t439 + t533 * t774;
t610 = (qJD(3) * t710 + t825) * t580;
t560 = qJD(3) * t761;
t523 = t531 * qJD(5);
t506 = (t768 - t799 / 0.2e1) * t586;
t492 = -0.2e1 * t538 * t766;
t467 = t856 / 0.2e1;
t429 = t440 * qJD(6);
t428 = t439 * qJD(6);
t402 = t496 - t504 / 0.2e1;
t390 = t408 * qJD(4);
t361 = t373 * qJD(3);
t280 = qJD(3) * t478 + qJD(5) * t401;
t259 = -t595 + t848;
t258 = t459 * t936 + t616 * t931 + (-t786 / 0.2e1 + t583 * t751) * t581;
t236 = t933 * t963 + t747;
t229 = t935 * t963 + t755;
t228 = t756 + t969 / 0.2e1;
t205 = t206 * qJD(6);
t200 = (t531 + t938) * t459;
t199 = (t644 + t941) * t459;
t174 = t906 + t858 / 0.2e1 - t671 * t588;
t172 = t905 - t880 / 0.2e1 + t671 * t584;
t161 = t162 * qJD(6);
t156 = t586 * t718 + (t583 * t721 + t750) * t581 + t593;
t132 = t652 + t711;
t129 = t410 * t752 + t765 + t855 / 0.2e1 + t584 * t761;
t127 = t598 - t740;
t125 = t126 * qJD(5);
t124 = -t594 + t849;
t121 = t909 / 0.2e1 + t635;
t107 = t694 * qJD(6);
t90 = t91 * qJD(6);
t88 = t902 / 0.2e1 + t649;
t87 = -t901 / 0.2e1 + t592;
t80 = t203 * t588;
t79 = t203 * t584;
t76 = t612 + t647;
t74 = t613 + t743;
t67 = t754 + t601;
t65 = t753 + t614;
t63 = t467 + t602;
t61 = t467 + t615;
t60 = t645 + t650;
t58 = t646 + t651;
t53 = t620 + t662;
t48 = t591 - t663;
t45 = t903 / 0.2e1 + t599;
t43 = -t904 / 0.2e1 + t600;
t40 = t920 / 0.2e1 + t203 * t763 + t757 + t907 / 0.2e1;
t39 = -t921 / 0.2e1 + t203 * t752 + t748 - t908 / 0.2e1;
t36 = t660 + t666;
t35 = t661 + t668;
t32 = pkin(11) * t383 + t796 + t152 * t935 + t890 / 0.2e1 + t865 / 0.2e1;
t30 = t584 * t793 + t797 + t152 * t932 + t868 / 0.2e1 - t887 / 0.2e1;
t25 = t590 + t965;
t23 = -t917 / 0.2e1 + t669;
t22 = t916 / 0.2e1 + t670;
t19 = t919 / 0.2e1 + t665;
t18 = -t918 / 0.2e1 + t667;
t15 = t605 + t608;
t13 = t606 - t607;
t12 = t152 * t763 + t758 + t863 / 0.2e1 - t622;
t11 = t152 * t752 + t749 - t885 / 0.2e1 - t624;
t4 = (t880 - t905) * t954 + (-t868 + t887) * t941 + t470 * t383 + t160 * t763 + t152 * t745 + pkin(11) * t746 + t794 + t201 * t935 - t621;
t3 = (t858 + t906) * t953 + (t865 + t890) * t940 + t470 * t756 + t160 * t442 / 0.2e1 + t440 * t973 + pkin(11) * t754 + t795 + t201 * t933 - t623;
t49 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t82 + qJD(3) * t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t891, -qJD(2) * t788, 0, 0, 0, 0, 0 (-t575 * t581 * t869 - t490 * t583) * qJD(2) + t259 * qJD(3) (-t491 * t583 + t575 * t790) * qJD(2) + t258 * qJD(3) (-t423 * t892 - t490 * t524) * qJD(2) + t87 * qJD(3) (t424 * t892 + t490 * t525) * qJD(2) + t88 * qJD(3) (-t423 * t525 + t424 * t524) * qJD(2) + t48 * qJD(3), t853 + (t362 * t423 + t363 * t424 + t490 * t505) * qJD(2) + t25 * qJD(3) + t124 * qJD(4), 0, 0, 0, 0, 0 (t246 * t892 + t490 * t410) * qJD(2) + t43 * qJD(3) + t60 * qJD(5) (t247 * t892 + t490 * t963) * qJD(2) + t45 * qJD(3) + t58 * qJD(5), 0, 0, 0, 0, 0 (t246 * t342 + t410 * t695) * qJD(2) + t13 * qJD(3) + t18 * qJD(5) + t36 * qJD(6) (t246 * t344 - t410 * t696) * qJD(2) + t15 * qJD(3) + t19 * qJD(5) + t35 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t259 - qJD(3) * t460, qJD(2) * t258 + qJD(3) * t459, qJD(2) * t87 - t460 * t829, qJD(2) * t88 + t460 * t830, t48 * qJD(2) - qJD(3) * t737, t854 + t25 * qJD(2) + (-t460 * pkin(3) - qJ(4) * t737) * qJD(3) + t156 * qJD(4), 0, 0, 0, 0, 0, qJD(2) * t43 + qJD(5) * t200 - t460 * t833, qJD(2) * t45 + qJD(5) * t199 + t460 * t832, 0, 0, 0, 0, 0, t13 * qJD(2) + (-t213 * t656 - t269 * t442) * qJD(3) + t22 * qJD(5) + t40 * qJD(6), t15 * qJD(2) + (t214 * t656 - t269 * t896) * qJD(3) + t23 * qJD(5) + t39 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t124 + qJD(3) * t156, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t60 + qJD(3) * t200 - qJD(5) * t204, qJD(2) * t58 + qJD(3) * t199 + qJD(5) * t203, 0, 0, 0, 0, 0, qJD(2) * t18 + qJD(3) * t22 + qJD(6) * t79 - t204 * t822, qJD(2) * t19 + qJD(3) * t23 + qJD(6) * t80 + t204 * t823; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t36 + qJD(3) * t40 + qJD(5) * t79 - qJD(6) * t149, qJD(2) * t35 + qJD(3) * t39 + qJD(5) * t80 + qJD(6) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256 * qJD(3), -t257 * qJD(3), t86 * qJD(3), t89 * qJD(3), t47 * qJD(3), qJD(3) * t26 - qJD(4) * t123 - t853, 0, 0, 0, 0, 0, qJD(3) * t44 + qJD(5) * t59, qJD(3) * t46 + qJD(5) * t57, 0, 0, 0, 0, 0, qJD(3) * t14 + qJD(5) * t17 - qJD(6) * t33, qJD(3) * t16 + qJD(5) * t20 - qJD(6) * t34; 0, 0, 0, 0, t727, t537 * qJD(3), t580 * t583 * t827, -t583 * t767, 0, t473 * qJD(3), -t472 * qJD(3), -t140 * qJD(3) + t525 * t773, t141 * qJD(3) + t524 * t773, qJD(3) * t85 + qJD(4) * t446, qJD(3) * t92 + qJD(4) * t176 (qJD(3) * t481 - t804) * t963, qJD(3) * t175 + qJD(5) * t188, -t266 * qJD(3) + t410 * t772, -t265 * qJD(3) + t772 * t963, -t727, -t41 * qJD(3) - t84 * qJD(5) + t773 * t963, t42 * qJD(3) + t83 * qJD(5) - t410 * t773 (qJD(3) * t427 - qJD(6) * t342 - t588 * t804) * t344, qJD(3) * t100 + qJD(5) * t99 + qJD(6) * t139, qJD(3) * t145 + qJD(5) * t105 - t342 * t821, qJD(3) * t144 - qJD(5) * t103 - t344 * t821 (qJD(3) * t480 + t803) * t410, qJD(3) * t7 - qJD(4) * t104 + qJD(5) * t5 + qJD(6) * t28, qJD(3) * t8 + qJD(4) * t106 + qJD(5) * t6 + qJD(6) * t27; 0, 0, 0, 0, t726, t800, t734 * t892, -t734 * t893, 0, -qJD(3) * t530 - t682, qJD(3) * t529 - t681, -t530 * t829 + t579 * t610 - t692, t530 * t830 + t582 * t610 - t691, qJD(3) * t693 - t390 + t699 (-t530 * pkin(3) + qJ(4) * t693) * qJD(3) + t127 * qJD(4) + t704, t481 * t975 - t966, t812 + (-t480 * t538 + t481 * t656) * qJD(3) + t125, qJD(5) * t373 + t538 * t767 - t840, -qJD(5) * t370 + t656 * t767 - t841, -t506 (-t470 * t893 - t479 * t656 + t480 * t573) * qJD(3) - t371 * qJD(4) + t76 * qJD(5) + t701 (-t471 * t893 + t479 * t538 + t481 * t573) * qJD(3) + t372 * qJD(4) + t74 * qJD(5) + t700, t121 * qJD(5) - t161 + (t775 + t839) * t427, t818 + t53 * qJD(5) - t90 + (-t426 * t588 - t909) * t832, t813 + (-t427 * t656 + t538 * t856) * qJD(3) + t65 * qJD(5) + t129 * qJD(6), t814 + (t426 * t656 - t442 * t480) * qJD(3) + t61 * qJD(5) + t132 * qJD(6), t480 * t967 + t801 + t966 (t201 * t442 - t288 * t480 + t426 * t470 - t656 * t93) * qJD(3) + t63 * qJD(4) + t3 * qJD(5) + t12 * qJD(6) + t708 (t201 * t896 - t289 * t480 + t427 * t470 + t656 * t94) * qJD(3) + t67 * qJD(4) + t4 * qJD(5) + t11 * qJD(6) + t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(2) * t525 + t830) * t892 (qJD(2) * t524 + t829) * t892, t673, qJD(3) * t127 + t683, 0, 0, 0, 0, 0, t402 * qJD(5) - t976, t977, 0, 0, 0, 0, 0, qJD(3) * t63 + qJD(5) * t236 + qJD(6) * t228 - t846, qJD(3) * t67 + qJD(5) * t229 + t844 + t978; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t655, t680, t410 * t709 + t361, t709 * t963 - t809, t560, qJD(3) * t76 + qJD(4) * t402 - qJD(5) * t160 - t697, qJD(3) * t74 + qJD(5) * t159 - t698, t121 * qJD(3) + t205 + (-t823 - t839) * t968, t53 * qJD(3) - t410 * t736 - t107 + t852, qJD(3) * t65 + t584 * t803 + t845 + t978, qJD(3) * t61 + t588 * t803 - t847, t626, t3 * qJD(3) + t236 * qJD(4) + (-t160 * t588 + t584 * t717) * qJD(5) + t32 * qJD(6) + t706, t4 * qJD(3) + t229 * qJD(4) + (t160 * t584 + t588 * t717) * qJD(5) + t30 * qJD(6) + t705; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t628, t637, t129 * qJD(3) + qJD(5) * t974 + t342 * t733, t132 * qJD(3) + t344 * t733, t280, qJD(3) * t12 + qJD(4) * t228 + qJD(5) * t32 - qJD(6) * t56 - t702, qJD(3) * t11 + qJD(4) * t974 + qJD(5) * t30 + qJD(6) * t55 - t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, t256 * qJD(2), t257 * qJD(2), -t86 * qJD(2), -t89 * qJD(2), -t47 * qJD(2), -qJD(2) * t26 - qJD(4) * t155 - t854, 0, 0, 0, 0, 0, -qJD(2) * t44 - qJD(5) * t197, -qJD(2) * t46 - qJD(5) * t198, 0, 0, 0, 0, 0, -qJD(2) * t14 + qJD(5) * t21 - qJD(6) * t37, -qJD(2) * t16 + qJD(5) * t24 - qJD(6) * t38; 0, 0, 0, 0, -t726, -t800, -t583 * t778, t835 * t893, 0, t682, t681, t692, t691, -t390 - t699, qJD(4) * t128 - t704, -t481 * t837 - t966, t125 - t812, -qJD(5) * t372 + t840, -qJD(5) * t371 + t841, t506, -qJD(4) * t370 + qJD(5) * t75 - t701, -qJD(4) * t373 + qJD(5) * t73 - t700, qJD(5) * t122 - t427 * t839 - t161, qJD(5) * t54 - t818 - t90, qJD(5) * t66 + qJD(6) * t130 - t813, qJD(5) * t62 + qJD(6) * t131 - t814, -t480 * t838 - t801 + t966, qJD(4) * t64 - qJD(5) * t1 - qJD(6) * t9 - t708, qJD(4) * t68 - qJD(5) * t2 - qJD(6) * t10 - t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t565 * qJD(4), t549 * qJD(4), t769, -t798 * qJD(5), 0, 0, 0, t573 * t824, t573 * t528, -t533 * t766 + t578 * t769, -qJD(6) * t452 - t656 * t715, t538 * t656 * t820 + qJD(5) * t367, -qJD(5) * t365 + t538 * t770, -t769, -qJD(4) * t366 + qJD(5) * t71 + qJD(6) * t180, -qJD(4) * t453 + qJD(5) * t72 + qJD(6) * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t672, -t636, 0, 0, 0, 0, 0, -t810, -t805, 0, 0, 0, 0, 0, t687, t685; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t654, t679, t528 - t806, -t808 - t824, qJD(2) * t762, -qJD(5) * t471 + t538 * t831 - t690, qJD(5) * t470 + t656 * t831 - t689, t817 - t428 - (-t771 - t776) * t656, -t656 * t961 + t492 + t895, t538 * t823 - t686, t522 - t688, t625 (t584 * t716 - t857) * qJD(5) + t174 * qJD(6) + t643 (t588 * t716 + t879) * qJD(5) + t172 * qJD(6) + t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t617, -t619, t442 * t731 + t816, t538 * t684 + t815, t523 - t836, qJD(5) * t174 - qJD(6) * t289 - t640, qJD(5) * t172 + qJD(6) * t288 + t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t123 + qJD(3) * t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t525 * t778, -t524 * t778, -t673, -qJD(3) * t128 - t683, 0, 0, 0, 0, 0, -t729 + t809 + t803, t361 + t728 - t804, 0, 0, 0, 0, 0, -qJD(3) * t64 + qJD(5) * t235 - qJD(6) * t227 + t846, -qJD(3) * t68 - qJD(5) * t230 + qJD(6) * t233 - t844; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t672, t636, 0, 0, 0, 0, 0, t810 + t824, t528 + t805, 0, 0, 0, 0, 0, t429 + t522 - t687, -qJD(5) * t442 - t685 + t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t975, -t967, 0, 0, 0, 0, 0, t653, t678; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t735 - t842, t684 + t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t59 + qJD(3) * t197, -qJD(2) * t57 + qJD(3) * t198, 0, 0, 0, 0, 0, -qJD(2) * t17 - qJD(3) * t21, -qJD(2) * t20 - qJD(3) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t655, -t680, t977, t976, t560, -qJD(3) * t75 - qJD(4) * t963 + t697, -qJD(3) * t73 + qJD(4) * t410 + t698, -qJD(3) * t122 + t588 * t780 + t205, -qJD(3) * t54 - t107 - t852, -qJD(3) * t66 + qJD(6) * t968 - t845, -qJD(3) * t62 - qJD(6) * t969 + t847, -t626, qJD(3) * t1 - qJD(4) * t235 + qJD(6) * t31 - t706, qJD(3) * t2 + qJD(4) * t230 + qJD(6) * t29 - t705; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t654, -t679, t806, t808, qJD(2) * t761, -t538 * t732 + t690, -t656 * t732 + t689, -t656 * t776 - t428 - t817, t656 * t714 + t492 - t895, -qJD(6) * t445 + t686, t429 + t688, -t625, -qJD(6) * t173 - t538 * t826 - t643, qJD(4) * t442 + qJD(6) * t171 - t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t975, t967, 0, 0, 0, 0, 0, -t653, -t678; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t766, t568 * qJD(6), 0, 0, 0, -pkin(5) * t820, -pkin(5) * t819; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t627, -t618, -t677 + t819, t735 - t843, -t674, -pkin(11) * t819 - t631, pkin(11) * t820 - t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t33 + qJD(3) * t37, qJD(2) * t34 + qJD(3) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t628, -t637, -qJD(3) * t130 - qJD(5) * t968 + t342 * t838, -qJD(3) * t131 + qJD(5) * t969 + t780, t280, qJD(3) * t9 + qJD(4) * t227 - qJD(5) * t31 + t702, qJD(3) * t10 - qJD(4) * t233 - qJD(5) * t29 + t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t617, t619, qJD(5) * t445 - t584 * t777 - t816, -qJD(5) * t440 - t656 * t775 - t815, t523 + t836, -qJD(4) * t440 + qJD(5) * t173 + t640, -qJD(5) * t171 - t656 * t826 - t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t802 + t842, -t656 * t828 - t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t627, t618, t677, -t802 + t843, t674, t631, t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t49;