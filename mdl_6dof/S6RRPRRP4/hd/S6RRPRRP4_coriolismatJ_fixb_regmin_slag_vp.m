% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:43
% EndTime: 2019-03-09 11:56:14
% DurationCPUTime: 17.49s
% Computational Cost: add. (22579->759), mult. (43599->948), div. (0->0), fcn. (49839->8), ass. (0->556)
t538 = sin(qJ(2));
t842 = cos(pkin(10));
t634 = t842 * t538;
t535 = sin(pkin(10));
t541 = cos(qJ(2));
t763 = t535 * t541;
t505 = t634 + t763;
t537 = sin(qJ(4));
t393 = t537 * t505;
t850 = -qJ(3) - pkin(7);
t515 = t850 * t538;
t516 = t850 * t541;
t425 = -t842 * t515 - t516 * t535;
t338 = pkin(4) * t393 + t425;
t536 = sin(qJ(5));
t539 = cos(qJ(5));
t540 = cos(qJ(4));
t752 = t539 * t540;
t369 = -t393 * t536 + t505 * t752;
t799 = t369 * qJ(6);
t753 = t539 * t537;
t759 = t536 * t540;
t512 = t753 + t759;
t283 = t512 * t505;
t857 = t283 * pkin(5);
t616 = -t799 + t857;
t159 = t338 + t616;
t510 = t536 * t537 - t752;
t870 = t510 / 0.2e1;
t525 = -pkin(2) * t842 - pkin(3);
t514 = -t540 * pkin(4) + t525;
t774 = t512 * qJ(6);
t855 = t510 * pkin(5);
t615 = -t774 + t855;
t371 = t514 + t615;
t942 = t371 / 0.2e1;
t503 = t535 * t538 - t541 * t842;
t523 = pkin(2) * t535 + pkin(8);
t851 = pkin(9) + t523;
t500 = t851 * t540;
t635 = t851 * t537;
t630 = t536 * t500 + t539 * t635;
t935 = t630 * t503;
t948 = -t935 / 0.2e1;
t950 = -t159 * t870 - t283 * t942 - t948;
t780 = t505 * t510;
t802 = t283 * t505;
t773 = t512 * t503;
t805 = t773 * t503;
t610 = t802 + t805;
t911 = qJD(1) * t610;
t580 = t759 / 0.2e1 + t753 / 0.2e1;
t868 = t512 / 0.2e1;
t560 = (t868 + t580) * t503;
t670 = qJD(4) + qJD(5);
t936 = t560 * t670;
t949 = -qJD(2) * t780 - t911 - t936;
t853 = t538 * pkin(2);
t396 = pkin(3) * t505 + pkin(8) * t503 + t853;
t377 = t540 * t396;
t395 = t540 * t503;
t756 = t425 * t537;
t190 = pkin(4) * t505 + pkin(9) * t395 + t377 + t756;
t189 = t536 * t190;
t390 = t537 * t503;
t376 = t537 * t396;
t937 = t425 * t540;
t727 = -t376 + t937;
t218 = pkin(9) * t390 - t727;
t215 = t539 * t218;
t745 = t189 / 0.2e1 + t215 / 0.2e1;
t871 = -t510 / 0.2e1;
t941 = -t514 / 0.2e1;
t896 = t283 * t941 + t338 * t871;
t58 = t948 - t745 - t896;
t947 = t937 / 0.2e1;
t940 = -t630 / 0.2e1;
t946 = t283 * t940;
t662 = -pkin(2) * t541 - pkin(1);
t388 = t503 * pkin(3) - t505 * pkin(8) + t662;
t914 = t535 * t515 - t842 * t516;
t786 = t914 * t537;
t261 = -t540 * t388 + t786;
t779 = t505 * t540;
t216 = -pkin(9) * t779 - t261;
t561 = t503 * pkin(4) + t216;
t188 = t539 * t561;
t785 = t914 * t540;
t262 = t388 * t537 + t785;
t217 = -pkin(9) * t393 + t262;
t761 = t536 * t217;
t125 = -t188 + t761;
t105 = -pkin(5) * t503 + t125;
t754 = t539 * t216;
t137 = t754 - t761;
t945 = t105 + t137;
t492 = t505 * qJ(6);
t944 = -t492 - t745;
t854 = t536 * pkin(4);
t214 = t539 * t217;
t762 = t536 * t216;
t917 = -t762 / 0.2e1 - t214 / 0.2e1;
t546 = -t503 * t854 + t917;
t886 = t214 / 0.2e1;
t618 = t886 + t762 / 0.2e1;
t543 = t546 + t618;
t847 = pkin(4) * qJD(4);
t934 = -qJD(1) * t543 + t536 * t847;
t647 = -t773 / 0.2e1;
t552 = t503 * t580 + t647;
t732 = t552 * t670;
t943 = t911 + t732;
t569 = t539 * t500 - t536 * t635;
t918 = t569 * t503;
t938 = -t918 / 0.2e1;
t559 = t159 * t868 + t369 * t942 + t938;
t755 = t539 * t190;
t760 = t536 * t218;
t744 = t755 / 0.2e1 - t760 / 0.2e1;
t869 = -t512 / 0.2e1;
t897 = t338 * t869 + t369 * t941;
t60 = t938 + t744 - t897;
t59 = t935 / 0.2e1 - t745 + t896;
t920 = -t569 / 0.2e1;
t939 = t630 / 0.2e1;
t705 = qJD(1) * t505;
t410 = t503 * t705;
t722 = -t763 / 0.2e1 - t634 / 0.2e1;
t916 = -qJD(4) * t722 + t410;
t223 = pkin(5) * t369 + qJ(6) * t283;
t668 = pkin(4) * t779;
t191 = t223 + t668;
t778 = t510 * qJ(6);
t859 = pkin(5) * t512;
t407 = t778 + t859;
t860 = pkin(4) * t537;
t399 = t407 + t860;
t876 = -t399 / 0.2e1;
t933 = -t191 * t869 - t369 * t876 + t950;
t846 = pkin(4) * qJD(5);
t665 = t536 * t846;
t932 = -t665 - t934;
t501 = t503 ^ 2;
t891 = t369 ^ 2;
t288 = t501 + t891;
t627 = t670 * t503;
t701 = qJD(2) * t512;
t931 = -t288 * qJD(1) - t369 * t701 - t627;
t136 = t214 + t762;
t545 = t546 + t917;
t903 = t552 * qJD(3);
t930 = -t136 * qJD(4) + qJD(5) * t545 + t903;
t554 = t536 * t561;
t126 = t554 + t214;
t929 = qJD(4) * t545 - t126 * qJD(5) + t903;
t901 = t560 * qJD(3);
t928 = qJD(5) * t543 + t901;
t927 = -qJD(4) * t543 + t901;
t498 = t512 * qJD(5);
t902 = t560 * qJD(1);
t926 = -t512 * qJD(4) - t498 - t902;
t798 = t369 * t505;
t447 = t536 * t390;
t368 = -t395 * t539 + t447;
t800 = t368 * t503;
t174 = -t798 + t800;
t924 = qJD(1) * t174;
t923 = qJD(3) * t174;
t57 = t918 / 0.2e1 + t744 + t897;
t922 = t670 * t630;
t921 = pkin(4) / 0.2e1;
t919 = t569 / 0.2e1;
t626 = t670 * t510;
t502 = t505 ^ 2;
t915 = -t502 - t501;
t913 = -qJD(5) * t722 + t916;
t419 = t283 / 0.2e1;
t621 = -t283 / 0.2e1 + t419;
t910 = qJD(2) * t621;
t908 = qJD(3) * t610;
t907 = qJD(3) * t780;
t900 = t569 * qJD(6);
t899 = t621 * qJD(3);
t898 = t780 * qJD(1);
t895 = t369 * t920 + t946;
t533 = t537 ^ 2;
t534 = t540 ^ 2;
t520 = t534 - t533;
t633 = 0.2e1 * t540 * t393;
t566 = qJD(1) * t633 - qJD(2) * t520;
t892 = t670 * t569;
t507 = t512 ^ 2;
t890 = -pkin(5) / 0.2e1;
t889 = -qJ(6) / 0.2e1;
t888 = qJ(6) / 0.2e1;
t782 = t503 * qJ(6);
t104 = t126 + t782;
t887 = t104 / 0.2e1;
t885 = -t371 / 0.2e1;
t875 = -t407 / 0.2e1;
t874 = t425 / 0.2e1;
t873 = t447 / 0.2e1;
t872 = -t505 / 0.2e1;
t524 = qJ(6) + t854;
t867 = -t524 / 0.2e1;
t866 = t524 / 0.2e1;
t852 = t539 * pkin(4);
t527 = -pkin(5) - t852;
t865 = -t527 / 0.2e1;
t864 = t527 / 0.2e1;
t863 = t536 / 0.2e1;
t862 = t539 / 0.2e1;
t861 = -t540 / 0.2e1;
t858 = t773 * pkin(5);
t856 = t505 * pkin(5);
t849 = t104 * t871 + t105 * t868;
t845 = qJD(2) * pkin(2);
t593 = t368 * t890 - t773 * t889;
t641 = t125 / 0.2e1 - t105 / 0.2e1;
t643 = t887 - t126 / 0.2e1;
t5 = -t510 * t641 + t512 * t643 + t593;
t844 = t5 * qJD(1);
t557 = t368 * t920 + t371 * t872 - t773 * t940;
t740 = t215 + t189;
t110 = t492 + t740;
t611 = t755 - t760;
t111 = -t611 - t856;
t587 = t110 * t868 + t111 * t870;
t29 = t557 + t587;
t841 = qJD(1) * t29;
t819 = t159 * t369;
t824 = t136 * t503;
t36 = t191 * t283 + t819 - t824;
t840 = qJD(1) * t36;
t820 = t159 * t283;
t823 = t137 * t503;
t37 = -t191 * t369 + t820 + t823;
t839 = qJD(1) * t37;
t827 = t125 * t503;
t38 = -t223 * t369 + t820 - t827;
t838 = qJD(1) * t38;
t826 = t126 * t503;
t39 = t223 * t283 + t819 - t826;
t837 = qJD(1) * t39;
t56 = t104 * t503 - t819;
t834 = qJD(1) * t56;
t811 = t338 * t369;
t61 = t283 * t668 + t811 - t824;
t833 = qJD(1) * t61;
t812 = t338 * t283;
t62 = t369 * t668 - t812 - t823;
t832 = qJD(1) * t62;
t65 = -t812 + t827;
t831 = qJD(1) * t65;
t66 = t811 - t826;
t830 = qJD(1) * t66;
t337 = -pkin(4) * t390 + t914;
t801 = t368 * qJ(6);
t158 = t337 - t801 - t858;
t11 = t104 * t110 + t105 * t111 + t158 * t159;
t829 = t11 * qJD(1);
t12 = -t104 * t125 + t105 * t126 + t159 * t223;
t828 = t12 * qJD(1);
t13 = t104 * t137 + t105 * t136 + t159 * t191;
t825 = t13 * qJD(1);
t14 = t104 * t773 + t105 * t368 - t110 * t283 + t111 * t369;
t822 = t14 * qJD(1);
t15 = (-t104 + t126) * t369 + (-t105 + t125) * t283;
t821 = t15 * qJD(1);
t590 = t858 / 0.2e1 + t801 / 0.2e1;
t16 = t510 * t643 + t512 * t641 + t590;
t817 = t16 * qJD(1);
t18 = (-t104 + t136) * t369 - t945 * t283;
t816 = t18 * qJD(1);
t21 = t104 * t505 + t110 * t503 - t158 * t369 - t159 * t368;
t815 = t21 * qJD(1);
t22 = -t105 * t505 - t111 * t503 + t158 * t283 - t159 * t773;
t814 = t22 * qJD(1);
t33 = t104 * t368 - t105 * t773 + t159 * t505;
t813 = t33 * qJD(1);
t34 = -t125 * t505 + t283 * t337 - t338 * t773 + t503 * t611;
t808 = t34 * qJD(1);
t35 = -t126 * t505 + t337 * t369 + t338 * t368 - t503 * t740;
t807 = t35 * qJD(1);
t806 = t773 * t369;
t804 = t773 * t512;
t803 = t283 * t368;
t796 = t371 * t510;
t795 = t371 * t512;
t791 = t630 * t505;
t787 = t569 * t505;
t542 = t782 / 0.2e1 + t503 * t866 + t554 / 0.2e1 + t886;
t46 = -t542 + t618;
t783 = t46 * qJD(1);
t781 = t503 * t539;
t777 = t510 * t368;
t776 = t510 * t503;
t775 = t510 * t514;
t770 = t524 * t369;
t769 = t524 * t510;
t768 = t524 * t512;
t767 = t527 * t283;
t766 = t527 * t510;
t765 = t527 * t512;
t758 = t537 * t283;
t757 = t537 * t369;
t71 = (-t261 + t786) * t505 + t377 * t503;
t751 = t71 * qJD(1);
t72 = (-t262 + t785) * t505 + (t727 - t937) * t503;
t750 = t72 * qJD(1);
t636 = t919 + t920;
t637 = t940 + t939;
t121 = -t510 * t636 + t512 * t637;
t749 = t121 * qJD(4);
t142 = t283 * t510 - t369 * t512;
t748 = t670 * t142;
t177 = -t283 * t868 + t369 * t871;
t746 = t670 * t177;
t646 = -t761 / 0.2e1;
t742 = t646 + t754 / 0.2e1;
t645 = t761 / 0.2e1;
t741 = t645 - t754 / 0.2e1;
t738 = t510 * qJD(6);
t267 = t873 + (-t752 / 0.2e1 + t870) * t503;
t491 = t503 * qJD(6);
t737 = t267 * qJD(3) + t491;
t644 = -t395 / 0.2e1;
t648 = -t776 / 0.2e1;
t269 = t539 * t644 + t648 + t873;
t736 = t670 * t269;
t735 = t269 * qJD(3) + t491;
t734 = t670 * t267;
t649 = t781 / 0.2e1;
t725 = -t447 / 0.2e1 + t540 * t649;
t650 = -t781 / 0.2e1;
t724 = t536 * t644 + t537 * t650;
t721 = qJ(6) * qJD(5);
t101 = -t777 - t804;
t720 = qJD(1) * t101;
t143 = -t803 - t806;
t719 = qJD(1) * t143;
t156 = t261 * t503 - t393 * t425;
t718 = qJD(1) * t156;
t157 = -t262 * t503 + t425 * t779;
t717 = qJD(1) * t157;
t171 = -t802 + t805;
t716 = qJD(1) * t171;
t173 = t798 + t800;
t714 = qJD(1) * t173;
t224 = t425 * t505 - t503 * t914;
t712 = qJD(1) * t224;
t669 = t502 - t501;
t327 = t669 * t537;
t711 = qJD(1) * t327;
t328 = t915 * t537;
t710 = qJD(1) * t328;
t329 = t669 * t540;
t709 = qJD(1) * t329;
t708 = qJD(1) * t369;
t398 = t915 * t540;
t707 = qJD(1) * t398;
t706 = qJD(1) * t503;
t704 = qJD(1) * t540;
t703 = qJD(1) * t541;
t702 = qJD(2) * t505;
t700 = qJD(2) * t514;
t699 = qJD(2) * t537;
t698 = qJD(2) * t538;
t697 = qJD(2) * t540;
t696 = qJD(2) * t541;
t695 = qJD(3) * t540;
t694 = qJD(4) * t137;
t693 = qJD(4) * t537;
t692 = qJD(4) * t540;
t691 = qJD(5) * t125;
t144 = -t803 + t806;
t689 = t144 * qJD(1);
t165 = t662 * t853;
t688 = t165 * qJD(1);
t264 = t647 + t724;
t685 = t264 * qJD(1);
t246 = t267 * qJD(1);
t268 = t648 + t725;
t248 = t268 * qJD(1);
t682 = t283 * qJD(1);
t285 = 0.2e1 * t419;
t681 = t285 * qJD(1);
t565 = -t535 * t503 / 0.2e1 + t842 * t872;
t343 = (-t538 / 0.2e1 + t565) * pkin(2);
t680 = t343 * qJD(1);
t679 = t283 * qJD(6);
t678 = t390 * qJD(1);
t677 = t393 * qJD(1);
t676 = t395 * qJD(1);
t675 = t915 * qJD(1);
t674 = t722 * qJD(1);
t497 = t512 * qJD(6);
t521 = -t538 ^ 2 + t541 ^ 2;
t672 = t521 * qJD(1);
t531 = t539 * t846;
t671 = t531 + qJD(6);
t667 = pkin(1) * t538 * qJD(1);
t666 = pkin(1) * t703;
t485 = t856 / 0.2e1;
t664 = t854 / 0.2e1;
t663 = t852 / 0.2e1;
t661 = t283 * t708;
t660 = t510 * t701;
t659 = t537 * t697;
t658 = t503 * t692;
t409 = t503 * t702;
t657 = t537 * t692;
t656 = t538 * t703;
t655 = t505 * t704;
t642 = t887 - t136 / 0.2e1;
t640 = t137 / 0.2e1 + t105 / 0.2e1;
t638 = -t376 / 0.2e1 + t947;
t595 = pkin(4) * t649 + t188 / 0.2e1 + t646;
t53 = t595 + t741;
t530 = t539 * t847;
t632 = -qJD(1) * t53 - t530;
t594 = pkin(4) * t650 - t188 / 0.2e1 + t645;
t54 = t594 + t742;
t631 = -qJD(1) * t54 + t530;
t629 = t104 * t869 + t105 * t871 + t895;
t625 = t668 / 0.2e1;
t624 = -qJD(4) - t706;
t620 = qJD(2) * t633;
t617 = -t856 / 0.2e1 - t744;
t585 = t159 * t875 + t223 * t885;
t592 = t110 * t888 + t111 * t890;
t3 = t569 * t641 + t630 * t643 + t585 + t592;
t98 = t371 * t407;
t614 = -t3 * qJD(1) + t98 * qJD(2);
t453 = qJD(5) - t624;
t581 = t368 * t864 - t773 * t867;
t7 = t283 * t637 + t369 * t636 + t510 * t640 + t512 * t642 + t581;
t612 = t7 * qJD(1);
t609 = -t503 * t525 - t505 * t523;
t583 = t368 * t866 - t773 * t864;
t19 = t510 * t642 - t512 * t640 + t583;
t608 = -t19 * qJD(1) + t121 * qJD(2);
t588 = t867 + t664 + t888;
t596 = t890 - t852 / 0.2e1 + t865;
t185 = t510 * t596 + t512 * t588;
t70 = t283 * t596 + t369 * t588;
t607 = qJD(1) * t70 + qJD(2) * t185;
t196 = t399 * t510 + t795;
t548 = t191 * t870 - t283 * t876 + t559;
t26 = (t864 + t890) * t505 + t548 - t744;
t606 = -qJD(1) * t26 - qJD(2) * t196;
t197 = -t399 * t512 + t796;
t570 = t492 / 0.2e1 + t505 * t866 + t745;
t23 = t570 + t933;
t605 = qJD(1) * t23 - qJD(2) * t197;
t200 = t407 * t510 + t795;
t547 = t223 * t870 - t283 * t875 + t559;
t27 = t547 - t744 - t856;
t604 = -qJD(1) * t27 - qJD(2) * t200;
t201 = -t407 * t512 + t796;
t553 = t223 * t869 + t369 * t875 - t950;
t32 = t553 + t944;
t603 = -qJD(1) * t32 - qJD(2) * t201;
t353 = t510 * t860 + t512 * t514;
t43 = (-t758 / 0.2e1 + (t510 * t861 + t862) * t505) * pkin(4) + t57;
t602 = qJD(1) * t43 - qJD(2) * t353;
t354 = t512 * t860 - t775;
t42 = (-t757 / 0.2e1 + (-t536 / 0.2e1 + t512 * t861) * t505) * pkin(4) + t58;
t601 = qJD(1) * t42 - qJD(2) * t354;
t600 = t624 * t540;
t166 = t283 ^ 2 - t891;
t67 = qJD(1) * t166 + qJD(2) * t142;
t355 = t510 ^ 2 - t507;
t116 = qJD(1) * t142 + qJD(2) * t355;
t591 = t136 * pkin(5) / 0.2e1 + t137 * t889;
t589 = pkin(5) * t919 + qJ(6) * t939;
t586 = t110 * t866 + t111 * t864;
t582 = t523 * t503 / 0.2e1 + t525 * t872;
t40 = t559 + t617;
t579 = qJD(1) * t40 + t371 * t701;
t578 = qJD(1) * t58 + t510 * t700;
t577 = qJD(1) * t57 - t512 * t700;
t576 = t505 * t600;
t551 = t582 * t537 + t947;
t146 = t551 - t638;
t575 = -qJD(1) * t146 - t525 * t697;
t563 = t582 * t540;
t148 = -t377 / 0.2e1 - t563 + (t874 - t425 / 0.2e1) * t537;
t574 = -qJD(1) * t148 - t525 * t699;
t163 = -t177 + t722;
t573 = qJD(1) * t163 + t660;
t131 = qJD(2) * t177 - t661;
t155 = -qJD(1) * t177 + t660;
t389 = (t533 / 0.2e1 - t534 / 0.2e1) * t505;
t571 = -qJD(1) * t389 + t659;
t568 = t502 * t537 * t704 + qJD(2) * t389;
t397 = t520 * t502;
t567 = qJD(1) * t397 + t620;
t544 = t104 * t939 + t136 * t940 + t159 * t876 + t191 * t885 + t920 * t945;
t1 = t544 + t586;
t89 = t371 * t399;
t564 = -t1 * qJD(1) + t89 * qJD(2) + t121 * qJD(3);
t550 = (t104 * t862 + t105 * t863) * pkin(4) + t125 * t867 + t126 * t864;
t10 = t550 + t591;
t186 = -t510 * t588 + t512 * t596;
t436 = (t524 * t539 + t527 * t536) * pkin(4);
t549 = (t569 * t862 + t630 * t863) * pkin(4) + t630 * t867 + t569 * t864;
t78 = t549 + t589;
t555 = t10 * qJD(1) + t78 * qJD(2) - t186 * qJD(3) + t436 * qJD(4);
t532 = qJ(6) * qJD(6);
t517 = t524 * qJD(6);
t479 = t722 * qJD(2);
t473 = t505 * t697;
t432 = t453 * qJ(6);
t381 = t390 * qJD(4);
t380 = t389 * qJD(4);
t342 = t853 / 0.2e1 + t565 * pkin(2);
t341 = -t678 - t693;
t305 = t369 * t497;
t271 = t773 / 0.2e1 + t724;
t270 = t776 / 0.2e1 + t725;
t263 = qJD(2) * t507 + t512 * t708;
t254 = t270 * qJD(3);
t247 = t268 * qJD(3);
t193 = -t248 + t626;
t192 = -t626 - t246;
t187 = -t769 / 0.2e1 + t512 * t663 + t765 / 0.2e1 + t510 * t664 - t778 / 0.2e1 - t859 / 0.2e1;
t184 = -t768 / 0.2e1 - t766 / 0.2e1 - t774 / 0.2e1 + t855 / 0.2e1 + (t512 * t863 + t539 * t871) * pkin(4);
t164 = t177 + t722;
t161 = qJD(2) * t560 + t369 * t706;
t160 = qJD(2) * t267 + t283 * t706;
t149 = t537 * t874 + t756 / 0.2e1 + t377 / 0.2e1 - t563;
t147 = t551 + t638;
t113 = qJD(2) * t552 - t369 * t453;
t112 = t269 * qJD(2) - t283 * t453;
t77 = t549 - t589;
t69 = -t770 / 0.2e1 - t767 / 0.2e1 - t799 / 0.2e1 + t857 / 0.2e1 + (-t283 * t862 + t369 * t863) * pkin(4);
t55 = t594 + t741;
t52 = t595 + t742;
t47 = t542 + t618;
t45 = t512 * t625 + t757 * t921 + t854 * t872 + t59;
t44 = t505 * t663 + t510 * t625 + t758 * t921 + t60;
t41 = -t559 + t617;
t31 = t553 - t944;
t30 = -t557 + t587;
t28 = 0.2e1 * t485 + t547 + t744;
t25 = t505 * t865 + t485 + t548 + t744;
t24 = t570 - t933;
t20 = t136 * t870 + t137 * t868 + t583 + t849;
t17 = t125 * t869 + t126 * t870 + t590 + t849;
t9 = t550 - t591;
t8 = t136 * t868 + t137 * t871 + t369 * t919 + t581 + t629 - t946;
t6 = t125 * t870 + t126 * t868 + t593 + t629 - t895;
t4 = t104 * t940 + t105 * t919 + t125 * t920 + t126 * t939 - t585 + t592;
t2 = -t544 + t586;
t48 = [0, 0, 0, t538 * t696, t521 * qJD(2), 0, 0, 0, -pkin(1) * t698, -pkin(1) * t696, -qJD(3) * t915, qJD(2) * t165 + qJD(3) * t224, -t409 * t534 - t502 * t657, -qJD(4) * t397 + t503 * t620, -t503 * t505 * t693 + qJD(2) * t329, -qJD(2) * t327 - t505 * t658, t409, qJD(2) * t71 - qJD(3) * t328 + qJD(4) * t157, qJD(2) * t72 - qJD(3) * t398 + qJD(4) * t156 (qJD(2) * t368 - t283 * t670) * t369, t144 * qJD(2) + t166 * t670, t173 * qJD(2) - t283 * t627, t171 * qJD(2) - t369 * t627, t409, qJD(2) * t34 + qJD(4) * t61 + qJD(5) * t66 + t908, qJD(2) * t35 + qJD(4) * t62 + qJD(5) * t65 - t923, qJD(2) * t22 + qJD(4) * t36 + qJD(5) * t39 - t369 * t679 + t908, qJD(2) * t14 + qJD(3) * t143 + qJD(4) * t18 + qJD(5) * t15 - t283 * t491, qJD(2) * t21 + qJD(4) * t37 + qJD(5) * t38 + qJD(6) * t288 + t923, qJD(2) * t11 + qJD(3) * t33 + qJD(4) * t13 + qJD(5) * t12 + qJD(6) * t56; 0, 0, 0, t656, t672, t696, -t698, 0, -pkin(7) * t696 - t667, pkin(7) * t698 - t666 (t503 * t842 - t505 * t535) * t845, t688 + (-t425 * t535 - t842 * t914) * t845 + t342 * qJD(3), -t380 + (-t534 * t705 - t659) * t503, t503 * t566 - 0.2e1 * t505 * t657, t505 * t699 + t709, t473 - t711, t916, t751 + (t537 * t609 - t785) * qJD(2) + t149 * qJD(4), t750 + (t540 * t609 + t786) * qJD(2) + t147 * qJD(4) (t701 + t708) * t368 + t746, t689 + (-t777 + t804) * qJD(2) + t748, t505 * t701 + t714 + t736, -t510 * t702 + t716 + t732, t913, t808 + (t337 * t510 - t514 * t773 - t791) * qJD(2) + t44 * qJD(4) + t60 * qJD(5), t807 + (t337 * t512 + t368 * t514 - t787) * qJD(2) + t899 + t45 * qJD(4) + t59 * qJD(5), t814 + (t158 * t510 - t371 * t773 - t791) * qJD(2) + t25 * qJD(4) + t28 * qJD(5) + t164 * qJD(6), t822 + (-t110 * t510 + t111 * t512 + t368 * t630 + t569 * t773) * qJD(2) + t8 * qJD(4) + t6 * qJD(5) + t269 * qJD(6), t815 + (-t158 * t512 - t368 * t371 + t787) * qJD(2) + t899 + t24 * qJD(4) + t31 * qJD(5) + t305, t829 + (t110 * t569 + t111 * t630 + t158 * t371) * qJD(2) + t30 * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t41 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t675, qJD(2) * t342 + t712, 0, 0, 0, 0, 0, -t710, -t707, 0, 0, 0, 0, 0, t943, t270 * t670 + t910 - t924, t943, t719, t924 + t736 + t910, t813 + t30 * qJD(2) + (t368 * t512 - t510 * t773) * qJD(3) + t20 * qJD(4) + t17 * qJD(5) + t271 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t568, -t567, t624 * t393, t576, -t479, qJD(2) * t149 - qJD(4) * t262 + t717, qJD(2) * t147 + qJD(4) * t261 + t718, t131, t67, t112, t113, -t479, qJD(2) * t44 + t833 + t930, qJD(2) * t45 + qJD(5) * t55 + t254 - t694 + t832, qJD(2) * t25 + t840 + t930, t816 + t8 * qJD(2) + (-t767 - t770) * qJD(4) + t69 * qJD(5) - t679, qJD(2) * t24 + qJD(5) * t52 + t694 + t735 + t839, t825 + t2 * qJD(2) + t20 * qJD(3) + (t136 * t527 + t137 * t524) * qJD(4) + t9 * qJD(5) + t47 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t67, t112, t113, -t479, qJD(2) * t60 + t830 + t929, qJD(2) * t59 + qJD(4) * t55 + t254 + t691 + t831, qJD(2) * t28 + t837 + t929, t6 * qJD(2) + t69 * qJD(4) + qJD(5) * t616 - t679 + t821, qJD(2) * t31 + qJD(4) * t52 - t691 + t735 + t838, t828 + t4 * qJD(2) + t17 * qJD(3) + t9 * qJD(4) + (-pkin(5) * t126 - qJ(6) * t125) * qJD(5) + t104 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t164 - t661, t112, -t931, qJD(2) * t41 + qJD(3) * t271 + qJD(4) * t47 + qJD(5) * t104 + t834; 0, 0, 0, -t656, -t672, 0, 0, 0, t667, t666, 0, qJD(3) * t343 - t688, t410 * t534 - t380, 0.2e1 * t537 * t576, qJD(4) * t395 - t709, -t381 + t711, -t916, qJD(4) * t148 - t505 * t695 - t751, qJD(3) * t393 + qJD(4) * t146 - t750, -t368 * t708 + t746, -t689 + t748, -t714 - t734, -t716 - t936, -t913, -qJD(4) * t43 - qJD(5) * t57 - t808 + t907, qJD(3) * t285 - qJD(4) * t42 - qJD(5) * t58 - t807, qJD(4) * t26 + qJD(5) * t27 - qJD(6) * t163 - t814 + t907, qJD(3) * t101 - qJD(4) * t7 - qJD(5) * t5 - qJD(6) * t267 - t822, -qJD(3) * t283 - qJD(4) * t23 + qJD(5) * t32 + t305 - t815, -qJD(3) * t29 - qJD(4) * t1 - qJD(5) * t3 - qJD(6) * t40 - t829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t657, t520 * qJD(4), 0, 0, 0, t525 * t693, t525 * t692, -t512 * t626, t670 * t355, 0, 0, 0, qJD(4) * t353 + t498 * t514, qJD(4) * t354 - qJD(5) * t775, qJD(4) * t196 + qJD(5) * t200 - t497 * t510, 0, qJD(4) * t197 + qJD(5) * t201 + qJD(6) * t507, qJD(4) * t89 + qJD(5) * t98 - t371 * t497; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t680, 0, 0, 0, 0, 0, -t655, t677, 0, 0, 0, 0, 0, t898, t681, t898, t720, -t682, t749 - t841; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t571, -t566, t676 + t692, t341, t674, -t523 * t692 - t574, t523 * t693 - t575, -t155, t116, t192, t926, t674, -t602 - t892, t922 - t601, -t606 - t892 (-t766 - t768) * qJD(4) + t184 * qJD(5) - t612 - t738, -t922 - t605 (-t524 * t630 + t527 * t569) * qJD(4) + t77 * qJD(5) + t900 + t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t116, t192, t926, t674, -t577 - t892, -t578 + t922, -t604 - t892, t184 * qJD(4) + qJD(5) * t615 - t738 - t844, -t922 - t603, t77 * qJD(4) + (-pkin(5) * t569 - qJ(6) * t630) * qJD(5) + t900 + t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t573, t192, t263, -t579 + t892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t675, -qJD(2) * t343 - t712, 0, 0, 0, 0, 0, -t381 + t473 + t710, -qJD(2) * t393 - t658 + t707, 0, 0, 0, 0, 0, t949, -t285 * qJD(2) - t268 * t670 + t924, t949, -qJD(2) * t101 - t719, qJD(2) * t283 - t734 - t924, qJD(2) * t29 - qJD(4) * t19 - qJD(5) * t16 - qJD(6) * t264 - t813; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t680, 0, 0, 0, 0, 0, t655, -t677, 0, 0, 0, 0, 0, -t898, -t681, -t898, -t720, t682, t749 + t841; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, t600, 0, 0, 0, 0, 0, t926, t193, t926, 0, t192 (t765 - t769) * qJD(4) + t187 * qJD(5) + t608 + t497; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t926, t193, t926, 0, t192, t187 * qJD(4) - qJD(5) * t407 + t497 - t817; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t512 * t670 - t685; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, t567, -qJD(2) * t395 + t410 * t537, qJD(2) * t390 + t503 * t655, -t479, -qJD(2) * t148 + qJD(3) * t390 - t717, -qJD(2) * t146 + t503 * t695 - t718, -t131, -t67, t160, t161, -t479, qJD(2) * t43 - t833 + t928, qJD(2) * t42 + qJD(5) * t54 + t247 - t832, -qJD(2) * t26 - t840 + t928, qJD(2) * t7 + qJD(5) * t70 - t816, qJD(2) * t23 + qJD(5) * t53 + t737 - t839, qJD(2) * t1 + qJD(3) * t19 + qJD(5) * t10 - qJD(6) * t46 - t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t571, t566, -t676, t678, -t674, t574, t575, t155, -t116, t246, t902, -t674, t602, t601, t606, qJD(5) * t185 + t612, t605, qJD(5) * t78 - t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t678, t503 * t704, 0, 0, 0, 0, 0, t902, t248, t902, 0, t246, -qJD(5) * t186 - t608; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t665, -t531, -t665, 0, t671, qJD(5) * t436 + t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t932, -t531 - t631, t932, t607, -t632 + t671, t517 + (-pkin(5) * t536 + qJ(6) * t539) * t846 + t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t453, t524 * t670 - t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t67, t160, t161, -t479, qJD(2) * t57 - t830 + t927, qJD(2) * t58 - qJD(4) * t54 + t247 - t831, -qJD(2) * t27 - t837 + t927, qJD(2) * t5 - qJD(4) * t70 - t821, -qJD(2) * t32 - qJD(4) * t53 + t737 - t838, qJ(6) * t491 + qJD(2) * t3 + qJD(3) * t16 - qJD(4) * t10 - t828; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t116, t246, t902, -t674, t577, t578, t604, -qJD(4) * t185 + t844, t603, -qJD(4) * t78 - t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t902, t248, t902, 0, t246, qJD(4) * t186 + t817; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t934, t631, t934, -t607, qJD(6) + t632, t532 - t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), t532; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t453, t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t163 + t661, t160, t931, qJD(2) * t40 + qJD(3) * t264 + qJD(4) * t46 - t503 * t721 - t834; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t573, t246, -t263, t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t685; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t453, -qJD(4) * t524 - t721 + t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t453, -t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t48;
