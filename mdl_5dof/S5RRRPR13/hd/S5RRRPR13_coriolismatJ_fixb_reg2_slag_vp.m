% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR13_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:39
% EndTime: 2019-12-31 21:47:15
% DurationCPUTime: 23.33s
% Computational Cost: add. (13407->867), mult. (32422->1155), div. (0->0), fcn. (33422->8), ass. (0->628)
t511 = sin(qJ(3));
t512 = sin(qJ(2));
t509 = sin(pkin(5));
t514 = cos(qJ(3));
t829 = t509 * t514;
t883 = cos(pkin(5));
t424 = t511 * t883 + t512 * t829;
t902 = t424 / 0.2e1;
t505 = t511 ^ 2;
t507 = t514 ^ 2;
t483 = t507 - t505;
t515 = cos(qJ(2));
t774 = qJD(2) * t515;
t706 = t509 * t774;
t494 = t883 * t514;
t830 = t509 * t512;
t422 = t511 * t830 - t494;
t199 = t422 * t514 + t424 * t511;
t957 = t199 * qJD(3);
t968 = t483 * t706 - t957;
t461 = t514 * t706;
t444 = t511 * t461;
t819 = t511 * t422;
t677 = -t819 / 0.2e1;
t838 = t424 * t514;
t689 = t838 / 0.2e1;
t564 = t689 + t677;
t926 = t564 * qJD(3);
t967 = t444 + t926;
t737 = t509 * qJD(3);
t714 = t515 * t737;
t501 = t509 ^ 2;
t508 = t515 ^ 2;
t832 = t501 * t508;
t324 = t422 * t830 - t511 * t832;
t958 = qJD(1) * t324;
t966 = t511 * t714 + t958;
t777 = qJD(2) * t509;
t707 = t512 * t777;
t965 = t514 * t707 - t958;
t918 = pkin(3) + pkin(9);
t964 = t918 * t422;
t804 = t918 * t514;
t510 = sin(qJ(5));
t917 = pkin(4) + pkin(8);
t941 = t917 * t514;
t833 = t941 * t510;
t513 = cos(qJ(5));
t449 = t941 * t513;
t550 = -t424 * t830 + t514 * t832;
t929 = qJD(1) * t550;
t963 = t514 * t714 - t929;
t962 = t511 * t707 - t929;
t779 = qJD(1) * t424;
t710 = t422 * t779;
t938 = -qJD(2) * t564 + t710;
t504 = t510 ^ 2;
t506 = t513 ^ 2;
t897 = t506 / 0.2e1;
t436 = (t897 - t504 / 0.2e1) * t514;
t769 = qJD(3) * t513;
t489 = t510 * t769;
t828 = t509 * t515;
t237 = t510 * t828;
t840 = t422 * t513;
t334 = t237 + t840;
t810 = t513 * t334;
t674 = t810 / 0.2e1;
t722 = t513 * t828;
t841 = t422 * t510;
t335 = t722 - t841;
t850 = t335 * t510;
t565 = t850 / 0.2e1 + t674;
t961 = -qJD(1) * t565 + qJD(2) * t436 + t489;
t736 = t514 * qJD(2);
t487 = t511 * t736;
t923 = qJD(1) * t564 + t487;
t778 = qJD(1) * t515;
t705 = t509 * t778;
t633 = t511 * t705;
t611 = t422 * t633;
t960 = t611 - t967;
t959 = -t926 - t611;
t809 = t513 * t335;
t851 = t334 * t510;
t562 = t809 / 0.2e1 - t851 / 0.2e1;
t945 = t562 * t514;
t956 = t945 * qJD(5);
t820 = t511 * qJ(4);
t887 = t514 * pkin(3);
t616 = -t820 - t887;
t457 = -pkin(2) + t616;
t720 = pkin(1) * t883;
t431 = pkin(7) * t830 - t515 * t720;
t398 = -pkin(2) * t883 + t431;
t839 = t424 * qJ(4);
t524 = t398 - t839;
t890 = t422 * pkin(3);
t214 = t524 + t890;
t866 = t214 * t511;
t939 = t457 * t902 + t866 / 0.2e1;
t955 = -qJD(1) * t199 + qJD(2) * t483;
t919 = t424 ^ 2;
t664 = t422 ^ 2 - t919;
t954 = qJD(1) * t664 - qJD(2) * t199;
t953 = qJD(2) * t324 - t424 * t714;
t780 = qJD(1) * t335;
t952 = -qJD(2) * t945 - qJD(3) * t565 + t334 * t780;
t775 = qJD(2) * t513;
t951 = -t507 * t510 * t775 - qJD(1) * t945 + qJD(3) * t436;
t771 = qJD(3) * t510;
t781 = qJD(1) * t334;
t950 = t424 * (-t771 + t781);
t165 = t524 + t964;
t891 = pkin(8) * t512;
t622 = -pkin(2) * t515 - t891;
t593 = -pkin(1) + t622;
t376 = t593 * t829;
t645 = t512 * t720;
t542 = pkin(8) * t883 + t645;
t888 = t424 * pkin(4);
t518 = t511 * t542 - t376 + t888 + (pkin(7) * t511 + t918) * t828;
t83 = t513 * t165 + t510 * t518;
t884 = t83 * t513;
t82 = t165 * t510 - t513 * t518;
t885 = t82 * t510;
t949 = t424 * (t884 + t885);
t463 = -qJD(3) + t705;
t948 = t463 * t424;
t947 = t510 * t919;
t946 = t513 * t919;
t734 = pkin(7) * t828;
t253 = t511 * (t542 + t734) - t376;
t186 = -t253 - t888;
t682 = t828 / 0.2e1;
t944 = (-t918 * t682 - t186 / 0.2e1) * t514;
t776 = qJD(2) * t511;
t943 = (t776 + t779) * t828;
t435 = -pkin(2) - t820 - t804;
t942 = t917 * t511;
t329 = t513 * t435 + t510 * t942;
t854 = t329 * t513;
t328 = t435 * t510 - t513 * t942;
t855 = t328 * t510;
t158 = (t854 + t855) * t511;
t632 = t514 * t705;
t347 = t424 * t632;
t940 = t347 + t444;
t895 = t511 / 0.2e1;
t575 = t329 * t902 + t83 * t895;
t684 = t830 / 0.2e1;
t412 = t511 * t684 - t494 / 0.2e1;
t931 = -qJD(5) * t412 - t938;
t735 = t514 * qJD(5);
t930 = t735 / 0.2e1 + t923;
t482 = t504 - t506;
t662 = qJD(3) * t482;
t927 = qJD(3) * t664;
t925 = t565 * qJD(5);
t857 = t253 * t514;
t217 = pkin(3) * t828 + t253;
t861 = t217 * t514;
t529 = t514 * t542;
t553 = t511 * t593;
t805 = t514 * t515;
t254 = t529 + (pkin(7) * t805 + t553) * t509;
t481 = qJ(4) * t828;
t216 = -t481 + t254;
t862 = t216 * t511;
t924 = -t857 / 0.2e1 + t861 / 0.2e1 - t862 / 0.2e1;
t921 = qJD(2) * t550 - t422 * t714;
t889 = t422 * pkin(4);
t187 = t254 - t889;
t814 = t513 * t187;
t842 = t422 * qJ(4);
t218 = t424 * t918 + t842;
t825 = t510 * t218;
t104 = t814 - t825;
t916 = -t104 / 0.2e1;
t812 = t513 * t218;
t827 = t510 * t187;
t105 = t812 + t827;
t915 = t105 / 0.2e1;
t407 = t511 * t431;
t432 = (pkin(2) * t512 - pkin(8) * t515) * t509;
t806 = t514 * t432;
t282 = t407 + t806;
t492 = pkin(3) * t830;
t258 = -t282 - t492;
t195 = (pkin(4) * t805 - pkin(9) * t512) * t509 + t258;
t813 = t513 * t195;
t817 = t511 * t515;
t727 = t509 * t817;
t583 = pkin(3) * t727 + t645;
t882 = qJ(4) * t514;
t615 = pkin(9) * t511 - t882;
t257 = (pkin(7) + t615) * t828 + t583;
t824 = t510 * t257;
t110 = t813 - t824;
t914 = -t110 / 0.2e1;
t811 = t513 * t257;
t826 = t510 * t195;
t111 = t811 + t826;
t913 = -t111 / 0.2e1;
t168 = -t481 + t187;
t912 = t168 / 0.2e1;
t911 = -t334 / 0.2e1;
t910 = t334 / 0.2e1;
t909 = -t335 / 0.2e1;
t908 = t335 / 0.2e1;
t451 = t511 * t722;
t382 = t510 * t830 - t451;
t907 = t382 / 0.2e1;
t815 = t512 * t513;
t383 = (t510 * t817 + t815) * t509;
t906 = t383 / 0.2e1;
t905 = -t407 / 0.2e1;
t904 = -t422 / 0.2e1;
t903 = -t424 / 0.2e1;
t452 = t514 * t237;
t901 = -t452 / 0.2e1;
t900 = -t457 / 0.2e1;
t899 = t941 / 0.2e1;
t898 = -t492 / 0.2e1;
t896 = -t511 / 0.2e1;
t894 = -t514 / 0.2e1;
t893 = t514 / 0.2e1;
t892 = t918 / 0.2e1;
t8 = -t104 * t82 + t105 * t83 + t168 * t186;
t886 = t8 * qJD(1);
t29 = -t168 * t828 - t949;
t881 = qJD(1) * t29;
t44 = t168 * t334 + t424 * t82;
t880 = qJD(1) * t44;
t45 = -t168 * t335 - t424 * t83;
t879 = qJD(1) * t45;
t295 = pkin(3) * t424 + t842;
t728 = t254 * t828;
t867 = t214 * t424;
t68 = -t295 * t422 - t728 - t867;
t878 = qJD(1) * t68;
t729 = t253 * t828;
t69 = t214 * t422 - t295 * t424 + t729;
t877 = qJD(1) * t69;
t96 = -t216 * t828 - t867;
t876 = qJD(1) * t96;
t875 = t104 * t513;
t874 = t105 * t510;
t408 = t511 * t432;
t409 = t514 * t431;
t283 = -t409 + t408;
t816 = t512 * qJ(4);
t215 = (-pkin(4) * t817 + t816) * t509 + t283;
t15 = -t110 * t82 + t111 * t83 + t168 * t215;
t873 = t15 * qJD(1);
t16 = t104 * t335 + t105 * t334 + t949;
t872 = t16 * qJD(1);
t871 = t168 * t510;
t870 = t168 * t513;
t17 = t110 * t335 + t111 * t334 - t382 * t83 + t383 * t82;
t869 = t17 * qJD(1);
t21 = -t186 * t334 + t82 * t422 + (t104 - t870) * t424;
t868 = t21 * qJD(1);
t865 = t214 * t514;
t864 = t215 * t510;
t863 = t215 * t513;
t22 = -t186 * t335 + t83 * t422 + (-t105 + t871) * t424;
t860 = t22 * qJD(1);
t725 = t509 * t805;
t23 = t110 * t424 + t168 * t382 - t215 * t334 - t725 * t82;
t859 = t23 * qJD(1);
t24 = -t111 * t424 + t168 * t383 - t215 * t335 - t725 * t83;
t858 = t24 * qJD(1);
t856 = t254 * t511;
t500 = pkin(3) * t511;
t437 = t500 + t615;
t823 = t510 * t437;
t330 = t449 - t823;
t853 = t330 * t513;
t808 = t513 * t437;
t331 = t808 + t833;
t852 = t331 * t510;
t849 = t335 * t511;
t848 = t382 * t510;
t847 = t383 * t513;
t846 = t398 * t514;
t40 = t214 * t295 - t216 * t253 + t217 * t254;
t845 = t40 * qJD(1);
t726 = t509 * t816;
t256 = -t726 - t283;
t303 = (pkin(7) - t882) * t828 + t583;
t41 = t214 * t303 - t216 * t256 + t217 * t258;
t844 = t41 * qJD(1);
t666 = t253 / 0.2e1 - t217 / 0.2e1;
t683 = -t828 / 0.2e1;
t696 = -t856 / 0.2e1;
t42 = t696 + (qJ(4) * t683 + t216 / 0.2e1) * t511 + (pkin(3) * t683 + t666) * t514;
t843 = t42 * qJD(1);
t270 = t424 * t510;
t46 = (-t216 + t254) * t424 + (-t217 + t253) * t422;
t836 = t46 * qJD(1);
t467 = t500 - t882;
t835 = t467 * t422;
t834 = t941 * t424;
t831 = t501 * t512;
t822 = t510 * t511;
t821 = t510 * t514;
t818 = t511 * t513;
t807 = t513 * t514;
t52 = t256 * t422 + t258 * t424 + (t861 - t862) * t828;
t803 = t52 * qJD(1);
t53 = -t303 * t424 + (t216 * t512 + (t256 - t865) * t515) * t509;
t802 = t53 * qJD(1);
t54 = -t303 * t422 + (t217 * t512 + (-t258 - t866) * t515) * t509;
t801 = t54 * qJD(1);
t55 = -t282 * t424 - t283 * t422 + (-t856 + t857) * t828;
t800 = t55 * qJD(1);
t433 = t645 + t734;
t60 = -t253 * t282 + t254 * t283 + t398 * t433;
t799 = t60 * qJD(1);
t86 = t253 * t830 + t282 * t828 - t398 * t727 - t422 * t433;
t798 = t86 * qJD(1);
t87 = t433 * t424 + (-t254 * t512 + (t283 + t846) * t515) * t509;
t797 = t87 * qJD(1);
t358 = t807 * t270;
t693 = -t849 / 0.2e1;
t694 = t334 * t895;
t89 = -t358 + (t693 + t907) * t513 + (t694 + t906) * t510;
t796 = t89 * qJD(1);
t121 = (-t810 + t850) * t424;
t795 = qJD(1) * t121;
t131 = -t398 * t422 - t729;
t794 = qJD(1) * t131;
t132 = -t398 * t424 - t728;
t793 = qJD(1) * t132;
t724 = t334 * t828;
t171 = -t724 - t947;
t792 = qJD(1) * t171;
t184 = t335 * t422 + t947;
t791 = qJD(1) * t184;
t185 = -t334 * t422 + t946;
t790 = qJD(1) * t185;
t230 = t422 * t725 + t424 * t727;
t787 = qJD(1) * t230;
t231 = t199 * t828;
t786 = qJD(1) * t231;
t773 = qJD(3) * t253;
t772 = qJD(3) * t254;
t770 = qJD(3) * t511;
t499 = qJD(3) * t514;
t768 = qJD(4) * t510;
t767 = qJD(4) * t511;
t766 = qJD(4) * t513;
t765 = qJD(4) * t514;
t764 = qJD(5) * t424;
t763 = qJD(5) * t510;
t762 = qJD(5) * t511;
t761 = qJD(5) * t513;
t760 = qJD(5) * t918;
t116 = t334 * t383 + t335 * t382;
t759 = t116 * qJD(1);
t118 = -t809 + t851;
t122 = t118 * t424;
t758 = t122 * qJD(1);
t392 = -t838 / 0.2e1;
t543 = t392 * t504 + t510 * t693;
t691 = -t847 / 0.2e1;
t135 = t691 + t543;
t757 = t135 * qJD(1);
t351 = -t848 / 0.2e1;
t609 = t392 * t506 + t511 * t674;
t136 = t351 + t609;
t756 = t136 * qJD(1);
t563 = t684 + t689;
t692 = t849 / 0.2e1;
t146 = t692 - t451 / 0.2e1 + t563 * t510;
t755 = t146 * qJD(1);
t680 = t822 / 0.2e1;
t541 = (t815 / 0.2e1 + t515 * t680) * t509;
t673 = -t807 / 0.2e1;
t549 = t424 * t673 + t694;
t148 = t541 - t549;
t754 = t148 * qJD(1);
t151 = -t382 * t424 + t514 * t724;
t753 = t151 * qJD(1);
t723 = t335 * t828;
t152 = t383 * t424 - t514 * t723;
t752 = t152 * qJD(1);
t378 = t424 * t822;
t255 = t378 / 0.2e1 + t424 * t680;
t212 = t514 * t722 - t255;
t751 = t212 * qJD(1);
t213 = t723 + t946;
t750 = t213 * qJD(1);
t379 = t424 * t818;
t307 = -t452 - t379;
t747 = t307 * qJD(1);
t345 = pkin(1) * t831 + t433 * t883;
t746 = t345 * qJD(1);
t346 = pkin(1) * t501 * t515 - t431 * t883;
t745 = t346 * qJD(1);
t410 = t422 * qJD(3);
t744 = t422 * qJD(4);
t741 = t436 * qJD(5);
t438 = (-t512 ^ 2 + t508) * t501;
t740 = t438 * qJD(1);
t738 = t505 * qJD(2);
t733 = pkin(8) * t770;
t718 = t510 * t779;
t717 = t513 * t779;
t716 = t510 * t736;
t715 = t513 * t736;
t713 = t510 * t499;
t712 = qJD(4) * t828;
t711 = t511 * t735;
t709 = t424 * t767;
t708 = t501 * t778;
t488 = t511 * t499;
t704 = t510 * t761;
t490 = t511 * t775;
t703 = qJ(4) * t910;
t702 = qJ(4) * t909;
t701 = t870 / 0.2e1;
t690 = t847 / 0.2e1;
t687 = -t834 / 0.2e1;
t686 = t834 / 0.2e1;
t685 = -t830 / 0.2e1;
t681 = -t826 / 0.2e1;
t679 = -t821 / 0.2e1;
t678 = t821 / 0.2e1;
t676 = -t818 / 0.2e1;
t675 = t813 / 0.2e1;
t672 = t807 / 0.2e1;
t671 = t805 / 0.2e1;
t670 = t382 * t892;
t669 = t383 * t892;
t668 = t424 * t892;
t667 = t912 - t187 / 0.2e1;
t665 = t408 / 0.2e1 - t409 / 0.2e1;
t663 = t883 * qJD(1);
t250 = t677 + t563;
t657 = qJD(1) * t250 + t487;
t486 = t510 * t776;
t656 = qJD(1) * t270 + t486;
t271 = t819 / 0.2e1 + t392;
t655 = qJD(1) * t271 - t487;
t650 = pkin(8) * t682;
t649 = -pkin(8) * t817 / 0.2e1;
t648 = qJD(5) + t779;
t647 = qJD(5) + t776;
t646 = (pkin(8) / 0.2e1 + pkin(4) / 0.2e1) * t514;
t638 = t514 * t489;
t637 = t507 * t704;
t636 = t774 * t831;
t635 = t512 * t708;
t634 = t422 * t705;
t631 = t510 * t715;
t630 = t942 / 0.2e1;
t629 = t514 * t683;
t628 = t509 * t671;
t627 = t941 * t683;
t626 = t513 * t682;
t624 = t509 * t663;
t623 = t883 * t777;
t621 = 0.2e1 * t631;
t620 = -t257 / 0.2e1 + t168 * t893;
t619 = t407 / 0.2e1 + t806 / 0.2e1;
t618 = t668 + t218 / 0.2e1;
t617 = t490 + t717;
t519 = t328 * t916 + t329 * t915 + t186 * t899 - t82 * t330 / 0.2e1 + t83 * t331 / 0.2e1 - t942 * t912;
t568 = t111 * t510 / 0.2e1 + t110 * t513 / 0.2e1;
t521 = -t568 * t918 + t215 * qJ(4) / 0.2e1;
t1 = -t519 + t521;
t95 = -t328 * t330 + t329 * t331 - t941 * t942;
t613 = -t1 * qJD(1) + t95 * qJD(2);
t547 = t914 - t575;
t576 = t328 * t903 + t82 * t896;
t548 = t913 + t576;
t567 = t330 * t909 + t331 * t911;
t3 = (t105 * t893 + t547 + t669) * t513 + (t104 * t894 + t548 + t670) * t510 + t567;
t90 = -t330 * t821 + t331 * t807 - t158;
t612 = -qJD(1) * t3 - qJD(2) * t90;
t610 = t514 * t626;
t608 = t839 - t964;
t120 = t329 * t514 - t941 * t822 + (t331 - t833) * t511;
t532 = t329 * t904 + t331 * t902 + t83 * t893;
t570 = qJ(4) * t906 + t863 / 0.2e1;
t9 = (-t908 * t917 + t915) * t511 + (t168 * t896 + t687 - t944) * t510 + t532 + t570;
t607 = -t9 * qJD(1) - t120 * qJD(2);
t606 = t874 + t875;
t605 = -t256 * t514 + t258 * t511;
t604 = -t282 * t511 + t283 * t514;
t602 = t852 + t853;
t601 = t663 + qJD(2);
t534 = t328 * t904 + t330 * t903 + t82 * t893;
t571 = qJ(4) * t907 + t864 / 0.2e1;
t10 = (t686 + t944) * t513 + (-t910 * t917 + t701 + t916) * t511 + t534 + t571;
t119 = -t328 * t514 + (t330 - 0.2e1 * t449) * t511;
t600 = -t10 * qJD(1) + t119 * qJD(2);
t14 = t510 * t548 + t513 * t547 + t627;
t599 = qJD(1) * t14 - qJD(2) * t158;
t232 = -t328 * t511 + t807 * t941;
t533 = t911 * t941 + t576;
t26 = t513 * t620 + t533 + t681;
t598 = -qJD(1) * t26 - qJD(2) * t232;
t233 = -t329 * t511 - t821 * t941;
t531 = t335 * t899 + t575;
t25 = t510 * t620 + t531 + t675;
t597 = -qJD(1) * t25 + qJD(2) * t233;
t343 = t457 * t514 + t467 * t511;
t520 = t865 / 0.2e1 + t295 * t895 + t422 * t900 + t467 * t902;
t56 = (t649 + t816) * t509 + t520 + t665;
t596 = qJD(1) * t56 + qJD(2) * t343;
t344 = -t457 * t511 + t467 * t514;
t592 = t650 - t432 / 0.2e1;
t58 = -t492 + t905 + t835 / 0.2e1 + (-t295 / 0.2e1 + t592) * t514 + t939;
t595 = qJD(1) * t58 - qJD(2) * t344;
t411 = (t504 + t506) * t514 * t511;
t352 = t848 / 0.2e1;
t522 = t334 * t676 + (t897 + t504 / 0.2e1) * t838 + t335 * t680;
t71 = t352 + t690 + t522;
t594 = qJD(1) * t71 + qJD(2) * t411;
t127 = t901 - t379 + (-t840 / 0.2e1 + t911) * t514;
t447 = t483 * t513;
t591 = -qJD(1) * t127 - qJD(2) * t447;
t128 = -t378 + (t626 - t841 / 0.2e1 + t908) * t514;
t445 = t483 * t510;
t590 = -qJD(1) * t128 - qJD(2) * t445;
t586 = qJD(1) * t237 - t771;
t584 = -t738 - t762;
t581 = -t410 + t634;
t306 = -t461 + t634;
t580 = pkin(2) * t903 + t398 * t895;
t579 = -t256 * qJ(4) / 0.2e1 - t258 * pkin(3) / 0.2e1;
t578 = t412 * qJD(1) - t736 / 0.2e1;
t577 = t509 * t622;
t523 = pkin(2) * t422 / 0.2e1 + t846 / 0.2e1 + t509 * t649;
t123 = t523 + t665;
t574 = pkin(2) * t736 - qJD(1) * t123;
t537 = t514 * t592 + t905;
t125 = t537 + t580;
t573 = pkin(2) * t776 - qJD(1) * t125;
t572 = qJD(5) * t629 - t347;
t569 = -t875 / 0.2e1 - t874 / 0.2e1;
t566 = -t853 / 0.2e1 - t852 / 0.2e1;
t561 = t511 * t892 - t882 / 0.2e1;
t528 = -t214 * t467 / 0.2e1 + t295 * t900 + pkin(8) * t696;
t30 = (t862 / 0.2e1 + t666 * t514) * pkin(8) + t528 + t579;
t560 = -t457 * t467 * qJD(2) + t30 * qJD(1);
t92 = t898 + t537 + t939;
t559 = qJD(1) * t92 + t457 * t776;
t145 = t646 + t566;
t517 = (pkin(7) * t671 + t553 / 0.2e1) * t509 - t481 - t889 / 0.2e1 + t529 / 0.2e1;
t36 = t517 + t569;
t503 = qJ(4) * qJD(3);
t558 = qJD(1) * t36 + qJD(2) * t145 + t503;
t557 = t509 * (-t457 * t515 + t891);
t556 = qJD(2) * t271 + t710;
t552 = (t769 - t780) * t424;
t551 = -t513 * t648 - t490;
t546 = pkin(8) * t629 - t619 - t939;
t545 = t437 / 0.2e1 + t561;
t113 = (t810 + t850) * t514;
t133 = t334 ^ 2 - t335 ^ 2;
t544 = qJD(1) * t133 - qJD(2) * t113 - qJD(3) * t118;
t198 = (t461 - t410) * t424;
t540 = (qJD(3) * t424 + t511 * t706) * t422;
t301 = t545 * t510;
t47 = t510 * t618 + t513 * t667 + t702;
t539 = -qJ(4) * t769 - qJD(1) * t47 - qJD(2) * t301;
t302 = t545 * t513;
t49 = -t510 * t667 + t513 * t618 + t703;
t538 = qJ(4) * t771 - qJD(1) * t49 - qJD(2) * t302;
t530 = qJD(3) * t616 + t765;
t446 = t482 * t507;
t527 = -qJD(1) * t113 - qJD(2) * t446 + 0.2e1 * t638;
t526 = -qJD(1) * t118 + t621 + t662;
t363 = t919 + t832;
t525 = qJD(1) * t363 + t424 * t776 - t714;
t502 = qJ(4) * qJD(4);
t498 = pkin(8) * t499;
t495 = t499 / 0.2e1;
t491 = t513 * t499;
t475 = qJD(2) * t684;
t474 = qJD(1) * t685;
t473 = qJD(1) * t684;
t468 = t483 * qJD(3);
t456 = 0.2e1 * t514 * t704;
t439 = qJ(4) * t705 - t503;
t428 = -t499 + t632;
t427 = t463 * t513;
t426 = t463 * t511;
t399 = (t708 - t737 / 0.2e1) * t512;
t348 = -t511 * t779 - t738;
t327 = t461 / 0.2e1 - t412 * qJD(3);
t260 = t271 * qJD(3);
t251 = t684 + t271;
t249 = -t833 - t808 / 0.2e1 + t561 * t513;
t248 = t449 - t823 / 0.2e1 + t561 * t510;
t149 = t541 + t549;
t147 = t424 * t678 + t692 + t510 * t685 + t451 / 0.2e1;
t144 = t646 - t566;
t137 = t352 + t609;
t134 = t690 + t543;
t130 = t334 * t893 + t422 * t672 + t379 + t901;
t129 = t335 * t894 + t422 * t678 + t378 + t610;
t126 = pkin(8) * t628 + t580 + t619;
t124 = t523 - t665;
t117 = t118 * qJD(5);
t112 = t113 * qJD(5);
t93 = t898 + t546;
t88 = -t358 - t513 * t382 / 0.2e1 - t383 * t510 / 0.2e1 - t562 * t511;
t70 = t351 + t691 + t522;
t59 = -t835 / 0.2e1 + t295 * t893 - t492 + t546;
t57 = t511 * t650 - t520 + t665 + t726;
t50 = t513 * t668 + t703 - t871 / 0.2e1 - t812 / 0.2e1 - t827 / 0.2e1;
t48 = t510 * t668 + t702 + t701 - t825 / 0.2e1 + t814 / 0.2e1;
t43 = t856 / 0.2e1 + (-t820 / 0.2e1 - t887 / 0.2e1) * t828 + t924;
t35 = t517 - t569;
t31 = pkin(8) * t924 - t528 + t579;
t28 = t168 * t679 - t824 / 0.2e1 + t675 - t531;
t27 = t168 * t673 - t811 / 0.2e1 + t681 - t533;
t13 = t627 + (-t884 / 0.2e1 - t885 / 0.2e1) * t511 + (-t854 / 0.2e1 - t855 / 0.2e1) * t424 + t568;
t12 = t105 * t896 + t168 * t680 + t186 * t679 + t335 * t630 - t532 + t570 + (-t629 * t918 + t686) * t510;
t11 = t104 * t895 + t168 * t676 + t186 * t672 + t334 * t630 + t513 * t687 - t610 * t918 - t534 + t571;
t4 = t105 * t673 + t328 * t270 / 0.2e1 + t104 * t678 + t82 * t680 + (t670 + t913) * t510 - t567 + (t669 + t914 + t575) * t513;
t2 = t519 + t521;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t636, t438 * qJD(2), t515 * t623, -t636, -t512 * t623, 0, -t345 * qJD(2), -t346 * qJD(2), 0, 0, t198, -qJD(2) * t231 + t927, -t921, t540, -t953, -t636, -qJD(2) * t86 - qJD(3) * t132, qJD(2) * t87 + qJD(3) * t131, qJD(2) * t55, qJD(2) * t60, -t636, t921, t953, t198, -qJD(2) * t230 + t927, t540, qJD(2) * t52 + qJD(3) * t46 + t422 * t712, qJD(2) * t54 + qJD(3) * t68 + t424 * t744, qJD(2) * t53 + qJD(3) * t69 + qJD(4) * t363, qJD(2) * t41 + qJD(3) * t40 + qJD(4) * t96, (-qJD(2) * t383 - qJD(5) * t334 - t424 * t771) * t335, qJD(2) * t116 + qJD(3) * t122 + qJD(5) * t133, qJD(2) * t152 + qJD(3) * t184 + t334 * t764, (-qJD(2) * t382 + qJD(5) * t335 + t424 * t769) * t334, qJD(2) * t151 + qJD(3) * t185 + t335 * t764, t198, qJD(2) * t23 + qJD(3) * t21 - qJD(4) * t171 + qJD(5) * t45, qJD(2) * t24 + qJD(3) * t22 + qJD(4) * t213 + qJD(5) * t44, qJD(2) * t17 + qJD(3) * t16 + qJD(4) * t121, qJD(2) * t15 + qJD(3) * t8 + qJD(4) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t635, t740, t601 * t828, -t635, -t601 * t830, 0, -qJD(2) * t433 - t746, qJD(2) * t431 - t745, 0, 0, t926 + t940, -t786 + t968, t962, t960, t965, -t399, -t798 + (-t433 * t514 + t511 * t577) * qJD(2) + t126 * qJD(3), t797 + (t433 * t511 + t514 * t577) * qJD(2) + t124 * qJD(3), qJD(2) * t604 + t800, t799 + (-t433 * pkin(2) + pkin(8) * t604) * qJD(2), -t399, -t962, -t965, -t260 + t940, -t787 + t968, t960, qJD(2) * t605 + qJD(3) * t43 + t803, t801 + (t303 * t514 + t511 * t557) * qJD(2) + t59 * qJD(3) + t251 * qJD(4), t802 + (-t303 * t511 + t514 * t557) * qJD(2) + t57 * qJD(3) + t709, t844 + (pkin(8) * t605 + t303 * t457) * qJD(2) + t31 * qJD(3) + t93 * qJD(4), t134 * qJD(3) + t956 + (-t716 - t780) * t383, t759 + t88 * qJD(3) - t112 + (-t847 + t848) * t736, t752 + (-t237 * t507 + t383 * t511) * qJD(2) + t129 * qJD(3) + t149 * qJD(5), t137 * qJD(3) - t956 + (t715 - t781) * t382, t753 + (-t382 * t511 - t507 * t722) * qJD(2) + t130 * qJD(3) + t147 * qJD(5), -t572 + t967, t859 + (t110 * t511 + t941 * t382 + (-t328 * t828 + t863) * t514) * qJD(2) + t11 * qJD(3) + t255 * qJD(4) + t28 * qJD(5), t858 + (-t111 * t511 + t941 * t383 + (-t329 * t828 - t864) * t514) * qJD(2) + t12 * qJD(3) + t513 * t709 + t27 * qJD(5), t869 + (t328 * t383 - t329 * t382 + (t110 * t510 - t111 * t513) * t514) * qJD(2) + t4 * qJD(3) + t70 * qJD(4), t873 + (-t110 * t328 + t111 * t329 + t215 * t941) * qJD(2) + t2 * qJD(3) + t13 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t938, t954, t581, t938, t948, t475, qJD(2) * t126 - t772 - t793, qJD(2) * t124 + t773 + t794, 0, 0, t475, -t463 * t422, -t948, -t556, t954, t938, t836 + t43 * qJD(2) + (-t839 + t890) * qJD(3) - t744, qJD(2) * t59 + t772 + t878, qJD(2) * t57 - t712 - t773 + t877, t845 + t31 * qJD(2) + (-pkin(3) * t254 - qJ(4) * t253) * qJD(3) + t216 * qJD(4), t134 * qJD(2) + t510 * t552 + t925, t88 * qJD(2) - t424 * t662 - t117 + t758, qJD(2) * t129 - t410 * t513 + t791, t137 * qJD(2) + t513 * t950 - t925, qJD(2) * t130 + t410 * t510 + t790, t931, t868 + t11 * qJD(2) + (t186 * t510 - t513 * t608) * qJD(3) - t334 * qJD(4) + t48 * qJD(5), t860 + t12 * qJD(2) + (t186 * t513 + t510 * t608) * qJD(3) - qJD(4) * t335 + t50 * qJD(5), qJD(2) * t4 - qJD(3) * t606 + t872, t886 + t2 * qJD(2) + (t186 * qJ(4) - t606 * t918) * qJD(3) + t35 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t581, qJD(2) * t251 + t710, t525, qJD(2) * t93 + qJD(3) * t216 + t876, 0, 0, 0, 0, 0, 0, qJD(2) * t255 - qJD(3) * t334 - t792, -qJD(3) * t335 + t424 * t490 + t750, qJD(2) * t70 + t795, qJD(2) * t13 + qJD(3) * t35 + t881; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t952, t544, t149 * qJD(2) + t334 * t648, t952, t147 * qJD(2) + t335 * t648, t327, qJD(2) * t28 + qJD(3) * t48 - qJD(5) * t83 + t879, qJD(2) * t27 + qJD(3) * t50 + qJD(5) * t82 + t880, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t635, -t740, -t515 * t624, t635, t512 * t624, 0, t746, t745, 0, 0, -t347 + t926, -t957 + t786, -t963, t959, t966, t399, qJD(3) * t125 + t798, qJD(3) * t123 - t797, -t800, -t799, t399, t963, -t966, -t260 - t347, -t957 + t787, t959, -qJD(3) * t42 - t514 * t712 - t803, -qJD(3) * t58 - qJD(4) * t250 - t801, -qJD(3) * t56 + t709 - t802, -qJD(3) * t30 - qJD(4) * t92 - t844, qJD(3) * t135 + t383 * t780 + t956, qJD(3) * t89 - t112 - t759, -qJD(3) * t128 - qJD(5) * t148 - t752, qJD(3) * t136 + t382 * t781 - t956, -qJD(3) * t127 + qJD(5) * t146 - t753, t926 + t572, -qJD(3) * t10 - qJD(4) * t212 - qJD(5) * t25 - t859, -qJD(3) * t9 - qJD(4) * t307 - qJD(5) * t26 - t858, -qJD(3) * t3 + qJD(4) * t71 - t869, -qJD(3) * t1 + qJD(4) * t14 - t873; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488, t468, 0, -t488, 0, 0, -pkin(2) * t770, -pkin(2) * t499, 0, 0, 0, 0, 0, t488, t468, -t488, 0, qJD(3) * t344 - t511 * t765, -qJD(3) * t343 + qJD(4) * t505, (qJD(3) * t467 - t767) * t457, -t488 * t504 + t637, -qJD(5) * t446 - 0.2e1 * t511 * t638, -qJD(3) * t445 - t513 * t711, -t488 * t506 - t637, -qJD(3) * t447 + t510 * t711, t488, qJD(3) * t119 + qJD(5) * t233 + t505 * t768, -qJD(3) * t120 - qJD(5) * t232 + t505 * t766, -qJD(3) * t90 + qJD(4) * t411, qJD(3) * t95 - qJD(4) * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t923, t955, -t428, -t923, t426, t474, -t498 - t573, -t574 + t733, 0, 0, t474, t428, -t426, -t655, t955, -t923, t530 - t843, t498 - t595, -t596 - t733, pkin(8) * t530 - t560, t757 - t741 + (-t504 * t736 + t489) * t511, t796 + t456 + (-0.2e1 * t631 - t662) * t511, t491 + t590, t756 + t741 + (-t506 * t736 - t489) * t511, t591 - t713, t930, (-t513 * t804 + (-qJ(4) * t513 - t510 * t917) * t511) * qJD(3) + t513 * t765 + t248 * qJD(5) + t600, (t510 * t804 + (qJ(4) * t510 - t513 * t917) * t511) * qJD(3) - t510 * t765 + t249 * qJD(5) + t607, -qJD(3) * t602 + t612, (-qJ(4) * t942 - t602 * t918) * qJD(3) + t144 * qJD(4) + t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t428, -t657, -t348, t498 - t559, 0, 0, 0, 0, 0, 0, t510 * t738 + t491 - t751, t513 * t738 - t713 - t747, t594, qJD(3) * t144 + t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t951, t527, -t647 * t807 - t754, t951, t647 * t821 + t755, qJD(1) * t629 + t495, qJD(3) * t248 - qJD(5) * t329 + t597, qJD(3) * t249 + qJD(5) * t328 + t598, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t938, -t954, -t306, -t938, -t943, t475, -qJD(2) * t125 + t793, -qJD(2) * t123 - t794, 0, 0, t475, t306, t943, t556, -t954, -t938, qJD(2) * t42 - t836, qJD(2) * t58 - t878, qJD(2) * t56 - t712 - t877, -qJ(4) * t712 + qJD(2) * t30 - t845, -qJD(2) * t135 + t335 * t718 + t925, -qJD(2) * t89 - t117 - t758, qJD(2) * t128 - t424 * t763 - t791, -qJD(2) * t136 - t334 * t717 - t925, qJD(2) * t127 - t424 * t761 - t790, -t931, qJD(2) * t10 - qJD(4) * t237 + qJD(5) * t47 - t868, qJD(2) * t9 + qJD(5) * t49 - t513 * t712 - t860, qJD(2) * t3 - t872, qJD(2) * t1 + qJD(4) * t36 - t886; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t923, -t955, t632, t923, -t633, t473, t573, t574, 0, 0, t473, -t632, t633, t655, -t955, t923, t843, t595, t596, t560, t487 * t504 - t741 - t757, t511 * t621 + t456 - t796, -t510 * t762 - t590, t487 * t506 + t741 - t756, -t511 * t761 - t591, -t930, qJD(5) * t301 - t600, qJD(5) * t302 - t607, -t612, qJD(4) * t145 - t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t502, -t704, t482 * qJD(5), 0, t704, 0, 0, qJ(4) * t761 + t768, -qJ(4) * t763 + t766, 0, t502; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, -t439, 0, 0, 0, 0, 0, 0, -t586, -t427, 0, t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t961, t526, -t510 * t648 - t486, t961, t551, t578, t510 * t760 - t539, t513 * t760 - t538, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306, qJD(2) * t250 - t710, -t525, qJ(4) * t714 + qJD(2) * t92 - t876, 0, 0, 0, 0, 0, 0, qJD(2) * t212 + qJD(3) * t237 - qJD(5) * t270 + t792, -t750 + t307 * qJD(2) + (t714 - t764) * t513, -qJD(2) * t71 - t795, -qJD(2) * t14 - qJD(3) * t36 - t881; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t632, t657, t348, t559, 0, 0, 0, 0, 0, 0, t510 * t584 + t751, t513 * t584 + t747, -t594, -qJD(3) * t145 - t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463, t439, 0, 0, 0, 0, 0, 0, t586, t427, 0, -t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t656 - t763, t551, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t952, -t544, t148 * qJD(2) - t950, -t952, -t146 * qJD(2) + t552, t327, qJD(2) * t25 - qJD(3) * t47 + qJD(4) * t270 - t879, qJD(2) * t26 - qJD(3) * t49 + t424 * t766 - t880, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t951, -t527, t754 + (t715 + t771) * t511, -t951, -t755 + (-t716 + t769) * t511, qJD(1) * t628 + t495, -qJD(3) * t301 + t510 * t767 - t597, -qJD(3) * t302 + t511 * t766 - t598, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t961, -t526, t486 + t718, -t961, t617, -t578, t539, t538, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t656, t617, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
