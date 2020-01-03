% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP11_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:35
% EndTime: 2019-12-31 22:19:20
% DurationCPUTime: 23.90s
% Computational Cost: add. (13832->952), mult. (34138->1282), div. (0->0), fcn. (35182->8), ass. (0->643)
t533 = cos(qJ(3));
t530 = sin(qJ(3));
t869 = cos(pkin(5));
t673 = t869 * t530;
t528 = sin(pkin(5));
t531 = sin(qJ(2));
t813 = t528 * t531;
t451 = t533 * t813 + t673;
t440 = t451 * qJ(5);
t723 = pkin(1) * t869;
t534 = cos(qJ(2));
t812 = t528 * t534;
t454 = pkin(7) * t812 + t531 * t723;
t424 = t869 * pkin(8) + t454;
t645 = -pkin(2) * t534 - pkin(8) * t531;
t612 = -pkin(1) + t645;
t576 = t612 * t528;
t251 = t424 * t530 - t533 * t576;
t532 = cos(qJ(4));
t241 = t532 * t251;
t449 = t530 * t813 - t869 * t533;
t879 = pkin(9) * t449;
t309 = pkin(3) * t451 + t879;
t529 = sin(qJ(4));
t287 = t529 * t309;
t787 = t241 / 0.2e1 - t287 / 0.2e1;
t925 = t787 - t440;
t225 = pkin(3) * t812 + t251;
t924 = (t251 / 0.2e1 - t225 / 0.2e1) * t529;
t274 = t532 * t449;
t684 = -t274 / 0.2e1;
t731 = t529 * t812;
t819 = t451 * t532;
t349 = -t731 + t819;
t837 = t349 * t533;
t923 = t530 * t684 + t837 / 0.2e1;
t807 = t529 * t530;
t692 = -t807 / 0.2e1;
t796 = t532 * t534;
t729 = t528 * t796;
t820 = t451 * t529;
t347 = t729 + t820;
t840 = t347 * t533;
t571 = t449 * t692 + t840 / 0.2e1;
t523 = t529 ^ 2;
t525 = t532 ^ 2;
t502 = t525 - t523;
t922 = qJD(3) * t502;
t887 = -t530 / 0.2e1;
t699 = t449 * t887;
t883 = t533 / 0.2e1;
t584 = t451 * t883 + t699;
t921 = t584 * qJD(3);
t453 = (pkin(2) * t531 - pkin(8) * t534) * t528;
t437 = t530 * t453;
t452 = pkin(7) * t813 - t534 * t723;
t438 = t533 * t452;
t782 = t437 - t438;
t259 = pkin(9) * t813 + t782;
t249 = t532 * t259;
t871 = t533 * pkin(9);
t873 = t530 * pkin(3);
t488 = -t871 + t873;
t310 = t488 * t812 + t454;
t296 = t529 * t310;
t784 = t249 / 0.2e1 + t296 / 0.2e1;
t802 = t532 * t310;
t809 = t529 * t259;
t785 = -t809 / 0.2e1 + t802 / 0.2e1;
t738 = t533 * qJD(2);
t507 = t530 * t738;
t667 = qJD(1) * t584 + t507;
t745 = t449 * qJD(1);
t582 = -qJD(2) * t584 + t451 * t745;
t920 = t349 ^ 2;
t919 = t449 ^ 2;
t918 = -qJ(5) / 0.2e1;
t174 = t249 + t296;
t730 = t530 * t812;
t662 = qJ(5) * t730;
t156 = t662 + t174;
t917 = t156 / 0.2e1;
t916 = -t347 / 0.2e1;
t915 = t347 / 0.2e1;
t914 = -t349 / 0.2e1;
t872 = t532 * pkin(8);
t874 = t529 * pkin(3);
t878 = pkin(9) * t530;
t370 = t529 * (-pkin(2) - t878) + (-qJ(5) + t872 - t874) * t533;
t913 = -t370 / 0.2e1;
t912 = t370 / 0.2e1;
t643 = -pkin(3) * t533 - t878;
t611 = -pkin(2) + t643;
t465 = t532 * t611;
t880 = pkin(8) * t529;
t371 = -t465 + (pkin(4) + t880) * t533;
t911 = t371 / 0.2e1;
t519 = t530 * qJ(5);
t473 = t529 * t488;
t805 = t530 * t532;
t649 = -pkin(8) * t805 + t473;
t373 = t519 + t649;
t910 = t373 / 0.2e1;
t520 = t530 * pkin(4);
t514 = pkin(8) * t807;
t799 = t532 * t488;
t670 = t514 + t799;
t378 = -t520 - t670;
t909 = t378 / 0.2e1;
t490 = t532 * t813;
t406 = t533 * t731 - t490;
t908 = -t406 / 0.2e1;
t793 = t533 * t534;
t804 = t531 * t529;
t407 = (t532 * t793 + t804) * t528;
t907 = -t407 / 0.2e1;
t806 = t529 * t533;
t736 = pkin(8) * t806;
t408 = -t465 + t736;
t906 = -t408 / 0.2e1;
t905 = t408 / 0.2e1;
t797 = t532 * t533;
t735 = pkin(8) * t797;
t409 = t529 * t611 + t735;
t904 = t409 / 0.2e1;
t512 = pkin(4) * t807;
t868 = qJ(5) * t532;
t428 = t512 + (pkin(8) - t868) * t530;
t903 = t428 / 0.2e1;
t487 = pkin(4) * t529 - t868;
t429 = (pkin(8) + t487) * t533;
t902 = t429 / 0.2e1;
t901 = -t449 / 0.2e1;
t900 = t449 / 0.2e1;
t899 = -t451 / 0.2e1;
t811 = t529 * qJ(5);
t882 = pkin(4) * t532;
t638 = t811 + t882;
t455 = t638 * t530;
t898 = -t455 / 0.2e1;
t897 = t455 / 0.2e1;
t460 = t473 / 0.2e1;
t478 = -pkin(3) - t638;
t896 = -t478 / 0.2e1;
t895 = t478 / 0.2e1;
t894 = -t487 / 0.2e1;
t893 = t487 / 0.2e1;
t892 = -t488 / 0.2e1;
t891 = t490 / 0.2e1;
t890 = -t514 / 0.2e1;
t889 = -t529 / 0.2e1;
t888 = t529 / 0.2e1;
t886 = t530 / 0.2e1;
t885 = t532 / 0.2e1;
t884 = -t533 / 0.2e1;
t881 = pkin(8) * t347;
t877 = t449 * pkin(4);
t876 = t451 * pkin(4);
t875 = t451 * pkin(9);
t870 = -qJD(4) / 0.2e1;
t639 = pkin(4) * t347 - qJ(5) * t349;
t121 = t225 + t639;
t203 = pkin(4) * t349 + qJ(5) * t347;
t795 = t533 * t424;
t226 = t795 + (-t534 * pkin(9) + t530 * t612) * t528;
t423 = -t869 * pkin(2) + t452;
t644 = t449 * pkin(3) - t875;
t543 = t423 + t644;
t119 = t226 * t529 - t532 * t543;
t862 = t119 * t449;
t45 = t121 * t347 - t203 * t349 - t862;
t867 = qJD(1) * t45;
t859 = t121 * t349;
t120 = t532 * t226 + t529 * t543;
t860 = t120 * t449;
t46 = t203 * t347 + t859 - t860;
t866 = qJD(1) * t46;
t822 = t449 * qJ(5);
t89 = t120 + t822;
t51 = t449 * t89 - t859;
t865 = qJD(1) * t51;
t56 = -t225 * t347 + t862;
t864 = qJD(1) * t56;
t57 = t225 * t349 - t860;
t863 = qJD(1) * t57;
t861 = t119 * t533;
t858 = t121 * t529;
t857 = t121 * t532;
t786 = -t241 + t287;
t138 = t440 + t786;
t803 = t532 * t309;
t810 = t529 * t251;
t626 = t803 + t810;
t139 = -t626 - t876;
t252 = t530 * t576 + t795;
t172 = -t449 * t487 + t252;
t90 = t119 - t877;
t15 = t121 * t172 + t138 * t89 + t139 * t90;
t856 = t15 * qJD(1);
t16 = -t119 * t89 + t120 * t90 + t121 * t203;
t855 = t16 * qJD(1);
t436 = t530 * t452;
t794 = t533 * t453;
t671 = t436 + t794;
t258 = -pkin(3) * t813 - t671;
t168 = pkin(4) * t406 - qJ(5) * t407 + t258;
t854 = t168 * t529;
t853 = t168 * t532;
t173 = t802 - t809;
t669 = pkin(4) * t730;
t157 = -t173 - t669;
t17 = t121 * t168 + t156 * t89 + t157 * t90;
t852 = t17 * qJD(1);
t20 = (t120 - t89) * t349 + (t119 - t90) * t347;
t851 = t20 * qJD(1);
t21 = -t138 * t347 + t139 * t349 + (t89 * t529 - t90 * t532) * t449;
t850 = t21 * qJD(1);
t22 = -t156 * t347 + t157 * t349 - t406 * t89 + t407 * t90;
t849 = t22 * qJD(1);
t848 = t225 * t532;
t847 = t252 * t532;
t846 = t258 * t529;
t845 = t258 * t532;
t31 = -t172 * t349 + t451 * t89 + (t138 + t857) * t449;
t844 = t31 * qJD(1);
t32 = t172 * t347 - t451 * t90 + (-t139 - t858) * t449;
t843 = t32 * qJD(1);
t33 = -t121 * t407 + t156 * t449 - t168 * t349 + t730 * t89;
t842 = t33 * qJD(1);
t34 = t121 * t406 - t157 * t449 + t168 * t347 - t730 * t90;
t841 = t34 * qJD(1);
t839 = t349 * t529;
t838 = t349 * t532;
t836 = t370 * t533;
t835 = t407 * t529;
t834 = t408 * t449;
t833 = t408 * t533;
t832 = t409 * t451;
t831 = t409 * t533;
t41 = -t119 * t451 + t252 * t347 + (t803 + (-t225 + t251) * t529) * t449;
t830 = t41 * qJD(1);
t42 = -t120 * t451 + t252 * t349 + (-t786 - t848) * t449;
t829 = t42 * qJD(1);
t828 = t423 * t533;
t827 = t428 * t349;
t826 = t428 * t529;
t825 = t428 * t532;
t43 = -t119 * t730 + t173 * t449 + t225 * t406 + t258 * t347;
t824 = t43 * qJD(1);
t44 = -t120 * t730 - t174 * t449 + t225 * t407 + t258 * t349;
t823 = t44 * qJD(1);
t821 = t449 * t533;
t818 = t455 * t529;
t817 = t478 * t349;
t816 = t478 * t529;
t522 = t528 ^ 2;
t527 = t534 ^ 2;
t815 = t522 * t527;
t814 = t522 * t531;
t808 = t529 * t347;
t270 = t529 * t449;
t801 = t532 * t347;
t800 = t532 * t406;
t524 = t530 ^ 2;
t798 = t532 * t524;
t66 = t251 * t813 - t423 * t730 - t454 * t449 + t671 * t812;
t792 = t66 * qJD(1);
t67 = t454 * t451 + (-t252 * t531 + (t782 + t828) * t534) * t528;
t791 = t67 * qJD(1);
t375 = t274 * t807;
t71 = t375 + (-t840 / 0.2e1 + t907) * t532 + (-t837 / 0.2e1 + t406 / 0.2e1) * t529;
t790 = t71 * qJD(1);
t278 = t803 / 0.2e1;
t789 = t810 / 0.2e1 + t278;
t697 = -t813 / 0.2e1;
t781 = t533 * t697 - t673 / 0.2e1;
t526 = t533 ^ 2;
t503 = t526 - t524;
t169 = -t251 * t812 - t423 * t449;
t780 = qJD(1) * t169;
t170 = -t252 * t812 - t423 * t451;
t779 = qJD(1) * t170;
t681 = t797 / 0.2e1;
t564 = t349 * t681 + t525 * t699;
t176 = -t835 / 0.2e1 + t564;
t778 = qJD(1) * t176;
t182 = t347 * t451 - t919 * t529;
t777 = qJD(1) * t182;
t183 = t349 * t451 - t919 * t532;
t776 = qJD(1) * t183;
t695 = -t812 / 0.2e1;
t656 = t529 * t695;
t683 = t274 / 0.2e1;
t184 = t891 + t530 * t683 + (t656 + t914) * t533;
t775 = qJD(1) * t184;
t680 = t796 / 0.2e1;
t560 = (t533 * t680 + t804 / 0.2e1) * t528;
t185 = t560 - t571;
t774 = qJD(1) * t185;
t773 = qJD(1) * t270;
t772 = qJD(1) * t274;
t771 = qJD(1) * t349;
t770 = qJD(1) * t534;
t769 = qJD(2) * t528;
t768 = qJD(2) * t534;
t767 = qJD(3) * t529;
t766 = qJD(3) * t530;
t765 = qJD(3) * t532;
t764 = qJD(3) * t533;
t763 = qJD(4) * t119;
t762 = qJD(4) * t408;
t761 = qJD(4) * t449;
t760 = qJD(4) * t529;
t759 = qJD(4) * t532;
t758 = qJD(4) * t533;
t757 = qJD(5) * t529;
t756 = qJD(5) * t533;
t145 = -t347 * t407 - t349 * t406;
t755 = t145 * qJD(1);
t625 = t801 + t839;
t155 = t625 * t449;
t754 = t155 * qJD(1);
t188 = -t347 * t730 - t406 * t449;
t753 = t188 * qJD(1);
t189 = t349 * t730 + t407 * t449;
t752 = t189 * qJD(1);
t209 = -t451 * t530 - t821;
t234 = t209 * t812;
t751 = t234 * qJD(1);
t329 = -t449 * t813 + t530 * t815;
t750 = t329 * qJD(1);
t330 = -t451 * t813 + t533 * t815;
t749 = t330 * qJD(1);
t748 = t347 * qJD(5);
t361 = pkin(1) * t814 + t454 * t869;
t747 = t361 * qJD(1);
t362 = t522 * pkin(1) * t534 - t452 * t869;
t746 = t362 * qJD(1);
t744 = t449 * qJD(3);
t439 = t449 * qJD(5);
t743 = t451 * qJD(3);
t464 = (-t531 ^ 2 + t527) * t522;
t742 = t464 * qJD(1);
t741 = t528 * qJD(3);
t740 = t530 * qJD(2);
t739 = t532 * qJD(5);
t737 = t890 - t520;
t734 = pkin(9) * t760;
t733 = pkin(9) * t759;
t431 = -t876 / 0.2e1;
t430 = t876 / 0.2e1;
t732 = t871 / 0.2e1;
t728 = t449 * t806;
t727 = t449 * t797;
t726 = -t119 / 0.2e1 + t90 / 0.2e1;
t725 = t89 / 0.2e1 - t120 / 0.2e1;
t724 = t430 + t789;
t722 = t347 * t771;
t721 = t349 * t745;
t719 = t528 * t770;
t718 = t528 * t768;
t717 = t532 * t740;
t716 = t534 * t741;
t715 = t529 * t758;
t714 = t532 * t758;
t713 = t529 * t739;
t712 = t522 * t770;
t711 = t531 * t769;
t710 = t529 * t759;
t709 = t349 * t757;
t708 = t530 * t757;
t506 = t529 * t765;
t707 = t530 * t764;
t115 = -t858 / 0.2e1;
t706 = t225 * t888;
t705 = t347 * t887;
t703 = t349 * t886;
t701 = t835 / 0.2e1;
t700 = t834 / 0.2e1;
t698 = t478 * t886;
t696 = t813 / 0.2e1;
t694 = t812 / 0.2e1;
t693 = t270 / 0.2e1;
t691 = t807 / 0.2e1;
t690 = -t806 / 0.2e1;
t689 = t806 / 0.2e1;
t688 = -t805 / 0.2e1;
t687 = t805 / 0.2e1;
t686 = -t803 / 0.2e1;
t685 = -t801 / 0.2e1;
t682 = -t797 / 0.2e1;
t678 = t912 - t409 / 0.2e1;
t677 = t906 + t911;
t676 = t437 / 0.2e1 - t438 / 0.2e1;
t675 = t890 - t520 / 0.2e1;
t674 = t523 / 0.2e1 - t525 / 0.2e1;
t672 = t869 * qJD(1);
t583 = t685 - t839 / 0.2e1;
t180 = (t695 - t583) * t530;
t485 = t529 * qJD(2) * t798;
t668 = qJD(1) * t180 + t485;
t666 = pkin(8) * t694;
t665 = pkin(3) * t914 + pkin(9) * t684;
t664 = -qJD(4) - t745;
t663 = -qJD(4) + t738;
t661 = t530 * t506;
t660 = t768 * t814;
t659 = t531 * t712;
t658 = t533 * t719;
t657 = t529 * t717;
t655 = t530 * t695;
t654 = t530 * t694;
t653 = t528 * t680;
t650 = t431 - t810 / 0.2e1;
t648 = -t879 / 0.2e1 - t203 / 0.2e1;
t647 = t528 * t672;
t646 = t869 * t769;
t642 = 0.2e1 * t657;
t641 = t487 * t886 + t903;
t640 = pkin(9) * t655;
t637 = t449 * t478 + t875;
t636 = -t478 * t533 + t878;
t635 = -qJD(3) + t719;
t536 = t121 * t902 + t138 * t912 + t139 * t911 + t172 * t903 + t89 * t910 + t90 * t909;
t588 = t156 * t885 + t157 * t888;
t546 = t588 * pkin(9) + t168 * t895;
t1 = -t536 + t546;
t78 = t370 * t373 + t371 * t378 + t428 * t429;
t634 = -t1 * qJD(1) + t78 * qJD(2);
t537 = t119 * t913 + t120 * t911 + t121 * t897 + t203 * t903 + t89 * t906 + t90 * t904;
t599 = qJ(5) * t917 - t157 * pkin(4) / 0.2e1;
t3 = -t537 + t599;
t92 = -t370 * t408 + t371 * t409 + t428 * t455;
t633 = -t3 * qJD(1) + t92 * qJD(2);
t586 = t347 * t910 + t378 * t914;
t593 = t370 * t901 + t89 * t883;
t5 = (pkin(9) * t908 + t139 * t887 + t371 * t900 + t90 * t884 + t917) * t532 + (pkin(9) * t407 / 0.2e1 + t157 / 0.2e1 + t138 * t886 + t593) * t529 + t586;
t96 = -t371 * t797 - t378 * t805 + (t373 * t530 + t836) * t529;
t632 = -t5 * qJD(1) - t96 * qJD(2);
t631 = t530 * t653;
t630 = t115 + pkin(9) * t683 - t817 / 0.2e1;
t118 = -t409 * t805 + (t370 * t532 + (t371 - t408) * t529) * t530;
t597 = pkin(4) * t907 + qJ(5) * t908;
t9 = t678 * t349 + t677 * t347 + (t529 * t726 + t532 * t725) * t530 + t597;
t629 = -t9 * qJD(1) - t118 * qJD(2);
t18 = (-t877 / 0.2e1 + t726) * t532 + (-t822 / 0.2e1 - t725) * t529;
t94 = (pkin(4) * t884 - t677) * t532 + (qJ(5) * t884 + t678) * t529;
t628 = -t18 * qJD(1) + t94 * qJD(2);
t627 = t138 * t532 + t139 * t529;
t624 = t373 * t532 + t378 * t529;
t623 = t428 * t533 + t429 * t530;
t622 = t672 + qJD(2);
t541 = t138 * t883 + t349 * t902 + t370 * t899 + t373 * t901 + t89 * t887;
t585 = t121 * t883 + t428 * t901;
t553 = t172 * t886 + t585;
t569 = pkin(9) * t631 - t854 / 0.2e1 + t407 * t896;
t11 = t532 * t553 + t541 + t569;
t140 = -t370 * t530 + t373 * t533 + t532 * t623;
t621 = -t11 * qJD(1) - t140 * qJD(2);
t540 = t139 * t883 + t347 * t902 + t371 * t899 + t378 * t901 + t90 * t887;
t587 = t853 / 0.2e1 + t406 * t896;
t14 = ((t172 / 0.2e1 + pkin(9) * t694) * t530 + t585) * t529 + t540 + t587;
t142 = -t371 * t530 + t378 * t533 + t529 * t623;
t620 = t14 * qJD(1) + t142 * qJD(2);
t195 = t831 + (t818 + t825) * t530;
t102 = t121 * t688;
t589 = t120 * t883 + t409 * t901;
t544 = -t827 / 0.2e1 + t347 * t898 - t589;
t29 = t102 + (pkin(4) * t812 + t203 * t889) * t530 + t544 + t785;
t619 = qJD(1) * t29 - qJD(2) * t195;
t196 = -t428 * t807 + t455 * t805 - t833;
t108 = -t861 / 0.2e1;
t554 = t349 * t897 + t428 * t916 + t700;
t27 = t108 + (qJ(5) * t812 + t203 * t885 + t115) * t530 + t554 + t784;
t618 = qJD(1) * t27 + qJD(2) * t196;
t213 = t408 * t530 + (t670 - 0.2e1 * t514) * t533;
t555 = t451 * t905 + t670 * t901;
t457 = t529 * t640;
t574 = t457 + pkin(3) * t908 - t845 / 0.2e1;
t23 = (t278 - t881 / 0.2e1 + t924) * t533 + (t119 / 0.2e1 + (pkin(8) * t900 - t252 / 0.2e1) * t529) * t530 + t555 + t574;
t617 = -t23 * qJD(1) - t213 * qJD(2);
t214 = t473 * t533 + (-t409 + t735) * t530;
t573 = t532 * t640 + pkin(3) * t907 + t846 / 0.2e1;
t601 = pkin(8) * t914 - t848 / 0.2e1;
t24 = t449 * t460 + t832 / 0.2e1 + (t601 + t787) * t533 + (t120 / 0.2e1 - t847 / 0.2e1) * t530 + t573;
t616 = -t24 * qJD(1) + t214 * qJD(2);
t222 = t428 * t805 + t836;
t101 = t121 * t687;
t556 = t827 / 0.2e1 + t593;
t575 = pkin(4) * t655 - t785;
t35 = t101 + t556 + t575;
t615 = -qJD(1) * t35 - qJD(2) * t222;
t333 = -t524 * t880 - t833;
t109 = t861 / 0.2e1;
t48 = -t834 / 0.2e1 + t109 + (t881 / 0.2e1 + t706) * t530 - t784;
t614 = qJD(1) * t48 - qJD(2) * t333;
t334 = -pkin(8) * t798 - t831;
t47 = t530 * t601 - t589 + t785;
t613 = qJD(1) * t47 + qJD(2) * t334;
t610 = t663 * t530;
t124 = -t727 + (t529 * t694 - t819 / 0.2e1 + t914) * t530;
t472 = t526 * t532 - t798;
t609 = -qJD(1) * t124 - qJD(2) * t472;
t125 = t728 + (t653 + t820 / 0.2e1 + t915) * t530;
t471 = t503 * t529;
t608 = -qJD(1) * t125 + qJD(2) * t471;
t235 = -t451 ^ 2 + t919;
t607 = qJD(1) * t235 + qJD(2) * t209;
t606 = qJD(1) * t209 + qJD(2) * t503;
t605 = -t738 + t745;
t604 = qJD(1) * t451 + t740;
t603 = t732 - t873 / 0.2e1;
t602 = pkin(2) * t899 + t423 * t886;
t600 = t138 * t918 + t139 * pkin(4) / 0.2e1;
t598 = pkin(4) * t909 + t373 * t918;
t445 = t533 * t696 + t673 / 0.2e1;
t596 = t445 * qJD(1) + t740 / 0.2e1;
t595 = t698 + t732;
t594 = t528 * t645;
t547 = pkin(2) * t900 + t828 / 0.2e1 + pkin(8) * t655;
t164 = t547 + t676;
t592 = pkin(2) * t738 - qJD(1) * t164;
t166 = -t436 / 0.2e1 + (t666 - t453 / 0.2e1) * t533 + t602;
t591 = pkin(2) * t740 - qJD(1) * t166;
t590 = t641 * t529;
t202 = -t808 / 0.2e1 + t838 / 0.2e1;
t581 = t717 + t771;
t177 = -t202 + t781;
t443 = (-0.1e1 / 0.2e1 + t674) * t530;
t580 = qJD(1) * t177 + qJD(2) * t443 - t506;
t456 = t674 * t530;
t579 = qJD(1) * t202 - qJD(2) * t456 + t506;
t194 = t583 * t530;
t578 = qJD(1) * t194 - qJD(3) * t456 - t485;
t577 = qJD(4) + t605;
t572 = t530 * t870 + t667;
t570 = t858 / 0.2e1 + t817 / 0.2e1 + t347 * t893;
t568 = t898 + t595;
t461 = -t473 / 0.2e1;
t493 = pkin(9) * t689;
t134 = t493 - t818 / 0.2e1 - t825 / 0.2e1 - t519 + t461 + (t816 / 0.2e1 + (t894 + pkin(8) / 0.2e1) * t532) * t530;
t354 = t478 * t532 + t487 * t529;
t539 = t648 * t529 - t857 / 0.2e1 + t347 * t895 + t349 * t894;
t40 = t539 + t925;
t567 = -qJD(1) * t40 - qJD(2) * t134 + qJD(3) * t354;
t137 = t590 + (t892 + t568) * t532 + t737;
t355 = t487 * t532 - t816;
t37 = t431 + (-t309 / 0.2e1 + t648) * t532 + t570 + t650;
t566 = -qJD(1) * t37 - qJD(2) * t137 + qJD(3) * t355;
t141 = (-t808 + t838) * t530;
t171 = t347 ^ 2 - t920;
t565 = qJD(1) * t171 - qJD(2) * t141 - qJD(3) * t625;
t494 = pkin(9) * t690;
t345 = pkin(3) * t691 + t460 + t494;
t551 = pkin(3) * t915 + t848 / 0.2e1 + pkin(9) * t693;
t59 = t551 - t787;
t563 = pkin(3) * t765 - qJD(1) * t59 - qJD(2) * t345;
t346 = (t892 + t603) * t532;
t61 = t665 + t686 - t924;
t562 = pkin(3) * t767 - qJD(1) * t61 - qJD(2) * t346;
t535 = (-t529 * t678 + t532 * t677) * pkin(9) + t428 * t893 + t455 * t895;
t55 = t535 + t598;
t548 = -t529 * t725 + t532 * t726;
t538 = t548 * pkin(9) + t121 * t893 + t203 * t895;
t8 = t538 + t600;
t561 = t478 * t487 * qJD(3) + t8 * qJD(1) + t55 * qJD(2);
t192 = t826 / 0.2e1 + (t892 + t595) * t532 + t675;
t53 = t630 + t724;
t559 = -qJD(1) * t53 + qJD(2) * t192 + t478 * t767;
t558 = -qJD(2) * t194 - qJD(3) * t202 + t722;
t557 = qJD(4) * t445 + t582;
t552 = -qJD(4) * t638 + t739;
t470 = t502 * t524;
t550 = qJD(1) * t141 + qJD(2) * t470 + 0.2e1 * t661;
t549 = qJD(1) * t625 + t642 - t922;
t221 = -t349 * t805 + t821;
t250 = t919 + t920;
t545 = qJD(1) * t250 - qJD(2) * t221 + t349 * t767 + t761;
t484 = t524 * t525 + t526;
t542 = -qJD(1) * t221 + qJD(2) * t484 + t661 - t758;
t516 = t766 / 0.2e1;
t489 = qJD(2) * t696;
t486 = t532 * t708;
t477 = -0.2e1 * t530 * t710;
t462 = t799 / 0.2e1;
t448 = t456 * qJD(4);
t444 = t523 * t887 + t525 * t886 + t887;
t425 = (t712 - t741 / 0.2e1) * t531;
t397 = t409 * qJD(4);
t344 = t577 * qJ(5);
t336 = qJD(2) * t654 + t445 * qJD(3);
t308 = t532 * t603 + t462 + t514;
t307 = t494 + t461 + (t874 / 0.2e1 + t872) * t530;
t248 = -qJD(3) * t523 - t529 * t581;
t233 = t532 * t738 - t772;
t218 = -t532 * t663 + t772;
t215 = t221 * qJD(5);
t204 = t209 * qJD(3);
t199 = t202 * qJD(4);
t193 = pkin(9) * t682 + t478 * t688 - t826 / 0.2e1 - t799 / 0.2e1 + t675;
t191 = t194 * qJD(4);
t187 = t533 * t656 + t891 + t923;
t186 = t560 + t571;
t181 = t349 * t692 + t530 * t685 + t655;
t178 = t202 + t781;
t175 = t701 + t564;
t167 = t533 * t666 + t436 / 0.2e1 + t794 / 0.2e1 + t602;
t165 = t547 - t676;
t154 = t529 * t610 - t774;
t153 = t774 + (-t529 * t740 + t765) * t533;
t149 = t625 * qJD(4);
t136 = t532 * t568 + t462 + t590 - t737;
t133 = t493 + t519 + pkin(8) * t688 + t460 - t641 * t532 + (t698 + t898) * t529;
t132 = t141 * qJD(4);
t127 = t451 * t687 + t529 * t654 + t703 + t727;
t126 = t451 * t692 + t631 + t705 - t728;
t117 = t120 * qJD(4);
t95 = t532 * t906 + t370 * t889 + t409 * t888 + t371 * t885 + (-t811 / 0.2e1 - t882 / 0.2e1) * t533;
t70 = t406 * t889 + t407 * t885 + t583 * t533 + t375;
t63 = qJD(2) * t185 - qJD(3) * t274 + t347 * t745;
t62 = t706 + t665 + t789;
t60 = t551 + t787;
t58 = qJD(2) * t186 + t347 * t664;
t54 = t535 - t598;
t52 = t686 + t630 + t650;
t50 = pkin(8) * t703 + t225 * t687 + t589 + t785;
t49 = pkin(8) * t705 + t225 * t692 + t108 + t700 - t784;
t39 = t539 - t925;
t38 = t532 * t648 + t430 + t570 + t724;
t36 = t102 - t556 + t575;
t30 = t203 * t691 + t101 - t544 + t669 + t785;
t28 = t121 * t691 + t203 * t688 + t109 - t554 + t662 + t784;
t26 = t649 * t901 - t832 / 0.2e1 + t786 * t883 + t120 * t887 + t252 * t687 + t225 * t681 + t573 + t923 * pkin(8);
t25 = t571 * pkin(8) + t119 * t887 + t225 * t689 + t252 * t691 + t626 * t884 - t555 + t574;
t19 = pkin(4) * t683 + qJ(5) * t693 + t548;
t13 = t529 * t553 + t457 + t540 - t587;
t12 = t121 * t682 + t172 * t688 + t428 * t683 - t541 + t569;
t10 = t119 * t691 + t120 * t687 + t347 * t905 + t371 * t916 + t89 * t688 + t90 * t692 + t597 + (t904 + t913) * t349;
t7 = t538 - t600;
t6 = t370 * t693 + t138 * t692 + t89 * t690 + t371 * t684 + t139 * t687 + t90 * t681 + (-t800 / 0.2e1 + t701) * pkin(9) - t586 + t588;
t4 = t537 + t599;
t2 = t536 + t546;
t64 = [0, 0, 0, t660, t464 * qJD(2), t534 * t646, -t531 * t646, 0, -t361 * qJD(2), -t362 * qJD(2), (t533 * t718 - t744) * t451, qJD(2) * t234 + qJD(3) * t235, -t330 * qJD(2) + t449 * t716, t329 * qJD(2) + t451 * t716, -t660, -qJD(2) * t66 - qJD(3) * t170, qJD(2) * t67 + qJD(3) * t169, (qJD(2) * t407 - qJD(4) * t347 - t532 * t744) * t349, qJD(2) * t145 + qJD(3) * t155 + qJD(4) * t171, qJD(2) * t189 + qJD(3) * t183 - t347 * t761, qJD(2) * t188 - qJD(3) * t182 - t349 * t761, (t530 * t718 + t743) * t449, qJD(2) * t43 + qJD(3) * t41 + qJD(4) * t57, qJD(2) * t44 + qJD(3) * t42 + qJD(4) * t56, qJD(2) * t34 + qJD(3) * t32 + qJD(4) * t46 - t349 * t748, qJD(2) * t22 + qJD(3) * t21 + qJD(4) * t20 - t347 * t439, qJD(2) * t33 + qJD(3) * t31 + qJD(4) * t45 + qJD(5) * t250, qJD(2) * t17 + qJD(3) * t15 + qJD(4) * t16 + qJD(5) * t51; 0, 0, 0, t659, t742, t622 * t812, -t622 * t813, 0, -qJD(2) * t454 - t747, qJD(2) * t452 - t746, t528 * t604 * t793 + t921, t503 * t718 + t204 + t751, t530 * t711 - t749, t533 * t711 + t750, -t425, -t792 + (-t454 * t533 + t530 * t594) * qJD(2) + t167 * qJD(3), t791 + (t454 * t530 + t533 * t594) * qJD(2) + t165 * qJD(3), qJD(3) * t175 + t407 * t581 + t191, t755 + t70 * qJD(3) - t132 + (-t800 - t835) * t740, t752 + (-t407 * t533 + t524 * t729) * qJD(2) + t127 * qJD(3) + t186 * qJD(4), t753 + (t406 * t533 - t524 * t731) * qJD(2) + t126 * qJD(3) + t187 * qJD(4), -t921 + (qJD(4) / 0.2e1 + t605) * t730, -t173 * t738 + t824 + t25 * qJD(3) + t50 * qJD(4) + (pkin(8) * t406 - t408 * t812 + t846) * t740, t174 * t738 + t823 + t26 * qJD(3) + t49 * qJD(4) + (pkin(8) * t407 - t409 * t812 + t845) * t740, t841 + (t157 * t533 + t428 * t406 + (-t371 * t812 + t854) * t530) * qJD(2) + t13 * qJD(3) + t30 * qJD(4) + t181 * qJD(5), t849 + (-t370 * t406 + t371 * t407 + (-t156 * t529 + t157 * t532) * t530) * qJD(2) + t6 * qJD(3) + t10 * qJD(4) + t186 * qJD(5), t842 + (-t156 * t533 - t428 * t407 + (t370 * t812 - t853) * t530) * qJD(2) + t12 * qJD(3) + t28 * qJD(4) - t215, t852 + (t156 * t370 + t157 * t371 + t168 * t428) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t36 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t582, t607, t635 * t449, t635 * t451, t489, qJD(2) * t167 - qJD(3) * t252 - t779, qJD(2) * t165 + qJD(3) * t251 + t780, qJD(2) * t175 + t199 + (-t767 - t771) * t274, t70 * qJD(2) - t502 * t744 - t149 + t754, qJD(2) * t127 + t529 * t743 + t776, qJD(2) * t126 + t532 * t743 - t777, t557, t830 + t25 * qJD(2) + (t529 * t644 - t847) * qJD(3) + t62 * qJD(4), t829 + t26 * qJD(2) + (t252 * t529 + t532 * t644) * qJD(3) + t60 * qJD(4), t843 + t13 * qJD(2) + (-t172 * t532 - t529 * t637) * qJD(3) + t38 * qJD(4) + t178 * qJD(5), t6 * qJD(2) + qJD(3) * t627 + t19 * qJD(4) + t850, t844 + t12 * qJD(2) + (-t172 * t529 + t532 * t637) * qJD(3) + t39 * qJD(4) + t709, t856 + t2 * qJD(2) + (pkin(9) * t627 + t172 * t478) * qJD(3) + t7 * qJD(4) + t52 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t558, t565, t58, qJD(2) * t187 + t349 * t664, t336, qJD(2) * t50 + qJD(3) * t62 - t117 + t863, qJD(2) * t49 + qJD(3) * t60 + t763 + t864, qJD(2) * t30 + qJD(3) * t38 - t117 + t866, t10 * qJD(2) + t19 * qJD(3) + qJD(4) * t639 - t748 + t851, qJD(2) * t28 + qJD(3) * t39 + t439 - t763 + t867, t855 + t4 * qJD(2) + t7 * qJD(3) + (-pkin(4) * t120 - qJ(5) * t119) * qJD(4) + t89 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t181 + qJD(3) * t178 - t722, t58, t545, qJD(2) * t36 + qJD(3) * t52 + qJD(4) * t89 + t865; 0, 0, 0, -t659, -t742, -t534 * t647, t531 * t647, 0, t747, t746, -t451 * t658 + t921, t204 - t751, -t533 * t716 + t749, t530 * t716 - t750, t425, qJD(3) * t166 + t792, qJD(3) * t164 - t791, qJD(3) * t176 - t407 * t771 + t191, qJD(3) * t71 - t132 - t755, -qJD(3) * t124 - qJD(4) * t185 - t752, -qJD(3) * t125 - qJD(4) * t184 - t753, -t921 + (-t745 + t870) * t730, -qJD(3) * t23 - qJD(4) * t47 - t824, -qJD(3) * t24 - qJD(4) * t48 - t823, qJD(3) * t14 - qJD(4) * t29 - qJD(5) * t180 - t841, -qJD(3) * t5 - qJD(4) * t9 - qJD(5) * t185 - t849, -qJD(3) * t11 - qJD(4) * t27 - t215 - t842, -qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t35 - t852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, t503 * qJD(3), 0, 0, 0, -pkin(2) * t766, -pkin(2) * t764, -t524 * t710 + t525 * t707, -qJD(4) * t470 - 0.2e1 * t533 * t661, -qJD(3) * t472 + t530 * t715, qJD(3) * t471 + t530 * t714, -t707, -qJD(3) * t213 - qJD(4) * t334, qJD(3) * t214 + qJD(4) * t333, qJD(3) * t142 + qJD(4) * t195 - t524 * t713, -qJD(3) * t96 - qJD(4) * t118 + t533 * t708, -qJD(3) * t140 - qJD(4) * t196 + qJD(5) * t484, qJD(3) * t78 + qJD(4) * t92 - qJD(5) * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t667, t606, -t635 * t533, t635 * t530, qJD(1) * t697, -pkin(8) * t764 - t591, pkin(8) * t766 - t592, t778 - t448 + (t525 * t740 + t506) * t533, t790 + t477 + (-0.2e1 * t657 + t922) * t533, t529 * t766 + t609, t530 * t765 + t608, -t572, (t529 * t643 - t735) * qJD(3) + t308 * qJD(4) + t617, (t532 * t643 + t736) * qJD(3) + t307 * qJD(4) + t616, (-t429 * t532 - t529 * t636) * qJD(3) + t136 * qJD(4) + t444 * qJD(5) + t620, qJD(3) * t624 + t95 * qJD(4) + t632, (-t429 * t529 + t532 * t636) * qJD(3) + t133 * qJD(4) + t486 + t621, (pkin(9) * t624 + t429 * t478) * qJD(3) + t54 * qJD(4) + t193 * qJD(5) + t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578, -t550, t154, t532 * t610 - t775, qJD(1) * t655 + t516, qJD(3) * t308 - t397 - t613, qJD(3) * t307 - t614 + t762, qJD(3) * t136 - t397 - t619, t95 * qJD(3) + (-qJ(5) * t805 + t512) * qJD(4) - t708 + t629, qJD(3) * t133 - t618 - t756 - t762, t54 * qJD(3) + (-pkin(4) * t409 - qJ(5) * t408) * qJD(4) + t370 * qJD(5) + t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t444 - t668, t154, t542, qJD(3) * t193 + qJD(4) * t370 + t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t582, -t607, -t605 * t812, -t604 * t812, t489, -qJD(2) * t166 + t779, -qJD(2) * t164 - t780, -qJD(2) * t176 + t532 * t721 + t199, -qJD(2) * t71 - t149 - t754, qJD(2) * t124 + qJD(4) * t274 - t776, qJD(2) * t125 - qJD(4) * t270 + t777, -t557, qJD(2) * t23 + qJD(4) * t61 - t830, qJD(2) * t24 + qJD(4) * t59 - t829, -qJD(2) * t14 + qJD(4) * t37 - qJD(5) * t177 - t843, qJD(2) * t5 + qJD(4) * t18 + qJD(5) * t274 - t850, qJD(2) * t11 + qJD(4) * t40 + t709 - t844, qJD(2) * t1 + qJD(4) * t8 + qJD(5) * t53 - t856; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t667, -t606, t658, -t530 * t719, qJD(1) * t696, t591, t592, -t507 * t525 - t448 - t778, t533 * t642 + t477 - t790, -t609 - t714, -t608 + t715, t572, qJD(4) * t346 - t617, qJD(4) * t345 - t616, qJD(4) * t137 - qJD(5) * t443 - t620, -qJD(4) * t94 - t533 * t739 - t632, qJD(4) * t134 + t486 - t621, qJD(4) * t55 - qJD(5) * t192 - t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t710, t502 * qJD(4), 0, 0, 0, -pkin(3) * t760, -pkin(3) * t759, -qJD(4) * t355 + t713, 0, -qJD(4) * t354 + qJD(5) * t523, (qJD(4) * t487 - t757) * t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579, -t549, t218, t529 * t663 - t773, -t596, -t562 - t733, -t563 + t734, -t566 - t733, t552 - t628, -t567 - t734, pkin(9) * t552 + t561; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t580, t218, -t248, -t559 + t733; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t558, -t565, t63, qJD(2) * t184 + qJD(3) * t270 + t721, t336, qJD(2) * t47 - qJD(3) * t61 - t863, qJD(2) * t48 - qJD(3) * t59 - t864, qJD(2) * t29 - qJD(3) * t37 - t866, qJD(2) * t9 - qJD(3) * t18 - t851, qJD(2) * t27 - qJD(3) * t40 + t439 - t867, qJ(5) * t439 + qJD(2) * t3 - qJD(3) * t8 - t855; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t578, t550, t153, t775 + (-t717 - t767) * t533, qJD(1) * t654 + t516, -qJD(3) * t346 + t613, -qJD(3) * t345 + t614, -qJD(3) * t137 + t619, qJD(3) * t94 - t629, -qJD(3) * t134 + t618 - t756, -qJ(5) * t756 - qJD(3) * t55 - t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, t549, t233, -t529 * t738 + t773, t596, t562, t563, t566, t628, t567, -t561; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t180 + qJD(3) * t177 + t722, t63, -t545, -qJ(5) * t761 + qJD(2) * t35 - qJD(3) * t53 - t865; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t443 + t668, t153, -t542, qJ(5) * t758 + qJD(3) * t192 - t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t580, t233, t248, t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t577, -t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t64;
