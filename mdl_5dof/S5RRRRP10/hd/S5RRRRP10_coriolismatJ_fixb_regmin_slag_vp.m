% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRP10
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
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:57
% EndTime: 2021-01-16 00:36:52
% DurationCPUTime: 21.55s
% Computational Cost: add. (14124->922), mult. (34905->1267), div. (0->0), fcn. (36152->8), ass. (0->625)
t531 = sin(qJ(2));
t863 = cos(pkin(5));
t722 = pkin(1) * t863;
t528 = sin(pkin(5));
t534 = cos(qJ(2));
t812 = t528 * t534;
t446 = pkin(7) * t812 + t531 * t722;
t420 = t863 * pkin(8) + t446;
t530 = sin(qJ(3));
t533 = cos(qJ(3));
t648 = -pkin(2) * t534 - pkin(8) * t531;
t619 = -pkin(1) + t648;
t578 = t619 * t528;
t236 = t530 * t420 - t533 * t578;
t532 = cos(qJ(4));
t227 = t532 * t236;
t813 = t528 * t531;
t441 = t530 * t813 - t863 * t533;
t670 = t863 * t530;
t443 = t533 * t813 + t670;
t303 = pkin(3) * t443 + pkin(9) * t441;
t529 = sin(qJ(4));
t273 = t529 * t303;
t790 = t227 / 0.2e1 - t273 / 0.2e1;
t876 = t441 * pkin(4);
t798 = t533 * t420;
t218 = t798 + (-t534 * pkin(9) + t530 * t619) * t528;
t444 = pkin(7) * t813 - t534 * t722;
t419 = -t863 * pkin(2) + t444;
t647 = t441 * pkin(3) - t443 * pkin(9);
t538 = t419 + t647;
t119 = t529 * t218 - t532 * t538;
t729 = t529 * t812;
t824 = t443 * t532;
t352 = -t729 + t824;
t98 = t352 * qJ(5) + t119;
t77 = -t98 + t876;
t866 = t77 + t98;
t217 = pkin(3) * t812 + t236;
t922 = (t236 / 0.2e1 - t217 / 0.2e1) * t529;
t808 = t529 * t530;
t686 = -t808 / 0.2e1;
t800 = t532 * t534;
t727 = t528 * t800;
t825 = t443 * t529;
t350 = t727 + t825;
t838 = t350 * t533;
t921 = t838 / 0.2e1 + t441 * t686;
t646 = -t533 * pkin(3) - t530 * pkin(9);
t618 = -pkin(2) + t646;
t457 = t532 * t618;
t804 = t532 * qJ(5);
t645 = -t530 * t804 + t457;
t880 = pkin(8) * t529;
t332 = (-pkin(4) - t880) * t533 + t645;
t807 = t529 * t533;
t739 = pkin(8) * t807;
t355 = -t645 + t739;
t920 = t332 + t355;
t318 = -t838 / 0.2e1;
t677 = t800 / 0.2e1;
t805 = t531 * t529;
t560 = (t533 * t677 + t805 / 0.2e1) * t528;
t685 = t808 / 0.2e1;
t536 = t441 * t685 + t318 + t560;
t919 = qJD(2) * t536;
t523 = t529 ^ 2;
t525 = t532 ^ 2;
t503 = t525 - t523;
t918 = t503 * qJD(3);
t883 = t533 / 0.2e1;
t888 = -t530 / 0.2e1;
t588 = t441 * t888 + t443 * t883;
t917 = t588 * qJD(3);
t425 = t530 * t444;
t445 = (pkin(2) * t531 - pkin(8) * t534) * t528;
t796 = t533 * t445;
t788 = t425 / 0.2e1 + t796 / 0.2e1;
t742 = t533 * qJD(2);
t507 = t530 * t742;
t665 = t588 * qJD(1) + t507;
t343 = t350 * qJD(4);
t537 = t560 + t921;
t752 = t441 * qJD(1);
t713 = t350 * t752;
t916 = qJD(2) * t537 - t343 - t713;
t584 = -qJD(2) * t588 + t443 * t752;
t915 = t352 ^ 2;
t914 = t441 ^ 2;
t913 = t77 / 0.2e1;
t258 = t532 * t441;
t274 = t532 * t303;
t811 = t529 * t236;
t668 = t274 + t811;
t875 = t443 * pkin(4);
t109 = qJ(5) * t258 + t668 + t875;
t912 = t109 / 0.2e1;
t868 = t533 * pkin(9);
t873 = t530 * pkin(3);
t483 = -t868 + t873;
t305 = t483 * t812 + t446;
t283 = t532 * t305;
t426 = t530 * t445;
t427 = t533 * t444;
t787 = t426 - t427;
t245 = pkin(9) * t813 + t787;
t810 = t529 * t245;
t159 = t283 - t810;
t728 = t530 * t812;
t666 = pkin(4) * t728;
t795 = t533 * t534;
t400 = (t532 * t795 + t805) * t528;
t831 = t400 * qJ(5);
t117 = t159 + t666 - t831;
t911 = t117 / 0.2e1;
t261 = t274 / 0.2e1;
t910 = -t332 / 0.2e1;
t909 = t332 / 0.2e1;
t333 = t352 * t529;
t697 = -t333 / 0.2e1;
t468 = t532 * t483;
t515 = pkin(8) * t808;
t785 = t515 + t468;
t799 = t533 * qJ(5);
t872 = t530 * pkin(4);
t341 = -t532 * t799 + t785 + t872;
t908 = t341 / 0.2e1;
t907 = t350 / 0.2e1;
t906 = -t352 / 0.2e1;
t801 = t532 * t533;
t737 = pkin(8) * t801;
t402 = t529 * t618 + t737;
t356 = -qJ(5) * t808 + t402;
t904 = t356 / 0.2e1;
t467 = t529 * t483;
t806 = t530 * t532;
t738 = pkin(8) * t806;
t651 = t467 - t738;
t359 = -t529 * t799 + t651;
t903 = t359 / 0.2e1;
t902 = -t400 / 0.2e1;
t901 = -t441 / 0.2e1;
t900 = t441 / 0.2e1;
t899 = t443 / 0.2e1;
t898 = t467 / 0.2e1;
t874 = t529 * pkin(4);
t723 = pkin(8) + t874;
t469 = t723 * t530;
t897 = -t469 / 0.2e1;
t896 = t469 / 0.2e1;
t470 = t723 * t533;
t895 = -t470 / 0.2e1;
t867 = qJ(5) + pkin(9);
t479 = t867 * t529;
t894 = -t479 / 0.2e1;
t480 = t867 * t532;
t893 = -t480 / 0.2e1;
t892 = t480 / 0.2e1;
t870 = t532 * pkin(4);
t519 = -pkin(3) - t870;
t891 = t519 / 0.2e1;
t890 = t523 / 0.2e1;
t889 = -t529 / 0.2e1;
t887 = t530 / 0.2e1;
t886 = -t532 / 0.2e1;
t885 = t532 / 0.2e1;
t884 = -t533 / 0.2e1;
t882 = pkin(4) * t352;
t881 = pkin(8) * t350;
t879 = t350 * pkin(4);
t486 = t532 * t813;
t399 = t533 * t729 - t486;
t878 = t399 * pkin(4);
t877 = t400 * pkin(4);
t871 = t530 * pkin(8);
t869 = t533 * pkin(8);
t865 = t77 * t532;
t120 = t532 * t218 + t529 * t538;
t99 = -t350 * qJ(5) + t120;
t864 = t99 * t529;
t15 = t866 * t350;
t862 = qJD(1) * t15;
t161 = t217 + t879;
t16 = t161 * t882 - t866 * t99;
t861 = qJD(1) * t16;
t43 = -t350 * t99 - t352 * t77;
t860 = qJD(1) * t43;
t50 = pkin(4) * t915 - t161 * t350 + t441 * t98;
t859 = qJD(1) * t50;
t51 = -t441 * t99 + (t161 + t879) * t352;
t858 = qJD(1) * t51;
t56 = t119 * t441 - t217 * t350;
t857 = qJD(1) * t56;
t57 = -t120 * t441 + t217 * t352;
t856 = qJD(1) * t57;
t855 = qJD(4) * t99;
t854 = t161 * t529;
t853 = t161 * t532;
t667 = t425 + t796;
t244 = -pkin(3) * t813 - t667;
t188 = t244 + t878;
t852 = t188 * t529;
t851 = t188 * t532;
t254 = t529 * t441;
t789 = t227 - t273;
t129 = qJ(5) * t254 - t789;
t572 = t530 * t578;
t237 = t572 + t798;
t190 = -pkin(4) * t254 + t237;
t19 = t109 * t77 + t129 * t99 + t161 * t190;
t850 = t19 * qJD(1);
t849 = t190 * t532;
t234 = t532 * t245;
t282 = t529 * t305;
t160 = t234 + t282;
t833 = t399 * qJ(5);
t130 = t160 - t833;
t20 = t117 * t77 + t130 * t99 + t161 * t188;
t848 = t20 * qJD(1);
t21 = -t109 * t352 - t129 * t350 + (t864 + t865) * t441;
t847 = t21 * qJD(1);
t846 = t217 * t532;
t22 = -t117 * t352 - t130 * t350 - t399 * t99 - t400 * t77;
t845 = t22 * qJD(1);
t844 = t244 * t529;
t843 = t244 * t532;
t27 = t190 * t350 + t77 * t443 + (t109 - t854) * t441;
t842 = t27 * qJD(1);
t841 = t332 * t532;
t34 = t117 * t441 + t161 * t399 + t188 * t350 + t728 * t77;
t840 = t34 * qJD(1);
t35 = t190 * t352 - t99 * t443 + (-t129 - t853) * t441;
t839 = t35 * qJD(1);
t837 = t352 * t532;
t836 = t356 * t533;
t36 = -t130 * t441 + t161 * t400 + t188 * t352 - t728 * t99;
t835 = t36 * qJD(1);
t39 = -t119 * t443 + t237 * t350 + (t274 + (-t217 + t236) * t529) * t441;
t834 = t39 * qJD(1);
t40 = -t120 * t443 + t237 * t352 + (t789 - t846) * t441;
t832 = t40 * qJD(1);
t830 = t400 * t529;
t829 = t402 * t443;
t41 = -t119 * t728 + t159 * t441 + t217 * t399 + t244 * t350;
t828 = t41 * qJD(1);
t827 = t419 * t533;
t42 = -t120 * t728 - t160 * t441 + t217 * t400 + t244 * t352;
t826 = t42 * qJD(1);
t823 = t469 * t532;
t822 = t470 * t532;
t821 = t479 * t530;
t820 = t479 * t533;
t819 = t480 * t530;
t818 = t480 * t533;
t817 = t519 * t529;
t816 = t519 * t532;
t522 = t528 ^ 2;
t527 = t534 ^ 2;
t815 = t522 * t527;
t814 = t522 * t531;
t809 = t529 * t350;
t803 = t532 * t350;
t524 = t530 ^ 2;
t802 = t532 * t524;
t797 = t533 * t441;
t82 = t236 * t813 - t419 * t728 - t446 * t441 + t667 * t812;
t794 = t82 * qJD(1);
t83 = t446 * t443 + (-t237 * t531 + (t787 + t827) * t534) * t528;
t793 = t83 * qJD(1);
t726 = t529 * t806;
t369 = t441 * t726;
t97 = t369 + (t318 + t902) * t532 + (t352 * t884 + t399 / 0.2e1) * t529;
t792 = t97 * qJD(1);
t791 = t811 / 0.2e1 + t261;
t477 = t529 * t507;
t744 = t532 * qJD(3);
t786 = t533 * t744 - t477;
t502 = t525 + t523;
t526 = t533 ^ 2;
t504 = t526 - t524;
t153 = -t237 * t812 - t419 * t443;
t784 = qJD(1) * t153;
t154 = -t236 * t812 - t419 * t441;
t783 = qJD(1) * t154;
t782 = qJD(1) * t534;
t781 = qJD(2) * t528;
t780 = qJD(2) * t534;
t779 = qJD(3) * t530;
t778 = qJD(3) * t533;
t777 = qJD(4) * t352;
t776 = qJD(4) * t356;
t775 = qJD(4) * t480;
t774 = qJD(4) * t529;
t521 = qJD(4) * t532;
t773 = qJD(4) * t533;
t772 = qJD(5) * t352;
t771 = qJD(5) * t441;
t134 = -t350 * t400 - t352 * t399;
t770 = t134 * qJD(1);
t139 = (t803 + t333) * t441;
t769 = t139 * qJD(1);
t678 = t801 / 0.2e1;
t692 = t525 * t888;
t564 = t352 * t678 + t441 * t692;
t163 = -t830 / 0.2e1 + t564;
t768 = t163 * qJD(1);
t164 = t350 * t443 - t914 * t529;
t767 = t164 * qJD(1);
t165 = t352 * t443 - t914 * t532;
t766 = t165 * qJD(1);
t688 = t812 / 0.2e1;
t606 = t529 * t688 + t906;
t680 = t258 / 0.2e1;
t653 = t530 * t680;
t169 = -t486 / 0.2e1 + t653 + t606 * t533;
t765 = t169 * qJD(1);
t689 = -t812 / 0.2e1;
t170 = t486 / 0.2e1 + t653 + (t529 * t689 + t906) * t533;
t764 = t170 * qJD(1);
t763 = t536 * qJD(1);
t762 = t537 * qJD(1);
t175 = -t350 * t728 - t399 * t441;
t761 = t175 * qJD(1);
t176 = t352 * t728 + t400 * t441;
t760 = t176 * qJD(1);
t201 = -t443 * t530 - t797;
t223 = t201 * t812;
t759 = t223 * qJD(1);
t758 = t254 * qJD(1);
t334 = -t441 * t813 + t530 * t815;
t757 = t334 * qJD(1);
t335 = -t443 * t813 + t533 * t815;
t756 = t335 * qJD(1);
t755 = t352 * qJD(1);
t360 = pkin(1) * t814 + t446 * t863;
t754 = t360 * qJD(1);
t361 = t522 * pkin(1) * t534 - t444 * t863;
t753 = t361 * qJD(1);
t751 = t441 * qJD(3);
t456 = (-t531 ^ 2 + t527) * t522;
t750 = t456 * qJD(1);
t749 = t528 * qJD(3);
t748 = t529 * qJD(3);
t747 = t529 * qJD(5);
t746 = t530 * qJD(2);
t745 = t530 * qJD(4);
t743 = t532 * qJD(5);
t455 = t468 / 0.2e1;
t741 = t455 + t515 / 0.2e1;
t740 = pkin(4) * t837;
t736 = t882 / 0.2e1;
t735 = t874 / 0.2e1;
t734 = pkin(4) * t884;
t733 = t869 / 0.2e1;
t732 = t98 / 0.2e1 + t913;
t731 = t469 * t806;
t730 = t529 * t802;
t725 = t529 * t797;
t724 = t532 * t797;
t720 = t528 * t782;
t719 = t528 * t780;
t718 = t534 * t749;
t717 = t529 * t745;
t716 = t529 * t773;
t715 = t532 * t745;
t714 = t532 * t773;
t235 = t352 * t752;
t712 = t522 * t782;
t711 = t531 * t781;
t710 = t529 * t521;
t709 = t529 * t744;
t708 = t533 * t747;
t707 = t530 * t778;
t706 = t530 * t744;
t705 = t532 * t746;
t704 = t530 * t743;
t703 = t532 * t742;
t702 = t533 * t743;
t701 = t854 / 0.2e1;
t700 = -t853 / 0.2e1;
t699 = t217 * t529 / 0.2e1;
t698 = t350 * t888;
t696 = t352 * t887;
t695 = t352 * t897;
t694 = t469 * t889;
t693 = t823 / 0.2e1;
t691 = -t813 / 0.2e1;
t690 = t813 / 0.2e1;
t687 = -t254 / 0.2e1;
t684 = t807 / 0.2e1;
t683 = -t806 / 0.2e1;
t682 = t806 / 0.2e1;
t681 = -t258 / 0.2e1;
t679 = -t801 / 0.2e1;
t676 = t799 / 0.2e1;
t675 = -t234 / 0.2e1 - t282 / 0.2e1;
t673 = t355 / 0.2e1 + t909;
t672 = -t426 / 0.2e1 + t427 / 0.2e1;
t671 = t890 - t525 / 0.2e1;
t669 = t863 * qJD(1);
t664 = pkin(4) * t682;
t663 = t780 * t814;
t662 = t531 * t712;
t661 = t533 * t720;
t660 = t529 * t705;
t659 = t529 * t706;
t658 = t161 * t682;
t657 = t530 * t689;
t656 = t530 * t688;
t655 = t528 * t677;
t652 = -pkin(4) * t726 + t820 / 0.2e1;
t650 = t528 * t669;
t649 = t863 * t781;
t644 = 0.2e1 * t660;
t643 = t283 / 0.2e1 - t810 / 0.2e1;
t642 = -t819 / 0.2e1 + t910;
t641 = t804 / 0.2e1 + t892;
t478 = t530 * t703;
t639 = t533 * t748 + t478;
t638 = t477 - t717;
t637 = -t478 + t715;
t636 = pkin(9) * t657;
t635 = -qJD(3) + t720;
t535 = t109 * t909 + t129 * t904 + t161 * t470 / 0.2e1 + t190 * t896 + t77 * t908 + t99 * t903;
t549 = t117 * t894 + t130 * t892 + t188 * t891;
t1 = -t535 + t549;
t100 = t332 * t341 + t356 * t359 + t469 * t470;
t634 = -t1 * qJD(1) + t100 * qJD(2);
t101 = (t332 * t533 + t341 * t530) * t532 + (t359 * t530 + t836) * t529;
t600 = t356 * t900 + t99 * t884;
t552 = t129 * t888 + t600;
t554 = t109 * t888 + t332 * t900 + t77 * t884;
t587 = t399 * t892 + t400 * t894;
t590 = t341 * t906 - t359 * t350 / 0.2e1;
t4 = (-t130 / 0.2e1 + t554) * t532 + (t911 + t552) * t529 + t587 + t590;
t633 = t4 * qJD(1) - t101 * qJD(2);
t104 = pkin(4) * t731 - t356 * t920;
t5 = t673 * t99 + t732 * t356 + (t161 * t683 + t695 + t911) * pkin(4);
t632 = -t5 * qJD(1) + t104 * qJD(2);
t11 = (t876 / 0.2e1 + t732) * t532;
t69 = (t734 + t673) * t532;
t631 = t11 * qJD(1) + t69 * qJD(2);
t630 = t669 + qJD(2);
t541 = t350 * t673 + t732 * t808;
t10 = t877 / 0.2e1 + t541;
t102 = t920 * t808;
t629 = -t10 * qJD(1) - t102 * qJD(2);
t539 = t109 * t883 + t341 * t901 + t350 * t895 + t443 * t910;
t592 = -t851 / 0.2e1 + t399 * t891;
t13 = (t479 * t689 - t77 / 0.2e1) * t530 + (t161 * t884 + t190 * t888 + t441 * t896) * t529 + t539 + t592;
t133 = -t332 * t530 + t341 * t533 - t469 * t807 - t470 * t808;
t628 = -t13 * qJD(1) - t133 * qJD(2);
t138 = (t359 + t823) * t533 + (-t356 + t822) * t530;
t589 = t352 * t895 + t356 * t899;
t593 = t852 / 0.2e1 + t400 * t891;
t17 = (-t129 / 0.2e1 + t700) * t533 + (t903 + t693) * t441 + (t480 * t689 + t99 / 0.2e1 - t849 / 0.2e1) * t530 + t589 + t593;
t627 = -t17 * qJD(1) + t138 * qJD(2);
t181 = (t356 * t529 + t841) * t530;
t563 = t878 / 0.2e1 + pkin(3) * t691 - t788;
t591 = t350 * t904 + t352 * t909;
t30 = (t864 / 0.2e1 + t865 / 0.2e1) * t530 + t563 + t591;
t626 = -t30 * qJD(1) - t181 * qJD(2);
t203 = -pkin(4) * t730 - t731 - t836;
t551 = t695 + t600;
t570 = -t831 / 0.2e1 + t643;
t586 = -t803 / 0.2e1 + t697;
t28 = (t700 + (t586 + t812) * pkin(4)) * t530 + t551 + t570;
t625 = -t28 * qJD(1) - t203 * qJD(2);
t401 = -t457 + t739;
t212 = t401 * t530 + (-t515 + t468) * t533;
t571 = t401 * t899 + t785 * t901;
t605 = -pkin(3) * t399 / 0.2e1 - t843 / 0.2e1;
t617 = pkin(9) * t689 - t237 / 0.2e1;
t23 = (t261 - t881 / 0.2e1 + t922) * t533 + (t119 / 0.2e1 + (pkin(8) * t900 + t617) * t529) * t530 + t571 + t605;
t624 = -t23 * qJD(1) - t212 * qJD(2);
t213 = t467 * t533 + (-t402 + t737) * t530;
t603 = pkin(8) * t906 - t846 / 0.2e1;
t604 = pkin(3) * t902 + t844 / 0.2e1;
t24 = t441 * t898 + t829 / 0.2e1 + (t603 + t790) * t533 + (t120 / 0.2e1 + t617 * t532) * t530 + t604;
t623 = -t24 * qJD(1) + t213 * qJD(2);
t214 = t525 * t524 * pkin(4) - t355 * t533 - t469 * t808;
t553 = t350 * t896 + t355 * t901 + t98 * t883;
t585 = t833 / 0.2e1 + t675;
t32 = (t701 - t740) * t530 + t553 + t585;
t622 = -t32 * qJD(1) + t214 * qJD(2);
t338 = -t401 * t533 - t524 * t880;
t596 = t119 * t883 + t401 * t901;
t47 = (t881 / 0.2e1 + t699) * t530 + t596 + t675;
t621 = qJD(1) * t47 - qJD(2) * t338;
t339 = -pkin(8) * t802 - t402 * t533;
t595 = t120 * t884 + t402 * t900;
t46 = t530 * t603 + t595 + t643;
t620 = qJD(1) * t46 + qJD(2) * t339;
t125 = -t724 + (-t824 / 0.2e1 + t606) * t530;
t466 = t532 * t526 - t802;
t616 = -t125 * qJD(1) - t466 * qJD(2);
t126 = t725 + (t655 + t825 / 0.2e1 + t907) * t530;
t465 = t504 * t529;
t615 = -t126 * qJD(1) + t465 * qJD(2);
t186 = (-t809 - t837) * t530;
t463 = t502 * t524;
t614 = qJD(1) * t186 - qJD(2) * t463;
t226 = -t443 ^ 2 + t914;
t613 = qJD(1) * t226 + qJD(2) * t201;
t612 = t201 * qJD(1) + t504 * qJD(2);
t611 = t742 - t752;
t610 = qJD(1) * t443 + t746;
t192 = t333 - t803;
t609 = t192 * qJD(1) + t502 * qJD(3);
t608 = t350 * qJD(1) - t744;
t607 = t868 / 0.2e1 - t873 / 0.2e1;
t431 = t533 * t690 + t670 / 0.2e1;
t602 = t431 * qJD(1) + t746 / 0.2e1;
t601 = t528 * t648;
t542 = pkin(2) * t899 + t419 * t888 + t689 * t869;
t147 = t542 + t788;
t599 = pkin(2) * t746 + t147 * qJD(1);
t543 = pkin(2) * t901 - t827 / 0.2e1 + pkin(8) * t656;
t148 = t543 + t672;
t598 = pkin(2) * t742 + t148 * qJD(1);
t597 = (t904 + t821 / 0.2e1) * t532;
t594 = -t854 / 0.2e1 + t519 * t906;
t191 = -t809 / 0.2e1 + t837 / 0.2e1;
t583 = -t258 * qJD(1) + t703;
t582 = t705 + t755;
t581 = qJD(2) * t170 + qJD(3) * t254 + t235;
t580 = qJD(2) * t169 + t235 + t777;
t579 = pkin(3) * t906 + pkin(9) * t681;
t577 = -t745 / 0.2e1 + t665;
t576 = t607 * t529;
t575 = t607 * t532;
t574 = t519 * t683 + t694;
t225 = pkin(4) * t817;
t54 = t673 * t480 + (t908 + t574) * pkin(4);
t7 = t732 * t480 + (t912 + t594) * pkin(4);
t573 = -t7 * qJD(1) - t54 * qJD(2) + t225 * qJD(3);
t124 = -t869 / 0.2e1 + t597 + (t734 + t642) * t529;
t363 = t479 * t529 + t480 * t532;
t540 = t798 / 0.2e1 + t572 / 0.2e1;
t550 = t350 * t892 + t352 * t894 + t99 * t886;
t37 = (-t876 / 0.2e1 + t913) * t529 + t540 + t550;
t569 = -t37 * qJD(1) + t124 * qJD(2) + t363 * qJD(3);
t151 = t694 - t641 * t533 + (-t816 / 0.2e1 + (0.1e1 - t671) * pkin(4)) * t530 + t741;
t437 = t529 * t870 - t817;
t44 = t641 * t441 + (t443 + t191) * pkin(4) + t594 + t791;
t568 = -t44 * qJD(1) - t151 * qJD(2) - t437 * qJD(3);
t454 = -t467 / 0.2e1;
t167 = t454 + (t871 / 0.2e1 + t897) * t532 + (t519 * t887 + t676) * t529 + t652;
t450 = t523 * pkin(4) + t816;
t548 = -pkin(4) * t333 + t350 * t891 + t700;
t52 = (qJ(5) * t889 + t894) * t441 + t548 + t790;
t567 = -t52 * qJD(1) - t167 * qJD(2) + t450 * qJD(3);
t132 = (-t809 + t837) * t530;
t136 = -t803 + 0.2e1 * t697;
t345 = t350 ^ 2;
t158 = t345 - t915;
t566 = qJD(1) * t158 - qJD(2) * t132 + qJD(3) * t136;
t202 = t345 + t915;
t565 = qJD(1) * t202 - qJD(2) * t186 + qJD(3) * t192;
t348 = t898 - t576;
t547 = pkin(3) * t907 + t846 / 0.2e1 + pkin(9) * t254 / 0.2e1;
t58 = t547 - t790;
t562 = pkin(3) * t744 - t58 * qJD(1) - t348 * qJD(2);
t349 = -t468 / 0.2e1 + t575;
t60 = -t274 / 0.2e1 - t922 + t579;
t561 = pkin(3) * t748 - t60 * qJD(1) - t349 * qJD(2);
t185 = t586 * t530;
t559 = -qJD(2) * t185 - qJD(3) * t191 + t350 * t755;
t449 = t671 * t530;
t558 = t191 * qJD(1) - t449 * qJD(2) + t709;
t557 = t529 * t746 + t608;
t556 = t582 + t748;
t555 = qJD(4) * t431 + t584;
t464 = t503 * t524;
t546 = t132 * qJD(1) + t464 * qJD(2) + 0.2e1 * t659;
t545 = -t136 * qJD(1) + t644 - t918;
t544 = -t185 * qJD(1) + qJD(2) * t730 + t449 * qJD(3);
t517 = t779 / 0.2e1;
t484 = qJD(2) * t690;
t473 = -0.2e1 * t530 * t710;
t440 = t449 * qJD(4);
t421 = (t712 - t749 / 0.2e1) * t531;
t342 = qJD(2) * t656 + t431 * qJD(3);
t302 = t515 + t455 + t575;
t301 = t454 - t576 + t738;
t238 = t556 * pkin(4);
t196 = t201 * qJD(3);
t189 = t192 * qJD(5);
t187 = t191 * qJD(4);
t182 = t186 * qJD(5);
t180 = t185 * qJD(4);
t168 = pkin(8) * t682 + t519 * t686 + t529 * t676 + t454 - t652 + t693;
t162 = t830 / 0.2e1 + t564;
t152 = t818 / 0.2e1 + t872 * t890 + pkin(4) * t692 + t872 + qJ(5) * t679 - t574 + t741;
t150 = -t542 + t788;
t149 = -t543 + t672;
t135 = t136 * qJD(4);
t131 = t132 * qJD(4);
t128 = t443 * t682 + t529 * t656 + t696 + t724;
t127 = t443 * t686 + t530 * t655 + t698 - t725;
t123 = pkin(4) * t684 + t529 * t642 + t597 + t733;
t96 = t399 * t889 + t400 * t885 + t586 * t533 + t369;
t70 = t355 * t886 - t841 / 0.2e1 + pkin(4) * t679;
t61 = t699 + t579 + t791;
t59 = t547 + t790;
t55 = pkin(4) * t908 + t469 * t735 + t519 * t664 + t893 * t920;
t53 = qJ(5) * t687 + t479 * t900 - t548 + t790;
t49 = pkin(8) * t696 + t217 * t682 - t595 + t643;
t48 = pkin(8) * t698 + t217 * t686 - t596 + t675;
t45 = t441 * t893 + t350 * t735 - t740 / 0.2e1 + t875 + qJ(5) * t680 - t594 + t791;
t38 = pkin(4) * t687 + t77 * t889 + t540 - t550;
t33 = t161 * t686 + t530 * t740 - t553 + t585;
t31 = t683 * t77 + t686 * t99 + t563 - t591;
t29 = t350 * t664 + t685 * t882 - t551 + t570 + t658 + t666;
t26 = t651 * t901 - t829 / 0.2e1 + t789 * t884 + t120 * t888 + t352 * t733 + t681 * t871 + t237 * t682 + t217 * t678 + t532 * t636 + t604;
t25 = pkin(8) * t921 + t119 * t888 + t217 * t684 + t237 * t685 + t529 * t636 + t668 * t884 - t571 + t605;
t18 = t129 * t883 + t161 * t678 + t190 * t682 + t359 * t901 + t469 * t681 + t480 * t657 + t99 * t888 - t589 + t593;
t14 = t161 * t684 + t190 * t685 + t469 * t687 + t479 * t657 + t77 * t887 - t539 + t592;
t12 = t98 * t886 - t865 / 0.2e1 + pkin(4) * t680;
t9 = -t877 / 0.2e1 + t541;
t8 = t519 * t736 + t866 * t893 + (t701 + t912) * pkin(4);
t6 = t469 * t736 + (-t355 / 0.2e1 + t910) * t99 - t866 * t356 / 0.2e1 + (t658 + t911) * pkin(4);
t3 = t117 * t889 + t130 * t885 + t552 * t529 + t554 * t532 - t587 + t590;
t2 = t535 + t549;
t62 = [0, 0, 0, t663, t456 * qJD(2), t534 * t649, -t531 * t649, 0, -t360 * qJD(2), -t361 * qJD(2), (t533 * t719 - t751) * t443, qJD(2) * t223 + qJD(3) * t226, -t335 * qJD(2) + t441 * t718, t334 * qJD(2) + t443 * t718, -t663, -qJD(2) * t82 - qJD(3) * t153, qJD(2) * t83 + qJD(3) * t154, (qJD(2) * t400 - t441 * t744 - t343) * t352, qJD(2) * t134 + qJD(3) * t139 + qJD(4) * t158, qJD(2) * t176 + qJD(3) * t165 - t343 * t441, qJD(2) * t175 - qJD(3) * t164 - t441 * t777, (t443 * qJD(3) + t530 * t719) * t441, qJD(2) * t41 + qJD(3) * t39 + qJD(4) * t57, qJD(2) * t42 + qJD(3) * t40 + qJD(4) * t56, qJD(2) * t34 + qJD(3) * t27 + qJD(4) * t51 - t352 * t771, qJD(2) * t36 + qJD(3) * t35 + qJD(4) * t50 + t350 * t771, qJD(2) * t22 + qJD(3) * t21 + qJD(4) * t15 + qJD(5) * t202, qJD(2) * t20 + qJD(3) * t19 + qJD(4) * t16 + qJD(5) * t43; 0, 0, 0, t662, t750, t630 * t812, -t630 * t813, 0, -qJD(2) * t446 - t754, qJD(2) * t444 - t753, t528 * t610 * t795 + t917, t504 * t719 + t196 + t759, t530 * t711 - t756, t533 * t711 + t757, -t421, -t794 + (-t446 * t533 + t530 * t601) * qJD(2) + t150 * qJD(3), t793 + (t446 * t530 + t533 * t601) * qJD(2) + t149 * qJD(3), t162 * qJD(3) + t400 * t582 + t180, t770 + t96 * qJD(3) - t131 + (-t399 * t532 - t830) * t746, t760 + (-t400 * t533 + t524 * t727) * qJD(2) + t128 * qJD(3) + t537 * qJD(4), t761 + (t399 * t533 - t524 * t729) * qJD(2) + t127 * qJD(3) - t169 * qJD(4), -t917 + (qJD(4) / 0.2e1 - t611) * t728, -t159 * t742 + t828 + t25 * qJD(3) + t49 * qJD(4) + (pkin(8) * t399 - t401 * t812 + t844) * t746, t160 * t742 + t826 + t26 * qJD(3) + t48 * qJD(4) + (pkin(8) * t400 - t402 * t812 + t843) * t746, t840 + (-t117 * t533 + t469 * t399 + (t332 * t812 + t852) * t530) * qJD(2) + t14 * qJD(3) + t29 * qJD(4) - t170 * qJD(5), t835 + (t130 * t533 + t469 * t400 + (-t356 * t812 + t851) * t530) * qJD(2) + t18 * qJD(3) + t33 * qJD(4) + t536 * qJD(5), t845 + (-t332 * t400 - t356 * t399 + (-t117 * t532 - t130 * t529) * t530) * qJD(2) + t3 * qJD(3) + t9 * qJD(4) - t182, t848 + (t117 * t332 + t130 * t356 + t188 * t469) * qJD(2) + t2 * qJD(3) + t6 * qJD(4) + t31 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t584, t613, t635 * t441, t635 * t443, t484, qJD(2) * t150 - qJD(3) * t237 - t784, qJD(2) * t149 + qJD(3) * t236 + t783, t162 * qJD(2) + t187 + (-t748 - t755) * t258, t96 * qJD(2) - t503 * t751 + t135 + t769, t128 * qJD(2) + t443 * t748 + t766, t127 * qJD(2) + t443 * t744 - t767, t555, t834 + t25 * qJD(2) + (-t237 * t532 + t529 * t647) * qJD(3) + t61 * qJD(4), t832 + t26 * qJD(2) + (t237 * t529 + t532 * t647) * qJD(3) + t59 * qJD(4), t842 + t14 * qJD(2) + (-t254 * t519 - t479 * t443 - t849) * qJD(3) + t45 * qJD(4) - t254 * qJD(5), t839 + t18 * qJD(2) + (t190 * t529 - t258 * t519 - t480 * t443) * qJD(3) + t53 * qJD(4) - t441 * t743, t847 + t3 * qJD(2) + (-t109 * t529 + t129 * t532 + (-t479 * t532 + t480 * t529) * t441) * qJD(3) + t12 * qJD(4) + t189, t850 + t2 * qJD(2) + (-t109 * t479 + t129 * t480 + t190 * t519) * qJD(3) + t8 * qJD(4) + t38 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t559, t566, t916, -t580, t342, qJD(2) * t49 + qJD(3) * t61 - qJD(4) * t120 + t856, qJD(2) * t48 + qJD(3) * t59 + qJD(4) * t119 + t857, qJD(2) * t29 + qJD(3) * t45 - t855 + t858, qJD(2) * t33 + qJD(3) * t53 + qJD(4) * t98 + t859, pkin(4) * t343 + qJD(2) * t9 + qJD(3) * t12 + t862, -pkin(4) * t855 + qJD(2) * t6 + qJD(3) * t8 + t861; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t581, t441 * t608 + t919, t565, qJD(2) * t31 + qJD(3) * t38 + t860; 0, 0, 0, -t662, -t750, -t534 * t650, t531 * t650, 0, t754, t753, -t443 * t661 + t917, t196 - t759, -t533 * t718 + t756, t530 * t718 - t757, t421, -qJD(3) * t147 + t794, -qJD(3) * t148 - t793, qJD(3) * t163 - t400 * t755 + t180, qJD(3) * t97 - t131 - t770, -qJD(3) * t125 - qJD(4) * t536 - t760, -qJD(3) * t126 - qJD(4) * t170 - t761, -t917 + (-t752 - qJD(4) / 0.2e1) * t728, -qJD(3) * t23 - qJD(4) * t46 - t828, -qJD(3) * t24 - qJD(4) * t47 - t826, -qJD(3) * t13 - qJD(4) * t28 - qJD(5) * t169 - t840, -qJD(3) * t17 - qJD(4) * t32 - qJD(5) * t537 - t835, qJD(3) * t4 + qJD(4) * t10 - t182 - t845, -qJD(3) * t1 - qJD(4) * t5 - qJD(5) * t30 - t848; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, t504 * qJD(3), 0, 0, 0, -pkin(2) * t779, -pkin(2) * t778, -t524 * t710 + t525 * t707, -t464 * qJD(4) - 0.2e1 * t533 * t659, -t466 * qJD(3) + t530 * t716, t465 * qJD(3) + t530 * t714, -t707, -qJD(3) * t212 - qJD(4) * t339, qJD(3) * t213 + qJD(4) * t338, -t133 * qJD(3) - t203 * qJD(4) + t530 * t702, t138 * qJD(3) + t214 * qJD(4) - t530 * t708, -qJD(3) * t101 + qJD(4) * t102 + qJD(5) * t463, qJD(3) * t100 + qJD(4) * t104 - qJD(5) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t665, t612, -t635 * t533, t635 * t530, qJD(1) * t691, -pkin(8) * t778 - t599, pkin(8) * t779 - t598, t768 - t440 + (t525 * t746 + t709) * t533, t792 + t473 + (-0.2e1 * t660 + t918) * t533, t530 * t748 + t616, t615 + t706, -t577, (t529 * t646 - t737) * qJD(3) + t302 * qJD(4) + t624, (t532 * t646 + t739) * qJD(3) + t301 * qJD(4) + t623, (t519 * t807 - t821 - t822) * qJD(3) + t152 * qJD(4) + t708 + t628, (t470 * t529 + t519 * t801 - t819) * qJD(3) + t168 * qJD(4) + t702 + t627, ((t359 + t820) * t532 + (-t341 - t818) * t529) * qJD(3) + t70 * qJD(4) + t633, (-t341 * t479 + t359 * t480 + t470 * t519) * qJD(3) + t55 * qJD(4) + t123 * qJD(5) + t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t544, -t546, t638 - t763, -t637 - t764, qJD(1) * t657 + t517, qJD(3) * t302 - qJD(4) * t402 - t620, qJD(3) * t301 + qJD(4) * t401 - t621, qJD(3) * t152 + t625 - t776, qJD(3) * t168 + qJD(4) * t355 + t622, pkin(4) * t717 + t70 * qJD(3) - t629, -pkin(4) * t776 + qJD(3) * t55 + t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t639 - t765, -t762 + t786, -t614, qJD(3) * t123 + t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t584, -t613, t611 * t812, -t610 * t812, t484, qJD(2) * t147 + t784, qJD(2) * t148 - t783, -t163 * qJD(2) + t235 * t532 + t187, -qJD(2) * t97 + t135 - t769, qJD(2) * t125 + qJD(4) * t258 - t766, qJD(2) * t126 - qJD(4) * t254 + t767, -t555, qJD(2) * t23 + qJD(4) * t60 - t834, qJD(2) * t24 + qJD(4) * t58 - t832, qJD(2) * t13 - qJD(4) * t44 - t842, qJD(2) * t17 - qJD(4) * t52 - t839, -qJD(2) * t4 - qJD(4) * t11 + t189 - t847, qJD(2) * t1 - qJD(4) * t7 - qJD(5) * t37 - t850; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t665, -t612, t661, -t530 * t720, qJD(1) * t690, t599, t598, -t507 * t525 - t440 - t768, t533 * t644 + t473 - t792, -t616 - t714, -t615 + t716, t577, qJD(4) * t349 - t624, qJD(4) * t348 - t623, -qJD(4) * t151 - t628, -qJD(4) * t167 - t627, -qJD(4) * t69 - t633, -qJD(4) * t54 + qJD(5) * t124 - t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t710, t503 * qJD(4), 0, 0, 0, -pkin(3) * t774, -pkin(3) * t521, -t437 * qJD(4), t450 * qJD(4), t502 * qJD(5), qJD(4) * t225 + qJD(5) * t363; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t558, -t545, t521 - t583, -t758 + (-qJD(4) + t742) * t529, -t602, -pkin(9) * t521 - t561, pkin(9) * t774 - t562, t568 - t775, qJD(4) * t479 + t567, -pkin(4) * t521 - t631, -pkin(4) * t775 + t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t609, t569; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t559, -t566, -qJD(3) * t258 + t713 + t919, t581, t342, qJD(2) * t46 - qJD(3) * t60 - t856, qJD(2) * t47 - qJD(3) * t58 - t857, qJD(2) * t28 + qJD(3) * t44 - t772 - t858, qJD(2) * t32 + qJD(3) * t52 + qJD(5) * t350 - t859, -qJD(2) * t10 + qJD(3) * t11 - t862, -pkin(4) * t772 + qJD(2) * t5 + qJD(3) * t7 - t861; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t544, t546, t763 + t786, -t639 + t764, qJD(1) * t656 + t517, -qJD(3) * t349 + t620, -qJD(3) * t348 + t621, t151 * qJD(3) - t625 - t704, t167 * qJD(3) + t530 * t747 - t622, qJD(3) * t69 + t629, -pkin(4) * t704 + t54 * qJD(3) - t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t558, t545, t583, -t529 * t742 + t758, t602, t561, t562, -t568 - t747, -t567 - t743, t631, -pkin(4) * t747 - t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t556, t557, 0, -t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t580, t916, -t565, pkin(4) * t777 + qJD(2) * t30 + qJD(3) * t37 - t860; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t637 + t765, t638 + t762, t614, pkin(4) * t715 - t124 * qJD(3) - t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t774, t521, -t609, pkin(4) * t774 - t569; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t556, -t557, 0, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t62;
