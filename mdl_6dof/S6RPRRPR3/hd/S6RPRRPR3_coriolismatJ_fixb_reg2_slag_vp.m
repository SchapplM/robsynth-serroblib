% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:22
% EndTime: 2019-03-09 05:07:49
% DurationCPUTime: 20.71s
% Computational Cost: add. (13825->782), mult. (26604->1004), div. (0->0), fcn. (25588->8), ass. (0->573)
t796 = qJD(4) - qJD(6);
t621 = sin(qJ(6));
t624 = cos(qJ(4));
t597 = -cos(pkin(10)) * pkin(1) - pkin(2);
t625 = cos(qJ(3));
t623 = sin(qJ(3));
t932 = t623 * pkin(8);
t711 = -t625 * pkin(3) - t932;
t646 = t711 + t597;
t482 = t624 * t646;
t881 = t624 * t623;
t794 = pkin(9) * t881;
t714 = -t482 - t794;
t680 = t621 * t714;
t596 = sin(pkin(10)) * pkin(1) + pkin(7);
t622 = sin(qJ(4));
t732 = t596 * t622 + pkin(4);
t723 = pkin(5) + t732;
t682 = t723 * t625;
t878 = t625 * t596;
t548 = t624 * t878;
t364 = t622 * t646 + t548;
t887 = t622 * t623;
t688 = pkin(9) * t887 + t364;
t879 = t625 * qJ(5);
t248 = t688 - t879;
t935 = cos(qJ(6));
t782 = t935 * t248;
t139 = t621 * t682 + t680 + t782;
t255 = t935 * t688;
t788 = t622 * t878;
t363 = -t482 + t788;
t286 = -t363 + t794;
t895 = t621 * t286;
t155 = -t255 + t895;
t986 = t139 + t155;
t594 = pkin(4) * t887;
t884 = t624 * qJ(5);
t730 = -t596 + t884;
t348 = -t594 + (-pkin(5) * t622 + t730) * t623;
t964 = pkin(8) - pkin(9);
t731 = t964 * t622;
t773 = t935 * t624;
t405 = -t621 * t731 - t964 * t773;
t536 = t621 * t622 + t773;
t478 = t536 * t623;
t937 = -t625 / 0.2e1;
t775 = t935 * t622;
t892 = t621 * t624;
t541 = t775 - t892;
t953 = t541 / 0.2e1;
t890 = t622 * qJ(5);
t965 = pkin(4) + pkin(5);
t969 = -t624 * t965 - t890;
t517 = pkin(3) - t969;
t956 = t517 / 0.2e1;
t640 = t348 * t953 - t405 * t937 + t478 * t956;
t985 = t796 * t405;
t724 = t964 * t892;
t974 = t935 * t731;
t980 = -t974 + t724;
t960 = -t980 / 0.2e1;
t774 = t935 * t623;
t477 = t621 * t881 - t622 * t774;
t864 = t477 * t796;
t862 = t478 * t796;
t798 = t625 * qJD(1);
t770 = t478 * t798;
t984 = t862 - t770;
t771 = t477 * t798;
t983 = t864 - t771;
t952 = t974 / 0.2e1;
t975 = t796 * t541;
t982 = t536 * t975;
t619 = t624 ^ 2;
t607 = t619 * t625;
t617 = t622 ^ 2;
t898 = t617 * t625;
t981 = t607 + t898;
t954 = t536 / 0.2e1;
t669 = t477 * t953 + t478 * t954;
t844 = qJD(3) * t541;
t979 = qJD(1) * t669 + t536 * t844;
t847 = qJD(1) * t478;
t978 = qJD(3) * t669 + t477 * t847;
t869 = t796 * t669;
t136 = t536 * t477 - t478 * t541;
t183 = t477 ^ 2 - t478 ^ 2;
t686 = qJD(1) * t183 + qJD(3) * t136;
t260 = t536 ^ 2 - t541 ^ 2;
t685 = qJD(1) * t136 + qJD(3) * t260;
t977 = t536 * t796;
t801 = t623 * qJD(3);
t755 = t624 * t801;
t620 = t625 ^ 2;
t618 = t623 ^ 2;
t888 = t622 * t618;
t540 = t620 * t622 - t888;
t806 = t540 * qJD(1);
t973 = t806 + t755;
t839 = qJD(4) * t625;
t765 = t622 * t839;
t972 = t806 - t765;
t931 = t625 * pkin(8);
t933 = t623 * pkin(3);
t568 = -t931 + t933;
t543 = t622 * t568;
t413 = -t596 * t881 + t543;
t912 = t413 * t623;
t547 = t596 * t887;
t883 = t624 * t568;
t412 = t547 + t883;
t913 = t412 * t623;
t919 = t364 * t625;
t920 = t363 * t625;
t60 = (t912 / 0.2e1 + t919 / 0.2e1) * t624 + (-t913 / 0.2e1 + t920 / 0.2e1) * t622 + (t618 / 0.2e1 - t620 / 0.2e1) * t596;
t605 = t623 * t625;
t420 = (-0.1e1 + t617 + t619) * t605;
t813 = t420 * qJD(2);
t971 = -t60 * qJD(1) - t813;
t315 = t364 - t879;
t316 = t625 * t732 - t482;
t934 = pkin(4) * t622;
t567 = -t884 + t934;
t417 = (t567 + t596) * t625;
t615 = t623 * pkin(4);
t380 = -t412 - t615;
t917 = t380 * t622;
t613 = t623 * qJ(5);
t377 = t413 + t613;
t918 = t377 * t624;
t938 = t624 / 0.2e1;
t942 = t622 / 0.2e1;
t416 = -t623 * t730 + t594;
t959 = t416 / 0.2e1;
t39 = (t315 * t938 - t417 / 0.2e1 + t316 * t942) * t625 + (t918 / 0.2e1 + t959 + t917 / 0.2e1) * t623;
t970 = t39 * qJD(1) + t813;
t968 = 0.2e1 * t622 * t881 * (qJD(4) + t798) + (-t607 + t898) * qJD(3);
t659 = t621 * t688;
t780 = t935 * t286;
t156 = t780 + t659;
t963 = -t156 / 0.2e1;
t962 = t364 / 0.2e1;
t961 = -t405 / 0.2e1;
t572 = t625 * t775;
t880 = t624 * t625;
t479 = t621 * t880 - t572;
t958 = t479 / 0.2e1;
t480 = t625 * t536;
t742 = t480 / 0.2e1;
t704 = t624 * pkin(4) + t890;
t513 = t704 * t623;
t957 = -t513 / 0.2e1;
t527 = -t622 * t965 + t884;
t955 = t527 / 0.2e1;
t951 = -t547 / 0.2e1;
t549 = qJ(5) * t621 + t935 * t965;
t950 = -t549 / 0.2e1;
t550 = t935 * qJ(5) - t621 * t965;
t949 = -t550 / 0.2e1;
t948 = t550 / 0.2e1;
t947 = t567 / 0.2e1;
t946 = -t568 / 0.2e1;
t945 = -t621 / 0.2e1;
t944 = t621 / 0.2e1;
t943 = -t622 / 0.2e1;
t941 = -t623 / 0.2e1;
t940 = t623 / 0.2e1;
t939 = -t624 / 0.2e1;
t936 = t625 / 0.2e1;
t739 = -t363 / 0.2e1 + t316 / 0.2e1;
t635 = t739 * t624 + (-t315 / 0.2e1 + t962) * t622;
t54 = t513 * t937 + t623 * t635;
t929 = qJD(1) * t54;
t446 = t969 * t623;
t636 = t682 + t714;
t195 = t935 * t636;
t897 = t621 * t248;
t138 = -t195 + t897;
t740 = t156 / 0.2e1 - t138 / 0.2e1;
t741 = t139 / 0.2e1 + t155 / 0.2e1;
t12 = t446 * t936 + t477 * t741 + t478 * t740;
t928 = t12 * qJD(1);
t256 = -t623 * pkin(5) - t547 - t615 + (-pkin(9) * t625 - t568) * t624;
t781 = t935 * t256;
t886 = t622 * t625;
t292 = pkin(9) * t886 + t377;
t894 = t621 * t292;
t149 = t781 - t894;
t779 = t935 * t292;
t896 = t621 * t256;
t150 = t779 + t896;
t14 = t138 * t480 - t139 * t479 - t149 * t478 - t150 * t477;
t927 = t14 * qJD(1);
t17 = t986 * t478 + (t138 - t156) * t477;
t926 = t17 * qJD(1);
t349 = (-t596 + t527) * t625;
t27 = t138 * t623 + t149 * t625 + t348 * t479 + t349 * t477;
t925 = t27 * qJD(1);
t28 = -t139 * t623 + t150 * t625 - t348 * t480 - t349 * t478;
t924 = t28 * qJD(1);
t923 = t315 * t625;
t922 = t348 * t477;
t921 = t348 * t478;
t915 = t980 * t621;
t914 = t980 * t625;
t911 = t416 * t622;
t910 = t416 * t624;
t709 = t255 / 0.2e1 - t895 / 0.2e1;
t42 = -t782 / 0.2e1 - t680 / 0.2e1 + (t723 * t945 + t948) * t625 + t709;
t909 = t42 * qJD(1);
t908 = t477 * t621;
t907 = t479 * t541;
t906 = t513 * t622;
t905 = t541 * t625;
t904 = t549 * t625;
t55 = -t155 * t625 + t446 * t477 - t921;
t903 = t55 * qJD(1);
t553 = -pkin(3) - t704;
t902 = t553 * t622;
t56 = t156 * t625 - t446 * t478 - t922;
t901 = t56 * qJD(1);
t744 = -t879 / 0.2e1;
t58 = t364 * t943 + (t744 + t315 / 0.2e1) * t622 + (pkin(4) * t937 - t739) * t624;
t900 = t58 * qJD(1);
t893 = t621 * t541;
t891 = t621 * t625;
t889 = t622 * t478;
t885 = t623 * t553;
t882 = t624 * t618;
t63 = ((t315 - t364) * t624 + (t316 - t363) * t622) * t623;
t877 = t63 * qJD(1);
t65 = -t316 * t880 - t380 * t881 + (t377 * t623 + t923) * t622;
t876 = t65 * qJD(1);
t68 = -t138 * t625 + t922;
t875 = t68 * qJD(1);
t69 = -t139 * t625 + t921;
t874 = t69 * qJD(1);
t70 = (t913 - t920) * t624 + (t912 + t919) * t622;
t873 = t70 * qJD(1);
t693 = t416 * t625 + t417 * t623;
t72 = -t315 * t623 + t377 * t625 + t624 * t693;
t872 = t72 * qJD(1);
t73 = -t316 * t623 + t380 * t625 + t622 * t693;
t871 = t73 * qJD(1);
t870 = t796 * t136;
t745 = t880 / 0.2e1;
t750 = t886 / 0.2e1;
t861 = t621 * t750 + t935 * t745;
t304 = t742 + t861;
t743 = -t480 / 0.2e1;
t783 = t625 * t935;
t716 = -t783 / 0.2e1;
t751 = -t886 / 0.2e1;
t860 = t621 * t751 + t624 * t716;
t305 = t743 + t860;
t867 = -t305 * qJD(4) - t304 * qJD(6);
t303 = t743 + t861;
t306 = t742 + t860;
t866 = -t303 * qJD(4) - t306 * qJD(6);
t746 = -t880 / 0.2e1;
t859 = t572 / 0.2e1 + t621 * t746;
t858 = -t572 / 0.2e1 + t621 * t745;
t857 = t981 * pkin(8);
t147 = t919 + (t906 + t910) * t623;
t856 = qJD(1) * t147;
t776 = t935 * t536;
t715 = t776 / 0.2e1;
t754 = -t905 / 0.2e1;
t777 = t935 * t480;
t154 = t625 * t715 + t777 / 0.2e1 + (t754 + t958) * t621;
t855 = qJD(1) * t154;
t748 = -t881 / 0.2e1;
t784 = t478 * t935;
t181 = (t784 / 0.2e1 + t908 / 0.2e1 + t748) * t625;
t854 = qJD(1) * t181;
t778 = t935 * t477;
t200 = t478 * t891 - t625 * t778;
t853 = qJD(1) * t200;
t201 = -t596 * t888 - t920;
t852 = qJD(1) * t201;
t202 = -t596 * t882 - t919;
t851 = qJD(1) * t202;
t238 = t623 * t477 - t625 * t479;
t850 = qJD(1) * t238;
t239 = -t478 * t623 + t480 * t625;
t849 = qJD(1) * t239;
t848 = qJD(1) * t477;
t785 = t405 * t935;
t373 = -t785 / 0.2e1;
t717 = t785 / 0.2e1;
t180 = t373 + t717;
t846 = qJD(3) * t180;
t845 = qJD(3) * t536;
t843 = qJD(3) * t622;
t842 = qJD(3) * t624;
t841 = qJD(4) * t622;
t840 = qJD(4) * t624;
t838 = qJD(5) * t180;
t837 = qJD(5) * t621;
t836 = qJD(5) * t622;
t835 = qJD(5) * t625;
t834 = qJD(6) * t517;
t148 = -t416 * t887 + t513 * t881 - t920;
t833 = t148 * qJD(1);
t158 = t363 * t623 + (t412 - 0.2e1 * t547) * t625;
t832 = t158 * qJD(1);
t159 = t413 * t625 + (-t364 + 0.2e1 * t548) * t623;
t831 = t159 * qJD(1);
t165 = -t480 * t477 - t478 * t479;
t830 = t165 * qJD(1);
t169 = t416 * t881 + t923;
t829 = t169 * qJD(1);
t299 = t754 + t858;
t828 = t299 * qJD(1);
t827 = t299 * qJD(3);
t753 = t905 / 0.2e1;
t300 = t753 + t858;
t826 = t300 * qJD(3);
t301 = t754 + t859;
t825 = t301 * qJD(3);
t302 = t753 + t859;
t824 = t302 * qJD(1);
t823 = t302 * qJD(3);
t822 = t303 * qJD(3);
t821 = t304 * qJD(1);
t820 = t304 * qJD(3);
t819 = t305 * qJD(1);
t818 = t305 * qJD(3);
t817 = t306 * qJD(3);
t816 = t363 * qJD(4);
t365 = t477 * t881 + t620 * t621;
t815 = t365 * qJD(1);
t366 = t478 * t881 + t620 * t935;
t814 = t366 * qJD(1);
t735 = t617 / 0.2e1 - t619 / 0.2e1;
t519 = t735 * t623;
t808 = t519 * qJD(4);
t542 = t620 * t624 - t882;
t805 = t542 * qJD(1);
t580 = t619 - t617;
t804 = t580 * qJD(4);
t581 = t620 - t618;
t803 = t581 * qJD(1);
t611 = t621 * qJD(4);
t802 = t623 * qJD(1);
t800 = t623 * qJD(4);
t799 = t624 * qJD(5);
t797 = t625 * qJD(3);
t795 = t951 - t615;
t793 = t935 / 0.2e1;
t792 = pkin(8) * t841;
t791 = pkin(8) * t840;
t790 = pkin(3) * t942;
t789 = t931 / 0.2e1;
t787 = t149 * t935;
t786 = t155 * t935;
t769 = t624 * t802;
t767 = t536 * t801;
t766 = t541 * t801;
t764 = t624 * t839;
t763 = t621 * t835;
t762 = t622 * t799;
t761 = t597 * t802;
t760 = t597 * t798;
t588 = t622 * t840;
t587 = t622 * t842;
t759 = t622 * t798;
t758 = t623 * t797;
t757 = t623 * t836;
t756 = t623 * t798;
t752 = t477 * t942;
t749 = t885 / 0.2e1;
t747 = t881 / 0.2e1;
t738 = t980 / 0.2e1 + t960;
t736 = t951 - t615 / 0.2e1;
t734 = qJD(6) * t935;
t612 = t935 * qJD(4);
t733 = t935 * qJD(5);
t487 = (-0.1e1 / 0.2e1 + t735) * t623;
t729 = qJD(1) * t487 - t587;
t418 = qJD(1) * t519 - t587;
t565 = t622 * qJD(1) * t882;
t395 = qJD(3) * t519 + t565;
t726 = t796 * t625;
t592 = -qJD(4) + t798;
t722 = t622 * t769;
t721 = t618 * t588;
t720 = t622 * t755;
t719 = t625 * t733;
t718 = t935 * t961;
t713 = t798 - qJD(4) / 0.2e1;
t710 = 0.2e1 * t720;
t708 = -t904 / 0.2e1 + t138 / 0.2e1;
t707 = t567 * t940 + t959;
t703 = -t553 * t625 + t932;
t18 = -t138 * t149 + t139 * t150 + t348 * t349;
t8 = t150 * t478 / 0.2e1 + t139 * t742 - t149 * t477 / 0.2e1 + t138 * t958 + t349 * t936 + t348 * t941;
t702 = t18 * qJD(1) + t8 * qJD(2);
t626 = -t477 * t738 - t536 * t740 + t541 * t741;
t667 = t479 * t949 + t549 * t742;
t5 = t626 - t667;
t701 = t5 * qJD(1);
t19 = t138 * t155 + t139 * t156 + t348 * t446;
t700 = t19 * qJD(1) + t12 * qJD(2);
t157 = t477 * t479 + t478 * t480 - t605;
t699 = t8 * qJD(1) + t157 * qJD(2);
t52 = t315 * t377 + t316 * t380 + t416 * t417;
t698 = t52 * qJD(1) + t39 * qJD(2);
t57 = -t315 * t363 + t316 * t364 + t416 * t513;
t697 = t57 * qJD(1) + t54 * qJD(2);
t78 = t596 ^ 2 * t605 - t363 * t412 + t364 * t413;
t696 = t78 * qJD(1) + t60 * qJD(2);
t695 = t917 + t918;
t694 = -t412 * t622 + t413 * t624;
t37 = t138 * t891 + t139 * t783 - t348 * t881;
t692 = -qJD(1) * t37 - qJD(2) * t181;
t388 = t553 * t624 + t567 * t622;
t523 = -t543 / 0.2e1;
t574 = pkin(8) * t750;
t93 = t574 - t906 / 0.2e1 - t910 / 0.2e1 - t613 + t523 + (t902 / 0.2e1 + (-t567 / 0.2e1 + t596 / 0.2e1) * t624) * t623;
t691 = -qJD(1) * t93 + qJD(3) * t388;
t389 = t567 * t624 - t902;
t676 = t749 + t789;
t653 = t957 + t676;
t673 = t707 * t622;
t95 = t673 + (t946 + t653) * t624 + t795;
t690 = -qJD(1) * t95 + qJD(3) * t389;
t392 = t952 - t974 / 0.2e1;
t632 = -t659 / 0.2e1 - t780 / 0.2e1;
t651 = t195 / 0.2e1 + t904 / 0.2e1 - t897 / 0.2e1;
t43 = t632 - t651;
t689 = -t43 * qJD(1) + t392 * qJD(3);
t687 = t592 * t623;
t681 = t789 - t933 / 0.2e1;
t534 = qJD(1) * t783 - t612;
t678 = -t934 / 0.2e1 + t884 / 0.2e1;
t677 = -t377 * qJ(5) / 0.2e1 + t380 * pkin(4) / 0.2e1;
t522 = t543 / 0.2e1;
t575 = pkin(8) * t751;
t382 = t623 * t790 + t522 + t575;
t675 = pkin(3) * t842 - qJD(1) * t382;
t383 = (t946 + t681) * t624;
t674 = pkin(3) * t843 - qJD(1) * t383;
t672 = t149 * t950 + t150 * t948;
t671 = t348 * t954 + t477 * t956;
t668 = t479 * t950 + t480 * t949;
t476 = t624 * t687;
t151 = t911 / 0.2e1 + (t946 + t676) * t624 + t736;
t666 = qJD(1) * t151 + t553 * t843;
t196 = -t889 / 0.2e1 + (t541 * t939 + t944) * t623;
t661 = -qJD(1) * t196 + t541 * t843;
t197 = t752 + (t536 * t938 + t793) * t623;
t660 = qJD(1) * t197 + t536 * t843;
t658 = t348 * t942 + t517 * t747;
t657 = t479 * t793 + t480 * t945;
t656 = -t896 / 0.2e1 - t779 / 0.2e1;
t655 = -t894 / 0.2e1 + t781 / 0.2e1;
t606 = t619 * t618;
t538 = t617 * t618 - t606;
t410 = -qJD(1) * t538 + t710;
t449 = -qJD(3) * t580 + 0.2e1 * t722;
t654 = t540 * qJD(3) + t623 * t764;
t627 = t138 * t961 - t405 * t963 - t348 * t527 / 0.2e1 - t446 * t517 / 0.2e1 + t986 * t960;
t2 = t627 + t672;
t631 = t478 * t738 + t527 * t936;
t36 = t631 + t668;
t71 = t517 * t527;
t652 = -t2 * qJD(1) + t36 * qJD(2) + t71 * qJD(3);
t639 = t139 * t793 + t550 * t716;
t21 = t786 / 0.2e1 + (t963 + t708) * t621 + t639;
t381 = t549 * t621 + t550 * t935;
t90 = t621 * t738 + t717 + t718;
t650 = -qJD(1) * t21 + qJD(3) * t90 - qJD(4) * t381;
t184 = -t517 * t541 + t527 * t536;
t629 = t446 * t954 + t477 * t955 - t640;
t634 = t549 * t941 + t655;
t29 = t629 + t634;
t649 = -qJD(1) * t29 + qJD(2) * t300 - qJD(3) * t184;
t185 = t517 * t536 + t527 * t541;
t630 = t446 * t953 + t478 * t955 + t937 * t980 + t671;
t633 = t550 * t941 + t656;
t32 = t630 + t633;
t648 = -qJD(1) * t32 + qJD(2) * t303 - qJD(3) * t185;
t647 = -t538 * qJD(4) + t625 * t710;
t187 = t750 + t657;
t26 = -t405 * t716 - t787 / 0.2e1 + (-t914 / 0.2e1 - t150 / 0.2e1) * t621 + t658;
t645 = -qJD(1) * t26 - qJD(2) * t187 - t517 * t843;
t46 = -t640 + t655;
t644 = qJD(1) * t46 + qJD(2) * t301 - t517 * t844;
t641 = t914 / 0.2e1 - t671;
t47 = -t641 + t656;
t643 = qJD(1) * t47 + qJD(2) * t306 + t517 * t845;
t628 = t635 * pkin(8) + t416 * t947 + t513 * t553 / 0.2e1;
t34 = t628 + t677;
t378 = (t947 + t678) * t625;
t642 = t553 * t567 * qJD(3) + t34 * qJD(1) - t378 * qJD(2);
t638 = -qJD(4) * t704 + t799;
t555 = t620 + t606;
t637 = t555 * qJD(1) + t720 - t839;
t602 = -t802 / 0.2e1;
t601 = t802 / 0.2e1;
t600 = -t801 / 0.2e1;
t599 = t801 / 0.2e1;
t591 = t624 * t797;
t590 = t624 * t798;
t589 = t622 * t801;
t566 = t624 * t757;
t554 = t592 * qJ(5);
t535 = -t590 + t840;
t533 = t592 * t622;
t532 = t621 * t798 - t611;
t529 = t713 * t623;
t526 = t883 / 0.2e1;
t516 = t981 * qJD(3);
t511 = -t589 + t764;
t510 = t622 * t797 + t624 * t800;
t509 = -t755 - t765;
t508 = t622 * t800 - t591;
t507 = qJD(3) * t617 + t722;
t488 = t617 * t941 + t619 * t940 + t941;
t486 = t734 + t534;
t485 = -t611 + (qJD(6) + t798) * t621;
t475 = (t769 + t843) * t625;
t474 = -t622 * t756 + t591;
t473 = t622 * t687;
t469 = (qJD(6) / 0.2e1 + t713) * t623;
t448 = t619 * t758 - t721;
t447 = t617 * t758 + t721;
t424 = t589 - t805;
t423 = -t764 + t805;
t422 = t622 * t476;
t411 = -t542 * qJD(3) + t623 * t765;
t398 = t420 * qJD(3);
t394 = -t619 * t756 - t808;
t393 = -t617 * t756 + t808;
t379 = t567 * t937 + t625 * t678;
t347 = t480 * t536;
t330 = t364 * qJD(4);
t290 = -t808 + (t619 * t802 + t587) * t625;
t289 = t808 + (t617 * t802 - t587) * t625;
t263 = 0.2e1 * t952 - t724;
t252 = t624 * t681 + t526 + t547;
t251 = t575 + t523 + (t596 * t624 + t790) * t623;
t199 = t889 / 0.2e1 + t541 * t747 + t621 * t940;
t198 = t752 + t536 * t747 - t774 / 0.2e1;
t190 = -t776 + t893;
t186 = t750 - t657;
t179 = t180 * qJD(6);
t175 = t784 + t908;
t174 = t478 * t621 - t778;
t153 = t479 * t945 - t777 / 0.2e1 + (t715 - t893 / 0.2e1) * t625;
t152 = pkin(8) * t746 + t553 * t748 - t911 / 0.2e1 - t883 / 0.2e1 + t736;
t94 = t624 * t653 + t526 + t673 - t795;
t92 = t574 + t613 + t596 * t748 + t522 - t707 * t624 + (t749 + t957) * t622;
t91 = t373 + t915 / 0.2e1 + t980 * t944 + t718;
t59 = t363 * t939 + t315 * t943 + t316 * t938 + pkin(4) * t746 + (t962 + t744) * t622;
t49 = t640 + t655;
t48 = t641 + t656;
t45 = t550 * t937 + t782 / 0.2e1 + t636 * t944 + t709;
t44 = t632 + t651;
t38 = qJD(3) * t60;
t35 = t631 - t668;
t33 = t628 - t677;
t31 = t630 - t633;
t30 = t629 - t634;
t25 = t150 * t944 + t787 / 0.2e1 + (t717 - t915 / 0.2e1) * t625 + t658;
t22 = qJD(3) * t39 + qJD(4) * t54;
t20 = t156 * t944 - t786 / 0.2e1 + t708 * t621 + t639;
t4 = t626 + t667;
t3 = -t627 + t672;
t1 = qJD(3) * t8 + qJD(4) * t12 - qJD(5) * t181;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t758, t581 * qJD(3), 0, -t758, 0, 0, t597 * t801, t597 * t797, 0, 0, t448, -t647, t411, t447, t654, -t758, -qJD(3) * t158 - qJD(4) * t202, qJD(3) * t159 + qJD(4) * t201, -qJD(3) * t70, qJD(3) * t78, t448, t411, t647, -t758, -t654, t447, qJD(3) * t73 + qJD(4) * t147 - t618 * t762, -t65 * qJD(3) - t63 * qJD(4) + t625 * t757, -qJD(3) * t72 - qJD(4) * t148 + qJD(5) * t555, qJD(3) * t52 + qJD(4) * t57 - qJD(5) * t169 (qJD(3) * t480 + t864) * t478, qJD(3) * t165 - t183 * t796, t239 * qJD(3) + t477 * t726 (qJD(3) * t479 - t862) * t477, t238 * qJD(3) + t478 * t726, -t758, qJD(3) * t27 + qJD(4) * t55 + qJD(5) * t365 + qJD(6) * t69, -qJD(3) * t28 - qJD(4) * t56 + qJD(5) * t366 - qJD(6) * t68, qJD(3) * t14 + qJD(4) * t17 - qJD(5) * t200, qJD(3) * t18 + qJD(4) * t19 - qJD(5) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t756, t803, t797, -t756, -t801, 0, -t596 * t797 + t761, t596 * t801 + t760, 0, 0, t290, -t968, t424, t289, t973, -t529, -t832 + (t622 * t711 - t548) * qJD(3) + t252 * qJD(4), t831 + (t624 * t711 + t788) * qJD(3) + t251 * qJD(4), qJD(3) * t694 - t873 (-pkin(3) * t878 + pkin(8) * t694) * qJD(3) + t696, t290, t424, t968, -t529, -t973, t289, t871 + (-t417 * t624 - t622 * t703) * qJD(3) + t94 * qJD(4) + t488 * qJD(5), qJD(3) * t695 + t59 * qJD(4) - t876, -t872 + (-t417 * t622 + t624 * t703) * qJD(3) + t92 * qJD(4) + t566 (pkin(8) * t695 + t417 * t553) * qJD(3) + t33 * qJD(4) + t152 * qJD(5) + t698 (t844 + t847) * t480 + t869, t830 + (-t347 - t907) * qJD(3) - t870, -t766 + t849 + t866 (t845 + t848) * t479 - t869, -qJD(4) * t301 - qJD(6) * t300 + t767 + t850, -t469, t925 + (t349 * t536 + t479 * t517 + t623 * t980) * qJD(3) + t30 * qJD(4) + t198 * qJD(5) + t49 * qJD(6), -t924 + (t349 * t541 - t405 * t623 + t480 * t517) * qJD(3) + t31 * qJD(4) + t199 * qJD(5) + t48 * qJD(6), t927 + (-t149 * t541 - t150 * t536 + t405 * t479 + t480 * t980) * qJD(3) + t4 * qJD(4) + t153 * qJD(5) (-t149 * t980 - t150 * t405 + t349 * t517) * qJD(3) + t3 * qJD(4) + t25 * qJD(5) + t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t395, -t410, t473, t395, t476, t599, qJD(3) * t252 - t330 - t851, qJD(3) * t251 + t816 + t852, 0, 0, -t395, t473, t410, t599, -t476, t395, qJD(3) * t94 - t330 + t856, -t877 + t59 * qJD(3) + (-qJ(5) * t881 + t594) * qJD(4) - t757, t92 * qJD(3) - t816 - t833 - t835, t33 * qJD(3) + (-pkin(4) * t364 - qJ(5) * t363) * qJD(4) + t315 * qJD(5) + t697, t978, -t686, -t822 - t983, -t978, -t825 - t984, t599, t30 * qJD(3) + t155 * qJD(4) + t45 * qJD(6) - t763 + t903, t31 * qJD(3) + t156 * qJD(4) + t44 * qJD(6) - t719 - t901, t926 + t4 * qJD(3) + (t549 * t477 + t478 * t550) * qJD(4) + t174 * qJD(5), t3 * qJD(3) + (t155 * t549 + t156 * t550) * qJD(4) + t20 * qJD(5) + t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t488 - t565, t473, t637, qJD(3) * t152 + qJD(4) * t315 - t829, 0, 0, 0, 0, 0, 0, t198 * qJD(3) - t611 * t625 + t815, t199 * qJD(3) - t612 * t625 + t814, qJD(3) * t153 + qJD(4) * t174 - t853, qJD(3) * t25 + qJD(4) * t20 + t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t978, t686, -t817 + t983, t978, -t826 + t984, t600, qJD(3) * t49 + qJD(4) * t45 - qJD(6) * t139 + t874, qJD(3) * t48 + qJD(4) * t44 + qJD(6) * t138 - t875, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t398, 0, 0, 0, 0, 0, 0, 0, 0, 0, t398, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t801, -t797, 0, 0, 0, 0, 0, 0, 0, 0, t509, -t511, t516 (t857 - t933) * qJD(3) - t971, 0, 0, 0, 0, 0, 0, t509, t516, t511 (t857 + t885) * qJD(3) + t379 * qJD(4) + t622 * t835 + t970, 0, 0, 0, 0, 0, 0, -qJD(4) * t302 - qJD(6) * t299 - t767, -t766 + t867 (-t347 + t907) * qJD(3) (-t405 * t480 + t479 * t980 - t517 * t623) * qJD(3) + t35 * qJD(4) + t186 * qJD(5) + t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t510, t508, 0, 0, 0, 0, 0, 0, 0, 0, -t510, 0, -t508, qJD(3) * t379 - qJD(4) * t513 + t623 * t799 + t929, 0, 0, 0, 0, 0, 0, -t823 - t862, -t818 + t864, 0, t928 + t35 * qJD(3) + (t477 * t550 - t478 * t549) * qJD(4) + t175 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t510, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t186 + qJD(4) * t175 - t854; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t827 + t862, -t820 - t864, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t756, -t803, 0, t756, 0, 0, -t761, -t760, 0, 0, t394, 0.2e1 * t422, t423, t393, -t972, t529, qJD(4) * t383 + t832, qJD(4) * t382 - t831, t873, -t696, t394, t423, -0.2e1 * t422, t529, t972, t393, qJD(4) * t95 - qJD(5) * t487 - t871, -t58 * qJD(4) - t625 * t799 + t876, qJD(4) * t93 + t566 + t872, qJD(4) * t34 - qJD(5) * t151 - t698, -t480 * t847 + t869, -t830 - t870, -t849 + t867, -t479 * t848 - t869, -qJD(4) * t299 - qJD(6) * t302 - t850, t469, qJD(4) * t29 + qJD(5) * t197 - qJD(6) * t46 - t925, qJD(4) * t32 - qJD(5) * t196 - qJD(6) * t47 + t924, qJD(4) * t5 + qJD(5) * t154 - t927, -qJD(4) * t2 + qJD(5) * t26 - t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t971, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t378 - t970, 0, 0, 0, 0, 0, 0, -qJD(4) * t300 - qJD(6) * t301, t866, 0, qJD(4) * t36 + qJD(5) * t187 - t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t588, t804, 0, -t588, 0, 0, -pkin(3) * t841, -pkin(3) * t840, 0, 0, t588, 0, -t804, 0, 0, -t588, -qJD(4) * t389 + t762, 0, -qJD(4) * t388 + qJD(5) * t617 (qJD(4) * t567 - t836) * t553, t982, -t796 * t260, 0, -t982, 0, 0, qJD(4) * t184 + t536 * t836 + t541 * t834, qJD(4) * t185 - t536 * t834 + t541 * t836, 0, qJD(4) * t71 + t517 * t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, -t449, t535, t418, t533, t602, -t674 - t791, -t675 + t792, 0, 0, -t418, t535, t449, t602, -t533, t418, -t690 - t791, t638 - t900, -t691 - t792, pkin(8) * t638 + t642, t979, -t685, -t819 - t977, -t979, -t975 - t828, t602, -t649 + t985, qJD(4) * t980 + qJD(6) * t263 - t648 (t549 * t536 + t550 * t541) * qJD(4) + t190 * qJD(5) + t701 (t405 * t549 + t550 * t980) * qJD(4) + t91 * qJD(5) + t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t729, t535, t507, -t666 + t791, 0, 0, 0, 0, 0, 0, t660, t661, qJD(4) * t190 + t855, qJD(4) * t91 + t179 - t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t979, t685, -t821 + t977, t979, t975 - t824, t601, -t644 - t985, qJD(4) * t263 + qJD(6) * t980 - t643, 0, t838; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t395, t410, t474, -t395, -t475, t599, -qJD(3) * t383 + t851, -qJD(3) * t382 - t852, 0, 0, t395, t474, -t410, t599, t475, -t395, -qJD(3) * t95 - t856, qJD(3) * t58 + t877, -t93 * qJD(3) + t833 - t835, -qJ(5) * t835 - t34 * qJD(3) - t697, -t978, t686, -t771 + t818, t978, -t770 + t827, t599, -t29 * qJD(3) - t42 * qJD(6) - t763 - t903, -t32 * qJD(3) - t43 * qJD(6) - t719 + t901, -qJD(3) * t5 - t926, qJD(3) * t2 + qJD(5) * t21 - t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t378 - t929, 0, 0, 0, 0, 0, 0, t826, t822, 0, -qJD(3) * t36 - t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, t449, t590, -t418, -t759, t601, t674, t675, 0, 0, t418, t590, -t449, t601, t759, -t418, t690, t900, t691, -t642, -t979, t685, t819, t979, t828, t601, t649, qJD(6) * t392 + t648, -t701, -qJD(5) * t90 - t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5), 0, 0, 0, 0, 0, 0, qJD(6) * t550 + t837, -t549 * qJD(6) + t733, 0, qJD(5) * t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t592, -t554, 0, 0, 0, 0, 0, 0, -t532, -t534, 0, -t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t550 * t796 - t909, -t549 * t796 + t689, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t487 + t565, t474, -t637, qJ(5) * t839 + t151 * qJD(3) + t829, 0, 0, 0, 0, 0, 0, -t197 * qJD(3) + t621 * t726 - t815, -t814 + t196 * qJD(3) + (t612 - t734) * t625, -qJD(3) * t154 + t853, -qJD(3) * t26 - qJD(4) * t21 - t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t187 + t854; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t729, t590, -t507, t666, 0, 0, 0, 0, 0, 0, -t660, -t661, -t855, qJD(4) * t90 + t179 + t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t592, t554, 0, 0, 0, 0, 0, 0, t485, t486, 0, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, -t486, 0, t846; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t978, -t686, t771 + t820, -t978, t770 + t823, t600, t46 * qJD(3) + t42 * qJD(4) + t763 - t874, t47 * qJD(3) + t43 * qJD(4) + t719 + t875, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t825, t817, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t979, -t685, t821, -t979, t824, t602, t644, -qJD(4) * t392 + t643, 0, -t838; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t550 - t837 + t909, t549 * qJD(4) - t689 - t733, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t532, t534, 0, -t846; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t6;