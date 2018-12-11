% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_6_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_6_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:25
% EndTime: 2018-12-10 18:38:38
% DurationCPUTime: 13.02s
% Computational Cost: add. (207561->344), mult. (217255->585), div. (983->12), fcn. (208653->33), ass. (0->253)
t886 = pkin(6) + qJ(2);
t856 = cos(t886) / 0.2e1;
t887 = pkin(6) - qJ(2);
t871 = cos(t887);
t825 = t871 / 0.2e1 + t856;
t917 = sin(qJ(2));
t920 = cos(qJ(1));
t859 = t920 * t917;
t918 = sin(qJ(1));
t809 = t825 * t918 + t859;
t868 = sin(t886);
t773 = -t868 / 0.2e1;
t869 = sin(t887);
t774 = t869 / 0.2e1;
t892 = t774 + t773;
t829 = t892 * qJD(2);
t919 = cos(qJ(2));
t858 = t918 * t919;
t735 = qJD(1) * t809 + qJD(2) * t858 - t829 * t920;
t852 = t868 / 0.2e1;
t853 = -t869 / 0.2e1;
t824 = t852 + t853;
t860 = t920 * t919;
t810 = -t824 * t918 + t860;
t755 = t825 * qJD(2);
t778 = t918 * t917;
t848 = qJD(2) * t778 - t755 * t920;
t736 = qJD(1) * t810 - t848;
t881 = pkin(7) + pkin(14);
t843 = sin(t881) / 0.2e1;
t882 = pkin(7) - pkin(14);
t861 = sin(t882);
t757 = t843 - t861 / 0.2e1;
t779 = cos(pkin(14));
t844 = cos(t882) / 0.2e1;
t862 = cos(t881);
t758 = t844 - t862 / 0.2e1;
t914 = sin(pkin(6));
t874 = t758 * t914;
t840 = t918 * t874;
t703 = qJD(1) * t840 - t735 * t757 + t736 * t779;
t821 = t843 + t861 / 0.2e1;
t816 = t821 * t914;
t812 = t918 * t816;
t822 = t844 + t862 / 0.2e1;
t911 = sin(pkin(14));
t704 = -qJD(1) * t812 + t735 * t822 + t736 * t911;
t808 = t824 * t920 + t858;
t813 = t920 * t816;
t922 = -t825 * t920 + t778;
t720 = t808 * t911 + t822 * t922 + t813;
t721 = t757 * t922 - t779 * t808 + t874 * t920;
t884 = pkin(8) + qJ(4);
t866 = sin(t884);
t850 = -t866 / 0.2e1;
t885 = pkin(8) - qJ(4);
t867 = sin(t885);
t851 = -t867 / 0.2e1;
t750 = (t850 + t851) * qJD(4);
t854 = cos(t884) / 0.2e1;
t870 = cos(t885);
t839 = t870 / 0.2e1 + t854;
t752 = t839 * qJD(4);
t849 = t866 / 0.2e1;
t759 = t849 + t851;
t761 = t854 - t870 / 0.2e1;
t786 = cos(qJ(4));
t916 = cos(pkin(7));
t846 = t916 * t914;
t831 = t918 * t846;
t913 = sin(pkin(7));
t815 = qJD(1) * t831 + t735 * t913;
t783 = sin(qJ(4));
t891 = qJD(4) * t783;
t832 = t920 * t846;
t921 = -t913 * t922 + t832;
t647 = -t703 * t786 + t704 * t759 + t720 * t752 - t721 * t891 - t750 * t921 + t761 * t815;
t686 = t720 * t759 + t721 * t786 - t761 * t921;
t782 = sin(qJ(5));
t785 = cos(qJ(5));
t912 = sin(pkin(8));
t915 = cos(pkin(8));
t798 = t720 * t912 - t915 * t921;
t668 = t686 * t782 + t785 * t798;
t801 = t704 * t912 + t815 * t915;
t626 = qJD(5) * t668 - t647 * t785 + t782 * t801;
t669 = t686 * t785 - t782 * t798;
t936 = qJD(5) * t669 + t647 * t782 + t785 * t801;
t772 = t867 / 0.2e1;
t751 = (t772 + t850) * qJD(4);
t753 = t761 * qJD(4);
t838 = t849 + t772;
t890 = qJD(4) * t786;
t931 = -t703 * t783 - t704 * t839 - t720 * t751 + t721 * t890 - t753 * t921 + t815 * t838;
t930 = -t720 * t839 + t721 * t783 - t838 * t921;
t760 = t852 + t774;
t762 = t856 - t871 / 0.2e1;
t780 = cos(pkin(6));
t731 = t760 * t822 + t762 * t911 + t780 * t821;
t732 = t757 * t760 + t758 * t780 - t762 * t779;
t754 = (t773 + t853) * qJD(2);
t756 = t762 * qJD(2);
t738 = t754 * t911 + t756 * t822;
t739 = -t754 * t779 + t756 * t757;
t745 = -t760 * t913 + t780 * t916;
t873 = t761 * t913;
t673 = t731 * t752 - t732 * t891 + t738 * t759 + t739 * t786 - t745 * t750 + t756 * t873;
t699 = t731 * t759 + t732 * t786 - t745 * t761;
t717 = -t731 * t912 + t745 * t915;
t681 = t699 * t785 + t717 * t782;
t845 = t915 * t913;
t724 = -t738 * t912 - t756 * t845;
t633 = qJD(5) * t681 + t673 * t782 - t724 * t785;
t680 = t699 * t782 - t717 * t785;
t678 = 0.1e1 / t680 ^ 2;
t924 = t633 * t678;
t677 = 0.1e1 / t680;
t698 = t731 * t839 - t732 * t783 + t745 * t838;
t897 = t668 * t678;
t833 = -t677 * t930 - t698 * t897;
t923 = t782 * t833;
t640 = atan2(t668, t680);
t635 = sin(t640);
t636 = cos(t640);
t622 = t635 * t668 + t636 * t680;
t619 = 0.1e1 / t622;
t799 = t809 * t822 + t810 * t911 - t812;
t804 = -t809 * t913 - t831;
t795 = t799 * t912 - t804 * t915;
t722 = -t757 * t809 + t779 * t810 + t840;
t796 = t722 * t786 - t759 * t799 + t761 * t804;
t671 = t782 * t795 + t785 * t796;
t687 = t722 * t783 + t799 * t839 + t804 * t838;
t781 = sin(qJ(6));
t784 = cos(qJ(6));
t655 = t671 * t784 + t687 * t781;
t649 = 0.1e1 / t655;
t620 = 0.1e1 / t622 ^ 2;
t650 = 0.1e1 / t655 ^ 2;
t670 = t782 * t796 - t785 * t795;
t664 = t670 ^ 2;
t618 = t620 * t664 + 0.1e1;
t811 = -qJD(2) * t860 - t829 * t918;
t823 = -qJD(2) * t859 - t755 * t918;
t702 = qJD(1) * t721 + t811 * t757 + t823 * t779;
t803 = -qJD(1) * t922 - t811;
t797 = (-qJD(1) * t808 + t823) * t911 + t803 * t822 - qJD(1) * t813;
t802 = -qJD(1) * t832 - t803 * t913;
t793 = t702 * t786 - t722 * t891 + t750 * t804 - t752 * t799 - t759 * t797 + t761 * t802;
t794 = t797 * t912 - t802 * t915;
t623 = qJD(5) * t671 + t782 * t793 - t785 * t794;
t904 = t620 * t670;
t663 = t668 ^ 2;
t639 = t663 * t678 + 0.1e1;
t637 = 0.1e1 / t639;
t837 = -t633 * t897 + t677 * t936;
t609 = t837 * t637;
t842 = -t635 * t680 + t636 * t668;
t603 = t609 * t842 + t633 * t636 + t635 * t936;
t621 = t619 * t620;
t909 = t603 * t621;
t910 = 0.2e1 * (t623 * t904 - t664 * t909) / t618 ^ 2;
t624 = -qJD(5) * t670 + t782 * t794 + t785 * t793;
t641 = t702 * t783 + t722 * t890 + t751 * t799 + t753 * t804 + t797 * t839 + t802 * t838;
t614 = qJD(6) * t655 + t624 * t781 - t641 * t784;
t654 = t671 * t781 - t687 * t784;
t648 = t654 ^ 2;
t630 = t648 * t650 + 0.1e1;
t900 = t650 * t654;
t888 = qJD(6) * t654;
t615 = t624 * t784 + t641 * t781 - t888;
t906 = t615 * t649 * t650;
t908 = 0.2e1 * (t614 * t900 - t648 * t906) / t630 ^ 2;
t903 = t677 * t924;
t907 = 0.2e1 * (-t663 * t903 + t897 * t936) / t639 ^ 2;
t905 = t620 * t623;
t902 = t635 * t670;
t901 = t636 * t670;
t899 = t654 * t784;
t898 = t668 * t677;
t896 = t687 * t782;
t895 = t687 * t785;
t893 = t781 * t649;
t889 = qJD(5) * t785;
t883 = 0.2e1 * t621 * t670;
t880 = t677 * t907;
t879 = t620 * t902;
t878 = t620 * t901;
t877 = t654 * t906;
t875 = t750 * t913;
t872 = t762 * t913;
t865 = t603 * t883;
t864 = 0.2e1 * t668 * t903;
t863 = 0.2e1 * t877;
t847 = qJD(6) * t895 + t793;
t653 = t669 * t784 + t781 * t930;
t652 = t669 * t781 - t784 * t930;
t747 = -t892 * t918 - t860;
t729 = t747 * t822 + t809 * t911;
t730 = t747 * t757 - t779 * t809;
t697 = t729 * t759 + t730 * t786 + t747 * t873;
t715 = -t729 * t912 - t747 * t845;
t676 = t697 * t785 + t715 * t782;
t827 = t838 * t913;
t696 = -t729 * t839 + t730 * t783 + t747 * t827;
t660 = t676 * t784 + t696 * t781;
t659 = t676 * t781 - t696 * t784;
t675 = t697 * t782 - t715 * t785;
t836 = t650 * t899 - t893;
t835 = t669 * t677 - t681 * t897;
t820 = -t892 * t920 + t858;
t727 = -t820 * t822 + t911 * t922;
t728 = -t757 * t820 - t779 * t922;
t695 = t727 * t759 + t728 * t786 - t820 * t873;
t819 = -t727 * t912 + t820 * t845;
t674 = t695 * t782 - t785 * t819;
t742 = -t760 * t911 + t762 * t822;
t743 = t757 * t762 + t760 * t779;
t708 = t742 * t759 + t743 * t786 + t761 * t872;
t726 = -t742 * t912 - t762 * t845;
t693 = t708 * t782 - t726 * t785;
t834 = -t674 * t677 - t693 * t897;
t828 = -t635 + (-t636 * t898 + t635) * t637;
t826 = qJD(5) * t896 + qJD(6) * t796 - t641 * t785;
t740 = t754 * t822 - t756 * t911;
t737 = qJD(1) * t747 + t848;
t734 = qJD(1) * t820 - t823;
t712 = t735 * t911 + t737 * t822;
t711 = t734 * t757 - t779 * t803;
t710 = t734 * t822 + t803 * t911;
t694 = -t710 * t912 - t734 * t845;
t672 = t731 * t751 - t732 * t890 + t738 * t839 - t739 * t783 + t745 * t753 - t756 * t827;
t662 = t781 * t796 - t784 * t895;
t661 = -t781 * t895 - t784 * t796;
t658 = t710 * t759 + t711 * t786 + t729 * t752 - t730 * t891 + t734 * t873 + t747 * t875;
t657 = t747 * t753 * t913 - t710 * t839 + t711 * t783 - t729 * t751 + t730 * t890 + t734 * t827;
t656 = ((t754 * t757 + t756 * t779) * t786 - t743 * t891 + t740 * t759 + t742 * t752 + t754 * t873 + t750 * t872) * t782 - (-t740 * t912 - t754 * t845) * t785 + (t708 * t785 + t726 * t782) * qJD(5);
t634 = -qJD(5) * t680 + t673 * t785 + t724 * t782;
t632 = ((-t735 * t779 + t737 * t757) * t786 - t728 * t891 + t712 * t759 + t727 * t752 + t737 * t873 - t820 * t875) * t782 - (-t712 * t912 - t737 * t845) * t785 + (t695 * t785 + t782 * t819) * qJD(5);
t631 = -qJD(5) * t675 + t658 * t785 + t694 * t782;
t628 = 0.1e1 / t630;
t616 = 0.1e1 / t618;
t613 = t637 * t923;
t612 = t834 * t637;
t611 = t835 * t637;
t608 = t828 * t670;
t606 = (-t635 * t930 + t636 * t698) * t782 + t842 * t613;
t605 = t612 * t842 - t635 * t674 + t636 * t693;
t604 = t611 * t842 + t635 * t669 + t636 * t681;
t601 = -t834 * t907 + (t693 * t864 - t632 * t677 + (t633 * t674 - t656 * t668 - t693 * t936) * t678) * t637;
t600 = -t835 * t907 + (t681 * t864 - t626 * t677 + (-t633 * t669 - t634 * t668 - t681 * t936) * t678) * t637;
t599 = -t907 * t923 + (t833 * t889 + (t698 * t864 - t931 * t677 + (t633 * t930 - t668 * t672 - t698 * t936) * t678) * t782) * t637;
t1 = [t670 * t880 + (-t623 * t677 + t670 * t924) * t637, t601, 0, t599, t600, 0; -t668 * t619 * t910 + (t936 * t619 + (-t603 * t668 - t608 * t623) * t620) * t616 + (t608 * t620 * t910 + (0.2e1 * t608 * t909 - (t609 * t637 * t898 - t907) * t879 - (t668 * t880 - t609 + (t609 - t837) * t637) * t878 - t828 * t905) * t616) * t670 (t605 * t904 - t619 * t675) * t910 + ((qJD(5) * t676 + t658 * t782 - t694 * t785) * t619 + t605 * t865 + (-t675 * t603 - t605 * t623 - (t601 * t668 + t612 * t936 + t656 + (-t612 * t680 - t674) * t609) * t901 - (-t601 * t680 - t612 * t633 - t632 + (-t612 * t668 - t693) * t609) * t902) * t620) * t616, 0 (t606 * t904 + t619 * t896) * t910 + (-t606 * t905 + (-t641 * t782 - t687 * t889) * t619 + (t606 * t883 + t620 * t896) * t603 - (t698 * t889 + t599 * t668 + t613 * t936 + t672 * t782 + (-t613 * t680 - t782 * t930) * t609) * t878 - (-t930 * t889 - t599 * t680 - t613 * t633 - t931 * t782 + (-t613 * t668 - t698 * t782) * t609) * t879) * t616 (t604 * t904 - t619 * t671) * t910 + (t604 * t865 + t624 * t619 + (-t671 * t603 - t604 * t623 - (t600 * t668 + t611 * t936 + t634 + (-t611 * t680 + t669) * t609) * t901 - (-t600 * t680 - t611 * t633 - t626 + (-t611 * t668 - t681) * t609) * t902) * t620) * t616, 0; (-t649 * t652 + t653 * t900) * t908 + ((qJD(6) * t653 - t626 * t781 - t784 * t931) * t649 + t653 * t863 + (-t652 * t615 - (-qJD(6) * t652 - t626 * t784 + t781 * t931) * t654 - t653 * t614) * t650) * t628 (-t649 * t659 + t660 * t900) * t908 + ((qJD(6) * t660 + t631 * t781 - t657 * t784) * t649 + t660 * t863 + (-t659 * t615 - (-qJD(6) * t659 + t631 * t784 + t657 * t781) * t654 - t660 * t614) * t650) * t628, 0 (-t649 * t661 + t662 * t900) * t908 + (t662 * t863 - t847 * t649 * t784 + t826 * t893 + (-t654 * t781 * t847 - t662 * t614 - t661 * t615 - t826 * t899) * t650) * t628, -t836 * t670 * t908 + (t836 * t623 + ((-qJD(6) * t649 - 0.2e1 * t877) * t784 + (t614 * t784 + (t615 - t888) * t781) * t650) * t670) * t628, -t908 + (0.2e1 * t614 * t650 * t628 + (-0.2e1 * t628 * t906 - t650 * t908) * t654) * t654;];
JaD_rot  = t1;
