% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:38
% EndTime: 2019-02-26 19:57:40
% DurationCPUTime: 2.16s
% Computational Cost: add. (3033->247), mult. (9956->443), div. (0->0), fcn. (11639->18), ass. (0->156)
t769 = sin(pkin(8));
t774 = cos(pkin(8));
t768 = sin(pkin(14));
t772 = cos(pkin(14));
t770 = sin(pkin(7));
t775 = cos(pkin(7));
t773 = cos(pkin(13));
t779 = sin(qJ(2));
t782 = cos(qJ(2));
t863 = sin(pkin(13));
t864 = cos(pkin(6));
t827 = t864 * t863;
t811 = t773 * t779 + t782 * t827;
t771 = sin(pkin(6));
t836 = t771 * t863;
t795 = t770 * t836 - t811 * t775;
t808 = -t773 * t782 + t779 * t827;
t788 = -t768 * t808 - t795 * t772;
t796 = t811 * t770 + t775 * t836;
t870 = t796 * t769 - t788 * t774;
t835 = t773 * t864;
t810 = t863 * t779 - t782 * t835;
t852 = t771 * t773;
t798 = -t770 * t852 - t810 * t775;
t809 = -t779 * t835 - t863 * t782;
t790 = -t768 * t809 - t798 * t772;
t799 = t810 * t770 - t775 * t852;
t869 = t799 * t769 - t790 * t774;
t847 = t775 * t782;
t818 = -t768 * t779 + t772 * t847;
t834 = t864 * t770;
t801 = t818 * t771 + t772 * t834;
t853 = t770 * t782;
t814 = -t771 * t853 + t864 * t775;
t867 = t814 * t769 + t801 * t774;
t865 = cos(qJ(4));
t872 = t870 * t865;
t871 = t869 * t865;
t868 = t867 * t865;
t776 = sin(qJ(6));
t780 = cos(qJ(6));
t816 = qJD(6) * (t776 * r_i_i_C(1) + t780 * r_i_i_C(2));
t866 = r_i_i_C(3) + pkin(12);
t757 = t810 * qJD(2);
t758 = t809 * qJD(2);
t851 = t772 * t775;
t734 = t757 * t851 - t758 * t768;
t862 = t734 * t769;
t759 = t811 * qJD(2);
t760 = t808 * qJD(2);
t737 = t759 * t851 - t760 * t768;
t861 = t737 * t769;
t846 = qJD(2) * t771;
t751 = t818 * t846;
t860 = t751 * t769;
t848 = t775 * t779;
t817 = t768 * t782 + t772 * t848;
t755 = t817 * t771;
t859 = t755 * t774;
t858 = t768 * t775;
t857 = t769 * t770;
t778 = sin(qJ(4));
t856 = t769 * t778;
t855 = t770 * t774;
t854 = t770 * t779;
t850 = t772 * t779;
t845 = qJD(4) * t778;
t844 = qJD(6) * t776;
t843 = qJD(6) * t780;
t842 = t770 * t856;
t841 = t774 * t854;
t840 = t774 * t865;
t839 = t770 * t846;
t838 = qJD(4) * t865;
t837 = pkin(10) * t774 + qJ(3);
t832 = t865 * t857;
t831 = t782 * t839;
t830 = t779 * t839;
t826 = t771 * t832;
t729 = t798 * t768 - t772 * t809;
t704 = t729 * t865 + t869 * t778;
t716 = t790 * t769 + t799 * t774;
t777 = sin(qJ(5));
t781 = cos(qJ(5));
t694 = t704 * t781 + t716 * t777;
t825 = -t704 * t777 + t716 * t781;
t730 = t795 * t768 - t772 * t808;
t706 = t730 * t865 + t870 * t778;
t717 = t788 * t769 + t796 * t774;
t696 = t706 * t781 + t717 * t777;
t824 = -t706 * t777 + t717 * t781;
t739 = t810 * t768 + t809 * t851;
t740 = -t810 * t772 + t809 * t858;
t710 = t740 * t865 + (t739 * t774 - t809 * t857) * t778;
t722 = -t739 * t769 - t809 * t855;
t697 = t710 * t781 + t722 * t777;
t741 = t811 * t768 + t808 * t851;
t742 = -t811 * t772 + t808 * t858;
t712 = t742 * t865 + (t741 * t774 - t808 * t857) * t778;
t723 = -t741 * t769 - t808 * t855;
t698 = t712 * t781 + t723 * t777;
t749 = t771 * t850 + (t771 * t847 + t834) * t768;
t715 = t749 * t865 + t867 * t778;
t728 = -t801 * t769 + t814 * t774;
t702 = t715 * t781 + t728 * t777;
t823 = -t715 * t777 + t728 * t781;
t756 = (-t768 * t848 + t772 * t782) * t771;
t725 = t756 * t865 + (t769 * t771 * t854 - t859) * t778;
t745 = t755 * t769 + t771 * t841;
t713 = t725 * t781 + t745 * t777;
t822 = r_i_i_C(1) * t780 - r_i_i_C(2) * t776 + pkin(5);
t820 = t757 * t768 + t758 * t851;
t819 = t759 * t768 + t760 * t851;
t815 = qJD(2) * t826;
t813 = t820 * t774;
t812 = t819 * t774;
t806 = qJD(2) * t859;
t804 = -t866 * t777 - t822 * t781 - pkin(4);
t803 = t739 * t840 - t740 * t778 - t809 * t832;
t802 = t741 * t840 - t742 * t778 - t808 * t832;
t800 = -t755 * t840 - t756 * t778 + t779 * t826;
t786 = t781 * t816 + (t822 * t777 - t866 * t781) * qJD(5);
t753 = qJD(2) * t756;
t752 = (-t768 * t847 - t850) * t846;
t744 = (t817 * t769 + t841) * t846;
t743 = t774 * t831 + t860;
t738 = t759 * t858 + t760 * t772;
t736 = -t759 * t772 + t760 * t858;
t735 = t757 * t858 + t758 * t772;
t733 = -t757 * t772 + t758 * t858;
t721 = -t759 * t855 - t861;
t720 = -t760 * t855 - t819 * t769;
t719 = -t757 * t855 - t862;
t718 = -t758 * t855 - t820 * t769;
t714 = t749 * t778 - t868;
t708 = t752 * t865 + (-t751 * t774 + t769 * t831) * t778 + t800 * qJD(4);
t707 = t725 * qJD(4) + t751 * t840 + t752 * t778 - t782 * t815;
t705 = t730 * t778 - t872;
t703 = t729 * t778 - t871;
t700 = t868 * qJD(4) - t749 * t845 + t753 * t865 - t778 * t806 + t830 * t856;
t699 = t749 * t838 + t753 * t778 - t779 * t815 + t865 * t806 + t867 * t845;
t692 = t738 * t865 + (t737 * t774 - t759 * t857) * t778 + t802 * qJD(4);
t691 = t712 * qJD(4) - t737 * t840 + t738 * t778 + t759 * t832;
t690 = t735 * t865 + (t734 * t774 - t757 * t857) * t778 + t803 * qJD(4);
t689 = t710 * qJD(4) - t734 * t840 + t735 * t778 + t757 * t832;
t688 = t872 * qJD(4) - t730 * t845 + t736 * t865 - t760 * t842 + t778 * t812;
t687 = t730 * t838 + t736 * t778 + t760 * t832 - t865 * t812 + t870 * t845;
t686 = t871 * qJD(4) - t729 * t845 + t733 * t865 - t758 * t842 + t778 * t813;
t685 = t729 * t838 + t733 * t778 + t758 * t832 - t865 * t813 + t869 * t845;
t684 = t708 * t781 + t743 * t777 + (-t725 * t777 + t745 * t781) * qJD(5);
t682 = t823 * qJD(5) + t700 * t781 + t744 * t777;
t680 = t692 * t781 + t721 * t777 + (-t712 * t777 + t723 * t781) * qJD(5);
t678 = t690 * t781 + t719 * t777 + (-t710 * t777 + t722 * t781) * qJD(5);
t676 = t824 * qJD(5) + t688 * t781 + t720 * t777;
t674 = t825 * qJD(5) + t686 * t781 + t718 * t777;
t1 = [0 (t680 * t780 + t691 * t776) * r_i_i_C(1) + (-t680 * t776 + t691 * t780) * r_i_i_C(2) + t680 * pkin(5) + t692 * pkin(4) + t691 * pkin(11) + t738 * pkin(3) - pkin(10) * t861 + t760 * pkin(2) + t866 * (t698 * qJD(5) + t692 * t777 - t721 * t781) + ((-t698 * t776 - t780 * t802) * r_i_i_C(1) + (-t698 * t780 + t776 * t802) * r_i_i_C(2)) * qJD(6) + (-qJD(3) * t808 - t837 * t759) * t770, -t760 * t770 (t688 * t776 + t706 * t843) * r_i_i_C(1) + (t688 * t780 - t706 * t844) * r_i_i_C(2) + t688 * pkin(11) + t804 * t687 + t786 * t705, t866 * t676 - t824 * t816 + t822 * (-qJD(5) * t696 - t688 * t777 + t720 * t781) (-t676 * t776 + t687 * t780) * r_i_i_C(1) + (-t676 * t780 - t687 * t776) * r_i_i_C(2) + ((-t696 * t780 - t705 * t776) * r_i_i_C(1) + (t696 * t776 - t705 * t780) * r_i_i_C(2)) * qJD(6); 0 (t678 * t780 + t689 * t776) * r_i_i_C(1) + (-t678 * t776 + t689 * t780) * r_i_i_C(2) + t678 * pkin(5) + t690 * pkin(4) + t689 * pkin(11) + t735 * pkin(3) - pkin(10) * t862 + t758 * pkin(2) + t866 * (t697 * qJD(5) + t690 * t777 - t719 * t781) + ((-t697 * t776 - t780 * t803) * r_i_i_C(1) + (-t697 * t780 + t776 * t803) * r_i_i_C(2)) * qJD(6) + (-qJD(3) * t809 - t837 * t757) * t770, -t758 * t770 (t686 * t776 + t704 * t843) * r_i_i_C(1) + (t686 * t780 - t704 * t844) * r_i_i_C(2) + t686 * pkin(11) + t804 * t685 + t786 * t703, t866 * t674 - t825 * t816 + t822 * (-qJD(5) * t694 - t686 * t777 + t718 * t781) (-t674 * t776 + t685 * t780) * r_i_i_C(1) + (-t674 * t780 - t685 * t776) * r_i_i_C(2) + ((-t694 * t780 - t703 * t776) * r_i_i_C(1) + (t694 * t776 - t703 * t780) * r_i_i_C(2)) * qJD(6); 0 (t684 * t780 + t707 * t776) * r_i_i_C(1) + (-t684 * t776 + t707 * t780) * r_i_i_C(2) + t684 * pkin(5) + t708 * pkin(4) + t707 * pkin(11) + t752 * pkin(3) + pkin(10) * t860 + t866 * (t713 * qJD(5) + t708 * t777 - t743 * t781) + ((-t713 * t776 - t780 * t800) * r_i_i_C(1) + (-t713 * t780 + t776 * t800) * r_i_i_C(2)) * qJD(6) + (qJD(3) * t854 + (-pkin(2) * t779 + t837 * t853) * qJD(2)) * t771, t830 (t700 * t776 + t715 * t843) * r_i_i_C(1) + (t700 * t780 - t715 * t844) * r_i_i_C(2) + t700 * pkin(11) + t804 * t699 + t786 * t714, t866 * t682 - t823 * t816 + t822 * (-qJD(5) * t702 - t700 * t777 + t744 * t781) (-t682 * t776 + t699 * t780) * r_i_i_C(1) + (-t682 * t780 - t699 * t776) * r_i_i_C(2) + ((-t702 * t780 - t714 * t776) * r_i_i_C(1) + (t702 * t776 - t714 * t780) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
