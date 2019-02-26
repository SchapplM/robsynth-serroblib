% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:31
% EndTime: 2019-02-26 22:51:34
% DurationCPUTime: 2.79s
% Computational Cost: add. (2908->254), mult. (7495->421), div. (0->0), fcn. (8135->16), ass. (0->159)
t741 = qJ(4) + qJ(5);
t739 = cos(t741);
t748 = sin(qJ(2));
t751 = cos(qJ(2));
t752 = cos(qJ(1));
t842 = cos(pkin(6));
t803 = t752 * t842;
t846 = sin(qJ(1));
t721 = t748 * t803 + t846 * t751;
t789 = t842 * t846;
t768 = t752 * t748 + t751 * t789;
t708 = t768 * qJD(1) + t721 * qJD(2);
t780 = t748 * t789;
t809 = t846 * t748;
t824 = t752 * t751;
t709 = -qJD(1) * t780 - qJD(2) * t809 + (qJD(2) * t842 + qJD(1)) * t824;
t720 = -t751 * t803 + t809;
t742 = sin(pkin(7));
t744 = cos(pkin(7));
t747 = sin(qJ(3));
t743 = sin(pkin(6));
t828 = t743 * t752;
t816 = t742 * t828;
t847 = cos(qJ(3));
t785 = t847 * t816;
t805 = t847 * qJD(3);
t791 = t744 * t805;
t804 = t846 * qJD(1);
t792 = t743 * t804;
t661 = (-qJD(3) * t721 - t708 * t744 + t742 * t792) * t747 - qJD(3) * t785 + t709 * t847 - t720 * t791;
t712 = -t720 * t742 + t744 * t828;
t740 = qJD(4) + qJD(5);
t801 = -t712 * t740 + t661;
t864 = t801 * t739;
t814 = t721 * t847;
t692 = (t720 * t744 + t816) * t747 - t814;
t738 = sin(t741);
t677 = t692 * t739 + t712 * t738;
t813 = t744 * t847;
t689 = t720 * t813 + t721 * t747 + t785;
t745 = sin(qJ(6));
t749 = cos(qJ(6));
t863 = -t677 * t745 - t689 * t749;
t862 = t677 * t749 - t689 * t745;
t839 = t708 * t742;
t696 = t744 * t792 + t839;
t798 = t692 * t740 + t696;
t861 = -t801 * t738 + t739 * t798;
t832 = t742 * t743;
t794 = t846 * t832;
t777 = t847 * t794;
t796 = t747 * t816;
t827 = t744 * t747;
t662 = qJD(1) * t777 - t708 * t813 - t709 * t747 + (t720 * t827 + t796 - t814) * qJD(3);
t860 = t662 * t745;
t859 = t662 * t749;
t848 = r_i_i_C(3) + pkin(13);
t857 = -t749 * r_i_i_C(1) - pkin(5);
t823 = qJD(6) * t745;
t821 = r_i_i_C(1) * t823;
t822 = qJD(6) * t749;
t852 = -t822 * r_i_i_C(2) - t821;
t782 = -t745 * r_i_i_C(2) - t857;
t769 = t780 - t824;
t694 = -t769 * t847 + (-t744 * t768 + t794) * t747;
t765 = t769 * qJD(2);
t755 = t712 * qJD(1) - t742 * t765;
t843 = pkin(4) * qJD(4);
t851 = pkin(4) * t755 - t694 * t843;
t850 = t747 * t769 + t777;
t707 = t721 * qJD(1) + t768 * qJD(2);
t759 = t720 * qJD(1) + t765;
t758 = t759 * t747;
t659 = qJD(1) * t796 + t850 * qJD(3) - t707 * t847 + t744 * t758 - t768 * t791;
t812 = t744 * t846;
t714 = t742 * t768 + t743 * t812;
t834 = t739 * t740;
t648 = t694 * t834 - t755 * t739 + (t714 * t740 + t659) * t738;
t836 = t738 * t740;
t649 = t659 * t739 - t694 * t836 + t714 * t834 + t755 * t738;
t678 = -t694 * t738 + t714 * t739;
t849 = (t648 * t745 - t678 * t822) * r_i_i_C(2) - t678 * t821 + t848 * t649 + t857 * t648;
t746 = sin(qJ(4));
t845 = pkin(4) * t746;
t840 = t707 * t742;
t838 = t709 * t742;
t835 = t738 * t742;
t833 = t740 * t742;
t831 = t742 * t746;
t750 = cos(qJ(4));
t830 = t742 * t750;
t829 = t743 * t751;
t826 = t747 * t748;
t825 = t747 * t751;
t820 = t746 * t843;
t818 = t714 * t843;
t817 = t748 * t832;
t815 = t743 * t826;
t811 = t847 * t748;
t810 = t847 * t751;
t808 = qJD(1) * t828;
t807 = qJD(2) * t832;
t802 = t842 * t742;
t773 = -t744 * t826 + t810;
t781 = t847 * t802;
t795 = t744 * t810;
t685 = qJD(3) * t781 + ((t795 - t826) * qJD(3) + t773 * qJD(2)) * t743;
t719 = -t742 * t829 + t842 * t744;
t799 = t719 * t740 + t685;
t797 = t830 * t843;
t793 = t748 * t807;
t790 = t747 * t802;
t756 = t759 * t847;
t774 = t747 * t768 + t769 * t813;
t670 = t774 * qJD(3) + t707 * t827 + t756;
t788 = -t769 * t833 + t670;
t775 = t720 * t747 - t721 * t813;
t672 = t775 * qJD(3) - t708 * t847 - t709 * t827;
t787 = t721 * t833 + t672;
t786 = t743 * t795;
t702 = -t720 * t847 - t721 * t827;
t784 = t702 * t740 - t838;
t704 = -t768 * t847 + t769 * t827;
t783 = t704 * t740 + t840;
t771 = t744 * t811 + t825;
t772 = t744 * t825 + t811;
t698 = (-t772 * qJD(2) - t771 * qJD(3)) * t743;
t779 = t740 * t817 + t698;
t718 = t773 * t743;
t770 = -t718 * t740 + t751 * t807;
t737 = t750 * pkin(4) + pkin(3);
t762 = -t848 * t738 - t782 * t739 - t737;
t651 = t692 * t836 + t696 * t738 + t864;
t761 = t852 * (t692 * t738 - t712 * t739) + t848 * t651 + t782 * t861;
t711 = t772 * t743 + t790;
t668 = -t711 * t836 + t738 * t793 + t799 * t739;
t760 = t852 * (-t711 * t738 + t719 * t739) + t848 * t668 + t782 * ((-t711 * t740 + t793) * t739 - t799 * t738);
t757 = t820 + (t745 * r_i_i_C(1) + t749 * r_i_i_C(2)) * t739 * qJD(6) + (t782 * t738 - t848 * t739) * t740;
t753 = -pkin(12) - pkin(11);
t717 = t771 * t743;
t710 = -t781 - t786 + t815;
t705 = t718 * t739 + t738 * t817;
t697 = -qJD(2) * t786 - t805 * t829 + (qJD(3) * t744 + qJD(2)) * t815;
t693 = t768 * t813 - t850;
t687 = t711 * t739 + t719 * t738;
t684 = qJD(3) * t790 + (t771 * qJD(2) + t772 * qJD(3)) * t743;
t681 = t704 * t739 - t769 * t835;
t680 = t702 * t739 + t721 * t835;
t679 = t694 * t739 + t714 * t738;
t674 = t770 * t738 + t779 * t739;
t671 = t702 * qJD(3) - t708 * t747 + t709 * t813;
t669 = t704 * qJD(3) - t707 * t813 + t758;
t658 = -qJD(1) * t785 + t694 * qJD(3) - t707 * t747 - t744 * t756;
t657 = -t784 * t738 + t787 * t739;
t655 = -t783 * t738 + t788 * t739;
t653 = -t798 * t738 - t864;
t641 = t649 * t749 + t658 * t745 + (-t679 * t745 + t693 * t749) * qJD(6);
t640 = -t649 * t745 + t658 * t749 + (-t679 * t749 - t693 * t745) * qJD(6);
t1 = [(t653 * t749 + t860) * r_i_i_C(1) + (-t653 * t745 + t859) * r_i_i_C(2) + t653 * pkin(5) - t661 * t737 - t662 * t753 - t709 * pkin(2) - pkin(10) * t839 + t848 * t861 + (t863 * r_i_i_C(1) - t862 * r_i_i_C(2)) * qJD(6) + (-t752 * pkin(1) + (-t846 * pkin(9) - pkin(10) * t812) * t743) * qJD(1) + (-t696 * t746 + (-t692 * t746 + t712 * t750) * qJD(4)) * pkin(4) (t655 * t749 + t669 * t745 + (-t681 * t745 - t749 * t774) * qJD(6)) * r_i_i_C(1) + (-t655 * t745 + t669 * t749 + (-t681 * t749 + t745 * t774) * qJD(6)) * r_i_i_C(2) + t655 * pkin(5) + t670 * t737 - t704 * t820 - t669 * t753 - t707 * pkin(4) * t831 - t769 * t797 + t759 * pkin(2) - pkin(10) * t840 + t848 * (t788 * t738 + t783 * t739) (t659 * t745 + t694 * t822) * r_i_i_C(1) + (t659 * t749 - t694 * t823) * r_i_i_C(2) - t659 * t753 + t762 * t658 + t757 * t693, -t659 * t845 - t746 * t818 + t851 * t750 + t849, t849, t640 * r_i_i_C(1) - t641 * r_i_i_C(2); -pkin(1) * t804 - t707 * pkin(2) + t649 * pkin(5) + pkin(9) * t808 + t641 * r_i_i_C(1) + t640 * r_i_i_C(2) - t658 * t753 + t659 * t737 + t750 * t818 + t851 * t746 + t848 * t648 + (-t742 * t759 + t744 * t808) * pkin(10) (t657 * t749 + t671 * t745) * r_i_i_C(1) + (-t657 * t745 + t671 * t749) * r_i_i_C(2) + t657 * pkin(5) + t672 * t737 - t671 * t753 - t708 * pkin(2) + pkin(10) * t838 + t848 * (t738 * t787 + t739 * t784) + ((-t680 * t745 - t749 * t775) * r_i_i_C(1) + (-t680 * t749 + t745 * t775) * r_i_i_C(2)) * qJD(6) + (t709 * t831 + (-t702 * t746 + t721 * t830) * qJD(4)) * pkin(4) (t661 * t745 - t692 * t822) * r_i_i_C(1) + (t661 * t749 + t692 * t823) * r_i_i_C(2) - t661 * t753 - t762 * t662 + t757 * t689 (-t661 * t746 + t696 * t750 + (t692 * t750 + t712 * t746) * qJD(4)) * pkin(4) + t761, t761 (-t651 * t745 - t859) * r_i_i_C(1) + (-t651 * t749 + t860) * r_i_i_C(2) + (t862 * r_i_i_C(1) + t863 * r_i_i_C(2)) * qJD(6); 0 (t674 * t749 - t697 * t745) * r_i_i_C(1) + (-t674 * t745 - t697 * t749) * r_i_i_C(2) + t674 * pkin(5) + t698 * t737 - t718 * t820 + t697 * t753 + t848 * (t738 * t779 - t739 * t770) + ((-t705 * t745 + t717 * t749) * r_i_i_C(1) + (-t705 * t749 - t717 * t745) * r_i_i_C(2)) * qJD(6) + (t748 * t797 + (-pkin(2) * t748 + (pkin(10) + t845) * t751 * t742) * qJD(2)) * t743 (t685 * t745 + t711 * t822) * r_i_i_C(1) + (t685 * t749 - t711 * t823) * r_i_i_C(2) - t685 * t753 + t762 * t684 + t757 * t710 (t750 * t793 - t685 * t746 + (-t711 * t750 - t719 * t746) * qJD(4)) * pkin(4) + t760, t760 (-t668 * t745 + t684 * t749) * r_i_i_C(1) + (-t668 * t749 - t684 * t745) * r_i_i_C(2) + ((-t687 * t749 - t710 * t745) * r_i_i_C(1) + (t687 * t745 - t710 * t749) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
