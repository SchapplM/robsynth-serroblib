% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:33
% EndTime: 2019-02-26 22:46:35
% DurationCPUTime: 2.91s
% Computational Cost: add. (3260->241), mult. (9997->391), div. (0->0), fcn. (11050->14), ass. (0->145)
t774 = sin(qJ(2));
t777 = cos(qJ(2));
t778 = cos(qJ(1));
t861 = cos(pkin(6));
t830 = t778 * t861;
t864 = sin(qJ(1));
t752 = t774 * t830 + t864 * t777;
t821 = t861 * t864;
t791 = t778 * t774 + t777 * t821;
t739 = t791 * qJD(1) + t752 * qJD(2);
t804 = t774 * t821;
t836 = t864 * t774;
t848 = t778 * t777;
t740 = -qJD(1) * t804 - qJD(2) * t836 + (qJD(2) * t861 + qJD(1)) * t848;
t751 = -t777 * t830 + t836;
t768 = sin(pkin(7));
t770 = cos(pkin(7));
t773 = sin(qJ(3));
t769 = sin(pkin(6));
t853 = t769 * t778;
t843 = t768 * t853;
t865 = cos(qJ(3));
t816 = t865 * t843;
t833 = qJD(3) * t865;
t823 = t770 * t833;
t832 = t864 * qJD(1);
t824 = t769 * t832;
t697 = (-qJD(3) * t752 - t739 * t770 + t768 * t824) * t773 - qJD(3) * t816 + t740 * t865 - t751 * t823;
t860 = t739 * t768;
t728 = t770 * t824 + t860;
t772 = sin(qJ(4));
t776 = cos(qJ(4));
t841 = t752 * t865;
t724 = (t751 * t770 + t843) * t773 - t841;
t743 = -t751 * t768 + t770 * t853;
t808 = t724 * t772 - t743 * t776;
t681 = -t808 * qJD(4) - t697 * t776 - t728 * t772;
t857 = t768 * t769;
t827 = t864 * t857;
t802 = t865 * t827;
t829 = t773 * t843;
t840 = t770 * t865;
t852 = t770 * t773;
t698 = qJD(1) * t802 - t739 * t840 - t740 * t773 + (t751 * t852 + t829 - t841) * qJD(3);
t771 = sin(qJ(5));
t775 = cos(qJ(5));
t710 = t724 * t776 + t743 * t772;
t721 = t751 * t840 + t752 * t773 + t816;
t881 = t710 * t775 - t721 * t771;
t886 = qJD(5) * t881 + t681 * t771 - t698 * t775;
t882 = t710 * t771 + t721 * t775;
t885 = qJD(5) * t882 - t681 * t775 - t698 * t771;
t880 = t710 * qJD(4) - t697 * t772 + t728 * t776;
t867 = r_i_i_C(1) + pkin(5);
t866 = r_i_i_C(2) + pkin(12);
t862 = r_i_i_C(3) + qJ(6);
t792 = t804 - t848;
t868 = t773 * t792 + t802;
t793 = t862 * t771 + t867 * t775 + pkin(4);
t863 = pkin(10) * t768;
t839 = t770 * t864;
t745 = t768 * t791 + t769 * t839;
t859 = t745 * t776;
t856 = t768 * t772;
t855 = t768 * t776;
t854 = t769 * t777;
t851 = t771 * t776;
t850 = t773 * t774;
t849 = t773 * t777;
t847 = qJD(2) * t769;
t846 = qJD(4) * t772;
t845 = qJD(5) * t776;
t844 = t774 * t857;
t842 = t769 * t850;
t838 = t865 * t774;
t837 = t865 * t777;
t834 = t768 * t847;
t831 = t768 * t861;
t828 = t770 * t837;
t826 = t774 * t834;
t825 = t777 * t834;
t822 = t773 * t831;
t738 = t752 * qJD(1) + t791 * qJD(2);
t786 = t792 * qJD(2);
t782 = t751 * qJD(1) + t786;
t781 = t782 * t773;
t695 = qJD(1) * t829 + qJD(3) * t868 - t738 * t865 + t770 * t781 - t791 * t823;
t725 = t791 * t840 - t868;
t820 = t725 * t845 + t695;
t819 = t721 * t845 + t697;
t797 = -t770 * t850 + t837;
t805 = t865 * t831;
t717 = qJD(3) * t805 + ((t828 - t850) * qJD(3) + t797 * qJD(2)) * t769;
t817 = t769 * t828;
t741 = -t805 - t817 + t842;
t818 = t741 * t845 + t717;
t726 = -t792 * t865 + (-t770 * t791 + t827) * t773;
t712 = t726 * t776 + t745 * t772;
t813 = t712 * t775 + t725 * t771;
t812 = -t712 * t771 + t725 * t775;
t734 = -t751 * t865 - t752 * t852;
t713 = t734 * t776 + t752 * t856;
t799 = t751 * t773 - t752 * t840;
t811 = -t713 * t771 - t775 * t799;
t736 = -t791 * t865 + t792 * t852;
t714 = t736 * t776 - t792 * t856;
t798 = t773 * t791 + t792 * t840;
t810 = -t714 * t771 - t775 * t798;
t796 = t770 * t849 + t838;
t742 = t796 * t769 + t822;
t750 = -t768 * t854 + t861 * t770;
t719 = t742 * t776 + t750 * t772;
t809 = t719 * t775 + t741 * t771;
t749 = t797 * t769;
t737 = t749 * t776 + t772 * t844;
t795 = t770 * t838 + t849;
t748 = t795 * t769;
t807 = -t737 * t771 + t748 * t775;
t806 = -t742 * t772 + t750 * t776;
t801 = -t776 * pkin(4) - t866 * t772 - pkin(3);
t794 = qJD(4) * (pkin(4) * t772 - t866 * t776);
t780 = t782 * t865;
t694 = -qJD(1) * t816 + t726 * qJD(3) - t738 * t773 - t770 * t780;
t790 = qJD(5) * t726 - t694 * t776 + t725 * t846;
t789 = -qJD(5) * t724 + t698 * t776 + t721 * t846;
t716 = qJD(3) * t822 + (t795 * qJD(2) + t796 * qJD(3)) * t769;
t788 = qJD(5) * t742 - t716 * t776 + t741 * t846;
t783 = qJD(6) * t771 + (-t867 * t771 + t862 * t775) * qJD(5);
t779 = t743 * qJD(1) - t768 * t786;
t730 = (-t796 * qJD(2) - t795 * qJD(3)) * t769;
t729 = -qJD(2) * t817 - t833 * t854 + (qJD(3) * t770 + qJD(2)) * t842;
t707 = t772 * t825 + t730 * t776 + (-t749 * t772 + t776 * t844) * qJD(4);
t705 = t799 * qJD(3) - t739 * t865 - t740 * t852;
t704 = t734 * qJD(3) - t739 * t773 + t740 * t840;
t703 = t798 * qJD(3) + t738 * t852 + t780;
t702 = t736 * qJD(3) - t738 * t840 + t781;
t701 = t806 * qJD(4) + t717 * t776 + t772 * t826;
t687 = t740 * t856 + t705 * t776 + (-t734 * t772 + t752 * t855) * qJD(4);
t685 = -t738 * t856 + t703 * t776 + (-t736 * t772 - t792 * t855) * qJD(4);
t682 = t809 * qJD(5) + t701 * t771 - t716 * t775;
t677 = qJD(4) * t859 + t695 * t776 - t726 * t846 + t779 * t772;
t676 = t712 * qJD(4) + t695 * t772 - t779 * t776;
t663 = t812 * qJD(5) + t677 * t775 + t694 * t771;
t662 = t813 * qJD(5) + t677 * t771 - t694 * t775;
t1 = [t882 * qJD(6) + t681 * pkin(4) - t697 * pkin(3) + t698 * pkin(11) - t740 * pkin(2) - pkin(10) * t860 + t866 * t880 - t867 * t885 + t862 * t886 + (-t778 * pkin(1) + (-t864 * pkin(9) - pkin(10) * t839) * t769) * qJD(1), t782 * pkin(2) + t703 * pkin(3) + t685 * pkin(4) + t702 * pkin(11) - t810 * qJD(6) - t738 * t863 + t866 * (t714 * qJD(4) + t703 * t772 + t738 * t855) + t867 * (t810 * qJD(5) + t685 * t775 + t702 * t771) + t862 * (t685 * t771 - t702 * t775 + (t714 * t775 - t771 * t798) * qJD(5)) -(t725 * t851 + t726 * t775) * qJD(6) + t695 * pkin(11) + t867 * (t820 * t771 + t790 * t775) + t862 * (t790 * t771 - t820 * t775) + t725 * t794 + t801 * t694, t866 * t677 + t783 * (-t726 * t772 + t859) - t793 * t676, t813 * qJD(6) - t867 * t662 + t862 * t663, t662; -pkin(1) * t832 - t738 * pkin(2) + t695 * pkin(3) + t677 * pkin(4) + t694 * pkin(11) - t812 * qJD(6) - t782 * t863 + (pkin(10) * t770 + pkin(9)) * qJD(1) * t853 + t866 * t676 + t867 * t663 + t862 * t662, -t811 * qJD(6) + t687 * pkin(4) + t705 * pkin(3) + t704 * pkin(11) - t739 * pkin(2) + t740 * t863 + t866 * (t713 * qJD(4) + t705 * t772 - t740 * t855) + t867 * (t811 * qJD(5) + t687 * t775 + t704 * t771) + t862 * (t687 * t771 - t704 * t775 + (t713 * t775 - t771 * t799) * qJD(5)) -(t721 * t851 - t724 * t775) * qJD(6) + t697 * pkin(11) + t867 * (t819 * t771 + t789 * t775) + t862 * (t789 * t771 - t819 * t775) + t721 * t794 - t801 * t698, -t681 * t866 + t783 * t808 + t793 * t880, -t881 * qJD(6) + t862 * t885 + t867 * t886, -t886; 0, -t807 * qJD(6) + t707 * pkin(4) + t730 * pkin(3) - t729 * pkin(11) + t866 * (t737 * qJD(4) + t730 * t772 - t776 * t825) + t867 * (t807 * qJD(5) + t707 * t775 - t729 * t771) + t862 * (t707 * t771 + t729 * t775 + (t737 * t775 + t748 * t771) * qJD(5)) + (-pkin(2) * t774 + t777 * t863) * t847 -(t741 * t851 + t742 * t775) * qJD(6) + t717 * pkin(11) + t867 * (t818 * t771 + t788 * t775) + t862 * (t788 * t771 - t818 * t775) + t741 * t794 + t801 * t716, t866 * t701 + t783 * t806 + t793 * (-t719 * qJD(4) - t717 * t772 + t776 * t826) t809 * qJD(6) + t862 * (t701 * t775 + t716 * t771 + (-t719 * t771 + t741 * t775) * qJD(5)) - t867 * t682, t682;];
JaD_transl  = t1;
