% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:56
% EndTime: 2019-02-26 22:53:00
% DurationCPUTime: 3.73s
% Computational Cost: add. (3094->275), mult. (9732->476), div. (0->0), fcn. (10685->16), ass. (0->164)
t838 = cos(qJ(1));
t919 = sin(qJ(2));
t920 = sin(qJ(1));
t878 = t920 * t919;
t837 = cos(qJ(2));
t917 = cos(pkin(6));
t882 = t837 * t917;
t814 = -t838 * t882 + t878;
t833 = sin(qJ(3));
t827 = sin(pkin(7));
t828 = sin(pkin(6));
t903 = t828 * t838;
t890 = t827 * t903;
t830 = cos(pkin(7));
t900 = t830 * t833;
t874 = t917 * t919;
t815 = t920 * t837 + t838 * t874;
t836 = cos(qJ(3));
t908 = t815 * t836;
t788 = t814 * t900 + t833 * t890 - t908;
t832 = sin(qJ(4));
t835 = cos(qJ(4));
t910 = t814 * t830;
t859 = t890 + t910;
t787 = t815 * t833 + t859 * t836;
t806 = -t814 * t827 + t830 * t903;
t826 = sin(pkin(8));
t829 = cos(pkin(8));
t870 = t787 * t829 + t806 * t826;
t759 = t788 * t835 + t870 * t832;
t772 = t787 * t826 - t806 * t829;
t831 = sin(qJ(5));
t834 = cos(qJ(5));
t940 = -t759 * t831 - t772 * t834;
t939 = t759 * t834 - t772 * t831;
t856 = t920 * t874;
t896 = t838 * t837;
t803 = -qJD(1) * t856 - qJD(2) * t878 + (qJD(2) * t917 + qJD(1)) * t896;
t840 = t838 * t919 + t920 * t882;
t802 = t840 * qJD(1) + t815 * qJD(2);
t889 = t828 * t920;
t877 = qJD(1) * t889;
t866 = t827 * t877;
t846 = -t802 * t830 + t866;
t753 = (t859 * t833 - t908) * qJD(3) - t803 * t833 + t846 * t836;
t914 = t802 * t827;
t792 = t830 * t877 + t914;
t873 = t753 * t829 + t792 * t826;
t938 = t759 * qJD(4) + t873 * t835;
t930 = -t788 * t832 + t870 * t835;
t937 = -qJD(4) * t930 + t873 * t832;
t743 = t753 * t826 - t792 * t829;
t936 = t743 * t831;
t935 = t743 * t834;
t885 = t919 * t836;
t898 = t833 * t837;
t855 = t830 * t898 + t885;
t881 = t917 * t827;
t875 = t833 * t881;
t805 = t855 * t828 + t875;
t886 = t919 * t833;
t897 = t836 * t837;
t854 = t830 * t897 - t886;
t804 = t854 * t828 + t836 * t881;
t904 = t828 * t837;
t813 = -t827 * t904 + t917 * t830;
t868 = t804 * t829 + t813 * t826;
t771 = t805 * t835 + t868 * t832;
t839 = t856 - t896;
t853 = t827 * t889 - t830 * t840;
t790 = t853 * t833 - t836 * t839;
t921 = r_i_i_C(3) + pkin(13);
t918 = t827 * pkin(11);
t852 = t830 * t886 - t897;
t793 = (-t854 * qJD(2) + t852 * qJD(3)) * t828;
t915 = t793 * t826;
t913 = t802 * t833;
t906 = t826 * t827;
t905 = t827 * t829;
t902 = t829 * t832;
t901 = t829 * t835;
t899 = t830 * t836;
t895 = qJD(2) * t828;
t894 = qJD(3) * t833;
t893 = qJD(3) * t836;
t892 = qJD(5) * t831;
t891 = qJD(5) * t834;
t888 = t828 * t919;
t887 = t830 * t920;
t884 = qJD(1) * t903;
t883 = t837 * t895;
t880 = t827 * t888;
t879 = t827 * t883;
t876 = qJD(2) * t888;
t789 = t833 * t839 + t853 * t836;
t808 = t827 * t840 + t828 * t887;
t869 = t789 * t829 + t808 * t826;
t867 = t827 * t876;
t865 = t834 * r_i_i_C(1) - t831 * r_i_i_C(2) + pkin(4);
t800 = t814 * qJD(1) + t839 * qJD(2);
t801 = t815 * qJD(1) + t840 * qJD(2);
t858 = t836 * t840 - t839 * t900;
t762 = t858 * qJD(3) - t800 * t833 + t801 * t899;
t748 = -t762 * t826 - t801 * t905;
t864 = -t762 * t829 + t801 * t906;
t860 = t814 * t836 + t815 * t900;
t764 = t860 * qJD(3) - t803 * t899 + t913;
t749 = -t764 * t826 + t803 * t905;
t863 = t764 * t829 + t803 * t906;
t768 = -t787 * t835 + t788 * t902;
t769 = t789 * t835 - t790 * t902;
t795 = t814 * t833 - t815 * t899;
t862 = t795 * t829 + t815 * t906;
t797 = t833 * t840 + t839 * t899;
t861 = t797 * t829 - t839 * t906;
t776 = t804 * t835 - t805 * t902;
t857 = qJD(5) * (-t831 * r_i_i_C(1) - t834 * r_i_i_C(2));
t851 = -t830 * t885 - t898;
t850 = -t800 * t827 + t830 * t884;
t849 = t800 * t830 + t827 * t884;
t809 = t851 * t828;
t848 = t809 * t829 + t826 * t880;
t847 = t793 * t829 + t826 * t879;
t845 = t850 * t826;
t781 = -qJD(3) * t875 + (t851 * qJD(2) - t855 * qJD(3)) * t828;
t844 = t781 * t829 + t826 * t867;
t842 = -t790 * t832 + t869 * t835;
t761 = t790 * t835 + t869 * t832;
t841 = -t805 * t832 + t868 * t835;
t766 = t862 * t832 - t835 * t860;
t767 = t861 * t832 - t835 * t858;
t751 = -qJD(3) * t790 + t801 * t833 + t849 * t836;
t741 = -t751 * t826 + t850 * t829;
t810 = t852 * t828;
t777 = -t810 * t835 + t848 * t832;
t819 = t890 * t893;
t799 = -t809 * t826 + t829 * t880;
t794 = (-t855 * qJD(2) + t851 * qJD(3)) * t828;
t784 = -t804 * t826 + t813 * t829;
t782 = -t876 * t900 - t828 * qJD(3) * t886 + (t883 + (t830 * t904 + t881) * qJD(3)) * t836;
t780 = t829 * t879 - t915;
t779 = -t797 * t826 - t839 * t905;
t778 = -t795 * t826 + t815 * t905;
t775 = -t781 * t826 + t829 * t867;
t774 = -t789 * t826 + t808 * t829;
t765 = t795 * qJD(3) - t802 * t836 - t803 * t900;
t763 = t797 * qJD(3) + t800 * t836 + t801 * t900;
t756 = t819 + (qJD(3) * t910 - t803) * t836 + (qJD(3) * t815 - t846) * t833;
t754 = t833 * t866 + t803 * t836 - t815 * t894 - t819 + (-t814 * t893 - t913) * t830;
t752 = t839 * t894 + t849 * t833 + (t853 * qJD(3) - t801) * t836;
t747 = t794 * t835 + t847 * t832 + (t810 * t832 + t848 * t835) * qJD(4);
t745 = -t782 * t902 + t781 * t835 + (-t804 * t832 - t805 * t901) * qJD(4);
t740 = t841 * qJD(4) + t782 * t835 + t844 * t832;
t738 = t765 * t835 + t863 * t832 + (t832 * t860 + t862 * t835) * qJD(4);
t736 = t763 * t835 - t864 * t832 + (t832 * t858 + t861 * t835) * qJD(4);
t734 = -t754 * t902 + t753 * t835 + (t787 * t832 + t788 * t901) * qJD(4);
t732 = -t752 * t902 + t751 * t835 + (-t789 * t832 - t790 * t901) * qJD(4);
t730 = t756 * t835 - t937;
t728 = t754 * t835 + t937;
t726 = t752 * t835 + (t751 * t829 + t845) * t832 + t842 * qJD(4);
t725 = t761 * qJD(4) - t751 * t901 + t752 * t832 - t835 * t845;
t724 = t726 * t834 + t741 * t831 + (-t761 * t831 + t774 * t834) * qJD(5);
t723 = -t726 * t831 + t741 * t834 + (-t761 * t834 - t774 * t831) * qJD(5);
t1 = [(t730 * t834 + t936) * r_i_i_C(1) + (-t730 * t831 + t935) * r_i_i_C(2) + t730 * pkin(4) + t756 * pkin(3) - t803 * pkin(2) - pkin(11) * t914 + t921 * (t756 * t832 + t938) + (r_i_i_C(1) * t940 - t939 * r_i_i_C(2)) * qJD(5) + t743 * pkin(12) + (-t838 * pkin(1) + (-t920 * pkin(10) - pkin(11) * t887) * t828) * qJD(1) (t736 * t834 + t748 * t831) * r_i_i_C(1) + (-t736 * t831 + t748 * t834) * r_i_i_C(2) + t736 * pkin(4) + t763 * pkin(3) + t800 * pkin(2) - t801 * t918 + t921 * (t767 * qJD(4) + t763 * t832 + t864 * t835) + ((-t767 * t831 + t779 * t834) * r_i_i_C(1) + (-t767 * t834 - t779 * t831) * r_i_i_C(2)) * qJD(5) + t748 * pkin(12) (t732 * t834 - t769 * t892) * r_i_i_C(1) + (-t732 * t831 - t769 * t891) * r_i_i_C(2) + t732 * pkin(4) + t751 * pkin(3) + t921 * (qJD(4) * t769 + t751 * t832 + t752 * t901) + ((t752 * t831 + t790 * t891) * r_i_i_C(1) + (t752 * t834 - t790 * t892) * r_i_i_C(2) + t752 * pkin(12)) * t826, -t865 * t725 + t921 * t726 + t842 * t857, t723 * r_i_i_C(1) - t724 * r_i_i_C(2), 0; -t801 * pkin(2) + t752 * pkin(3) + t726 * pkin(4) + t724 * r_i_i_C(1) + t723 * r_i_i_C(2) + t921 * t725 + (-t920 * pkin(1) + pkin(10) * t903) * qJD(1) + t741 * pkin(12) + t850 * pkin(11) (t738 * t834 + t749 * t831) * r_i_i_C(1) + (-t738 * t831 + t749 * t834) * r_i_i_C(2) + t738 * pkin(4) + t765 * pkin(3) - t802 * pkin(2) + t803 * t918 + t921 * (t766 * qJD(4) + t765 * t832 - t863 * t835) + ((-t766 * t831 + t778 * t834) * r_i_i_C(1) + (-t766 * t834 - t778 * t831) * r_i_i_C(2)) * qJD(5) + t749 * pkin(12) (t734 * t834 - t768 * t892) * r_i_i_C(1) + (-t734 * t831 - t768 * t891) * r_i_i_C(2) + t734 * pkin(4) + t753 * pkin(3) + t921 * (qJD(4) * t768 + t753 * t832 + t754 * t901) + ((t754 * t831 - t788 * t891) * r_i_i_C(1) + (t754 * t834 + t788 * t892) * r_i_i_C(2) + t754 * pkin(12)) * t826, t921 * t728 - t930 * t857 + t865 * (-t754 * t832 + t938) (-t728 * t831 - t935) * r_i_i_C(1) + (-t728 * t834 + t936) * r_i_i_C(2) + (t939 * r_i_i_C(1) + r_i_i_C(2) * t940) * qJD(5), 0; 0 (t747 * t834 + t780 * t831) * r_i_i_C(1) + (-t747 * t831 + t780 * t834) * r_i_i_C(2) + t747 * pkin(4) + t794 * pkin(3) - pkin(12) * t915 + t921 * (qJD(4) * t777 + t794 * t832 - t835 * t847) + ((-t777 * t831 + t799 * t834) * r_i_i_C(1) + (-t777 * t834 - t799 * t831) * r_i_i_C(2)) * qJD(5) + (-t919 * pkin(2) + (pkin(12) * t829 + pkin(11)) * t837 * t827) * t895 (t745 * t834 - t776 * t892) * r_i_i_C(1) + (-t745 * t831 - t776 * t891) * r_i_i_C(2) + t745 * pkin(4) + t781 * pkin(3) + t921 * (qJD(4) * t776 + t781 * t832 + t782 * t901) + ((t782 * t831 + t805 * t891) * r_i_i_C(1) + (t782 * t834 - t805 * t892) * r_i_i_C(2) + t782 * pkin(12)) * t826, t921 * t740 + t841 * t857 + t865 * (-qJD(4) * t771 - t782 * t832 + t844 * t835) (-t740 * t831 + t775 * t834) * r_i_i_C(1) + (-t740 * t834 - t775 * t831) * r_i_i_C(2) + ((-t771 * t834 - t784 * t831) * r_i_i_C(1) + (t771 * t831 - t784 * t834) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
