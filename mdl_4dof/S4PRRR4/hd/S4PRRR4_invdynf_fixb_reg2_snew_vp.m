% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:50
% EndTime: 2019-12-31 16:32:51
% DurationCPUTime: 0.96s
% Computational Cost: add. (2733->130), mult. (5682->193), div. (0->0), fcn. (3935->8), ass. (0->112)
t853 = qJD(3) + qJD(4);
t888 = qJD(4) + t853;
t861 = sin(qJ(4));
t864 = cos(qJ(4));
t865 = cos(qJ(3));
t862 = sin(qJ(3));
t881 = qJD(2) * t862;
t823 = -t864 * t865 * qJD(2) + t861 * t881;
t887 = t823 ^ 2;
t825 = (t861 * t865 + t862 * t864) * qJD(2);
t886 = t825 ^ 2;
t885 = t853 ^ 2;
t884 = t825 * t823;
t855 = t865 ^ 2;
t868 = qJD(2) ^ 2;
t883 = t855 * t868;
t858 = sin(pkin(7));
t859 = cos(pkin(7));
t839 = -t859 * g(1) - t858 * g(2);
t863 = sin(qJ(2));
t866 = cos(qJ(2));
t872 = t858 * g(1) - t859 * g(2);
t816 = t866 * t839 + t863 * t872;
t812 = -t868 * pkin(2) + qJDD(2) * pkin(5) + t816;
t856 = -g(3) + qJDD(1);
t804 = t865 * t812 + t862 * t856;
t854 = t862 ^ 2;
t882 = t854 + t855;
t880 = qJD(4) - t853;
t879 = qJD(2) * qJD(3);
t878 = t862 * qJDD(2);
t845 = t862 * t868 * t865;
t840 = qJDD(3) + t845;
t877 = -qJDD(3) - qJDD(4);
t876 = t862 * t879;
t875 = t865 * t879;
t801 = -t862 * t812 + t865 * t856;
t833 = t875 + t878;
t851 = t865 * qJDD(2);
t871 = -t851 + t876;
t874 = -t861 * t833 - t864 * t871;
t836 = t866 * qJDD(2) - t863 * t868;
t837 = -t863 * qJDD(2) - t866 * t868;
t873 = -t858 * t836 + t859 * t837;
t815 = -t863 * t839 + t866 * t872;
t870 = t859 * t836 + t858 * t837;
t811 = -qJDD(2) * pkin(2) - t868 * pkin(5) - t815;
t869 = -t864 * t833 + t861 * t871;
t867 = qJD(3) ^ 2;
t844 = -t867 - t883;
t843 = -t854 * t868 - t867;
t842 = qJD(3) * pkin(3) - pkin(6) * t881;
t841 = -qJDD(3) + t845;
t838 = t882 * t868;
t835 = t882 * qJDD(2);
t834 = t851 - 0.2e1 * t876;
t832 = 0.2e1 * t875 + t878;
t821 = -t885 - t886;
t820 = t865 * t841 - t862 * t843;
t819 = -t862 * t840 + t865 * t844;
t818 = t862 * t841 + t865 * t843;
t817 = t865 * t840 + t862 * t844;
t814 = t866 * t835 - t863 * t838;
t813 = t863 * t835 + t866 * t838;
t810 = t877 - t884;
t809 = -t877 - t884;
t807 = -t885 - t887;
t806 = t866 * t820 + t863 * t832;
t805 = t866 * t819 - t863 * t834;
t803 = t863 * t820 - t866 * t832;
t802 = t863 * t819 + t866 * t834;
t800 = -t886 - t887;
t799 = -t863 * t815 + t866 * t816;
t798 = t866 * t815 + t863 * t816;
t797 = t871 * pkin(3) - pkin(6) * t883 + t842 * t881 + t811;
t796 = t864 * t810 - t861 * t821;
t795 = t861 * t810 + t864 * t821;
t794 = t880 * t823 + t869;
t793 = -t888 * t823 - t869;
t792 = -t880 * t825 + t874;
t791 = t888 * t825 - t874;
t790 = -pkin(3) * t883 - t871 * pkin(6) - qJD(3) * t842 + t804;
t789 = (-t833 + t875) * pkin(6) + t840 * pkin(3) + t801;
t788 = t864 * t807 - t861 * t809;
t787 = t861 * t807 + t864 * t809;
t786 = -t862 * t801 + t865 * t804;
t785 = t865 * t801 + t862 * t804;
t784 = t866 * t786 + t863 * t811;
t783 = t863 * t786 - t866 * t811;
t782 = -t862 * t795 + t865 * t796;
t781 = t865 * t795 + t862 * t796;
t780 = t864 * t792 - t861 * t794;
t779 = t861 * t792 + t864 * t794;
t778 = t861 * t789 + t864 * t790;
t777 = t864 * t789 - t861 * t790;
t776 = -t862 * t787 + t865 * t788;
t775 = t865 * t787 + t862 * t788;
t774 = t866 * t782 + t863 * t793;
t773 = t863 * t782 - t866 * t793;
t772 = t866 * t776 + t863 * t791;
t771 = t863 * t776 - t866 * t791;
t770 = -t862 * t779 + t865 * t780;
t769 = t865 * t779 + t862 * t780;
t768 = -t861 * t777 + t864 * t778;
t767 = t864 * t777 + t861 * t778;
t766 = t866 * t770 + t863 * t800;
t765 = t863 * t770 - t866 * t800;
t764 = -t862 * t767 + t865 * t768;
t763 = t865 * t767 + t862 * t768;
t762 = t866 * t764 + t863 * t797;
t761 = t863 * t764 - t866 * t797;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t859 * t839 - t858 * t872, 0, 0, 0, 0, 0, 0, t873, -t870, 0, -t858 * t798 + t859 * t799, 0, 0, 0, 0, 0, 0, -t858 * t802 + t859 * t805, -t858 * t803 + t859 * t806, -t858 * t813 + t859 * t814, -t858 * t783 + t859 * t784, 0, 0, 0, 0, 0, 0, -t858 * t771 + t859 * t772, -t858 * t773 + t859 * t774, -t858 * t765 + t859 * t766, -t858 * t761 + t859 * t762; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t858 * t839 + t859 * t872, 0, 0, 0, 0, 0, 0, t870, t873, 0, t859 * t798 + t858 * t799, 0, 0, 0, 0, 0, 0, t859 * t802 + t858 * t805, t859 * t803 + t858 * t806, t859 * t813 + t858 * t814, t859 * t783 + t858 * t784, 0, 0, 0, 0, 0, 0, t859 * t771 + t858 * t772, t859 * t773 + t858 * t774, t859 * t765 + t858 * t766, t859 * t761 + t858 * t762; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, 0, 0, 0, 0, 0, 0, t817, t818, 0, t785, 0, 0, 0, 0, 0, 0, t775, t781, t769, t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t839, 0, 0, 0, 0, 0, 0, t837, -t836, 0, t799, 0, 0, 0, 0, 0, 0, t805, t806, t814, t784, 0, 0, 0, 0, 0, 0, t772, t774, t766, t762; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t872, 0, 0, 0, 0, 0, 0, t836, t837, 0, t798, 0, 0, 0, 0, 0, 0, t802, t803, t813, t783, 0, 0, 0, 0, 0, 0, t771, t773, t765, t761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, 0, 0, 0, 0, 0, 0, t817, t818, 0, t785, 0, 0, 0, 0, 0, 0, t775, t781, t769, t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t868, -qJDD(2), 0, t816, 0, 0, 0, 0, 0, 0, t819, t820, t835, t786, 0, 0, 0, 0, 0, 0, t776, t782, t770, t764; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t868, 0, t815, 0, 0, 0, 0, 0, 0, t834, -t832, t838, -t811, 0, 0, 0, 0, 0, 0, -t791, -t793, -t800, -t797; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, 0, 0, 0, 0, 0, 0, t817, t818, 0, t785, 0, 0, 0, 0, 0, 0, t775, t781, t769, t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t844, t841, t851, t804, 0, 0, 0, 0, 0, 0, t788, t796, t780, t768; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t840, t843, -t878, t801, 0, 0, 0, 0, 0, 0, t787, t795, t779, t767; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t834, t832, -t838, t811, 0, 0, 0, 0, 0, 0, t791, t793, t800, t797; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t807, t810, t792, t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t809, t821, t794, t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t791, t793, t800, t797;];
f_new_reg = t1;
