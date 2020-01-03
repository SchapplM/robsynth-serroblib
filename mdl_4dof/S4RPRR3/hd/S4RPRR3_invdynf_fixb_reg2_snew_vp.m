% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:30
% EndTime: 2019-12-31 16:49:31
% DurationCPUTime: 1.05s
% Computational Cost: add. (3030->140), mult. (6131->199), div. (0->0), fcn. (3943->8), ass. (0->116)
t870 = qJD(3) + qJD(4);
t905 = qJD(4) + t870;
t878 = sin(qJ(4));
t881 = cos(qJ(4));
t882 = cos(qJ(3));
t879 = sin(qJ(3));
t898 = qJD(1) * t879;
t839 = -t881 * t882 * qJD(1) + t878 * t898;
t904 = t839 ^ 2;
t841 = (t878 * t882 + t879 * t881) * qJD(1);
t903 = t841 ^ 2;
t902 = t870 ^ 2;
t901 = t841 * t839;
t872 = t882 ^ 2;
t885 = qJD(1) ^ 2;
t900 = t872 * t885;
t880 = sin(qJ(1));
t883 = cos(qJ(1));
t859 = -t883 * g(1) - t880 * g(2);
t845 = -t885 * pkin(1) + t859;
t875 = sin(pkin(7));
t876 = cos(pkin(7));
t858 = t880 * g(1) - t883 * g(2);
t887 = qJDD(1) * pkin(1) + t858;
t827 = t876 * t845 + t875 * t887;
t822 = -t885 * pkin(2) + qJDD(1) * pkin(5) + t827;
t873 = -g(3) + qJDD(2);
t814 = t882 * t822 + t879 * t873;
t871 = t879 ^ 2;
t899 = t871 + t872;
t897 = qJD(4) - t870;
t896 = qJD(1) * qJD(3);
t895 = t879 * qJDD(1);
t863 = t882 * t885 * t879;
t855 = qJDD(3) + t863;
t894 = -qJDD(3) - qJDD(4);
t893 = t879 * t896;
t892 = t882 * t896;
t813 = -t879 * t822 + t882 * t873;
t826 = -t875 * t845 + t876 * t887;
t847 = t892 + t895;
t869 = t882 * qJDD(1);
t889 = -t869 + t893;
t891 = -t878 * t847 - t881 * t889;
t849 = -t875 * qJDD(1) - t876 * t885;
t850 = t876 * qJDD(1) - t875 * t885;
t890 = t883 * t849 - t880 * t850;
t888 = t880 * t849 + t883 * t850;
t821 = -qJDD(1) * pkin(2) - t885 * pkin(5) - t826;
t886 = -t881 * t847 + t878 * t889;
t884 = qJD(3) ^ 2;
t861 = -t884 - t900;
t860 = -t871 * t885 - t884;
t857 = qJD(3) * pkin(3) - pkin(6) * t898;
t856 = -qJDD(3) + t863;
t854 = t899 * t885;
t853 = -t880 * qJDD(1) - t883 * t885;
t852 = t883 * qJDD(1) - t880 * t885;
t851 = t899 * qJDD(1);
t848 = t869 - 0.2e1 * t893;
t846 = 0.2e1 * t892 + t895;
t834 = -t902 - t903;
t833 = t882 * t856 - t879 * t860;
t832 = -t879 * t855 + t882 * t861;
t831 = t879 * t856 + t882 * t860;
t830 = t882 * t855 + t879 * t861;
t829 = t876 * t851 - t875 * t854;
t828 = t875 * t851 + t876 * t854;
t825 = t894 - t901;
t824 = -t894 - t901;
t823 = -t902 - t904;
t819 = t876 * t833 + t875 * t846;
t818 = t876 * t832 - t875 * t848;
t817 = t875 * t833 - t876 * t846;
t816 = t875 * t832 + t876 * t848;
t815 = -t903 - t904;
t812 = t881 * t825 - t878 * t834;
t811 = t878 * t825 + t881 * t834;
t810 = t897 * t839 + t886;
t809 = -t905 * t839 - t886;
t808 = -t897 * t841 + t891;
t807 = t905 * t841 - t891;
t806 = t889 * pkin(3) - pkin(6) * t900 + t857 * t898 + t821;
t805 = -t875 * t826 + t876 * t827;
t804 = t876 * t826 + t875 * t827;
t803 = t881 * t823 - t878 * t824;
t802 = t878 * t823 + t881 * t824;
t801 = -pkin(3) * t900 - t889 * pkin(6) - qJD(3) * t857 + t814;
t800 = (-t847 + t892) * pkin(6) + t855 * pkin(3) + t813;
t799 = -t879 * t813 + t882 * t814;
t798 = t882 * t813 + t879 * t814;
t797 = -t879 * t811 + t882 * t812;
t796 = t882 * t811 + t879 * t812;
t795 = t881 * t808 - t878 * t810;
t794 = t878 * t808 + t881 * t810;
t793 = t876 * t799 + t875 * t821;
t792 = t875 * t799 - t876 * t821;
t791 = -t879 * t802 + t882 * t803;
t790 = t882 * t802 + t879 * t803;
t789 = t878 * t800 + t881 * t801;
t788 = t881 * t800 - t878 * t801;
t787 = t876 * t797 + t875 * t809;
t786 = t875 * t797 - t876 * t809;
t785 = t876 * t791 + t875 * t807;
t784 = t875 * t791 - t876 * t807;
t783 = -t879 * t794 + t882 * t795;
t782 = t882 * t794 + t879 * t795;
t781 = -t878 * t788 + t881 * t789;
t780 = t881 * t788 + t878 * t789;
t779 = t876 * t783 + t875 * t815;
t778 = t875 * t783 - t876 * t815;
t777 = -t879 * t780 + t882 * t781;
t776 = t882 * t780 + t879 * t781;
t775 = t876 * t777 + t875 * t806;
t774 = t875 * t777 - t876 * t806;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t853, -t852, 0, -t880 * t858 + t883 * t859, 0, 0, 0, 0, 0, 0, t890, -t888, 0, -t880 * t804 + t883 * t805, 0, 0, 0, 0, 0, 0, -t880 * t816 + t883 * t818, -t880 * t817 + t883 * t819, -t880 * t828 + t883 * t829, -t880 * t792 + t883 * t793, 0, 0, 0, 0, 0, 0, -t880 * t784 + t883 * t785, -t880 * t786 + t883 * t787, -t880 * t778 + t883 * t779, -t880 * t774 + t883 * t775; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t852, t853, 0, t883 * t858 + t880 * t859, 0, 0, 0, 0, 0, 0, t888, t890, 0, t883 * t804 + t880 * t805, 0, 0, 0, 0, 0, 0, t883 * t816 + t880 * t818, t883 * t817 + t880 * t819, t883 * t828 + t880 * t829, t883 * t792 + t880 * t793, 0, 0, 0, 0, 0, 0, t883 * t784 + t880 * t785, t883 * t786 + t880 * t787, t883 * t778 + t880 * t779, t883 * t774 + t880 * t775; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t873, 0, 0, 0, 0, 0, 0, t830, t831, 0, t798, 0, 0, 0, 0, 0, 0, t790, t796, t782, t776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t885, -qJDD(1), 0, t859, 0, 0, 0, 0, 0, 0, t849, -t850, 0, t805, 0, 0, 0, 0, 0, 0, t818, t819, t829, t793, 0, 0, 0, 0, 0, 0, t785, t787, t779, t775; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t885, 0, t858, 0, 0, 0, 0, 0, 0, t850, t849, 0, t804, 0, 0, 0, 0, 0, 0, t816, t817, t828, t792, 0, 0, 0, 0, 0, 0, t784, t786, t778, t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t873, 0, 0, 0, 0, 0, 0, t830, t831, 0, t798, 0, 0, 0, 0, 0, 0, t790, t796, t782, t776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t885, -qJDD(1), 0, t827, 0, 0, 0, 0, 0, 0, t832, t833, t851, t799, 0, 0, 0, 0, 0, 0, t791, t797, t783, t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t885, 0, t826, 0, 0, 0, 0, 0, 0, t848, -t846, t854, -t821, 0, 0, 0, 0, 0, 0, -t807, -t809, -t815, -t806; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t873, 0, 0, 0, 0, 0, 0, t830, t831, 0, t798, 0, 0, 0, 0, 0, 0, t790, t796, t782, t776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t861, t856, t869, t814, 0, 0, 0, 0, 0, 0, t803, t812, t795, t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t855, t860, -t895, t813, 0, 0, 0, 0, 0, 0, t802, t811, t794, t780; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t848, t846, -t854, t821, 0, 0, 0, 0, 0, 0, t807, t809, t815, t806; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, t825, t808, t789; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t824, t834, t810, t788; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t807, t809, t815, t806;];
f_new_reg = t1;
