% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:47
% EndTime: 2019-12-31 16:50:48
% DurationCPUTime: 1.05s
% Computational Cost: add. (2752->144), mult. (5402->200), div. (0->0), fcn. (3434->8), ass. (0->116)
t891 = cos(qJ(3));
t907 = t891 * qJD(1);
t868 = -qJD(4) + t907;
t914 = qJD(4) - t868;
t887 = sin(qJ(4));
t890 = cos(qJ(4));
t888 = sin(qJ(3));
t908 = qJD(1) * t888;
t847 = -t890 * qJD(3) + t887 * t908;
t913 = t847 ^ 2;
t849 = t887 * qJD(3) + t890 * t908;
t912 = t849 ^ 2;
t911 = t868 ^ 2;
t910 = t849 * t847;
t889 = sin(qJ(1));
t892 = cos(qJ(1));
t864 = -t892 * g(1) - t889 * g(2);
t894 = qJD(1) ^ 2;
t850 = -t894 * pkin(1) + t864;
t883 = sin(pkin(7));
t884 = cos(pkin(7));
t863 = t889 * g(1) - t892 * g(2);
t896 = qJDD(1) * pkin(1) + t863;
t834 = t884 * t850 + t883 * t896;
t831 = -t894 * pkin(2) + qJDD(1) * pkin(5) + t834;
t880 = -g(3) + qJDD(2);
t821 = t891 * t831 + t888 * t880;
t878 = t888 ^ 2;
t879 = t891 ^ 2;
t909 = t878 + t879;
t906 = qJD(4) + t868;
t905 = qJD(1) * qJD(3);
t904 = t888 * qJDD(1);
t903 = t888 * t905;
t902 = t891 * t905;
t833 = -t883 * t850 + t884 * t896;
t855 = -t883 * qJDD(1) - t884 * t894;
t856 = t884 * qJDD(1) - t883 * t894;
t901 = t892 * t855 - t889 * t856;
t853 = t902 + t904;
t900 = t890 * qJDD(3) - t887 * t853;
t875 = t891 * qJDD(1);
t899 = -t875 + 0.2e1 * t903;
t898 = t889 * t855 + t892 * t856;
t897 = -t887 * qJDD(3) - t890 * t853;
t830 = -qJDD(1) * pkin(2) - t894 * pkin(5) - t833;
t895 = -qJDD(4) + t875 - t903;
t893 = qJD(3) ^ 2;
t872 = t891 * t880;
t867 = t891 * t894 * t888;
t866 = -t879 * t894 - t893;
t865 = -t878 * t894 - t893;
t862 = -qJDD(3) + t867;
t861 = qJDD(3) + t867;
t860 = t909 * t894;
t859 = -t889 * qJDD(1) - t892 * t894;
t858 = t892 * qJDD(1) - t889 * t894;
t857 = t909 * qJDD(1);
t852 = 0.2e1 * t902 + t904;
t851 = (-pkin(3) * t891 - pkin(6) * t888) * qJD(1);
t841 = t891 * t862 - t888 * t865;
t840 = -t888 * t861 + t891 * t866;
t839 = t888 * t862 + t891 * t865;
t838 = t891 * t861 + t888 * t866;
t837 = -t911 - t912;
t836 = t884 * t857 - t883 * t860;
t835 = t883 * t857 + t884 * t860;
t832 = -t911 - t913;
t829 = -t895 - t910;
t828 = t895 - t910;
t827 = -t912 - t913;
t825 = t884 * t841 + t883 * t852;
t824 = t884 * t840 + t883 * t899;
t823 = t883 * t841 - t884 * t852;
t822 = t883 * t840 - t884 * t899;
t820 = -t888 * t831 + t872;
t819 = t906 * t847 + t897;
t818 = -t914 * t847 - t897;
t817 = -t906 * t849 + t900;
t816 = t914 * t849 - t900;
t815 = -t893 * pkin(3) + qJDD(3) * pkin(6) + t851 * t907 + t821;
t814 = -t872 - qJDD(3) * pkin(3) - t893 * pkin(6) + (qJD(1) * t851 + t831) * t888;
t813 = -t883 * t833 + t884 * t834;
t812 = t884 * t833 + t883 * t834;
t811 = t890 * t828 - t887 * t837;
t810 = t887 * t828 + t890 * t837;
t809 = (-t853 - t902) * pkin(6) + t899 * pkin(3) + t830;
t808 = -t887 * t829 + t890 * t832;
t807 = t890 * t829 + t887 * t832;
t806 = -t888 * t820 + t891 * t821;
t805 = t891 * t820 + t888 * t821;
t804 = t890 * t817 - t887 * t819;
t803 = t887 * t817 + t890 * t819;
t802 = t891 * t811 + t888 * t818;
t801 = t888 * t811 - t891 * t818;
t800 = t891 * t808 + t888 * t816;
t799 = t888 * t808 - t891 * t816;
t798 = t884 * t806 + t883 * t830;
t797 = t883 * t806 - t884 * t830;
t796 = t887 * t809 + t890 * t815;
t795 = t890 * t809 - t887 * t815;
t794 = t891 * t804 + t888 * t827;
t793 = t888 * t804 - t891 * t827;
t792 = t884 * t802 + t883 * t810;
t791 = t883 * t802 - t884 * t810;
t790 = t884 * t800 + t883 * t807;
t789 = t883 * t800 - t884 * t807;
t788 = t884 * t794 + t883 * t803;
t787 = t883 * t794 - t884 * t803;
t786 = -t887 * t795 + t890 * t796;
t785 = t890 * t795 + t887 * t796;
t784 = t891 * t786 + t888 * t814;
t783 = t888 * t786 - t891 * t814;
t782 = t884 * t784 + t883 * t785;
t781 = t883 * t784 - t884 * t785;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t859, -t858, 0, -t889 * t863 + t892 * t864, 0, 0, 0, 0, 0, 0, t901, -t898, 0, -t889 * t812 + t892 * t813, 0, 0, 0, 0, 0, 0, -t889 * t822 + t892 * t824, -t889 * t823 + t892 * t825, -t889 * t835 + t892 * t836, -t889 * t797 + t892 * t798, 0, 0, 0, 0, 0, 0, -t889 * t789 + t892 * t790, -t889 * t791 + t892 * t792, -t889 * t787 + t892 * t788, -t889 * t781 + t892 * t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t858, t859, 0, t892 * t863 + t889 * t864, 0, 0, 0, 0, 0, 0, t898, t901, 0, t892 * t812 + t889 * t813, 0, 0, 0, 0, 0, 0, t892 * t822 + t889 * t824, t892 * t823 + t889 * t825, t892 * t835 + t889 * t836, t892 * t797 + t889 * t798, 0, 0, 0, 0, 0, 0, t892 * t789 + t889 * t790, t892 * t791 + t889 * t792, t892 * t787 + t889 * t788, t892 * t781 + t889 * t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t880, 0, 0, 0, 0, 0, 0, t838, t839, 0, t805, 0, 0, 0, 0, 0, 0, t799, t801, t793, t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t894, -qJDD(1), 0, t864, 0, 0, 0, 0, 0, 0, t855, -t856, 0, t813, 0, 0, 0, 0, 0, 0, t824, t825, t836, t798, 0, 0, 0, 0, 0, 0, t790, t792, t788, t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t894, 0, t863, 0, 0, 0, 0, 0, 0, t856, t855, 0, t812, 0, 0, 0, 0, 0, 0, t822, t823, t835, t797, 0, 0, 0, 0, 0, 0, t789, t791, t787, t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t880, 0, 0, 0, 0, 0, 0, t838, t839, 0, t805, 0, 0, 0, 0, 0, 0, t799, t801, t793, t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t894, -qJDD(1), 0, t834, 0, 0, 0, 0, 0, 0, t840, t841, t857, t806, 0, 0, 0, 0, 0, 0, t800, t802, t794, t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t894, 0, t833, 0, 0, 0, 0, 0, 0, -t899, -t852, t860, -t830, 0, 0, 0, 0, 0, 0, -t807, -t810, -t803, -t785; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t880, 0, 0, 0, 0, 0, 0, t838, t839, 0, t805, 0, 0, 0, 0, 0, 0, t799, t801, t793, t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t866, t862, t875, t821, 0, 0, 0, 0, 0, 0, t808, t811, t804, t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t861, t865, -t904, t820, 0, 0, 0, 0, 0, 0, -t816, -t818, -t827, -t814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t899, t852, -t860, t830, 0, 0, 0, 0, 0, 0, t807, t810, t803, t785; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t832, t828, t817, t796; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t829, t837, t819, t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t816, t818, t827, t814;];
f_new_reg = t1;
