% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR8
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:35
% EndTime: 2019-12-31 17:08:36
% DurationCPUTime: 1.05s
% Computational Cost: add. (2033->157), mult. (4440->175), div. (0->0), fcn. (2547->6), ass. (0->106)
t885 = sin(qJ(2));
t888 = cos(qJ(2));
t891 = qJD(1) ^ 2;
t907 = t888 * t891;
t868 = t885 * t907;
t861 = qJDD(2) - t868;
t881 = t885 ^ 2;
t890 = qJD(2) ^ 2;
t865 = t881 * t891 + t890;
t836 = t888 * t861 - t885 * t865;
t901 = qJD(1) * qJD(2);
t897 = t888 * t901;
t900 = t885 * qJDD(1);
t852 = 0.2e1 * t897 + t900;
t886 = sin(qJ(1));
t889 = cos(qJ(1));
t918 = t886 * t836 + t889 * t852;
t917 = t889 * t836 - t886 * t852;
t877 = qJD(2) - qJD(4);
t916 = qJD(4) - t877;
t884 = sin(qJ(4));
t887 = cos(qJ(4));
t844 = (-t884 * t885 - t887 * t888) * qJD(1);
t915 = t844 ^ 2;
t903 = qJD(1) * t888;
t904 = qJD(1) * t885;
t846 = -t884 * t903 + t887 * t904;
t914 = t846 ^ 2;
t913 = t877 ^ 2;
t912 = 2 * qJD(3);
t911 = t888 * g(3);
t910 = t846 * t844;
t882 = t888 ^ 2;
t909 = t882 * t891;
t905 = t881 + t882;
t902 = qJD(4) + t877;
t899 = t888 * qJDD(1);
t898 = -qJDD(2) + qJDD(4);
t873 = t885 * t901;
t863 = t886 * g(1) - t889 * g(2);
t864 = -t889 * g(1) - t886 * g(2);
t848 = -t891 * pkin(1) + qJDD(1) * pkin(5) + t864;
t840 = -t885 * g(3) + t888 * t848;
t853 = t897 + t900;
t854 = -t873 + t899;
t896 = -t884 * t853 - t887 * t854;
t851 = (-pkin(2) * t888 - qJ(3) * t885) * qJD(1);
t895 = qJD(1) * t851 + t848;
t894 = -t887 * t853 + t884 * t854;
t833 = t885 * t861 + t888 * t865;
t893 = -qJDD(2) * pkin(2) - t890 * qJ(3) + qJDD(3) + t911;
t847 = qJDD(1) * pkin(1) + t891 * pkin(5) + t863;
t822 = -t890 * pkin(2) + qJDD(2) * qJ(3) + (qJD(2) * t912) + t851 * t903 + t840;
t892 = t847 + (t854 - t873) * pkin(2);
t866 = -t890 - t909;
t862 = -(qJD(2) * pkin(3)) - pkin(6) * t904;
t860 = qJDD(2) + t868;
t859 = t905 * t891;
t858 = -t886 * qJDD(1) - t889 * t891;
t857 = t889 * qJDD(1) - t886 * t891;
t856 = t905 * qJDD(1);
t855 = -0.2e1 * t873 + t899;
t839 = -t885 * t848 - t911;
t838 = -t913 - t914;
t835 = -t885 * t860 + t888 * t866;
t832 = t888 * t860 + t885 * t866;
t831 = t889 * t856 - t886 * t859;
t830 = t886 * t856 + t889 * t859;
t829 = t898 + t910;
t828 = -t898 + t910;
t827 = -t913 - t915;
t826 = t889 * t835 - t886 * t855;
t825 = t886 * t835 + t889 * t855;
t824 = -t914 - t915;
t823 = t895 * t885 + t893;
t821 = -t885 * t839 + t888 * t840;
t820 = t888 * t839 + t885 * t840;
t819 = t904 * t912 + (t853 + t897) * qJ(3) + t892;
t818 = t887 * t828 - t884 * t838;
t817 = t884 * t828 + t887 * t838;
t816 = -t902 * t844 + t894;
t815 = t916 * t844 - t894;
t814 = -t902 * t846 + t896;
t813 = t916 * t846 - t896;
t812 = t887 * t827 - t884 * t829;
t811 = t884 * t827 + t887 * t829;
t810 = -qJDD(2) * pkin(3) + (-t853 + t897) * pkin(6) + (-pkin(3) * t907 + t895) * t885 + t893;
t809 = -pkin(3) * t909 - t854 * pkin(6) + qJD(2) * t862 + t822;
t808 = t853 * qJ(3) + t854 * pkin(3) - pkin(6) * t909 + (qJ(3) * qJD(2) * t888 + (t912 + t862) * t885) * qJD(1) + t892;
t807 = t888 * t822 + t885 * t823;
t806 = t885 * t822 - t888 * t823;
t805 = t885 * t817 + t888 * t818;
t804 = -t888 * t817 + t885 * t818;
t803 = t887 * t814 - t884 * t816;
t802 = t884 * t814 + t887 * t816;
t801 = t885 * t811 + t888 * t812;
t800 = -t888 * t811 + t885 * t812;
t799 = t887 * t809 + t884 * t810;
t798 = -t884 * t809 + t887 * t810;
t797 = t885 * t802 + t888 * t803;
t796 = -t888 * t802 + t885 * t803;
t795 = -t884 * t798 + t887 * t799;
t794 = t887 * t798 + t884 * t799;
t793 = t885 * t794 + t888 * t795;
t792 = -t888 * t794 + t885 * t795;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t858, -t857, 0, -t886 * t863 + t889 * t864, 0, 0, 0, 0, 0, 0, t826, -t917, t831, t889 * t821 - t886 * t847, 0, 0, 0, 0, 0, 0, t826, t831, t917, t889 * t807 - t886 * t819, 0, 0, 0, 0, 0, 0, t889 * t801 - t886 * t813, t889 * t805 - t886 * t815, t889 * t797 - t886 * t824, t889 * t793 - t886 * t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t857, t858, 0, t889 * t863 + t886 * t864, 0, 0, 0, 0, 0, 0, t825, -t918, t830, t886 * t821 + t889 * t847, 0, 0, 0, 0, 0, 0, t825, t830, t918, t886 * t807 + t889 * t819, 0, 0, 0, 0, 0, 0, t886 * t801 + t889 * t813, t886 * t805 + t889 * t815, t886 * t797 + t889 * t824, t886 * t793 + t889 * t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t832, -t833, 0, t820, 0, 0, 0, 0, 0, 0, t832, 0, t833, t806, 0, 0, 0, 0, 0, 0, t800, t804, t796, t792; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t891, -qJDD(1), 0, t864, 0, 0, 0, 0, 0, 0, t835, -t836, t856, t821, 0, 0, 0, 0, 0, 0, t835, t856, t836, t807, 0, 0, 0, 0, 0, 0, t801, t805, t797, t793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t891, 0, t863, 0, 0, 0, 0, 0, 0, t855, -t852, t859, t847, 0, 0, 0, 0, 0, 0, t855, t859, t852, t819, 0, 0, 0, 0, 0, 0, t813, t815, t824, t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t832, -t833, 0, t820, 0, 0, 0, 0, 0, 0, t832, 0, t833, t806, 0, 0, 0, 0, 0, 0, t800, t804, t796, t792; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t866, -t861, t899, t840, 0, 0, 0, 0, 0, 0, t866, t899, t861, t822, 0, 0, 0, 0, 0, 0, t812, t818, t803, t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t860, -t865, -t900, t839, 0, 0, 0, 0, 0, 0, t860, -t900, t865, -t823, 0, 0, 0, 0, 0, 0, -t811, -t817, -t802, -t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t855, t852, -t859, -t847, 0, 0, 0, 0, 0, 0, -t855, -t859, -t852, -t819, 0, 0, 0, 0, 0, 0, -t813, -t815, -t824, -t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t866, t899, t861, t822, 0, 0, 0, 0, 0, 0, t812, t818, t803, t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t855, -t859, -t852, -t819, 0, 0, 0, 0, 0, 0, -t813, -t815, -t824, -t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t860, t900, -t865, t823, 0, 0, 0, 0, 0, 0, t811, t817, t802, t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t827, t828, t814, t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t829, t838, t816, t798; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t813, t815, t824, t808;];
f_new_reg = t1;
