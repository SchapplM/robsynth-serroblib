% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRRP8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:44
% EndTime: 2019-12-31 18:47:45
% DurationCPUTime: 1.43s
% Computational Cost: add. (2612->166), mult. (3759->142), div. (0->0), fcn. (1635->6), ass. (0->80)
t883 = -qJD(1) + qJD(3);
t881 = t883 ^ 2;
t889 = sin(qJ(4));
t892 = cos(qJ(4));
t869 = t892 * t881 * t889;
t863 = qJDD(4) - t869;
t886 = t889 ^ 2;
t895 = qJD(4) ^ 2;
t866 = t886 * t881 + t895;
t843 = t892 * t863 - t889 * t866;
t904 = qJD(4) * t883;
t882 = qJDD(1) - qJDD(3);
t909 = t889 * t882;
t854 = 0.2e1 * t892 * t904 - t909;
t890 = sin(qJ(3));
t893 = cos(qJ(3));
t830 = t890 * t843 + t893 * t854;
t833 = t893 * t843 - t890 * t854;
t891 = sin(qJ(1));
t894 = cos(qJ(1));
t915 = t894 * t830 - t891 * t833;
t914 = t891 * t830 + t894 * t833;
t900 = t890 * t881 + t893 * t882;
t901 = -t893 * t881 + t890 * t882;
t913 = t891 * t901 + t894 * t900;
t912 = t891 * t900 - t894 * t901;
t911 = -pkin(1) - pkin(2);
t910 = (-pkin(4) * t892 - qJ(5) * t889) * t881;
t907 = t892 * t882;
t896 = qJD(1) ^ 2;
t871 = -t894 * g(1) - t891 * g(2);
t899 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t871;
t848 = t911 * t896 + t899;
t870 = t891 * g(1) - t894 * g(2);
t898 = -t896 * qJ(2) + qJDD(2) - t870;
t897 = t911 * qJDD(1) + t898;
t836 = t893 * t848 + t890 * t897;
t828 = -t881 * pkin(3) - t882 * pkin(7) + t836;
t825 = t889 * g(3) + t892 * t828;
t887 = t892 ^ 2;
t905 = t886 + t887;
t903 = t889 * t904;
t835 = -t890 * t848 + t893 * t897;
t840 = t889 * t863 + t892 * t866;
t827 = t882 * pkin(3) - t881 * pkin(7) - t835;
t880 = t892 * g(3);
t867 = -t887 * t881 - t895;
t865 = t894 * qJDD(1) - t891 * t896;
t864 = t891 * qJDD(1) + t894 * t896;
t862 = qJDD(4) + t869;
t861 = t905 * t881;
t856 = t905 * t882;
t855 = -0.2e1 * t903 - t907;
t852 = qJDD(1) * pkin(1) - t898;
t849 = -t896 * pkin(1) + t899;
t842 = -t889 * t862 + t892 * t867;
t839 = t892 * t862 + t889 * t867;
t838 = -t893 * t856 - t890 * t861;
t837 = -t890 * t856 + t893 * t861;
t832 = t893 * t842 - t890 * t855;
t829 = t890 * t842 + t893 * t855;
t824 = -t889 * t828 + t880;
t823 = t891 * t837 + t894 * t838;
t822 = -t894 * t837 + t891 * t838;
t821 = qJDD(5) - t880 - t895 * qJ(5) - qJDD(4) * pkin(4) + (t828 + t910) * t889;
t820 = -t895 * pkin(4) + qJDD(4) * qJ(5) + (2 * qJD(5) * qJD(4)) + t892 * t910 + t825;
t819 = -t890 * t835 + t893 * t836;
t818 = t893 * t835 + t890 * t836;
t817 = t891 * t829 + t894 * t832;
t816 = -t894 * t829 + t891 * t832;
t815 = -(-t903 - t907) * pkin(4) + (pkin(4) * qJD(4) - (2 * qJD(5))) * t889 * t883 + t827 - t854 * qJ(5);
t814 = -t889 * t824 + t892 * t825;
t813 = t892 * t824 + t889 * t825;
t812 = t892 * t820 + t889 * t821;
t811 = t889 * t820 - t892 * t821;
t810 = t893 * t814 + t890 * t827;
t809 = t890 * t814 - t893 * t827;
t808 = t893 * t812 + t890 * t815;
t807 = t890 * t812 - t893 * t815;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t864, -t865, 0, -t891 * t870 + t894 * t871, 0, 0, 0, 0, 0, 0, -t864, 0, t865, t894 * t849 - t891 * t852, 0, 0, 0, 0, 0, 0, -t912, t913, 0, t891 * t818 + t894 * t819, 0, 0, 0, 0, 0, 0, t817, -t914, t823, t891 * t809 + t894 * t810, 0, 0, 0, 0, 0, 0, t817, t823, t914, t891 * t807 + t894 * t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t865, -t864, 0, t894 * t870 + t891 * t871, 0, 0, 0, 0, 0, 0, t865, 0, t864, t891 * t849 + t894 * t852, 0, 0, 0, 0, 0, 0, t913, t912, 0, -t894 * t818 + t891 * t819, 0, 0, 0, 0, 0, 0, t816, t915, t822, -t894 * t809 + t891 * t810, 0, 0, 0, 0, 0, 0, t816, t822, -t915, -t894 * t807 + t891 * t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t839, t840, 0, -t813, 0, 0, 0, 0, 0, 0, -t839, 0, -t840, -t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t896, -qJDD(1), 0, t871, 0, 0, 0, 0, 0, 0, -t896, 0, qJDD(1), t849, 0, 0, 0, 0, 0, 0, t901, t900, 0, t819, 0, 0, 0, 0, 0, 0, t832, -t833, t838, t810, 0, 0, 0, 0, 0, 0, t832, t838, t833, t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t896, 0, t870, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t896, t852, 0, 0, 0, 0, 0, 0, t900, -t901, 0, -t818, 0, 0, 0, 0, 0, 0, -t829, t830, -t837, -t809, 0, 0, 0, 0, 0, 0, -t829, -t837, -t830, -t807; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t839, t840, 0, -t813, 0, 0, 0, 0, 0, 0, -t839, 0, -t840, -t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t896, 0, qJDD(1), t849, 0, 0, 0, 0, 0, 0, t901, t900, 0, t819, 0, 0, 0, 0, 0, 0, t832, -t833, t838, t810, 0, 0, 0, 0, 0, 0, t832, t838, t833, t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t839, t840, 0, -t813, 0, 0, 0, 0, 0, 0, -t839, 0, -t840, -t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t896, -t852, 0, 0, 0, 0, 0, 0, -t900, t901, 0, t818, 0, 0, 0, 0, 0, 0, t829, -t830, t837, t809, 0, 0, 0, 0, 0, 0, t829, t837, t830, t807; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t881, t882, 0, t836, 0, 0, 0, 0, 0, 0, t842, -t843, -t856, t814, 0, 0, 0, 0, 0, 0, t842, -t856, t843, t812; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t882, -t881, 0, t835, 0, 0, 0, 0, 0, 0, t855, -t854, t861, -t827, 0, 0, 0, 0, 0, 0, t855, t861, t854, -t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, t839, -t840, 0, t813, 0, 0, 0, 0, 0, 0, t839, 0, t840, t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t867, -t863, -t907, t825, 0, 0, 0, 0, 0, 0, t867, -t907, t863, t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t862, -t866, t909, t824, 0, 0, 0, 0, 0, 0, t862, t909, t866, -t821; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t855, t854, -t861, t827, 0, 0, 0, 0, 0, 0, -t855, -t861, -t854, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t867, -t907, t863, t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t855, -t861, -t854, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t862, -t909, -t866, t821;];
f_new_reg = t1;
