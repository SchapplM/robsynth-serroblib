% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:07:43
% EndTime: 2019-05-04 19:07:44
% DurationCPUTime: 1.41s
% Computational Cost: add. (1523->162), mult. (5086->175), div. (0->0), fcn. (3166->6), ass. (0->100)
t937 = sin(pkin(6));
t938 = sin(pkin(4));
t940 = cos(pkin(4));
t986 = qJD(1) ^ 2;
t971 = t940 * t986;
t963 = t938 * t971;
t918 = t937 * t963;
t939 = cos(pkin(6));
t967 = qJDD(1) * t938;
t962 = t939 * t967;
t906 = -t918 + t962;
t972 = t986 * t938 ^ 2;
t957 = t939 * t937 * t972;
t966 = t940 * qJDD(1);
t909 = t957 + t966;
t987 = t939 ^ 2;
t927 = t987 * t972;
t935 = t940 ^ 2 * t986;
t914 = -t927 - t935;
t950 = t909 * t939 + t914 * t937;
t890 = t938 * t906 + t950 * t940;
t898 = t937 * t909 - t939 * t914;
t941 = sin(qJ(1));
t942 = cos(qJ(1));
t993 = t941 * t890 + t942 * t898;
t877 = t942 * t890 - t941 * t898;
t964 = t939 * t971;
t968 = qJDD(1) * t937;
t903 = (t964 + t968) * t938;
t910 = -t957 + t966;
t926 = t937 ^ 2 * t972;
t913 = -t926 - t935;
t949 = t910 * t937 - t913 * t939;
t892 = t938 * t903 + t949 * t940;
t899 = t939 * t910 + t937 * t913;
t992 = t941 * t892 - t942 * t899;
t878 = t942 * t892 + t941 * t899;
t888 = -t940 * t903 + t949 * t938;
t886 = -t940 * t906 + t950 * t938;
t985 = 2 * qJD(2);
t984 = 2 * qJD(3);
t983 = -2 * qJD(4);
t982 = qJ(2) * t938;
t981 = qJ(3) * t937;
t980 = t937 * t938;
t928 = g(3) * t980;
t975 = -pkin(2) * t935 - t928;
t974 = qJD(1) * t938;
t973 = qJD(1) * t940;
t955 = -pkin(2) * t939 - t981;
t970 = t955 * t974 + t985;
t969 = qJ(3) * qJDD(1);
t965 = qJ(3) * t939 * t940;
t961 = -t941 * g(1) + t942 * g(2);
t960 = -t940 * g(3) + qJDD(2);
t959 = -0.2e1 * t937 * t974;
t958 = t970 * qJD(1);
t947 = qJDD(1) * pkin(1) - t961;
t944 = t986 * t982 + t947;
t943 = t940 * t944;
t956 = (g(3) * t938 - t943) * t939;
t925 = -t942 * g(1) - t941 * g(2);
t908 = -t986 * pkin(1) + qJ(2) * t967 + t925;
t945 = -pkin(2) * t966 - qJ(3) * t935 + qJDD(3) + t956;
t871 = t973 * t983 - t939 * pkin(3) * t963 + (t908 + (pkin(3) * qJDD(1) + t958) * t938) * t937 + t945 - t909 * qJ(4);
t912 = (pkin(3) * t980 - qJ(4) * t940) * qJD(1);
t872 = qJDD(4) + (t937 * t947 + t969) * t940 + (pkin(3) * t967 + t908) * t939 + ((t984 + t912) * t940 + (qJ(2) * t937 * t973 + (-qJ(4) * t939 * t974 + t970) * t939) * t938) * qJD(1) + t975;
t954 = -t871 * t939 + t872 * t937;
t946 = t938 * t958 + t908;
t881 = (qJD(1) * t984 + t937 * t944 + t969) * t940 + t946 * t939 + t975;
t882 = t946 * t937 + t945;
t953 = t881 * t937 - t882 * t939;
t894 = qJD(2) * t959 - t937 * t908 - t956;
t895 = t937 * t943 - t928 + (t974 * t985 + t908) * t939;
t952 = t894 * t939 + t895 * t937;
t904 = (-t964 + t968) * t938;
t905 = t918 + t962;
t951 = -t904 * t939 + t905 * t937;
t948 = pkin(2) * t918 + qJD(3) * t959 + t960;
t920 = -t941 * qJDD(1) - t942 * t986;
t919 = t942 * qJDD(1) - t941 * t986;
t911 = -t926 - t927;
t901 = -t938 * t944 + t960;
t896 = t937 * t904 + t939 * t905;
t885 = -t938 * t911 + t951 * t940;
t884 = t940 * t911 + t951 * t938;
t883 = ((-pkin(1) + t955) * qJDD(1) + (-t965 - t982) * t986 + t961) * t938 + t948;
t876 = -t941 * t885 + t942 * t896;
t875 = t942 * t885 + t941 * t896;
t874 = ((-t981 - pkin(1) + (-pkin(2) - qJ(4)) * t939) * qJDD(1) + (t939 * t983 - t937 * t912 + (-t965 + (-pkin(3) * t987 - qJ(2)) * t938) * qJD(1)) * qJD(1) + t961) * t938 + t948;
t873 = -t937 * t894 + t939 * t895;
t870 = -t938 * t901 + t952 * t940;
t869 = t940 * t901 + t952 * t938;
t868 = t939 * t881 + t937 * t882;
t867 = t937 * t871 + t939 * t872;
t866 = -t938 * t883 + t953 * t940;
t865 = t940 * t883 + t953 * t938;
t864 = -t938 * t874 + t954 * t940;
t863 = t940 * t874 + t954 * t938;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t920, -t919, 0, t942 * t925 + t941 * t961, 0, 0, 0, 0, 0, 0, -t993, t992, t876, -t941 * t870 + t942 * t873, 0, 0, 0, 0, 0, 0, t876, t993, -t992, -t941 * t866 + t942 * t868, 0, 0, 0, 0, 0, 0, t876, -t992, -t993, -t941 * t864 + t942 * t867; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t919, t920, 0, t941 * t925 - t942 * t961, 0, 0, 0, 0, 0, 0, t877, -t878, t875, t942 * t870 + t941 * t873, 0, 0, 0, 0, 0, 0, t875, -t877, t878, t942 * t866 + t941 * t868, 0, 0, 0, 0, 0, 0, t875, t878, t877, t942 * t864 + t941 * t867; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t886, -t888, t884, t869, 0, 0, 0, 0, 0, 0, t884, -t886, t888, t865, 0, 0, 0, 0, 0, 0, t884, t888, t886, t863; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t986, -qJDD(1), 0, t925, 0, 0, 0, 0, 0, 0, -t898, -t899, t896, t873, 0, 0, 0, 0, 0, 0, t896, t898, t899, t868, 0, 0, 0, 0, 0, 0, t896, t899, -t898, t867; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t986, 0, -t961, 0, 0, 0, 0, 0, 0, t890, -t892, t885, t870, 0, 0, 0, 0, 0, 0, t885, -t890, t892, t866, 0, 0, 0, 0, 0, 0, t885, t892, t890, t864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t886, -t888, t884, t869, 0, 0, 0, 0, 0, 0, t884, -t886, t888, t865, 0, 0, 0, 0, 0, 0, t884, t888, t886, t863; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t914, -t910, t905, t895, 0, 0, 0, 0, 0, 0, t905, -t914, t910, t881, 0, 0, 0, 0, 0, 0, t905, t910, t914, t872; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t909, t913, -t904, t894, 0, 0, 0, 0, 0, 0, -t904, -t909, -t913, -t882, 0, 0, 0, 0, 0, 0, -t904, -t913, t909, -t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t906, t903, t911, t901, 0, 0, 0, 0, 0, 0, t911, t906, -t903, t883, 0, 0, 0, 0, 0, 0, t911, -t903, -t906, t874; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t911, t906, -t903, t883, 0, 0, 0, 0, 0, 0, t911, -t903, -t906, t874; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t905, t914, -t910, -t881, 0, 0, 0, 0, 0, 0, -t905, -t910, -t914, -t872; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t904, t909, t913, t882, 0, 0, 0, 0, 0, 0, t904, t913, -t909, t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t911, -t903, -t906, t874; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t904, t913, -t909, t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t905, t910, t914, t872;];
f_new_reg  = t1;
