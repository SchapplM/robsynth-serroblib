% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 14:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:43:56
% EndTime: 2019-05-06 14:44:08
% DurationCPUTime: 10.66s
% Computational Cost: add. (158170->365), mult. (350961->447), div. (0->0), fcn. (232191->10), ass. (0->147)
t1035 = Ifges(3,1) + Ifges(4,1);
t1026 = Ifges(3,4) - Ifges(4,5);
t1025 = (Ifges(3,5) + Ifges(4,4));
t1034 = Ifges(3,2) + Ifges(4,3);
t1024 = (Ifges(3,6) - Ifges(4,6));
t1033 = Ifges(3,3) + Ifges(4,2);
t1028 = 2 * qJD(5);
t993 = cos(qJ(2));
t996 = qJD(1) ^ 2;
t1022 = t993 ^ 2 * t996;
t1016 = qJD(1) * t993;
t1029 = 2 * qJD(3);
t990 = sin(qJ(1));
t994 = cos(qJ(1));
t965 = -g(1) * t994 - g(2) * t990;
t942 = -pkin(1) * t996 + qJDD(1) * pkin(7) + t965;
t989 = sin(qJ(2));
t923 = -g(3) * t989 + t942 * t993;
t951 = (-pkin(2) * t993 - qJ(3) * t989) * qJD(1);
t995 = qJD(2) ^ 2;
t899 = -pkin(2) * t995 + qJDD(2) * qJ(3) + (qJD(2) * t1029) + t1016 * t951 + t923;
t1015 = qJD(1) * qJD(2);
t1013 = t989 * t1015;
t955 = qJDD(1) * t993 - t1013;
t1017 = qJD(1) * t989;
t963 = -(qJD(2) * pkin(3)) - pkin(8) * t1017;
t894 = -pkin(3) * t1022 - pkin(8) * t955 + qJD(2) * t963 + t899;
t1014 = t993 * t1015;
t922 = -g(3) * t993 - t942 * t989;
t906 = -qJDD(2) * pkin(2) - t995 * qJ(3) + t951 * t1017 + qJDD(3) - t922;
t954 = qJDD(1) * t989 + t1014;
t895 = (-t954 + t1014) * pkin(8) + (-t989 * t993 * t996 - qJDD(2)) * pkin(3) + t906;
t988 = sin(qJ(4));
t992 = cos(qJ(4));
t863 = -t988 * t894 + t895 * t992;
t939 = (-t988 * t989 - t992 * t993) * qJD(1);
t908 = qJD(4) * t939 + t954 * t992 - t955 * t988;
t940 = (-t988 * t993 + t989 * t992) * qJD(1);
t976 = -qJDD(2) + qJDD(4);
t977 = -qJD(2) + qJD(4);
t856 = (t939 * t977 - t908) * qJ(5) + (t939 * t940 + t976) * pkin(4) + t863;
t864 = t894 * t992 + t895 * t988;
t907 = -qJD(4) * t940 - t954 * t988 - t955 * t992;
t925 = pkin(4) * t977 - qJ(5) * t940;
t932 = t939 ^ 2;
t858 = -pkin(4) * t932 + qJ(5) * t907 - t925 * t977 + t864;
t984 = sin(pkin(10));
t985 = cos(pkin(10));
t919 = t939 * t985 - t940 * t984;
t853 = t1028 * t919 + t856 * t984 + t858 * t985;
t920 = t939 * t984 + t940 * t985;
t893 = -pkin(5) * t919 - pkin(9) * t920;
t975 = t977 ^ 2;
t850 = -pkin(5) * t975 + pkin(9) * t976 + t893 * t919 + t853;
t964 = g(1) * t990 - g(2) * t994;
t941 = -qJDD(1) * pkin(1) - pkin(7) * t996 - t964;
t1006 = -t955 * pkin(2) + t941 + (-t1014 - t954) * qJ(3);
t881 = -pkin(2) * t1013 + t955 * pkin(3) - pkin(8) * t1022 - t1006 + (t1029 + t963) * t1017;
t860 = -t907 * pkin(4) - t932 * qJ(5) + t925 * t940 + qJDD(5) + t881;
t877 = t907 * t985 - t908 * t984;
t878 = t907 * t984 + t908 * t985;
t854 = t860 + (-t919 * t977 - t878) * pkin(9) + (t920 * t977 - t877) * pkin(5);
t987 = sin(qJ(6));
t991 = cos(qJ(6));
t847 = -t850 * t987 + t854 * t991;
t904 = -t920 * t987 + t977 * t991;
t866 = qJD(6) * t904 + t878 * t991 + t976 * t987;
t876 = qJDD(6) - t877;
t905 = t920 * t991 + t977 * t987;
t879 = -mrSges(7,1) * t904 + mrSges(7,2) * t905;
t912 = qJD(6) - t919;
t882 = -mrSges(7,2) * t912 + mrSges(7,3) * t904;
t843 = m(7) * t847 + mrSges(7,1) * t876 - mrSges(7,3) * t866 - t879 * t905 + t882 * t912;
t848 = t850 * t991 + t854 * t987;
t865 = -qJD(6) * t905 - t878 * t987 + t976 * t991;
t883 = mrSges(7,1) * t912 - mrSges(7,3) * t905;
t844 = m(7) * t848 - mrSges(7,2) * t876 + mrSges(7,3) * t865 + t879 * t904 - t883 * t912;
t1008 = -t843 * t987 + t844 * t991;
t892 = -mrSges(6,1) * t919 + mrSges(6,2) * t920;
t910 = mrSges(6,1) * t977 - mrSges(6,3) * t920;
t829 = m(6) * t853 - mrSges(6,2) * t976 + mrSges(6,3) * t877 + t892 * t919 - t910 * t977 + t1008;
t1007 = -t985 * t856 + t984 * t858;
t849 = -t976 * pkin(5) - t975 * pkin(9) + (t1028 + t893) * t920 + t1007;
t1003 = -m(7) * t849 + mrSges(7,1) * t865 - mrSges(7,2) * t866 + t882 * t904 - t883 * t905;
t852 = -0.2e1 * qJD(5) * t920 - t1007;
t909 = -mrSges(6,2) * t977 + mrSges(6,3) * t919;
t839 = m(6) * t852 + mrSges(6,1) * t976 - mrSges(6,3) * t878 - t892 * t920 + t909 * t977 + t1003;
t822 = t829 * t984 + t839 * t985;
t921 = -mrSges(5,1) * t939 + mrSges(5,2) * t940;
t924 = -mrSges(5,2) * t977 + mrSges(5,3) * t939;
t819 = m(5) * t863 + mrSges(5,1) * t976 - mrSges(5,3) * t908 - t921 * t940 + t924 * t977 + t822;
t1009 = t829 * t985 - t839 * t984;
t926 = mrSges(5,1) * t977 - mrSges(5,3) * t940;
t820 = m(5) * t864 - mrSges(5,2) * t976 + mrSges(5,3) * t907 + t921 * t939 - t926 * t977 + t1009;
t1010 = -t988 * t819 + t820 * t992;
t952 = (-mrSges(4,1) * t993 - mrSges(4,3) * t989) * qJD(1);
t960 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t1017;
t1005 = m(4) * t899 + qJDD(2) * mrSges(4,3) + qJD(2) * t960 + t1016 * t952 + t1010;
t1018 = (t1025 * qJD(2)) + (t1026 * t993 + t1035 * t989) * qJD(1);
t1019 = -(t1024 * qJD(2)) + (-t1026 * t989 - t1034 * t993) * qJD(1);
t815 = t992 * t819 + t988 * t820;
t962 = mrSges(4,2) * t1016 + qJD(2) * mrSges(4,3);
t1002 = -m(4) * t906 + qJDD(2) * mrSges(4,1) + qJD(2) * t962 - t815;
t814 = t954 * mrSges(4,2) + t1017 * t952 - t1002;
t867 = Ifges(7,5) * t905 + Ifges(7,6) * t904 + Ifges(7,3) * t912;
t869 = Ifges(7,1) * t905 + Ifges(7,4) * t904 + Ifges(7,5) * t912;
t836 = -mrSges(7,1) * t849 + mrSges(7,3) * t848 + Ifges(7,4) * t866 + Ifges(7,2) * t865 + Ifges(7,6) * t876 - t867 * t905 + t869 * t912;
t868 = Ifges(7,4) * t905 + Ifges(7,2) * t904 + Ifges(7,6) * t912;
t837 = mrSges(7,2) * t849 - mrSges(7,3) * t847 + Ifges(7,1) * t866 + Ifges(7,4) * t865 + Ifges(7,5) * t876 + t867 * t904 - t868 * t912;
t885 = Ifges(6,4) * t920 + Ifges(6,2) * t919 + Ifges(6,6) * t977;
t886 = Ifges(6,1) * t920 + Ifges(6,4) * t919 + Ifges(6,5) * t977;
t914 = Ifges(5,4) * t940 + Ifges(5,2) * t939 + Ifges(5,6) * t977;
t915 = Ifges(5,1) * t940 + Ifges(5,4) * t939 + Ifges(5,5) * t977;
t999 = -mrSges(5,1) * t863 - mrSges(6,1) * t852 + mrSges(5,2) * t864 + mrSges(6,2) * t853 - pkin(4) * t822 - pkin(5) * t1003 - pkin(9) * t1008 - t991 * t836 - t987 * t837 - t920 * t885 + t919 * t886 + t939 * t915 - Ifges(6,6) * t877 - Ifges(6,5) * t878 - t940 * t914 - Ifges(5,6) * t907 - Ifges(5,5) * t908 + (-Ifges(5,3) - Ifges(6,3)) * t976;
t1032 = -(t1018 * t993 + t1019 * t989) * qJD(1) + t1033 * qJDD(2) + t1024 * t955 + t1025 * t954 + mrSges(3,1) * t922 - mrSges(4,1) * t906 - mrSges(3,2) * t923 + mrSges(4,3) * t899 - pkin(2) * t814 - pkin(3) * t815 + qJ(3) * (mrSges(4,2) * t955 + t1005) + t999;
t1027 = mrSges(3,3) + mrSges(4,2);
t953 = (-mrSges(3,1) * t993 + mrSges(3,2) * t989) * qJD(1);
t959 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1017;
t811 = m(3) * t923 - qJDD(2) * mrSges(3,2) - qJD(2) * t959 + t1016 * t953 + t1027 * t955 + t1005;
t961 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1016;
t812 = m(3) * t922 + qJDD(2) * mrSges(3,1) + qJD(2) * t961 - t1027 * t954 + (-t952 - t953) * t1017 + t1002;
t1011 = t811 * t993 - t812 * t989;
t804 = m(2) * t965 - mrSges(2,1) * t996 - qJDD(1) * mrSges(2,2) + t1011;
t832 = t843 * t991 + t844 * t987;
t830 = m(6) * t860 - mrSges(6,1) * t877 + mrSges(6,2) * t878 - t909 * t919 + t910 * t920 + t832;
t1001 = -m(5) * t881 + mrSges(5,1) * t907 - mrSges(5,2) * t908 + t924 * t939 - t926 * t940 - t830;
t896 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t1017 + t1006;
t826 = m(4) * t896 - mrSges(4,1) * t955 - mrSges(4,3) * t954 - t1016 * t962 - t1017 * t960 + t1001;
t998 = -m(3) * t941 + mrSges(3,1) * t955 - mrSges(3,2) * t954 + t1016 * t961 - t1017 * t959 - t826;
t824 = m(2) * t964 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t996 + t998;
t1021 = t804 * t990 + t824 * t994;
t806 = t811 * t989 + t812 * t993;
t1020 = t1033 * qJD(2) + (t1024 * t993 + t1025 * t989) * qJD(1);
t1012 = t804 * t994 - t824 * t990;
t884 = Ifges(6,5) * t920 + Ifges(6,6) * t919 + Ifges(6,3) * t977;
t816 = mrSges(6,2) * t860 - mrSges(6,3) * t852 + Ifges(6,1) * t878 + Ifges(6,4) * t877 + Ifges(6,5) * t976 - pkin(9) * t832 - t836 * t987 + t837 * t991 + t884 * t919 - t885 * t977;
t1000 = mrSges(7,1) * t847 - mrSges(7,2) * t848 + Ifges(7,5) * t866 + Ifges(7,6) * t865 + Ifges(7,3) * t876 + t868 * t905 - t869 * t904;
t817 = -mrSges(6,1) * t860 + mrSges(6,3) * t853 + Ifges(6,4) * t878 + Ifges(6,2) * t877 + Ifges(6,6) * t976 - pkin(5) * t832 - t884 * t920 + t886 * t977 - t1000;
t913 = Ifges(5,5) * t940 + Ifges(5,6) * t939 + Ifges(5,3) * t977;
t801 = -mrSges(5,1) * t881 + mrSges(5,3) * t864 + Ifges(5,4) * t908 + Ifges(5,2) * t907 + Ifges(5,6) * t976 - pkin(4) * t830 + qJ(5) * t1009 + t984 * t816 + t985 * t817 - t940 * t913 + t977 * t915;
t807 = mrSges(5,2) * t881 - mrSges(5,3) * t863 + Ifges(5,1) * t908 + Ifges(5,4) * t907 + Ifges(5,5) * t976 - qJ(5) * t822 + t816 * t985 - t817 * t984 + t913 * t939 - t914 * t977;
t798 = -mrSges(3,1) * t941 - mrSges(4,1) * t896 + mrSges(4,2) * t899 + mrSges(3,3) * t923 - pkin(2) * t826 - pkin(3) * t1001 - pkin(8) * t1010 + qJD(2) * t1018 + qJDD(2) * t1024 - t1017 * t1020 + t1026 * t954 + t1034 * t955 - t992 * t801 - t988 * t807;
t800 = mrSges(3,2) * t941 + mrSges(4,2) * t906 - mrSges(3,3) * t922 - mrSges(4,3) * t896 - pkin(8) * t815 - qJ(3) * t826 + qJD(2) * t1019 + qJDD(2) * t1025 + t1016 * t1020 + t1026 * t955 + t1035 * t954 - t988 * t801 + t992 * t807;
t1004 = mrSges(2,1) * t964 - mrSges(2,2) * t965 + Ifges(2,3) * qJDD(1) + pkin(1) * t998 + pkin(7) * t1011 + t798 * t993 + t800 * t989;
t796 = mrSges(2,1) * g(3) + mrSges(2,3) * t965 + t996 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t806 - t1032;
t795 = -mrSges(2,2) * g(3) - mrSges(2,3) * t964 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t996 - pkin(7) * t806 - t798 * t989 + t800 * t993;
t1 = [-m(1) * g(1) + t1012; -m(1) * g(2) + t1021; (-m(1) - m(2)) * g(3) + t806; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1021 + t795 * t994 - t796 * t990; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1012 + t990 * t795 + t994 * t796; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1004; t1004; t1032; t814; -t999; t830; t1000;];
tauJB  = t1;
