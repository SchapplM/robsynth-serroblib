% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-05-06 12:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:42:11
% EndTime: 2019-05-06 12:42:17
% DurationCPUTime: 4.17s
% Computational Cost: add. (30425->328), mult. (61261->365), div. (0->0), fcn. (33015->6), ass. (0->131)
t1020 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t1064 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t1063 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t991 = cos(qJ(2));
t1024 = qJD(1) * t991;
t989 = sin(qJ(1));
t992 = cos(qJ(1));
t970 = -t992 * g(1) - t989 * g(2);
t994 = qJD(1) ^ 2;
t939 = -t994 * pkin(1) + qJDD(1) * pkin(7) + t970;
t988 = sin(qJ(2));
t926 = -t988 * g(3) + t991 * t939;
t955 = (-pkin(2) * t991 - qJ(3) * t988) * qJD(1);
t993 = qJD(2) ^ 2;
t1008 = t993 * pkin(2) - qJDD(2) * qJ(3) - t955 * t1024 - t926;
t1022 = qJD(1) * qJD(2);
t1016 = t988 * t1022;
t959 = t991 * qJDD(1) - t1016;
t1023 = t988 * qJD(1);
t968 = pkin(3) * t1023 - qJD(2) * pkin(8);
t985 = t991 ^ 2;
t1001 = t985 * t994 * pkin(8) - t959 * pkin(3) - qJD(2) * t968 + t1008;
t1042 = 2 * qJD(5);
t987 = sin(qJ(4));
t990 = cos(qJ(4));
t953 = t987 * qJD(2) + t990 * t1024;
t974 = qJD(4) + t1023;
t1033 = t953 * t974;
t911 = -t953 * qJD(4) + t990 * qJDD(2) - t987 * t959;
t1048 = (-t911 + t1033) * qJ(5);
t954 = t990 * qJD(2) - t987 * t1024;
t910 = t954 * qJD(4) + t987 * qJDD(2) + t990 * t959;
t921 = -t974 * pkin(5) - t954 * qJ(6);
t951 = t953 ^ 2;
t1021 = qJD(3) * qJD(2);
t978 = -0.2e1 * t1021;
t862 = -t951 * qJ(6) + qJDD(6) + t978 + (-pkin(4) - pkin(5)) * t910 - t1048 + (-pkin(4) * t974 + t1042 + t921) * t954 + t1001;
t918 = t974 * mrSges(7,2) + t953 * mrSges(7,3);
t1009 = m(7) * t862 - t910 * mrSges(7,1) - t953 * t918;
t1043 = -0.2e1 * t954;
t877 = -t1001 + 0.2e1 * t1021;
t868 = qJD(5) * t1043 + t1048 + (t954 * t974 + t910) * pkin(4) + t877;
t920 = -t953 * mrSges(6,2) + t974 * mrSges(6,3);
t922 = -t974 * mrSges(7,1) - t954 * mrSges(7,3);
t924 = -t974 * mrSges(6,1) + t954 * mrSges(6,2);
t855 = m(6) * t868 + t910 * mrSges(6,1) + t953 * t920 - t1009 - (t922 + t924) * t954 - (mrSges(7,2) + mrSges(6,3)) * t911;
t1062 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t1053 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t919 = -t974 * mrSges(5,2) - t953 * mrSges(5,3);
t923 = t974 * mrSges(5,1) - t954 * mrSges(5,3);
t1061 = -m(5) * t877 - t910 * mrSges(5,1) - t911 * mrSges(5,2) - t953 * t919 - t954 * t923 - t855;
t1060 = Ifges(3,1) + Ifges(4,2);
t1037 = Ifges(3,4) + Ifges(4,6);
t1036 = Ifges(3,5) - Ifges(4,4);
t1058 = Ifges(3,2) + Ifges(4,3);
t1035 = Ifges(3,6) - Ifges(4,5);
t1055 = Ifges(3,3) + Ifges(4,1);
t1054 = t1020 * t954 + t1063 * t974 + t1064 * t953;
t1052 = -Ifges(6,2) - Ifges(5,3) - Ifges(7,3);
t1050 = -t1020 * t953 + t1053 * t974 + t1062 * t954;
t1025 = t1036 * qJD(2) + (t1037 * t991 + t1060 * t988) * qJD(1);
t1026 = t1035 * qJD(2) + (t1037 * t988 + t1058 * t991) * qJD(1);
t1018 = t1052 * t974 - t1053 * t954 + t1063 * t953;
t969 = t989 * g(1) - t992 * g(2);
t1007 = -qJDD(1) * pkin(1) - t969;
t1015 = t991 * t1022;
t958 = t988 * qJDD(1) + t1015;
t1000 = pkin(2) * t1016 - 0.2e1 * qJD(3) * t1023 + (-t958 - t1015) * qJ(3) + t1007;
t874 = -t968 * t1023 + (-pkin(3) * t985 - pkin(7)) * t994 + (-pkin(2) - pkin(8)) * t959 + t1000;
t925 = -t991 * g(3) - t988 * t939;
t884 = -qJDD(2) * pkin(2) - t993 * qJ(3) + t955 * t1023 + qJDD(3) - t925;
t878 = (-t988 * t991 * t994 - qJDD(2)) * pkin(8) + (t958 - t1015) * pkin(3) + t884;
t871 = t990 * t874 + t987 * t878;
t914 = t953 * pkin(4) - t954 * qJ(5);
t952 = qJDD(4) + t958;
t971 = t974 ^ 2;
t865 = -t971 * pkin(4) + t952 * qJ(5) + t974 * t1042 - t953 * t914 + t871;
t861 = -t951 * pkin(5) + t910 * qJ(6) + 0.2e1 * qJD(6) * t953 + t974 * t921 + t865;
t916 = -t953 * mrSges(7,1) + t954 * mrSges(7,2);
t1019 = m(7) * t861 + t910 * mrSges(7,3) + t953 * t916;
t857 = t911 * mrSges(7,2) + t954 * t922 + t1009;
t834 = -mrSges(5,1) * t877 + mrSges(5,3) * t871 - mrSges(6,1) * t868 + mrSges(6,2) * t865 + mrSges(7,1) * t862 - mrSges(7,3) * t861 + pkin(5) * t857 - qJ(6) * t1019 - pkin(4) * t855 + (-qJ(6) * t922 + t1050) * t974 + t1018 * t954 + (-qJ(6) * mrSges(7,2) + t1063) * t952 + t1020 * t911 + t1064 * t910;
t870 = -t987 * t874 + t990 * t878;
t866 = -t952 * pkin(4) - t971 * qJ(5) + t954 * t914 + qJDD(5) - t870;
t858 = qJD(6) * t1043 + (-t911 - t1033) * qJ(6) + (t953 * t954 - t952) * pkin(5) + t866;
t1010 = -m(7) * t858 + t911 * mrSges(7,3) + t954 * t916;
t856 = -t952 * mrSges(7,1) - t974 * t918 - t1010;
t835 = mrSges(5,2) * t877 + mrSges(6,2) * t866 + mrSges(7,2) * t862 - mrSges(5,3) * t870 - mrSges(6,3) * t868 - mrSges(7,3) * t858 - qJ(5) * t855 - qJ(6) * t856 + t1018 * t953 - t1020 * t910 + t1053 * t952 - t1054 * t974 + t1062 * t911;
t915 = t953 * mrSges(6,1) - t954 * mrSges(6,3);
t1029 = -t953 * mrSges(5,1) - t954 * mrSges(5,2) - t915;
t1038 = -mrSges(5,3) - mrSges(6,2);
t1046 = -m(6) * t866 + t952 * mrSges(6,1) + t974 * t920;
t849 = m(5) * t870 + (t918 + t919) * t974 + t1029 * t954 + (mrSges(5,1) + mrSges(7,1)) * t952 + t1038 * t911 + t1010 + t1046;
t1006 = m(6) * t865 + t952 * mrSges(6,3) + t974 * t924 + t1019;
t850 = m(5) * t871 + (t922 - t923) * t974 + t1029 * t953 + (-mrSges(5,2) + mrSges(7,2)) * t952 + t1038 * t910 + t1006;
t843 = t990 * t849 + t987 * t850;
t1003 = -m(4) * t884 - t958 * mrSges(4,1) - t843;
t956 = (mrSges(4,2) * t991 - mrSges(4,3) * t988) * qJD(1);
t965 = -mrSges(4,1) * t1024 - qJD(2) * mrSges(4,3);
t842 = qJDD(2) * mrSges(4,2) + qJD(2) * t965 + t956 * t1023 - t1003;
t883 = t978 + t1008;
t966 = mrSges(4,1) * t1023 + qJD(2) * mrSges(4,2);
t996 = -m(4) * t883 + qJDD(2) * mrSges(4,3) + qJD(2) * t966 + t956 * t1024 - t1061;
t1049 = -(t1025 * t991 - t1026 * t988) * qJD(1) + t1055 * qJDD(2) + t1035 * t959 + t1036 * t958 + mrSges(3,1) * t925 - mrSges(3,2) * t926 + mrSges(4,2) * t884 - mrSges(4,3) * t883 - pkin(2) * t842 - pkin(8) * t843 + qJ(3) * (t959 * mrSges(4,1) + t996) - t987 * t834 + t990 * t835;
t1040 = t994 * pkin(7);
t957 = (-mrSges(3,1) * t991 + mrSges(3,2) * t988) * qJD(1);
t964 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1024;
t840 = m(3) * t925 - t958 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t964 - t965) * qJD(2) + (-t956 - t957) * t1023 + t1003;
t963 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1023;
t846 = t996 + t957 * t1024 + (mrSges(3,3) + mrSges(4,1)) * t959 - qJDD(2) * mrSges(3,2) - qJD(2) * t963 + m(3) * t926;
t1012 = -t988 * t840 + t991 * t846;
t831 = m(2) * t970 - t994 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1012;
t1031 = -t987 * t849 + t990 * t850;
t879 = -t959 * pkin(2) + t1000 - t1040;
t1005 = -m(4) * t879 - t959 * mrSges(4,2) + t966 * t1023 - t1031;
t938 = t1007 - t1040;
t998 = -m(3) * t938 + t964 * t1024 + t959 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t958 + (-t963 * t988 - t965 * t991) * qJD(1) + t1005;
t837 = m(2) * t969 + qJDD(1) * mrSges(2,1) - t994 * mrSges(2,2) + t998;
t1032 = t989 * t831 + t992 * t837;
t833 = t991 * t840 + t988 * t846;
t1027 = t1055 * qJD(2) + (t1035 * t991 + t1036 * t988) * qJD(1);
t1013 = t992 * t831 - t989 * t837;
t841 = -t958 * mrSges(4,3) + t965 * t1024 - t1005;
t826 = -mrSges(3,1) * t938 - mrSges(4,1) * t883 + mrSges(4,2) * t879 + mrSges(3,3) * t926 - pkin(2) * t841 - pkin(3) * t1061 - pkin(8) * t1031 + t1025 * qJD(2) + t1035 * qJDD(2) - t1027 * t1023 + t1037 * t958 + t1058 * t959 - t990 * t834 - t987 * t835;
t853 = t911 * mrSges(6,2) + t954 * t915 - t1046 + t856;
t995 = -mrSges(6,1) * t866 - mrSges(7,1) * t858 - mrSges(5,2) * t871 - pkin(5) * t856 - pkin(4) * t853 + qJ(5) * (t974 * t922 + t1006) + mrSges(7,2) * t861 + mrSges(6,3) * t865 + mrSges(5,1) * t870 + t1054 * t954 + (-qJ(5) * t915 + t1050) * t953 + (mrSges(7,2) * qJ(5) - t1052) * t952 + t1053 * t911 + (-mrSges(6,2) * qJ(5) - t1063) * t910;
t828 = mrSges(4,1) * t884 + mrSges(3,2) * t938 - mrSges(3,3) * t925 - mrSges(4,3) * t879 + pkin(3) * t843 - qJ(3) * t841 - t1026 * qJD(2) + t1036 * qJDD(2) + t1027 * t1024 + t1037 * t959 + t1060 * t958 + t995;
t1002 = mrSges(2,1) * t969 - mrSges(2,2) * t970 + Ifges(2,3) * qJDD(1) + pkin(1) * t998 + pkin(7) * t1012 + t991 * t826 + t988 * t828;
t824 = mrSges(2,1) * g(3) + mrSges(2,3) * t970 + t994 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t833 - t1049;
t823 = -mrSges(2,2) * g(3) - mrSges(2,3) * t969 + Ifges(2,5) * qJDD(1) - t994 * Ifges(2,6) - pkin(7) * t833 - t988 * t826 + t991 * t828;
t1 = [-m(1) * g(1) + t1013; -m(1) * g(2) + t1032; (-m(1) - m(2)) * g(3) + t833; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1032 + t992 * t823 - t989 * t824; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1013 + t989 * t823 + t992 * t824; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1002; t1002; t1049; t842; t995; t853; t857;];
tauJB  = t1;
