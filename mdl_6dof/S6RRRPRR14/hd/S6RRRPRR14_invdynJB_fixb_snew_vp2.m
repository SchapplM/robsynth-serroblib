% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 16:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:35:24
% EndTime: 2019-05-07 16:35:50
% DurationCPUTime: 24.76s
% Computational Cost: add. (419089->382), mult. (897237->475), div. (0->0), fcn. (699843->12), ass. (0->160)
t1065 = Ifges(4,1) + Ifges(5,2);
t1060 = Ifges(4,5) - Ifges(5,4);
t1064 = -Ifges(4,2) - Ifges(5,3);
t1059 = Ifges(4,6) - Ifges(5,5);
t1058 = -Ifges(5,6) - Ifges(4,4);
t1063 = Ifges(4,3) + Ifges(5,1);
t1011 = sin(pkin(6));
t1016 = sin(qJ(2));
t1020 = cos(qJ(2));
t995 = (qJD(1) * qJD(2) * t1016 - qJDD(1) * t1020) * t1011;
t1043 = t1011 * t1020;
t1039 = qJD(1) * t1043;
t1001 = -qJD(3) + t1039;
t1000 = t1001 ^ 2;
t1015 = sin(qJ(3));
t1061 = cos(qJ(3));
t1012 = cos(pkin(6));
t1007 = t1012 * qJD(1) + qJD(2);
t1005 = t1007 ^ 2;
t1006 = t1012 * qJDD(1) + qJDD(2);
t1046 = qJD(1) * t1020;
t1042 = t1012 * t1016;
t1017 = sin(qJ(1));
t1021 = cos(qJ(1));
t1003 = t1017 * g(1) - t1021 * g(2);
t1022 = qJD(1) ^ 2;
t1056 = pkin(8) * t1011;
t990 = qJDD(1) * pkin(1) + t1022 * t1056 + t1003;
t1004 = -t1021 * g(1) - t1017 * g(2);
t991 = -t1022 * pkin(1) + qJDD(1) * t1056 + t1004;
t1049 = t1020 * t991 + t990 * t1042;
t1047 = qJD(1) * t1011;
t993 = (-pkin(2) * t1020 - pkin(9) * t1016) * t1047;
t927 = -t1005 * pkin(2) + t1006 * pkin(9) + (-g(3) * t1016 + t993 * t1046) * t1011 + t1049;
t1055 = t1012 * g(3);
t994 = (qJD(2) * t1046 + qJDD(1) * t1016) * t1011;
t928 = t995 * pkin(2) - t994 * pkin(9) - t1055 + (-t990 + (pkin(2) * t1016 - pkin(9) * t1020) * t1007 * qJD(1)) * t1011;
t905 = t1015 * t928 + t1061 * t927;
t1044 = t1011 * t1016;
t1040 = qJD(1) * t1044;
t981 = -t1061 * t1007 + t1015 * t1040;
t982 = t1015 * t1007 + t1061 * t1040;
t957 = t981 * pkin(3) - t982 * qJ(4);
t987 = qJDD(3) + t995;
t900 = t1000 * pkin(3) - t987 * qJ(4) + 0.2e1 * qJD(4) * t1001 + t981 * t957 - t905;
t1014 = sin(qJ(5));
t1019 = cos(qJ(5));
t953 = t982 * qJD(3) - t1061 * t1006 + t1015 * t994;
t969 = t982 * pkin(4) + t1001 * pkin(10);
t980 = t981 ^ 2;
t893 = -t953 * pkin(4) - t980 * pkin(10) - t1001 * t969 - t900;
t964 = -t1019 * t1001 + t1014 * t981;
t915 = -t964 * qJD(5) - t1014 * t987 + t1019 * t953;
t979 = qJD(5) + t982;
t937 = t979 * pkin(5) - t964 * pkin(11);
t963 = t1014 * t1001 + t1019 * t981;
t962 = t963 ^ 2;
t887 = -t915 * pkin(5) - t962 * pkin(11) + t964 * t937 + t893;
t1013 = sin(qJ(6));
t1018 = cos(qJ(6));
t916 = t963 * qJD(5) + t1014 * t953 + t1019 * t987;
t930 = t1013 * t963 + t1018 * t964;
t898 = -t930 * qJD(6) - t1013 * t916 + t1018 * t915;
t929 = -t1013 * t964 + t1018 * t963;
t899 = t929 * qJD(6) + t1013 * t915 + t1018 * t916;
t977 = qJD(6) + t979;
t917 = -t977 * mrSges(7,2) + t929 * mrSges(7,3);
t918 = t977 * mrSges(7,1) - t930 * mrSges(7,3);
t1032 = m(7) * t887 - t898 * mrSges(7,1) + t899 * mrSges(7,2) - t929 * t917 + t930 * t918;
t935 = -t979 * mrSges(6,2) + t963 * mrSges(6,3);
t936 = t979 * mrSges(6,1) - t964 * mrSges(6,3);
t1028 = -m(6) * t893 + t915 * mrSges(6,1) - t916 * mrSges(6,2) + t963 * t935 - t964 * t936 - t1032;
t968 = t982 * mrSges(5,1) - t1001 * mrSges(5,2);
t1025 = -m(5) * t900 + t987 * mrSges(5,3) - t1001 * t968 - t1028;
t1050 = -t1060 * t1001 + t1058 * t981 + t1065 * t982;
t1051 = -t1059 * t1001 - t1058 * t982 + t1064 * t981;
t1048 = t1001 * t981;
t904 = -t1015 * t927 + t1061 * t928;
t901 = -t987 * pkin(3) - t1000 * qJ(4) + t982 * t957 + qJDD(4) - t904;
t954 = -t981 * qJD(3) + t1015 * t1006 + t1061 * t994;
t890 = (t981 * t982 - t987) * pkin(10) + (t954 - t1048) * pkin(4) + t901;
t1041 = t1012 * t1020;
t955 = -g(3) * t1043 - t1016 * t991 + t990 * t1041;
t926 = -t1006 * pkin(2) - t1005 * pkin(9) + t993 * t1040 - t955;
t1026 = (-t954 - t1048) * qJ(4) + t926 + (-t1001 * pkin(3) - 0.2e1 * qJD(4)) * t982;
t894 = -t980 * pkin(4) - t982 * t969 + (pkin(3) + pkin(10)) * t953 + t1026;
t884 = -t1014 * t894 + t1019 * t890;
t950 = qJDD(5) + t954;
t881 = (t963 * t979 - t916) * pkin(11) + (t963 * t964 + t950) * pkin(5) + t884;
t885 = t1014 * t890 + t1019 * t894;
t882 = -t962 * pkin(5) + t915 * pkin(11) - t979 * t937 + t885;
t879 = -t1013 * t882 + t1018 * t881;
t910 = -mrSges(7,1) * t929 + mrSges(7,2) * t930;
t940 = qJDD(6) + t950;
t876 = m(7) * t879 + t940 * mrSges(7,1) - t899 * mrSges(7,3) - t930 * t910 + t977 * t917;
t880 = t1013 * t881 + t1018 * t882;
t877 = m(7) * t880 - t940 * mrSges(7,2) + t898 * mrSges(7,3) + t929 * t910 - t977 * t918;
t1038 = -t1013 * t876 + t1018 * t877;
t906 = Ifges(7,5) * t930 + Ifges(7,6) * t929 + Ifges(7,3) * t977;
t908 = Ifges(7,1) * t930 + Ifges(7,4) * t929 + Ifges(7,5) * t977;
t867 = -mrSges(7,1) * t887 + mrSges(7,3) * t880 + Ifges(7,4) * t899 + Ifges(7,2) * t898 + Ifges(7,6) * t940 - t930 * t906 + t977 * t908;
t907 = Ifges(7,4) * t930 + Ifges(7,2) * t929 + Ifges(7,6) * t977;
t868 = mrSges(7,2) * t887 - mrSges(7,3) * t879 + Ifges(7,1) * t899 + Ifges(7,4) * t898 + Ifges(7,5) * t940 + t929 * t906 - t977 * t907;
t919 = Ifges(6,5) * t964 + Ifges(6,6) * t963 + Ifges(6,3) * t979;
t921 = Ifges(6,1) * t964 + Ifges(6,4) * t963 + Ifges(6,5) * t979;
t851 = -mrSges(6,1) * t893 + mrSges(6,3) * t885 + Ifges(6,4) * t916 + Ifges(6,2) * t915 + Ifges(6,6) * t950 - pkin(5) * t1032 + pkin(11) * t1038 + t1013 * t868 + t1018 * t867 - t964 * t919 + t979 * t921;
t866 = t1013 * t877 + t1018 * t876;
t920 = Ifges(6,4) * t964 + Ifges(6,2) * t963 + Ifges(6,6) * t979;
t852 = mrSges(6,2) * t893 - mrSges(6,3) * t884 + Ifges(6,1) * t916 + Ifges(6,4) * t915 + Ifges(6,5) * t950 - pkin(11) * t866 - t1013 * t867 + t1018 * t868 + t963 * t919 - t979 * t920;
t931 = -mrSges(6,1) * t963 + mrSges(6,2) * t964;
t863 = m(6) * t884 + t950 * mrSges(6,1) - t916 * mrSges(6,3) - t964 * t931 + t979 * t935 + t866;
t864 = m(6) * t885 - t950 * mrSges(6,2) + t915 * mrSges(6,3) + t963 * t931 - t979 * t936 + t1038;
t860 = t1014 * t864 + t1019 * t863;
t959 = -t981 * mrSges(5,2) - t982 * mrSges(5,3);
t1031 = -m(5) * t901 - t954 * mrSges(5,1) - t982 * t959 - t860;
t967 = t981 * mrSges(5,1) + t1001 * mrSges(5,3);
t859 = t987 * mrSges(5,2) - t1001 * t967 - t1031;
t1062 = t1050 * t981 + t1051 * t982 + t1063 * t987 - t1059 * t953 + t1060 * t954 + mrSges(4,1) * t904 - mrSges(4,2) * t905 + mrSges(5,2) * t901 - mrSges(5,3) * t900 - pkin(3) * t859 - pkin(10) * t860 + qJ(4) * (-t953 * mrSges(5,1) - t981 * t959 + t1025) - t1014 * t851 + t1019 * t852;
t958 = t981 * mrSges(4,1) + t982 * mrSges(4,2);
t965 = t1001 * mrSges(4,2) - t981 * mrSges(4,3);
t857 = m(4) * t904 - t954 * mrSges(4,3) - t982 * t958 + (mrSges(4,1) - mrSges(5,2)) * t987 + (-t965 + t967) * t1001 + t1031;
t966 = -t1001 * mrSges(4,1) - t982 * mrSges(4,3);
t871 = (-t958 - t959) * t981 + (-mrSges(4,3) - mrSges(5,1)) * t953 + t1001 * t966 - t987 * mrSges(4,2) + m(4) * t905 + t1025;
t1037 = -t1015 * t857 + t1061 * t871;
t956 = -g(3) * t1044 + t1049;
t988 = t1007 * mrSges(3,1) - mrSges(3,3) * t1040;
t992 = (-mrSges(3,1) * t1020 + mrSges(3,2) * t1016) * t1047;
t847 = m(3) * t956 - t1006 * mrSges(3,2) - t995 * mrSges(3,3) - t1007 * t988 + t992 * t1039 + t1037;
t850 = t1015 * t871 + t1061 * t857;
t974 = -t1011 * t990 - t1055;
t989 = -t1007 * mrSges(3,2) + mrSges(3,3) * t1039;
t849 = m(3) * t974 + t995 * mrSges(3,1) + t994 * mrSges(3,2) + (t1016 * t988 - t1020 * t989) * t1047 + t850;
t1053 = -t1014 * t863 + t1019 * t864;
t902 = t953 * pkin(3) + t1026;
t1035 = -m(5) * t902 + t953 * mrSges(5,2) + t981 * t967 - t1053;
t1027 = -m(4) * t926 - t953 * mrSges(4,1) - t981 * t965 + (-t966 + t968) * t982 + (-mrSges(4,2) + mrSges(5,3)) * t954 + t1035;
t855 = m(3) * t955 + t1006 * mrSges(3,1) - t994 * mrSges(3,3) + t1007 * t989 - t992 * t1040 + t1027;
t837 = -t1011 * t849 + t855 * t1041 + t847 * t1042;
t834 = m(2) * t1003 + qJDD(1) * mrSges(2,1) - t1022 * mrSges(2,2) + t837;
t843 = -t1016 * t855 + t1020 * t847;
t841 = m(2) * t1004 - t1022 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t843;
t1054 = t1017 * t841 + t1021 * t834;
t1052 = t1063 * t1001 + t1059 * t981 - t1060 * t982;
t836 = t1012 * t849 + t855 * t1043 + t847 * t1044;
t1036 = -t1017 * t834 + t1021 * t841;
t858 = -t954 * mrSges(5,3) - t982 * t968 - t1035;
t832 = -mrSges(4,1) * t926 - mrSges(5,1) * t900 + mrSges(5,2) * t902 + mrSges(4,3) * t905 - pkin(3) * t858 - pkin(4) * t1028 - pkin(10) * t1053 - t1050 * t1001 - t1014 * t852 - t1019 * t851 + t1052 * t982 - t1058 * t954 + t1059 * t987 + t1064 * t953;
t1029 = mrSges(7,1) * t879 - mrSges(7,2) * t880 + Ifges(7,5) * t899 + Ifges(7,6) * t898 + Ifges(7,3) * t940 + t930 * t907 - t929 * t908;
t1024 = mrSges(6,1) * t884 - mrSges(6,2) * t885 + Ifges(6,5) * t916 + Ifges(6,6) * t915 + Ifges(6,3) * t950 + pkin(5) * t866 + t964 * t920 - t963 * t921 + t1029;
t838 = mrSges(5,1) * t901 + mrSges(4,2) * t926 - mrSges(4,3) * t904 - mrSges(5,3) * t902 + pkin(4) * t860 - qJ(4) * t858 + t1051 * t1001 + t1052 * t981 + t1058 * t953 + t1060 * t987 + t1065 * t954 + t1024;
t972 = Ifges(3,6) * t1007 + (Ifges(3,4) * t1016 + Ifges(3,2) * t1020) * t1047;
t973 = Ifges(3,5) * t1007 + (Ifges(3,1) * t1016 + Ifges(3,4) * t1020) * t1047;
t827 = Ifges(3,5) * t994 - Ifges(3,6) * t995 + Ifges(3,3) * t1006 + mrSges(3,1) * t955 - mrSges(3,2) * t956 + t1015 * t838 + t1061 * t832 + pkin(2) * t1027 + pkin(9) * t1037 + (t1016 * t972 - t1020 * t973) * t1047;
t971 = Ifges(3,3) * t1007 + (Ifges(3,5) * t1016 + Ifges(3,6) * t1020) * t1047;
t829 = mrSges(3,2) * t974 - mrSges(3,3) * t955 + Ifges(3,1) * t994 - Ifges(3,4) * t995 + Ifges(3,5) * t1006 - pkin(9) * t850 - t1007 * t972 - t1015 * t832 + t971 * t1039 + t1061 * t838;
t831 = -mrSges(3,1) * t974 + mrSges(3,3) * t956 + Ifges(3,4) * t994 - Ifges(3,2) * t995 + Ifges(3,6) * t1006 - pkin(2) * t850 + t1007 * t973 - t971 * t1040 - t1062;
t1030 = mrSges(2,1) * t1003 - mrSges(2,2) * t1004 + Ifges(2,3) * qJDD(1) + pkin(1) * t837 + t1012 * t827 + t831 * t1043 + t829 * t1044 + t843 * t1056;
t825 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1003 + Ifges(2,5) * qJDD(1) - t1022 * Ifges(2,6) - t1016 * t831 + t1020 * t829 + (-t1011 * t836 - t1012 * t837) * pkin(8);
t824 = mrSges(2,1) * g(3) + mrSges(2,3) * t1004 + t1022 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t836 - t1011 * t827 + (pkin(8) * t843 + t1016 * t829 + t1020 * t831) * t1012;
t1 = [-m(1) * g(1) + t1036; -m(1) * g(2) + t1054; (-m(1) - m(2)) * g(3) + t836; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1054 - t1017 * t824 + t1021 * t825; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1036 + t1017 * t825 + t1021 * t824; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1030; t1030; t827; t1062; t859; t1024; t1029;];
tauJB  = t1;
