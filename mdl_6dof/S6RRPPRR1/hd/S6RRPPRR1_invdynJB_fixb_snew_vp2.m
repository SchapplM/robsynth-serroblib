% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:31:52
% EndTime: 2019-05-06 09:32:05
% DurationCPUTime: 10.97s
% Computational Cost: add. (150203->365), mult. (350699->444), div. (0->0), fcn. (244920->10), ass. (0->147)
t1035 = -2 * qJD(3);
t1034 = Ifges(4,1) + Ifges(5,1);
t1028 = Ifges(4,4) - Ifges(5,5);
t1027 = Ifges(4,5) + Ifges(5,4);
t1033 = -Ifges(4,2) - Ifges(5,3);
t1032 = -Ifges(5,2) - Ifges(4,3);
t1026 = Ifges(4,6) - Ifges(5,6);
t1024 = cos(pkin(10));
t1012 = qJD(1) * qJD(2);
t987 = sin(qJ(1));
t991 = cos(qJ(1));
t966 = -g(1) * t991 - g(2) * t987;
t993 = qJD(1) ^ 2;
t953 = -pkin(1) * t993 + qJDD(1) * pkin(7) + t966;
t986 = sin(qJ(2));
t1022 = t986 * t953;
t990 = cos(qJ(2));
t959 = qJDD(1) * t986 + t990 * t1012;
t899 = qJDD(2) * pkin(2) - t959 * qJ(3) - t1022 + (pkin(2) * t986 * t993 + qJ(3) * t1012 - g(3)) * t990;
t1023 = t990 ^ 2 * t993;
t934 = -g(3) * t986 + t990 * t953;
t960 = qJDD(1) * t990 - t986 * t1012;
t1016 = qJD(1) * t986;
t962 = qJD(2) * pkin(2) - qJ(3) * t1016;
t900 = -pkin(2) * t1023 + qJ(3) * t960 - qJD(2) * t962 + t934;
t982 = sin(pkin(10));
t947 = (t1024 * t986 + t982 * t990) * qJD(1);
t876 = t1024 * t899 + t947 * t1035 - t982 * t900;
t1015 = qJD(1) * t990;
t946 = -t1024 * t1015 + t982 * t1016;
t1014 = qJD(2) * t946;
t921 = pkin(3) * t946 - qJ(4) * t947;
t992 = qJD(2) ^ 2;
t866 = -qJDD(2) * pkin(3) - t992 * qJ(4) + t947 * t921 + qJDD(4) - t876;
t930 = t1024 * t959 + t982 * t960;
t859 = (-t930 - t1014) * pkin(8) + (t946 * t947 - qJDD(2)) * pkin(4) + t866;
t1030 = 2 * qJD(4);
t877 = t1024 * t900 + t946 * t1035 + t982 * t899;
t865 = -pkin(3) * t992 + qJDD(2) * qJ(4) + qJD(2) * t1030 - t946 * t921 + t877;
t929 = -t1024 * t960 + t982 * t959;
t939 = -qJD(2) * pkin(4) - pkin(8) * t947;
t945 = t946 ^ 2;
t861 = -pkin(4) * t945 + pkin(8) * t929 + qJD(2) * t939 + t865;
t985 = sin(qJ(5));
t989 = cos(qJ(5));
t857 = t985 * t859 + t989 * t861;
t915 = t946 * t989 - t947 * t985;
t916 = t946 * t985 + t947 * t989;
t895 = -pkin(5) * t915 - pkin(9) * t916;
t975 = -qJD(2) + qJD(5);
t973 = t975 ^ 2;
t974 = -qJDD(2) + qJDD(5);
t853 = -pkin(5) * t973 + pkin(9) * t974 + t895 * t915 + t857;
t1025 = pkin(3) * qJD(2);
t965 = t987 * g(1) - t991 * g(2);
t952 = -qJDD(1) * pkin(1) - t993 * pkin(7) - t965;
t903 = -t960 * pkin(2) - qJ(3) * t1023 + t962 * t1016 + qJDD(3) + t952;
t998 = t929 * pkin(3) + t903 + (t1014 - t930) * qJ(4);
t863 = -t929 * pkin(4) - t945 * pkin(8) - t998 + (-t1025 + t1030 + t939) * t947;
t884 = -qJD(5) * t916 + t929 * t989 - t930 * t985;
t885 = qJD(5) * t915 + t929 * t985 + t930 * t989;
t854 = t863 + (t916 * t975 - t884) * pkin(5) + (-t915 * t975 - t885) * pkin(9);
t984 = sin(qJ(6));
t988 = cos(qJ(6));
t850 = -t853 * t984 + t854 * t988;
t904 = -t916 * t984 + t975 * t988;
t869 = qJD(6) * t904 + t885 * t988 + t974 * t984;
t883 = qJDD(6) - t884;
t905 = t916 * t988 + t975 * t984;
t886 = -mrSges(7,1) * t904 + mrSges(7,2) * t905;
t908 = qJD(6) - t915;
t887 = -mrSges(7,2) * t908 + mrSges(7,3) * t904;
t846 = m(7) * t850 + mrSges(7,1) * t883 - mrSges(7,3) * t869 - t886 * t905 + t887 * t908;
t851 = t853 * t988 + t854 * t984;
t868 = -qJD(6) * t905 - t885 * t984 + t974 * t988;
t888 = mrSges(7,1) * t908 - mrSges(7,3) * t905;
t847 = m(7) * t851 - mrSges(7,2) * t883 + mrSges(7,3) * t868 + t886 * t904 - t888 * t908;
t1006 = -t846 * t984 + t988 * t847;
t894 = -mrSges(6,1) * t915 + mrSges(6,2) * t916;
t907 = mrSges(6,1) * t975 - mrSges(6,3) * t916;
t834 = m(6) * t857 - mrSges(6,2) * t974 + mrSges(6,3) * t884 + t894 * t915 - t907 * t975 + t1006;
t856 = t859 * t989 - t861 * t985;
t852 = -pkin(5) * t974 - pkin(9) * t973 + t895 * t916 - t856;
t1000 = -m(7) * t852 + t868 * mrSges(7,1) - mrSges(7,2) * t869 + t904 * t887 - t888 * t905;
t906 = -mrSges(6,2) * t975 + mrSges(6,3) * t915;
t842 = m(6) * t856 + mrSges(6,1) * t974 - mrSges(6,3) * t885 - t894 * t916 + t906 * t975 + t1000;
t1007 = t989 * t834 - t985 * t842;
t937 = -qJD(2) * mrSges(5,1) + mrSges(5,2) * t947;
t1003 = m(5) * t865 + qJDD(2) * mrSges(5,3) + qJD(2) * t937 + t1007;
t1018 = t1027 * qJD(2) - t1028 * t946 + t1034 * t947;
t1019 = t1026 * qJD(2) + t1028 * t947 + t1033 * t946;
t922 = mrSges(5,1) * t946 - mrSges(5,3) * t947;
t1017 = -mrSges(4,1) * t946 - mrSges(4,2) * t947 - t922;
t1029 = -mrSges(4,3) - mrSges(5,2);
t936 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t947;
t823 = m(4) * t877 - qJDD(2) * mrSges(4,2) - qJD(2) * t936 + t1017 * t946 + t1029 * t929 + t1003;
t935 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t946;
t827 = t985 * t834 + t989 * t842;
t938 = -mrSges(5,2) * t946 + qJD(2) * mrSges(5,3);
t999 = -m(5) * t866 + qJDD(2) * mrSges(5,1) + qJD(2) * t938 - t827;
t824 = m(4) * t876 + qJDD(2) * mrSges(4,1) + qJD(2) * t935 + t1017 * t947 + t1029 * t930 + t999;
t817 = t1024 * t824 + t982 * t823;
t826 = t930 * mrSges(5,2) + t947 * t922 - t999;
t933 = -t990 * g(3) - t1022;
t949 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t986 + Ifges(3,2) * t990) * qJD(1);
t950 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t986 + Ifges(3,4) * t990) * qJD(1);
t872 = Ifges(7,5) * t905 + Ifges(7,6) * t904 + Ifges(7,3) * t908;
t874 = Ifges(7,1) * t905 + Ifges(7,4) * t904 + Ifges(7,5) * t908;
t840 = -mrSges(7,1) * t852 + mrSges(7,3) * t851 + Ifges(7,4) * t869 + Ifges(7,2) * t868 + Ifges(7,6) * t883 - t872 * t905 + t874 * t908;
t873 = Ifges(7,4) * t905 + Ifges(7,2) * t904 + Ifges(7,6) * t908;
t841 = mrSges(7,2) * t852 - mrSges(7,3) * t850 + Ifges(7,1) * t869 + Ifges(7,4) * t868 + Ifges(7,5) * t883 + t872 * t904 - t873 * t908;
t890 = Ifges(6,4) * t916 + Ifges(6,2) * t915 + Ifges(6,6) * t975;
t891 = Ifges(6,1) * t916 + Ifges(6,4) * t915 + Ifges(6,5) * t975;
t997 = -mrSges(6,1) * t856 + mrSges(6,2) * t857 - Ifges(6,5) * t885 - Ifges(6,6) * t884 - Ifges(6,3) * t974 - pkin(5) * t1000 - pkin(9) * t1006 - t988 * t840 - t984 * t841 - t916 * t890 + t915 * t891;
t1031 = (t949 * t986 - t950 * t990) * qJD(1) + (Ifges(3,3) - t1032) * qJDD(2) + t1018 * t946 + t1019 * t947 - t1026 * t929 + t1027 * t930 + mrSges(3,1) * t933 + mrSges(4,1) * t876 - mrSges(5,1) * t866 - mrSges(3,2) * t934 - mrSges(4,2) * t877 + mrSges(5,3) * t865 + Ifges(3,5) * t959 + Ifges(3,6) * t960 + pkin(2) * t817 - pkin(3) * t826 - pkin(4) * t827 + qJ(4) * (-t929 * mrSges(5,2) - t946 * t922 + t1003) + t997;
t958 = (-mrSges(3,1) * t990 + mrSges(3,2) * t986) * qJD(1);
t964 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1015;
t815 = m(3) * t933 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t959 + qJD(2) * t964 - t958 * t1016 + t817;
t1008 = t1024 * t823 - t824 * t982;
t963 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1016;
t816 = m(3) * t934 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t960 - qJD(2) * t963 + t958 * t1015 + t1008;
t1009 = -t815 * t986 + t990 * t816;
t809 = m(2) * t966 - mrSges(2,1) * t993 - qJDD(1) * mrSges(2,2) + t1009;
t836 = t988 * t846 + t984 * t847;
t1002 = -m(6) * t863 + t884 * mrSges(6,1) - t885 * mrSges(6,2) + t915 * t906 - t916 * t907 - t836;
t871 = (-(2 * qJD(4)) + t1025) * t947 + t998;
t832 = m(5) * t871 + t929 * mrSges(5,1) - t930 * mrSges(5,3) - t947 * t937 + t946 * t938 + t1002;
t831 = m(4) * t903 + t929 * mrSges(4,1) + t930 * mrSges(4,2) + t946 * t935 + t947 * t936 + t832;
t995 = -m(3) * t952 + t960 * mrSges(3,1) - t959 * mrSges(3,2) + t964 * t1015 - t963 * t1016 - t831;
t829 = m(2) * t965 + qJDD(1) * mrSges(2,1) - t993 * mrSges(2,2) + t995;
t1021 = t987 * t809 + t991 * t829;
t811 = t990 * t815 + t986 * t816;
t1020 = t1032 * qJD(2) + t1026 * t946 - t1027 * t947;
t1010 = t991 * t809 - t829 * t987;
t889 = Ifges(6,5) * t916 + Ifges(6,6) * t915 + Ifges(6,3) * t975;
t818 = mrSges(6,2) * t863 - mrSges(6,3) * t856 + Ifges(6,1) * t885 + Ifges(6,4) * t884 + Ifges(6,5) * t974 - pkin(9) * t836 - t840 * t984 + t841 * t988 + t889 * t915 - t890 * t975;
t996 = mrSges(7,1) * t850 - mrSges(7,2) * t851 + Ifges(7,5) * t869 + Ifges(7,6) * t868 + Ifges(7,3) * t883 + t873 * t905 - t874 * t904;
t819 = -mrSges(6,1) * t863 + mrSges(6,3) * t857 + Ifges(6,4) * t885 + Ifges(6,2) * t884 + Ifges(6,6) * t974 - pkin(5) * t836 - t889 * t916 + t891 * t975 - t996;
t805 = -mrSges(4,1) * t903 - mrSges(5,1) * t871 + mrSges(5,2) * t865 + mrSges(4,3) * t877 - pkin(3) * t832 - pkin(4) * t1002 - pkin(8) * t1007 + t1018 * qJD(2) + t1026 * qJDD(2) + t1020 * t947 + t1028 * t930 + t1033 * t929 - t985 * t818 - t989 * t819;
t806 = mrSges(4,2) * t903 + mrSges(5,2) * t866 - mrSges(4,3) * t876 - mrSges(5,3) * t871 - pkin(8) * t827 - qJ(4) * t832 - t1019 * qJD(2) + t1027 * qJDD(2) + t1020 * t946 - t1028 * t929 + t1034 * t930 + t989 * t818 - t985 * t819;
t948 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t986 + Ifges(3,6) * t990) * qJD(1);
t801 = -mrSges(3,1) * t952 + mrSges(3,3) * t934 + Ifges(3,4) * t959 + Ifges(3,2) * t960 + Ifges(3,6) * qJDD(2) - pkin(2) * t831 + qJ(3) * t1008 + qJD(2) * t950 - t1016 * t948 + t1024 * t805 + t982 * t806;
t803 = mrSges(3,2) * t952 - mrSges(3,3) * t933 + Ifges(3,1) * t959 + Ifges(3,4) * t960 + Ifges(3,5) * qJDD(2) - qJ(3) * t817 - qJD(2) * t949 + t1015 * t948 + t1024 * t806 - t982 * t805;
t1001 = mrSges(2,1) * t965 - mrSges(2,2) * t966 + Ifges(2,3) * qJDD(1) + pkin(1) * t995 + pkin(7) * t1009 + t990 * t801 + t986 * t803;
t804 = mrSges(2,1) * g(3) + mrSges(2,3) * t966 + t993 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t811 - t1031;
t799 = -mrSges(2,2) * g(3) - mrSges(2,3) * t965 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t993 - pkin(7) * t811 - t801 * t986 + t803 * t990;
t1 = [-m(1) * g(1) + t1010; -m(1) * g(2) + t1021; (-m(1) - m(2)) * g(3) + t811; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1021 + t991 * t799 - t987 * t804; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1010 + t987 * t799 + t991 * t804; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1001; t1001; t1031; t831; t826; -t997; t996;];
tauJB  = t1;
