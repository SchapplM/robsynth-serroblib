% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR8
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
% Datum: 2019-05-06 11:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:18:41
% EndTime: 2019-05-06 11:18:53
% DurationCPUTime: 11.68s
% Computational Cost: add. (173139->362), mult. (381664->439), div. (0->0), fcn. (258286->10), ass. (0->146)
t1037 = -2 * qJD(3);
t1036 = -2 * qJD(4);
t1035 = Ifges(4,1) + Ifges(5,1);
t1028 = Ifges(4,4) - Ifges(5,5);
t1034 = Ifges(4,5) + Ifges(5,4);
t1033 = Ifges(4,2) + Ifges(5,3);
t1032 = Ifges(5,2) + Ifges(4,3);
t1031 = Ifges(4,6) - Ifges(5,6);
t988 = sin(qJ(2));
t1018 = qJD(1) * t988;
t1025 = cos(pkin(10));
t985 = sin(pkin(10));
t957 = -t1025 * qJD(2) + t985 * t1018;
t992 = cos(qJ(2));
t978 = t992 * qJD(1);
t1014 = t957 * t978;
t989 = sin(qJ(1));
t993 = cos(qJ(1));
t971 = -t993 * g(1) - t989 * g(2);
t995 = qJD(1) ^ 2;
t951 = -t995 * pkin(1) + qJDD(1) * pkin(7) + t971;
t927 = -t992 * g(3) - t988 * t951;
t963 = (-pkin(2) * t992 - qJ(3) * t988) * qJD(1);
t994 = qJD(2) ^ 2;
t905 = -qJDD(2) * pkin(2) - t994 * qJ(3) + t963 * t1018 + qJDD(3) - t927;
t1015 = qJD(1) * qJD(2);
t1013 = t992 * t1015;
t965 = t988 * qJDD(1) + t1013;
t936 = -t1025 * qJDD(2) + t985 * t965;
t937 = t985 * qJDD(2) + t1025 * t965;
t958 = t985 * qJD(2) + t1025 * t1018;
t883 = t905 + (-t1014 - t937) * qJ(4) + (-t958 * t978 + t936) * pkin(3) + t958 * t1036;
t970 = t989 * g(1) - t993 * g(2);
t950 = -qJDD(1) * pkin(1) - t995 * pkin(7) - t970;
t975 = t988 * t1015;
t966 = t992 * qJDD(1) - t975;
t902 = (-t965 - t1013) * qJ(3) + (-t966 + t975) * pkin(2) + t950;
t928 = -t988 * g(3) + t992 * t951;
t906 = -t994 * pkin(2) + qJDD(2) * qJ(3) + t963 * t978 + t928;
t884 = t1025 * t902 + t958 * t1037 - t985 * t906;
t1024 = t992 ^ 2 * t995;
t924 = t957 * pkin(3) - t958 * qJ(4);
t877 = t966 * pkin(3) - qJ(4) * t1024 + t958 * t924 + qJDD(4) - t884;
t864 = (-t937 + t1014) * pkin(8) + (t957 * t958 + t966) * pkin(4) + t877;
t885 = t1025 * t906 + t957 * t1037 + t985 * t902;
t876 = -pkin(3) * t1024 - t966 * qJ(4) + t978 * t1036 - t957 * t924 + t885;
t938 = pkin(4) * t978 - t958 * pkin(8);
t955 = t957 ^ 2;
t866 = -t955 * pkin(4) + t936 * pkin(8) - t938 * t978 + t876;
t987 = sin(qJ(5));
t991 = cos(qJ(5));
t858 = t991 * t864 - t987 * t866;
t922 = t991 * t957 - t987 * t958;
t893 = t922 * qJD(5) + t987 * t936 + t991 * t937;
t923 = t987 * t957 + t991 * t958;
t962 = qJDD(5) + t966;
t973 = t978 + qJD(5);
t855 = (t922 * t973 - t893) * pkin(9) + (t922 * t923 + t962) * pkin(5) + t858;
t859 = t987 * t864 + t991 * t866;
t892 = -t923 * qJD(5) + t991 * t936 - t987 * t937;
t911 = t973 * pkin(5) - t923 * pkin(9);
t921 = t922 ^ 2;
t856 = -t921 * pkin(5) + t892 * pkin(9) - t973 * t911 + t859;
t986 = sin(qJ(6));
t990 = cos(qJ(6));
t853 = t990 * t855 - t986 * t856;
t898 = t990 * t922 - t986 * t923;
t871 = t898 * qJD(6) + t986 * t892 + t990 * t893;
t899 = t986 * t922 + t990 * t923;
t882 = -t898 * mrSges(7,1) + mrSges(7,2) * t899;
t972 = qJD(6) + t973;
t887 = -t972 * mrSges(7,2) + t898 * mrSges(7,3);
t953 = qJDD(6) + t962;
t850 = m(7) * t853 + t953 * mrSges(7,1) - t871 * mrSges(7,3) - t899 * t882 + t972 * t887;
t854 = t986 * t855 + t990 * t856;
t870 = -t899 * qJD(6) + t990 * t892 - t986 * t893;
t888 = t972 * mrSges(7,1) - t899 * mrSges(7,3);
t851 = m(7) * t854 - t953 * mrSges(7,2) + t870 * mrSges(7,3) + t898 * t882 - t972 * t888;
t840 = t990 * t850 + t986 * t851;
t900 = -t922 * mrSges(6,1) + t923 * mrSges(6,2);
t907 = -t973 * mrSges(6,2) + t922 * mrSges(6,3);
t837 = m(6) * t858 + t962 * mrSges(6,1) - t893 * mrSges(6,3) - t923 * t900 + t973 * t907 + t840;
t1009 = -t986 * t850 + t990 * t851;
t908 = t973 * mrSges(6,1) - t923 * mrSges(6,3);
t838 = m(6) * t859 - t962 * mrSges(6,2) + t892 * mrSges(6,3) + t922 * t900 - t973 * t908 + t1009;
t1010 = -t987 * t837 + t991 * t838;
t1020 = t1028 * t957 + t1034 * t978 - t1035 * t958;
t1021 = t1031 * t957 + t1032 * t978 - t1034 * t958;
t874 = -t936 * pkin(4) - t955 * pkin(8) + t958 * t938 - t883;
t861 = -t892 * pkin(5) - t921 * pkin(9) + t923 * t911 + t874;
t1004 = m(7) * t861 - t870 * mrSges(7,1) + t871 * mrSges(7,2) - t898 * t887 + t899 * t888;
t878 = Ifges(7,5) * t899 + Ifges(7,6) * t898 + Ifges(7,3) * t972;
t880 = Ifges(7,1) * t899 + Ifges(7,4) * t898 + Ifges(7,5) * t972;
t841 = -mrSges(7,1) * t861 + mrSges(7,3) * t854 + Ifges(7,4) * t871 + Ifges(7,2) * t870 + Ifges(7,6) * t953 - t899 * t878 + t972 * t880;
t879 = Ifges(7,4) * t899 + Ifges(7,2) * t898 + Ifges(7,6) * t972;
t842 = mrSges(7,2) * t861 - mrSges(7,3) * t853 + Ifges(7,1) * t871 + Ifges(7,4) * t870 + Ifges(7,5) * t953 + t898 * t878 - t972 * t879;
t894 = Ifges(6,5) * t923 + Ifges(6,6) * t922 + Ifges(6,3) * t973;
t896 = Ifges(6,1) * t923 + Ifges(6,4) * t922 + Ifges(6,5) * t973;
t829 = -mrSges(6,1) * t874 + mrSges(6,3) * t859 + Ifges(6,4) * t893 + Ifges(6,2) * t892 + Ifges(6,6) * t962 - pkin(5) * t1004 + pkin(9) * t1009 + t990 * t841 + t986 * t842 - t923 * t894 + t973 * t896;
t895 = Ifges(6,4) * t923 + Ifges(6,2) * t922 + Ifges(6,6) * t973;
t830 = mrSges(6,2) * t874 - mrSges(6,3) * t858 + Ifges(6,1) * t893 + Ifges(6,4) * t892 + Ifges(6,5) * t962 - pkin(9) * t840 - t986 * t841 + t990 * t842 + t922 * t894 - t973 * t895;
t933 = -t957 * mrSges(5,2) - mrSges(5,3) * t978;
t935 = mrSges(5,1) * t978 + t958 * mrSges(5,2);
t999 = -m(6) * t874 + t892 * mrSges(6,1) - t893 * mrSges(6,2) + t922 * t907 - t923 * t908 - t1004;
t846 = m(5) * t883 + t936 * mrSges(5,1) - t937 * mrSges(5,3) + t957 * t933 - t958 * t935 + t999;
t814 = -mrSges(4,1) * t905 - mrSges(5,1) * t883 + mrSges(5,2) * t876 + mrSges(4,3) * t885 - pkin(3) * t846 - pkin(4) * t999 - pkin(8) * t1010 + t1020 * t978 + t1021 * t958 + t1028 * t937 - t1031 * t966 - t1033 * t936 - t991 * t829 - t987 * t830;
t1022 = -t1028 * t958 + t1031 * t978 + t1033 * t957;
t835 = t991 * t837 + t987 * t838;
t815 = mrSges(4,2) * t905 + mrSges(5,2) * t877 - mrSges(4,3) * t884 - mrSges(5,3) * t883 - pkin(8) * t835 - qJ(4) * t846 + t1021 * t957 - t1022 * t978 - t1028 * t936 - t1034 * t966 + t1035 * t937 - t987 * t829 + t991 * t830;
t1006 = m(5) * t876 - t966 * mrSges(5,3) + t1010;
t925 = t957 * mrSges(5,1) - t958 * mrSges(5,3);
t1019 = -t957 * mrSges(4,1) - t958 * mrSges(4,2) - t925;
t1029 = -mrSges(4,3) - mrSges(5,2);
t934 = -mrSges(4,1) * t978 - t958 * mrSges(4,3);
t832 = m(4) * t885 + t966 * mrSges(4,2) + t1019 * t957 + t1029 * t936 + (t934 - t935) * t978 + t1006;
t1001 = -m(5) * t877 - t966 * mrSges(5,1) - t933 * t978 - t835;
t1005 = mrSges(4,2) * t978 - t957 * mrSges(4,3);
t833 = m(4) * t884 - t966 * mrSges(4,1) - t1005 * t978 + t1019 * t958 + t1029 * t937 + t1001;
t828 = t1025 * t832 - t985 * t833;
t845 = m(4) * t905 + t936 * mrSges(4,1) + t937 * mrSges(4,2) + t957 * t1005 + t958 * t934 + t846;
t948 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t988 + Ifges(3,2) * t992) * qJD(1);
t949 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t988 + Ifges(3,4) * t992) * qJD(1);
t1030 = mrSges(3,1) * t927 - mrSges(3,2) * t928 + Ifges(3,5) * t965 + Ifges(3,6) * t966 + Ifges(3,3) * qJDD(2) - pkin(2) * t845 + qJ(3) * t828 + (t948 * t988 - t949 * t992) * qJD(1) + t1025 * t814 + t985 * t815;
t964 = (-mrSges(3,1) * t992 + mrSges(3,2) * t988) * qJD(1);
t968 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1018;
t826 = m(3) * t928 - qJDD(2) * mrSges(3,2) + t966 * mrSges(3,3) - qJD(2) * t968 + t964 * t978 + t828;
t969 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t978;
t844 = m(3) * t927 + qJDD(2) * mrSges(3,1) - t965 * mrSges(3,3) + qJD(2) * t969 - t964 * t1018 - t845;
t1011 = t992 * t826 - t988 * t844;
t818 = m(2) * t971 - t995 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1011;
t827 = t1025 * t833 + t985 * t832;
t998 = -m(3) * t950 + t966 * mrSges(3,1) - t965 * mrSges(3,2) - t968 * t1018 + t969 * t978 - t827;
t822 = m(2) * t970 + qJDD(1) * mrSges(2,1) - t995 * mrSges(2,2) + t998;
t1023 = t989 * t818 + t993 * t822;
t820 = t988 * t826 + t992 * t844;
t1012 = t993 * t818 - t989 * t822;
t947 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t988 + Ifges(3,6) * t992) * qJD(1);
t811 = mrSges(3,2) * t950 - mrSges(3,3) * t927 + Ifges(3,1) * t965 + Ifges(3,4) * t966 + Ifges(3,5) * qJDD(2) - qJ(3) * t827 - qJD(2) * t948 + t1025 * t815 - t985 * t814 + t947 * t978;
t834 = t937 * mrSges(5,2) + t958 * t925 - t1001;
t1002 = mrSges(7,1) * t853 - mrSges(7,2) * t854 + Ifges(7,5) * t871 + Ifges(7,6) * t870 + Ifges(7,3) * t953 + t899 * t879 - t898 * t880;
t997 = mrSges(6,1) * t858 - mrSges(6,2) * t859 + Ifges(6,5) * t893 + Ifges(6,6) * t892 + Ifges(6,3) * t962 + pkin(5) * t840 + t923 * t895 - t922 * t896 + t1002;
t813 = Ifges(3,6) * qJDD(2) + (qJ(4) * t925 + t1020) * t957 + Ifges(3,4) * t965 + qJD(2) * t949 - mrSges(3,1) * t950 + mrSges(3,3) * t928 + (Ifges(3,2) + t1032) * t966 + mrSges(4,2) * t885 - mrSges(5,3) * t876 + mrSges(5,1) * t877 - mrSges(4,1) * t884 + pkin(4) * t835 + pkin(3) * t834 + (qJ(4) * mrSges(5,2) + t1031) * t936 + t1022 * t958 - pkin(2) * t827 + t997 - t1034 * t937 - qJ(4) * (-t935 * t978 + t1006) - t947 * t1018;
t1003 = mrSges(2,1) * t970 - mrSges(2,2) * t971 + Ifges(2,3) * qJDD(1) + pkin(1) * t998 + pkin(7) * t1011 + t988 * t811 + t992 * t813;
t809 = mrSges(2,1) * g(3) + mrSges(2,3) * t971 + t995 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t820 - t1030;
t808 = -mrSges(2,2) * g(3) - mrSges(2,3) * t970 + Ifges(2,5) * qJDD(1) - t995 * Ifges(2,6) - pkin(7) * t820 + t992 * t811 - t988 * t813;
t1 = [-m(1) * g(1) + t1012; -m(1) * g(2) + t1023; (-m(1) - m(2)) * g(3) + t820; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1023 + t993 * t808 - t989 * t809; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1012 + t989 * t808 + t993 * t809; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1003; t1003; t1030; t845; t834; t997; t1002;];
tauJB  = t1;
