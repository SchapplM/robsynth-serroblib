% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:45:11
% EndTime: 2019-05-07 07:45:28
% DurationCPUTime: 8.29s
% Computational Cost: add. (87915->343), mult. (180945->402), div. (0->0), fcn. (121303->8), ass. (0->139)
t1036 = Ifges(6,1) + Ifges(7,1);
t1015 = Ifges(6,4) - Ifges(7,5);
t1030 = Ifges(7,4) + Ifges(6,5);
t1035 = Ifges(4,2) + Ifges(5,3);
t1034 = Ifges(6,2) + Ifges(7,3);
t1028 = Ifges(4,6) - Ifges(5,5);
t1014 = -Ifges(5,6) - Ifges(4,4);
t1027 = Ifges(6,6) - Ifges(7,6);
t1029 = Ifges(4,5) - Ifges(5,4);
t1031 = Ifges(4,1) + Ifges(5,2);
t980 = cos(qJ(2));
t1004 = qJD(1) * t980;
t978 = sin(qJ(2));
t1005 = qJD(1) * t978;
t1019 = cos(qJ(3));
t977 = sin(qJ(3));
t947 = -t1019 * t1004 + t977 * t1005;
t948 = (t1019 * t978 + t977 * t980) * qJD(1);
t973 = qJD(2) + qJD(3);
t1006 = t1014 * t947 + t1029 * t973 + t1031 * t948;
t1018 = cos(qJ(5));
t1022 = t1014 * t948 - t1028 * t973 + t1035 * t947;
t1026 = Ifges(4,3) + Ifges(5,1);
t934 = t947 * mrSges(5,1) - t973 * mrSges(5,3);
t972 = qJDD(2) + qJDD(3);
t1002 = qJD(1) * qJD(2);
t956 = t978 * qJDD(1) + t980 * t1002;
t957 = t980 * qJDD(1) - t978 * t1002;
t907 = t948 * qJD(3) - t1019 * t957 + t977 * t956;
t936 = t948 * pkin(4) - t973 * pkin(9);
t943 = t947 ^ 2;
t1013 = t947 * t973;
t1032 = -2 * qJD(4);
t908 = -t947 * qJD(3) + t1019 * t956 + t977 * t957;
t961 = qJD(2) * pkin(2) - pkin(8) * t1005;
t975 = t980 ^ 2;
t982 = qJD(1) ^ 2;
t979 = sin(qJ(1));
t981 = cos(qJ(1));
t962 = t979 * g(1) - t981 * g(2);
t994 = -qJDD(1) * pkin(1) - t962;
t909 = -t957 * pkin(2) + t961 * t1005 + (-pkin(8) * t975 - pkin(7)) * t982 + t994;
t987 = (-t908 + t1013) * qJ(4) + t909 + (t973 * pkin(3) + t1032) * t948;
t847 = -t943 * pkin(4) - t948 * t936 + (pkin(3) + pkin(9)) * t907 + t987;
t963 = -t981 * g(1) - t979 * g(2);
t950 = -t982 * pkin(1) + qJDD(1) * pkin(7) + t963;
t1012 = t978 * t950;
t1017 = pkin(2) * t982;
t888 = qJDD(2) * pkin(2) - t956 * pkin(8) - t1012 + (pkin(8) * t1002 + t978 * t1017 - g(3)) * t980;
t931 = -t978 * g(3) + t980 * t950;
t889 = t957 * pkin(8) - qJD(2) * t961 - t975 * t1017 + t931;
t858 = t1019 * t888 - t977 * t889;
t923 = t947 * pkin(3) - t948 * qJ(4);
t971 = t973 ^ 2;
t856 = -t972 * pkin(3) - t971 * qJ(4) + t948 * t923 + qJDD(4) - t858;
t849 = (t947 * t948 - t972) * pkin(9) + (t908 + t1013) * pkin(4) + t856;
t976 = sin(qJ(5));
t843 = t1018 * t847 + t976 * t849;
t928 = -t1018 * t947 + t976 * t973;
t929 = t1018 * t973 + t976 * t947;
t883 = t928 * pkin(5) - t929 * qJ(6);
t906 = qJDD(5) + t908;
t942 = qJD(5) + t948;
t941 = t942 ^ 2;
t839 = -t941 * pkin(5) + t906 * qJ(6) + 0.2e1 * qJD(6) * t942 - t928 * t883 + t843;
t913 = -t942 * mrSges(7,1) + t929 * mrSges(7,2);
t1001 = m(7) * t839 + t906 * mrSges(7,3) + t942 * t913;
t884 = t928 * mrSges(7,1) - t929 * mrSges(7,3);
t1008 = -t928 * mrSges(6,1) - t929 * mrSges(6,2) - t884;
t1016 = -mrSges(6,3) - mrSges(7,2);
t868 = t929 * qJD(5) - t1018 * t907 + t976 * t972;
t912 = t942 * mrSges(6,1) - t929 * mrSges(6,3);
t828 = m(6) * t843 - t906 * mrSges(6,2) + t1008 * t928 + t1016 * t868 - t942 * t912 + t1001;
t842 = t1018 * t849 - t976 * t847;
t869 = -t928 * qJD(5) + t1018 * t972 + t976 * t907;
t911 = -t942 * mrSges(6,2) - t928 * mrSges(6,3);
t840 = -t906 * pkin(5) - t941 * qJ(6) + t929 * t883 + qJDD(6) - t842;
t910 = -t928 * mrSges(7,2) + t942 * mrSges(7,3);
t996 = -m(7) * t840 + t906 * mrSges(7,1) + t942 * t910;
t830 = m(6) * t842 + t906 * mrSges(6,1) + t1008 * t929 + t1016 * t869 + t942 * t911 + t996;
t823 = t1018 * t830 + t976 * t828;
t925 = -t947 * mrSges(5,2) - t948 * mrSges(5,3);
t991 = m(5) * t856 + t908 * mrSges(5,1) + t948 * t925 + t823;
t819 = t972 * mrSges(5,2) + t973 * t934 + t991;
t1025 = Ifges(6,3) + Ifges(7,2);
t1009 = -t1025 * t942 + t1027 * t928 - t1030 * t929;
t1023 = -t1015 * t928 + t1030 * t942 + t1036 * t929;
t859 = t1019 * t889 + t977 * t888;
t854 = t971 * pkin(3) - t972 * qJ(4) + t1032 * t973 + t947 * t923 - t859;
t851 = -t907 * pkin(4) - t943 * pkin(9) + t973 * t936 - t854;
t845 = -0.2e1 * qJD(6) * t929 + (t928 * t942 - t869) * qJ(6) + (t929 * t942 + t868) * pkin(5) + t851;
t836 = m(7) * t845 + t868 * mrSges(7,1) - t869 * mrSges(7,3) + t928 * t910 - t929 * t913;
t820 = -mrSges(6,1) * t851 - mrSges(7,1) * t845 + mrSges(7,2) * t839 + mrSges(6,3) * t843 - pkin(5) * t836 + t1009 * t929 + t1015 * t869 + t1023 * t942 + t1027 * t906 - t1034 * t868;
t1024 = -t1015 * t929 - t1027 * t942 + t1034 * t928;
t822 = mrSges(6,2) * t851 + mrSges(7,2) * t840 - mrSges(6,3) * t842 - mrSges(7,3) * t845 - qJ(6) * t836 + t1009 * t928 - t1015 * t868 + t1024 * t942 + t1030 * t906 + t1036 * t869;
t935 = t948 * mrSges(5,1) + t973 * mrSges(5,2);
t990 = m(6) * t851 + t868 * mrSges(6,1) + t869 * mrSges(6,2) + t928 * t911 + t929 * t912 + t836;
t988 = -m(5) * t854 + t972 * mrSges(5,3) + t973 * t935 + t990;
t1033 = -mrSges(4,2) * t859 - mrSges(5,3) * t854 - pkin(3) * t819 - pkin(9) * t823 + t1018 * t822 - t976 * t820 + qJ(4) * (-t947 * t925 + t988) + mrSges(5,2) * t856 + mrSges(4,1) * t858 + t1026 * t972 - t1022 * t948 + t1029 * t908 + (-qJ(4) * mrSges(5,1) - t1028) * t907 + t1006 * t947;
t924 = t947 * mrSges(4,1) + t948 * mrSges(4,2);
t932 = -t973 * mrSges(4,2) - t947 * mrSges(4,3);
t816 = m(4) * t858 - t908 * mrSges(4,3) - t948 * t924 + (t932 - t934) * t973 + (mrSges(4,1) - mrSges(5,2)) * t972 - t991;
t933 = t973 * mrSges(4,1) - t948 * mrSges(4,3);
t826 = t988 + (-t924 - t925) * t947 + (-mrSges(4,3) - mrSges(5,1)) * t907 + m(4) * t859 - t972 * mrSges(4,2) - t973 * t933;
t811 = t1019 * t816 + t977 * t826;
t930 = -t980 * g(3) - t1012;
t945 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t978 + Ifges(3,2) * t980) * qJD(1);
t946 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t978 + Ifges(3,4) * t980) * qJD(1);
t1021 = mrSges(3,1) * t930 - mrSges(3,2) * t931 + Ifges(3,5) * t956 + Ifges(3,6) * t957 + Ifges(3,3) * qJDD(2) + pkin(2) * t811 + (t978 * t945 - t980 * t946) * qJD(1) + t1033;
t955 = (-mrSges(3,1) * t980 + mrSges(3,2) * t978) * qJD(1);
t960 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1004;
t809 = m(3) * t930 + qJDD(2) * mrSges(3,1) - t956 * mrSges(3,3) + qJD(2) * t960 - t955 * t1005 + t811;
t959 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1005;
t997 = t1019 * t826 - t977 * t816;
t810 = m(3) * t931 - qJDD(2) * mrSges(3,2) + t957 * mrSges(3,3) - qJD(2) * t959 + t955 * t1004 + t997;
t998 = -t978 * t809 + t980 * t810;
t803 = m(2) * t963 - t982 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t998;
t949 = -t982 * pkin(7) + t994;
t1010 = t1018 * t828 - t976 * t830;
t853 = t907 * pkin(3) + t987;
t817 = m(5) * t853 - t907 * mrSges(5,2) - t908 * mrSges(5,3) - t947 * t934 - t948 * t935 + t1010;
t989 = m(4) * t909 + t907 * mrSges(4,1) + t908 * mrSges(4,2) + t947 * t932 + t948 * t933 + t817;
t984 = -m(3) * t949 + t957 * mrSges(3,1) - t956 * mrSges(3,2) + t960 * t1004 - t959 * t1005 - t989;
t813 = m(2) * t962 + qJDD(1) * mrSges(2,1) - t982 * mrSges(2,2) + t984;
t1011 = t979 * t803 + t981 * t813;
t805 = t980 * t809 + t978 * t810;
t1007 = -t1026 * t973 + t1028 * t947 - t1029 * t948;
t999 = t981 * t803 - t979 * t813;
t799 = -mrSges(4,1) * t909 - mrSges(5,1) * t854 + mrSges(5,2) * t853 + mrSges(4,3) * t859 - pkin(3) * t817 + pkin(4) * t990 - pkin(9) * t1010 + t1006 * t973 + t1007 * t948 - t1014 * t908 - t1018 * t820 + t1028 * t972 - t1035 * t907 - t976 * t822;
t835 = t869 * mrSges(7,2) + t929 * t884 - t996;
t986 = mrSges(6,1) * t842 - mrSges(7,1) * t840 - mrSges(6,2) * t843 + mrSges(7,3) * t839 - pkin(5) * t835 + qJ(6) * t1001 - t1024 * t929 + (-qJ(6) * t884 + t1023) * t928 + t1025 * t906 + t1030 * t869 + (-qJ(6) * mrSges(7,2) - t1027) * t868;
t800 = mrSges(5,1) * t856 + mrSges(4,2) * t909 - mrSges(4,3) * t858 - mrSges(5,3) * t853 + pkin(4) * t823 - qJ(4) * t817 + t1007 * t947 + t1014 * t907 + t1022 * t973 + t1029 * t972 + t1031 * t908 + t986;
t944 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t978 + Ifges(3,6) * t980) * qJD(1);
t795 = -mrSges(3,1) * t949 + mrSges(3,3) * t931 + Ifges(3,4) * t956 + Ifges(3,2) * t957 + Ifges(3,6) * qJDD(2) - pkin(2) * t989 + pkin(8) * t997 + qJD(2) * t946 - t944 * t1005 + t1019 * t799 + t977 * t800;
t797 = mrSges(3,2) * t949 - mrSges(3,3) * t930 + Ifges(3,1) * t956 + Ifges(3,4) * t957 + Ifges(3,5) * qJDD(2) - pkin(8) * t811 - qJD(2) * t945 + t944 * t1004 + t1019 * t800 - t977 * t799;
t992 = mrSges(2,1) * t962 - mrSges(2,2) * t963 + Ifges(2,3) * qJDD(1) + pkin(1) * t984 + pkin(7) * t998 + t980 * t795 + t978 * t797;
t798 = mrSges(2,1) * g(3) + mrSges(2,3) * t963 + t982 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t805 - t1021;
t793 = -mrSges(2,2) * g(3) - mrSges(2,3) * t962 + Ifges(2,5) * qJDD(1) - t982 * Ifges(2,6) - pkin(7) * t805 - t978 * t795 + t980 * t797;
t1 = [-m(1) * g(1) + t999; -m(1) * g(2) + t1011; (-m(1) - m(2)) * g(3) + t805; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1011 + t981 * t793 - t979 * t798; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t999 + t979 * t793 + t981 * t798; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t992; t992; t1021; t1033; t819; t986; t835;];
tauJB  = t1;
