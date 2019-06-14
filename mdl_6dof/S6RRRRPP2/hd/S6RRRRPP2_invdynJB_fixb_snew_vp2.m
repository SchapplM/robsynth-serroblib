% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:01:57
% EndTime: 2019-05-07 18:02:08
% DurationCPUTime: 7.00s
% Computational Cost: add. (101246->343), mult. (202829->405), div. (0->0), fcn. (139752->8), ass. (0->134)
t1002 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t981 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t980 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t1001 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t979 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t1000 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t951 = sin(qJ(3));
t952 = sin(qJ(2));
t954 = cos(qJ(2));
t993 = cos(qJ(3));
t924 = (t951 * t954 + t952 * t993) * qJD(1);
t982 = qJD(1) * qJD(2);
t932 = qJDD(1) * t952 + t954 * t982;
t933 = qJDD(1) * t954 - t952 * t982;
t893 = -qJD(3) * t924 - t932 * t951 + t993 * t933;
t983 = qJD(1) * t954;
t984 = qJD(1) * t952;
t923 = -t951 * t984 + t993 * t983;
t894 = t923 * qJD(3) + t932 * t993 + t951 * t933;
t937 = qJD(2) * pkin(2) - pkin(8) * t984;
t949 = t954 ^ 2;
t956 = qJD(1) ^ 2;
t953 = sin(qJ(1));
t955 = cos(qJ(1));
t938 = t953 * g(1) - t955 * g(2);
t966 = -qJDD(1) * pkin(1) - t938;
t895 = -t933 * pkin(2) + t937 * t984 + (-pkin(8) * t949 - pkin(7)) * t956 + t966;
t947 = qJD(2) + qJD(3);
t839 = (-t923 * t947 - t894) * pkin(9) + (t924 * t947 - t893) * pkin(3) + t895;
t939 = -g(1) * t955 - g(2) * t953;
t926 = -pkin(1) * t956 + qJDD(1) * pkin(7) + t939;
t987 = t952 * t926;
t990 = pkin(2) * t956;
t883 = qJDD(2) * pkin(2) - t932 * pkin(8) - t987 + (pkin(8) * t982 + t952 * t990 - g(3)) * t954;
t913 = -g(3) * t952 + t954 * t926;
t884 = pkin(8) * t933 - qJD(2) * t937 - t949 * t990 + t913;
t846 = t951 * t883 + t993 * t884;
t908 = -pkin(3) * t923 - pkin(9) * t924;
t945 = t947 ^ 2;
t946 = qJDD(2) + qJDD(3);
t843 = -pkin(3) * t945 + pkin(9) * t946 + t908 * t923 + t846;
t950 = sin(qJ(4));
t992 = cos(qJ(4));
t836 = t839 * t992 - t950 * t843;
t910 = t950 * t924 - t947 * t992;
t911 = t924 * t992 + t950 * t947;
t878 = pkin(4) * t910 - qJ(5) * t911;
t892 = qJDD(4) - t893;
t919 = qJD(4) - t923;
t918 = t919 ^ 2;
t834 = -t892 * pkin(4) - t918 * qJ(5) + t911 * t878 + qJDD(5) - t836;
t896 = -mrSges(6,2) * t910 + mrSges(6,3) * t919;
t999 = -m(6) * t834 + t892 * mrSges(6,1) + t919 * t896;
t854 = -t910 * qJD(4) + t894 * t992 + t950 * t946;
t845 = t993 * t883 - t951 * t884;
t964 = t946 * pkin(3) + t945 * pkin(9) - t924 * t908 + t845;
t988 = t910 * t919;
t998 = (-t854 + t988) * qJ(5) - t964;
t907 = -mrSges(4,1) * t923 + mrSges(4,2) * t924;
t915 = mrSges(4,1) * t947 - mrSges(4,3) * t924;
t897 = mrSges(7,2) * t919 + mrSges(7,3) * t910;
t898 = -mrSges(5,2) * t919 - mrSges(5,3) * t910;
t995 = -0.2e1 * t911;
t827 = qJD(6) * t995 + (-t854 - t988) * qJ(6) + (t910 * t911 - t892) * pkin(5) + t834;
t880 = -mrSges(7,1) * t910 + mrSges(7,2) * t911;
t968 = -m(7) * t827 + t854 * mrSges(7,3) + t911 * t880;
t879 = mrSges(6,1) * t910 - mrSges(6,3) * t911;
t985 = -mrSges(5,1) * t910 - mrSges(5,2) * t911 - t879;
t989 = -mrSges(5,3) - mrSges(6,2);
t814 = m(5) * t836 + (t897 + t898) * t919 + t985 * t911 + (mrSges(5,1) + mrSges(7,1)) * t892 + t989 * t854 + t968 + t999;
t837 = t950 * t839 + t992 * t843;
t853 = t911 * qJD(4) + t950 * t894 - t946 * t992;
t900 = -mrSges(7,1) * t919 - mrSges(7,3) * t911;
t901 = mrSges(5,1) * t919 - mrSges(5,3) * t911;
t994 = 2 * qJD(5);
t833 = -pkin(4) * t918 + t892 * qJ(5) - t910 * t878 + t919 * t994 + t837;
t902 = -mrSges(6,1) * t919 + mrSges(6,2) * t911;
t899 = -pkin(5) * t919 - qJ(6) * t911;
t909 = t910 ^ 2;
t829 = -pkin(5) * t909 + t853 * qJ(6) + 0.2e1 * qJD(6) * t910 + t899 * t919 + t833;
t977 = m(7) * t829 + t853 * mrSges(7,3) + t910 * t880;
t965 = m(6) * t833 + t892 * mrSges(6,3) + t919 * t902 + t977;
t819 = m(5) * t837 + (t900 - t901) * t919 + t985 * t910 + (-mrSges(5,2) + mrSges(7,2)) * t892 + t989 * t853 + t965;
t970 = -t814 * t950 + t992 * t819;
t809 = m(4) * t846 - mrSges(4,2) * t946 + mrSges(4,3) * t893 + t907 * t923 - t915 * t947 + t970;
t914 = -mrSges(4,2) * t947 + mrSges(4,3) * t923;
t831 = -t909 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t853 + (-pkin(4) * t919 + t899 + t994) * t911 - t998;
t826 = m(7) * t831 - t853 * mrSges(7,1) + t854 * mrSges(7,2) - t910 * t897 + t911 * t900;
t835 = qJD(5) * t995 + (t911 * t919 + t853) * pkin(4) + t998;
t824 = m(6) * t835 + t853 * mrSges(6,1) - t854 * mrSges(6,3) + t896 * t910 - t911 * t902 - t826;
t958 = m(5) * t964 - t853 * mrSges(5,1) - t854 * mrSges(5,2) - t910 * t898 - t901 * t911 - t824;
t816 = m(4) * t845 + mrSges(4,1) * t946 - mrSges(4,3) * t894 - t907 * t924 + t914 * t947 + t958;
t801 = t951 * t809 + t993 * t816;
t912 = -t954 * g(3) - t987;
t921 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t952 + Ifges(3,2) * t954) * qJD(1);
t922 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t952 + Ifges(3,4) * t954) * qJD(1);
t974 = t1002 * t911 - t981 * t910 + t980 * t919;
t976 = -t1000 * t919 + t979 * t910 - t980 * t911;
t795 = mrSges(5,1) * t964 + mrSges(5,3) * t837 - mrSges(6,1) * t835 + mrSges(6,2) * t833 + mrSges(7,1) * t831 - mrSges(7,3) * t829 + pkin(5) * t826 - qJ(6) * t977 - pkin(4) * t824 + (-qJ(6) * t900 + t974) * t919 + t976 * t911 + (-mrSges(7,2) * qJ(6) + t979) * t892 + t981 * t854 + t1001 * t853;
t825 = -t892 * mrSges(7,1) - t919 * t897 - t968;
t975 = t1001 * t910 + t981 * t911 + t979 * t919;
t803 = -mrSges(5,2) * t964 + mrSges(6,2) * t834 + mrSges(7,2) * t831 - mrSges(5,3) * t836 - mrSges(6,3) * t835 - mrSges(7,3) * t827 - qJ(5) * t824 - qJ(6) * t825 + t1002 * t854 - t981 * t853 + t980 * t892 + t976 * t910 - t975 * t919;
t904 = Ifges(4,4) * t924 + Ifges(4,2) * t923 + Ifges(4,6) * t947;
t905 = Ifges(4,1) * t924 + Ifges(4,4) * t923 + Ifges(4,5) * t947;
t961 = -mrSges(4,1) * t845 + mrSges(4,2) * t846 - Ifges(4,5) * t894 - Ifges(4,6) * t893 - Ifges(4,3) * t946 - pkin(3) * t958 - pkin(9) * t970 - t992 * t795 - t950 * t803 - t924 * t904 + t923 * t905;
t997 = mrSges(3,1) * t912 - mrSges(3,2) * t913 + Ifges(3,5) * t932 + Ifges(3,6) * t933 + Ifges(3,3) * qJDD(2) + pkin(2) * t801 + (t921 * t952 - t922 * t954) * qJD(1) - t961;
t822 = t854 * mrSges(6,2) + t911 * t879 + t825 - t999;
t996 = -t853 * t979 + t854 * t980 + t1000 * t892 + t910 * t974 + t911 * t975 + mrSges(5,1) * t836 - mrSges(6,1) * t834 - mrSges(7,1) * t827 - mrSges(5,2) * t837 + mrSges(7,2) * t829 + mrSges(6,3) * t833 - pkin(4) * t822 - pkin(5) * t825 + qJ(5) * (-t853 * mrSges(6,2) + t892 * mrSges(7,2) - t910 * t879 + t919 * t900 + t965);
t931 = (-mrSges(3,1) * t954 + mrSges(3,2) * t952) * qJD(1);
t936 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t983;
t799 = m(3) * t912 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t932 + qJD(2) * t936 - t931 * t984 + t801;
t935 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t984;
t971 = t993 * t809 - t816 * t951;
t800 = m(3) * t913 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t933 - qJD(2) * t935 + t931 * t983 + t971;
t972 = -t799 * t952 + t954 * t800;
t790 = m(2) * t939 - mrSges(2,1) * t956 - qJDD(1) * mrSges(2,2) + t972;
t925 = -t956 * pkin(7) + t966;
t811 = t992 * t814 + t950 * t819;
t962 = m(4) * t895 - t893 * mrSges(4,1) + mrSges(4,2) * t894 - t923 * t914 + t915 * t924 + t811;
t959 = -m(3) * t925 + t933 * mrSges(3,1) - mrSges(3,2) * t932 - t935 * t984 + t936 * t983 - t962;
t805 = m(2) * t938 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t956 + t959;
t986 = t953 * t790 + t955 * t805;
t792 = t954 * t799 + t952 * t800;
t973 = t955 * t790 - t805 * t953;
t903 = Ifges(4,5) * t924 + Ifges(4,6) * t923 + Ifges(4,3) * t947;
t787 = mrSges(4,2) * t895 - mrSges(4,3) * t845 + Ifges(4,1) * t894 + Ifges(4,4) * t893 + Ifges(4,5) * t946 - pkin(9) * t811 - t950 * t795 + t803 * t992 + t923 * t903 - t947 * t904;
t793 = -mrSges(4,1) * t895 + mrSges(4,3) * t846 + Ifges(4,4) * t894 + Ifges(4,2) * t893 + Ifges(4,6) * t946 - pkin(3) * t811 - t924 * t903 + t947 * t905 - t996;
t920 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t952 + Ifges(3,6) * t954) * qJD(1);
t783 = -mrSges(3,1) * t925 + mrSges(3,3) * t913 + Ifges(3,4) * t932 + Ifges(3,2) * t933 + Ifges(3,6) * qJDD(2) - pkin(2) * t962 + pkin(8) * t971 + qJD(2) * t922 + t951 * t787 + t793 * t993 - t920 * t984;
t786 = mrSges(3,2) * t925 - mrSges(3,3) * t912 + Ifges(3,1) * t932 + Ifges(3,4) * t933 + Ifges(3,5) * qJDD(2) - pkin(8) * t801 - qJD(2) * t921 + t787 * t993 - t951 * t793 + t920 * t983;
t963 = mrSges(2,1) * t938 - mrSges(2,2) * t939 + Ifges(2,3) * qJDD(1) + pkin(1) * t959 + pkin(7) * t972 + t954 * t783 + t952 * t786;
t784 = mrSges(2,1) * g(3) + mrSges(2,3) * t939 + t956 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t792 - t997;
t781 = -mrSges(2,2) * g(3) - mrSges(2,3) * t938 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t956 - pkin(7) * t792 - t783 * t952 + t786 * t954;
t1 = [-m(1) * g(1) + t973; -m(1) * g(2) + t986; (-m(1) - m(2)) * g(3) + t792; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t986 + t955 * t781 - t953 * t784; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t973 + t953 * t781 + t955 * t784; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t963; t963; t997; -t961; t996; t822; t826;];
tauJB  = t1;
