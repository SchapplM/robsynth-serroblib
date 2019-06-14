% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:59:05
% EndTime: 2019-05-08 04:59:33
% DurationCPUTime: 21.16s
% Computational Cost: add. (344442->364), mult. (706244->441), div. (0->0), fcn. (516462->10), ass. (0->145)
t990 = Ifges(6,4) + Ifges(7,4);
t1000 = Ifges(6,2) + Ifges(7,2);
t995 = Ifges(6,6) + Ifges(7,6);
t957 = sin(qJ(3));
t962 = cos(qJ(3));
t958 = sin(qJ(2));
t986 = qJD(1) * t958;
t938 = qJD(2) * t957 + t962 * t986;
t963 = cos(qJ(2));
t984 = qJD(1) * qJD(2);
t981 = t963 * t984;
t941 = qJDD(1) * t958 + t981;
t908 = -qJD(3) * t938 + qJDD(2) * t962 - t941 * t957;
t937 = qJD(2) * t962 - t957 * t986;
t909 = qJD(3) * t937 + qJDD(2) * t957 + t941 * t962;
t956 = sin(qJ(4));
t961 = cos(qJ(4));
t912 = t937 * t956 + t938 * t961;
t874 = -qJD(4) * t912 + t908 * t961 - t909 * t956;
t911 = t937 * t961 - t938 * t956;
t875 = qJD(4) * t911 + t908 * t956 + t909 * t961;
t955 = sin(qJ(5));
t960 = cos(qJ(5));
t890 = t911 * t960 - t912 * t955;
t848 = qJD(5) * t890 + t874 * t955 + t875 * t960;
t891 = t911 * t955 + t912 * t960;
t867 = -mrSges(7,1) * t890 + mrSges(7,2) * t891;
t959 = sin(qJ(1));
t964 = cos(qJ(1));
t947 = g(1) * t959 - t964 * g(2);
t966 = qJD(1) ^ 2;
t930 = -qJDD(1) * pkin(1) - pkin(7) * t966 - t947;
t951 = t958 * t984;
t942 = qJDD(1) * t963 - t951;
t895 = (-t941 - t981) * pkin(8) + (-t942 + t951) * pkin(2) + t930;
t948 = -g(1) * t964 - g(2) * t959;
t931 = -pkin(1) * t966 + qJDD(1) * pkin(7) + t948;
t917 = -g(3) * t958 + t963 * t931;
t940 = (-pkin(2) * t963 - pkin(8) * t958) * qJD(1);
t965 = qJD(2) ^ 2;
t985 = qJD(1) * t963;
t898 = -pkin(2) * t965 + qJDD(2) * pkin(8) + t940 * t985 + t917;
t876 = t962 * t895 - t898 * t957;
t936 = qJDD(3) - t942;
t950 = qJD(3) - t985;
t854 = (t937 * t950 - t909) * pkin(9) + (t937 * t938 + t936) * pkin(3) + t876;
t877 = t957 * t895 + t962 * t898;
t918 = pkin(3) * t950 - pkin(9) * t938;
t935 = t937 ^ 2;
t856 = -pkin(3) * t935 + pkin(9) * t908 - t918 * t950 + t877;
t834 = t961 * t854 - t856 * t956;
t932 = qJDD(4) + t936;
t949 = qJD(4) + t950;
t830 = (t911 * t949 - t875) * pkin(10) + (t911 * t912 + t932) * pkin(4) + t834;
t835 = t956 * t854 + t961 * t856;
t901 = pkin(4) * t949 - pkin(10) * t912;
t910 = t911 ^ 2;
t832 = -pkin(4) * t910 + pkin(10) * t874 - t901 * t949 + t835;
t824 = t960 * t830 - t832 * t955;
t926 = qJDD(5) + t932;
t944 = qJD(5) + t949;
t819 = -0.2e1 * qJD(6) * t891 + (t890 * t944 - t848) * qJ(6) + (t890 * t891 + t926) * pkin(5) + t824;
t879 = -mrSges(7,2) * t944 + mrSges(7,3) * t890;
t983 = m(7) * t819 + t926 * mrSges(7,1) + t944 * t879;
t816 = -mrSges(7,3) * t848 - t867 * t891 + t983;
t825 = t955 * t830 + t960 * t832;
t847 = -qJD(5) * t891 + t874 * t960 - t875 * t955;
t881 = pkin(5) * t944 - qJ(6) * t891;
t889 = t890 ^ 2;
t822 = -pkin(5) * t889 + qJ(6) * t847 + 0.2e1 * qJD(6) * t890 - t881 * t944 + t825;
t996 = Ifges(6,5) + Ifges(7,5);
t997 = Ifges(6,1) + Ifges(7,1);
t987 = t990 * t890 + t891 * t997 + t996 * t944;
t993 = t1000 * t890 + t990 * t891 + t995 * t944;
t994 = Ifges(6,3) + Ifges(7,3);
t999 = mrSges(6,1) * t824 + mrSges(7,1) * t819 - mrSges(6,2) * t825 - mrSges(7,2) * t822 + pkin(5) * t816 + t847 * t995 + t848 * t996 - t987 * t890 + t891 * t993 + t926 * t994;
t868 = -mrSges(6,1) * t890 + mrSges(6,2) * t891;
t880 = -mrSges(6,2) * t944 + mrSges(6,3) * t890;
t808 = m(6) * t824 + mrSges(6,1) * t926 + t880 * t944 + (-t867 - t868) * t891 + (-mrSges(6,3) - mrSges(7,3)) * t848 + t983;
t882 = mrSges(7,1) * t944 - mrSges(7,3) * t891;
t883 = mrSges(6,1) * t944 - mrSges(6,3) * t891;
t982 = m(7) * t822 + t847 * mrSges(7,3) + t890 * t867;
t811 = m(6) * t825 + mrSges(6,3) * t847 + t868 * t890 + (-t882 - t883) * t944 + (-mrSges(6,2) - mrSges(7,2)) * t926 + t982;
t806 = t960 * t808 + t955 * t811;
t885 = Ifges(5,4) * t912 + Ifges(5,2) * t911 + Ifges(5,6) * t949;
t886 = Ifges(5,1) * t912 + Ifges(5,4) * t911 + Ifges(5,5) * t949;
t998 = mrSges(5,1) * t834 - mrSges(5,2) * t835 + Ifges(5,5) * t875 + Ifges(5,6) * t874 + Ifges(5,3) * t932 + pkin(4) * t806 + t912 * t885 - t911 * t886 + t999;
t892 = -mrSges(5,1) * t911 + mrSges(5,2) * t912;
t899 = -mrSges(5,2) * t949 + mrSges(5,3) * t911;
t802 = m(5) * t834 + mrSges(5,1) * t932 - mrSges(5,3) * t875 - t892 * t912 + t899 * t949 + t806;
t900 = mrSges(5,1) * t949 - mrSges(5,3) * t912;
t976 = -t808 * t955 + t960 * t811;
t803 = m(5) * t835 - mrSges(5,2) * t932 + mrSges(5,3) * t874 + t892 * t911 - t900 * t949 + t976;
t797 = t961 * t802 + t956 * t803;
t903 = Ifges(4,4) * t938 + Ifges(4,2) * t937 + Ifges(4,6) * t950;
t904 = Ifges(4,1) * t938 + Ifges(4,4) * t937 + Ifges(4,5) * t950;
t992 = mrSges(4,1) * t876 - mrSges(4,2) * t877 + Ifges(4,5) * t909 + Ifges(4,6) * t908 + Ifges(4,3) * t936 + pkin(3) * t797 + t938 * t903 - t937 * t904 + t998;
t916 = -t963 * g(3) - t958 * t931;
t897 = -qJDD(2) * pkin(2) - pkin(8) * t965 + t940 * t986 - t916;
t869 = -pkin(3) * t908 - pkin(9) * t935 + t938 * t918 + t897;
t837 = -pkin(4) * t874 - pkin(10) * t910 + t912 * t901 + t869;
t827 = -pkin(5) * t847 - qJ(6) * t889 + t881 * t891 + qJDD(6) + t837;
t820 = m(7) * t827 - t847 * mrSges(7,1) + t848 * mrSges(7,2) - t890 * t879 + t891 * t882;
t988 = -t890 * t995 - t891 * t996 - t944 * t994;
t798 = -mrSges(6,1) * t837 + mrSges(6,3) * t825 - mrSges(7,1) * t827 + mrSges(7,3) * t822 - pkin(5) * t820 + qJ(6) * t982 + (-qJ(6) * t882 + t987) * t944 + (-mrSges(7,2) * qJ(6) + t995) * t926 + t988 * t891 + t990 * t848 + t1000 * t847;
t804 = mrSges(6,2) * t837 + mrSges(7,2) * t827 - mrSges(6,3) * t824 - mrSges(7,3) * t819 - qJ(6) * t816 + t990 * t847 + t848 * t997 - t988 * t890 + t996 * t926 - t993 * t944;
t884 = Ifges(5,5) * t912 + Ifges(5,6) * t911 + Ifges(5,3) * t949;
t973 = m(6) * t837 - t847 * mrSges(6,1) + t848 * mrSges(6,2) - t890 * t880 + t891 * t883 + t820;
t792 = -mrSges(5,1) * t869 + mrSges(5,3) * t835 + Ifges(5,4) * t875 + Ifges(5,2) * t874 + Ifges(5,6) * t932 - pkin(4) * t973 + pkin(10) * t976 + t960 * t798 + t955 * t804 - t912 * t884 + t949 * t886;
t793 = mrSges(5,2) * t869 - mrSges(5,3) * t834 + Ifges(5,1) * t875 + Ifges(5,4) * t874 + Ifges(5,5) * t932 - pkin(10) * t806 - t798 * t955 + t804 * t960 + t884 * t911 - t885 * t949;
t902 = Ifges(4,5) * t938 + Ifges(4,6) * t937 + Ifges(4,3) * t950;
t970 = m(5) * t869 - t874 * mrSges(5,1) + t875 * mrSges(5,2) - t911 * t899 + t912 * t900 + t973;
t977 = -t802 * t956 + t961 * t803;
t775 = -mrSges(4,1) * t897 + mrSges(4,3) * t877 + Ifges(4,4) * t909 + Ifges(4,2) * t908 + Ifges(4,6) * t936 - pkin(3) * t970 + pkin(9) * t977 + t961 * t792 + t956 * t793 - t938 * t902 + t950 * t904;
t776 = mrSges(4,2) * t897 - mrSges(4,3) * t876 + Ifges(4,1) * t909 + Ifges(4,4) * t908 + Ifges(4,5) * t936 - pkin(9) * t797 - t792 * t956 + t793 * t961 + t902 * t937 - t903 * t950;
t913 = -mrSges(4,1) * t937 + mrSges(4,2) * t938;
t914 = -mrSges(4,2) * t950 + mrSges(4,3) * t937;
t795 = m(4) * t876 + mrSges(4,1) * t936 - mrSges(4,3) * t909 - t913 * t938 + t914 * t950 + t797;
t915 = mrSges(4,1) * t950 - mrSges(4,3) * t938;
t796 = m(4) * t877 - mrSges(4,2) * t936 + mrSges(4,3) * t908 + t913 * t937 - t915 * t950 + t977;
t791 = -t795 * t957 + t962 * t796;
t814 = -m(4) * t897 + t908 * mrSges(4,1) - t909 * mrSges(4,2) + t937 * t914 - t938 * t915 - t970;
t928 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t958 + Ifges(3,2) * t963) * qJD(1);
t929 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t958 + Ifges(3,4) * t963) * qJD(1);
t991 = mrSges(3,1) * t916 - mrSges(3,2) * t917 + Ifges(3,5) * t941 + Ifges(3,6) * t942 + Ifges(3,3) * qJDD(2) + pkin(2) * t814 + pkin(8) * t791 + t962 * t775 + t957 * t776 + (t928 * t958 - t929 * t963) * qJD(1);
t939 = (-mrSges(3,1) * t963 + mrSges(3,2) * t958) * qJD(1);
t945 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t986;
t789 = m(3) * t917 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t942 - qJD(2) * t945 + t939 * t985 + t791;
t946 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t985;
t813 = m(3) * t916 + qJDD(2) * mrSges(3,1) - t941 * mrSges(3,3) + qJD(2) * t946 - t939 * t986 + t814;
t978 = t963 * t789 - t813 * t958;
t781 = m(2) * t948 - mrSges(2,1) * t966 - qJDD(1) * mrSges(2,2) + t978;
t790 = t795 * t962 + t796 * t957;
t972 = -m(3) * t930 + t942 * mrSges(3,1) - mrSges(3,2) * t941 - t945 * t986 + t946 * t985 - t790;
t785 = m(2) * t947 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t966 + t972;
t989 = t959 * t781 + t964 * t785;
t783 = t958 * t789 + t963 * t813;
t979 = t964 * t781 - t785 * t959;
t927 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t958 + Ifges(3,6) * t963) * qJD(1);
t774 = mrSges(3,2) * t930 - mrSges(3,3) * t916 + Ifges(3,1) * t941 + Ifges(3,4) * t942 + Ifges(3,5) * qJDD(2) - pkin(8) * t790 - qJD(2) * t928 - t775 * t957 + t776 * t962 + t927 * t985;
t778 = -mrSges(3,1) * t930 + mrSges(3,3) * t917 + Ifges(3,4) * t941 + Ifges(3,2) * t942 + Ifges(3,6) * qJDD(2) - pkin(2) * t790 + qJD(2) * t929 - t927 * t986 - t992;
t974 = mrSges(2,1) * t947 - mrSges(2,2) * t948 + Ifges(2,3) * qJDD(1) + pkin(1) * t972 + pkin(7) * t978 + t958 * t774 + t963 * t778;
t772 = mrSges(2,1) * g(3) + mrSges(2,3) * t948 + t966 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t783 - t991;
t771 = -mrSges(2,2) * g(3) - mrSges(2,3) * t947 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t966 - pkin(7) * t783 + t774 * t963 - t778 * t958;
t1 = [-m(1) * g(1) + t979; -m(1) * g(2) + t989; (-m(1) - m(2)) * g(3) + t783; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t989 + t964 * t771 - t959 * t772; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t979 + t959 * t771 + t964 * t772; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t974; t974; t991; t992; t998; t999; t820;];
tauJB  = t1;
