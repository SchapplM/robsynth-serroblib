% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:23:34
% EndTime: 2019-05-06 12:23:45
% DurationCPUTime: 7.44s
% Computational Cost: add. (81278->344), mult. (183485->407), div. (0->0), fcn. (124047->8), ass. (0->135)
t995 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t971 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t970 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t994 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t969 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t993 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t945 = sin(qJ(2));
t947 = cos(qJ(2));
t972 = qJD(1) * qJD(2);
t929 = t945 * qJDD(1) + t947 * t972;
t930 = t947 * qJDD(1) - t945 * t972;
t943 = sin(pkin(9));
t981 = cos(pkin(9));
t902 = t929 * t981 + t943 * t930;
t918 = (t943 * t947 + t945 * t981) * qJD(1);
t944 = sin(qJ(4));
t986 = cos(qJ(4));
t904 = -qJD(2) * t986 + t944 * t918;
t868 = -t904 * qJD(4) + t944 * qJDD(2) + t902 * t986;
t974 = qJD(1) * t947;
t975 = qJD(1) * t945;
t917 = -t943 * t975 + t981 * t974;
t916 = qJD(4) - t917;
t980 = t904 * t916;
t992 = (-t868 + t980) * qJ(5);
t946 = sin(qJ(1));
t948 = cos(qJ(1));
t936 = -t948 * g(1) - t946 * g(2);
t950 = qJD(1) ^ 2;
t924 = -t950 * pkin(1) + qJDD(1) * pkin(7) + t936;
t979 = t945 * t924;
t984 = pkin(2) * t950;
t875 = qJDD(2) * pkin(2) - t929 * qJ(3) - t979 + (qJ(3) * t972 + t945 * t984 - g(3)) * t947;
t907 = -t945 * g(3) + t947 * t924;
t932 = qJD(2) * pkin(2) - qJ(3) * t975;
t942 = t947 ^ 2;
t876 = t930 * qJ(3) - qJD(2) * t932 - t942 * t984 + t907;
t843 = 0.2e1 * qJD(3) * t917 + t943 * t875 + t981 * t876;
t895 = -t917 * pkin(3) - t918 * pkin(8);
t949 = qJD(2) ^ 2;
t839 = -pkin(3) * t949 + qJDD(2) * pkin(8) + t895 * t917 + t843;
t935 = t946 * g(1) - t948 * g(2);
t957 = -qJDD(1) * pkin(1) - t935;
t882 = -t930 * pkin(2) + qJDD(3) + t932 * t975 + (-qJ(3) * t942 - pkin(7)) * t950 + t957;
t901 = -t943 * t929 + t981 * t930;
t841 = (-qJD(2) * t917 - t902) * pkin(8) + (qJD(2) * t918 - t901) * pkin(3) + t882;
t834 = -t944 * t839 + t841 * t986;
t905 = t944 * qJD(2) + t918 * t986;
t877 = t904 * pkin(4) - t905 * qJ(5);
t900 = qJDD(4) - t901;
t915 = t916 ^ 2;
t832 = -t900 * pkin(4) - t915 * qJ(5) + t905 * t877 + qJDD(5) - t834;
t883 = -t904 * mrSges(6,2) + t916 * mrSges(6,3);
t991 = -m(6) * t832 + t900 * mrSges(6,1) + t916 * t883;
t867 = t905 * qJD(4) - qJDD(2) * t986 + t944 * t902;
t886 = -t916 * pkin(5) - t905 * qJ(6);
t903 = t904 ^ 2;
t973 = qJD(3) * t918;
t912 = -0.2e1 * t973;
t977 = t981 * t875 - t943 * t876;
t955 = qJDD(2) * pkin(3) + t949 * pkin(8) - t918 * t895 + t977;
t987 = 2 * qJD(5);
t829 = -t903 * qJ(6) + qJDD(6) + t912 + (-pkin(4) - pkin(5)) * t867 - t992 + (-pkin(4) * t916 + t886 + t987) * t905 + t955;
t884 = t916 * mrSges(7,2) + t904 * mrSges(7,3);
t887 = -t916 * mrSges(7,1) - t905 * mrSges(7,3);
t824 = m(7) * t829 - t867 * mrSges(7,1) + t868 * mrSges(7,2) - t904 * t884 + t905 * t887;
t838 = 0.2e1 * t973 - t955;
t988 = -0.2e1 * t905;
t833 = qJD(5) * t988 + t992 + (t905 * t916 + t867) * pkin(4) + t838;
t889 = -t916 * mrSges(6,1) + t905 * mrSges(6,2);
t822 = m(6) * t833 + t867 * mrSges(6,1) - t868 * mrSges(6,3) + t904 * t883 - t905 * t889 - t824;
t835 = t986 * t839 + t944 * t841;
t831 = -pkin(4) * t915 + t900 * qJ(5) - t904 * t877 + t916 * t987 + t835;
t827 = -pkin(5) * t903 + qJ(6) * t867 + 0.2e1 * qJD(6) * t904 + t886 * t916 + t831;
t964 = -t971 * t904 + t995 * t905 + t970 * t916;
t966 = t969 * t904 - t970 * t905 - t993 * t916;
t879 = -t904 * mrSges(7,1) + t905 * mrSges(7,2);
t967 = m(7) * t827 + t867 * mrSges(7,3) + t904 * t879;
t795 = -mrSges(5,1) * t838 + mrSges(5,3) * t835 - mrSges(6,1) * t833 + mrSges(6,2) * t831 + mrSges(7,1) * t829 - mrSges(7,3) * t827 + pkin(5) * t824 - qJ(6) * t967 - pkin(4) * t822 + (-qJ(6) * t887 + t964) * t916 + t966 * t905 + (-mrSges(7,2) * qJ(6) + t969) * t900 + t971 * t868 + t994 * t867;
t885 = -t916 * mrSges(5,2) - t904 * mrSges(5,3);
t825 = qJD(6) * t988 + (-t868 - t980) * qJ(6) + (t904 * t905 - t900) * pkin(5) + t832;
t959 = -m(7) * t825 + t868 * mrSges(7,3) + t905 * t879;
t878 = t904 * mrSges(6,1) - t905 * mrSges(6,3);
t976 = -t904 * mrSges(5,1) - t905 * mrSges(5,2) - t878;
t983 = -mrSges(5,3) - mrSges(6,2);
t817 = m(5) * t834 + (t884 + t885) * t916 + t976 * t905 + (mrSges(5,1) + mrSges(7,1)) * t900 + t983 * t868 + t959 + t991;
t888 = t916 * mrSges(5,1) - t905 * mrSges(5,3);
t956 = m(6) * t831 + t900 * mrSges(6,3) + t916 * t889 + t967;
t818 = m(5) * t835 + (t887 - t888) * t916 + t976 * t904 + (-mrSges(5,2) + mrSges(7,2)) * t900 + t983 * t867 + t956;
t811 = -t817 * t944 + t986 * t818;
t894 = -t917 * mrSges(4,1) + t918 * mrSges(4,2);
t909 = qJD(2) * mrSges(4,1) - t918 * mrSges(4,3);
t808 = m(4) * t843 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t901 - qJD(2) * t909 + t894 * t917 + t811;
t819 = -m(5) * t838 - t867 * mrSges(5,1) - t868 * mrSges(5,2) - t904 * t885 - t905 * t888 - t822;
t842 = t912 + t977;
t908 = -qJD(2) * mrSges(4,2) + t917 * mrSges(4,3);
t813 = m(4) * t842 + qJDD(2) * mrSges(4,1) - t902 * mrSges(4,3) + qJD(2) * t908 - t918 * t894 + t819;
t801 = t943 * t808 + t981 * t813;
t823 = -t900 * mrSges(7,1) - t916 * t884 - t959;
t965 = t994 * t904 + t971 * t905 + t969 * t916;
t802 = mrSges(5,2) * t838 + mrSges(6,2) * t832 + mrSges(7,2) * t829 - mrSges(5,3) * t834 - mrSges(6,3) * t833 - mrSges(7,3) * t825 - qJ(5) * t822 - qJ(6) * t823 - t971 * t867 + t995 * t868 + t970 * t900 + t966 * t904 - t965 * t916;
t891 = Ifges(4,4) * t918 + Ifges(4,2) * t917 + Ifges(4,6) * qJD(2);
t892 = Ifges(4,1) * t918 + Ifges(4,4) * t917 + Ifges(4,5) * qJD(2);
t906 = -t947 * g(3) - t979;
t920 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t945 + Ifges(3,2) * t947) * qJD(1);
t921 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t945 + Ifges(3,4) * t947) * qJD(1);
t990 = (t920 * t945 - t921 * t947) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t906 + mrSges(4,1) * t842 - mrSges(3,2) * t907 - mrSges(4,2) * t843 + Ifges(3,5) * t929 + Ifges(4,5) * t902 + Ifges(3,6) * t930 + Ifges(4,6) * t901 + pkin(2) * t801 + pkin(3) * t819 + pkin(8) * t811 + t795 * t986 + t944 * t802 + t918 * t891 - t917 * t892;
t820 = t868 * mrSges(6,2) + t905 * t878 + t823 - t991;
t989 = -t867 * t969 + t868 * t970 + t993 * t900 + t904 * t964 + t905 * t965 + mrSges(5,1) * t834 - mrSges(6,1) * t832 - mrSges(7,1) * t825 - mrSges(5,2) * t835 + mrSges(7,2) * t827 + mrSges(6,3) * t831 - pkin(4) * t820 - pkin(5) * t823 + qJ(5) * (-t867 * mrSges(6,2) + t900 * mrSges(7,2) - t904 * t878 + t916 * t887 + t956);
t928 = (-mrSges(3,1) * t947 + mrSges(3,2) * t945) * qJD(1);
t934 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t974;
t799 = m(3) * t906 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t929 + qJD(2) * t934 - t928 * t975 + t801;
t933 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t975;
t961 = t981 * t808 - t813 * t943;
t800 = m(3) * t907 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t930 - qJD(2) * t933 + t928 * t974 + t961;
t962 = -t799 * t945 + t947 * t800;
t791 = m(2) * t936 - mrSges(2,1) * t950 - qJDD(1) * mrSges(2,2) + t962;
t810 = t986 * t817 + t944 * t818;
t809 = m(4) * t882 - t901 * mrSges(4,1) + mrSges(4,2) * t902 - t917 * t908 + t909 * t918 + t810;
t923 = -t950 * pkin(7) + t957;
t952 = -m(3) * t923 + t930 * mrSges(3,1) - mrSges(3,2) * t929 - t933 * t975 + t934 * t974 - t809;
t804 = m(2) * t935 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t950 + t952;
t978 = t946 * t791 + t948 * t804;
t793 = t947 * t799 + t945 * t800;
t963 = t948 * t791 - t804 * t946;
t890 = Ifges(4,5) * t918 + Ifges(4,6) * t917 + Ifges(4,3) * qJD(2);
t788 = mrSges(4,2) * t882 - mrSges(4,3) * t842 + Ifges(4,1) * t902 + Ifges(4,4) * t901 + Ifges(4,5) * qJDD(2) - pkin(8) * t810 - qJD(2) * t891 - t944 * t795 + t802 * t986 + t917 * t890;
t794 = -mrSges(4,1) * t882 + mrSges(4,3) * t843 + Ifges(4,4) * t902 + Ifges(4,2) * t901 + Ifges(4,6) * qJDD(2) - pkin(3) * t810 + qJD(2) * t892 - t918 * t890 - t989;
t919 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t945 + Ifges(3,6) * t947) * qJD(1);
t784 = -mrSges(3,1) * t923 + mrSges(3,3) * t907 + Ifges(3,4) * t929 + Ifges(3,2) * t930 + Ifges(3,6) * qJDD(2) - pkin(2) * t809 + qJ(3) * t961 + qJD(2) * t921 + t943 * t788 + t794 * t981 - t919 * t975;
t787 = mrSges(3,2) * t923 - mrSges(3,3) * t906 + Ifges(3,1) * t929 + Ifges(3,4) * t930 + Ifges(3,5) * qJDD(2) - qJ(3) * t801 - qJD(2) * t920 + t788 * t981 - t943 * t794 + t919 * t974;
t954 = mrSges(2,1) * t935 - mrSges(2,2) * t936 + Ifges(2,3) * qJDD(1) + pkin(1) * t952 + pkin(7) * t962 + t947 * t784 + t945 * t787;
t785 = mrSges(2,1) * g(3) + mrSges(2,3) * t936 + t950 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t793 - t990;
t782 = -mrSges(2,2) * g(3) - mrSges(2,3) * t935 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t950 - pkin(7) * t793 - t784 * t945 + t787 * t947;
t1 = [-m(1) * g(1) + t963; -m(1) * g(2) + t978; (-m(1) - m(2)) * g(3) + t793; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t978 + t948 * t782 - t946 * t785; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t963 + t946 * t782 + t948 * t785; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t954; t954; t990; t809; t989; t820; t824;];
tauJB  = t1;
