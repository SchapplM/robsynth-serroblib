% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP4
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
% Datum: 2019-05-08 04:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:48:35
% EndTime: 2019-05-08 04:49:02
% DurationCPUTime: 17.38s
% Computational Cost: add. (280851->366), mult. (558631->446), div. (0->0), fcn. (404803->10), ass. (0->147)
t1004 = Ifges(6,1) + Ifges(7,1);
t993 = Ifges(6,4) - Ifges(7,5);
t1002 = Ifges(7,4) + Ifges(6,5);
t1003 = Ifges(6,2) + Ifges(7,3);
t1001 = Ifges(6,6) - Ifges(7,6);
t1000 = -Ifges(6,3) - Ifges(7,2);
t960 = sin(qJ(3));
t961 = sin(qJ(2));
t964 = cos(qJ(3));
t965 = cos(qJ(2));
t934 = (t960 * t965 + t961 * t964) * qJD(1);
t955 = qJD(2) + qJD(3);
t959 = sin(qJ(4));
t963 = cos(qJ(4));
t919 = -t934 * t959 + t955 * t963;
t920 = t934 * t963 + t955 * t959;
t958 = sin(qJ(5));
t996 = cos(qJ(5));
t886 = -t919 * t996 + t958 * t920;
t887 = t958 * t919 + t920 * t996;
t987 = qJD(1) * t965;
t988 = qJD(1) * t961;
t933 = -t960 * t988 + t964 * t987;
t929 = qJD(4) - t933;
t927 = qJD(5) + t929;
t999 = -t1001 * t927 + t1003 * t886 - t993 * t887;
t998 = t1002 * t927 + t1004 * t887 - t993 * t886;
t986 = qJD(1) * qJD(2);
t942 = qJDD(1) * t961 + t965 * t986;
t962 = sin(qJ(1));
t966 = cos(qJ(1));
t949 = -g(1) * t966 - g(2) * t962;
t967 = qJD(1) ^ 2;
t936 = -pkin(1) * t967 + qJDD(1) * pkin(7) + t949;
t992 = t961 * t936;
t995 = pkin(2) * t967;
t893 = qJDD(2) * pkin(2) - t942 * pkin(8) - t992 + (pkin(8) * t986 + t961 * t995 - g(3)) * t965;
t922 = -g(3) * t961 + t965 * t936;
t943 = qJDD(1) * t965 - t961 * t986;
t947 = qJD(2) * pkin(2) - pkin(8) * t988;
t957 = t965 ^ 2;
t894 = pkin(8) * t943 - qJD(2) * t947 - t957 * t995 + t922;
t870 = t960 * t893 + t964 * t894;
t905 = -t934 * qJD(3) - t960 * t942 + t943 * t964;
t915 = -mrSges(4,1) * t933 + mrSges(4,2) * t934;
t924 = mrSges(4,1) * t955 - mrSges(4,3) * t934;
t954 = qJDD(2) + qJDD(3);
t906 = qJD(3) * t933 + t942 * t964 + t943 * t960;
t948 = t962 * g(1) - t966 * g(2);
t977 = -qJDD(1) * pkin(1) - t948;
t907 = -t943 * pkin(2) + t947 * t988 + (-pkin(8) * t957 - pkin(7)) * t967 + t977;
t852 = (-t933 * t955 - t906) * pkin(9) + (t934 * t955 - t905) * pkin(3) + t907;
t916 = -pkin(3) * t933 - pkin(9) * t934;
t953 = t955 ^ 2;
t855 = -pkin(3) * t953 + pkin(9) * t954 + t916 * t933 + t870;
t836 = t963 * t852 - t959 * t855;
t876 = qJD(4) * t919 + t906 * t963 + t954 * t959;
t904 = qJDD(4) - t905;
t833 = (t919 * t929 - t876) * pkin(10) + (t919 * t920 + t904) * pkin(4) + t836;
t837 = t959 * t852 + t963 * t855;
t875 = -qJD(4) * t920 - t906 * t959 + t954 * t963;
t910 = pkin(4) * t929 - pkin(10) * t920;
t918 = t919 ^ 2;
t835 = -pkin(4) * t918 + pkin(10) * t875 - t910 * t929 + t837;
t829 = t958 * t833 + t835 * t996;
t846 = t887 * qJD(5) - t875 * t996 + t958 * t876;
t879 = mrSges(6,1) * t927 - mrSges(6,3) * t887;
t899 = qJDD(5) + t904;
t865 = pkin(5) * t886 - qJ(6) * t887;
t925 = t927 ^ 2;
t825 = -pkin(5) * t925 + qJ(6) * t899 + 0.2e1 * qJD(6) * t927 - t865 * t886 + t829;
t880 = -mrSges(7,1) * t927 + mrSges(7,2) * t887;
t985 = m(7) * t825 + t899 * mrSges(7,3) + t927 * t880;
t866 = mrSges(7,1) * t886 - mrSges(7,3) * t887;
t989 = -mrSges(6,1) * t886 - mrSges(6,2) * t887 - t866;
t994 = -mrSges(6,3) - mrSges(7,2);
t811 = m(6) * t829 - t899 * mrSges(6,2) + t846 * t994 - t927 * t879 + t886 * t989 + t985;
t828 = t833 * t996 - t958 * t835;
t847 = -t886 * qJD(5) + t958 * t875 + t876 * t996;
t878 = -mrSges(6,2) * t927 - mrSges(6,3) * t886;
t826 = -t899 * pkin(5) - t925 * qJ(6) + t887 * t865 + qJDD(6) - t828;
t877 = -mrSges(7,2) * t886 + mrSges(7,3) * t927;
t979 = -m(7) * t826 + t899 * mrSges(7,1) + t927 * t877;
t813 = m(6) * t828 + t899 * mrSges(6,1) + t847 * t994 + t927 * t878 + t887 * t989 + t979;
t808 = t958 * t811 + t813 * t996;
t891 = -mrSges(5,1) * t919 + mrSges(5,2) * t920;
t908 = -mrSges(5,2) * t929 + mrSges(5,3) * t919;
t804 = m(5) * t836 + mrSges(5,1) * t904 - mrSges(5,3) * t876 - t891 * t920 + t908 * t929 + t808;
t909 = mrSges(5,1) * t929 - mrSges(5,3) * t920;
t980 = t811 * t996 - t813 * t958;
t805 = m(5) * t837 - mrSges(5,2) * t904 + mrSges(5,3) * t875 + t891 * t919 - t909 * t929 + t980;
t981 = -t804 * t959 + t963 * t805;
t797 = m(4) * t870 - mrSges(4,2) * t954 + mrSges(4,3) * t905 + t915 * t933 - t924 * t955 + t981;
t869 = t964 * t893 - t960 * t894;
t923 = -mrSges(4,2) * t955 + mrSges(4,3) * t933;
t854 = -t954 * pkin(3) - t953 * pkin(9) + t934 * t916 - t869;
t838 = -t875 * pkin(4) - t918 * pkin(10) + t920 * t910 + t854;
t831 = -0.2e1 * qJD(6) * t887 + (t886 * t927 - t847) * qJ(6) + (t887 * t927 + t846) * pkin(5) + t838;
t822 = m(7) * t831 + t846 * mrSges(7,1) - t847 * mrSges(7,3) + t886 * t877 - t887 * t880;
t973 = m(6) * t838 + t846 * mrSges(6,1) + mrSges(6,2) * t847 + t886 * t878 + t879 * t887 + t822;
t970 = -m(5) * t854 + t875 * mrSges(5,1) - mrSges(5,2) * t876 + t919 * t908 - t909 * t920 - t973;
t815 = m(4) * t869 + mrSges(4,1) * t954 - mrSges(4,3) * t906 - t915 * t934 + t923 * t955 + t970;
t791 = t960 * t797 + t964 * t815;
t921 = -t965 * g(3) - t992;
t931 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t961 + Ifges(3,2) * t965) * qJD(1);
t932 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t961 + Ifges(3,4) * t965) * qJD(1);
t990 = t1000 * t927 + t1001 * t886 - t1002 * t887;
t806 = -mrSges(6,1) * t838 - mrSges(7,1) * t831 + mrSges(7,2) * t825 + mrSges(6,3) * t829 - pkin(5) * t822 + t1001 * t899 - t1003 * t846 + t993 * t847 + t990 * t887 + t998 * t927;
t807 = mrSges(6,2) * t838 + mrSges(7,2) * t826 - mrSges(6,3) * t828 - mrSges(7,3) * t831 - qJ(6) * t822 + t1002 * t899 + t1004 * t847 - t993 * t846 + t990 * t886 + t999 * t927;
t881 = Ifges(5,5) * t920 + Ifges(5,6) * t919 + Ifges(5,3) * t929;
t883 = Ifges(5,1) * t920 + Ifges(5,4) * t919 + Ifges(5,5) * t929;
t783 = -mrSges(5,1) * t854 + mrSges(5,3) * t837 + Ifges(5,4) * t876 + Ifges(5,2) * t875 + Ifges(5,6) * t904 - pkin(4) * t973 + pkin(10) * t980 + t806 * t996 + t958 * t807 - t920 * t881 + t929 * t883;
t882 = Ifges(5,4) * t920 + Ifges(5,2) * t919 + Ifges(5,6) * t929;
t785 = mrSges(5,2) * t854 - mrSges(5,3) * t836 + Ifges(5,1) * t876 + Ifges(5,4) * t875 + Ifges(5,5) * t904 - pkin(10) * t808 - t958 * t806 + t807 * t996 + t919 * t881 - t929 * t882;
t912 = Ifges(4,4) * t934 + Ifges(4,2) * t933 + Ifges(4,6) * t955;
t913 = Ifges(4,1) * t934 + Ifges(4,4) * t933 + Ifges(4,5) * t955;
t974 = -mrSges(4,1) * t869 + mrSges(4,2) * t870 - Ifges(4,5) * t906 - Ifges(4,6) * t905 - Ifges(4,3) * t954 - pkin(3) * t970 - pkin(9) * t981 - t963 * t783 - t959 * t785 - t934 * t912 + t933 * t913;
t997 = mrSges(3,1) * t921 - mrSges(3,2) * t922 + Ifges(3,5) * t942 + Ifges(3,6) * t943 + Ifges(3,3) * qJDD(2) + pkin(2) * t791 + (t931 * t961 - t932 * t965) * qJD(1) - t974;
t941 = (-mrSges(3,1) * t965 + mrSges(3,2) * t961) * qJD(1);
t946 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t987;
t789 = m(3) * t921 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t942 + qJD(2) * t946 - t941 * t988 + t791;
t945 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t988;
t982 = t964 * t797 - t815 * t960;
t790 = m(3) * t922 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t943 - qJD(2) * t945 + t941 * t987 + t982;
t983 = -t789 * t961 + t965 * t790;
t778 = m(2) * t949 - mrSges(2,1) * t967 - qJDD(1) * mrSges(2,2) + t983;
t935 = -t967 * pkin(7) + t977;
t799 = t963 * t804 + t959 * t805;
t975 = m(4) * t907 - t905 * mrSges(4,1) + mrSges(4,2) * t906 - t933 * t923 + t924 * t934 + t799;
t972 = -m(3) * t935 + t943 * mrSges(3,1) - mrSges(3,2) * t942 - t945 * t988 + t946 * t987 - t975;
t793 = m(2) * t948 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t967 + t972;
t991 = t962 * t778 + t966 * t793;
t780 = t965 * t789 + t961 * t790;
t984 = t966 * t778 - t793 * t962;
t911 = Ifges(4,5) * t934 + Ifges(4,6) * t933 + Ifges(4,3) * t955;
t775 = mrSges(4,2) * t907 - mrSges(4,3) * t869 + Ifges(4,1) * t906 + Ifges(4,4) * t905 + Ifges(4,5) * t954 - pkin(9) * t799 - t783 * t959 + t785 * t963 + t911 * t933 - t912 * t955;
t821 = t847 * mrSges(7,2) + t887 * t866 - t979;
t971 = -mrSges(6,1) * t828 + mrSges(7,1) * t826 + mrSges(6,2) * t829 - mrSges(7,3) * t825 + pkin(5) * t821 - qJ(6) * t985 + t1000 * t899 + t999 * t887 + (qJ(6) * t866 - t998) * t886 - t1002 * t847 + (mrSges(7,2) * qJ(6) + t1001) * t846;
t968 = mrSges(5,1) * t836 - mrSges(5,2) * t837 + Ifges(5,5) * t876 + Ifges(5,6) * t875 + Ifges(5,3) * t904 + pkin(4) * t808 + t920 * t882 - t919 * t883 - t971;
t781 = -mrSges(4,1) * t907 + mrSges(4,3) * t870 + Ifges(4,4) * t906 + Ifges(4,2) * t905 + Ifges(4,6) * t954 - pkin(3) * t799 - t934 * t911 + t955 * t913 - t968;
t930 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t961 + Ifges(3,6) * t965) * qJD(1);
t771 = -mrSges(3,1) * t935 + mrSges(3,3) * t922 + Ifges(3,4) * t942 + Ifges(3,2) * t943 + Ifges(3,6) * qJDD(2) - pkin(2) * t975 + pkin(8) * t982 + qJD(2) * t932 + t960 * t775 + t964 * t781 - t930 * t988;
t774 = mrSges(3,2) * t935 - mrSges(3,3) * t921 + Ifges(3,1) * t942 + Ifges(3,4) * t943 + Ifges(3,5) * qJDD(2) - pkin(8) * t791 - qJD(2) * t931 + t775 * t964 - t781 * t960 + t930 * t987;
t976 = mrSges(2,1) * t948 - mrSges(2,2) * t949 + Ifges(2,3) * qJDD(1) + pkin(1) * t972 + pkin(7) * t983 + t965 * t771 + t961 * t774;
t772 = mrSges(2,1) * g(3) + mrSges(2,3) * t949 + t967 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t780 - t997;
t769 = -mrSges(2,2) * g(3) - mrSges(2,3) * t948 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t967 - pkin(7) * t780 - t771 * t961 + t774 * t965;
t1 = [-m(1) * g(1) + t984; -m(1) * g(2) + t991; (-m(1) - m(2)) * g(3) + t780; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t991 + t966 * t769 - t962 * t772; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t984 + t962 * t769 + t966 * t772; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t976; t976; t997; -t974; t968; -t971; t821;];
tauJB  = t1;
