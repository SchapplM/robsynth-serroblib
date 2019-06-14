% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:40:45
% EndTime: 2019-05-07 10:41:11
% DurationCPUTime: 13.60s
% Computational Cost: add. (202157->367), mult. (415579->444), div. (0->0), fcn. (290629->10), ass. (0->149)
t1022 = Ifges(4,4) + Ifges(5,6);
t1034 = Ifges(4,2) + Ifges(5,3);
t1029 = Ifges(4,6) - Ifges(5,5);
t1030 = Ifges(4,5) - Ifges(5,4);
t1031 = Ifges(4,1) + Ifges(5,2);
t989 = cos(qJ(2));
t1014 = qJD(1) * t989;
t985 = sin(qJ(2));
t1015 = qJD(1) * t985;
t1024 = cos(qJ(3));
t984 = sin(qJ(3));
t953 = -t1024 * t1014 + t984 * t1015;
t954 = (t1024 * t985 + t984 * t989) * qJD(1);
t979 = qJD(2) + qJD(3);
t1016 = -t1022 * t953 + t1030 * t979 + t1031 * t954;
t1027 = -t1022 * t954 - t1029 * t979 + t1034 * t953;
t1028 = Ifges(4,3) + Ifges(5,1);
t1032 = -2 * qJD(4);
t1012 = qJD(1) * qJD(2);
t986 = sin(qJ(1));
t990 = cos(qJ(1));
t969 = -t990 * g(1) - t986 * g(2);
t991 = qJD(1) ^ 2;
t956 = -t991 * pkin(1) + qJDD(1) * pkin(7) + t969;
t1020 = t985 * t956;
t1023 = pkin(2) * t991;
t962 = t985 * qJDD(1) + t989 * t1012;
t897 = qJDD(2) * pkin(2) - t962 * pkin(8) - t1020 + (pkin(8) * t1012 + t985 * t1023 - g(3)) * t989;
t938 = -t985 * g(3) + t989 * t956;
t963 = t989 * qJDD(1) - t985 * t1012;
t967 = qJD(2) * pkin(2) - pkin(8) * t1015;
t981 = t989 ^ 2;
t898 = t963 * pkin(8) - qJD(2) * t967 - t981 * t1023 + t938;
t876 = t1024 * t898 + t984 * t897;
t929 = t953 * pkin(3) - t954 * qJ(4);
t977 = t979 ^ 2;
t978 = qJDD(2) + qJDD(3);
t866 = t977 * pkin(3) - t978 * qJ(4) + t979 * t1032 + t953 * t929 - t876;
t914 = t954 * qJD(3) - t1024 * t963 + t984 * t962;
t943 = t954 * pkin(4) - t979 * pkin(9);
t949 = t953 ^ 2;
t857 = -t914 * pkin(4) - t949 * pkin(9) + t979 * t943 - t866;
t983 = sin(qJ(5));
t988 = cos(qJ(5));
t936 = t983 * t953 + t988 * t979;
t881 = -t936 * qJD(5) + t988 * t914 - t983 * t978;
t948 = qJD(5) + t954;
t919 = t948 * pkin(5) - t936 * pkin(10);
t935 = t988 * t953 - t983 * t979;
t934 = t935 ^ 2;
t850 = -t881 * pkin(5) - t934 * pkin(10) + t936 * t919 + t857;
t882 = t935 * qJD(5) + t983 * t914 + t988 * t978;
t982 = sin(qJ(6));
t987 = cos(qJ(6));
t890 = t982 * t935 + t987 * t936;
t861 = -t890 * qJD(6) + t987 * t881 - t982 * t882;
t889 = t987 * t935 - t982 * t936;
t862 = t889 * qJD(6) + t982 * t881 + t987 * t882;
t946 = qJD(6) + t948;
t883 = -t946 * mrSges(7,2) + t889 * mrSges(7,3);
t884 = t946 * mrSges(7,1) - t890 * mrSges(7,3);
t1003 = m(7) * t850 - t861 * mrSges(7,1) + t862 * mrSges(7,2) - t889 * t883 + t890 * t884;
t1021 = t953 * t979;
t915 = -t953 * qJD(3) + t1024 * t962 + t984 * t963;
t968 = t986 * g(1) - t990 * g(2);
t1005 = -qJDD(1) * pkin(1) - t968;
t916 = -t963 * pkin(2) + t967 * t1015 + (-pkin(8) * t981 - pkin(7)) * t991 + t1005;
t995 = (-t915 + t1021) * qJ(4) + t916 + (t979 * pkin(3) + t1032) * t954;
t852 = -t949 * pkin(4) - t954 * t943 + (pkin(3) + pkin(9)) * t914 + t995;
t875 = t1024 * t897 - t984 * t898;
t868 = -t978 * pkin(3) - t977 * qJ(4) + t954 * t929 + qJDD(4) - t875;
t855 = (t953 * t954 - t978) * pkin(9) + (t915 + t1021) * pkin(4) + t868;
t847 = -t983 * t852 + t988 * t855;
t913 = qJDD(5) + t915;
t844 = (t935 * t948 - t882) * pkin(10) + (t935 * t936 + t913) * pkin(5) + t847;
t848 = t988 * t852 + t983 * t855;
t845 = -t934 * pkin(5) + t881 * pkin(10) - t948 * t919 + t848;
t842 = t987 * t844 - t982 * t845;
t873 = -t889 * mrSges(7,1) + t890 * mrSges(7,2);
t901 = qJDD(6) + t913;
t837 = m(7) * t842 + t901 * mrSges(7,1) - t862 * mrSges(7,3) - t890 * t873 + t946 * t883;
t843 = t982 * t844 + t987 * t845;
t838 = m(7) * t843 - t901 * mrSges(7,2) + t861 * mrSges(7,3) + t889 * t873 - t946 * t884;
t1007 = -t982 * t837 + t987 * t838;
t869 = Ifges(7,5) * t890 + Ifges(7,6) * t889 + Ifges(7,3) * t946;
t871 = Ifges(7,1) * t890 + Ifges(7,4) * t889 + Ifges(7,5) * t946;
t829 = -mrSges(7,1) * t850 + mrSges(7,3) * t843 + Ifges(7,4) * t862 + Ifges(7,2) * t861 + Ifges(7,6) * t901 - t890 * t869 + t946 * t871;
t870 = Ifges(7,4) * t890 + Ifges(7,2) * t889 + Ifges(7,6) * t946;
t830 = mrSges(7,2) * t850 - mrSges(7,3) * t842 + Ifges(7,1) * t862 + Ifges(7,4) * t861 + Ifges(7,5) * t901 + t889 * t869 - t946 * t870;
t885 = Ifges(6,5) * t936 + Ifges(6,6) * t935 + Ifges(6,3) * t948;
t887 = Ifges(6,1) * t936 + Ifges(6,4) * t935 + Ifges(6,5) * t948;
t811 = -mrSges(6,1) * t857 + mrSges(6,3) * t848 + Ifges(6,4) * t882 + Ifges(6,2) * t881 + Ifges(6,6) * t913 - pkin(5) * t1003 + pkin(10) * t1007 + t987 * t829 + t982 * t830 - t936 * t885 + t948 * t887;
t828 = t987 * t837 + t982 * t838;
t886 = Ifges(6,4) * t936 + Ifges(6,2) * t935 + Ifges(6,6) * t948;
t813 = mrSges(6,2) * t857 - mrSges(6,3) * t847 + Ifges(6,1) * t882 + Ifges(6,4) * t881 + Ifges(6,5) * t913 - pkin(10) * t828 - t982 * t829 + t987 * t830 + t935 * t885 - t948 * t886;
t894 = -t935 * mrSges(6,1) + t936 * mrSges(6,2);
t917 = -t948 * mrSges(6,2) + t935 * mrSges(6,3);
t825 = m(6) * t847 + t913 * mrSges(6,1) - t882 * mrSges(6,3) - t936 * t894 + t948 * t917 + t828;
t918 = t948 * mrSges(6,1) - t936 * mrSges(6,3);
t826 = m(6) * t848 - t913 * mrSges(6,2) + t881 * mrSges(6,3) + t935 * t894 - t948 * t918 + t1007;
t822 = t988 * t825 + t983 * t826;
t931 = -t953 * mrSges(5,2) - t954 * mrSges(5,3);
t1001 = -m(5) * t868 - t915 * mrSges(5,1) - t954 * t931 - t822;
t941 = t953 * mrSges(5,1) - t979 * mrSges(5,3);
t821 = t978 * mrSges(5,2) + t979 * t941 - t1001;
t942 = t954 * mrSges(5,1) + t979 * mrSges(5,2);
t999 = -m(6) * t857 + t881 * mrSges(6,1) - t882 * mrSges(6,2) + t935 * t917 - t936 * t918 - t1003;
t997 = -m(5) * t866 + t978 * mrSges(5,3) + t979 * t942 - t999;
t1033 = t1016 * t953 - mrSges(4,2) * t876 - mrSges(5,3) * t866 - pkin(3) * t821 - pkin(9) * t822 - t983 * t811 + t988 * t813 + qJ(4) * (-t953 * t931 + t997) + mrSges(5,2) * t868 + mrSges(4,1) * t875 + t1028 * t978 - t1027 * t954 + t1030 * t915 + (-qJ(4) * mrSges(5,1) - t1029) * t914;
t930 = t953 * mrSges(4,1) + t954 * mrSges(4,2);
t939 = -t979 * mrSges(4,2) - t953 * mrSges(4,3);
t818 = m(4) * t875 - t915 * mrSges(4,3) - t954 * t930 + (t939 - t941) * t979 + (mrSges(4,1) - mrSges(5,2)) * t978 + t1001;
t940 = t979 * mrSges(4,1) - t954 * mrSges(4,3);
t833 = (-t930 - t931) * t953 + (-mrSges(4,3) - mrSges(5,1)) * t914 - t979 * t940 - t978 * mrSges(4,2) + m(4) * t876 + t997;
t810 = t1024 * t818 + t984 * t833;
t937 = -t989 * g(3) - t1020;
t951 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t985 + Ifges(3,2) * t989) * qJD(1);
t952 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t985 + Ifges(3,4) * t989) * qJD(1);
t1026 = mrSges(3,1) * t937 - mrSges(3,2) * t938 + Ifges(3,5) * t962 + Ifges(3,6) * t963 + Ifges(3,3) * qJDD(2) + pkin(2) * t810 + (t985 * t951 - t989 * t952) * qJD(1) + t1033;
t961 = (-mrSges(3,1) * t989 + mrSges(3,2) * t985) * qJD(1);
t966 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1014;
t808 = m(3) * t937 + qJDD(2) * mrSges(3,1) - t962 * mrSges(3,3) + qJD(2) * t966 - t961 * t1015 + t810;
t1008 = t1024 * t833 - t984 * t818;
t965 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1015;
t809 = m(3) * t938 - qJDD(2) * mrSges(3,2) + t963 * mrSges(3,3) - qJD(2) * t965 + t961 * t1014 + t1008;
t1009 = -t985 * t808 + t989 * t809;
t802 = m(2) * t969 - t991 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1009;
t955 = -t991 * pkin(7) + t1005;
t1018 = -t983 * t825 + t988 * t826;
t864 = t914 * pkin(3) + t995;
t819 = m(5) * t864 - t914 * mrSges(5,2) - t915 * mrSges(5,3) - t953 * t941 - t954 * t942 + t1018;
t998 = m(4) * t916 + t914 * mrSges(4,1) + t915 * mrSges(4,2) + t953 * t939 + t954 * t940 + t819;
t993 = -m(3) * t955 + t963 * mrSges(3,1) - t962 * mrSges(3,2) + t966 * t1014 - t965 * t1015 - t998;
t815 = m(2) * t968 + qJDD(1) * mrSges(2,1) - t991 * mrSges(2,2) + t993;
t1019 = t986 * t802 + t990 * t815;
t804 = t989 * t808 + t985 * t809;
t1017 = -t1028 * t979 + t1029 * t953 - t1030 * t954;
t1010 = t990 * t802 - t986 * t815;
t798 = -mrSges(4,1) * t916 - mrSges(5,1) * t866 + mrSges(5,2) * t864 + mrSges(4,3) * t876 - pkin(3) * t819 - pkin(4) * t999 - pkin(9) * t1018 + t1016 * t979 + t1017 * t954 + t1022 * t915 + t1029 * t978 - t1034 * t914 - t988 * t811 - t983 * t813;
t1000 = mrSges(7,1) * t842 - mrSges(7,2) * t843 + Ifges(7,5) * t862 + Ifges(7,6) * t861 + Ifges(7,3) * t901 + t890 * t870 - t889 * t871;
t996 = mrSges(6,1) * t847 - mrSges(6,2) * t848 + Ifges(6,5) * t882 + Ifges(6,6) * t881 + Ifges(6,3) * t913 + pkin(5) * t828 + t936 * t886 - t935 * t887 + t1000;
t799 = mrSges(5,1) * t868 + mrSges(4,2) * t916 - mrSges(4,3) * t875 - mrSges(5,3) * t864 + pkin(4) * t822 - qJ(4) * t819 + t1017 * t953 - t1022 * t914 + t1027 * t979 + t1030 * t978 + t1031 * t915 + t996;
t950 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t985 + Ifges(3,6) * t989) * qJD(1);
t794 = -mrSges(3,1) * t955 + mrSges(3,3) * t938 + Ifges(3,4) * t962 + Ifges(3,2) * t963 + Ifges(3,6) * qJDD(2) - pkin(2) * t998 + pkin(8) * t1008 + qJD(2) * t952 - t950 * t1015 + t1024 * t798 + t984 * t799;
t796 = mrSges(3,2) * t955 - mrSges(3,3) * t937 + Ifges(3,1) * t962 + Ifges(3,4) * t963 + Ifges(3,5) * qJDD(2) - pkin(8) * t810 - qJD(2) * t951 + t950 * t1014 + t1024 * t799 - t984 * t798;
t1002 = mrSges(2,1) * t968 - mrSges(2,2) * t969 + Ifges(2,3) * qJDD(1) + pkin(1) * t993 + pkin(7) * t1009 + t989 * t794 + t985 * t796;
t797 = mrSges(2,1) * g(3) + mrSges(2,3) * t969 + t991 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t804 - t1026;
t792 = -mrSges(2,2) * g(3) - mrSges(2,3) * t968 + Ifges(2,5) * qJDD(1) - t991 * Ifges(2,6) - pkin(7) * t804 - t985 * t794 + t989 * t796;
t1 = [-m(1) * g(1) + t1010; -m(1) * g(2) + t1019; (-m(1) - m(2)) * g(3) + t804; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1019 + t990 * t792 - t986 * t797; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1010 + t986 * t792 + t990 * t797; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1002; t1002; t1026; t1033; t821; t996; t1000;];
tauJB  = t1;
