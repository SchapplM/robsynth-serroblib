% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-05-06 12:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:36:05
% EndTime: 2019-05-06 12:36:17
% DurationCPUTime: 7.72s
% Computational Cost: add. (81877->345), mult. (171727->405), div. (0->0), fcn. (103251->8), ass. (0->139)
t1019 = Ifges(6,1) + Ifges(7,1);
t1001 = Ifges(6,4) - Ifges(7,5);
t1015 = Ifges(7,4) + Ifges(6,5);
t1018 = Ifges(6,2) + Ifges(7,3);
t1012 = Ifges(6,6) - Ifges(7,6);
t1017 = -2 * qJD(3);
t1016 = Ifges(3,1) + Ifges(4,2);
t1002 = Ifges(3,4) + Ifges(4,6);
t1000 = Ifges(3,5) - Ifges(4,4);
t1014 = Ifges(3,2) + Ifges(4,3);
t1013 = -Ifges(7,2) - Ifges(6,3);
t999 = Ifges(3,6) - Ifges(4,5);
t1011 = Ifges(3,3) + Ifges(4,1);
t959 = sin(qJ(4));
t962 = cos(qJ(4));
t963 = cos(qJ(2));
t989 = qJD(1) * t963;
t929 = -t959 * qJD(2) - t962 * t989;
t930 = t962 * qJD(2) - t959 * t989;
t958 = sin(pkin(9));
t997 = cos(pkin(9));
t895 = -t997 * t929 + t958 * t930;
t896 = t958 * t929 + t997 * t930;
t960 = sin(qJ(2));
t988 = t960 * qJD(1);
t948 = qJD(4) + t988;
t1010 = -t1001 * t896 - t1012 * t948 + t1018 * t895;
t1009 = -t1001 * t895 + t1015 * t948 + t1019 * t896;
t961 = sin(qJ(1));
t964 = cos(qJ(1));
t945 = -t964 * g(1) - t961 * g(2);
t966 = qJD(1) ^ 2;
t916 = -t966 * pkin(1) + qJDD(1) * pkin(7) + t945;
t902 = -t960 * g(3) + t963 * t916;
t931 = (-pkin(2) * t963 - qJ(3) * t960) * qJD(1);
t965 = qJD(2) ^ 2;
t877 = t965 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t1017 - t931 * t989 - t902;
t987 = qJD(1) * qJD(2);
t984 = t960 * t987;
t935 = t963 * qJDD(1) - t984;
t943 = pkin(3) * t988 - qJD(2) * pkin(8);
t957 = t963 ^ 2;
t858 = -t957 * t966 * pkin(8) + t935 * pkin(3) + qJD(2) * t943 - t877;
t893 = -t930 * qJD(4) - t959 * qJDD(2) - t962 * t935;
t899 = t948 * pkin(4) - t930 * qJ(5);
t927 = t929 ^ 2;
t844 = -t893 * pkin(4) - t927 * qJ(5) + t930 * t899 + qJDD(5) + t858;
t894 = t929 * qJD(4) + t962 * qJDD(2) - t959 * t935;
t866 = -t997 * t893 + t958 * t894;
t867 = t958 * t893 + t997 * t894;
t835 = -0.2e1 * qJD(6) * t896 + (t896 * t948 + t866) * pkin(5) + t844 + (t895 * t948 - t867) * qJ(6);
t880 = -t895 * mrSges(7,2) + t948 * mrSges(7,3);
t882 = -t948 * mrSges(7,1) + t896 * mrSges(7,2);
t826 = m(7) * t835 + t866 * mrSges(7,1) - t867 * mrSges(7,3) + t895 * t880 - t896 * t882;
t1005 = -2 * qJD(5);
t983 = t963 * t987;
t934 = t960 * qJDD(1) + t983;
t944 = t961 * g(1) - t964 * g(2);
t978 = -qJDD(1) * pkin(1) - t944;
t972 = pkin(2) * t984 + t988 * t1017 + (-t934 - t983) * qJ(3) + t978;
t848 = -t943 * t988 + (-pkin(3) * t957 - pkin(7)) * t966 + (-pkin(2) - pkin(8)) * t935 + t972;
t901 = -t963 * g(3) - t960 * t916;
t878 = -qJDD(2) * pkin(2) - t965 * qJ(3) + t931 * t988 + qJDD(3) - t901;
t859 = (-t960 * t963 * t966 - qJDD(2)) * pkin(8) + (t934 - t983) * pkin(3) + t878;
t841 = -t959 * t848 + t962 * t859;
t928 = qJDD(4) + t934;
t837 = (t929 * t948 - t894) * qJ(5) + (t929 * t930 + t928) * pkin(4) + t841;
t842 = t962 * t848 + t959 * t859;
t839 = -t927 * pkin(4) + t893 * qJ(5) - t948 * t899 + t842;
t833 = t895 * t1005 + t958 * t837 + t997 * t839;
t870 = t895 * pkin(5) - t896 * qJ(6);
t946 = t948 ^ 2;
t829 = -t946 * pkin(5) + t928 * qJ(6) + 0.2e1 * qJD(6) * t948 - t895 * t870 + t833;
t994 = t1012 * t895 + t1013 * t948 - t1015 * t896;
t811 = -mrSges(6,1) * t844 - mrSges(7,1) * t835 + mrSges(7,2) * t829 + mrSges(6,3) * t833 - pkin(5) * t826 + t1001 * t867 + t1009 * t948 + t1012 * t928 - t1018 * t866 + t994 * t896;
t976 = t997 * t837 - t958 * t839;
t830 = -t928 * pkin(5) - t946 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t870) * t896 - t976;
t832 = t896 * t1005 + t976;
t812 = mrSges(6,2) * t844 + mrSges(7,2) * t830 - mrSges(6,3) * t832 - mrSges(7,3) * t835 - qJ(6) * t826 - t1001 * t866 + t1010 * t948 + t1015 * t928 + t1019 * t867 + t994 * t895;
t879 = -t948 * mrSges(6,2) - t895 * mrSges(6,3);
t881 = t948 * mrSges(6,1) - t896 * mrSges(6,3);
t822 = m(6) * t844 + t866 * mrSges(6,1) + t867 * mrSges(6,2) + t895 * t879 + t896 * t881 + t826;
t884 = Ifges(5,5) * t930 + Ifges(5,6) * t929 + Ifges(5,3) * t948;
t886 = Ifges(5,1) * t930 + Ifges(5,4) * t929 + Ifges(5,5) * t948;
t1003 = -mrSges(6,3) - mrSges(7,2);
t985 = m(7) * t829 + t928 * mrSges(7,3) + t948 * t882;
t871 = t895 * mrSges(7,1) - t896 * mrSges(7,3);
t993 = -t895 * mrSges(6,1) - t896 * mrSges(6,2) - t871;
t815 = m(6) * t833 - t928 * mrSges(6,2) + t1003 * t866 - t948 * t881 + t993 * t895 + t985;
t979 = -m(7) * t830 + t928 * mrSges(7,1) + t948 * t880;
t817 = m(6) * t832 + t928 * mrSges(6,1) + t1003 * t867 + t948 * t879 + t993 * t896 + t979;
t980 = t997 * t815 - t958 * t817;
t790 = -mrSges(5,1) * t858 + mrSges(5,3) * t842 + Ifges(5,4) * t894 + Ifges(5,2) * t893 + Ifges(5,6) * t928 - pkin(4) * t822 + qJ(5) * t980 + t997 * t811 + t958 * t812 - t930 * t884 + t948 * t886;
t810 = t958 * t815 + t997 * t817;
t885 = Ifges(5,4) * t930 + Ifges(5,2) * t929 + Ifges(5,6) * t948;
t791 = mrSges(5,2) * t858 - mrSges(5,3) * t841 + Ifges(5,1) * t894 + Ifges(5,4) * t893 + Ifges(5,5) * t928 - qJ(5) * t810 - t958 * t811 + t997 * t812 + t929 * t884 - t948 * t885;
t932 = (mrSges(4,2) * t963 - mrSges(4,3) * t960) * qJD(1);
t941 = -mrSges(4,1) * t989 - qJD(2) * mrSges(4,3);
t897 = -t929 * mrSges(5,1) + t930 * mrSges(5,2);
t898 = -t948 * mrSges(5,2) + t929 * mrSges(5,3);
t807 = m(5) * t841 + t928 * mrSges(5,1) - t894 * mrSges(5,3) - t930 * t897 + t948 * t898 + t810;
t900 = t948 * mrSges(5,1) - t930 * mrSges(5,3);
t808 = m(5) * t842 - t928 * mrSges(5,2) + t893 * mrSges(5,3) + t929 * t897 - t948 * t900 + t980;
t804 = t962 * t807 + t959 * t808;
t974 = -m(4) * t878 - t934 * mrSges(4,1) - t804;
t803 = qJDD(2) * mrSges(4,2) + qJD(2) * t941 + t932 * t988 - t974;
t942 = mrSges(4,1) * t988 + qJD(2) * mrSges(4,2);
t969 = -m(5) * t858 + t893 * mrSges(5,1) - t894 * mrSges(5,2) + t929 * t898 - t930 * t900 - t822;
t968 = -m(4) * t877 + qJDD(2) * mrSges(4,3) + qJD(2) * t942 + t932 * t989 - t969;
t990 = t1000 * qJD(2) + (t1002 * t963 + t1016 * t960) * qJD(1);
t991 = t999 * qJD(2) + (t1002 * t960 + t1014 * t963) * qJD(1);
t1008 = (t991 * t960 - t990 * t963) * qJD(1) + t1011 * qJDD(2) + t1000 * t934 + t999 * t935 + mrSges(3,1) * t901 - mrSges(3,2) * t902 + mrSges(4,2) * t878 - mrSges(4,3) * t877 - pkin(2) * t803 - pkin(8) * t804 + qJ(3) * (t935 * mrSges(4,1) + t968) - t959 * t790 + t962 * t791;
t1004 = t966 * pkin(7);
t933 = (-mrSges(3,1) * t963 + mrSges(3,2) * t960) * qJD(1);
t940 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t989;
t801 = m(3) * t901 - t934 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t940 - t941) * qJD(2) + (-t932 - t933) * t988 + t974;
t939 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t988;
t820 = t968 + (mrSges(3,3) + mrSges(4,1)) * t935 - qJDD(2) * mrSges(3,2) + t933 * t989 - qJD(2) * t939 + m(3) * t902;
t981 = -t960 * t801 + t963 * t820;
t794 = m(2) * t945 - t966 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t981;
t915 = t978 - t1004;
t873 = -t935 * pkin(2) - t1004 + t972;
t995 = -t959 * t807 + t962 * t808;
t977 = -m(4) * t873 - t935 * mrSges(4,2) + t942 * t988 - t995;
t971 = -m(3) * t915 + t940 * t989 + t935 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t934 + (-t939 * t960 - t941 * t963) * qJD(1) + t977;
t798 = m(2) * t944 + qJDD(1) * mrSges(2,1) - t966 * mrSges(2,2) + t971;
t996 = t961 * t794 + t964 * t798;
t796 = t963 * t801 + t960 * t820;
t992 = t1011 * qJD(2) + (t1000 * t960 + t999 * t963) * qJD(1);
t982 = t964 * t794 - t961 * t798;
t802 = -t934 * mrSges(4,3) + t941 * t989 - t977;
t787 = -mrSges(3,1) * t915 - mrSges(4,1) * t877 + mrSges(4,2) * t873 + mrSges(3,3) * t902 - pkin(2) * t802 - pkin(3) * t969 - pkin(8) * t995 + t990 * qJD(2) + t999 * qJDD(2) + t1002 * t934 + t1014 * t935 - t962 * t790 - t959 * t791 - t992 * t988;
t825 = t867 * mrSges(7,2) + t896 * t871 - t979;
t967 = mrSges(5,1) * t841 + mrSges(6,1) * t832 - mrSges(7,1) * t830 - mrSges(5,2) * t842 - mrSges(6,2) * t833 + mrSges(7,3) * t829 + Ifges(5,5) * t894 + Ifges(5,6) * t893 + pkin(4) * t810 - pkin(5) * t825 + qJ(6) * t985 + t930 * t885 - t929 * t886 - t1010 * t896 + (-qJ(6) * t871 + t1009) * t895 + t1015 * t867 + (-mrSges(7,2) * qJ(6) - t1012) * t866 + (Ifges(5,3) - t1013) * t928;
t789 = mrSges(4,1) * t878 + mrSges(3,2) * t915 - mrSges(3,3) * t901 - mrSges(4,3) * t873 + pkin(3) * t804 - qJ(3) * t802 - t991 * qJD(2) + t1000 * qJDD(2) + t1002 * t935 + t1016 * t934 + t992 * t989 + t967;
t973 = mrSges(2,1) * t944 - mrSges(2,2) * t945 + Ifges(2,3) * qJDD(1) + pkin(1) * t971 + pkin(7) * t981 + t963 * t787 + t960 * t789;
t785 = mrSges(2,1) * g(3) + mrSges(2,3) * t945 + t966 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t796 - t1008;
t784 = -mrSges(2,2) * g(3) - mrSges(2,3) * t944 + Ifges(2,5) * qJDD(1) - t966 * Ifges(2,6) - pkin(7) * t796 - t960 * t787 + t963 * t789;
t1 = [-m(1) * g(1) + t982; -m(1) * g(2) + t996; (-m(1) - m(2)) * g(3) + t796; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t996 + t964 * t784 - t961 * t785; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t982 + t961 * t784 + t964 * t785; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t973; t973; t1008; t803; t967; t822; t825;];
tauJB  = t1;
