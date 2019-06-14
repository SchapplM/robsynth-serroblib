% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP5
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
% Datum: 2019-05-07 18:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:24:36
% EndTime: 2019-05-07 18:24:51
% DurationCPUTime: 8.54s
% Computational Cost: add. (107669->340), mult. (214142->401), div. (0->0), fcn. (147258->8), ass. (0->132)
t1004 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t994 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t978 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t951 = sin(qJ(3));
t954 = cos(qJ(3));
t952 = sin(qJ(2));
t981 = qJD(1) * t952;
t930 = t954 * qJD(2) - t951 * t981;
t931 = t951 * qJD(2) + t954 * t981;
t950 = sin(qJ(4));
t988 = cos(qJ(4));
t903 = -t988 * t930 + t950 * t931;
t904 = t950 * t930 + t988 * t931;
t955 = cos(qJ(2));
t980 = t955 * qJD(1);
t945 = qJD(3) - t980;
t943 = -qJD(4) - t945;
t1003 = t1004 * t903 + t978 * t904 - t994 * t943;
t1002 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t996 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t887 = -t943 * mrSges(7,2) + t903 * mrSges(7,3);
t979 = qJD(1) * qJD(2);
t946 = t952 * t979;
t935 = t955 * qJDD(1) - t946;
t929 = qJDD(3) - t935;
t925 = -qJDD(4) - t929;
t953 = sin(qJ(1));
t956 = cos(qJ(1));
t940 = t953 * g(1) - t956 * g(2);
t958 = qJD(1) ^ 2;
t923 = -qJDD(1) * pkin(1) - t958 * pkin(7) - t940;
t975 = t955 * t979;
t934 = t952 * qJDD(1) + t975;
t880 = (-t934 - t975) * pkin(8) + (-t935 + t946) * pkin(2) + t923;
t941 = -t956 * g(1) - t953 * g(2);
t924 = -t958 * pkin(1) + qJDD(1) * pkin(7) + t941;
t909 = -t952 * g(3) + t955 * t924;
t933 = (-pkin(2) * t955 - pkin(8) * t952) * qJD(1);
t957 = qJD(2) ^ 2;
t885 = -t957 * pkin(2) + qJDD(2) * pkin(8) + t933 * t980 + t909;
t854 = t954 * t880 - t951 * t885;
t901 = t930 * qJD(3) + t951 * qJDD(2) + t954 * t934;
t833 = (t930 * t945 - t901) * pkin(9) + (t930 * t931 + t929) * pkin(3) + t854;
t855 = t951 * t880 + t954 * t885;
t900 = -t931 * qJD(3) + t954 * qJDD(2) - t951 * t934;
t910 = t945 * pkin(3) - t931 * pkin(9);
t928 = t930 ^ 2;
t836 = -t928 * pkin(3) + t900 * pkin(9) - t945 * t910 + t855;
t830 = t988 * t833 - t950 * t836;
t874 = t903 * pkin(4) - t904 * qJ(5);
t942 = t943 ^ 2;
t826 = t925 * pkin(4) - t942 * qJ(5) + t904 * t874 + qJDD(5) - t830;
t853 = -t903 * qJD(4) + t950 * t900 + t988 * t901;
t985 = t903 * t943;
t817 = -0.2e1 * qJD(6) * t904 + (-t853 + t985) * qJ(6) + (t903 * t904 + t925) * pkin(5) + t826;
t876 = -t903 * mrSges(7,1) + t904 * mrSges(7,2);
t969 = -m(7) * t817 + t853 * mrSges(7,3) + t904 * t876;
t815 = t925 * mrSges(7,1) + t943 * t887 - t969;
t875 = t903 * mrSges(6,1) - t904 * mrSges(6,3);
t886 = -t903 * mrSges(6,2) - t943 * mrSges(6,3);
t992 = -m(6) * t826 - t925 * mrSges(6,1) - t943 * t886;
t812 = t853 * mrSges(6,2) + t904 * t875 + t815 - t992;
t831 = t950 * t833 + t988 * t836;
t989 = -2 * qJD(5);
t825 = -t942 * pkin(4) - t925 * qJ(5) - t903 * t874 + t943 * t989 + t831;
t852 = t904 * qJD(4) - t988 * t900 + t950 * t901;
t889 = t943 * pkin(5) - t904 * qJ(6);
t902 = t903 ^ 2;
t820 = -t902 * pkin(5) + t852 * qJ(6) + 0.2e1 * qJD(6) * t903 - t943 * t889 + t825;
t890 = t943 * mrSges(7,1) - t904 * mrSges(7,3);
t892 = t943 * mrSges(6,1) + t904 * mrSges(6,2);
t977 = m(7) * t820 + t852 * mrSges(7,3) + t903 * t876;
t967 = m(6) * t825 - t925 * mrSges(6,3) - t943 * t892 + t977;
t993 = t1002 * t904 - t978 * t903 - t996 * t943;
t995 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t1001 = -mrSges(6,1) * t826 - mrSges(7,1) * t817 - mrSges(5,2) * t831 - pkin(4) * t812 - pkin(5) * t815 + qJ(5) * (-t943 * t890 + t967) + mrSges(7,2) * t820 + mrSges(6,3) * t825 + mrSges(5,1) * t830 + (-qJ(5) * mrSges(7,2) - t995) * t925 + (-qJ(5) * t875 + t993) * t903 + t996 * t853 + (-qJ(5) * mrSges(6,2) - t994) * t852 + t1003 * t904;
t888 = t943 * mrSges(5,2) - t903 * mrSges(5,3);
t982 = -t903 * mrSges(5,1) - t904 * mrSges(5,2) - t875;
t986 = -mrSges(5,3) - mrSges(6,2);
t806 = m(5) * t830 + (-t887 - t888) * t943 + (-mrSges(5,1) - mrSges(7,1)) * t925 + t982 * t904 + t986 * t853 + t969 + t992;
t891 = -t943 * mrSges(5,1) - t904 * mrSges(5,3);
t809 = m(5) * t831 + (-t890 + t891) * t943 + (mrSges(5,2) - mrSges(7,2)) * t925 + t982 * t903 + t986 * t852 + t967;
t801 = t988 * t806 + t950 * t809;
t895 = Ifges(4,4) * t931 + Ifges(4,2) * t930 + Ifges(4,6) * t945;
t896 = Ifges(4,1) * t931 + Ifges(4,4) * t930 + Ifges(4,5) * t945;
t997 = mrSges(4,1) * t854 - mrSges(4,2) * t855 + Ifges(4,5) * t901 + Ifges(4,6) * t900 + Ifges(4,3) * t929 + pkin(3) * t801 + t931 * t895 - t930 * t896 + t1001;
t908 = -t955 * g(3) - t952 * t924;
t966 = qJDD(2) * pkin(2) + t957 * pkin(8) - t933 * t981 + t908;
t964 = t900 * pkin(3) + t928 * pkin(9) - t931 * t910 + t966;
t991 = -(t853 + t985) * qJ(5) - t964;
t822 = -t902 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t852 + (pkin(4) * t943 + (2 * qJD(5)) + t889) * t904 - t991;
t816 = m(7) * t822 - t852 * mrSges(7,1) + t853 * mrSges(7,2) - t903 * t887 + t904 * t890;
t828 = t904 * t989 + (-t904 * t943 + t852) * pkin(4) + t991;
t810 = m(6) * t828 + t852 * mrSges(6,1) - t853 * mrSges(6,3) + t903 * t886 - t904 * t892 - t816;
t976 = t903 * t994 - t904 * t996 + t943 * t995;
t796 = mrSges(5,1) * t964 + mrSges(5,3) * t831 - mrSges(6,1) * t828 + mrSges(6,2) * t825 + mrSges(7,1) * t822 - mrSges(7,3) * t820 + pkin(5) * t816 - qJ(6) * t977 - pkin(4) * t810 + (qJ(6) * t890 - t993) * t943 + (qJ(6) * mrSges(7,2) - t994) * t925 + t976 * t904 + t978 * t853 + t1004 * t852;
t797 = -mrSges(5,2) * t964 + mrSges(6,2) * t826 + mrSges(7,2) * t822 - mrSges(5,3) * t830 - mrSges(6,3) * t828 - mrSges(7,3) * t817 - qJ(5) * t810 - qJ(6) * t815 + t1002 * t853 + t1003 * t943 - t978 * t852 + t976 * t903 - t996 * t925;
t894 = Ifges(4,5) * t931 + Ifges(4,6) * t930 + Ifges(4,3) * t945;
t962 = -m(5) * t964 + t852 * mrSges(5,1) + t853 * mrSges(5,2) + t903 * t888 + t904 * t891 + t810;
t971 = -t950 * t806 + t988 * t809;
t781 = mrSges(4,1) * t966 + mrSges(4,3) * t855 + Ifges(4,4) * t901 + Ifges(4,2) * t900 + Ifges(4,6) * t929 - pkin(3) * t962 + pkin(9) * t971 + t988 * t796 + t950 * t797 - t931 * t894 + t945 * t896;
t782 = -mrSges(4,2) * t966 - mrSges(4,3) * t854 + Ifges(4,1) * t901 + Ifges(4,4) * t900 + Ifges(4,5) * t929 - pkin(9) * t801 - t950 * t796 + t988 * t797 + t930 * t894 - t945 * t895;
t905 = -t930 * mrSges(4,1) + t931 * mrSges(4,2);
t906 = -t945 * mrSges(4,2) + t930 * mrSges(4,3);
t799 = m(4) * t854 + t929 * mrSges(4,1) - t901 * mrSges(4,3) - t931 * t905 + t945 * t906 + t801;
t907 = t945 * mrSges(4,1) - t931 * mrSges(4,3);
t800 = m(4) * t855 - t929 * mrSges(4,2) + t900 * mrSges(4,3) + t930 * t905 - t945 * t907 + t971;
t795 = -t951 * t799 + t954 * t800;
t804 = m(4) * t966 + t900 * mrSges(4,1) - t901 * mrSges(4,2) + t930 * t906 - t931 * t907 - t962;
t921 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t952 + Ifges(3,2) * t955) * qJD(1);
t922 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t952 + Ifges(3,4) * t955) * qJD(1);
t990 = mrSges(3,1) * t908 - mrSges(3,2) * t909 + Ifges(3,5) * t934 + Ifges(3,6) * t935 + Ifges(3,3) * qJDD(2) + pkin(2) * t804 + pkin(8) * t795 + t954 * t781 + t951 * t782 + (t952 * t921 - t955 * t922) * qJD(1);
t932 = (-mrSges(3,1) * t955 + mrSges(3,2) * t952) * qJD(1);
t938 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t981;
t793 = m(3) * t909 - qJDD(2) * mrSges(3,2) + t935 * mrSges(3,3) - qJD(2) * t938 + t932 * t980 + t795;
t939 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t980;
t803 = m(3) * t908 + qJDD(2) * mrSges(3,1) - t934 * mrSges(3,3) + qJD(2) * t939 - t932 * t981 + t804;
t972 = t955 * t793 - t952 * t803;
t785 = m(2) * t941 - t958 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t972;
t794 = t954 * t799 + t951 * t800;
t963 = -m(3) * t923 + t935 * mrSges(3,1) - t934 * mrSges(3,2) - t938 * t981 + t939 * t980 - t794;
t789 = m(2) * t940 + qJDD(1) * mrSges(2,1) - t958 * mrSges(2,2) + t963;
t984 = t953 * t785 + t956 * t789;
t787 = t952 * t793 + t955 * t803;
t973 = t956 * t785 - t953 * t789;
t920 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t952 + Ifges(3,6) * t955) * qJD(1);
t778 = mrSges(3,2) * t923 - mrSges(3,3) * t908 + Ifges(3,1) * t934 + Ifges(3,4) * t935 + Ifges(3,5) * qJDD(2) - pkin(8) * t794 - qJD(2) * t921 - t951 * t781 + t954 * t782 + t920 * t980;
t780 = -mrSges(3,1) * t923 + mrSges(3,3) * t909 + Ifges(3,4) * t934 + Ifges(3,2) * t935 + Ifges(3,6) * qJDD(2) - pkin(2) * t794 + qJD(2) * t922 - t920 * t981 - t997;
t965 = mrSges(2,1) * t940 - mrSges(2,2) * t941 + Ifges(2,3) * qJDD(1) + pkin(1) * t963 + pkin(7) * t972 + t952 * t778 + t955 * t780;
t776 = mrSges(2,1) * g(3) + mrSges(2,3) * t941 + t958 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t787 - t990;
t775 = -mrSges(2,2) * g(3) - mrSges(2,3) * t940 + Ifges(2,5) * qJDD(1) - t958 * Ifges(2,6) - pkin(7) * t787 + t955 * t778 - t952 * t780;
t1 = [-m(1) * g(1) + t973; -m(1) * g(2) + t984; (-m(1) - m(2)) * g(3) + t787; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t984 + t956 * t775 - t953 * t776; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t973 + t953 * t775 + t956 * t776; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t965; t965; t990; t997; t1001; t812; t816;];
tauJB  = t1;
