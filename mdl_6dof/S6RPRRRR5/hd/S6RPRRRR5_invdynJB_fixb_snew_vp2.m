% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:31:26
% EndTime: 2019-05-06 03:31:56
% DurationCPUTime: 30.64s
% Computational Cost: add. (500253->367), mult. (1196346->456), div. (0->0), fcn. (956008->12), ass. (0->156)
t962 = qJD(1) ^ 2;
t951 = cos(pkin(11));
t994 = pkin(2) * t951;
t950 = sin(pkin(11));
t993 = mrSges(3,2) * t950;
t947 = t951 ^ 2;
t992 = t947 * t962;
t956 = sin(qJ(1));
t961 = cos(qJ(1));
t936 = -g(1) * t961 - g(2) * t956;
t931 = -pkin(1) * t962 + qJDD(1) * qJ(2) + t936;
t988 = qJD(1) * qJD(2);
t986 = -t951 * g(3) - 0.2e1 * t950 * t988;
t902 = (-pkin(7) * qJDD(1) + t962 * t994 - t931) * t950 + t986;
t921 = -g(3) * t950 + (t931 + 0.2e1 * t988) * t951;
t987 = qJDD(1) * t951;
t903 = -pkin(2) * t992 + pkin(7) * t987 + t921;
t955 = sin(qJ(3));
t960 = cos(qJ(3));
t882 = t960 * t902 - t955 * t903;
t975 = t950 * t960 + t951 * t955;
t974 = -t950 * t955 + t951 * t960;
t929 = t974 * qJD(1);
t989 = t929 * qJD(3);
t919 = t975 * qJDD(1) + t989;
t930 = t975 * qJD(1);
t853 = (-t919 + t989) * pkin(8) + (t929 * t930 + qJDD(3)) * pkin(3) + t882;
t883 = t955 * t902 + t960 * t903;
t918 = -t930 * qJD(3) + t974 * qJDD(1);
t924 = qJD(3) * pkin(3) - pkin(8) * t930;
t928 = t929 ^ 2;
t859 = -pkin(3) * t928 + pkin(8) * t918 - qJD(3) * t924 + t883;
t954 = sin(qJ(4));
t959 = cos(qJ(4));
t846 = t954 * t853 + t959 * t859;
t911 = t929 * t954 + t930 * t959;
t877 = -t911 * qJD(4) + t918 * t959 - t954 * t919;
t910 = t929 * t959 - t954 * t930;
t892 = -mrSges(5,1) * t910 + mrSges(5,2) * t911;
t948 = qJD(3) + qJD(4);
t900 = mrSges(5,1) * t948 - mrSges(5,3) * t911;
t945 = qJDD(3) + qJDD(4);
t893 = -pkin(4) * t910 - pkin(9) * t911;
t944 = t948 ^ 2;
t834 = -pkin(4) * t944 + pkin(9) * t945 + t893 * t910 + t846;
t946 = t950 ^ 2;
t935 = t956 * g(1) - t961 * g(2);
t979 = qJDD(2) - t935;
t917 = (-pkin(1) - t994) * qJDD(1) + (-qJ(2) + (-t946 - t947) * pkin(7)) * t962 + t979;
t868 = -pkin(3) * t918 - pkin(8) * t928 + t930 * t924 + t917;
t878 = qJD(4) * t910 + t918 * t954 + t919 * t959;
t842 = (-t910 * t948 - t878) * pkin(9) + (t911 * t948 - t877) * pkin(4) + t868;
t953 = sin(qJ(5));
t958 = cos(qJ(5));
t829 = -t953 * t834 + t958 * t842;
t895 = -t911 * t953 + t948 * t958;
t856 = qJD(5) * t895 + t878 * t958 + t945 * t953;
t876 = qJDD(5) - t877;
t896 = t911 * t958 + t948 * t953;
t906 = qJD(5) - t910;
t827 = (t895 * t906 - t856) * pkin(10) + (t895 * t896 + t876) * pkin(5) + t829;
t830 = t958 * t834 + t953 * t842;
t855 = -qJD(5) * t896 - t878 * t953 + t945 * t958;
t886 = pkin(5) * t906 - pkin(10) * t896;
t894 = t895 ^ 2;
t828 = -pkin(5) * t894 + pkin(10) * t855 - t886 * t906 + t830;
t952 = sin(qJ(6));
t957 = cos(qJ(6));
t825 = t827 * t957 - t828 * t952;
t879 = t895 * t957 - t896 * t952;
t839 = qJD(6) * t879 + t855 * t952 + t856 * t957;
t880 = t895 * t952 + t896 * t957;
t851 = -mrSges(7,1) * t879 + mrSges(7,2) * t880;
t904 = qJD(6) + t906;
t860 = -mrSges(7,2) * t904 + mrSges(7,3) * t879;
t871 = qJDD(6) + t876;
t820 = m(7) * t825 + mrSges(7,1) * t871 - mrSges(7,3) * t839 - t851 * t880 + t860 * t904;
t826 = t827 * t952 + t828 * t957;
t838 = -qJD(6) * t880 + t855 * t957 - t856 * t952;
t861 = mrSges(7,1) * t904 - mrSges(7,3) * t880;
t821 = m(7) * t826 - mrSges(7,2) * t871 + mrSges(7,3) * t838 + t851 * t879 - t861 * t904;
t812 = t957 * t820 + t952 * t821;
t881 = -mrSges(6,1) * t895 + mrSges(6,2) * t896;
t884 = -mrSges(6,2) * t906 + mrSges(6,3) * t895;
t810 = m(6) * t829 + mrSges(6,1) * t876 - mrSges(6,3) * t856 - t881 * t896 + t884 * t906 + t812;
t885 = mrSges(6,1) * t906 - mrSges(6,3) * t896;
t980 = -t820 * t952 + t957 * t821;
t811 = m(6) * t830 - mrSges(6,2) * t876 + mrSges(6,3) * t855 + t881 * t895 - t885 * t906 + t980;
t981 = -t810 * t953 + t958 * t811;
t803 = m(5) * t846 - mrSges(5,2) * t945 + mrSges(5,3) * t877 + t892 * t910 - t900 * t948 + t981;
t845 = t853 * t959 - t954 * t859;
t899 = -mrSges(5,2) * t948 + mrSges(5,3) * t910;
t833 = -pkin(4) * t945 - pkin(9) * t944 + t911 * t893 - t845;
t831 = -pkin(5) * t855 - pkin(10) * t894 + t886 * t896 + t833;
t970 = m(7) * t831 - t838 * mrSges(7,1) + mrSges(7,2) * t839 - t879 * t860 + t861 * t880;
t966 = -m(6) * t833 + t855 * mrSges(6,1) - mrSges(6,2) * t856 + t895 * t884 - t885 * t896 - t970;
t816 = m(5) * t845 + mrSges(5,1) * t945 - mrSges(5,3) * t878 - t892 * t911 + t899 * t948 + t966;
t792 = t954 * t803 + t959 * t816;
t915 = -mrSges(4,1) * t929 + mrSges(4,2) * t930;
t922 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t929;
t790 = m(4) * t882 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t919 + qJD(3) * t922 - t915 * t930 + t792;
t923 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t930;
t982 = t959 * t803 - t816 * t954;
t791 = m(4) * t883 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t918 - qJD(3) * t923 + t915 * t929 + t982;
t785 = t960 * t790 + t955 * t791;
t920 = -t950 * t931 + t986;
t973 = mrSges(3,3) * qJDD(1) + t962 * (-mrSges(3,1) * t951 + t993);
t783 = m(3) * t920 - t973 * t950 + t785;
t983 = -t955 * t790 + t960 * t791;
t784 = m(3) * t921 + t973 * t951 + t983;
t984 = -t783 * t950 + t951 * t784;
t776 = m(2) * t936 - mrSges(2,1) * t962 - qJDD(1) * mrSges(2,2) + t984;
t927 = -qJDD(1) * pkin(1) - t962 * qJ(2) + t979;
t805 = t958 * t810 + t953 * t811;
t972 = m(5) * t868 - t877 * mrSges(5,1) + t878 * mrSges(5,2) - t910 * t899 + t911 * t900 + t805;
t967 = m(4) * t917 - t918 * mrSges(4,1) + mrSges(4,2) * t919 - t929 * t922 + t923 * t930 + t972;
t965 = -m(3) * t927 + mrSges(3,1) * t987 - t967 + (t946 * t962 + t992) * mrSges(3,3);
t798 = t965 + m(2) * t935 - mrSges(2,2) * t962 + (mrSges(2,1) - t993) * qJDD(1);
t991 = t956 * t776 + t961 * t798;
t778 = t951 * t783 + t950 * t784;
t976 = Ifges(3,5) * t950 + Ifges(3,6) * t951;
t990 = t962 * t976;
t985 = t961 * t776 - t798 * t956;
t978 = Ifges(3,1) * t950 + Ifges(3,4) * t951;
t977 = Ifges(3,4) * t950 + Ifges(3,2) * t951;
t847 = Ifges(7,5) * t880 + Ifges(7,6) * t879 + Ifges(7,3) * t904;
t849 = Ifges(7,1) * t880 + Ifges(7,4) * t879 + Ifges(7,5) * t904;
t813 = -mrSges(7,1) * t831 + mrSges(7,3) * t826 + Ifges(7,4) * t839 + Ifges(7,2) * t838 + Ifges(7,6) * t871 - t847 * t880 + t849 * t904;
t848 = Ifges(7,4) * t880 + Ifges(7,2) * t879 + Ifges(7,6) * t904;
t814 = mrSges(7,2) * t831 - mrSges(7,3) * t825 + Ifges(7,1) * t839 + Ifges(7,4) * t838 + Ifges(7,5) * t871 + t847 * t879 - t848 * t904;
t862 = Ifges(6,5) * t896 + Ifges(6,6) * t895 + Ifges(6,3) * t906;
t864 = Ifges(6,1) * t896 + Ifges(6,4) * t895 + Ifges(6,5) * t906;
t794 = -mrSges(6,1) * t833 + mrSges(6,3) * t830 + Ifges(6,4) * t856 + Ifges(6,2) * t855 + Ifges(6,6) * t876 - pkin(5) * t970 + pkin(10) * t980 + t957 * t813 + t952 * t814 - t896 * t862 + t906 * t864;
t863 = Ifges(6,4) * t896 + Ifges(6,2) * t895 + Ifges(6,6) * t906;
t796 = mrSges(6,2) * t833 - mrSges(6,3) * t829 + Ifges(6,1) * t856 + Ifges(6,4) * t855 + Ifges(6,5) * t876 - pkin(10) * t812 - t813 * t952 + t814 * t957 + t862 * t895 - t863 * t906;
t887 = Ifges(5,5) * t911 + Ifges(5,6) * t910 + Ifges(5,3) * t948;
t888 = Ifges(5,4) * t911 + Ifges(5,2) * t910 + Ifges(5,6) * t948;
t779 = mrSges(5,2) * t868 - mrSges(5,3) * t845 + Ifges(5,1) * t878 + Ifges(5,4) * t877 + Ifges(5,5) * t945 - pkin(9) * t805 - t794 * t953 + t796 * t958 + t887 * t910 - t888 * t948;
t889 = Ifges(5,1) * t911 + Ifges(5,4) * t910 + Ifges(5,5) * t948;
t969 = -mrSges(7,1) * t825 + mrSges(7,2) * t826 - Ifges(7,5) * t839 - Ifges(7,6) * t838 - Ifges(7,3) * t871 - t880 * t848 + t879 * t849;
t964 = mrSges(6,1) * t829 - mrSges(6,2) * t830 + Ifges(6,5) * t856 + Ifges(6,6) * t855 + Ifges(6,3) * t876 + pkin(5) * t812 + t896 * t863 - t895 * t864 - t969;
t786 = -mrSges(5,1) * t868 + mrSges(5,3) * t846 + Ifges(5,4) * t878 + Ifges(5,2) * t877 + Ifges(5,6) * t945 - pkin(4) * t805 - t911 * t887 + t948 * t889 - t964;
t907 = Ifges(4,5) * t930 + Ifges(4,6) * t929 + Ifges(4,3) * qJD(3);
t909 = Ifges(4,1) * t930 + Ifges(4,4) * t929 + Ifges(4,5) * qJD(3);
t772 = -mrSges(4,1) * t917 + mrSges(4,3) * t883 + Ifges(4,4) * t919 + Ifges(4,2) * t918 + Ifges(4,6) * qJDD(3) - pkin(3) * t972 + pkin(8) * t982 + qJD(3) * t909 + t954 * t779 + t959 * t786 - t930 * t907;
t908 = Ifges(4,4) * t930 + Ifges(4,2) * t929 + Ifges(4,6) * qJD(3);
t773 = mrSges(4,2) * t917 - mrSges(4,3) * t882 + Ifges(4,1) * t919 + Ifges(4,4) * t918 + Ifges(4,5) * qJDD(3) - pkin(8) * t792 - qJD(3) * t908 + t779 * t959 - t786 * t954 + t907 * t929;
t768 = -mrSges(3,1) * t927 + mrSges(3,3) * t921 - pkin(2) * t967 + pkin(7) * t983 + t977 * qJDD(1) + t960 * t772 + t955 * t773 - t950 * t990;
t770 = mrSges(3,2) * t927 - mrSges(3,3) * t920 - pkin(7) * t785 + t978 * qJDD(1) - t772 * t955 + t773 * t960 + t951 * t990;
t800 = qJDD(1) * t993 - t965;
t971 = mrSges(2,1) * t935 - mrSges(2,2) * t936 + Ifges(2,3) * qJDD(1) - pkin(1) * t800 + qJ(2) * t984 + t951 * t768 + t950 * t770;
t968 = -mrSges(5,1) * t845 + mrSges(5,2) * t846 - Ifges(5,5) * t878 - Ifges(5,6) * t877 - Ifges(5,3) * t945 - pkin(4) * t966 - pkin(9) * t981 - t958 * t794 - t953 * t796 - t911 * t888 + t910 * t889;
t963 = mrSges(4,1) * t882 - mrSges(4,2) * t883 + Ifges(4,5) * t919 + Ifges(4,6) * t918 + Ifges(4,3) * qJDD(3) + pkin(3) * t792 + t930 * t908 - t929 * t909 - t968;
t771 = -pkin(2) * t785 - t963 + mrSges(2,1) * g(3) + (Ifges(2,6) - t976) * qJDD(1) + mrSges(2,3) * t936 - mrSges(3,1) * t920 + mrSges(3,2) * t921 - pkin(1) * t778 + (-t950 * t977 + t951 * t978 + Ifges(2,5)) * t962;
t766 = -mrSges(2,2) * g(3) - mrSges(2,3) * t935 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t962 - qJ(2) * t778 - t768 * t950 + t770 * t951;
t1 = [-m(1) * g(1) + t985; -m(1) * g(2) + t991; (-m(1) - m(2)) * g(3) + t778; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t991 + t961 * t766 - t956 * t771; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t985 + t956 * t766 + t961 * t771; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t971; t971; t800; t963; -t968; t964; -t969;];
tauJB  = t1;
