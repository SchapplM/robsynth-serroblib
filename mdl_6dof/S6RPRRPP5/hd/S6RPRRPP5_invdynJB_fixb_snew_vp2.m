% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:36:46
% EndTime: 2019-05-05 21:36:56
% DurationCPUTime: 7.02s
% Computational Cost: add. (73785->320), mult. (172788->376), div. (0->0), fcn. (123884->8), ass. (0->136)
t992 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t968 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t967 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t991 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t966 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t990 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t940 = qJD(1) ^ 2;
t937 = sin(qJ(1));
t938 = cos(qJ(1));
t919 = -t938 * g(1) - t937 * g(2);
t914 = -t940 * pkin(1) + qJDD(1) * qJ(2) + t919;
t933 = sin(pkin(9));
t934 = cos(pkin(9));
t971 = qJD(1) * qJD(2);
t959 = -t934 * g(3) - 0.2e1 * t933 * t971;
t981 = pkin(2) * t934;
t884 = (-pkin(7) * qJDD(1) + t940 * t981 - t914) * t933 + t959;
t901 = -t933 * g(3) + (t914 + 0.2e1 * t971) * t934;
t969 = qJDD(1) * t934;
t931 = t934 ^ 2;
t977 = t931 * t940;
t885 = -pkin(2) * t977 + pkin(7) * t969 + t901;
t936 = sin(qJ(3));
t984 = cos(qJ(3));
t840 = t936 * t884 + t984 * t885;
t960 = t934 * t984;
t972 = t933 * qJD(1);
t912 = qJD(1) * t960 - t936 * t972;
t947 = t984 * t933 + t934 * t936;
t913 = t947 * qJD(1);
t896 = -t912 * pkin(3) - t913 * pkin(8);
t939 = qJD(3) ^ 2;
t836 = -t939 * pkin(3) + qJDD(3) * pkin(8) + t912 * t896 + t840;
t930 = t933 ^ 2;
t918 = t937 * g(1) - t938 * g(2);
t953 = qJDD(2) - t918;
t897 = (-pkin(1) - t981) * qJDD(1) + (-qJ(2) + (-t930 - t931) * pkin(7)) * t940 + t953;
t970 = qJDD(1) * t933;
t973 = t913 * qJD(3);
t898 = qJDD(1) * t960 - t936 * t970 - t973;
t974 = t912 * qJD(3);
t899 = qJDD(1) * t947 + t974;
t838 = (-t899 - t974) * pkin(8) + (-t898 + t973) * pkin(3) + t897;
t935 = sin(qJ(4));
t983 = cos(qJ(4));
t831 = -t935 * t836 + t983 * t838;
t903 = -t983 * qJD(3) + t935 * t913;
t904 = t935 * qJD(3) + t983 * t913;
t868 = t903 * pkin(4) - t904 * qJ(5);
t895 = qJDD(4) - t898;
t910 = qJD(4) - t912;
t909 = t910 ^ 2;
t829 = -t895 * pkin(4) - t909 * qJ(5) + t904 * t868 + qJDD(5) - t831;
t875 = -t903 * mrSges(6,2) + t910 * mrSges(6,3);
t989 = -m(6) * t829 + t895 * mrSges(6,1) + t910 * t875;
t865 = -t903 * qJD(4) + t935 * qJDD(3) + t983 * t899;
t839 = t984 * t884 - t936 * t885;
t946 = qJDD(3) * pkin(3) + t939 * pkin(8) - t913 * t896 + t839;
t978 = t903 * t910;
t988 = (-t865 + t978) * qJ(5) - t946;
t876 = t910 * mrSges(7,2) + t903 * mrSges(7,3);
t986 = -0.2e1 * t904;
t822 = qJD(6) * t986 + (-t865 - t978) * qJ(6) + (t903 * t904 - t895) * pkin(5) + t829;
t870 = -t903 * mrSges(7,1) + t904 * mrSges(7,2);
t954 = -m(7) * t822 + t865 * mrSges(7,3) + t904 * t870;
t820 = -t895 * mrSges(7,1) - t910 * t876 - t954;
t869 = t903 * mrSges(6,1) - t904 * mrSges(6,3);
t817 = t865 * mrSges(6,2) + t904 * t869 + t820 - t989;
t832 = t983 * t836 + t935 * t838;
t985 = 2 * qJD(5);
t828 = -t909 * pkin(4) + t895 * qJ(5) - t903 * t868 + t910 * t985 + t832;
t864 = t904 * qJD(4) - t983 * qJDD(3) + t935 * t899;
t878 = -t910 * pkin(5) - t904 * qJ(6);
t902 = t903 ^ 2;
t824 = -t902 * pkin(5) + t864 * qJ(6) + 0.2e1 * qJD(6) * t903 + t910 * t878 + t828;
t879 = -t910 * mrSges(7,1) - t904 * mrSges(7,3);
t881 = -t910 * mrSges(6,1) + t904 * mrSges(6,2);
t964 = m(7) * t824 + t864 * mrSges(7,3) + t903 * t870;
t948 = m(6) * t828 + t895 * mrSges(6,3) + t910 * t881 + t964;
t961 = -t968 * t903 + t904 * t992 + t967 * t910;
t962 = t903 * t991 + t904 * t968 + t910 * t966;
t987 = -t966 * t864 + t967 * t865 + t990 * t895 + t961 * t903 + t962 * t904 + mrSges(5,1) * t831 - mrSges(6,1) * t829 - mrSges(7,1) * t822 - mrSges(5,2) * t832 + mrSges(7,2) * t824 + mrSges(6,3) * t828 - pkin(4) * t817 - pkin(5) * t820 + qJ(5) * (-t864 * mrSges(6,2) + t895 * mrSges(7,2) - t903 * t869 + t910 * t879 + t948);
t980 = -mrSges(5,3) - mrSges(6,2);
t979 = mrSges(3,2) * t933;
t877 = -t910 * mrSges(5,2) - t903 * mrSges(5,3);
t975 = -t903 * mrSges(5,1) - t904 * mrSges(5,2) - t869;
t812 = m(5) * t831 + (t876 + t877) * t910 + t975 * t904 + (mrSges(5,1) + mrSges(7,1)) * t895 + t980 * t865 + t954 + t989;
t880 = t910 * mrSges(5,1) - t904 * mrSges(5,3);
t815 = m(5) * t832 + (t879 - t880) * t910 + t975 * t903 + (-mrSges(5,2) + mrSges(7,2)) * t895 + t980 * t864 + t948;
t808 = -t935 * t812 + t983 * t815;
t893 = -t912 * mrSges(4,1) + t913 * mrSges(4,2);
t906 = qJD(3) * mrSges(4,1) - t913 * mrSges(4,3);
t806 = m(4) * t840 - qJDD(3) * mrSges(4,2) + t898 * mrSges(4,3) - qJD(3) * t906 + t912 * t893 + t808;
t826 = -t902 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t864 + (-pkin(4) * t910 + t878 + t985) * t904 - t988;
t821 = m(7) * t826 - t864 * mrSges(7,1) + t865 * mrSges(7,2) - t903 * t876 + t904 * t879;
t830 = qJD(5) * t986 + (t904 * t910 + t864) * pkin(4) + t988;
t819 = m(6) * t830 + t864 * mrSges(6,1) - t865 * mrSges(6,3) + t903 * t875 - t904 * t881 - t821;
t816 = m(5) * t946 - t864 * mrSges(5,1) - t865 * mrSges(5,2) - t903 * t877 - t904 * t880 - t819;
t905 = -qJD(3) * mrSges(4,2) + t912 * mrSges(4,3);
t810 = m(4) * t839 + qJDD(3) * mrSges(4,1) - t899 * mrSges(4,3) + qJD(3) * t905 - t913 * t893 + t816;
t798 = t936 * t806 + t984 * t810;
t900 = -t933 * t914 + t959;
t949 = mrSges(3,3) * qJDD(1) + t940 * (-mrSges(3,1) * t934 + t979);
t796 = m(3) * t900 - t933 * t949 + t798;
t956 = t984 * t806 - t936 * t810;
t797 = m(3) * t901 + t934 * t949 + t956;
t957 = -t933 * t796 + t934 * t797;
t788 = m(2) * t919 - t940 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t957;
t911 = -qJDD(1) * pkin(1) - t940 * qJ(2) + t953;
t807 = t983 * t812 + t935 * t815;
t944 = m(4) * t897 - t898 * mrSges(4,1) + t899 * mrSges(4,2) - t912 * t905 + t913 * t906 + t807;
t943 = -m(3) * t911 + mrSges(3,1) * t969 - t944 + (t930 * t940 + t977) * mrSges(3,3);
t801 = t943 - t940 * mrSges(2,2) + m(2) * t918 + (mrSges(2,1) - t979) * qJDD(1);
t976 = t937 * t788 + t938 * t801;
t790 = t934 * t796 + t933 * t797;
t963 = t903 * t966 - t904 * t967 - t910 * t990;
t958 = t938 * t788 - t937 * t801;
t952 = Ifges(3,1) * t933 + Ifges(3,4) * t934;
t951 = Ifges(3,4) * t933 + Ifges(3,2) * t934;
t950 = Ifges(3,5) * t933 + Ifges(3,6) * t934;
t792 = mrSges(5,1) * t946 + mrSges(5,3) * t832 - mrSges(6,1) * t830 + mrSges(6,2) * t828 + mrSges(7,1) * t826 - mrSges(7,3) * t824 + pkin(5) * t821 - qJ(6) * t964 - pkin(4) * t819 + (-qJ(6) * t879 + t961) * t910 + t963 * t904 + (-qJ(6) * mrSges(7,2) + t966) * t895 + t968 * t865 + t991 * t864;
t799 = -mrSges(5,2) * t946 + mrSges(6,2) * t829 + mrSges(7,2) * t826 - mrSges(5,3) * t831 - mrSges(6,3) * t830 - mrSges(7,3) * t822 - qJ(5) * t819 - qJ(6) * t820 - t968 * t864 + t865 * t992 + t967 * t895 + t963 * t903 - t962 * t910;
t886 = Ifges(4,5) * t913 + Ifges(4,6) * t912 + Ifges(4,3) * qJD(3);
t887 = Ifges(4,4) * t913 + Ifges(4,2) * t912 + Ifges(4,6) * qJD(3);
t785 = mrSges(4,2) * t897 - mrSges(4,3) * t839 + Ifges(4,1) * t899 + Ifges(4,4) * t898 + Ifges(4,5) * qJDD(3) - pkin(8) * t807 - qJD(3) * t887 - t935 * t792 + t983 * t799 + t912 * t886;
t888 = Ifges(4,1) * t913 + Ifges(4,4) * t912 + Ifges(4,5) * qJD(3);
t791 = -mrSges(4,1) * t897 + mrSges(4,3) * t840 + Ifges(4,4) * t899 + Ifges(4,2) * t898 + Ifges(4,6) * qJDD(3) - pkin(3) * t807 + qJD(3) * t888 - t913 * t886 - t987;
t916 = t950 * qJD(1);
t781 = -mrSges(3,1) * t911 + mrSges(3,3) * t901 - pkin(2) * t944 + pkin(7) * t956 + t951 * qJDD(1) + t936 * t785 + t984 * t791 - t916 * t972;
t784 = t934 * qJD(1) * t916 + mrSges(3,2) * t911 - mrSges(3,3) * t900 - pkin(7) * t798 + t952 * qJDD(1) + t984 * t785 - t936 * t791;
t803 = mrSges(3,2) * t970 - t943;
t945 = mrSges(2,1) * t918 - mrSges(2,2) * t919 + Ifges(2,3) * qJDD(1) - pkin(1) * t803 + qJ(2) * t957 + t934 * t781 + t933 * t784;
t941 = mrSges(4,1) * t839 - mrSges(4,2) * t840 + Ifges(4,5) * t899 + Ifges(4,6) * t898 + Ifges(4,3) * qJDD(3) + pkin(3) * t816 + pkin(8) * t808 + t983 * t792 + t935 * t799 + t913 * t887 - t912 * t888;
t782 = -t941 + (Ifges(2,6) - t950) * qJDD(1) + mrSges(2,3) * t919 - mrSges(3,1) * t900 + mrSges(3,2) * t901 + mrSges(2,1) * g(3) - pkin(2) * t798 - pkin(1) * t790 + (-t933 * t951 + t934 * t952 + Ifges(2,5)) * t940;
t779 = -mrSges(2,2) * g(3) - mrSges(2,3) * t918 + Ifges(2,5) * qJDD(1) - t940 * Ifges(2,6) - qJ(2) * t790 - t933 * t781 + t934 * t784;
t1 = [-m(1) * g(1) + t958; -m(1) * g(2) + t976; (-m(1) - m(2)) * g(3) + t790; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t976 + t938 * t779 - t937 * t782; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t958 + t937 * t779 + t938 * t782; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t945; t945; t803; t941; t987; t817; t821;];
tauJB  = t1;
