% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:38:15
% EndTime: 2019-05-07 07:38:32
% DurationCPUTime: 15.64s
% Computational Cost: add. (257789->365), mult. (532731->447), div. (0->0), fcn. (381862->10), ass. (0->144)
t974 = Ifges(6,1) + Ifges(7,1);
t965 = Ifges(6,4) - Ifges(7,5);
t964 = -Ifges(6,5) - Ifges(7,4);
t973 = Ifges(6,2) + Ifges(7,3);
t963 = Ifges(6,6) - Ifges(7,6);
t972 = -Ifges(6,3) - Ifges(7,2);
t931 = sin(qJ(2));
t933 = cos(qJ(2));
t953 = qJD(1) * qJD(2);
t911 = t931 * qJDD(1) + t933 * t953;
t932 = sin(qJ(1));
t934 = cos(qJ(1));
t918 = -t934 * g(1) - t932 * g(2);
t935 = qJD(1) ^ 2;
t905 = -t935 * pkin(1) + qJDD(1) * pkin(7) + t918;
t961 = t931 * t905;
t967 = pkin(2) * t935;
t866 = qJDD(2) * pkin(2) - t911 * pkin(8) - t961 + (pkin(8) * t953 + t931 * t967 - g(3)) * t933;
t893 = -t931 * g(3) + t933 * t905;
t912 = t933 * qJDD(1) - t931 * t953;
t955 = qJD(1) * t931;
t916 = qJD(2) * pkin(2) - pkin(8) * t955;
t926 = t933 ^ 2;
t867 = t912 * pkin(8) - qJD(2) * t916 - t926 * t967 + t893;
t930 = sin(qJ(3));
t969 = cos(qJ(3));
t843 = t930 * t866 + t867 * t969;
t903 = (t930 * t933 + t931 * t969) * qJD(1);
t875 = t903 * qJD(3) + t930 * t911 - t912 * t969;
t954 = qJD(1) * t933;
t902 = t930 * t955 - t954 * t969;
t885 = t902 * mrSges(4,1) + t903 * mrSges(4,2);
t924 = qJD(2) + qJD(3);
t895 = t924 * mrSges(4,1) - t903 * mrSges(4,3);
t923 = qJDD(2) + qJDD(3);
t876 = -t902 * qJD(3) + t911 * t969 + t930 * t912;
t917 = t932 * g(1) - t934 * g(2);
t943 = -qJDD(1) * pkin(1) - t917;
t877 = -t912 * pkin(2) + t916 * t955 + (-pkin(8) * t926 - pkin(7)) * t935 + t943;
t826 = (t902 * t924 - t876) * qJ(4) + (t903 * t924 + t875) * pkin(3) + t877;
t884 = t902 * pkin(3) - t903 * qJ(4);
t922 = t924 ^ 2;
t829 = -t922 * pkin(3) + t923 * qJ(4) - t902 * t884 + t843;
t927 = sin(pkin(10));
t928 = cos(pkin(10));
t891 = t928 * t903 + t927 * t924;
t817 = -0.2e1 * qJD(4) * t891 + t928 * t826 - t927 * t829;
t860 = t928 * t876 + t927 * t923;
t890 = -t927 * t903 + t928 * t924;
t814 = (t902 * t890 - t860) * pkin(9) + (t890 * t891 + t875) * pkin(4) + t817;
t818 = 0.2e1 * qJD(4) * t890 + t927 * t826 + t928 * t829;
t859 = -t927 * t876 + t928 * t923;
t879 = t902 * pkin(4) - t891 * pkin(9);
t889 = t890 ^ 2;
t816 = -t889 * pkin(4) + t859 * pkin(9) - t902 * t879 + t818;
t929 = sin(qJ(5));
t968 = cos(qJ(5));
t811 = t929 * t814 + t968 * t816;
t858 = t929 * t890 + t891 * t968;
t824 = t858 * qJD(5) - t859 * t968 + t929 * t860;
t898 = qJD(5) + t902;
t848 = t898 * mrSges(6,1) - t858 * mrSges(6,3);
t857 = -t890 * t968 + t929 * t891;
t874 = qJDD(5) + t875;
t838 = t857 * pkin(5) - t858 * qJ(6);
t897 = t898 ^ 2;
t808 = -t897 * pkin(5) + t874 * qJ(6) + 0.2e1 * qJD(6) * t898 - t857 * t838 + t811;
t849 = -t898 * mrSges(7,1) + t858 * mrSges(7,2);
t952 = m(7) * t808 + t874 * mrSges(7,3) + t898 * t849;
t839 = t857 * mrSges(7,1) - t858 * mrSges(7,3);
t956 = -t857 * mrSges(6,1) - t858 * mrSges(6,2) - t839;
t966 = -mrSges(6,3) - mrSges(7,2);
t797 = m(6) * t811 - t874 * mrSges(6,2) + t824 * t966 - t898 * t848 + t857 * t956 + t952;
t810 = t814 * t968 - t929 * t816;
t825 = -t857 * qJD(5) + t929 * t859 + t860 * t968;
t847 = -t898 * mrSges(6,2) - t857 * mrSges(6,3);
t809 = -t874 * pkin(5) - t897 * qJ(6) + t858 * t838 + qJDD(6) - t810;
t846 = -t857 * mrSges(7,2) + t898 * mrSges(7,3);
t946 = -m(7) * t809 + t874 * mrSges(7,1) + t898 * t846;
t799 = m(6) * t810 + t874 * mrSges(6,1) + t825 * t966 + t898 * t847 + t858 * t956 + t946;
t792 = t929 * t797 + t968 * t799;
t862 = -t890 * mrSges(5,1) + t891 * mrSges(5,2);
t945 = -t902 * mrSges(5,2) + t890 * mrSges(5,3);
t790 = m(5) * t817 + t875 * mrSges(5,1) - t860 * mrSges(5,3) - t891 * t862 + t902 * t945 + t792;
t878 = t902 * mrSges(5,1) - t891 * mrSges(5,3);
t947 = t968 * t797 - t929 * t799;
t791 = m(5) * t818 - t875 * mrSges(5,2) + t859 * mrSges(5,3) + t890 * t862 - t902 * t878 + t947;
t948 = -t927 * t790 + t928 * t791;
t783 = m(4) * t843 - t923 * mrSges(4,2) - t875 * mrSges(4,3) - t902 * t885 - t924 * t895 + t948;
t842 = t866 * t969 - t930 * t867;
t828 = -t923 * pkin(3) - t922 * qJ(4) + t903 * t884 + qJDD(4) - t842;
t819 = -t859 * pkin(4) - t889 * pkin(9) + t891 * t879 + t828;
t812 = -0.2e1 * qJD(6) * t858 + (t857 * t898 - t825) * qJ(6) + (t858 * t898 + t824) * pkin(5) + t819;
t806 = m(7) * t812 + t824 * mrSges(7,1) - t825 * mrSges(7,3) + t857 * t846 - t858 * t849;
t938 = m(6) * t819 + t824 * mrSges(6,1) + t825 * mrSges(6,2) + t857 * t847 + t858 * t848 + t806;
t803 = m(5) * t828 - t859 * mrSges(5,1) + t860 * mrSges(5,2) + t891 * t878 - t890 * t945 + t938;
t894 = -t924 * mrSges(4,2) - t902 * mrSges(4,3);
t801 = m(4) * t842 + t923 * mrSges(4,1) - t876 * mrSges(4,3) - t903 * t885 + t924 * t894 - t803;
t777 = t930 * t783 + t969 * t801;
t892 = -t933 * g(3) - t961;
t900 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t931 + Ifges(3,2) * t933) * qJD(1);
t901 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t931 + Ifges(3,4) * t933) * qJD(1);
t957 = -t965 * t857 + t974 * t858 - t964 * t898;
t958 = t963 * t857 + t964 * t858 + t972 * t898;
t793 = -mrSges(6,1) * t819 - mrSges(7,1) * t812 + mrSges(7,2) * t808 + mrSges(6,3) * t811 - pkin(5) * t806 - t973 * t824 + t965 * t825 + t958 * t858 + t963 * t874 + t957 * t898;
t959 = t973 * t857 - t965 * t858 - t963 * t898;
t794 = mrSges(6,2) * t819 + mrSges(7,2) * t809 - mrSges(6,3) * t810 - mrSges(7,3) * t812 - qJ(6) * t806 - t965 * t824 + t974 * t825 + t958 * t857 - t964 * t874 + t959 * t898;
t850 = Ifges(5,5) * t891 + Ifges(5,6) * t890 + Ifges(5,3) * t902;
t852 = Ifges(5,1) * t891 + Ifges(5,4) * t890 + Ifges(5,5) * t902;
t769 = -mrSges(5,1) * t828 + mrSges(5,3) * t818 + Ifges(5,4) * t860 + Ifges(5,2) * t859 + Ifges(5,6) * t875 - pkin(4) * t938 + pkin(9) * t947 + t793 * t968 + t929 * t794 - t891 * t850 + t902 * t852;
t851 = Ifges(5,4) * t891 + Ifges(5,2) * t890 + Ifges(5,6) * t902;
t771 = mrSges(5,2) * t828 - mrSges(5,3) * t817 + Ifges(5,1) * t860 + Ifges(5,4) * t859 + Ifges(5,5) * t875 - pkin(9) * t792 - t929 * t793 + t794 * t968 + t890 * t850 - t902 * t851;
t881 = Ifges(4,4) * t903 - Ifges(4,2) * t902 + Ifges(4,6) * t924;
t882 = Ifges(4,1) * t903 - Ifges(4,4) * t902 + Ifges(4,5) * t924;
t940 = -mrSges(4,1) * t842 + mrSges(4,2) * t843 - Ifges(4,5) * t876 + Ifges(4,6) * t875 - Ifges(4,3) * t923 + pkin(3) * t803 - qJ(4) * t948 - t928 * t769 - t927 * t771 - t903 * t881 - t902 * t882;
t971 = mrSges(3,1) * t892 - mrSges(3,2) * t893 + Ifges(3,5) * t911 + Ifges(3,6) * t912 + Ifges(3,3) * qJDD(2) + pkin(2) * t777 + (t931 * t900 - t933 * t901) * qJD(1) - t940;
t805 = t825 * mrSges(7,2) + t858 * t839 - t946;
t970 = t963 * t824 + t964 * t825 - t957 * t857 + t959 * t858 + t972 * t874 - mrSges(6,1) * t810 + mrSges(7,1) * t809 + mrSges(6,2) * t811 - mrSges(7,3) * t808 + pkin(5) * t805 - qJ(6) * (-t824 * mrSges(7,2) - t857 * t839 + t952);
t910 = (-mrSges(3,1) * t933 + mrSges(3,2) * t931) * qJD(1);
t915 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t954;
t775 = m(3) * t892 + qJDD(2) * mrSges(3,1) - t911 * mrSges(3,3) + qJD(2) * t915 - t910 * t955 + t777;
t914 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t955;
t949 = t783 * t969 - t930 * t801;
t776 = m(3) * t893 - qJDD(2) * mrSges(3,2) + t912 * mrSges(3,3) - qJD(2) * t914 + t910 * t954 + t949;
t950 = -t931 * t775 + t933 * t776;
t764 = m(2) * t918 - t935 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t950;
t904 = -t935 * pkin(7) + t943;
t785 = t928 * t790 + t927 * t791;
t941 = m(4) * t877 + t875 * mrSges(4,1) + t876 * mrSges(4,2) + t902 * t894 + t903 * t895 + t785;
t937 = -m(3) * t904 + t912 * mrSges(3,1) - t911 * mrSges(3,2) - t914 * t955 + t915 * t954 - t941;
t779 = m(2) * t917 + qJDD(1) * mrSges(2,1) - t935 * mrSges(2,2) + t937;
t960 = t932 * t764 + t934 * t779;
t766 = t933 * t775 + t931 * t776;
t951 = t934 * t764 - t932 * t779;
t880 = Ifges(4,5) * t903 - Ifges(4,6) * t902 + Ifges(4,3) * t924;
t761 = mrSges(4,2) * t877 - mrSges(4,3) * t842 + Ifges(4,1) * t876 - Ifges(4,4) * t875 + Ifges(4,5) * t923 - qJ(4) * t785 - t927 * t769 + t928 * t771 - t902 * t880 - t924 * t881;
t767 = t970 + Ifges(4,6) * t923 + t924 * t882 - t903 * t880 - t891 * t851 + t890 * t852 + Ifges(4,4) * t876 - mrSges(4,1) * t877 - Ifges(5,6) * t859 - Ifges(5,5) * t860 + mrSges(4,3) * t843 - mrSges(5,1) * t817 + mrSges(5,2) * t818 - pkin(4) * t792 + (-Ifges(4,2) - Ifges(5,3)) * t875 - pkin(3) * t785;
t899 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t931 + Ifges(3,6) * t933) * qJD(1);
t757 = -mrSges(3,1) * t904 + mrSges(3,3) * t893 + Ifges(3,4) * t911 + Ifges(3,2) * t912 + Ifges(3,6) * qJDD(2) - pkin(2) * t941 + pkin(8) * t949 + qJD(2) * t901 + t930 * t761 + t767 * t969 - t899 * t955;
t760 = mrSges(3,2) * t904 - mrSges(3,3) * t892 + Ifges(3,1) * t911 + Ifges(3,4) * t912 + Ifges(3,5) * qJDD(2) - pkin(8) * t777 - qJD(2) * t900 + t761 * t969 - t930 * t767 + t899 * t954;
t942 = mrSges(2,1) * t917 - mrSges(2,2) * t918 + Ifges(2,3) * qJDD(1) + pkin(1) * t937 + pkin(7) * t950 + t933 * t757 + t931 * t760;
t758 = mrSges(2,1) * g(3) + mrSges(2,3) * t918 + t935 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t766 - t971;
t755 = -mrSges(2,2) * g(3) - mrSges(2,3) * t917 + Ifges(2,5) * qJDD(1) - t935 * Ifges(2,6) - pkin(7) * t766 - t931 * t757 + t933 * t760;
t1 = [-m(1) * g(1) + t951; -m(1) * g(2) + t960; (-m(1) - m(2)) * g(3) + t766; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t960 + t934 * t755 - t932 * t758; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t951 + t932 * t755 + t934 * t758; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t942; t942; t971; -t940; t803; -t970; t805;];
tauJB  = t1;
