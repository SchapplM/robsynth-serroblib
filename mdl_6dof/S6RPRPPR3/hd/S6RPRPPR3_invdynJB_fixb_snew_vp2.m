% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-05-05 16:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:41:47
% EndTime: 2019-05-05 16:41:51
% DurationCPUTime: 3.09s
% Computational Cost: add. (27957->306), mult. (55420->357), div. (0->0), fcn. (26369->8), ass. (0->128)
t982 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t958 = Ifges(4,4) - Ifges(5,5) + Ifges(6,4);
t957 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t981 = Ifges(4,2) + Ifges(5,3) + Ifges(6,1);
t956 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t980 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t928 = sin(qJ(3));
t931 = cos(qJ(3));
t885 = (-mrSges(5,1) * t931 - mrSges(5,3) * t928) * qJD(1);
t960 = qJD(1) * qJD(3);
t951 = t931 * t960;
t890 = qJDD(1) * t928 + t951;
t961 = qJD(1) * t931;
t900 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t961;
t950 = t928 * t960;
t891 = qJDD(1) * t931 - t950;
t962 = qJD(1) * t928;
t896 = -qJD(3) * pkin(4) - qJ(5) * t962;
t929 = sin(qJ(1));
t932 = cos(qJ(1));
t903 = t929 * g(1) - g(2) * t932;
t883 = qJDD(1) * pkin(1) + t903;
t904 = -g(1) * t932 - g(2) * t929;
t934 = qJD(1) ^ 2;
t888 = -pkin(1) * t934 + t904;
t924 = sin(pkin(9));
t925 = cos(pkin(9));
t844 = t925 * t883 - t924 * t888;
t832 = -qJDD(1) * pkin(2) - t934 * pkin(7) - t844;
t944 = -t891 * pkin(3) + t832 + (-t890 - t951) * qJ(4);
t965 = t931 ^ 2 * t934;
t972 = 2 * qJD(4);
t938 = -qJ(5) * t965 + qJDD(5) - t944 + (t896 + t972) * t962;
t970 = pkin(4) + pkin(8);
t971 = -pkin(3) - pkin(8);
t815 = t938 + t970 * t891 + t890 * pkin(5) + (pkin(5) * t931 + t971 * t928) * t960;
t889 = (pkin(5) * t928 + pkin(8) * t931) * qJD(1);
t933 = qJD(3) ^ 2;
t845 = t924 * t883 + t925 * t888;
t833 = -pkin(2) * t934 + qJDD(1) * pkin(7) + t845;
t922 = -g(3) + qJDD(2);
t828 = -t928 * t833 + t931 * t922;
t884 = (-pkin(3) * t931 - qJ(4) * t928) * qJD(1);
t946 = t884 * t962 + qJDD(4) - t828;
t964 = t931 * t934;
t959 = qJD(1) * qJD(5);
t977 = -0.2e1 * t928 * t959 + (-t890 + t951) * qJ(5);
t818 = (-pkin(5) - qJ(4)) * t933 + (-pkin(4) * t964 - qJD(1) * t889) * t928 + (-pkin(3) - t970) * qJDD(3) + t946 + t977;
t927 = sin(qJ(6));
t930 = cos(qJ(6));
t812 = t815 * t930 - t818 * t927;
t881 = -qJD(3) * t930 + t927 * t961;
t842 = qJD(6) * t881 - qJDD(3) * t927 - t891 * t930;
t882 = -qJD(3) * t927 - t930 * t961;
t846 = -mrSges(7,1) * t881 + mrSges(7,2) * t882;
t907 = qJD(6) + t962;
t847 = -mrSges(7,2) * t907 + mrSges(7,3) * t881;
t879 = qJDD(6) + t890;
t809 = m(7) * t812 + mrSges(7,1) * t879 - mrSges(7,3) * t842 - t846 * t882 + t847 * t907;
t813 = t815 * t927 + t818 * t930;
t841 = -qJD(6) * t882 - qJDD(3) * t930 + t891 * t927;
t848 = mrSges(7,1) * t907 - mrSges(7,3) * t882;
t810 = m(7) * t813 - mrSges(7,2) * t879 + mrSges(7,3) * t841 + t846 * t881 - t848 * t907;
t799 = -t927 * t809 + t930 * t810;
t826 = -qJDD(3) * pkin(3) - t933 * qJ(4) + t946;
t822 = (-t928 * t964 - qJDD(3)) * pkin(4) + t826 + t977;
t887 = (mrSges(6,1) * t928 - mrSges(6,2) * t931) * qJD(1);
t945 = -m(6) * t822 + t887 * t962 - t799;
t941 = qJDD(3) * mrSges(6,2) + qJD(3) * t900 - t945;
t902 = mrSges(5,2) * t961 + qJD(3) * mrSges(5,3);
t976 = m(5) * t826 - qJDD(3) * mrSges(5,1) - qJD(3) * t902;
t795 = t885 * t962 + (mrSges(5,2) - mrSges(6,3)) * t890 + t941 + t976;
t797 = -t890 * mrSges(6,3) + t941;
t913 = qJD(3) * t972;
t829 = t931 * t833 + t928 * t922;
t978 = qJDD(3) * qJ(4) + t884 * t961 + t829;
t943 = pkin(4) * t965 + t891 * qJ(5) - t978;
t817 = qJDD(3) * pkin(5) + qJD(3) * t896 + t913 + t971 * t933 + (-0.2e1 * qJD(5) - t889) * t961 - t943;
t834 = Ifges(7,5) * t882 + Ifges(7,6) * t881 + Ifges(7,3) * t907;
t836 = Ifges(7,1) * t882 + Ifges(7,4) * t881 + Ifges(7,5) * t907;
t803 = -mrSges(7,1) * t817 + mrSges(7,3) * t813 + Ifges(7,4) * t842 + Ifges(7,2) * t841 + Ifges(7,6) * t879 - t834 * t882 + t836 * t907;
t835 = Ifges(7,4) * t882 + Ifges(7,2) * t881 + Ifges(7,6) * t907;
t804 = mrSges(7,2) * t817 - mrSges(7,3) * t812 + Ifges(7,1) * t842 + Ifges(7,4) * t841 + Ifges(7,5) * t879 + t834 * t881 - t835 * t907;
t814 = -m(7) * t817 + t841 * mrSges(7,1) - t842 * mrSges(7,2) + t881 * t847 - t882 * t848;
t968 = t933 * pkin(3);
t973 = -2 * qJD(4);
t821 = 0.2e1 * t931 * t959 + t968 + (t973 - t896) * qJD(3) + t943;
t825 = t913 - t968 + t978;
t899 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t962;
t897 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t962;
t940 = -m(6) * t821 + qJDD(3) * mrSges(6,1) - t891 * mrSges(6,3) + qJD(3) * t897 - t814;
t937 = m(5) * t825 + qJDD(3) * mrSges(5,3) + qJD(3) * t899 + t885 * t961 + t940;
t952 = -t957 * qJD(3) + (-t982 * t928 - t958 * t931) * qJD(1);
t953 = -t956 * qJD(3) + (-t958 * t928 - t981 * t931) * qJD(1);
t979 = -(t953 * t928 - t952 * t931) * qJD(1) + t980 * qJDD(3) + t957 * t890 + t956 * t891 + mrSges(4,1) * t828 - mrSges(5,1) * t826 - mrSges(6,1) * t821 - mrSges(4,2) * t829 + mrSges(6,2) * t822 + mrSges(5,3) * t825 - pkin(3) * t795 - pkin(4) * t797 - pkin(5) * t814 - pkin(8) * t799 + qJ(4) * (t891 * mrSges(5,2) - t887 * t961 + t937) - t930 * t803 - t927 * t804;
t967 = mrSges(4,3) + mrSges(5,2);
t886 = (-mrSges(4,1) * t931 + mrSges(4,2) * t928) * qJD(1);
t901 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t961;
t793 = m(4) * t828 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t900 + t901) * qJD(3) + (-t885 - t886) * t962 + (mrSges(6,3) - t967) * t890 + t945 - t976;
t898 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t962;
t802 = t937 + t967 * t891 - qJDD(3) * mrSges(4,2) + (t886 - t887) * t961 - qJD(3) * t898 + m(4) * t829;
t947 = -t793 * t928 + t931 * t802;
t785 = m(3) * t845 - mrSges(3,1) * t934 - qJDD(1) * mrSges(3,2) + t947;
t798 = t930 * t809 + t927 * t810;
t820 = -pkin(3) * t950 + t891 * pkin(4) + t938;
t796 = m(6) * t820 + t890 * mrSges(6,1) - t891 * mrSges(6,2) + t897 * t962 - t900 * t961 + t798;
t823 = (pkin(3) * qJD(3) + t973) * t962 + t944;
t794 = m(5) * t823 - mrSges(5,1) * t891 - t890 * mrSges(5,3) - t899 * t962 - t902 * t961 - t796;
t936 = -m(4) * t832 + t891 * mrSges(4,1) - mrSges(4,2) * t890 - t898 * t962 + t901 * t961 - t794;
t790 = m(3) * t844 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t934 + t936;
t782 = t924 * t785 + t925 * t790;
t779 = m(2) * t903 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t934 + t782;
t948 = t925 * t785 - t790 * t924;
t780 = m(2) * t904 - mrSges(2,1) * t934 - qJDD(1) * mrSges(2,2) + t948;
t963 = t932 * t779 + t929 * t780;
t788 = t931 * t793 + t928 * t802;
t954 = -t980 * qJD(3) + (-t957 * t928 - t956 * t931) * qJD(1);
t786 = m(3) * t922 + t788;
t949 = -t779 * t929 + t932 * t780;
t942 = mrSges(7,1) * t812 - mrSges(7,2) * t813 + Ifges(7,5) * t842 + Ifges(7,6) * t841 + Ifges(7,3) * t879 + t882 * t835 - t881 * t836;
t773 = -t930 * t804 - qJ(5) * t940 + t927 * t803 + mrSges(4,3) * t829 - mrSges(4,1) * t832 - mrSges(6,2) * t820 + mrSges(6,3) * t821 - mrSges(5,1) * t823 + mrSges(5,2) * t825 + pkin(8) * t798 + pkin(4) * t796 - pkin(3) * t794 + t981 * t891 + t958 * t890 + t956 * qJDD(3) - t952 * qJD(3) + (qJ(5) * t887 * t931 + t954 * t928) * qJD(1);
t775 = mrSges(6,1) * t820 + mrSges(4,2) * t832 + mrSges(5,2) * t826 - mrSges(4,3) * t828 - mrSges(5,3) * t823 - mrSges(6,3) * t822 + pkin(5) * t798 - qJ(4) * t794 - qJ(5) * t797 + t953 * qJD(3) + t957 * qJDD(3) + t982 * t890 + t958 * t891 - t954 * t961 + t942;
t939 = mrSges(2,1) * t903 + mrSges(3,1) * t844 - mrSges(2,2) * t904 - mrSges(3,2) * t845 + pkin(1) * t782 + pkin(2) * t936 + pkin(7) * t947 + t931 * t773 + t928 * t775 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t771 = -mrSges(3,1) * t922 + mrSges(3,3) * t845 + t934 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t788 - t979;
t770 = mrSges(3,2) * t922 - mrSges(3,3) * t844 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t934 - pkin(7) * t788 - t773 * t928 + t775 * t931;
t769 = -mrSges(2,2) * g(3) - mrSges(2,3) * t903 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t934 - qJ(2) * t782 + t770 * t925 - t771 * t924;
t768 = mrSges(2,1) * g(3) + mrSges(2,3) * t904 + t934 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t786 + qJ(2) * t948 + t924 * t770 + t925 * t771;
t1 = [-m(1) * g(1) + t949; -m(1) * g(2) + t963; (-m(1) - m(2)) * g(3) + t786; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t963 - t929 * t768 + t932 * t769; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t949 + t932 * t768 + t929 * t769; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t939; t939; t786; t979; t795; t796; t942;];
tauJB  = t1;
