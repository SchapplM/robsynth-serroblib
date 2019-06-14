% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:22:29
% EndTime: 2019-05-05 17:22:33
% DurationCPUTime: 2.88s
% Computational Cost: add. (16136->327), mult. (32147->356), div. (0->0), fcn. (14259->6), ass. (0->144)
t999 = Ifges(4,1) + Ifges(5,1);
t998 = Ifges(5,4) + Ifges(4,5);
t997 = -Ifges(6,1) - Ifges(5,3);
t996 = Ifges(6,4) - Ifges(5,5);
t995 = -Ifges(5,6) + Ifges(6,5);
t925 = sin(qJ(3));
t928 = cos(qJ(3));
t994 = t998 * qJD(3) + (t999 * t928 + (-Ifges(4,4) + Ifges(5,5)) * t925) * qJD(1);
t993 = Ifges(6,6) + t998;
t992 = Ifges(4,6) + t995;
t991 = -Ifges(6,3) - Ifges(4,3) - Ifges(5,2);
t957 = qJD(2) * qJD(1);
t903 = 0.2e1 * t957;
t958 = qJD(1) * qJD(3);
t882 = qJDD(1) * t925 + t928 * t958;
t901 = t925 * t958;
t883 = qJDD(1) * t928 - t901;
t931 = qJD(1) ^ 2;
t926 = sin(qJ(1));
t929 = cos(qJ(1));
t896 = -t929 * g(1) - t926 * g(2);
t952 = t931 * pkin(1) - qJDD(1) * qJ(2) - t896;
t949 = -t931 * pkin(7) - t952;
t946 = t882 * pkin(3) - t883 * qJ(4) + t949;
t950 = pkin(3) * t928 + qJ(4) * t925;
t960 = qJD(4) * t928;
t816 = t903 + (qJD(3) * t950 - 0.2e1 * t960) * qJD(1) + t946;
t961 = qJD(1) * t928;
t891 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t961;
t893 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t961;
t962 = qJD(1) * t925;
t894 = -mrSges(5,2) * t962 + qJD(3) * mrSges(5,3);
t888 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t962;
t969 = t888 * t925;
t990 = m(5) * t816 + t882 * mrSges(5,1) - t883 * mrSges(5,3) - qJD(1) * ((t891 + t893) * t928 + t969) + t894 * t962;
t895 = t926 * g(1) - t929 * g(2);
t945 = -t931 * qJ(2) + qJDD(2) - t895;
t840 = (-pkin(1) - pkin(7)) * qJDD(1) + t945;
t834 = -t928 * g(3) + t925 * t840;
t989 = -qJDD(3) * qJ(4) - t834;
t878 = (pkin(3) * t925 - qJ(4) * t928) * qJD(1);
t988 = t878 * t961 + qJDD(4);
t967 = t928 * t840;
t979 = t925 * g(3);
t833 = t967 + t979;
t930 = qJD(3) ^ 2;
t966 = t930 * qJ(4);
t820 = -qJDD(3) * pkin(3) - t833 - t966 + t988;
t987 = m(5) * t820 - qJDD(3) * mrSges(5,1) - qJD(3) * t894;
t839 = t903 + t949;
t889 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t962;
t892 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t961;
t984 = m(4) * t839 + t883 * mrSges(4,2) + t889 * t962 + t892 * t961;
t983 = -2 * qJD(5);
t982 = -pkin(3) - pkin(8);
t981 = -pkin(4) - pkin(8);
t978 = t930 * pkin(3);
t977 = mrSges(2,1) - mrSges(3,2);
t976 = mrSges(4,1) - mrSges(6,2);
t975 = -mrSges(4,3) - mrSges(5,2);
t974 = Ifges(2,5) - Ifges(3,4);
t973 = -Ifges(2,6) + Ifges(3,5);
t972 = -pkin(5) - qJ(4);
t971 = t882 * mrSges(6,2);
t968 = t925 ^ 2 * t931;
t890 = -qJD(3) * pkin(4) - qJ(5) * t961;
t904 = -0.2e1 * t957;
t937 = -qJ(5) * t968 + 0.2e1 * qJD(1) * t960 + t890 * t961 + qJDD(5) + t904 - t946;
t807 = t883 * pkin(5) + t937 + t981 * t882 + (t925 * t972 + t928 * t982) * t958;
t881 = (pkin(5) * t928 - pkin(8) * t925) * qJD(1);
t939 = t928 * t931 * t925 * pkin(4) + t961 * t983 + (-t883 - t901) * qJ(5) - t979 + t988;
t810 = t972 * t930 + (-qJD(1) * t881 - t840) * t928 + (-pkin(3) + t981) * qJDD(3) + t939;
t924 = sin(qJ(6));
t927 = cos(qJ(6));
t804 = t807 * t927 - t810 * t924;
t875 = -qJD(3) * t927 - t924 * t962;
t830 = qJD(6) * t875 - qJDD(3) * t924 + t882 * t927;
t876 = -qJD(3) * t924 + t927 * t962;
t832 = -mrSges(7,1) * t875 + mrSges(7,2) * t876;
t898 = qJD(6) + t961;
t835 = -mrSges(7,2) * t898 + mrSges(7,3) * t875;
t873 = qJDD(6) + t883;
t801 = m(7) * t804 + mrSges(7,1) * t873 - mrSges(7,3) * t830 - t832 * t876 + t835 * t898;
t805 = t807 * t924 + t810 * t927;
t829 = -qJD(6) * t876 - qJDD(3) * t927 - t882 * t924;
t836 = mrSges(7,1) * t898 - mrSges(7,3) * t876;
t802 = m(7) * t805 - mrSges(7,2) * t873 + mrSges(7,3) * t829 + t832 * t875 - t836 * t898;
t791 = -t924 * t801 + t927 * t802;
t814 = -t966 - t967 + (-pkin(3) - pkin(4)) * qJDD(3) + t939;
t877 = (mrSges(6,1) * t928 + mrSges(6,2) * t925) * qJD(1);
t948 = -m(6) * t814 + t877 * t961 - t791;
t879 = (mrSges(5,1) * t925 - mrSges(5,3) * t928) * qJD(1);
t953 = qJD(1) * (-t879 - (mrSges(4,1) * t925 + mrSges(4,2) * t928) * qJD(1));
t784 = m(4) * t833 + t976 * qJDD(3) + (-t888 + t889) * qJD(3) + t928 * t953 + (mrSges(6,3) + t975) * t883 + t948 - t987;
t902 = 0.2e1 * qJD(4) * qJD(3);
t819 = -t878 * t962 + t902 - t978 - t989;
t941 = pkin(4) * t968 - t882 * qJ(5) + t989;
t959 = t983 + t878;
t809 = qJDD(3) * pkin(5) + qJD(3) * t890 + t902 + t982 * t930 + (t881 - t959) * t962 - t941;
t806 = -m(7) * t809 + t829 * mrSges(7,1) - t830 * mrSges(7,2) + t875 * t835 - t876 * t836;
t813 = t978 + (-0.2e1 * qJD(4) - t890) * qJD(3) + t959 * t962 + t941;
t938 = -m(6) * t813 + qJDD(3) * mrSges(6,1) + t882 * mrSges(6,3) + qJD(3) * t891 + t877 * t962 - t806;
t935 = m(5) * t819 + qJDD(3) * mrSges(5,3) + qJD(3) * t893 + t938;
t793 = m(4) * t834 - qJDD(3) * mrSges(4,2) - qJD(3) * t892 + t882 * t975 + t925 * t953 + t935;
t779 = t928 * t784 + t925 * t793;
t846 = -qJDD(1) * pkin(1) + t945;
t943 = -m(3) * t846 + t931 * mrSges(3,3) - t779;
t775 = m(2) * t895 - t931 * mrSges(2,2) + qJDD(1) * t977 + t943;
t844 = t904 + t952;
t790 = t927 * t801 + t924 * t802;
t811 = -t882 * pkin(4) - t950 * t958 + t937;
t947 = m(6) * t811 + t883 * mrSges(6,1) + t790;
t934 = -t947 + t990;
t933 = -m(3) * t844 + t931 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t882 * t976 + t934 + t984;
t782 = m(2) * t896 - t931 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t933;
t965 = t929 * t775 + t926 * t782;
t964 = t995 * qJD(3) + (t925 * t997 + t996 * t928) * qJD(1);
t956 = -Ifges(4,4) - t996;
t955 = -t775 * t926 + t929 * t782;
t954 = -t925 * t784 + t928 * t793;
t951 = qJD(1) * (t991 * qJD(3) + (t925 * t992 - t928 * t993) * qJD(1));
t944 = t947 + t971;
t823 = Ifges(7,4) * t876 + Ifges(7,2) * t875 + Ifges(7,6) * t898;
t824 = Ifges(7,1) * t876 + Ifges(7,4) * t875 + Ifges(7,5) * t898;
t942 = mrSges(7,1) * t804 - mrSges(7,2) * t805 + Ifges(7,5) * t830 + Ifges(7,6) * t829 + Ifges(7,3) * t873 + t876 * t823 - t875 * t824;
t940 = qJDD(3) * mrSges(6,2) + qJD(3) * t888 - t948;
t785 = t934 - t971;
t788 = (t891 * t928 + t969) * qJD(1) + t944;
t822 = Ifges(7,5) * t876 + Ifges(7,6) * t875 + Ifges(7,3) * t898;
t794 = -mrSges(7,1) * t809 + mrSges(7,3) * t805 + Ifges(7,4) * t830 + Ifges(7,2) * t829 + Ifges(7,6) * t873 - t822 * t876 + t824 * t898;
t795 = mrSges(7,2) * t809 - mrSges(7,3) * t804 + Ifges(7,1) * t830 + Ifges(7,4) * t829 + Ifges(7,5) * t873 + t822 * t875 - t823 * t898;
t853 = -Ifges(6,6) * qJD(3) + (Ifges(6,4) * t925 - Ifges(6,2) * t928) * qJD(1);
t771 = -t927 * t795 - qJ(5) * t938 + t924 * t794 + pkin(4) * t788 + pkin(8) * t790 - pkin(3) * t785 + mrSges(4,3) * t834 - mrSges(4,1) * t839 - mrSges(6,2) * t811 + mrSges(6,3) * t813 - mrSges(5,1) * t816 + mrSges(5,2) * t819 - t956 * t883 + (-Ifges(4,2) + t997) * t882 + t992 * qJDD(3) + (-t853 + t994) * qJD(3) + t928 * t951;
t789 = -t883 * mrSges(6,3) + t940;
t855 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t928 - Ifges(4,2) * t925) * qJD(1);
t773 = t942 + t993 * qJDD(3) - qJ(5) * t789 + pkin(5) * t790 - qJ(4) * t785 + t956 * t882 + (Ifges(6,2) + t999) * t883 + t925 * t951 + (-t855 - t964) * qJD(3) - mrSges(4,3) * t833 + mrSges(4,2) * t839 + mrSges(6,1) * t811 - mrSges(6,3) * t814 - mrSges(5,3) * t816 + mrSges(5,2) * t820;
t777 = qJDD(1) * mrSges(3,2) - t943;
t936 = mrSges(2,1) * t895 - mrSges(2,2) * t896 + mrSges(3,2) * t846 - mrSges(3,3) * t844 - pkin(1) * t777 - pkin(7) * t779 + qJ(2) * t933 - t771 * t925 + t928 * t773 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t787 = t879 * t961 + (mrSges(5,2) - mrSges(6,3)) * t883 + t940 + t987;
t932 = t855 * t961 + (-t925 * t853 + t928 * t964) * qJD(1) + qJ(4) * t935 - t927 * t794 - t924 * t795 - pkin(8) * t791 - pkin(4) * t789 - pkin(3) * t787 + mrSges(4,1) * t833 - mrSges(4,2) * t834 - mrSges(6,1) * t813 + mrSges(6,2) * t814 + mrSges(5,3) * t819 - mrSges(5,1) * t820 - pkin(5) * t806 + (-qJ(4) * t879 + t994) * t962 + t993 * t883 + (-mrSges(5,2) * qJ(4) - t992) * t882 - t991 * qJDD(3);
t778 = -m(3) * g(3) + t954;
t770 = t932 - mrSges(2,3) * t895 + pkin(2) * t779 - qJ(2) * t778 + t973 * t931 + mrSges(3,1) * t846 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t974 * qJDD(1);
t769 = -mrSges(3,1) * t844 + mrSges(2,3) * t896 - pkin(1) * t778 - pkin(7) * t954 + t977 * g(3) - t973 * qJDD(1) - t928 * t771 - t925 * t773 + t974 * t931 + (t882 * mrSges(4,1) - t944 + t984 + t990) * pkin(2);
t1 = [-m(1) * g(1) + t955; -m(1) * g(2) + t965; (-m(1) - m(2) - m(3)) * g(3) + t954; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t965 - t926 * t769 + t929 * t770; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t955 + t929 * t769 + t926 * t770; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t936; t936; t777; t932; t787; t788; t942;];
tauJB  = t1;
