% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 11:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:47:37
% EndTime: 2019-05-06 11:47:50
% DurationCPUTime: 12.53s
% Computational Cost: add. (195564->368), mult. (430576->447), div. (0->0), fcn. (274094->10), ass. (0->147)
t1002 = -2 * qJD(3);
t1001 = Ifges(3,1) + Ifges(4,2);
t994 = Ifges(3,4) + Ifges(4,6);
t993 = Ifges(3,5) - Ifges(4,4);
t1000 = Ifges(3,2) + Ifges(4,3);
t992 = Ifges(3,6) - Ifges(4,5);
t999 = Ifges(3,3) + Ifges(4,1);
t960 = cos(qJ(2));
t956 = sin(qJ(2));
t984 = qJD(1) * qJD(2);
t981 = t956 * t984;
t928 = qJDD(1) * t960 - t981;
t944 = t956 * qJD(1);
t934 = pkin(3) * t944 - qJD(2) * qJ(4);
t951 = t960 ^ 2;
t963 = qJD(1) ^ 2;
t982 = t960 * t984;
t927 = qJDD(1) * t956 + t982;
t957 = sin(qJ(1));
t961 = cos(qJ(1));
t937 = t957 * g(1) - t961 * g(2);
t976 = -qJDD(1) * pkin(1) - t937;
t969 = pkin(2) * t981 + t944 * t1002 + (-t927 - t982) * qJ(3) + t976;
t854 = -t934 * t944 + (-pkin(3) * t951 - pkin(7)) * t963 + (-pkin(2) - qJ(4)) * t928 + t969;
t938 = -g(1) * t961 - g(2) * t957;
t912 = -pkin(1) * t963 + qJDD(1) * pkin(7) + t938;
t889 = -t960 * g(3) - t956 * t912;
t924 = (-pkin(2) * t960 - qJ(3) * t956) * qJD(1);
t962 = qJD(2) ^ 2;
t876 = -qJDD(2) * pkin(2) - t962 * qJ(3) + t924 * t944 + qJDD(3) - t889;
t869 = (-t956 * t960 * t963 - qJDD(2)) * qJ(4) + (t927 - t982) * pkin(3) + t876;
t952 = sin(pkin(10));
t953 = cos(pkin(10));
t985 = qJD(1) * t960;
t919 = qJD(2) * t953 - t952 * t985;
t843 = -0.2e1 * qJD(4) * t919 - t952 * t854 + t953 * t869;
t895 = qJDD(2) * t953 - t928 * t952;
t918 = -qJD(2) * t952 - t953 * t985;
t834 = (t918 * t944 - t895) * pkin(8) + (t918 * t919 + t927) * pkin(4) + t843;
t844 = 0.2e1 * qJD(4) * t918 + t953 * t854 + t952 * t869;
t894 = -qJDD(2) * t952 - t928 * t953;
t896 = pkin(4) * t944 - pkin(8) * t919;
t917 = t918 ^ 2;
t836 = -pkin(4) * t917 + pkin(8) * t894 - t896 * t944 + t844;
t955 = sin(qJ(5));
t959 = cos(qJ(5));
t828 = t959 * t834 - t955 * t836;
t886 = t918 * t959 - t919 * t955;
t861 = qJD(5) * t886 + t894 * t955 + t895 * t959;
t887 = t918 * t955 + t919 * t959;
t923 = qJDD(5) + t927;
t941 = t944 + qJD(5);
t825 = (t886 * t941 - t861) * pkin(9) + (t886 * t887 + t923) * pkin(5) + t828;
t829 = t955 * t834 + t959 * t836;
t860 = -qJD(5) * t887 + t894 * t959 - t895 * t955;
t879 = pkin(5) * t941 - pkin(9) * t887;
t885 = t886 ^ 2;
t826 = -pkin(5) * t885 + pkin(9) * t860 - t879 * t941 + t829;
t954 = sin(qJ(6));
t958 = cos(qJ(6));
t824 = t825 * t954 + t826 * t958;
t890 = -g(3) * t956 + t960 * t912;
t875 = pkin(2) * t962 - qJDD(2) * qJ(3) + qJD(2) * t1002 - t924 * t985 - t890;
t865 = -qJ(4) * t951 * t963 + pkin(3) * t928 + qJD(2) * t934 + qJDD(4) - t875;
t850 = -pkin(4) * t894 - pkin(8) * t917 + t919 * t896 + t865;
t831 = -pkin(5) * t860 - pkin(9) * t885 + t879 * t887 + t850;
t872 = t886 * t954 + t887 * t958;
t840 = -qJD(6) * t872 + t860 * t958 - t861 * t954;
t871 = t886 * t958 - t887 * t954;
t841 = qJD(6) * t871 + t860 * t954 + t861 * t958;
t939 = qJD(6) + t941;
t845 = Ifges(7,5) * t872 + Ifges(7,6) * t871 + Ifges(7,3) * t939;
t847 = Ifges(7,1) * t872 + Ifges(7,4) * t871 + Ifges(7,5) * t939;
t914 = qJDD(6) + t923;
t810 = -mrSges(7,1) * t831 + mrSges(7,3) * t824 + Ifges(7,4) * t841 + Ifges(7,2) * t840 + Ifges(7,6) * t914 - t845 * t872 + t847 * t939;
t823 = t825 * t958 - t826 * t954;
t846 = Ifges(7,4) * t872 + Ifges(7,2) * t871 + Ifges(7,6) * t939;
t811 = mrSges(7,2) * t831 - mrSges(7,3) * t823 + Ifges(7,1) * t841 + Ifges(7,4) * t840 + Ifges(7,5) * t914 + t845 * t871 - t846 * t939;
t866 = Ifges(6,5) * t887 + Ifges(6,6) * t886 + Ifges(6,3) * t941;
t868 = Ifges(6,1) * t887 + Ifges(6,4) * t886 + Ifges(6,5) * t941;
t855 = -mrSges(7,2) * t939 + mrSges(7,3) * t871;
t856 = mrSges(7,1) * t939 - mrSges(7,3) * t872;
t973 = m(7) * t831 - t840 * mrSges(7,1) + t841 * mrSges(7,2) - t871 * t855 + t872 * t856;
t851 = -mrSges(7,1) * t871 + mrSges(7,2) * t872;
t815 = m(7) * t823 + mrSges(7,1) * t914 - mrSges(7,3) * t841 - t851 * t872 + t855 * t939;
t816 = m(7) * t824 - mrSges(7,2) * t914 + mrSges(7,3) * t840 + t851 * t871 - t856 * t939;
t977 = -t815 * t954 + t958 * t816;
t796 = -mrSges(6,1) * t850 + mrSges(6,3) * t829 + Ifges(6,4) * t861 + Ifges(6,2) * t860 + Ifges(6,6) * t923 - pkin(5) * t973 + pkin(9) * t977 + t958 * t810 + t954 * t811 - t887 * t866 + t941 * t868;
t809 = t958 * t815 + t954 * t816;
t867 = Ifges(6,4) * t887 + Ifges(6,2) * t886 + Ifges(6,6) * t941;
t797 = mrSges(6,2) * t850 - mrSges(6,3) * t828 + Ifges(6,1) * t861 + Ifges(6,4) * t860 + Ifges(6,5) * t923 - pkin(9) * t809 - t810 * t954 + t811 * t958 + t866 * t886 - t867 * t941;
t880 = Ifges(5,5) * t919 + Ifges(5,6) * t918 + Ifges(5,3) * t944;
t882 = Ifges(5,1) * t919 + Ifges(5,4) * t918 + Ifges(5,5) * t944;
t877 = -mrSges(6,2) * t941 + mrSges(6,3) * t886;
t878 = mrSges(6,1) * t941 - mrSges(6,3) * t887;
t968 = m(6) * t850 - t860 * mrSges(6,1) + t861 * mrSges(6,2) - t886 * t877 + t887 * t878 + t973;
t873 = -mrSges(6,1) * t886 + mrSges(6,2) * t887;
t806 = m(6) * t828 + mrSges(6,1) * t923 - mrSges(6,3) * t861 - t873 * t887 + t877 * t941 + t809;
t807 = m(6) * t829 - mrSges(6,2) * t923 + mrSges(6,3) * t860 + t873 * t886 - t878 * t941 + t977;
t978 = -t806 * t955 + t959 * t807;
t781 = -mrSges(5,1) * t865 + mrSges(5,3) * t844 + Ifges(5,4) * t895 + Ifges(5,2) * t894 + Ifges(5,6) * t927 - pkin(4) * t968 + pkin(8) * t978 + t959 * t796 + t955 * t797 - t919 * t880 + t882 * t944;
t802 = t959 * t806 + t955 * t807;
t881 = Ifges(5,4) * t919 + Ifges(5,2) * t918 + Ifges(5,6) * t944;
t782 = mrSges(5,2) * t865 - mrSges(5,3) * t843 + Ifges(5,1) * t895 + Ifges(5,4) * t894 + Ifges(5,5) * t927 - pkin(8) * t802 - t796 * t955 + t797 * t959 + t880 * t918 - t881 * t944;
t925 = (mrSges(4,2) * t960 - mrSges(4,3) * t956) * qJD(1);
t935 = -mrSges(4,1) * t985 - qJD(2) * mrSges(4,3);
t888 = -mrSges(5,1) * t918 + mrSges(5,2) * t919;
t892 = -mrSges(5,2) * t944 + mrSges(5,3) * t918;
t800 = m(5) * t843 + mrSges(5,1) * t927 - mrSges(5,3) * t895 - t888 * t919 + t892 * t944 + t802;
t893 = mrSges(5,1) * t944 - mrSges(5,3) * t919;
t801 = m(5) * t844 - mrSges(5,2) * t927 + mrSges(5,3) * t894 + t888 * t918 - t893 * t944 + t978;
t795 = t953 * t800 + t952 * t801;
t972 = -m(4) * t876 - t927 * mrSges(4,1) - t795;
t794 = qJDD(2) * mrSges(4,2) + qJD(2) * t935 + t925 * t944 - t972;
t821 = m(5) * t865 - t894 * mrSges(5,1) + t895 * mrSges(5,2) - t918 * t892 + t919 * t893 + t968;
t936 = mrSges(4,1) * t944 + qJD(2) * mrSges(4,2);
t964 = -m(4) * t875 + qJDD(2) * mrSges(4,3) + qJD(2) * t936 + t925 * t985 + t821;
t986 = t993 * qJD(2) + (t1001 * t956 + t994 * t960) * qJD(1);
t987 = t992 * qJD(2) + (t1000 * t960 + t994 * t956) * qJD(1);
t998 = (t956 * t987 - t960 * t986) * qJD(1) + t999 * qJDD(2) + t993 * t927 + t992 * t928 + mrSges(3,1) * t889 - mrSges(3,2) * t890 + mrSges(4,2) * t876 - mrSges(4,3) * t875 - pkin(2) * t794 + qJ(3) * (t928 * mrSges(4,1) + t964) - qJ(4) * t795 - t952 * t781 + t953 * t782;
t995 = t963 * pkin(7);
t926 = (-mrSges(3,1) * t960 + mrSges(3,2) * t956) * qJD(1);
t933 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t985;
t792 = m(3) * t889 - t927 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t933 - t935) * qJD(2) + (-t925 - t926) * t944 + t972;
t932 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t944;
t819 = -qJD(2) * t932 + t926 * t985 - qJDD(2) * mrSges(3,2) + m(3) * t890 + t964 + (mrSges(3,3) + mrSges(4,1)) * t928;
t979 = -t792 * t956 + t960 * t819;
t785 = m(2) * t938 - mrSges(2,1) * t963 - qJDD(1) * mrSges(2,2) + t979;
t911 = t976 - t995;
t874 = -t928 * pkin(2) + t969 - t995;
t989 = -t952 * t800 + t953 * t801;
t975 = -m(4) * t874 - t928 * mrSges(4,2) + t936 * t944 - t989;
t966 = -m(3) * t911 + t933 * t985 + t928 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t927 + (-t932 * t956 - t935 * t960) * qJD(1) + t975;
t789 = m(2) * t937 + qJDD(1) * mrSges(2,1) - t963 * mrSges(2,2) + t966;
t990 = t957 * t785 + t961 * t789;
t787 = t960 * t792 + t956 * t819;
t988 = t999 * qJD(2) + (t993 * t956 + t992 * t960) * qJD(1);
t980 = t961 * t785 - t789 * t957;
t793 = -t927 * mrSges(4,3) + t935 * t985 - t975;
t778 = -mrSges(3,1) * t911 - mrSges(4,1) * t875 + mrSges(4,2) * t874 + mrSges(3,3) * t890 - pkin(2) * t793 + pkin(3) * t821 - qJ(4) * t989 + t986 * qJD(2) + t992 * qJDD(2) + t1000 * t928 - t953 * t781 - t952 * t782 + t994 * t927 - t988 * t944;
t970 = mrSges(7,1) * t823 - mrSges(7,2) * t824 + Ifges(7,5) * t841 + Ifges(7,6) * t840 + Ifges(7,3) * t914 + t872 * t846 - t871 * t847;
t967 = mrSges(6,1) * t828 - mrSges(6,2) * t829 + Ifges(6,5) * t861 + Ifges(6,6) * t860 + Ifges(6,3) * t923 + pkin(5) * t809 + t887 * t867 - t886 * t868 + t970;
t780 = mrSges(5,1) * t843 - mrSges(5,2) * t844 + t919 * t881 - qJ(3) * t793 + mrSges(3,2) * t911 - t918 * t882 + t988 * t985 - mrSges(4,3) * t874 + mrSges(4,1) * t876 + pkin(3) * t795 - mrSges(3,3) * t889 + Ifges(5,6) * t894 + Ifges(5,5) * t895 + t993 * qJDD(2) + t994 * t928 + (Ifges(5,3) + t1001) * t927 - t987 * qJD(2) + pkin(4) * t802 + t967;
t971 = mrSges(2,1) * t937 - mrSges(2,2) * t938 + Ifges(2,3) * qJDD(1) + pkin(1) * t966 + pkin(7) * t979 + t960 * t778 + t956 * t780;
t776 = mrSges(2,1) * g(3) + mrSges(2,3) * t938 + t963 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t787 - t998;
t775 = -mrSges(2,2) * g(3) - mrSges(2,3) * t937 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t963 - pkin(7) * t787 - t778 * t956 + t780 * t960;
t1 = [-m(1) * g(1) + t980; -m(1) * g(2) + t990; (-m(1) - m(2)) * g(3) + t787; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t990 + t961 * t775 - t957 * t776; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t980 + t957 * t775 + t961 * t776; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t971; t971; t998; t794; t821; t967; t970;];
tauJB  = t1;
