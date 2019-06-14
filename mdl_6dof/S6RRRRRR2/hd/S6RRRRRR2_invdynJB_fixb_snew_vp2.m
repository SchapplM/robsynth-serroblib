% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 08:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:27:06
% EndTime: 2019-05-08 08:27:48
% DurationCPUTime: 39.49s
% Computational Cost: add. (640881->391), mult. (1381274->489), div. (0->0), fcn. (1054430->12), ass. (0->157)
t965 = sin(qJ(2));
t971 = cos(qJ(2));
t993 = qJD(1) * qJD(2);
t943 = qJDD(1) * t965 + t971 * t993;
t966 = sin(qJ(1));
t972 = cos(qJ(1));
t950 = -g(1) * t972 - g(2) * t966;
t973 = qJD(1) ^ 2;
t938 = -pkin(1) * t973 + qJDD(1) * pkin(7) + t950;
t997 = t965 * t938;
t998 = pkin(2) * t973;
t899 = qJDD(2) * pkin(2) - t943 * pkin(8) - t997 + (pkin(8) * t993 + t965 * t998 - g(3)) * t971;
t925 = -g(3) * t965 + t971 * t938;
t944 = qJDD(1) * t971 - t965 * t993;
t995 = qJD(1) * t965;
t948 = qJD(2) * pkin(2) - pkin(8) * t995;
t960 = t971 ^ 2;
t900 = pkin(8) * t944 - qJD(2) * t948 - t960 * t998 + t925;
t964 = sin(qJ(3));
t970 = cos(qJ(3));
t880 = t970 * t899 - t964 * t900;
t935 = (-t964 * t965 + t970 * t971) * qJD(1);
t909 = qJD(3) * t935 + t943 * t970 + t944 * t964;
t936 = (t964 * t971 + t965 * t970) * qJD(1);
t957 = qJDD(2) + qJDD(3);
t958 = qJD(2) + qJD(3);
t852 = (t935 * t958 - t909) * pkin(9) + (t935 * t936 + t957) * pkin(3) + t880;
t881 = t964 * t899 + t970 * t900;
t908 = -qJD(3) * t936 - t943 * t964 + t944 * t970;
t928 = pkin(3) * t958 - pkin(9) * t936;
t931 = t935 ^ 2;
t857 = -pkin(3) * t931 + pkin(9) * t908 - t928 * t958 + t881;
t963 = sin(qJ(4));
t969 = cos(qJ(4));
t845 = t963 * t852 + t969 * t857;
t922 = t935 * t963 + t936 * t969;
t875 = -t922 * qJD(4) + t908 * t969 - t963 * t909;
t921 = t935 * t969 - t963 * t936;
t893 = -mrSges(5,1) * t921 + mrSges(5,2) * t922;
t955 = qJD(4) + t958;
t912 = mrSges(5,1) * t955 - mrSges(5,3) * t922;
t954 = qJDD(4) + t957;
t894 = -pkin(4) * t921 - pkin(10) * t922;
t953 = t955 ^ 2;
t834 = -pkin(4) * t953 + pkin(10) * t954 + t894 * t921 + t845;
t949 = t966 * g(1) - t972 * g(2);
t985 = -qJDD(1) * pkin(1) - t949;
t910 = -t944 * pkin(2) + t948 * t995 + (-pkin(8) * t960 - pkin(7)) * t973 + t985;
t861 = -t908 * pkin(3) - t931 * pkin(9) + t936 * t928 + t910;
t876 = qJD(4) * t921 + t908 * t963 + t909 * t969;
t842 = (-t921 * t955 - t876) * pkin(10) + (t922 * t955 - t875) * pkin(4) + t861;
t962 = sin(qJ(5));
t968 = cos(qJ(5));
t829 = -t962 * t834 + t968 * t842;
t903 = -t922 * t962 + t955 * t968;
t859 = qJD(5) * t903 + t876 * t968 + t954 * t962;
t873 = qJDD(5) - t875;
t904 = t922 * t968 + t955 * t962;
t918 = qJD(5) - t921;
t827 = (t903 * t918 - t859) * pkin(11) + (t903 * t904 + t873) * pkin(5) + t829;
t830 = t968 * t834 + t962 * t842;
t858 = -qJD(5) * t904 - t876 * t962 + t954 * t968;
t887 = pkin(5) * t918 - pkin(11) * t904;
t902 = t903 ^ 2;
t828 = -pkin(5) * t902 + pkin(11) * t858 - t887 * t918 + t830;
t961 = sin(qJ(6));
t967 = cos(qJ(6));
t825 = t827 * t967 - t828 * t961;
t882 = t903 * t967 - t904 * t961;
t841 = qJD(6) * t882 + t858 * t961 + t859 * t967;
t883 = t903 * t961 + t904 * t967;
t853 = -mrSges(7,1) * t882 + mrSges(7,2) * t883;
t913 = qJD(6) + t918;
t862 = -mrSges(7,2) * t913 + mrSges(7,3) * t882;
t865 = qJDD(6) + t873;
t820 = m(7) * t825 + mrSges(7,1) * t865 - mrSges(7,3) * t841 - t853 * t883 + t862 * t913;
t826 = t827 * t961 + t828 * t967;
t840 = -qJD(6) * t883 + t858 * t967 - t859 * t961;
t863 = mrSges(7,1) * t913 - mrSges(7,3) * t883;
t821 = m(7) * t826 - mrSges(7,2) * t865 + mrSges(7,3) * t840 + t853 * t882 - t863 * t913;
t812 = t967 * t820 + t961 * t821;
t884 = -mrSges(6,1) * t903 + mrSges(6,2) * t904;
t885 = -mrSges(6,2) * t918 + mrSges(6,3) * t903;
t810 = m(6) * t829 + mrSges(6,1) * t873 - mrSges(6,3) * t859 - t884 * t904 + t885 * t918 + t812;
t886 = mrSges(6,1) * t918 - mrSges(6,3) * t904;
t987 = -t820 * t961 + t967 * t821;
t811 = m(6) * t830 - mrSges(6,2) * t873 + mrSges(6,3) * t858 + t884 * t903 - t886 * t918 + t987;
t988 = -t810 * t962 + t968 * t811;
t803 = m(5) * t845 - mrSges(5,2) * t954 + mrSges(5,3) * t875 + t893 * t921 - t912 * t955 + t988;
t844 = t852 * t969 - t963 * t857;
t911 = -mrSges(5,2) * t955 + mrSges(5,3) * t921;
t833 = -pkin(4) * t954 - pkin(10) * t953 + t922 * t894 - t844;
t831 = -pkin(5) * t858 - pkin(11) * t902 + t887 * t904 + t833;
t982 = m(7) * t831 - t840 * mrSges(7,1) + mrSges(7,2) * t841 - t882 * t862 + t863 * t883;
t978 = -m(6) * t833 + t858 * mrSges(6,1) - mrSges(6,2) * t859 + t903 * t885 - t886 * t904 - t982;
t816 = m(5) * t844 + mrSges(5,1) * t954 - mrSges(5,3) * t876 - t893 * t922 + t911 * t955 + t978;
t793 = t963 * t803 + t969 * t816;
t923 = -mrSges(4,1) * t935 + mrSges(4,2) * t936;
t926 = -mrSges(4,2) * t958 + mrSges(4,3) * t935;
t790 = m(4) * t880 + mrSges(4,1) * t957 - mrSges(4,3) * t909 - t923 * t936 + t926 * t958 + t793;
t927 = mrSges(4,1) * t958 - mrSges(4,3) * t936;
t989 = t969 * t803 - t816 * t963;
t791 = m(4) * t881 - mrSges(4,2) * t957 + mrSges(4,3) * t908 + t923 * t935 - t927 * t958 + t989;
t785 = t970 * t790 + t964 * t791;
t924 = -t971 * g(3) - t997;
t933 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t965 + Ifges(3,2) * t971) * qJD(1);
t934 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t965 + Ifges(3,4) * t971) * qJD(1);
t915 = Ifges(4,4) * t936 + Ifges(4,2) * t935 + Ifges(4,6) * t958;
t916 = Ifges(4,1) * t936 + Ifges(4,4) * t935 + Ifges(4,5) * t958;
t847 = Ifges(7,5) * t883 + Ifges(7,6) * t882 + Ifges(7,3) * t913;
t849 = Ifges(7,1) * t883 + Ifges(7,4) * t882 + Ifges(7,5) * t913;
t813 = -mrSges(7,1) * t831 + mrSges(7,3) * t826 + Ifges(7,4) * t841 + Ifges(7,2) * t840 + Ifges(7,6) * t865 - t847 * t883 + t849 * t913;
t848 = Ifges(7,4) * t883 + Ifges(7,2) * t882 + Ifges(7,6) * t913;
t814 = mrSges(7,2) * t831 - mrSges(7,3) * t825 + Ifges(7,1) * t841 + Ifges(7,4) * t840 + Ifges(7,5) * t865 + t847 * t882 - t848 * t913;
t866 = Ifges(6,5) * t904 + Ifges(6,6) * t903 + Ifges(6,3) * t918;
t868 = Ifges(6,1) * t904 + Ifges(6,4) * t903 + Ifges(6,5) * t918;
t795 = -mrSges(6,1) * t833 + mrSges(6,3) * t830 + Ifges(6,4) * t859 + Ifges(6,2) * t858 + Ifges(6,6) * t873 - pkin(5) * t982 + pkin(11) * t987 + t967 * t813 + t961 * t814 - t904 * t866 + t918 * t868;
t867 = Ifges(6,4) * t904 + Ifges(6,2) * t903 + Ifges(6,6) * t918;
t797 = mrSges(6,2) * t833 - mrSges(6,3) * t829 + Ifges(6,1) * t859 + Ifges(6,4) * t858 + Ifges(6,5) * t873 - pkin(11) * t812 - t813 * t961 + t814 * t967 + t866 * t903 - t867 * t918;
t889 = Ifges(5,4) * t922 + Ifges(5,2) * t921 + Ifges(5,6) * t955;
t890 = Ifges(5,1) * t922 + Ifges(5,4) * t921 + Ifges(5,5) * t955;
t980 = -mrSges(5,1) * t844 + mrSges(5,2) * t845 - Ifges(5,5) * t876 - Ifges(5,6) * t875 - Ifges(5,3) * t954 - pkin(4) * t978 - pkin(10) * t988 - t968 * t795 - t962 * t797 - t922 * t889 + t921 * t890;
t977 = -mrSges(4,1) * t880 + mrSges(4,2) * t881 - Ifges(4,5) * t909 - Ifges(4,6) * t908 - Ifges(4,3) * t957 - pkin(3) * t793 - t936 * t915 + t935 * t916 + t980;
t999 = mrSges(3,1) * t924 - mrSges(3,2) * t925 + Ifges(3,5) * t943 + Ifges(3,6) * t944 + Ifges(3,3) * qJDD(2) + pkin(2) * t785 + (t933 * t965 - t934 * t971) * qJD(1) - t977;
t942 = (-mrSges(3,1) * t971 + mrSges(3,2) * t965) * qJD(1);
t994 = qJD(1) * t971;
t947 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t994;
t783 = m(3) * t924 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t943 + qJD(2) * t947 - t942 * t995 + t785;
t946 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t995;
t990 = -t790 * t964 + t970 * t791;
t784 = m(3) * t925 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t944 - qJD(2) * t946 + t942 * t994 + t990;
t991 = -t783 * t965 + t971 * t784;
t776 = m(2) * t950 - mrSges(2,1) * t973 - qJDD(1) * mrSges(2,2) + t991;
t937 = -t973 * pkin(7) + t985;
t805 = t968 * t810 + t962 * t811;
t984 = m(5) * t861 - t875 * mrSges(5,1) + t876 * mrSges(5,2) - t921 * t911 + t922 * t912 + t805;
t979 = m(4) * t910 - t908 * mrSges(4,1) + mrSges(4,2) * t909 - t935 * t926 + t927 * t936 + t984;
t976 = -m(3) * t937 + t944 * mrSges(3,1) - mrSges(3,2) * t943 - t946 * t995 + t947 * t994 - t979;
t799 = m(2) * t949 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t973 + t976;
t996 = t966 * t776 + t972 * t799;
t778 = t971 * t783 + t965 * t784;
t992 = t972 * t776 - t799 * t966;
t888 = Ifges(5,5) * t922 + Ifges(5,6) * t921 + Ifges(5,3) * t955;
t779 = mrSges(5,2) * t861 - mrSges(5,3) * t844 + Ifges(5,1) * t876 + Ifges(5,4) * t875 + Ifges(5,5) * t954 - pkin(10) * t805 - t795 * t962 + t797 * t968 + t888 * t921 - t889 * t955;
t981 = -mrSges(7,1) * t825 + mrSges(7,2) * t826 - Ifges(7,5) * t841 - Ifges(7,6) * t840 - Ifges(7,3) * t865 - t883 * t848 + t882 * t849;
t975 = mrSges(6,1) * t829 - mrSges(6,2) * t830 + Ifges(6,5) * t859 + Ifges(6,6) * t858 + Ifges(6,3) * t873 + pkin(5) * t812 + t904 * t867 - t903 * t868 - t981;
t786 = -mrSges(5,1) * t861 + mrSges(5,3) * t845 + Ifges(5,4) * t876 + Ifges(5,2) * t875 + Ifges(5,6) * t954 - pkin(4) * t805 - t922 * t888 + t955 * t890 - t975;
t914 = Ifges(4,5) * t936 + Ifges(4,6) * t935 + Ifges(4,3) * t958;
t772 = -mrSges(4,1) * t910 + mrSges(4,3) * t881 + Ifges(4,4) * t909 + Ifges(4,2) * t908 + Ifges(4,6) * t957 - pkin(3) * t984 + pkin(9) * t989 + t963 * t779 + t969 * t786 - t936 * t914 + t958 * t916;
t773 = mrSges(4,2) * t910 - mrSges(4,3) * t880 + Ifges(4,1) * t909 + Ifges(4,4) * t908 + Ifges(4,5) * t957 - pkin(9) * t793 + t779 * t969 - t786 * t963 + t914 * t935 - t915 * t958;
t932 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t965 + Ifges(3,6) * t971) * qJD(1);
t768 = -mrSges(3,1) * t937 + mrSges(3,3) * t925 + Ifges(3,4) * t943 + Ifges(3,2) * t944 + Ifges(3,6) * qJDD(2) - pkin(2) * t979 + pkin(8) * t990 + qJD(2) * t934 + t970 * t772 + t964 * t773 - t932 * t995;
t770 = mrSges(3,2) * t937 - mrSges(3,3) * t924 + Ifges(3,1) * t943 + Ifges(3,4) * t944 + Ifges(3,5) * qJDD(2) - pkin(8) * t785 - qJD(2) * t933 - t772 * t964 + t773 * t970 + t932 * t994;
t983 = mrSges(2,1) * t949 - mrSges(2,2) * t950 + Ifges(2,3) * qJDD(1) + pkin(1) * t976 + pkin(7) * t991 + t971 * t768 + t965 * t770;
t771 = mrSges(2,1) * g(3) + mrSges(2,3) * t950 + t973 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t778 - t999;
t766 = -mrSges(2,2) * g(3) - mrSges(2,3) * t949 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t973 - pkin(7) * t778 - t768 * t965 + t770 * t971;
t1 = [-m(1) * g(1) + t992; -m(1) * g(2) + t996; (-m(1) - m(2)) * g(3) + t778; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t996 + t972 * t766 - t966 * t771; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t992 + t966 * t766 + t972 * t771; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t983; t983; t999; -t977; -t980; t975; -t981;];
tauJB  = t1;
