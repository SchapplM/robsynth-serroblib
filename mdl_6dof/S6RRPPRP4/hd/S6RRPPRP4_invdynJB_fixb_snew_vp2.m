% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-05-06 09:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:19:19
% EndTime: 2019-05-06 09:19:30
% DurationCPUTime: 7.30s
% Computational Cost: add. (70947->338), mult. (154775->397), div. (0->0), fcn. (99511->8), ass. (0->136)
t1038 = Ifges(6,1) + Ifges(7,1);
t1019 = Ifges(6,4) - Ifges(7,5);
t1033 = Ifges(7,4) + Ifges(6,5);
t1037 = Ifges(6,2) + Ifges(7,3);
t1028 = Ifges(6,6) - Ifges(7,6);
t1036 = -2 * qJD(3);
t1035 = -2 * qJD(4);
t1034 = Ifges(4,1) + Ifges(5,1);
t1020 = Ifges(4,4) - Ifges(5,5);
t1032 = Ifges(4,5) + Ifges(5,4);
t1031 = Ifges(4,2) + Ifges(5,3);
t1030 = Ifges(5,2) + Ifges(4,3);
t1029 = Ifges(4,6) - Ifges(5,6);
t1027 = Ifges(6,3) + Ifges(7,2);
t1023 = cos(qJ(5));
t979 = sin(qJ(2));
t1007 = qJD(1) * t979;
t1016 = cos(pkin(9));
t977 = sin(pkin(9));
t950 = -qJD(2) * t1016 + t1007 * t977;
t951 = t977 * qJD(2) + t1007 * t1016;
t978 = sin(qJ(5));
t912 = -t1023 * t950 + t951 * t978;
t913 = t1023 * t951 + t978 * t950;
t981 = cos(qJ(2));
t1006 = qJD(1) * t981;
t966 = qJD(5) + t1006;
t1026 = -t1019 * t913 - t1028 * t966 + t1037 * t912;
t1025 = -t1019 * t912 + t1033 * t966 + t1038 * t913;
t1001 = t950 * t1006;
t980 = sin(qJ(1));
t982 = cos(qJ(1));
t964 = -g(1) * t982 - g(2) * t980;
t984 = qJD(1) ^ 2;
t942 = -pkin(1) * t984 + qJDD(1) * pkin(7) + t964;
t918 = -t981 * g(3) - t979 * t942;
t956 = (-pkin(2) * t981 - qJ(3) * t979) * qJD(1);
t983 = qJD(2) ^ 2;
t895 = -qJDD(2) * pkin(2) - t983 * qJ(3) + t956 * t1007 + qJDD(3) - t918;
t1003 = qJD(1) * qJD(2);
t1000 = t981 * t1003;
t958 = qJDD(1) * t979 + t1000;
t927 = -qJDD(2) * t1016 + t958 * t977;
t928 = t977 * qJDD(2) + t1016 * t958;
t865 = t895 - (t1001 + t928) * qJ(4) - (t1006 * t951 - t927) * pkin(3) + t951 * t1035;
t963 = t980 * g(1) - t982 * g(2);
t941 = -qJDD(1) * pkin(1) - t984 * pkin(7) - t963;
t968 = t979 * t1003;
t959 = qJDD(1) * t981 - t968;
t890 = (-t958 - t1000) * qJ(3) + (-t959 + t968) * pkin(2) + t941;
t919 = -g(3) * t979 + t981 * t942;
t896 = -pkin(2) * t983 + qJDD(2) * qJ(3) + t1006 * t956 + t919;
t866 = t1016 * t890 + t1036 * t951 - t977 * t896;
t1009 = t1006 * t1032 + t1020 * t950 - t1034 * t951;
t1010 = t1006 * t1030 + t1029 * t950 - t1032 * t951;
t1013 = -t1027 * t966 + t1028 * t912 - t1033 * t913;
t929 = pkin(4) * t1006 - pkin(8) * t951;
t948 = t950 ^ 2;
t861 = -t927 * pkin(4) - t948 * pkin(8) + t951 * t929 - t865;
t876 = qJD(5) * t913 - t1023 * t927 + t928 * t978;
t877 = -t912 * qJD(5) + t1023 * t928 + t978 * t927;
t855 = -0.2e1 * qJD(6) * t913 + t861 + (t912 * t966 - t877) * qJ(6) + (t913 * t966 + t876) * pkin(5);
t899 = -mrSges(7,1) * t966 + mrSges(7,2) * t913;
t900 = -mrSges(7,2) * t912 + mrSges(7,3) * t966;
t846 = m(7) * t855 + t876 * mrSges(7,1) - t877 * mrSges(7,3) - t913 * t899 + t912 * t900;
t1015 = t981 ^ 2 * t984;
t914 = pkin(3) * t950 - qJ(4) * t951;
t864 = t959 * pkin(3) - qJ(4) * t1015 + t951 * t914 + qJDD(4) - t866;
t857 = (-t928 + t1001) * pkin(8) + (t950 * t951 + t959) * pkin(4) + t864;
t867 = t1016 * t896 + t1036 * t950 + t977 * t890;
t863 = -pkin(3) * t1015 - t959 * qJ(4) + t1006 * t1035 - t950 * t914 + t867;
t859 = -pkin(4) * t948 + pkin(8) * t927 - t1006 * t929 + t863;
t853 = t1023 * t859 + t978 * t857;
t886 = pkin(5) * t912 - qJ(6) * t913;
t955 = qJDD(5) + t959;
t965 = t966 ^ 2;
t849 = -pkin(5) * t965 + qJ(6) * t955 + 0.2e1 * qJD(6) * t966 - t886 * t912 + t853;
t833 = -mrSges(6,1) * t861 - mrSges(7,1) * t855 + mrSges(7,2) * t849 + mrSges(6,3) * t853 - pkin(5) * t846 + t1013 * t913 + t1019 * t877 + t1025 * t966 + t1028 * t955 - t1037 * t876;
t852 = t1023 * t857 - t978 * t859;
t850 = -t955 * pkin(5) - t965 * qJ(6) + t913 * t886 + qJDD(6) - t852;
t834 = mrSges(6,2) * t861 + mrSges(7,2) * t850 - mrSges(6,3) * t852 - mrSges(7,3) * t855 - qJ(6) * t846 + t1013 * t912 - t1019 * t876 + t1026 * t966 + t1033 * t955 + t1038 * t877;
t924 = -mrSges(5,2) * t950 - mrSges(5,3) * t1006;
t926 = mrSges(5,1) * t1006 + mrSges(5,2) * t951;
t897 = -mrSges(6,2) * t966 - mrSges(6,3) * t912;
t898 = mrSges(6,1) * t966 - mrSges(6,3) * t913;
t988 = m(6) * t861 + t876 * mrSges(6,1) + t877 * mrSges(6,2) + t912 * t897 + t913 * t898 + t846;
t842 = m(5) * t865 + t927 * mrSges(5,1) - t928 * mrSges(5,3) + t950 * t924 - t951 * t926 - t988;
t1002 = m(7) * t849 + t955 * mrSges(7,3) + t966 * t899;
t887 = mrSges(7,1) * t912 - mrSges(7,3) * t913;
t1012 = -mrSges(6,1) * t912 - mrSges(6,2) * t913 - t887;
t1021 = -mrSges(6,3) - mrSges(7,2);
t840 = m(6) * t853 - t955 * mrSges(6,2) + t1012 * t912 + t1021 * t876 - t966 * t898 + t1002;
t996 = -m(7) * t850 + t955 * mrSges(7,1) + t966 * t900;
t841 = m(6) * t852 + t955 * mrSges(6,1) + t1012 * t913 + t1021 * t877 + t966 * t897 + t996;
t997 = t1023 * t840 - t978 * t841;
t814 = -mrSges(4,1) * t895 - mrSges(5,1) * t865 + mrSges(5,2) * t863 + mrSges(4,3) * t867 - pkin(3) * t842 + pkin(4) * t988 - pkin(8) * t997 + t1009 * t1006 + t1010 * t951 + t1020 * t928 - t1023 * t833 - t1029 * t959 - t1031 * t927 - t978 * t834;
t1011 = t1006 * t1029 - t1020 * t951 + t1031 * t950;
t835 = t1023 * t841 + t978 * t840;
t815 = mrSges(4,2) * t895 + mrSges(5,2) * t864 - mrSges(4,3) * t866 - mrSges(5,3) * t865 - pkin(8) * t835 - qJ(4) * t842 - t1011 * t1006 + t1010 * t950 - t1020 * t927 + t1023 * t834 - t1032 * t959 + t1034 * t928 - t978 * t833;
t915 = mrSges(5,1) * t950 - mrSges(5,3) * t951;
t1008 = -mrSges(4,1) * t950 - mrSges(4,2) * t951 - t915;
t1022 = -mrSges(4,3) - mrSges(5,2);
t925 = -mrSges(4,1) * t1006 - mrSges(4,3) * t951;
t993 = m(5) * t863 - t959 * mrSges(5,3) + t997;
t830 = m(4) * t867 + t959 * mrSges(4,2) + t1008 * t950 + t1022 * t927 + (t925 - t926) * t1006 + t993;
t990 = m(5) * t864 + t959 * mrSges(5,1) + t1006 * t924 + t835;
t992 = mrSges(4,2) * t1006 - mrSges(4,3) * t950;
t831 = m(4) * t866 - t959 * mrSges(4,1) - t1006 * t992 + t1008 * t951 + t1022 * t928 - t990;
t828 = t1016 * t830 - t831 * t977;
t838 = m(4) * t895 + t927 * mrSges(4,1) + t928 * mrSges(4,2) + t951 * t925 + t950 * t992 + t842;
t939 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t979 + Ifges(3,2) * t981) * qJD(1);
t940 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t979 + Ifges(3,4) * t981) * qJD(1);
t1024 = mrSges(3,1) * t918 - mrSges(3,2) * t919 + Ifges(3,5) * t958 + Ifges(3,6) * t959 + Ifges(3,3) * qJDD(2) - pkin(2) * t838 + qJ(3) * t828 + (t939 * t979 - t940 * t981) * qJD(1) + t1016 * t814 + t977 * t815;
t957 = (-mrSges(3,1) * t981 + mrSges(3,2) * t979) * qJD(1);
t961 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1007;
t826 = m(3) * t919 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t959 - qJD(2) * t961 + t1006 * t957 + t828;
t962 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1006;
t837 = m(3) * t918 + qJDD(2) * mrSges(3,1) - t958 * mrSges(3,3) + qJD(2) * t962 - t1007 * t957 - t838;
t998 = t981 * t826 - t837 * t979;
t818 = m(2) * t964 - mrSges(2,1) * t984 - qJDD(1) * mrSges(2,2) + t998;
t827 = t1016 * t831 + t977 * t830;
t987 = -m(3) * t941 + t959 * mrSges(3,1) - t958 * mrSges(3,2) + t962 * t1006 - t1007 * t961 - t827;
t822 = m(2) * t963 + qJDD(1) * mrSges(2,1) - t984 * mrSges(2,2) + t987;
t1014 = t980 * t818 + t982 * t822;
t820 = t979 * t826 + t981 * t837;
t999 = t982 * t818 - t822 * t980;
t938 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t979 + Ifges(3,6) * t981) * qJD(1);
t811 = mrSges(3,2) * t941 - mrSges(3,3) * t918 + Ifges(3,1) * t958 + Ifges(3,4) * t959 + Ifges(3,5) * qJDD(2) - qJ(3) * t827 - qJD(2) * t939 + t1006 * t938 + t1016 * t815 - t977 * t814;
t832 = t928 * mrSges(5,2) + t951 * t915 + t990;
t845 = t877 * mrSges(7,2) + t913 * t887 - t996;
t986 = mrSges(6,1) * t852 - mrSges(7,1) * t850 - mrSges(6,2) * t853 + mrSges(7,3) * t849 - pkin(5) * t845 + qJ(6) * t1002 + t1027 * t955 - t1026 * t913 + (-qJ(6) * t887 + t1025) * t912 + t1033 * t877 + (-mrSges(7,2) * qJ(6) - t1028) * t876;
t813 = Ifges(3,4) * t958 + qJD(2) * t940 - mrSges(3,1) * t941 + mrSges(3,3) * t919 - mrSges(5,3) * t863 + mrSges(5,1) * t864 - mrSges(4,1) * t866 + mrSges(4,2) * t867 + pkin(4) * t835 + pkin(3) * t832 - t938 * t1007 - pkin(2) * t827 + t986 + Ifges(3,6) * qJDD(2) - qJ(4) * (-t1006 * t926 + t993) + (Ifges(3,2) + t1030) * t959 + t1011 * t951 + (mrSges(5,2) * qJ(4) + t1029) * t927 - t1032 * t928 + (qJ(4) * t915 + t1009) * t950;
t991 = mrSges(2,1) * t963 - mrSges(2,2) * t964 + Ifges(2,3) * qJDD(1) + pkin(1) * t987 + pkin(7) * t998 + t979 * t811 + t981 * t813;
t809 = mrSges(2,1) * g(3) + mrSges(2,3) * t964 + t984 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t820 - t1024;
t808 = -mrSges(2,2) * g(3) - mrSges(2,3) * t963 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t984 - pkin(7) * t820 + t811 * t981 - t813 * t979;
t1 = [-m(1) * g(1) + t999; -m(1) * g(2) + t1014; (-m(1) - m(2)) * g(3) + t820; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1014 + t982 * t808 - t980 * t809; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t999 + t980 * t808 + t982 * t809; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t991; t991; t1024; t838; t832; t986; t845;];
tauJB  = t1;
