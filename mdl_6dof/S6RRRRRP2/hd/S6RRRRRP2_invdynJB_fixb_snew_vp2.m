% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:31:01
% EndTime: 2019-05-08 04:31:24
% DurationCPUTime: 17.93s
% Computational Cost: add. (281289->366), mult. (609201->447), div. (0->0), fcn. (452999->10), ass. (0->147)
t984 = Ifges(6,1) + Ifges(7,1);
t976 = Ifges(6,4) - Ifges(7,5);
t975 = -Ifges(6,5) - Ifges(7,4);
t983 = Ifges(6,2) + Ifges(7,3);
t974 = Ifges(6,6) - Ifges(7,6);
t982 = -Ifges(6,3) - Ifges(7,2);
t939 = sin(qJ(2));
t943 = cos(qJ(2));
t964 = qJD(1) * qJD(2);
t916 = qJDD(1) * t939 + t943 * t964;
t940 = sin(qJ(1));
t944 = cos(qJ(1));
t923 = -g(1) * t944 - g(2) * t940;
t945 = qJD(1) ^ 2;
t911 = -pkin(1) * t945 + qJDD(1) * pkin(7) + t923;
t972 = t939 * t911;
t978 = pkin(2) * t945;
t876 = qJDD(2) * pkin(2) - t916 * pkin(8) - t972 + (pkin(8) * t964 + t939 * t978 - g(3)) * t943;
t899 = -g(3) * t939 + t943 * t911;
t917 = qJDD(1) * t943 - t939 * t964;
t966 = qJD(1) * t939;
t921 = qJD(2) * pkin(2) - pkin(8) * t966;
t935 = t943 ^ 2;
t877 = pkin(8) * t917 - qJD(2) * t921 - t935 * t978 + t899;
t938 = sin(qJ(3));
t942 = cos(qJ(3));
t856 = t942 * t876 - t938 * t877;
t908 = (-t938 * t939 + t942 * t943) * qJD(1);
t884 = qJD(3) * t908 + t916 * t942 + t917 * t938;
t909 = (t938 * t943 + t939 * t942) * qJD(1);
t932 = qJDD(2) + qJDD(3);
t933 = qJD(2) + qJD(3);
t824 = (t908 * t933 - t884) * pkin(9) + (t908 * t909 + t932) * pkin(3) + t856;
t857 = t938 * t876 + t942 * t877;
t883 = -qJD(3) * t909 - t916 * t938 + t917 * t942;
t902 = pkin(3) * t933 - pkin(9) * t909;
t904 = t908 ^ 2;
t829 = -pkin(3) * t904 + pkin(9) * t883 - t902 * t933 + t857;
t937 = sin(qJ(4));
t941 = cos(qJ(4));
t822 = t937 * t824 + t941 * t829;
t896 = t908 * t937 + t909 * t941;
t849 = -qJD(4) * t896 + t883 * t941 - t884 * t937;
t895 = t908 * t941 - t909 * t937;
t870 = -mrSges(5,1) * t895 + mrSges(5,2) * t896;
t930 = qJD(4) + t933;
t887 = mrSges(5,1) * t930 - mrSges(5,3) * t896;
t929 = qJDD(4) + t932;
t871 = -pkin(4) * t895 - pkin(10) * t896;
t928 = t930 ^ 2;
t817 = -pkin(4) * t928 + pkin(10) * t929 + t871 * t895 + t822;
t922 = t940 * g(1) - t944 * g(2);
t955 = -qJDD(1) * pkin(1) - t922;
t885 = -t917 * pkin(2) + t921 * t966 + (-pkin(8) * t935 - pkin(7)) * t945 + t955;
t833 = -t883 * pkin(3) - t904 * pkin(9) + t909 * t902 + t885;
t850 = qJD(4) * t895 + t883 * t937 + t884 * t941;
t819 = (-t895 * t930 - t850) * pkin(10) + (t896 * t930 - t849) * pkin(4) + t833;
t936 = sin(qJ(5));
t979 = cos(qJ(5));
t814 = t979 * t817 + t936 * t819;
t879 = t896 * t979 + t936 * t930;
t830 = qJD(5) * t879 + t850 * t936 - t929 * t979;
t846 = qJDD(5) - t849;
t892 = qJD(5) - t895;
t863 = mrSges(6,1) * t892 - mrSges(6,3) * t879;
t878 = t896 * t936 - t930 * t979;
t858 = pkin(5) * t878 - qJ(6) * t879;
t891 = t892 ^ 2;
t810 = -pkin(5) * t891 + qJ(6) * t846 + 0.2e1 * qJD(6) * t892 - t858 * t878 + t814;
t864 = -mrSges(7,1) * t892 + mrSges(7,2) * t879;
t963 = m(7) * t810 + t846 * mrSges(7,3) + t892 * t864;
t859 = mrSges(7,1) * t878 - mrSges(7,3) * t879;
t967 = -mrSges(6,1) * t878 - mrSges(6,2) * t879 - t859;
t977 = -mrSges(6,3) - mrSges(7,2);
t801 = m(6) * t814 - t846 * mrSges(6,2) + t830 * t977 - t892 * t863 + t878 * t967 + t963;
t813 = -t936 * t817 + t819 * t979;
t831 = -t878 * qJD(5) + t850 * t979 + t936 * t929;
t862 = -mrSges(6,2) * t892 - mrSges(6,3) * t878;
t811 = -t846 * pkin(5) - t891 * qJ(6) + t879 * t858 + qJDD(6) - t813;
t861 = -mrSges(7,2) * t878 + mrSges(7,3) * t892;
t957 = -m(7) * t811 + t846 * mrSges(7,1) + t892 * t861;
t803 = m(6) * t813 + t846 * mrSges(6,1) + t831 * t977 + t892 * t862 + t879 * t967 + t957;
t958 = t979 * t801 - t803 * t936;
t789 = m(5) * t822 - mrSges(5,2) * t929 + mrSges(5,3) * t849 + t870 * t895 - t887 * t930 + t958;
t821 = t941 * t824 - t937 * t829;
t886 = -mrSges(5,2) * t930 + mrSges(5,3) * t895;
t816 = -t929 * pkin(4) - t928 * pkin(10) + t896 * t871 - t821;
t812 = -0.2e1 * qJD(6) * t879 + (t878 * t892 - t831) * qJ(6) + (t879 * t892 + t830) * pkin(5) + t816;
t808 = m(7) * t812 + mrSges(7,1) * t830 - t831 * mrSges(7,3) + t861 * t878 - t879 * t864;
t949 = -m(6) * t816 - t830 * mrSges(6,1) - mrSges(6,2) * t831 - t878 * t862 - t863 * t879 - t808;
t798 = m(5) * t821 + mrSges(5,1) * t929 - mrSges(5,3) * t850 - t870 * t896 + t886 * t930 + t949;
t783 = t937 * t789 + t941 * t798;
t897 = -mrSges(4,1) * t908 + mrSges(4,2) * t909;
t900 = -mrSges(4,2) * t933 + mrSges(4,3) * t908;
t780 = m(4) * t856 + mrSges(4,1) * t932 - mrSges(4,3) * t884 - t897 * t909 + t900 * t933 + t783;
t901 = mrSges(4,1) * t933 - mrSges(4,3) * t909;
t959 = t941 * t789 - t798 * t937;
t781 = m(4) * t857 - mrSges(4,2) * t932 + mrSges(4,3) * t883 + t897 * t908 - t901 * t933 + t959;
t774 = t942 * t780 + t938 * t781;
t898 = -t943 * g(3) - t972;
t906 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t939 + Ifges(3,2) * t943) * qJD(1);
t907 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t939 + Ifges(3,4) * t943) * qJD(1);
t889 = Ifges(4,4) * t909 + Ifges(4,2) * t908 + Ifges(4,6) * t933;
t890 = Ifges(4,1) * t909 + Ifges(4,4) * t908 + Ifges(4,5) * t933;
t968 = t976 * t878 - t879 * t984 + t975 * t892;
t970 = t878 * t974 + t879 * t975 + t892 * t982;
t791 = -mrSges(6,1) * t816 - mrSges(7,1) * t812 + mrSges(7,2) * t810 + mrSges(6,3) * t814 - pkin(5) * t808 - t830 * t983 + t976 * t831 + t974 * t846 + t970 * t879 - t968 * t892;
t969 = t878 * t983 - t879 * t976 - t892 * t974;
t793 = mrSges(6,2) * t816 + mrSges(7,2) * t811 - mrSges(6,3) * t813 - mrSges(7,3) * t812 - qJ(6) * t808 - t976 * t830 + t831 * t984 - t975 * t846 + t970 * t878 + t969 * t892;
t866 = Ifges(5,4) * t896 + Ifges(5,2) * t895 + Ifges(5,6) * t930;
t867 = Ifges(5,1) * t896 + Ifges(5,4) * t895 + Ifges(5,5) * t930;
t952 = -mrSges(5,1) * t821 + mrSges(5,2) * t822 - Ifges(5,5) * t850 - Ifges(5,6) * t849 - Ifges(5,3) * t929 - pkin(4) * t949 - pkin(10) * t958 - t979 * t791 - t936 * t793 - t896 * t866 + t895 * t867;
t948 = -mrSges(4,1) * t856 + mrSges(4,2) * t857 - Ifges(4,5) * t884 - Ifges(4,6) * t883 - Ifges(4,3) * t932 - pkin(3) * t783 - t909 * t889 + t908 * t890 + t952;
t981 = mrSges(3,1) * t898 - mrSges(3,2) * t899 + Ifges(3,5) * t916 + Ifges(3,6) * t917 + Ifges(3,3) * qJDD(2) + pkin(2) * t774 + (t906 * t939 - t907 * t943) * qJD(1) - t948;
t807 = t831 * mrSges(7,2) + t879 * t859 - t957;
t980 = -t830 * t974 - t831 * t975 - t982 * t846 - t878 * t968 - t879 * t969 + mrSges(6,1) * t813 - mrSges(7,1) * t811 - mrSges(6,2) * t814 + mrSges(7,3) * t810 - pkin(5) * t807 + qJ(6) * (-t830 * mrSges(7,2) - t878 * t859 + t963);
t915 = (-mrSges(3,1) * t943 + mrSges(3,2) * t939) * qJD(1);
t965 = qJD(1) * t943;
t920 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t965;
t772 = m(3) * t898 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t916 + qJD(2) * t920 - t915 * t966 + t774;
t919 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t966;
t960 = -t780 * t938 + t942 * t781;
t773 = m(3) * t899 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t917 - qJD(2) * t919 + t915 * t965 + t960;
t961 = -t772 * t939 + t943 * t773;
t766 = m(2) * t923 - mrSges(2,1) * t945 - qJDD(1) * mrSges(2,2) + t961;
t910 = -t945 * pkin(7) + t955;
t795 = t936 * t801 + t979 * t803;
t954 = m(5) * t833 - t849 * mrSges(5,1) + t850 * mrSges(5,2) - t895 * t886 + t896 * t887 + t795;
t950 = m(4) * t885 - t883 * mrSges(4,1) + mrSges(4,2) * t884 - t908 * t900 + t901 * t909 + t954;
t947 = -m(3) * t910 + t917 * mrSges(3,1) - mrSges(3,2) * t916 - t919 * t966 + t920 * t965 - t950;
t785 = m(2) * t922 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t945 + t947;
t971 = t940 * t766 + t944 * t785;
t768 = t943 * t772 + t939 * t773;
t962 = t944 * t766 - t785 * t940;
t865 = Ifges(5,5) * t896 + Ifges(5,6) * t895 + Ifges(5,3) * t930;
t775 = mrSges(5,2) * t833 - mrSges(5,3) * t821 + Ifges(5,1) * t850 + Ifges(5,4) * t849 + Ifges(5,5) * t929 - pkin(10) * t795 - t936 * t791 + t793 * t979 + t895 * t865 - t930 * t866;
t776 = -mrSges(5,1) * t833 + mrSges(5,3) * t822 + Ifges(5,4) * t850 + Ifges(5,2) * t849 + Ifges(5,6) * t929 - pkin(4) * t795 - t896 * t865 + t930 * t867 - t980;
t888 = Ifges(4,5) * t909 + Ifges(4,6) * t908 + Ifges(4,3) * t933;
t762 = -mrSges(4,1) * t885 + mrSges(4,3) * t857 + Ifges(4,4) * t884 + Ifges(4,2) * t883 + Ifges(4,6) * t932 - pkin(3) * t954 + pkin(9) * t959 + t937 * t775 + t941 * t776 - t909 * t888 + t933 * t890;
t763 = mrSges(4,2) * t885 - mrSges(4,3) * t856 + Ifges(4,1) * t884 + Ifges(4,4) * t883 + Ifges(4,5) * t932 - pkin(9) * t783 + t775 * t941 - t776 * t937 + t888 * t908 - t889 * t933;
t905 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t939 + Ifges(3,6) * t943) * qJD(1);
t758 = -mrSges(3,1) * t910 + mrSges(3,3) * t899 + Ifges(3,4) * t916 + Ifges(3,2) * t917 + Ifges(3,6) * qJDD(2) - pkin(2) * t950 + pkin(8) * t960 + qJD(2) * t907 + t942 * t762 + t938 * t763 - t905 * t966;
t760 = mrSges(3,2) * t910 - mrSges(3,3) * t898 + Ifges(3,1) * t916 + Ifges(3,4) * t917 + Ifges(3,5) * qJDD(2) - pkin(8) * t774 - qJD(2) * t906 - t762 * t938 + t763 * t942 + t905 * t965;
t953 = mrSges(2,1) * t922 - mrSges(2,2) * t923 + Ifges(2,3) * qJDD(1) + pkin(1) * t947 + pkin(7) * t961 + t943 * t758 + t939 * t760;
t761 = mrSges(2,1) * g(3) + mrSges(2,3) * t923 + t945 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t768 - t981;
t756 = -mrSges(2,2) * g(3) - mrSges(2,3) * t922 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t945 - pkin(7) * t768 - t758 * t939 + t760 * t943;
t1 = [-m(1) * g(1) + t962; -m(1) * g(2) + t971; (-m(1) - m(2)) * g(3) + t768; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t971 + t944 * t756 - t940 * t761; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t962 + t940 * t756 + t944 * t761; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t953; t953; t981; -t948; -t952; t980; t807;];
tauJB  = t1;
