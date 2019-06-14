% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:17:06
% EndTime: 2019-05-06 12:17:21
% DurationCPUTime: 14.44s
% Computational Cost: add. (213360->363), mult. (492544->448), div. (0->0), fcn. (349473->10), ass. (0->144)
t970 = -2 * qJD(3);
t969 = Ifges(6,1) + Ifges(7,1);
t961 = Ifges(6,4) - Ifges(7,5);
t960 = Ifges(6,5) + Ifges(7,4);
t968 = -Ifges(6,2) - Ifges(7,3);
t967 = -Ifges(7,2) - Ifges(6,3);
t959 = Ifges(6,6) - Ifges(7,6);
t924 = sin(qJ(2));
t927 = cos(qJ(2));
t947 = qJD(1) * qJD(2);
t909 = t924 * qJDD(1) + t927 * t947;
t925 = sin(qJ(1));
t928 = cos(qJ(1));
t916 = -t928 * g(1) - t925 * g(2);
t930 = qJD(1) ^ 2;
t904 = -t930 * pkin(1) + qJDD(1) * pkin(7) + t916;
t956 = t924 * t904;
t963 = pkin(2) * t930;
t865 = qJDD(2) * pkin(2) - t909 * qJ(3) - t956 + (qJ(3) * t947 + t924 * t963 - g(3)) * t927;
t890 = -t924 * g(3) + t927 * t904;
t910 = t927 * qJDD(1) - t924 * t947;
t950 = qJD(1) * t924;
t912 = qJD(2) * pkin(2) - qJ(3) * t950;
t919 = t927 ^ 2;
t866 = t910 * qJ(3) - qJD(2) * t912 - t919 * t963 + t890;
t921 = sin(pkin(9));
t922 = cos(pkin(9));
t899 = (t921 * t927 + t922 * t924) * qJD(1);
t839 = t922 * t865 - t921 * t866 + t899 * t970;
t898 = (t921 * t924 - t922 * t927) * qJD(1);
t878 = t898 * pkin(3) - t899 * pkin(8);
t929 = qJD(2) ^ 2;
t820 = -qJDD(2) * pkin(3) - t929 * pkin(8) + t899 * t878 - t839;
t884 = t922 * t909 + t921 * t910;
t923 = sin(qJ(4));
t926 = cos(qJ(4));
t888 = t923 * qJD(2) + t926 * t899;
t857 = -t888 * qJD(4) + t926 * qJDD(2) - t923 * t884;
t897 = qJD(4) + t898;
t871 = t897 * pkin(4) - t888 * qJ(5);
t887 = t926 * qJD(2) - t923 * t899;
t886 = t887 ^ 2;
t816 = -t857 * pkin(4) - t886 * qJ(5) + t888 * t871 + qJDD(5) + t820;
t858 = t887 * qJD(4) + t923 * qJDD(2) + t926 * t884;
t920 = sin(pkin(10));
t957 = cos(pkin(10));
t828 = -t957 * t857 + t920 * t858;
t829 = t920 * t857 + t957 * t858;
t863 = -t957 * t887 + t920 * t888;
t864 = t920 * t887 + t957 * t888;
t811 = -0.2e1 * qJD(6) * t864 + (t863 * t897 - t829) * qJ(6) + (t864 * t897 + t828) * pkin(5) + t816;
t846 = -t863 * mrSges(7,2) + t897 * mrSges(7,3);
t849 = -t897 * mrSges(7,1) + t864 * mrSges(7,2);
t804 = m(7) * t811 + t828 * mrSges(7,1) - t829 * mrSges(7,3) + t863 * t846 - t864 * t849;
t840 = t921 * t865 + t922 * t866 + t898 * t970;
t821 = -t929 * pkin(3) + qJDD(2) * pkin(8) - t898 * t878 + t840;
t915 = t925 * g(1) - t928 * g(2);
t936 = -qJDD(1) * pkin(1) - t915;
t869 = -t910 * pkin(2) + qJDD(3) + t912 * t950 + (-qJ(3) * t919 - pkin(7)) * t930 + t936;
t883 = -t921 * t909 + t922 * t910;
t824 = (qJD(2) * t898 - t884) * pkin(8) + (qJD(2) * t899 - t883) * pkin(3) + t869;
t817 = -t923 * t821 + t926 * t824;
t882 = qJDD(4) - t883;
t813 = (t887 * t897 - t858) * qJ(5) + (t887 * t888 + t882) * pkin(4) + t817;
t818 = t926 * t821 + t923 * t824;
t815 = -t886 * pkin(4) + t857 * qJ(5) - t897 * t871 + t818;
t964 = -2 * qJD(5);
t809 = t920 * t813 + t957 * t815 + t863 * t964;
t841 = t863 * pkin(5) - t864 * qJ(6);
t896 = t897 ^ 2;
t806 = -t896 * pkin(5) + t882 * qJ(6) + 0.2e1 * qJD(6) * t897 - t863 * t841 + t809;
t952 = -t863 * t961 + t864 * t969 + t897 * t960;
t954 = t863 * t959 - t864 * t960 + t897 * t967;
t791 = -mrSges(6,1) * t816 - mrSges(7,1) * t811 + mrSges(7,2) * t806 + mrSges(6,3) * t809 - pkin(5) * t804 + t828 * t968 + t961 * t829 + t954 * t864 + t959 * t882 + t952 * t897;
t935 = t957 * t813 - t920 * t815;
t807 = -t882 * pkin(5) - t896 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t841) * t864 - t935;
t808 = t864 * t964 + t935;
t953 = t863 * t968 + t864 * t961 + t897 * t959;
t792 = mrSges(6,2) * t816 + mrSges(7,2) * t807 - mrSges(6,3) * t808 - mrSges(7,3) * t811 - qJ(6) * t804 - t961 * t828 + t829 * t969 + t954 * t863 + t960 * t882 - t953 * t897;
t847 = -t897 * mrSges(6,2) - t863 * mrSges(6,3);
t848 = t897 * mrSges(6,1) - t864 * mrSges(6,3);
t801 = m(6) * t816 + t828 * mrSges(6,1) + t829 * mrSges(6,2) + t863 * t847 + t864 * t848 + t804;
t850 = Ifges(5,5) * t888 + Ifges(5,6) * t887 + Ifges(5,3) * t897;
t852 = Ifges(5,1) * t888 + Ifges(5,4) * t887 + Ifges(5,5) * t897;
t945 = m(7) * t806 + t882 * mrSges(7,3) + t897 * t849;
t842 = t863 * mrSges(7,1) - t864 * mrSges(7,3);
t951 = -t863 * mrSges(6,1) - t864 * mrSges(6,2) - t842;
t962 = -mrSges(6,3) - mrSges(7,2);
t795 = m(6) * t809 - t882 * mrSges(6,2) + t962 * t828 - t897 * t848 + t951 * t863 + t945;
t939 = -m(7) * t807 + t882 * mrSges(7,1) + t897 * t846;
t797 = m(6) * t808 + t882 * mrSges(6,1) + t962 * t829 + t897 * t847 + t951 * t864 + t939;
t941 = t957 * t795 - t920 * t797;
t768 = -mrSges(5,1) * t820 + mrSges(5,3) * t818 + Ifges(5,4) * t858 + Ifges(5,2) * t857 + Ifges(5,6) * t882 - pkin(4) * t801 + qJ(5) * t941 + t957 * t791 + t920 * t792 - t888 * t850 + t897 * t852;
t790 = t920 * t795 + t957 * t797;
t851 = Ifges(5,4) * t888 + Ifges(5,2) * t887 + Ifges(5,6) * t897;
t769 = mrSges(5,2) * t820 - mrSges(5,3) * t817 + Ifges(5,1) * t858 + Ifges(5,4) * t857 + Ifges(5,5) * t882 - qJ(5) * t790 - t920 * t791 + t957 * t792 + t887 * t850 - t897 * t851;
t867 = -t887 * mrSges(5,1) + t888 * mrSges(5,2);
t870 = -t897 * mrSges(5,2) + t887 * mrSges(5,3);
t788 = m(5) * t817 + t882 * mrSges(5,1) - t858 * mrSges(5,3) - t888 * t867 + t897 * t870 + t790;
t872 = t897 * mrSges(5,1) - t888 * mrSges(5,3);
t789 = m(5) * t818 - t882 * mrSges(5,2) + t857 * mrSges(5,3) + t887 * t867 - t897 * t872 + t941;
t784 = -t923 * t788 + t926 * t789;
t877 = t898 * mrSges(4,1) + t899 * mrSges(4,2);
t892 = qJD(2) * mrSges(4,1) - t899 * mrSges(4,3);
t781 = m(4) * t840 - qJDD(2) * mrSges(4,2) + t883 * mrSges(4,3) - qJD(2) * t892 - t898 * t877 + t784;
t800 = -m(5) * t820 + t857 * mrSges(5,1) - t858 * mrSges(5,2) + t887 * t870 - t888 * t872 - t801;
t891 = -qJD(2) * mrSges(4,2) - t898 * mrSges(4,3);
t799 = m(4) * t839 + qJDD(2) * mrSges(4,1) - t884 * mrSges(4,3) + qJD(2) * t891 - t899 * t877 + t800;
t775 = t921 * t781 + t922 * t799;
t874 = Ifges(4,4) * t899 - Ifges(4,2) * t898 + Ifges(4,6) * qJD(2);
t875 = Ifges(4,1) * t899 - Ifges(4,4) * t898 + Ifges(4,5) * qJD(2);
t889 = -t927 * g(3) - t956;
t901 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t924 + Ifges(3,2) * t927) * qJD(1);
t902 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t924 + Ifges(3,4) * t927) * qJD(1);
t966 = (t924 * t901 - t927 * t902) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t889 + mrSges(4,1) * t839 - mrSges(3,2) * t890 - mrSges(4,2) * t840 + Ifges(3,5) * t909 + Ifges(4,5) * t884 + Ifges(3,6) * t910 + Ifges(4,6) * t883 + pkin(2) * t775 + pkin(3) * t800 + pkin(8) * t784 + t926 * t768 + t923 * t769 + t899 * t874 + t898 * t875;
t803 = t829 * mrSges(7,2) + t864 * t842 - t939;
t965 = -t959 * t828 + t960 * t829 + t952 * t863 + t953 * t864 + (Ifges(5,3) - t967) * t882 + mrSges(5,1) * t817 + mrSges(6,1) * t808 - mrSges(7,1) * t807 - mrSges(5,2) * t818 - mrSges(6,2) * t809 + mrSges(7,3) * t806 + Ifges(5,5) * t858 + Ifges(5,6) * t857 + pkin(4) * t790 - pkin(5) * t803 + qJ(6) * (-t828 * mrSges(7,2) - t863 * t842 + t945) + t888 * t851 - t887 * t852;
t908 = (-mrSges(3,1) * t927 + mrSges(3,2) * t924) * qJD(1);
t949 = qJD(1) * t927;
t914 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t949;
t773 = m(3) * t889 + qJDD(2) * mrSges(3,1) - t909 * mrSges(3,3) + qJD(2) * t914 - t908 * t950 + t775;
t913 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t950;
t942 = t922 * t781 - t921 * t799;
t774 = m(3) * t890 - qJDD(2) * mrSges(3,2) + t910 * mrSges(3,3) - qJD(2) * t913 + t908 * t949 + t942;
t943 = -t924 * t773 + t927 * t774;
t764 = m(2) * t916 - t930 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t943;
t783 = t926 * t788 + t923 * t789;
t782 = m(4) * t869 - t883 * mrSges(4,1) + t884 * mrSges(4,2) + t898 * t891 + t899 * t892 + t783;
t903 = -t930 * pkin(7) + t936;
t933 = -m(3) * t903 + t910 * mrSges(3,1) - t909 * mrSges(3,2) - t913 * t950 + t914 * t949 - t782;
t777 = m(2) * t915 + qJDD(1) * mrSges(2,1) - t930 * mrSges(2,2) + t933;
t955 = t925 * t764 + t928 * t777;
t766 = t927 * t773 + t924 * t774;
t944 = t928 * t764 - t925 * t777;
t873 = Ifges(4,5) * t899 - Ifges(4,6) * t898 + Ifges(4,3) * qJD(2);
t761 = mrSges(4,2) * t869 - mrSges(4,3) * t839 + Ifges(4,1) * t884 + Ifges(4,4) * t883 + Ifges(4,5) * qJDD(2) - pkin(8) * t783 - qJD(2) * t874 - t923 * t768 + t926 * t769 - t898 * t873;
t767 = -mrSges(4,1) * t869 + mrSges(4,3) * t840 + Ifges(4,4) * t884 + Ifges(4,2) * t883 + Ifges(4,6) * qJDD(2) - pkin(3) * t783 + qJD(2) * t875 - t899 * t873 - t965;
t900 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t924 + Ifges(3,6) * t927) * qJD(1);
t757 = -mrSges(3,1) * t903 + mrSges(3,3) * t890 + Ifges(3,4) * t909 + Ifges(3,2) * t910 + Ifges(3,6) * qJDD(2) - pkin(2) * t782 + qJ(3) * t942 + qJD(2) * t902 + t921 * t761 + t922 * t767 - t900 * t950;
t760 = mrSges(3,2) * t903 - mrSges(3,3) * t889 + Ifges(3,1) * t909 + Ifges(3,4) * t910 + Ifges(3,5) * qJDD(2) - qJ(3) * t775 - qJD(2) * t901 + t922 * t761 - t921 * t767 + t900 * t949;
t934 = mrSges(2,1) * t915 - mrSges(2,2) * t916 + Ifges(2,3) * qJDD(1) + pkin(1) * t933 + pkin(7) * t943 + t927 * t757 + t924 * t760;
t758 = mrSges(2,1) * g(3) + mrSges(2,3) * t916 + t930 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t766 - t966;
t755 = -mrSges(2,2) * g(3) - mrSges(2,3) * t915 + Ifges(2,5) * qJDD(1) - t930 * Ifges(2,6) - pkin(7) * t766 - t924 * t757 + t927 * t760;
t1 = [-m(1) * g(1) + t944; -m(1) * g(2) + t955; (-m(1) - m(2)) * g(3) + t766; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t955 + t928 * t755 - t925 * t758; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t944 + t925 * t755 + t928 * t758; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t934; t934; t966; t782; t965; t801; t803;];
tauJB  = t1;
