% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:16:15
% EndTime: 2019-05-06 13:16:52
% DurationCPUTime: 34.44s
% Computational Cost: add. (533576->387), mult. (1241047->490), div. (0->0), fcn. (907295->12), ass. (0->151)
t943 = sin(qJ(2));
t947 = cos(qJ(2));
t965 = qJD(1) * qJD(2);
t925 = qJDD(1) * t943 + t947 * t965;
t944 = sin(qJ(1));
t948 = cos(qJ(1));
t932 = -g(1) * t948 - g(2) * t944;
t950 = qJD(1) ^ 2;
t920 = -pkin(1) * t950 + qJDD(1) * pkin(7) + t932;
t970 = t943 * t920;
t973 = pkin(2) * t950;
t881 = qJDD(2) * pkin(2) - t925 * qJ(3) - t970 + (qJ(3) * t965 + t943 * t973 - g(3)) * t947;
t905 = -g(3) * t943 + t947 * t920;
t926 = qJDD(1) * t947 - t943 * t965;
t968 = qJD(1) * t943;
t928 = qJD(2) * pkin(2) - qJ(3) * t968;
t936 = t947 ^ 2;
t882 = qJ(3) * t926 - qJD(2) * t928 - t936 * t973 + t905;
t938 = sin(pkin(10));
t940 = cos(pkin(10));
t914 = (t938 * t947 + t940 * t943) * qJD(1);
t857 = -0.2e1 * qJD(3) * t914 + t881 * t940 - t938 * t882;
t967 = qJD(1) * t947;
t913 = -t938 * t968 + t940 * t967;
t858 = 0.2e1 * qJD(3) * t913 + t938 * t881 + t940 * t882;
t895 = -pkin(3) * t913 - pkin(8) * t914;
t949 = qJD(2) ^ 2;
t844 = -pkin(3) * t949 + qJDD(2) * pkin(8) + t895 * t913 + t858;
t931 = t944 * g(1) - t948 * g(2);
t956 = -qJDD(1) * pkin(1) - t931;
t885 = -t926 * pkin(2) + qJDD(3) + t928 * t968 + (-qJ(3) * t936 - pkin(7)) * t950 + t956;
t899 = -t938 * t925 + t926 * t940;
t900 = t925 * t940 + t926 * t938;
t847 = (-qJD(2) * t913 - t900) * pkin(8) + (qJD(2) * t914 - t899) * pkin(3) + t885;
t942 = sin(qJ(4));
t946 = cos(qJ(4));
t834 = -t942 * t844 + t946 * t847;
t902 = qJD(2) * t946 - t914 * t942;
t872 = qJD(4) * t902 + qJDD(2) * t942 + t900 * t946;
t898 = qJDD(4) - t899;
t903 = qJD(2) * t942 + t914 * t946;
t912 = qJD(4) - t913;
t823 = (t902 * t912 - t872) * qJ(5) + (t902 * t903 + t898) * pkin(4) + t834;
t835 = t946 * t844 + t942 * t847;
t871 = -qJD(4) * t903 + qJDD(2) * t946 - t900 * t942;
t887 = pkin(4) * t912 - qJ(5) * t903;
t901 = t902 ^ 2;
t825 = -pkin(4) * t901 + qJ(5) * t871 - t887 * t912 + t835;
t937 = sin(pkin(11));
t939 = cos(pkin(11));
t880 = t902 * t937 + t903 * t939;
t817 = -0.2e1 * qJD(5) * t880 + t939 * t823 - t937 * t825;
t852 = t871 * t937 + t872 * t939;
t879 = t902 * t939 - t903 * t937;
t815 = (t879 * t912 - t852) * pkin(9) + (t879 * t880 + t898) * pkin(5) + t817;
t818 = 0.2e1 * qJD(5) * t879 + t937 * t823 + t939 * t825;
t851 = t871 * t939 - t872 * t937;
t864 = pkin(5) * t912 - pkin(9) * t880;
t876 = t879 ^ 2;
t816 = -pkin(5) * t876 + pkin(9) * t851 - t864 * t912 + t818;
t941 = sin(qJ(6));
t945 = cos(qJ(6));
t813 = t815 * t945 - t816 * t941;
t859 = t879 * t945 - t880 * t941;
t831 = qJD(6) * t859 + t851 * t941 + t852 * t945;
t860 = t879 * t941 + t880 * t945;
t841 = -mrSges(7,1) * t859 + mrSges(7,2) * t860;
t908 = qJD(6) + t912;
t848 = -mrSges(7,2) * t908 + mrSges(7,3) * t859;
t896 = qJDD(6) + t898;
t806 = m(7) * t813 + mrSges(7,1) * t896 - mrSges(7,3) * t831 - t841 * t860 + t848 * t908;
t814 = t815 * t941 + t816 * t945;
t830 = -qJD(6) * t860 + t851 * t945 - t852 * t941;
t849 = mrSges(7,1) * t908 - mrSges(7,3) * t860;
t807 = m(7) * t814 - mrSges(7,2) * t896 + mrSges(7,3) * t830 + t841 * t859 - t849 * t908;
t800 = t945 * t806 + t941 * t807;
t861 = -mrSges(6,1) * t879 + mrSges(6,2) * t880;
t862 = -mrSges(6,2) * t912 + mrSges(6,3) * t879;
t798 = m(6) * t817 + mrSges(6,1) * t898 - mrSges(6,3) * t852 - t861 * t880 + t862 * t912 + t800;
t863 = mrSges(6,1) * t912 - mrSges(6,3) * t880;
t960 = -t806 * t941 + t945 * t807;
t799 = m(6) * t818 - mrSges(6,2) * t898 + mrSges(6,3) * t851 + t861 * t879 - t863 * t912 + t960;
t794 = t939 * t798 + t937 * t799;
t854 = Ifges(6,4) * t880 + Ifges(6,2) * t879 + Ifges(6,6) * t912;
t855 = Ifges(6,1) * t880 + Ifges(6,4) * t879 + Ifges(6,5) * t912;
t866 = Ifges(5,4) * t903 + Ifges(5,2) * t902 + Ifges(5,6) * t912;
t867 = Ifges(5,1) * t903 + Ifges(5,4) * t902 + Ifges(5,5) * t912;
t837 = Ifges(7,4) * t860 + Ifges(7,2) * t859 + Ifges(7,6) * t908;
t838 = Ifges(7,1) * t860 + Ifges(7,4) * t859 + Ifges(7,5) * t908;
t954 = -mrSges(7,1) * t813 + mrSges(7,2) * t814 - Ifges(7,5) * t831 - Ifges(7,6) * t830 - Ifges(7,3) * t896 - t860 * t837 + t859 * t838;
t975 = mrSges(5,1) * t834 + mrSges(6,1) * t817 - mrSges(5,2) * t835 - mrSges(6,2) * t818 + Ifges(5,5) * t872 + Ifges(6,5) * t852 + Ifges(5,6) * t871 + Ifges(6,6) * t851 + pkin(4) * t794 + pkin(5) * t800 + t880 * t854 - t879 * t855 + t903 * t866 - t902 * t867 + (Ifges(5,3) + Ifges(6,3)) * t898 - t954;
t843 = -qJDD(2) * pkin(3) - pkin(8) * t949 + t914 * t895 - t857;
t833 = -pkin(4) * t871 - qJ(5) * t901 + t903 * t887 + qJDD(5) + t843;
t820 = -pkin(5) * t851 - pkin(9) * t876 + t864 * t880 + t833;
t836 = Ifges(7,5) * t860 + Ifges(7,6) * t859 + Ifges(7,3) * t908;
t801 = -mrSges(7,1) * t820 + mrSges(7,3) * t814 + Ifges(7,4) * t831 + Ifges(7,2) * t830 + Ifges(7,6) * t896 - t836 * t860 + t838 * t908;
t802 = mrSges(7,2) * t820 - mrSges(7,3) * t813 + Ifges(7,1) * t831 + Ifges(7,4) * t830 + Ifges(7,5) * t896 + t836 * t859 - t837 * t908;
t853 = Ifges(6,5) * t880 + Ifges(6,6) * t879 + Ifges(6,3) * t912;
t958 = m(7) * t820 - t830 * mrSges(7,1) + t831 * mrSges(7,2) - t859 * t848 + t860 * t849;
t787 = -mrSges(6,1) * t833 + mrSges(6,3) * t818 + Ifges(6,4) * t852 + Ifges(6,2) * t851 + Ifges(6,6) * t898 - pkin(5) * t958 + pkin(9) * t960 + t945 * t801 + t941 * t802 - t880 * t853 + t912 * t855;
t788 = mrSges(6,2) * t833 - mrSges(6,3) * t817 + Ifges(6,1) * t852 + Ifges(6,4) * t851 + Ifges(6,5) * t898 - pkin(9) * t800 - t801 * t941 + t802 * t945 + t853 * t879 - t854 * t912;
t811 = m(6) * t833 - t851 * mrSges(6,1) + mrSges(6,2) * t852 - t879 * t862 + t863 * t880 + t958;
t865 = Ifges(5,5) * t903 + Ifges(5,6) * t902 + Ifges(5,3) * t912;
t961 = -t798 * t937 + t939 * t799;
t770 = -mrSges(5,1) * t843 + mrSges(5,3) * t835 + Ifges(5,4) * t872 + Ifges(5,2) * t871 + Ifges(5,6) * t898 - pkin(4) * t811 + qJ(5) * t961 + t939 * t787 + t937 * t788 - t903 * t865 + t912 * t867;
t771 = mrSges(5,2) * t843 - mrSges(5,3) * t834 + Ifges(5,1) * t872 + Ifges(5,4) * t871 + Ifges(5,5) * t898 - qJ(5) * t794 - t787 * t937 + t788 * t939 + t865 * t902 - t866 * t912;
t883 = -mrSges(5,1) * t902 + mrSges(5,2) * t903;
t886 = -mrSges(5,2) * t912 + mrSges(5,3) * t902;
t792 = m(5) * t834 + mrSges(5,1) * t898 - mrSges(5,3) * t872 - t883 * t903 + t886 * t912 + t794;
t888 = mrSges(5,1) * t912 - mrSges(5,3) * t903;
t793 = m(5) * t835 - mrSges(5,2) * t898 + mrSges(5,3) * t871 + t883 * t902 - t888 * t912 + t961;
t786 = -t792 * t942 + t946 * t793;
t893 = -mrSges(4,1) * t913 + mrSges(4,2) * t914;
t907 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t914;
t783 = m(4) * t858 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t899 - qJD(2) * t907 + t893 * t913 + t786;
t810 = -m(5) * t843 + t871 * mrSges(5,1) - mrSges(5,2) * t872 + t902 * t886 - t888 * t903 - t811;
t906 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t913;
t809 = m(4) * t857 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t900 + qJD(2) * t906 - t893 * t914 + t810;
t777 = t938 * t783 + t940 * t809;
t890 = Ifges(4,4) * t914 + Ifges(4,2) * t913 + Ifges(4,6) * qJD(2);
t891 = Ifges(4,1) * t914 + Ifges(4,4) * t913 + Ifges(4,5) * qJD(2);
t904 = -t947 * g(3) - t970;
t916 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t943 + Ifges(3,2) * t947) * qJD(1);
t917 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t943 + Ifges(3,4) * t947) * qJD(1);
t974 = (t916 * t943 - t917 * t947) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t904 + mrSges(4,1) * t857 - mrSges(3,2) * t905 - mrSges(4,2) * t858 + Ifges(3,5) * t925 + Ifges(4,5) * t900 + Ifges(3,6) * t926 + Ifges(4,6) * t899 + pkin(2) * t777 + pkin(3) * t810 + pkin(8) * t786 + t946 * t770 + t942 * t771 + t914 * t890 - t913 * t891;
t924 = (-mrSges(3,1) * t947 + mrSges(3,2) * t943) * qJD(1);
t930 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t967;
t775 = m(3) * t904 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t925 + qJD(2) * t930 - t924 * t968 + t777;
t929 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t968;
t962 = t940 * t783 - t809 * t938;
t776 = m(3) * t905 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t926 - qJD(2) * t929 + t924 * t967 + t962;
t963 = -t775 * t943 + t947 * t776;
t766 = m(2) * t932 - mrSges(2,1) * t950 - qJDD(1) * mrSges(2,2) + t963;
t785 = t946 * t792 + t942 * t793;
t784 = m(4) * t885 - t899 * mrSges(4,1) + mrSges(4,2) * t900 - t913 * t906 + t907 * t914 + t785;
t919 = -t950 * pkin(7) + t956;
t953 = -m(3) * t919 + t926 * mrSges(3,1) - mrSges(3,2) * t925 - t929 * t968 + t930 * t967 - t784;
t779 = m(2) * t931 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t950 + t953;
t969 = t944 * t766 + t948 * t779;
t768 = t947 * t775 + t943 * t776;
t964 = t948 * t766 - t779 * t944;
t889 = Ifges(4,5) * t914 + Ifges(4,6) * t913 + Ifges(4,3) * qJD(2);
t763 = mrSges(4,2) * t885 - mrSges(4,3) * t857 + Ifges(4,1) * t900 + Ifges(4,4) * t899 + Ifges(4,5) * qJDD(2) - pkin(8) * t785 - qJD(2) * t890 - t770 * t942 + t771 * t946 + t889 * t913;
t769 = -mrSges(4,1) * t885 + mrSges(4,3) * t858 + Ifges(4,4) * t900 + Ifges(4,2) * t899 + Ifges(4,6) * qJDD(2) - pkin(3) * t785 + qJD(2) * t891 - t914 * t889 - t975;
t915 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t943 + Ifges(3,6) * t947) * qJD(1);
t759 = -mrSges(3,1) * t919 + mrSges(3,3) * t905 + Ifges(3,4) * t925 + Ifges(3,2) * t926 + Ifges(3,6) * qJDD(2) - pkin(2) * t784 + qJ(3) * t962 + qJD(2) * t917 + t938 * t763 + t940 * t769 - t915 * t968;
t762 = mrSges(3,2) * t919 - mrSges(3,3) * t904 + Ifges(3,1) * t925 + Ifges(3,4) * t926 + Ifges(3,5) * qJDD(2) - qJ(3) * t777 - qJD(2) * t916 + t763 * t940 - t769 * t938 + t915 * t967;
t955 = mrSges(2,1) * t931 - mrSges(2,2) * t932 + Ifges(2,3) * qJDD(1) + pkin(1) * t953 + pkin(7) * t963 + t947 * t759 + t943 * t762;
t760 = mrSges(2,1) * g(3) + mrSges(2,3) * t932 + t950 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t768 - t974;
t757 = -mrSges(2,2) * g(3) - mrSges(2,3) * t931 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t950 - pkin(7) * t768 - t759 * t943 + t762 * t947;
t1 = [-m(1) * g(1) + t964; -m(1) * g(2) + t969; (-m(1) - m(2)) * g(3) + t768; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t969 + t948 * t757 - t944 * t760; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t964 + t944 * t757 + t948 * t760; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t955; t955; t974; t784; t975; t811; -t954;];
tauJB  = t1;
