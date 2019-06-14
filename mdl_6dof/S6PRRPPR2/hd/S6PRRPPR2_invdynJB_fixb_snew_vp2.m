% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-05-05 02:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:44:43
% EndTime: 2019-05-05 02:44:51
% DurationCPUTime: 8.28s
% Computational Cost: add. (123181->324), mult. (267108->400), div. (0->0), fcn. (182184->12), ass. (0->141)
t938 = -2 * qJD(4);
t937 = Ifges(5,1) + Ifges(6,2);
t936 = -Ifges(6,1) - Ifges(5,3);
t932 = Ifges(5,4) + Ifges(6,6);
t931 = Ifges(5,5) - Ifges(6,4);
t935 = Ifges(5,2) + Ifges(6,3);
t930 = Ifges(5,6) - Ifges(6,5);
t884 = sin(pkin(10));
t886 = cos(pkin(10));
t872 = g(1) * t884 - g(2) * t886;
t873 = -g(1) * t886 - g(2) * t884;
t882 = -g(3) + qJDD(1);
t893 = cos(qJ(2));
t887 = cos(pkin(6));
t890 = sin(qJ(2));
t926 = t887 * t890;
t885 = sin(pkin(6));
t927 = t885 * t890;
t825 = t872 * t926 + t893 * t873 + t882 * t927;
t895 = qJD(2) ^ 2;
t817 = -pkin(2) * t895 + qJDD(2) * pkin(8) + t825;
t850 = -t872 * t885 + t882 * t887;
t889 = sin(qJ(3));
t892 = cos(qJ(3));
t799 = -t817 * t889 + t892 * t850;
t915 = qJD(2) * qJD(3);
t913 = t892 * t915;
t870 = qJDD(2) * t889 + t913;
t795 = (-t870 + t913) * qJ(4) + (t889 * t892 * t895 + qJDD(3)) * pkin(3) + t799;
t800 = t892 * t817 + t889 * t850;
t871 = qJDD(2) * t892 - t889 * t915;
t919 = qJD(2) * t889;
t874 = qJD(3) * pkin(3) - qJ(4) * t919;
t881 = t892 ^ 2;
t796 = -pkin(3) * t881 * t895 + qJ(4) * t871 - qJD(3) * t874 + t800;
t883 = sin(pkin(11));
t929 = cos(pkin(11));
t858 = (t883 * t892 + t889 * t929) * qJD(2);
t788 = t795 * t929 - t883 * t796 + t858 * t938;
t824 = -t890 * t873 + (t872 * t887 + t882 * t885) * t893;
t918 = qJD(2) * t892;
t857 = t883 * t919 - t918 * t929;
t830 = mrSges(5,1) * t857 + mrSges(5,2) * t858;
t839 = t870 * t929 + t883 * t871;
t845 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t857;
t847 = mrSges(6,1) * t857 - qJD(3) * mrSges(6,3);
t829 = pkin(4) * t857 - qJ(5) * t858;
t894 = qJD(3) ^ 2;
t786 = -qJDD(3) * pkin(4) - t894 * qJ(5) + t858 * t829 + qJDD(5) - t788;
t917 = qJD(3) * t857;
t782 = (t857 * t858 - qJDD(3)) * pkin(9) + (t839 + t917) * pkin(5) + t786;
t838 = t870 * t883 - t871 * t929;
t849 = pkin(5) * t858 - qJD(3) * pkin(9);
t856 = t857 ^ 2;
t900 = -qJDD(2) * pkin(2) - t824;
t798 = -pkin(3) * t871 + qJDD(4) + t874 * t919 + (-qJ(4) * t881 - pkin(8)) * t895 + t900;
t933 = -2 * qJD(5);
t897 = (-t839 + t917) * qJ(5) + t798 + (pkin(4) * qJD(3) + t933) * t858;
t787 = -pkin(5) * t856 - t849 * t858 + (pkin(4) + pkin(9)) * t838 + t897;
t888 = sin(qJ(6));
t891 = cos(qJ(6));
t780 = t782 * t891 - t787 * t888;
t840 = -qJD(3) * t888 + t857 * t891;
t809 = qJD(6) * t840 + qJDD(3) * t891 + t838 * t888;
t841 = qJD(3) * t891 + t857 * t888;
t810 = -mrSges(7,1) * t840 + mrSges(7,2) * t841;
t855 = qJD(6) + t858;
t814 = -mrSges(7,2) * t855 + mrSges(7,3) * t840;
t837 = qJDD(6) + t839;
t777 = m(7) * t780 + mrSges(7,1) * t837 - mrSges(7,3) * t809 - t810 * t841 + t814 * t855;
t781 = t782 * t888 + t787 * t891;
t808 = -qJD(6) * t841 - qJDD(3) * t888 + t838 * t891;
t815 = mrSges(7,1) * t855 - mrSges(7,3) * t841;
t778 = m(7) * t781 - mrSges(7,2) * t837 + t808 * mrSges(7,3) + t810 * t840 - t815 * t855;
t768 = t777 * t891 + t778 * t888;
t831 = -mrSges(6,2) * t857 - mrSges(6,3) * t858;
t902 = -m(6) * t786 - t839 * mrSges(6,1) - t858 * t831 - t768;
t761 = m(5) * t788 - mrSges(5,3) * t839 - t830 * t858 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t845 - t847) * qJD(3) + t902;
t853 = t857 * t938;
t923 = t883 * t795 + t929 * t796;
t789 = t853 + t923;
t846 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t858;
t905 = pkin(4) * t894 - qJDD(3) * qJ(5) - t923;
t785 = qJD(3) * t933 + ((2 * qJD(4)) + t829) * t857 + t905;
t848 = mrSges(6,1) * t858 + qJD(3) * mrSges(6,2);
t784 = -pkin(5) * t838 - pkin(9) * t856 - t829 * t857 + t853 + ((2 * qJD(5)) + t849) * qJD(3) - t905;
t903 = -m(7) * t784 + t808 * mrSges(7,1) - t809 * mrSges(7,2) + t814 * t840 - t841 * t815;
t899 = -m(6) * t785 + qJDD(3) * mrSges(6,3) + qJD(3) * t848 - t903;
t773 = m(5) * t789 - qJDD(3) * mrSges(5,2) - qJD(3) * t846 + (-t830 - t831) * t857 + (-mrSges(5,3) - mrSges(6,1)) * t838 + t899;
t759 = t929 * t761 + t883 * t773;
t766 = qJDD(3) * mrSges(6,2) + qJD(3) * t847 - t902;
t801 = Ifges(7,5) * t841 + Ifges(7,6) * t840 + Ifges(7,3) * t855;
t803 = Ifges(7,1) * t841 + Ifges(7,4) * t840 + Ifges(7,5) * t855;
t771 = -mrSges(7,1) * t784 + mrSges(7,3) * t781 + Ifges(7,4) * t809 + Ifges(7,2) * t808 + Ifges(7,6) * t837 - t801 * t841 + t803 * t855;
t802 = Ifges(7,4) * t841 + Ifges(7,2) * t840 + Ifges(7,6) * t855;
t772 = mrSges(7,2) * t784 - mrSges(7,3) * t780 + Ifges(7,1) * t809 + Ifges(7,4) * t808 + Ifges(7,5) * t837 + t801 * t840 - t802 * t855;
t861 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t889 + Ifges(4,2) * t892) * qJD(2);
t862 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t889 + Ifges(4,4) * t892) * qJD(2);
t920 = t931 * qJD(3) - t932 * t857 + t937 * t858;
t921 = -t930 * qJD(3) + t935 * t857 - t932 * t858;
t934 = (t889 * t861 - t892 * t862) * qJD(2) + (Ifges(4,3) - t936) * qJDD(3) - t838 * t930 + t839 * t931 + t857 * t920 - t858 * t921 + mrSges(4,1) * t799 + mrSges(5,1) * t788 - mrSges(4,2) * t800 - mrSges(5,2) * t789 + mrSges(6,2) * t786 - mrSges(6,3) * t785 + Ifges(4,5) * t870 + Ifges(4,6) * t871 + pkin(3) * t759 - pkin(4) * t766 - pkin(9) * t768 + qJ(5) * (-mrSges(6,1) * t838 - t831 * t857 + t899) - t888 * t771 + t891 * t772;
t791 = pkin(4) * t838 + t897;
t924 = -t888 * t777 + t891 * t778;
t767 = m(6) * t791 - t838 * mrSges(6,2) - t839 * mrSges(6,3) - t857 * t847 - t858 * t848 + t924;
t765 = m(5) * t798 + t838 * mrSges(5,1) + mrSges(5,2) * t839 + t857 * t845 + t846 * t858 + t767;
t816 = -pkin(8) * t895 + t900;
t875 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t919;
t876 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t918;
t898 = -m(4) * t816 + t871 * mrSges(4,1) - mrSges(4,2) * t870 - t875 * t919 + t876 * t918 - t765;
t764 = m(3) * t824 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t895 + t898;
t928 = t764 * t893;
t869 = (-mrSges(4,1) * t892 + mrSges(4,2) * t889) * qJD(2);
t757 = m(4) * t799 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t870 + qJD(3) * t876 - t869 * t919 + t759;
t910 = -t761 * t883 + t929 * t773;
t758 = m(4) * t800 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t871 - qJD(3) * t875 + t869 * t918 + t910;
t911 = -t757 * t889 + t892 * t758;
t749 = m(3) * t825 - mrSges(3,1) * t895 - qJDD(2) * mrSges(3,2) + t911;
t752 = t892 * t757 + t889 * t758;
t751 = m(3) * t850 + t752;
t740 = t749 * t926 - t751 * t885 + t887 * t928;
t738 = m(2) * t872 + t740;
t744 = t893 * t749 - t764 * t890;
t743 = m(2) * t873 + t744;
t925 = t886 * t738 + t884 * t743;
t922 = t936 * qJD(3) + t930 * t857 - t931 * t858;
t739 = t749 * t927 + t887 * t751 + t885 * t928;
t912 = -t738 * t884 + t886 * t743;
t909 = m(2) * t882 + t739;
t745 = -mrSges(5,1) * t798 - mrSges(6,1) * t785 + mrSges(6,2) * t791 + mrSges(5,3) * t789 - pkin(4) * t767 - pkin(5) * t903 - pkin(9) * t924 + t920 * qJD(3) + t930 * qJDD(3) - t891 * t771 - t888 * t772 - t935 * t838 + t932 * t839 + t922 * t858;
t901 = mrSges(7,1) * t780 - mrSges(7,2) * t781 + Ifges(7,5) * t809 + Ifges(7,6) * t808 + Ifges(7,3) * t837 + t841 * t802 - t840 * t803;
t753 = mrSges(6,1) * t786 + mrSges(5,2) * t798 - mrSges(5,3) * t788 - mrSges(6,3) * t791 + pkin(5) * t768 - qJ(5) * t767 + t921 * qJD(3) + t931 * qJDD(3) - t932 * t838 + t937 * t839 + t922 * t857 + t901;
t860 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t889 + Ifges(4,6) * t892) * qJD(2);
t734 = -mrSges(4,1) * t816 + mrSges(4,3) * t800 + Ifges(4,4) * t870 + Ifges(4,2) * t871 + Ifges(4,6) * qJDD(3) - pkin(3) * t765 + qJ(4) * t910 + qJD(3) * t862 + t745 * t929 + t883 * t753 - t860 * t919;
t736 = mrSges(4,2) * t816 - mrSges(4,3) * t799 + Ifges(4,1) * t870 + Ifges(4,4) * t871 + Ifges(4,5) * qJDD(3) - qJ(4) * t759 - qJD(3) * t861 - t883 * t745 + t753 * t929 + t860 * t918;
t733 = mrSges(3,2) * t850 - mrSges(3,3) * t824 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t895 - pkin(8) * t752 - t734 * t889 + t736 * t892;
t735 = -mrSges(3,1) * t850 + mrSges(3,3) * t825 + t895 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t752 - t934;
t904 = pkin(7) * t744 + t733 * t890 + t735 * t893;
t732 = mrSges(3,1) * t824 - mrSges(3,2) * t825 + Ifges(3,3) * qJDD(2) + pkin(2) * t898 + pkin(8) * t911 + t892 * t734 + t889 * t736;
t731 = mrSges(2,2) * t882 - mrSges(2,3) * t872 + t733 * t893 - t735 * t890 + (-t739 * t885 - t740 * t887) * pkin(7);
t730 = -mrSges(2,1) * t882 + mrSges(2,3) * t873 - pkin(1) * t739 - t732 * t885 + t887 * t904;
t1 = [-m(1) * g(1) + t912; -m(1) * g(2) + t925; -m(1) * g(3) + t909; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t925 - t884 * t730 + t886 * t731; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t912 + t886 * t730 + t884 * t731; -mrSges(1,1) * g(2) + mrSges(2,1) * t872 + mrSges(1,2) * g(1) - mrSges(2,2) * t873 + pkin(1) * t740 + t732 * t887 + t885 * t904; t909; t732; t934; t765; t766; t901;];
tauJB  = t1;
