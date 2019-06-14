% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPPR4
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
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:09:05
% EndTime: 2019-05-05 03:09:13
% DurationCPUTime: 8.13s
% Computational Cost: add. (121115->320), mult. (252456->395), div. (0->0), fcn. (168870->12), ass. (0->139)
t936 = -2 * qJD(4);
t935 = Ifges(5,1) + Ifges(6,1);
t926 = Ifges(5,4) - Ifges(6,5);
t934 = Ifges(5,5) + Ifges(6,4);
t933 = Ifges(5,2) + Ifges(6,3);
t932 = Ifges(6,2) + Ifges(5,3);
t931 = Ifges(5,6) - Ifges(6,6);
t879 = sin(pkin(10));
t881 = cos(pkin(10));
t865 = t879 * g(1) - t881 * g(2);
t866 = -t881 * g(1) - t879 * g(2);
t877 = -g(3) + qJDD(1);
t888 = cos(qJ(2));
t882 = cos(pkin(6));
t885 = sin(qJ(2));
t919 = t882 * t885;
t880 = sin(pkin(6));
t920 = t880 * t885;
t817 = t865 * t919 + t888 * t866 + t877 * t920;
t890 = qJD(2) ^ 2;
t811 = -t890 * pkin(2) + qJDD(2) * pkin(8) + t817;
t841 = -t880 * t865 + t882 * t877;
t884 = sin(qJ(3));
t887 = cos(qJ(3));
t806 = t887 * t811 + t884 * t841;
t861 = (-pkin(3) * t887 - qJ(4) * t884) * qJD(2);
t889 = qJD(3) ^ 2;
t911 = t887 * qJD(2);
t793 = -t889 * pkin(3) + qJDD(3) * qJ(4) + t861 * t911 + t806;
t816 = -t885 * t866 + (t865 * t882 + t877 * t880) * t888;
t810 = -qJDD(2) * pkin(2) - t890 * pkin(8) - t816;
t910 = qJD(2) * qJD(3);
t908 = t887 * t910;
t863 = t884 * qJDD(2) + t908;
t872 = t884 * t910;
t864 = t887 * qJDD(2) - t872;
t798 = (-t863 - t908) * qJ(4) + (-t864 + t872) * pkin(3) + t810;
t878 = sin(pkin(11));
t913 = qJD(2) * t884;
t923 = cos(pkin(11));
t854 = t878 * qJD(3) + t923 * t913;
t786 = -t878 * t793 + t923 * t798 + t854 * t936;
t839 = t878 * qJDD(3) + t923 * t863;
t805 = -t884 * t811 + t887 * t841;
t895 = qJDD(3) * pkin(3) + t889 * qJ(4) - t861 * t913 - qJDD(4) + t805;
t853 = -t923 * qJD(3) + t878 * t913;
t909 = t853 * t911;
t930 = -(t839 + t909) * qJ(5) - t895;
t828 = t853 * pkin(4) - t854 * qJ(5);
t921 = t887 ^ 2 * t890;
t784 = t864 * pkin(4) - qJ(5) * t921 + t854 * t828 + qJDD(5) - t786;
t779 = (-t839 + t909) * pkin(9) + (t853 * t854 + t864) * pkin(5) + t784;
t787 = t923 * t793 + t878 * t798 + t853 * t936;
t928 = -2 * qJD(5);
t782 = -pkin(4) * t921 - t864 * qJ(5) - t853 * t828 + t911 * t928 + t787;
t838 = -t923 * qJDD(3) + t878 * t863;
t840 = pkin(5) * t911 - t854 * pkin(9);
t851 = t853 ^ 2;
t780 = -t851 * pkin(5) + t838 * pkin(9) - t840 * t911 + t782;
t883 = sin(qJ(6));
t886 = cos(qJ(6));
t778 = t883 * t779 + t886 * t780;
t785 = -t851 * pkin(9) + (-pkin(4) - pkin(5)) * t838 + (pkin(4) * t911 + (2 * qJD(5)) + t840) * t854 - t930;
t827 = t883 * t853 + t886 * t854;
t799 = -t827 * qJD(6) + t886 * t838 - t883 * t839;
t826 = t886 * t853 - t883 * t854;
t800 = t826 * qJD(6) + t883 * t838 + t886 * t839;
t870 = qJD(6) + t911;
t801 = Ifges(7,5) * t827 + Ifges(7,6) * t826 + Ifges(7,3) * t870;
t803 = Ifges(7,1) * t827 + Ifges(7,4) * t826 + Ifges(7,5) * t870;
t858 = qJDD(6) + t864;
t767 = -mrSges(7,1) * t785 + mrSges(7,3) * t778 + Ifges(7,4) * t800 + Ifges(7,2) * t799 + Ifges(7,6) * t858 - t827 * t801 + t870 * t803;
t777 = t886 * t779 - t883 * t780;
t802 = Ifges(7,4) * t827 + Ifges(7,2) * t826 + Ifges(7,6) * t870;
t768 = mrSges(7,2) * t785 - mrSges(7,3) * t777 + Ifges(7,1) * t800 + Ifges(7,4) * t799 + Ifges(7,5) * t858 + t826 * t801 - t870 * t802;
t788 = t854 * t928 + (-t854 * t911 + t838) * pkin(4) + t930;
t835 = -t853 * mrSges(6,2) - mrSges(6,3) * t911;
t837 = mrSges(6,1) * t911 + t854 * mrSges(6,2);
t812 = -t870 * mrSges(7,2) + t826 * mrSges(7,3);
t813 = t870 * mrSges(7,1) - t827 * mrSges(7,3);
t896 = -m(7) * t785 + t799 * mrSges(7,1) - t800 * mrSges(7,2) + t826 * t812 - t827 * t813;
t775 = m(6) * t788 + t838 * mrSges(6,1) - t839 * mrSges(6,3) + t853 * t835 - t854 * t837 + t896;
t807 = -t826 * mrSges(7,1) + t827 * mrSges(7,2);
t773 = m(7) * t777 + t858 * mrSges(7,1) - t800 * mrSges(7,3) - t827 * t807 + t870 * t812;
t774 = m(7) * t778 - t858 * mrSges(7,2) + t799 * mrSges(7,3) + t826 * t807 - t870 * t813;
t905 = -t883 * t773 + t886 * t774;
t915 = t853 * t926 - t854 * t935 + t911 * t934;
t917 = t853 * t931 - t854 * t934 + t911 * t932;
t745 = mrSges(5,1) * t895 - mrSges(6,1) * t788 + mrSges(6,2) * t782 + mrSges(5,3) * t787 - pkin(4) * t775 - pkin(5) * t896 - pkin(9) * t905 - t886 * t767 - t883 * t768 - t838 * t933 + t926 * t839 + t917 * t854 - t864 * t931 + t915 * t911;
t766 = t886 * t773 + t883 * t774;
t916 = t853 * t933 - t854 * t926 + t911 * t931;
t750 = -mrSges(5,2) * t895 + mrSges(6,2) * t784 - mrSges(5,3) * t786 - mrSges(6,3) * t788 - pkin(9) * t766 - qJ(5) * t775 - t883 * t767 + t886 * t768 - t926 * t838 + t839 * t935 + t917 * t853 - t864 * t934 - t916 * t911;
t836 = -mrSges(5,1) * t911 - t854 * mrSges(5,3);
t899 = m(6) * t782 - t864 * mrSges(6,3) + t905;
t829 = t853 * mrSges(6,1) - t854 * mrSges(6,3);
t914 = -t853 * mrSges(5,1) - t854 * mrSges(5,2) - t829;
t927 = -mrSges(5,3) - mrSges(6,2);
t763 = m(5) * t787 + t864 * mrSges(5,2) + t914 * t853 + t927 * t838 + (t836 - t837) * t911 + t899;
t893 = -m(6) * t784 - t864 * mrSges(6,1) - t835 * t911 - t766;
t898 = mrSges(5,2) * t911 - t853 * mrSges(5,3);
t764 = m(5) * t786 - t864 * mrSges(5,1) + t927 * t839 + t914 * t854 - t898 * t911 + t893;
t761 = t923 * t763 - t878 * t764;
t771 = -m(5) * t895 + t838 * mrSges(5,1) + t839 * mrSges(5,2) + t854 * t836 + t853 * t898 + t775;
t847 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t884 + Ifges(4,2) * t887) * qJD(2);
t848 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t884 + Ifges(4,4) * t887) * qJD(2);
t929 = mrSges(4,1) * t805 - mrSges(4,2) * t806 + Ifges(4,5) * t863 + Ifges(4,6) * t864 + Ifges(4,3) * qJDD(3) - pkin(3) * t771 + qJ(4) * t761 + (t847 * t884 - t848 * t887) * qJD(2) + t923 * t745 + t878 * t750;
t760 = t878 * t763 + t923 * t764;
t867 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t913;
t868 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t911;
t892 = -m(4) * t810 + t864 * mrSges(4,1) - t863 * mrSges(4,2) - t867 * t913 + t868 * t911 - t760;
t756 = m(3) * t816 + qJDD(2) * mrSges(3,1) - t890 * mrSges(3,2) + t892;
t922 = t756 * t888;
t862 = (-mrSges(4,1) * t887 + mrSges(4,2) * t884) * qJD(2);
t759 = m(4) * t806 - qJDD(3) * mrSges(4,2) + t864 * mrSges(4,3) - qJD(3) * t867 + t862 * t911 + t761;
t770 = m(4) * t805 + qJDD(3) * mrSges(4,1) - t863 * mrSges(4,3) + qJD(3) * t868 - t862 * t913 - t771;
t906 = t887 * t759 - t884 * t770;
t749 = m(3) * t817 - t890 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t906;
t753 = t884 * t759 + t887 * t770;
t752 = m(3) * t841 + t753;
t739 = t749 * t919 - t880 * t752 + t882 * t922;
t737 = m(2) * t865 + t739;
t744 = t888 * t749 - t885 * t756;
t743 = m(2) * t866 + t744;
t918 = t881 * t737 + t879 * t743;
t738 = t749 * t920 + t882 * t752 + t880 * t922;
t907 = -t879 * t737 + t881 * t743;
t904 = m(2) * t877 + t738;
t846 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t884 + Ifges(4,6) * t887) * qJD(2);
t735 = mrSges(4,2) * t810 - mrSges(4,3) * t805 + Ifges(4,1) * t863 + Ifges(4,4) * t864 + Ifges(4,5) * qJDD(3) - qJ(4) * t760 - qJD(3) * t847 - t878 * t745 + t923 * t750 + t846 * t911;
t765 = t839 * mrSges(6,2) + t854 * t829 - t893;
t894 = mrSges(7,1) * t777 - mrSges(7,2) * t778 + Ifges(7,5) * t800 + Ifges(7,6) * t799 + Ifges(7,3) * t858 + t827 * t802 - t826 * t803;
t740 = Ifges(4,6) * qJDD(3) + t894 + (Ifges(4,2) + t932) * t864 + t916 * t854 + (qJ(5) * t829 + t915) * t853 - t934 * t839 + (qJ(5) * mrSges(6,2) + t931) * t838 + Ifges(4,4) * t863 + qJD(3) * t848 + mrSges(4,3) * t806 - mrSges(4,1) * t810 + mrSges(6,1) * t784 - mrSges(5,1) * t786 + mrSges(5,2) * t787 - mrSges(6,3) * t782 + pkin(4) * t765 + pkin(5) * t766 - t846 * t913 - pkin(3) * t760 - qJ(5) * (-t837 * t911 + t899);
t733 = mrSges(3,2) * t841 - mrSges(3,3) * t816 + Ifges(3,5) * qJDD(2) - t890 * Ifges(3,6) - pkin(8) * t753 + t887 * t735 - t884 * t740;
t734 = -mrSges(3,1) * t841 + mrSges(3,3) * t817 + t890 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t753 - t929;
t897 = pkin(7) * t744 + t733 * t885 + t734 * t888;
t732 = mrSges(3,1) * t816 - mrSges(3,2) * t817 + Ifges(3,3) * qJDD(2) + pkin(2) * t892 + pkin(8) * t906 + t884 * t735 + t887 * t740;
t731 = mrSges(2,2) * t877 - mrSges(2,3) * t865 + t888 * t733 - t885 * t734 + (-t738 * t880 - t739 * t882) * pkin(7);
t730 = -mrSges(2,1) * t877 + mrSges(2,3) * t866 - pkin(1) * t738 - t880 * t732 + t897 * t882;
t1 = [-m(1) * g(1) + t907; -m(1) * g(2) + t918; -m(1) * g(3) + t904; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t918 - t879 * t730 + t881 * t731; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t907 + t881 * t730 + t879 * t731; -mrSges(1,1) * g(2) + mrSges(2,1) * t865 + mrSges(1,2) * g(1) - mrSges(2,2) * t866 + pkin(1) * t739 + t882 * t732 + t897 * t880; t904; t732; t929; t771; t765; t894;];
tauJB  = t1;
