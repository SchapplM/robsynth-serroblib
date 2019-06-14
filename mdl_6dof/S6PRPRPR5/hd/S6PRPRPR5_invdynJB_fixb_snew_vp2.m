% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 23:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:57:33
% EndTime: 2019-05-04 22:57:41
% DurationCPUTime: 7.58s
% Computational Cost: add. (112749->301), mult. (252046->371), div. (0->0), fcn. (181538->12), ass. (0->143)
t939 = Ifges(5,1) + Ifges(6,2);
t932 = Ifges(5,4) + Ifges(6,6);
t931 = Ifges(5,5) - Ifges(6,4);
t938 = -Ifges(5,2) - Ifges(6,3);
t930 = Ifges(5,6) - Ifges(6,5);
t937 = Ifges(5,3) + Ifges(6,1);
t889 = qJD(2) ^ 2;
t878 = sin(pkin(10));
t881 = cos(pkin(10));
t862 = g(1) * t878 - g(2) * t881;
t863 = -g(1) * t881 - g(2) * t878;
t876 = -g(3) + qJDD(1);
t879 = sin(pkin(6));
t882 = cos(pkin(6));
t885 = sin(qJ(2));
t887 = cos(qJ(2));
t824 = -t885 * t863 + (t862 * t882 + t876 * t879) * t887;
t884 = sin(qJ(4));
t880 = cos(pkin(11));
t934 = cos(qJ(4));
t911 = t880 * t934;
t877 = sin(pkin(11));
t917 = qJD(2) * t877;
t855 = -qJD(2) * t911 + t884 * t917;
t847 = mrSges(6,1) * t855 - qJD(4) * mrSges(6,3);
t924 = t882 * t885;
t925 = t879 * t885;
t825 = t862 * t924 + t887 * t863 + t876 * t925;
t817 = -pkin(2) * t889 + qJDD(2) * qJ(3) + t825;
t850 = -t862 * t879 + t876 * t882;
t914 = qJD(2) * qJD(3);
t918 = t880 * t850 - 0.2e1 * t877 * t914;
t933 = pkin(3) * t880;
t796 = (-pkin(8) * qJDD(2) + t889 * t933 - t817) * t877 + t918;
t799 = t877 * t850 + (t817 + 0.2e1 * t914) * t880;
t912 = qJDD(2) * t880;
t874 = t880 ^ 2;
t926 = t874 * t889;
t797 = -pkin(3) * t926 + pkin(8) * t912 + t799;
t788 = t796 * t934 - t884 * t797;
t901 = t877 * t934 + t880 * t884;
t856 = t901 * qJD(2);
t830 = pkin(4) * t855 - qJ(5) * t856;
t888 = qJD(4) ^ 2;
t787 = -qJDD(4) * pkin(4) - t888 * qJ(5) + t856 * t830 + qJDD(5) - t788;
t916 = t855 * qJD(4);
t839 = qJDD(2) * t901 - t916;
t782 = (t855 * t856 - qJDD(4)) * pkin(9) + (t839 + t916) * pkin(5) + t787;
t913 = qJDD(2) * t877;
t915 = t856 * qJD(4);
t838 = -qJDD(2) * t911 + t884 * t913 + t915;
t849 = pkin(5) * t856 - qJD(4) * pkin(9);
t854 = t855 ^ 2;
t873 = t877 ^ 2;
t897 = qJDD(3) - t824;
t809 = (-pkin(2) - t933) * qJDD(2) + (-qJ(3) + (-t873 - t874) * pkin(8)) * t889 + t897;
t935 = -2 * qJD(5);
t890 = pkin(4) * t915 + t856 * t935 + (-t839 + t916) * qJ(5) + t809;
t785 = -t854 * pkin(5) - t856 * t849 + (pkin(4) + pkin(9)) * t838 + t890;
t883 = sin(qJ(6));
t886 = cos(qJ(6));
t780 = t782 * t886 - t785 * t883;
t840 = -qJD(4) * t883 + t855 * t886;
t808 = qJD(6) * t840 + qJDD(4) * t886 + t838 * t883;
t841 = qJD(4) * t886 + t855 * t883;
t810 = -mrSges(7,1) * t840 + mrSges(7,2) * t841;
t853 = qJD(6) + t856;
t815 = -mrSges(7,2) * t853 + mrSges(7,3) * t840;
t837 = qJDD(6) + t839;
t777 = m(7) * t780 + mrSges(7,1) * t837 - t808 * mrSges(7,3) - t810 * t841 + t815 * t853;
t781 = t782 * t883 + t785 * t886;
t807 = -qJD(6) * t841 - qJDD(4) * t883 + t838 * t886;
t816 = mrSges(7,1) * t853 - mrSges(7,3) * t841;
t778 = m(7) * t781 - mrSges(7,2) * t837 + t807 * mrSges(7,3) + t810 * t840 - t816 * t853;
t768 = t886 * t777 + t883 * t778;
t832 = -mrSges(6,2) * t855 - mrSges(6,3) * t856;
t898 = -m(6) * t787 - t839 * mrSges(6,1) - t856 * t832 - t768;
t766 = qJDD(4) * mrSges(6,2) + qJD(4) * t847 - t898;
t789 = t884 * t796 + t934 * t797;
t896 = -t888 * pkin(4) + qJDD(4) * qJ(5) - t855 * t830 + t789;
t784 = -t838 * pkin(5) - t854 * pkin(9) + ((2 * qJD(5)) + t849) * qJD(4) + t896;
t800 = Ifges(7,5) * t841 + Ifges(7,6) * t840 + Ifges(7,3) * t853;
t802 = Ifges(7,1) * t841 + Ifges(7,4) * t840 + Ifges(7,5) * t853;
t769 = -mrSges(7,1) * t784 + mrSges(7,3) * t781 + Ifges(7,4) * t808 + Ifges(7,2) * t807 + Ifges(7,6) * t837 - t800 * t841 + t802 * t853;
t801 = Ifges(7,4) * t841 + Ifges(7,2) * t840 + Ifges(7,6) * t853;
t770 = mrSges(7,2) * t784 - mrSges(7,3) * t780 + Ifges(7,1) * t808 + Ifges(7,4) * t807 + Ifges(7,5) * t837 + t800 * t840 - t801 * t853;
t786 = qJD(4) * t935 - t896;
t848 = mrSges(6,1) * t856 + qJD(4) * mrSges(6,2);
t899 = -m(7) * t784 + t807 * mrSges(7,1) - t808 * mrSges(7,2) + t840 * t815 - t841 * t816;
t894 = -m(6) * t786 + qJDD(4) * mrSges(6,3) + qJD(4) * t848 - t899;
t919 = t931 * qJD(4) - t932 * t855 + t939 * t856;
t920 = t930 * qJD(4) + t938 * t855 + t932 * t856;
t936 = t937 * qJDD(4) - t930 * t838 + t931 * t839 + t919 * t855 + t920 * t856 + mrSges(5,1) * t788 - mrSges(5,2) * t789 + mrSges(6,2) * t787 - mrSges(6,3) * t786 - pkin(4) * t766 - pkin(9) * t768 + qJ(5) * (-t838 * mrSges(6,1) - t855 * t832 + t894) - t883 * t769 + t886 * t770;
t928 = mrSges(4,2) * t877;
t814 = -qJDD(2) * pkin(2) - t889 * qJ(3) + t897;
t791 = t838 * pkin(4) + t890;
t922 = -t883 * t777 + t886 * t778;
t767 = m(6) * t791 - t838 * mrSges(6,2) - t839 * mrSges(6,3) - t855 * t847 - t856 * t848 + t922;
t845 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t855;
t846 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t856;
t893 = m(5) * t809 + t838 * mrSges(5,1) + t839 * mrSges(5,2) + t855 * t845 + t856 * t846 + t767;
t892 = -m(4) * t814 + mrSges(4,1) * t912 - t893 + (t873 * t889 + t926) * mrSges(4,3);
t764 = t892 + (mrSges(3,1) - t928) * qJDD(2) - t889 * mrSges(3,2) + m(3) * t824;
t927 = t764 * t887;
t831 = mrSges(5,1) * t855 + mrSges(5,2) * t856;
t761 = m(5) * t788 - t839 * mrSges(5,3) - t856 * t831 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t845 - t847) * qJD(4) + t898;
t773 = m(5) * t789 - qJDD(4) * mrSges(5,2) - qJD(4) * t846 + (-t831 - t832) * t855 + (-mrSges(5,3) - mrSges(6,1)) * t838 + t894;
t759 = t934 * t761 + t884 * t773;
t798 = -t817 * t877 + t918;
t902 = mrSges(4,3) * qJDD(2) + t889 * (-mrSges(4,1) * t880 + t928);
t757 = m(4) * t798 - t877 * t902 + t759;
t908 = -t884 * t761 + t934 * t773;
t758 = m(4) * t799 + t880 * t902 + t908;
t909 = -t757 * t877 + t880 * t758;
t749 = m(3) * t825 - mrSges(3,1) * t889 - qJDD(2) * mrSges(3,2) + t909;
t752 = t880 * t757 + t877 * t758;
t751 = m(3) * t850 + t752;
t740 = t749 * t924 - t751 * t879 + t882 * t927;
t738 = m(2) * t862 + t740;
t744 = t887 * t749 - t764 * t885;
t743 = m(2) * t863 + t744;
t923 = t881 * t738 + t878 * t743;
t921 = -t937 * qJD(4) + t930 * t855 - t931 * t856;
t739 = t749 * t925 + t882 * t751 + t879 * t927;
t910 = -t738 * t878 + t881 * t743;
t907 = m(2) * t876 + t739;
t906 = Ifges(4,1) * t877 + Ifges(4,4) * t880;
t905 = Ifges(4,4) * t877 + Ifges(4,2) * t880;
t904 = Ifges(4,5) * t877 + Ifges(4,6) * t880;
t745 = -mrSges(5,1) * t809 - mrSges(6,1) * t786 + mrSges(6,2) * t791 + mrSges(5,3) * t789 - pkin(4) * t767 - pkin(5) * t899 - pkin(9) * t922 + t919 * qJD(4) + t930 * qJDD(4) - t886 * t769 - t883 * t770 + t938 * t838 + t932 * t839 + t921 * t856;
t895 = mrSges(7,1) * t780 - mrSges(7,2) * t781 + Ifges(7,5) * t808 + Ifges(7,6) * t807 + Ifges(7,3) * t837 + t841 * t801 - t840 * t802;
t753 = mrSges(6,1) * t787 + mrSges(5,2) * t809 - mrSges(5,3) * t788 - mrSges(6,3) * t791 + pkin(5) * t768 - qJ(5) * t767 - t920 * qJD(4) + t931 * qJDD(4) - t932 * t838 + t939 * t839 + t921 * t855 + t895;
t861 = t904 * qJD(2);
t734 = -mrSges(4,1) * t814 + mrSges(4,3) * t799 - pkin(3) * t893 + pkin(8) * t908 + qJDD(2) * t905 + t745 * t934 + t884 * t753 - t861 * t917;
t736 = t880 * qJD(2) * t861 + mrSges(4,2) * t814 - mrSges(4,3) * t798 - pkin(8) * t759 + qJDD(2) * t906 - t884 * t745 + t753 * t934;
t733 = mrSges(3,2) * t850 - mrSges(3,3) * t824 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t889 - qJ(3) * t752 - t734 * t877 + t736 * t880;
t735 = (Ifges(3,6) - t904) * qJDD(2) - mrSges(3,1) * t850 + mrSges(3,3) * t825 - mrSges(4,1) * t798 + mrSges(4,2) * t799 - pkin(3) * t759 - pkin(2) * t752 + (-t877 * t905 + t880 * t906 + Ifges(3,5)) * t889 - t936;
t900 = pkin(7) * t744 + t733 * t885 + t735 * t887;
t765 = mrSges(4,2) * t913 - t892;
t732 = mrSges(3,1) * t824 - mrSges(3,2) * t825 + Ifges(3,3) * qJDD(2) - pkin(2) * t765 + qJ(3) * t909 + t880 * t734 + t877 * t736;
t731 = mrSges(2,2) * t876 - mrSges(2,3) * t862 + t887 * t733 - t885 * t735 + (-t739 * t879 - t740 * t882) * pkin(7);
t730 = -mrSges(2,1) * t876 + mrSges(2,3) * t863 - pkin(1) * t739 - t879 * t732 + t882 * t900;
t1 = [-m(1) * g(1) + t910; -m(1) * g(2) + t923; -m(1) * g(3) + t907; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t923 - t878 * t730 + t881 * t731; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t910 + t881 * t730 + t878 * t731; -mrSges(1,1) * g(2) + mrSges(2,1) * t862 + mrSges(1,2) * g(1) - mrSges(2,2) * t863 + pkin(1) * t740 + t882 * t732 + t879 * t900; t907; t732; t765; t936; t766; t895;];
tauJB  = t1;
