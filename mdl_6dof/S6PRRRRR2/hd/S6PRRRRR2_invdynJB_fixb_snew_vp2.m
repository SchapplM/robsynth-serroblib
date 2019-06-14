% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:43:02
% EndTime: 2019-05-05 10:43:26
% DurationCPUTime: 23.20s
% Computational Cost: add. (422352->344), mult. (825754->441), div. (0->0), fcn. (607582->14), ass. (0->148)
t889 = sin(pkin(12));
t891 = cos(pkin(12));
t875 = g(1) * t889 - g(2) * t891;
t876 = -g(1) * t891 - g(2) * t889;
t888 = -g(3) + qJDD(1);
t890 = sin(pkin(6));
t892 = cos(pkin(6));
t897 = sin(qJ(2));
t902 = cos(qJ(2));
t844 = -t897 * t876 + (t875 * t892 + t888 * t890) * t902;
t927 = t892 * t897;
t928 = t890 * t897;
t845 = t875 * t927 + t902 * t876 + t888 * t928;
t903 = qJD(2) ^ 2;
t840 = -pkin(2) * t903 + qJDD(2) * pkin(8) + t845;
t856 = -t875 * t890 + t888 * t892;
t896 = sin(qJ(3));
t901 = cos(qJ(3));
t820 = -t896 * t840 + t901 * t856;
t923 = qJD(2) * qJD(3);
t922 = t901 * t923;
t873 = qJDD(2) * t896 + t922;
t808 = (-t873 + t922) * pkin(9) + (t896 * t901 * t903 + qJDD(3)) * pkin(3) + t820;
t821 = t901 * t840 + t896 * t856;
t874 = qJDD(2) * t901 - t896 * t923;
t925 = qJD(2) * t896;
t880 = qJD(3) * pkin(3) - pkin(9) * t925;
t887 = t901 ^ 2;
t809 = -pkin(3) * t887 * t903 + pkin(9) * t874 - qJD(3) * t880 + t821;
t895 = sin(qJ(4));
t900 = cos(qJ(4));
t798 = t895 * t808 + t900 * t809;
t865 = (t895 * t901 + t896 * t900) * qJD(2);
t834 = -t865 * qJD(4) - t895 * t873 + t874 * t900;
t924 = qJD(2) * t901;
t864 = -t895 * t925 + t900 * t924;
t847 = -mrSges(5,1) * t864 + mrSges(5,2) * t865;
t886 = qJD(3) + qJD(4);
t855 = mrSges(5,1) * t886 - mrSges(5,3) * t865;
t885 = qJDD(3) + qJDD(4);
t848 = -pkin(4) * t864 - pkin(10) * t865;
t884 = t886 ^ 2;
t787 = -pkin(4) * t884 + pkin(10) * t885 + t848 * t864 + t798;
t910 = -qJDD(2) * pkin(2) - t844;
t815 = -t874 * pkin(3) + t880 * t925 + (-pkin(9) * t887 - pkin(8)) * t903 + t910;
t835 = qJD(4) * t864 + t873 * t900 + t874 * t895;
t795 = (-t864 * t886 - t835) * pkin(10) + (t865 * t886 - t834) * pkin(4) + t815;
t894 = sin(qJ(5));
t899 = cos(qJ(5));
t782 = -t894 * t787 + t899 * t795;
t850 = -t865 * t894 + t886 * t899;
t812 = qJD(5) * t850 + t835 * t899 + t885 * t894;
t832 = qJDD(5) - t834;
t851 = t865 * t899 + t886 * t894;
t859 = qJD(5) - t864;
t780 = (t850 * t859 - t812) * pkin(11) + (t850 * t851 + t832) * pkin(5) + t782;
t783 = t899 * t787 + t894 * t795;
t811 = -qJD(5) * t851 - t835 * t894 + t885 * t899;
t838 = pkin(5) * t859 - pkin(11) * t851;
t849 = t850 ^ 2;
t781 = -pkin(5) * t849 + pkin(11) * t811 - t838 * t859 + t783;
t893 = sin(qJ(6));
t898 = cos(qJ(6));
t778 = t780 * t898 - t781 * t893;
t822 = t850 * t898 - t851 * t893;
t792 = qJD(6) * t822 + t811 * t893 + t812 * t898;
t823 = t850 * t893 + t851 * t898;
t804 = -mrSges(7,1) * t822 + mrSges(7,2) * t823;
t857 = qJD(6) + t859;
t813 = -mrSges(7,2) * t857 + mrSges(7,3) * t822;
t827 = qJDD(6) + t832;
t773 = m(7) * t778 + mrSges(7,1) * t827 - t792 * mrSges(7,3) - t804 * t823 + t813 * t857;
t779 = t780 * t893 + t781 * t898;
t791 = -qJD(6) * t823 + t811 * t898 - t812 * t893;
t814 = mrSges(7,1) * t857 - mrSges(7,3) * t823;
t774 = m(7) * t779 - mrSges(7,2) * t827 + t791 * mrSges(7,3) + t804 * t822 - t814 * t857;
t765 = t898 * t773 + t893 * t774;
t824 = -mrSges(6,1) * t850 + mrSges(6,2) * t851;
t836 = -mrSges(6,2) * t859 + mrSges(6,3) * t850;
t763 = m(6) * t782 + mrSges(6,1) * t832 - mrSges(6,3) * t812 - t824 * t851 + t836 * t859 + t765;
t837 = mrSges(6,1) * t859 - mrSges(6,3) * t851;
t917 = -t773 * t893 + t898 * t774;
t764 = m(6) * t783 - mrSges(6,2) * t832 + mrSges(6,3) * t811 + t824 * t850 - t837 * t859 + t917;
t918 = -t763 * t894 + t899 * t764;
t756 = m(5) * t798 - mrSges(5,2) * t885 + mrSges(5,3) * t834 + t847 * t864 - t855 * t886 + t918;
t797 = t808 * t900 - t895 * t809;
t854 = -mrSges(5,2) * t886 + mrSges(5,3) * t864;
t786 = -pkin(4) * t885 - pkin(10) * t884 + t865 * t848 - t797;
t784 = -pkin(5) * t811 - pkin(11) * t849 + t838 * t851 + t786;
t912 = m(7) * t784 - t791 * mrSges(7,1) + t792 * mrSges(7,2) - t822 * t813 + t814 * t823;
t907 = -m(6) * t786 + t811 * mrSges(6,1) - mrSges(6,2) * t812 + t850 * t836 - t837 * t851 - t912;
t769 = m(5) * t797 + mrSges(5,1) * t885 - mrSges(5,3) * t835 - t847 * t865 + t854 * t886 + t907;
t746 = t895 * t756 + t900 * t769;
t862 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t896 + Ifges(4,2) * t901) * qJD(2);
t863 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t896 + Ifges(4,4) * t901) * qJD(2);
t800 = Ifges(7,5) * t823 + Ifges(7,6) * t822 + Ifges(7,3) * t857;
t802 = Ifges(7,1) * t823 + Ifges(7,4) * t822 + Ifges(7,5) * t857;
t766 = -mrSges(7,1) * t784 + mrSges(7,3) * t779 + Ifges(7,4) * t792 + Ifges(7,2) * t791 + Ifges(7,6) * t827 - t800 * t823 + t802 * t857;
t801 = Ifges(7,4) * t823 + Ifges(7,2) * t822 + Ifges(7,6) * t857;
t767 = mrSges(7,2) * t784 - mrSges(7,3) * t778 + Ifges(7,1) * t792 + Ifges(7,4) * t791 + Ifges(7,5) * t827 + t800 * t822 - t801 * t857;
t816 = Ifges(6,5) * t851 + Ifges(6,6) * t850 + Ifges(6,3) * t859;
t818 = Ifges(6,1) * t851 + Ifges(6,4) * t850 + Ifges(6,5) * t859;
t748 = -mrSges(6,1) * t786 + mrSges(6,3) * t783 + Ifges(6,4) * t812 + Ifges(6,2) * t811 + Ifges(6,6) * t832 - pkin(5) * t912 + pkin(11) * t917 + t898 * t766 + t893 * t767 - t851 * t816 + t859 * t818;
t817 = Ifges(6,4) * t851 + Ifges(6,2) * t850 + Ifges(6,6) * t859;
t750 = mrSges(6,2) * t786 - mrSges(6,3) * t782 + Ifges(6,1) * t812 + Ifges(6,4) * t811 + Ifges(6,5) * t832 - pkin(11) * t765 - t766 * t893 + t767 * t898 + t816 * t850 - t817 * t859;
t842 = Ifges(5,4) * t865 + Ifges(5,2) * t864 + Ifges(5,6) * t886;
t843 = Ifges(5,1) * t865 + Ifges(5,4) * t864 + Ifges(5,5) * t886;
t908 = -mrSges(5,1) * t797 + mrSges(5,2) * t798 - Ifges(5,5) * t835 - Ifges(5,6) * t834 - Ifges(5,3) * t885 - pkin(4) * t907 - pkin(10) * t918 - t899 * t748 - t894 * t750 - t865 * t842 + t864 * t843;
t930 = mrSges(4,1) * t820 - mrSges(4,2) * t821 + Ifges(4,5) * t873 + Ifges(4,6) * t874 + Ifges(4,3) * qJDD(3) + pkin(3) * t746 + (t862 * t896 - t863 * t901) * qJD(2) - t908;
t839 = -t903 * pkin(8) + t910;
t877 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t925;
t878 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t924;
t758 = t899 * t763 + t894 * t764;
t909 = m(5) * t815 - t834 * mrSges(5,1) + mrSges(5,2) * t835 - t864 * t854 + t855 * t865 + t758;
t906 = -m(4) * t839 + t874 * mrSges(4,1) - mrSges(4,2) * t873 - t877 * t925 + t878 * t924 - t909;
t753 = m(3) * t844 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t903 + t906;
t929 = t753 * t902;
t872 = (-mrSges(4,1) * t901 + mrSges(4,2) * t896) * qJD(2);
t744 = m(4) * t820 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t873 + qJD(3) * t878 - t872 * t925 + t746;
t919 = t900 * t756 - t769 * t895;
t745 = m(4) * t821 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t874 - qJD(3) * t877 + t872 * t924 + t919;
t920 = -t744 * t896 + t901 * t745;
t736 = m(3) * t845 - mrSges(3,1) * t903 - qJDD(2) * mrSges(3,2) + t920;
t739 = t901 * t744 + t896 * t745;
t738 = m(3) * t856 + t739;
t727 = t736 * t927 - t738 * t890 + t892 * t929;
t725 = m(2) * t875 + t727;
t731 = t902 * t736 - t753 * t897;
t730 = m(2) * t876 + t731;
t926 = t891 * t725 + t889 * t730;
t726 = t736 * t928 + t892 * t738 + t890 * t929;
t921 = -t725 * t889 + t891 * t730;
t916 = m(2) * t888 + t726;
t841 = Ifges(5,5) * t865 + Ifges(5,6) * t864 + Ifges(5,3) * t886;
t732 = mrSges(5,2) * t815 - mrSges(5,3) * t797 + Ifges(5,1) * t835 + Ifges(5,4) * t834 + Ifges(5,5) * t885 - pkin(10) * t758 - t748 * t894 + t750 * t899 + t841 * t864 - t842 * t886;
t911 = -mrSges(7,1) * t778 + mrSges(7,2) * t779 - Ifges(7,5) * t792 - Ifges(7,6) * t791 - Ifges(7,3) * t827 - t823 * t801 + t822 * t802;
t904 = mrSges(6,1) * t782 - mrSges(6,2) * t783 + Ifges(6,5) * t812 + Ifges(6,6) * t811 + Ifges(6,3) * t832 + pkin(5) * t765 + t851 * t817 - t850 * t818 - t911;
t740 = -mrSges(5,1) * t815 + mrSges(5,3) * t798 + Ifges(5,4) * t835 + Ifges(5,2) * t834 + Ifges(5,6) * t885 - pkin(4) * t758 - t865 * t841 + t886 * t843 - t904;
t861 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t896 + Ifges(4,6) * t901) * qJD(2);
t721 = -mrSges(4,1) * t839 + mrSges(4,3) * t821 + Ifges(4,4) * t873 + Ifges(4,2) * t874 + Ifges(4,6) * qJDD(3) - pkin(3) * t909 + pkin(9) * t919 + qJD(3) * t863 + t895 * t732 + t900 * t740 - t861 * t925;
t723 = mrSges(4,2) * t839 - mrSges(4,3) * t820 + Ifges(4,1) * t873 + Ifges(4,4) * t874 + Ifges(4,5) * qJDD(3) - pkin(9) * t746 - qJD(3) * t862 + t732 * t900 - t740 * t895 + t861 * t924;
t720 = mrSges(3,2) * t856 - mrSges(3,3) * t844 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t903 - pkin(8) * t739 - t721 * t896 + t723 * t901;
t722 = -mrSges(3,1) * t856 + mrSges(3,3) * t845 + t903 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t739 - t930;
t913 = pkin(7) * t731 + t720 * t897 + t722 * t902;
t719 = mrSges(3,1) * t844 - mrSges(3,2) * t845 + Ifges(3,3) * qJDD(2) + pkin(2) * t906 + pkin(8) * t920 + t901 * t721 + t896 * t723;
t718 = mrSges(2,2) * t888 - mrSges(2,3) * t875 + t902 * t720 - t897 * t722 + (-t726 * t890 - t727 * t892) * pkin(7);
t717 = -mrSges(2,1) * t888 + mrSges(2,3) * t876 - pkin(1) * t726 - t890 * t719 + t892 * t913;
t1 = [-m(1) * g(1) + t921; -m(1) * g(2) + t926; -m(1) * g(3) + t916; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t926 - t889 * t717 + t891 * t718; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t921 + t891 * t717 + t889 * t718; -mrSges(1,1) * g(2) + mrSges(2,1) * t875 + mrSges(1,2) * g(1) - mrSges(2,2) * t876 + pkin(1) * t727 + t892 * t719 + t890 * t913; t916; t719; t930; -t908; t904; -t911;];
tauJB  = t1;
