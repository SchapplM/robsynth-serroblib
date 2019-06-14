% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-05-05 14:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:37:55
% EndTime: 2019-05-05 14:37:59
% DurationCPUTime: 4.33s
% Computational Cost: add. (40616->295), mult. (92649->343), div. (0->0), fcn. (60646->8), ass. (0->133)
t908 = Ifges(5,4) + Ifges(6,6);
t923 = -Ifges(5,2) - Ifges(6,3);
t920 = Ifges(5,6) - Ifges(6,5);
t922 = Ifges(5,1) + Ifges(6,2);
t921 = Ifges(5,5) - Ifges(6,4);
t919 = Ifges(5,3) + Ifges(6,1);
t858 = sin(pkin(9));
t859 = cos(pkin(9));
t863 = cos(qJ(4));
t913 = sin(qJ(4));
t881 = t858 * t863 + t913 * t859;
t830 = t881 * qJD(1);
t892 = t858 * t913;
t897 = t859 * qJD(1);
t831 = -qJD(1) * t892 + t863 * t897;
t918 = t920 * qJD(4) + t923 * t830 + t908 * t831;
t861 = sin(qJ(1));
t864 = cos(qJ(1));
t835 = t861 * g(1) - t864 * g(2);
t866 = qJD(1) ^ 2;
t879 = -t866 * qJ(2) + qJDD(2) - t835;
t905 = -pkin(1) - qJ(3);
t917 = -0.2e1 * qJD(1) * qJD(3) + t905 * qJDD(1) + t879;
t836 = -t864 * g(1) - t861 * g(2);
t916 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t836;
t823 = t866 * pkin(1) - t916;
t915 = -m(3) * t823 + t866 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t914 = -2 * qJD(5);
t911 = pkin(3) * t866;
t910 = mrSges(2,1) - mrSges(3,2);
t909 = mrSges(5,1) - mrSges(6,2);
t907 = Ifges(2,5) - Ifges(3,4);
t906 = -Ifges(2,6) + Ifges(3,5);
t806 = t858 * g(3) + t859 * t917;
t780 = (-pkin(7) * qJDD(1) - t858 * t911) * t859 + t806;
t807 = -t859 * g(3) + t858 * t917;
t847 = t858 ^ 2;
t895 = qJDD(1) * t858;
t781 = -pkin(7) * t895 - t847 * t911 + t807;
t760 = t863 * t780 - t913 * t781;
t796 = t830 * mrSges(5,1) + t831 * mrSges(5,2);
t894 = qJDD(1) * t859;
t899 = t830 * qJD(4);
t809 = -qJDD(1) * t892 + t863 * t894 - t899;
t898 = t831 * qJD(4);
t808 = t881 * qJDD(1) + t898;
t821 = t831 * pkin(5) - qJD(4) * pkin(8);
t829 = t830 ^ 2;
t878 = qJDD(3) + t916;
t900 = -t859 ^ 2 - t847;
t785 = pkin(3) * t895 + (t900 * pkin(7) + t905) * t866 + t878;
t868 = pkin(4) * t898 + t831 * t914 + (-t809 + t899) * qJ(5) + t785;
t750 = -t829 * pkin(5) - t831 * t821 + (pkin(4) + pkin(8)) * t808 + t868;
t795 = t830 * pkin(4) - t831 * qJ(5);
t865 = qJD(4) ^ 2;
t758 = -qJDD(4) * pkin(4) - t865 * qJ(5) + t831 * t795 + qJDD(5) - t760;
t751 = (t830 * t831 - qJDD(4)) * pkin(8) + (t809 + t899) * pkin(5) + t758;
t860 = sin(qJ(6));
t862 = cos(qJ(6));
t748 = -t860 * t750 + t862 * t751;
t811 = -t860 * qJD(4) + t862 * t830;
t771 = t811 * qJD(6) + t862 * qJDD(4) + t860 * t808;
t812 = t862 * qJD(4) + t860 * t830;
t775 = -t811 * mrSges(7,1) + t812 * mrSges(7,2);
t827 = qJD(6) + t831;
t782 = -t827 * mrSges(7,2) + t811 * mrSges(7,3);
t805 = qJDD(6) + t809;
t745 = m(7) * t748 + t805 * mrSges(7,1) - t771 * mrSges(7,3) - t812 * t775 + t827 * t782;
t749 = t862 * t750 + t860 * t751;
t770 = -t812 * qJD(6) - t860 * qJDD(4) + t862 * t808;
t783 = t827 * mrSges(7,1) - t812 * mrSges(7,3);
t746 = m(7) * t749 - t805 * mrSges(7,2) + t770 * mrSges(7,3) + t811 * t775 - t827 * t783;
t736 = t862 * t745 + t860 * t746;
t797 = -t830 * mrSges(6,2) - t831 * mrSges(6,3);
t875 = -m(6) * t758 - t809 * mrSges(6,1) - t831 * t797 - t736;
t819 = t830 * mrSges(6,1) - qJD(4) * mrSges(6,3);
t901 = -qJD(4) * mrSges(5,2) - t830 * mrSges(5,3) - t819;
t732 = m(5) * t760 - t809 * mrSges(5,3) + t901 * qJD(4) + t909 * qJDD(4) - t831 * t796 + t875;
t761 = t913 * t780 + t863 * t781;
t818 = qJD(4) * mrSges(5,1) - t831 * mrSges(5,3);
t874 = -t865 * pkin(4) + qJDD(4) * qJ(5) - t830 * t795 + t761;
t756 = qJD(4) * t914 - t874;
t820 = t831 * mrSges(6,1) + qJD(4) * mrSges(6,2);
t753 = -t808 * pkin(5) - t829 * pkin(8) + ((2 * qJD(5)) + t821) * qJD(4) + t874;
t877 = -m(7) * t753 + t770 * mrSges(7,1) - t771 * mrSges(7,2) + t811 * t782 - t812 * t783;
t872 = -m(6) * t756 + qJDD(4) * mrSges(6,3) + qJD(4) * t820 - t877;
t742 = m(5) * t761 - qJDD(4) * mrSges(5,2) - qJD(4) * t818 + (-t796 - t797) * t830 + (-mrSges(5,3) - mrSges(6,1)) * t808 + t872;
t725 = t863 * t732 + t913 * t742;
t883 = -qJDD(1) * mrSges(4,3) - t866 * (mrSges(4,1) * t858 + mrSges(4,2) * t859);
t723 = m(4) * t806 + t883 * t859 + t725;
t887 = -t913 * t732 + t863 * t742;
t724 = m(4) * t807 + t883 * t858 + t887;
t720 = t859 * t723 + t858 * t724;
t828 = -qJDD(1) * pkin(1) + t879;
t876 = -m(3) * t828 + t866 * mrSges(3,3) - t720;
t716 = m(2) * t835 - t866 * mrSges(2,2) + t910 * qJDD(1) + t876;
t815 = t905 * t866 + t878;
t755 = t808 * pkin(4) + t868;
t888 = -t860 * t745 + t862 * t746;
t880 = m(6) * t755 - t809 * mrSges(6,3) - t831 * t820 + t888;
t870 = m(5) * t785 + t809 * mrSges(5,2) + t909 * t808 + t831 * t818 + t901 * t830 + t880;
t869 = m(4) * t815 + mrSges(4,1) * t895 + mrSges(4,2) * t894 + t870;
t891 = t900 * mrSges(4,3);
t728 = t869 + m(2) * t836 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) + t891) * t866 + t915;
t904 = t864 * t716 + t861 * t728;
t903 = -qJD(4) * t919 + t830 * t920 - t831 * t921;
t902 = t921 * qJD(4) - t908 * t830 + t831 * t922;
t890 = -t861 * t716 + t864 * t728;
t889 = -t858 * t723 + t859 * t724;
t886 = Ifges(4,1) * t859 - Ifges(4,4) * t858;
t885 = Ifges(4,4) * t859 - Ifges(4,2) * t858;
t884 = Ifges(4,5) * t859 - Ifges(4,6) * t858;
t764 = Ifges(7,4) * t812 + Ifges(7,2) * t811 + Ifges(7,6) * t827;
t765 = Ifges(7,1) * t812 + Ifges(7,4) * t811 + Ifges(7,5) * t827;
t873 = mrSges(7,1) * t748 - mrSges(7,2) * t749 + Ifges(7,5) * t771 + Ifges(7,6) * t770 + Ifges(7,3) * t805 + t812 * t764 - t811 * t765;
t733 = -t808 * mrSges(6,2) - t830 * t819 + t880;
t763 = Ifges(7,5) * t812 + Ifges(7,6) * t811 + Ifges(7,3) * t827;
t738 = -mrSges(7,1) * t753 + mrSges(7,3) * t749 + Ifges(7,4) * t771 + Ifges(7,2) * t770 + Ifges(7,6) * t805 - t812 * t763 + t827 * t765;
t739 = mrSges(7,2) * t753 - mrSges(7,3) * t748 + Ifges(7,1) * t771 + Ifges(7,4) * t770 + Ifges(7,5) * t805 + t811 * t763 - t827 * t764;
t714 = -mrSges(5,1) * t785 - mrSges(6,1) * t756 + mrSges(6,2) * t755 + mrSges(5,3) * t761 - pkin(4) * t733 - pkin(5) * t877 - pkin(8) * t888 + t902 * qJD(4) + t920 * qJDD(4) - t862 * t738 - t860 * t739 + t923 * t808 + t908 * t809 + t903 * t831;
t721 = mrSges(6,1) * t758 + mrSges(5,2) * t785 - mrSges(5,3) * t760 - mrSges(6,3) * t755 + pkin(5) * t736 - qJ(5) * t733 - t918 * qJD(4) + t921 * qJDD(4) - t908 * t808 + t809 * t922 + t903 * t830 + t873;
t833 = t884 * qJD(1);
t711 = -mrSges(4,1) * t815 + mrSges(4,3) * t807 - pkin(3) * t870 + pkin(7) * t887 + t885 * qJDD(1) + t863 * t714 + t913 * t721 - t833 * t897;
t713 = -t858 * qJD(1) * t833 + mrSges(4,2) * t815 - mrSges(4,3) * t806 - pkin(7) * t725 + t886 * qJDD(1) - t913 * t714 + t863 * t721;
t718 = qJDD(1) * mrSges(3,2) - t876;
t730 = t866 * t891 + t869;
t871 = -mrSges(2,2) * t836 - mrSges(3,3) * t823 - pkin(1) * t718 - qJ(3) * t720 - t858 * t711 + t859 * t713 + qJ(2) * (t730 + t915) + mrSges(3,2) * t828 + mrSges(2,1) * t835 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t735 = qJDD(4) * mrSges(6,2) + qJD(4) * t819 - t875;
t867 = -mrSges(5,2) * t761 - mrSges(6,3) * t756 - pkin(4) * t735 - pkin(8) * t736 - t860 * t738 + t862 * t739 + t902 * t830 + qJ(5) * (-t830 * t797 + t872) + mrSges(6,2) * t758 + mrSges(5,1) * t760 + t918 * t831 + t921 * t809 + (-mrSges(6,1) * qJ(5) - t920) * t808 + t919 * qJDD(4);
t719 = -m(3) * g(3) + t889;
t710 = t867 - mrSges(2,3) * t835 + mrSges(3,1) * t828 + mrSges(4,1) * t806 - mrSges(4,2) * t807 + pkin(3) * t725 + pkin(2) * t720 - qJ(2) * t719 + (t884 + t907) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t858 * t886 + t859 * t885 + t906) * t866;
t709 = -mrSges(3,1) * t823 + mrSges(2,3) * t836 - pkin(1) * t719 + pkin(2) * t730 + t910 * g(3) - qJ(3) * t889 - t906 * qJDD(1) - t859 * t711 - t858 * t713 + t907 * t866;
t1 = [-m(1) * g(1) + t890; -m(1) * g(2) + t904; (-m(1) - m(2) - m(3)) * g(3) + t889; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t904 - t861 * t709 + t864 * t710; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t890 + t864 * t709 + t861 * t710; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t871; t871; t718; t730; t867; t735; t873;];
tauJB  = t1;
