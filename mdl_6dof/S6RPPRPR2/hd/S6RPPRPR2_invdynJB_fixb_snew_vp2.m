% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-05-05 14:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:02:11
% EndTime: 2019-05-05 14:02:17
% DurationCPUTime: 6.38s
% Computational Cost: add. (69974->303), mult. (157314->362), div. (0->0), fcn. (104208->10), ass. (0->136)
t907 = Ifges(5,1) + Ifges(6,2);
t900 = Ifges(5,4) + Ifges(6,6);
t899 = Ifges(5,5) - Ifges(6,4);
t906 = -Ifges(5,2) - Ifges(6,3);
t898 = Ifges(5,6) - Ifges(6,5);
t905 = Ifges(5,3) + Ifges(6,1);
t861 = qJD(1) ^ 2;
t856 = sin(qJ(4));
t853 = cos(pkin(10));
t902 = cos(qJ(4));
t882 = t853 * t902;
t851 = sin(pkin(10));
t888 = qJD(1) * t851;
t822 = -qJD(1) * t882 + t856 * t888;
t815 = mrSges(6,1) * t822 - qJD(4) * mrSges(6,3);
t857 = sin(qJ(1));
t859 = cos(qJ(1));
t832 = t857 * g(1) - g(2) * t859;
t829 = qJDD(1) * pkin(1) + t832;
t833 = -g(1) * t859 - g(2) * t857;
t830 = -pkin(1) * t861 + t833;
t852 = sin(pkin(9));
t854 = cos(pkin(9));
t809 = t852 * t829 + t854 * t830;
t792 = -pkin(2) * t861 + qJDD(1) * qJ(3) + t809;
t850 = -g(3) + qJDD(2);
t885 = qJD(1) * qJD(3);
t889 = t853 * t850 - 0.2e1 * t851 * t885;
t901 = pkin(3) * t853;
t772 = (-pkin(7) * qJDD(1) + t861 * t901 - t792) * t851 + t889;
t778 = t851 * t850 + (t792 + 0.2e1 * t885) * t853;
t883 = qJDD(1) * t853;
t846 = t853 ^ 2;
t895 = t846 * t861;
t775 = -pkin(3) * t895 + pkin(7) * t883 + t778;
t759 = t772 * t902 - t856 * t775;
t872 = t851 * t902 + t853 * t856;
t823 = t872 * qJD(1);
t797 = pkin(4) * t822 - qJ(5) * t823;
t860 = qJD(4) ^ 2;
t756 = -qJDD(4) * pkin(4) - t860 * qJ(5) + t823 * t797 + qJDD(5) - t759;
t887 = t822 * qJD(4);
t806 = qJDD(1) * t872 - t887;
t751 = (t822 * t823 - qJDD(4)) * pkin(8) + (t806 + t887) * pkin(5) + t756;
t884 = qJDD(1) * t851;
t886 = t823 * qJD(4);
t805 = -qJDD(1) * t882 + t856 * t884 + t886;
t817 = pkin(5) * t823 - qJD(4) * pkin(8);
t821 = t822 ^ 2;
t845 = t851 ^ 2;
t808 = t854 * t829 - t852 * t830;
t874 = qJDD(3) - t808;
t776 = (-pkin(2) - t901) * qJDD(1) + (-qJ(3) + (-t845 - t846) * pkin(7)) * t861 + t874;
t903 = -2 * qJD(5);
t862 = pkin(4) * t886 + t823 * t903 + (-t806 + t887) * qJ(5) + t776;
t754 = -t821 * pkin(5) - t823 * t817 + (pkin(4) + pkin(8)) * t805 + t862;
t855 = sin(qJ(6));
t858 = cos(qJ(6));
t749 = t751 * t858 - t754 * t855;
t810 = -qJD(4) * t855 + t822 * t858;
t774 = qJD(6) * t810 + qJDD(4) * t858 + t805 * t855;
t811 = qJD(4) * t858 + t822 * t855;
t779 = -mrSges(7,1) * t810 + mrSges(7,2) * t811;
t820 = qJD(6) + t823;
t782 = -mrSges(7,2) * t820 + mrSges(7,3) * t810;
t804 = qJDD(6) + t806;
t746 = m(7) * t749 + mrSges(7,1) * t804 - mrSges(7,3) * t774 - t779 * t811 + t782 * t820;
t750 = t751 * t855 + t754 * t858;
t773 = -qJD(6) * t811 - qJDD(4) * t855 + t805 * t858;
t783 = mrSges(7,1) * t820 - mrSges(7,3) * t811;
t747 = m(7) * t750 - mrSges(7,2) * t804 + mrSges(7,3) * t773 + t779 * t810 - t783 * t820;
t737 = t858 * t746 + t855 * t747;
t799 = -mrSges(6,2) * t822 - mrSges(6,3) * t823;
t870 = -m(6) * t756 - t806 * mrSges(6,1) - t823 * t799 - t737;
t736 = qJDD(4) * mrSges(6,2) + qJD(4) * t815 - t870;
t760 = t856 * t772 + t902 * t775;
t869 = -t860 * pkin(4) + qJDD(4) * qJ(5) - t822 * t797 + t760;
t753 = -t805 * pkin(5) - t821 * pkin(8) + ((2 * qJD(5)) + t817) * qJD(4) + t869;
t763 = Ifges(7,5) * t811 + Ifges(7,6) * t810 + Ifges(7,3) * t820;
t765 = Ifges(7,1) * t811 + Ifges(7,4) * t810 + Ifges(7,5) * t820;
t738 = -mrSges(7,1) * t753 + mrSges(7,3) * t750 + Ifges(7,4) * t774 + Ifges(7,2) * t773 + Ifges(7,6) * t804 - t763 * t811 + t765 * t820;
t764 = Ifges(7,4) * t811 + Ifges(7,2) * t810 + Ifges(7,6) * t820;
t739 = mrSges(7,2) * t753 - mrSges(7,3) * t749 + Ifges(7,1) * t774 + Ifges(7,4) * t773 + Ifges(7,5) * t804 + t763 * t810 - t764 * t820;
t755 = qJD(4) * t903 - t869;
t816 = mrSges(6,1) * t823 + qJD(4) * mrSges(6,2);
t871 = -m(7) * t753 + t773 * mrSges(7,1) - t774 * mrSges(7,2) + t810 * t782 - t811 * t783;
t867 = -m(6) * t755 + qJDD(4) * mrSges(6,3) + qJD(4) * t816 - t871;
t890 = t899 * qJD(4) - t900 * t822 + t907 * t823;
t891 = t898 * qJD(4) + t906 * t822 + t900 * t823;
t904 = t905 * qJDD(4) - t898 * t805 + t899 * t806 + t890 * t822 + t891 * t823 + mrSges(5,1) * t759 - mrSges(5,2) * t760 + mrSges(6,2) * t756 - mrSges(6,3) * t755 - pkin(4) * t736 - pkin(8) * t737 + qJ(5) * (-t805 * mrSges(6,1) - t822 * t799 + t867) - t855 * t738 + t858 * t739;
t896 = mrSges(4,2) * t851;
t798 = mrSges(5,1) * t822 + mrSges(5,2) * t823;
t813 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t822;
t734 = m(5) * t759 - t806 * mrSges(5,3) - t823 * t798 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t813 - t815) * qJD(4) + t870;
t814 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t823;
t742 = m(5) * t760 - qJDD(4) * mrSges(5,2) - qJD(4) * t814 + (-t798 - t799) * t822 + (-mrSges(5,3) - mrSges(6,1)) * t805 + t867;
t728 = t902 * t734 + t856 * t742;
t777 = -t792 * t851 + t889;
t873 = mrSges(4,3) * qJDD(1) + t861 * (-mrSges(4,1) * t853 + t896);
t726 = m(4) * t777 - t851 * t873 + t728;
t878 = -t856 * t734 + t902 * t742;
t727 = m(4) * t778 + t853 * t873 + t878;
t879 = -t726 * t851 + t853 * t727;
t718 = m(3) * t809 - mrSges(3,1) * t861 - qJDD(1) * mrSges(3,2) + t879;
t785 = -qJDD(1) * pkin(2) - t861 * qJ(3) + t874;
t758 = t805 * pkin(4) + t862;
t893 = -t855 * t746 + t858 * t747;
t735 = m(6) * t758 - t805 * mrSges(6,2) - t806 * mrSges(6,3) - t822 * t815 - t823 * t816 + t893;
t865 = m(5) * t776 + t805 * mrSges(5,1) + t806 * mrSges(5,2) + t822 * t813 + t823 * t814 + t735;
t864 = -m(4) * t785 + mrSges(4,1) * t883 - t865 + (t845 * t861 + t895) * mrSges(4,3);
t730 = t864 + (mrSges(3,1) - t896) * qJDD(1) - t861 * mrSges(3,2) + m(3) * t808;
t714 = t852 * t718 + t854 * t730;
t711 = m(2) * t832 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t861 + t714;
t880 = t854 * t718 - t730 * t852;
t712 = m(2) * t833 - mrSges(2,1) * t861 - qJDD(1) * mrSges(2,2) + t880;
t894 = t859 * t711 + t857 * t712;
t721 = t853 * t726 + t851 * t727;
t892 = -t905 * qJD(4) + t898 * t822 - t899 * t823;
t719 = m(3) * t850 + t721;
t881 = -t711 * t857 + t859 * t712;
t877 = Ifges(4,1) * t851 + Ifges(4,4) * t853;
t876 = Ifges(4,4) * t851 + Ifges(4,2) * t853;
t875 = Ifges(4,5) * t851 + Ifges(4,6) * t853;
t868 = mrSges(7,1) * t749 - mrSges(7,2) * t750 + Ifges(7,5) * t774 + Ifges(7,6) * t773 + Ifges(7,3) * t804 + t811 * t764 - t810 * t765;
t715 = -mrSges(5,1) * t776 - mrSges(6,1) * t755 + mrSges(6,2) * t758 + mrSges(5,3) * t760 - pkin(4) * t735 - pkin(5) * t871 - pkin(8) * t893 + t890 * qJD(4) + t898 * qJDD(4) - t858 * t738 - t855 * t739 + t906 * t805 + t900 * t806 + t892 * t823;
t722 = mrSges(6,1) * t756 + mrSges(5,2) * t776 - mrSges(5,3) * t759 - mrSges(6,3) * t758 + pkin(5) * t737 - qJ(5) * t735 - t891 * qJD(4) + t899 * qJDD(4) - t900 * t805 + t907 * t806 + t892 * t822 + t868;
t828 = t875 * qJD(1);
t704 = -mrSges(4,1) * t785 + mrSges(4,3) * t778 - pkin(3) * t865 + pkin(7) * t878 + qJDD(1) * t876 + t715 * t902 + t856 * t722 - t828 * t888;
t707 = t853 * qJD(1) * t828 + mrSges(4,2) * t785 - mrSges(4,3) * t777 - pkin(7) * t728 + qJDD(1) * t877 - t856 * t715 + t722 * t902;
t732 = mrSges(4,2) * t884 - t864;
t866 = mrSges(2,1) * t832 + mrSges(3,1) * t808 - mrSges(2,2) * t833 - mrSges(3,2) * t809 + pkin(1) * t714 - pkin(2) * t732 + qJ(3) * t879 + t853 * t704 + t851 * t707 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t705 = (Ifges(3,6) - t875) * qJDD(1) - mrSges(3,1) * t850 + mrSges(3,3) * t809 - mrSges(4,1) * t777 + mrSges(4,2) * t778 - pkin(3) * t728 - pkin(2) * t721 + (-t851 * t876 + t853 * t877 + Ifges(3,5)) * t861 - t904;
t702 = mrSges(3,2) * t850 - mrSges(3,3) * t808 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t861 - qJ(3) * t721 - t704 * t851 + t707 * t853;
t701 = -mrSges(2,2) * g(3) - mrSges(2,3) * t832 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t861 - qJ(2) * t714 + t702 * t854 - t705 * t852;
t700 = mrSges(2,1) * g(3) + mrSges(2,3) * t833 + t861 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t719 + qJ(2) * t880 + t852 * t702 + t854 * t705;
t1 = [-m(1) * g(1) + t881; -m(1) * g(2) + t894; (-m(1) - m(2)) * g(3) + t719; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t894 - t857 * t700 + t859 * t701; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t881 + t859 * t700 + t857 * t701; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t866; t866; t719; t732; t904; t736; t868;];
tauJB  = t1;
