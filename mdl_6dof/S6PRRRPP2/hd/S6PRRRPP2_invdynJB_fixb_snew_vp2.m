% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:42:33
% EndTime: 2019-05-05 06:42:42
% DurationCPUTime: 5.64s
% Computational Cost: add. (69471->295), mult. (131328->355), div. (0->0), fcn. (85247->10), ass. (0->124)
t903 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t883 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t882 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t902 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t881 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t901 = Ifges(7,3) + Ifges(5,3) + Ifges(6,2);
t851 = sin(pkin(10));
t853 = cos(pkin(10));
t837 = g(1) * t851 - g(2) * t853;
t838 = -g(1) * t853 - g(2) * t851;
t850 = -g(3) + qJDD(1);
t859 = cos(qJ(2));
t854 = cos(pkin(6));
t857 = sin(qJ(2));
t889 = t854 * t857;
t852 = sin(pkin(6));
t890 = t852 * t857;
t779 = t837 * t889 + t859 * t838 + t850 * t890;
t861 = qJD(2) ^ 2;
t773 = -pkin(2) * t861 + qJDD(2) * pkin(8) + t779;
t816 = -t837 * t852 + t850 * t854;
t856 = sin(qJ(3));
t858 = cos(qJ(3));
t769 = t858 * t773 + t856 * t816;
t834 = (-pkin(3) * t858 - pkin(9) * t856) * qJD(2);
t860 = qJD(3) ^ 2;
t885 = qJD(2) * t858;
t765 = -pkin(3) * t860 + qJDD(3) * pkin(9) + t834 * t885 + t769;
t778 = -t857 * t838 + (t837 * t854 + t850 * t852) * t859;
t772 = -qJDD(2) * pkin(2) - t861 * pkin(8) - t778;
t884 = qJD(2) * qJD(3);
t875 = t858 * t884;
t835 = qJDD(2) * t856 + t875;
t846 = t856 * t884;
t836 = qJDD(2) * t858 - t846;
t767 = (-t835 - t875) * pkin(9) + (-t836 + t846) * pkin(3) + t772;
t855 = sin(qJ(4));
t895 = cos(qJ(4));
t760 = -t855 * t765 + t767 * t895;
t886 = qJD(2) * t856;
t831 = -qJD(3) * t895 + t855 * t886;
t832 = t855 * qJD(3) + t886 * t895;
t803 = pkin(4) * t831 - qJ(5) * t832;
t828 = -qJDD(4) + t836;
t844 = -qJD(4) + t885;
t843 = t844 ^ 2;
t758 = t828 * pkin(4) - t843 * qJ(5) + t832 * t803 + qJDD(5) - t760;
t815 = -mrSges(6,2) * t831 - mrSges(6,3) * t844;
t900 = -m(6) * t758 - t828 * mrSges(6,1) - t844 * t815;
t800 = -t831 * qJD(4) + t855 * qJDD(3) + t835 * t895;
t768 = -t856 * t773 + t858 * t816;
t865 = qJDD(3) * pkin(3) + t860 * pkin(9) - t834 * t886 + t768;
t891 = t831 * t844;
t899 = -(t800 + t891) * qJ(5) - t865;
t799 = qJD(4) * t832 - qJDD(3) * t895 + t835 * t855;
t811 = pkin(5) * t844 - qJ(6) * t832;
t827 = t831 ^ 2;
t755 = -t827 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t799 + (pkin(4) * t844 + (2 * qJD(5)) + t811) * t832 - t899;
t809 = -mrSges(7,2) * t844 + mrSges(7,3) * t831;
t812 = mrSges(7,1) * t844 - mrSges(7,3) * t832;
t750 = m(7) * t755 - t799 * mrSges(7,1) + t800 * mrSges(7,2) - t831 * t809 + t832 * t812;
t896 = -2 * qJD(5);
t759 = t832 * t896 + (-t832 * t844 + t799) * pkin(4) + t899;
t814 = mrSges(6,1) * t844 + mrSges(6,2) * t832;
t748 = m(6) * t759 + t799 * mrSges(6,1) - t800 * mrSges(6,3) - t832 * t814 + t831 * t815 - t750;
t761 = t895 * t765 + t855 * t767;
t757 = -pkin(4) * t843 - t828 * qJ(5) - t831 * t803 + t844 * t896 + t761;
t753 = -pkin(5) * t827 + qJ(6) * t799 + 0.2e1 * qJD(6) * t831 - t811 * t844 + t757;
t876 = -t883 * t831 + t903 * t832 - t882 * t844;
t878 = t881 * t831 - t882 * t832 + t901 * t844;
t805 = -mrSges(7,1) * t831 + mrSges(7,2) * t832;
t879 = m(7) * t753 + t799 * mrSges(7,3) + t831 * t805;
t723 = mrSges(5,1) * t865 + mrSges(5,3) * t761 - mrSges(6,1) * t759 + mrSges(6,2) * t757 + mrSges(7,1) * t755 - mrSges(7,3) * t753 + pkin(5) * t750 - qJ(6) * t879 - pkin(4) * t748 + (qJ(6) * t812 - t876) * t844 + t878 * t832 + (mrSges(7,2) * qJ(6) - t881) * t828 + t883 * t800 + t902 * t799;
t751 = -0.2e1 * qJD(6) * t832 + (-t800 + t891) * qJ(6) + (t831 * t832 + t828) * pkin(5) + t758;
t870 = -m(7) * t751 + t800 * mrSges(7,3) + t832 * t805;
t749 = t828 * mrSges(7,1) + t844 * t809 - t870;
t877 = t902 * t831 + t883 * t832 - t881 * t844;
t731 = -mrSges(5,2) * t865 + mrSges(6,2) * t758 + mrSges(7,2) * t755 - mrSges(5,3) * t760 - mrSges(6,3) * t759 - mrSges(7,3) * t751 - qJ(5) * t748 - qJ(6) * t749 - t883 * t799 + t903 * t800 - t882 * t828 + t878 * t831 + t877 * t844;
t810 = mrSges(5,2) * t844 - mrSges(5,3) * t831;
t804 = mrSges(6,1) * t831 - mrSges(6,3) * t832;
t887 = -mrSges(5,1) * t831 - mrSges(5,2) * t832 - t804;
t893 = -mrSges(5,3) - mrSges(6,2);
t743 = m(5) * t760 + (-t809 - t810) * t844 + t887 * t832 + (-mrSges(5,1) - mrSges(7,1)) * t828 + t893 * t800 + t870 + t900;
t813 = -mrSges(5,1) * t844 - mrSges(5,3) * t832;
t867 = m(6) * t757 - t828 * mrSges(6,3) - t844 * t814 + t879;
t744 = m(5) * t761 + (-t812 + t813) * t844 + t887 * t831 + (mrSges(5,2) - mrSges(7,2)) * t828 + t893 * t799 + t867;
t739 = -t743 * t855 + t895 * t744;
t745 = m(5) * t865 - t799 * mrSges(5,1) - t800 * mrSges(5,2) - t831 * t810 - t832 * t813 - t748;
t820 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t856 + Ifges(4,2) * t858) * qJD(2);
t821 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t856 + Ifges(4,4) * t858) * qJD(2);
t898 = mrSges(4,1) * t768 - mrSges(4,2) * t769 + Ifges(4,5) * t835 + Ifges(4,6) * t836 + Ifges(4,3) * qJDD(3) + pkin(3) * t745 + pkin(9) * t739 + (t820 * t856 - t821 * t858) * qJD(2) + t723 * t895 + t855 * t731;
t746 = t800 * mrSges(6,2) + t832 * t804 + t749 - t900;
t897 = -t799 * t881 + t800 * t882 - t901 * t828 + t831 * t876 + t832 * t877 + mrSges(5,1) * t760 - mrSges(6,1) * t758 - mrSges(7,1) * t751 - mrSges(5,2) * t761 + mrSges(7,2) * t753 + mrSges(6,3) * t757 - pkin(4) * t746 - pkin(5) * t749 + qJ(5) * (-t799 * mrSges(6,2) - t828 * mrSges(7,2) - t831 * t804 - t844 * t812 + t867);
t738 = t743 * t895 + t855 * t744;
t839 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t886;
t840 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t885;
t864 = -m(4) * t772 + t836 * mrSges(4,1) - t835 * mrSges(4,2) - t839 * t886 + t840 * t885 - t738;
t734 = m(3) * t778 + qJDD(2) * mrSges(3,1) - t861 * mrSges(3,2) + t864;
t892 = t734 * t859;
t833 = (-mrSges(4,1) * t858 + mrSges(4,2) * t856) * qJD(2);
t737 = m(4) * t769 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t836 - qJD(3) * t839 + t833 * t885 + t739;
t741 = m(4) * t768 + qJDD(3) * mrSges(4,1) - t835 * mrSges(4,3) + qJD(3) * t840 - t833 * t886 + t745;
t873 = t858 * t737 - t741 * t856;
t727 = m(3) * t779 - mrSges(3,1) * t861 - qJDD(2) * mrSges(3,2) + t873;
t730 = t856 * t737 + t858 * t741;
t729 = m(3) * t816 + t730;
t716 = t727 * t889 - t729 * t852 + t854 * t892;
t714 = m(2) * t837 + t716;
t722 = t859 * t727 - t734 * t857;
t721 = m(2) * t838 + t722;
t888 = t853 * t714 + t851 * t721;
t715 = t727 * t890 + t854 * t729 + t852 * t892;
t874 = -t714 * t851 + t853 * t721;
t871 = m(2) * t850 + t715;
t819 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t856 + Ifges(4,6) * t858) * qJD(2);
t717 = mrSges(4,2) * t772 - mrSges(4,3) * t768 + Ifges(4,1) * t835 + Ifges(4,4) * t836 + Ifges(4,5) * qJDD(3) - pkin(9) * t738 - qJD(3) * t820 - t855 * t723 + t731 * t895 + t819 * t885;
t718 = -mrSges(4,1) * t772 + mrSges(4,3) * t769 + Ifges(4,4) * t835 + Ifges(4,2) * t836 + Ifges(4,6) * qJDD(3) - pkin(3) * t738 + qJD(3) * t821 - t819 * t886 - t897;
t711 = mrSges(3,2) * t816 - mrSges(3,3) * t778 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t861 - pkin(8) * t730 + t717 * t858 - t718 * t856;
t712 = -mrSges(3,1) * t816 + mrSges(3,3) * t779 + t861 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t730 - t898;
t866 = pkin(7) * t722 + t711 * t857 + t712 * t859;
t710 = mrSges(3,1) * t778 - mrSges(3,2) * t779 + Ifges(3,3) * qJDD(2) + pkin(2) * t864 + pkin(8) * t873 + t856 * t717 + t858 * t718;
t709 = mrSges(2,2) * t850 - mrSges(2,3) * t837 + t859 * t711 - t857 * t712 + (-t715 * t852 - t716 * t854) * pkin(7);
t708 = -mrSges(2,1) * t850 + mrSges(2,3) * t838 - pkin(1) * t715 - t852 * t710 + t854 * t866;
t1 = [-m(1) * g(1) + t874; -m(1) * g(2) + t888; -m(1) * g(3) + t871; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t888 - t851 * t708 + t853 * t709; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t874 + t853 * t708 + t851 * t709; -mrSges(1,1) * g(2) + mrSges(2,1) * t837 + mrSges(1,2) * g(1) - mrSges(2,2) * t838 + pkin(1) * t716 + t854 * t710 + t852 * t866; t871; t710; t898; t897; t746; t750;];
tauJB  = t1;
