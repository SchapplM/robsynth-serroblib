% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:21:03
% EndTime: 2019-05-05 21:21:10
% DurationCPUTime: 4.13s
% Computational Cost: add. (44717->297), mult. (84327->346), div. (0->0), fcn. (49219->8), ass. (0->117)
t872 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t855 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t854 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t871 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t853 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t870 = Ifges(7,3) + Ifges(5,3) + Ifges(6,2);
t830 = sin(qJ(1));
t832 = cos(qJ(1));
t811 = g(1) * t830 - g(2) * t832;
t801 = qJDD(1) * pkin(1) + t811;
t812 = -g(1) * t832 - g(2) * t830;
t834 = qJD(1) ^ 2;
t803 = -pkin(1) * t834 + t812;
t826 = sin(pkin(9));
t827 = cos(pkin(9));
t770 = t801 * t827 - t803 * t826;
t746 = -qJDD(1) * pkin(2) - pkin(7) * t834 - t770;
t829 = sin(qJ(3));
t831 = cos(qJ(3));
t856 = qJD(1) * qJD(3);
t847 = t831 * t856;
t805 = qJDD(1) * t829 + t847;
t818 = t829 * t856;
t806 = qJDD(1) * t831 - t818;
t733 = (-t805 - t847) * pkin(8) + (-t806 + t818) * pkin(3) + t746;
t771 = t801 * t826 + t803 * t827;
t747 = -pkin(2) * t834 + qJDD(1) * pkin(7) + t771;
t825 = -g(3) + qJDD(2);
t739 = t747 * t831 + t825 * t829;
t804 = (-pkin(3) * t831 - pkin(8) * t829) * qJD(1);
t833 = qJD(3) ^ 2;
t857 = qJD(1) * t831;
t737 = -pkin(3) * t833 + qJDD(3) * pkin(8) + t804 * t857 + t739;
t828 = sin(qJ(4));
t864 = cos(qJ(4));
t730 = t733 * t864 - t737 * t828;
t858 = qJD(1) * t829;
t799 = -qJD(3) * t864 + t828 * t858;
t800 = qJD(3) * t828 + t858 * t864;
t774 = pkin(4) * t799 - qJ(5) * t800;
t798 = -qJDD(4) + t806;
t814 = -qJD(4) + t857;
t813 = t814 ^ 2;
t728 = t798 * pkin(4) - t813 * qJ(5) + t774 * t800 + qJDD(5) - t730;
t784 = -mrSges(6,2) * t799 - mrSges(6,3) * t814;
t869 = -m(6) * t728 - mrSges(6,1) * t798 - t784 * t814;
t768 = -qJD(4) * t799 + qJDD(3) * t828 + t805 * t864;
t738 = -t747 * t829 + t825 * t831;
t839 = qJDD(3) * pkin(3) + pkin(8) * t833 - t804 * t858 + t738;
t861 = t799 * t814;
t868 = -(t768 + t861) * qJ(5) - t839;
t767 = qJD(4) * t800 - qJDD(3) * t864 + t805 * t828;
t780 = pkin(5) * t814 - qJ(6) * t800;
t797 = t799 ^ 2;
t725 = -qJ(6) * t797 + qJDD(6) + (-pkin(4) - pkin(5)) * t767 + (pkin(4) * t814 + (2 * qJD(5)) + t780) * t800 - t868;
t778 = -mrSges(7,2) * t814 + mrSges(7,3) * t799;
t781 = mrSges(7,1) * t814 - mrSges(7,3) * t800;
t720 = m(7) * t725 - mrSges(7,1) * t767 + mrSges(7,2) * t768 - t778 * t799 + t781 * t800;
t865 = -2 * qJD(5);
t729 = t800 * t865 + (-t800 * t814 + t767) * pkin(4) + t868;
t783 = mrSges(6,1) * t814 + mrSges(6,2) * t800;
t718 = m(6) * t729 + mrSges(6,1) * t767 - mrSges(6,3) * t768 - t783 * t800 + t784 * t799 - t720;
t731 = t733 * t828 + t737 * t864;
t727 = -pkin(4) * t813 - qJ(5) * t798 - t774 * t799 + t814 * t865 + t731;
t723 = -pkin(5) * t797 + qJ(6) * t767 + 0.2e1 * qJD(6) * t799 - t780 * t814 + t727;
t848 = -t799 * t855 + t800 * t872 - t814 * t854;
t850 = t799 * t853 - t800 * t854 + t814 * t870;
t776 = -mrSges(7,1) * t799 + mrSges(7,2) * t800;
t851 = m(7) * t723 + mrSges(7,3) * t767 + t776 * t799;
t694 = mrSges(5,1) * t839 + mrSges(5,3) * t731 - mrSges(6,1) * t729 + mrSges(6,2) * t727 + mrSges(7,1) * t725 - mrSges(7,3) * t723 + pkin(5) * t720 - qJ(6) * t851 - pkin(4) * t718 + (qJ(6) * t781 - t848) * t814 + t850 * t800 + (mrSges(7,2) * qJ(6) - t853) * t798 + t855 * t768 + t871 * t767;
t721 = -0.2e1 * qJD(6) * t800 + (-t768 + t861) * qJ(6) + (t799 * t800 + t798) * pkin(5) + t728;
t842 = -m(7) * t721 + mrSges(7,3) * t768 + t776 * t800;
t719 = mrSges(7,1) * t798 + t778 * t814 - t842;
t849 = t799 * t871 + t800 * t855 - t814 * t853;
t701 = -mrSges(5,2) * t839 + mrSges(6,2) * t728 + mrSges(7,2) * t725 - mrSges(5,3) * t730 - mrSges(6,3) * t729 - mrSges(7,3) * t721 - qJ(5) * t718 - qJ(6) * t719 - t767 * t855 + t768 * t872 - t798 * t854 + t799 * t850 + t814 * t849;
t779 = mrSges(5,2) * t814 - mrSges(5,3) * t799;
t775 = mrSges(6,1) * t799 - mrSges(6,3) * t800;
t859 = -mrSges(5,1) * t799 - mrSges(5,2) * t800 - t775;
t862 = -mrSges(5,3) - mrSges(6,2);
t713 = m(5) * t730 + (-t778 - t779) * t814 + t859 * t800 + (-mrSges(5,1) - mrSges(7,1)) * t798 + t862 * t768 + t842 + t869;
t782 = -mrSges(5,1) * t814 - mrSges(5,3) * t800;
t840 = m(6) * t727 - mrSges(6,3) * t798 - t783 * t814 + t851;
t714 = m(5) * t731 + (-t781 + t782) * t814 + t859 * t799 + (mrSges(5,2) - mrSges(7,2)) * t798 + t862 * t767 + t840;
t709 = -t713 * t828 + t714 * t864;
t715 = m(5) * t839 - mrSges(5,1) * t767 - mrSges(5,2) * t768 - t779 * t799 - t782 * t800 - t718;
t790 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t829 + Ifges(4,2) * t831) * qJD(1);
t791 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t829 + Ifges(4,4) * t831) * qJD(1);
t867 = mrSges(4,1) * t738 - mrSges(4,2) * t739 + Ifges(4,5) * t805 + Ifges(4,6) * t806 + Ifges(4,3) * qJDD(3) + pkin(3) * t715 + pkin(8) * t709 + (t790 * t829 - t791 * t831) * qJD(1) + t694 * t864 + t828 * t701;
t716 = mrSges(6,2) * t768 + t775 * t800 + t719 - t869;
t866 = -t767 * t853 + t768 * t854 - t870 * t798 + t799 * t848 + t800 * t849 + mrSges(5,1) * t730 - mrSges(6,1) * t728 - mrSges(7,1) * t721 - mrSges(5,2) * t731 + mrSges(7,2) * t723 + mrSges(6,3) * t727 - pkin(4) * t716 - pkin(5) * t719 + qJ(5) * (-mrSges(6,2) * t767 - t798 * mrSges(7,2) - t775 * t799 - t781 * t814 + t840);
t802 = (-mrSges(4,1) * t831 + mrSges(4,2) * t829) * qJD(1);
t808 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t858;
t707 = m(4) * t739 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t806 - qJD(3) * t808 + t802 * t857 + t709;
t809 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t857;
t711 = m(4) * t738 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t805 + qJD(3) * t809 - t802 * t858 + t715;
t844 = t707 * t831 - t711 * t829;
t697 = m(3) * t771 - mrSges(3,1) * t834 - qJDD(1) * mrSges(3,2) + t844;
t708 = t713 * t864 + t714 * t828;
t837 = -m(4) * t746 + mrSges(4,1) * t806 - mrSges(4,2) * t805 - t808 * t858 + t809 * t857 - t708;
t703 = m(3) * t770 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t834 + t837;
t693 = t697 * t826 + t703 * t827;
t690 = m(2) * t811 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t834 + t693;
t845 = t697 * t827 - t703 * t826;
t691 = m(2) * t812 - mrSges(2,1) * t834 - qJDD(1) * mrSges(2,2) + t845;
t860 = t690 * t832 + t691 * t830;
t700 = t707 * t829 + t711 * t831;
t698 = m(3) * t825 + t700;
t846 = -t690 * t830 + t691 * t832;
t789 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t829 + Ifges(4,6) * t831) * qJD(1);
t684 = mrSges(4,2) * t746 - mrSges(4,3) * t738 + Ifges(4,1) * t805 + Ifges(4,4) * t806 + Ifges(4,5) * qJDD(3) - pkin(8) * t708 - qJD(3) * t790 - t694 * t828 + t701 * t864 + t789 * t857;
t686 = -mrSges(4,1) * t746 + mrSges(4,3) * t739 + Ifges(4,4) * t805 + Ifges(4,2) * t806 + Ifges(4,6) * qJDD(3) - pkin(3) * t708 + qJD(3) * t791 - t789 * t858 - t866;
t838 = mrSges(2,1) * t811 + mrSges(3,1) * t770 - mrSges(2,2) * t812 - mrSges(3,2) * t771 + pkin(1) * t693 + pkin(2) * t837 + pkin(7) * t844 + t829 * t684 + t831 * t686 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t682 = -mrSges(3,1) * t825 + mrSges(3,3) * t771 + t834 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t700 - t867;
t681 = mrSges(3,2) * t825 - mrSges(3,3) * t770 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t834 - pkin(7) * t700 + t684 * t831 - t686 * t829;
t680 = -mrSges(2,2) * g(3) - mrSges(2,3) * t811 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t834 - qJ(2) * t693 + t681 * t827 - t682 * t826;
t679 = mrSges(2,1) * g(3) + mrSges(2,3) * t812 + t834 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t698 + qJ(2) * t845 + t826 * t681 + t827 * t682;
t1 = [-m(1) * g(1) + t846; -m(1) * g(2) + t860; (-m(1) - m(2)) * g(3) + t698; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t860 - t679 * t830 + t680 * t832; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t846 + t832 * t679 + t830 * t680; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t838; t838; t698; t867; t866; t716; t720;];
tauJB  = t1;
