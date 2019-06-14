% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:51:05
% EndTime: 2019-05-05 21:51:10
% DurationCPUTime: 2.79s
% Computational Cost: add. (25898->293), mult. (48755->326), div. (0->0), fcn. (27295->6), ass. (0->119)
t856 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t835 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t834 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t855 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t833 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t854 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t803 = sin(qJ(1));
t805 = cos(qJ(1));
t787 = -t805 * g(1) - t803 * g(2);
t853 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t787;
t807 = qJD(1) ^ 2;
t849 = (-pkin(1) - pkin(7));
t754 = (t849 * t807) - t853;
t802 = sin(qJ(3));
t804 = cos(qJ(3));
t837 = qJD(1) * qJD(3);
t826 = t804 * t837;
t781 = -t802 * qJDD(1) - t826;
t827 = t802 * t837;
t782 = t804 * qJDD(1) - t827;
t708 = (-t782 + t827) * pkin(8) + (-t781 + t826) * pkin(3) + t754;
t786 = t803 * g(1) - t805 * g(2);
t819 = -t807 * qJ(2) + qJDD(2) - t786;
t756 = t849 * qJDD(1) + t819;
t744 = -t804 * g(3) + t802 * t756;
t780 = (pkin(3) * t802 - pkin(8) * t804) * qJD(1);
t806 = qJD(3) ^ 2;
t838 = t802 * qJD(1);
t712 = -t806 * pkin(3) + qJDD(3) * pkin(8) - t780 * t838 + t744;
t801 = sin(qJ(4));
t848 = cos(qJ(4));
t705 = t848 * t708 - t801 * t712;
t839 = qJD(1) * t804;
t777 = -t848 * qJD(3) + t801 * t839;
t778 = t801 * qJD(3) + t848 * t839;
t740 = t777 * pkin(4) - t778 * qJ(5);
t776 = qJDD(4) - t781;
t789 = qJD(4) + t838;
t788 = t789 ^ 2;
t703 = -t776 * pkin(4) - t788 * qJ(5) + t778 * t740 + qJDD(5) - t705;
t734 = -t777 * qJD(4) + t801 * qJDD(3) + t848 * t782;
t742 = -t777 * mrSges(6,2) - t778 * mrSges(6,3);
t852 = -m(6) * t703 - t734 * mrSges(6,1) - t778 * t742;
t739 = -t778 * mrSges(7,2) + t777 * mrSges(7,3);
t842 = t777 * t789;
t697 = -0.2e1 * qJD(6) * t789 + (t777 * t778 - t776) * qJ(6) + (t734 + t842) * pkin(5) + t703;
t750 = -t777 * mrSges(7,1) + t789 * mrSges(7,2);
t822 = -m(7) * t697 + t776 * mrSges(7,3) + t789 * t750;
t695 = t734 * mrSges(7,1) + t778 * t739 - t822;
t749 = t777 * mrSges(6,1) - t789 * mrSges(6,3);
t693 = t776 * mrSges(6,2) + t789 * t749 + t695 - t852;
t733 = t778 * qJD(4) - t848 * qJDD(3) + t801 * t782;
t747 = t778 * pkin(5) - t789 * qJ(6);
t775 = t777 ^ 2;
t706 = t801 * t708 + t848 * t712;
t814 = -t788 * pkin(4) + t776 * qJ(5) - t777 * t740 + t706;
t701 = -t733 * pkin(5) - t775 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t747) * t789 + t814;
t850 = -2 * qJD(5);
t702 = t789 * t850 - t814;
t751 = t778 * mrSges(6,1) + t789 * mrSges(6,2);
t748 = t778 * mrSges(7,1) - t789 * mrSges(7,3);
t831 = m(7) * t701 + t776 * mrSges(7,2) + t789 * t748;
t818 = -m(6) * t702 + t776 * mrSges(6,3) + t789 * t751 + t831;
t828 = -t835 * t777 + t856 * t778 + t834 * t789;
t829 = t855 * t777 + t835 * t778 + t833 * t789;
t840 = -t739 - t742;
t851 = -t833 * t733 + t834 * t734 + t854 * t776 + t828 * t777 + t829 * t778 + mrSges(5,1) * t705 - mrSges(5,2) * t706 + mrSges(6,2) * t703 + mrSges(7,2) * t701 - mrSges(6,3) * t702 - mrSges(7,3) * t697 - pkin(4) * t693 + qJ(5) * (t840 * t777 + (-mrSges(6,1) - mrSges(7,1)) * t733 + t818) - qJ(6) * t695;
t846 = mrSges(2,1) - mrSges(3,2);
t845 = -mrSges(7,1) - mrSges(5,3);
t844 = -Ifges(3,4) + Ifges(2,5);
t843 = (Ifges(3,5) - Ifges(2,6));
t779 = (mrSges(4,1) * t802 + mrSges(4,2) * t804) * qJD(1);
t785 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t839;
t741 = t777 * mrSges(5,1) + t778 * mrSges(5,2);
t745 = -t789 * mrSges(5,2) - t777 * mrSges(5,3);
t687 = m(5) * t705 + (t745 - t749) * t789 + (-t739 - t741) * t778 + (mrSges(5,1) - mrSges(6,2)) * t776 + t845 * t734 + t822 + t852;
t746 = t789 * mrSges(5,1) - t778 * mrSges(5,3);
t690 = m(5) * t706 - t776 * mrSges(5,2) - t789 * t746 + (-t741 + t840) * t777 + (-mrSges(6,1) + t845) * t733 + t818;
t823 = -t801 * t687 + t848 * t690;
t682 = m(4) * t744 - qJDD(3) * mrSges(4,2) + t781 * mrSges(4,3) - qJD(3) * t785 - t779 * t838 + t823;
t743 = t802 * g(3) + t804 * t756;
t784 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t838;
t711 = -qJDD(3) * pkin(3) - t806 * pkin(8) + t780 * t839 - t743;
t810 = (-t734 + t842) * qJ(5) + t711 + (t789 * pkin(4) + t850) * t778;
t704 = t733 * pkin(4) + t810;
t699 = -t775 * pkin(5) + 0.2e1 * qJD(6) * t777 - t778 * t747 + (pkin(4) + qJ(6)) * t733 + t810;
t821 = m(7) * t699 - t734 * mrSges(7,2) + t733 * mrSges(7,3) - t778 * t748 + t777 * t750;
t815 = -m(6) * t704 + t733 * mrSges(6,2) + t777 * t749 - t821;
t809 = -m(5) * t711 - t733 * mrSges(5,1) - t777 * t745 + (-t746 + t751) * t778 + (-mrSges(5,2) + mrSges(6,3)) * t734 + t815;
t685 = m(4) * t743 + qJDD(3) * mrSges(4,1) - t782 * mrSges(4,3) + qJD(3) * t784 - t779 * t839 + t809;
t674 = t802 * t682 + t804 * t685;
t761 = -qJDD(1) * pkin(1) + t819;
t817 = -m(3) * t761 + (t807 * mrSges(3,3)) - t674;
t670 = m(2) * t786 - (t807 * mrSges(2,2)) + t846 * qJDD(1) + t817;
t759 = t807 * pkin(1) + t853;
t684 = t848 * t687 + t801 * t690;
t816 = -m(4) * t754 + t781 * mrSges(4,1) - t782 * mrSges(4,2) - t784 * t838 - t785 * t839 - t684;
t812 = -m(3) * t759 + (t807 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t816;
t679 = m(2) * t787 - (t807 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t812;
t841 = t805 * t670 + t803 * t679;
t830 = t833 * t777 - t834 * t778 - t854 * t789;
t825 = -t803 * t670 + t805 * t679;
t824 = t804 * t682 - t802 * t685;
t692 = -t734 * mrSges(6,3) - t778 * t751 - t815;
t696 = -t733 * mrSges(7,1) - t777 * t739 + t831;
t668 = -mrSges(5,1) * t711 - mrSges(6,1) * t702 + mrSges(7,1) * t701 + mrSges(6,2) * t704 + mrSges(5,3) * t706 - mrSges(7,3) * t699 - pkin(4) * t692 + pkin(5) * t696 - qJ(6) * t821 + t855 * t733 + t835 * t734 + t833 * t776 + t830 * t778 + t828 * t789;
t676 = mrSges(6,1) * t703 + mrSges(7,1) * t697 + mrSges(5,2) * t711 - mrSges(7,2) * t699 - mrSges(5,3) * t705 - mrSges(6,3) * t704 + pkin(5) * t695 - qJ(5) * t692 - t835 * t733 + t856 * t734 + t834 * t776 + t830 * t777 - t829 * t789;
t765 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t804 - Ifges(4,2) * t802) * qJD(1);
t766 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t804 - Ifges(4,4) * t802) * qJD(1);
t813 = mrSges(4,1) * t743 - mrSges(4,2) * t744 + Ifges(4,5) * t782 + Ifges(4,6) * t781 + Ifges(4,3) * qJDD(3) + pkin(3) * t809 + pkin(8) * t823 + t848 * t668 + t801 * t676 + t765 * t839 + t766 * t838;
t764 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t804 - Ifges(4,6) * t802) * qJD(1);
t665 = mrSges(4,2) * t754 - mrSges(4,3) * t743 + Ifges(4,1) * t782 + Ifges(4,4) * t781 + Ifges(4,5) * qJDD(3) - pkin(8) * t684 - qJD(3) * t765 - t801 * t668 + t848 * t676 - t764 * t838;
t666 = -mrSges(4,1) * t754 + mrSges(4,3) * t744 + Ifges(4,4) * t782 + Ifges(4,2) * t781 + Ifges(4,6) * qJDD(3) - pkin(3) * t684 + qJD(3) * t766 - t764 * t839 - t851;
t672 = qJDD(1) * mrSges(3,2) - t817;
t811 = mrSges(2,1) * t786 - mrSges(2,2) * t787 + mrSges(3,2) * t761 - mrSges(3,3) * t759 - pkin(1) * t672 - pkin(7) * t674 + qJ(2) * t812 + t804 * t665 - t802 * t666 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t673 = -m(3) * g(3) + t824;
t663 = -mrSges(2,3) * t786 + mrSges(3,1) * t761 + pkin(2) * t674 - qJ(2) * t673 + t813 + (t843 * t807) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t844 * qJDD(1);
t662 = -mrSges(3,1) * t759 + mrSges(2,3) * t787 - pkin(1) * t673 - pkin(2) * t816 - pkin(7) * t824 + t846 * g(3) - t843 * qJDD(1) - t802 * t665 - t804 * t666 + t844 * t807;
t1 = [-m(1) * g(1) + t825; -m(1) * g(2) + t841; (-m(1) - m(2) - m(3)) * g(3) + t824; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t841 - t803 * t662 + t805 * t663; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t825 + t805 * t662 + t803 * t663; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t811; t811; t672; t813; t851; t693; t696;];
tauJB  = t1;
