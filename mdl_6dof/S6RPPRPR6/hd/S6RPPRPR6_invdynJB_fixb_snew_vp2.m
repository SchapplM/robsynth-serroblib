% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-05-05 14:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:26:40
% EndTime: 2019-05-05 14:26:43
% DurationCPUTime: 1.76s
% Computational Cost: add. (13139->281), mult. (24977->315), div. (0->0), fcn. (10889->6), ass. (0->114)
t825 = 2 * qJD(1);
t824 = Ifges(5,5) - Ifges(6,4);
t823 = Ifges(5,6) - Ifges(6,5);
t822 = (Ifges(5,3) + Ifges(6,1));
t775 = sin(qJ(1));
t778 = cos(qJ(1));
t746 = -t778 * g(1) - t775 * g(2);
t821 = qJDD(1) * qJ(2) + (qJD(2) * t825) + t746;
t745 = t775 * g(1) - t778 * g(2);
t780 = qJD(1) ^ 2;
t708 = -qJDD(1) * pkin(1) - t780 * qJ(2) + qJDD(2) - t745;
t793 = qJDD(1) * qJ(3) + (qJD(3) * t825) - t708;
t820 = -2 * qJD(5);
t819 = -m(3) - m(4);
t818 = pkin(4) + pkin(8);
t817 = pkin(8) * t780;
t774 = sin(qJ(4));
t816 = t774 * g(3);
t815 = t780 * pkin(7);
t814 = mrSges(2,1) - mrSges(3,2);
t813 = Ifges(5,4) + Ifges(6,6);
t812 = mrSges(4,3) * t780;
t701 = qJDD(3) + (-pkin(1) - qJ(3)) * t780 + t821;
t697 = -qJDD(1) * pkin(7) + t701;
t777 = cos(qJ(4));
t811 = t777 * t697;
t706 = pkin(1) * t780 - t821;
t691 = t811 + t816;
t806 = qJD(1) * qJD(4);
t751 = t774 * t806;
t737 = qJDD(1) * t777 - t751;
t808 = qJD(1) * t774;
t740 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t808;
t742 = mrSges(6,1) * t808 - (qJD(4) * mrSges(6,3));
t801 = t777 * t806;
t736 = qJDD(1) * t774 + t801;
t807 = qJD(1) * t777;
t744 = pkin(5) * t807 - (qJD(4) * pkin(8));
t768 = t774 ^ 2;
t784 = pkin(4) * t801 + t807 * t820 + t793 + (-t737 + t751) * qJ(5);
t671 = t784 + t818 * t736 + (-pkin(5) * t768 - pkin(7)) * t780 - t744 * t807;
t733 = (pkin(4) * t774 - qJ(5) * t777) * qJD(1);
t779 = qJD(4) ^ 2;
t792 = -t779 * qJ(5) + t733 * t807 + qJDD(5) - t811;
t674 = t737 * pkin(5) - t818 * qJDD(4) + (pkin(5) * t806 + t777 * t817 - g(3)) * t774 + t792;
t773 = sin(qJ(6));
t776 = cos(qJ(6));
t669 = -t671 * t773 + t674 * t776;
t731 = -qJD(4) * t773 + t776 * t808;
t690 = qJD(6) * t731 + qJDD(4) * t776 + t736 * t773;
t732 = qJD(4) * t776 + t773 * t808;
t694 = -mrSges(7,1) * t731 + mrSges(7,2) * t732;
t749 = qJD(6) + t807;
t702 = -mrSges(7,2) * t749 + mrSges(7,3) * t731;
t730 = qJDD(6) + t737;
t666 = m(7) * t669 + mrSges(7,1) * t730 - mrSges(7,3) * t690 - t694 * t732 + t702 * t749;
t670 = t671 * t776 + t674 * t773;
t689 = -qJD(6) * t732 - qJDD(4) * t773 + t736 * t776;
t703 = mrSges(7,1) * t749 - mrSges(7,3) * t732;
t667 = m(7) * t670 - mrSges(7,2) * t730 + mrSges(7,3) * t689 + t694 * t731 - t703 * t749;
t656 = t776 * t666 + t773 * t667;
t679 = -qJDD(4) * pkin(4) + t792 - t816;
t789 = -m(6) * t679 - t737 * mrSges(6,1) - t656;
t734 = (-mrSges(6,2) * t774 - mrSges(6,3) * t777) * qJD(1);
t797 = qJD(1) * (-t734 - (mrSges(5,1) * t774 + mrSges(5,2) * t777) * qJD(1));
t652 = m(5) * t691 - t737 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t740 - t742) * qJD(4) + t777 * t797 + t789;
t692 = -t777 * g(3) + t774 * t697;
t741 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t807;
t786 = -t779 * pkin(4) + qJDD(4) * qJ(5) - t733 * t808 + t692;
t677 = (qJD(4) * t820) - t786;
t743 = mrSges(6,1) * t807 + qJD(4) * mrSges(6,2);
t673 = -t768 * t817 - t736 * pkin(5) + ((2 * qJD(5)) + t744) * qJD(4) + t786;
t791 = -m(7) * t673 + t689 * mrSges(7,1) - t690 * mrSges(7,2) + t731 * t702 - t732 * t703;
t785 = -m(6) * t677 + qJDD(4) * mrSges(6,3) + qJD(4) * t743 - t791;
t662 = m(5) * t692 - qJDD(4) * mrSges(5,2) - qJD(4) * t741 + (-mrSges(5,3) - mrSges(6,1)) * t736 + t774 * t797 + t785;
t645 = t777 * t652 + t774 * t662;
t796 = m(4) * t701 + qJDD(1) * mrSges(4,2) + t645;
t790 = -m(3) * t706 + (t780 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t796;
t641 = m(2) * t746 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t780) + t790;
t676 = t736 * pkin(4) + t784 - t815;
t809 = -t773 * t666 + t776 * t667;
t653 = m(6) * t676 - t736 * mrSges(6,2) - t737 * mrSges(6,3) - t742 * t808 - t743 * t807 + t809;
t696 = t793 - t815;
t787 = -m(5) * t696 - t736 * mrSges(5,1) - t737 * mrSges(5,2) - t740 * t808 - t741 * t807 - t653;
t650 = -m(4) * t793 - t780 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t787;
t783 = -m(3) * t708 + t780 * mrSges(3,3) - t650;
t647 = m(2) * t745 - t780 * mrSges(2,2) + t814 * qJDD(1) + t783;
t810 = t775 * t641 + t778 * t647;
t803 = (Ifges(2,5) - Ifges(3,4) + Ifges(4,5));
t802 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t800 = t778 * t641 - t647 * t775;
t799 = -t774 * t652 + t777 * t662;
t798 = qJD(1) * (-(t822 * qJD(4)) + (t823 * t774 - t824 * t777) * qJD(1));
t682 = Ifges(7,4) * t732 + Ifges(7,2) * t731 + Ifges(7,6) * t749;
t683 = Ifges(7,1) * t732 + Ifges(7,4) * t731 + Ifges(7,5) * t749;
t788 = mrSges(7,1) * t669 - mrSges(7,2) * t670 + Ifges(7,5) * t690 + Ifges(7,6) * t689 + Ifges(7,3) * t730 + t732 * t682 - t731 * t683;
t681 = Ifges(7,5) * t732 + Ifges(7,6) * t731 + Ifges(7,3) * t749;
t658 = -mrSges(7,1) * t673 + mrSges(7,3) * t670 + Ifges(7,4) * t690 + Ifges(7,2) * t689 + Ifges(7,6) * t730 - t681 * t732 + t683 * t749;
t659 = mrSges(7,2) * t673 - mrSges(7,3) * t669 + Ifges(7,1) * t690 + Ifges(7,4) * t689 + Ifges(7,5) * t730 + t681 * t731 - t682 * t749;
t712 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t777 - Ifges(5,4) * t774) * qJD(1);
t714 = (Ifges(6,4) * qJD(4)) + (-Ifges(6,2) * t777 + Ifges(6,6) * t774) * qJD(1);
t636 = -mrSges(5,1) * t696 + mrSges(5,3) * t692 - mrSges(6,1) * t677 + mrSges(6,2) * t676 - t773 * t659 - t776 * t658 - pkin(5) * t791 - pkin(8) * t809 - pkin(4) * t653 + t813 * t737 + (-Ifges(5,2) - Ifges(6,3)) * t736 + t823 * qJDD(4) + (t712 - t714) * qJD(4) + t777 * t798;
t711 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t777 - Ifges(5,2) * t774) * qJD(1);
t713 = (Ifges(6,5) * qJD(4)) + (-Ifges(6,6) * t777 + Ifges(6,3) * t774) * qJD(1);
t638 = (-t711 + t713) * qJD(4) - t813 * t736 + (Ifges(5,1) + Ifges(6,2)) * t737 + t788 + t774 * t798 - mrSges(5,3) * t691 + mrSges(5,2) * t696 - mrSges(6,3) * t676 + mrSges(6,1) * t679 + pkin(5) * t656 - qJ(5) * t653 + t824 * qJDD(4);
t649 = qJDD(1) * mrSges(3,2) - t783;
t782 = -mrSges(2,2) * t746 - mrSges(3,3) * t706 + mrSges(4,3) * t793 - pkin(7) * t645 - qJ(3) * t650 - t636 * t774 + t777 * t638 + qJ(2) * (t790 - t812) - pkin(1) * t649 + mrSges(4,2) * t701 + mrSges(3,2) * t708 + mrSges(2,1) * t745 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t655 = qJDD(4) * mrSges(6,2) + qJD(4) * t742 + t734 * t807 - t789;
t781 = -mrSges(5,2) * t692 - mrSges(6,3) * t677 - pkin(8) * t656 - t773 * t658 - pkin(4) * t655 + t776 * t659 + qJ(5) * (-t734 * t808 + t785) + mrSges(6,2) * t679 + mrSges(5,1) * t691 + t712 * t808 + t711 * t807 + (-t713 * t777 - t714 * t774) * qJD(1) + t824 * t737 + (-mrSges(6,1) * qJ(5) - t823) * t736 + t822 * qJDD(4);
t644 = t819 * g(3) + t799;
t643 = t796 - t812;
t635 = -qJ(3) * t799 + t781 + (t803 * t780) + mrSges(2,3) * t746 - mrSges(3,1) * t706 + mrSges(4,1) * t701 + (m(4) * qJ(3) + mrSges(4,3) + t814) * g(3) + pkin(2) * t643 - pkin(1) * t644 + pkin(3) * t645 + t802 * qJDD(1);
t634 = -qJ(2) * t644 - mrSges(2,3) * t745 + pkin(2) * t650 + mrSges(3,1) * t708 + t777 * t636 + pkin(3) * t787 + pkin(7) * t799 + t774 * t638 - mrSges(4,1) * t793 - t802 * t780 + t803 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t800; -m(1) * g(2) + t810; (-m(1) - m(2) + t819) * g(3) + t799; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t810 + t778 * t634 - t775 * t635; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t800 + t775 * t634 + t778 * t635; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t782; t782; t649; t643; t781; t655; t788;];
tauJB  = t1;
