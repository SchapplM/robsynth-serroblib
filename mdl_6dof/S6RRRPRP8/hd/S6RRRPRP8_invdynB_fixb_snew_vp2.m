% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 08:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:29:46
% EndTime: 2019-05-07 08:29:58
% DurationCPUTime: 5.85s
% Computational Cost: add. (65296->342), mult. (128848->398), div. (0->0), fcn. (85019->8), ass. (0->134)
t802 = Ifges(4,1) + Ifges(5,1);
t801 = Ifges(6,1) + Ifges(7,1);
t792 = Ifges(4,4) - Ifges(5,5);
t791 = Ifges(6,4) + Ifges(7,4);
t790 = Ifges(4,5) + Ifges(5,4);
t789 = Ifges(6,5) + Ifges(7,5);
t800 = Ifges(4,2) + Ifges(5,3);
t799 = Ifges(6,2) + Ifges(7,2);
t788 = Ifges(4,6) - Ifges(5,6);
t787 = Ifges(6,6) + Ifges(7,6);
t798 = -Ifges(4,3) - Ifges(5,2);
t797 = Ifges(6,3) + Ifges(7,3);
t796 = 2 * qJD(4);
t795 = cos(qJ(3));
t756 = cos(qJ(2));
t776 = qJD(1) * t756;
t743 = qJD(3) - t776;
t794 = pkin(3) * t743;
t793 = -mrSges(4,3) - mrSges(5,2);
t752 = sin(qJ(3));
t753 = sin(qJ(2));
t777 = qJD(1) * t753;
t730 = -qJD(2) * t795 + t752 * t777;
t786 = t730 * t743;
t754 = sin(qJ(1));
t757 = cos(qJ(1));
t740 = -g(1) * t757 - g(2) * t754;
t759 = qJD(1) ^ 2;
t721 = -pkin(1) * t759 + qJDD(1) * pkin(7) + t740;
t709 = -g(3) * t753 + t756 * t721;
t732 = (-mrSges(3,1) * t756 + mrSges(3,2) * t753) * qJD(1);
t775 = qJD(1) * qJD(2);
t772 = t753 * t775;
t735 = t756 * qJDD(1) - t772;
t736 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t777;
t739 = g(1) * t754 - t757 * g(2);
t720 = -qJDD(1) * pkin(1) - pkin(7) * t759 - t739;
t771 = t756 * t775;
t734 = qJDD(1) * t753 + t771;
t668 = (-t734 - t771) * pkin(8) + (-t735 + t772) * pkin(2) + t720;
t733 = (-pkin(2) * t756 - pkin(8) * t753) * qJD(1);
t758 = qJD(2) ^ 2;
t673 = -pkin(2) * t758 + qJDD(2) * pkin(8) + t733 * t776 + t709;
t653 = t752 * t668 + t795 * t673;
t731 = t752 * qJD(2) + t777 * t795;
t694 = qJD(3) * t731 - qJDD(2) * t795 + t734 * t752;
t705 = mrSges(4,1) * t743 - mrSges(4,3) * t731;
t729 = qJDD(3) - t735;
t701 = pkin(3) * t730 - qJ(4) * t731;
t742 = t743 ^ 2;
t642 = -pkin(3) * t742 + t729 * qJ(4) - t730 * t701 + t743 * t796 + t653;
t706 = -mrSges(5,1) * t743 + mrSges(5,2) * t731;
t652 = t668 * t795 - t752 * t673;
t644 = -t729 * pkin(3) - t742 * qJ(4) + t731 * t701 + qJDD(4) - t652;
t695 = -t730 * qJD(3) + t752 * qJDD(2) + t734 * t795;
t636 = (-t695 - t786) * pkin(9) + (t730 * t731 - t729) * pkin(4) + t644;
t710 = -pkin(4) * t743 - pkin(9) * t731;
t728 = t730 ^ 2;
t639 = -pkin(4) * t728 + pkin(9) * t694 + t710 * t743 + t642;
t751 = sin(qJ(5));
t755 = cos(qJ(5));
t630 = t755 * t636 - t639 * t751;
t697 = t730 * t755 - t731 * t751;
t651 = qJD(5) * t697 + t694 * t751 + t695 * t755;
t698 = t730 * t751 + t731 * t755;
t665 = -mrSges(7,1) * t697 + mrSges(7,2) * t698;
t666 = -mrSges(6,1) * t697 + mrSges(6,2) * t698;
t741 = qJD(5) - t743;
t675 = -mrSges(6,2) * t741 + mrSges(6,3) * t697;
t725 = qJDD(5) - t729;
t627 = -0.2e1 * qJD(6) * t698 + (t697 * t741 - t651) * qJ(6) + (t697 * t698 + t725) * pkin(5) + t630;
t674 = -mrSges(7,2) * t741 + mrSges(7,3) * t697;
t774 = m(7) * t627 + t725 * mrSges(7,1) + t741 * t674;
t621 = m(6) * t630 + mrSges(6,1) * t725 + t675 * t741 + (-t665 - t666) * t698 + (-mrSges(6,3) - mrSges(7,3)) * t651 + t774;
t631 = t751 * t636 + t755 * t639;
t650 = -qJD(5) * t698 + t694 * t755 - t695 * t751;
t677 = mrSges(7,1) * t741 - mrSges(7,3) * t698;
t678 = mrSges(6,1) * t741 - mrSges(6,3) * t698;
t676 = pkin(5) * t741 - qJ(6) * t698;
t696 = t697 ^ 2;
t629 = -pkin(5) * t696 + qJ(6) * t650 + 0.2e1 * qJD(6) * t697 - t676 * t741 + t631;
t773 = m(7) * t629 + t650 * mrSges(7,3) + t697 * t665;
t623 = m(6) * t631 + mrSges(6,3) * t650 + t666 * t697 + (-t677 - t678) * t741 + (-mrSges(6,2) - mrSges(7,2)) * t725 + t773;
t767 = -t621 * t751 + t755 * t623;
t765 = m(5) * t642 + t729 * mrSges(5,3) + t743 * t706 + t767;
t702 = mrSges(5,1) * t730 - mrSges(5,3) * t731;
t778 = -mrSges(4,1) * t730 - mrSges(4,2) * t731 - t702;
t615 = m(4) * t653 - mrSges(4,2) * t729 + t694 * t793 - t705 * t743 + t730 * t778 + t765;
t704 = -mrSges(4,2) * t743 - mrSges(4,3) * t730;
t618 = t621 * t755 + t623 * t751;
t707 = -mrSges(5,2) * t730 + mrSges(5,3) * t743;
t763 = -m(5) * t644 + t729 * mrSges(5,1) + t743 * t707 - t618;
t616 = m(4) * t652 + mrSges(4,1) * t729 + t695 * t793 + t704 * t743 + t731 * t778 + t763;
t768 = t795 * t615 - t616 * t752;
t611 = m(3) * t709 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t735 - qJD(2) * t736 + t732 * t776 + t768;
t708 = -t756 * g(3) - t753 * t721;
t737 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t776;
t672 = -qJDD(2) * pkin(2) - t758 * pkin(8) + t733 * t777 - t708;
t764 = t694 * pkin(3) + t672 + (-t695 + t786) * qJ(4);
t643 = (-(2 * qJD(4)) + t794) * t731 + t764;
t640 = -pkin(4) * t694 - pkin(9) * t728 - t764 + (t710 - t794 + t796) * t731;
t633 = -pkin(5) * t650 - qJ(6) * t696 + t676 * t698 + qJDD(6) + t640;
t766 = m(7) * t633 - t650 * mrSges(7,1) + t651 * mrSges(7,2) - t697 * t674 + t698 * t677;
t762 = -m(6) * t640 + t650 * mrSges(6,1) - t651 * mrSges(6,2) + t697 * t675 - t698 * t678 - t766;
t624 = m(5) * t643 + t694 * mrSges(5,1) - t695 * mrSges(5,3) - t731 * t706 + t730 * t707 + t762;
t760 = -m(4) * t672 - t694 * mrSges(4,1) - t695 * mrSges(4,2) - t730 * t704 - t731 * t705 - t624;
t620 = m(3) * t708 + qJDD(2) * mrSges(3,1) - t734 * mrSges(3,3) + qJD(2) * t737 - t732 * t777 + t760;
t769 = t756 * t611 - t620 * t753;
t605 = m(2) * t740 - mrSges(2,1) * t759 - qJDD(1) * mrSges(2,2) + t769;
t612 = t752 * t615 + t616 * t795;
t761 = -m(3) * t720 + t735 * mrSges(3,1) - t734 * mrSges(3,2) - t736 * t777 + t737 * t776 - t612;
t608 = m(2) * t739 + qJDD(1) * mrSges(2,1) - t759 * mrSges(2,2) + t761;
t785 = t754 * t605 + t757 * t608;
t606 = t753 * t611 + t756 * t620;
t784 = -t787 * t697 - t789 * t698 - t797 * t741;
t783 = t799 * t697 + t791 * t698 + t787 * t741;
t782 = -t791 * t697 - t801 * t698 - t789 * t741;
t781 = t800 * t730 - t792 * t731 - t788 * t743;
t780 = t788 * t730 - t790 * t731 + t798 * t743;
t779 = -t792 * t730 + t802 * t731 + t790 * t743;
t770 = t757 * t605 - t608 * t754;
t719 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t753 + Ifges(3,4) * t756) * qJD(1);
t718 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t753 + Ifges(3,2) * t756) * qJD(1);
t717 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t753 + Ifges(3,6) * t756) * qJD(1);
t625 = -mrSges(7,3) * t651 - t665 * t698 + t774;
t617 = mrSges(6,2) * t640 + mrSges(7,2) * t633 - mrSges(6,3) * t630 - mrSges(7,3) * t627 - qJ(6) * t625 + t791 * t650 + t801 * t651 - t784 * t697 + t789 * t725 - t783 * t741;
t613 = -mrSges(6,1) * t640 + mrSges(6,3) * t631 - mrSges(7,1) * t633 + mrSges(7,3) * t629 - pkin(5) * t766 + qJ(6) * t773 + (-qJ(6) * t677 - t782) * t741 + (-mrSges(7,2) * qJ(6) + t787) * t725 + t784 * t698 + t791 * t651 + t799 * t650;
t602 = mrSges(4,2) * t672 + mrSges(5,2) * t644 - mrSges(4,3) * t652 - mrSges(5,3) * t643 - pkin(9) * t618 - qJ(4) * t624 - t613 * t751 + t617 * t755 - t792 * t694 + t802 * t695 + t790 * t729 + t780 * t730 + t781 * t743;
t601 = -mrSges(4,1) * t672 - mrSges(5,1) * t643 + mrSges(5,2) * t642 + mrSges(4,3) * t653 - pkin(3) * t624 - pkin(4) * t762 - pkin(9) * t767 - t755 * t613 - t751 * t617 - t800 * t694 + t792 * t695 + t788 * t729 + t780 * t731 + t779 * t743;
t600 = -qJ(4) * t765 - pkin(3) * t763 - t717 * t777 + (qJ(4) * t702 - t779) * t730 + (pkin(3) * t702 + t781) * t731 + t782 * t697 + t783 * t698 + t787 * t650 + (qJ(4) * mrSges(5,2) + t788) * t694 + t789 * t651 + (pkin(3) * mrSges(5,2) - t790) * t695 - mrSges(5,3) * t642 + mrSges(5,1) * t644 + mrSges(6,1) * t630 - mrSges(6,2) * t631 - mrSges(7,2) * t629 + pkin(5) * t625 + mrSges(7,1) * t627 + pkin(4) * t618 - pkin(2) * t612 + Ifges(3,6) * qJDD(2) - mrSges(4,1) * t652 + mrSges(4,2) * t653 + t797 * t725 + t798 * t729 + Ifges(3,4) * t734 + Ifges(3,2) * t735 + mrSges(3,3) * t709 + qJD(2) * t719 - mrSges(3,1) * t720;
t599 = mrSges(3,2) * t720 - mrSges(3,3) * t708 + Ifges(3,1) * t734 + Ifges(3,4) * t735 + Ifges(3,5) * qJDD(2) - pkin(8) * t612 - qJD(2) * t718 - t752 * t601 + t602 * t795 + t717 * t776;
t598 = Ifges(2,6) * qJDD(1) + t759 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t740 - Ifges(3,5) * t734 - Ifges(3,6) * t735 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t708 + mrSges(3,2) * t709 - t752 * t602 - t795 * t601 - pkin(2) * t760 - pkin(8) * t768 - pkin(1) * t606 + (-t718 * t753 + t719 * t756) * qJD(1);
t597 = -mrSges(2,2) * g(3) - mrSges(2,3) * t739 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t759 - pkin(7) * t606 + t599 * t756 - t600 * t753;
t1 = [-m(1) * g(1) + t770; -m(1) * g(2) + t785; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t785 + t757 * t597 - t754 * t598; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t770 + t754 * t597 + t757 * t598; -mrSges(1,1) * g(2) + mrSges(2,1) * t739 + mrSges(1,2) * g(1) - mrSges(2,2) * t740 + Ifges(2,3) * qJDD(1) + pkin(1) * t761 + pkin(7) * t769 + t753 * t599 + t756 * t600;];
tauB  = t1;
