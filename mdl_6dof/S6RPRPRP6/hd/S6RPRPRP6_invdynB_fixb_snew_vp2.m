% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:52:50
% EndTime: 2019-05-05 17:52:59
% DurationCPUTime: 5.45s
% Computational Cost: add. (50510->329), mult. (121114->375), div. (0->0), fcn. (84656->8), ass. (0->137)
t805 = -2 * qJD(4);
t804 = Ifges(4,1) + Ifges(5,2);
t803 = Ifges(6,1) + Ifges(7,1);
t793 = Ifges(6,4) + Ifges(7,4);
t792 = Ifges(4,5) - Ifges(5,4);
t791 = Ifges(6,5) + Ifges(7,5);
t802 = Ifges(4,2) + Ifges(5,3);
t801 = Ifges(6,2) + Ifges(7,2);
t790 = Ifges(4,6) - Ifges(5,5);
t789 = -Ifges(5,6) - Ifges(4,4);
t788 = Ifges(6,6) + Ifges(7,6);
t800 = -Ifges(4,3) - Ifges(5,1);
t799 = Ifges(6,3) + Ifges(7,3);
t742 = sin(qJ(3));
t740 = cos(pkin(9));
t796 = cos(qJ(3));
t765 = t740 * t796;
t739 = sin(pkin(9));
t770 = qJDD(1) * t739;
t754 = t739 * t796 + t740 * t742;
t720 = t754 * qJD(1);
t774 = qJD(3) * t720;
t699 = -qJDD(1) * t765 + t742 * t770 + t774;
t775 = qJD(1) * t739;
t719 = -qJD(1) * t765 + t742 * t775;
t741 = sin(qJ(5));
t744 = cos(qJ(5));
t705 = qJD(3) * t744 + t719 * t741;
t662 = -qJD(5) * t705 - qJDD(3) * t741 + t699 * t744;
t798 = (mrSges(6,1) + mrSges(7,1)) * t662;
t747 = qJD(1) ^ 2;
t743 = sin(qJ(1));
t745 = cos(qJ(1));
t725 = -g(1) * t745 - g(2) * t743;
t721 = -pkin(1) * t747 + qJDD(1) * qJ(2) + t725;
t772 = qJD(1) * qJD(2);
t764 = -t740 * g(3) - 0.2e1 * t739 * t772;
t795 = pkin(2) * t740;
t678 = (-pkin(7) * qJDD(1) + t747 * t795 - t721) * t739 + t764;
t702 = -g(3) * t739 + (t721 + 0.2e1 * t772) * t740;
t769 = qJDD(1) * t740;
t736 = t740 ^ 2;
t785 = t736 * t747;
t679 = -pkin(2) * t785 + pkin(7) * t769 + t702;
t649 = t742 * t678 + t796 * t679;
t690 = pkin(3) * t719 - qJ(4) * t720;
t746 = qJD(3) ^ 2;
t646 = pkin(3) * t746 - qJDD(3) * qJ(4) + qJD(3) * t805 + t719 * t690 - t649;
t704 = -qJD(3) * t741 + t719 * t744;
t716 = qJD(5) + t720;
t671 = -mrSges(7,2) * t716 + mrSges(7,3) * t704;
t672 = -mrSges(6,2) * t716 + mrSges(6,3) * t704;
t797 = -t798 - (t671 + t672) * t704;
t787 = mrSges(3,2) * t739;
t786 = t704 * t671;
t648 = t678 * t796 - t742 * t679;
t691 = mrSges(4,1) * t719 + mrSges(4,2) * t720;
t773 = t719 * qJD(3);
t700 = qJDD(1) * t754 - t773;
t708 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t719;
t710 = mrSges(5,1) * t719 - qJD(3) * mrSges(5,3);
t712 = pkin(4) * t720 - qJD(3) * pkin(8);
t718 = t719 ^ 2;
t735 = t739 ^ 2;
t724 = t743 * g(1) - t745 * g(2);
t760 = qJDD(2) - t724;
t698 = (-pkin(1) - t795) * qJDD(1) + (-qJ(2) + (-t735 - t736) * pkin(7)) * t747 + t760;
t748 = pkin(3) * t774 + t720 * t805 + (-t700 + t773) * qJ(4) + t698;
t638 = -t718 * pkin(4) - t720 * t712 + (pkin(3) + pkin(8)) * t699 + t748;
t647 = -qJDD(3) * pkin(3) - t746 * qJ(4) + t720 * t690 + qJDD(4) - t648;
t641 = (t719 * t720 - qJDD(3)) * pkin(8) + (t700 + t773) * pkin(4) + t647;
t633 = -t741 * t638 + t744 * t641;
t663 = qJD(5) * t704 + qJDD(3) * t744 + t699 * t741;
t665 = -mrSges(7,1) * t704 + mrSges(7,2) * t705;
t666 = -mrSges(6,1) * t704 + mrSges(6,2) * t705;
t697 = qJDD(5) + t700;
t630 = -0.2e1 * qJD(6) * t705 + (t704 * t716 - t663) * qJ(6) + (t704 * t705 + t697) * pkin(5) + t633;
t768 = m(7) * t630 + t697 * mrSges(7,1) + t716 * t671;
t622 = m(6) * t633 + t697 * mrSges(6,1) + t716 * t672 + (-t665 - t666) * t705 + (-mrSges(6,3) - mrSges(7,3)) * t663 + t768;
t634 = t744 * t638 + t741 * t641;
t674 = mrSges(7,1) * t716 - mrSges(7,3) * t705;
t675 = mrSges(6,1) * t716 - mrSges(6,3) * t705;
t673 = pkin(5) * t716 - qJ(6) * t705;
t703 = t704 ^ 2;
t632 = -pkin(5) * t703 + qJ(6) * t662 + 0.2e1 * qJD(6) * t704 - t673 * t716 + t634;
t767 = m(7) * t632 + t662 * mrSges(7,3) + t704 * t665;
t624 = m(6) * t634 + t662 * mrSges(6,3) + t704 * t666 + (-t674 - t675) * t716 + (-mrSges(6,2) - mrSges(7,2)) * t697 + t767;
t620 = t744 * t622 + t741 * t624;
t692 = -mrSges(5,2) * t719 - mrSges(5,3) * t720;
t752 = -m(5) * t647 - t700 * mrSges(5,1) - t720 * t692 - t620;
t617 = m(4) * t648 - t700 * mrSges(4,3) - t720 * t691 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t708 - t710) * qJD(3) + t752;
t709 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t720;
t711 = mrSges(5,1) * t720 + qJD(3) * mrSges(5,2);
t643 = -pkin(4) * t699 - pkin(8) * t718 + qJD(3) * t712 - t646;
t636 = -pkin(5) * t662 - qJ(6) * t703 + t673 * t705 + qJDD(6) + t643;
t766 = m(7) * t636 + t663 * mrSges(7,2) + t705 * t674;
t755 = -m(6) * t643 - t663 * mrSges(6,2) - t705 * t675 - t766;
t751 = -m(5) * t646 + qJDD(3) * mrSges(5,3) + qJD(3) * t711 - t755;
t627 = t751 + (-mrSges(4,3) - mrSges(5,1)) * t699 + (-t691 - t692) * t719 - qJDD(3) * mrSges(4,2) - qJD(3) * t709 + m(4) * t649 + t797;
t612 = t796 * t617 + t742 * t627;
t701 = -t739 * t721 + t764;
t756 = mrSges(3,3) * qJDD(1) + t747 * (-mrSges(3,1) * t740 + t787);
t610 = m(3) * t701 - t739 * t756 + t612;
t761 = -t742 * t617 + t796 * t627;
t611 = m(3) * t702 + t740 * t756 + t761;
t762 = -t610 * t739 + t740 * t611;
t605 = m(2) * t725 - mrSges(2,1) * t747 - qJDD(1) * mrSges(2,2) + t762;
t717 = -qJDD(1) * pkin(1) - t747 * qJ(2) + t760;
t645 = t699 * pkin(3) + t748;
t783 = -t741 * t622 + t744 * t624;
t618 = m(5) * t645 - t699 * mrSges(5,2) - t700 * mrSges(5,3) - t719 * t710 - t720 * t711 + t783;
t750 = m(4) * t698 + t699 * mrSges(4,1) + t700 * mrSges(4,2) + t719 * t708 + t720 * t709 + t618;
t749 = -m(3) * t717 + mrSges(3,1) * t769 - t750 + (t735 * t747 + t785) * mrSges(3,3);
t615 = t749 - t747 * mrSges(2,2) + m(2) * t724 + (mrSges(2,1) - t787) * qJDD(1);
t784 = t743 * t605 + t745 * t615;
t606 = t740 * t610 + t739 * t611;
t782 = -t788 * t704 - t791 * t705 - t799 * t716;
t781 = t801 * t704 + t793 * t705 + t788 * t716;
t780 = -t793 * t704 - t803 * t705 - t791 * t716;
t778 = -t790 * qJD(3) + t802 * t719 + t789 * t720;
t777 = t800 * qJD(3) + t790 * t719 - t792 * t720;
t776 = t792 * qJD(3) + t789 * t719 + t804 * t720;
t763 = t745 * t605 - t615 * t743;
t759 = Ifges(3,1) * t739 + Ifges(3,4) * t740;
t758 = Ifges(3,4) * t739 + Ifges(3,2) * t740;
t757 = Ifges(3,5) * t739 + Ifges(3,6) * t740;
t723 = t757 * qJD(1);
t628 = -t663 * mrSges(7,3) - t705 * t665 + t768;
t619 = mrSges(6,2) * t643 + mrSges(7,2) * t636 - mrSges(6,3) * t633 - mrSges(7,3) * t630 - qJ(6) * t628 + t793 * t662 + t803 * t663 + t791 * t697 - t782 * t704 - t781 * t716;
t613 = -mrSges(6,1) * t643 + mrSges(6,3) * t634 - mrSges(7,1) * t636 + mrSges(7,3) * t632 - pkin(5) * (t766 - t786) + qJ(6) * t767 + (-qJ(6) * t674 - t780) * t716 + t782 * t705 + (-mrSges(7,2) * qJ(6) + t788) * t697 + t793 * t663 + (mrSges(7,1) * pkin(5) + t801) * t662;
t602 = mrSges(5,1) * t647 + mrSges(6,1) * t633 + mrSges(7,1) * t630 + mrSges(4,2) * t698 - mrSges(6,2) * t634 - mrSges(7,2) * t632 - mrSges(4,3) * t648 - mrSges(5,3) * t645 + pkin(4) * t620 + pkin(5) * t628 - qJ(4) * t618 + t777 * t719 + t781 * t705 + t780 * t704 + t804 * t700 + t789 * t699 + t799 * t697 + t791 * t663 + t788 * t662 + t792 * qJDD(3) + t778 * qJD(3);
t601 = -mrSges(4,1) * t698 + mrSges(4,3) * t649 - mrSges(5,1) * t646 + mrSges(5,2) * t645 - t741 * t619 - t744 * t613 - pkin(4) * (t755 - t797) - pkin(8) * t783 - pkin(3) * t618 + t777 * t720 - t789 * t700 - t802 * t699 + t790 * qJDD(3) + t776 * qJD(3);
t600 = mrSges(2,1) * g(3) + t741 * t613 - pkin(3) * (-qJD(3) * t710 + t752) - t744 * t619 - qJ(4) * (-t704 * t672 + t751 - t786 - t798) + mrSges(2,3) * t725 - mrSges(3,1) * t701 + mrSges(3,2) * t702 + mrSges(5,3) * t646 - mrSges(5,2) * t647 - mrSges(4,1) * t648 + mrSges(4,2) * t649 + pkin(8) * t620 - pkin(2) * t612 - pkin(1) * t606 + t778 * t720 + (qJ(4) * t692 - t776) * t719 - t792 * t700 + (mrSges(5,1) * qJ(4) + t790) * t699 + (mrSges(5,2) * pkin(3) + t800) * qJDD(3) + (Ifges(2,6) - t757) * qJDD(1) + (-t739 * t758 + t740 * t759 + Ifges(2,5)) * t747;
t599 = t740 * qJD(1) * t723 + mrSges(3,2) * t717 - mrSges(3,3) * t701 - pkin(7) * t612 + qJDD(1) * t759 - t742 * t601 + t602 * t796;
t598 = -mrSges(3,1) * t717 + mrSges(3,3) * t702 - pkin(2) * t750 + pkin(7) * t761 + qJDD(1) * t758 + t601 * t796 + t742 * t602 - t723 * t775;
t597 = -mrSges(2,2) * g(3) - mrSges(2,3) * t724 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t747 - qJ(2) * t606 - t598 * t739 + t599 * t740;
t1 = [-m(1) * g(1) + t763; -m(1) * g(2) + t784; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t784 + t745 * t597 - t743 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t763 + t743 * t597 + t745 * t600; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t724 - mrSges(2,2) * t725 + t739 * t599 + t740 * t598 + pkin(1) * (-mrSges(3,2) * t770 + t749) + qJ(2) * t762;];
tauB  = t1;
