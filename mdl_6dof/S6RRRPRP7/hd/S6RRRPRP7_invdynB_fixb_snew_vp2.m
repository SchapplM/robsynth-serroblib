% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:14:59
% EndTime: 2019-05-07 08:15:37
% DurationCPUTime: 24.96s
% Computational Cost: add. (403997->374), mult. (890330->475), div. (0->0), fcn. (710473->12), ass. (0->153)
t797 = Ifges(6,1) + Ifges(7,1);
t789 = Ifges(6,4) - Ifges(7,5);
t796 = -Ifges(6,5) - Ifges(7,4);
t795 = Ifges(6,2) + Ifges(7,3);
t787 = Ifges(6,6) - Ifges(7,6);
t794 = -Ifges(6,3) - Ifges(7,2);
t748 = sin(pkin(6));
t753 = sin(qJ(2));
t756 = cos(qJ(2));
t773 = qJD(1) * qJD(2);
t736 = (-qJDD(1) * t756 + t753 * t773) * t748;
t776 = qJD(1) * t748;
t734 = (-pkin(2) * t756 - pkin(9) * t753) * t776;
t750 = cos(pkin(6));
t744 = qJD(1) * t750 + qJD(2);
t742 = t744 ^ 2;
t743 = qJDD(1) * t750 + qJDD(2);
t775 = qJD(1) * t756;
t754 = sin(qJ(1));
t757 = cos(qJ(1));
t740 = t754 * g(1) - g(2) * t757;
t758 = qJD(1) ^ 2;
t792 = pkin(8) * t748;
t731 = qJDD(1) * pkin(1) + t758 * t792 + t740;
t741 = -g(1) * t757 - g(2) * t754;
t732 = -pkin(1) * t758 + qJDD(1) * t792 + t741;
t784 = t750 * t753;
t777 = t731 * t784 + t756 * t732;
t690 = -pkin(2) * t742 + pkin(9) * t743 + (-g(3) * t753 + t734 * t775) * t748 + t777;
t735 = (qJDD(1) * t753 + t756 * t773) * t748;
t791 = g(3) * t750;
t691 = pkin(2) * t736 - pkin(9) * t735 - t791 + (-t731 + (pkin(2) * t753 - pkin(9) * t756) * t744 * qJD(1)) * t748;
t752 = sin(qJ(3));
t755 = cos(qJ(3));
t654 = -t690 * t752 + t755 * t691;
t771 = t753 * t776;
t724 = t744 * t755 - t752 * t771;
t703 = qJD(3) * t724 + t735 * t755 + t743 * t752;
t725 = t744 * t752 + t755 * t771;
t728 = qJDD(3) + t736;
t770 = t748 * t775;
t739 = qJD(3) - t770;
t645 = (t724 * t739 - t703) * qJ(4) + (t724 * t725 + t728) * pkin(3) + t654;
t655 = t755 * t690 + t752 * t691;
t702 = -qJD(3) * t725 - t735 * t752 + t743 * t755;
t714 = pkin(3) * t739 - qJ(4) * t725;
t723 = t724 ^ 2;
t648 = -pkin(3) * t723 + qJ(4) * t702 - t714 * t739 + t655;
t747 = sin(pkin(11));
t749 = cos(pkin(11));
t711 = t724 * t747 + t725 * t749;
t640 = -0.2e1 * qJD(4) * t711 + t645 * t749 - t747 * t648;
t793 = cos(qJ(5));
t790 = -mrSges(6,3) - mrSges(7,2);
t786 = t748 * t753;
t785 = t748 * t756;
t783 = t750 * t756;
t705 = -g(3) * t786 + t777;
t729 = mrSges(3,1) * t744 - mrSges(3,3) * t771;
t733 = (-mrSges(3,1) * t756 + mrSges(3,2) * t753) * t776;
t710 = t724 * t749 - t725 * t747;
t641 = 0.2e1 * qJD(4) * t710 + t747 * t645 + t749 * t648;
t678 = t702 * t749 - t703 * t747;
t684 = -mrSges(5,1) * t710 + mrSges(5,2) * t711;
t695 = mrSges(5,1) * t739 - mrSges(5,3) * t711;
t685 = -pkin(4) * t710 - pkin(10) * t711;
t738 = t739 ^ 2;
t639 = -pkin(4) * t738 + pkin(10) * t728 + t685 * t710 + t641;
t704 = -g(3) * t785 + t731 * t783 - t753 * t732;
t689 = -pkin(2) * t743 - pkin(9) * t742 + t734 * t771 - t704;
t649 = -pkin(3) * t702 - qJ(4) * t723 + t725 * t714 + qJDD(4) + t689;
t679 = t702 * t747 + t703 * t749;
t643 = (-t710 * t739 - t679) * pkin(10) + (t711 * t739 - t678) * pkin(4) + t649;
t751 = sin(qJ(5));
t636 = t793 * t639 + t751 * t643;
t693 = t711 * t793 + t751 * t739;
t652 = qJD(5) * t693 + t679 * t751 - t728 * t793;
t709 = qJD(5) - t710;
t672 = mrSges(6,1) * t709 - mrSges(6,3) * t693;
t677 = qJDD(5) - t678;
t692 = t711 * t751 - t739 * t793;
t667 = pkin(5) * t692 - qJ(6) * t693;
t708 = t709 ^ 2;
t632 = -pkin(5) * t708 + qJ(6) * t677 + 0.2e1 * qJD(6) * t709 - t667 * t692 + t636;
t673 = -mrSges(7,1) * t709 + mrSges(7,2) * t693;
t772 = m(7) * t632 + t677 * mrSges(7,3) + t709 * t673;
t668 = mrSges(7,1) * t692 - mrSges(7,3) * t693;
t778 = -mrSges(6,1) * t692 - mrSges(6,2) * t693 - t668;
t627 = m(6) * t636 - mrSges(6,2) * t677 + t652 * t790 - t672 * t709 + t692 * t778 + t772;
t635 = -t751 * t639 + t643 * t793;
t653 = -t692 * qJD(5) + t679 * t793 + t751 * t728;
t671 = -mrSges(6,2) * t709 - mrSges(6,3) * t692;
t633 = -t677 * pkin(5) - t708 * qJ(6) + t693 * t667 + qJDD(6) - t635;
t670 = -mrSges(7,2) * t692 + mrSges(7,3) * t709;
t764 = -m(7) * t633 + t677 * mrSges(7,1) + t709 * t670;
t629 = m(6) * t635 + mrSges(6,1) * t677 + t653 * t790 + t671 * t709 + t693 * t778 + t764;
t766 = t793 * t627 - t751 * t629;
t619 = m(5) * t641 - mrSges(5,2) * t728 + mrSges(5,3) * t678 + t684 * t710 - t695 * t739 + t766;
t694 = -mrSges(5,2) * t739 + mrSges(5,3) * t710;
t638 = -pkin(4) * t728 - pkin(10) * t738 + t711 * t685 - t640;
t634 = -0.2e1 * qJD(6) * t693 + (t692 * t709 - t653) * qJ(6) + (t693 * t709 + t652) * pkin(5) + t638;
t630 = m(7) * t634 + t652 * mrSges(7,1) - t653 * mrSges(7,3) + t692 * t670 - t693 * t673;
t760 = -m(6) * t638 - t652 * mrSges(6,1) - t653 * mrSges(6,2) - t692 * t671 - t693 * t672 - t630;
t624 = m(5) * t640 + mrSges(5,1) * t728 - mrSges(5,3) * t679 - t684 * t711 + t694 * t739 + t760;
t613 = t747 * t619 + t749 * t624;
t712 = -mrSges(4,1) * t724 + mrSges(4,2) * t725;
t713 = -mrSges(4,2) * t739 + mrSges(4,3) * t724;
t611 = m(4) * t654 + mrSges(4,1) * t728 - mrSges(4,3) * t703 - t712 * t725 + t713 * t739 + t613;
t715 = mrSges(4,1) * t739 - mrSges(4,3) * t725;
t767 = t749 * t619 - t624 * t747;
t612 = m(4) * t655 - mrSges(4,2) * t728 + mrSges(4,3) * t702 + t712 * t724 - t715 * t739 + t767;
t768 = -t611 * t752 + t755 * t612;
t602 = m(3) * t705 - mrSges(3,2) * t743 - mrSges(3,3) * t736 - t729 * t744 + t733 * t770 + t768;
t605 = t755 * t611 + t752 * t612;
t719 = -t731 * t748 - t791;
t730 = -mrSges(3,2) * t744 + mrSges(3,3) * t770;
t604 = m(3) * t719 + mrSges(3,1) * t736 + mrSges(3,2) * t735 + (t729 * t753 - t730 * t756) * t776 + t605;
t622 = t751 * t627 + t793 * t629;
t761 = m(5) * t649 - t678 * mrSges(5,1) + mrSges(5,2) * t679 - t710 * t694 + t695 * t711 + t622;
t759 = -m(4) * t689 + t702 * mrSges(4,1) - mrSges(4,2) * t703 + t724 * t713 - t715 * t725 - t761;
t616 = m(3) * t704 + mrSges(3,1) * t743 - mrSges(3,3) * t735 + t730 * t744 - t733 * t771 + t759;
t593 = t602 * t784 - t604 * t748 + t616 * t783;
t591 = m(2) * t740 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t758 + t593;
t598 = t756 * t602 - t616 * t753;
t597 = m(2) * t741 - mrSges(2,1) * t758 - qJDD(1) * mrSges(2,2) + t598;
t782 = t757 * t591 + t754 * t597;
t781 = t795 * t692 - t789 * t693 - t787 * t709;
t780 = t787 * t692 + t796 * t693 + t794 * t709;
t779 = -t789 * t692 + t797 * t693 - t796 * t709;
t592 = t602 * t786 + t750 * t604 + t616 * t785;
t769 = -t591 * t754 + t757 * t597;
t620 = -mrSges(6,1) * t638 - mrSges(7,1) * t634 + mrSges(7,2) * t632 + mrSges(6,3) * t636 - pkin(5) * t630 - t795 * t652 + t789 * t653 + t787 * t677 + t780 * t693 + t779 * t709;
t621 = mrSges(6,2) * t638 + mrSges(7,2) * t633 - mrSges(6,3) * t635 - mrSges(7,3) * t634 - qJ(6) * t630 - t789 * t652 + t797 * t653 - t677 * t796 + t780 * t692 + t781 * t709;
t680 = Ifges(5,5) * t711 + Ifges(5,6) * t710 + Ifges(5,3) * t739;
t681 = Ifges(5,4) * t711 + Ifges(5,2) * t710 + Ifges(5,6) * t739;
t606 = mrSges(5,2) * t649 - mrSges(5,3) * t640 + Ifges(5,1) * t679 + Ifges(5,4) * t678 + Ifges(5,5) * t728 - pkin(10) * t622 - t751 * t620 + t621 * t793 + t710 * t680 - t739 * t681;
t682 = Ifges(5,1) * t711 + Ifges(5,4) * t710 + Ifges(5,5) * t739;
t607 = Ifges(5,4) * t679 + Ifges(5,2) * t678 + Ifges(5,6) * t728 - t711 * t680 + t739 * t682 - mrSges(5,1) * t649 + mrSges(5,3) * t641 - mrSges(6,1) * t635 + mrSges(6,2) * t636 + mrSges(7,1) * t633 - mrSges(7,3) * t632 - pkin(5) * t764 - qJ(6) * t772 - pkin(4) * t622 + (pkin(5) * t668 + t781) * t693 + (qJ(6) * t668 - t779) * t692 + t794 * t677 + (mrSges(7,2) * pkin(5) + t796) * t653 + (mrSges(7,2) * qJ(6) + t787) * t652;
t696 = Ifges(4,5) * t725 + Ifges(4,6) * t724 + Ifges(4,3) * t739;
t698 = Ifges(4,1) * t725 + Ifges(4,4) * t724 + Ifges(4,5) * t739;
t589 = -mrSges(4,1) * t689 + mrSges(4,3) * t655 + Ifges(4,4) * t703 + Ifges(4,2) * t702 + Ifges(4,6) * t728 - pkin(3) * t761 + qJ(4) * t767 + t747 * t606 + t749 * t607 - t725 * t696 + t739 * t698;
t697 = Ifges(4,4) * t725 + Ifges(4,2) * t724 + Ifges(4,6) * t739;
t594 = mrSges(4,2) * t689 - mrSges(4,3) * t654 + Ifges(4,1) * t703 + Ifges(4,4) * t702 + Ifges(4,5) * t728 - qJ(4) * t613 + t606 * t749 - t607 * t747 + t696 * t724 - t697 * t739;
t716 = Ifges(3,3) * t744 + (Ifges(3,5) * t753 + Ifges(3,6) * t756) * t776;
t717 = Ifges(3,6) * t744 + (Ifges(3,4) * t753 + Ifges(3,2) * t756) * t776;
t587 = mrSges(3,2) * t719 - mrSges(3,3) * t704 + Ifges(3,1) * t735 - Ifges(3,4) * t736 + Ifges(3,5) * t743 - pkin(9) * t605 - t589 * t752 + t594 * t755 + t716 * t770 - t717 * t744;
t718 = Ifges(3,5) * t744 + (Ifges(3,1) * t753 + Ifges(3,4) * t756) * t776;
t588 = -t716 * t771 - t793 * t620 - pkin(4) * t760 + (-Ifges(4,3) - Ifges(5,3)) * t728 - pkin(3) * t613 - pkin(2) * t605 - pkin(10) * t766 - t751 * t621 + Ifges(3,6) * t743 + t744 * t718 + Ifges(3,4) * t735 - Ifges(3,2) * t736 - t725 * t697 - mrSges(3,1) * t719 + t724 * t698 + t710 * t682 - t711 * t681 + mrSges(3,3) * t705 - Ifges(4,6) * t702 - Ifges(4,5) * t703 - Ifges(5,6) * t678 - Ifges(5,5) * t679 - mrSges(4,1) * t654 + mrSges(4,2) * t655 - mrSges(5,1) * t640 + mrSges(5,2) * t641;
t762 = pkin(8) * t598 + t587 * t753 + t588 * t756;
t586 = Ifges(3,5) * t735 - Ifges(3,6) * t736 + Ifges(3,3) * t743 + mrSges(3,1) * t704 - mrSges(3,2) * t705 + t752 * t594 + t755 * t589 + pkin(2) * t759 + pkin(9) * t768 + (t717 * t753 - t718 * t756) * t776;
t585 = -mrSges(2,2) * g(3) - mrSges(2,3) * t740 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t758 + t587 * t756 - t588 * t753 + (-t592 * t748 - t593 * t750) * pkin(8);
t584 = mrSges(2,1) * g(3) + mrSges(2,3) * t741 + Ifges(2,5) * t758 + Ifges(2,6) * qJDD(1) - pkin(1) * t592 - t586 * t748 + t750 * t762;
t1 = [-m(1) * g(1) + t769; -m(1) * g(2) + t782; (-m(1) - m(2)) * g(3) + t592; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t782 - t754 * t584 + t757 * t585; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t769 + t757 * t584 + t754 * t585; -mrSges(1,1) * g(2) + mrSges(2,1) * t740 + mrSges(1,2) * g(1) - mrSges(2,2) * t741 + Ifges(2,3) * qJDD(1) + pkin(1) * t593 + t586 * t750 + t748 * t762;];
tauB  = t1;
