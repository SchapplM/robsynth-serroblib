% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:39:08
% EndTime: 2019-05-08 04:39:32
% DurationCPUTime: 13.70s
% Computational Cost: add. (218325->366), mult. (434918->446), div. (0->0), fcn. (316428->10), ass. (0->141)
t776 = Ifges(6,1) + Ifges(7,1);
t772 = Ifges(6,4) + Ifges(7,4);
t771 = Ifges(6,5) + Ifges(7,5);
t775 = Ifges(6,2) + Ifges(7,2);
t770 = -Ifges(6,6) - Ifges(7,6);
t774 = -Ifges(6,3) - Ifges(7,3);
t748 = qJD(1) ^ 2;
t773 = pkin(2) * t748;
t742 = sin(qJ(1));
t747 = cos(qJ(1));
t731 = -g(1) * t747 - g(2) * t742;
t719 = -pkin(1) * t748 + qJDD(1) * pkin(7) + t731;
t741 = sin(qJ(2));
t769 = t741 * t719;
t746 = cos(qJ(2));
t762 = qJD(1) * qJD(2);
t725 = qJDD(1) * t741 + t746 * t762;
t683 = qJDD(2) * pkin(2) - t725 * pkin(8) - t769 + (pkin(8) * t762 + t741 * t773 - g(3)) * t746;
t706 = -g(3) * t741 + t746 * t719;
t726 = qJDD(1) * t746 - t741 * t762;
t764 = qJD(1) * t741;
t729 = qJD(2) * pkin(2) - pkin(8) * t764;
t737 = t746 ^ 2;
t684 = pkin(8) * t726 - qJD(2) * t729 - t737 * t773 + t706;
t740 = sin(qJ(3));
t745 = cos(qJ(3));
t661 = t740 * t683 + t745 * t684;
t717 = (t740 * t746 + t741 * t745) * qJD(1);
t690 = -t717 * qJD(3) - t740 * t725 + t726 * t745;
t763 = qJD(1) * t746;
t716 = -t740 * t764 + t745 * t763;
t700 = -mrSges(4,1) * t716 + mrSges(4,2) * t717;
t736 = qJD(2) + qJD(3);
t708 = mrSges(4,1) * t736 - mrSges(4,3) * t717;
t735 = qJDD(2) + qJDD(3);
t691 = qJD(3) * t716 + t725 * t745 + t726 * t740;
t730 = t742 * g(1) - t747 * g(2);
t753 = -qJDD(1) * pkin(1) - t730;
t692 = -t726 * pkin(2) + t729 * t764 + (-pkin(8) * t737 - pkin(7)) * t748 + t753;
t644 = (-t716 * t736 - t691) * pkin(9) + (t717 * t736 - t690) * pkin(3) + t692;
t701 = -pkin(3) * t716 - pkin(9) * t717;
t734 = t736 ^ 2;
t648 = -pkin(3) * t734 + pkin(9) * t735 + t701 * t716 + t661;
t739 = sin(qJ(4));
t744 = cos(qJ(4));
t633 = t744 * t644 - t739 * t648;
t703 = -t717 * t739 + t736 * t744;
t665 = qJD(4) * t703 + t691 * t744 + t735 * t739;
t689 = qJDD(4) - t690;
t704 = t717 * t744 + t736 * t739;
t712 = qJD(4) - t716;
t630 = (t703 * t712 - t665) * pkin(10) + (t703 * t704 + t689) * pkin(4) + t633;
t634 = t739 * t644 + t744 * t648;
t664 = -qJD(4) * t704 - t691 * t739 + t735 * t744;
t695 = pkin(4) * t712 - pkin(10) * t704;
t702 = t703 ^ 2;
t632 = -pkin(4) * t702 + pkin(10) * t664 - t695 * t712 + t634;
t738 = sin(qJ(5));
t743 = cos(qJ(5));
t624 = t743 * t630 - t738 * t632;
t677 = t703 * t743 - t704 * t738;
t641 = qJD(5) * t677 + t664 * t738 + t665 * t743;
t678 = t703 * t738 + t704 * t743;
t658 = -mrSges(7,1) * t677 + mrSges(7,2) * t678;
t659 = -mrSges(6,1) * t677 + mrSges(6,2) * t678;
t710 = qJD(5) + t712;
t667 = -mrSges(6,2) * t710 + mrSges(6,3) * t677;
t686 = qJDD(5) + t689;
t621 = -0.2e1 * qJD(6) * t678 + (t677 * t710 - t641) * qJ(6) + (t677 * t678 + t686) * pkin(5) + t624;
t666 = -mrSges(7,2) * t710 + mrSges(7,3) * t677;
t761 = m(7) * t621 + t686 * mrSges(7,1) + t710 * t666;
t613 = m(6) * t624 + t686 * mrSges(6,1) + t710 * t667 + (-t658 - t659) * t678 + (-mrSges(6,3) - mrSges(7,3)) * t641 + t761;
t625 = t738 * t630 + t743 * t632;
t640 = -qJD(5) * t678 + t664 * t743 - t665 * t738;
t669 = mrSges(7,1) * t710 - mrSges(7,3) * t678;
t670 = mrSges(6,1) * t710 - mrSges(6,3) * t678;
t668 = pkin(5) * t710 - qJ(6) * t678;
t676 = t677 ^ 2;
t623 = -pkin(5) * t676 + qJ(6) * t640 + 0.2e1 * qJD(6) * t677 - t668 * t710 + t625;
t760 = m(7) * t623 + t640 * mrSges(7,3) + t677 * t658;
t616 = m(6) * t625 + t640 * mrSges(6,3) + t677 * t659 + (-t669 - t670) * t710 + (-mrSges(6,2) - mrSges(7,2)) * t686 + t760;
t611 = t743 * t613 + t738 * t616;
t682 = -mrSges(5,1) * t703 + mrSges(5,2) * t704;
t693 = -mrSges(5,2) * t712 + mrSges(5,3) * t703;
t608 = m(5) * t633 + mrSges(5,1) * t689 - mrSges(5,3) * t665 - t682 * t704 + t693 * t712 + t611;
t694 = mrSges(5,1) * t712 - mrSges(5,3) * t704;
t755 = -t613 * t738 + t743 * t616;
t609 = m(5) * t634 - mrSges(5,2) * t689 + mrSges(5,3) * t664 + t682 * t703 - t694 * t712 + t755;
t756 = -t608 * t739 + t744 * t609;
t602 = m(4) * t661 - mrSges(4,2) * t735 + mrSges(4,3) * t690 + t700 * t716 - t708 * t736 + t756;
t660 = t683 * t745 - t740 * t684;
t707 = -mrSges(4,2) * t736 + mrSges(4,3) * t716;
t647 = -pkin(3) * t735 - pkin(9) * t734 + t717 * t701 - t660;
t635 = -pkin(4) * t664 - pkin(10) * t702 + t704 * t695 + t647;
t627 = -pkin(5) * t640 - qJ(6) * t676 + t668 * t678 + qJDD(6) + t635;
t754 = m(7) * t627 - t640 * mrSges(7,1) + t641 * mrSges(7,2) - t677 * t666 + t678 * t669;
t751 = m(6) * t635 - t640 * mrSges(6,1) + mrSges(6,2) * t641 - t677 * t667 + t670 * t678 + t754;
t749 = -m(5) * t647 + t664 * mrSges(5,1) - mrSges(5,2) * t665 + t703 * t693 - t694 * t704 - t751;
t618 = m(4) * t660 + mrSges(4,1) * t735 - mrSges(4,3) * t691 - t700 * t717 + t707 * t736 + t749;
t597 = t740 * t602 + t745 * t618;
t705 = -t746 * g(3) - t769;
t724 = (-mrSges(3,1) * t746 + mrSges(3,2) * t741) * qJD(1);
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t763;
t595 = m(3) * t705 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t725 + qJD(2) * t728 - t724 * t764 + t597;
t727 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t764;
t757 = t745 * t602 - t618 * t740;
t596 = m(3) * t706 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t726 - qJD(2) * t727 + t724 * t763 + t757;
t758 = -t595 * t741 + t746 * t596;
t587 = m(2) * t731 - mrSges(2,1) * t748 - qJDD(1) * mrSges(2,2) + t758;
t718 = -t748 * pkin(7) + t753;
t603 = t744 * t608 + t739 * t609;
t752 = m(4) * t692 - t690 * mrSges(4,1) + mrSges(4,2) * t691 - t716 * t707 + t708 * t717 + t603;
t750 = -m(3) * t718 + t726 * mrSges(3,1) - mrSges(3,2) * t725 - t727 * t764 + t728 * t763 - t752;
t599 = m(2) * t730 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t748 + t750;
t768 = t742 * t587 + t747 * t599;
t588 = t746 * t595 + t741 * t596;
t767 = t770 * t677 - t771 * t678 + t774 * t710;
t766 = -t775 * t677 - t772 * t678 + t770 * t710;
t765 = t772 * t677 + t776 * t678 + t771 * t710;
t759 = t747 * t587 - t599 * t742;
t715 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t741 + Ifges(3,4) * t746) * qJD(1);
t714 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t741 + Ifges(3,2) * t746) * qJD(1);
t713 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t741 + Ifges(3,6) * t746) * qJD(1);
t698 = Ifges(4,1) * t717 + Ifges(4,4) * t716 + Ifges(4,5) * t736;
t697 = Ifges(4,4) * t717 + Ifges(4,2) * t716 + Ifges(4,6) * t736;
t696 = Ifges(4,5) * t717 + Ifges(4,6) * t716 + Ifges(4,3) * t736;
t673 = Ifges(5,1) * t704 + Ifges(5,4) * t703 + Ifges(5,5) * t712;
t672 = Ifges(5,4) * t704 + Ifges(5,2) * t703 + Ifges(5,6) * t712;
t671 = Ifges(5,5) * t704 + Ifges(5,6) * t703 + Ifges(5,3) * t712;
t619 = -t641 * mrSges(7,3) - t678 * t658 + t761;
t610 = mrSges(6,2) * t635 + mrSges(7,2) * t627 - mrSges(6,3) * t624 - mrSges(7,3) * t621 - qJ(6) * t619 + t772 * t640 + t776 * t641 - t767 * t677 + t771 * t686 + t766 * t710;
t604 = -mrSges(6,1) * t635 + mrSges(6,3) * t625 - mrSges(7,1) * t627 + mrSges(7,3) * t623 - pkin(5) * t754 + qJ(6) * t760 + (-qJ(6) * t669 + t765) * t710 + (-mrSges(7,2) * qJ(6) - t770) * t686 + t767 * t678 + t772 * t641 + t775 * t640;
t591 = mrSges(5,2) * t647 - mrSges(5,3) * t633 + Ifges(5,1) * t665 + Ifges(5,4) * t664 + Ifges(5,5) * t689 - pkin(10) * t611 - t604 * t738 + t610 * t743 + t671 * t703 - t672 * t712;
t590 = -mrSges(5,1) * t647 + mrSges(5,3) * t634 + Ifges(5,4) * t665 + Ifges(5,2) * t664 + Ifges(5,6) * t689 - pkin(4) * t751 + pkin(10) * t755 + t743 * t604 + t738 * t610 - t704 * t671 + t712 * t673;
t589 = t774 * t686 + t736 * t698 + Ifges(4,6) * t735 - t717 * t696 + t703 * t673 - t704 * t672 - mrSges(4,1) * t692 - Ifges(5,3) * t689 + Ifges(4,2) * t690 + Ifges(4,4) * t691 + mrSges(4,3) * t661 - Ifges(5,6) * t664 - Ifges(5,5) * t665 + t770 * t640 - t771 * t641 - mrSges(5,1) * t633 + mrSges(5,2) * t634 + mrSges(6,2) * t625 - mrSges(6,1) * t624 + mrSges(7,2) * t623 - mrSges(7,1) * t621 - pkin(5) * t619 - pkin(4) * t611 - pkin(3) * t603 + t765 * t677 + t766 * t678;
t584 = mrSges(4,2) * t692 - mrSges(4,3) * t660 + Ifges(4,1) * t691 + Ifges(4,4) * t690 + Ifges(4,5) * t735 - pkin(9) * t603 - t590 * t739 + t591 * t744 + t696 * t716 - t697 * t736;
t583 = mrSges(3,2) * t718 - mrSges(3,3) * t705 + Ifges(3,1) * t725 + Ifges(3,4) * t726 + Ifges(3,5) * qJDD(2) - pkin(8) * t597 - qJD(2) * t714 + t584 * t745 - t589 * t740 + t713 * t763;
t582 = -pkin(1) * t588 + mrSges(2,3) * t731 - pkin(2) * t597 - Ifges(3,5) * t725 - Ifges(3,6) * t726 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t705 + mrSges(3,2) * t706 - pkin(3) * t749 - pkin(9) * t756 + mrSges(4,2) * t661 - t739 * t591 - t744 * t590 - Ifges(4,5) * t691 - Ifges(4,6) * t690 - Ifges(4,3) * t735 - mrSges(4,1) * t660 + mrSges(2,1) * g(3) + t748 * Ifges(2,5) - t717 * t697 + t716 * t698 + Ifges(2,6) * qJDD(1) + (-t714 * t741 + t715 * t746) * qJD(1);
t581 = -mrSges(3,1) * t718 + mrSges(3,3) * t706 + Ifges(3,4) * t725 + Ifges(3,2) * t726 + Ifges(3,6) * qJDD(2) - pkin(2) * t752 + pkin(8) * t757 + qJD(2) * t715 + t740 * t584 + t745 * t589 - t713 * t764;
t580 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t748 - pkin(7) * t588 - t581 * t741 + t583 * t746;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t588; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t768 + t747 * t580 - t742 * t582; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t759 + t742 * t580 + t747 * t582; -mrSges(1,1) * g(2) + mrSges(2,1) * t730 + mrSges(1,2) * g(1) - mrSges(2,2) * t731 + Ifges(2,3) * qJDD(1) + pkin(1) * t750 + pkin(7) * t758 + t746 * t581 + t741 * t583;];
tauB  = t1;
