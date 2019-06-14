% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-05-05 17:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:58:15
% EndTime: 2019-05-05 16:58:23
% DurationCPUTime: 7.96s
% Computational Cost: add. (106218->344), mult. (263816->417), div. (0->0), fcn. (190440->10), ass. (0->142)
t779 = -2 * qJD(4);
t778 = Ifges(4,1) + Ifges(5,2);
t773 = Ifges(4,5) - Ifges(5,4);
t777 = Ifges(4,2) + Ifges(5,3);
t772 = Ifges(4,6) - Ifges(5,5);
t771 = -Ifges(5,6) - Ifges(4,4);
t776 = (-Ifges(4,3) - Ifges(5,1));
t736 = qJD(1) ^ 2;
t732 = sin(qJ(1));
t734 = cos(qJ(1));
t712 = -g(1) * t734 - g(2) * t732;
t708 = -pkin(1) * t736 + qJDD(1) * qJ(2) + t712;
t727 = sin(pkin(9));
t729 = cos(pkin(9));
t760 = qJD(1) * qJD(2);
t755 = -g(3) * t729 - 0.2e1 * t727 * t760;
t774 = pkin(2) * t729;
t658 = (-pkin(7) * qJDD(1) + t736 * t774 - t708) * t727 + t755;
t687 = -g(3) * t727 + (t708 + 0.2e1 * t760) * t729;
t757 = qJDD(1) * t729;
t723 = t729 ^ 2;
t769 = t723 * t736;
t665 = -pkin(2) * t769 + pkin(7) * t757 + t687;
t731 = sin(qJ(3));
t775 = cos(qJ(3));
t642 = t731 * t658 + t775 * t665;
t756 = t729 * t775;
t763 = qJD(1) * t727;
t706 = -qJD(1) * t756 + t731 * t763;
t745 = t727 * t775 + t729 * t731;
t707 = t745 * qJD(1);
t675 = pkin(3) * t706 - qJ(4) * t707;
t735 = qJD(3) ^ 2;
t631 = pkin(3) * t735 - qJDD(3) * qJ(4) + (qJD(3) * t779) + t706 * t675 - t642;
t770 = mrSges(3,2) * t727;
t641 = t658 * t775 - t731 * t665;
t676 = mrSges(4,1) * t706 + mrSges(4,2) * t707;
t761 = t706 * qJD(3);
t685 = qJDD(1) * t745 - t761;
t695 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t706;
t698 = mrSges(5,1) * t706 - (qJD(3) * mrSges(5,3));
t758 = qJDD(1) * t727;
t762 = qJD(3) * t707;
t684 = -qJDD(1) * t756 + t731 * t758 + t762;
t697 = pkin(4) * t707 - (qJD(3) * qJ(5));
t705 = t706 ^ 2;
t722 = t727 ^ 2;
t711 = g(1) * t732 - t734 * g(2);
t750 = qJDD(2) - t711;
t683 = (-pkin(1) - t774) * qJDD(1) + (-qJ(2) + (-t722 - t723) * pkin(7)) * t736 + t750;
t737 = pkin(3) * t762 + t707 * t779 + (-t685 + t761) * qJ(4) + t683;
t623 = -pkin(4) * t705 - t697 * t707 + (pkin(3) + qJ(5)) * t684 + t737;
t633 = -qJDD(3) * pkin(3) - t735 * qJ(4) + t707 * t675 + qJDD(4) - t641;
t626 = (t706 * t707 - qJDD(3)) * qJ(5) + (t685 + t761) * pkin(4) + t633;
t726 = sin(pkin(10));
t728 = cos(pkin(10));
t692 = qJD(3) * t728 + t706 * t726;
t618 = -0.2e1 * qJD(5) * t692 - t623 * t726 + t728 * t626;
t664 = qJDD(3) * t728 + t684 * t726;
t691 = -qJD(3) * t726 + t706 * t728;
t616 = (t691 * t707 - t664) * pkin(8) + (t691 * t692 + t685) * pkin(5) + t618;
t619 = 0.2e1 * qJD(5) * t691 + t728 * t623 + t726 * t626;
t662 = pkin(5) * t707 - pkin(8) * t692;
t663 = -qJDD(3) * t726 + t684 * t728;
t690 = t691 ^ 2;
t617 = -pkin(5) * t690 + pkin(8) * t663 - t662 * t707 + t619;
t730 = sin(qJ(6));
t733 = cos(qJ(6));
t614 = t616 * t733 - t617 * t730;
t649 = t691 * t733 - t692 * t730;
t635 = qJD(6) * t649 + t663 * t730 + t664 * t733;
t650 = t691 * t730 + t692 * t733;
t640 = -mrSges(7,1) * t649 + mrSges(7,2) * t650;
t703 = qJD(6) + t707;
t643 = -mrSges(7,2) * t703 + mrSges(7,3) * t649;
t682 = qJDD(6) + t685;
t612 = m(7) * t614 + mrSges(7,1) * t682 - mrSges(7,3) * t635 - t640 * t650 + t643 * t703;
t615 = t616 * t730 + t617 * t733;
t634 = -qJD(6) * t650 + t663 * t733 - t664 * t730;
t644 = mrSges(7,1) * t703 - mrSges(7,3) * t650;
t613 = m(7) * t615 - mrSges(7,2) * t682 + mrSges(7,3) * t634 + t640 * t649 - t644 * t703;
t603 = t733 * t612 + t730 * t613;
t651 = -mrSges(6,1) * t691 + mrSges(6,2) * t692;
t660 = -mrSges(6,2) * t707 + mrSges(6,3) * t691;
t601 = m(6) * t618 + mrSges(6,1) * t685 - mrSges(6,3) * t664 - t651 * t692 + t660 * t707 + t603;
t661 = mrSges(6,1) * t707 - mrSges(6,3) * t692;
t751 = -t612 * t730 + t733 * t613;
t602 = m(6) * t619 - mrSges(6,2) * t685 + mrSges(6,3) * t663 + t651 * t691 - t661 * t707 + t751;
t598 = t601 * t728 + t602 * t726;
t677 = -mrSges(5,2) * t706 - mrSges(5,3) * t707;
t742 = -m(5) * t633 - t685 * mrSges(5,1) - t707 * t677 - t598;
t596 = m(4) * t641 - mrSges(4,3) * t685 - t676 * t707 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t695 - t698) * qJD(3) + t742;
t696 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t707;
t699 = mrSges(5,1) * t707 + qJD(3) * mrSges(5,2);
t628 = -pkin(4) * t684 - qJ(5) * t705 + qJD(3) * t697 + qJDD(5) - t631;
t621 = -pkin(5) * t663 - pkin(8) * t690 + t662 * t692 + t628;
t743 = m(7) * t621 - mrSges(7,1) * t634 + t635 * mrSges(7,2) - t643 * t649 + t650 * t644;
t741 = -m(6) * t628 + mrSges(6,1) * t663 - t664 * mrSges(6,2) + t660 * t691 - t692 * t661 - t743;
t739 = -m(5) * t631 + qJDD(3) * mrSges(5,3) + qJD(3) * t699 - t741;
t608 = -qJDD(3) * mrSges(4,2) + (-t676 - t677) * t706 + (-mrSges(4,3) - mrSges(5,1)) * t684 + m(4) * t642 - qJD(3) * t696 + t739;
t590 = t775 * t596 + t731 * t608;
t686 = -t708 * t727 + t755;
t746 = qJDD(1) * mrSges(3,3) + t736 * (-mrSges(3,1) * t729 + t770);
t588 = m(3) * t686 - t727 * t746 + t590;
t752 = -t596 * t731 + t775 * t608;
t589 = m(3) * t687 + t729 * t746 + t752;
t753 = -t588 * t727 + t729 * t589;
t583 = m(2) * t712 - mrSges(2,1) * t736 - qJDD(1) * mrSges(2,2) + t753;
t704 = -qJDD(1) * pkin(1) - qJ(2) * t736 + t750;
t630 = pkin(3) * t684 + t737;
t767 = -t726 * t601 + t728 * t602;
t597 = m(5) * t630 - t684 * mrSges(5,2) - t685 * mrSges(5,3) - t706 * t698 - t707 * t699 + t767;
t740 = m(4) * t683 + t684 * mrSges(4,1) + mrSges(4,2) * t685 + t706 * t695 + t696 * t707 + t597;
t738 = -m(3) * t704 + mrSges(3,1) * t757 - t740 + (t722 * t736 + t769) * mrSges(3,3);
t594 = (mrSges(2,1) - t770) * qJDD(1) + m(2) * t711 - mrSges(2,2) * t736 + t738;
t768 = t732 * t583 + t734 * t594;
t584 = t729 * t588 + t727 * t589;
t766 = (t776 * qJD(3)) + t772 * t706 - t773 * t707;
t765 = -t772 * qJD(3) + t777 * t706 + t771 * t707;
t764 = t773 * qJD(3) + t771 * t706 + t778 * t707;
t754 = t734 * t583 - t594 * t732;
t749 = Ifges(3,1) * t727 + Ifges(3,4) * t729;
t748 = Ifges(3,4) * t727 + Ifges(3,2) * t729;
t747 = Ifges(3,5) * t727 + Ifges(3,6) * t729;
t710 = t747 * qJD(1);
t647 = Ifges(6,1) * t692 + Ifges(6,4) * t691 + Ifges(6,5) * t707;
t646 = Ifges(6,4) * t692 + Ifges(6,2) * t691 + Ifges(6,6) * t707;
t645 = Ifges(6,5) * t692 + Ifges(6,6) * t691 + Ifges(6,3) * t707;
t638 = Ifges(7,1) * t650 + Ifges(7,4) * t649 + Ifges(7,5) * t703;
t637 = Ifges(7,4) * t650 + Ifges(7,2) * t649 + Ifges(7,6) * t703;
t636 = Ifges(7,5) * t650 + Ifges(7,6) * t649 + Ifges(7,3) * t703;
t605 = mrSges(7,2) * t621 - mrSges(7,3) * t614 + Ifges(7,1) * t635 + Ifges(7,4) * t634 + Ifges(7,5) * t682 + t636 * t649 - t637 * t703;
t604 = -mrSges(7,1) * t621 + mrSges(7,3) * t615 + Ifges(7,4) * t635 + Ifges(7,2) * t634 + Ifges(7,6) * t682 - t636 * t650 + t638 * t703;
t592 = mrSges(6,2) * t628 - mrSges(6,3) * t618 + Ifges(6,1) * t664 + Ifges(6,4) * t663 + Ifges(6,5) * t685 - pkin(8) * t603 - t604 * t730 + t605 * t733 + t645 * t691 - t646 * t707;
t591 = -mrSges(6,1) * t628 + mrSges(6,3) * t619 + Ifges(6,4) * t664 + Ifges(6,2) * t663 + Ifges(6,6) * t685 - pkin(5) * t743 + pkin(8) * t751 + t733 * t604 + t730 * t605 - t692 * t645 + t707 * t647;
t580 = t771 * t684 + t773 * qJDD(3) + t765 * qJD(3) + t766 * t706 + (Ifges(6,3) + t778) * t685 + t692 * t646 - t691 * t647 + Ifges(7,3) * t682 + mrSges(4,2) * t683 + Ifges(6,6) * t663 + Ifges(6,5) * t664 - t649 * t638 + t650 * t637 - mrSges(4,3) * t641 + mrSges(5,1) * t633 + Ifges(7,6) * t634 + Ifges(7,5) * t635 - mrSges(5,3) * t630 + mrSges(6,1) * t618 - mrSges(6,2) * t619 - mrSges(7,2) * t615 + mrSges(7,1) * t614 + pkin(5) * t603 + pkin(4) * t598 - qJ(4) * t597;
t579 = -mrSges(4,1) * t683 - mrSges(5,1) * t631 + mrSges(5,2) * t630 + mrSges(4,3) * t642 - pkin(3) * t597 - pkin(4) * t741 - qJ(5) * t767 + t764 * qJD(3) + t772 * qJDD(3) - t728 * t591 - t726 * t592 - t777 * t684 - t771 * t685 + t766 * t707;
t578 = -pkin(3) * (-qJD(3) * t698 + t742) + mrSges(2,1) * g(3) - t728 * t592 - qJ(4) * t739 + t726 * t591 + mrSges(2,3) * t712 - mrSges(3,1) * t686 + mrSges(3,2) * t687 - mrSges(4,1) * t641 + mrSges(4,2) * t642 + mrSges(5,3) * t631 - mrSges(5,2) * t633 + qJ(5) * t598 - pkin(2) * t590 - pkin(1) * t584 + t765 * t707 + (qJ(4) * t677 - t764) * t706 - t773 * t685 + (qJ(4) * mrSges(5,1) + t772) * t684 + (mrSges(5,2) * pkin(3) + t776) * qJDD(3) + (Ifges(2,6) - t747) * qJDD(1) + (-t727 * t748 + t729 * t749 + Ifges(2,5)) * t736;
t577 = t729 * qJD(1) * t710 + mrSges(3,2) * t704 - mrSges(3,3) * t686 - pkin(7) * t590 + qJDD(1) * t749 - t731 * t579 + t580 * t775;
t576 = -mrSges(3,1) * t704 + mrSges(3,3) * t687 - pkin(2) * t740 + pkin(7) * t752 + qJDD(1) * t748 + t579 * t775 + t731 * t580 - t710 * t763;
t575 = -mrSges(2,2) * g(3) - mrSges(2,3) * t711 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t736 - qJ(2) * t584 - t576 * t727 + t577 * t729;
t1 = [-m(1) * g(1) + t754; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t584; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t768 + t734 * t575 - t732 * t578; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t754 + t732 * t575 + t734 * t578; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t711 - mrSges(2,2) * t712 + t727 * t577 + t729 * t576 + pkin(1) * (-mrSges(3,2) * t758 + t738) + qJ(2) * t753;];
tauB  = t1;
