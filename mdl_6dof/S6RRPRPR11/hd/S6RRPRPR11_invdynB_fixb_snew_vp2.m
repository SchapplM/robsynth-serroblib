% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:00:30
% EndTime: 2019-05-06 16:00:45
% DurationCPUTime: 10.23s
% Computational Cost: add. (154397->368), mult. (327693->447), div. (0->0), fcn. (208586->10), ass. (0->142)
t780 = -2 * qJD(3);
t779 = Ifges(3,1) + Ifges(4,2);
t775 = Ifges(3,4) + Ifges(4,6);
t774 = Ifges(3,5) - Ifges(4,4);
t778 = Ifges(3,2) + Ifges(4,3);
t773 = Ifges(3,6) - Ifges(4,5);
t777 = (Ifges(3,3) + Ifges(4,1));
t742 = sin(qJ(1));
t746 = cos(qJ(1));
t724 = -g(1) * t746 - g(2) * t742;
t748 = qJD(1) ^ 2;
t700 = -pkin(1) * t748 + qJDD(1) * pkin(7) + t724;
t741 = sin(qJ(2));
t745 = cos(qJ(2));
t687 = -g(3) * t741 + t745 * t700;
t711 = (-pkin(2) * t745 - qJ(3) * t741) * qJD(1);
t747 = qJD(2) ^ 2;
t767 = qJD(1) * t745;
t664 = pkin(2) * t747 - qJDD(2) * qJ(3) + (qJD(2) * t780) - t711 * t767 - t687;
t776 = t748 * pkin(7);
t686 = -t745 * g(3) - t741 * t700;
t712 = (mrSges(4,2) * t745 - mrSges(4,3) * t741) * qJD(1);
t713 = (-mrSges(3,1) * t745 + mrSges(3,2) * t741) * qJD(1);
t766 = qJD(1) * qJD(2);
t764 = t745 * t766;
t714 = qJDD(1) * t741 + t764;
t719 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t767;
t720 = -mrSges(4,1) * t767 - qJD(2) * mrSges(4,3);
t730 = t741 * qJD(1);
t763 = t741 * t766;
t715 = qJDD(1) * t745 - t763;
t722 = pkin(3) * t730 - (qJD(2) * pkin(8));
t736 = t745 ^ 2;
t723 = t742 * g(1) - t746 * g(2);
t758 = -qJDD(1) * pkin(1) - t723;
t753 = pkin(2) * t763 + t730 * t780 + (-t714 - t764) * qJ(3) + t758;
t645 = -t722 * t730 + (-pkin(3) * t736 - pkin(7)) * t748 + (-pkin(2) - pkin(8)) * t715 + t753;
t665 = -qJDD(2) * pkin(2) - t747 * qJ(3) + t711 * t730 + qJDD(3) - t686;
t653 = (-t741 * t745 * t748 - qJDD(2)) * pkin(8) + (t714 - t764) * pkin(3) + t665;
t740 = sin(qJ(4));
t744 = cos(qJ(4));
t635 = -t740 * t645 + t744 * t653;
t709 = -qJD(2) * t740 - t744 * t767;
t678 = qJD(4) * t709 + qJDD(2) * t744 - t715 * t740;
t708 = qJDD(4) + t714;
t710 = qJD(2) * t744 - t740 * t767;
t727 = t730 + qJD(4);
t629 = (t709 * t727 - t678) * qJ(5) + (t709 * t710 + t708) * pkin(4) + t635;
t636 = t744 * t645 + t740 * t653;
t677 = -qJD(4) * t710 - qJDD(2) * t740 - t715 * t744;
t684 = pkin(4) * t727 - qJ(5) * t710;
t707 = t709 ^ 2;
t631 = -pkin(4) * t707 + qJ(5) * t677 - t684 * t727 + t636;
t737 = sin(pkin(10));
t738 = cos(pkin(10));
t681 = t709 * t737 + t710 * t738;
t623 = -0.2e1 * qJD(5) * t681 + t738 * t629 - t737 * t631;
t658 = t677 * t737 + t678 * t738;
t680 = t709 * t738 - t710 * t737;
t621 = (t680 * t727 - t658) * pkin(9) + (t680 * t681 + t708) * pkin(5) + t623;
t624 = 0.2e1 * qJD(5) * t680 + t737 * t629 + t738 * t631;
t657 = t677 * t738 - t678 * t737;
t668 = pkin(5) * t727 - pkin(9) * t681;
t679 = t680 ^ 2;
t622 = -pkin(5) * t679 + pkin(9) * t657 - t668 * t727 + t624;
t739 = sin(qJ(6));
t743 = cos(qJ(6));
t619 = t621 * t743 - t622 * t739;
t660 = t680 * t743 - t681 * t739;
t634 = qJD(6) * t660 + t657 * t739 + t658 * t743;
t661 = t680 * t739 + t681 * t743;
t643 = -mrSges(7,1) * t660 + mrSges(7,2) * t661;
t725 = qJD(6) + t727;
t646 = -mrSges(7,2) * t725 + mrSges(7,3) * t660;
t701 = qJDD(6) + t708;
t614 = m(7) * t619 + mrSges(7,1) * t701 - mrSges(7,3) * t634 - t643 * t661 + t646 * t725;
t620 = t621 * t739 + t622 * t743;
t633 = -qJD(6) * t661 + t657 * t743 - t658 * t739;
t647 = mrSges(7,1) * t725 - mrSges(7,3) * t661;
t615 = m(7) * t620 - mrSges(7,2) * t701 + mrSges(7,3) * t633 + t643 * t660 - t647 * t725;
t608 = t743 * t614 + t739 * t615;
t662 = -mrSges(6,1) * t680 + mrSges(6,2) * t681;
t666 = -mrSges(6,2) * t727 + mrSges(6,3) * t680;
t606 = m(6) * t623 + mrSges(6,1) * t708 - mrSges(6,3) * t658 - t662 * t681 + t666 * t727 + t608;
t667 = mrSges(6,1) * t727 - mrSges(6,3) * t681;
t759 = -t614 * t739 + t743 * t615;
t607 = m(6) * t624 - mrSges(6,2) * t708 + mrSges(6,3) * t657 + t662 * t680 - t667 * t727 + t759;
t602 = t738 * t606 + t737 * t607;
t682 = -mrSges(5,1) * t709 + mrSges(5,2) * t710;
t683 = -mrSges(5,2) * t727 + mrSges(5,3) * t709;
t600 = m(5) * t635 + mrSges(5,1) * t708 - mrSges(5,3) * t678 - t682 * t710 + t683 * t727 + t602;
t685 = mrSges(5,1) * t727 - mrSges(5,3) * t710;
t760 = -t606 * t737 + t738 * t607;
t601 = m(5) * t636 - mrSges(5,2) * t708 + mrSges(5,3) * t677 + t682 * t709 - t685 * t727 + t760;
t595 = t744 * t600 + t740 * t601;
t754 = -m(4) * t665 - t714 * mrSges(4,1) - t595;
t593 = m(3) * t686 - t714 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t719 - t720) * qJD(2) + (-t712 - t713) * t730 + t754;
t718 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t730;
t721 = mrSges(4,1) * t730 + qJD(2) * mrSges(4,2);
t652 = -pkin(8) * t736 * t748 + pkin(3) * t715 + qJD(2) * t722 - t664;
t638 = -pkin(4) * t677 - qJ(5) * t707 + t710 * t684 + qJDD(5) + t652;
t626 = -pkin(5) * t657 - pkin(9) * t679 + t668 * t681 + t638;
t755 = m(7) * t626 - t633 * mrSges(7,1) + t634 * mrSges(7,2) - t660 * t646 + t661 * t647;
t752 = m(6) * t638 - t657 * mrSges(6,1) + t658 * mrSges(6,2) - t680 * t666 + t681 * t667 + t755;
t750 = -m(5) * t652 + t677 * mrSges(5,1) - t678 * mrSges(5,2) + t709 * t683 - t710 * t685 - t752;
t749 = -m(4) * t664 + qJDD(2) * mrSges(4,3) + qJD(2) * t721 + t712 * t767 - t750;
t618 = t713 * t767 + (mrSges(3,3) + mrSges(4,1)) * t715 - qJDD(2) * mrSges(3,2) + t749 - qJD(2) * t718 + m(3) * t687;
t761 = -t593 * t741 + t745 * t618;
t588 = m(2) * t724 - mrSges(2,1) * t748 - qJDD(1) * mrSges(2,2) + t761;
t699 = t758 - t776;
t663 = -t715 * pkin(2) + t753 - t776;
t771 = -t740 * t600 + t744 * t601;
t757 = -m(4) * t663 - t715 * mrSges(4,2) + t721 * t730 - t771;
t751 = -m(3) * t699 + t719 * t767 + t715 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t714 + (-t718 * t741 - t720 * t745) * qJD(1) + t757;
t591 = m(2) * t723 + qJDD(1) * mrSges(2,1) - t748 * mrSges(2,2) + t751;
t772 = t742 * t588 + t746 * t591;
t589 = t745 * t593 + t741 * t618;
t770 = (t777 * qJD(2)) + (t774 * t741 + t773 * t745) * qJD(1);
t769 = t774 * qJD(2) + (t779 * t741 + t775 * t745) * qJD(1);
t768 = -t773 * qJD(2) + (-t775 * t741 - t778 * t745) * qJD(1);
t762 = t746 * t588 - t591 * t742;
t671 = Ifges(5,1) * t710 + Ifges(5,4) * t709 + Ifges(5,5) * t727;
t670 = Ifges(5,4) * t710 + Ifges(5,2) * t709 + Ifges(5,6) * t727;
t669 = Ifges(5,5) * t710 + Ifges(5,6) * t709 + Ifges(5,3) * t727;
t656 = Ifges(6,1) * t681 + Ifges(6,4) * t680 + Ifges(6,5) * t727;
t655 = Ifges(6,4) * t681 + Ifges(6,2) * t680 + Ifges(6,6) * t727;
t654 = Ifges(6,5) * t681 + Ifges(6,6) * t680 + Ifges(6,3) * t727;
t641 = Ifges(7,1) * t661 + Ifges(7,4) * t660 + Ifges(7,5) * t725;
t640 = Ifges(7,4) * t661 + Ifges(7,2) * t660 + Ifges(7,6) * t725;
t639 = Ifges(7,5) * t661 + Ifges(7,6) * t660 + Ifges(7,3) * t725;
t610 = mrSges(7,2) * t626 - mrSges(7,3) * t619 + Ifges(7,1) * t634 + Ifges(7,4) * t633 + Ifges(7,5) * t701 + t639 * t660 - t640 * t725;
t609 = -mrSges(7,1) * t626 + mrSges(7,3) * t620 + Ifges(7,4) * t634 + Ifges(7,2) * t633 + Ifges(7,6) * t701 - t639 * t661 + t641 * t725;
t597 = mrSges(6,2) * t638 - mrSges(6,3) * t623 + Ifges(6,1) * t658 + Ifges(6,4) * t657 + Ifges(6,5) * t708 - pkin(9) * t608 - t609 * t739 + t610 * t743 + t654 * t680 - t655 * t727;
t596 = -mrSges(6,1) * t638 + mrSges(6,3) * t624 + Ifges(6,4) * t658 + Ifges(6,2) * t657 + Ifges(6,6) * t708 - pkin(5) * t755 + pkin(9) * t759 + t743 * t609 + t739 * t610 - t681 * t654 + t727 * t656;
t594 = -t714 * mrSges(4,3) + t720 * t767 - t757;
t585 = mrSges(5,2) * t652 - mrSges(5,3) * t635 + Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t708 - qJ(5) * t602 - t596 * t737 + t597 * t738 + t669 * t709 - t670 * t727;
t584 = -mrSges(5,1) * t652 + mrSges(5,3) * t636 + Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t708 - pkin(4) * t752 + qJ(5) * t760 + t738 * t596 + t737 * t597 - t710 * t669 + t727 * t671;
t583 = (Ifges(6,3) + Ifges(5,3)) * t708 + pkin(3) * t595 + t779 * t714 - qJ(3) * t594 - t709 * t671 + t710 * t670 + mrSges(3,2) * t699 + Ifges(7,3) * t701 + t681 * t655 - mrSges(3,3) * t686 + Ifges(5,6) * t677 + Ifges(5,5) * t678 - t680 * t656 + mrSges(4,1) * t665 + Ifges(6,5) * t658 - t660 * t641 + t661 * t640 - mrSges(4,3) * t663 + Ifges(6,6) * t657 + Ifges(7,5) * t634 + mrSges(5,1) * t635 - mrSges(5,2) * t636 + Ifges(7,6) * t633 + mrSges(6,1) * t623 - mrSges(6,2) * t624 - mrSges(7,2) * t620 + mrSges(7,1) * t619 + pkin(5) * t608 + pkin(4) * t602 + t768 * qJD(2) + t770 * t767 + t774 * qJDD(2) + t775 * t715;
t582 = -mrSges(3,1) * t699 - mrSges(4,1) * t664 + mrSges(4,2) * t663 + mrSges(3,3) * t687 - pkin(2) * t594 - pkin(3) * t750 - pkin(8) * t771 + t769 * qJD(2) + t773 * qJDD(2) - t744 * t584 - t740 * t585 + t775 * t714 + t778 * t715 - t770 * t730;
t581 = -pkin(1) * t589 + mrSges(2,3) * t724 - pkin(2) * (-qJD(2) * t720 + t754) - qJ(3) * t749 - mrSges(3,1) * t686 + mrSges(3,2) * t687 - t744 * t585 + t740 * t584 + pkin(8) * t595 - mrSges(4,2) * t665 + mrSges(4,3) * t664 + mrSges(2,1) * g(3) + t748 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-mrSges(4,1) * qJ(3) - t773) * t715 - t774 * t714 + (mrSges(4,2) * pkin(2) - t777) * qJDD(2) + (t769 * t745 + (pkin(2) * t712 + t768) * t741) * qJD(1);
t580 = -mrSges(2,2) * g(3) - mrSges(2,3) * t723 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t748 - pkin(7) * t589 - t582 * t741 + t583 * t745;
t1 = [-m(1) * g(1) + t762; -m(1) * g(2) + t772; (-m(1) - m(2)) * g(3) + t589; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t772 + t746 * t580 - t742 * t581; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t762 + t742 * t580 + t746 * t581; -mrSges(1,1) * g(2) + mrSges(2,1) * t723 + mrSges(1,2) * g(1) - mrSges(2,2) * t724 + Ifges(2,3) * qJDD(1) + pkin(1) * t751 + pkin(7) * t761 + t745 * t582 + t741 * t583;];
tauB  = t1;
