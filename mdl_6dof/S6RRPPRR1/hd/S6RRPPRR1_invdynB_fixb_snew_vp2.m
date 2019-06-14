% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:31:41
% EndTime: 2019-05-06 09:31:52
% DurationCPUTime: 8.32s
% Computational Cost: add. (115339->364), mult. (269240->444), div. (0->0), fcn. (187916->10), ass. (0->142)
t798 = -2 * qJD(3);
t797 = Ifges(4,1) + Ifges(5,1);
t792 = Ifges(4,4) - Ifges(5,5);
t791 = Ifges(4,5) + Ifges(5,4);
t796 = Ifges(4,2) + Ifges(5,3);
t795 = -Ifges(5,2) - Ifges(4,3);
t790 = Ifges(4,6) - Ifges(5,6);
t755 = sin(qJ(2));
t759 = cos(qJ(2));
t776 = qJD(1) * qJD(2);
t731 = qJDD(1) * t755 + t759 * t776;
t762 = qJD(1) ^ 2;
t756 = sin(qJ(1));
t760 = cos(qJ(1));
t737 = -g(1) * t760 - g(2) * t756;
t725 = -pkin(1) * t762 + qJDD(1) * pkin(7) + t737;
t786 = t755 * t725;
t671 = qJDD(2) * pkin(2) - t731 * qJ(3) - t786 + (pkin(2) * t755 * t762 + qJ(3) * t776 - g(3)) * t759;
t706 = -g(3) * t755 + t759 * t725;
t732 = qJDD(1) * t759 - t755 * t776;
t780 = qJD(1) * t755;
t733 = qJD(2) * pkin(2) - qJ(3) * t780;
t787 = t759 ^ 2 * t762;
t672 = -pkin(2) * t787 + qJ(3) * t732 - qJD(2) * t733 + t706;
t751 = sin(pkin(10));
t788 = cos(pkin(10));
t719 = (t751 * t759 + t788 * t755) * qJD(1);
t651 = t788 * t671 - t751 * t672 + t719 * t798;
t794 = 2 * qJD(4);
t793 = -mrSges(4,3) - mrSges(5,2);
t789 = pkin(3) * qJD(2);
t779 = qJD(1) * t759;
t718 = t751 * t780 - t788 * t779;
t652 = t751 * t671 + t788 * t672 + t718 * t798;
t701 = t731 * t751 - t788 * t732;
t708 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t719;
t693 = pkin(3) * t718 - qJ(4) * t719;
t761 = qJD(2) ^ 2;
t640 = -pkin(3) * t761 + qJDD(2) * qJ(4) + qJD(2) * t794 - t718 * t693 + t652;
t709 = -qJD(2) * mrSges(5,1) + mrSges(5,2) * t719;
t641 = -qJDD(2) * pkin(3) - t761 * qJ(4) + t719 * t693 + qJDD(4) - t651;
t702 = t788 * t731 + t751 * t732;
t778 = qJD(2) * t718;
t634 = (-t702 - t778) * pkin(8) + (t718 * t719 - qJDD(2)) * pkin(4) + t641;
t711 = -qJD(2) * pkin(4) - pkin(8) * t719;
t717 = t718 ^ 2;
t636 = -pkin(4) * t717 + pkin(8) * t701 + qJD(2) * t711 + t640;
t754 = sin(qJ(5));
t758 = cos(qJ(5));
t632 = t754 * t634 + t758 * t636;
t688 = t718 * t754 + t719 * t758;
t656 = -qJD(5) * t688 + t701 * t758 - t702 * t754;
t687 = t718 * t758 - t719 * t754;
t666 = -mrSges(6,1) * t687 + mrSges(6,2) * t688;
t745 = -qJD(2) + qJD(5);
t679 = mrSges(6,1) * t745 - mrSges(6,3) * t688;
t744 = -qJDD(2) + qJDD(5);
t667 = -pkin(5) * t687 - pkin(9) * t688;
t743 = t745 ^ 2;
t629 = -pkin(5) * t743 + pkin(9) * t744 + t667 * t687 + t632;
t736 = t756 * g(1) - t760 * g(2);
t724 = -qJDD(1) * pkin(1) - t762 * pkin(7) - t736;
t675 = -t732 * pkin(2) - qJ(3) * t787 + t733 * t780 + qJDD(3) + t724;
t765 = t701 * pkin(3) + t675 + (-t702 + t778) * qJ(4);
t638 = -t701 * pkin(4) - t717 * pkin(8) - t765 + (t711 - t789 + t794) * t719;
t657 = qJD(5) * t687 + t701 * t754 + t702 * t758;
t630 = (-t687 * t745 - t657) * pkin(9) + t638 + (t688 * t745 - t656) * pkin(5);
t753 = sin(qJ(6));
t757 = cos(qJ(6));
t626 = -t629 * t753 + t630 * t757;
t676 = -t688 * t753 + t745 * t757;
t644 = qJD(6) * t676 + t657 * t757 + t744 * t753;
t655 = qJDD(6) - t656;
t677 = t688 * t757 + t745 * t753;
t658 = -mrSges(7,1) * t676 + mrSges(7,2) * t677;
t680 = qJD(6) - t687;
t659 = -mrSges(7,2) * t680 + mrSges(7,3) * t676;
t624 = m(7) * t626 + mrSges(7,1) * t655 - mrSges(7,3) * t644 - t658 * t677 + t659 * t680;
t627 = t629 * t757 + t630 * t753;
t643 = -qJD(6) * t677 - t657 * t753 + t744 * t757;
t660 = mrSges(7,1) * t680 - mrSges(7,3) * t677;
t625 = m(7) * t627 - mrSges(7,2) * t655 + mrSges(7,3) * t643 + t658 * t676 - t660 * t680;
t771 = -t624 * t753 + t757 * t625;
t616 = m(6) * t632 - mrSges(6,2) * t744 + mrSges(6,3) * t656 + t666 * t687 - t679 * t745 + t771;
t631 = t634 * t758 - t636 * t754;
t678 = -mrSges(6,2) * t745 + mrSges(6,3) * t687;
t628 = -pkin(5) * t744 - pkin(9) * t743 + t667 * t688 - t631;
t767 = -m(7) * t628 + t643 * mrSges(7,1) - mrSges(7,2) * t644 + t676 * t659 - t660 * t677;
t620 = m(6) * t631 + mrSges(6,1) * t744 - mrSges(6,3) * t657 - t666 * t688 + t678 * t745 + t767;
t772 = t758 * t616 - t754 * t620;
t769 = m(5) * t640 + qJDD(2) * mrSges(5,3) + qJD(2) * t709 + t772;
t694 = mrSges(5,1) * t718 - mrSges(5,3) * t719;
t781 = -mrSges(4,1) * t718 - mrSges(4,2) * t719 - t694;
t609 = m(4) * t652 - qJDD(2) * mrSges(4,2) - qJD(2) * t708 + t793 * t701 + t781 * t718 + t769;
t707 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t718;
t611 = t754 * t616 + t758 * t620;
t710 = -mrSges(5,2) * t718 + qJD(2) * mrSges(5,3);
t766 = -m(5) * t641 + qJDD(2) * mrSges(5,1) + qJD(2) * t710 - t611;
t610 = m(4) * t651 + qJDD(2) * mrSges(4,1) + qJD(2) * t707 + t793 * t702 + t781 * t719 + t766;
t603 = t751 * t609 + t788 * t610;
t705 = -t759 * g(3) - t786;
t730 = (-mrSges(3,1) * t759 + mrSges(3,2) * t755) * qJD(1);
t735 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t779;
t601 = m(3) * t705 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t731 + qJD(2) * t735 - t730 * t780 + t603;
t734 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t780;
t773 = t788 * t609 - t610 * t751;
t602 = m(3) * t706 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t732 - qJD(2) * t734 + t730 * t779 + t773;
t774 = -t601 * t755 + t759 * t602;
t596 = m(2) * t737 - mrSges(2,1) * t762 - qJDD(1) * mrSges(2,2) + t774;
t646 = (-(2 * qJD(4)) + t789) * t719 + t765;
t617 = t757 * t624 + t753 * t625;
t768 = -m(6) * t638 + t656 * mrSges(6,1) - t657 * mrSges(6,2) + t687 * t678 - t688 * t679 - t617;
t614 = m(5) * t646 + t701 * mrSges(5,1) - t702 * mrSges(5,3) - t719 * t709 + t718 * t710 + t768;
t764 = m(4) * t675 + t701 * mrSges(4,1) + t702 * mrSges(4,2) + t718 * t707 + t719 * t708 + t614;
t763 = -m(3) * t724 + t732 * mrSges(3,1) - t731 * mrSges(3,2) - t734 * t780 + t735 * t779 - t764;
t613 = m(2) * t736 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t763;
t785 = t756 * t596 + t760 * t613;
t597 = t759 * t601 + t755 * t602;
t784 = -qJD(2) * t790 + t718 * t796 - t719 * t792;
t783 = qJD(2) * t795 + t718 * t790 - t719 * t791;
t782 = qJD(2) * t791 - t718 * t792 + t719 * t797;
t775 = t760 * t596 - t613 * t756;
t722 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t755 + Ifges(3,4) * t759) * qJD(1);
t721 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t755 + Ifges(3,2) * t759) * qJD(1);
t720 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t755 + Ifges(3,6) * t759) * qJD(1);
t663 = Ifges(6,1) * t688 + Ifges(6,4) * t687 + Ifges(6,5) * t745;
t662 = Ifges(6,4) * t688 + Ifges(6,2) * t687 + Ifges(6,6) * t745;
t661 = Ifges(6,5) * t688 + Ifges(6,6) * t687 + Ifges(6,3) * t745;
t649 = Ifges(7,1) * t677 + Ifges(7,4) * t676 + Ifges(7,5) * t680;
t648 = Ifges(7,4) * t677 + Ifges(7,2) * t676 + Ifges(7,6) * t680;
t647 = Ifges(7,5) * t677 + Ifges(7,6) * t676 + Ifges(7,3) * t680;
t619 = mrSges(7,2) * t628 - mrSges(7,3) * t626 + Ifges(7,1) * t644 + Ifges(7,4) * t643 + Ifges(7,5) * t655 + t647 * t676 - t648 * t680;
t618 = -mrSges(7,1) * t628 + mrSges(7,3) * t627 + Ifges(7,4) * t644 + Ifges(7,2) * t643 + Ifges(7,6) * t655 - t647 * t677 + t649 * t680;
t605 = -mrSges(6,1) * t638 - mrSges(7,1) * t626 + mrSges(7,2) * t627 + mrSges(6,3) * t632 + Ifges(6,4) * t657 - Ifges(7,5) * t644 + Ifges(6,2) * t656 + Ifges(6,6) * t744 - Ifges(7,6) * t643 - Ifges(7,3) * t655 - pkin(5) * t617 - t648 * t677 + t649 * t676 - t661 * t688 + t663 * t745;
t604 = mrSges(6,2) * t638 - mrSges(6,3) * t631 + Ifges(6,1) * t657 + Ifges(6,4) * t656 + Ifges(6,5) * t744 - pkin(9) * t617 - t618 * t753 + t619 * t757 + t661 * t687 - t662 * t745;
t593 = mrSges(4,2) * t675 + mrSges(5,2) * t641 - mrSges(4,3) * t651 - mrSges(5,3) * t646 - pkin(8) * t611 - qJ(4) * t614 + t784 * qJD(2) + t791 * qJDD(2) + t758 * t604 - t754 * t605 - t792 * t701 + t702 * t797 + t783 * t718;
t592 = -mrSges(4,1) * t675 - mrSges(5,1) * t646 + mrSges(5,2) * t640 + mrSges(4,3) * t652 - pkin(3) * t614 - pkin(4) * t768 - pkin(8) * t772 + t782 * qJD(2) + t790 * qJDD(2) - t754 * t604 - t758 * t605 - t701 * t796 + t792 * t702 + t783 * t719;
t591 = (mrSges(5,2) * qJ(4) + t790) * t701 + (mrSges(5,2) * pkin(3) - t791) * t702 - qJ(4) * t769 - pkin(3) * t766 + pkin(5) * t767 + (qJ(4) * t694 - t782) * t718 + (pkin(3) * t694 + t784) * t719 + pkin(9) * t771 + (-Ifges(3,3) + t795) * qJDD(2) + Ifges(2,6) * qJDD(1) - pkin(1) * t597 + mrSges(2,1) * g(3) + (-t721 * t755 + t722 * t759) * qJD(1) - pkin(2) * t603 + t762 * Ifges(2,5) + t757 * t618 + t753 * t619 + Ifges(6,3) * t744 - Ifges(3,6) * t732 + mrSges(2,3) * t737 - Ifges(3,5) * t731 - mrSges(3,1) * t705 + mrSges(3,2) * t706 - t687 * t663 + t688 * t662 - mrSges(4,1) * t651 + mrSges(4,2) * t652 + Ifges(6,6) * t656 + Ifges(6,5) * t657 - mrSges(5,3) * t640 + mrSges(5,1) * t641 + mrSges(6,1) * t631 - mrSges(6,2) * t632 + pkin(4) * t611;
t590 = mrSges(3,2) * t724 - mrSges(3,3) * t705 + Ifges(3,1) * t731 + Ifges(3,4) * t732 + Ifges(3,5) * qJDD(2) - qJ(3) * t603 - qJD(2) * t721 - t751 * t592 + t788 * t593 + t720 * t779;
t589 = -mrSges(3,1) * t724 + mrSges(3,3) * t706 + Ifges(3,4) * t731 + Ifges(3,2) * t732 + Ifges(3,6) * qJDD(2) - pkin(2) * t764 + qJ(3) * t773 + qJD(2) * t722 + t788 * t592 + t751 * t593 - t720 * t780;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t736 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t762 - pkin(7) * t597 - t589 * t755 + t590 * t759;
t1 = [-m(1) * g(1) + t775; -m(1) * g(2) + t785; (-m(1) - m(2)) * g(3) + t597; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t785 + t760 * t588 - t756 * t591; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t775 + t756 * t588 + t760 * t591; -mrSges(1,1) * g(2) + mrSges(2,1) * t736 + mrSges(1,2) * g(1) - mrSges(2,2) * t737 + Ifges(2,3) * qJDD(1) + pkin(1) * t763 + pkin(7) * t774 + t759 * t589 + t755 * t590;];
tauB  = t1;
