% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-05-06 08:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:11:05
% EndTime: 2019-05-06 08:11:14
% DurationCPUTime: 8.60s
% Computational Cost: add. (111023->368), mult. (263796->448), div. (0->0), fcn. (182496->10), ass. (0->145)
t796 = -2 * qJD(4);
t795 = Ifges(5,1) + Ifges(6,1);
t787 = Ifges(5,4) - Ifges(6,5);
t794 = -Ifges(5,5) - Ifges(6,4);
t793 = Ifges(5,2) + Ifges(6,3);
t792 = -Ifges(6,2) - Ifges(5,3);
t785 = Ifges(5,6) - Ifges(6,6);
t748 = sin(qJ(2));
t751 = cos(qJ(2));
t770 = qJD(1) * qJD(2);
t733 = t748 * qJDD(1) + t751 * t770;
t734 = t751 * qJDD(1) - t748 * t770;
t746 = sin(pkin(9));
t784 = cos(pkin(9));
t704 = t784 * t733 + t746 * t734;
t745 = sin(pkin(10));
t783 = cos(pkin(10));
t691 = t745 * qJDD(2) + t783 * t704;
t722 = (t746 * t751 + t784 * t748) * qJD(1);
t708 = -t783 * qJD(2) + t745 * t722;
t773 = qJD(1) * t751;
t774 = qJD(1) * t748;
t721 = t746 * t774 - t784 * t773;
t782 = t708 * t721;
t791 = (-t691 + t782) * qJ(5);
t749 = sin(qJ(1));
t752 = cos(qJ(1));
t739 = -t752 * g(1) - t749 * g(2);
t754 = qJD(1) ^ 2;
t728 = -t754 * pkin(1) + qJDD(1) * pkin(7) + t739;
t781 = t748 * t728;
t789 = pkin(2) * t754;
t676 = qJDD(2) * pkin(2) - t733 * qJ(3) - t781 + (qJ(3) * t770 + t748 * t789 - g(3)) * t751;
t711 = -t748 * g(3) + t751 * t728;
t735 = qJD(2) * pkin(2) - qJ(3) * t774;
t744 = t751 ^ 2;
t680 = t734 * qJ(3) - qJD(2) * t735 - t744 * t789 + t711;
t654 = -0.2e1 * qJD(3) * t721 + t746 * t676 + t784 * t680;
t696 = t721 * pkin(3) - t722 * qJ(4);
t753 = qJD(2) ^ 2;
t641 = -t753 * pkin(3) + qJDD(2) * qJ(4) - t721 * t696 + t654;
t738 = t749 * g(1) - t752 * g(2);
t761 = -qJDD(1) * pkin(1) - t738;
t682 = -t734 * pkin(2) + qJDD(3) + t735 * t774 + (-qJ(3) * t744 - pkin(7)) * t754 + t761;
t703 = t746 * t733 - t784 * t734;
t643 = (qJD(2) * t721 - t704) * qJ(4) + (qJD(2) * t722 + t703) * pkin(3) + t682;
t709 = t745 * qJD(2) + t783 * t722;
t636 = -t745 * t641 + t783 * t643 + t709 * t796;
t790 = 2 * qJD(5);
t788 = -mrSges(5,3) - mrSges(6,2);
t697 = t721 * mrSges(4,1) + t722 * mrSges(4,2);
t713 = qJD(2) * mrSges(4,1) - t722 * mrSges(4,3);
t637 = t783 * t641 + t745 * t643 + t708 * t796;
t685 = t721 * mrSges(5,1) - t709 * mrSges(5,3);
t690 = -t783 * qJDD(2) + t745 * t704;
t677 = t708 * pkin(4) - t709 * qJ(5);
t720 = t721 ^ 2;
t633 = -t720 * pkin(4) + t703 * qJ(5) - t708 * t677 + t721 * t790 + t637;
t686 = -t721 * mrSges(6,1) + t709 * mrSges(6,2);
t634 = -t703 * pkin(4) - t720 * qJ(5) + t709 * t677 + qJDD(5) - t636;
t628 = (-t691 - t782) * pkin(8) + (t708 * t709 - t703) * pkin(5) + t634;
t687 = -t721 * pkin(5) - t709 * pkin(8);
t707 = t708 ^ 2;
t629 = -t707 * pkin(5) + t690 * pkin(8) + t721 * t687 + t633;
t747 = sin(qJ(6));
t750 = cos(qJ(6));
t626 = t750 * t628 - t747 * t629;
t672 = t750 * t708 - t747 * t709;
t647 = t672 * qJD(6) + t747 * t690 + t750 * t691;
t673 = t747 * t708 + t750 * t709;
t655 = -t672 * mrSges(7,1) + t673 * mrSges(7,2);
t719 = qJD(6) - t721;
t658 = -t719 * mrSges(7,2) + t672 * mrSges(7,3);
t702 = qJDD(6) - t703;
t623 = m(7) * t626 + t702 * mrSges(7,1) - t647 * mrSges(7,3) - t673 * t655 + t719 * t658;
t627 = t747 * t628 + t750 * t629;
t646 = -t673 * qJD(6) + t750 * t690 - t747 * t691;
t659 = t719 * mrSges(7,1) - t673 * mrSges(7,3);
t624 = m(7) * t627 - t702 * mrSges(7,2) + t646 * mrSges(7,3) + t672 * t655 - t719 * t659;
t765 = -t747 * t623 + t750 * t624;
t760 = m(6) * t633 + t703 * mrSges(6,3) + t721 * t686 + t765;
t678 = t708 * mrSges(6,1) - t709 * mrSges(6,3);
t775 = -t708 * mrSges(5,1) - t709 * mrSges(5,2) - t678;
t614 = m(5) * t637 - t703 * mrSges(5,2) - t721 * t685 + t788 * t690 + t775 * t708 + t760;
t684 = -t721 * mrSges(5,2) - t708 * mrSges(5,3);
t617 = t750 * t623 + t747 * t624;
t683 = -t708 * mrSges(6,2) + t721 * mrSges(6,3);
t759 = -m(6) * t634 + t703 * mrSges(6,1) + t721 * t683 - t617;
t616 = m(5) * t636 + t703 * mrSges(5,1) + t721 * t684 + t788 * t691 + t775 * t709 + t759;
t766 = t783 * t614 - t745 * t616;
t610 = m(4) * t654 - qJDD(2) * mrSges(4,2) - t703 * mrSges(4,3) - qJD(2) * t713 - t721 * t697 + t766;
t772 = qJD(3) * t722;
t716 = -0.2e1 * t772;
t776 = t784 * t676 - t746 * t680;
t653 = t716 + t776;
t712 = -qJD(2) * mrSges(4,2) - t721 * mrSges(4,3);
t758 = qJDD(2) * pkin(3) + t753 * qJ(4) - t722 * t696 - qJDD(4) + t776;
t640 = 0.2e1 * t772 - t758;
t635 = -0.2e1 * qJD(5) * t709 + t791 + (t709 * t721 + t690) * pkin(4) + t640;
t632 = -t707 * pkin(8) + t716 + (-pkin(4) - pkin(5)) * t690 - t791 + (-pkin(4) * t721 + t687 + t790) * t709 + t758;
t762 = -m(7) * t632 + t646 * mrSges(7,1) - t647 * mrSges(7,2) + t672 * t658 - t673 * t659;
t625 = m(6) * t635 + t690 * mrSges(6,1) - t691 * mrSges(6,3) + t708 * t683 - t709 * t686 + t762;
t755 = -m(5) * t640 - t690 * mrSges(5,1) - t691 * mrSges(5,2) - t708 * t684 - t709 * t685 - t625;
t621 = m(4) * t653 + qJDD(2) * mrSges(4,1) - t704 * mrSges(4,3) + qJD(2) * t712 - t722 * t697 + t755;
t605 = t746 * t610 + t784 * t621;
t710 = -t751 * g(3) - t781;
t732 = (-mrSges(3,1) * t751 + mrSges(3,2) * t748) * qJD(1);
t737 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t773;
t602 = m(3) * t710 + qJDD(2) * mrSges(3,1) - t733 * mrSges(3,3) + qJD(2) * t737 - t732 * t774 + t605;
t736 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t774;
t767 = t784 * t610 - t746 * t621;
t603 = m(3) * t711 - qJDD(2) * mrSges(3,2) + t734 * mrSges(3,3) - qJD(2) * t736 + t732 * t773 + t767;
t768 = -t748 * t602 + t751 * t603;
t596 = m(2) * t739 - t754 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t768;
t727 = -t754 * pkin(7) + t761;
t611 = t745 * t614 + t783 * t616;
t757 = m(4) * t682 + t703 * mrSges(4,1) + t704 * mrSges(4,2) + t721 * t712 + t722 * t713 + t611;
t756 = -m(3) * t727 + t734 * mrSges(3,1) - t733 * mrSges(3,2) - t736 * t774 + t737 * t773 - t757;
t607 = m(2) * t738 + qJDD(1) * mrSges(2,1) - t754 * mrSges(2,2) + t756;
t780 = t749 * t596 + t752 * t607;
t597 = t751 * t602 + t748 * t603;
t779 = t793 * t708 - t787 * t709 - t785 * t721;
t778 = t785 * t708 + t794 * t709 + t792 * t721;
t777 = -t787 * t708 + t795 * t709 - t794 * t721;
t769 = t752 * t596 - t749 * t607;
t725 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t748 + Ifges(3,4) * t751) * qJD(1);
t724 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t748 + Ifges(3,2) * t751) * qJD(1);
t723 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t748 + Ifges(3,6) * t751) * qJD(1);
t694 = Ifges(4,1) * t722 - Ifges(4,4) * t721 + Ifges(4,5) * qJD(2);
t693 = Ifges(4,4) * t722 - Ifges(4,2) * t721 + Ifges(4,6) * qJD(2);
t692 = Ifges(4,5) * t722 - Ifges(4,6) * t721 + Ifges(4,3) * qJD(2);
t650 = Ifges(7,1) * t673 + Ifges(7,4) * t672 + Ifges(7,5) * t719;
t649 = Ifges(7,4) * t673 + Ifges(7,2) * t672 + Ifges(7,6) * t719;
t648 = Ifges(7,5) * t673 + Ifges(7,6) * t672 + Ifges(7,3) * t719;
t619 = mrSges(7,2) * t632 - mrSges(7,3) * t626 + Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t702 + t672 * t648 - t719 * t649;
t618 = -mrSges(7,1) * t632 + mrSges(7,3) * t627 + Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t702 - t673 * t648 + t719 * t650;
t604 = mrSges(5,2) * t640 + mrSges(6,2) * t634 - mrSges(5,3) * t636 - mrSges(6,3) * t635 - pkin(8) * t617 - qJ(5) * t625 - t747 * t618 + t750 * t619 - t787 * t690 + t795 * t691 - t703 * t794 + t778 * t708 + t779 * t721;
t598 = -mrSges(5,1) * t640 - mrSges(6,1) * t635 + mrSges(6,2) * t633 + mrSges(5,3) * t637 - pkin(4) * t625 - pkin(5) * t762 - pkin(8) * t765 - t750 * t618 - t747 * t619 - t793 * t690 + t787 * t691 + t785 * t703 + t778 * t709 + t777 * t721;
t593 = (-Ifges(4,2) + t792) * t703 - pkin(4) * t759 + (qJ(5) * t678 - t777) * t708 + (pkin(4) * t678 + t779) * t709 + (qJ(5) * mrSges(6,2) + t785) * t690 + (pkin(4) * mrSges(6,2) + t794) * t691 + Ifges(4,6) * qJDD(2) - t722 * t692 + Ifges(4,4) * t704 - qJ(5) * t760 + qJD(2) * t694 + Ifges(7,3) * t702 - mrSges(4,1) * t682 - t672 * t650 + t673 * t649 + mrSges(4,3) * t654 + Ifges(7,6) * t646 + Ifges(7,5) * t647 - mrSges(5,1) * t636 + mrSges(5,2) * t637 - mrSges(6,3) * t633 + mrSges(6,1) * t634 + mrSges(7,1) * t626 - mrSges(7,2) * t627 + pkin(5) * t617 - pkin(3) * t611;
t592 = mrSges(4,2) * t682 - mrSges(4,3) * t653 + Ifges(4,1) * t704 - Ifges(4,4) * t703 + Ifges(4,5) * qJDD(2) - qJ(4) * t611 - qJD(2) * t693 - t745 * t598 + t783 * t604 - t721 * t692;
t591 = -pkin(1) * t597 + mrSges(2,3) * t739 - pkin(2) * t605 - Ifges(3,5) * t733 - Ifges(3,6) * t734 - mrSges(3,1) * t710 + mrSges(3,2) * t711 - qJ(4) * t766 - t745 * t604 - t783 * t598 - pkin(3) * t755 - Ifges(4,5) * t704 + Ifges(4,6) * t703 - mrSges(4,1) * t653 + mrSges(4,2) * t654 + mrSges(2,1) * g(3) + t754 * Ifges(2,5) - t722 * t693 - t721 * t694 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t748 * t724 + t751 * t725) * qJD(1);
t590 = mrSges(3,2) * t727 - mrSges(3,3) * t710 + Ifges(3,1) * t733 + Ifges(3,4) * t734 + Ifges(3,5) * qJDD(2) - qJ(3) * t605 - qJD(2) * t724 + t784 * t592 - t746 * t593 + t723 * t773;
t589 = -mrSges(3,1) * t727 + mrSges(3,3) * t711 + Ifges(3,4) * t733 + Ifges(3,2) * t734 + Ifges(3,6) * qJDD(2) - pkin(2) * t757 + qJ(3) * t767 + qJD(2) * t725 + t746 * t592 + t784 * t593 - t723 * t774;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t738 + Ifges(2,5) * qJDD(1) - t754 * Ifges(2,6) - pkin(7) * t597 - t748 * t589 + t751 * t590;
t1 = [-m(1) * g(1) + t769; -m(1) * g(2) + t780; (-m(1) - m(2)) * g(3) + t597; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t780 + t752 * t588 - t749 * t591; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t769 + t749 * t588 + t752 * t591; -mrSges(1,1) * g(2) + mrSges(2,1) * t738 + mrSges(1,2) * g(1) - mrSges(2,2) * t739 + Ifges(2,3) * qJDD(1) + pkin(1) * t756 + pkin(7) * t768 + t751 * t589 + t748 * t590;];
tauB  = t1;
