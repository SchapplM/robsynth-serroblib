% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP5
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
% Datum: 2019-05-08 05:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:58:31
% EndTime: 2019-05-08 04:59:05
% DurationCPUTime: 16.58s
% Computational Cost: add. (256491->364), mult. (525775->441), div. (0->0), fcn. (384291->10), ass. (0->140)
t770 = Ifges(6,1) + Ifges(7,1);
t767 = Ifges(6,4) + Ifges(7,4);
t766 = Ifges(6,5) + Ifges(7,5);
t769 = Ifges(6,2) + Ifges(7,2);
t765 = -Ifges(6,6) - Ifges(7,6);
t768 = -Ifges(6,3) - Ifges(7,3);
t737 = sin(qJ(1));
t742 = cos(qJ(1));
t727 = -g(1) * t742 - g(2) * t737;
t744 = qJD(1) ^ 2;
t711 = -pkin(1) * t744 + qJDD(1) * pkin(7) + t727;
t736 = sin(qJ(2));
t741 = cos(qJ(2));
t700 = -g(3) * t736 + t741 * t711;
t719 = (-mrSges(3,1) * t741 + mrSges(3,2) * t736) * qJD(1);
t758 = qJD(1) * qJD(2);
t730 = t736 * t758;
t722 = qJDD(1) * t741 - t730;
t760 = qJD(1) * t736;
t724 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t760;
t726 = t737 * g(1) - t742 * g(2);
t710 = -qJDD(1) * pkin(1) - t744 * pkin(7) - t726;
t755 = t741 * t758;
t721 = qJDD(1) * t736 + t755;
t678 = (-t721 - t755) * pkin(8) + (-t722 + t730) * pkin(2) + t710;
t720 = (-pkin(2) * t741 - pkin(8) * t736) * qJD(1);
t743 = qJD(2) ^ 2;
t759 = qJD(1) * t741;
t681 = -pkin(2) * t743 + qJDD(2) * pkin(8) + t720 * t759 + t700;
t735 = sin(qJ(3));
t740 = cos(qJ(3));
t659 = t740 * t678 - t735 * t681;
t717 = qJD(2) * t740 - t735 * t760;
t692 = qJD(3) * t717 + qJDD(2) * t735 + t721 * t740;
t716 = qJDD(3) - t722;
t718 = qJD(2) * t735 + t740 * t760;
t729 = qJD(3) - t759;
t640 = (t717 * t729 - t692) * pkin(9) + (t717 * t718 + t716) * pkin(3) + t659;
t660 = t735 * t678 + t740 * t681;
t691 = -qJD(3) * t718 + qJDD(2) * t740 - t721 * t735;
t701 = pkin(3) * t729 - pkin(9) * t718;
t715 = t717 ^ 2;
t642 = -pkin(3) * t715 + pkin(9) * t691 - t701 * t729 + t660;
t734 = sin(qJ(4));
t739 = cos(qJ(4));
t626 = t739 * t640 - t734 * t642;
t694 = t717 * t739 - t718 * t734;
t658 = qJD(4) * t694 + t691 * t734 + t692 * t739;
t695 = t717 * t734 + t718 * t739;
t712 = qJDD(4) + t716;
t728 = qJD(4) + t729;
t623 = (t694 * t728 - t658) * pkin(10) + (t694 * t695 + t712) * pkin(4) + t626;
t627 = t734 * t640 + t739 * t642;
t657 = -qJD(4) * t695 + t691 * t739 - t692 * t734;
t684 = pkin(4) * t728 - pkin(10) * t695;
t693 = t694 ^ 2;
t625 = -pkin(4) * t693 + pkin(10) * t657 - t684 * t728 + t627;
t733 = sin(qJ(5));
t738 = cos(qJ(5));
t617 = t738 * t623 - t733 * t625;
t673 = t694 * t738 - t695 * t733;
t636 = qJD(5) * t673 + t657 * t733 + t658 * t738;
t674 = t694 * t733 + t695 * t738;
t653 = -mrSges(7,1) * t673 + mrSges(7,2) * t674;
t654 = -mrSges(6,1) * t673 + mrSges(6,2) * t674;
t723 = qJD(5) + t728;
t663 = -mrSges(6,2) * t723 + mrSges(6,3) * t673;
t706 = qJDD(5) + t712;
t614 = -0.2e1 * qJD(6) * t674 + (t673 * t723 - t636) * qJ(6) + (t673 * t674 + t706) * pkin(5) + t617;
t662 = -mrSges(7,2) * t723 + mrSges(7,3) * t673;
t757 = m(7) * t614 + t706 * mrSges(7,1) + t723 * t662;
t606 = m(6) * t617 + t706 * mrSges(6,1) + t723 * t663 + (-t653 - t654) * t674 + (-mrSges(6,3) - mrSges(7,3)) * t636 + t757;
t618 = t733 * t623 + t738 * t625;
t635 = -qJD(5) * t674 + t657 * t738 - t658 * t733;
t665 = mrSges(7,1) * t723 - mrSges(7,3) * t674;
t666 = mrSges(6,1) * t723 - mrSges(6,3) * t674;
t664 = pkin(5) * t723 - qJ(6) * t674;
t672 = t673 ^ 2;
t616 = -pkin(5) * t672 + qJ(6) * t635 + 0.2e1 * qJD(6) * t673 - t664 * t723 + t618;
t756 = m(7) * t616 + t635 * mrSges(7,3) + t673 * t653;
t609 = m(6) * t618 + t635 * mrSges(6,3) + t673 * t654 + (-t665 - t666) * t723 + (-mrSges(6,2) - mrSges(7,2)) * t706 + t756;
t604 = t738 * t606 + t733 * t609;
t675 = -mrSges(5,1) * t694 + mrSges(5,2) * t695;
t682 = -mrSges(5,2) * t728 + mrSges(5,3) * t694;
t601 = m(5) * t626 + mrSges(5,1) * t712 - mrSges(5,3) * t658 - t675 * t695 + t682 * t728 + t604;
t683 = mrSges(5,1) * t728 - mrSges(5,3) * t695;
t750 = -t606 * t733 + t738 * t609;
t602 = m(5) * t627 - mrSges(5,2) * t712 + mrSges(5,3) * t657 + t675 * t694 - t683 * t728 + t750;
t596 = t739 * t601 + t734 * t602;
t696 = -mrSges(4,1) * t717 + mrSges(4,2) * t718;
t697 = -mrSges(4,2) * t729 + mrSges(4,3) * t717;
t594 = m(4) * t659 + mrSges(4,1) * t716 - mrSges(4,3) * t692 - t696 * t718 + t697 * t729 + t596;
t698 = mrSges(4,1) * t729 - mrSges(4,3) * t718;
t751 = -t601 * t734 + t739 * t602;
t595 = m(4) * t660 - mrSges(4,2) * t716 + mrSges(4,3) * t691 + t696 * t717 - t698 * t729 + t751;
t752 = -t594 * t735 + t740 * t595;
t589 = m(3) * t700 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t722 - qJD(2) * t724 + t719 * t759 + t752;
t699 = -t741 * g(3) - t736 * t711;
t725 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t759;
t680 = -qJDD(2) * pkin(2) - pkin(8) * t743 + t720 * t760 - t699;
t655 = -pkin(3) * t691 - pkin(9) * t715 + t718 * t701 + t680;
t629 = -pkin(4) * t657 - pkin(10) * t693 + t695 * t684 + t655;
t620 = -pkin(5) * t635 - qJ(6) * t672 + t664 * t674 + qJDD(6) + t629;
t749 = m(7) * t620 - t635 * mrSges(7,1) + t636 * mrSges(7,2) - t673 * t662 + t674 * t665;
t748 = m(6) * t629 - t635 * mrSges(6,1) + t636 * mrSges(6,2) - t673 * t663 + t674 * t666 + t749;
t746 = m(5) * t655 - t657 * mrSges(5,1) + t658 * mrSges(5,2) - t694 * t682 + t695 * t683 + t748;
t745 = -m(4) * t680 + t691 * mrSges(4,1) - t692 * mrSges(4,2) + t717 * t697 - t718 * t698 - t746;
t611 = m(3) * t699 + qJDD(2) * mrSges(3,1) - t721 * mrSges(3,3) + qJD(2) * t725 - t719 * t760 + t745;
t753 = t741 * t589 - t611 * t736;
t583 = m(2) * t727 - mrSges(2,1) * t744 - qJDD(1) * mrSges(2,2) + t753;
t590 = t594 * t740 + t595 * t735;
t747 = -m(3) * t710 + t722 * mrSges(3,1) - mrSges(3,2) * t721 - t724 * t760 + t725 * t759 - t590;
t586 = m(2) * t726 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t744 + t747;
t764 = t737 * t583 + t742 * t586;
t584 = t736 * t589 + t741 * t611;
t763 = t765 * t673 - t766 * t674 + t768 * t723;
t762 = -t769 * t673 - t767 * t674 + t765 * t723;
t761 = t767 * t673 + t770 * t674 + t766 * t723;
t754 = t742 * t583 - t586 * t737;
t709 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t736 + Ifges(3,4) * t741) * qJD(1);
t708 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t736 + Ifges(3,2) * t741) * qJD(1);
t707 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t736 + Ifges(3,6) * t741) * qJD(1);
t687 = Ifges(4,1) * t718 + Ifges(4,4) * t717 + Ifges(4,5) * t729;
t686 = Ifges(4,4) * t718 + Ifges(4,2) * t717 + Ifges(4,6) * t729;
t685 = Ifges(4,5) * t718 + Ifges(4,6) * t717 + Ifges(4,3) * t729;
t669 = Ifges(5,1) * t695 + Ifges(5,4) * t694 + Ifges(5,5) * t728;
t668 = Ifges(5,4) * t695 + Ifges(5,2) * t694 + Ifges(5,6) * t728;
t667 = Ifges(5,5) * t695 + Ifges(5,6) * t694 + Ifges(5,3) * t728;
t612 = -t636 * mrSges(7,3) - t674 * t653 + t757;
t603 = mrSges(6,2) * t629 + mrSges(7,2) * t620 - mrSges(6,3) * t617 - mrSges(7,3) * t614 - qJ(6) * t612 + t767 * t635 + t770 * t636 - t763 * t673 + t766 * t706 + t762 * t723;
t597 = -mrSges(6,1) * t629 + mrSges(6,3) * t618 - mrSges(7,1) * t620 + mrSges(7,3) * t616 - pkin(5) * t749 + qJ(6) * t756 + (-qJ(6) * t665 + t761) * t723 + (-mrSges(7,2) * qJ(6) - t765) * t706 + t763 * t674 + t767 * t636 + t769 * t635;
t592 = mrSges(5,2) * t655 - mrSges(5,3) * t626 + Ifges(5,1) * t658 + Ifges(5,4) * t657 + Ifges(5,5) * t712 - pkin(10) * t604 - t597 * t733 + t603 * t738 + t667 * t694 - t668 * t728;
t591 = -mrSges(5,1) * t655 + mrSges(5,3) * t627 + Ifges(5,4) * t658 + Ifges(5,2) * t657 + Ifges(5,6) * t712 - pkin(4) * t748 + pkin(10) * t750 + t738 * t597 + t733 * t603 - t695 * t667 + t728 * t669;
t580 = t765 * t635 - t766 * t636 + t761 * t673 + t762 * t674 - t707 * t760 + t768 * t706 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t721 + Ifges(3,2) * t722 - Ifges(4,3) * t716 + t717 * t687 - t718 * t686 + mrSges(3,3) * t700 + qJD(2) * t709 - mrSges(3,1) * t710 - Ifges(5,3) * t712 + t694 * t669 - t695 * t668 - Ifges(4,6) * t691 - Ifges(4,5) * t692 - Ifges(5,6) * t657 - Ifges(5,5) * t658 - mrSges(4,1) * t659 + mrSges(4,2) * t660 - mrSges(5,1) * t626 + mrSges(5,2) * t627 + mrSges(6,2) * t618 + mrSges(7,2) * t616 - mrSges(6,1) * t617 - mrSges(7,1) * t614 - pkin(5) * t612 - pkin(4) * t604 - pkin(3) * t596 - pkin(2) * t590;
t579 = mrSges(4,2) * t680 - mrSges(4,3) * t659 + Ifges(4,1) * t692 + Ifges(4,4) * t691 + Ifges(4,5) * t716 - pkin(9) * t596 - t591 * t734 + t592 * t739 + t685 * t717 - t686 * t729;
t578 = -mrSges(4,1) * t680 + mrSges(4,3) * t660 + Ifges(4,4) * t692 + Ifges(4,2) * t691 + Ifges(4,6) * t716 - pkin(3) * t746 + pkin(9) * t751 + t739 * t591 + t734 * t592 - t718 * t685 + t729 * t687;
t577 = mrSges(3,2) * t710 - mrSges(3,3) * t699 + Ifges(3,1) * t721 + Ifges(3,4) * t722 + Ifges(3,5) * qJDD(2) - pkin(8) * t590 - qJD(2) * t708 - t578 * t735 + t579 * t740 + t707 * t759;
t576 = Ifges(2,6) * qJDD(1) + t744 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t727 - Ifges(3,5) * t721 - Ifges(3,6) * t722 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t699 + mrSges(3,2) * t700 - t735 * t579 - t740 * t578 - pkin(2) * t745 - pkin(8) * t752 - pkin(1) * t584 + (-t708 * t736 + t709 * t741) * qJD(1);
t575 = -mrSges(2,2) * g(3) - mrSges(2,3) * t726 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t744 - pkin(7) * t584 + t577 * t741 - t580 * t736;
t1 = [-m(1) * g(1) + t754; -m(1) * g(2) + t764; (-m(1) - m(2)) * g(3) + t584; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t764 + t742 * t575 - t737 * t576; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t754 + t737 * t575 + t742 * t576; -mrSges(1,1) * g(2) + mrSges(2,1) * t726 + mrSges(1,2) * g(1) - mrSges(2,2) * t727 + Ifges(2,3) * qJDD(1) + pkin(1) * t747 + pkin(7) * t753 + t736 * t577 + t741 * t580;];
tauB  = t1;
