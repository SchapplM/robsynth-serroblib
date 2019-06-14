% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:16:14
% EndTime: 2019-05-05 06:16:29
% DurationCPUTime: 14.80s
% Computational Cost: add. (242144->339), mult. (519937->432), div. (0->0), fcn. (384402->14), ass. (0->155)
t787 = -2 * qJD(4);
t786 = Ifges(4,1) + Ifges(5,2);
t779 = Ifges(4,4) + Ifges(5,6);
t778 = Ifges(4,5) - Ifges(5,4);
t785 = Ifges(4,2) + Ifges(5,3);
t777 = Ifges(4,6) - Ifges(5,5);
t784 = Ifges(4,3) + Ifges(5,1);
t729 = sin(pkin(12));
t732 = cos(pkin(12));
t715 = t729 * g(1) - t732 * g(2);
t728 = -g(3) + qJDD(1);
t731 = sin(pkin(6));
t734 = cos(pkin(6));
t694 = -t731 * t715 + t734 * t728;
t730 = sin(pkin(7));
t741 = cos(qJ(3));
t716 = -t732 * g(1) - t729 * g(2);
t738 = sin(qJ(2));
t742 = cos(qJ(2));
t768 = t734 * t742;
t771 = t731 * t742;
t671 = t715 * t768 - t738 * t716 + t728 * t771;
t743 = qJD(2) ^ 2;
t667 = t743 * t730 * pkin(9) + qJDD(2) * pkin(2) + t671;
t733 = cos(pkin(7));
t775 = t667 * t733;
t783 = t741 * (t694 * t730 + t775);
t725 = t733 * qJD(2) + qJD(3);
t737 = sin(qJ(3));
t763 = qJD(2) * t730;
t757 = t737 * t763;
t782 = (pkin(3) * t725 + t787) * t757;
t769 = t734 * t738;
t772 = t731 * t738;
t672 = t715 * t769 + t742 * t716 + t728 * t772;
t760 = qJDD(2) * t730;
t668 = -t743 * pkin(2) + pkin(9) * t760 + t672;
t770 = t733 * t737;
t773 = t730 * t737;
t644 = t667 * t770 + t741 * t668 + t694 * t773;
t702 = (-pkin(3) * t741 - qJ(4) * t737) * t763;
t723 = t725 ^ 2;
t724 = t733 * qJDD(2) + qJDD(3);
t762 = qJD(2) * t741;
t758 = t730 * t762;
t639 = t723 * pkin(3) - t724 * qJ(4) - t702 * t758 + t725 * t787 - t644;
t781 = -pkin(3) - pkin(10);
t780 = mrSges(4,1) - mrSges(5,2);
t665 = t737 * t668;
t643 = -t665 + t783;
t699 = -t725 * mrSges(4,2) + mrSges(4,3) * t758;
t700 = -mrSges(5,1) * t758 - t725 * mrSges(5,3);
t703 = (mrSges(5,2) * t741 - mrSges(5,3) * t737) * t763;
t704 = (-mrSges(4,1) * t741 + mrSges(4,2) * t737) * t763;
t706 = (qJD(3) * t762 + qJDD(2) * t737) * t730;
t752 = -t723 * qJ(4) + t702 * t757 + qJDD(4) + t665;
t774 = t730 ^ 2 * t743;
t636 = t706 * pkin(4) + t781 * t724 + (-pkin(10) * t737 * t774 - t775 + (-pkin(4) * qJD(2) * t725 - t694) * t730) * t741 + t752;
t690 = t733 * t694;
t705 = pkin(4) * t757 - t725 * pkin(10);
t707 = -qJD(3) * t757 + t741 * t760;
t759 = t741 ^ 2 * t774;
t638 = -pkin(4) * t759 - t706 * qJ(4) + t690 + t781 * t707 + (-t667 + (-qJ(4) * t725 * t741 - t705 * t737) * qJD(2)) * t730 + t782;
t736 = sin(qJ(5));
t740 = cos(qJ(5));
t631 = t736 * t636 + t740 * t638;
t693 = t740 * t725 - t736 * t758;
t662 = -t693 * qJD(5) - t740 * t707 - t736 * t724;
t692 = -t736 * t725 - t740 * t758;
t669 = -t692 * mrSges(6,1) + t693 * mrSges(6,2);
t714 = qJD(5) + t757;
t676 = t714 * mrSges(6,1) - t693 * mrSges(6,3);
t697 = qJDD(5) + t706;
t670 = -t692 * pkin(5) - t693 * pkin(11);
t712 = t714 ^ 2;
t629 = -t712 * pkin(5) + t697 * pkin(11) + t692 * t670 + t631;
t634 = t707 * pkin(4) - pkin(10) * t759 + t725 * t705 - t639;
t663 = t692 * qJD(5) - t736 * t707 + t740 * t724;
t632 = (-t692 * t714 - t663) * pkin(11) + (t693 * t714 - t662) * pkin(5) + t634;
t735 = sin(qJ(6));
t739 = cos(qJ(6));
t626 = -t735 * t629 + t739 * t632;
t673 = -t735 * t693 + t739 * t714;
t647 = t673 * qJD(6) + t739 * t663 + t735 * t697;
t674 = t739 * t693 + t735 * t714;
t652 = -t673 * mrSges(7,1) + t674 * mrSges(7,2);
t691 = qJD(6) - t692;
t654 = -t691 * mrSges(7,2) + t673 * mrSges(7,3);
t660 = qJDD(6) - t662;
t624 = m(7) * t626 + t660 * mrSges(7,1) - t647 * mrSges(7,3) - t674 * t652 + t691 * t654;
t627 = t739 * t629 + t735 * t632;
t646 = -t674 * qJD(6) - t735 * t663 + t739 * t697;
t655 = t691 * mrSges(7,1) - t674 * mrSges(7,3);
t625 = m(7) * t627 - t660 * mrSges(7,2) + t646 * mrSges(7,3) + t673 * t652 - t691 * t655;
t754 = -t735 * t624 + t739 * t625;
t616 = m(6) * t631 - t697 * mrSges(6,2) + t662 * mrSges(6,3) + t692 * t669 - t714 * t676 + t754;
t630 = t740 * t636 - t736 * t638;
t675 = -t714 * mrSges(6,2) + t692 * mrSges(6,3);
t628 = -t697 * pkin(5) - t712 * pkin(11) + t693 * t670 - t630;
t746 = -m(7) * t628 + t646 * mrSges(7,1) - t647 * mrSges(7,2) + t673 * t654 - t674 * t655;
t620 = m(6) * t630 + t697 * mrSges(6,1) - t663 * mrSges(6,3) - t693 * t669 + t714 * t675 + t746;
t610 = t736 * t616 + t740 * t620;
t640 = -t724 * pkin(3) + t752 - t783;
t747 = -m(5) * t640 - t706 * mrSges(5,1) - t610;
t606 = m(4) * t643 - t706 * mrSges(4,3) + (t699 - t700) * t725 + t780 * t724 + (-t703 - t704) * t757 + t747;
t776 = t606 * t741;
t653 = -t730 * t667 + t690;
t698 = t725 * mrSges(4,1) - mrSges(4,3) * t757;
t701 = mrSges(5,1) * t757 + t725 * mrSges(5,2);
t642 = -t707 * pkin(3) + (-t725 * t758 - t706) * qJ(4) + t653 + t782;
t755 = t740 * t616 - t736 * t620;
t751 = m(5) * t642 - t706 * mrSges(5,3) + t700 * t758 + t755;
t608 = m(4) * t653 + t706 * mrSges(4,2) - t780 * t707 + (-t699 * t741 + (t698 - t701) * t737) * t763 + t751;
t617 = t739 * t624 + t735 * t625;
t745 = -m(6) * t634 + t662 * mrSges(6,1) - t663 * mrSges(6,2) + t692 * t675 - t693 * t676 - t617;
t744 = -m(5) * t639 + t724 * mrSges(5,3) + t725 * t701 + t703 * t758 - t745;
t614 = t744 + (mrSges(4,3) + mrSges(5,1)) * t707 + m(4) * t644 - t724 * mrSges(4,2) - t725 * t698 + t704 * t758;
t596 = -t730 * t608 + t614 * t770 + t733 * t776;
t592 = m(3) * t671 + qJDD(2) * mrSges(3,1) - t743 * mrSges(3,2) + t596;
t595 = t733 * t608 + t614 * t773 + t730 * t776;
t594 = m(3) * t694 + t595;
t601 = -t737 * t606 + t741 * t614;
t600 = m(3) * t672 - t743 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t601;
t582 = t592 * t768 - t731 * t594 + t600 * t769;
t580 = m(2) * t715 + t582;
t588 = -t738 * t592 + t742 * t600;
t587 = m(2) * t716 + t588;
t767 = t732 * t580 + t729 * t587;
t766 = (t737 * t778 + t741 * t777) * t763 + t784 * t725;
t765 = (-t737 * t779 - t741 * t785) * t763 - t777 * t725;
t764 = (t737 * t786 + t741 * t779) * t763 + t778 * t725;
t581 = t592 * t771 + t734 * t594 + t600 * t772;
t756 = -t729 * t580 + t732 * t587;
t648 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t691;
t650 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t691;
t618 = -mrSges(7,1) * t628 + mrSges(7,3) * t627 + Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t660 - t674 * t648 + t691 * t650;
t649 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t691;
t619 = mrSges(7,2) * t628 - mrSges(7,3) * t626 + Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t660 + t673 * t648 - t691 * t649;
t656 = Ifges(6,5) * t693 + Ifges(6,6) * t692 + Ifges(6,3) * t714;
t657 = Ifges(6,4) * t693 + Ifges(6,2) * t692 + Ifges(6,6) * t714;
t602 = mrSges(6,2) * t634 - mrSges(6,3) * t630 + Ifges(6,1) * t663 + Ifges(6,4) * t662 + Ifges(6,5) * t697 - pkin(11) * t617 - t735 * t618 + t739 * t619 + t692 * t656 - t714 * t657;
t658 = Ifges(6,1) * t693 + Ifges(6,4) * t692 + Ifges(6,5) * t714;
t603 = -mrSges(6,1) * t634 - mrSges(7,1) * t626 + mrSges(7,2) * t627 + mrSges(6,3) * t631 + Ifges(6,4) * t663 - Ifges(7,5) * t647 + Ifges(6,2) * t662 + Ifges(6,6) * t697 - Ifges(7,6) * t646 - Ifges(7,3) * t660 - pkin(5) * t617 - t674 * t649 + t673 * t650 - t693 * t656 + t714 * t658;
t583 = mrSges(4,1) * t643 - mrSges(4,2) * t644 + mrSges(5,2) * t640 - mrSges(5,3) * t639 + t740 * t602 - t736 * t603 - pkin(10) * t610 + pkin(3) * (-t725 * t700 + t747) + qJ(4) * t744 + (-pkin(3) * mrSges(5,2) + t784) * t724 + (qJ(4) * mrSges(5,1) + t777) * t707 + t778 * t706 + (-t764 * t741 + (-pkin(3) * t703 - t765) * t737) * t763;
t609 = t707 * mrSges(5,2) - t701 * t757 + t751;
t584 = -mrSges(4,1) * t653 - mrSges(5,1) * t639 + mrSges(5,2) * t642 + mrSges(4,3) * t644 - pkin(3) * t609 - pkin(4) * t745 - pkin(10) * t755 - t736 * t602 - t740 * t603 + t779 * t706 + t707 * t785 + t777 * t724 + t764 * t725 - t766 * t757;
t589 = mrSges(5,1) * t640 - mrSges(5,3) * t642 + pkin(5) * t746 - qJ(4) * t609 + mrSges(6,1) * t630 - mrSges(6,2) * t631 - mrSges(4,3) * t643 + mrSges(4,2) * t653 - t692 * t658 + t693 * t657 + t735 * t619 + pkin(11) * t754 + t739 * t618 + pkin(4) * t610 + Ifges(6,6) * t662 + Ifges(6,5) * t663 + Ifges(6,3) * t697 + t765 * t725 + t778 * t724 + t779 * t707 + t786 * t706 + t766 * t758;
t748 = pkin(9) * t601 + t584 * t741 + t589 * t737;
t577 = -mrSges(3,1) * t694 + mrSges(3,3) * t672 + t743 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t595 - t730 * t583 + t733 * t748;
t578 = mrSges(3,2) * t694 - mrSges(3,3) * t671 + Ifges(3,5) * qJDD(2) - t743 * Ifges(3,6) - t737 * t584 + t741 * t589 + (-t595 * t730 - t596 * t733) * pkin(9);
t749 = pkin(8) * t588 + t577 * t742 + t578 * t738;
t576 = mrSges(3,1) * t671 - mrSges(3,2) * t672 + Ifges(3,3) * qJDD(2) + pkin(2) * t596 + t733 * t583 + t730 * t748;
t575 = mrSges(2,2) * t728 - mrSges(2,3) * t715 - t738 * t577 + t742 * t578 + (-t581 * t731 - t582 * t734) * pkin(8);
t574 = -mrSges(2,1) * t728 + mrSges(2,3) * t716 - pkin(1) * t581 - t731 * t576 + t734 * t749;
t1 = [-m(1) * g(1) + t756; -m(1) * g(2) + t767; -m(1) * g(3) + m(2) * t728 + t581; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t767 - t729 * t574 + t732 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t756 + t732 * t574 + t729 * t575; -mrSges(1,1) * g(2) + mrSges(2,1) * t715 + mrSges(1,2) * g(1) - mrSges(2,2) * t716 + pkin(1) * t582 + t734 * t576 + t731 * t749;];
tauB  = t1;
