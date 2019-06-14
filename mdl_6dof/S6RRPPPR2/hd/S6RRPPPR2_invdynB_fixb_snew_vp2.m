% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-05-06 08:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:22:47
% EndTime: 2019-05-06 08:22:57
% DurationCPUTime: 8.58s
% Computational Cost: add. (114841->367), mult. (273546->447), div. (0->0), fcn. (185882->10), ass. (0->143)
t786 = -2 * qJD(4);
t785 = -Ifges(5,1) - Ifges(4,3);
t780 = Ifges(4,5) - Ifges(5,4);
t784 = Ifges(4,2) + Ifges(5,3);
t783 = Ifges(5,2) + Ifges(4,1);
t779 = Ifges(4,6) - Ifges(5,5);
t778 = -Ifges(5,6) - Ifges(4,4);
t739 = sin(pkin(9));
t745 = cos(qJ(2));
t768 = qJD(1) * t745;
t742 = sin(qJ(2));
t769 = qJD(1) * t742;
t777 = cos(pkin(9));
t713 = t739 * t769 - t777 * t768;
t714 = (t739 * t745 + t777 * t742) * qJD(1);
t682 = t713 * pkin(3) - t714 * qJ(4);
t747 = qJD(2) ^ 2;
t764 = qJD(1) * qJD(2);
t725 = t742 * qJDD(1) + t745 * t764;
t743 = sin(qJ(1));
t746 = cos(qJ(1));
t731 = -t746 * g(1) - t743 * g(2);
t748 = qJD(1) ^ 2;
t720 = -t748 * pkin(1) + qJDD(1) * pkin(7) + t731;
t776 = t742 * t720;
t781 = pkin(2) * t748;
t660 = qJDD(2) * pkin(2) - t725 * qJ(3) - t776 + (qJ(3) * t764 + t742 * t781 - g(3)) * t745;
t700 = -t742 * g(3) + t745 * t720;
t726 = t745 * qJDD(1) - t742 * t764;
t727 = qJD(2) * pkin(2) - qJ(3) * t769;
t737 = t745 ^ 2;
t662 = t726 * qJ(3) - qJD(2) * t727 - t737 * t781 + t700;
t773 = t739 * t660 + t777 * t662;
t782 = t747 * pkin(3) - qJDD(2) * qJ(4) + qJD(2) * t786 + t713 * t682 - t773;
t646 = -0.2e1 * qJD(3) * t714 + t777 * t660 - t739 * t662;
t683 = t713 * mrSges(4,1) + t714 * mrSges(4,2);
t691 = t777 * t725 + t739 * t726;
t701 = -qJD(2) * mrSges(4,2) - t713 * mrSges(4,3);
t704 = t713 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t636 = -qJDD(2) * pkin(3) - t747 * qJ(4) + t714 * t682 + qJDD(4) - t646;
t767 = qJD(2) * t713;
t630 = (t713 * t714 - qJDD(2)) * qJ(5) + (t691 + t767) * pkin(4) + t636;
t690 = t739 * t725 - t777 * t726;
t703 = t714 * pkin(4) - qJD(2) * qJ(5);
t712 = t713 ^ 2;
t730 = t743 * g(1) - t746 * g(2);
t757 = -qJDD(1) * pkin(1) - t730;
t666 = -t726 * pkin(2) + qJDD(3) + t727 * t769 + (-qJ(3) * t737 - pkin(7)) * t748 + t757;
t750 = (-t691 + t767) * qJ(4) + t666 + (qJD(2) * pkin(3) + t786) * t714;
t634 = -t712 * pkin(4) - t714 * t703 + (pkin(3) + qJ(5)) * t690 + t750;
t738 = sin(pkin(10));
t740 = cos(pkin(10));
t696 = t740 * qJD(2) + t738 * t713;
t624 = -0.2e1 * qJD(5) * t696 + t740 * t630 - t738 * t634;
t672 = t740 * qJDD(2) + t738 * t690;
t695 = -t738 * qJD(2) + t740 * t713;
t622 = (t695 * t714 - t672) * pkin(8) + (t695 * t696 + t691) * pkin(5) + t624;
t625 = 0.2e1 * qJD(5) * t695 + t738 * t630 + t740 * t634;
t669 = t714 * pkin(5) - t696 * pkin(8);
t671 = -t738 * qJDD(2) + t740 * t690;
t694 = t695 ^ 2;
t623 = -t694 * pkin(5) + t671 * pkin(8) - t714 * t669 + t625;
t741 = sin(qJ(6));
t744 = cos(qJ(6));
t620 = t744 * t622 - t741 * t623;
t656 = t744 * t695 - t741 * t696;
t641 = t656 * qJD(6) + t741 * t671 + t744 * t672;
t657 = t741 * t695 + t744 * t696;
t648 = -t656 * mrSges(7,1) + t657 * mrSges(7,2);
t711 = qJD(6) + t714;
t649 = -t711 * mrSges(7,2) + t656 * mrSges(7,3);
t689 = qJDD(6) + t691;
t618 = m(7) * t620 + t689 * mrSges(7,1) - t641 * mrSges(7,3) - t657 * t648 + t711 * t649;
t621 = t741 * t622 + t744 * t623;
t640 = -t657 * qJD(6) + t744 * t671 - t741 * t672;
t650 = t711 * mrSges(7,1) - t657 * mrSges(7,3);
t619 = m(7) * t621 - t689 * mrSges(7,2) + t640 * mrSges(7,3) + t656 * t648 - t711 * t650;
t609 = t744 * t618 + t741 * t619;
t661 = -t695 * mrSges(6,1) + t696 * mrSges(6,2);
t667 = -t714 * mrSges(6,2) + t695 * mrSges(6,3);
t607 = m(6) * t624 + t691 * mrSges(6,1) - t672 * mrSges(6,3) - t696 * t661 + t714 * t667 + t609;
t668 = t714 * mrSges(6,1) - t696 * mrSges(6,3);
t759 = -t741 * t618 + t744 * t619;
t608 = m(6) * t625 - t691 * mrSges(6,2) + t671 * mrSges(6,3) + t695 * t661 - t714 * t668 + t759;
t604 = t740 * t607 + t738 * t608;
t684 = -t713 * mrSges(5,2) - t714 * mrSges(5,3);
t754 = -m(5) * t636 - t691 * mrSges(5,1) - t714 * t684 - t604;
t602 = m(4) * t646 - t691 * mrSges(4,3) - t714 * t683 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t701 - t704) * qJD(2) + t754;
t766 = qJD(3) * t713;
t708 = -0.2e1 * t766;
t647 = t708 + t773;
t702 = qJD(2) * mrSges(4,1) - t714 * mrSges(4,3);
t635 = 0.2e1 * t766 + t782;
t705 = t714 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t632 = -t690 * pkin(4) - t712 * qJ(5) + qJD(2) * t703 + qJDD(5) + t708 - t782;
t627 = -t671 * pkin(5) - t694 * pkin(8) + t696 * t669 + t632;
t755 = m(7) * t627 - t640 * mrSges(7,1) + t641 * mrSges(7,2) - t656 * t649 + t657 * t650;
t753 = -m(6) * t632 + t671 * mrSges(6,1) - t672 * mrSges(6,2) + t695 * t667 - t696 * t668 - t755;
t751 = -m(5) * t635 + qJDD(2) * mrSges(5,3) + qJD(2) * t705 - t753;
t614 = (-mrSges(4,3) - mrSges(5,1)) * t690 + (-t683 - t684) * t713 - qJDD(2) * mrSges(4,2) + t751 - qJD(2) * t702 + m(4) * t647;
t596 = t777 * t602 + t739 * t614;
t699 = -t745 * g(3) - t776;
t724 = (-mrSges(3,1) * t745 + mrSges(3,2) * t742) * qJD(1);
t729 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t768;
t594 = m(3) * t699 + qJDD(2) * mrSges(3,1) - t725 * mrSges(3,3) + qJD(2) * t729 - t724 * t769 + t596;
t728 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t769;
t760 = -t739 * t602 + t777 * t614;
t595 = m(3) * t700 - qJDD(2) * mrSges(3,2) + t726 * mrSges(3,3) - qJD(2) * t728 + t724 * t768 + t760;
t761 = -t742 * t594 + t745 * t595;
t589 = m(2) * t731 - t748 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t761;
t719 = -t748 * pkin(7) + t757;
t638 = t690 * pkin(3) + t750;
t774 = -t738 * t607 + t740 * t608;
t603 = m(5) * t638 - t690 * mrSges(5,2) - t691 * mrSges(5,3) - t713 * t704 - t714 * t705 + t774;
t752 = m(4) * t666 + t690 * mrSges(4,1) + t691 * mrSges(4,2) + t713 * t701 + t714 * t702 + t603;
t749 = -m(3) * t719 + t726 * mrSges(3,1) - t725 * mrSges(3,2) - t728 * t769 + t729 * t768 - t752;
t600 = m(2) * t730 + qJDD(1) * mrSges(2,1) - t748 * mrSges(2,2) + t749;
t775 = t743 * t589 + t746 * t600;
t590 = t745 * t594 + t742 * t595;
t772 = -t779 * qJD(2) + t784 * t713 + t778 * t714;
t771 = t785 * qJD(2) + t779 * t713 - t780 * t714;
t770 = t780 * qJD(2) + t778 * t713 + t783 * t714;
t762 = t746 * t589 - t743 * t600;
t717 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t742 + Ifges(3,4) * t745) * qJD(1);
t716 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t742 + Ifges(3,2) * t745) * qJD(1);
t715 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t742 + Ifges(3,6) * t745) * qJD(1);
t653 = Ifges(6,1) * t696 + Ifges(6,4) * t695 + Ifges(6,5) * t714;
t652 = Ifges(6,4) * t696 + Ifges(6,2) * t695 + Ifges(6,6) * t714;
t651 = Ifges(6,5) * t696 + Ifges(6,6) * t695 + Ifges(6,3) * t714;
t644 = Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t711;
t643 = Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t711;
t642 = Ifges(7,5) * t657 + Ifges(7,6) * t656 + Ifges(7,3) * t711;
t611 = mrSges(7,2) * t627 - mrSges(7,3) * t620 + Ifges(7,1) * t641 + Ifges(7,4) * t640 + Ifges(7,5) * t689 + t656 * t642 - t711 * t643;
t610 = -mrSges(7,1) * t627 + mrSges(7,3) * t621 + Ifges(7,4) * t641 + Ifges(7,2) * t640 + Ifges(7,6) * t689 - t657 * t642 + t711 * t644;
t598 = mrSges(6,2) * t632 - mrSges(6,3) * t624 + Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t691 - pkin(8) * t609 - t741 * t610 + t744 * t611 + t695 * t651 - t714 * t652;
t597 = -mrSges(6,1) * t632 + mrSges(6,3) * t625 + Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t691 - pkin(5) * t755 + pkin(8) * t759 + t744 * t610 + t741 * t611 - t696 * t651 + t714 * t653;
t586 = pkin(5) * t609 - qJ(4) * t603 + (Ifges(6,3) + t783) * t691 + t778 * t690 + t780 * qJDD(2) + t771 * t713 + t772 * qJD(2) + pkin(4) * t604 - t695 * t653 + t696 * t652 + Ifges(7,3) * t689 + Ifges(6,6) * t671 + Ifges(6,5) * t672 + mrSges(4,2) * t666 - t656 * t644 + t657 * t643 - mrSges(4,3) * t646 - mrSges(5,3) * t638 + Ifges(7,6) * t640 + Ifges(7,5) * t641 + mrSges(5,1) * t636 - mrSges(6,2) * t625 + mrSges(6,1) * t624 - mrSges(7,2) * t621 + mrSges(7,1) * t620;
t585 = -mrSges(4,1) * t666 - mrSges(5,1) * t635 + mrSges(5,2) * t638 + mrSges(4,3) * t647 - pkin(3) * t603 - pkin(4) * t753 - qJ(5) * t774 + t770 * qJD(2) + t779 * qJDD(2) - t740 * t597 - t738 * t598 - t784 * t690 - t778 * t691 + t771 * t714;
t584 = -pkin(1) * t590 + mrSges(2,1) * g(3) + (pkin(3) * mrSges(5,2) - Ifges(3,3) + t785) * qJDD(2) + (qJ(4) * mrSges(5,1) + t779) * t690 - t780 * t691 + (qJ(4) * t684 - t770) * t713 + t772 * t714 - qJ(4) * t751 - pkin(3) * (-qJD(2) * t704 + t754) + qJ(5) * t604 + (-t742 * t716 + t745 * t717) * qJD(1) - pkin(2) * t596 + Ifges(2,6) * qJDD(1) + t748 * Ifges(2,5) + t738 * t597 - t740 * t598 + mrSges(2,3) * t731 - Ifges(3,5) * t725 - Ifges(3,6) * t726 - mrSges(3,1) * t699 + mrSges(3,2) * t700 - mrSges(4,1) * t646 + mrSges(4,2) * t647 + mrSges(5,3) * t635 - mrSges(5,2) * t636;
t583 = mrSges(3,2) * t719 - mrSges(3,3) * t699 + Ifges(3,1) * t725 + Ifges(3,4) * t726 + Ifges(3,5) * qJDD(2) - qJ(3) * t596 - qJD(2) * t716 - t739 * t585 + t777 * t586 + t715 * t768;
t582 = -mrSges(3,1) * t719 + mrSges(3,3) * t700 + Ifges(3,4) * t725 + Ifges(3,2) * t726 + Ifges(3,6) * qJDD(2) - pkin(2) * t752 + qJ(3) * t760 + qJD(2) * t717 + t777 * t585 + t739 * t586 - t715 * t769;
t581 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - t748 * Ifges(2,6) - pkin(7) * t590 - t742 * t582 + t745 * t583;
t1 = [-m(1) * g(1) + t762; -m(1) * g(2) + t775; (-m(1) - m(2)) * g(3) + t590; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t775 + t746 * t581 - t743 * t584; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t762 + t743 * t581 + t746 * t584; -mrSges(1,1) * g(2) + mrSges(2,1) * t730 + mrSges(1,2) * g(1) - mrSges(2,2) * t731 + Ifges(2,3) * qJDD(1) + pkin(1) * t749 + pkin(7) * t761 + t745 * t582 + t742 * t583;];
tauB  = t1;
