% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:23:26
% EndTime: 2019-05-06 17:23:44
% DurationCPUTime: 14.38s
% Computational Cost: add. (183865->363), mult. (425622->448), div. (0->0), fcn. (312732->10), ass. (0->139)
t748 = Ifges(6,1) + Ifges(7,1);
t741 = Ifges(6,4) - Ifges(7,5);
t747 = -Ifges(6,5) - Ifges(7,4);
t746 = Ifges(6,2) + Ifges(7,3);
t739 = Ifges(6,6) - Ifges(7,6);
t745 = -Ifges(6,3) - Ifges(7,2);
t744 = cos(qJ(5));
t717 = qJD(1) ^ 2;
t743 = pkin(2) * t717;
t742 = -mrSges(6,3) - mrSges(7,2);
t713 = sin(qJ(1));
t716 = cos(qJ(1));
t700 = -t716 * g(1) - t713 * g(2);
t689 = -t717 * pkin(1) + qJDD(1) * pkin(7) + t700;
t712 = sin(qJ(2));
t738 = t712 * t689;
t715 = cos(qJ(2));
t730 = qJD(1) * qJD(2);
t694 = t712 * qJDD(1) + t715 * t730;
t654 = qJDD(2) * pkin(2) - t694 * qJ(3) - t738 + (qJ(3) * t730 + t712 * t743 - g(3)) * t715;
t675 = -t712 * g(3) + t715 * t689;
t695 = t715 * qJDD(1) - t712 * t730;
t732 = qJD(1) * t712;
t696 = qJD(2) * pkin(2) - qJ(3) * t732;
t707 = t715 ^ 2;
t655 = t695 * qJ(3) - qJD(2) * t696 - t707 * t743 + t675;
t708 = sin(pkin(10));
t709 = cos(pkin(10));
t684 = (t708 * t715 + t709 * t712) * qJD(1);
t627 = -0.2e1 * qJD(3) * t684 + t709 * t654 - t708 * t655;
t673 = t709 * t694 + t708 * t695;
t683 = (-t708 * t712 + t709 * t715) * qJD(1);
t607 = (qJD(2) * t683 - t673) * pkin(8) + (t683 * t684 + qJDD(2)) * pkin(3) + t627;
t628 = 0.2e1 * qJD(3) * t683 + t708 * t654 + t709 * t655;
t672 = -t708 * t694 + t709 * t695;
t678 = qJD(2) * pkin(3) - t684 * pkin(8);
t682 = t683 ^ 2;
t610 = -t682 * pkin(3) + t672 * pkin(8) - qJD(2) * t678 + t628;
t711 = sin(qJ(4));
t714 = cos(qJ(4));
t605 = t711 * t607 + t714 * t610;
t667 = t711 * t683 + t714 * t684;
t635 = -t667 * qJD(4) + t714 * t672 - t711 * t673;
t666 = t714 * t683 - t711 * t684;
t649 = -t666 * mrSges(5,1) + t667 * mrSges(5,2);
t706 = qJD(2) + qJD(4);
t660 = t706 * mrSges(5,1) - t667 * mrSges(5,3);
t705 = qJDD(2) + qJDD(4);
t650 = -t666 * pkin(4) - t667 * pkin(9);
t704 = t706 ^ 2;
t601 = -t704 * pkin(4) + t705 * pkin(9) + t666 * t650 + t605;
t699 = t713 * g(1) - t716 * g(2);
t722 = -qJDD(1) * pkin(1) - t699;
t656 = -t695 * pkin(2) + qJDD(3) + t696 * t732 + (-qJ(3) * t707 - pkin(7)) * t717 + t722;
t624 = -t672 * pkin(3) - t682 * pkin(8) + t684 * t678 + t656;
t636 = t666 * qJD(4) + t711 * t672 + t714 * t673;
t603 = (-t666 * t706 - t636) * pkin(9) + (t667 * t706 - t635) * pkin(4) + t624;
t710 = sin(qJ(5));
t598 = t744 * t601 + t710 * t603;
t658 = t744 * t667 + t710 * t706;
t613 = t658 * qJD(5) + t710 * t636 - t744 * t705;
t634 = qJDD(5) - t635;
t662 = qJD(5) - t666;
t642 = t662 * mrSges(6,1) - t658 * mrSges(6,3);
t657 = t710 * t667 - t744 * t706;
t637 = t657 * pkin(5) - t658 * qJ(6);
t661 = t662 ^ 2;
t594 = -t661 * pkin(5) + t634 * qJ(6) + 0.2e1 * qJD(6) * t662 - t657 * t637 + t598;
t643 = -t662 * mrSges(7,1) + t658 * mrSges(7,2);
t729 = m(7) * t594 + t634 * mrSges(7,3) + t662 * t643;
t638 = t657 * mrSges(7,1) - t658 * mrSges(7,3);
t733 = -t657 * mrSges(6,1) - t658 * mrSges(6,2) - t638;
t589 = m(6) * t598 - t634 * mrSges(6,2) + t742 * t613 - t662 * t642 + t733 * t657 + t729;
t597 = -t710 * t601 + t744 * t603;
t614 = -t657 * qJD(5) + t744 * t636 + t710 * t705;
t641 = -t662 * mrSges(6,2) - t657 * mrSges(6,3);
t595 = -t634 * pkin(5) - t661 * qJ(6) + t658 * t637 + qJDD(6) - t597;
t640 = -t657 * mrSges(7,2) + t662 * mrSges(7,3);
t723 = -m(7) * t595 + t634 * mrSges(7,1) + t662 * t640;
t591 = m(6) * t597 + t634 * mrSges(6,1) + t742 * t614 + t662 * t641 + t733 * t658 + t723;
t724 = t744 * t589 - t710 * t591;
t581 = m(5) * t605 - t705 * mrSges(5,2) + t635 * mrSges(5,3) + t666 * t649 - t706 * t660 + t724;
t604 = t714 * t607 - t711 * t610;
t659 = -t706 * mrSges(5,2) + t666 * mrSges(5,3);
t600 = -t705 * pkin(4) - t704 * pkin(9) + t667 * t650 - t604;
t596 = -0.2e1 * qJD(6) * t658 + (t657 * t662 - t614) * qJ(6) + (t658 * t662 + t613) * pkin(5) + t600;
t592 = m(7) * t596 + t613 * mrSges(7,1) - t614 * mrSges(7,3) + t657 * t640 - t658 * t643;
t719 = -m(6) * t600 - t613 * mrSges(6,1) - t614 * mrSges(6,2) - t657 * t641 - t658 * t642 - t592;
t586 = m(5) * t604 + t705 * mrSges(5,1) - t636 * mrSges(5,3) - t667 * t649 + t706 * t659 + t719;
t576 = t711 * t581 + t714 * t586;
t670 = -t683 * mrSges(4,1) + t684 * mrSges(4,2);
t676 = -qJD(2) * mrSges(4,2) + t683 * mrSges(4,3);
t574 = m(4) * t627 + qJDD(2) * mrSges(4,1) - t673 * mrSges(4,3) + qJD(2) * t676 - t684 * t670 + t576;
t677 = qJD(2) * mrSges(4,1) - t684 * mrSges(4,3);
t725 = t714 * t581 - t711 * t586;
t575 = m(4) * t628 - qJDD(2) * mrSges(4,2) + t672 * mrSges(4,3) - qJD(2) * t677 + t683 * t670 + t725;
t568 = t709 * t574 + t708 * t575;
t674 = -t715 * g(3) - t738;
t693 = (-mrSges(3,1) * t715 + mrSges(3,2) * t712) * qJD(1);
t731 = qJD(1) * t715;
t698 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t731;
t566 = m(3) * t674 + qJDD(2) * mrSges(3,1) - t694 * mrSges(3,3) + qJD(2) * t698 - t693 * t732 + t568;
t697 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t732;
t726 = -t708 * t574 + t709 * t575;
t567 = m(3) * t675 - qJDD(2) * mrSges(3,2) + t695 * mrSges(3,3) - qJD(2) * t697 + t693 * t731 + t726;
t727 = -t712 * t566 + t715 * t567;
t561 = m(2) * t700 - t717 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t727;
t688 = -t717 * pkin(7) + t722;
t584 = t710 * t589 + t744 * t591;
t721 = m(5) * t624 - t635 * mrSges(5,1) + t636 * mrSges(5,2) - t666 * t659 + t667 * t660 + t584;
t720 = m(4) * t656 - t672 * mrSges(4,1) + t673 * mrSges(4,2) - t683 * t676 + t684 * t677 + t721;
t718 = -m(3) * t688 + t695 * mrSges(3,1) - t694 * mrSges(3,2) - t697 * t732 + t698 * t731 - t720;
t578 = m(2) * t699 + qJDD(1) * mrSges(2,1) - t717 * mrSges(2,2) + t718;
t737 = t713 * t561 + t716 * t578;
t562 = t715 * t566 + t712 * t567;
t736 = t657 * t746 - t658 * t741 - t662 * t739;
t735 = t657 * t739 + t658 * t747 + t662 * t745;
t734 = -t741 * t657 + t658 * t748 - t747 * t662;
t728 = t716 * t561 - t713 * t578;
t687 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t712 + Ifges(3,4) * t715) * qJD(1);
t686 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t712 + Ifges(3,2) * t715) * qJD(1);
t685 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t712 + Ifges(3,6) * t715) * qJD(1);
t665 = Ifges(4,1) * t684 + Ifges(4,4) * t683 + Ifges(4,5) * qJD(2);
t664 = Ifges(4,4) * t684 + Ifges(4,2) * t683 + Ifges(4,6) * qJD(2);
t663 = Ifges(4,5) * t684 + Ifges(4,6) * t683 + Ifges(4,3) * qJD(2);
t646 = Ifges(5,1) * t667 + Ifges(5,4) * t666 + Ifges(5,5) * t706;
t645 = Ifges(5,4) * t667 + Ifges(5,2) * t666 + Ifges(5,6) * t706;
t644 = Ifges(5,5) * t667 + Ifges(5,6) * t666 + Ifges(5,3) * t706;
t583 = mrSges(6,2) * t600 + mrSges(7,2) * t595 - mrSges(6,3) * t597 - mrSges(7,3) * t596 - qJ(6) * t592 - t741 * t613 + t614 * t748 - t634 * t747 + t735 * t657 + t736 * t662;
t582 = -mrSges(6,1) * t600 - mrSges(7,1) * t596 + mrSges(7,2) * t594 + mrSges(6,3) * t598 - pkin(5) * t592 - t613 * t746 + t741 * t614 + t739 * t634 + t735 * t658 + t734 * t662;
t570 = Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * t705 - t667 * t644 + t706 * t646 - mrSges(5,1) * t624 + mrSges(5,3) * t605 - mrSges(6,1) * t597 + mrSges(6,2) * t598 + mrSges(7,1) * t595 - mrSges(7,3) * t594 - pkin(5) * t723 - qJ(6) * t729 - pkin(4) * t584 + (pkin(5) * t638 + t736) * t658 + (qJ(6) * t638 - t734) * t657 + t745 * t634 + (pkin(5) * mrSges(7,2) + t747) * t614 + (qJ(6) * mrSges(7,2) + t739) * t613;
t569 = mrSges(5,2) * t624 - mrSges(5,3) * t604 + Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * t705 - pkin(9) * t584 - t710 * t582 + t744 * t583 + t666 * t644 - t706 * t645;
t558 = mrSges(4,2) * t656 - mrSges(4,3) * t627 + Ifges(4,1) * t673 + Ifges(4,4) * t672 + Ifges(4,5) * qJDD(2) - pkin(8) * t576 - qJD(2) * t664 + t714 * t569 - t711 * t570 + t683 * t663;
t557 = -mrSges(4,1) * t656 + mrSges(4,3) * t628 + Ifges(4,4) * t673 + Ifges(4,2) * t672 + Ifges(4,6) * qJDD(2) - pkin(3) * t721 + pkin(8) * t725 + qJD(2) * t665 + t711 * t569 + t714 * t570 - t684 * t663;
t556 = -t744 * t582 - pkin(1) * t562 + mrSges(2,1) * g(3) - pkin(4) * t719 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t712 * t686 + t715 * t687) * qJD(1) - pkin(9) * t724 - pkin(2) * t568 + t717 * Ifges(2,5) - t710 * t583 - Ifges(5,3) * t705 - Ifges(3,5) * t694 - Ifges(3,6) * t695 + mrSges(2,3) * t700 + t683 * t665 - t684 * t664 - Ifges(4,6) * t672 - Ifges(4,5) * t673 - mrSges(3,1) * t674 + mrSges(3,2) * t675 - t667 * t645 + t666 * t646 - Ifges(5,6) * t635 - Ifges(5,5) * t636 - mrSges(4,1) * t627 + mrSges(4,2) * t628 - mrSges(5,1) * t604 + mrSges(5,2) * t605 + Ifges(2,6) * qJDD(1) - pkin(3) * t576;
t555 = mrSges(3,2) * t688 - mrSges(3,3) * t674 + Ifges(3,1) * t694 + Ifges(3,4) * t695 + Ifges(3,5) * qJDD(2) - qJ(3) * t568 - qJD(2) * t686 - t708 * t557 + t709 * t558 + t685 * t731;
t554 = -mrSges(3,1) * t688 + mrSges(3,3) * t675 + Ifges(3,4) * t694 + Ifges(3,2) * t695 + Ifges(3,6) * qJDD(2) - pkin(2) * t720 + qJ(3) * t726 + qJD(2) * t687 + t709 * t557 + t708 * t558 - t685 * t732;
t553 = -mrSges(2,2) * g(3) - mrSges(2,3) * t699 + Ifges(2,5) * qJDD(1) - t717 * Ifges(2,6) - pkin(7) * t562 - t712 * t554 + t715 * t555;
t1 = [-m(1) * g(1) + t728; -m(1) * g(2) + t737; (-m(1) - m(2)) * g(3) + t562; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t737 + t716 * t553 - t713 * t556; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t728 + t713 * t553 + t716 * t556; -mrSges(1,1) * g(2) + mrSges(2,1) * t699 + mrSges(1,2) * g(1) - mrSges(2,2) * t700 + Ifges(2,3) * qJDD(1) + pkin(1) * t718 + pkin(7) * t727 + t715 * t554 + t712 * t555;];
tauB  = t1;
