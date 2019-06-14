% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:41:20
% EndTime: 2019-05-05 05:42:03
% DurationCPUTime: 42.16s
% Computational Cost: add. (724755->352), mult. (1572881->468), div. (0->0), fcn. (1274518->16), ass. (0->154)
t692 = sin(pkin(12));
t696 = cos(pkin(12));
t682 = g(1) * t692 - g(2) * t696;
t683 = -g(1) * t696 - g(2) * t692;
t690 = -g(3) + qJDD(1);
t702 = sin(qJ(2));
t698 = cos(pkin(6));
t706 = cos(qJ(2));
t724 = t698 * t706;
t694 = sin(pkin(6));
t727 = t694 * t706;
t650 = t682 * t724 - t683 * t702 + t690 * t727;
t693 = sin(pkin(7));
t707 = qJD(2) ^ 2;
t648 = pkin(9) * t693 * t707 + qJDD(2) * pkin(2) + t650;
t725 = t698 * t702;
t728 = t694 * t702;
t651 = t682 * t725 + t706 * t683 + t690 * t728;
t720 = qJDD(2) * t693;
t649 = -pkin(2) * t707 + pkin(9) * t720 + t651;
t669 = -t682 * t694 + t690 * t698;
t697 = cos(pkin(7));
t701 = sin(qJ(3));
t705 = cos(qJ(3));
t617 = -t701 * t649 + (t648 * t697 + t669 * t693) * t705;
t689 = qJD(2) * t697 + qJD(3);
t721 = qJD(2) * t705;
t718 = t693 * t721;
t672 = -mrSges(4,2) * t689 + mrSges(4,3) * t718;
t722 = qJD(2) * t693;
t674 = (-mrSges(4,1) * t705 + mrSges(4,2) * t701) * t722;
t675 = (qJD(3) * t721 + qJDD(2) * t701) * t693;
t688 = qJDD(2) * t697 + qJDD(3);
t673 = (-pkin(3) * t705 - qJ(4) * t701) * t722;
t687 = t689 ^ 2;
t719 = t701 * t722;
t611 = -t688 * pkin(3) - t687 * qJ(4) + t673 * t719 + qJDD(4) - t617;
t691 = sin(pkin(13));
t695 = cos(pkin(13));
t667 = t689 * t695 - t691 * t719;
t653 = mrSges(5,2) * t718 + mrSges(5,3) * t667;
t668 = t689 * t691 + t695 * t719;
t654 = -mrSges(5,1) * t718 - mrSges(5,3) * t668;
t655 = -t675 * t691 + t688 * t695;
t656 = t675 * t695 + t688 * t691;
t726 = t697 * t701;
t729 = t693 * t701;
t618 = t648 * t726 + t705 * t649 + t669 * t729;
t612 = -pkin(3) * t687 + qJ(4) * t688 + t673 * t718 + t618;
t665 = t697 * t669;
t676 = -qJD(3) * t719 + t705 * t720;
t615 = -t676 * pkin(3) - t675 * qJ(4) + t665 + (-t648 + (pkin(3) * t701 - qJ(4) * t705) * t689 * qJD(2)) * t693;
t600 = -0.2e1 * qJD(4) * t668 - t691 * t612 + t695 * t615;
t597 = (-t667 * t718 - t656) * pkin(10) + (t667 * t668 - t676) * pkin(4) + t600;
t601 = 0.2e1 * qJD(4) * t667 + t695 * t612 + t691 * t615;
t657 = -pkin(4) * t718 - pkin(10) * t668;
t666 = t667 ^ 2;
t599 = -pkin(4) * t666 + pkin(10) * t655 + t657 * t718 + t601;
t700 = sin(qJ(5));
t704 = cos(qJ(5));
t594 = t700 * t597 + t704 * t599;
t643 = t667 * t704 - t668 * t700;
t644 = t667 * t700 + t668 * t704;
t630 = -pkin(5) * t643 - pkin(11) * t644;
t670 = qJDD(5) - t676;
t681 = qJD(5) - t718;
t680 = t681 ^ 2;
t592 = -pkin(5) * t680 + pkin(11) * t670 + t630 * t643 + t594;
t602 = -t655 * pkin(4) - t666 * pkin(10) + t668 * t657 + t611;
t621 = -qJD(5) * t644 + t655 * t704 - t656 * t700;
t622 = qJD(5) * t643 + t655 * t700 + t656 * t704;
t595 = (-t643 * t681 - t622) * pkin(11) + (t644 * t681 - t621) * pkin(5) + t602;
t699 = sin(qJ(6));
t703 = cos(qJ(6));
t589 = -t592 * t699 + t595 * t703;
t632 = -t644 * t699 + t681 * t703;
t605 = qJD(6) * t632 + t622 * t703 + t670 * t699;
t633 = t644 * t703 + t681 * t699;
t616 = -mrSges(7,1) * t632 + mrSges(7,2) * t633;
t620 = qJDD(6) - t621;
t642 = qJD(6) - t643;
t623 = -mrSges(7,2) * t642 + mrSges(7,3) * t632;
t587 = m(7) * t589 + mrSges(7,1) * t620 - mrSges(7,3) * t605 - t616 * t633 + t623 * t642;
t590 = t592 * t703 + t595 * t699;
t604 = -qJD(6) * t633 - t622 * t699 + t670 * t703;
t624 = mrSges(7,1) * t642 - mrSges(7,3) * t633;
t588 = m(7) * t590 - mrSges(7,2) * t620 + mrSges(7,3) * t604 + t616 * t632 - t624 * t642;
t579 = t703 * t587 + t699 * t588;
t634 = -mrSges(6,2) * t681 + mrSges(6,3) * t643;
t635 = mrSges(6,1) * t681 - mrSges(6,3) * t644;
t709 = m(6) * t602 - t621 * mrSges(6,1) + mrSges(6,2) * t622 - t643 * t634 + t635 * t644 + t579;
t708 = -m(5) * t611 + t655 * mrSges(5,1) - mrSges(5,2) * t656 + t667 * t653 - t654 * t668 - t709;
t575 = m(4) * t617 + mrSges(4,1) * t688 - mrSges(4,3) * t675 + t672 * t689 - t674 * t719 + t708;
t730 = t575 * t705;
t671 = mrSges(4,1) * t689 - mrSges(4,3) * t719;
t629 = -mrSges(6,1) * t643 + mrSges(6,2) * t644;
t714 = -t587 * t699 + t703 * t588;
t578 = m(6) * t594 - mrSges(6,2) * t670 + mrSges(6,3) * t621 + t629 * t643 - t635 * t681 + t714;
t593 = t597 * t704 - t599 * t700;
t591 = -pkin(5) * t670 - pkin(11) * t680 + t630 * t644 - t593;
t710 = -m(7) * t591 + t604 * mrSges(7,1) - mrSges(7,2) * t605 + t632 * t623 - t624 * t633;
t583 = m(6) * t593 + mrSges(6,1) * t670 - mrSges(6,3) * t622 - t629 * t644 + t634 * t681 + t710;
t572 = t700 * t578 + t704 * t583;
t647 = -mrSges(5,1) * t667 + mrSges(5,2) * t668;
t570 = m(5) * t600 - mrSges(5,1) * t676 - mrSges(5,3) * t656 - t647 * t668 - t653 * t718 + t572;
t715 = t704 * t578 - t583 * t700;
t571 = m(5) * t601 + mrSges(5,2) * t676 + mrSges(5,3) * t655 + t647 * t667 + t654 * t718 + t715;
t716 = -t570 * t691 + t695 * t571;
t561 = m(4) * t618 - mrSges(4,2) * t688 + mrSges(4,3) * t676 - t671 * t689 + t674 * t718 + t716;
t564 = t695 * t570 + t691 * t571;
t631 = -t693 * t648 + t665;
t563 = m(4) * t631 - t676 * mrSges(4,1) + t675 * mrSges(4,2) + (t671 * t701 - t672 * t705) * t722 + t564;
t550 = t561 * t726 - t563 * t693 + t697 * t730;
t546 = m(3) * t650 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t707 + t550;
t549 = t561 * t729 + t697 * t563 + t693 * t730;
t548 = m(3) * t669 + t549;
t557 = t705 * t561 - t575 * t701;
t556 = m(3) * t651 - mrSges(3,1) * t707 - qJDD(2) * mrSges(3,2) + t557;
t536 = t546 * t724 - t548 * t694 + t556 * t725;
t534 = m(2) * t682 + t536;
t542 = -t546 * t702 + t706 * t556;
t541 = m(2) * t683 + t542;
t723 = t696 * t534 + t692 * t541;
t535 = t546 * t727 + t698 * t548 + t556 * t728;
t717 = -t534 * t692 + t696 * t541;
t606 = Ifges(7,5) * t633 + Ifges(7,6) * t632 + Ifges(7,3) * t642;
t608 = Ifges(7,1) * t633 + Ifges(7,4) * t632 + Ifges(7,5) * t642;
t580 = -mrSges(7,1) * t591 + mrSges(7,3) * t590 + Ifges(7,4) * t605 + Ifges(7,2) * t604 + Ifges(7,6) * t620 - t606 * t633 + t608 * t642;
t607 = Ifges(7,4) * t633 + Ifges(7,2) * t632 + Ifges(7,6) * t642;
t581 = mrSges(7,2) * t591 - mrSges(7,3) * t589 + Ifges(7,1) * t605 + Ifges(7,4) * t604 + Ifges(7,5) * t620 + t606 * t632 - t607 * t642;
t625 = Ifges(6,5) * t644 + Ifges(6,6) * t643 + Ifges(6,3) * t681;
t626 = Ifges(6,4) * t644 + Ifges(6,2) * t643 + Ifges(6,6) * t681;
t565 = mrSges(6,2) * t602 - mrSges(6,3) * t593 + Ifges(6,1) * t622 + Ifges(6,4) * t621 + Ifges(6,5) * t670 - pkin(11) * t579 - t580 * t699 + t581 * t703 + t625 * t643 - t626 * t681;
t627 = Ifges(6,1) * t644 + Ifges(6,4) * t643 + Ifges(6,5) * t681;
t566 = -mrSges(6,1) * t602 - mrSges(7,1) * t589 + mrSges(7,2) * t590 + mrSges(6,3) * t594 + Ifges(6,4) * t622 - Ifges(7,5) * t605 + Ifges(6,2) * t621 + Ifges(6,6) * t670 - Ifges(7,6) * t604 - Ifges(7,3) * t620 - pkin(5) * t579 - t607 * t633 + t608 * t632 - t625 * t644 + t627 * t681;
t636 = Ifges(5,5) * t668 + Ifges(5,6) * t667 - Ifges(5,3) * t718;
t638 = Ifges(5,1) * t668 + Ifges(5,4) * t667 - Ifges(5,5) * t718;
t551 = -mrSges(5,1) * t611 + mrSges(5,3) * t601 + Ifges(5,4) * t656 + Ifges(5,2) * t655 - Ifges(5,6) * t676 - pkin(4) * t709 + pkin(10) * t715 + t700 * t565 + t704 * t566 - t668 * t636 - t638 * t718;
t637 = Ifges(5,4) * t668 + Ifges(5,2) * t667 - Ifges(5,6) * t718;
t552 = mrSges(5,2) * t611 - mrSges(5,3) * t600 + Ifges(5,1) * t656 + Ifges(5,4) * t655 - Ifges(5,5) * t676 - pkin(10) * t572 + t565 * t704 - t566 * t700 + t636 * t667 + t637 * t718;
t660 = Ifges(4,6) * t689 + (Ifges(4,4) * t701 + Ifges(4,2) * t705) * t722;
t661 = Ifges(4,5) * t689 + (Ifges(4,1) * t701 + Ifges(4,4) * t705) * t722;
t537 = Ifges(4,5) * t675 + Ifges(4,6) * t676 + Ifges(4,3) * t688 + mrSges(4,1) * t617 - mrSges(4,2) * t618 + t691 * t552 + t695 * t551 + pkin(3) * t708 + qJ(4) * t716 + (t660 * t701 - t661 * t705) * t722;
t659 = Ifges(4,3) * t689 + (Ifges(4,5) * t701 + Ifges(4,6) * t705) * t722;
t538 = mrSges(4,2) * t631 - mrSges(4,3) * t617 + Ifges(4,1) * t675 + Ifges(4,4) * t676 + Ifges(4,5) * t688 - qJ(4) * t564 - t551 * t691 + t552 * t695 + t659 * t718 - t660 * t689;
t543 = -Ifges(6,3) * t670 + Ifges(4,4) * t675 - t659 * t719 - mrSges(5,1) * t600 + mrSges(5,2) * t601 + Ifges(4,6) * t688 + mrSges(4,3) * t618 - Ifges(6,6) * t621 - Ifges(6,5) * t622 - pkin(11) * t714 - Ifges(5,6) * t655 - Ifges(5,5) * t656 - pkin(3) * t564 - pkin(5) * t710 + t689 * t661 - pkin(4) * t572 - mrSges(6,1) * t593 + mrSges(6,2) * t594 + t643 * t627 - t644 * t626 - mrSges(4,1) * t631 - t699 * t581 - t703 * t580 + (Ifges(4,2) + Ifges(5,3)) * t676 + t667 * t638 - t668 * t637;
t711 = pkin(9) * t557 + t538 * t701 + t543 * t705;
t531 = -mrSges(3,1) * t669 + mrSges(3,3) * t651 + t707 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t549 - t693 * t537 + t697 * t711;
t532 = mrSges(3,2) * t669 - mrSges(3,3) * t650 + Ifges(3,5) * qJDD(2) - t707 * Ifges(3,6) + t705 * t538 - t701 * t543 + (-t549 * t693 - t550 * t697) * pkin(9);
t712 = pkin(8) * t542 + t531 * t706 + t532 * t702;
t530 = mrSges(3,1) * t650 - mrSges(3,2) * t651 + Ifges(3,3) * qJDD(2) + pkin(2) * t550 + t697 * t537 + t693 * t711;
t529 = mrSges(2,2) * t690 - mrSges(2,3) * t682 - t702 * t531 + t706 * t532 + (-t535 * t694 - t536 * t698) * pkin(8);
t528 = -mrSges(2,1) * t690 + mrSges(2,3) * t683 - pkin(1) * t535 - t694 * t530 + t698 * t712;
t1 = [-m(1) * g(1) + t717; -m(1) * g(2) + t723; -m(1) * g(3) + m(2) * t690 + t535; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t723 - t692 * t528 + t696 * t529; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t717 + t696 * t528 + t692 * t529; -mrSges(1,1) * g(2) + mrSges(2,1) * t682 + mrSges(1,2) * g(1) - mrSges(2,2) * t683 + pkin(1) * t536 + t698 * t530 + t694 * t712;];
tauB  = t1;
