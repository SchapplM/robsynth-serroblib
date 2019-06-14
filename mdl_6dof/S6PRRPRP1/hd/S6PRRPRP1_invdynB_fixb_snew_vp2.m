% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:36:58
% EndTime: 2019-05-05 03:37:11
% DurationCPUTime: 9.45s
% Computational Cost: add. (144885->319), mult. (310421->400), div. (0->0), fcn. (218607->12), ass. (0->133)
t715 = -2 * qJD(4);
t714 = Ifges(6,1) + Ifges(7,1);
t709 = Ifges(6,4) + Ifges(7,4);
t708 = Ifges(6,5) + Ifges(7,5);
t713 = Ifges(6,2) + Ifges(7,2);
t712 = Ifges(6,6) + Ifges(7,6);
t711 = Ifges(6,3) + Ifges(7,3);
t666 = sin(pkin(10));
t669 = cos(pkin(10));
t657 = g(1) * t666 - g(2) * t669;
t658 = -g(1) * t669 - g(2) * t666;
t664 = -g(3) + qJDD(1);
t676 = cos(qJ(2));
t670 = cos(pkin(6));
t673 = sin(qJ(2));
t704 = t670 * t673;
t667 = sin(pkin(6));
t705 = t667 * t673;
t623 = t657 * t704 + t676 * t658 + t664 * t705;
t678 = qJD(2) ^ 2;
t618 = -pkin(2) * t678 + qJDD(2) * pkin(8) + t623;
t639 = -t657 * t667 + t664 * t670;
t672 = sin(qJ(3));
t675 = cos(qJ(3));
t592 = -t672 * t618 + t675 * t639;
t695 = qJD(2) * qJD(3);
t692 = t675 * t695;
t655 = qJDD(2) * t672 + t692;
t588 = (-t655 + t692) * qJ(4) + (t672 * t675 * t678 + qJDD(3)) * pkin(3) + t592;
t593 = t675 * t618 + t672 * t639;
t656 = qJDD(2) * t675 - t672 * t695;
t698 = qJD(2) * t672;
t659 = qJD(3) * pkin(3) - qJ(4) * t698;
t663 = t675 ^ 2;
t589 = -pkin(3) * t663 * t678 + qJ(4) * t656 - qJD(3) * t659 + t593;
t665 = sin(pkin(11));
t668 = cos(pkin(11));
t644 = (t665 * t675 + t668 * t672) * qJD(2);
t580 = t588 * t668 - t665 * t589 + t644 * t715;
t643 = (t665 * t672 - t668 * t675) * qJD(2);
t622 = -t673 * t658 + (t657 * t670 + t664 * t667) * t676;
t710 = -mrSges(6,2) - mrSges(7,2);
t682 = -qJDD(2) * pkin(2) - t622;
t617 = -t678 * pkin(8) + t682;
t660 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t698;
t697 = qJD(2) * t675;
t661 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t697;
t581 = t665 * t588 + t668 * t589 + t643 * t715;
t626 = pkin(4) * t643 - pkin(9) * t644;
t677 = qJD(3) ^ 2;
t579 = -pkin(4) * t677 + qJDD(3) * pkin(9) - t626 * t643 + t581;
t590 = -t656 * pkin(3) + qJDD(4) + t659 * t698 + (-qJ(4) * t663 - pkin(8)) * t678 + t682;
t630 = -t655 * t665 + t656 * t668;
t631 = t655 * t668 + t656 * t665;
t584 = (qJD(3) * t643 - t631) * pkin(9) + (qJD(3) * t644 - t630) * pkin(4) + t590;
t671 = sin(qJ(5));
t674 = cos(qJ(5));
t574 = -t671 * t579 + t674 * t584;
t633 = qJD(3) * t674 - t644 * t671;
t606 = qJD(5) * t633 + qJDD(3) * t671 + t631 * t674;
t634 = qJD(3) * t671 + t644 * t674;
t608 = -mrSges(7,1) * t633 + mrSges(7,2) * t634;
t609 = -mrSges(6,1) * t633 + mrSges(6,2) * t634;
t642 = qJD(5) + t643;
t613 = -mrSges(6,2) * t642 + mrSges(6,3) * t633;
t629 = qJDD(5) - t630;
t571 = -0.2e1 * qJD(6) * t634 + (t633 * t642 - t606) * qJ(6) + (t633 * t634 + t629) * pkin(5) + t574;
t612 = -mrSges(7,2) * t642 + mrSges(7,3) * t633;
t694 = m(7) * t571 + t629 * mrSges(7,1) + t642 * t612;
t563 = m(6) * t574 + t629 * mrSges(6,1) + t642 * t613 + (-t608 - t609) * t634 + (-mrSges(6,3) - mrSges(7,3)) * t606 + t694;
t575 = t674 * t579 + t671 * t584;
t605 = -qJD(5) * t634 + qJDD(3) * t674 - t631 * t671;
t614 = pkin(5) * t642 - qJ(6) * t634;
t632 = t633 ^ 2;
t573 = -pkin(5) * t632 + qJ(6) * t605 + 0.2e1 * qJD(6) * t633 - t614 * t642 + t575;
t693 = m(7) * t573 + t605 * mrSges(7,3) + t633 * t608;
t615 = mrSges(7,1) * t642 - mrSges(7,3) * t634;
t699 = -mrSges(6,1) * t642 + mrSges(6,3) * t634 - t615;
t568 = m(6) * t575 + t605 * mrSges(6,3) + t633 * t609 + t629 * t710 + t699 * t642 + t693;
t561 = t674 * t563 + t671 * t568;
t637 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t643;
t638 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t644;
t681 = m(5) * t590 - t630 * mrSges(5,1) + mrSges(5,2) * t631 + t643 * t637 + t638 * t644 + t561;
t679 = -m(4) * t617 + t656 * mrSges(4,1) - mrSges(4,2) * t655 - t660 * t698 + t661 * t697 - t681;
t556 = m(3) * t622 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t678 + t679;
t706 = t556 * t676;
t625 = mrSges(5,1) * t643 + mrSges(5,2) * t644;
t688 = -t563 * t671 + t674 * t568;
t559 = m(5) * t581 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t630 - qJD(3) * t638 - t625 * t643 + t688;
t578 = -qJDD(3) * pkin(4) - pkin(9) * t677 + t644 * t626 - t580;
t576 = -pkin(5) * t605 - qJ(6) * t632 + t614 * t634 + qJDD(6) + t578;
t686 = m(7) * t576 - t605 * mrSges(7,1) - t633 * t612;
t680 = -m(6) * t578 + t605 * mrSges(6,1) + t710 * t606 + t633 * t613 + t699 * t634 - t686;
t565 = m(5) * t580 + qJDD(3) * mrSges(5,1) - t631 * mrSges(5,3) + qJD(3) * t637 - t644 * t625 + t680;
t552 = t665 * t559 + t668 * t565;
t654 = (-mrSges(4,1) * t675 + mrSges(4,2) * t672) * qJD(2);
t550 = m(4) * t592 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t655 + qJD(3) * t661 - t654 * t698 + t552;
t689 = t668 * t559 - t565 * t665;
t551 = m(4) * t593 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t656 - qJD(3) * t660 + t654 * t697 + t689;
t690 = -t550 * t672 + t675 * t551;
t542 = m(3) * t623 - mrSges(3,1) * t678 - qJDD(2) * mrSges(3,2) + t690;
t545 = t675 * t550 + t672 * t551;
t544 = m(3) * t639 + t545;
t532 = t542 * t704 - t544 * t667 + t670 * t706;
t530 = m(2) * t657 + t532;
t537 = t676 * t542 - t556 * t673;
t536 = m(2) * t658 + t537;
t703 = t669 * t530 + t666 * t536;
t702 = t633 * t712 + t634 * t708 + t642 * t711;
t701 = -t633 * t713 - t634 * t709 - t642 * t712;
t700 = t633 * t709 + t634 * t714 + t642 * t708;
t531 = t542 * t705 + t670 * t544 + t667 * t706;
t691 = -t530 * t666 + t669 * t536;
t553 = -mrSges(6,1) * t578 + mrSges(6,3) * t575 - mrSges(7,1) * t576 + mrSges(7,3) * t573 - pkin(5) * t686 + qJ(6) * t693 + (-qJ(6) * t615 + t700) * t642 + (-pkin(5) * t615 - t702) * t634 + (-mrSges(7,2) * qJ(6) + t712) * t629 + (-mrSges(7,2) * pkin(5) + t709) * t606 + t713 * t605;
t569 = -t606 * mrSges(7,3) - t634 * t608 + t694;
t560 = mrSges(6,2) * t578 + mrSges(7,2) * t576 - mrSges(6,3) * t574 - mrSges(7,3) * t571 - qJ(6) * t569 + t709 * t605 + t606 * t714 + t708 * t629 + t702 * t633 + t701 * t642;
t619 = Ifges(5,5) * t644 - Ifges(5,6) * t643 + Ifges(5,3) * qJD(3);
t620 = Ifges(5,4) * t644 - Ifges(5,2) * t643 + Ifges(5,6) * qJD(3);
t538 = mrSges(5,2) * t590 - mrSges(5,3) * t580 + Ifges(5,1) * t631 + Ifges(5,4) * t630 + Ifges(5,5) * qJDD(3) - pkin(9) * t561 - qJD(3) * t620 - t553 * t671 + t560 * t674 - t619 * t643;
t621 = Ifges(5,1) * t644 - Ifges(5,4) * t643 + Ifges(5,5) * qJD(3);
t546 = -mrSges(5,1) * t590 - mrSges(6,1) * t574 - mrSges(7,1) * t571 + mrSges(6,2) * t575 + mrSges(7,2) * t573 + mrSges(5,3) * t581 + Ifges(5,4) * t631 + Ifges(5,2) * t630 + Ifges(5,6) * qJDD(3) - pkin(4) * t561 - pkin(5) * t569 + qJD(3) * t621 - t644 * t619 + t701 * t634 + t700 * t633 - t711 * t629 - t708 * t606 - t712 * t605;
t646 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t672 + Ifges(4,6) * t675) * qJD(2);
t648 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t672 + Ifges(4,4) * t675) * qJD(2);
t528 = -mrSges(4,1) * t617 + mrSges(4,3) * t593 + Ifges(4,4) * t655 + Ifges(4,2) * t656 + Ifges(4,6) * qJDD(3) - pkin(3) * t681 + qJ(4) * t689 + qJD(3) * t648 + t665 * t538 + t668 * t546 - t646 * t698;
t647 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t672 + Ifges(4,2) * t675) * qJD(2);
t533 = mrSges(4,2) * t617 - mrSges(4,3) * t592 + Ifges(4,1) * t655 + Ifges(4,4) * t656 + Ifges(4,5) * qJDD(3) - qJ(4) * t552 - qJD(3) * t647 + t538 * t668 - t546 * t665 + t646 * t697;
t526 = mrSges(3,2) * t639 - mrSges(3,3) * t622 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t678 - pkin(8) * t545 - t528 * t672 + t533 * t675;
t527 = -t644 * t620 - t643 * t621 + t678 * Ifges(3,5) - pkin(2) * t545 + mrSges(3,3) * t623 - mrSges(3,1) * t639 + Ifges(3,6) * qJDD(2) - pkin(3) * t552 - Ifges(4,5) * t655 - Ifges(4,6) * t656 - mrSges(4,1) * t592 + mrSges(4,2) * t593 - t671 * t560 - t674 * t553 - pkin(4) * t680 - pkin(9) * t688 - Ifges(5,5) * t631 - Ifges(5,6) * t630 - mrSges(5,1) * t580 + mrSges(5,2) * t581 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t647 * t672 + t648 * t675) * qJD(2);
t683 = pkin(7) * t537 + t526 * t673 + t527 * t676;
t525 = mrSges(3,1) * t622 - mrSges(3,2) * t623 + Ifges(3,3) * qJDD(2) + pkin(2) * t679 + pkin(8) * t690 + t675 * t528 + t672 * t533;
t524 = mrSges(2,2) * t664 - mrSges(2,3) * t657 + t676 * t526 - t673 * t527 + (-t531 * t667 - t532 * t670) * pkin(7);
t523 = -mrSges(2,1) * t664 + mrSges(2,3) * t658 - pkin(1) * t531 - t667 * t525 + t670 * t683;
t1 = [-m(1) * g(1) + t691; -m(1) * g(2) + t703; -m(1) * g(3) + m(2) * t664 + t531; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t703 - t666 * t523 + t669 * t524; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t691 + t669 * t523 + t666 * t524; -mrSges(1,1) * g(2) + mrSges(2,1) * t657 + mrSges(1,2) * g(1) - mrSges(2,2) * t658 + pkin(1) * t532 + t670 * t525 + t667 * t683;];
tauB  = t1;
