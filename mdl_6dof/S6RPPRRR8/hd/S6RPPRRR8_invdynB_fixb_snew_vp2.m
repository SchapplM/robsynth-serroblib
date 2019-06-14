% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:14:17
% EndTime: 2019-05-05 16:14:24
% DurationCPUTime: 6.63s
% Computational Cost: add. (93656->315), mult. (212325->385), div. (0->0), fcn. (151226->10), ass. (0->133)
t663 = sin(qJ(1));
t667 = cos(qJ(1));
t638 = t663 * g(1) - t667 * g(2);
t669 = qJD(1) ^ 2;
t677 = -t669 * qJ(2) + qJDD(2) - t638;
t699 = -pkin(1) - qJ(3);
t705 = -(2 * qJD(1) * qJD(3)) + t699 * qJDD(1) + t677;
t658 = sin(pkin(10));
t652 = t658 ^ 2;
t659 = cos(pkin(10));
t697 = t659 ^ 2 + t652;
t689 = t697 * mrSges(4,3);
t639 = -t667 * g(1) - t663 * g(2);
t704 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t639;
t703 = pkin(3) * t669;
t702 = mrSges(2,1) - mrSges(3,2);
t701 = -Ifges(3,4) + Ifges(2,5);
t700 = -Ifges(2,6) + Ifges(3,5);
t614 = t658 * g(3) + t659 * t705;
t599 = (-pkin(7) * qJDD(1) - t658 * t703) * t659 + t614;
t615 = -t659 * g(3) + t658 * t705;
t692 = qJDD(1) * t658;
t600 = -pkin(7) * t692 - t652 * t703 + t615;
t662 = sin(qJ(4));
t666 = cos(qJ(4));
t582 = t662 * t599 + t666 * t600;
t694 = t659 * qJD(1);
t695 = t658 * qJD(1);
t634 = -t662 * t694 - t666 * t695;
t680 = -t658 * t662 + t659 * t666;
t635 = t680 * qJD(1);
t611 = -t634 * mrSges(5,1) + t635 * mrSges(5,2);
t631 = t635 * qJD(4);
t691 = qJDD(1) * t659;
t617 = -t662 * t691 - t666 * t692 - t631;
t627 = qJD(4) * mrSges(5,1) - t635 * mrSges(5,3);
t616 = -t634 * pkin(4) - t635 * pkin(8);
t668 = qJD(4) ^ 2;
t572 = -t668 * pkin(4) + qJDD(4) * pkin(8) + t634 * t616 + t582;
t676 = qJDD(3) + t704;
t605 = pkin(3) * t692 + (-t697 * pkin(7) + t699) * t669 + t676;
t696 = t634 * qJD(4);
t618 = qJDD(1) * t680 + t696;
t575 = (-t618 - t696) * pkin(8) + (-t617 + t631) * pkin(4) + t605;
t661 = sin(qJ(5));
t665 = cos(qJ(5));
t565 = -t661 * t572 + t665 * t575;
t621 = t665 * qJD(4) - t661 * t635;
t591 = t621 * qJD(5) + t661 * qJDD(4) + t665 * t618;
t613 = qJDD(5) - t617;
t622 = t661 * qJD(4) + t665 * t635;
t632 = qJD(5) - t634;
t562 = (t621 * t632 - t591) * pkin(9) + (t621 * t622 + t613) * pkin(5) + t565;
t566 = t665 * t572 + t661 * t575;
t590 = -t622 * qJD(5) + t665 * qJDD(4) - t661 * t618;
t603 = t632 * pkin(5) - t622 * pkin(9);
t620 = t621 ^ 2;
t563 = -t620 * pkin(5) + t590 * pkin(9) - t632 * t603 + t566;
t660 = sin(qJ(6));
t664 = cos(qJ(6));
t560 = t664 * t562 - t660 * t563;
t592 = t664 * t621 - t660 * t622;
t569 = t592 * qJD(6) + t660 * t590 + t664 * t591;
t593 = t660 * t621 + t664 * t622;
t580 = -t592 * mrSges(7,1) + t593 * mrSges(7,2);
t630 = qJD(6) + t632;
t583 = -t630 * mrSges(7,2) + t592 * mrSges(7,3);
t610 = qJDD(6) + t613;
t558 = m(7) * t560 + t610 * mrSges(7,1) - t569 * mrSges(7,3) - t593 * t580 + t630 * t583;
t561 = t660 * t562 + t664 * t563;
t568 = -t593 * qJD(6) + t664 * t590 - t660 * t591;
t584 = t630 * mrSges(7,1) - t593 * mrSges(7,3);
t559 = m(7) * t561 - t610 * mrSges(7,2) + t568 * mrSges(7,3) + t592 * t580 - t630 * t584;
t550 = t664 * t558 + t660 * t559;
t595 = -t621 * mrSges(6,1) + t622 * mrSges(6,2);
t601 = -t632 * mrSges(6,2) + t621 * mrSges(6,3);
t548 = m(6) * t565 + t613 * mrSges(6,1) - t591 * mrSges(6,3) - t622 * t595 + t632 * t601 + t550;
t602 = t632 * mrSges(6,1) - t622 * mrSges(6,3);
t684 = -t660 * t558 + t664 * t559;
t549 = m(6) * t566 - t613 * mrSges(6,2) + t590 * mrSges(6,3) + t621 * t595 - t632 * t602 + t684;
t685 = -t661 * t548 + t665 * t549;
t543 = m(5) * t582 - qJDD(4) * mrSges(5,2) + t617 * mrSges(5,3) - qJD(4) * t627 + t634 * t611 + t685;
t581 = t666 * t599 - t662 * t600;
t626 = -qJD(4) * mrSges(5,2) + t634 * mrSges(5,3);
t571 = -qJDD(4) * pkin(4) - t668 * pkin(8) + t635 * t616 - t581;
t564 = -t590 * pkin(5) - t620 * pkin(9) + t622 * t603 + t571;
t674 = m(7) * t564 - t568 * mrSges(7,1) + t569 * mrSges(7,2) - t592 * t583 + t593 * t584;
t670 = -m(6) * t571 + t590 * mrSges(6,1) - t591 * mrSges(6,2) + t621 * t601 - t622 * t602 - t674;
t554 = m(5) * t581 + qJDD(4) * mrSges(5,1) - t618 * mrSges(5,3) + qJD(4) * t626 - t635 * t611 + t670;
t536 = t662 * t543 + t666 * t554;
t679 = -qJDD(1) * mrSges(4,3) - t669 * (mrSges(4,1) * t658 + mrSges(4,2) * t659);
t533 = m(4) * t614 + t659 * t679 + t536;
t686 = t666 * t543 - t662 * t554;
t534 = m(4) * t615 + t658 * t679 + t686;
t530 = t659 * t533 + t658 * t534;
t633 = -qJDD(1) * pkin(1) + t677;
t675 = -m(3) * t633 + t669 * mrSges(3,3) - t530;
t528 = m(2) * t638 - t669 * mrSges(2,2) + t702 * qJDD(1) + t675;
t629 = t669 * pkin(1) - t704;
t625 = t699 * t669 + t676;
t544 = t665 * t548 + t661 * t549;
t673 = m(5) * t605 - t617 * mrSges(5,1) + t618 * mrSges(5,2) - t634 * t626 + t635 * t627 + t544;
t672 = -m(4) * t625 - mrSges(4,1) * t692 - mrSges(4,2) * t691 - t673;
t671 = -m(3) * t629 + t669 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t672;
t540 = t671 + m(2) * t639 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t689) * t669;
t698 = t667 * t528 + t663 * t540;
t688 = -t663 * t528 + t667 * t540;
t687 = -t658 * t533 + t659 * t534;
t683 = Ifges(4,1) * t659 - Ifges(4,4) * t658;
t682 = Ifges(4,4) * t659 - Ifges(4,2) * t658;
t681 = Ifges(4,5) * t659 - Ifges(4,6) * t658;
t637 = t681 * qJD(1);
t608 = Ifges(5,1) * t635 + Ifges(5,4) * t634 + Ifges(5,5) * qJD(4);
t607 = Ifges(5,4) * t635 + Ifges(5,2) * t634 + Ifges(5,6) * qJD(4);
t606 = Ifges(5,5) * t635 + Ifges(5,6) * t634 + Ifges(5,3) * qJD(4);
t587 = Ifges(6,1) * t622 + Ifges(6,4) * t621 + Ifges(6,5) * t632;
t586 = Ifges(6,4) * t622 + Ifges(6,2) * t621 + Ifges(6,6) * t632;
t585 = Ifges(6,5) * t622 + Ifges(6,6) * t621 + Ifges(6,3) * t632;
t578 = Ifges(7,1) * t593 + Ifges(7,4) * t592 + Ifges(7,5) * t630;
t577 = Ifges(7,4) * t593 + Ifges(7,2) * t592 + Ifges(7,6) * t630;
t576 = Ifges(7,5) * t593 + Ifges(7,6) * t592 + Ifges(7,3) * t630;
t552 = mrSges(7,2) * t564 - mrSges(7,3) * t560 + Ifges(7,1) * t569 + Ifges(7,4) * t568 + Ifges(7,5) * t610 + t592 * t576 - t630 * t577;
t551 = -mrSges(7,1) * t564 + mrSges(7,3) * t561 + Ifges(7,4) * t569 + Ifges(7,2) * t568 + Ifges(7,6) * t610 - t593 * t576 + t630 * t578;
t537 = mrSges(6,2) * t571 - mrSges(6,3) * t565 + Ifges(6,1) * t591 + Ifges(6,4) * t590 + Ifges(6,5) * t613 - pkin(9) * t550 - t660 * t551 + t664 * t552 + t621 * t585 - t632 * t586;
t535 = -mrSges(6,1) * t571 + mrSges(6,3) * t566 + Ifges(6,4) * t591 + Ifges(6,2) * t590 + Ifges(6,6) * t613 - pkin(5) * t674 + pkin(9) * t684 + t664 * t551 + t660 * t552 - t622 * t585 + t632 * t587;
t531 = Ifges(5,4) * t618 + Ifges(5,2) * t617 + Ifges(5,6) * qJDD(4) - t635 * t606 + qJD(4) * t608 - mrSges(5,1) * t605 + mrSges(5,3) * t582 - Ifges(6,5) * t591 - Ifges(6,6) * t590 - Ifges(6,3) * t613 - t622 * t586 + t621 * t587 - mrSges(6,1) * t565 + mrSges(6,2) * t566 - Ifges(7,5) * t569 - Ifges(7,6) * t568 - Ifges(7,3) * t610 - t593 * t577 + t592 * t578 - mrSges(7,1) * t560 + mrSges(7,2) * t561 - pkin(5) * t550 - pkin(4) * t544;
t529 = -m(3) * g(3) + t687;
t526 = mrSges(5,2) * t605 - mrSges(5,3) * t581 + Ifges(5,1) * t618 + Ifges(5,4) * t617 + Ifges(5,5) * qJDD(4) - pkin(8) * t544 - qJD(4) * t607 - t661 * t535 + t665 * t537 + t634 * t606;
t525 = mrSges(4,2) * t625 - mrSges(4,3) * t614 - pkin(7) * t536 + t683 * qJDD(1) + t666 * t526 - t662 * t531 - t637 * t695;
t524 = -mrSges(4,1) * t625 + mrSges(4,3) * t615 - pkin(3) * t673 + pkin(7) * t686 + t682 * qJDD(1) + t662 * t526 + t666 * t531 - t637 * t694;
t523 = -qJ(2) * t529 + Ifges(5,3) * qJDD(4) + t665 * t535 + pkin(8) * t685 + t661 * t537 + t635 * t607 - mrSges(2,3) * t638 + mrSges(3,1) * t633 - t634 * t608 + pkin(4) * t670 + mrSges(4,1) * t614 - mrSges(4,2) * t615 + Ifges(5,6) * t617 + Ifges(5,5) * t618 + mrSges(5,1) * t581 - mrSges(5,2) * t582 + pkin(2) * t530 + pkin(3) * t536 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t681 + t701) * qJDD(1) + (t658 * t683 + t659 * t682 + t700) * t669;
t522 = mrSges(2,3) * t639 - mrSges(3,1) * t629 - t658 * t525 - t659 * t524 - pkin(2) * t672 - qJ(3) * t687 - pkin(1) * t529 - t700 * qJDD(1) + t702 * g(3) + (-pkin(2) * t689 + t701) * t669;
t1 = [-m(1) * g(1) + t688; -m(1) * g(2) + t698; (-m(1) - m(2) - m(3)) * g(3) + t687; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t698 - t663 * t522 + t667 * t523; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t688 + t667 * t522 + t663 * t523; pkin(1) * t675 + qJ(2) * (-t669 * t689 + t671) + t659 * t525 - t658 * t524 - qJ(3) * t530 + mrSges(2,1) * t638 - mrSges(2,2) * t639 + mrSges(3,2) * t633 - mrSges(3,3) * t629 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
