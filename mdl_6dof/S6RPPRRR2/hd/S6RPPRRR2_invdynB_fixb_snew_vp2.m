% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:20:53
% EndTime: 2019-05-05 15:21:05
% DurationCPUTime: 11.27s
% Computational Cost: add. (175612->323), mult. (392906->403), div. (0->0), fcn. (281659->12), ass. (0->138)
t663 = qJD(1) ^ 2;
t652 = cos(pkin(11));
t690 = pkin(3) * t652;
t650 = sin(pkin(11));
t689 = mrSges(4,2) * t650;
t648 = t652 ^ 2;
t688 = t648 * t663;
t657 = sin(qJ(1));
t661 = cos(qJ(1));
t634 = t657 * g(1) - g(2) * t661;
t632 = qJDD(1) * pkin(1) + t634;
t635 = -g(1) * t661 - g(2) * t657;
t633 = -pkin(1) * t663 + t635;
t651 = sin(pkin(10));
t653 = cos(pkin(10));
t616 = t651 * t632 + t653 * t633;
t606 = -pkin(2) * t663 + qJDD(1) * qJ(3) + t616;
t649 = -g(3) + qJDD(2);
t682 = qJD(1) * qJD(3);
t686 = t652 * t649 - 0.2e1 * t650 * t682;
t587 = (-pkin(7) * qJDD(1) + t663 * t690 - t606) * t650 + t686;
t595 = t650 * t649 + (t606 + 0.2e1 * t682) * t652;
t681 = qJDD(1) * t652;
t590 = -pkin(3) * t688 + pkin(7) * t681 + t595;
t656 = sin(qJ(4));
t660 = cos(qJ(4));
t571 = t656 * t587 + t660 * t590;
t684 = qJD(1) * t652;
t685 = qJD(1) * t650;
t625 = -t656 * t685 + t660 * t684;
t669 = t650 * t660 + t652 * t656;
t626 = t669 * qJD(1);
t609 = -mrSges(5,1) * t625 + mrSges(5,2) * t626;
t623 = t626 * qJD(4);
t613 = -t656 * t650 * qJDD(1) + t660 * t681 - t623;
t621 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t626;
t612 = -pkin(4) * t625 - pkin(8) * t626;
t662 = qJD(4) ^ 2;
t566 = -pkin(4) * t662 + qJDD(4) * pkin(8) + t612 * t625 + t571;
t647 = t650 ^ 2;
t615 = t653 * t632 - t651 * t633;
t670 = qJDD(3) - t615;
t591 = (-pkin(2) - t690) * qJDD(1) + (-qJ(3) + (-t647 - t648) * pkin(7)) * t663 + t670;
t683 = t625 * qJD(4);
t614 = qJDD(1) * t669 + t683;
t569 = (-t614 - t683) * pkin(8) + (-t613 + t623) * pkin(4) + t591;
t655 = sin(qJ(5));
t659 = cos(qJ(5));
t558 = -t655 * t566 + t659 * t569;
t618 = qJD(4) * t659 - t626 * t655;
t589 = qJD(5) * t618 + qJDD(4) * t655 + t614 * t659;
t611 = qJDD(5) - t613;
t619 = qJD(4) * t655 + t626 * t659;
t624 = qJD(5) - t625;
t556 = (t618 * t624 - t589) * pkin(9) + (t618 * t619 + t611) * pkin(5) + t558;
t559 = t659 * t566 + t655 * t569;
t588 = -qJD(5) * t619 + qJDD(4) * t659 - t614 * t655;
t600 = pkin(5) * t624 - pkin(9) * t619;
t617 = t618 ^ 2;
t557 = -pkin(5) * t617 + pkin(9) * t588 - t600 * t624 + t559;
t654 = sin(qJ(6));
t658 = cos(qJ(6));
t554 = t556 * t658 - t557 * t654;
t592 = t618 * t658 - t619 * t654;
t564 = qJD(6) * t592 + t588 * t654 + t589 * t658;
t593 = t618 * t654 + t619 * t658;
t576 = -mrSges(7,1) * t592 + mrSges(7,2) * t593;
t622 = qJD(6) + t624;
t577 = -mrSges(7,2) * t622 + mrSges(7,3) * t592;
t608 = qJDD(6) + t611;
t552 = m(7) * t554 + mrSges(7,1) * t608 - mrSges(7,3) * t564 - t576 * t593 + t577 * t622;
t555 = t556 * t654 + t557 * t658;
t563 = -qJD(6) * t593 + t588 * t658 - t589 * t654;
t578 = mrSges(7,1) * t622 - mrSges(7,3) * t593;
t553 = m(7) * t555 - mrSges(7,2) * t608 + mrSges(7,3) * t563 + t576 * t592 - t578 * t622;
t544 = t658 * t552 + t654 * t553;
t596 = -mrSges(6,1) * t618 + mrSges(6,2) * t619;
t598 = -mrSges(6,2) * t624 + mrSges(6,3) * t618;
t542 = m(6) * t558 + mrSges(6,1) * t611 - mrSges(6,3) * t589 - t596 * t619 + t598 * t624 + t544;
t599 = mrSges(6,1) * t624 - mrSges(6,3) * t619;
t674 = -t552 * t654 + t658 * t553;
t543 = m(6) * t559 - mrSges(6,2) * t611 + mrSges(6,3) * t588 + t596 * t618 - t599 * t624 + t674;
t675 = -t542 * t655 + t659 * t543;
t537 = m(5) * t571 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t613 - qJD(4) * t621 + t609 * t625 + t675;
t570 = t587 * t660 - t656 * t590;
t620 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t625;
t565 = -qJDD(4) * pkin(4) - pkin(8) * t662 + t626 * t612 - t570;
t560 = -pkin(5) * t588 - pkin(9) * t617 + t600 * t619 + t565;
t667 = m(7) * t560 - t563 * mrSges(7,1) + mrSges(7,2) * t564 - t592 * t577 + t578 * t593;
t664 = -m(6) * t565 + t588 * mrSges(6,1) - mrSges(6,2) * t589 + t618 * t598 - t599 * t619 - t667;
t548 = m(5) * t570 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t614 + qJD(4) * t620 - t609 * t626 + t664;
t530 = t656 * t537 + t660 * t548;
t594 = -t606 * t650 + t686;
t668 = mrSges(4,3) * qJDD(1) + t663 * (-mrSges(4,1) * t652 + t689);
t528 = m(4) * t594 - t650 * t668 + t530;
t676 = t660 * t537 - t656 * t548;
t529 = m(4) * t595 + t652 * t668 + t676;
t677 = -t528 * t650 + t652 * t529;
t522 = m(3) * t616 - mrSges(3,1) * t663 - qJDD(1) * mrSges(3,2) + t677;
t602 = -qJDD(1) * pkin(2) - t663 * qJ(3) + t670;
t538 = t659 * t542 + t655 * t543;
t666 = m(5) * t591 - t613 * mrSges(5,1) + t614 * mrSges(5,2) - t625 * t620 + t626 * t621 + t538;
t665 = -m(4) * t602 + mrSges(4,1) * t681 - t666 + (t647 * t663 + t688) * mrSges(4,3);
t534 = -t663 * mrSges(3,2) + m(3) * t615 + (mrSges(3,1) - t689) * qJDD(1) + t665;
t518 = t651 * t522 + t653 * t534;
t516 = m(2) * t634 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t663 + t518;
t678 = t653 * t522 - t534 * t651;
t517 = m(2) * t635 - mrSges(2,1) * t663 - qJDD(1) * mrSges(2,2) + t678;
t687 = t661 * t516 + t657 * t517;
t523 = t652 * t528 + t650 * t529;
t680 = m(3) * t649 + t523;
t679 = -t516 * t657 + t661 * t517;
t673 = Ifges(4,1) * t650 + Ifges(4,4) * t652;
t672 = Ifges(4,4) * t650 + Ifges(4,2) * t652;
t671 = Ifges(4,5) * t650 + Ifges(4,6) * t652;
t631 = t671 * qJD(1);
t605 = Ifges(5,1) * t626 + Ifges(5,4) * t625 + Ifges(5,5) * qJD(4);
t604 = Ifges(5,4) * t626 + Ifges(5,2) * t625 + Ifges(5,6) * qJD(4);
t603 = Ifges(5,5) * t626 + Ifges(5,6) * t625 + Ifges(5,3) * qJD(4);
t582 = Ifges(6,1) * t619 + Ifges(6,4) * t618 + Ifges(6,5) * t624;
t581 = Ifges(6,4) * t619 + Ifges(6,2) * t618 + Ifges(6,6) * t624;
t580 = Ifges(6,5) * t619 + Ifges(6,6) * t618 + Ifges(6,3) * t624;
t574 = Ifges(7,1) * t593 + Ifges(7,4) * t592 + Ifges(7,5) * t622;
t573 = Ifges(7,4) * t593 + Ifges(7,2) * t592 + Ifges(7,6) * t622;
t572 = Ifges(7,5) * t593 + Ifges(7,6) * t592 + Ifges(7,3) * t622;
t546 = mrSges(7,2) * t560 - mrSges(7,3) * t554 + Ifges(7,1) * t564 + Ifges(7,4) * t563 + Ifges(7,5) * t608 + t572 * t592 - t573 * t622;
t545 = -mrSges(7,1) * t560 + mrSges(7,3) * t555 + Ifges(7,4) * t564 + Ifges(7,2) * t563 + Ifges(7,6) * t608 - t572 * t593 + t574 * t622;
t532 = mrSges(6,2) * t565 - mrSges(6,3) * t558 + Ifges(6,1) * t589 + Ifges(6,4) * t588 + Ifges(6,5) * t611 - pkin(9) * t544 - t545 * t654 + t546 * t658 + t580 * t618 - t581 * t624;
t531 = -mrSges(6,1) * t565 + mrSges(6,3) * t559 + Ifges(6,4) * t589 + Ifges(6,2) * t588 + Ifges(6,6) * t611 - pkin(5) * t667 + pkin(9) * t674 + t658 * t545 + t654 * t546 - t619 * t580 + t624 * t582;
t524 = Ifges(5,4) * t614 + Ifges(5,2) * t613 + Ifges(5,6) * qJDD(4) - t626 * t603 + qJD(4) * t605 - mrSges(5,1) * t591 + mrSges(5,3) * t571 - Ifges(6,5) * t589 - Ifges(6,6) * t588 - Ifges(6,3) * t611 - t619 * t581 + t618 * t582 - mrSges(6,1) * t558 + mrSges(6,2) * t559 - Ifges(7,5) * t564 - Ifges(7,6) * t563 - Ifges(7,3) * t608 - t593 * t573 + t592 * t574 - mrSges(7,1) * t554 + mrSges(7,2) * t555 - pkin(5) * t544 - pkin(4) * t538;
t519 = mrSges(5,2) * t591 - mrSges(5,3) * t570 + Ifges(5,1) * t614 + Ifges(5,4) * t613 + Ifges(5,5) * qJDD(4) - pkin(8) * t538 - qJD(4) * t604 - t531 * t655 + t532 * t659 + t603 * t625;
t512 = mrSges(4,2) * t602 - mrSges(4,3) * t594 - pkin(7) * t530 + qJDD(1) * t673 + t660 * t519 - t656 * t524 + t631 * t684;
t511 = -pkin(2) * t523 - mrSges(3,1) * t649 + mrSges(3,3) * t616 - pkin(3) * t530 - mrSges(4,1) * t594 + mrSges(4,2) * t595 - pkin(4) * t664 - pkin(8) * t675 - t655 * t532 - t659 * t531 - mrSges(5,1) * t570 + mrSges(5,2) * t571 - Ifges(5,5) * t614 - Ifges(5,6) * t613 - Ifges(5,3) * qJDD(4) - t626 * t604 + t625 * t605 + (Ifges(3,6) - t671) * qJDD(1) + (-t650 * t672 + t652 * t673 + Ifges(3,5)) * t663;
t510 = -mrSges(4,1) * t602 + mrSges(4,3) * t595 - pkin(3) * t666 + pkin(7) * t676 + qJDD(1) * t672 + t656 * t519 + t660 * t524 - t631 * t685;
t509 = mrSges(3,2) * t649 - mrSges(3,3) * t615 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t663 - qJ(3) * t523 - t510 * t650 + t512 * t652;
t508 = -mrSges(2,2) * g(3) - mrSges(2,3) * t634 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t663 - qJ(2) * t518 + t509 * t653 - t511 * t651;
t507 = mrSges(2,1) * g(3) + mrSges(2,3) * t635 + t663 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t680 + qJ(2) * t678 + t651 * t509 + t653 * t511;
t1 = [-m(1) * g(1) + t679; -m(1) * g(2) + t687; (-m(1) - m(2)) * g(3) + t680; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t687 - t657 * t507 + t661 * t508; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t679 + t661 * t507 + t657 * t508; pkin(1) * t518 + mrSges(2,1) * t634 - mrSges(2,2) * t635 + t650 * t512 + t652 * t510 + pkin(2) * t665 + qJ(3) * t677 + mrSges(3,1) * t615 - mrSges(3,2) * t616 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t689 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
