% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:27:06
% EndTime: 2019-05-06 04:27:19
% DurationCPUTime: 8.68s
% Computational Cost: add. (142956->338), mult. (283467->412), div. (0->0), fcn. (191647->10), ass. (0->133)
t647 = sin(qJ(1));
t652 = cos(qJ(1));
t631 = -t652 * g(1) - t647 * g(2);
t677 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t631;
t676 = -pkin(1) - pkin(7);
t675 = mrSges(2,1) - mrSges(3,2);
t674 = -Ifges(3,4) + Ifges(2,5);
t673 = (Ifges(3,5) - Ifges(2,6));
t630 = t647 * g(1) - t652 * g(2);
t654 = qJD(1) ^ 2;
t660 = -t654 * qJ(2) + qJDD(2) - t630;
t607 = qJDD(1) * t676 + t660;
t646 = sin(qJ(3));
t651 = cos(qJ(3));
t600 = -t651 * g(3) + t646 * t607;
t623 = (mrSges(4,1) * t646 + mrSges(4,2) * t651) * qJD(1);
t670 = qJD(1) * qJD(3);
t634 = t651 * t670;
t625 = -t646 * qJDD(1) - t634;
t671 = qJD(1) * t651;
t629 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t671;
t636 = t646 * qJD(1);
t606 = t654 * t676 - t677;
t668 = t646 * t670;
t626 = t651 * qJDD(1) - t668;
t580 = (-t626 + t668) * pkin(8) + (-t625 + t634) * pkin(3) + t606;
t624 = (pkin(3) * t646 - pkin(8) * t651) * qJD(1);
t653 = qJD(3) ^ 2;
t583 = -t653 * pkin(3) + qJDD(3) * pkin(8) - t624 * t636 + t600;
t645 = sin(qJ(4));
t650 = cos(qJ(4));
t564 = t650 * t580 - t645 * t583;
t621 = t650 * qJD(3) - t645 * t671;
t594 = t621 * qJD(4) + t645 * qJDD(3) + t650 * t626;
t620 = qJDD(4) - t625;
t622 = t645 * qJD(3) + t650 * t671;
t633 = t636 + qJD(4);
t554 = (t621 * t633 - t594) * pkin(9) + (t621 * t622 + t620) * pkin(4) + t564;
t565 = t645 * t580 + t650 * t583;
t593 = -t622 * qJD(4) + t650 * qJDD(3) - t645 * t626;
t605 = t633 * pkin(4) - t622 * pkin(9);
t619 = t621 ^ 2;
t556 = -t619 * pkin(4) + t593 * pkin(9) - t633 * t605 + t565;
t644 = sin(qJ(5));
t649 = cos(qJ(5));
t544 = t649 * t554 - t644 * t556;
t596 = t649 * t621 - t644 * t622;
t568 = t596 * qJD(5) + t644 * t593 + t649 * t594;
t597 = t644 * t621 + t649 * t622;
t615 = qJDD(5) + t620;
t632 = qJD(5) + t633;
t542 = (t596 * t632 - t568) * pkin(10) + (t596 * t597 + t615) * pkin(5) + t544;
t545 = t644 * t554 + t649 * t556;
t567 = -t597 * qJD(5) + t649 * t593 - t644 * t594;
t586 = t632 * pkin(5) - t597 * pkin(10);
t595 = t596 ^ 2;
t543 = -t595 * pkin(5) + t567 * pkin(10) - t632 * t586 + t545;
t643 = sin(qJ(6));
t648 = cos(qJ(6));
t540 = t648 * t542 - t643 * t543;
t575 = t648 * t596 - t643 * t597;
t551 = t575 * qJD(6) + t643 * t567 + t648 * t568;
t576 = t643 * t596 + t648 * t597;
t562 = -mrSges(7,1) * t575 + mrSges(7,2) * t576;
t627 = qJD(6) + t632;
t569 = -t627 * mrSges(7,2) + t575 * mrSges(7,3);
t611 = qJDD(6) + t615;
t538 = m(7) * t540 + t611 * mrSges(7,1) - t551 * mrSges(7,3) - t576 * t562 + t627 * t569;
t541 = t643 * t542 + t648 * t543;
t550 = -t576 * qJD(6) + t648 * t567 - t643 * t568;
t570 = t627 * mrSges(7,1) - t576 * mrSges(7,3);
t539 = m(7) * t541 - t611 * mrSges(7,2) + t550 * mrSges(7,3) + t575 * t562 - t627 * t570;
t531 = t648 * t538 + t643 * t539;
t577 = -t596 * mrSges(6,1) + mrSges(6,2) * t597;
t584 = -t632 * mrSges(6,2) + t596 * mrSges(6,3);
t529 = m(6) * t544 + t615 * mrSges(6,1) - t568 * mrSges(6,3) - t597 * t577 + t632 * t584 + t531;
t585 = t632 * mrSges(6,1) - t597 * mrSges(6,3);
t663 = -t643 * t538 + t648 * t539;
t530 = m(6) * t545 - t615 * mrSges(6,2) + t567 * mrSges(6,3) + t596 * t577 - t632 * t585 + t663;
t525 = t649 * t529 + t644 * t530;
t598 = -t621 * mrSges(5,1) + t622 * mrSges(5,2);
t601 = -t633 * mrSges(5,2) + t621 * mrSges(5,3);
t523 = m(5) * t564 + t620 * mrSges(5,1) - t594 * mrSges(5,3) - t622 * t598 + t633 * t601 + t525;
t602 = t633 * mrSges(5,1) - t622 * mrSges(5,3);
t664 = -t644 * t529 + t649 * t530;
t524 = m(5) * t565 - t620 * mrSges(5,2) + t593 * mrSges(5,3) + t621 * t598 - t633 * t602 + t664;
t665 = -t645 * t523 + t650 * t524;
t516 = m(4) * t600 - qJDD(3) * mrSges(4,2) + t625 * mrSges(4,3) - qJD(3) * t629 - t623 * t636 + t665;
t599 = t646 * g(3) + t651 * t607;
t628 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t636;
t582 = -qJDD(3) * pkin(3) - t653 * pkin(8) + t624 * t671 - t599;
t563 = -t593 * pkin(4) - t619 * pkin(9) + t622 * t605 + t582;
t547 = -t567 * pkin(5) - t595 * pkin(10) + t597 * t586 + t563;
t662 = m(7) * t547 - t550 * mrSges(7,1) + t551 * mrSges(7,2) - t575 * t569 + t576 * t570;
t656 = m(6) * t563 - t567 * mrSges(6,1) + t568 * mrSges(6,2) - t596 * t584 + t597 * t585 + t662;
t655 = -m(5) * t582 + t593 * mrSges(5,1) - t594 * mrSges(5,2) + t621 * t601 - t622 * t602 - t656;
t534 = m(4) * t599 + qJDD(3) * mrSges(4,1) - t626 * mrSges(4,3) + qJD(3) * t628 - t623 * t671 + t655;
t511 = t646 * t516 + t651 * t534;
t609 = -qJDD(1) * pkin(1) + t660;
t659 = -m(3) * t609 + (t654 * mrSges(3,3)) - t511;
t509 = m(2) * t630 - (t654 * mrSges(2,2)) + qJDD(1) * t675 + t659;
t608 = t654 * pkin(1) + t677;
t517 = t650 * t523 + t645 * t524;
t658 = -m(4) * t606 + t625 * mrSges(4,1) - t626 * mrSges(4,2) - t628 * t636 - t629 * t671 - t517;
t657 = -m(3) * t608 + (t654 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t658;
t514 = m(2) * t631 - (t654 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t657;
t672 = t652 * t509 + t647 * t514;
t667 = -t647 * t509 + t652 * t514;
t666 = t651 * t516 - t646 * t534;
t614 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t651 - Ifges(4,4) * t646) * qJD(1);
t613 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t651 - Ifges(4,2) * t646) * qJD(1);
t612 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t651 - Ifges(4,6) * t646) * qJD(1);
t589 = Ifges(5,1) * t622 + Ifges(5,4) * t621 + Ifges(5,5) * t633;
t588 = Ifges(5,4) * t622 + Ifges(5,2) * t621 + Ifges(5,6) * t633;
t587 = Ifges(5,5) * t622 + Ifges(5,6) * t621 + Ifges(5,3) * t633;
t573 = Ifges(6,1) * t597 + Ifges(6,4) * t596 + Ifges(6,5) * t632;
t572 = Ifges(6,4) * t597 + Ifges(6,2) * t596 + Ifges(6,6) * t632;
t571 = Ifges(6,5) * t597 + Ifges(6,6) * t596 + Ifges(6,3) * t632;
t559 = Ifges(7,1) * t576 + Ifges(7,4) * t575 + Ifges(7,5) * t627;
t558 = Ifges(7,4) * t576 + Ifges(7,2) * t575 + Ifges(7,6) * t627;
t557 = Ifges(7,5) * t576 + Ifges(7,6) * t575 + Ifges(7,3) * t627;
t533 = mrSges(7,2) * t547 - mrSges(7,3) * t540 + Ifges(7,1) * t551 + Ifges(7,4) * t550 + Ifges(7,5) * t611 + t575 * t557 - t627 * t558;
t532 = -mrSges(7,1) * t547 + mrSges(7,3) * t541 + Ifges(7,4) * t551 + Ifges(7,2) * t550 + Ifges(7,6) * t611 - t576 * t557 + t627 * t559;
t519 = mrSges(6,2) * t563 - mrSges(6,3) * t544 + Ifges(6,1) * t568 + Ifges(6,4) * t567 + Ifges(6,5) * t615 - pkin(10) * t531 - t643 * t532 + t648 * t533 + t596 * t571 - t632 * t572;
t518 = -mrSges(6,1) * t563 + mrSges(6,3) * t545 + Ifges(6,4) * t568 + Ifges(6,2) * t567 + Ifges(6,6) * t615 - pkin(5) * t662 + pkin(10) * t663 + t648 * t532 + t643 * t533 - t597 * t571 + t632 * t573;
t510 = -m(3) * g(3) + t666;
t507 = mrSges(5,2) * t582 - mrSges(5,3) * t564 + Ifges(5,1) * t594 + Ifges(5,4) * t593 + Ifges(5,5) * t620 - pkin(9) * t525 - t644 * t518 + t649 * t519 + t621 * t587 - t633 * t588;
t506 = -mrSges(5,1) * t582 + mrSges(5,3) * t565 + Ifges(5,4) * t594 + Ifges(5,2) * t593 + Ifges(5,6) * t620 - pkin(4) * t656 + pkin(9) * t664 + t649 * t518 + t644 * t519 - t622 * t587 + t633 * t589;
t505 = Ifges(4,6) * qJDD(3) + Ifges(4,4) * t626 - Ifges(5,3) * t620 + t621 * t589 - t622 * t588 + Ifges(4,2) * t625 - mrSges(4,1) * t606 - Ifges(7,3) * t611 + qJD(3) * t614 - Ifges(6,3) * t615 + t596 * t573 - t597 * t572 + mrSges(4,3) * t600 - Ifges(5,6) * t593 - Ifges(5,5) * t594 + t575 * t559 - t576 * t558 + mrSges(5,2) * t565 - Ifges(6,6) * t567 - Ifges(6,5) * t568 - mrSges(5,1) * t564 - Ifges(7,6) * t550 - Ifges(7,5) * t551 - mrSges(6,1) * t544 + mrSges(6,2) * t545 + mrSges(7,2) * t541 - mrSges(7,1) * t540 - pkin(5) * t531 - t612 * t671 - pkin(4) * t525 - pkin(3) * t517;
t504 = mrSges(4,2) * t606 - mrSges(4,3) * t599 + Ifges(4,1) * t626 + Ifges(4,4) * t625 + Ifges(4,5) * qJDD(3) - pkin(8) * t517 - qJD(3) * t613 - t645 * t506 + t650 * t507 - t612 * t636;
t503 = -qJ(2) * t510 - mrSges(2,3) * t630 + pkin(2) * t511 + mrSges(3,1) * t609 + t650 * t506 + pkin(3) * t655 + pkin(8) * t665 + t645 * t507 + Ifges(4,5) * t626 + Ifges(4,6) * t625 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t599 - mrSges(4,2) * t600 + (t673 * t654) + t674 * qJDD(1) + (t651 * t613 + t646 * t614) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t502 = -mrSges(3,1) * t608 + mrSges(2,3) * t631 - pkin(1) * t510 - pkin(2) * t658 - pkin(7) * t666 + g(3) * t675 - qJDD(1) * t673 - t646 * t504 - t651 * t505 + t654 * t674;
t1 = [-m(1) * g(1) + t667; -m(1) * g(2) + t672; (-m(1) - m(2) - m(3)) * g(3) + t666; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t672 - t647 * t502 + t652 * t503; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t667 + t652 * t502 + t647 * t503; qJ(2) * t657 + pkin(1) * t659 + t651 * t504 - t646 * t505 - pkin(7) * t511 + mrSges(2,1) * t630 - mrSges(2,2) * t631 - mrSges(3,3) * t608 + mrSges(3,2) * t609 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
