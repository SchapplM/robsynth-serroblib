% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:14:01
% EndTime: 2019-05-05 18:14:14
% DurationCPUTime: 12.68s
% Computational Cost: add. (204822->343), mult. (451295->434), div. (0->0), fcn. (314798->12), ass. (0->134)
t654 = sin(qJ(1));
t658 = cos(qJ(1));
t636 = t654 * g(1) - g(2) * t658;
t628 = qJDD(1) * pkin(1) + t636;
t637 = -g(1) * t658 - g(2) * t654;
t659 = qJD(1) ^ 2;
t630 = -pkin(1) * t659 + t637;
t648 = sin(pkin(10));
t650 = cos(pkin(10));
t608 = t648 * t628 + t650 * t630;
t602 = -pkin(2) * t659 + qJDD(1) * pkin(7) + t608;
t646 = -g(3) + qJDD(2);
t653 = sin(qJ(3));
t657 = cos(qJ(3));
t590 = -t653 * t602 + t657 * t646;
t673 = qJD(1) * qJD(3);
t671 = t657 * t673;
t631 = qJDD(1) * t653 + t671;
t585 = (-t631 + t671) * qJ(4) + (t653 * t657 * t659 + qJDD(3)) * pkin(3) + t590;
t591 = t657 * t602 + t653 * t646;
t632 = qJDD(1) * t657 - t653 * t673;
t675 = qJD(1) * t653;
t633 = qJD(3) * pkin(3) - qJ(4) * t675;
t645 = t657 ^ 2;
t586 = -pkin(3) * t645 * t659 + qJ(4) * t632 - qJD(3) * t633 + t591;
t647 = sin(pkin(11));
t649 = cos(pkin(11));
t618 = (t647 * t657 + t649 * t653) * qJD(1);
t556 = -0.2e1 * qJD(4) * t618 + t649 * t585 - t647 * t586;
t610 = t631 * t649 + t632 * t647;
t617 = (-t647 * t653 + t649 * t657) * qJD(1);
t553 = (qJD(3) * t617 - t610) * pkin(8) + (t617 * t618 + qJDD(3)) * pkin(4) + t556;
t557 = 0.2e1 * qJD(4) * t617 + t647 * t585 + t649 * t586;
t609 = -t631 * t647 + t632 * t649;
t613 = qJD(3) * pkin(4) - pkin(8) * t618;
t616 = t617 ^ 2;
t555 = -pkin(4) * t616 + pkin(8) * t609 - qJD(3) * t613 + t557;
t652 = sin(qJ(5));
t656 = cos(qJ(5));
t550 = t652 * t553 + t656 * t555;
t600 = t617 * t652 + t618 * t656;
t570 = -qJD(5) * t600 + t609 * t656 - t610 * t652;
t599 = t617 * t656 - t618 * t652;
t583 = -mrSges(6,1) * t599 + mrSges(6,2) * t600;
t644 = qJD(3) + qJD(5);
t593 = mrSges(6,1) * t644 - mrSges(6,3) * t600;
t643 = qJDD(3) + qJDD(5);
t584 = -pkin(5) * t599 - pkin(9) * t600;
t642 = t644 ^ 2;
t548 = -pkin(5) * t642 + pkin(9) * t643 + t584 * t599 + t550;
t607 = t650 * t628 - t648 * t630;
t664 = -qJDD(1) * pkin(2) - t607;
t587 = -t632 * pkin(3) + qJDD(4) + t633 * t675 + (-qJ(4) * t645 - pkin(7)) * t659 + t664;
t562 = -t609 * pkin(4) - t616 * pkin(8) + t618 * t613 + t587;
t571 = qJD(5) * t599 + t609 * t652 + t610 * t656;
t551 = (-t599 * t644 - t571) * pkin(9) + (t600 * t644 - t570) * pkin(5) + t562;
t651 = sin(qJ(6));
t655 = cos(qJ(6));
t545 = -t548 * t651 + t551 * t655;
t588 = -t600 * t651 + t644 * t655;
t560 = qJD(6) * t588 + t571 * t655 + t643 * t651;
t569 = qJDD(6) - t570;
t589 = t600 * t655 + t644 * t651;
t572 = -mrSges(7,1) * t588 + mrSges(7,2) * t589;
t595 = qJD(6) - t599;
t573 = -mrSges(7,2) * t595 + mrSges(7,3) * t588;
t543 = m(7) * t545 + mrSges(7,1) * t569 - mrSges(7,3) * t560 - t572 * t589 + t573 * t595;
t546 = t548 * t655 + t551 * t651;
t559 = -qJD(6) * t589 - t571 * t651 + t643 * t655;
t574 = mrSges(7,1) * t595 - mrSges(7,3) * t589;
t544 = m(7) * t546 - mrSges(7,2) * t569 + mrSges(7,3) * t559 + t572 * t588 - t574 * t595;
t665 = -t543 * t651 + t655 * t544;
t534 = m(6) * t550 - mrSges(6,2) * t643 + mrSges(6,3) * t570 + t583 * t599 - t593 * t644 + t665;
t549 = t553 * t656 - t555 * t652;
t592 = -mrSges(6,2) * t644 + mrSges(6,3) * t599;
t547 = -pkin(5) * t643 - pkin(9) * t642 + t584 * t600 - t549;
t662 = -m(7) * t547 + t559 * mrSges(7,1) - mrSges(7,2) * t560 + t588 * t573 - t574 * t589;
t539 = m(6) * t549 + mrSges(6,1) * t643 - mrSges(6,3) * t571 - t583 * t600 + t592 * t644 + t662;
t529 = t652 * t534 + t656 * t539;
t605 = -mrSges(5,1) * t617 + mrSges(5,2) * t618;
t611 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t617;
t527 = m(5) * t556 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t610 + qJD(3) * t611 - t605 * t618 + t529;
t612 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t618;
t666 = t656 * t534 - t539 * t652;
t528 = m(5) * t557 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t609 - qJD(3) * t612 + t605 * t617 + t666;
t521 = t649 * t527 + t647 * t528;
t629 = (-mrSges(4,1) * t657 + mrSges(4,2) * t653) * qJD(1);
t674 = qJD(1) * t657;
t635 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t674;
t519 = m(4) * t590 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t631 + qJD(3) * t635 - t629 * t675 + t521;
t634 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t675;
t667 = -t527 * t647 + t649 * t528;
t520 = m(4) * t591 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t632 - qJD(3) * t634 + t629 * t674 + t667;
t668 = -t519 * t653 + t657 * t520;
t513 = m(3) * t608 - mrSges(3,1) * t659 - qJDD(1) * mrSges(3,2) + t668;
t601 = -t659 * pkin(7) + t664;
t535 = t655 * t543 + t651 * t544;
t663 = m(6) * t562 - t570 * mrSges(6,1) + t571 * mrSges(6,2) - t599 * t592 + t600 * t593 + t535;
t661 = m(5) * t587 - t609 * mrSges(5,1) + mrSges(5,2) * t610 - t617 * t611 + t612 * t618 + t663;
t660 = -m(4) * t601 + t632 * mrSges(4,1) - mrSges(4,2) * t631 - t634 * t675 + t635 * t674 - t661;
t531 = m(3) * t607 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t659 + t660;
t509 = t648 * t513 + t650 * t531;
t507 = m(2) * t636 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t659 + t509;
t669 = t650 * t513 - t531 * t648;
t508 = m(2) * t637 - mrSges(2,1) * t659 - qJDD(1) * mrSges(2,2) + t669;
t676 = t658 * t507 + t654 * t508;
t514 = t657 * t519 + t653 * t520;
t672 = m(3) * t646 + t514;
t670 = -t507 * t654 + t658 * t508;
t624 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t653 + Ifges(4,4) * t657) * qJD(1);
t623 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t653 + Ifges(4,2) * t657) * qJD(1);
t622 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t653 + Ifges(4,6) * t657) * qJD(1);
t598 = Ifges(5,1) * t618 + Ifges(5,4) * t617 + Ifges(5,5) * qJD(3);
t597 = Ifges(5,4) * t618 + Ifges(5,2) * t617 + Ifges(5,6) * qJD(3);
t596 = Ifges(5,5) * t618 + Ifges(5,6) * t617 + Ifges(5,3) * qJD(3);
t577 = Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t644;
t576 = Ifges(6,4) * t600 + Ifges(6,2) * t599 + Ifges(6,6) * t644;
t575 = Ifges(6,5) * t600 + Ifges(6,6) * t599 + Ifges(6,3) * t644;
t565 = Ifges(7,1) * t589 + Ifges(7,4) * t588 + Ifges(7,5) * t595;
t564 = Ifges(7,4) * t589 + Ifges(7,2) * t588 + Ifges(7,6) * t595;
t563 = Ifges(7,5) * t589 + Ifges(7,6) * t588 + Ifges(7,3) * t595;
t537 = mrSges(7,2) * t547 - mrSges(7,3) * t545 + Ifges(7,1) * t560 + Ifges(7,4) * t559 + Ifges(7,5) * t569 + t563 * t588 - t564 * t595;
t536 = -mrSges(7,1) * t547 + mrSges(7,3) * t546 + Ifges(7,4) * t560 + Ifges(7,2) * t559 + Ifges(7,6) * t569 - t563 * t589 + t565 * t595;
t523 = -mrSges(6,1) * t562 - mrSges(7,1) * t545 + mrSges(7,2) * t546 + mrSges(6,3) * t550 + Ifges(6,4) * t571 - Ifges(7,5) * t560 + Ifges(6,2) * t570 + Ifges(6,6) * t643 - Ifges(7,6) * t559 - Ifges(7,3) * t569 - pkin(5) * t535 - t564 * t589 + t565 * t588 - t575 * t600 + t577 * t644;
t522 = mrSges(6,2) * t562 - mrSges(6,3) * t549 + Ifges(6,1) * t571 + Ifges(6,4) * t570 + Ifges(6,5) * t643 - pkin(9) * t535 - t536 * t651 + t537 * t655 + t575 * t599 - t576 * t644;
t515 = mrSges(5,2) * t587 - mrSges(5,3) * t556 + Ifges(5,1) * t610 + Ifges(5,4) * t609 + Ifges(5,5) * qJDD(3) - pkin(8) * t529 - qJD(3) * t597 + t522 * t656 - t523 * t652 + t596 * t617;
t510 = -mrSges(5,1) * t587 + mrSges(5,3) * t557 + Ifges(5,4) * t610 + Ifges(5,2) * t609 + Ifges(5,6) * qJDD(3) - pkin(4) * t663 + pkin(8) * t666 + qJD(3) * t598 + t652 * t522 + t656 * t523 - t618 * t596;
t503 = Ifges(3,6) * qJDD(1) - pkin(2) * t514 - pkin(9) * t665 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t623 * t653 + t624 * t657) * qJD(1) + t659 * Ifges(3,5) - t651 * t537 - t655 * t536 - Ifges(6,3) * t643 - mrSges(3,1) * t646 - Ifges(4,6) * t632 + t617 * t598 - t618 * t597 - Ifges(4,5) * t631 + mrSges(3,3) * t608 - Ifges(5,6) * t609 - Ifges(5,5) * t610 + t599 * t577 - t600 * t576 - mrSges(4,1) * t590 + mrSges(4,2) * t591 - Ifges(6,6) * t570 - Ifges(6,5) * t571 - mrSges(5,1) * t556 + mrSges(5,2) * t557 - mrSges(6,1) * t549 + mrSges(6,2) * t550 - pkin(4) * t529 - pkin(3) * t521 - pkin(5) * t662;
t502 = mrSges(4,2) * t601 - mrSges(4,3) * t590 + Ifges(4,1) * t631 + Ifges(4,4) * t632 + Ifges(4,5) * qJDD(3) - qJ(4) * t521 - qJD(3) * t623 - t510 * t647 + t515 * t649 + t622 * t674;
t501 = -mrSges(4,1) * t601 + mrSges(4,3) * t591 + Ifges(4,4) * t631 + Ifges(4,2) * t632 + Ifges(4,6) * qJDD(3) - pkin(3) * t661 + qJ(4) * t667 + qJD(3) * t624 + t649 * t510 + t647 * t515 - t622 * t675;
t500 = mrSges(3,2) * t646 - mrSges(3,3) * t607 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t659 - pkin(7) * t514 - t501 * t653 + t502 * t657;
t499 = -mrSges(2,2) * g(3) - mrSges(2,3) * t636 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t659 - qJ(2) * t509 + t500 * t650 - t503 * t648;
t498 = mrSges(2,1) * g(3) + mrSges(2,3) * t637 + t659 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t672 + qJ(2) * t669 + t648 * t500 + t650 * t503;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t676; (-m(1) - m(2)) * g(3) + t672; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t676 - t654 * t498 + t658 * t499; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t670 + t658 * t498 + t654 * t499; pkin(1) * t509 - mrSges(2,2) * t637 + mrSges(2,1) * t636 + t653 * t502 + t657 * t501 + pkin(2) * t660 + pkin(7) * t668 + mrSges(3,1) * t607 - mrSges(3,2) * t608 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
