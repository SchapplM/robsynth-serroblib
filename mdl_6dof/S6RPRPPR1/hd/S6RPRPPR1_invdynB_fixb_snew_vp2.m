% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:24:22
% EndTime: 2019-05-05 16:24:34
% DurationCPUTime: 11.92s
% Computational Cost: add. (178880->342), mult. (399735->433), div. (0->0), fcn. (268609->12), ass. (0->133)
t679 = -2 * qJD(4);
t653 = sin(qJ(1));
t656 = cos(qJ(1));
t637 = g(1) * t653 - g(2) * t656;
t629 = qJDD(1) * pkin(1) + t637;
t638 = -g(1) * t656 - g(2) * t653;
t658 = qJD(1) ^ 2;
t631 = -pkin(1) * t658 + t638;
t648 = sin(pkin(9));
t650 = cos(pkin(9));
t604 = t629 * t648 + t631 * t650;
t597 = -pkin(2) * t658 + qJDD(1) * pkin(7) + t604;
t645 = -g(3) + qJDD(2);
t652 = sin(qJ(3));
t655 = cos(qJ(3));
t584 = -t652 * t597 + t645 * t655;
t673 = qJD(1) * qJD(3);
t671 = t655 * t673;
t632 = qJDD(1) * t652 + t671;
t571 = (-t632 + t671) * qJ(4) + (t652 * t655 * t658 + qJDD(3)) * pkin(3) + t584;
t585 = t597 * t655 + t645 * t652;
t633 = qJDD(1) * t655 - t652 * t673;
t676 = qJD(1) * t652;
t634 = qJD(3) * pkin(3) - qJ(4) * t676;
t644 = t655 ^ 2;
t574 = -pkin(3) * t644 * t658 + qJ(4) * t633 - qJD(3) * t634 + t585;
t647 = sin(pkin(10));
t678 = cos(pkin(10));
t618 = (t647 * t655 + t652 * t678) * qJD(1);
t555 = t571 * t678 - t574 * t647 + t618 * t679;
t675 = qJD(1) * t655;
t617 = t647 * t676 - t675 * t678;
t556 = t571 * t647 + t574 * t678 + t617 * t679;
t600 = mrSges(5,1) * t617 + mrSges(5,2) * t618;
t605 = t632 * t647 - t633 * t678;
t613 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t618;
t599 = pkin(4) * t617 - qJ(5) * t618;
t657 = qJD(3) ^ 2;
t554 = -pkin(4) * t657 + qJDD(3) * qJ(5) - t599 * t617 + t556;
t603 = t650 * t629 - t631 * t648;
t663 = -qJDD(1) * pkin(2) - t603;
t575 = -t633 * pkin(3) + qJDD(4) + t634 * t676 + (-qJ(4) * t644 - pkin(7)) * t658 + t663;
t606 = t632 * t678 + t633 * t647;
t559 = (qJD(3) * t617 - t606) * qJ(5) + (qJD(3) * t618 + t605) * pkin(4) + t575;
t646 = sin(pkin(11));
t649 = cos(pkin(11));
t611 = qJD(3) * t646 + t618 * t649;
t549 = -0.2e1 * qJD(5) * t611 - t646 * t554 + t559 * t649;
t592 = qJDD(3) * t646 + t606 * t649;
t610 = qJD(3) * t649 - t618 * t646;
t547 = (t610 * t617 - t592) * pkin(8) + (t610 * t611 + t605) * pkin(5) + t549;
t550 = 0.2e1 * qJD(5) * t610 + t554 * t649 + t559 * t646;
t588 = pkin(5) * t617 - pkin(8) * t611;
t591 = qJDD(3) * t649 - t606 * t646;
t609 = t610 ^ 2;
t548 = -pkin(5) * t609 + pkin(8) * t591 - t588 * t617 + t550;
t651 = sin(qJ(6));
t654 = cos(qJ(6));
t545 = t547 * t654 - t548 * t651;
t580 = t610 * t654 - t611 * t651;
t562 = qJD(6) * t580 + t591 * t651 + t592 * t654;
t581 = t610 * t651 + t611 * t654;
t567 = -mrSges(7,1) * t580 + mrSges(7,2) * t581;
t616 = qJD(6) + t617;
t572 = -mrSges(7,2) * t616 + mrSges(7,3) * t580;
t602 = qJDD(6) + t605;
t543 = m(7) * t545 + mrSges(7,1) * t602 - mrSges(7,3) * t562 - t567 * t581 + t572 * t616;
t546 = t547 * t651 + t548 * t654;
t561 = -qJD(6) * t581 + t591 * t654 - t592 * t651;
t573 = mrSges(7,1) * t616 - mrSges(7,3) * t581;
t544 = m(7) * t546 - mrSges(7,2) * t602 + mrSges(7,3) * t561 + t567 * t580 - t573 * t616;
t535 = t543 * t654 + t544 * t651;
t582 = -mrSges(6,1) * t610 + mrSges(6,2) * t611;
t586 = -mrSges(6,2) * t617 + mrSges(6,3) * t610;
t533 = m(6) * t549 + mrSges(6,1) * t605 - mrSges(6,3) * t592 - t582 * t611 + t586 * t617 + t535;
t587 = mrSges(6,1) * t617 - mrSges(6,3) * t611;
t665 = -t543 * t651 + t544 * t654;
t534 = m(6) * t550 - mrSges(6,2) * t605 + mrSges(6,3) * t591 + t582 * t610 - t587 * t617 + t665;
t666 = -t533 * t646 + t534 * t649;
t528 = m(5) * t556 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t605 - qJD(3) * t613 - t600 * t617 + t666;
t612 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t617;
t553 = -qJDD(3) * pkin(4) - qJ(5) * t657 + t599 * t618 + qJDD(5) - t555;
t551 = -pkin(5) * t591 - pkin(8) * t609 + t588 * t611 + t553;
t662 = m(7) * t551 - mrSges(7,1) * t561 + mrSges(7,2) * t562 - t572 * t580 + t573 * t581;
t660 = -m(6) * t553 + mrSges(6,1) * t591 - mrSges(6,2) * t592 + t586 * t610 - t587 * t611 - t662;
t539 = m(5) * t555 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t606 + qJD(3) * t612 - t600 * t618 + t660;
t521 = t528 * t647 + t539 * t678;
t630 = (-mrSges(4,1) * t655 + mrSges(4,2) * t652) * qJD(1);
t636 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t675;
t519 = m(4) * t584 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t632 + qJD(3) * t636 - t630 * t676 + t521;
t635 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t676;
t667 = t528 * t678 - t539 * t647;
t520 = m(4) * t585 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t633 - qJD(3) * t635 + t630 * t675 + t667;
t668 = -t519 * t652 + t520 * t655;
t513 = m(3) * t604 - mrSges(3,1) * t658 - qJDD(1) * mrSges(3,2) + t668;
t596 = -t658 * pkin(7) + t663;
t529 = t533 * t649 + t534 * t646;
t661 = m(5) * t575 + mrSges(5,1) * t605 + mrSges(5,2) * t606 + t612 * t617 + t613 * t618 + t529;
t659 = -m(4) * t596 + mrSges(4,1) * t633 - mrSges(4,2) * t632 - t635 * t676 + t636 * t675 - t661;
t525 = m(3) * t603 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t658 + t659;
t509 = t513 * t648 + t525 * t650;
t507 = m(2) * t637 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t658 + t509;
t669 = t513 * t650 - t525 * t648;
t508 = m(2) * t638 - mrSges(2,1) * t658 - qJDD(1) * mrSges(2,2) + t669;
t677 = t507 * t656 + t508 * t653;
t514 = t519 * t655 + t520 * t652;
t672 = m(3) * t645 + t514;
t670 = -t507 * t653 + t508 * t656;
t624 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t652 + Ifges(4,4) * t655) * qJD(1);
t623 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t652 + Ifges(4,2) * t655) * qJD(1);
t622 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t652 + Ifges(4,6) * t655) * qJD(1);
t595 = Ifges(5,1) * t618 - Ifges(5,4) * t617 + Ifges(5,5) * qJD(3);
t594 = Ifges(5,4) * t618 - Ifges(5,2) * t617 + Ifges(5,6) * qJD(3);
t593 = Ifges(5,5) * t618 - Ifges(5,6) * t617 + Ifges(5,3) * qJD(3);
t578 = Ifges(6,1) * t611 + Ifges(6,4) * t610 + Ifges(6,5) * t617;
t577 = Ifges(6,4) * t611 + Ifges(6,2) * t610 + Ifges(6,6) * t617;
t576 = Ifges(6,5) * t611 + Ifges(6,6) * t610 + Ifges(6,3) * t617;
t565 = Ifges(7,1) * t581 + Ifges(7,4) * t580 + Ifges(7,5) * t616;
t564 = Ifges(7,4) * t581 + Ifges(7,2) * t580 + Ifges(7,6) * t616;
t563 = Ifges(7,5) * t581 + Ifges(7,6) * t580 + Ifges(7,3) * t616;
t537 = mrSges(7,2) * t551 - mrSges(7,3) * t545 + Ifges(7,1) * t562 + Ifges(7,4) * t561 + Ifges(7,5) * t602 + t563 * t580 - t564 * t616;
t536 = -mrSges(7,1) * t551 + mrSges(7,3) * t546 + Ifges(7,4) * t562 + Ifges(7,2) * t561 + Ifges(7,6) * t602 - t563 * t581 + t565 * t616;
t523 = mrSges(6,2) * t553 - mrSges(6,3) * t549 + Ifges(6,1) * t592 + Ifges(6,4) * t591 + Ifges(6,5) * t605 - pkin(8) * t535 - t536 * t651 + t537 * t654 + t576 * t610 - t577 * t617;
t522 = -mrSges(6,1) * t553 + mrSges(6,3) * t550 + Ifges(6,4) * t592 + Ifges(6,2) * t591 + Ifges(6,6) * t605 - pkin(5) * t662 + pkin(8) * t665 + t654 * t536 + t651 * t537 - t611 * t576 + t617 * t578;
t515 = Ifges(5,4) * t606 + Ifges(5,6) * qJDD(3) - t618 * t593 + qJD(3) * t595 - mrSges(5,1) * t575 + mrSges(5,3) * t556 - Ifges(6,5) * t592 - Ifges(6,6) * t591 - t611 * t577 + t610 * t578 - mrSges(6,1) * t549 + mrSges(6,2) * t550 - Ifges(7,5) * t562 - Ifges(7,6) * t561 - Ifges(7,3) * t602 - t581 * t564 + t580 * t565 - mrSges(7,1) * t545 + mrSges(7,2) * t546 - pkin(5) * t535 - pkin(4) * t529 + (-Ifges(5,2) - Ifges(6,3)) * t605;
t510 = mrSges(5,2) * t575 - mrSges(5,3) * t555 + Ifges(5,1) * t606 - Ifges(5,4) * t605 + Ifges(5,5) * qJDD(3) - qJ(5) * t529 - qJD(3) * t594 - t522 * t646 + t523 * t649 - t593 * t617;
t503 = mrSges(4,2) * t596 - mrSges(4,3) * t584 + Ifges(4,1) * t632 + Ifges(4,4) * t633 + Ifges(4,5) * qJDD(3) - qJ(4) * t521 - qJD(3) * t623 + t510 * t678 - t515 * t647 + t622 * t675;
t502 = -pkin(2) * t514 - mrSges(3,1) * t645 + mrSges(3,3) * t604 - pkin(3) * t521 - Ifges(4,5) * t632 - Ifges(4,6) * t633 - mrSges(4,1) * t584 + mrSges(4,2) * t585 - t646 * t523 - t649 * t522 - pkin(4) * t660 - qJ(5) * t666 - Ifges(5,5) * t606 + Ifges(5,6) * t605 - mrSges(5,1) * t555 + mrSges(5,2) * t556 + t658 * Ifges(3,5) - t618 * t594 - t617 * t595 + Ifges(3,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t623 * t652 + t624 * t655) * qJD(1);
t501 = -mrSges(4,1) * t596 + mrSges(4,3) * t585 + Ifges(4,4) * t632 + Ifges(4,2) * t633 + Ifges(4,6) * qJDD(3) - pkin(3) * t661 + qJ(4) * t667 + qJD(3) * t624 + t647 * t510 + t515 * t678 - t622 * t676;
t500 = mrSges(3,2) * t645 - mrSges(3,3) * t603 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t658 - pkin(7) * t514 - t501 * t652 + t503 * t655;
t499 = -mrSges(2,2) * g(3) - mrSges(2,3) * t637 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t658 - qJ(2) * t509 + t500 * t650 - t502 * t648;
t498 = mrSges(2,1) * g(3) + mrSges(2,3) * t638 + t658 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t672 + qJ(2) * t669 + t648 * t500 + t650 * t502;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t677; (-m(1) - m(2)) * g(3) + t672; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t677 - t498 * t653 + t499 * t656; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t670 + t656 * t498 + t653 * t499; pkin(1) * t509 + mrSges(2,1) * t637 - mrSges(2,2) * t638 + t655 * t501 + pkin(2) * t659 + pkin(7) * t668 + t652 * t503 + mrSges(3,1) * t603 - mrSges(3,2) * t604 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
