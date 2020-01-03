% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:14
% EndTime: 2019-12-31 18:19:18
% DurationCPUTime: 4.25s
% Computational Cost: add. (47664->270), mult. (101010->344), div. (0->0), fcn. (62811->10), ass. (0->111)
t642 = sin(pkin(9));
t644 = cos(pkin(9));
t647 = sin(qJ(3));
t650 = cos(qJ(3));
t612 = (t642 * t647 - t644 * t650) * qJD(1);
t648 = sin(qJ(1));
t651 = cos(qJ(1));
t632 = t648 * g(1) - t651 * g(2);
t623 = qJDD(1) * pkin(1) + t632;
t633 = -t651 * g(1) - t648 * g(2);
t653 = qJD(1) ^ 2;
t625 = -t653 * pkin(1) + t633;
t643 = sin(pkin(8));
t645 = cos(pkin(8));
t603 = t643 * t623 + t645 * t625;
t595 = -t653 * pkin(2) + qJDD(1) * pkin(6) + t603;
t641 = -g(3) + qJDD(2);
t586 = -t647 * t595 + t650 * t641;
t667 = qJD(1) * qJD(3);
t666 = t650 * t667;
t626 = t647 * qJDD(1) + t666;
t575 = (-t626 + t666) * qJ(4) + (t647 * t650 * t653 + qJDD(3)) * pkin(3) + t586;
t587 = t650 * t595 + t647 * t641;
t627 = t650 * qJDD(1) - t647 * t667;
t669 = qJD(1) * t647;
t629 = qJD(3) * pkin(3) - qJ(4) * t669;
t640 = t650 ^ 2;
t576 = -t640 * t653 * pkin(3) + t627 * qJ(4) - qJD(3) * t629 + t587;
t672 = 2 * qJD(4);
t571 = t642 * t575 + t644 * t576 - t612 * t672;
t613 = (t642 * t650 + t644 * t647) * qJD(1);
t598 = t612 * pkin(4) - t613 * pkin(7);
t652 = qJD(3) ^ 2;
t569 = -t652 * pkin(4) + qJDD(3) * pkin(7) - t612 * t598 + t571;
t602 = t645 * t623 - t643 * t625;
t658 = -qJDD(1) * pkin(2) - t602;
t577 = -t627 * pkin(3) + qJDD(4) + t629 * t669 + (-qJ(4) * t640 - pkin(6)) * t653 + t658;
t604 = -t642 * t626 + t644 * t627;
t605 = t644 * t626 + t642 * t627;
t572 = (qJD(3) * t612 - t605) * pkin(7) + (qJD(3) * t613 - t604) * pkin(4) + t577;
t646 = sin(qJ(5));
t649 = cos(qJ(5));
t566 = -t646 * t569 + t649 * t572;
t606 = t649 * qJD(3) - t646 * t613;
t584 = t606 * qJD(5) + t646 * qJDD(3) + t649 * t605;
t607 = t646 * qJD(3) + t649 * t613;
t585 = -t606 * mrSges(6,1) + t607 * mrSges(6,2);
t611 = qJD(5) + t612;
t588 = -t611 * mrSges(6,2) + t606 * mrSges(6,3);
t601 = qJDD(5) - t604;
t563 = m(6) * t566 + t601 * mrSges(6,1) - t584 * mrSges(6,3) - t607 * t585 + t611 * t588;
t567 = t649 * t569 + t646 * t572;
t583 = -t607 * qJD(5) + t649 * qJDD(3) - t646 * t605;
t589 = t611 * mrSges(6,1) - t607 * mrSges(6,3);
t564 = m(6) * t567 - t601 * mrSges(6,2) + t583 * mrSges(6,3) + t606 * t585 - t611 * t589;
t555 = -t646 * t563 + t649 * t564;
t597 = t612 * mrSges(5,1) + t613 * mrSges(5,2);
t609 = qJD(3) * mrSges(5,1) - t613 * mrSges(5,3);
t552 = m(5) * t571 - qJDD(3) * mrSges(5,2) + t604 * mrSges(5,3) - qJD(3) * t609 - t612 * t597 + t555;
t661 = -t644 * t575 + t642 * t576;
t568 = -qJDD(3) * pkin(4) - t652 * pkin(7) + (t672 + t598) * t613 + t661;
t565 = -m(6) * t568 + t583 * mrSges(6,1) - t584 * mrSges(6,2) + t606 * t588 - t607 * t589;
t570 = -0.2e1 * qJD(4) * t613 - t661;
t608 = -qJD(3) * mrSges(5,2) - t612 * mrSges(5,3);
t559 = m(5) * t570 + qJDD(3) * mrSges(5,1) - t605 * mrSges(5,3) + qJD(3) * t608 - t613 * t597 + t565;
t546 = t642 * t552 + t644 * t559;
t578 = Ifges(6,5) * t607 + Ifges(6,6) * t606 + Ifges(6,3) * t611;
t580 = Ifges(6,1) * t607 + Ifges(6,4) * t606 + Ifges(6,5) * t611;
t556 = -mrSges(6,1) * t568 + mrSges(6,3) * t567 + Ifges(6,4) * t584 + Ifges(6,2) * t583 + Ifges(6,6) * t601 - t607 * t578 + t611 * t580;
t579 = Ifges(6,4) * t607 + Ifges(6,2) * t606 + Ifges(6,6) * t611;
t557 = mrSges(6,2) * t568 - mrSges(6,3) * t566 + Ifges(6,1) * t584 + Ifges(6,4) * t583 + Ifges(6,5) * t601 + t606 * t578 - t611 * t579;
t592 = Ifges(5,4) * t613 - Ifges(5,2) * t612 + Ifges(5,6) * qJD(3);
t593 = Ifges(5,1) * t613 - Ifges(5,4) * t612 + Ifges(5,5) * qJD(3);
t618 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t647 + Ifges(4,2) * t650) * qJD(1);
t619 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t647 + Ifges(4,4) * t650) * qJD(1);
t673 = (t647 * t618 - t650 * t619) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t586 + mrSges(5,1) * t570 - mrSges(4,2) * t587 - mrSges(5,2) * t571 + Ifges(4,5) * t626 + Ifges(5,5) * t605 + Ifges(4,6) * t627 + Ifges(5,6) * t604 + pkin(3) * t546 + pkin(4) * t565 + pkin(7) * t555 + t649 * t556 + t646 * t557 + t613 * t592 + t612 * t593;
t624 = (-mrSges(4,1) * t650 + mrSges(4,2) * t647) * qJD(1);
t668 = qJD(1) * t650;
t631 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t668;
t544 = m(4) * t586 + qJDD(3) * mrSges(4,1) - t626 * mrSges(4,3) + qJD(3) * t631 - t624 * t669 + t546;
t630 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t669;
t662 = t644 * t552 - t642 * t559;
t545 = m(4) * t587 - qJDD(3) * mrSges(4,2) + t627 * mrSges(4,3) - qJD(3) * t630 + t624 * t668 + t662;
t663 = -t647 * t544 + t650 * t545;
t535 = m(3) * t603 - t653 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t663;
t554 = t649 * t563 + t646 * t564;
t553 = m(5) * t577 - t604 * mrSges(5,1) + t605 * mrSges(5,2) + t612 * t608 + t613 * t609 + t554;
t594 = -t653 * pkin(6) + t658;
t655 = -m(4) * t594 + t627 * mrSges(4,1) - t626 * mrSges(4,2) - t630 * t669 + t631 * t668 - t553;
t548 = m(3) * t602 + qJDD(1) * mrSges(3,1) - t653 * mrSges(3,2) + t655;
t532 = t643 * t535 + t645 * t548;
t529 = m(2) * t632 + qJDD(1) * mrSges(2,1) - t653 * mrSges(2,2) + t532;
t664 = t645 * t535 - t643 * t548;
t530 = m(2) * t633 - t653 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t664;
t670 = t651 * t529 + t648 * t530;
t538 = t650 * t544 + t647 * t545;
t536 = m(3) * t641 + t538;
t665 = -t648 * t529 + t651 * t530;
t591 = Ifges(5,5) * t613 - Ifges(5,6) * t612 + Ifges(5,3) * qJD(3);
t539 = mrSges(5,2) * t577 - mrSges(5,3) * t570 + Ifges(5,1) * t605 + Ifges(5,4) * t604 + Ifges(5,5) * qJDD(3) - pkin(7) * t554 - qJD(3) * t592 - t646 * t556 + t649 * t557 - t612 * t591;
t656 = mrSges(6,1) * t566 - mrSges(6,2) * t567 + Ifges(6,5) * t584 + Ifges(6,6) * t583 + Ifges(6,3) * t601 + t607 * t579 - t606 * t580;
t540 = -mrSges(5,1) * t577 + mrSges(5,3) * t571 + Ifges(5,4) * t605 + Ifges(5,2) * t604 + Ifges(5,6) * qJDD(3) - pkin(4) * t554 + qJD(3) * t593 - t613 * t591 - t656;
t617 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t647 + Ifges(4,6) * t650) * qJD(1);
t523 = -mrSges(4,1) * t594 + mrSges(4,3) * t587 + Ifges(4,4) * t626 + Ifges(4,2) * t627 + Ifges(4,6) * qJDD(3) - pkin(3) * t553 + qJ(4) * t662 + qJD(3) * t619 + t642 * t539 + t644 * t540 - t617 * t669;
t525 = mrSges(4,2) * t594 - mrSges(4,3) * t586 + Ifges(4,1) * t626 + Ifges(4,4) * t627 + Ifges(4,5) * qJDD(3) - qJ(4) * t546 - qJD(3) * t618 + t644 * t539 - t642 * t540 + t617 * t668;
t657 = mrSges(2,1) * t632 + mrSges(3,1) * t602 - mrSges(2,2) * t633 - mrSges(3,2) * t603 + pkin(1) * t532 + pkin(2) * t655 + pkin(6) * t663 + t650 * t523 + t647 * t525 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t521 = -mrSges(3,1) * t641 + mrSges(3,3) * t603 + t653 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t538 - t673;
t520 = mrSges(3,2) * t641 - mrSges(3,3) * t602 + Ifges(3,5) * qJDD(1) - t653 * Ifges(3,6) - pkin(6) * t538 - t647 * t523 + t650 * t525;
t519 = -mrSges(2,2) * g(3) - mrSges(2,3) * t632 + Ifges(2,5) * qJDD(1) - t653 * Ifges(2,6) - qJ(2) * t532 + t645 * t520 - t643 * t521;
t518 = mrSges(2,1) * g(3) + mrSges(2,3) * t633 + t653 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t536 + qJ(2) * t664 + t643 * t520 + t645 * t521;
t1 = [-m(1) * g(1) + t665; -m(1) * g(2) + t670; (-m(1) - m(2)) * g(3) + t536; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t670 - t648 * t518 + t651 * t519; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t665 + t651 * t518 + t648 * t519; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t657; t657; t536; t673; t553; t656;];
tauJB = t1;
