% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:47
% EndTime: 2019-12-31 18:44:50
% DurationCPUTime: 2.41s
% Computational Cost: add. (24589->245), mult. (46031->296), div. (0->0), fcn. (26181->8), ass. (0->102)
t666 = Ifges(5,1) + Ifges(6,1);
t659 = Ifges(5,4) - Ifges(6,5);
t658 = -Ifges(5,5) - Ifges(6,4);
t665 = Ifges(5,2) + Ifges(6,3);
t657 = Ifges(5,6) - Ifges(6,6);
t664 = -Ifges(5,3) - Ifges(6,2);
t631 = sin(qJ(1));
t633 = cos(qJ(1));
t615 = t631 * g(1) - t633 * g(2);
t606 = qJDD(1) * pkin(1) + t615;
t616 = -t633 * g(1) - t631 * g(2);
t635 = qJD(1) ^ 2;
t608 = -t635 * pkin(1) + t616;
t627 = sin(pkin(8));
t628 = cos(pkin(8));
t582 = t627 * t606 + t628 * t608;
t567 = -t635 * pkin(2) + qJDD(1) * pkin(6) + t582;
t626 = -g(3) + qJDD(2);
t630 = sin(qJ(3));
t632 = cos(qJ(3));
t562 = -t630 * t567 + t632 * t626;
t609 = (-pkin(3) * t632 - pkin(7) * t630) * qJD(1);
t634 = qJD(3) ^ 2;
t650 = qJD(1) * t630;
t560 = -qJDD(3) * pkin(3) - t634 * pkin(7) + t609 * t650 - t562;
t629 = sin(qJ(4));
t661 = cos(qJ(4));
t605 = t629 * qJD(3) + t661 * t650;
t648 = qJD(1) * qJD(3);
t645 = t632 * t648;
t610 = t630 * qJDD(1) + t645;
t578 = t605 * qJD(4) - t661 * qJDD(3) + t629 * t610;
t604 = -t661 * qJD(3) + t629 * t650;
t579 = -t604 * qJD(4) + t629 * qJDD(3) + t661 * t610;
t649 = t632 * qJD(1);
t618 = qJD(4) - t649;
t554 = -0.2e1 * qJD(5) * t605 + (t604 * t618 - t579) * qJ(5) + (t605 * t618 + t578) * pkin(4) + t560;
t590 = -t618 * mrSges(6,1) + t605 * mrSges(6,2);
t591 = -t604 * mrSges(6,2) + t618 * mrSges(6,3);
t550 = m(6) * t554 + t578 * mrSges(6,1) - t579 * mrSges(6,3) - t605 * t590 + t604 * t591;
t581 = t628 * t606 - t627 * t608;
t566 = -qJDD(1) * pkin(2) - t635 * pkin(6) - t581;
t646 = t630 * t648;
t611 = t632 * qJDD(1) - t646;
t558 = (-t610 - t645) * pkin(7) + (-t611 + t646) * pkin(3) + t566;
t563 = t632 * t567 + t630 * t626;
t561 = -t634 * pkin(3) + qJDD(3) * pkin(7) + t609 * t649 + t563;
t556 = t629 * t558 + t661 * t561;
t585 = t604 * pkin(4) - t605 * qJ(5);
t603 = qJDD(4) - t611;
t617 = t618 ^ 2;
t552 = -t617 * pkin(4) + t603 * qJ(5) + 0.2e1 * qJD(5) * t618 - t604 * t585 + t556;
t652 = t659 * t604 - t666 * t605 + t658 * t618;
t654 = t657 * t604 + t658 * t605 + t664 * t618;
t534 = -mrSges(5,1) * t560 - mrSges(6,1) * t554 + mrSges(6,2) * t552 + mrSges(5,3) * t556 - pkin(4) * t550 - t665 * t578 + t659 * t579 + t657 * t603 + t654 * t605 - t652 * t618;
t555 = t661 * t558 - t629 * t561;
t553 = -t603 * pkin(4) - t617 * qJ(5) + t605 * t585 + qJDD(5) - t555;
t653 = t665 * t604 - t659 * t605 - t657 * t618;
t535 = mrSges(5,2) * t560 + mrSges(6,2) * t553 - mrSges(5,3) * t555 - mrSges(6,3) * t554 - qJ(5) * t550 - t659 * t578 + t666 * t579 - t658 * t603 + t654 * t604 + t653 * t618;
t589 = t618 * mrSges(5,1) - t605 * mrSges(5,3);
t647 = m(6) * t552 + t603 * mrSges(6,3) + t618 * t590;
t586 = t604 * mrSges(6,1) - t605 * mrSges(6,3);
t651 = -t604 * mrSges(5,1) - t605 * mrSges(5,2) - t586;
t660 = -mrSges(5,3) - mrSges(6,2);
t545 = m(5) * t556 - t603 * mrSges(5,2) + t660 * t578 - t618 * t589 + t651 * t604 + t647;
t588 = -t618 * mrSges(5,2) - t604 * mrSges(5,3);
t641 = -m(6) * t553 + t603 * mrSges(6,1) + t618 * t591;
t546 = m(5) * t555 + t603 * mrSges(5,1) + t660 * t579 + t618 * t588 + t651 * t605 + t641;
t541 = t661 * t545 - t629 * t546;
t547 = -m(5) * t560 - t578 * mrSges(5,1) - t579 * mrSges(5,2) - t604 * t588 - t605 * t589 - t550;
t597 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t630 + Ifges(4,2) * t632) * qJD(1);
t598 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t630 + Ifges(4,4) * t632) * qJD(1);
t663 = mrSges(4,1) * t562 - mrSges(4,2) * t563 + Ifges(4,5) * t610 + Ifges(4,6) * t611 + Ifges(4,3) * qJDD(3) + pkin(3) * t547 + pkin(7) * t541 + (t597 * t630 - t598 * t632) * qJD(1) + t661 * t534 + t629 * t535;
t549 = t579 * mrSges(6,2) + t605 * t586 - t641;
t662 = -t657 * t578 - t658 * t579 - t664 * t603 - t652 * t604 - t653 * t605 + mrSges(5,1) * t555 - mrSges(6,1) * t553 - mrSges(5,2) * t556 + mrSges(6,3) * t552 - pkin(4) * t549 + qJ(5) * (-t578 * mrSges(6,2) - t604 * t586 + t647);
t607 = (-mrSges(4,1) * t632 + mrSges(4,2) * t630) * qJD(1);
t613 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t650;
t539 = m(4) * t563 - qJDD(3) * mrSges(4,2) + t611 * mrSges(4,3) - qJD(3) * t613 + t607 * t649 + t541;
t614 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t649;
t543 = m(4) * t562 + qJDD(3) * mrSges(4,1) - t610 * mrSges(4,3) + qJD(3) * t614 - t607 * t650 + t547;
t642 = t632 * t539 - t630 * t543;
t528 = m(3) * t582 - t635 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t642;
t540 = t629 * t545 + t661 * t546;
t637 = -m(4) * t566 + t611 * mrSges(4,1) - t610 * mrSges(4,2) - t613 * t650 + t614 * t649 - t540;
t533 = m(3) * t581 + qJDD(1) * mrSges(3,1) - t635 * mrSges(3,2) + t637;
t523 = t627 * t528 + t628 * t533;
t520 = m(2) * t615 + qJDD(1) * mrSges(2,1) - t635 * mrSges(2,2) + t523;
t643 = t628 * t528 - t627 * t533;
t521 = m(2) * t616 - t635 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t643;
t655 = t633 * t520 + t631 * t521;
t531 = t630 * t539 + t632 * t543;
t529 = m(3) * t626 + t531;
t644 = -t631 * t520 + t633 * t521;
t596 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t630 + Ifges(4,6) * t632) * qJD(1);
t516 = mrSges(4,2) * t566 - mrSges(4,3) * t562 + Ifges(4,1) * t610 + Ifges(4,4) * t611 + Ifges(4,5) * qJDD(3) - pkin(7) * t540 - qJD(3) * t597 - t629 * t534 + t661 * t535 + t596 * t649;
t525 = -mrSges(4,1) * t566 + mrSges(4,3) * t563 + Ifges(4,4) * t610 + Ifges(4,2) * t611 + Ifges(4,6) * qJDD(3) - pkin(3) * t540 + qJD(3) * t598 - t596 * t650 - t662;
t639 = mrSges(2,1) * t615 + mrSges(3,1) * t581 - mrSges(2,2) * t616 - mrSges(3,2) * t582 + pkin(1) * t523 + pkin(2) * t637 + pkin(6) * t642 + t630 * t516 + t632 * t525 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t514 = -mrSges(3,1) * t626 + mrSges(3,3) * t582 + t635 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t531 - t663;
t513 = mrSges(3,2) * t626 - mrSges(3,3) * t581 + Ifges(3,5) * qJDD(1) - t635 * Ifges(3,6) - pkin(6) * t531 + t632 * t516 - t630 * t525;
t512 = -mrSges(2,2) * g(3) - mrSges(2,3) * t615 + Ifges(2,5) * qJDD(1) - t635 * Ifges(2,6) - qJ(2) * t523 + t628 * t513 - t627 * t514;
t511 = mrSges(2,1) * g(3) + mrSges(2,3) * t616 + t635 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t529 + qJ(2) * t643 + t627 * t513 + t628 * t514;
t1 = [-m(1) * g(1) + t644; -m(1) * g(2) + t655; (-m(1) - m(2)) * g(3) + t529; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t655 - t631 * t511 + t633 * t512; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t644 + t633 * t511 + t631 * t512; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t639; t639; t529; t663; t662; t549;];
tauJB = t1;
