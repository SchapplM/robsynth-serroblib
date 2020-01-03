% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:23
% EndTime: 2019-12-31 18:42:26
% DurationCPUTime: 2.48s
% Computational Cost: add. (25201->247), mult. (47601->296), div. (0->0), fcn. (27269->8), ass. (0->103)
t668 = Ifges(5,1) + Ifges(6,1);
t662 = Ifges(5,4) + Ifges(6,4);
t661 = Ifges(5,5) + Ifges(6,5);
t667 = Ifges(5,2) + Ifges(6,2);
t660 = Ifges(5,6) + Ifges(6,6);
t666 = Ifges(5,3) + Ifges(6,3);
t630 = sin(qJ(4));
t633 = cos(qJ(4));
t631 = sin(qJ(3));
t653 = qJD(1) * t631;
t608 = qJD(3) * t633 - t630 * t653;
t634 = cos(qJ(3));
t651 = qJD(1) * qJD(3);
t647 = t634 * t651;
t614 = qJDD(1) * t631 + t647;
t584 = qJD(4) * t608 + qJDD(3) * t630 + t614 * t633;
t609 = qJD(3) * t630 + t633 * t653;
t652 = qJD(1) * t634;
t621 = qJD(4) - t652;
t594 = mrSges(6,1) * t621 - mrSges(6,3) * t609;
t632 = sin(qJ(1));
t635 = cos(qJ(1));
t619 = t632 * g(1) - g(2) * t635;
t610 = qJDD(1) * pkin(1) + t619;
t620 = -g(1) * t635 - g(2) * t632;
t637 = qJD(1) ^ 2;
t612 = -pkin(1) * t637 + t620;
t628 = sin(pkin(8));
t629 = cos(pkin(8));
t587 = t628 * t610 + t629 * t612;
t571 = -pkin(2) * t637 + qJDD(1) * pkin(6) + t587;
t627 = -g(3) + qJDD(2);
t565 = -t631 * t571 + t627 * t634;
t613 = (-pkin(3) * t634 - pkin(7) * t631) * qJD(1);
t636 = qJD(3) ^ 2;
t563 = -qJDD(3) * pkin(3) - pkin(7) * t636 + t613 * t653 - t565;
t583 = -qJD(4) * t609 + qJDD(3) * t633 - t614 * t630;
t593 = pkin(4) * t621 - qJ(5) * t609;
t606 = t608 ^ 2;
t556 = -t583 * pkin(4) - qJ(5) * t606 + t593 * t609 + qJDD(5) + t563;
t591 = -mrSges(6,2) * t621 + mrSges(6,3) * t608;
t643 = -m(6) * t556 + t583 * mrSges(6,1) + t608 * t591;
t551 = mrSges(6,2) * t584 + t594 * t609 - t643;
t586 = t629 * t610 - t628 * t612;
t570 = -qJDD(1) * pkin(2) - t637 * pkin(6) - t586;
t648 = t631 * t651;
t615 = qJDD(1) * t634 - t648;
t561 = (-t614 - t647) * pkin(7) + (-t615 + t648) * pkin(3) + t570;
t566 = t634 * t571 + t631 * t627;
t564 = -pkin(3) * t636 + qJDD(3) * pkin(7) + t613 * t652 + t566;
t558 = t630 * t561 + t633 * t564;
t555 = -pkin(4) * t606 + t583 * qJ(5) + 0.2e1 * qJD(5) * t608 - t593 * t621 + t558;
t607 = qJDD(4) - t615;
t589 = -mrSges(6,1) * t608 + mrSges(6,2) * t609;
t649 = m(6) * t555 + t583 * mrSges(6,3) + t608 * t589;
t655 = t662 * t608 + t668 * t609 + t661 * t621;
t657 = -t660 * t608 - t661 * t609 - t666 * t621;
t534 = -mrSges(5,1) * t563 + mrSges(5,3) * t558 - mrSges(6,1) * t556 + mrSges(6,3) * t555 - pkin(4) * t551 + qJ(5) * t649 + (-qJ(5) * t594 + t655) * t621 + t657 * t609 + (-mrSges(6,2) * qJ(5) + t660) * t607 + t662 * t584 + t667 * t583;
t557 = t633 * t561 - t630 * t564;
t553 = -0.2e1 * qJD(5) * t609 + (t608 * t621 - t584) * qJ(5) + (t608 * t609 + t607) * pkin(4) + t557;
t650 = m(6) * t553 + t607 * mrSges(6,1) + t621 * t591;
t550 = -mrSges(6,3) * t584 - t589 * t609 + t650;
t656 = -t667 * t608 - t662 * t609 - t660 * t621;
t541 = mrSges(5,2) * t563 + mrSges(6,2) * t556 - mrSges(5,3) * t557 - mrSges(6,3) * t553 - qJ(5) * t550 + t662 * t583 + t668 * t584 + t661 * t607 - t657 * t608 + t656 * t621;
t590 = -mrSges(5,1) * t608 + mrSges(5,2) * t609;
t592 = -mrSges(5,2) * t621 + mrSges(5,3) * t608;
t544 = m(5) * t557 + mrSges(5,1) * t607 + t592 * t621 + (-t589 - t590) * t609 + (-mrSges(5,3) - mrSges(6,3)) * t584 + t650;
t654 = -mrSges(5,1) * t621 + mrSges(5,3) * t609 - t594;
t663 = -mrSges(5,2) - mrSges(6,2);
t546 = m(5) * t558 + t583 * mrSges(5,3) + t590 * t608 + t663 * t607 + t654 * t621 + t649;
t543 = -t544 * t630 + t633 * t546;
t549 = -m(5) * t563 + t583 * mrSges(5,1) + t663 * t584 + t608 * t592 + t654 * t609 + t643;
t601 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t631 + Ifges(4,2) * t634) * qJD(1);
t602 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t631 + Ifges(4,4) * t634) * qJD(1);
t665 = mrSges(4,1) * t565 - mrSges(4,2) * t566 + Ifges(4,5) * t614 + Ifges(4,6) * t615 + Ifges(4,3) * qJDD(3) + pkin(3) * t549 + pkin(7) * t543 + t633 * t534 + t630 * t541 + (t601 * t631 - t602 * t634) * qJD(1);
t664 = mrSges(5,1) * t557 + mrSges(6,1) * t553 - mrSges(5,2) * t558 - mrSges(6,2) * t555 + pkin(4) * t550 + t660 * t583 + t661 * t584 + t666 * t607 - t655 * t608 - t656 * t609;
t611 = (-mrSges(4,1) * t634 + mrSges(4,2) * t631) * qJD(1);
t617 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t653;
t540 = m(4) * t566 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t615 - qJD(3) * t617 + t611 * t652 + t543;
t618 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t652;
t548 = m(4) * t565 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t614 + qJD(3) * t618 - t611 * t653 + t549;
t644 = t634 * t540 - t548 * t631;
t530 = m(3) * t587 - mrSges(3,1) * t637 - qJDD(1) * mrSges(3,2) + t644;
t542 = t544 * t633 + t546 * t630;
t639 = -m(4) * t570 + t615 * mrSges(4,1) - mrSges(4,2) * t614 - t617 * t653 + t618 * t652 - t542;
t536 = m(3) * t586 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t637 + t639;
t525 = t628 * t530 + t629 * t536;
t522 = m(2) * t619 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t637 + t525;
t645 = t629 * t530 - t628 * t536;
t523 = m(2) * t620 - mrSges(2,1) * t637 - qJDD(1) * mrSges(2,2) + t645;
t658 = t635 * t522 + t632 * t523;
t533 = t631 * t540 + t634 * t548;
t531 = m(3) * t627 + t533;
t646 = -t522 * t632 + t635 * t523;
t600 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t631 + Ifges(4,6) * t634) * qJD(1);
t518 = mrSges(4,2) * t570 - mrSges(4,3) * t565 + Ifges(4,1) * t614 + Ifges(4,4) * t615 + Ifges(4,5) * qJDD(3) - pkin(7) * t542 - qJD(3) * t601 - t534 * t630 + t541 * t633 + t600 * t652;
t527 = -mrSges(4,1) * t570 + mrSges(4,3) * t566 + Ifges(4,4) * t614 + Ifges(4,2) * t615 + Ifges(4,6) * qJDD(3) - pkin(3) * t542 + qJD(3) * t602 - t600 * t653 - t664;
t640 = mrSges(2,1) * t619 + mrSges(3,1) * t586 - mrSges(2,2) * t620 - mrSges(3,2) * t587 + pkin(1) * t525 + pkin(2) * t639 + pkin(6) * t644 + t631 * t518 + t634 * t527 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t516 = -mrSges(3,1) * t627 + mrSges(3,3) * t587 + t637 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t533 - t665;
t515 = mrSges(3,2) * t627 - mrSges(3,3) * t586 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t637 - pkin(6) * t533 + t518 * t634 - t527 * t631;
t514 = -mrSges(2,2) * g(3) - mrSges(2,3) * t619 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t637 - qJ(2) * t525 + t515 * t629 - t516 * t628;
t513 = mrSges(2,1) * g(3) + mrSges(2,3) * t620 + t637 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t531 + qJ(2) * t645 + t628 * t515 + t629 * t516;
t1 = [-m(1) * g(1) + t646; -m(1) * g(2) + t658; (-m(1) - m(2)) * g(3) + t531; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t658 - t632 * t513 + t635 * t514; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t646 + t635 * t513 + t632 * t514; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t640; t640; t531; t665; t664; t551;];
tauJB = t1;
