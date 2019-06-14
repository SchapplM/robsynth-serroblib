% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:18:10
% EndTime: 2019-05-06 01:18:21
% DurationCPUTime: 6.78s
% Computational Cost: add. (106167->319), mult. (204187->386), div. (0->0), fcn. (133422->10), ass. (0->126)
t669 = Ifges(6,1) + Ifges(7,1);
t664 = Ifges(6,4) - Ifges(7,5);
t663 = Ifges(7,4) + Ifges(6,5);
t668 = Ifges(6,2) + Ifges(7,3);
t662 = Ifges(6,6) - Ifges(7,6);
t667 = -Ifges(6,3) - Ifges(7,2);
t666 = cos(qJ(5));
t665 = -mrSges(6,3) - mrSges(7,2);
t636 = sin(qJ(1));
t639 = cos(qJ(1));
t621 = t636 * g(1) - g(2) * t639;
t613 = qJDD(1) * pkin(1) + t621;
t622 = -g(1) * t639 - g(2) * t636;
t641 = qJD(1) ^ 2;
t615 = -pkin(1) * t641 + t622;
t631 = sin(pkin(10));
t632 = cos(pkin(10));
t591 = t631 * t613 + t632 * t615;
t579 = -pkin(2) * t641 + qJDD(1) * pkin(7) + t591;
t630 = -g(3) + qJDD(2);
t635 = sin(qJ(3));
t638 = cos(qJ(3));
t570 = t638 * t579 + t635 * t630;
t614 = (-mrSges(4,1) * t638 + mrSges(4,2) * t635) * qJD(1);
t654 = qJD(1) * qJD(3);
t627 = t635 * t654;
t618 = qJDD(1) * t638 - t627;
t656 = qJD(1) * t635;
t619 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t656;
t590 = t632 * t613 - t631 * t615;
t578 = -qJDD(1) * pkin(2) - t641 * pkin(7) - t590;
t651 = t638 * t654;
t617 = qJDD(1) * t635 + t651;
t557 = (-t617 - t651) * pkin(8) + (-t618 + t627) * pkin(3) + t578;
t616 = (-pkin(3) * t638 - pkin(8) * t635) * qJD(1);
t640 = qJD(3) ^ 2;
t655 = qJD(1) * t638;
t563 = -pkin(3) * t640 + qJDD(3) * pkin(8) + t616 * t655 + t570;
t634 = sin(qJ(4));
t637 = cos(qJ(4));
t542 = t637 * t557 - t634 * t563;
t611 = qJD(3) * t637 - t634 * t656;
t587 = qJD(4) * t611 + qJDD(3) * t634 + t617 * t637;
t610 = qJDD(4) - t618;
t612 = qJD(3) * t634 + t637 * t656;
t625 = qJD(4) - t655;
t538 = (t611 * t625 - t587) * pkin(9) + (t611 * t612 + t610) * pkin(4) + t542;
t543 = t634 * t557 + t637 * t563;
t586 = -qJD(4) * t612 + qJDD(3) * t637 - t617 * t634;
t595 = pkin(4) * t625 - pkin(9) * t612;
t609 = t611 ^ 2;
t540 = -pkin(4) * t609 + pkin(9) * t586 - t595 * t625 + t543;
t633 = sin(qJ(5));
t536 = t633 * t538 + t666 * t540;
t589 = t633 * t611 + t666 * t612;
t547 = qJD(5) * t589 - t666 * t586 + t587 * t633;
t624 = qJD(5) + t625;
t573 = mrSges(6,1) * t624 - mrSges(6,3) * t589;
t588 = -t666 * t611 + t612 * t633;
t606 = qJDD(5) + t610;
t564 = pkin(5) * t588 - qJ(6) * t589;
t623 = t624 ^ 2;
t531 = -pkin(5) * t623 + qJ(6) * t606 + 0.2e1 * qJD(6) * t624 - t564 * t588 + t536;
t574 = -mrSges(7,1) * t624 + mrSges(7,2) * t589;
t653 = m(7) * t531 + t606 * mrSges(7,3) + t624 * t574;
t565 = mrSges(7,1) * t588 - mrSges(7,3) * t589;
t657 = -mrSges(6,1) * t588 - mrSges(6,2) * t589 - t565;
t526 = m(6) * t536 - t606 * mrSges(6,2) + t665 * t547 - t624 * t573 + t657 * t588 + t653;
t535 = t666 * t538 - t633 * t540;
t548 = -t588 * qJD(5) + t633 * t586 + t666 * t587;
t572 = -mrSges(6,2) * t624 - mrSges(6,3) * t588;
t532 = -t606 * pkin(5) - t623 * qJ(6) + t589 * t564 + qJDD(6) - t535;
t571 = -mrSges(7,2) * t588 + mrSges(7,3) * t624;
t645 = -m(7) * t532 + t606 * mrSges(7,1) + t624 * t571;
t528 = m(6) * t535 + t606 * mrSges(6,1) + t665 * t548 + t624 * t572 + t657 * t589 + t645;
t521 = t633 * t526 + t666 * t528;
t592 = -mrSges(5,1) * t611 + mrSges(5,2) * t612;
t593 = -mrSges(5,2) * t625 + mrSges(5,3) * t611;
t517 = m(5) * t542 + mrSges(5,1) * t610 - mrSges(5,3) * t587 - t592 * t612 + t593 * t625 + t521;
t594 = mrSges(5,1) * t625 - mrSges(5,3) * t612;
t646 = t666 * t526 - t528 * t633;
t518 = m(5) * t543 - mrSges(5,2) * t610 + mrSges(5,3) * t586 + t592 * t611 - t594 * t625 + t646;
t647 = -t517 * t634 + t637 * t518;
t514 = m(4) * t570 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t618 - qJD(3) * t619 + t614 * t655 + t647;
t569 = -t635 * t579 + t638 * t630;
t620 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t655;
t562 = -qJDD(3) * pkin(3) - t640 * pkin(8) + t616 * t656 - t569;
t541 = -t586 * pkin(4) - t609 * pkin(9) + t612 * t595 + t562;
t534 = -0.2e1 * qJD(6) * t589 + (t588 * t624 - t548) * qJ(6) + (t589 * t624 + t547) * pkin(5) + t541;
t529 = m(7) * t534 + t547 * mrSges(7,1) - t548 * mrSges(7,3) + t588 * t571 - t589 * t574;
t644 = m(6) * t541 + t547 * mrSges(6,1) + t548 * mrSges(6,2) + t588 * t572 + t589 * t573 + t529;
t642 = -m(5) * t562 + t586 * mrSges(5,1) - t587 * mrSges(5,2) + t611 * t593 - t612 * t594 - t644;
t523 = m(4) * t569 + qJDD(3) * mrSges(4,1) - t617 * mrSges(4,3) + qJD(3) * t620 - t614 * t656 + t642;
t648 = t638 * t514 - t523 * t635;
t508 = m(3) * t591 - mrSges(3,1) * t641 - qJDD(1) * mrSges(3,2) + t648;
t515 = t517 * t637 + t518 * t634;
t643 = -m(4) * t578 + t618 * mrSges(4,1) - mrSges(4,2) * t617 - t619 * t656 + t620 * t655 - t515;
t511 = m(3) * t590 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t641 + t643;
t502 = t631 * t508 + t632 * t511;
t500 = m(2) * t621 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t641 + t502;
t649 = t632 * t508 - t511 * t631;
t501 = m(2) * t622 - mrSges(2,1) * t641 - qJDD(1) * mrSges(2,2) + t649;
t661 = t639 * t500 + t636 * t501;
t509 = t635 * t514 + t638 * t523;
t660 = t588 * t668 - t589 * t664 - t624 * t662;
t659 = t588 * t662 - t589 * t663 + t624 * t667;
t658 = -t664 * t588 + t589 * t669 + t663 * t624;
t652 = m(3) * t630 + t509;
t650 = -t500 * t636 + t639 * t501;
t605 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t635 + Ifges(4,4) * t638) * qJD(1);
t604 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t635 + Ifges(4,2) * t638) * qJD(1);
t603 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t635 + Ifges(4,6) * t638) * qJD(1);
t582 = Ifges(5,1) * t612 + Ifges(5,4) * t611 + Ifges(5,5) * t625;
t581 = Ifges(5,4) * t612 + Ifges(5,2) * t611 + Ifges(5,6) * t625;
t580 = Ifges(5,5) * t612 + Ifges(5,6) * t611 + Ifges(5,3) * t625;
t520 = mrSges(6,2) * t541 + mrSges(7,2) * t532 - mrSges(6,3) * t535 - mrSges(7,3) * t534 - qJ(6) * t529 - t664 * t547 + t548 * t669 + t659 * t588 + t663 * t606 + t660 * t624;
t519 = -mrSges(6,1) * t541 - mrSges(7,1) * t534 + mrSges(7,2) * t531 + mrSges(6,3) * t536 - pkin(5) * t529 - t547 * t668 + t664 * t548 + t659 * t589 + t662 * t606 + t658 * t624;
t505 = mrSges(5,2) * t562 - mrSges(5,3) * t542 + Ifges(5,1) * t587 + Ifges(5,4) * t586 + Ifges(5,5) * t610 - pkin(9) * t521 - t633 * t519 + t666 * t520 + t611 * t580 - t625 * t581;
t504 = -mrSges(5,1) * t562 + mrSges(5,3) * t543 + Ifges(5,4) * t587 + Ifges(5,2) * t586 + Ifges(5,6) * t610 - pkin(4) * t644 + pkin(9) * t646 + t666 * t519 + t633 * t520 - t612 * t580 + t625 * t582;
t503 = (mrSges(7,2) * qJ(6) + t662) * t547 + (mrSges(7,2) * pkin(5) - t663) * t548 - t603 * t656 + (qJ(6) * t565 - t658) * t588 + (pkin(5) * t565 + t660) * t589 - qJ(6) * t653 + t667 * t606 + Ifges(4,6) * qJDD(3) - pkin(5) * t645 + t611 * t582 - t612 * t581 + Ifges(4,4) * t617 + Ifges(4,2) * t618 + qJD(3) * t605 - Ifges(5,3) * t610 - Ifges(5,6) * t586 - Ifges(5,5) * t587 + mrSges(4,3) * t570 - mrSges(4,1) * t578 - mrSges(5,1) * t542 + mrSges(5,2) * t543 + mrSges(7,1) * t532 - mrSges(6,1) * t535 + mrSges(6,2) * t536 - mrSges(7,3) * t531 - pkin(4) * t521 - pkin(3) * t515;
t496 = mrSges(4,2) * t578 - mrSges(4,3) * t569 + Ifges(4,1) * t617 + Ifges(4,4) * t618 + Ifges(4,5) * qJDD(3) - pkin(8) * t515 - qJD(3) * t604 - t504 * t634 + t505 * t637 + t603 * t655;
t495 = Ifges(3,6) * qJDD(1) + t641 * Ifges(3,5) - mrSges(3,1) * t630 + mrSges(3,3) * t591 - Ifges(4,5) * t617 - Ifges(4,6) * t618 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t569 + mrSges(4,2) * t570 - t634 * t505 - t637 * t504 - pkin(3) * t642 - pkin(8) * t647 - pkin(2) * t509 + (-t604 * t635 + t605 * t638) * qJD(1);
t494 = mrSges(3,2) * t630 - mrSges(3,3) * t590 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t641 - pkin(7) * t509 + t496 * t638 - t503 * t635;
t493 = -mrSges(2,2) * g(3) - mrSges(2,3) * t621 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t641 - qJ(2) * t502 + t494 * t632 - t495 * t631;
t492 = mrSges(2,1) * g(3) + mrSges(2,3) * t622 + t641 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t652 + qJ(2) * t649 + t631 * t494 + t632 * t495;
t1 = [-m(1) * g(1) + t650; -m(1) * g(2) + t661; (-m(1) - m(2)) * g(3) + t652; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t661 - t636 * t492 + t639 * t493; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t650 + t639 * t492 + t636 * t493; pkin(1) * t502 + mrSges(2,1) * t621 - mrSges(2,2) * t622 + t638 * t503 + pkin(2) * t643 + pkin(7) * t648 + t635 * t496 + mrSges(3,1) * t590 - mrSges(3,2) * t591 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
