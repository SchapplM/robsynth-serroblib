% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:20:18
% EndTime: 2019-05-05 19:20:25
% DurationCPUTime: 6.81s
% Computational Cost: add. (101807->338), mult. (223193->416), div. (0->0), fcn. (150282->10), ass. (0->130)
t655 = sin(qJ(1));
t659 = cos(qJ(1));
t637 = t655 * g(1) - t659 * g(2);
t661 = qJD(1) ^ 2;
t668 = -t661 * qJ(2) + qJDD(2) - t637;
t687 = -pkin(1) - pkin(7);
t614 = qJDD(1) * t687 + t668;
t654 = sin(qJ(3));
t658 = cos(qJ(3));
t603 = t654 * g(3) + t658 * t614;
t679 = qJD(1) * qJD(3);
t677 = t654 * t679;
t633 = qJDD(1) * t658 - t677;
t581 = (-t633 - t677) * qJ(4) + (-t654 * t658 * t661 + qJDD(3)) * pkin(3) + t603;
t604 = -g(3) * t658 + t654 * t614;
t632 = -qJDD(1) * t654 - t658 * t679;
t681 = qJD(1) * t658;
t635 = qJD(3) * pkin(3) - qJ(4) * t681;
t647 = t654 ^ 2;
t582 = -pkin(3) * t647 * t661 + qJ(4) * t632 - qJD(3) * t635 + t604;
t650 = sin(pkin(10));
t651 = cos(pkin(10));
t682 = qJD(1) * t654;
t622 = -t650 * t682 + t651 * t681;
t562 = -0.2e1 * qJD(4) * t622 + t581 * t651 - t650 * t582;
t638 = -t659 * g(1) - t655 * g(2);
t669 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t638;
t686 = mrSges(2,1) - mrSges(3,2);
t685 = Ifges(2,5) - Ifges(3,4);
t684 = -Ifges(2,6) + Ifges(3,5);
t621 = -t650 * t681 - t651 * t682;
t563 = 0.2e1 * qJD(4) * t621 + t650 * t581 + t651 * t582;
t596 = -mrSges(5,1) * t621 + mrSges(5,2) * t622;
t601 = t632 * t651 - t650 * t633;
t613 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t622;
t597 = -pkin(4) * t621 - pkin(8) * t622;
t660 = qJD(3) ^ 2;
t555 = -pkin(4) * t660 + qJDD(3) * pkin(8) + t597 * t621 + t563;
t584 = -t632 * pkin(3) + qJDD(4) + t635 * t681 + (-qJ(4) * t647 + t687) * t661 + t669;
t602 = t632 * t650 + t633 * t651;
t561 = (-qJD(3) * t621 - t602) * pkin(8) + (qJD(3) * t622 - t601) * pkin(4) + t584;
t653 = sin(qJ(5));
t657 = cos(qJ(5));
t550 = -t653 * t555 + t657 * t561;
t606 = qJD(3) * t657 - t622 * t653;
t577 = qJD(5) * t606 + qJDD(3) * t653 + t602 * t657;
t600 = qJDD(5) - t601;
t607 = qJD(3) * t653 + t622 * t657;
t619 = qJD(5) - t621;
t548 = (t606 * t619 - t577) * pkin(9) + (t606 * t607 + t600) * pkin(5) + t550;
t551 = t657 * t555 + t653 * t561;
t576 = -qJD(5) * t607 + qJDD(3) * t657 - t602 * t653;
t591 = pkin(5) * t619 - pkin(9) * t607;
t605 = t606 ^ 2;
t549 = -pkin(5) * t605 + pkin(9) * t576 - t591 * t619 + t551;
t652 = sin(qJ(6));
t656 = cos(qJ(6));
t546 = t548 * t656 - t549 * t652;
t585 = t606 * t656 - t607 * t652;
t558 = qJD(6) * t585 + t576 * t652 + t577 * t656;
t586 = t606 * t652 + t607 * t656;
t568 = -mrSges(7,1) * t585 + mrSges(7,2) * t586;
t616 = qJD(6) + t619;
t569 = -mrSges(7,2) * t616 + mrSges(7,3) * t585;
t598 = qJDD(6) + t600;
t544 = m(7) * t546 + mrSges(7,1) * t598 - mrSges(7,3) * t558 - t568 * t586 + t569 * t616;
t547 = t548 * t652 + t549 * t656;
t557 = -qJD(6) * t586 + t576 * t656 - t577 * t652;
t570 = mrSges(7,1) * t616 - mrSges(7,3) * t586;
t545 = m(7) * t547 - mrSges(7,2) * t598 + mrSges(7,3) * t557 + t568 * t585 - t570 * t616;
t536 = t656 * t544 + t652 * t545;
t587 = -mrSges(6,1) * t606 + mrSges(6,2) * t607;
t589 = -mrSges(6,2) * t619 + mrSges(6,3) * t606;
t534 = m(6) * t550 + mrSges(6,1) * t600 - mrSges(6,3) * t577 - t587 * t607 + t589 * t619 + t536;
t590 = mrSges(6,1) * t619 - mrSges(6,3) * t607;
t672 = -t544 * t652 + t656 * t545;
t535 = m(6) * t551 - mrSges(6,2) * t600 + mrSges(6,3) * t576 + t587 * t606 - t590 * t619 + t672;
t673 = -t534 * t653 + t657 * t535;
t529 = m(5) * t563 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t601 - qJD(3) * t613 + t596 * t621 + t673;
t612 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t621;
t554 = -qJDD(3) * pkin(4) - pkin(8) * t660 + t622 * t597 - t562;
t552 = -pkin(5) * t576 - pkin(9) * t605 + t591 * t607 + t554;
t666 = m(7) * t552 - t557 * mrSges(7,1) + mrSges(7,2) * t558 - t585 * t569 + t570 * t586;
t663 = -m(6) * t554 + t576 * mrSges(6,1) - mrSges(6,2) * t577 + t606 * t589 - t590 * t607 - t666;
t540 = m(5) * t562 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t602 + qJD(3) * t612 - t596 * t622 + t663;
t521 = t650 * t529 + t651 * t540;
t631 = (mrSges(4,1) * t654 + mrSges(4,2) * t658) * qJD(1);
t634 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t682;
t519 = m(4) * t603 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t633 + qJD(3) * t634 - t631 * t681 + t521;
t636 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t681;
t674 = t651 * t529 - t540 * t650;
t520 = m(4) * t604 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t632 - qJD(3) * t636 - t631 * t682 + t674;
t516 = t658 * t519 + t654 * t520;
t620 = -qJDD(1) * pkin(1) + t668;
t667 = -m(3) * t620 + t661 * mrSges(3,3) - t516;
t514 = m(2) * t637 - t661 * mrSges(2,2) + qJDD(1) * t686 + t667;
t615 = t661 * pkin(1) - t669;
t611 = t661 * t687 + t669;
t530 = t657 * t534 + t653 * t535;
t665 = m(5) * t584 - mrSges(5,1) * t601 + t602 * mrSges(5,2) - t612 * t621 + t622 * t613 + t530;
t664 = -m(4) * t611 + mrSges(4,1) * t632 - t633 * mrSges(4,2) - t634 * t682 - t636 * t681 - t665;
t662 = -m(3) * t615 + t661 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t664;
t526 = m(2) * t638 - mrSges(2,1) * t661 - qJDD(1) * mrSges(2,2) + t662;
t683 = t659 * t514 + t655 * t526;
t676 = -t514 * t655 + t659 * t526;
t675 = -t654 * t519 + t658 * t520;
t625 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t658 - Ifges(4,4) * t654) * qJD(1);
t624 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t658 - Ifges(4,2) * t654) * qJD(1);
t623 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t658 - Ifges(4,6) * t654) * qJD(1);
t594 = Ifges(5,1) * t622 + Ifges(5,4) * t621 + Ifges(5,5) * qJD(3);
t593 = Ifges(5,4) * t622 + Ifges(5,2) * t621 + Ifges(5,6) * qJD(3);
t592 = Ifges(5,5) * t622 + Ifges(5,6) * t621 + Ifges(5,3) * qJD(3);
t573 = Ifges(6,1) * t607 + Ifges(6,4) * t606 + Ifges(6,5) * t619;
t572 = Ifges(6,4) * t607 + Ifges(6,2) * t606 + Ifges(6,6) * t619;
t571 = Ifges(6,5) * t607 + Ifges(6,6) * t606 + Ifges(6,3) * t619;
t566 = Ifges(7,1) * t586 + Ifges(7,4) * t585 + Ifges(7,5) * t616;
t565 = Ifges(7,4) * t586 + Ifges(7,2) * t585 + Ifges(7,6) * t616;
t564 = Ifges(7,5) * t586 + Ifges(7,6) * t585 + Ifges(7,3) * t616;
t538 = mrSges(7,2) * t552 - mrSges(7,3) * t546 + Ifges(7,1) * t558 + Ifges(7,4) * t557 + Ifges(7,5) * t598 + t564 * t585 - t565 * t616;
t537 = -mrSges(7,1) * t552 + mrSges(7,3) * t547 + Ifges(7,4) * t558 + Ifges(7,2) * t557 + Ifges(7,6) * t598 - t564 * t586 + t566 * t616;
t523 = mrSges(6,2) * t554 - mrSges(6,3) * t550 + Ifges(6,1) * t577 + Ifges(6,4) * t576 + Ifges(6,5) * t600 - pkin(9) * t536 - t537 * t652 + t538 * t656 + t571 * t606 - t572 * t619;
t522 = -mrSges(6,1) * t554 + mrSges(6,3) * t551 + Ifges(6,4) * t577 + Ifges(6,2) * t576 + Ifges(6,6) * t600 - pkin(5) * t666 + pkin(9) * t672 + t656 * t537 + t652 * t538 - t607 * t571 + t619 * t573;
t517 = Ifges(5,4) * t602 + Ifges(5,2) * t601 + Ifges(5,6) * qJDD(3) - t622 * t592 + qJD(3) * t594 - mrSges(5,1) * t584 + mrSges(5,3) * t563 - Ifges(6,5) * t577 - Ifges(6,6) * t576 - Ifges(6,3) * t600 - t607 * t572 + t606 * t573 - mrSges(6,1) * t550 + mrSges(6,2) * t551 - Ifges(7,5) * t558 - Ifges(7,6) * t557 - Ifges(7,3) * t598 - t586 * t565 + t585 * t566 - mrSges(7,1) * t546 + mrSges(7,2) * t547 - pkin(5) * t536 - pkin(4) * t530;
t515 = -m(3) * g(3) + t675;
t512 = mrSges(5,2) * t584 - mrSges(5,3) * t562 + Ifges(5,1) * t602 + Ifges(5,4) * t601 + Ifges(5,5) * qJDD(3) - pkin(8) * t530 - qJD(3) * t593 - t522 * t653 + t523 * t657 + t592 * t621;
t511 = mrSges(4,2) * t611 - mrSges(4,3) * t603 + Ifges(4,1) * t633 + Ifges(4,4) * t632 + Ifges(4,5) * qJDD(3) - qJ(4) * t521 - qJD(3) * t624 + t512 * t651 - t517 * t650 - t623 * t682;
t510 = -mrSges(4,1) * t611 + mrSges(4,3) * t604 + Ifges(4,4) * t633 + Ifges(4,2) * t632 + Ifges(4,6) * qJDD(3) - pkin(3) * t665 + qJ(4) * t674 + qJD(3) * t625 + t650 * t512 + t651 * t517 - t623 * t681;
t509 = pkin(3) * t521 + pkin(2) * t516 - qJ(2) * t515 + t657 * t522 + t653 * t523 + pkin(8) * t673 + Ifges(4,6) * t632 + Ifges(4,5) * t633 - mrSges(2,3) * t637 + mrSges(3,1) * t620 - t621 * t594 + t622 * t593 + pkin(4) * t663 + Ifges(5,6) * t601 + Ifges(5,5) * t602 + mrSges(4,1) * t603 - mrSges(4,2) * t604 + mrSges(5,1) * t562 - mrSges(5,2) * t563 + t684 * t661 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t685 * qJDD(1) + (t624 * t658 + t625 * t654) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t508 = -mrSges(3,1) * t615 + mrSges(2,3) * t638 - pkin(1) * t515 - pkin(2) * t664 - pkin(7) * t675 + g(3) * t686 - qJDD(1) * t684 - t658 * t510 - t654 * t511 + t661 * t685;
t1 = [-m(1) * g(1) + t676; -m(1) * g(2) + t683; (-m(1) - m(2) - m(3)) * g(3) + t675; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t683 - t655 * t508 + t659 * t509; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t676 + t659 * t508 + t655 * t509; pkin(1) * t667 + qJ(2) * t662 - t654 * t510 - pkin(7) * t516 + mrSges(2,1) * t637 - mrSges(2,2) * t638 + t658 * t511 + mrSges(3,2) * t620 - mrSges(3,3) * t615 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
