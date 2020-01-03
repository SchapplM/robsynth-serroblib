% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:10
% EndTime: 2019-12-31 18:04:13
% DurationCPUTime: 2.87s
% Computational Cost: add. (23577->252), mult. (56733->317), div. (0->0), fcn. (37719->8), ass. (0->104)
t570 = cos(pkin(8));
t611 = (Ifges(3,6) - Ifges(4,6)) * t570;
t573 = sin(qJ(1));
t576 = cos(qJ(1));
t545 = -t576 * g(1) - t573 * g(2);
t577 = qJD(1) ^ 2;
t610 = -t577 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t545;
t609 = Ifges(4,4) + Ifges(3,5);
t569 = sin(pkin(8));
t564 = t569 ^ 2;
t565 = t570 ^ 2;
t601 = t565 * t577;
t607 = t564 * t577 + t601;
t544 = t573 * g(1) - t576 * g(2);
t534 = -qJDD(1) * pkin(1) - t577 * qJ(2) + qJDD(2) - t544;
t592 = qJDD(1) * t570;
t593 = qJDD(1) * t569;
t596 = t569 * qJD(1);
t521 = -pkin(2) * t592 - qJ(3) * t593 - 0.2e1 * qJD(3) * t596 + t534;
t524 = -t570 * g(3) - t610 * t569;
t606 = Ifges(3,1) + Ifges(4,1);
t605 = Ifges(3,4) - Ifges(4,5);
t604 = Ifges(3,2) + Ifges(4,3);
t603 = mrSges(3,2) * t569;
t540 = (-mrSges(4,1) * t570 - mrSges(4,3) * t569) * qJD(1);
t541 = (-mrSges(3,1) * t570 + t603) * qJD(1);
t539 = (-pkin(2) * t570 - qJ(3) * t569) * qJD(1);
t506 = t539 * t596 + qJDD(3) - t524;
t501 = (-pkin(3) * t570 * t577 - pkin(6) * qJDD(1)) * t569 + t506;
t525 = -t569 * g(3) + t610 * t570;
t595 = t570 * qJD(1);
t508 = t539 * t595 + t525;
t503 = -pkin(3) * t601 - pkin(6) * t592 + t508;
t572 = sin(qJ(4));
t575 = cos(qJ(4));
t487 = t575 * t501 - t572 * t503;
t584 = t569 * t575 - t570 * t572;
t583 = -t569 * t572 - t570 * t575;
t536 = t583 * qJD(1);
t597 = t536 * qJD(4);
t523 = t584 * qJDD(1) + t597;
t537 = t584 * qJD(1);
t483 = (-t523 + t597) * pkin(7) + (t536 * t537 + qJDD(4)) * pkin(4) + t487;
t488 = t572 * t501 + t575 * t503;
t522 = -t537 * qJD(4) + t583 * qJDD(1);
t528 = qJD(4) * pkin(4) - t537 * pkin(7);
t535 = t536 ^ 2;
t484 = -t535 * pkin(4) + t522 * pkin(7) - qJD(4) * t528 + t488;
t571 = sin(qJ(5));
t574 = cos(qJ(5));
t481 = t574 * t483 - t571 * t484;
t514 = t574 * t536 - t571 * t537;
t492 = t514 * qJD(5) + t571 * t522 + t574 * t523;
t515 = t571 * t536 + t574 * t537;
t498 = -t514 * mrSges(6,1) + t515 * mrSges(6,2);
t566 = qJD(4) + qJD(5);
t509 = -t566 * mrSges(6,2) + t514 * mrSges(6,3);
t563 = qJDD(4) + qJDD(5);
t479 = m(6) * t481 + t563 * mrSges(6,1) - t492 * mrSges(6,3) - t515 * t498 + t566 * t509;
t482 = t571 * t483 + t574 * t484;
t491 = -t515 * qJD(5) + t574 * t522 - t571 * t523;
t510 = t566 * mrSges(6,1) - t515 * mrSges(6,3);
t480 = m(6) * t482 - t563 * mrSges(6,2) + t491 * mrSges(6,3) + t514 * t498 - t566 * t510;
t470 = t574 * t479 + t571 * t480;
t518 = -t536 * mrSges(5,1) + t537 * mrSges(5,2);
t526 = -qJD(4) * mrSges(5,2) + t536 * mrSges(5,3);
t468 = m(5) * t487 + qJDD(4) * mrSges(5,1) - t523 * mrSges(5,3) + qJD(4) * t526 - t537 * t518 + t470;
t527 = qJD(4) * mrSges(5,1) - t537 * mrSges(5,3);
t586 = -t571 * t479 + t574 * t480;
t469 = m(5) * t488 - qJDD(4) * mrSges(5,2) + t522 * mrSges(5,3) - qJD(4) * t527 + t536 * t518 + t586;
t466 = t575 * t468 + t572 * t469;
t580 = -m(4) * t506 - t466;
t464 = m(3) * t524 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t540 - t541) * qJD(1)) * t569 + t580;
t587 = -t572 * t468 + t575 * t469;
t582 = m(4) * t508 + mrSges(4,2) * t592 + t540 * t595 + t587;
t465 = m(3) * t525 + (qJDD(1) * mrSges(3,3) + qJD(1) * t541) * t570 + t582;
t588 = -t569 * t464 + t570 * t465;
t457 = m(2) * t545 - t577 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t588;
t505 = pkin(3) * t592 + (-t564 - t565) * t577 * pkin(6) - t521;
t486 = -t522 * pkin(4) - t535 * pkin(7) + t537 * t528 + t505;
t585 = m(6) * t486 - t491 * mrSges(6,1) + t492 * mrSges(6,2) - t514 * t509 + t515 * t510;
t579 = -m(5) * t505 + t522 * mrSges(5,1) - t523 * mrSges(5,2) + t536 * t526 - t537 * t527 - t585;
t475 = m(4) * t521 - mrSges(4,1) * t592 - t607 * mrSges(4,2) - mrSges(4,3) * t593 + t579;
t578 = -m(3) * t534 + mrSges(3,1) * t592 + t607 * mrSges(3,3) - t475;
t474 = t578 - t577 * mrSges(2,2) + m(2) * t544 + (mrSges(2,1) - t603) * qJDD(1);
t600 = t573 * t457 + t576 * t474;
t458 = t570 * t464 + t569 * t465;
t598 = (t609 * t569 + t611) * qJD(1);
t589 = t576 * t457 - t573 * t474;
t513 = Ifges(5,1) * t537 + Ifges(5,4) * t536 + Ifges(5,5) * qJD(4);
t512 = Ifges(5,4) * t537 + Ifges(5,2) * t536 + Ifges(5,6) * qJD(4);
t511 = Ifges(5,5) * t537 + Ifges(5,6) * t536 + Ifges(5,3) * qJD(4);
t495 = Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t566;
t494 = Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t566;
t493 = Ifges(6,5) * t515 + Ifges(6,6) * t514 + Ifges(6,3) * t566;
t472 = mrSges(6,2) * t486 - mrSges(6,3) * t481 + Ifges(6,1) * t492 + Ifges(6,4) * t491 + Ifges(6,5) * t563 + t514 * t493 - t566 * t494;
t471 = -mrSges(6,1) * t486 + mrSges(6,3) * t482 + Ifges(6,4) * t492 + Ifges(6,2) * t491 + Ifges(6,6) * t563 - t515 * t493 + t566 * t495;
t460 = mrSges(5,2) * t505 - mrSges(5,3) * t487 + Ifges(5,1) * t523 + Ifges(5,4) * t522 + Ifges(5,5) * qJDD(4) - pkin(7) * t470 - qJD(4) * t512 - t571 * t471 + t574 * t472 + t536 * t511;
t459 = -mrSges(5,1) * t505 + mrSges(5,3) * t488 + Ifges(5,4) * t523 + Ifges(5,2) * t522 + Ifges(5,6) * qJDD(4) - pkin(4) * t585 + pkin(7) * t586 + qJD(4) * t513 + t574 * t471 + t571 * t472 - t537 * t511;
t454 = mrSges(3,2) * t534 + mrSges(4,2) * t506 - mrSges(3,3) * t524 - mrSges(4,3) * t521 - pkin(6) * t466 - qJ(3) * t475 - t572 * t459 + t575 * t460 + t598 * t595 + (t606 * t569 + t605 * t570) * qJDD(1);
t453 = -mrSges(3,1) * t534 + mrSges(3,3) * t525 - mrSges(4,1) * t521 + mrSges(4,2) * t508 - t572 * t460 - t575 * t459 - pkin(3) * t579 - pkin(6) * t587 - pkin(2) * t475 - t598 * t596 + (t605 * t569 + t604 * t570) * qJDD(1);
t452 = -t514 * t495 + t515 * t494 + Ifges(6,6) * t491 + Ifges(6,5) * t492 + Ifges(5,6) * t522 + Ifges(5,5) * t523 - mrSges(3,1) * t524 + mrSges(3,2) * t525 + pkin(3) * t466 - mrSges(6,2) * t482 - pkin(2) * t580 + Ifges(6,3) * t563 + (Ifges(2,6) - t611 + (mrSges(4,2) * pkin(2) - t609) * t569) * qJDD(1) - pkin(1) * t458 + Ifges(5,3) * qJDD(4) + t577 * Ifges(2,5) + mrSges(2,1) * g(3) - t536 * t513 + t537 * t512 + mrSges(2,3) * t545 + mrSges(5,1) * t487 - mrSges(5,2) * t488 + mrSges(4,1) * t506 - mrSges(4,3) * t508 + pkin(4) * t470 - qJ(3) * t582 + (t605 * t565 * qJD(1) + (pkin(2) * t540 - t605 * t596 + (-t604 + t606) * t595) * t569) * qJD(1) + mrSges(6,1) * t481;
t451 = -mrSges(2,2) * g(3) - mrSges(2,3) * t544 + Ifges(2,5) * qJDD(1) - t577 * Ifges(2,6) - qJ(2) * t458 - t569 * t453 + t570 * t454;
t1 = [-m(1) * g(1) + t589; -m(1) * g(2) + t600; (-m(1) - m(2)) * g(3) + t458; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t600 + t576 * t451 - t573 * t452; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t589 + t573 * t451 + t576 * t452; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t544 - mrSges(2,2) * t545 + t569 * t454 + t570 * t453 + pkin(1) * (-mrSges(3,2) * t593 + t578) + qJ(2) * t588;];
tauB = t1;
