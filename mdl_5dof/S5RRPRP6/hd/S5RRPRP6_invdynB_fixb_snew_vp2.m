% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:57:05
% EndTime: 2019-12-31 19:57:11
% DurationCPUTime: 3.94s
% Computational Cost: add. (35721->291), mult. (80792->357), div. (0->0), fcn. (53186->8), ass. (0->113)
t603 = -2 * qJD(3);
t602 = Ifges(5,1) + Ifges(6,1);
t596 = Ifges(5,4) + Ifges(6,4);
t595 = Ifges(5,5) + Ifges(6,5);
t601 = Ifges(5,2) + Ifges(6,2);
t600 = Ifges(5,6) + Ifges(6,6);
t599 = Ifges(5,3) + Ifges(6,3);
t564 = sin(qJ(2));
t567 = cos(qJ(2));
t584 = qJD(1) * qJD(2);
t552 = qJDD(1) * t564 + t567 * t584;
t565 = sin(qJ(1));
t568 = cos(qJ(1));
t558 = -g(1) * t568 - g(2) * t565;
t570 = qJD(1) ^ 2;
t547 = -pkin(1) * t570 + qJDD(1) * pkin(6) + t558;
t593 = t547 * t564;
t598 = pkin(2) * t570;
t508 = qJDD(2) * pkin(2) - qJ(3) * t552 - t593 + (qJ(3) * t584 + t564 * t598 - g(3)) * t567;
t534 = -g(3) * t564 + t567 * t547;
t553 = qJDD(1) * t567 - t564 * t584;
t587 = qJD(1) * t564;
t554 = qJD(2) * pkin(2) - qJ(3) * t587;
t560 = t567 ^ 2;
t509 = qJ(3) * t553 - qJD(2) * t554 - t560 * t598 + t534;
t561 = sin(pkin(8));
t562 = cos(pkin(8));
t542 = (t561 * t567 + t562 * t564) * qJD(1);
t488 = t508 * t562 - t561 * t509 + t542 * t603;
t541 = (t561 * t564 - t562 * t567) * qJD(1);
t597 = -mrSges(5,2) - mrSges(6,2);
t489 = t561 * t508 + t562 * t509 + t541 * t603;
t523 = mrSges(4,1) * t541 + mrSges(4,2) * t542;
t528 = -t552 * t561 + t553 * t562;
t536 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t542;
t524 = pkin(3) * t541 - pkin(7) * t542;
t569 = qJD(2) ^ 2;
t484 = -pkin(3) * t569 + qJDD(2) * pkin(7) - t524 * t541 + t489;
t557 = g(1) * t565 - t568 * g(2);
t574 = -qJDD(1) * pkin(1) - t557;
t513 = -pkin(2) * t553 + qJDD(3) + t554 * t587 + (-qJ(3) * t560 - pkin(6)) * t570 + t574;
t529 = t552 * t562 + t553 * t561;
t487 = (qJD(2) * t541 - t529) * pkin(7) + (qJD(2) * t542 - t528) * pkin(3) + t513;
t563 = sin(qJ(4));
t566 = cos(qJ(4));
t480 = -t484 * t563 + t566 * t487;
t531 = qJD(2) * t566 - t542 * t563;
t503 = qJD(4) * t531 + qJDD(2) * t563 + t529 * t566;
t532 = qJD(2) * t563 + t542 * t566;
t510 = -mrSges(6,1) * t531 + mrSges(6,2) * t532;
t511 = -mrSges(5,1) * t531 + mrSges(5,2) * t532;
t540 = qJD(4) + t541;
t515 = -mrSges(5,2) * t540 + mrSges(5,3) * t531;
t527 = qJDD(4) - t528;
t476 = -0.2e1 * qJD(5) * t532 + (t531 * t540 - t503) * qJ(5) + (t531 * t532 + t527) * pkin(4) + t480;
t514 = -mrSges(6,2) * t540 + mrSges(6,3) * t531;
t583 = m(6) * t476 + t527 * mrSges(6,1) + t540 * t514;
t468 = m(5) * t480 + mrSges(5,1) * t527 + t515 * t540 + (-t510 - t511) * t532 + (-mrSges(5,3) - mrSges(6,3)) * t503 + t583;
t481 = t566 * t484 + t563 * t487;
t502 = -qJD(4) * t532 + qJDD(2) * t566 - t529 * t563;
t516 = pkin(4) * t540 - qJ(5) * t532;
t530 = t531 ^ 2;
t478 = -pkin(4) * t530 + qJ(5) * t502 + 0.2e1 * qJD(5) * t531 - t516 * t540 + t481;
t582 = m(6) * t478 + t502 * mrSges(6,3) + t531 * t510;
t517 = mrSges(6,1) * t540 - mrSges(6,3) * t532;
t588 = -mrSges(5,1) * t540 + mrSges(5,3) * t532 - t517;
t471 = m(5) * t481 + mrSges(5,3) * t502 + t511 * t531 + t597 * t527 + t588 * t540 + t582;
t578 = -t563 * t468 + t566 * t471;
t464 = m(4) * t489 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t528 - qJD(2) * t536 - t523 * t541 + t578;
t535 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t541;
t483 = -qJDD(2) * pkin(3) - pkin(7) * t569 + t542 * t524 - t488;
t479 = -pkin(4) * t502 - qJ(5) * t530 + t516 * t532 + qJDD(5) + t483;
t576 = m(6) * t479 - t502 * mrSges(6,1) - t531 * t514;
t572 = -m(5) * t483 + t502 * mrSges(5,1) + t597 * t503 + t531 * t515 + t588 * t532 - t576;
t473 = m(4) * t488 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t529 + qJD(2) * t535 - t523 * t542 + t572;
t458 = t561 * t464 + t562 * t473;
t533 = -g(3) * t567 - t593;
t551 = (-mrSges(3,1) * t567 + mrSges(3,2) * t564) * qJD(1);
t586 = qJD(1) * t567;
t556 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t586;
t456 = m(3) * t533 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t552 + qJD(2) * t556 - t551 * t587 + t458;
t555 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t587;
t579 = t562 * t464 - t473 * t561;
t457 = m(3) * t534 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t553 - qJD(2) * t555 + t551 * t586 + t579;
t580 = -t456 * t564 + t567 * t457;
t450 = m(2) * t558 - mrSges(2,1) * t570 - qJDD(1) * mrSges(2,2) + t580;
t546 = -pkin(6) * t570 + t574;
t466 = t566 * t468 + t563 * t471;
t573 = m(4) * t513 - t528 * mrSges(4,1) + t529 * mrSges(4,2) + t541 * t535 + t536 * t542 + t466;
t571 = -m(3) * t546 + t553 * mrSges(3,1) - mrSges(3,2) * t552 - t555 * t587 + t556 * t586 - t573;
t461 = m(2) * t557 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t570 + t571;
t592 = t565 * t450 + t568 * t461;
t451 = t567 * t456 + t564 * t457;
t591 = t600 * t531 + t595 * t532 + t599 * t540;
t590 = -t601 * t531 - t596 * t532 - t600 * t540;
t589 = t596 * t531 + t602 * t532 + t595 * t540;
t581 = t568 * t450 - t461 * t565;
t545 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t564 + Ifges(3,4) * t567) * qJD(1);
t544 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t564 + Ifges(3,2) * t567) * qJD(1);
t543 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t564 + Ifges(3,6) * t567) * qJD(1);
t521 = Ifges(4,1) * t542 - Ifges(4,4) * t541 + Ifges(4,5) * qJD(2);
t520 = Ifges(4,4) * t542 - Ifges(4,2) * t541 + Ifges(4,6) * qJD(2);
t519 = Ifges(4,5) * t542 - Ifges(4,6) * t541 + Ifges(4,3) * qJD(2);
t474 = -mrSges(6,3) * t503 - t510 * t532 + t583;
t465 = mrSges(5,2) * t483 + mrSges(6,2) * t479 - mrSges(5,3) * t480 - mrSges(6,3) * t476 - qJ(5) * t474 + t596 * t502 + t602 * t503 + t595 * t527 + t591 * t531 + t590 * t540;
t459 = -mrSges(5,1) * t483 + mrSges(5,3) * t481 - mrSges(6,1) * t479 + mrSges(6,3) * t478 - pkin(4) * t576 + qJ(5) * t582 + (-qJ(5) * t517 + t589) * t540 + (-pkin(4) * t517 - t591) * t532 + (-mrSges(6,2) * qJ(5) + t600) * t527 + (-mrSges(6,2) * pkin(4) + t596) * t503 + t601 * t502;
t452 = -mrSges(4,1) * t513 - mrSges(5,1) * t480 - mrSges(6,1) * t476 + mrSges(5,2) * t481 + mrSges(6,2) * t478 + mrSges(4,3) * t489 + Ifges(4,4) * t529 + Ifges(4,2) * t528 + Ifges(4,6) * qJDD(2) - pkin(3) * t466 - pkin(4) * t474 + qJD(2) * t521 - t542 * t519 + t590 * t532 + t589 * t531 - t599 * t527 - t595 * t503 - t600 * t502;
t447 = mrSges(4,2) * t513 - mrSges(4,3) * t488 + Ifges(4,1) * t529 + Ifges(4,4) * t528 + Ifges(4,5) * qJDD(2) - pkin(7) * t466 - qJD(2) * t520 - t459 * t563 + t465 * t566 - t519 * t541;
t446 = mrSges(3,2) * t546 - mrSges(3,3) * t533 + Ifges(3,1) * t552 + Ifges(3,4) * t553 + Ifges(3,5) * qJDD(2) - qJ(3) * t458 - qJD(2) * t544 + t447 * t562 - t452 * t561 + t543 * t586;
t445 = -mrSges(3,1) * t546 + mrSges(3,3) * t534 + Ifges(3,4) * t552 + Ifges(3,2) * t553 + Ifges(3,6) * qJDD(2) - pkin(2) * t573 + qJ(3) * t579 + qJD(2) * t545 + t561 * t447 + t562 * t452 - t543 * t587;
t444 = -pkin(1) * t451 + mrSges(2,3) * t558 - pkin(2) * t458 - Ifges(3,5) * t552 - Ifges(3,6) * t553 - mrSges(3,1) * t533 + mrSges(3,2) * t534 - t563 * t465 - t566 * t459 - pkin(3) * t572 - pkin(7) * t578 - Ifges(4,5) * t529 - Ifges(4,6) * t528 - mrSges(4,1) * t488 + mrSges(4,2) * t489 + mrSges(2,1) * g(3) - t541 * t521 - t542 * t520 + t570 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t544 * t564 + t545 * t567) * qJD(1);
t443 = -mrSges(2,2) * g(3) - mrSges(2,3) * t557 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t570 - pkin(6) * t451 - t445 * t564 + t446 * t567;
t1 = [-m(1) * g(1) + t581; -m(1) * g(2) + t592; (-m(1) - m(2)) * g(3) + t451; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t592 + t568 * t443 - t565 * t444; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t581 + t565 * t443 + t568 * t444; -mrSges(1,1) * g(2) + mrSges(2,1) * t557 + mrSges(1,2) * g(1) - mrSges(2,2) * t558 + Ifges(2,3) * qJDD(1) + pkin(1) * t571 + pkin(6) * t580 + t567 * t445 + t564 * t446;];
tauB = t1;
