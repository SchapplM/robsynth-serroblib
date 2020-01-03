% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:28:26
% EndTime: 2019-12-31 22:28:37
% DurationCPUTime: 7.49s
% Computational Cost: add. (114983->312), mult. (233456->393), div. (0->0), fcn. (163843->10), ass. (0->122)
t571 = sin(qJ(1));
t576 = cos(qJ(1));
t561 = -t576 * g(1) - t571 * g(2);
t578 = qJD(1) ^ 2;
t545 = -t578 * pkin(1) + qJDD(1) * pkin(6) + t561;
t570 = sin(qJ(2));
t575 = cos(qJ(2));
t535 = -t570 * g(3) + t575 * t545;
t553 = (-mrSges(3,1) * t575 + mrSges(3,2) * t570) * qJD(1);
t589 = qJD(1) * qJD(2);
t564 = t570 * t589;
t556 = t575 * qJDD(1) - t564;
t591 = qJD(1) * t570;
t558 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t591;
t560 = t571 * g(1) - t576 * g(2);
t544 = -qJDD(1) * pkin(1) - t578 * pkin(6) - t560;
t588 = t575 * t589;
t555 = t570 * qJDD(1) + t588;
t513 = (-t555 - t588) * pkin(7) + (-t556 + t564) * pkin(2) + t544;
t554 = (-pkin(2) * t575 - pkin(7) * t570) * qJD(1);
t577 = qJD(2) ^ 2;
t590 = t575 * qJD(1);
t516 = -t577 * pkin(2) + qJDD(2) * pkin(7) + t554 * t590 + t535;
t569 = sin(qJ(3));
t574 = cos(qJ(3));
t500 = t574 * t513 - t569 * t516;
t551 = t574 * qJD(2) - t569 * t591;
t527 = t551 * qJD(3) + t569 * qJDD(2) + t574 * t555;
t550 = qJDD(3) - t556;
t552 = t569 * qJD(2) + t574 * t591;
t563 = qJD(3) - t590;
t487 = (t551 * t563 - t527) * pkin(8) + (t551 * t552 + t550) * pkin(3) + t500;
t501 = t569 * t513 + t574 * t516;
t526 = -t552 * qJD(3) + t574 * qJDD(2) - t569 * t555;
t536 = t563 * pkin(3) - t552 * pkin(8);
t549 = t551 ^ 2;
t489 = -t549 * pkin(3) + t526 * pkin(8) - t563 * t536 + t501;
t568 = sin(qJ(4));
t573 = cos(qJ(4));
t477 = t573 * t487 - t568 * t489;
t529 = t573 * t551 - t568 * t552;
t499 = t529 * qJD(4) + t568 * t526 + t573 * t527;
t530 = t568 * t551 + t573 * t552;
t546 = qJDD(4) + t550;
t562 = qJD(4) + t563;
t475 = (t529 * t562 - t499) * pkin(9) + (t529 * t530 + t546) * pkin(4) + t477;
t478 = t568 * t487 + t573 * t489;
t498 = -t530 * qJD(4) + t573 * t526 - t568 * t527;
t519 = t562 * pkin(4) - t530 * pkin(9);
t528 = t529 ^ 2;
t476 = -t528 * pkin(4) + t498 * pkin(9) - t562 * t519 + t478;
t567 = sin(qJ(5));
t572 = cos(qJ(5));
t473 = t572 * t475 - t567 * t476;
t508 = t572 * t529 - t567 * t530;
t484 = t508 * qJD(5) + t567 * t498 + t572 * t499;
t509 = t567 * t529 + t572 * t530;
t495 = -t508 * mrSges(6,1) + t509 * mrSges(6,2);
t557 = qJD(5) + t562;
t502 = -t557 * mrSges(6,2) + t508 * mrSges(6,3);
t540 = qJDD(5) + t546;
t471 = m(6) * t473 + t540 * mrSges(6,1) - t484 * mrSges(6,3) - t509 * t495 + t557 * t502;
t474 = t567 * t475 + t572 * t476;
t483 = -t509 * qJD(5) + t572 * t498 - t567 * t499;
t503 = t557 * mrSges(6,1) - t509 * mrSges(6,3);
t472 = m(6) * t474 - t540 * mrSges(6,2) + t483 * mrSges(6,3) + t508 * t495 - t557 * t503;
t463 = t572 * t471 + t567 * t472;
t510 = -t529 * mrSges(5,1) + t530 * mrSges(5,2);
t517 = -t562 * mrSges(5,2) + t529 * mrSges(5,3);
t461 = m(5) * t477 + t546 * mrSges(5,1) - t499 * mrSges(5,3) - t530 * t510 + t562 * t517 + t463;
t518 = t562 * mrSges(5,1) - t530 * mrSges(5,3);
t583 = -t567 * t471 + t572 * t472;
t462 = m(5) * t478 - t546 * mrSges(5,2) + t498 * mrSges(5,3) + t529 * t510 - t562 * t518 + t583;
t457 = t573 * t461 + t568 * t462;
t531 = -t551 * mrSges(4,1) + t552 * mrSges(4,2);
t532 = -t563 * mrSges(4,2) + t551 * mrSges(4,3);
t455 = m(4) * t500 + t550 * mrSges(4,1) - t527 * mrSges(4,3) - t552 * t531 + t563 * t532 + t457;
t533 = t563 * mrSges(4,1) - t552 * mrSges(4,3);
t584 = -t568 * t461 + t573 * t462;
t456 = m(4) * t501 - t550 * mrSges(4,2) + t526 * mrSges(4,3) + t551 * t531 - t563 * t533 + t584;
t585 = -t569 * t455 + t574 * t456;
t450 = m(3) * t535 - qJDD(2) * mrSges(3,2) + t556 * mrSges(3,3) - qJD(2) * t558 + t553 * t590 + t585;
t534 = -t575 * g(3) - t570 * t545;
t559 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t590;
t515 = -qJDD(2) * pkin(2) - t577 * pkin(7) + t554 * t591 - t534;
t496 = -t526 * pkin(3) - t549 * pkin(8) + t552 * t536 + t515;
t480 = -t498 * pkin(4) - t528 * pkin(9) + t530 * t519 + t496;
t582 = m(6) * t480 - t483 * mrSges(6,1) + t484 * mrSges(6,2) - t508 * t502 + t509 * t503;
t581 = m(5) * t496 - t498 * mrSges(5,1) + t499 * mrSges(5,2) - t529 * t517 + t530 * t518 + t582;
t579 = -m(4) * t515 + t526 * mrSges(4,1) - t527 * mrSges(4,2) + t551 * t532 - t552 * t533 - t581;
t467 = m(3) * t534 + qJDD(2) * mrSges(3,1) - t555 * mrSges(3,3) + qJD(2) * t559 - t553 * t591 + t579;
t586 = t575 * t450 - t570 * t467;
t444 = m(2) * t561 - t578 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t586;
t451 = t574 * t455 + t569 * t456;
t580 = -m(3) * t544 + t556 * mrSges(3,1) - t555 * mrSges(3,2) - t558 * t591 + t559 * t590 - t451;
t447 = m(2) * t560 + qJDD(1) * mrSges(2,1) - t578 * mrSges(2,2) + t580;
t592 = t571 * t444 + t576 * t447;
t445 = t570 * t450 + t575 * t467;
t587 = t576 * t444 - t571 * t447;
t543 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t570 + Ifges(3,4) * t575) * qJD(1);
t542 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t570 + Ifges(3,2) * t575) * qJD(1);
t541 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t570 + Ifges(3,6) * t575) * qJD(1);
t522 = Ifges(4,1) * t552 + Ifges(4,4) * t551 + Ifges(4,5) * t563;
t521 = Ifges(4,4) * t552 + Ifges(4,2) * t551 + Ifges(4,6) * t563;
t520 = Ifges(4,5) * t552 + Ifges(4,6) * t551 + Ifges(4,3) * t563;
t506 = Ifges(5,1) * t530 + Ifges(5,4) * t529 + Ifges(5,5) * t562;
t505 = Ifges(5,4) * t530 + Ifges(5,2) * t529 + Ifges(5,6) * t562;
t504 = Ifges(5,5) * t530 + Ifges(5,6) * t529 + Ifges(5,3) * t562;
t492 = Ifges(6,1) * t509 + Ifges(6,4) * t508 + Ifges(6,5) * t557;
t491 = Ifges(6,4) * t509 + Ifges(6,2) * t508 + Ifges(6,6) * t557;
t490 = Ifges(6,5) * t509 + Ifges(6,6) * t508 + Ifges(6,3) * t557;
t465 = mrSges(6,2) * t480 - mrSges(6,3) * t473 + Ifges(6,1) * t484 + Ifges(6,4) * t483 + Ifges(6,5) * t540 + t508 * t490 - t557 * t491;
t464 = -mrSges(6,1) * t480 + mrSges(6,3) * t474 + Ifges(6,4) * t484 + Ifges(6,2) * t483 + Ifges(6,6) * t540 - t509 * t490 + t557 * t492;
t453 = mrSges(5,2) * t496 - mrSges(5,3) * t477 + Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t546 - pkin(9) * t463 - t567 * t464 + t572 * t465 + t529 * t504 - t562 * t505;
t452 = -mrSges(5,1) * t496 + mrSges(5,3) * t478 + Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t546 - pkin(4) * t582 + pkin(9) * t583 + t572 * t464 + t567 * t465 - t530 * t504 + t562 * t506;
t441 = mrSges(4,2) * t515 - mrSges(4,3) * t500 + Ifges(4,1) * t527 + Ifges(4,4) * t526 + Ifges(4,5) * t550 - pkin(8) * t457 - t568 * t452 + t573 * t453 + t551 * t520 - t563 * t521;
t440 = -mrSges(4,1) * t515 + mrSges(4,3) * t501 + Ifges(4,4) * t527 + Ifges(4,2) * t526 + Ifges(4,6) * t550 - pkin(3) * t581 + pkin(8) * t584 + t573 * t452 + t568 * t453 - t552 * t520 + t563 * t522;
t439 = -t541 * t591 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t555 + Ifges(3,2) * t556 - Ifges(5,3) * t546 - Ifges(4,3) * t550 + t551 * t522 - t552 * t521 - t530 * t505 + mrSges(3,3) * t535 - Ifges(6,3) * t540 + qJD(2) * t543 - mrSges(3,1) * t544 - Ifges(4,6) * t526 - Ifges(4,5) * t527 + t529 * t506 + mrSges(4,2) * t501 + t508 * t492 - t509 * t491 - Ifges(5,6) * t498 - Ifges(5,5) * t499 - mrSges(4,1) * t500 - Ifges(6,5) * t484 + mrSges(5,2) * t478 - Ifges(6,6) * t483 - mrSges(5,1) * t477 - mrSges(6,1) * t473 + mrSges(6,2) * t474 - pkin(4) * t463 - pkin(3) * t457 - pkin(2) * t451;
t438 = mrSges(3,2) * t544 - mrSges(3,3) * t534 + Ifges(3,1) * t555 + Ifges(3,4) * t556 + Ifges(3,5) * qJDD(2) - pkin(7) * t451 - qJD(2) * t542 - t569 * t440 + t574 * t441 + t541 * t590;
t437 = Ifges(2,6) * qJDD(1) + t578 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t561 - Ifges(3,5) * t555 - Ifges(3,6) * t556 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t534 + mrSges(3,2) * t535 - t569 * t441 - t574 * t440 - pkin(2) * t579 - pkin(7) * t585 - pkin(1) * t445 + (-t570 * t542 + t575 * t543) * qJD(1);
t436 = -mrSges(2,2) * g(3) - mrSges(2,3) * t560 + Ifges(2,5) * qJDD(1) - t578 * Ifges(2,6) - pkin(6) * t445 + t575 * t438 - t570 * t439;
t1 = [-m(1) * g(1) + t587; -m(1) * g(2) + t592; (-m(1) - m(2)) * g(3) + t445; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t592 + t576 * t436 - t571 * t437; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t587 + t571 * t436 + t576 * t437; -mrSges(1,1) * g(2) + mrSges(2,1) * t560 + mrSges(1,2) * g(1) - mrSges(2,2) * t561 + Ifges(2,3) * qJDD(1) + pkin(1) * t580 + pkin(6) * t586 + t570 * t438 + t575 * t439;];
tauB = t1;
