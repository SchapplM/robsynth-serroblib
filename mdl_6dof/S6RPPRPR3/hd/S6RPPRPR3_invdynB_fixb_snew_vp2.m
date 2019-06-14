% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:09:13
% EndTime: 2019-05-05 14:09:18
% DurationCPUTime: 4.39s
% Computational Cost: add. (51987->296), mult. (105556->363), div. (0->0), fcn. (63707->10), ass. (0->119)
t565 = sin(qJ(1));
t568 = cos(qJ(1));
t545 = t565 * g(1) - t568 * g(2);
t537 = qJDD(1) * pkin(1) + t545;
t546 = -t568 * g(1) - t565 * g(2);
t570 = qJD(1) ^ 2;
t539 = -t570 * pkin(1) + t546;
t560 = sin(pkin(9));
t562 = cos(pkin(9));
t516 = t560 * t537 + t562 * t539;
t580 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t516;
t559 = sin(pkin(10));
t561 = cos(pkin(10));
t564 = sin(qJ(4));
t567 = cos(qJ(4));
t525 = (t559 * t567 + t561 * t564) * qJD(1);
t597 = 2 * qJD(5);
t596 = -pkin(2) - pkin(7);
t595 = pkin(4) * t570;
t594 = mrSges(3,1) - mrSges(4,2);
t593 = Ifges(3,5) - Ifges(4,4);
t592 = -Ifges(3,6) + Ifges(4,5);
t515 = t562 * t537 - t560 * t539;
t576 = -t570 * qJ(3) + qJDD(3) - t515;
t504 = qJDD(1) * t596 + t576;
t499 = t567 * t504;
t588 = qJD(1) * qJD(4);
t541 = t567 * qJDD(1) - t564 * t588;
t556 = -g(3) + qJDD(2);
t484 = qJDD(4) * pkin(4) - t541 * qJ(5) + t499 + (-qJ(5) * t588 - t567 * t595 - t556) * t564;
t496 = t564 * t504 + t567 * t556;
t540 = -t564 * qJDD(1) - t567 * t588;
t589 = qJD(1) * t567;
t543 = qJD(4) * pkin(4) - qJ(5) * t589;
t555 = t564 ^ 2;
t485 = t540 * qJ(5) - qJD(4) * t543 - t555 * t595 + t496;
t480 = t559 * t484 + t561 * t485 - t525 * t597;
t590 = qJD(1) * t564;
t526 = -t559 * t590 + t561 * t589;
t511 = t525 * mrSges(6,1) + t526 * mrSges(6,2);
t517 = t561 * t540 - t559 * t541;
t522 = qJD(4) * mrSges(6,1) - t526 * mrSges(6,3);
t512 = t525 * pkin(5) - t526 * pkin(8);
t569 = qJD(4) ^ 2;
t478 = -t569 * pkin(5) + qJDD(4) * pkin(8) - t525 * t512 + t480;
t487 = -t540 * pkin(4) + qJDD(5) + t543 * t589 + (-qJ(5) * t555 + t596) * t570 + t580;
t518 = t559 * t540 + t561 * t541;
t481 = (qJD(4) * t525 - t518) * pkin(8) + (qJD(4) * t526 - t517) * pkin(5) + t487;
t563 = sin(qJ(6));
t566 = cos(qJ(6));
t475 = -t563 * t478 + t566 * t481;
t519 = t566 * qJD(4) - t563 * t526;
t494 = t519 * qJD(6) + t563 * qJDD(4) + t566 * t518;
t520 = t563 * qJD(4) + t566 * t526;
t497 = -t519 * mrSges(7,1) + t520 * mrSges(7,2);
t524 = qJD(6) + t525;
t501 = -t524 * mrSges(7,2) + t519 * mrSges(7,3);
t514 = qJDD(6) - t517;
t473 = m(7) * t475 + t514 * mrSges(7,1) - t494 * mrSges(7,3) - t520 * t497 + t524 * t501;
t476 = t566 * t478 + t563 * t481;
t493 = -t520 * qJD(6) + t566 * qJDD(4) - t563 * t518;
t502 = t524 * mrSges(7,1) - t520 * mrSges(7,3);
t474 = m(7) * t476 - t514 * mrSges(7,2) + t493 * mrSges(7,3) + t519 * t497 - t524 * t502;
t581 = -t563 * t473 + t566 * t474;
t464 = m(6) * t480 - qJDD(4) * mrSges(6,2) + t517 * mrSges(6,3) - qJD(4) * t522 - t525 * t511 + t581;
t579 = -t561 * t484 + t559 * t485;
t479 = -0.2e1 * qJD(5) * t526 - t579;
t521 = -qJD(4) * mrSges(6,2) - t525 * mrSges(6,3);
t477 = -qJDD(4) * pkin(5) - t569 * pkin(8) + (t597 + t512) * t526 + t579;
t574 = -m(7) * t477 + t493 * mrSges(7,1) - t494 * mrSges(7,2) + t519 * t501 - t520 * t502;
t469 = m(6) * t479 + qJDD(4) * mrSges(6,1) - t518 * mrSges(6,3) + qJD(4) * t521 - t526 * t511 + t574;
t458 = t559 * t464 + t561 * t469;
t495 = -t564 * t556 + t499;
t538 = (mrSges(5,1) * t564 + mrSges(5,2) * t567) * qJD(1);
t542 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t590;
t456 = m(5) * t495 + qJDD(4) * mrSges(5,1) - t541 * mrSges(5,3) + qJD(4) * t542 - t538 * t589 + t458;
t544 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t589;
t582 = t561 * t464 - t559 * t469;
t457 = m(5) * t496 - qJDD(4) * mrSges(5,2) + t540 * mrSges(5,3) - qJD(4) * t544 - t538 * t590 + t582;
t452 = t567 * t456 + t564 * t457;
t506 = -qJDD(1) * pkin(2) + t576;
t575 = -m(4) * t506 + t570 * mrSges(4,3) - t452;
t450 = m(3) * t515 - t570 * mrSges(3,2) + qJDD(1) * t594 + t575;
t505 = t570 * pkin(2) - t580;
t503 = t570 * t596 + t580;
t465 = t566 * t473 + t563 * t474;
t573 = m(6) * t487 - t517 * mrSges(6,1) + t518 * mrSges(6,2) + t525 * t521 + t526 * t522 + t465;
t572 = -m(5) * t503 + t540 * mrSges(5,1) - t541 * mrSges(5,2) - t542 * t590 - t544 * t589 - t573;
t571 = -m(4) * t505 + t570 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t572;
t461 = m(3) * t516 - t570 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t571;
t448 = t562 * t450 + t560 * t461;
t446 = m(2) * t545 + qJDD(1) * mrSges(2,1) - t570 * mrSges(2,2) + t448;
t584 = -t560 * t450 + t562 * t461;
t447 = m(2) * t546 - t570 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t584;
t591 = t568 * t446 + t565 * t447;
t585 = -t565 * t446 + t568 * t447;
t583 = -t564 * t456 + t567 * t457;
t451 = m(4) * t556 + t583;
t577 = m(3) * t556 + t451;
t532 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t567 - Ifges(5,4) * t564) * qJD(1);
t531 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t567 - Ifges(5,2) * t564) * qJD(1);
t530 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t567 - Ifges(5,6) * t564) * qJD(1);
t509 = Ifges(6,1) * t526 - Ifges(6,4) * t525 + (Ifges(6,5) * qJD(4));
t508 = Ifges(6,4) * t526 - Ifges(6,2) * t525 + (Ifges(6,6) * qJD(4));
t507 = Ifges(6,5) * t526 - Ifges(6,6) * t525 + (Ifges(6,3) * qJD(4));
t490 = Ifges(7,1) * t520 + Ifges(7,4) * t519 + Ifges(7,5) * t524;
t489 = Ifges(7,4) * t520 + Ifges(7,2) * t519 + Ifges(7,6) * t524;
t488 = Ifges(7,5) * t520 + Ifges(7,6) * t519 + Ifges(7,3) * t524;
t467 = mrSges(7,2) * t477 - mrSges(7,3) * t475 + Ifges(7,1) * t494 + Ifges(7,4) * t493 + Ifges(7,5) * t514 + t519 * t488 - t524 * t489;
t466 = -mrSges(7,1) * t477 + mrSges(7,3) * t476 + Ifges(7,4) * t494 + Ifges(7,2) * t493 + Ifges(7,6) * t514 - t520 * t488 + t524 * t490;
t454 = -mrSges(6,1) * t487 - mrSges(7,1) * t475 + mrSges(7,2) * t476 + mrSges(6,3) * t480 + Ifges(6,4) * t518 - Ifges(7,5) * t494 + Ifges(6,2) * t517 + Ifges(6,6) * qJDD(4) - Ifges(7,6) * t493 - Ifges(7,3) * t514 - pkin(5) * t465 + qJD(4) * t509 - t520 * t489 + t519 * t490 - t526 * t507;
t453 = mrSges(6,2) * t487 - mrSges(6,3) * t479 + Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * qJDD(4) - pkin(8) * t465 - qJD(4) * t508 - t563 * t466 + t566 * t467 - t525 * t507;
t442 = mrSges(5,2) * t503 - mrSges(5,3) * t495 + Ifges(5,1) * t541 + Ifges(5,4) * t540 + Ifges(5,5) * qJDD(4) - qJ(5) * t458 - qJD(4) * t531 + t561 * t453 - t559 * t454 - t530 * t590;
t441 = -mrSges(5,1) * t503 + mrSges(5,3) * t496 + Ifges(5,4) * t541 + Ifges(5,2) * t540 + Ifges(5,6) * qJDD(4) - pkin(4) * t573 + qJ(5) * t582 + qJD(4) * t532 + t559 * t453 + t561 * t454 - t530 * t589;
t440 = pkin(4) * t458 + pkin(3) * t452 - qJ(3) * t451 + mrSges(6,1) * t479 - mrSges(6,2) * t480 + mrSges(5,1) * t495 - mrSges(5,2) * t496 + mrSges(4,1) * t506 - mrSges(3,3) * t515 + Ifges(6,6) * t517 + Ifges(6,5) * t518 + pkin(5) * t574 + t525 * t509 + t526 * t508 + Ifges(5,6) * t540 + Ifges(5,5) * t541 + t563 * t467 + pkin(8) * t581 + t566 * t466 + t592 * t570 + (mrSges(3,2) - mrSges(4,3)) * t556 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + t593 * qJDD(1) + (t567 * t531 + t564 * t532) * qJD(1);
t439 = -mrSges(4,1) * t505 + mrSges(3,3) * t516 - pkin(2) * t451 - pkin(3) * t572 - pkin(7) * t583 - qJDD(1) * t592 - t567 * t441 - t564 * t442 - t556 * t594 + t570 * t593;
t438 = -mrSges(2,2) * g(3) - mrSges(2,3) * t545 + Ifges(2,5) * qJDD(1) - t570 * Ifges(2,6) - qJ(2) * t448 - t560 * t439 + t562 * t440;
t437 = mrSges(2,1) * g(3) + mrSges(2,3) * t546 + t570 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t577 + qJ(2) * t584 + t562 * t439 + t560 * t440;
t1 = [-m(1) * g(1) + t585; -m(1) * g(2) + t591; (-m(1) - m(2)) * g(3) + t577; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t591 - t565 * t437 + t568 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t585 + t568 * t437 + t565 * t438; pkin(1) * t448 + mrSges(2,1) * t545 - mrSges(2,2) * t546 + qJ(3) * t571 + pkin(2) * t575 - mrSges(3,2) * t516 + t567 * t442 - t564 * t441 - pkin(7) * t452 + mrSges(3,1) * t515 + mrSges(4,2) * t506 - mrSges(4,3) * t505 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
