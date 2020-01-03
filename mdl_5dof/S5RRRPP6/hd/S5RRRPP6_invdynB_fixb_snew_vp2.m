% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:21
% EndTime: 2019-12-31 21:00:27
% DurationCPUTime: 3.91s
% Computational Cost: add. (41669->288), mult. (85416->352), div. (0->0), fcn. (56037->8), ass. (0->112)
t588 = Ifges(5,1) + Ifges(6,1);
t583 = Ifges(5,4) - Ifges(6,5);
t582 = Ifges(5,5) + Ifges(6,4);
t587 = Ifges(5,2) + Ifges(6,3);
t586 = -Ifges(6,2) - Ifges(5,3);
t581 = Ifges(5,6) - Ifges(6,6);
t585 = -2 * qJD(4);
t584 = -mrSges(5,3) - mrSges(6,2);
t580 = cos(pkin(8));
t554 = sin(qJ(1));
t557 = cos(qJ(1));
t546 = -t557 * g(1) - t554 * g(2);
t559 = qJD(1) ^ 2;
t530 = -t559 * pkin(1) + qJDD(1) * pkin(6) + t546;
t553 = sin(qJ(2));
t556 = cos(qJ(2));
t521 = -t553 * g(3) + t556 * t530;
t539 = (-mrSges(3,1) * t556 + mrSges(3,2) * t553) * qJD(1);
t572 = qJD(1) * qJD(2);
t570 = t553 * t572;
t542 = t556 * qJDD(1) - t570;
t574 = qJD(1) * t553;
t543 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t574;
t545 = t554 * g(1) - t557 * g(2);
t529 = -qJDD(1) * pkin(1) - t559 * pkin(6) - t545;
t569 = t556 * t572;
t541 = t553 * qJDD(1) + t569;
t494 = (-t541 - t569) * pkin(7) + (-t542 + t570) * pkin(2) + t529;
t540 = (-pkin(2) * t556 - pkin(7) * t553) * qJD(1);
t558 = qJD(2) ^ 2;
t573 = t556 * qJD(1);
t499 = -t558 * pkin(2) + qJDD(2) * pkin(7) + t540 * t573 + t521;
t552 = sin(qJ(3));
t555 = cos(qJ(3));
t473 = t555 * t494 - t552 * t499;
t537 = t555 * qJD(2) - t552 * t574;
t513 = t537 * qJD(3) + t552 * qJDD(2) + t555 * t541;
t536 = qJDD(3) - t542;
t538 = t552 * qJD(2) + t555 * t574;
t548 = qJD(3) - t573;
t469 = (t537 * t548 - t513) * qJ(4) + (t537 * t538 + t536) * pkin(3) + t473;
t474 = t552 * t494 + t555 * t499;
t512 = -t538 * qJD(3) + t555 * qJDD(2) - t552 * t541;
t518 = t548 * pkin(3) - t538 * qJ(4);
t535 = t537 ^ 2;
t471 = -t535 * pkin(3) + t512 * qJ(4) - t548 * t518 + t474;
t551 = sin(pkin(8));
t514 = -t580 * t537 + t551 * t538;
t465 = t551 * t469 + t580 * t471 + t514 * t585;
t484 = -t580 * t512 + t551 * t513;
t515 = t551 * t537 + t580 * t538;
t501 = t548 * mrSges(5,1) - t515 * mrSges(5,3);
t489 = t514 * pkin(4) - t515 * qJ(5);
t547 = t548 ^ 2;
t462 = -t547 * pkin(4) + t536 * qJ(5) + 0.2e1 * qJD(5) * t548 - t514 * t489 + t465;
t502 = -t548 * mrSges(6,1) + t515 * mrSges(6,2);
t571 = m(6) * t462 + t536 * mrSges(6,3) + t548 * t502;
t490 = t514 * mrSges(6,1) - t515 * mrSges(6,3);
t575 = -t514 * mrSges(5,1) - t515 * mrSges(5,2) - t490;
t455 = m(5) * t465 - t536 * mrSges(5,2) + t584 * t484 - t548 * t501 + t575 * t514 + t571;
t563 = t580 * t469 - t551 * t471;
t464 = t515 * t585 + t563;
t485 = t551 * t512 + t580 * t513;
t500 = -t548 * mrSges(5,2) - t514 * mrSges(5,3);
t463 = -t536 * pkin(4) - t547 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t489) * t515 - t563;
t503 = -t514 * mrSges(6,2) + t548 * mrSges(6,3);
t564 = -m(6) * t463 + t536 * mrSges(6,1) + t548 * t503;
t457 = m(5) * t464 + t536 * mrSges(5,1) + t584 * t485 + t548 * t500 + t575 * t515 + t564;
t452 = t551 * t455 + t580 * t457;
t516 = -t537 * mrSges(4,1) + t538 * mrSges(4,2);
t517 = -t548 * mrSges(4,2) + t537 * mrSges(4,3);
t448 = m(4) * t473 + t536 * mrSges(4,1) - t513 * mrSges(4,3) - t538 * t516 + t548 * t517 + t452;
t519 = t548 * mrSges(4,1) - t538 * mrSges(4,3);
t565 = t580 * t455 - t551 * t457;
t449 = m(4) * t474 - t536 * mrSges(4,2) + t512 * mrSges(4,3) + t537 * t516 - t548 * t519 + t565;
t566 = -t552 * t448 + t555 * t449;
t445 = m(3) * t521 - qJDD(2) * mrSges(3,2) + t542 * mrSges(3,3) - qJD(2) * t543 + t539 * t573 + t566;
t520 = -t556 * g(3) - t553 * t530;
t544 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t573;
t498 = -qJDD(2) * pkin(2) - t558 * pkin(7) + t540 * t574 - t520;
t472 = -t512 * pkin(3) - t535 * qJ(4) + t538 * t518 + qJDD(4) + t498;
t467 = -0.2e1 * qJD(5) * t515 + (t514 * t548 - t485) * qJ(5) + (t515 * t548 + t484) * pkin(4) + t472;
t460 = m(6) * t467 + t484 * mrSges(6,1) - t485 * mrSges(6,3) - t515 * t502 + t514 * t503;
t562 = m(5) * t472 + t484 * mrSges(5,1) + t485 * mrSges(5,2) + t514 * t500 + t515 * t501 + t460;
t560 = -m(4) * t498 + t512 * mrSges(4,1) - t513 * mrSges(4,2) + t537 * t517 - t538 * t519 - t562;
t459 = m(3) * t520 + qJDD(2) * mrSges(3,1) - t541 * mrSges(3,3) + qJD(2) * t544 - t539 * t574 + t560;
t567 = t556 * t445 - t553 * t459;
t439 = m(2) * t546 - t559 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t567;
t446 = t555 * t448 + t552 * t449;
t561 = -m(3) * t529 + t542 * mrSges(3,1) - t541 * mrSges(3,2) - t543 * t574 + t544 * t573 - t446;
t442 = m(2) * t545 + qJDD(1) * mrSges(2,1) - t559 * mrSges(2,2) + t561;
t579 = t554 * t439 + t557 * t442;
t440 = t553 * t445 + t556 * t459;
t578 = t587 * t514 - t583 * t515 - t581 * t548;
t577 = t581 * t514 - t582 * t515 + t586 * t548;
t576 = -t583 * t514 + t588 * t515 + t582 * t548;
t568 = t557 * t439 - t554 * t442;
t528 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t553 + Ifges(3,4) * t556) * qJD(1);
t527 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t553 + Ifges(3,2) * t556) * qJD(1);
t526 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t553 + Ifges(3,6) * t556) * qJD(1);
t507 = Ifges(4,1) * t538 + Ifges(4,4) * t537 + Ifges(4,5) * t548;
t506 = Ifges(4,4) * t538 + Ifges(4,2) * t537 + Ifges(4,6) * t548;
t505 = Ifges(4,5) * t538 + Ifges(4,6) * t537 + Ifges(4,3) * t548;
t451 = mrSges(5,2) * t472 + mrSges(6,2) * t463 - mrSges(5,3) * t464 - mrSges(6,3) * t467 - qJ(5) * t460 - t583 * t484 + t588 * t485 + t577 * t514 + t582 * t536 + t578 * t548;
t450 = -mrSges(5,1) * t472 - mrSges(6,1) * t467 + mrSges(6,2) * t462 + mrSges(5,3) * t465 - pkin(4) * t460 - t587 * t484 + t583 * t485 + t577 * t515 + t581 * t536 + t576 * t548;
t436 = mrSges(4,2) * t498 - mrSges(4,3) * t473 + Ifges(4,1) * t513 + Ifges(4,4) * t512 + Ifges(4,5) * t536 - qJ(4) * t452 - t551 * t450 + t580 * t451 + t537 * t505 - t548 * t506;
t435 = -mrSges(4,1) * t498 + mrSges(4,3) * t474 + Ifges(4,4) * t513 + Ifges(4,2) * t512 + Ifges(4,6) * t536 - pkin(3) * t562 + qJ(4) * t565 + t580 * t450 + t551 * t451 - t538 * t505 + t548 * t507;
t434 = -t526 * t574 + (qJ(5) * mrSges(6,2) + t581) * t484 + (pkin(4) * mrSges(6,2) - t582) * t485 + (-Ifges(4,3) + t586) * t536 + t537 * t507 - t538 * t506 + Ifges(3,4) * t541 + Ifges(3,2) * t542 + qJD(2) * t528 - mrSges(3,1) * t529 + mrSges(3,3) * t521 - Ifges(4,5) * t513 - Ifges(4,6) * t512 - mrSges(4,1) * t473 + mrSges(4,2) * t474 - mrSges(5,1) * t464 + mrSges(5,2) * t465 - mrSges(6,3) * t462 + mrSges(6,1) * t463 - pkin(3) * t452 - pkin(2) * t446 + Ifges(3,6) * qJDD(2) + (qJ(5) * t490 - t576) * t514 + (pkin(4) * t490 + t578) * t515 - qJ(5) * t571 - pkin(4) * t564;
t433 = mrSges(3,2) * t529 - mrSges(3,3) * t520 + Ifges(3,1) * t541 + Ifges(3,4) * t542 + Ifges(3,5) * qJDD(2) - pkin(7) * t446 - qJD(2) * t527 - t552 * t435 + t555 * t436 + t526 * t573;
t432 = Ifges(2,6) * qJDD(1) + t559 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t546 - Ifges(3,5) * t541 - Ifges(3,6) * t542 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t520 + mrSges(3,2) * t521 - t552 * t436 - t555 * t435 - pkin(2) * t560 - pkin(7) * t566 - pkin(1) * t440 + (-t553 * t527 + t556 * t528) * qJD(1);
t431 = -mrSges(2,2) * g(3) - mrSges(2,3) * t545 + Ifges(2,5) * qJDD(1) - t559 * Ifges(2,6) - pkin(6) * t440 + t556 * t433 - t553 * t434;
t1 = [-m(1) * g(1) + t568; -m(1) * g(2) + t579; (-m(1) - m(2)) * g(3) + t440; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t579 + t557 * t431 - t554 * t432; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t568 + t554 * t431 + t557 * t432; -mrSges(1,1) * g(2) + mrSges(2,1) * t545 + mrSges(1,2) * g(1) - mrSges(2,2) * t546 + Ifges(2,3) * qJDD(1) + pkin(1) * t561 + pkin(6) * t567 + t553 * t433 + t556 * t434;];
tauB = t1;
