% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP9
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:06:02
% EndTime: 2019-12-31 20:06:08
% DurationCPUTime: 4.10s
% Computational Cost: add. (38778->287), mult. (83170->352), div. (0->0), fcn. (54413->8), ass. (0->110)
t586 = Ifges(5,1) + Ifges(6,1);
t580 = Ifges(5,4) - Ifges(6,5);
t585 = -Ifges(5,5) - Ifges(6,4);
t584 = Ifges(5,2) + Ifges(6,3);
t578 = Ifges(5,6) - Ifges(6,6);
t583 = -Ifges(5,3) - Ifges(6,2);
t582 = cos(qJ(4));
t581 = -mrSges(5,3) - mrSges(6,2);
t555 = sin(qJ(1));
t557 = cos(qJ(1));
t545 = -t557 * g(1) - t555 * g(2);
t559 = qJD(1) ^ 2;
t529 = -t559 * pkin(1) + qJDD(1) * pkin(6) + t545;
t554 = sin(qJ(2));
t556 = cos(qJ(2));
t512 = -t554 * g(3) + t556 * t529;
t539 = (-mrSges(3,1) * t556 + mrSges(3,2) * t554) * qJD(1);
t570 = qJD(1) * qJD(2);
t548 = t554 * t570;
t541 = t556 * qJDD(1) - t548;
t572 = qJD(1) * t554;
t542 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t572;
t544 = t555 * g(1) - t557 * g(2);
t528 = -qJDD(1) * pkin(1) - t559 * pkin(6) - t544;
t568 = t556 * t570;
t540 = t554 * qJDD(1) + t568;
t492 = (-t540 - t568) * qJ(3) + (-t541 + t548) * pkin(2) + t528;
t538 = (-pkin(2) * t556 - qJ(3) * t554) * qJD(1);
t558 = qJD(2) ^ 2;
t571 = t556 * qJD(1);
t497 = -t558 * pkin(2) + qJDD(2) * qJ(3) + t538 * t571 + t512;
t551 = sin(pkin(8));
t552 = cos(pkin(8));
t535 = t551 * qJD(2) + t552 * t572;
t470 = -0.2e1 * qJD(3) * t535 + t552 * t492 - t551 * t497;
t517 = t551 * qJDD(2) + t552 * t540;
t534 = t552 * qJD(2) - t551 * t572;
t467 = (-t534 * t571 - t517) * pkin(7) + (t534 * t535 - t541) * pkin(3) + t470;
t471 = 0.2e1 * qJD(3) * t534 + t551 * t492 + t552 * t497;
t516 = t552 * qJDD(2) - t551 * t540;
t518 = -pkin(3) * t571 - t535 * pkin(7);
t533 = t534 ^ 2;
t469 = -t533 * pkin(3) + t516 * pkin(7) + t518 * t571 + t471;
t553 = sin(qJ(4));
t465 = t553 * t467 + t582 * t469;
t508 = t553 * t534 + t582 * t535;
t476 = t508 * qJD(4) - t582 * t516 + t553 * t517;
t547 = qJD(4) - t571;
t499 = t547 * mrSges(5,1) - t508 * mrSges(5,3);
t507 = -t582 * t534 + t553 * t535;
t537 = qJDD(4) - t541;
t487 = t507 * pkin(4) - t508 * qJ(5);
t546 = t547 ^ 2;
t460 = -t546 * pkin(4) + t537 * qJ(5) + 0.2e1 * qJD(5) * t547 - t507 * t487 + t465;
t500 = -t547 * mrSges(6,1) + t508 * mrSges(6,2);
t569 = m(6) * t460 + t537 * mrSges(6,3) + t547 * t500;
t488 = t507 * mrSges(6,1) - t508 * mrSges(6,3);
t573 = -t507 * mrSges(5,1) - t508 * mrSges(5,2) - t488;
t453 = m(5) * t465 - t537 * mrSges(5,2) + t581 * t476 - t547 * t499 + t573 * t507 + t569;
t464 = t582 * t467 - t553 * t469;
t477 = -t507 * qJD(4) + t553 * t516 + t582 * t517;
t498 = -t547 * mrSges(5,2) - t507 * mrSges(5,3);
t461 = -t537 * pkin(4) - t546 * qJ(5) + t508 * t487 + qJDD(5) - t464;
t501 = -t507 * mrSges(6,2) + t547 * mrSges(6,3);
t563 = -m(6) * t461 + t537 * mrSges(6,1) + t547 * t501;
t455 = m(5) * t464 + t537 * mrSges(5,1) + t581 * t477 + t547 * t498 + t573 * t508 + t563;
t450 = t553 * t453 + t582 * t455;
t509 = -t534 * mrSges(4,1) + t535 * mrSges(4,2);
t514 = mrSges(4,2) * t571 + t534 * mrSges(4,3);
t446 = m(4) * t470 - t541 * mrSges(4,1) - t517 * mrSges(4,3) - t535 * t509 - t514 * t571 + t450;
t515 = -mrSges(4,1) * t571 - t535 * mrSges(4,3);
t564 = t582 * t453 - t553 * t455;
t447 = m(4) * t471 + t541 * mrSges(4,2) + t516 * mrSges(4,3) + t534 * t509 + t515 * t571 + t564;
t565 = -t551 * t446 + t552 * t447;
t443 = m(3) * t512 - qJDD(2) * mrSges(3,2) + t541 * mrSges(3,3) - qJD(2) * t542 + t539 * t571 + t565;
t511 = -t556 * g(3) - t554 * t529;
t543 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t571;
t496 = -qJDD(2) * pkin(2) - t558 * qJ(3) + t538 * t572 + qJDD(3) - t511;
t472 = -t516 * pkin(3) - t533 * pkin(7) + t535 * t518 + t496;
t463 = -0.2e1 * qJD(5) * t508 + (t507 * t547 - t477) * qJ(5) + (t508 * t547 + t476) * pkin(4) + t472;
t458 = m(6) * t463 + t476 * mrSges(6,1) - t477 * mrSges(6,3) - t508 * t500 + t507 * t501;
t562 = m(5) * t472 + t476 * mrSges(5,1) + t477 * mrSges(5,2) + t507 * t498 + t508 * t499 + t458;
t560 = -m(4) * t496 + t516 * mrSges(4,1) - t517 * mrSges(4,2) + t534 * t514 - t535 * t515 - t562;
t457 = m(3) * t511 + qJDD(2) * mrSges(3,1) - t540 * mrSges(3,3) + qJD(2) * t543 - t539 * t572 + t560;
t566 = t556 * t443 - t554 * t457;
t437 = m(2) * t545 - t559 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t566;
t444 = t552 * t446 + t551 * t447;
t561 = -m(3) * t528 + t541 * mrSges(3,1) - t540 * mrSges(3,2) - t542 * t572 + t543 * t571 - t444;
t440 = m(2) * t544 + qJDD(1) * mrSges(2,1) - t559 * mrSges(2,2) + t561;
t577 = t555 * t437 + t557 * t440;
t438 = t554 * t443 + t556 * t457;
t576 = t584 * t507 - t580 * t508 - t578 * t547;
t575 = t578 * t507 + t585 * t508 + t583 * t547;
t574 = -t580 * t507 + t586 * t508 - t585 * t547;
t567 = t557 * t437 - t555 * t440;
t527 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t554 + Ifges(3,4) * t556) * qJD(1);
t526 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t554 + Ifges(3,2) * t556) * qJD(1);
t525 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t554 + Ifges(3,6) * t556) * qJD(1);
t504 = Ifges(4,1) * t535 + Ifges(4,4) * t534 - Ifges(4,5) * t571;
t503 = Ifges(4,4) * t535 + Ifges(4,2) * t534 - Ifges(4,6) * t571;
t502 = Ifges(4,5) * t535 + Ifges(4,6) * t534 - Ifges(4,3) * t571;
t449 = mrSges(5,2) * t472 + mrSges(6,2) * t461 - mrSges(5,3) * t464 - mrSges(6,3) * t463 - qJ(5) * t458 - t580 * t476 + t586 * t477 + t575 * t507 - t537 * t585 + t576 * t547;
t448 = -mrSges(5,1) * t472 - mrSges(6,1) * t463 + mrSges(6,2) * t460 + mrSges(5,3) * t465 - pkin(4) * t458 - t584 * t476 + t580 * t477 + t575 * t508 + t578 * t537 + t574 * t547;
t434 = mrSges(4,2) * t496 - mrSges(4,3) * t470 + Ifges(4,1) * t517 + Ifges(4,4) * t516 - Ifges(4,5) * t541 - pkin(7) * t450 - t553 * t448 + t582 * t449 + t534 * t502 + t503 * t571;
t433 = -mrSges(4,1) * t496 + mrSges(4,3) * t471 + Ifges(4,4) * t517 + Ifges(4,2) * t516 - Ifges(4,6) * t541 - pkin(3) * t562 + pkin(7) * t564 + t582 * t448 + t553 * t449 - t535 * t502 - t504 * t571;
t432 = t583 * t537 - pkin(4) * t563 + mrSges(3,3) * t512 + t534 * t504 - t535 * t503 + (Ifges(4,3) + Ifges(3,2)) * t541 - Ifges(4,6) * t516 - Ifges(4,5) * t517 + qJD(2) * t527 - mrSges(3,1) * t528 - pkin(2) * t444 - mrSges(6,3) * t460 + mrSges(6,1) * t461 + Ifges(3,6) * qJDD(2) - pkin(3) * t450 + Ifges(3,4) * t540 - mrSges(5,1) * t464 + mrSges(5,2) * t465 - mrSges(4,1) * t470 + mrSges(4,2) * t471 + (pkin(4) * mrSges(6,2) + t585) * t477 - qJ(5) * t569 - t525 * t572 + (qJ(5) * t488 - t574) * t507 + (pkin(4) * t488 + t576) * t508 + (qJ(5) * mrSges(6,2) + t578) * t476;
t431 = mrSges(3,2) * t528 - mrSges(3,3) * t511 + Ifges(3,1) * t540 + Ifges(3,4) * t541 + Ifges(3,5) * qJDD(2) - qJ(3) * t444 - qJD(2) * t526 - t551 * t433 + t552 * t434 + t525 * t571;
t430 = Ifges(2,6) * qJDD(1) + t559 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t545 - Ifges(3,5) * t540 - Ifges(3,6) * t541 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t511 + mrSges(3,2) * t512 - t551 * t434 - t552 * t433 - pkin(2) * t560 - qJ(3) * t565 - pkin(1) * t438 + (-t554 * t526 + t556 * t527) * qJD(1);
t429 = -mrSges(2,2) * g(3) - mrSges(2,3) * t544 + Ifges(2,5) * qJDD(1) - t559 * Ifges(2,6) - pkin(6) * t438 + t556 * t431 - t554 * t432;
t1 = [-m(1) * g(1) + t567; -m(1) * g(2) + t577; (-m(1) - m(2)) * g(3) + t438; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t577 + t557 * t429 - t555 * t430; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t567 + t555 * t429 + t557 * t430; -mrSges(1,1) * g(2) + mrSges(2,1) * t544 + mrSges(1,2) * g(1) - mrSges(2,2) * t545 + Ifges(2,3) * qJDD(1) + pkin(1) * t561 + pkin(6) * t566 + t554 * t431 + t556 * t432;];
tauB = t1;
