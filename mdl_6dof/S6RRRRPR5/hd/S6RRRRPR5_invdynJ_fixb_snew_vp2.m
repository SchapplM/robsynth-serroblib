% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 20:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:26:26
% EndTime: 2019-05-07 20:26:38
% DurationCPUTime: 5.01s
% Computational Cost: add. (52090->330), mult. (103953->402), div. (0->0), fcn. (74125->10), ass. (0->134)
t586 = Ifges(5,1) + Ifges(6,1);
t577 = Ifges(5,4) - Ifges(6,5);
t576 = Ifges(5,5) + Ifges(6,4);
t585 = -Ifges(5,2) - Ifges(6,3);
t575 = Ifges(5,6) - Ifges(6,6);
t584 = Ifges(5,3) + Ifges(6,2);
t541 = sin(qJ(3));
t545 = cos(qJ(3));
t546 = cos(qJ(2));
t566 = qJD(1) * t546;
t542 = sin(qJ(2));
t567 = qJD(1) * t542;
t518 = -t541 * t567 + t545 * t566;
t565 = qJD(1) * qJD(2);
t525 = qJDD(1) * t542 + t546 * t565;
t526 = qJDD(1) * t546 - t542 * t565;
t491 = qJD(3) * t518 + t525 * t545 + t526 * t541;
t519 = (t541 * t546 + t542 * t545) * qJD(1);
t537 = qJD(2) + qJD(3);
t540 = sin(qJ(4));
t580 = cos(qJ(4));
t504 = t519 * t540 - t537 * t580;
t536 = qJDD(2) + qJDD(3);
t454 = -t504 * qJD(4) + t491 * t580 + t540 * t536;
t548 = qJD(1) ^ 2;
t543 = sin(qJ(1));
t547 = cos(qJ(1));
t559 = -g(1) * t547 - g(2) * t543;
t521 = -pkin(1) * t548 + qJDD(1) * pkin(7) + t559;
t572 = t542 * t521;
t579 = pkin(2) * t548;
t479 = qJDD(2) * pkin(2) - t525 * pkin(8) - t572 + (pkin(8) * t565 + t542 * t579 - g(3)) * t546;
t507 = -g(3) * t542 + t546 * t521;
t529 = qJD(2) * pkin(2) - pkin(8) * t567;
t538 = t546 ^ 2;
t480 = pkin(8) * t526 - qJD(2) * t529 - t538 * t579 + t507;
t449 = t545 * t479 - t541 * t480;
t502 = -pkin(3) * t518 - pkin(9) * t519;
t535 = t537 ^ 2;
t555 = t536 * pkin(3) + t535 * pkin(9) - t519 * t502 + t449;
t514 = qJD(4) - t518;
t573 = t504 * t514;
t583 = (-t454 + t573) * qJ(5) - t555;
t505 = t519 * t580 + t540 * t537;
t476 = mrSges(6,1) * t504 - mrSges(6,3) * t505;
t490 = -qJD(3) * t519 - t525 * t541 + t545 * t526;
t564 = t543 * g(1) - t547 * g(2);
t557 = -qJDD(1) * pkin(1) - t564;
t492 = -t526 * pkin(2) + t529 * t567 + (-pkin(8) * t538 - pkin(7)) * t548 + t557;
t437 = (-t518 * t537 - t491) * pkin(9) + (t519 * t537 - t490) * pkin(3) + t492;
t450 = t541 * t479 + t545 * t480;
t441 = -pkin(3) * t535 + pkin(9) * t536 + t502 * t518 + t450;
t427 = t437 * t580 - t540 * t441;
t475 = pkin(4) * t504 - qJ(5) * t505;
t489 = qJDD(4) - t490;
t513 = t514 ^ 2;
t425 = -t489 * pkin(4) - t513 * qJ(5) + t505 * t475 + qJDD(5) - t427;
t419 = (-t454 - t573) * pkin(10) + (t504 * t505 - t489) * pkin(5) + t425;
t428 = t540 * t437 + t580 * t441;
t581 = 2 * qJD(5);
t424 = -pkin(4) * t513 + t489 * qJ(5) - t504 * t475 + t514 * t581 + t428;
t453 = qJD(4) * t505 + t491 * t540 - t536 * t580;
t497 = -pkin(5) * t514 - pkin(10) * t505;
t503 = t504 ^ 2;
t420 = -pkin(5) * t503 + pkin(10) * t453 + t497 * t514 + t424;
t539 = sin(qJ(6));
t544 = cos(qJ(6));
t417 = t419 * t544 - t420 * t539;
t469 = t504 * t544 - t505 * t539;
t434 = qJD(6) * t469 + t453 * t539 + t454 * t544;
t470 = t504 * t539 + t505 * t544;
t447 = -mrSges(7,1) * t469 + mrSges(7,2) * t470;
t512 = qJD(6) - t514;
t457 = -mrSges(7,2) * t512 + mrSges(7,3) * t469;
t485 = qJDD(6) - t489;
t414 = m(7) * t417 + mrSges(7,1) * t485 - mrSges(7,3) * t434 - t447 * t470 + t457 * t512;
t418 = t419 * t539 + t420 * t544;
t433 = -qJD(6) * t470 + t453 * t544 - t454 * t539;
t458 = mrSges(7,1) * t512 - mrSges(7,3) * t470;
t415 = m(7) * t418 - mrSges(7,2) * t485 + mrSges(7,3) * t433 + t447 * t469 - t458 * t512;
t406 = t544 * t414 + t539 * t415;
t493 = -mrSges(6,2) * t504 + mrSges(6,3) * t514;
t554 = -m(6) * t425 + t489 * mrSges(6,1) + t514 * t493 - t406;
t405 = t454 * mrSges(6,2) + t505 * t476 - t554;
t443 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t512;
t444 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t512;
t553 = -mrSges(7,1) * t417 + mrSges(7,2) * t418 - Ifges(7,5) * t434 - Ifges(7,6) * t433 - Ifges(7,3) * t485 - t470 * t443 + t469 * t444;
t496 = -mrSges(6,1) * t514 + mrSges(6,2) * t505;
t561 = -t539 * t414 + t544 * t415;
t556 = m(6) * t424 + t489 * mrSges(6,3) + t514 * t496 + t561;
t569 = -t577 * t504 + t586 * t505 + t576 * t514;
t570 = t585 * t504 + t577 * t505 + t575 * t514;
t582 = -t453 * t575 + t454 * t576 + t584 * t489 + t504 * t569 + t505 * t570 + mrSges(5,1) * t427 - mrSges(6,1) * t425 - mrSges(5,2) * t428 + mrSges(6,3) * t424 - pkin(4) * t405 - pkin(5) * t406 + qJ(5) * (-t453 * mrSges(6,2) - t504 * t476 + t556) + t553;
t578 = -mrSges(5,3) - mrSges(6,2);
t501 = -mrSges(4,1) * t518 + mrSges(4,2) * t519;
t509 = mrSges(4,1) * t537 - mrSges(4,3) * t519;
t495 = mrSges(5,1) * t514 - mrSges(5,3) * t505;
t568 = -mrSges(5,1) * t504 - mrSges(5,2) * t505 - t476;
t401 = m(5) * t428 - t489 * mrSges(5,2) + t453 * t578 - t514 * t495 + t504 * t568 + t556;
t494 = -mrSges(5,2) * t514 - mrSges(5,3) * t504;
t403 = m(5) * t427 + t489 * mrSges(5,1) + t454 * t578 + t514 * t494 + t505 * t568 + t554;
t562 = t580 * t401 - t403 * t540;
t396 = m(4) * t450 - mrSges(4,2) * t536 + mrSges(4,3) * t490 + t501 * t518 - t509 * t537 + t562;
t508 = -mrSges(4,2) * t537 + mrSges(4,3) * t518;
t426 = -0.2e1 * qJD(5) * t505 + (t505 * t514 + t453) * pkin(4) + t583;
t422 = -t503 * pkin(10) + (-pkin(4) - pkin(5)) * t453 + (-pkin(4) * t514 + t497 + t581) * t505 - t583;
t558 = -m(7) * t422 + t433 * mrSges(7,1) - t434 * mrSges(7,2) + t469 * t457 - t470 * t458;
t412 = m(6) * t426 + mrSges(6,1) * t453 - t454 * mrSges(6,3) + t493 * t504 - t505 * t496 + t558;
t550 = m(5) * t555 - t453 * mrSges(5,1) - mrSges(5,2) * t454 - t504 * t494 - t495 * t505 - t412;
t410 = m(4) * t449 + mrSges(4,1) * t536 - mrSges(4,3) * t491 - t501 * t519 + t508 * t537 + t550;
t393 = t541 * t396 + t545 * t410;
t398 = t540 * t401 + t580 * t403;
t571 = t575 * t504 - t576 * t505 - t584 * t514;
t563 = t545 * t396 - t410 * t541;
t442 = Ifges(7,5) * t470 + Ifges(7,6) * t469 + Ifges(7,3) * t512;
t407 = -mrSges(7,1) * t422 + mrSges(7,3) * t418 + Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t485 - t442 * t470 + t444 * t512;
t408 = mrSges(7,2) * t422 - mrSges(7,3) * t417 + Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t485 + t442 * t469 - t443 * t512;
t390 = mrSges(5,1) * t555 - mrSges(6,1) * t426 + mrSges(6,2) * t424 + mrSges(5,3) * t428 - pkin(4) * t412 - pkin(5) * t558 - pkin(10) * t561 - t544 * t407 - t539 * t408 + t585 * t453 + t577 * t454 + t575 * t489 + t571 * t505 + t569 * t514;
t392 = -mrSges(5,2) * t555 + mrSges(6,2) * t425 - mrSges(5,3) * t427 - mrSges(6,3) * t426 - pkin(10) * t406 - qJ(5) * t412 - t539 * t407 + t544 * t408 - t577 * t453 + t586 * t454 + t576 * t489 + t571 * t504 - t570 * t514;
t499 = Ifges(4,4) * t519 + Ifges(4,2) * t518 + Ifges(4,6) * t537;
t500 = Ifges(4,1) * t519 + Ifges(4,4) * t518 + Ifges(4,5) * t537;
t552 = mrSges(4,1) * t449 - mrSges(4,2) * t450 + Ifges(4,5) * t491 + Ifges(4,6) * t490 + Ifges(4,3) * t536 + pkin(3) * t550 + pkin(9) * t562 + t580 * t390 + t540 * t392 + t519 * t499 - t500 * t518;
t551 = m(4) * t492 - t490 * mrSges(4,1) + t491 * mrSges(4,2) - t518 * t508 + t519 * t509 + t398;
t528 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t566;
t527 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t567;
t524 = (-mrSges(3,1) * t546 + mrSges(3,2) * t542) * qJD(1);
t520 = -t548 * pkin(7) + t557;
t517 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t542 + Ifges(3,4) * t546) * qJD(1);
t516 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t542 + Ifges(3,2) * t546) * qJD(1);
t506 = -t546 * g(3) - t572;
t498 = Ifges(4,5) * t519 + Ifges(4,6) * t518 + Ifges(4,3) * t537;
t388 = -mrSges(4,1) * t492 + mrSges(4,3) * t450 + Ifges(4,4) * t491 + Ifges(4,2) * t490 + Ifges(4,6) * t536 - pkin(3) * t398 - t519 * t498 + t537 * t500 - t582;
t387 = mrSges(4,2) * t492 - mrSges(4,3) * t449 + Ifges(4,1) * t491 + Ifges(4,4) * t490 + Ifges(4,5) * t536 - pkin(9) * t398 - t540 * t390 + t392 * t580 + t518 * t498 - t537 * t499;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t564 - mrSges(2,2) * t559 + t542 * (mrSges(3,2) * t520 - mrSges(3,3) * t506 + Ifges(3,1) * t525 + Ifges(3,4) * t526 + Ifges(3,5) * qJDD(2) - pkin(8) * t393 - qJD(2) * t516 + t545 * t387 - t541 * t388) + t546 * (-mrSges(3,1) * t520 + mrSges(3,3) * t507 + Ifges(3,4) * t525 + Ifges(3,2) * t526 + Ifges(3,6) * qJDD(2) - pkin(2) * t551 + pkin(8) * t563 + qJD(2) * t517 + t541 * t387 + t545 * t388) + pkin(1) * (-m(3) * t520 + t526 * mrSges(3,1) - t525 * mrSges(3,2) + (-t527 * t542 + t528 * t546) * qJD(1) - t551) + pkin(7) * (t546 * (m(3) * t507 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t526 - qJD(2) * t527 + t524 * t566 + t563) - t542 * (m(3) * t506 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t525 + qJD(2) * t528 - t524 * t567 + t393)); t552 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t506 - mrSges(3,2) * t507 + Ifges(3,5) * t525 + Ifges(3,6) * t526 + pkin(2) * t393 + (t542 * t516 - t546 * t517) * qJD(1); t552; t582; t405; -t553;];
tauJ  = t1;
