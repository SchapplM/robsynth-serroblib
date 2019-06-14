% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 14:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:21:32
% EndTime: 2019-05-06 14:21:43
% DurationCPUTime: 6.53s
% Computational Cost: add. (59949->343), mult. (157567->428), div. (0->0), fcn. (122788->12), ass. (0->145)
t597 = Ifges(5,1) + Ifges(6,2);
t590 = Ifges(5,4) + Ifges(6,6);
t589 = Ifges(5,5) - Ifges(6,4);
t596 = -Ifges(5,2) - Ifges(6,3);
t588 = Ifges(5,6) - Ifges(6,5);
t595 = Ifges(5,3) + Ifges(6,1);
t545 = sin(pkin(6));
t550 = sin(qJ(2));
t553 = cos(qJ(2));
t575 = qJD(1) * qJD(2);
t533 = (qJDD(1) * t550 + t553 * t575) * t545;
t547 = cos(pkin(6));
t539 = qJDD(1) * t547 + qJDD(2);
t540 = qJD(1) * t547 + qJD(2);
t555 = qJD(1) ^ 2;
t551 = sin(qJ(1));
t554 = cos(qJ(1));
t566 = -g(1) * t554 - g(2) * t551;
t591 = pkin(8) * t545;
t531 = -pkin(1) * t555 + qJDD(1) * t591 + t566;
t571 = g(1) * t551 - g(2) * t554;
t530 = qJDD(1) * pkin(1) + t555 * t591 + t571;
t585 = t530 * t547;
t567 = -t550 * t531 + t553 * t585;
t584 = t545 ^ 2 * t555;
t459 = t539 * pkin(2) - t533 * qJ(3) + (pkin(2) * t550 * t584 + (qJ(3) * qJD(1) * t540 - g(3)) * t545) * t553 + t567;
t583 = t545 * t550;
t498 = -g(3) * t583 + t531 * t553 + t550 * t585;
t577 = qJD(1) * t545;
t573 = t550 * t577;
t527 = pkin(2) * t540 - qJ(3) * t573;
t534 = (qJDD(1) * t553 - t550 * t575) * t545;
t574 = t553 ^ 2 * t584;
t462 = -pkin(2) * t574 + qJ(3) * t534 - t527 * t540 + t498;
t544 = sin(pkin(11));
t546 = cos(pkin(11));
t525 = (t544 * t553 + t546 * t550) * t577;
t439 = -0.2e1 * qJD(3) * t525 + t546 * t459 - t462 * t544;
t549 = sin(qJ(4));
t592 = cos(qJ(4));
t509 = t525 * t549 - t540 * t592;
t572 = t553 * t577;
t524 = -t544 * t573 + t546 * t572;
t523 = qJD(4) - t524;
t487 = mrSges(6,1) * t509 - mrSges(6,3) * t523;
t504 = -t533 * t544 + t534 * t546;
t503 = qJDD(4) - t504;
t440 = 0.2e1 * qJD(3) * t524 + t459 * t544 + t462 * t546;
t500 = -pkin(3) * t524 - pkin(9) * t525;
t538 = t540 ^ 2;
t438 = -pkin(3) * t538 + pkin(9) * t539 + t500 * t524 + t440;
t516 = -t547 * g(3) - t545 * t530;
t480 = -t534 * pkin(2) - qJ(3) * t574 + t527 * t573 + qJDD(3) + t516;
t505 = t533 * t546 + t534 * t544;
t443 = (-t524 * t540 - t505) * pkin(9) + (t525 * t540 - t504) * pkin(3) + t480;
t433 = -t438 * t549 + t443 * t592;
t510 = t525 * t592 + t540 * t549;
t481 = pkin(4) * t509 - qJ(5) * t510;
t522 = t523 ^ 2;
t431 = -t503 * pkin(4) - t522 * qJ(5) + t481 * t510 + qJDD(5) - t433;
t477 = -qJD(4) * t509 + t505 * t592 + t539 * t549;
t586 = t509 * t523;
t426 = (t509 * t510 - t503) * pkin(10) + (t477 + t586) * pkin(5) + t431;
t476 = qJD(4) * t510 + t505 * t549 - t539 * t592;
t491 = pkin(5) * t510 - pkin(10) * t523;
t508 = t509 ^ 2;
t437 = -t539 * pkin(3) - t538 * pkin(9) + t500 * t525 - t439;
t593 = -2 * qJD(5);
t558 = (-t477 + t586) * qJ(5) + t437 + (pkin(4) * t523 + t593) * t510;
t429 = -t508 * pkin(5) - t510 * t491 + (pkin(4) + pkin(10)) * t476 + t558;
t548 = sin(qJ(6));
t552 = cos(qJ(6));
t424 = t426 * t552 - t429 * t548;
t485 = t509 * t552 - t523 * t548;
t448 = qJD(6) * t485 + t476 * t548 + t503 * t552;
t486 = t509 * t548 + t523 * t552;
t453 = -mrSges(7,1) * t485 + mrSges(7,2) * t486;
t507 = qJD(6) + t510;
t460 = -mrSges(7,2) * t507 + mrSges(7,3) * t485;
t475 = qJDD(6) + t477;
t421 = m(7) * t424 + mrSges(7,1) * t475 - mrSges(7,3) * t448 - t453 * t486 + t460 * t507;
t425 = t426 * t548 + t429 * t552;
t447 = -qJD(6) * t486 + t476 * t552 - t503 * t548;
t461 = mrSges(7,1) * t507 - mrSges(7,3) * t486;
t422 = m(7) * t425 - mrSges(7,2) * t475 + mrSges(7,3) * t447 + t453 * t485 - t461 * t507;
t412 = t552 * t421 + t548 * t422;
t483 = -mrSges(6,2) * t509 - mrSges(6,3) * t510;
t563 = -m(6) * t431 - mrSges(6,1) * t477 - t483 * t510 - t412;
t410 = t503 * mrSges(6,2) + t523 * t487 - t563;
t434 = t438 * t592 + t443 * t549;
t562 = -t522 * pkin(4) + t503 * qJ(5) - t509 * t481 + t434;
t428 = -t476 * pkin(5) - t508 * pkin(10) + ((2 * qJD(5)) + t491) * t523 + t562;
t449 = Ifges(7,5) * t486 + Ifges(7,6) * t485 + Ifges(7,3) * t507;
t451 = Ifges(7,1) * t486 + Ifges(7,4) * t485 + Ifges(7,5) * t507;
t413 = -mrSges(7,1) * t428 + mrSges(7,3) * t425 + Ifges(7,4) * t448 + Ifges(7,2) * t447 + Ifges(7,6) * t475 - t449 * t486 + t451 * t507;
t450 = Ifges(7,4) * t486 + Ifges(7,2) * t485 + Ifges(7,6) * t507;
t414 = mrSges(7,2) * t428 - mrSges(7,3) * t424 + Ifges(7,1) * t448 + Ifges(7,4) * t447 + Ifges(7,5) * t475 + t449 * t485 - t450 * t507;
t430 = t523 * t593 - t562;
t488 = mrSges(6,1) * t510 + mrSges(6,2) * t523;
t564 = -m(7) * t428 + t447 * mrSges(7,1) - mrSges(7,2) * t448 + t485 * t460 - t461 * t486;
t559 = -m(6) * t430 + mrSges(6,3) * t503 + t488 * t523 - t564;
t578 = -t509 * t590 + t510 * t597 + t523 * t589;
t579 = t509 * t596 + t510 * t590 + t523 * t588;
t594 = -t476 * t588 + t477 * t589 + t595 * t503 + t509 * t578 + t510 * t579 + mrSges(5,1) * t433 - mrSges(5,2) * t434 + mrSges(6,2) * t431 - mrSges(6,3) * t430 - pkin(4) * t410 - pkin(10) * t412 + qJ(5) * (-t476 * mrSges(6,1) - t509 * t483 + t559) - t548 * t413 + t552 * t414;
t582 = t545 * t553;
t499 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t512 = mrSges(4,1) * t540 - mrSges(4,3) * t525;
t482 = mrSges(5,1) * t509 + mrSges(5,2) * t510;
t489 = -mrSges(5,2) * t523 - mrSges(5,3) * t509;
t409 = m(5) * t433 - t477 * mrSges(5,3) - t510 * t482 + (-t487 + t489) * t523 + (mrSges(5,1) - mrSges(6,2)) * t503 + t563;
t490 = mrSges(5,1) * t523 - mrSges(5,3) * t510;
t417 = m(5) * t434 - t503 * mrSges(5,2) - t523 * t490 + (-t482 - t483) * t509 + (-mrSges(5,3) - mrSges(6,1)) * t476 + t559;
t569 = -t409 * t549 + t417 * t592;
t403 = m(4) * t440 - mrSges(4,2) * t539 + mrSges(4,3) * t504 + t499 * t524 - t512 * t540 + t569;
t511 = -mrSges(4,2) * t540 + mrSges(4,3) * t524;
t432 = t476 * pkin(4) + t558;
t581 = -t421 * t548 + t422 * t552;
t565 = -m(6) * t432 + mrSges(6,2) * t476 + t487 * t509 - t581;
t557 = -m(5) * t437 - t509 * t489 - t476 * mrSges(5,1) + (t488 - t490) * t510 + (-mrSges(5,2) + mrSges(6,3)) * t477 + t565;
t407 = m(4) * t439 + t539 * mrSges(4,1) - t505 * mrSges(4,3) - t525 * t499 + t540 * t511 + t557;
t400 = t403 * t544 + t407 * t546;
t405 = t409 * t592 + t417 * t549;
t580 = t509 * t588 - t510 * t589 - t523 * t595;
t570 = t403 * t546 - t407 * t544;
t561 = -m(4) * t480 + t504 * mrSges(4,1) - mrSges(4,2) * t505 + t524 * t511 - t512 * t525 - t405;
t560 = mrSges(7,1) * t424 - mrSges(7,2) * t425 + Ifges(7,5) * t448 + Ifges(7,6) * t447 + Ifges(7,3) * t475 + t450 * t486 - t485 * t451;
t532 = (-mrSges(3,1) * t553 + mrSges(3,2) * t550) * t577;
t529 = -mrSges(3,2) * t540 + mrSges(3,3) * t572;
t528 = mrSges(3,1) * t540 - mrSges(3,3) * t573;
t515 = Ifges(3,5) * t540 + (Ifges(3,1) * t550 + Ifges(3,4) * t553) * t577;
t514 = Ifges(3,6) * t540 + (Ifges(3,4) * t550 + Ifges(3,2) * t553) * t577;
t513 = Ifges(3,3) * t540 + (Ifges(3,5) * t550 + Ifges(3,6) * t553) * t577;
t497 = -g(3) * t582 + t567;
t495 = Ifges(4,1) * t525 + Ifges(4,4) * t524 + Ifges(4,5) * t540;
t494 = Ifges(4,4) * t525 + Ifges(4,2) * t524 + Ifges(4,6) * t540;
t493 = Ifges(4,5) * t525 + Ifges(4,6) * t524 + Ifges(4,3) * t540;
t411 = -t477 * mrSges(6,3) - t510 * t488 - t565;
t399 = m(3) * t498 - mrSges(3,2) * t539 + mrSges(3,3) * t534 - t528 * t540 + t532 * t572 + t570;
t398 = m(3) * t497 + mrSges(3,1) * t539 - mrSges(3,3) * t533 + t529 * t540 - t532 * t573 + t400;
t397 = mrSges(6,1) * t431 + mrSges(5,2) * t437 - mrSges(5,3) * t433 - mrSges(6,3) * t432 + pkin(5) * t412 - qJ(5) * t411 - t476 * t590 + t477 * t597 + t503 * t589 + t509 * t580 - t523 * t579 + t560;
t396 = -mrSges(5,1) * t437 - mrSges(6,1) * t430 + mrSges(6,2) * t432 + mrSges(5,3) * t434 - pkin(4) * t411 - pkin(5) * t564 - pkin(10) * t581 - t552 * t413 - t548 * t414 + t476 * t596 + t477 * t590 + t503 * t588 + t510 * t580 + t523 * t578;
t395 = -mrSges(4,1) * t480 + mrSges(4,3) * t440 + Ifges(4,4) * t505 + Ifges(4,2) * t504 + Ifges(4,6) * t539 - pkin(3) * t405 - t525 * t493 + t540 * t495 - t594;
t394 = mrSges(4,2) * t480 - mrSges(4,3) * t439 + Ifges(4,1) * t505 + Ifges(4,4) * t504 + Ifges(4,5) * t539 - pkin(9) * t405 - t396 * t549 + t397 * t592 + t493 * t524 - t494 * t540;
t393 = Ifges(3,5) * t533 + Ifges(3,6) * t534 + mrSges(3,1) * t497 - mrSges(3,2) * t498 + Ifges(4,5) * t505 + Ifges(4,6) * t504 + t525 * t494 - t524 * t495 + mrSges(4,1) * t439 - mrSges(4,2) * t440 + t549 * t397 + t592 * t396 + pkin(3) * t557 + pkin(9) * t569 + pkin(2) * t400 + (Ifges(3,3) + Ifges(4,3)) * t539 + (t514 * t550 - t515 * t553) * t577;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t571 - mrSges(2,2) * t566 + (mrSges(3,2) * t516 - mrSges(3,3) * t497 + Ifges(3,1) * t533 + Ifges(3,4) * t534 + Ifges(3,5) * t539 - qJ(3) * t400 + t394 * t546 - t395 * t544 + t513 * t572 - t514 * t540) * t583 + (-mrSges(3,1) * t516 + mrSges(3,3) * t498 + Ifges(3,4) * t533 + Ifges(3,2) * t534 + Ifges(3,6) * t539 + pkin(2) * t561 + qJ(3) * t570 + t544 * t394 + t546 * t395 - t513 * t573 + t540 * t515) * t582 + t547 * t393 + pkin(1) * ((t398 * t553 + t399 * t550) * t547 + (-m(3) * t516 + t534 * mrSges(3,1) - t533 * mrSges(3,2) + (-t528 * t550 + t529 * t553) * t577 + t561) * t545) + (-t398 * t550 + t399 * t553) * t591; t393; -t561; t594; t410; t560;];
tauJ  = t1;
