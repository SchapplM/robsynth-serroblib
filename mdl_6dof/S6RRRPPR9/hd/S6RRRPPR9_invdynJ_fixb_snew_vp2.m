% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:33:00
% EndTime: 2019-05-07 06:33:11
% DurationCPUTime: 7.44s
% Computational Cost: add. (82535->340), mult. (180544->429), div. (0->0), fcn. (141669->12), ass. (0->146)
t603 = -2 * qJD(4);
t602 = Ifges(5,1) + Ifges(6,1);
t593 = Ifges(5,4) - Ifges(6,5);
t592 = -Ifges(5,5) - Ifges(6,4);
t601 = Ifges(5,2) + Ifges(6,3);
t600 = -Ifges(6,2) - Ifges(5,3);
t591 = Ifges(5,6) - Ifges(6,6);
t550 = sin(pkin(6));
t554 = sin(qJ(2));
t557 = cos(qJ(2));
t577 = qJD(1) * qJD(2);
t538 = (-qJDD(1) * t557 + t554 * t577) * t550;
t580 = qJD(1) * t550;
t536 = (-pkin(2) * t557 - pkin(9) * t554) * t580;
t551 = cos(pkin(6));
t546 = qJD(1) * t551 + qJD(2);
t544 = t546 ^ 2;
t545 = qJDD(1) * t551 + qJDD(2);
t579 = qJD(1) * t557;
t559 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t558 = cos(qJ(1));
t570 = -g(1) * t558 - g(2) * t555;
t596 = pkin(8) * t550;
t534 = -pkin(1) * t559 + qJDD(1) * t596 + t570;
t574 = t555 * g(1) - g(2) * t558;
t533 = qJDD(1) * pkin(1) + t559 * t596 + t574;
t588 = t533 * t551;
t581 = t557 * t534 + t554 * t588;
t475 = -pkin(2) * t544 + pkin(9) * t545 + (-g(3) * t554 + t536 * t579) * t550 + t581;
t537 = (qJDD(1) * t554 + t557 * t577) * t550;
t595 = g(3) * t551;
t476 = pkin(2) * t538 - pkin(9) * t537 - t595 + (-t533 + (pkin(2) * t554 - pkin(9) * t557) * t546 * qJD(1)) * t550;
t553 = sin(qJ(3));
t597 = cos(qJ(3));
t452 = t597 * t475 + t553 * t476;
t576 = t554 * t580;
t526 = -t597 * t546 + t553 * t576;
t527 = t553 * t546 + t597 * t576;
t507 = pkin(3) * t526 - qJ(4) * t527;
t530 = qJDD(3) + t538;
t575 = t550 * t579;
t543 = qJD(3) - t575;
t542 = t543 ^ 2;
t442 = -pkin(3) * t542 + qJ(4) * t530 - t507 * t526 + t452;
t586 = t550 * t557;
t505 = -g(3) * t586 - t554 * t534 + t557 * t588;
t474 = -pkin(2) * t545 - pkin(9) * t544 + t536 * t576 - t505;
t503 = qJD(3) * t527 + t537 * t553 - t597 * t545;
t504 = -t526 * qJD(3) + t597 * t537 + t553 * t545;
t445 = (t526 * t543 - t504) * qJ(4) + (t527 * t543 + t503) * pkin(3) + t474;
t549 = sin(pkin(11));
t590 = cos(pkin(11));
t513 = t590 * t527 + t549 * t543;
t437 = -t549 * t442 + t590 * t445 + t513 * t603;
t484 = t590 * t504 + t549 * t530;
t451 = -t553 * t475 + t597 * t476;
t564 = pkin(3) * t530 + qJ(4) * t542 - t527 * t507 - qJDD(4) + t451;
t512 = t527 * t549 - t590 * t543;
t589 = t512 * t526;
t599 = (-t484 + t589) * qJ(5) - t564;
t598 = 2 * qJD(5);
t594 = -mrSges(5,3) - mrSges(6,2);
t587 = t550 * t554;
t438 = t590 * t442 + t549 * t445 + t512 * t603;
t483 = t504 * t549 - t590 * t530;
t490 = mrSges(5,1) * t526 - mrSges(5,3) * t513;
t485 = pkin(4) * t512 - qJ(5) * t513;
t525 = t526 ^ 2;
t434 = -pkin(4) * t525 + t503 * qJ(5) - t512 * t485 + t526 * t598 + t438;
t491 = -mrSges(6,1) * t526 + mrSges(6,2) * t513;
t435 = -t503 * pkin(4) - t525 * qJ(5) + t513 * t485 + qJDD(5) - t437;
t429 = (-t484 - t589) * pkin(10) + (t512 * t513 - t503) * pkin(5) + t435;
t492 = -pkin(5) * t526 - pkin(10) * t513;
t511 = t512 ^ 2;
t430 = -pkin(5) * t511 + pkin(10) * t483 + t492 * t526 + t434;
t552 = sin(qJ(6));
t556 = cos(qJ(6));
t427 = t429 * t556 - t430 * t552;
t479 = t512 * t556 - t513 * t552;
t450 = qJD(6) * t479 + t483 * t552 + t484 * t556;
t480 = t512 * t552 + t513 * t556;
t457 = -mrSges(7,1) * t479 + mrSges(7,2) * t480;
t523 = qJD(6) - t526;
t460 = -mrSges(7,2) * t523 + mrSges(7,3) * t479;
t501 = qJDD(6) - t503;
t423 = m(7) * t427 + mrSges(7,1) * t501 - mrSges(7,3) * t450 - t457 * t480 + t460 * t523;
t428 = t429 * t552 + t430 * t556;
t449 = -qJD(6) * t480 + t483 * t556 - t484 * t552;
t461 = mrSges(7,1) * t523 - mrSges(7,3) * t480;
t424 = m(7) * t428 - mrSges(7,2) * t501 + mrSges(7,3) * t449 + t457 * t479 - t461 * t523;
t572 = -t423 * t552 + t556 * t424;
t567 = m(6) * t434 + t503 * mrSges(6,3) + t526 * t491 + t572;
t486 = mrSges(6,1) * t512 - mrSges(6,3) * t513;
t582 = -mrSges(5,1) * t512 - mrSges(5,2) * t513 - t486;
t413 = m(5) * t438 - mrSges(5,2) * t503 + t594 * t483 - t490 * t526 + t582 * t512 + t567;
t416 = t556 * t423 + t552 * t424;
t489 = -mrSges(6,2) * t512 + mrSges(6,3) * t526;
t563 = -m(6) * t435 + t503 * mrSges(6,1) + t526 * t489 - t416;
t569 = -mrSges(5,2) * t526 - mrSges(5,3) * t512;
t414 = m(5) * t437 + t503 * mrSges(5,1) + t594 * t484 + t582 * t513 + t526 * t569 + t563;
t411 = t590 * t413 - t414 * t549;
t508 = mrSges(4,1) * t526 + mrSges(4,2) * t527;
t515 = mrSges(4,1) * t543 - mrSges(4,3) * t527;
t409 = m(4) * t452 - mrSges(4,2) * t530 - mrSges(4,3) * t503 - t508 * t526 - t515 * t543 + t411;
t436 = -0.2e1 * qJD(5) * t513 + (t513 * t526 + t483) * pkin(4) + t599;
t432 = -pkin(10) * t511 + (-pkin(4) - pkin(5)) * t483 + (-pkin(4) * t526 + t492 + t598) * t513 - t599;
t565 = -m(7) * t432 + t449 * mrSges(7,1) - t450 * mrSges(7,2) + t479 * t460 - t480 * t461;
t425 = m(6) * t436 + t483 * mrSges(6,1) - t484 * mrSges(6,3) + t512 * t489 - t513 * t491 + t565;
t421 = -m(5) * t564 + t483 * mrSges(5,1) + t484 * mrSges(5,2) + t513 * t490 + t512 * t569 + t425;
t514 = -mrSges(4,2) * t543 - mrSges(4,3) * t526;
t420 = m(4) * t451 + t530 * mrSges(4,1) - t504 * mrSges(4,3) - t527 * t508 + t543 * t514 - t421;
t405 = t553 * t409 + t597 * t420;
t585 = t512 * t601 - t513 * t593 - t526 * t591;
t584 = t512 * t591 + t513 * t592 + t526 * t600;
t583 = -t512 * t593 + t513 * t602 - t526 * t592;
t573 = t597 * t409 - t420 * t553;
t410 = t549 * t413 + t590 * t414;
t454 = Ifges(7,4) * t480 + Ifges(7,2) * t479 + Ifges(7,6) * t523;
t455 = Ifges(7,1) * t480 + Ifges(7,4) * t479 + Ifges(7,5) * t523;
t562 = mrSges(7,1) * t427 - mrSges(7,2) * t428 + Ifges(7,5) * t450 + Ifges(7,6) * t449 + Ifges(7,3) * t501 + t480 * t454 - t479 * t455;
t561 = -m(4) * t474 - t503 * mrSges(4,1) - t504 * mrSges(4,2) - t526 * t514 - t527 * t515 - t410;
t453 = Ifges(7,5) * t480 + Ifges(7,6) * t479 + Ifges(7,3) * t523;
t417 = -mrSges(7,1) * t432 + mrSges(7,3) * t428 + Ifges(7,4) * t450 + Ifges(7,2) * t449 + Ifges(7,6) * t501 - t453 * t480 + t455 * t523;
t418 = mrSges(7,2) * t432 - mrSges(7,3) * t427 + Ifges(7,1) * t450 + Ifges(7,4) * t449 + Ifges(7,5) * t501 + t453 * t479 - t454 * t523;
t402 = mrSges(5,1) * t564 - mrSges(6,1) * t436 + mrSges(6,2) * t434 + mrSges(5,3) * t438 - pkin(4) * t425 - pkin(5) * t565 - pkin(10) * t572 - t556 * t417 - t552 * t418 - t483 * t601 + t593 * t484 + t591 * t503 + t584 * t513 + t583 * t526;
t404 = -mrSges(5,2) * t564 + mrSges(6,2) * t435 - mrSges(5,3) * t437 - mrSges(6,3) * t436 - pkin(10) * t416 - qJ(5) * t425 - t417 * t552 + t418 * t556 - t593 * t483 + t484 * t602 - t592 * t503 + t584 * t512 + t585 * t526;
t495 = Ifges(4,4) * t527 - Ifges(4,2) * t526 + Ifges(4,6) * t543;
t496 = Ifges(4,1) * t527 - Ifges(4,4) * t526 + Ifges(4,5) * t543;
t560 = mrSges(4,1) * t451 - mrSges(4,2) * t452 + Ifges(4,5) * t504 - Ifges(4,6) * t503 + Ifges(4,3) * t530 - pkin(3) * t421 + qJ(4) * t411 + t590 * t402 + t549 * t404 + t527 * t495 + t526 * t496;
t535 = (-mrSges(3,1) * t557 + mrSges(3,2) * t554) * t580;
t532 = -mrSges(3,2) * t546 + mrSges(3,3) * t575;
t531 = mrSges(3,1) * t546 - mrSges(3,3) * t576;
t519 = -t533 * t550 - t595;
t518 = Ifges(3,5) * t546 + (Ifges(3,1) * t554 + Ifges(3,4) * t557) * t580;
t517 = Ifges(3,6) * t546 + (Ifges(3,4) * t554 + Ifges(3,2) * t557) * t580;
t516 = Ifges(3,3) * t546 + (Ifges(3,5) * t554 + Ifges(3,6) * t557) * t580;
t506 = -g(3) * t587 + t581;
t494 = Ifges(4,5) * t527 - Ifges(4,6) * t526 + Ifges(4,3) * t543;
t415 = t484 * mrSges(6,2) + t513 * t486 - t563;
t406 = m(3) * t505 + t545 * mrSges(3,1) - t537 * mrSges(3,3) + t546 * t532 - t535 * t576 + t561;
t403 = m(3) * t506 - mrSges(3,2) * t545 - mrSges(3,3) * t538 - t531 * t546 + t535 * t575 + t573;
t401 = t562 - pkin(3) * t410 + t543 * t496 + pkin(5) * t416 + pkin(4) * t415 + Ifges(4,6) * t530 - t527 * t494 + Ifges(4,4) * t504 - mrSges(4,1) * t474 + mrSges(4,3) * t452 - mrSges(5,1) * t437 + mrSges(5,2) * t438 - mrSges(6,3) * t434 + mrSges(6,1) * t435 + (mrSges(6,2) * qJ(5) + t591) * t483 + t592 * t484 + (-Ifges(4,2) + t600) * t503 + (qJ(5) * t486 - t583) * t512 + t585 * t513 - qJ(5) * t567;
t400 = mrSges(4,2) * t474 - mrSges(4,3) * t451 + Ifges(4,1) * t504 - Ifges(4,4) * t503 + Ifges(4,5) * t530 - qJ(4) * t410 - t549 * t402 + t590 * t404 - t526 * t494 - t543 * t495;
t399 = Ifges(3,5) * t537 - Ifges(3,6) * t538 + Ifges(3,3) * t545 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + t553 * t400 + t597 * t401 + pkin(2) * t561 + pkin(9) * t573 + (t517 * t554 - t518 * t557) * t580;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t574 - mrSges(2,2) * t570 + (mrSges(3,2) * t519 - mrSges(3,3) * t505 + Ifges(3,1) * t537 - Ifges(3,4) * t538 + Ifges(3,5) * t545 - pkin(9) * t405 + t597 * t400 - t553 * t401 + t516 * t575 - t546 * t517) * t587 + (-mrSges(3,1) * t519 + mrSges(3,3) * t506 + Ifges(3,4) * t537 - Ifges(3,2) * t538 + Ifges(3,6) * t545 - pkin(2) * t405 - t516 * t576 + t546 * t518 - t560) * t586 + t551 * t399 + pkin(1) * ((t403 * t554 + t406 * t557) * t551 + (-m(3) * t519 - t538 * mrSges(3,1) - t537 * mrSges(3,2) + (-t531 * t554 + t532 * t557) * t580 - t405) * t550) + (t403 * t557 - t406 * t554) * t596; t399; t560; t421; t415; t562;];
tauJ  = t1;
