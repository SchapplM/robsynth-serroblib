% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 01:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:20:43
% EndTime: 2019-05-08 01:20:59
% DurationCPUTime: 6.65s
% Computational Cost: add. (91844->342), mult. (194236->429), div. (0->0), fcn. (153522->12), ass. (0->148)
t602 = Ifges(5,1) + Ifges(6,1);
t592 = Ifges(5,4) - Ifges(6,5);
t591 = Ifges(5,5) + Ifges(6,4);
t601 = -Ifges(5,2) - Ifges(6,3);
t590 = Ifges(5,6) - Ifges(6,6);
t600 = Ifges(5,3) + Ifges(6,2);
t549 = sin(pkin(6));
t554 = sin(qJ(2));
t558 = cos(qJ(2));
t577 = qJD(1) * qJD(2);
t538 = (-qJDD(1) * t558 + t554 * t577) * t549;
t550 = cos(pkin(6));
t546 = qJD(1) * t550 + qJD(2);
t553 = sin(qJ(3));
t557 = cos(qJ(3));
t579 = qJD(1) * t549;
t576 = t554 * t579;
t526 = t557 * t546 - t553 * t576;
t537 = (qJDD(1) * t554 + t558 * t577) * t549;
t545 = qJDD(1) * t550 + qJDD(2);
t505 = qJD(3) * t526 + t537 * t557 + t545 * t553;
t527 = t546 * t553 + t557 * t576;
t578 = qJD(1) * t558;
t575 = t549 * t578;
t543 = qJD(3) - t575;
t552 = sin(qJ(4));
t596 = cos(qJ(4));
t511 = t527 * t552 - t543 * t596;
t530 = qJDD(3) + t538;
t460 = -t511 * qJD(4) + t505 * t596 + t552 * t530;
t536 = (-pkin(2) * t558 - pkin(9) * t554) * t579;
t544 = t546 ^ 2;
t560 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t559 = cos(qJ(1));
t570 = -g(1) * t559 - g(2) * t555;
t595 = pkin(8) * t549;
t534 = -pkin(1) * t560 + qJDD(1) * t595 + t570;
t574 = t555 * g(1) - g(2) * t559;
t533 = qJDD(1) * pkin(1) + t560 * t595 + t574;
t587 = t533 * t550;
t580 = t558 * t534 + t554 * t587;
t478 = -t544 * pkin(2) + t545 * pkin(9) + (-g(3) * t554 + t536 * t578) * t549 + t580;
t594 = t550 * g(3);
t479 = t538 * pkin(2) - t537 * pkin(9) - t594 + (-t533 + (pkin(2) * t554 - pkin(9) * t558) * t546 * qJD(1)) * t549;
t449 = -t553 * t478 + t557 * t479;
t509 = -pkin(3) * t526 - pkin(10) * t527;
t542 = t543 ^ 2;
t566 = t530 * pkin(3) + t542 * pkin(10) - t527 * t509 + t449;
t524 = qJD(4) - t526;
t588 = t511 * t524;
t599 = (-t460 + t588) * qJ(5) - t566;
t512 = t527 * t596 + t552 * t543;
t485 = mrSges(6,1) * t511 - mrSges(6,3) * t512;
t450 = t557 * t478 + t553 * t479;
t445 = -pkin(3) * t542 + pkin(10) * t530 + t509 * t526 + t450;
t585 = t549 * t558;
t506 = -g(3) * t585 - t554 * t534 + t558 * t587;
t477 = -t545 * pkin(2) - t544 * pkin(9) + t536 * t576 - t506;
t504 = -qJD(3) * t527 - t537 * t553 + t557 * t545;
t448 = (-t526 * t543 - t505) * pkin(10) + (t527 * t543 - t504) * pkin(3) + t477;
t434 = -t552 * t445 + t448 * t596;
t484 = pkin(4) * t511 - qJ(5) * t512;
t502 = qJDD(4) - t504;
t523 = t524 ^ 2;
t432 = -t502 * pkin(4) - t523 * qJ(5) + t512 * t484 + qJDD(5) - t434;
t426 = (-t460 - t588) * pkin(11) + (t511 * t512 - t502) * pkin(5) + t432;
t435 = t596 * t445 + t552 * t448;
t597 = 2 * qJD(5);
t431 = -pkin(4) * t523 + t502 * qJ(5) - t511 * t484 + t524 * t597 + t435;
t459 = qJD(4) * t512 + t505 * t552 - t530 * t596;
t492 = -pkin(5) * t524 - pkin(11) * t512;
t510 = t511 ^ 2;
t427 = -pkin(5) * t510 + pkin(11) * t459 + t492 * t524 + t431;
t551 = sin(qJ(6));
t556 = cos(qJ(6));
t424 = t426 * t556 - t427 * t551;
t482 = t511 * t556 - t512 * t551;
t441 = qJD(6) * t482 + t459 * t551 + t460 * t556;
t483 = t511 * t551 + t512 * t556;
t456 = -mrSges(7,1) * t482 + mrSges(7,2) * t483;
t522 = qJD(6) - t524;
t463 = -mrSges(7,2) * t522 + mrSges(7,3) * t482;
t497 = qJDD(6) - t502;
t421 = m(7) * t424 + mrSges(7,1) * t497 - mrSges(7,3) * t441 - t456 * t483 + t463 * t522;
t425 = t426 * t551 + t427 * t556;
t440 = -qJD(6) * t483 + t459 * t556 - t460 * t551;
t464 = mrSges(7,1) * t522 - mrSges(7,3) * t483;
t422 = m(7) * t425 - mrSges(7,2) * t497 + mrSges(7,3) * t440 + t456 * t482 - t464 * t522;
t413 = t556 * t421 + t551 * t422;
t488 = -mrSges(6,2) * t511 + mrSges(6,3) * t524;
t565 = -m(6) * t432 + t502 * mrSges(6,1) + t524 * t488 - t413;
t412 = t460 * mrSges(6,2) + t512 * t485 - t565;
t452 = Ifges(7,4) * t483 + Ifges(7,2) * t482 + Ifges(7,6) * t522;
t453 = Ifges(7,1) * t483 + Ifges(7,4) * t482 + Ifges(7,5) * t522;
t564 = -mrSges(7,1) * t424 + mrSges(7,2) * t425 - Ifges(7,5) * t441 - Ifges(7,6) * t440 - Ifges(7,3) * t497 - t483 * t452 + t482 * t453;
t491 = -mrSges(6,1) * t524 + mrSges(6,2) * t512;
t572 = -t551 * t421 + t556 * t422;
t568 = m(6) * t431 + t502 * mrSges(6,3) + t524 * t491 + t572;
t582 = -t592 * t511 + t602 * t512 + t591 * t524;
t583 = t601 * t511 + t592 * t512 + t590 * t524;
t598 = -t459 * t590 + t460 * t591 + t600 * t502 + t511 * t582 + t512 * t583 + mrSges(5,1) * t434 - mrSges(6,1) * t432 - mrSges(5,2) * t435 + mrSges(6,3) * t431 - pkin(4) * t412 - pkin(5) * t413 + qJ(5) * (-t459 * mrSges(6,2) - t511 * t485 + t568) + t564;
t593 = -mrSges(5,3) - mrSges(6,2);
t586 = t549 * t554;
t490 = mrSges(5,1) * t524 - mrSges(5,3) * t512;
t581 = -mrSges(5,1) * t511 - mrSges(5,2) * t512 - t485;
t409 = m(5) * t435 - t502 * mrSges(5,2) + t459 * t593 - t524 * t490 + t511 * t581 + t568;
t489 = -mrSges(5,2) * t524 - mrSges(5,3) * t511;
t410 = m(5) * t434 + t502 * mrSges(5,1) + t460 * t593 + t524 * t489 + t512 * t581 + t565;
t407 = t596 * t409 - t410 * t552;
t508 = -mrSges(4,1) * t526 + mrSges(4,2) * t527;
t514 = mrSges(4,1) * t543 - mrSges(4,3) * t527;
t405 = m(4) * t450 - mrSges(4,2) * t530 + mrSges(4,3) * t504 + t508 * t526 - t514 * t543 + t407;
t433 = -0.2e1 * qJD(5) * t512 + (t512 * t524 + t459) * pkin(4) + t599;
t429 = -t510 * pkin(11) + (-pkin(4) - pkin(5)) * t459 + (-pkin(4) * t524 + t492 + t597) * t512 - t599;
t569 = -m(7) * t429 + t440 * mrSges(7,1) - t441 * mrSges(7,2) + t482 * t463 - t483 * t464;
t419 = m(6) * t433 + mrSges(6,1) * t459 - t460 * mrSges(6,3) + t488 * t511 - t512 * t491 + t569;
t418 = m(5) * t566 - t459 * mrSges(5,1) - mrSges(5,2) * t460 - t511 * t489 - t490 * t512 - t419;
t513 = -mrSges(4,2) * t543 + mrSges(4,3) * t526;
t417 = m(4) * t449 + mrSges(4,1) * t530 - mrSges(4,3) * t505 - t508 * t527 + t513 * t543 + t418;
t401 = t553 * t405 + t557 * t417;
t584 = t590 * t511 - t591 * t512 - t600 * t524;
t573 = t557 * t405 - t417 * t553;
t406 = t552 * t409 + t410 * t596;
t563 = -m(4) * t477 + t504 * mrSges(4,1) - t505 * mrSges(4,2) + t526 * t513 - t527 * t514 - t406;
t451 = Ifges(7,5) * t483 + Ifges(7,6) * t482 + Ifges(7,3) * t522;
t414 = -mrSges(7,1) * t429 + mrSges(7,3) * t425 + Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t497 - t451 * t483 + t453 * t522;
t415 = mrSges(7,2) * t429 - mrSges(7,3) * t424 + Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t497 + t451 * t482 - t452 * t522;
t398 = mrSges(5,1) * t566 - mrSges(6,1) * t433 + mrSges(6,2) * t431 + mrSges(5,3) * t435 - pkin(4) * t419 - pkin(5) * t569 - pkin(11) * t572 - t556 * t414 - t551 * t415 + t601 * t459 + t592 * t460 + t590 * t502 + t584 * t512 + t582 * t524;
t399 = -mrSges(5,2) * t566 + mrSges(6,2) * t432 - mrSges(5,3) * t434 - mrSges(6,3) * t433 - pkin(11) * t413 - qJ(5) * t419 - t551 * t414 + t556 * t415 - t592 * t459 + t602 * t460 + t591 * t502 + t584 * t511 - t583 * t524;
t499 = Ifges(4,4) * t527 + Ifges(4,2) * t526 + Ifges(4,6) * t543;
t500 = Ifges(4,1) * t527 + Ifges(4,4) * t526 + Ifges(4,5) * t543;
t561 = mrSges(4,1) * t449 - mrSges(4,2) * t450 + Ifges(4,5) * t505 + Ifges(4,6) * t504 + Ifges(4,3) * t530 + pkin(3) * t418 + pkin(10) * t407 + t398 * t596 + t552 * t399 + t527 * t499 - t526 * t500;
t535 = (-mrSges(3,1) * t558 + mrSges(3,2) * t554) * t579;
t532 = -mrSges(3,2) * t546 + mrSges(3,3) * t575;
t531 = mrSges(3,1) * t546 - mrSges(3,3) * t576;
t518 = -t549 * t533 - t594;
t517 = Ifges(3,5) * t546 + (Ifges(3,1) * t554 + Ifges(3,4) * t558) * t579;
t516 = Ifges(3,6) * t546 + (Ifges(3,4) * t554 + Ifges(3,2) * t558) * t579;
t515 = Ifges(3,3) * t546 + (Ifges(3,5) * t554 + Ifges(3,6) * t558) * t579;
t507 = -g(3) * t586 + t580;
t498 = Ifges(4,5) * t527 + Ifges(4,6) * t526 + Ifges(4,3) * t543;
t402 = m(3) * t506 + t545 * mrSges(3,1) - t537 * mrSges(3,3) + t546 * t532 - t535 * t576 + t563;
t400 = m(3) * t507 - mrSges(3,2) * t545 - mrSges(3,3) * t538 - t531 * t546 + t535 * t575 + t573;
t397 = -mrSges(4,1) * t477 + mrSges(4,3) * t450 + Ifges(4,4) * t505 + Ifges(4,2) * t504 + Ifges(4,6) * t530 - pkin(3) * t406 - t527 * t498 + t543 * t500 - t598;
t396 = mrSges(4,2) * t477 - mrSges(4,3) * t449 + Ifges(4,1) * t505 + Ifges(4,4) * t504 + Ifges(4,5) * t530 - pkin(10) * t406 - t552 * t398 + t399 * t596 + t526 * t498 - t543 * t499;
t395 = Ifges(3,5) * t537 - Ifges(3,6) * t538 + Ifges(3,3) * t545 + mrSges(3,1) * t506 - mrSges(3,2) * t507 + t553 * t396 + t557 * t397 + pkin(2) * t563 + pkin(9) * t573 + (t516 * t554 - t517 * t558) * t579;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t574 - mrSges(2,2) * t570 + (mrSges(3,2) * t518 - mrSges(3,3) * t506 + Ifges(3,1) * t537 - Ifges(3,4) * t538 + Ifges(3,5) * t545 - pkin(9) * t401 + t396 * t557 - t397 * t553 + t515 * t575 - t516 * t546) * t586 + (-mrSges(3,1) * t518 + mrSges(3,3) * t507 + Ifges(3,4) * t537 - Ifges(3,2) * t538 + Ifges(3,6) * t545 - pkin(2) * t401 - t515 * t576 + t546 * t517 - t561) * t585 + t550 * t395 + pkin(1) * ((t400 * t554 + t402 * t558) * t550 + (-m(3) * t518 - t538 * mrSges(3,1) - t537 * mrSges(3,2) + (-t531 * t554 + t532 * t558) * t579 - t401) * t549) + (t400 * t558 - t402 * t554) * t595; t395; t561; t598; t412; -t564;];
tauJ  = t1;
