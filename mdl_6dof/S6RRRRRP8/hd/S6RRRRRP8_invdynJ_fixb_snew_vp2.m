% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:39:52
% EndTime: 2019-05-08 05:40:10
% DurationCPUTime: 8.66s
% Computational Cost: add. (124779->339), mult. (265090->428), div. (0->0), fcn. (213917->12), ass. (0->145)
t586 = Ifges(6,1) + Ifges(7,1);
t578 = Ifges(6,4) - Ifges(7,5);
t577 = -Ifges(6,5) - Ifges(7,4);
t585 = Ifges(6,2) + Ifges(7,3);
t576 = Ifges(6,6) - Ifges(7,6);
t584 = -Ifges(6,3) - Ifges(7,2);
t538 = cos(pkin(6));
t534 = qJD(1) * t538 + qJD(2);
t541 = sin(qJ(3));
t545 = cos(qJ(3));
t542 = sin(qJ(2));
t537 = sin(pkin(6));
t566 = qJD(1) * t537;
t562 = t542 * t566;
t514 = t534 * t541 + t545 * t562;
t546 = cos(qJ(2));
t565 = qJD(1) * t546;
t525 = (qJD(2) * t565 + qJDD(1) * t542) * t537;
t533 = qJDD(1) * t538 + qJDD(2);
t492 = -qJD(3) * t514 - t525 * t541 + t533 * t545;
t513 = t534 * t545 - t541 * t562;
t493 = qJD(3) * t513 + t525 * t545 + t533 * t541;
t540 = sin(qJ(4));
t544 = cos(qJ(4));
t498 = t513 * t544 - t514 * t540;
t460 = qJD(4) * t498 + t492 * t540 + t493 * t544;
t499 = t513 * t540 + t514 * t544;
t561 = t537 * t565;
t530 = qJD(3) - t561;
t529 = qJD(4) + t530;
t539 = sin(qJ(5));
t582 = cos(qJ(5));
t482 = t539 * t499 - t582 * t529;
t564 = qJDD(1) * t537;
t526 = -qJD(2) * t562 + t546 * t564;
t518 = qJDD(3) - t526;
t517 = qJDD(4) + t518;
t439 = -t482 * qJD(5) + t582 * t460 + t539 * t517;
t483 = t582 * t499 + t539 * t529;
t464 = mrSges(7,1) * t482 - mrSges(7,3) * t483;
t524 = (-pkin(2) * t546 - pkin(9) * t542) * t566;
t532 = t534 ^ 2;
t548 = qJD(1) ^ 2;
t543 = sin(qJ(1));
t547 = cos(qJ(1));
t555 = -g(1) * t547 - g(2) * t543;
t522 = -pkin(1) * t548 + pkin(8) * t564 + t555;
t560 = t543 * g(1) - g(2) * t547;
t581 = pkin(8) * t537;
t521 = qJDD(1) * pkin(1) + t548 * t581 + t560;
t574 = t521 * t538;
t567 = t546 * t522 + t542 * t574;
t480 = -t532 * pkin(2) + t533 * pkin(9) + (-g(3) * t542 + t524 * t565) * t537 + t567;
t580 = t538 * g(3);
t481 = -t526 * pkin(2) - t525 * pkin(9) - t580 + (-t521 + (pkin(2) * t542 - pkin(9) * t546) * t534 * qJD(1)) * t537;
t441 = -t541 * t480 + t545 * t481;
t432 = (t513 * t530 - t493) * pkin(10) + (t513 * t514 + t518) * pkin(3) + t441;
t442 = t545 * t480 + t541 * t481;
t503 = pkin(3) * t530 - pkin(10) * t514;
t512 = t513 ^ 2;
t435 = -pkin(3) * t512 + pkin(10) * t492 - t503 * t530 + t442;
t430 = t540 * t432 + t544 * t435;
t475 = -pkin(4) * t498 - pkin(11) * t499;
t528 = t529 ^ 2;
t425 = -pkin(4) * t528 + pkin(11) * t517 + t475 * t498 + t430;
t572 = t537 * t546;
t494 = -g(3) * t572 - t542 * t522 + t546 * t574;
t479 = -t533 * pkin(2) - t532 * pkin(9) + t524 * t562 - t494;
t440 = -t492 * pkin(3) - t512 * pkin(10) + t514 * t503 + t479;
t459 = -qJD(4) * t499 + t492 * t544 - t493 * t540;
t427 = (-t498 * t529 - t460) * pkin(11) + (t499 * t529 - t459) * pkin(4) + t440;
t421 = -t539 * t425 + t582 * t427;
t458 = qJDD(5) - t459;
t463 = pkin(5) * t482 - qJ(6) * t483;
t497 = qJD(5) - t498;
t496 = t497 ^ 2;
t419 = -t458 * pkin(5) - t496 * qJ(6) + t483 * t463 + qJDD(6) - t421;
t466 = -mrSges(7,2) * t482 + mrSges(7,3) * t497;
t556 = -m(7) * t419 + t458 * mrSges(7,1) + t497 * t466;
t415 = t439 * mrSges(7,2) + t483 * t464 - t556;
t422 = t582 * t425 + t539 * t427;
t418 = -pkin(5) * t496 + qJ(6) * t458 + 0.2e1 * qJD(6) * t497 - t463 * t482 + t422;
t438 = t483 * qJD(5) + t539 * t460 - t582 * t517;
t469 = -mrSges(7,1) * t497 + mrSges(7,2) * t483;
t563 = m(7) * t418 + t458 * mrSges(7,3) + t497 * t469;
t569 = t578 * t482 - t483 * t586 + t577 * t497;
t570 = t482 * t585 - t483 * t578 - t497 * t576;
t583 = -t576 * t438 - t577 * t439 - t584 * t458 - t569 * t482 - t570 * t483 + mrSges(6,1) * t421 - mrSges(7,1) * t419 - mrSges(6,2) * t422 + mrSges(7,3) * t418 - pkin(5) * t415 + qJ(6) * (-t438 * mrSges(7,2) - t482 * t464 + t563);
t579 = -mrSges(6,3) - mrSges(7,2);
t573 = t537 * t542;
t474 = -mrSges(5,1) * t498 + mrSges(5,2) * t499;
t485 = mrSges(5,1) * t529 - mrSges(5,3) * t499;
t468 = mrSges(6,1) * t497 - mrSges(6,3) * t483;
t568 = -mrSges(6,1) * t482 - mrSges(6,2) * t483 - t464;
t409 = m(6) * t422 - t458 * mrSges(6,2) + t579 * t438 - t497 * t468 + t568 * t482 + t563;
t467 = -mrSges(6,2) * t497 - mrSges(6,3) * t482;
t411 = m(6) * t421 + t458 * mrSges(6,1) + t579 * t439 + t497 * t467 + t568 * t483 + t556;
t557 = t582 * t409 - t411 * t539;
t397 = m(5) * t430 - mrSges(5,2) * t517 + mrSges(5,3) * t459 + t474 * t498 - t485 * t529 + t557;
t429 = t544 * t432 - t540 * t435;
t484 = -mrSges(5,2) * t529 + mrSges(5,3) * t498;
t424 = -t517 * pkin(4) - t528 * pkin(11) + t499 * t475 - t429;
t420 = -0.2e1 * qJD(6) * t483 + (t482 * t497 - t439) * qJ(6) + (t483 * t497 + t438) * pkin(5) + t424;
t416 = m(7) * t420 + mrSges(7,1) * t438 - t439 * mrSges(7,3) + t466 * t482 - t483 * t469;
t551 = -m(6) * t424 - t438 * mrSges(6,1) - mrSges(6,2) * t439 - t482 * t467 - t468 * t483 - t416;
t406 = m(5) * t429 + mrSges(5,1) * t517 - mrSges(5,3) * t460 - t474 * t499 + t484 * t529 + t551;
t393 = t540 * t397 + t544 * t406;
t500 = -mrSges(4,1) * t513 + mrSges(4,2) * t514;
t501 = -mrSges(4,2) * t530 + mrSges(4,3) * t513;
t391 = m(4) * t441 + mrSges(4,1) * t518 - mrSges(4,3) * t493 - t500 * t514 + t501 * t530 + t393;
t502 = mrSges(4,1) * t530 - mrSges(4,3) * t514;
t558 = t544 * t397 - t406 * t540;
t392 = m(4) * t442 - mrSges(4,2) * t518 + mrSges(4,3) * t492 + t500 * t513 - t502 * t530 + t558;
t385 = t545 * t391 + t541 * t392;
t403 = t539 * t409 + t582 * t411;
t571 = t482 * t576 + t483 * t577 + t497 * t584;
t559 = -t391 * t541 + t545 * t392;
t554 = m(5) * t440 - t459 * mrSges(5,1) + mrSges(5,2) * t460 - t498 * t484 + t485 * t499 + t403;
t399 = -mrSges(6,1) * t424 - mrSges(7,1) * t420 + mrSges(7,2) * t418 + mrSges(6,3) * t422 - pkin(5) * t416 - t438 * t585 + t578 * t439 + t576 * t458 + t571 * t483 - t569 * t497;
t401 = mrSges(6,2) * t424 + mrSges(7,2) * t419 - mrSges(6,3) * t421 - mrSges(7,3) * t420 - qJ(6) * t416 - t578 * t438 + t439 * t586 - t577 * t458 + t571 * t482 + t570 * t497;
t471 = Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t529;
t472 = Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t529;
t553 = -mrSges(5,1) * t429 + mrSges(5,2) * t430 - Ifges(5,5) * t460 - Ifges(5,6) * t459 - Ifges(5,3) * t517 - pkin(4) * t551 - pkin(11) * t557 - t582 * t399 - t539 * t401 - t499 * t471 + t498 * t472;
t550 = -m(4) * t479 + t492 * mrSges(4,1) - mrSges(4,2) * t493 + t513 * t501 - t502 * t514 - t554;
t487 = Ifges(4,4) * t514 + Ifges(4,2) * t513 + Ifges(4,6) * t530;
t488 = Ifges(4,1) * t514 + Ifges(4,4) * t513 + Ifges(4,5) * t530;
t549 = mrSges(4,1) * t441 - mrSges(4,2) * t442 + Ifges(4,5) * t493 + Ifges(4,6) * t492 + Ifges(4,3) * t518 + pkin(3) * t393 + t514 * t487 - t513 * t488 - t553;
t523 = (-mrSges(3,1) * t546 + mrSges(3,2) * t542) * t566;
t520 = -mrSges(3,2) * t534 + mrSges(3,3) * t561;
t519 = mrSges(3,1) * t534 - mrSges(3,3) * t562;
t507 = -t537 * t521 - t580;
t506 = Ifges(3,5) * t534 + (Ifges(3,1) * t542 + Ifges(3,4) * t546) * t566;
t505 = Ifges(3,6) * t534 + (Ifges(3,4) * t542 + Ifges(3,2) * t546) * t566;
t504 = Ifges(3,3) * t534 + (Ifges(3,5) * t542 + Ifges(3,6) * t546) * t566;
t495 = -g(3) * t573 + t567;
t486 = Ifges(4,5) * t514 + Ifges(4,6) * t513 + Ifges(4,3) * t530;
t470 = Ifges(5,5) * t499 + Ifges(5,6) * t498 + Ifges(5,3) * t529;
t394 = m(3) * t494 + mrSges(3,1) * t533 - mrSges(3,3) * t525 + t520 * t534 - t523 * t562 + t550;
t387 = -mrSges(5,1) * t440 + mrSges(5,3) * t430 + Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * t517 - pkin(4) * t403 - t499 * t470 + t529 * t472 - t583;
t386 = mrSges(5,2) * t440 - mrSges(5,3) * t429 + Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * t517 - pkin(11) * t403 - t539 * t399 + t582 * t401 + t498 * t470 - t529 * t471;
t384 = m(3) * t495 - mrSges(3,2) * t533 + mrSges(3,3) * t526 - t519 * t534 + t523 * t561 + t559;
t383 = mrSges(4,2) * t479 - mrSges(4,3) * t441 + Ifges(4,1) * t493 + Ifges(4,4) * t492 + Ifges(4,5) * t518 - pkin(10) * t393 + t386 * t544 - t387 * t540 + t486 * t513 - t487 * t530;
t382 = -mrSges(4,1) * t479 + mrSges(4,3) * t442 + Ifges(4,4) * t493 + Ifges(4,2) * t492 + Ifges(4,6) * t518 - pkin(3) * t554 + pkin(10) * t558 + t540 * t386 + t544 * t387 - t514 * t486 + t530 * t488;
t381 = Ifges(3,5) * t525 + Ifges(3,6) * t526 + Ifges(3,3) * t533 + mrSges(3,1) * t494 - mrSges(3,2) * t495 + t541 * t383 + t545 * t382 + pkin(2) * t550 + pkin(9) * t559 + (t505 * t542 - t506 * t546) * t566;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t560 - mrSges(2,2) * t555 + (mrSges(3,2) * t507 - mrSges(3,3) * t494 + Ifges(3,1) * t525 + Ifges(3,4) * t526 + Ifges(3,5) * t533 - pkin(9) * t385 - t382 * t541 + t383 * t545 + t504 * t561 - t505 * t534) * t573 + (-mrSges(3,1) * t507 + mrSges(3,3) * t495 + Ifges(3,4) * t525 + Ifges(3,2) * t526 + Ifges(3,6) * t533 - pkin(2) * t385 - t504 * t562 + t534 * t506 - t549) * t572 + t538 * t381 + pkin(1) * ((t384 * t542 + t394 * t546) * t538 + (-m(3) * t507 + t526 * mrSges(3,1) - t525 * mrSges(3,2) + (-t519 * t542 + t520 * t546) * t566 - t385) * t537) + (t384 * t546 - t394 * t542) * t581; t381; t549; -t553; t583; t415;];
tauJ  = t1;
