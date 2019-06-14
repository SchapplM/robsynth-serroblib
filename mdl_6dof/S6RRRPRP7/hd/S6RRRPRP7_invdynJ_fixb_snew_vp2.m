% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:14:44
% EndTime: 2019-05-07 08:14:59
% DurationCPUTime: 7.98s
% Computational Cost: add. (110720->337), mult. (243884->430), div. (0->0), fcn. (194833->12), ass. (0->142)
t583 = Ifges(6,1) + Ifges(7,1);
t574 = Ifges(6,4) - Ifges(7,5);
t573 = -Ifges(6,5) - Ifges(7,4);
t582 = Ifges(6,2) + Ifges(7,3);
t572 = Ifges(6,6) - Ifges(7,6);
t581 = -Ifges(6,3) - Ifges(7,2);
t534 = sin(pkin(6));
t539 = sin(qJ(2));
t542 = cos(qJ(2));
t558 = qJD(1) * qJD(2);
t524 = (-qJDD(1) * t542 + t539 * t558) * t534;
t561 = qJD(1) * t534;
t522 = (-pkin(2) * t542 - pkin(9) * t539) * t561;
t536 = cos(pkin(6));
t530 = qJD(1) * t536 + qJD(2);
t528 = t530 ^ 2;
t529 = qJDD(1) * t536 + qJDD(2);
t560 = qJD(1) * t542;
t544 = qJD(1) ^ 2;
t540 = sin(qJ(1));
t543 = cos(qJ(1));
t549 = -g(1) * t543 - g(2) * t540;
t577 = pkin(8) * t534;
t520 = -pkin(1) * t544 + qJDD(1) * t577 + t549;
t554 = t540 * g(1) - g(2) * t543;
t519 = qJDD(1) * pkin(1) + t544 * t577 + t554;
t569 = t519 * t536;
t562 = t542 * t520 + t539 * t569;
t478 = -t528 * pkin(2) + t529 * pkin(9) + (-g(3) * t539 + t522 * t560) * t534 + t562;
t523 = (qJDD(1) * t539 + t542 * t558) * t534;
t576 = t536 * g(3);
t479 = t524 * pkin(2) - t523 * pkin(9) - t576 + (-t519 + (pkin(2) * t539 - pkin(9) * t542) * t530 * qJD(1)) * t534;
t538 = sin(qJ(3));
t541 = cos(qJ(3));
t442 = -t538 * t478 + t541 * t479;
t556 = t539 * t561;
t512 = t530 * t541 - t538 * t556;
t491 = qJD(3) * t512 + t523 * t541 + t529 * t538;
t513 = t530 * t538 + t541 * t556;
t516 = qJDD(3) + t524;
t555 = t534 * t560;
t527 = qJD(3) - t555;
t433 = (t512 * t527 - t491) * qJ(4) + (t512 * t513 + t516) * pkin(3) + t442;
t443 = t541 * t478 + t538 * t479;
t490 = -qJD(3) * t513 - t523 * t538 + t529 * t541;
t502 = pkin(3) * t527 - qJ(4) * t513;
t511 = t512 ^ 2;
t436 = -pkin(3) * t511 + qJ(4) * t490 - t502 * t527 + t443;
t533 = sin(pkin(11));
t535 = cos(pkin(11));
t499 = t512 * t533 + t513 * t535;
t428 = -0.2e1 * qJD(4) * t499 + t535 * t433 - t533 * t436;
t498 = t512 * t535 - t513 * t533;
t429 = 0.2e1 * qJD(4) * t498 + t533 * t433 + t535 * t436;
t473 = -pkin(4) * t498 - pkin(10) * t499;
t526 = t527 ^ 2;
t427 = -pkin(4) * t526 + pkin(10) * t516 + t473 * t498 + t429;
t567 = t534 * t542;
t492 = -g(3) * t567 - t539 * t520 + t542 * t569;
t477 = -t529 * pkin(2) - t528 * pkin(9) + t522 * t556 - t492;
t437 = -t490 * pkin(3) - t511 * qJ(4) + t513 * t502 + qJDD(4) + t477;
t466 = t490 * t535 - t491 * t533;
t467 = t490 * t533 + t491 * t535;
t431 = (-t498 * t527 - t467) * pkin(10) + (t499 * t527 - t466) * pkin(4) + t437;
t537 = sin(qJ(5));
t578 = cos(qJ(5));
t424 = t427 * t578 + t537 * t431;
t481 = t499 * t578 + t537 * t527;
t440 = qJD(5) * t481 + t467 * t537 - t516 * t578;
t497 = qJD(5) - t498;
t460 = mrSges(6,1) * t497 - mrSges(6,3) * t481;
t465 = qJDD(5) - t466;
t480 = t499 * t537 - t527 * t578;
t455 = pkin(5) * t480 - qJ(6) * t481;
t496 = t497 ^ 2;
t420 = -pkin(5) * t496 + qJ(6) * t465 + 0.2e1 * qJD(6) * t497 - t455 * t480 + t424;
t461 = -mrSges(7,1) * t497 + mrSges(7,2) * t481;
t557 = m(7) * t420 + t465 * mrSges(7,3) + t497 * t461;
t456 = mrSges(7,1) * t480 - mrSges(7,3) * t481;
t563 = -mrSges(6,1) * t480 - mrSges(6,2) * t481 - t456;
t575 = -mrSges(6,3) - mrSges(7,2);
t412 = m(6) * t424 - t465 * mrSges(6,2) + t440 * t575 - t497 * t460 + t480 * t563 + t557;
t423 = -t537 * t427 + t431 * t578;
t441 = -t480 * qJD(5) + t467 * t578 + t537 * t516;
t459 = -mrSges(6,2) * t497 - mrSges(6,3) * t480;
t421 = -t465 * pkin(5) - t496 * qJ(6) + t481 * t455 + qJDD(6) - t423;
t458 = -mrSges(7,2) * t480 + mrSges(7,3) * t497;
t550 = -m(7) * t421 + t465 * mrSges(7,1) + t497 * t458;
t414 = m(6) * t423 + t465 * mrSges(6,1) + t441 * t575 + t497 * t459 + t481 * t563 + t550;
t407 = t412 * t578 - t414 * t537;
t472 = -mrSges(5,1) * t498 + mrSges(5,2) * t499;
t483 = mrSges(5,1) * t527 - mrSges(5,3) * t499;
t402 = m(5) * t429 - mrSges(5,2) * t516 + mrSges(5,3) * t466 + t472 * t498 - t483 * t527 + t407;
t426 = -t516 * pkin(4) - t526 * pkin(10) + t499 * t473 - t428;
t422 = -0.2e1 * qJD(6) * t481 + (t480 * t497 - t441) * qJ(6) + (t481 * t497 + t440) * pkin(5) + t426;
t418 = m(7) * t422 + mrSges(7,1) * t440 - t441 * mrSges(7,3) + t458 * t480 - t481 * t461;
t415 = -m(6) * t426 - t440 * mrSges(6,1) - mrSges(6,2) * t441 - t480 * t459 - t460 * t481 - t418;
t482 = -mrSges(5,2) * t527 + mrSges(5,3) * t498;
t409 = m(5) * t428 + mrSges(5,1) * t516 - mrSges(5,3) * t467 - t472 * t499 + t482 * t527 + t415;
t398 = t533 * t402 + t535 * t409;
t564 = t574 * t480 - t583 * t481 + t573 * t497;
t566 = t572 * t480 + t573 * t481 + t581 * t497;
t403 = -mrSges(6,1) * t426 - mrSges(7,1) * t422 + mrSges(7,2) * t420 + mrSges(6,3) * t424 - pkin(5) * t418 - t440 * t582 + t574 * t441 + t572 * t465 + t566 * t481 - t564 * t497;
t565 = t582 * t480 - t574 * t481 - t572 * t497;
t404 = mrSges(6,2) * t426 + mrSges(7,2) * t421 - mrSges(6,3) * t423 - mrSges(7,3) * t422 - qJ(6) * t418 - t574 * t440 + t583 * t441 - t573 * t465 + t566 * t480 + t565 * t497;
t469 = Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t527;
t470 = Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t527;
t485 = Ifges(4,4) * t513 + Ifges(4,2) * t512 + Ifges(4,6) * t527;
t486 = Ifges(4,1) * t513 + Ifges(4,4) * t512 + Ifges(4,5) * t527;
t580 = (Ifges(4,3) + Ifges(5,3)) * t516 + mrSges(4,1) * t442 + mrSges(5,1) * t428 + Ifges(4,5) * t491 + Ifges(5,5) * t467 + Ifges(4,6) * t490 + Ifges(5,6) * t466 + pkin(3) * t398 + pkin(4) * t415 + pkin(10) * t407 + t499 * t469 + t513 * t485 + t537 * t404 + t578 * t403 - mrSges(4,2) * t443 - mrSges(5,2) * t429 - t498 * t470 - t512 * t486;
t417 = t441 * mrSges(7,2) + t481 * t456 - t550;
t579 = -t440 * t572 - t441 * t573 - t581 * t465 - t480 * t564 - t481 * t565 + mrSges(6,1) * t423 - mrSges(7,1) * t421 - mrSges(6,2) * t424 + mrSges(7,3) * t420 - pkin(5) * t417 + qJ(6) * (-t440 * mrSges(7,2) - t480 * t456 + t557);
t568 = t534 * t539;
t500 = -mrSges(4,1) * t512 + mrSges(4,2) * t513;
t501 = -mrSges(4,2) * t527 + mrSges(4,3) * t512;
t396 = m(4) * t442 + mrSges(4,1) * t516 - mrSges(4,3) * t491 - t500 * t513 + t501 * t527 + t398;
t503 = mrSges(4,1) * t527 - mrSges(4,3) * t513;
t552 = t535 * t402 - t409 * t533;
t397 = m(4) * t443 - mrSges(4,2) * t516 + mrSges(4,3) * t490 + t500 * t512 - t503 * t527 + t552;
t390 = t541 * t396 + t538 * t397;
t406 = t537 * t412 + t414 * t578;
t553 = -t396 * t538 + t541 * t397;
t405 = m(5) * t437 - t466 * mrSges(5,1) + mrSges(5,2) * t467 - t498 * t482 + t483 * t499 + t406;
t546 = -m(4) * t477 + t490 * mrSges(4,1) - mrSges(4,2) * t491 + t512 * t501 - t503 * t513 - t405;
t521 = (-mrSges(3,1) * t542 + mrSges(3,2) * t539) * t561;
t518 = -mrSges(3,2) * t530 + mrSges(3,3) * t555;
t517 = mrSges(3,1) * t530 - mrSges(3,3) * t556;
t507 = -t534 * t519 - t576;
t506 = Ifges(3,5) * t530 + (Ifges(3,1) * t539 + Ifges(3,4) * t542) * t561;
t505 = Ifges(3,6) * t530 + (Ifges(3,4) * t539 + Ifges(3,2) * t542) * t561;
t504 = Ifges(3,3) * t530 + (Ifges(3,5) * t539 + Ifges(3,6) * t542) * t561;
t493 = -g(3) * t568 + t562;
t484 = Ifges(4,5) * t513 + Ifges(4,6) * t512 + Ifges(4,3) * t527;
t468 = Ifges(5,5) * t499 + Ifges(5,6) * t498 + Ifges(5,3) * t527;
t399 = m(3) * t492 + mrSges(3,1) * t529 - mrSges(3,3) * t523 + t518 * t530 - t521 * t556 + t546;
t392 = -mrSges(5,1) * t437 + mrSges(5,3) * t429 + Ifges(5,4) * t467 + Ifges(5,2) * t466 + Ifges(5,6) * t516 - pkin(4) * t406 - t499 * t468 + t527 * t470 - t579;
t391 = mrSges(5,2) * t437 - mrSges(5,3) * t428 + Ifges(5,1) * t467 + Ifges(5,4) * t466 + Ifges(5,5) * t516 - pkin(10) * t406 - t537 * t403 + t404 * t578 + t498 * t468 - t527 * t469;
t389 = m(3) * t493 - mrSges(3,2) * t529 - mrSges(3,3) * t524 - t517 * t530 + t521 * t555 + t553;
t388 = mrSges(4,2) * t477 - mrSges(4,3) * t442 + Ifges(4,1) * t491 + Ifges(4,4) * t490 + Ifges(4,5) * t516 - qJ(4) * t398 + t391 * t535 - t392 * t533 + t484 * t512 - t485 * t527;
t387 = -mrSges(4,1) * t477 + mrSges(4,3) * t443 + Ifges(4,4) * t491 + Ifges(4,2) * t490 + Ifges(4,6) * t516 - pkin(3) * t405 + qJ(4) * t552 + t533 * t391 + t535 * t392 - t513 * t484 + t527 * t486;
t386 = Ifges(3,5) * t523 - Ifges(3,6) * t524 + Ifges(3,3) * t529 + mrSges(3,1) * t492 - mrSges(3,2) * t493 + t538 * t388 + t541 * t387 + pkin(2) * t546 + pkin(9) * t553 + (t505 * t539 - t506 * t542) * t561;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t554 - mrSges(2,2) * t549 + (mrSges(3,2) * t507 - mrSges(3,3) * t492 + Ifges(3,1) * t523 - Ifges(3,4) * t524 + Ifges(3,5) * t529 - pkin(9) * t390 - t387 * t538 + t388 * t541 + t504 * t555 - t505 * t530) * t568 + (-mrSges(3,1) * t507 + mrSges(3,3) * t493 + Ifges(3,4) * t523 - Ifges(3,2) * t524 + Ifges(3,6) * t529 - pkin(2) * t390 - t504 * t556 + t530 * t506 - t580) * t567 + t536 * t386 + pkin(1) * ((t389 * t539 + t399 * t542) * t536 + (-m(3) * t507 - t524 * mrSges(3,1) - t523 * mrSges(3,2) + (-t517 * t539 + t518 * t542) * t561 - t390) * t534) + (t389 * t542 - t399 * t539) * t577; t386; t580; t405; t579; t417;];
tauJ  = t1;
