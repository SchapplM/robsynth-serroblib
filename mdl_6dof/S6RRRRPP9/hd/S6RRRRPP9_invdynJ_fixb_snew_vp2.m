% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:12:39
% EndTime: 2019-05-07 19:12:49
% DurationCPUTime: 4.59s
% Computational Cost: add. (45737->317), mult. (96897->386), div. (0->0), fcn. (74972->10), ass. (0->135)
t573 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t553 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t552 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t572 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t551 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t571 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t521 = sin(pkin(6));
t525 = sin(qJ(2));
t528 = cos(qJ(2));
t554 = qJD(1) * qJD(2);
t511 = (-qJDD(1) * t528 + t525 * t554) * t521;
t556 = qJD(1) * t521;
t509 = (-pkin(2) * t528 - pkin(9) * t525) * t556;
t522 = cos(pkin(6));
t518 = qJD(1) * t522 + qJD(2);
t516 = t518 ^ 2;
t517 = qJDD(1) * t522 + qJDD(2);
t555 = qJD(1) * t528;
t530 = qJD(1) ^ 2;
t526 = sin(qJ(1));
t529 = cos(qJ(1));
t540 = -g(1) * t529 - g(2) * t526;
t565 = pkin(8) * t521;
t507 = -pkin(1) * t530 + qJDD(1) * t565 + t540;
t543 = t526 * g(1) - g(2) * t529;
t506 = qJDD(1) * pkin(1) + t530 * t565 + t543;
t561 = t506 * t522;
t557 = t528 * t507 + t525 * t561;
t451 = -t516 * pkin(2) + t517 * pkin(9) + (-g(3) * t525 + t509 * t555) * t521 + t557;
t510 = (qJDD(1) * t525 + t528 * t554) * t521;
t564 = t522 * g(3);
t452 = t511 * pkin(2) - t510 * pkin(9) - t564 + (-t506 + (pkin(2) * t525 - pkin(9) * t528) * t518 * qJD(1)) * t521;
t524 = sin(qJ(3));
t527 = cos(qJ(3));
t424 = t527 * t451 + t524 * t452;
t545 = t525 * t556;
t499 = t527 * t518 - t524 * t545;
t500 = t518 * t524 + t527 * t545;
t483 = -pkin(3) * t499 - pkin(10) * t500;
t503 = qJDD(3) + t511;
t544 = t521 * t555;
t515 = qJD(3) - t544;
t514 = t515 ^ 2;
t420 = -pkin(3) * t514 + pkin(10) * t503 + t483 * t499 + t424;
t559 = t521 * t528;
t480 = -g(3) * t559 - t525 * t507 + t528 * t561;
t450 = -t517 * pkin(2) - t516 * pkin(9) + t509 * t545 - t480;
t478 = -qJD(3) * t500 - t510 * t524 + t517 * t527;
t479 = qJD(3) * t499 + t510 * t527 + t517 * t524;
t422 = (-t499 * t515 - t479) * pkin(10) + (t500 * t515 - t478) * pkin(3) + t450;
t523 = sin(qJ(4));
t567 = cos(qJ(4));
t415 = -t523 * t420 + t422 * t567;
t486 = t523 * t500 - t515 * t567;
t487 = t500 * t567 + t523 * t515;
t457 = pkin(4) * t486 - qJ(5) * t487;
t476 = qJDD(4) - t478;
t497 = qJD(4) - t499;
t496 = t497 ^ 2;
t413 = -t476 * pkin(4) - t496 * qJ(5) + t487 * t457 + qJDD(5) - t415;
t431 = -t486 * qJD(4) + t479 * t567 + t523 * t503;
t459 = -mrSges(6,2) * t486 - mrSges(6,3) * t487;
t570 = -m(6) * t413 - t431 * mrSges(6,1) - t487 * t459;
t456 = -mrSges(7,2) * t487 + mrSges(7,3) * t486;
t562 = t486 * t497;
t407 = -0.2e1 * qJD(6) * t497 + (t486 * t487 - t476) * qJ(6) + (t431 + t562) * pkin(5) + t413;
t464 = -mrSges(7,1) * t486 + mrSges(7,2) * t497;
t541 = -m(7) * t407 + t476 * mrSges(7,3) + t497 * t464;
t405 = t431 * mrSges(7,1) + t487 * t456 - t541;
t463 = mrSges(6,1) * t486 - mrSges(6,3) * t497;
t402 = t476 * mrSges(6,2) + t497 * t463 + t405 - t570;
t430 = t487 * qJD(4) + t523 * t479 - t503 * t567;
t461 = pkin(5) * t487 - qJ(6) * t497;
t485 = t486 ^ 2;
t416 = t567 * t420 + t523 * t422;
t535 = -t496 * pkin(4) + t476 * qJ(5) - t486 * t457 + t416;
t409 = -t430 * pkin(5) - t485 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t461) * t497 + t535;
t568 = -2 * qJD(5);
t412 = t497 * t568 - t535;
t465 = mrSges(6,1) * t487 + mrSges(6,2) * t497;
t462 = mrSges(7,1) * t487 - mrSges(7,3) * t497;
t549 = m(7) * t409 + t476 * mrSges(7,2) + t497 * t462;
t537 = -m(6) * t412 + t476 * mrSges(6,3) + t497 * t465 + t549;
t546 = -t553 * t486 + t573 * t487 + t552 * t497;
t547 = t572 * t486 + t553 * t487 + t551 * t497;
t558 = -t456 - t459;
t569 = -t430 * t551 + t431 * t552 + t571 * t476 + t486 * t546 + t487 * t547 + mrSges(5,1) * t415 - mrSges(5,2) * t416 + mrSges(6,2) * t413 + mrSges(7,2) * t409 - mrSges(6,3) * t412 - mrSges(7,3) * t407 - pkin(4) * t402 + qJ(5) * (t558 * t486 + (-mrSges(6,1) - mrSges(7,1)) * t430 + t537) - qJ(6) * t405;
t563 = -mrSges(7,1) - mrSges(5,3);
t560 = t521 * t525;
t458 = mrSges(5,1) * t486 + mrSges(5,2) * t487;
t466 = -mrSges(5,2) * t497 - mrSges(5,3) * t486;
t398 = m(5) * t415 + (-t463 + t466) * t497 + (-t456 - t458) * t487 + (mrSges(5,1) - mrSges(6,2)) * t476 + t563 * t431 + t541 + t570;
t467 = mrSges(5,1) * t497 - mrSges(5,3) * t487;
t400 = m(5) * t416 - t476 * mrSges(5,2) - t497 * t467 + (-t458 + t558) * t486 + (-mrSges(6,1) + t563) * t430 + t537;
t395 = -t398 * t523 + t567 * t400;
t482 = -mrSges(4,1) * t499 + mrSges(4,2) * t500;
t489 = mrSges(4,1) * t515 - mrSges(4,3) * t500;
t393 = m(4) * t424 - mrSges(4,2) * t503 + mrSges(4,3) * t478 + t482 * t499 - t489 * t515 + t395;
t423 = -t524 * t451 + t527 * t452;
t419 = -t503 * pkin(3) - t514 * pkin(10) + t500 * t483 - t423;
t533 = (-t431 + t562) * qJ(5) + t419 + (t497 * pkin(4) + t568) * t487;
t414 = t430 * pkin(4) + t533;
t411 = -t485 * pkin(5) + 0.2e1 * qJD(6) * t486 - t487 * t461 + (pkin(4) + qJ(6)) * t430 + t533;
t539 = m(7) * t411 - t431 * mrSges(7,2) + t430 * mrSges(7,3) - t487 * t462 + t486 * t464;
t536 = -m(6) * t414 + t430 * mrSges(6,2) + t486 * t463 - t539;
t401 = -m(5) * t419 - t430 * mrSges(5,1) - t486 * t466 + (t465 - t467) * t487 + (-mrSges(5,2) + mrSges(6,3)) * t431 + t536;
t488 = -mrSges(4,2) * t515 + mrSges(4,3) * t499;
t397 = m(4) * t423 + t503 * mrSges(4,1) - t479 * mrSges(4,3) - t500 * t482 + t515 * t488 + t401;
t388 = t524 * t393 + t527 * t397;
t548 = t551 * t486 - t552 * t487 - t571 * t497;
t542 = t527 * t393 - t397 * t524;
t394 = t398 * t567 + t523 * t400;
t534 = -m(4) * t450 + t478 * mrSges(4,1) - t479 * mrSges(4,2) + t499 * t488 - t500 * t489 - t394;
t404 = -t431 * mrSges(6,3) - t487 * t465 - t536;
t406 = -t430 * mrSges(7,1) - t486 * t456 + t549;
t386 = -mrSges(5,1) * t419 - mrSges(6,1) * t412 + mrSges(7,1) * t409 + mrSges(6,2) * t414 + mrSges(5,3) * t416 - mrSges(7,3) * t411 - pkin(4) * t404 + pkin(5) * t406 - qJ(6) * t539 + t572 * t430 + t553 * t431 + t551 * t476 + t548 * t487 + t546 * t497;
t389 = mrSges(6,1) * t413 + mrSges(7,1) * t407 + mrSges(5,2) * t419 - mrSges(7,2) * t411 - mrSges(5,3) * t415 - mrSges(6,3) * t414 + pkin(5) * t405 - qJ(5) * t404 - t553 * t430 + t573 * t431 + t552 * t476 + t548 * t486 - t547 * t497;
t473 = Ifges(4,4) * t500 + Ifges(4,2) * t499 + Ifges(4,6) * t515;
t474 = Ifges(4,1) * t500 + Ifges(4,4) * t499 + Ifges(4,5) * t515;
t531 = mrSges(4,1) * t423 - mrSges(4,2) * t424 + Ifges(4,5) * t479 + Ifges(4,6) * t478 + Ifges(4,3) * t503 + pkin(3) * t401 + pkin(10) * t395 + t386 * t567 + t523 * t389 + t500 * t473 - t499 * t474;
t508 = (-mrSges(3,1) * t528 + mrSges(3,2) * t525) * t556;
t505 = -mrSges(3,2) * t518 + mrSges(3,3) * t544;
t504 = mrSges(3,1) * t518 - mrSges(3,3) * t545;
t493 = -t521 * t506 - t564;
t492 = Ifges(3,5) * t518 + (Ifges(3,1) * t525 + Ifges(3,4) * t528) * t556;
t491 = Ifges(3,6) * t518 + (Ifges(3,4) * t525 + Ifges(3,2) * t528) * t556;
t490 = Ifges(3,3) * t518 + (Ifges(3,5) * t525 + Ifges(3,6) * t528) * t556;
t481 = -g(3) * t560 + t557;
t472 = Ifges(4,5) * t500 + Ifges(4,6) * t499 + Ifges(4,3) * t515;
t390 = m(3) * t480 + t517 * mrSges(3,1) - t510 * mrSges(3,3) + t518 * t505 - t508 * t545 + t534;
t387 = m(3) * t481 - mrSges(3,2) * t517 - mrSges(3,3) * t511 - t504 * t518 + t508 * t544 + t542;
t385 = -mrSges(4,1) * t450 + mrSges(4,3) * t424 + Ifges(4,4) * t479 + Ifges(4,2) * t478 + Ifges(4,6) * t503 - pkin(3) * t394 - t500 * t472 + t515 * t474 - t569;
t384 = mrSges(4,2) * t450 - mrSges(4,3) * t423 + Ifges(4,1) * t479 + Ifges(4,4) * t478 + Ifges(4,5) * t503 - pkin(10) * t394 - t523 * t386 + t389 * t567 + t499 * t472 - t515 * t473;
t383 = Ifges(3,5) * t510 - Ifges(3,6) * t511 + Ifges(3,3) * t517 + mrSges(3,1) * t480 - mrSges(3,2) * t481 + t524 * t384 + t527 * t385 + pkin(2) * t534 + pkin(9) * t542 + (t491 * t525 - t492 * t528) * t556;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t543 - mrSges(2,2) * t540 + (mrSges(3,2) * t493 - mrSges(3,3) * t480 + Ifges(3,1) * t510 - Ifges(3,4) * t511 + Ifges(3,5) * t517 - pkin(9) * t388 + t384 * t527 - t385 * t524 + t490 * t544 - t491 * t518) * t560 + (-mrSges(3,1) * t493 + mrSges(3,3) * t481 + Ifges(3,4) * t510 - Ifges(3,2) * t511 + Ifges(3,6) * t517 - pkin(2) * t388 - t490 * t545 + t518 * t492 - t531) * t559 + t522 * t383 + pkin(1) * ((t387 * t525 + t390 * t528) * t522 + (-m(3) * t493 - t511 * mrSges(3,1) - t510 * mrSges(3,2) + (-t504 * t525 + t505 * t528) * t556 - t388) * t521) + (t387 * t528 - t390 * t525) * t565; t383; t531; t569; t402; t406;];
tauJ  = t1;
