% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-06 22:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:16:58
% EndTime: 2019-05-06 22:17:07
% DurationCPUTime: 4.75s
% Computational Cost: add. (45682->335), mult. (92685->408), div. (0->0), fcn. (60702->10), ass. (0->134)
t591 = Ifges(3,1) + Ifges(4,1);
t586 = Ifges(3,4) - Ifges(4,5);
t585 = Ifges(3,5) + Ifges(4,4);
t590 = Ifges(3,2) + Ifges(4,3);
t584 = Ifges(3,6) - Ifges(4,6);
t589 = Ifges(3,3) + Ifges(4,2);
t588 = 2 * qJD(3);
t587 = mrSges(3,3) + mrSges(4,2);
t558 = cos(qJ(2));
t561 = qJD(1) ^ 2;
t583 = t558 ^ 2 * t561;
t553 = sin(qJ(2));
t576 = qJD(1) * qJD(2);
t574 = t553 * t576;
t524 = qJDD(1) * t558 - t574;
t578 = qJD(1) * t553;
t531 = -qJD(2) * pkin(3) - pkin(8) * t578;
t554 = sin(qJ(1));
t559 = cos(qJ(1));
t579 = t554 * g(1) - t559 * g(2);
t514 = -qJDD(1) * pkin(1) - t561 * pkin(7) - t579;
t575 = t558 * t576;
t523 = qJDD(1) * t553 + t575;
t570 = -t524 * pkin(2) + t514 + (-t523 - t575) * qJ(3);
t453 = -pkin(2) * t574 + pkin(3) * t524 - pkin(8) * t583 - t570 + (t531 + t588) * t578;
t552 = sin(qJ(4));
t557 = cos(qJ(4));
t513 = (-t558 * t552 + t553 * t557) * qJD(1);
t480 = -t513 * qJD(4) - t552 * t523 - t524 * t557;
t577 = qJD(1) * t558;
t512 = -t552 * t578 - t557 * t577;
t481 = qJD(4) * t512 + t523 * t557 - t524 * t552;
t543 = -qJD(2) + qJD(4);
t441 = (-t512 * t543 - t481) * pkin(9) + (t513 * t543 - t480) * pkin(4) + t453;
t571 = -g(1) * t559 - g(2) * t554;
t515 = -pkin(1) * t561 + qJDD(1) * pkin(7) + t571;
t495 = -g(3) * t553 + t558 * t515;
t520 = (-pkin(2) * t558 - qJ(3) * t553) * qJD(1);
t560 = qJD(2) ^ 2;
t476 = -pkin(2) * t560 + qJDD(2) * qJ(3) + qJD(2) * t588 + t520 * t577 + t495;
t462 = -pkin(3) * t583 - pkin(8) * t524 + qJD(2) * t531 + t476;
t494 = -t558 * g(3) - t553 * t515;
t479 = -qJDD(2) * pkin(2) - qJ(3) * t560 + t520 * t578 + qJDD(3) - t494;
t463 = (-t523 + t575) * pkin(8) + (-t553 * t558 * t561 - qJDD(2)) * pkin(3) + t479;
t447 = t557 * t462 + t552 * t463;
t490 = -pkin(4) * t512 - pkin(9) * t513;
t541 = t543 ^ 2;
t542 = -qJDD(2) + qJDD(4);
t444 = -pkin(4) * t541 + pkin(9) * t542 + t490 * t512 + t447;
t551 = sin(qJ(5));
t556 = cos(qJ(5));
t431 = t556 * t441 - t444 * t551;
t492 = -t513 * t551 + t543 * t556;
t456 = qJD(5) * t492 + t481 * t556 + t542 * t551;
t478 = qJDD(5) - t480;
t493 = t513 * t556 + t543 * t551;
t505 = qJD(5) - t512;
t429 = (t492 * t505 - t456) * pkin(10) + (t492 * t493 + t478) * pkin(5) + t431;
t432 = t551 * t441 + t556 * t444;
t455 = -qJD(5) * t493 - t481 * t551 + t542 * t556;
t484 = pkin(5) * t505 - pkin(10) * t493;
t491 = t492 ^ 2;
t430 = -pkin(5) * t491 + pkin(10) * t455 - t484 * t505 + t432;
t550 = sin(qJ(6));
t555 = cos(qJ(6));
t427 = t429 * t555 - t430 * t550;
t469 = t492 * t555 - t493 * t550;
t438 = qJD(6) * t469 + t455 * t550 + t456 * t555;
t470 = t492 * t550 + t493 * t555;
t452 = -mrSges(7,1) * t469 + mrSges(7,2) * t470;
t500 = qJD(6) + t505;
t460 = -mrSges(7,2) * t500 + mrSges(7,3) * t469;
t474 = qJDD(6) + t478;
t423 = m(7) * t427 + mrSges(7,1) * t474 - mrSges(7,3) * t438 - t452 * t470 + t460 * t500;
t428 = t429 * t550 + t430 * t555;
t437 = -qJD(6) * t470 + t455 * t555 - t456 * t550;
t461 = mrSges(7,1) * t500 - mrSges(7,3) * t470;
t424 = m(7) * t428 - mrSges(7,2) * t474 + mrSges(7,3) * t437 + t452 * t469 - t461 * t500;
t416 = t555 * t423 + t550 * t424;
t471 = -mrSges(6,1) * t492 + mrSges(6,2) * t493;
t482 = -mrSges(6,2) * t505 + mrSges(6,3) * t492;
t414 = m(6) * t431 + mrSges(6,1) * t478 - mrSges(6,3) * t456 - t471 * t493 + t482 * t505 + t416;
t483 = mrSges(6,1) * t505 - mrSges(6,3) * t493;
t572 = -t423 * t550 + t555 * t424;
t415 = m(6) * t432 - mrSges(6,2) * t478 + mrSges(6,3) * t455 + t471 * t492 - t483 * t505 + t572;
t409 = t556 * t414 + t551 * t415;
t582 = t589 * qJD(2) + (t553 * t585 + t558 * t584) * qJD(1);
t581 = -t584 * qJD(2) + (-t553 * t586 - t590 * t558) * qJD(1);
t580 = t585 * qJD(2) + (t553 * t591 + t558 * t586) * qJD(1);
t410 = -t414 * t551 + t556 * t415;
t489 = -mrSges(5,1) * t512 + mrSges(5,2) * t513;
t497 = mrSges(5,1) * t543 - mrSges(5,3) * t513;
t408 = m(5) * t447 - mrSges(5,2) * t542 + mrSges(5,3) * t480 + t489 * t512 - t497 * t543 + t410;
t446 = -t552 * t462 + t463 * t557;
t443 = -pkin(4) * t542 - pkin(9) * t541 + t513 * t490 - t446;
t433 = -pkin(5) * t455 - pkin(10) * t491 + t484 * t493 + t443;
t568 = m(7) * t433 - t437 * mrSges(7,1) + mrSges(7,2) * t438 - t469 * t460 + t461 * t470;
t425 = -m(6) * t443 + t455 * mrSges(6,1) - mrSges(6,2) * t456 + t492 * t482 - t483 * t493 - t568;
t496 = -mrSges(5,2) * t543 + mrSges(5,3) * t512;
t419 = m(5) * t446 + mrSges(5,1) * t542 - mrSges(5,3) * t481 - t489 * t513 + t496 * t543 + t425;
t573 = t557 * t408 - t552 * t419;
t403 = t552 * t408 + t557 * t419;
t521 = (-t558 * mrSges(4,1) - t553 * mrSges(4,3)) * qJD(1);
t528 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t578;
t569 = m(4) * t476 + qJDD(2) * mrSges(4,3) + qJD(2) * t528 + t521 * t577 + t573;
t530 = mrSges(4,2) * t577 + qJD(2) * mrSges(4,3);
t567 = m(4) * t479 - qJDD(2) * mrSges(4,1) - qJD(2) * t530 + t403;
t449 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t500;
t450 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t500;
t566 = -mrSges(7,1) * t427 + mrSges(7,2) * t428 - Ifges(7,5) * t438 - Ifges(7,6) * t437 - Ifges(7,3) * t474 - t470 * t449 + t469 * t450;
t565 = m(5) * t453 - t480 * mrSges(5,1) + mrSges(5,2) * t481 - t512 * t496 + t497 * t513 + t409;
t468 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t578 + t570;
t564 = m(4) * t468 - t565;
t448 = Ifges(7,5) * t470 + Ifges(7,6) * t469 + Ifges(7,3) * t500;
t417 = -mrSges(7,1) * t433 + mrSges(7,3) * t428 + Ifges(7,4) * t438 + Ifges(7,2) * t437 + Ifges(7,6) * t474 - t448 * t470 + t450 * t500;
t418 = mrSges(7,2) * t433 - mrSges(7,3) * t427 + Ifges(7,1) * t438 + Ifges(7,4) * t437 + Ifges(7,5) * t474 + t448 * t469 - t449 * t500;
t464 = Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t505;
t466 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t505;
t404 = -mrSges(6,1) * t443 + mrSges(6,3) * t432 + Ifges(6,4) * t456 + Ifges(6,2) * t455 + Ifges(6,6) * t478 - pkin(5) * t568 + pkin(10) * t572 + t555 * t417 + t550 * t418 - t493 * t464 + t505 * t466;
t465 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t505;
t405 = mrSges(6,2) * t443 - mrSges(6,3) * t431 + Ifges(6,1) * t456 + Ifges(6,4) * t455 + Ifges(6,5) * t478 - pkin(10) * t416 - t417 * t550 + t418 * t555 + t464 * t492 - t465 * t505;
t486 = Ifges(5,4) * t513 + Ifges(5,2) * t512 + Ifges(5,6) * t543;
t487 = Ifges(5,1) * t513 + Ifges(5,4) * t512 + Ifges(5,5) * t543;
t563 = mrSges(5,1) * t446 - mrSges(5,2) * t447 + Ifges(5,5) * t481 + Ifges(5,6) * t480 + Ifges(5,3) * t542 + pkin(4) * t425 + pkin(9) * t410 + t556 * t404 + t551 * t405 + t513 * t486 - t512 * t487;
t562 = mrSges(6,1) * t431 - mrSges(6,2) * t432 + Ifges(6,5) * t456 + Ifges(6,6) * t455 + Ifges(6,3) * t478 + pkin(5) * t416 + t493 * t465 - t492 * t466 - t566;
t529 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t577;
t527 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t578;
t522 = (-t558 * mrSges(3,1) + t553 * mrSges(3,2)) * qJD(1);
t485 = Ifges(5,5) * t513 + Ifges(5,6) * t512 + Ifges(5,3) * t543;
t406 = -mrSges(4,1) * t524 - mrSges(4,3) * t523 + (-t528 * t553 - t530 * t558) * qJD(1) + t564;
t402 = t523 * mrSges(4,2) + t521 * t578 + t567;
t401 = -mrSges(5,1) * t453 + mrSges(5,3) * t447 + Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * t542 - pkin(4) * t409 - t513 * t485 + t543 * t487 - t562;
t400 = mrSges(5,2) * t453 - mrSges(5,3) * t446 + Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * t542 - pkin(9) * t409 - t404 * t551 + t405 * t556 + t485 * t512 - t486 * t543;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t579 - mrSges(2,2) * t571 + t553 * (mrSges(3,2) * t514 + mrSges(4,2) * t479 - mrSges(3,3) * t494 - mrSges(4,3) * t468 - pkin(8) * t403 - qJ(3) * t406 + t581 * qJD(2) + t585 * qJDD(2) + t557 * t400 - t552 * t401 + t591 * t523 + t586 * t524 + t582 * t577) + t558 * (-mrSges(3,1) * t514 - mrSges(4,1) * t468 + mrSges(4,2) * t476 + mrSges(3,3) * t495 - pkin(2) * t406 + pkin(3) * t565 - pkin(8) * t573 + t580 * qJD(2) + t584 * qJDD(2) - t552 * t400 - t557 * t401 + t586 * t523 + t590 * t524 - t582 * t578) + pkin(1) * (-m(3) * t514 + (mrSges(3,1) + mrSges(4,1)) * t524 + (-mrSges(3,2) + mrSges(4,3)) * t523 + ((t529 + t530) * t558 + (-t527 + t528) * t553) * qJD(1) - t564) + pkin(7) * (t558 * (m(3) * t495 - qJDD(2) * mrSges(3,2) - qJD(2) * t527 + t522 * t577 + t587 * t524 + t569) + (-m(3) * t494 - qJDD(2) * mrSges(3,1) - qJD(2) * t529 + t587 * t523 + (t521 + t522) * t578 + t567) * t553); (-t581 * t553 - t580 * t558) * qJD(1) + qJ(3) * t569 + (qJ(3) * mrSges(4,2) + t584) * t524 + t585 * t523 + mrSges(3,1) * t494 - mrSges(3,2) * t495 + mrSges(4,3) * t476 - mrSges(4,1) * t479 + t589 * qJDD(2) - pkin(2) * t402 - pkin(3) * t403 - t563; t402; t563; t562; -t566;];
tauJ  = t1;
