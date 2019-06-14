% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 09:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:31:35
% EndTime: 2019-05-06 09:31:41
% DurationCPUTime: 4.01s
% Computational Cost: add. (34864->328), mult. (81459->400), div. (0->0), fcn. (57004->10), ass. (0->132)
t587 = -2 * qJD(3);
t586 = Ifges(4,1) + Ifges(5,1);
t580 = Ifges(4,4) - Ifges(5,5);
t579 = Ifges(4,5) + Ifges(5,4);
t585 = Ifges(4,2) + Ifges(5,3);
t584 = -Ifges(5,2) - Ifges(4,3);
t578 = Ifges(4,6) - Ifges(5,6);
t543 = sin(qJ(2));
t547 = cos(qJ(2));
t563 = qJD(1) * qJD(2);
t521 = qJDD(1) * t543 + t547 * t563;
t550 = qJD(1) ^ 2;
t544 = sin(qJ(1));
t548 = cos(qJ(1));
t560 = -g(1) * t548 - g(2) * t544;
t517 = -pkin(1) * t550 + qJDD(1) * pkin(7) + t560;
t575 = t517 * t543;
t467 = qJDD(2) * pkin(2) - qJ(3) * t521 - t575 + (pkin(2) * t543 * t550 + qJ(3) * t563 - g(3)) * t547;
t498 = -g(3) * t543 + t547 * t517;
t522 = qJDD(1) * t547 - t543 * t563;
t567 = qJD(1) * t543;
t523 = qJD(2) * pkin(2) - qJ(3) * t567;
t574 = t547 ^ 2 * t550;
t468 = -pkin(2) * t574 + qJ(3) * t522 - qJD(2) * t523 + t498;
t539 = sin(pkin(10));
t576 = cos(pkin(10));
t511 = (t539 * t547 + t576 * t543) * qJD(1);
t449 = t576 * t467 - t539 * t468 + t511 * t587;
t568 = t544 * g(1) - t548 * g(2);
t516 = -qJDD(1) * pkin(1) - t550 * pkin(7) - t568;
t471 = -t522 * pkin(2) - qJ(3) * t574 + t523 * t567 + qJDD(3) + t516;
t493 = t521 * t539 - t576 * t522;
t494 = t576 * t521 + t539 * t522;
t566 = qJD(1) * t547;
t510 = t539 * t567 - t576 * t566;
t499 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t510;
t500 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t511;
t501 = -qJD(2) * mrSges(5,1) + mrSges(5,2) * t511;
t565 = qJD(2) * t510;
t556 = t493 * pkin(3) + t471 + (-t494 + t565) * qJ(4);
t577 = pkin(3) * qJD(2);
t444 = (-(2 * qJD(4)) + t577) * t511 + t556;
t502 = -mrSges(5,2) * t510 + qJD(2) * mrSges(5,3);
t487 = pkin(3) * t510 - qJ(4) * t511;
t549 = qJD(2) ^ 2;
t439 = -qJDD(2) * pkin(3) - t549 * qJ(4) + t511 * t487 + qJDD(4) - t449;
t433 = (-t494 - t565) * pkin(8) + (t510 * t511 - qJDD(2)) * pkin(4) + t439;
t450 = t539 * t467 + t576 * t468 + t510 * t587;
t583 = 2 * qJD(4);
t438 = -pkin(3) * t549 + qJDD(2) * qJ(4) + qJD(2) * t583 - t510 * t487 + t450;
t503 = -qJD(2) * pkin(4) - pkin(8) * t511;
t509 = t510 ^ 2;
t435 = -pkin(4) * t509 + pkin(8) * t493 + qJD(2) * t503 + t438;
t542 = sin(qJ(5));
t546 = cos(qJ(5));
t431 = t542 * t433 + t546 * t435;
t483 = t510 * t546 - t511 * t542;
t484 = t510 * t542 + t511 * t546;
t463 = -pkin(5) * t483 - pkin(9) * t484;
t533 = -qJD(2) + qJD(5);
t531 = t533 ^ 2;
t532 = -qJDD(2) + qJDD(5);
t428 = -pkin(5) * t531 + pkin(9) * t532 + t463 * t483 + t431;
t436 = -pkin(4) * t493 - pkin(8) * t509 - t556 + (t503 - t577 + t583) * t511;
t453 = -qJD(5) * t484 + t493 * t546 - t494 * t542;
t454 = qJD(5) * t483 + t493 * t542 + t494 * t546;
t429 = (t484 * t533 - t453) * pkin(5) + (-t483 * t533 - t454) * pkin(9) + t436;
t541 = sin(qJ(6));
t545 = cos(qJ(6));
t425 = -t428 * t541 + t429 * t545;
t472 = -t484 * t541 + t533 * t545;
t442 = qJD(6) * t472 + t454 * t545 + t532 * t541;
t452 = qJDD(6) - t453;
t473 = t484 * t545 + t533 * t541;
t455 = -mrSges(7,1) * t472 + mrSges(7,2) * t473;
t476 = qJD(6) - t483;
t456 = -mrSges(7,2) * t476 + mrSges(7,3) * t472;
t422 = m(7) * t425 + mrSges(7,1) * t452 - mrSges(7,3) * t442 - t455 * t473 + t456 * t476;
t426 = t428 * t545 + t429 * t541;
t441 = -qJD(6) * t473 - t454 * t541 + t532 * t545;
t457 = mrSges(7,1) * t476 - mrSges(7,3) * t473;
t423 = m(7) * t426 - mrSges(7,2) * t452 + mrSges(7,3) * t441 + t455 * t472 - t457 * t476;
t414 = t545 * t422 + t541 * t423;
t474 = -mrSges(6,2) * t533 + mrSges(6,3) * t483;
t475 = mrSges(6,1) * t533 - mrSges(6,3) * t484;
t555 = m(6) * t436 - t453 * mrSges(6,1) + t454 * mrSges(6,2) - t483 * t474 + t484 * t475 + t414;
t554 = -m(5) * t444 - t493 * mrSges(5,1) - t510 * t502 + t555;
t410 = m(4) * t471 + t493 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t494 + t510 * t499 + (t500 - t501) * t511 - t554;
t581 = -mrSges(4,3) - mrSges(5,2);
t415 = -t422 * t541 + t545 * t423;
t462 = -mrSges(6,1) * t483 + mrSges(6,2) * t484;
t413 = m(6) * t431 - mrSges(6,2) * t532 + mrSges(6,3) * t453 + t462 * t483 - t475 * t533 + t415;
t430 = t433 * t546 - t435 * t542;
t427 = -pkin(5) * t532 - pkin(9) * t531 + t463 * t484 - t430;
t424 = -m(7) * t427 + t441 * mrSges(7,1) - mrSges(7,2) * t442 + t472 * t456 - t457 * t473;
t418 = m(6) * t430 + mrSges(6,1) * t532 - mrSges(6,3) * t454 - t462 * t484 + t474 * t533 + t424;
t561 = t546 * t413 - t542 * t418;
t558 = m(5) * t438 + qJDD(2) * mrSges(5,3) + qJD(2) * t501 + t561;
t488 = mrSges(5,1) * t510 - mrSges(5,3) * t511;
t570 = -mrSges(4,1) * t510 - mrSges(4,2) * t511 - t488;
t406 = m(4) * t450 - qJDD(2) * mrSges(4,2) - qJD(2) * t500 + t581 * t493 + t570 * t510 + t558;
t409 = t413 * t542 + t418 * t546;
t557 = -m(5) * t439 + qJDD(2) * mrSges(5,1) + qJD(2) * t502 - t409;
t407 = m(4) * t449 + qJDD(2) * mrSges(4,1) + qJD(2) * t499 + t581 * t494 + t570 * t511 + t557;
t400 = t539 * t406 + t576 * t407;
t573 = -t578 * qJD(2) + t585 * t510 - t580 * t511;
t572 = t584 * qJD(2) + t578 * t510 - t579 * t511;
t571 = t579 * qJD(2) - t580 * t510 + t586 * t511;
t562 = t576 * t406 - t407 * t539;
t446 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t476;
t447 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t476;
t553 = mrSges(7,1) * t425 - mrSges(7,2) * t426 + Ifges(7,5) * t442 + Ifges(7,6) * t441 + Ifges(7,3) * t452 + t446 * t473 - t472 * t447;
t445 = Ifges(7,5) * t473 + Ifges(7,6) * t472 + Ifges(7,3) * t476;
t416 = -mrSges(7,1) * t427 + mrSges(7,3) * t426 + Ifges(7,4) * t442 + Ifges(7,2) * t441 + Ifges(7,6) * t452 - t445 * t473 + t447 * t476;
t417 = mrSges(7,2) * t427 - mrSges(7,3) * t425 + Ifges(7,1) * t442 + Ifges(7,4) * t441 + Ifges(7,5) * t452 + t445 * t472 - t446 * t476;
t459 = Ifges(6,4) * t484 + Ifges(6,2) * t483 + Ifges(6,6) * t533;
t460 = Ifges(6,1) * t484 + Ifges(6,4) * t483 + Ifges(6,5) * t533;
t551 = mrSges(6,1) * t430 - mrSges(6,2) * t431 + Ifges(6,5) * t454 + Ifges(6,6) * t453 + Ifges(6,3) * t532 + pkin(5) * t424 + pkin(9) * t415 + t545 * t416 + t541 * t417 + t484 * t459 - t483 * t460;
t525 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t566;
t524 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t567;
t520 = (-mrSges(3,1) * t547 + mrSges(3,2) * t543) * qJD(1);
t514 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t543 + Ifges(3,4) * t547) * qJD(1);
t513 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t543 + Ifges(3,2) * t547) * qJD(1);
t497 = -g(3) * t547 - t575;
t458 = Ifges(6,5) * t484 + Ifges(6,6) * t483 + Ifges(6,3) * t533;
t411 = -t494 * mrSges(5,3) - t511 * t501 - t554;
t408 = mrSges(5,2) * t494 + t488 * t511 - t557;
t402 = -mrSges(6,1) * t436 + mrSges(6,3) * t431 + Ifges(6,4) * t454 + Ifges(6,2) * t453 + Ifges(6,6) * t532 - pkin(5) * t414 - t458 * t484 + t460 * t533 - t553;
t401 = mrSges(6,2) * t436 - mrSges(6,3) * t430 + Ifges(6,1) * t454 + Ifges(6,4) * t453 + Ifges(6,5) * t532 - pkin(9) * t414 - t416 * t541 + t417 * t545 + t458 * t483 - t459 * t533;
t399 = mrSges(4,2) * t471 + mrSges(5,2) * t439 - mrSges(4,3) * t449 - mrSges(5,3) * t444 - pkin(8) * t409 - qJ(4) * t411 + t573 * qJD(2) + t579 * qJDD(2) + t401 * t546 - t402 * t542 - t580 * t493 + t586 * t494 + t572 * t510;
t398 = -mrSges(4,1) * t471 - mrSges(5,1) * t444 + mrSges(5,2) * t438 + mrSges(4,3) * t450 - pkin(3) * t411 + pkin(4) * t555 - pkin(8) * t561 + t571 * qJD(2) + t578 * qJDD(2) - t542 * t401 - t546 * t402 - t585 * t493 + t580 * t494 + t572 * t511;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t568 - mrSges(2,2) * t560 + t543 * (mrSges(3,2) * t516 - mrSges(3,3) * t497 + Ifges(3,1) * t521 + Ifges(3,4) * t522 + Ifges(3,5) * qJDD(2) - qJ(3) * t400 - qJD(2) * t513 - t539 * t398 + t576 * t399) + t547 * (-mrSges(3,1) * t516 + mrSges(3,3) * t498 + Ifges(3,4) * t521 + Ifges(3,2) * t522 + Ifges(3,6) * qJDD(2) - pkin(2) * t410 + qJ(3) * t562 + qJD(2) * t514 + t576 * t398 + t539 * t399) + pkin(1) * (t522 * mrSges(3,1) - m(3) * t516 - t521 * mrSges(3,2) + (-t524 * t543 + t525 * t547) * qJD(1) - t410) + pkin(7) * (t547 * (m(3) * t498 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t522 - qJD(2) * t524 + t520 * t566 + t562) - t543 * (m(3) * t497 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t521 + qJD(2) * t525 - t520 * t567 + t400)); -t551 + (Ifges(3,3) - t584) * qJDD(2) + Ifges(3,5) * t521 + Ifges(3,6) * t522 + mrSges(3,1) * t497 - mrSges(3,2) * t498 + (-mrSges(5,2) * qJ(4) - t578) * t493 + t579 * t494 - t573 * t511 + mrSges(4,1) * t449 - mrSges(4,2) * t450 + mrSges(5,3) * t438 - mrSges(5,1) * t439 + (t543 * t513 - t547 * t514) * qJD(1) - pkin(4) * t409 - pkin(3) * t408 + pkin(2) * t400 + qJ(4) * t558 + (-qJ(4) * t488 + t571) * t510; t410; t408; t551; t553;];
tauJ  = t1;
