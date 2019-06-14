% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:20:34
% EndTime: 2019-05-06 10:20:40
% DurationCPUTime: 5.72s
% Computational Cost: add. (52638->341), mult. (141622->428), div. (0->0), fcn. (108403->12), ass. (0->149)
t614 = -2 * qJD(4);
t613 = Ifges(4,1) + Ifges(5,2);
t612 = -Ifges(5,1) - Ifges(4,3);
t607 = Ifges(4,5) - Ifges(5,4);
t611 = Ifges(4,2) + Ifges(5,3);
t606 = Ifges(4,6) - Ifges(5,5);
t605 = -Ifges(5,6) - Ifges(4,4);
t557 = sin(pkin(11));
t566 = cos(qJ(2));
t558 = sin(pkin(6));
t593 = qJD(1) * t558;
t586 = t566 * t593;
t562 = sin(qJ(2));
t587 = t562 * t593;
t604 = cos(pkin(11));
t535 = t557 * t587 - t604 * t586;
t536 = (t557 * t566 + t604 * t562) * t593;
t501 = pkin(3) * t535 - qJ(4) * t536;
t559 = cos(pkin(6));
t553 = qJD(1) * t559 + qJD(2);
t551 = t553 ^ 2;
t552 = qJDD(1) * t559 + qJDD(2);
t589 = qJD(1) * qJD(2);
t545 = (qJDD(1) * t562 + t566 * t589) * t558;
t568 = qJD(1) ^ 2;
t563 = sin(qJ(1));
t567 = cos(qJ(1));
t580 = -g(1) * t567 - g(2) * t563;
t609 = pkin(8) * t558;
t543 = -pkin(1) * t568 + qJDD(1) * t609 + t580;
t585 = t563 * g(1) - g(2) * t567;
t542 = qJDD(1) * pkin(1) + t568 * t609 + t585;
t602 = t542 * t559;
t581 = -t543 * t562 + t566 * t602;
t601 = t558 ^ 2 * t568;
t466 = pkin(2) * t552 - qJ(3) * t545 + (pkin(2) * t562 * t601 + (qJ(3) * qJD(1) * t553 - g(3)) * t558) * t566 + t581;
t600 = t558 * t562;
t500 = -g(3) * t600 + t566 * t543 + t562 * t602;
t539 = pkin(2) * t553 - qJ(3) * t587;
t546 = (qJDD(1) * t566 - t562 * t589) * t558;
t588 = t566 ^ 2 * t601;
t469 = -pkin(2) * t588 + qJ(3) * t546 - t539 * t553 + t500;
t598 = t557 * t466 + t604 * t469;
t610 = pkin(3) * t551 - t552 * qJ(4) + t535 * t501 + t553 * t614 - t598;
t452 = -0.2e1 * qJD(3) * t536 + t604 * t466 - t557 * t469;
t525 = -g(3) * t559 - t542 * t558;
t481 = -pkin(2) * t546 - qJ(3) * t588 + t539 * t587 + qJDD(3) + t525;
t509 = t545 * t557 - t604 * t546;
t510 = t604 * t545 + t557 * t546;
t517 = mrSges(4,1) * t553 - mrSges(4,3) * t536;
t603 = t535 * t553;
t569 = (-t510 + t603) * qJ(4) + t481 + (pkin(3) * t553 + t614) * t536;
t451 = pkin(3) * t509 + t569;
t519 = mrSges(5,1) * t536 + mrSges(5,2) * t553;
t449 = -t552 * pkin(3) - t551 * qJ(4) + t536 * t501 + qJDD(4) - t452;
t443 = (t535 * t536 - t552) * pkin(9) + (t510 + t603) * pkin(4) + t449;
t520 = pkin(4) * t536 - pkin(9) * t553;
t534 = t535 ^ 2;
t447 = -pkin(4) * t534 - t520 * t536 + (pkin(3) + pkin(9)) * t509 + t569;
t561 = sin(qJ(5));
t565 = cos(qJ(5));
t440 = t561 * t443 + t565 * t447;
t515 = t535 * t561 + t553 * t565;
t478 = -qJD(5) * t515 + t509 * t565 - t552 * t561;
t514 = t535 * t565 - t553 * t561;
t482 = -mrSges(6,1) * t514 + mrSges(6,2) * t515;
t533 = qJD(5) + t536;
t490 = mrSges(6,1) * t533 - mrSges(6,3) * t515;
t508 = qJDD(5) + t510;
t483 = -pkin(5) * t514 - pkin(10) * t515;
t532 = t533 ^ 2;
t437 = -pkin(5) * t532 + pkin(10) * t508 + t483 * t514 + t440;
t592 = qJD(3) * t535;
t528 = -0.2e1 * t592;
t445 = -pkin(4) * t509 - pkin(9) * t534 + t553 * t520 + t528 - t610;
t479 = qJD(5) * t514 + t509 * t561 + t552 * t565;
t441 = (-t514 * t533 - t479) * pkin(10) + (t515 * t533 - t478) * pkin(5) + t445;
t560 = sin(qJ(6));
t564 = cos(qJ(6));
t434 = -t437 * t560 + t441 * t564;
t487 = -t515 * t560 + t533 * t564;
t456 = qJD(6) * t487 + t479 * t564 + t508 * t560;
t488 = t515 * t564 + t533 * t560;
t462 = -mrSges(7,1) * t487 + mrSges(7,2) * t488;
t513 = qJD(6) - t514;
t467 = -mrSges(7,2) * t513 + mrSges(7,3) * t487;
t477 = qJDD(6) - t478;
t431 = m(7) * t434 + mrSges(7,1) * t477 - mrSges(7,3) * t456 - t462 * t488 + t467 * t513;
t435 = t437 * t564 + t441 * t560;
t455 = -qJD(6) * t488 - t479 * t560 + t508 * t564;
t468 = mrSges(7,1) * t513 - mrSges(7,3) * t488;
t432 = m(7) * t435 - mrSges(7,2) * t477 + mrSges(7,3) * t455 + t462 * t487 - t468 * t513;
t582 = -t431 * t560 + t564 * t432;
t420 = m(6) * t440 - mrSges(6,2) * t508 + mrSges(6,3) * t478 + t482 * t514 - t490 * t533 + t582;
t439 = t443 * t565 - t447 * t561;
t489 = -mrSges(6,2) * t533 + mrSges(6,3) * t514;
t436 = -pkin(5) * t508 - pkin(10) * t532 + t483 * t515 - t439;
t576 = -m(7) * t436 + t455 * mrSges(7,1) - mrSges(7,2) * t456 + t487 * t467 - t468 * t488;
t427 = m(6) * t439 + mrSges(6,1) * t508 - mrSges(6,3) * t479 - t482 * t515 + t489 * t533 + t576;
t583 = t565 * t420 - t561 * t427;
t577 = -m(5) * t451 + t510 * mrSges(5,3) + t536 * t519 - t583;
t518 = mrSges(5,1) * t535 - mrSges(5,3) * t553;
t594 = -mrSges(4,2) * t553 - mrSges(4,3) * t535 - t518;
t608 = mrSges(4,1) - mrSges(5,2);
t412 = m(4) * t481 + t510 * mrSges(4,2) + t608 * t509 + t536 * t517 + t594 * t535 - t577;
t599 = t558 * t566;
t502 = mrSges(4,1) * t535 + mrSges(4,2) * t536;
t415 = t420 * t561 + t427 * t565;
t503 = -mrSges(5,2) * t535 - mrSges(5,3) * t536;
t575 = -m(5) * t449 - t510 * mrSges(5,1) - t536 * t503 - t415;
t411 = m(4) * t452 - mrSges(4,3) * t510 - t502 * t536 + t608 * t552 + t594 * t553 + t575;
t453 = t528 + t598;
t448 = 0.2e1 * t592 + t610;
t422 = t564 * t431 + t560 * t432;
t573 = -m(6) * t445 + mrSges(6,1) * t478 - t479 * mrSges(6,2) + t489 * t514 - t515 * t490 - t422;
t571 = -m(5) * t448 + t552 * mrSges(5,3) + t553 * t519 - t573;
t418 = (-mrSges(4,3) - mrSges(5,1)) * t509 + (-t502 - t503) * t535 + t571 + m(4) * t453 - mrSges(4,2) * t552 - t517 * t553;
t407 = t604 * t411 + t557 * t418;
t597 = t606 * t535 - t607 * t536 + t612 * t553;
t596 = t611 * t535 + t605 * t536 - t606 * t553;
t595 = t605 * t535 + t613 * t536 + t607 * t553;
t584 = -t411 * t557 + t604 * t418;
t457 = Ifges(7,5) * t488 + Ifges(7,6) * t487 + Ifges(7,3) * t513;
t459 = Ifges(7,1) * t488 + Ifges(7,4) * t487 + Ifges(7,5) * t513;
t425 = -mrSges(7,1) * t436 + mrSges(7,3) * t435 + Ifges(7,4) * t456 + Ifges(7,2) * t455 + Ifges(7,6) * t477 - t457 * t488 + t459 * t513;
t458 = Ifges(7,4) * t488 + Ifges(7,2) * t487 + Ifges(7,6) * t513;
t426 = mrSges(7,2) * t436 - mrSges(7,3) * t434 + Ifges(7,1) * t456 + Ifges(7,4) * t455 + Ifges(7,5) * t477 + t457 * t487 - t458 * t513;
t471 = Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t533;
t472 = Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t533;
t572 = mrSges(6,1) * t439 - mrSges(6,2) * t440 + Ifges(6,5) * t479 + Ifges(6,6) * t478 + Ifges(6,3) * t508 + pkin(5) * t576 + pkin(10) * t582 + t564 * t425 + t560 * t426 + t515 * t471 - t514 * t472;
t570 = mrSges(7,1) * t434 - mrSges(7,2) * t435 + Ifges(7,5) * t456 + Ifges(7,6) * t455 + Ifges(7,3) * t477 + t458 * t488 - t459 * t487;
t544 = (-mrSges(3,1) * t566 + mrSges(3,2) * t562) * t593;
t541 = -mrSges(3,2) * t553 + mrSges(3,3) * t586;
t540 = mrSges(3,1) * t553 - mrSges(3,3) * t587;
t524 = Ifges(3,5) * t553 + (Ifges(3,1) * t562 + Ifges(3,4) * t566) * t593;
t523 = Ifges(3,6) * t553 + (Ifges(3,4) * t562 + Ifges(3,2) * t566) * t593;
t522 = Ifges(3,3) * t553 + (Ifges(3,5) * t562 + Ifges(3,6) * t566) * t593;
t499 = -g(3) * t599 + t581;
t470 = Ifges(6,5) * t515 + Ifges(6,6) * t514 + Ifges(6,3) * t533;
t414 = -t509 * mrSges(5,2) - t535 * t518 - t577;
t413 = mrSges(5,2) * t552 + t518 * t553 - t575;
t409 = -mrSges(6,1) * t445 + mrSges(6,3) * t440 + Ifges(6,4) * t479 + Ifges(6,2) * t478 + Ifges(6,6) * t508 - pkin(5) * t422 - t470 * t515 + t472 * t533 - t570;
t408 = mrSges(6,2) * t445 - mrSges(6,3) * t439 + Ifges(6,1) * t479 + Ifges(6,4) * t478 + Ifges(6,5) * t508 - pkin(10) * t422 - t425 * t560 + t426 * t564 + t470 * t514 - t471 * t533;
t406 = m(3) * t500 - mrSges(3,2) * t552 + mrSges(3,3) * t546 - t540 * t553 + t544 * t586 + t584;
t405 = m(3) * t499 + mrSges(3,1) * t552 - mrSges(3,3) * t545 + t541 * t553 - t544 * t587 + t407;
t404 = mrSges(5,1) * t449 + mrSges(4,2) * t481 - mrSges(4,3) * t452 - mrSges(5,3) * t451 + pkin(4) * t415 - qJ(4) * t414 + t605 * t509 + t613 * t510 + t597 * t535 + t607 * t552 + t596 * t553 + t572;
t403 = -mrSges(4,1) * t481 - mrSges(5,1) * t448 + mrSges(5,2) * t451 + mrSges(4,3) * t453 - pkin(3) * t414 - pkin(4) * t573 - pkin(9) * t583 - t561 * t408 - t565 * t409 - t611 * t509 - t605 * t510 + t597 * t536 + t606 * t552 + t595 * t553;
t402 = t565 * t408 - t561 * t409 + Ifges(3,6) * t546 + qJ(4) * t571 + Ifges(3,5) * t545 + mrSges(3,1) * t499 - mrSges(3,2) * t500 - mrSges(5,3) * t448 + mrSges(5,2) * t449 + mrSges(4,1) * t452 - mrSges(4,2) * t453 - pkin(9) * t415 - pkin(3) * t413 + pkin(2) * t407 - t596 * t536 + (-qJ(4) * t503 + t595) * t535 + t607 * t510 + (-mrSges(5,1) * qJ(4) - t606) * t509 + (t523 * t562 - t524 * t566) * t593 + (Ifges(3,3) - t612) * t552;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t585 - mrSges(2,2) * t580 + (mrSges(3,2) * t525 - mrSges(3,3) * t499 + Ifges(3,1) * t545 + Ifges(3,4) * t546 + Ifges(3,5) * t552 - qJ(3) * t407 - t557 * t403 + t604 * t404 + t522 * t586 - t553 * t523) * t600 + (-mrSges(3,1) * t525 + mrSges(3,3) * t500 + Ifges(3,4) * t545 + Ifges(3,2) * t546 + Ifges(3,6) * t552 - pkin(2) * t412 + qJ(3) * t584 + t604 * t403 + t557 * t404 - t522 * t587 + t553 * t524) * t599 + t559 * t402 + pkin(1) * ((t405 * t566 + t406 * t562) * t559 + (-m(3) * t525 + t546 * mrSges(3,1) - t545 * mrSges(3,2) + (-t540 * t562 + t541 * t566) * t593 - t412) * t558) + (-t405 * t562 + t406 * t566) * t609; t402; t412; t413; t572; t570;];
tauJ  = t1;
