% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:53:17
% EndTime: 2019-05-06 18:53:26
% DurationCPUTime: 3.88s
% Computational Cost: add. (23534->320), mult. (48369->367), div. (0->0), fcn. (29765->8), ass. (0->128)
t588 = Ifges(6,4) + Ifges(7,4);
t604 = Ifges(6,2) + Ifges(7,2);
t598 = Ifges(6,6) + Ifges(7,6);
t603 = -2 * qJD(3);
t602 = Ifges(3,1) + Ifges(4,2);
t601 = Ifges(6,1) + Ifges(7,1);
t589 = Ifges(3,4) + Ifges(4,6);
t587 = Ifges(3,5) - Ifges(4,4);
t600 = Ifges(6,5) + Ifges(7,5);
t599 = Ifges(3,2) + Ifges(4,3);
t586 = Ifges(3,6) - Ifges(4,5);
t597 = Ifges(3,3) + Ifges(4,1);
t596 = Ifges(6,3) + Ifges(7,3);
t544 = sin(qJ(4));
t548 = cos(qJ(4));
t549 = cos(qJ(2));
t576 = qJD(1) * t549;
t517 = -qJD(2) * t544 - t548 * t576;
t518 = qJD(2) * t548 - t544 * t576;
t543 = sin(qJ(5));
t547 = cos(qJ(5));
t488 = t517 * t547 - t518 * t543;
t489 = t517 * t543 + t518 * t547;
t545 = sin(qJ(2));
t536 = t545 * qJD(1);
t533 = t536 + qJD(4);
t531 = qJD(5) + t533;
t595 = t604 * t488 + t588 * t489 + t598 * t531;
t575 = qJD(1) * qJD(2);
t568 = t545 * t575;
t523 = qJDD(1) * t549 - t568;
t485 = -qJD(4) * t518 - qJDD(2) * t544 - t523 * t548;
t486 = qJD(4) * t517 + qJDD(2) * t548 - t523 * t544;
t447 = -qJD(5) * t489 + t485 * t547 - t486 * t543;
t472 = -mrSges(7,2) * t531 + mrSges(7,3) * t488;
t473 = -mrSges(6,2) * t531 + mrSges(6,3) * t488;
t594 = -(t472 + t473) * t488 - (mrSges(6,1) + mrSges(7,1)) * t447;
t552 = qJD(1) ^ 2;
t546 = sin(qJ(1));
t550 = cos(qJ(1));
t563 = -g(1) * t550 - g(2) * t546;
t511 = -pkin(1) * t552 + qJDD(1) * pkin(7) + t563;
t494 = -g(3) * t545 + t549 * t511;
t519 = (-pkin(2) * t549 - qJ(3) * t545) * qJD(1);
t551 = qJD(2) ^ 2;
t470 = pkin(2) * t551 - qJDD(2) * qJ(3) + qJD(2) * t603 - t519 * t576 - t494;
t593 = -t447 * mrSges(7,1) - t488 * t472;
t592 = pkin(7) * t552;
t591 = mrSges(3,1) - mrSges(4,2);
t530 = pkin(3) * t536 - qJD(2) * pkin(8);
t542 = t549 ^ 2;
t569 = t549 * t575;
t522 = qJDD(1) * t545 + t569;
t567 = g(1) * t546 - t550 * g(2);
t561 = -qJDD(1) * pkin(1) - t567;
t556 = pkin(2) * t568 + t536 * t603 + (-t522 - t569) * qJ(3) + t561;
t451 = -t530 * t536 + (-pkin(3) * t542 - pkin(7)) * t552 + (-pkin(2) - pkin(8)) * t523 + t556;
t493 = -t549 * g(3) - t545 * t511;
t471 = -qJDD(2) * pkin(2) - qJ(3) * t551 + t519 * t536 + qJDD(3) - t493;
t456 = (-t545 * t549 * t552 - qJDD(2)) * pkin(8) + (t522 - t569) * pkin(3) + t471;
t434 = -t451 * t544 + t548 * t456;
t516 = qJDD(4) + t522;
t430 = (t517 * t533 - t486) * pkin(9) + (t517 * t518 + t516) * pkin(4) + t434;
t435 = t548 * t451 + t544 * t456;
t495 = pkin(4) * t533 - pkin(9) * t518;
t515 = t517 ^ 2;
t432 = -pkin(4) * t515 + pkin(9) * t485 - t495 * t533 + t435;
t424 = t547 * t430 - t432 * t543;
t448 = qJD(5) * t488 + t485 * t543 + t486 * t547;
t465 = -mrSges(7,1) * t488 + mrSges(7,2) * t489;
t466 = -mrSges(6,1) * t488 + mrSges(6,2) * t489;
t513 = qJDD(5) + t516;
t419 = -0.2e1 * qJD(6) * t489 + (t488 * t531 - t448) * qJ(6) + (t488 * t489 + t513) * pkin(5) + t424;
t573 = m(7) * t419 + t513 * mrSges(7,1) + t531 * t472;
t411 = m(6) * t424 + mrSges(6,1) * t513 + t473 * t531 + (-t465 - t466) * t489 + (-mrSges(6,3) - mrSges(7,3)) * t448 + t573;
t425 = t543 * t430 + t547 * t432;
t475 = mrSges(7,1) * t531 - mrSges(7,3) * t489;
t476 = mrSges(6,1) * t531 - mrSges(6,3) * t489;
t474 = pkin(5) * t531 - qJ(6) * t489;
t487 = t488 ^ 2;
t421 = -pkin(5) * t487 + qJ(6) * t447 + 0.2e1 * qJD(6) * t488 - t474 * t531 + t425;
t572 = m(7) * t421 + t447 * mrSges(7,3) + t488 * t465;
t414 = m(6) * t425 + mrSges(6,3) * t447 + t466 * t488 + (-t475 - t476) * t531 + (-mrSges(6,2) - mrSges(7,2)) * t513 + t572;
t409 = t547 * t411 + t543 * t414;
t583 = -t598 * t488 - t600 * t489 - t596 * t531;
t582 = -t588 * t488 - t601 * t489 - t600 * t531;
t580 = t597 * qJD(2) + (t545 * t587 + t549 * t586) * qJD(1);
t579 = t586 * qJD(2) + (t545 * t589 + t549 * t599) * qJD(1);
t578 = t587 * qJD(2) + (t545 * t602 + t549 * t589) * qJD(1);
t528 = -mrSges(4,1) * t576 - qJD(2) * mrSges(4,3);
t577 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t576 - t528;
t455 = -pkin(8) * t542 * t552 + pkin(3) * t523 + qJD(2) * t530 - t470;
t437 = -pkin(4) * t485 - pkin(9) * t515 + t518 * t495 + t455;
t427 = -pkin(5) * t447 - qJ(6) * t487 + t474 * t489 + qJDD(6) + t437;
t571 = m(7) * t427 + t448 * mrSges(7,2) + t489 * t475;
t490 = -mrSges(5,1) * t517 + mrSges(5,2) * t518;
t491 = -mrSges(5,2) * t533 + mrSges(5,3) * t517;
t405 = m(5) * t434 + mrSges(5,1) * t516 - mrSges(5,3) * t486 - t490 * t518 + t491 * t533 + t409;
t492 = mrSges(5,1) * t533 - mrSges(5,3) * t518;
t564 = -t411 * t543 + t547 * t414;
t406 = m(5) * t435 - mrSges(5,2) * t516 + mrSges(5,3) * t485 + t490 * t517 - t492 * t533 + t564;
t565 = -t405 * t544 + t548 * t406;
t402 = t548 * t405 + t544 * t406;
t467 = -pkin(2) * t523 + t556 - t592;
t562 = m(4) * t467 + t565;
t560 = m(6) * t437 + t448 * mrSges(6,2) + t489 * t476 + t571;
t558 = m(4) * t471 + t522 * mrSges(4,1) + t402;
t557 = -m(5) * t455 + t485 * mrSges(5,1) - t486 * mrSges(5,2) + t517 * t491 - t518 * t492 - t560;
t520 = (mrSges(4,2) * t549 - mrSges(4,3) * t545) * qJD(1);
t529 = mrSges(4,1) * t536 + qJD(2) * mrSges(4,2);
t555 = -m(4) * t470 + qJDD(2) * mrSges(4,3) + qJD(2) * t529 + t520 * t576 - t557;
t416 = -mrSges(7,3) * t448 - t465 * t489 + t573;
t554 = mrSges(6,1) * t424 + mrSges(7,1) * t419 - mrSges(6,2) * t425 - mrSges(7,2) * t421 + pkin(5) * t416 + t598 * t447 + t600 * t448 + t582 * t488 + t595 * t489 + t596 * t513;
t478 = Ifges(5,4) * t518 + Ifges(5,2) * t517 + Ifges(5,6) * t533;
t479 = Ifges(5,1) * t518 + Ifges(5,4) * t517 + Ifges(5,5) * t533;
t553 = mrSges(5,1) * t434 - mrSges(5,2) * t435 + Ifges(5,5) * t486 + Ifges(5,6) * t485 + Ifges(5,3) * t516 + pkin(4) * t409 + t518 * t478 - t517 * t479 + t554;
t526 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t536;
t521 = (-mrSges(3,1) * t549 + mrSges(3,2) * t545) * qJD(1);
t510 = t561 - t592;
t477 = Ifges(5,5) * t518 + Ifges(5,6) * t517 + Ifges(5,3) * t533;
t422 = t571 + t593;
t407 = mrSges(6,2) * t437 + mrSges(7,2) * t427 - mrSges(6,3) * t424 - mrSges(7,3) * t419 - qJ(6) * t416 + t588 * t447 + t601 * t448 - t583 * t488 + t600 * t513 - t595 * t531;
t403 = -mrSges(6,1) * t437 + mrSges(6,3) * t425 - mrSges(7,1) * t427 + mrSges(7,3) * t421 - pkin(5) * t422 + qJ(6) * t572 + (-qJ(6) * t475 - t582) * t531 + (-mrSges(7,2) * qJ(6) + t598) * t513 + t583 * t489 + t588 * t448 + t604 * t447;
t401 = qJDD(2) * mrSges(4,2) + qJD(2) * t528 + t520 * t536 + t558;
t400 = mrSges(4,2) * t523 - mrSges(4,3) * t522 + (t528 * t549 - t529 * t545) * qJD(1) + t562;
t399 = mrSges(5,2) * t455 - mrSges(5,3) * t434 + Ifges(5,1) * t486 + Ifges(5,4) * t485 + Ifges(5,5) * t516 - pkin(9) * t409 - t403 * t543 + t407 * t547 + t477 * t517 - t478 * t533;
t398 = Ifges(5,4) * t486 + Ifges(5,2) * t485 + Ifges(5,6) * t516 - t518 * t477 + t533 * t479 - mrSges(5,1) * t455 + mrSges(5,3) * t435 + t543 * t407 + t547 * t403 - pkin(4) * (t560 + t594) + pkin(9) * t564;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t567 - mrSges(2,2) * t563 + t545 * (mrSges(4,1) * t471 + mrSges(3,2) * t510 - mrSges(3,3) * t493 - mrSges(4,3) * t467 + pkin(3) * t402 - qJ(3) * t400 - t579 * qJD(2) + t587 * qJDD(2) + t602 * t522 + t589 * t523 + t580 * t576 + t553) + t549 * (-mrSges(3,1) * t510 + mrSges(3,3) * t494 - mrSges(4,1) * t470 + mrSges(4,2) * t467 - t544 * t399 - t548 * t398 - pkin(3) * (t557 - t594) - pkin(8) * t565 - pkin(2) * t400 + t599 * t523 + t589 * t522 + t586 * qJDD(2) + t578 * qJD(2) - t580 * t536) + pkin(1) * (-m(3) * t510 + t591 * t523 + (-mrSges(3,2) + mrSges(4,3)) * t522 + (t577 * t549 + (-t526 + t529) * t545) * qJD(1) - t562) + pkin(7) * (t549 * (t555 + t521 * t576 - qJDD(2) * mrSges(3,2) + (mrSges(4,1) + mrSges(3,3)) * t523 - qJD(2) * t526 + m(3) * t494 + t594) + (-m(3) * t493 + t522 * mrSges(3,3) - t591 * qJDD(2) - t577 * qJD(2) + (t520 + t521) * t536 + t558) * t545); mrSges(3,1) * t493 - mrSges(3,2) * t494 + mrSges(4,2) * t471 - mrSges(4,3) * t470 + t548 * t399 - t544 * t398 - pkin(8) * t402 - pkin(2) * t401 + qJ(3) * (-t447 * mrSges(6,1) - t488 * t473 + t555 + t593) + (mrSges(4,1) * qJ(3) + t586) * t523 + t587 * t522 + t597 * qJDD(2) + (t579 * t545 - t578 * t549) * qJD(1); t401; t553; t554; t422;];
tauJ  = t1;
