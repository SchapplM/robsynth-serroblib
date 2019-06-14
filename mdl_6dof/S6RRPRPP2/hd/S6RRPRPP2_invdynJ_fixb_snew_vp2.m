% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:23:17
% EndTime: 2019-05-06 12:23:22
% DurationCPUTime: 2.64s
% Computational Cost: add. (19151->308), mult. (43237->363), div. (0->0), fcn. (29271->8), ass. (0->122)
t583 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t559 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t558 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t582 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t557 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t581 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t534 = sin(qJ(2));
t536 = cos(qJ(2));
t561 = qJD(1) * qJD(2);
t523 = qJDD(1) * t534 + t536 * t561;
t532 = sin(pkin(9));
t560 = qJDD(1) * t536;
t543 = -t534 * t561 + t560;
t571 = cos(pkin(9));
t498 = t571 * t523 + t532 * t543;
t514 = (t532 * t536 + t571 * t534) * qJD(1);
t533 = sin(qJ(4));
t575 = cos(qJ(4));
t500 = -t575 * qJD(2) + t514 * t533;
t466 = -t500 * qJD(4) + t533 * qJDD(2) + t575 * t498;
t563 = qJD(1) * t536;
t564 = qJD(1) * t534;
t513 = -t532 * t564 + t571 * t563;
t512 = qJD(4) - t513;
t568 = t500 * t512;
t580 = (-t466 + t568) * qJ(5);
t539 = qJD(1) ^ 2;
t535 = sin(qJ(1));
t537 = cos(qJ(1));
t546 = -g(1) * t537 - g(2) * t535;
t520 = -pkin(1) * t539 + qJDD(1) * pkin(7) + t546;
t567 = t534 * t520;
t573 = pkin(2) * t539;
t473 = qJDD(2) * pkin(2) - t523 * qJ(3) - t567 + (qJ(3) * t561 + t534 * t573 - g(3)) * t536;
t503 = -t534 * g(3) + t536 * t520;
t524 = qJD(2) * pkin(2) - qJ(3) * t564;
t531 = t536 ^ 2;
t474 = t543 * qJ(3) - qJD(2) * t524 - t531 * t573 + t503;
t441 = 0.2e1 * qJD(3) * t513 + t532 * t473 + t571 * t474;
t492 = -pkin(3) * t513 - pkin(8) * t514;
t538 = qJD(2) ^ 2;
t437 = -pkin(3) * t538 + qJDD(2) * pkin(8) + t492 * t513 + t441;
t551 = t535 * g(1) - t537 * g(2);
t545 = -qJDD(1) * pkin(1) - t551;
t480 = t524 * t564 - t543 * pkin(2) + qJDD(3) + (-qJ(3) * t531 - pkin(7)) * t539 + t545;
t497 = -t523 * t532 + t571 * t543;
t439 = (-qJD(2) * t513 - t498) * pkin(8) + (qJD(2) * t514 - t497) * pkin(3) + t480;
t432 = -t533 * t437 + t575 * t439;
t501 = t533 * qJD(2) + t575 * t514;
t475 = pkin(4) * t500 - qJ(5) * t501;
t496 = qJDD(4) - t497;
t511 = t512 ^ 2;
t430 = -t496 * pkin(4) - t511 * qJ(5) + t501 * t475 + qJDD(5) - t432;
t481 = -mrSges(6,2) * t500 + mrSges(6,3) * t512;
t579 = -m(6) * t430 + t496 * mrSges(6,1) + t512 * t481;
t482 = mrSges(7,2) * t512 + mrSges(7,3) * t500;
t577 = -0.2e1 * t501;
t423 = qJD(6) * t577 + (-t466 - t568) * qJ(6) + (t500 * t501 - t496) * pkin(5) + t430;
t477 = -mrSges(7,1) * t500 + mrSges(7,2) * t501;
t547 = -m(7) * t423 + t466 * mrSges(7,3) + t501 * t477;
t421 = -mrSges(7,1) * t496 - t482 * t512 - t547;
t476 = mrSges(6,1) * t500 - mrSges(6,3) * t501;
t418 = mrSges(6,2) * t466 + t476 * t501 + t421 - t579;
t433 = t575 * t437 + t533 * t439;
t576 = 2 * qJD(5);
t429 = -pkin(4) * t511 + t496 * qJ(5) - t500 * t475 + t512 * t576 + t433;
t465 = qJD(4) * t501 - t575 * qJDD(2) + t498 * t533;
t484 = -pkin(5) * t512 - qJ(6) * t501;
t499 = t500 ^ 2;
t425 = -pkin(5) * t499 + qJ(6) * t465 + 0.2e1 * qJD(6) * t500 + t484 * t512 + t429;
t485 = -mrSges(7,1) * t512 - mrSges(7,3) * t501;
t487 = -mrSges(6,1) * t512 + mrSges(6,2) * t501;
t555 = m(7) * t425 + t465 * mrSges(7,3) + t500 * t477;
t544 = m(6) * t429 + t496 * mrSges(6,3) + t512 * t487 + t555;
t552 = -t559 * t500 + t583 * t501 + t558 * t512;
t553 = t582 * t500 + t559 * t501 + t557 * t512;
t578 = -t557 * t465 + t558 * t466 + t581 * t496 + t552 * t500 + t553 * t501 + mrSges(5,1) * t432 - mrSges(6,1) * t430 - mrSges(7,1) * t423 - mrSges(5,2) * t433 + mrSges(7,2) * t425 + mrSges(6,3) * t429 - pkin(4) * t418 - pkin(5) * t421 + qJ(5) * (-mrSges(6,2) * t465 + mrSges(7,2) * t496 - t476 * t500 + t485 * t512 + t544);
t572 = -mrSges(5,3) - mrSges(6,2);
t570 = Ifges(3,6) * qJD(2);
t569 = qJD(2) * mrSges(3,1);
t491 = -mrSges(4,1) * t513 + mrSges(4,2) * t514;
t505 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t514;
t483 = -mrSges(5,2) * t512 - mrSges(5,3) * t500;
t565 = -mrSges(5,1) * t500 - mrSges(5,2) * t501 - t476;
t416 = m(5) * t432 + (t482 + t483) * t512 + t565 * t501 + (mrSges(5,1) + mrSges(7,1)) * t496 + t572 * t466 + t547 + t579;
t486 = mrSges(5,1) * t512 - mrSges(5,3) * t501;
t417 = m(5) * t433 + (t485 - t486) * t512 + t565 * t500 + (-mrSges(5,2) + mrSges(7,2)) * t496 + t572 * t465 + t544;
t549 = -t416 * t533 + t575 * t417;
t408 = m(4) * t441 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t497 - qJD(2) * t505 + t491 * t513 + t549;
t562 = qJD(3) * t514;
t508 = -0.2e1 * t562;
t566 = t571 * t473 - t532 * t474;
t440 = t508 + t566;
t504 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t513;
t542 = qJDD(2) * pkin(3) + t538 * pkin(8) - t514 * t492 + t566;
t427 = -qJ(6) * t499 + qJDD(6) + t508 + (-pkin(4) - pkin(5)) * t465 - t580 + (-pkin(4) * t512 + t484 + t576) * t501 + t542;
t422 = m(7) * t427 - t465 * mrSges(7,1) + t466 * mrSges(7,2) - t500 * t482 + t501 * t485;
t436 = 0.2e1 * t562 - t542;
t431 = qJD(5) * t577 + t580 + (t501 * t512 + t465) * pkin(4) + t436;
t420 = m(6) * t431 + mrSges(6,1) * t465 - t466 * mrSges(6,3) + t481 * t500 - t501 * t487 - t422;
t540 = -m(5) * t436 - t465 * mrSges(5,1) - mrSges(5,2) * t466 - t500 * t483 - t486 * t501 - t420;
t412 = m(4) * t440 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t498 + qJD(2) * t504 - t491 * t514 + t540;
t404 = t532 * t408 + t571 * t412;
t410 = t575 * t416 + t533 * t417;
t554 = t557 * t500 - t558 * t501 - t581 * t512;
t550 = t571 * t408 - t532 * t412;
t409 = m(4) * t480 - t497 * mrSges(4,1) + t498 * mrSges(4,2) - t513 * t504 + t514 * t505 + t410;
t526 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t563;
t525 = -mrSges(3,3) * t564 + t569;
t522 = (-mrSges(3,1) * t536 + mrSges(3,2) * t534) * qJD(1);
t519 = -t539 * pkin(7) + t545;
t517 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t534 + Ifges(3,4) * t536) * qJD(1);
t516 = t570 + (Ifges(3,4) * t534 + Ifges(3,2) * t536) * qJD(1);
t502 = -t536 * g(3) - t567;
t490 = Ifges(4,1) * t514 + Ifges(4,4) * t513 + Ifges(4,5) * qJD(2);
t489 = Ifges(4,4) * t514 + Ifges(4,2) * t513 + Ifges(4,6) * qJD(2);
t488 = Ifges(4,5) * t514 + Ifges(4,6) * t513 + Ifges(4,3) * qJD(2);
t405 = mrSges(5,2) * t436 + mrSges(6,2) * t430 + mrSges(7,2) * t427 - mrSges(5,3) * t432 - mrSges(6,3) * t431 - mrSges(7,3) * t423 - qJ(5) * t420 - qJ(6) * t421 - t559 * t465 + t583 * t466 + t558 * t496 + t554 * t500 - t553 * t512;
t403 = -mrSges(5,1) * t436 + mrSges(5,3) * t433 - mrSges(6,1) * t431 + mrSges(6,2) * t429 + mrSges(7,1) * t427 - mrSges(7,3) * t425 + pkin(5) * t422 - qJ(6) * t555 - pkin(4) * t420 + (-qJ(6) * t485 + t552) * t512 + t554 * t501 + (-mrSges(7,2) * qJ(6) + t557) * t496 + t559 * t466 + t582 * t465;
t402 = -mrSges(4,1) * t480 + mrSges(4,3) * t441 + Ifges(4,4) * t498 + Ifges(4,2) * t497 + Ifges(4,6) * qJDD(2) - pkin(3) * t410 + qJD(2) * t490 - t514 * t488 - t578;
t401 = mrSges(4,2) * t480 - mrSges(4,3) * t440 + Ifges(4,1) * t498 + Ifges(4,4) * t497 + Ifges(4,5) * qJDD(2) - pkin(8) * t410 - qJD(2) * t489 - t533 * t403 + t575 * t405 + t513 * t488;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t551 - mrSges(2,2) * t546 + t534 * (mrSges(3,2) * t519 - mrSges(3,3) * t502 + Ifges(3,1) * t523 + Ifges(3,4) * t543 + Ifges(3,5) * qJDD(2) - qJ(3) * t404 - qJD(2) * t516 + t571 * t401 - t532 * t402) + t536 * (-mrSges(3,1) * t519 + mrSges(3,3) * t503 + Ifges(3,4) * t523 + Ifges(3,2) * t543 + Ifges(3,6) * qJDD(2) - pkin(2) * t409 + qJ(3) * t550 + qJD(2) * t517 + t532 * t401 + t571 * t402) + pkin(1) * (mrSges(3,1) * t560 - m(3) * t519 - t523 * mrSges(3,2) + (t536 * t526 + (-t525 - t569) * t534) * qJD(1) - t409) + pkin(7) * (t536 * (m(3) * t503 - qJDD(2) * mrSges(3,2) + t543 * mrSges(3,3) - qJD(2) * t525 + t522 * t563 + t550) - t534 * (m(3) * t502 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t523 + qJD(2) * t526 - t522 * t564 + t404)); Ifges(3,5) * t523 + Ifges(3,6) * t560 + mrSges(3,1) * t502 - mrSges(3,2) * t503 + Ifges(4,5) * t498 + Ifges(4,6) * t497 + t514 * t489 - t513 * t490 + mrSges(4,1) * t440 - mrSges(4,2) * t441 + t533 * t405 + t575 * t403 + pkin(3) * t540 + pkin(8) * t549 + pkin(2) * t404 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (-t536 * t517 + (t516 - t570) * t534) * qJD(1); t409; t578; t418; t422;];
tauJ  = t1;
