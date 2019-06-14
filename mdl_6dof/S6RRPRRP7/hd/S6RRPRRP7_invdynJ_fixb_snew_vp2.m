% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP7
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
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:12:55
% EndTime: 2019-05-06 18:13:00
% DurationCPUTime: 2.73s
% Computational Cost: add. (19197->310), mult. (39064->367), div. (0->0), fcn. (24295->8), ass. (0->124)
t583 = Ifges(3,1) + Ifges(4,1);
t582 = Ifges(6,1) + Ifges(7,1);
t572 = Ifges(3,4) - Ifges(4,5);
t571 = Ifges(6,4) - Ifges(7,5);
t570 = Ifges(3,5) + Ifges(4,4);
t569 = -Ifges(6,5) - Ifges(7,4);
t581 = Ifges(3,2) + Ifges(4,3);
t580 = Ifges(6,2) + Ifges(7,3);
t568 = Ifges(3,6) - Ifges(4,6);
t567 = Ifges(6,6) - Ifges(7,6);
t579 = Ifges(3,3) + Ifges(4,2);
t578 = -Ifges(6,3) - Ifges(7,2);
t532 = sin(qJ(4));
t533 = sin(qJ(2));
t535 = cos(qJ(4));
t536 = cos(qJ(2));
t494 = (t532 * t533 + t535 * t536) * qJD(1);
t554 = qJD(1) * qJD(2);
t552 = t536 * t554;
t504 = qJDD(1) * t533 + t552;
t551 = t533 * t554;
t505 = qJDD(1) * t536 - t551;
t464 = -qJD(4) * t494 + t504 * t535 - t505 * t532;
t495 = (-t532 * t536 + t533 * t535) * qJD(1);
t524 = -qJD(2) + qJD(4);
t531 = sin(qJ(5));
t575 = cos(qJ(5));
t475 = t495 * t531 - t524 * t575;
t523 = -qJDD(2) + qJDD(4);
t435 = -qJD(5) * t475 + t464 * t575 + t523 * t531;
t476 = t495 * t575 + t524 * t531;
t453 = mrSges(7,1) * t475 - mrSges(7,3) * t476;
t556 = qJD(1) * t533;
t512 = -qJD(2) * pkin(3) - pkin(8) * t556;
t539 = qJD(1) ^ 2;
t534 = sin(qJ(1));
t537 = cos(qJ(1));
t557 = g(1) * t534 - g(2) * t537;
t496 = -qJDD(1) * pkin(1) - pkin(7) * t539 - t557;
t546 = -t505 * pkin(2) + t496 + (-t504 - t552) * qJ(3);
t565 = t536 ^ 2 * t539;
t576 = 2 * qJD(3);
t431 = -pkin(2) * t551 + t505 * pkin(3) - pkin(8) * t565 - t546 + (t512 + t576) * t556;
t463 = -qJD(4) * t495 - t504 * t532 - t505 * t535;
t425 = t431 + (t495 * t524 - t463) * pkin(4) + (t494 * t524 - t464) * pkin(9);
t548 = -g(1) * t537 - g(2) * t534;
t497 = -pkin(1) * t539 + qJDD(1) * pkin(7) + t548;
t478 = -g(3) * t533 + t497 * t536;
t501 = (-pkin(2) * t536 - qJ(3) * t533) * qJD(1);
t538 = qJD(2) ^ 2;
t555 = qJD(1) * t536;
t459 = -pkin(2) * t538 + qJDD(2) * qJ(3) + qJD(2) * t576 + t501 * t555 + t478;
t439 = -pkin(3) * t565 - pkin(8) * t505 + qJD(2) * t512 + t459;
t477 = -g(3) * t536 - t497 * t533;
t462 = -qJDD(2) * pkin(2) - t538 * qJ(3) + t501 * t556 + qJDD(3) - t477;
t440 = (-t504 + t552) * pkin(8) + (-t533 * t536 * t539 - qJDD(2)) * pkin(3) + t462;
t430 = t439 * t535 + t440 * t532;
t474 = pkin(4) * t494 - pkin(9) * t495;
t522 = t524 ^ 2;
t428 = -pkin(4) * t522 + pkin(9) * t523 - t474 * t494 + t430;
t422 = t425 * t575 - t428 * t531;
t452 = pkin(5) * t475 - qJ(6) * t476;
t461 = qJDD(5) - t463;
t487 = qJD(5) + t494;
t483 = t487 ^ 2;
t420 = -pkin(5) * t461 - qJ(6) * t483 + t452 * t476 + qJDD(6) - t422;
t465 = -mrSges(7,2) * t475 + mrSges(7,3) * t487;
t549 = -m(7) * t420 + mrSges(7,1) * t461 + t465 * t487;
t416 = t435 * mrSges(7,2) + t476 * t453 - t549;
t423 = t425 * t531 + t428 * t575;
t419 = -pkin(5) * t483 + qJ(6) * t461 + 0.2e1 * qJD(6) * t487 - t452 * t475 + t423;
t434 = qJD(5) * t476 + t464 * t531 - t523 * t575;
t468 = -mrSges(7,1) * t487 + mrSges(7,2) * t476;
t553 = m(7) * t419 + mrSges(7,3) * t461 + t468 * t487;
t562 = t475 * t571 - t476 * t582 + t487 * t569;
t563 = t475 * t580 - t476 * t571 - t487 * t567;
t577 = -t434 * t567 - t435 * t569 - t578 * t461 - t475 * t562 - t476 * t563 + mrSges(6,1) * t422 - mrSges(7,1) * t420 - mrSges(6,2) * t423 + mrSges(7,3) * t419 - pkin(5) * t416 + qJ(6) * (-t434 * mrSges(7,2) - t475 * t453 + t553);
t574 = mrSges(3,3) + mrSges(4,2);
t573 = -mrSges(6,3) - mrSges(7,2);
t467 = mrSges(6,1) * t487 - mrSges(6,3) * t476;
t561 = -mrSges(6,1) * t475 - mrSges(6,2) * t476 - t453;
t411 = m(6) * t423 - t461 * mrSges(6,2) + t434 * t573 - t487 * t467 + t475 * t561 + t553;
t466 = -mrSges(6,2) * t487 - mrSges(6,3) * t475;
t413 = m(6) * t422 + t461 * mrSges(6,1) + t435 * t573 + t487 * t466 + t476 * t561 + t549;
t406 = t411 * t531 + t413 * t575;
t564 = t475 * t567 + t476 * t569 + t487 * t578;
t560 = t579 * qJD(2) + (t533 * t570 + t536 * t568) * qJD(1);
t559 = -t568 * qJD(2) + (-t533 * t572 - t581 * t536) * qJD(1);
t558 = t570 * qJD(2) + (t533 * t583 + t536 * t572) * qJD(1);
t407 = t411 * t575 - t413 * t531;
t473 = mrSges(5,1) * t494 + mrSges(5,2) * t495;
t480 = mrSges(5,1) * t524 - mrSges(5,3) * t495;
t403 = m(5) * t430 - mrSges(5,2) * t523 + mrSges(5,3) * t463 - t473 * t494 - t480 * t524 + t407;
t429 = -t439 * t532 + t535 * t440;
t427 = -t523 * pkin(4) - t522 * pkin(9) + t474 * t495 - t429;
t421 = -0.2e1 * qJD(6) * t476 + (t475 * t487 - t435) * qJ(6) + (t476 * t487 + t434) * pkin(5) + t427;
t417 = m(7) * t421 + mrSges(7,1) * t434 - mrSges(7,3) * t435 + t465 * t475 - t468 * t476;
t414 = -m(6) * t427 - mrSges(6,1) * t434 - mrSges(6,2) * t435 - t466 * t475 - t467 * t476 - t417;
t479 = -mrSges(5,2) * t524 - mrSges(5,3) * t494;
t408 = m(5) * t429 + mrSges(5,1) * t523 - mrSges(5,3) * t464 - t473 * t495 + t479 * t524 + t414;
t550 = t403 * t535 - t408 * t532;
t400 = t532 * t403 + t535 * t408;
t502 = (-mrSges(4,1) * t536 - mrSges(4,3) * t533) * qJD(1);
t509 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t556;
t545 = m(4) * t459 + qJDD(2) * mrSges(4,3) + qJD(2) * t509 + t502 * t555 + t550;
t511 = mrSges(4,2) * t555 + qJD(2) * mrSges(4,3);
t544 = m(4) * t462 - qJDD(2) * mrSges(4,1) - qJD(2) * t511 + t400;
t543 = m(5) * t431 - mrSges(5,1) * t463 + t464 * mrSges(5,2) + t479 * t494 + t495 * t480 + t406;
t451 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t556 + t546;
t542 = m(4) * t451 - t543;
t404 = -mrSges(6,1) * t427 - mrSges(7,1) * t421 + mrSges(7,2) * t419 + mrSges(6,3) * t423 - pkin(5) * t417 - t434 * t580 + t435 * t571 + t461 * t567 + t476 * t564 - t487 * t562;
t405 = mrSges(6,2) * t427 + mrSges(7,2) * t420 - mrSges(6,3) * t422 - mrSges(7,3) * t421 - qJ(6) * t417 - t434 * t571 + t435 * t582 - t461 * t569 + t475 * t564 + t487 * t563;
t470 = Ifges(5,4) * t495 - Ifges(5,2) * t494 + Ifges(5,6) * t524;
t471 = Ifges(5,1) * t495 - Ifges(5,4) * t494 + Ifges(5,5) * t524;
t540 = mrSges(5,1) * t429 - mrSges(5,2) * t430 + Ifges(5,5) * t464 + Ifges(5,6) * t463 + Ifges(5,3) * t523 + pkin(4) * t414 + pkin(9) * t407 + t404 * t575 + t531 * t405 + t495 * t470 + t494 * t471;
t510 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t555;
t508 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t556;
t503 = (-mrSges(3,1) * t536 + mrSges(3,2) * t533) * qJD(1);
t469 = Ifges(5,5) * t495 - Ifges(5,6) * t494 + Ifges(5,3) * t524;
t401 = -t505 * mrSges(4,1) - t504 * mrSges(4,3) + (-t509 * t533 - t511 * t536) * qJD(1) + t542;
t399 = t504 * mrSges(4,2) + t502 * t556 + t544;
t398 = -mrSges(5,1) * t431 + mrSges(5,3) * t430 + Ifges(5,4) * t464 + Ifges(5,2) * t463 + Ifges(5,6) * t523 - pkin(4) * t406 - t495 * t469 + t524 * t471 - t577;
t397 = mrSges(5,2) * t431 - mrSges(5,3) * t429 + Ifges(5,1) * t464 + Ifges(5,4) * t463 + Ifges(5,5) * t523 - pkin(9) * t406 - t404 * t531 + t405 * t575 - t469 * t494 - t470 * t524;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t548 + t533 * (mrSges(3,2) * t496 + mrSges(4,2) * t462 - mrSges(3,3) * t477 - mrSges(4,3) * t451 - pkin(8) * t400 - qJ(3) * t401 + t559 * qJD(2) + t570 * qJDD(2) + t535 * t397 - t532 * t398 + t583 * t504 + t572 * t505 + t560 * t555) + t536 * (-mrSges(3,1) * t496 - mrSges(4,1) * t451 + mrSges(4,2) * t459 + mrSges(3,3) * t478 - pkin(2) * t401 + pkin(3) * t543 - pkin(8) * t550 + t558 * qJD(2) + t568 * qJDD(2) - t532 * t397 - t535 * t398 + t572 * t504 + t581 * t505 - t560 * t556) + pkin(1) * (-m(3) * t496 + (mrSges(3,1) + mrSges(4,1)) * t505 + (-mrSges(3,2) + mrSges(4,3)) * t504 + ((t510 + t511) * t536 + (-t508 + t509) * t533) * qJD(1) - t542) + pkin(7) * (t536 * (m(3) * t478 - qJDD(2) * mrSges(3,2) - qJD(2) * t508 + t503 * t555 + t505 * t574 + t545) + (-m(3) * t477 - qJDD(2) * mrSges(3,1) - qJD(2) * t510 + t574 * t504 + (t502 + t503) * t556 + t544) * t533); -t540 + t570 * t504 + mrSges(3,1) * t477 - mrSges(3,2) * t478 - mrSges(4,1) * t462 + mrSges(4,3) * t459 - pkin(3) * t400 - pkin(2) * t399 + qJ(3) * t545 + (mrSges(4,2) * qJ(3) + t568) * t505 + (-t559 * t533 - t558 * t536) * qJD(1) + t579 * qJDD(2); t399; t540; t577; t416;];
tauJ  = t1;
