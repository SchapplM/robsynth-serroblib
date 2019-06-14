% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:48:15
% EndTime: 2019-05-05 16:48:20
% DurationCPUTime: 3.52s
% Computational Cost: add. (29366->301), mult. (72741->370), div. (0->0), fcn. (53448->10), ass. (0->130)
t590 = -2 * qJD(4);
t589 = Ifges(5,1) + Ifges(6,1);
t581 = Ifges(5,4) - Ifges(6,5);
t580 = -Ifges(5,5) - Ifges(6,4);
t588 = Ifges(5,2) + Ifges(6,3);
t587 = -Ifges(6,2) - Ifges(5,3);
t579 = Ifges(5,6) - Ifges(6,6);
t545 = qJD(1) ^ 2;
t541 = sin(qJ(1));
t543 = cos(qJ(1));
t557 = -g(1) * t543 - g(2) * t541;
t524 = -pkin(1) * t545 + qJDD(1) * qJ(2) + t557;
t537 = sin(pkin(9));
t538 = cos(pkin(9));
t565 = qJD(1) * qJD(2);
t562 = -g(3) * t538 - 0.2e1 * t537 * t565;
t577 = pkin(7) * qJDD(1);
t583 = pkin(2) * t545;
t484 = (t538 * t583 - t524 - t577) * t537 + t562;
t509 = -g(3) * t537 + (t524 + 0.2e1 * t565) * t538;
t535 = t538 ^ 2;
t494 = -t535 * t583 + t538 * t577 + t509;
t540 = sin(qJ(3));
t584 = cos(qJ(3));
t459 = t540 * t484 + t584 * t494;
t564 = t538 * t584;
t569 = qJD(1) * t537;
t522 = -qJD(1) * t564 + t540 * t569;
t551 = t537 * t584 + t538 * t540;
t523 = t551 * qJD(1);
t500 = pkin(3) * t522 - qJ(4) * t523;
t544 = qJD(3) ^ 2;
t449 = -pkin(3) * t544 + qJDD(3) * qJ(4) - t500 * t522 + t459;
t563 = g(1) * t541 - t543 * g(2);
t556 = qJDD(2) - t563;
t570 = -t537 ^ 2 - t535;
t505 = (-pkin(2) * t538 - pkin(1)) * qJDD(1) + (pkin(7) * t570 - qJ(2)) * t545 + t556;
t568 = qJD(3) * t523;
t506 = t568 + (t537 * t540 - t564) * qJDD(1);
t566 = t522 * qJD(3);
t507 = qJDD(1) * t551 - t566;
t451 = (-t507 + t566) * qJ(4) + (t506 + t568) * pkin(3) + t505;
t536 = sin(pkin(10));
t578 = cos(pkin(10));
t514 = t536 * qJD(3) + t523 * t578;
t437 = -t536 * t449 + t451 * t578 + t514 * t590;
t493 = t536 * qJDD(3) + t507 * t578;
t458 = t584 * t484 - t540 * t494;
t548 = qJDD(3) * pkin(3) + qJ(4) * t544 - t523 * t500 - qJDD(4) + t458;
t513 = -qJD(3) * t578 + t523 * t536;
t576 = t513 * t522;
t586 = (-t493 + t576) * qJ(5) - t548;
t585 = 2 * qJD(5);
t582 = -mrSges(5,3) - mrSges(6,2);
t501 = mrSges(4,1) * t522 + mrSges(4,2) * t523;
t516 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t523;
t438 = t578 * t449 + t536 * t451 + t513 * t590;
t489 = mrSges(5,1) * t522 - mrSges(5,3) * t514;
t492 = -qJDD(3) * t578 + t507 * t536;
t476 = pkin(4) * t513 - qJ(5) * t514;
t521 = t522 ^ 2;
t434 = -pkin(4) * t521 + t506 * qJ(5) - t513 * t476 + t522 * t585 + t438;
t490 = -mrSges(6,1) * t522 + mrSges(6,2) * t514;
t435 = -t506 * pkin(4) - t521 * qJ(5) + t514 * t476 + qJDD(5) - t437;
t429 = (-t493 - t576) * pkin(8) + (t513 * t514 - t506) * pkin(5) + t435;
t491 = -pkin(5) * t522 - pkin(8) * t514;
t512 = t513 ^ 2;
t430 = -pkin(5) * t512 + pkin(8) * t492 + t491 * t522 + t434;
t539 = sin(qJ(6));
t542 = cos(qJ(6));
t427 = t429 * t542 - t430 * t539;
t474 = t513 * t542 - t514 * t539;
t448 = qJD(6) * t474 + t492 * t539 + t493 * t542;
t475 = t513 * t539 + t514 * t542;
t457 = -mrSges(7,1) * t474 + mrSges(7,2) * t475;
t519 = qJD(6) - t522;
t462 = -mrSges(7,2) * t519 + mrSges(7,3) * t474;
t504 = qJDD(6) - t506;
t423 = m(7) * t427 + mrSges(7,1) * t504 - mrSges(7,3) * t448 - t457 * t475 + t462 * t519;
t428 = t429 * t539 + t430 * t542;
t447 = -qJD(6) * t475 + t492 * t542 - t493 * t539;
t463 = mrSges(7,1) * t519 - mrSges(7,3) * t475;
t424 = m(7) * t428 - mrSges(7,2) * t504 + mrSges(7,3) * t447 + t457 * t474 - t463 * t519;
t559 = -t423 * t539 + t542 * t424;
t550 = m(6) * t434 + t506 * mrSges(6,3) + t522 * t490 + t559;
t477 = mrSges(6,1) * t513 - mrSges(6,3) * t514;
t571 = -mrSges(5,1) * t513 - mrSges(5,2) * t514 - t477;
t412 = m(5) * t438 - mrSges(5,2) * t506 - t489 * t522 + t492 * t582 + t513 * t571 + t550;
t488 = -mrSges(5,2) * t522 - mrSges(5,3) * t513;
t416 = t423 * t542 + t424 * t539;
t487 = -mrSges(6,2) * t513 + mrSges(6,3) * t522;
t549 = -m(6) * t435 + t506 * mrSges(6,1) + t522 * t487 - t416;
t414 = m(5) * t437 + mrSges(5,1) * t506 + t488 * t522 + t493 * t582 + t514 * t571 + t549;
t560 = t578 * t412 - t414 * t536;
t408 = m(4) * t459 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t506 - qJD(3) * t516 - t501 * t522 + t560;
t436 = -0.2e1 * qJD(5) * t514 + (t514 * t522 + t492) * pkin(4) + t586;
t432 = -pkin(8) * t512 + (-pkin(4) - pkin(5)) * t492 + (-pkin(4) * t522 + t491 + t585) * t514 - t586;
t553 = -m(7) * t432 + t447 * mrSges(7,1) - t448 * mrSges(7,2) + t474 * t462 - t475 * t463;
t425 = m(6) * t436 + mrSges(6,1) * t492 - t493 * mrSges(6,3) + t487 * t513 - t514 * t490 + t553;
t421 = -m(5) * t548 + t492 * mrSges(5,1) + mrSges(5,2) * t493 + t513 * t488 + t489 * t514 + t425;
t515 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t522;
t420 = m(4) * t458 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t507 + qJD(3) * t515 - t501 * t523 - t421;
t575 = t540 * t408 + t584 * t420;
t409 = t536 * t412 + t578 * t414;
t574 = t513 * t588 - t514 * t581 - t522 * t579;
t573 = t513 * t579 + t514 * t580 + t522 * t587;
t572 = -t513 * t581 + t514 * t589 - t522 * t580;
t561 = t584 * t408 - t540 * t420;
t555 = -mrSges(3,1) * t538 + mrSges(3,2) * t537;
t552 = mrSges(3,3) * qJDD(1) + t545 * t555;
t453 = Ifges(7,4) * t475 + Ifges(7,2) * t474 + Ifges(7,6) * t519;
t454 = Ifges(7,1) * t475 + Ifges(7,4) * t474 + Ifges(7,5) * t519;
t547 = mrSges(7,1) * t427 - mrSges(7,2) * t428 + Ifges(7,5) * t448 + Ifges(7,6) * t447 + Ifges(7,3) * t504 + t475 * t453 - t474 * t454;
t546 = m(4) * t505 + mrSges(4,1) * t506 + mrSges(4,2) * t507 + t515 * t522 + t516 * t523 + t409;
t526 = (Ifges(3,5) * t537 + Ifges(3,6) * t538) * qJD(1);
t520 = -qJDD(1) * pkin(1) - qJ(2) * t545 + t556;
t508 = -t524 * t537 + t562;
t497 = Ifges(4,1) * t523 - Ifges(4,4) * t522 + Ifges(4,5) * qJD(3);
t496 = Ifges(4,4) * t523 - Ifges(4,2) * t522 + Ifges(4,6) * qJD(3);
t495 = Ifges(4,5) * t523 - Ifges(4,6) * t522 + Ifges(4,3) * qJD(3);
t452 = Ifges(7,5) * t475 + Ifges(7,6) * t474 + Ifges(7,3) * t519;
t418 = mrSges(7,2) * t432 - mrSges(7,3) * t427 + Ifges(7,1) * t448 + Ifges(7,4) * t447 + Ifges(7,5) * t504 + t452 * t474 - t453 * t519;
t417 = -mrSges(7,1) * t432 + mrSges(7,3) * t428 + Ifges(7,4) * t448 + Ifges(7,2) * t447 + Ifges(7,6) * t504 - t452 * t475 + t454 * t519;
t415 = mrSges(6,2) * t493 + t477 * t514 - t549;
t405 = mrSges(3,3) * t545 * t570 + m(3) * t520 + qJDD(1) * t555 + t546;
t404 = -mrSges(5,2) * t548 + mrSges(6,2) * t435 - mrSges(5,3) * t437 - mrSges(6,3) * t436 - pkin(8) * t416 - qJ(5) * t425 - t417 * t539 + t418 * t542 - t581 * t492 + t493 * t589 - t580 * t506 + t573 * t513 + t574 * t522;
t403 = mrSges(5,1) * t548 - mrSges(6,1) * t436 + mrSges(6,2) * t434 + mrSges(5,3) * t438 - pkin(4) * t425 - pkin(5) * t553 - pkin(8) * t559 - t542 * t417 - t539 * t418 - t492 * t588 + t581 * t493 + t579 * t506 + t573 * t514 + t572 * t522;
t402 = -t523 * t495 - mrSges(4,1) * t505 + Ifges(4,4) * t507 + qJD(3) * t497 + mrSges(4,3) * t459 - mrSges(5,1) * t437 + mrSges(5,2) * t438 - mrSges(6,3) * t434 + mrSges(6,1) * t435 + Ifges(4,6) * qJDD(3) + pkin(5) * t416 + pkin(4) * t415 + (qJ(5) * t477 - t572) * t513 + (mrSges(6,2) * qJ(5) + t579) * t492 + t574 * t514 - pkin(3) * t409 + t547 + (-Ifges(4,2) + t587) * t506 + t580 * t493 - qJ(5) * t550;
t401 = mrSges(4,2) * t505 - mrSges(4,3) * t458 + Ifges(4,1) * t507 - Ifges(4,4) * t506 + Ifges(4,5) * qJDD(3) - qJ(4) * t409 - qJD(3) * t496 - t536 * t403 + t404 * t578 - t522 * t495;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t563 - mrSges(2,2) * t557 + t537 * (t538 * qJD(1) * t526 + mrSges(3,2) * t520 - mrSges(3,3) * t508 + t584 * t401 - t540 * t402 - pkin(7) * t575 + (Ifges(3,1) * t537 + Ifges(3,4) * t538) * qJDD(1)) + t538 * (-t526 * t569 - mrSges(3,1) * t520 + mrSges(3,3) * t509 + t540 * t401 + t584 * t402 - pkin(2) * t546 + pkin(7) * t561 + (Ifges(3,4) * t537 + Ifges(3,2) * t538) * qJDD(1)) - pkin(1) * t405 + qJ(2) * ((m(3) * t509 + t538 * t552 + t561) * t538 + (-m(3) * t508 + t537 * t552 - t575) * t537); t405; mrSges(4,1) * t458 - mrSges(4,2) * t459 + Ifges(4,5) * t507 - Ifges(4,6) * t506 + Ifges(4,3) * qJDD(3) - pkin(3) * t421 + qJ(4) * t560 + t403 * t578 + t536 * t404 + t523 * t496 + t522 * t497; t421; t415; t547;];
tauJ  = t1;
