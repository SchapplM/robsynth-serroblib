% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:44:49
% EndTime: 2019-05-07 07:44:58
% DurationCPUTime: 3.16s
% Computational Cost: add. (22422->305), mult. (46169->358), div. (0->0), fcn. (31025->8), ass. (0->125)
t591 = Ifges(6,1) + Ifges(7,1);
t571 = Ifges(6,4) - Ifges(7,5);
t586 = Ifges(7,4) + Ifges(6,5);
t590 = Ifges(4,2) + Ifges(5,3);
t589 = Ifges(6,2) + Ifges(7,3);
t584 = Ifges(4,6) - Ifges(5,5);
t570 = -Ifges(5,6) - Ifges(4,4);
t583 = Ifges(6,6) - Ifges(7,6);
t588 = -2 * qJD(4);
t587 = Ifges(4,1) + Ifges(5,2);
t585 = Ifges(4,5) - Ifges(5,4);
t582 = Ifges(4,3) + Ifges(5,1);
t581 = Ifges(6,3) + Ifges(7,2);
t537 = sin(qJ(3));
t540 = cos(qJ(2));
t561 = qJD(1) * t540;
t538 = sin(qJ(2));
t562 = qJD(1) * t538;
t576 = cos(qJ(3));
t513 = t537 * t562 - t561 * t576;
t534 = qJD(2) + qJD(3);
t536 = sin(qJ(5));
t575 = cos(qJ(5));
t494 = -t513 * t575 + t534 * t536;
t495 = t536 * t513 + t534 * t575;
t514 = (t537 * t540 + t538 * t576) * qJD(1);
t508 = qJD(5) + t514;
t580 = t494 * t589 - t495 * t571 - t508 * t583;
t579 = -t571 * t494 + t495 * t591 + t586 * t508;
t578 = t513 * t590 + t514 * t570 - t534 * t584;
t559 = qJD(1) * qJD(2);
t520 = qJDD(1) * t538 + t540 * t559;
t542 = qJD(1) ^ 2;
t539 = sin(qJ(1));
t541 = cos(qJ(1));
t553 = -g(1) * t541 - g(2) * t539;
t516 = -pkin(1) * t542 + qJDD(1) * pkin(7) + t553;
t568 = t516 * t538;
t574 = pkin(2) * t542;
t458 = qJDD(2) * pkin(2) - pkin(8) * t520 - t568 + (pkin(8) * t559 + t538 * t574 - g(3)) * t540;
t497 = -g(3) * t538 + t540 * t516;
t521 = qJDD(1) * t540 - t538 * t559;
t524 = qJD(2) * pkin(2) - pkin(8) * t562;
t535 = t540 ^ 2;
t459 = pkin(8) * t521 - qJD(2) * t524 - t535 * t574 + t497;
t429 = t537 * t458 + t576 * t459;
t489 = pkin(3) * t513 - qJ(4) * t514;
t532 = t534 ^ 2;
t533 = qJDD(2) + qJDD(3);
t424 = pkin(3) * t532 - t533 * qJ(4) + t513 * t489 + t534 * t588 - t429;
t475 = qJD(3) * t514 + t520 * t537 - t521 * t576;
t476 = -t513 * qJD(3) + t520 * t576 + t537 * t521;
t557 = g(1) * t539 - t541 * g(2);
t552 = -qJDD(1) * pkin(1) - t557;
t477 = -pkin(2) * t521 + t524 * t562 + (-pkin(8) * t535 - pkin(7)) * t542 + t552;
t499 = mrSges(4,1) * t534 - mrSges(4,3) * t514;
t569 = t513 * t534;
t545 = (-t476 + t569) * qJ(4) + t477 + (t534 * pkin(3) + t588) * t514;
t423 = pkin(3) * t475 + t545;
t501 = mrSges(5,1) * t514 + mrSges(5,2) * t534;
t502 = pkin(4) * t514 - pkin(9) * t534;
t509 = t513 ^ 2;
t417 = -pkin(4) * t509 - t502 * t514 + (pkin(3) + pkin(9)) * t475 + t545;
t428 = t458 * t576 - t537 * t459;
t426 = -t533 * pkin(3) - t532 * qJ(4) + t514 * t489 + qJDD(4) - t428;
t419 = (t513 * t514 - t533) * pkin(9) + (t476 + t569) * pkin(4) + t426;
t413 = t575 * t417 + t536 * t419;
t438 = qJD(5) * t495 - t475 * t575 + t533 * t536;
t474 = qJDD(5) + t476;
t480 = mrSges(6,1) * t508 - mrSges(6,3) * t495;
t453 = pkin(5) * t494 - qJ(6) * t495;
t507 = t508 ^ 2;
t409 = -pkin(5) * t507 + qJ(6) * t474 + 0.2e1 * qJD(6) * t508 - t453 * t494 + t413;
t481 = -mrSges(7,1) * t508 + mrSges(7,2) * t495;
t558 = m(7) * t409 + t474 * mrSges(7,3) + t508 * t481;
t454 = mrSges(7,1) * t494 - mrSges(7,3) * t495;
t566 = -mrSges(6,1) * t494 - mrSges(6,2) * t495 - t454;
t572 = -mrSges(6,3) - mrSges(7,2);
t400 = m(6) * t413 - mrSges(6,2) * t474 + t438 * t572 - t480 * t508 + t494 * t566 + t558;
t412 = -t536 * t417 + t419 * t575;
t439 = -t494 * qJD(5) + t536 * t475 + t533 * t575;
t479 = -mrSges(6,2) * t508 - mrSges(6,3) * t494;
t410 = -t474 * pkin(5) - t507 * qJ(6) + t495 * t453 + qJDD(6) - t412;
t478 = -mrSges(7,2) * t494 + mrSges(7,3) * t508;
t554 = -m(7) * t410 + t474 * mrSges(7,1) + t508 * t478;
t401 = m(6) * t412 + mrSges(6,1) * t474 + t439 * t572 + t479 * t508 + t495 * t566 + t554;
t555 = t575 * t400 - t401 * t536;
t550 = -m(5) * t423 + t476 * mrSges(5,3) + t514 * t501 - t555;
t500 = mrSges(5,1) * t513 - mrSges(5,3) * t534;
t563 = -mrSges(4,2) * t534 - mrSges(4,3) * t513 - t500;
t573 = mrSges(4,1) - mrSges(5,2);
t577 = m(4) * t477 + mrSges(4,2) * t476 + t573 * t475 + t499 * t514 + t563 * t513 - t550;
t490 = mrSges(4,1) * t513 + mrSges(4,2) * t514;
t395 = t536 * t400 + t401 * t575;
t491 = -mrSges(5,2) * t513 - mrSges(5,3) * t514;
t549 = m(5) * t426 + t476 * mrSges(5,1) + t514 * t491 + t395;
t388 = m(4) * t428 - t476 * mrSges(4,3) - t514 * t490 + t533 * t573 + t534 * t563 - t549;
t421 = -pkin(4) * t475 - pkin(9) * t509 + t534 * t502 - t424;
t415 = -0.2e1 * qJD(6) * t495 + (t494 * t508 - t439) * qJ(6) + (t495 * t508 + t438) * pkin(5) + t421;
t406 = m(7) * t415 + t438 * mrSges(7,1) - mrSges(7,3) * t439 + t494 * t478 - t481 * t495;
t547 = m(6) * t421 + mrSges(6,1) * t438 + t439 * mrSges(6,2) + t479 * t494 + t495 * t480 + t406;
t546 = -m(5) * t424 + t533 * mrSges(5,3) + t534 * t501 + t547;
t398 = (-t490 - t491) * t513 + t546 + m(4) * t429 - mrSges(4,2) * t533 - t499 * t534 + (-mrSges(4,3) - mrSges(5,1)) * t475;
t386 = t576 * t388 + t537 * t398;
t567 = t583 * t494 - t586 * t495 - t581 * t508;
t565 = t584 * t513 - t585 * t514 - t582 * t534;
t564 = t570 * t513 + t587 * t514 + t585 * t534;
t556 = -t388 * t537 + t576 * t398;
t405 = mrSges(7,2) * t439 + t454 * t495 - t554;
t544 = mrSges(6,1) * t412 - mrSges(7,1) * t410 - mrSges(6,2) * t413 + mrSges(7,3) * t409 - pkin(5) * t405 + qJ(6) * t558 - t580 * t495 + (-qJ(6) * t454 + t579) * t494 + t581 * t474 + t586 * t439 + (-qJ(6) * mrSges(7,2) - t583) * t438;
t391 = t533 * mrSges(5,2) + t534 * t500 + t549;
t392 = -mrSges(6,1) * t421 - mrSges(7,1) * t415 + mrSges(7,2) * t409 + mrSges(6,3) * t413 - pkin(5) * t406 - t438 * t589 + t571 * t439 + t583 * t474 + t567 * t495 + t579 * t508;
t394 = mrSges(6,2) * t421 + mrSges(7,2) * t410 - mrSges(6,3) * t412 - mrSges(7,3) * t415 - qJ(6) * t406 - t571 * t438 + t439 * t591 + t586 * t474 + t567 * t494 + t580 * t508;
t543 = -mrSges(4,2) * t429 - mrSges(5,3) * t424 - pkin(3) * t391 - pkin(9) * t395 - t536 * t392 + t575 * t394 + t513 * t564 + qJ(4) * (-t491 * t513 + t546) + mrSges(5,2) * t426 + mrSges(4,1) * t428 + t582 * t533 - t578 * t514 + t585 * t476 + (-qJ(4) * mrSges(5,1) - t584) * t475;
t523 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t561;
t522 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t562;
t519 = (-mrSges(3,1) * t540 + mrSges(3,2) * t538) * qJD(1);
t515 = -pkin(7) * t542 + t552;
t512 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t538 + Ifges(3,4) * t540) * qJD(1);
t511 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t538 + Ifges(3,2) * t540) * qJD(1);
t496 = -g(3) * t540 - t568;
t389 = -mrSges(5,2) * t475 - t500 * t513 - t550;
t385 = mrSges(5,1) * t426 + mrSges(4,2) * t477 - mrSges(4,3) * t428 - mrSges(5,3) * t423 + pkin(4) * t395 - qJ(4) * t389 + t570 * t475 + t587 * t476 + t565 * t513 + t585 * t533 + t578 * t534 + t544;
t384 = -mrSges(4,1) * t477 - mrSges(5,1) * t424 + mrSges(5,2) * t423 + mrSges(4,3) * t429 - pkin(3) * t389 + pkin(4) * t547 - pkin(9) * t555 - t575 * t392 - t536 * t394 - t475 * t590 - t570 * t476 + t565 * t514 + t584 * t533 + t564 * t534;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t553 + t538 * (mrSges(3,2) * t515 - mrSges(3,3) * t496 + Ifges(3,1) * t520 + Ifges(3,4) * t521 + Ifges(3,5) * qJDD(2) - pkin(8) * t386 - qJD(2) * t511 - t537 * t384 + t385 * t576) + t540 * (-mrSges(3,1) * t515 + mrSges(3,3) * t497 + Ifges(3,4) * t520 + Ifges(3,2) * t521 + Ifges(3,6) * qJDD(2) - pkin(2) * t577 + pkin(8) * t556 + qJD(2) * t512 + t576 * t384 + t537 * t385) + pkin(1) * (-m(3) * t515 + mrSges(3,1) * t521 - mrSges(3,2) * t520 + (-t522 * t538 + t523 * t540) * qJD(1) - t577) + pkin(7) * (t540 * (m(3) * t497 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t521 - qJD(2) * t522 + t519 * t561 + t556) - t538 * (m(3) * t496 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t520 + qJD(2) * t523 - t519 * t562 + t386)); Ifges(3,3) * qJDD(2) + Ifges(3,6) * t521 + Ifges(3,5) * t520 - mrSges(3,2) * t497 + mrSges(3,1) * t496 + pkin(2) * t386 + t543 + (t538 * t511 - t540 * t512) * qJD(1); t543; t391; t544; t405;];
tauJ  = t1;
