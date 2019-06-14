% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-05-06 08:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:22:40
% EndTime: 2019-05-06 08:22:47
% DurationCPUTime: 4.17s
% Computational Cost: add. (35792->329), mult. (85329->403), div. (0->0), fcn. (58082->10), ass. (0->133)
t580 = -2 * qJD(4);
t579 = -Ifges(5,1) - Ifges(4,3);
t573 = Ifges(4,5) - Ifges(5,4);
t578 = Ifges(4,2) + Ifges(5,3);
t577 = Ifges(5,2) + Ifges(4,1);
t572 = Ifges(4,6) - Ifges(5,5);
t571 = -Ifges(5,6) - Ifges(4,4);
t532 = sin(pkin(9));
t538 = cos(qJ(2));
t562 = qJD(1) * t538;
t535 = sin(qJ(2));
t563 = qJD(1) * t535;
t570 = cos(pkin(9));
t510 = t532 * t563 - t570 * t562;
t511 = (t532 * t538 + t570 * t535) * qJD(1);
t480 = pkin(3) * t510 - qJ(4) * t511;
t540 = qJD(2) ^ 2;
t558 = qJD(1) * qJD(2);
t520 = qJDD(1) * t535 + t538 * t558;
t541 = qJD(1) ^ 2;
t536 = sin(qJ(1));
t539 = cos(qJ(1));
t552 = -g(1) * t539 - g(2) * t536;
t517 = -pkin(1) * t541 + qJDD(1) * pkin(7) + t552;
t569 = t535 * t517;
t575 = pkin(2) * t541;
t460 = qJDD(2) * pkin(2) - t520 * qJ(3) - t569 + (qJ(3) * t558 + t535 * t575 - g(3)) * t538;
t497 = -g(3) * t535 + t538 * t517;
t521 = qJDD(1) * t538 - t535 * t558;
t522 = qJD(2) * pkin(2) - qJ(3) * t563;
t530 = t538 ^ 2;
t462 = qJ(3) * t521 - qJD(2) * t522 - t530 * t575 + t497;
t568 = t532 * t460 + t570 * t462;
t576 = pkin(3) * t540 - qJDD(2) * qJ(4) + qJD(2) * t580 + t510 * t480 - t568;
t446 = -0.2e1 * qJD(3) * t511 + t570 * t460 - t532 * t462;
t556 = t536 * g(1) - t539 * g(2);
t550 = -qJDD(1) * pkin(1) - t556;
t466 = -t521 * pkin(2) + qJDD(3) + t522 * t563 + (-qJ(3) * t530 - pkin(7)) * t541 + t550;
t487 = t520 * t532 - t570 * t521;
t488 = t570 * t520 + t532 * t521;
t499 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t511;
t561 = qJD(2) * t510;
t542 = (-t488 + t561) * qJ(4) + t466 + (pkin(3) * qJD(2) + t580) * t511;
t436 = t487 * pkin(3) + t542;
t502 = mrSges(5,1) * t511 + qJD(2) * mrSges(5,2);
t434 = -qJDD(2) * pkin(3) - t540 * qJ(4) + t511 * t480 + qJDD(4) - t446;
t427 = (t510 * t511 - qJDD(2)) * qJ(5) + (t488 + t561) * pkin(4) + t434;
t500 = pkin(4) * t511 - qJD(2) * qJ(5);
t509 = t510 ^ 2;
t431 = -t509 * pkin(4) - t511 * t500 + (pkin(3) + qJ(5)) * t487 + t542;
t531 = sin(pkin(10));
t533 = cos(pkin(10));
t493 = qJD(2) * t533 + t510 * t531;
t421 = -0.2e1 * qJD(5) * t493 + t533 * t427 - t531 * t431;
t472 = qJDD(2) * t533 + t487 * t531;
t492 = -qJD(2) * t531 + t510 * t533;
t419 = (t492 * t511 - t472) * pkin(8) + (t492 * t493 + t488) * pkin(5) + t421;
t422 = 0.2e1 * qJD(5) * t492 + t531 * t427 + t533 * t431;
t469 = pkin(5) * t511 - pkin(8) * t493;
t471 = -qJDD(2) * t531 + t487 * t533;
t491 = t492 ^ 2;
t420 = -pkin(5) * t491 + pkin(8) * t471 - t469 * t511 + t422;
t534 = sin(qJ(6));
t537 = cos(qJ(6));
t417 = t419 * t537 - t420 * t534;
t456 = t492 * t537 - t493 * t534;
t441 = qJD(6) * t456 + t471 * t534 + t472 * t537;
t457 = t492 * t534 + t493 * t537;
t448 = -mrSges(7,1) * t456 + mrSges(7,2) * t457;
t508 = qJD(6) + t511;
t449 = -mrSges(7,2) * t508 + mrSges(7,3) * t456;
t486 = qJDD(6) + t488;
t413 = m(7) * t417 + mrSges(7,1) * t486 - mrSges(7,3) * t441 - t448 * t457 + t449 * t508;
t418 = t419 * t534 + t420 * t537;
t440 = -qJD(6) * t457 + t471 * t537 - t472 * t534;
t450 = mrSges(7,1) * t508 - mrSges(7,3) * t457;
t414 = m(7) * t418 - mrSges(7,2) * t486 + mrSges(7,3) * t440 + t448 * t456 - t450 * t508;
t404 = t537 * t413 + t534 * t414;
t461 = -mrSges(6,1) * t492 + mrSges(6,2) * t493;
t467 = -mrSges(6,2) * t511 + mrSges(6,3) * t492;
t402 = m(6) * t421 + mrSges(6,1) * t488 - mrSges(6,3) * t472 - t461 * t493 + t467 * t511 + t404;
t468 = mrSges(6,1) * t511 - mrSges(6,3) * t493;
t553 = -t413 * t534 + t537 * t414;
t403 = m(6) * t422 - mrSges(6,2) * t488 + mrSges(6,3) * t471 + t461 * t492 - t468 * t511 + t553;
t554 = -t531 * t402 + t533 * t403;
t548 = -m(5) * t436 + t488 * mrSges(5,3) + t511 * t502 - t554;
t501 = mrSges(5,1) * t510 - qJD(2) * mrSges(5,3);
t564 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t510 - t501;
t574 = mrSges(4,1) - mrSges(5,2);
t397 = m(4) * t466 + t488 * mrSges(4,2) + t574 * t487 + t511 * t499 + t564 * t510 - t548;
t481 = mrSges(4,1) * t510 + mrSges(4,2) * t511;
t400 = t533 * t402 + t531 * t403;
t482 = -mrSges(5,2) * t510 - mrSges(5,3) * t511;
t546 = -m(5) * t434 - t488 * mrSges(5,1) - t511 * t482 - t400;
t396 = m(4) * t446 - t488 * mrSges(4,3) + t564 * qJD(2) + t574 * qJDD(2) - t511 * t481 + t546;
t560 = qJD(3) * t510;
t505 = -0.2e1 * t560;
t447 = t505 + t568;
t429 = -pkin(4) * t487 - qJ(5) * t509 + qJD(2) * t500 + qJDD(5) + t505 - t576;
t424 = -pkin(5) * t471 - pkin(8) * t491 + t469 * t493 + t429;
t547 = m(7) * t424 - t440 * mrSges(7,1) + t441 * mrSges(7,2) - t456 * t449 + t457 * t450;
t415 = m(6) * t429 - t471 * mrSges(6,1) + t472 * mrSges(6,2) - t492 * t467 + t493 * t468 + t547;
t432 = 0.2e1 * t560 + t576;
t543 = -m(5) * t432 + qJDD(2) * mrSges(5,3) + qJD(2) * t502 + t415;
t409 = t543 - qJDD(2) * mrSges(4,2) + (-mrSges(4,3) - mrSges(5,1)) * t487 - qJD(2) * t499 + m(4) * t447 + (-t481 - t482) * t510;
t392 = t570 * t396 + t532 * t409;
t567 = qJD(2) * t579 + t510 * t572 - t511 * t573;
t566 = -qJD(2) * t572 + t510 * t578 + t511 * t571;
t565 = qJD(2) * t573 + t510 * t571 + t511 * t577;
t555 = -t396 * t532 + t570 * t409;
t443 = Ifges(7,4) * t457 + Ifges(7,2) * t456 + Ifges(7,6) * t508;
t444 = Ifges(7,1) * t457 + Ifges(7,4) * t456 + Ifges(7,5) * t508;
t545 = mrSges(7,1) * t417 - mrSges(7,2) * t418 + Ifges(7,5) * t441 + Ifges(7,6) * t440 + Ifges(7,3) * t486 + t457 * t443 - t456 * t444;
t524 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t562;
t523 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t563;
t519 = (-mrSges(3,1) * t538 + mrSges(3,2) * t535) * qJD(1);
t516 = -t541 * pkin(7) + t550;
t514 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t535 + Ifges(3,4) * t538) * qJD(1);
t513 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t535 + Ifges(3,2) * t538) * qJD(1);
t496 = -t538 * g(3) - t569;
t453 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t511;
t452 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t511;
t451 = Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t511;
t442 = Ifges(7,5) * t457 + Ifges(7,6) * t456 + Ifges(7,3) * t508;
t406 = mrSges(7,2) * t424 - mrSges(7,3) * t417 + Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t486 + t442 * t456 - t443 * t508;
t405 = -mrSges(7,1) * t424 + mrSges(7,3) * t418 + Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t486 - t442 * t457 + t444 * t508;
t399 = qJDD(2) * mrSges(5,2) + qJD(2) * t501 - t546;
t398 = -t487 * mrSges(5,2) - t510 * t501 - t548;
t394 = mrSges(6,2) * t429 - mrSges(6,3) * t421 + Ifges(6,1) * t472 + Ifges(6,4) * t471 + Ifges(6,5) * t488 - pkin(8) * t404 - t405 * t534 + t406 * t537 + t451 * t492 - t452 * t511;
t393 = -mrSges(6,1) * t429 + mrSges(6,3) * t422 + Ifges(6,4) * t472 + Ifges(6,2) * t471 + Ifges(6,6) * t488 - pkin(5) * t547 + pkin(8) * t553 + t537 * t405 + t534 * t406 - t493 * t451 + t511 * t453;
t391 = t567 * t510 + t571 * t487 + (Ifges(6,3) + t577) * t488 - t492 * t453 + t493 * t452 + Ifges(6,6) * t471 + Ifges(6,5) * t472 + mrSges(4,2) * t466 - mrSges(4,3) * t446 + mrSges(5,1) * t434 - mrSges(5,3) * t436 - mrSges(6,2) * t422 + mrSges(6,1) * t421 + pkin(5) * t404 + pkin(4) * t400 + t573 * qJDD(2) - qJ(4) * t398 + t545 + t566 * qJD(2);
t390 = -mrSges(4,1) * t466 - mrSges(5,1) * t432 + mrSges(5,2) * t436 + mrSges(4,3) * t447 - pkin(3) * t398 + pkin(4) * t415 - qJ(5) * t554 + t565 * qJD(2) + t572 * qJDD(2) - t533 * t393 - t531 * t394 - t487 * t578 - t571 * t488 + t567 * t511;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t556 - mrSges(2,2) * t552 + t535 * (mrSges(3,2) * t516 - mrSges(3,3) * t496 + Ifges(3,1) * t520 + Ifges(3,4) * t521 + Ifges(3,5) * qJDD(2) - qJ(3) * t392 - qJD(2) * t513 - t532 * t390 + t570 * t391) + t538 * (-mrSges(3,1) * t516 + mrSges(3,3) * t497 + Ifges(3,4) * t520 + Ifges(3,2) * t521 + Ifges(3,6) * qJDD(2) - pkin(2) * t397 + qJ(3) * t555 + qJD(2) * t514 + t570 * t390 + t532 * t391) + pkin(1) * (-m(3) * t516 + t521 * mrSges(3,1) - t520 * mrSges(3,2) + (-t523 * t535 + t524 * t538) * qJD(1) - t397) + pkin(7) * (t538 * (m(3) * t497 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t521 - qJD(2) * t523 + t519 * t562 + t555) - t535 * (m(3) * t496 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t520 + qJD(2) * t524 - t519 * t563 + t392)); pkin(2) * t392 + t533 * t394 + qJ(4) * t543 - t531 * t393 + Ifges(3,5) * t520 + Ifges(3,6) * t521 + mrSges(3,1) * t496 - mrSges(3,2) * t497 + mrSges(4,1) * t446 - mrSges(4,2) * t447 - mrSges(5,3) * t432 + mrSges(5,2) * t434 - qJ(5) * t400 - pkin(3) * t399 - t566 * t511 + (-qJ(4) * t482 + t565) * t510 + t573 * t488 + (-mrSges(5,1) * qJ(4) - t572) * t487 + (t535 * t513 - t538 * t514) * qJD(1) + (Ifges(3,3) - t579) * qJDD(2); t397; t399; t415; t545;];
tauJ  = t1;
