% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:08:41
% EndTime: 2019-05-07 20:08:56
% DurationCPUTime: 10.98s
% Computational Cost: add. (178917->350), mult. (361197->445), div. (0->0), fcn. (267487->12), ass. (0->138)
t536 = sin(qJ(3));
t537 = sin(qJ(2));
t541 = cos(qJ(3));
t542 = cos(qJ(2));
t514 = (t536 * t542 + t537 * t541) * qJD(1);
t558 = qJD(1) * qJD(2);
t520 = qJDD(1) * t537 + t542 * t558;
t521 = qJDD(1) * t542 - t537 * t558;
t488 = -t514 * qJD(3) - t536 * t520 + t521 * t541;
t559 = qJD(1) * t542;
t560 = qJD(1) * t537;
t513 = -t536 * t560 + t541 * t559;
t489 = qJD(3) * t513 + t520 * t541 + t521 * t536;
t524 = qJD(2) * pkin(2) - pkin(8) * t560;
t531 = t542 ^ 2;
t544 = qJD(1) ^ 2;
t538 = sin(qJ(1));
t543 = cos(qJ(1));
t557 = t538 * g(1) - t543 * g(2);
t550 = -qJDD(1) * pkin(1) - t557;
t490 = -t521 * pkin(2) + t524 * t560 + (-pkin(8) * t531 - pkin(7)) * t544 + t550;
t530 = qJD(2) + qJD(3);
t440 = (-t513 * t530 - t489) * pkin(9) + (t514 * t530 - t488) * pkin(3) + t490;
t552 = -g(1) * t543 - g(2) * t538;
t516 = -pkin(1) * t544 + qJDD(1) * pkin(7) + t552;
t561 = t537 * t516;
t563 = pkin(2) * t544;
t480 = qJDD(2) * pkin(2) - t520 * pkin(8) - t561 + (pkin(8) * t558 + t537 * t563 - g(3)) * t542;
t503 = -g(3) * t537 + t542 * t516;
t481 = pkin(8) * t521 - qJD(2) * t524 - t531 * t563 + t503;
t458 = t536 * t480 + t541 * t481;
t498 = -pkin(3) * t513 - pkin(9) * t514;
t528 = t530 ^ 2;
t529 = qJDD(2) + qJDD(3);
t446 = -pkin(3) * t528 + pkin(9) * t529 + t498 * t513 + t458;
t535 = sin(qJ(4));
t540 = cos(qJ(4));
t429 = t540 * t440 - t535 * t446;
t500 = -t514 * t535 + t530 * t540;
t461 = qJD(4) * t500 + t489 * t540 + t529 * t535;
t487 = qJDD(4) - t488;
t501 = t514 * t540 + t530 * t535;
t509 = qJD(4) - t513;
t419 = (t500 * t509 - t461) * qJ(5) + (t500 * t501 + t487) * pkin(4) + t429;
t430 = t535 * t440 + t540 * t446;
t460 = -qJD(4) * t501 - t489 * t535 + t529 * t540;
t492 = pkin(4) * t509 - qJ(5) * t501;
t499 = t500 ^ 2;
t421 = -pkin(4) * t499 + qJ(5) * t460 - t492 * t509 + t430;
t532 = sin(pkin(11));
t533 = cos(pkin(11));
t474 = t500 * t532 + t501 * t533;
t413 = -0.2e1 * qJD(5) * t474 + t533 * t419 - t532 * t421;
t443 = t460 * t532 + t461 * t533;
t473 = t500 * t533 - t501 * t532;
t411 = (t473 * t509 - t443) * pkin(10) + (t473 * t474 + t487) * pkin(5) + t413;
t414 = 0.2e1 * qJD(5) * t473 + t532 * t419 + t533 * t421;
t442 = t460 * t533 - t461 * t532;
t464 = pkin(5) * t509 - pkin(10) * t474;
t472 = t473 ^ 2;
t412 = -pkin(5) * t472 + pkin(10) * t442 - t464 * t509 + t414;
t534 = sin(qJ(6));
t539 = cos(qJ(6));
t409 = t411 * t539 - t412 * t534;
t453 = t473 * t539 - t474 * t534;
t427 = qJD(6) * t453 + t442 * t534 + t443 * t539;
t454 = t473 * t534 + t474 * t539;
t437 = -mrSges(7,1) * t453 + mrSges(7,2) * t454;
t507 = qJD(6) + t509;
t447 = -mrSges(7,2) * t507 + mrSges(7,3) * t453;
t483 = qJDD(6) + t487;
t402 = m(7) * t409 + mrSges(7,1) * t483 - mrSges(7,3) * t427 - t437 * t454 + t447 * t507;
t410 = t411 * t534 + t412 * t539;
t426 = -qJD(6) * t454 + t442 * t539 - t443 * t534;
t448 = mrSges(7,1) * t507 - mrSges(7,3) * t454;
t403 = m(7) * t410 - mrSges(7,2) * t483 + mrSges(7,3) * t426 + t437 * t453 - t448 * t507;
t396 = t539 * t402 + t534 * t403;
t455 = -mrSges(6,1) * t473 + mrSges(6,2) * t474;
t462 = -mrSges(6,2) * t509 + mrSges(6,3) * t473;
t394 = m(6) * t413 + mrSges(6,1) * t487 - mrSges(6,3) * t443 - t455 * t474 + t462 * t509 + t396;
t463 = mrSges(6,1) * t509 - mrSges(6,3) * t474;
t553 = -t402 * t534 + t539 * t403;
t395 = m(6) * t414 - mrSges(6,2) * t487 + mrSges(6,3) * t442 + t455 * t473 - t463 * t509 + t553;
t390 = t533 * t394 + t532 * t395;
t450 = Ifges(6,4) * t474 + Ifges(6,2) * t473 + Ifges(6,6) * t509;
t451 = Ifges(6,1) * t474 + Ifges(6,4) * t473 + Ifges(6,5) * t509;
t466 = Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t509;
t467 = Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t509;
t433 = Ifges(7,4) * t454 + Ifges(7,2) * t453 + Ifges(7,6) * t507;
t434 = Ifges(7,1) * t454 + Ifges(7,4) * t453 + Ifges(7,5) * t507;
t549 = -mrSges(7,1) * t409 + mrSges(7,2) * t410 - Ifges(7,5) * t427 - Ifges(7,6) * t426 - Ifges(7,3) * t483 - t454 * t433 + t453 * t434;
t564 = mrSges(5,1) * t429 + mrSges(6,1) * t413 - mrSges(5,2) * t430 - mrSges(6,2) * t414 + Ifges(5,5) * t461 + Ifges(6,5) * t443 + Ifges(5,6) * t460 + Ifges(6,6) * t442 + pkin(4) * t390 + pkin(5) * t396 + t474 * t450 - t473 * t451 + t501 * t466 - t500 * t467 - (-Ifges(6,3) - Ifges(5,3)) * t487 - t549;
t497 = -mrSges(4,1) * t513 + mrSges(4,2) * t514;
t505 = mrSges(4,1) * t530 - mrSges(4,3) * t514;
t478 = -mrSges(5,1) * t500 + mrSges(5,2) * t501;
t491 = -mrSges(5,2) * t509 + mrSges(5,3) * t500;
t388 = m(5) * t429 + mrSges(5,1) * t487 - mrSges(5,3) * t461 - t478 * t501 + t491 * t509 + t390;
t493 = mrSges(5,1) * t509 - mrSges(5,3) * t501;
t554 = -t394 * t532 + t533 * t395;
t389 = m(5) * t430 - mrSges(5,2) * t487 + mrSges(5,3) * t460 + t478 * t500 - t493 * t509 + t554;
t555 = -t388 * t535 + t540 * t389;
t380 = m(4) * t458 - mrSges(4,2) * t529 + mrSges(4,3) * t488 + t497 * t513 - t505 * t530 + t555;
t457 = t480 * t541 - t536 * t481;
t504 = -mrSges(4,2) * t530 + mrSges(4,3) * t513;
t445 = -pkin(3) * t529 - pkin(9) * t528 + t514 * t498 - t457;
t431 = -pkin(4) * t460 - qJ(5) * t499 + t501 * t492 + qJDD(5) + t445;
t416 = -pkin(5) * t442 - pkin(10) * t472 + t464 * t474 + t431;
t551 = m(7) * t416 - t426 * mrSges(7,1) + t427 * mrSges(7,2) - t453 * t447 + t454 * t448;
t407 = m(6) * t431 - t442 * mrSges(6,1) + mrSges(6,2) * t443 - t473 * t462 + t463 * t474 + t551;
t546 = -m(5) * t445 + t460 * mrSges(5,1) - mrSges(5,2) * t461 + t500 * t491 - t493 * t501 - t407;
t405 = m(4) * t457 + mrSges(4,1) * t529 - mrSges(4,3) * t489 - t497 * t514 + t504 * t530 + t546;
t377 = t536 * t380 + t541 * t405;
t382 = t540 * t388 + t535 * t389;
t556 = t541 * t380 - t405 * t536;
t432 = Ifges(7,5) * t454 + Ifges(7,6) * t453 + Ifges(7,3) * t507;
t397 = -mrSges(7,1) * t416 + mrSges(7,3) * t410 + Ifges(7,4) * t427 + Ifges(7,2) * t426 + Ifges(7,6) * t483 - t432 * t454 + t434 * t507;
t398 = mrSges(7,2) * t416 - mrSges(7,3) * t409 + Ifges(7,1) * t427 + Ifges(7,4) * t426 + Ifges(7,5) * t483 + t432 * t453 - t433 * t507;
t449 = Ifges(6,5) * t474 + Ifges(6,6) * t473 + Ifges(6,3) * t509;
t383 = -mrSges(6,1) * t431 + mrSges(6,3) * t414 + Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t487 - pkin(5) * t551 + pkin(10) * t553 + t539 * t397 + t534 * t398 - t474 * t449 + t509 * t451;
t384 = mrSges(6,2) * t431 - mrSges(6,3) * t413 + Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t487 - pkin(10) * t396 - t397 * t534 + t398 * t539 + t449 * t473 - t450 * t509;
t465 = Ifges(5,5) * t501 + Ifges(5,6) * t500 + Ifges(5,3) * t509;
t374 = -mrSges(5,1) * t445 + mrSges(5,3) * t430 + Ifges(5,4) * t461 + Ifges(5,2) * t460 + Ifges(5,6) * t487 - pkin(4) * t407 + qJ(5) * t554 + t533 * t383 + t532 * t384 - t501 * t465 + t509 * t467;
t376 = mrSges(5,2) * t445 - mrSges(5,3) * t429 + Ifges(5,1) * t461 + Ifges(5,4) * t460 + Ifges(5,5) * t487 - qJ(5) * t390 - t383 * t532 + t384 * t533 + t465 * t500 - t466 * t509;
t495 = Ifges(4,4) * t514 + Ifges(4,2) * t513 + Ifges(4,6) * t530;
t496 = Ifges(4,1) * t514 + Ifges(4,4) * t513 + Ifges(4,5) * t530;
t548 = mrSges(4,1) * t457 - mrSges(4,2) * t458 + Ifges(4,5) * t489 + Ifges(4,6) * t488 + Ifges(4,3) * t529 + pkin(3) * t546 + pkin(9) * t555 + t540 * t374 + t535 * t376 + t514 * t495 - t513 * t496;
t547 = m(4) * t490 - t488 * mrSges(4,1) + t489 * mrSges(4,2) - t513 * t504 + t514 * t505 + t382;
t523 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t559;
t522 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t560;
t519 = (-mrSges(3,1) * t542 + mrSges(3,2) * t537) * qJD(1);
t515 = -t544 * pkin(7) + t550;
t512 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t537 + Ifges(3,4) * t542) * qJD(1);
t511 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t537 + Ifges(3,2) * t542) * qJD(1);
t502 = -t542 * g(3) - t561;
t494 = Ifges(4,5) * t514 + Ifges(4,6) * t513 + Ifges(4,3) * t530;
t372 = -mrSges(4,1) * t490 + mrSges(4,3) * t458 + Ifges(4,4) * t489 + Ifges(4,2) * t488 + Ifges(4,6) * t529 - pkin(3) * t382 - t514 * t494 + t530 * t496 - t564;
t371 = mrSges(4,2) * t490 - mrSges(4,3) * t457 + Ifges(4,1) * t489 + Ifges(4,4) * t488 + Ifges(4,5) * t529 - pkin(9) * t382 - t374 * t535 + t376 * t540 + t494 * t513 - t495 * t530;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t552 + t537 * (mrSges(3,2) * t515 - mrSges(3,3) * t502 + Ifges(3,1) * t520 + Ifges(3,4) * t521 + Ifges(3,5) * qJDD(2) - pkin(8) * t377 - qJD(2) * t511 + t541 * t371 - t536 * t372) + t542 * (-mrSges(3,1) * t515 + mrSges(3,3) * t503 + Ifges(3,4) * t520 + Ifges(3,2) * t521 + Ifges(3,6) * qJDD(2) - pkin(2) * t547 + pkin(8) * t556 + qJD(2) * t512 + t536 * t371 + t541 * t372) + pkin(1) * (-m(3) * t515 + t521 * mrSges(3,1) - t520 * mrSges(3,2) + (-t522 * t537 + t523 * t542) * qJD(1) - t547) + pkin(7) * (t542 * (m(3) * t503 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t521 - qJD(2) * t522 + t519 * t559 + t556) - t537 * (m(3) * t502 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t520 + qJD(2) * t523 - t519 * t560 + t377)); pkin(2) * t377 + Ifges(3,3) * qJDD(2) + Ifges(3,5) * t520 + Ifges(3,6) * t521 + mrSges(3,1) * t502 - mrSges(3,2) * t503 + (t537 * t511 - t542 * t512) * qJD(1) + t548; t548; t564; t407; -t549;];
tauJ  = t1;
