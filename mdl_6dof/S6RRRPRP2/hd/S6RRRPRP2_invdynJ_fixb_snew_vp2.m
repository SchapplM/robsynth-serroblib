% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:31:13
% EndTime: 2019-05-07 07:31:21
% DurationCPUTime: 5.90s
% Computational Cost: add. (59737->326), mult. (134225->404), div. (0->0), fcn. (98316->10), ass. (0->128)
t565 = Ifges(6,1) + Ifges(7,1);
t558 = Ifges(6,4) - Ifges(7,5);
t557 = -Ifges(6,5) - Ifges(7,4);
t564 = Ifges(6,2) + Ifges(7,3);
t556 = Ifges(6,6) - Ifges(7,6);
t563 = -Ifges(6,3) - Ifges(7,2);
t526 = sin(qJ(2));
t529 = cos(qJ(2));
t546 = qJD(1) * qJD(2);
t508 = qJDD(1) * t526 + t529 * t546;
t531 = qJD(1) ^ 2;
t527 = sin(qJ(1));
t530 = cos(qJ(1));
t538 = -g(1) * t530 - g(2) * t527;
t505 = -pkin(1) * t531 + qJDD(1) * pkin(7) + t538;
t554 = t505 * t526;
t560 = pkin(2) * t531;
t470 = qJDD(2) * pkin(2) - pkin(8) * t508 - t554 + (pkin(8) * t546 + t526 * t560 - g(3)) * t529;
t493 = -g(3) * t526 + t529 * t505;
t509 = qJDD(1) * t529 - t526 * t546;
t549 = qJD(1) * t526;
t512 = qJD(2) * pkin(2) - pkin(8) * t549;
t521 = t529 ^ 2;
t471 = pkin(8) * t509 - qJD(2) * t512 - t521 * t560 + t493;
t525 = sin(qJ(3));
t528 = cos(qJ(3));
t443 = t528 * t470 - t471 * t525;
t502 = (-t525 * t526 + t528 * t529) * qJD(1);
t477 = qJD(3) * t502 + t508 * t528 + t509 * t525;
t503 = (t525 * t529 + t526 * t528) * qJD(1);
t519 = qJDD(2) + qJDD(3);
t520 = qJD(2) + qJD(3);
t420 = (t502 * t520 - t477) * qJ(4) + (t502 * t503 + t519) * pkin(3) + t443;
t444 = t525 * t470 + t528 * t471;
t476 = -qJD(3) * t503 - t508 * t525 + t509 * t528;
t495 = pkin(3) * t520 - qJ(4) * t503;
t498 = t502 ^ 2;
t423 = -pkin(3) * t498 + qJ(4) * t476 - t495 * t520 + t444;
t522 = sin(pkin(10));
t523 = cos(pkin(10));
t490 = t502 * t522 + t503 * t523;
t415 = -0.2e1 * qJD(4) * t490 + t420 * t523 - t522 * t423;
t452 = t476 * t522 + t477 * t523;
t524 = sin(qJ(5));
t561 = cos(qJ(5));
t474 = t490 * t524 - t520 * t561;
t427 = -t474 * qJD(5) + t452 * t561 + t524 * t519;
t475 = t490 * t561 + t524 * t520;
t454 = mrSges(7,1) * t474 - mrSges(7,3) * t475;
t489 = t502 * t523 - t503 * t522;
t416 = 0.2e1 * qJD(4) * t489 + t522 * t420 + t523 * t423;
t465 = -pkin(4) * t489 - pkin(9) * t490;
t518 = t520 ^ 2;
t413 = -pkin(4) * t518 + pkin(9) * t519 + t465 * t489 + t416;
t544 = g(1) * t527 - t530 * g(2);
t537 = -qJDD(1) * pkin(1) - t544;
t478 = -pkin(2) * t509 + t512 * t549 + (-pkin(8) * t521 - pkin(7)) * t531 + t537;
t429 = -pkin(3) * t476 - qJ(4) * t498 + t503 * t495 + qJDD(4) + t478;
t451 = t476 * t523 - t477 * t522;
t418 = (-t489 * t520 - t452) * pkin(9) + (t490 * t520 - t451) * pkin(4) + t429;
t409 = -t524 * t413 + t418 * t561;
t450 = qJDD(5) - t451;
t453 = pkin(5) * t474 - qJ(6) * t475;
t484 = qJD(5) - t489;
t483 = t484 ^ 2;
t407 = -t450 * pkin(5) - t483 * qJ(6) + t475 * t453 + qJDD(6) - t409;
t456 = -mrSges(7,2) * t474 + mrSges(7,3) * t484;
t539 = -m(7) * t407 + t450 * mrSges(7,1) + t484 * t456;
t403 = mrSges(7,2) * t427 + t454 * t475 - t539;
t410 = t413 * t561 + t524 * t418;
t406 = -pkin(5) * t483 + qJ(6) * t450 + 0.2e1 * qJD(6) * t484 - t453 * t474 + t410;
t426 = qJD(5) * t475 + t452 * t524 - t519 * t561;
t459 = -mrSges(7,1) * t484 + mrSges(7,2) * t475;
t545 = m(7) * t406 + t450 * mrSges(7,3) + t484 * t459;
t551 = t558 * t474 - t565 * t475 + t557 * t484;
t552 = t564 * t474 - t558 * t475 - t556 * t484;
t562 = -t426 * t556 - t427 * t557 - t563 * t450 - t474 * t551 - t475 * t552 + mrSges(6,1) * t409 - mrSges(7,1) * t407 - mrSges(6,2) * t410 + mrSges(7,3) * t406 - pkin(5) * t403 + qJ(6) * (-mrSges(7,2) * t426 - t454 * t474 + t545);
t559 = -mrSges(6,3) - mrSges(7,2);
t464 = -mrSges(5,1) * t489 + mrSges(5,2) * t490;
t480 = mrSges(5,1) * t520 - mrSges(5,3) * t490;
t458 = mrSges(6,1) * t484 - mrSges(6,3) * t475;
t550 = -mrSges(6,1) * t474 - mrSges(6,2) * t475 - t454;
t398 = m(6) * t410 - mrSges(6,2) * t450 + t426 * t559 - t458 * t484 + t474 * t550 + t545;
t457 = -mrSges(6,2) * t484 - mrSges(6,3) * t474;
t400 = m(6) * t409 + mrSges(6,1) * t450 + t427 * t559 + t457 * t484 + t475 * t550 + t539;
t541 = t398 * t561 - t400 * t524;
t386 = m(5) * t416 - mrSges(5,2) * t519 + mrSges(5,3) * t451 + t464 * t489 - t480 * t520 + t541;
t479 = -mrSges(5,2) * t520 + mrSges(5,3) * t489;
t412 = -pkin(4) * t519 - pkin(9) * t518 + t490 * t465 - t415;
t408 = -0.2e1 * qJD(6) * t475 + (t474 * t484 - t427) * qJ(6) + (t475 * t484 + t426) * pkin(5) + t412;
t404 = m(7) * t408 + mrSges(7,1) * t426 - t427 * mrSges(7,3) + t456 * t474 - t475 * t459;
t534 = -m(6) * t412 - t426 * mrSges(6,1) - mrSges(6,2) * t427 - t474 * t457 - t458 * t475 - t404;
t395 = m(5) * t415 + mrSges(5,1) * t519 - mrSges(5,3) * t452 - t464 * t490 + t479 * t520 + t534;
t383 = t522 * t386 + t523 * t395;
t491 = -mrSges(4,1) * t502 + mrSges(4,2) * t503;
t494 = -mrSges(4,2) * t520 + mrSges(4,3) * t502;
t380 = m(4) * t443 + mrSges(4,1) * t519 - mrSges(4,3) * t477 - t491 * t503 + t494 * t520 + t383;
t496 = mrSges(4,1) * t520 - mrSges(4,3) * t503;
t542 = t523 * t386 - t395 * t522;
t381 = m(4) * t444 - mrSges(4,2) * t519 + mrSges(4,3) * t476 + t491 * t502 - t496 * t520 + t542;
t374 = t528 * t380 + t525 * t381;
t393 = t524 * t398 + t400 * t561;
t553 = t556 * t474 + t557 * t475 + t563 * t484;
t548 = qJD(1) * t529;
t543 = -t380 * t525 + t528 * t381;
t536 = -m(5) * t429 + mrSges(5,1) * t451 - t452 * mrSges(5,2) + t479 * t489 - t490 * t480 - t393;
t533 = m(4) * t478 - mrSges(4,1) * t476 + mrSges(4,2) * t477 - t494 * t502 + t496 * t503 - t536;
t388 = -mrSges(6,1) * t412 - mrSges(7,1) * t408 + mrSges(7,2) * t406 + mrSges(6,3) * t410 - pkin(5) * t404 - t564 * t426 + t558 * t427 + t556 * t450 + t553 * t475 - t551 * t484;
t391 = mrSges(6,2) * t412 + mrSges(7,2) * t407 - mrSges(6,3) * t409 - mrSges(7,3) * t408 - qJ(6) * t404 - t558 * t426 + t565 * t427 - t557 * t450 + t553 * t474 + t552 * t484;
t461 = Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * t520;
t462 = Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t520;
t486 = Ifges(4,4) * t503 + Ifges(4,2) * t502 + Ifges(4,6) * t520;
t487 = Ifges(4,1) * t503 + Ifges(4,4) * t502 + Ifges(4,5) * t520;
t532 = mrSges(4,1) * t443 + mrSges(5,1) * t415 - mrSges(4,2) * t444 - mrSges(5,2) * t416 + pkin(3) * t383 + pkin(4) * t534 + pkin(9) * t541 + t561 * t388 + t524 * t391 + t490 * t461 - t489 * t462 - t502 * t487 + Ifges(5,6) * t451 + Ifges(5,5) * t452 + t503 * t486 + Ifges(4,6) * t476 + Ifges(4,5) * t477 + (Ifges(5,3) + Ifges(4,3)) * t519;
t511 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t548;
t510 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t549;
t507 = (-mrSges(3,1) * t529 + mrSges(3,2) * t526) * qJD(1);
t504 = -pkin(7) * t531 + t537;
t501 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t526 + Ifges(3,4) * t529) * qJD(1);
t500 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t526 + Ifges(3,2) * t529) * qJD(1);
t492 = -g(3) * t529 - t554;
t485 = Ifges(4,5) * t503 + Ifges(4,6) * t502 + Ifges(4,3) * t520;
t460 = Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * t520;
t376 = -mrSges(5,1) * t429 + mrSges(5,3) * t416 + Ifges(5,4) * t452 + Ifges(5,2) * t451 + Ifges(5,6) * t519 - pkin(4) * t393 - t490 * t460 + t520 * t462 - t562;
t375 = mrSges(5,2) * t429 - mrSges(5,3) * t415 + Ifges(5,1) * t452 + Ifges(5,4) * t451 + Ifges(5,5) * t519 - pkin(9) * t393 - t524 * t388 + t391 * t561 + t489 * t460 - t520 * t461;
t373 = mrSges(4,2) * t478 - mrSges(4,3) * t443 + Ifges(4,1) * t477 + Ifges(4,4) * t476 + Ifges(4,5) * t519 - qJ(4) * t383 + t375 * t523 - t376 * t522 + t485 * t502 - t486 * t520;
t372 = -mrSges(4,1) * t478 + mrSges(4,3) * t444 + Ifges(4,4) * t477 + Ifges(4,2) * t476 + Ifges(4,6) * t519 + pkin(3) * t536 + qJ(4) * t542 + t522 * t375 + t523 * t376 - t503 * t485 + t520 * t487;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t544 - mrSges(2,2) * t538 + t526 * (mrSges(3,2) * t504 - mrSges(3,3) * t492 + Ifges(3,1) * t508 + Ifges(3,4) * t509 + Ifges(3,5) * qJDD(2) - pkin(8) * t374 - qJD(2) * t500 - t525 * t372 + t528 * t373) + t529 * (-mrSges(3,1) * t504 + mrSges(3,3) * t493 + Ifges(3,4) * t508 + Ifges(3,2) * t509 + Ifges(3,6) * qJDD(2) - pkin(2) * t533 + pkin(8) * t543 + qJD(2) * t501 + t528 * t372 + t525 * t373) + pkin(1) * (-t533 + (-t510 * t526 + t511 * t529) * qJD(1) - m(3) * t504 + mrSges(3,1) * t509 - mrSges(3,2) * t508) + pkin(7) * (t529 * (m(3) * t493 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t509 - qJD(2) * t510 + t507 * t548 + t543) - t526 * (m(3) * t492 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t508 + qJD(2) * t511 - t507 * t549 + t374)); Ifges(3,5) * t508 + Ifges(3,6) * t509 + mrSges(3,1) * t492 - mrSges(3,2) * t493 + Ifges(3,3) * qJDD(2) + t532 + pkin(2) * t374 + (t526 * t500 - t529 * t501) * qJD(1); t532; -t536; t562; t403;];
tauJ  = t1;
