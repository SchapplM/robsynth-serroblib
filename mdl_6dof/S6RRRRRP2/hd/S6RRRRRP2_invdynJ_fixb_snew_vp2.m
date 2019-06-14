% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:30:29
% EndTime: 2019-05-08 04:30:41
% DurationCPUTime: 6.43s
% Computational Cost: add. (68436->326), mult. (148209->403), div. (0->0), fcn. (110371->10), ass. (0->131)
t559 = Ifges(6,1) + Ifges(7,1);
t552 = Ifges(6,4) - Ifges(7,5);
t551 = -Ifges(6,5) - Ifges(7,4);
t558 = Ifges(6,2) + Ifges(7,3);
t550 = Ifges(6,6) - Ifges(7,6);
t557 = -Ifges(6,3) - Ifges(7,2);
t519 = sin(qJ(3));
t520 = sin(qJ(2));
t523 = cos(qJ(3));
t524 = cos(qJ(2));
t496 = (t519 * t524 + t520 * t523) * qJD(1);
t541 = qJD(1) * qJD(2);
t501 = qJDD(1) * t520 + t524 * t541;
t502 = qJDD(1) * t524 - t520 * t541;
t471 = -qJD(3) * t496 - t501 * t519 + t502 * t523;
t495 = (-t519 * t520 + t523 * t524) * qJD(1);
t472 = qJD(3) * t495 + t501 * t523 + t502 * t519;
t518 = sin(qJ(4));
t522 = cos(qJ(4));
t482 = t495 * t522 - t496 * t518;
t440 = qJD(4) * t482 + t471 * t518 + t472 * t522;
t483 = t495 * t518 + t496 * t522;
t515 = qJD(2) + qJD(3);
t512 = qJD(4) + t515;
t517 = sin(qJ(5));
t555 = cos(qJ(5));
t467 = t483 * t517 - t512 * t555;
t514 = qJDD(2) + qJDD(3);
t511 = qJDD(4) + t514;
t422 = -t467 * qJD(5) + t440 * t555 + t517 * t511;
t468 = t483 * t555 + t517 * t512;
t449 = mrSges(7,1) * t467 - mrSges(7,3) * t468;
t526 = qJD(1) ^ 2;
t521 = sin(qJ(1));
t525 = cos(qJ(1));
t534 = -g(1) * t525 - g(2) * t521;
t498 = -pkin(1) * t526 + qJDD(1) * pkin(7) + t534;
t548 = t520 * t498;
t554 = pkin(2) * t526;
t465 = qJDD(2) * pkin(2) - t501 * pkin(8) - t548 + (pkin(8) * t541 + t520 * t554 - g(3)) * t524;
t486 = -g(3) * t520 + t524 * t498;
t543 = qJD(1) * t520;
t505 = qJD(2) * pkin(2) - pkin(8) * t543;
t516 = t524 ^ 2;
t466 = pkin(8) * t502 - qJD(2) * t505 - t516 * t554 + t486;
t446 = t523 * t465 - t519 * t466;
t415 = (t495 * t515 - t472) * pkin(9) + (t495 * t496 + t514) * pkin(3) + t446;
t447 = t519 * t465 + t523 * t466;
t489 = pkin(3) * t515 - pkin(9) * t496;
t491 = t495 ^ 2;
t420 = -pkin(3) * t491 + pkin(9) * t471 - t489 * t515 + t447;
t413 = t518 * t415 + t522 * t420;
t460 = -pkin(4) * t482 - pkin(10) * t483;
t510 = t512 ^ 2;
t408 = -pkin(4) * t510 + pkin(10) * t511 + t460 * t482 + t413;
t539 = t521 * g(1) - t525 * g(2);
t533 = -qJDD(1) * pkin(1) - t539;
t473 = -t502 * pkin(2) + t505 * t543 + (-pkin(8) * t516 - pkin(7)) * t526 + t533;
t424 = -t471 * pkin(3) - t491 * pkin(9) + t496 * t489 + t473;
t439 = -qJD(4) * t483 + t471 * t522 - t472 * t518;
t410 = (-t482 * t512 - t440) * pkin(10) + (t483 * t512 - t439) * pkin(4) + t424;
t404 = -t517 * t408 + t410 * t555;
t436 = qJDD(5) - t439;
t448 = pkin(5) * t467 - qJ(6) * t468;
t480 = qJD(5) - t482;
t479 = t480 ^ 2;
t402 = -t436 * pkin(5) - t479 * qJ(6) + t468 * t448 + qJDD(6) - t404;
t451 = -mrSges(7,2) * t467 + mrSges(7,3) * t480;
t535 = -m(7) * t402 + t436 * mrSges(7,1) + t480 * t451;
t398 = t422 * mrSges(7,2) + t468 * t449 - t535;
t405 = t408 * t555 + t517 * t410;
t401 = -pkin(5) * t479 + qJ(6) * t436 + 0.2e1 * qJD(6) * t480 - t448 * t467 + t405;
t421 = qJD(5) * t468 + t440 * t517 - t511 * t555;
t454 = -mrSges(7,1) * t480 + mrSges(7,2) * t468;
t540 = m(7) * t401 + t436 * mrSges(7,3) + t480 * t454;
t545 = t552 * t467 - t559 * t468 + t551 * t480;
t546 = t558 * t467 - t552 * t468 - t550 * t480;
t556 = -t421 * t550 - t422 * t551 - t557 * t436 - t467 * t545 - t468 * t546 + mrSges(6,1) * t404 - mrSges(7,1) * t402 - mrSges(6,2) * t405 + mrSges(7,3) * t401 - pkin(5) * t398 + qJ(6) * (-t421 * mrSges(7,2) - t467 * t449 + t540);
t553 = -mrSges(6,3) - mrSges(7,2);
t459 = -mrSges(5,1) * t482 + mrSges(5,2) * t483;
t475 = mrSges(5,1) * t512 - mrSges(5,3) * t483;
t453 = mrSges(6,1) * t480 - mrSges(6,3) * t468;
t544 = -mrSges(6,1) * t467 - mrSges(6,2) * t468 - t449;
t393 = m(6) * t405 - t436 * mrSges(6,2) + t421 * t553 - t480 * t453 + t467 * t544 + t540;
t452 = -mrSges(6,2) * t480 - mrSges(6,3) * t467;
t395 = m(6) * t404 + t436 * mrSges(6,1) + t422 * t553 + t480 * t452 + t468 * t544 + t535;
t536 = t393 * t555 - t395 * t517;
t382 = m(5) * t413 - mrSges(5,2) * t511 + mrSges(5,3) * t439 + t459 * t482 - t475 * t512 + t536;
t412 = t522 * t415 - t518 * t420;
t474 = -mrSges(5,2) * t512 + mrSges(5,3) * t482;
t407 = -t511 * pkin(4) - t510 * pkin(10) + t483 * t460 - t412;
t403 = -0.2e1 * qJD(6) * t468 + (t467 * t480 - t422) * qJ(6) + (t468 * t480 + t421) * pkin(5) + t407;
t399 = m(7) * t403 + mrSges(7,1) * t421 - t422 * mrSges(7,3) + t451 * t467 - t468 * t454;
t529 = -m(6) * t407 - t421 * mrSges(6,1) - mrSges(6,2) * t422 - t467 * t452 - t453 * t468 - t399;
t390 = m(5) * t412 + mrSges(5,1) * t511 - mrSges(5,3) * t440 - t459 * t483 + t474 * t512 + t529;
t379 = t518 * t382 + t522 * t390;
t484 = -mrSges(4,1) * t495 + mrSges(4,2) * t496;
t487 = -mrSges(4,2) * t515 + mrSges(4,3) * t495;
t376 = m(4) * t446 + mrSges(4,1) * t514 - mrSges(4,3) * t472 - t484 * t496 + t487 * t515 + t379;
t488 = mrSges(4,1) * t515 - mrSges(4,3) * t496;
t537 = t522 * t382 - t390 * t518;
t377 = m(4) * t447 - mrSges(4,2) * t514 + mrSges(4,3) * t471 + t484 * t495 - t488 * t515 + t537;
t370 = t523 * t376 + t519 * t377;
t388 = t517 * t393 + t395 * t555;
t547 = t550 * t467 + t551 * t468 + t557 * t480;
t542 = qJD(1) * t524;
t538 = -t376 * t519 + t523 * t377;
t532 = -m(5) * t424 + mrSges(5,1) * t439 - t440 * mrSges(5,2) + t474 * t482 - t483 * t475 - t388;
t384 = -mrSges(6,1) * t407 - mrSges(7,1) * t403 + mrSges(7,2) * t401 + mrSges(6,3) * t405 - pkin(5) * t399 - t558 * t421 + t552 * t422 + t550 * t436 + t547 * t468 - t545 * t480;
t386 = mrSges(6,2) * t407 + mrSges(7,2) * t402 - mrSges(6,3) * t404 - mrSges(7,3) * t403 - qJ(6) * t399 - t552 * t421 + t559 * t422 - t551 * t436 + t547 * t467 + t546 * t480;
t456 = Ifges(5,4) * t483 + Ifges(5,2) * t482 + Ifges(5,6) * t512;
t457 = Ifges(5,1) * t483 + Ifges(5,4) * t482 + Ifges(5,5) * t512;
t531 = mrSges(5,1) * t412 - mrSges(5,2) * t413 + Ifges(5,5) * t440 + Ifges(5,6) * t439 + Ifges(5,3) * t511 + pkin(4) * t529 + pkin(10) * t536 + t384 * t555 + t517 * t386 + t483 * t456 - t482 * t457;
t528 = m(4) * t473 - mrSges(4,1) * t471 + mrSges(4,2) * t472 - t487 * t495 + t488 * t496 - t532;
t477 = Ifges(4,4) * t496 + Ifges(4,2) * t495 + Ifges(4,6) * t515;
t478 = Ifges(4,1) * t496 + Ifges(4,4) * t495 + Ifges(4,5) * t515;
t527 = mrSges(4,1) * t446 - mrSges(4,2) * t447 + Ifges(4,5) * t472 + Ifges(4,6) * t471 + Ifges(4,3) * t514 + pkin(3) * t379 + t496 * t477 - t495 * t478 + t531;
t504 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t542;
t503 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t543;
t500 = (-mrSges(3,1) * t524 + mrSges(3,2) * t520) * qJD(1);
t497 = -t526 * pkin(7) + t533;
t494 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t520 + Ifges(3,4) * t524) * qJD(1);
t493 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t520 + Ifges(3,2) * t524) * qJD(1);
t485 = -t524 * g(3) - t548;
t476 = Ifges(4,5) * t496 + Ifges(4,6) * t495 + Ifges(4,3) * t515;
t455 = Ifges(5,5) * t483 + Ifges(5,6) * t482 + Ifges(5,3) * t512;
t372 = -mrSges(5,1) * t424 + mrSges(5,3) * t413 + Ifges(5,4) * t440 + Ifges(5,2) * t439 + Ifges(5,6) * t511 - pkin(4) * t388 - t483 * t455 + t512 * t457 - t556;
t371 = mrSges(5,2) * t424 - mrSges(5,3) * t412 + Ifges(5,1) * t440 + Ifges(5,4) * t439 + Ifges(5,5) * t511 - pkin(10) * t388 - t517 * t384 + t386 * t555 + t482 * t455 - t512 * t456;
t369 = mrSges(4,2) * t473 - mrSges(4,3) * t446 + Ifges(4,1) * t472 + Ifges(4,4) * t471 + Ifges(4,5) * t514 - pkin(9) * t379 + t371 * t522 - t372 * t518 + t476 * t495 - t477 * t515;
t368 = -mrSges(4,1) * t473 + mrSges(4,3) * t447 + Ifges(4,4) * t472 + Ifges(4,2) * t471 + Ifges(4,6) * t514 + pkin(3) * t532 + pkin(9) * t537 + t518 * t371 + t522 * t372 - t496 * t476 + t515 * t478;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t539 - mrSges(2,2) * t534 + t520 * (mrSges(3,2) * t497 - mrSges(3,3) * t485 + Ifges(3,1) * t501 + Ifges(3,4) * t502 + Ifges(3,5) * qJDD(2) - pkin(8) * t370 - qJD(2) * t493 - t519 * t368 + t523 * t369) + t524 * (-mrSges(3,1) * t497 + mrSges(3,3) * t486 + Ifges(3,4) * t501 + Ifges(3,2) * t502 + Ifges(3,6) * qJDD(2) - pkin(2) * t528 + pkin(8) * t538 + qJD(2) * t494 + t523 * t368 + t519 * t369) + pkin(1) * (-t528 + (-t503 * t520 + t504 * t524) * qJD(1) - m(3) * t497 + mrSges(3,1) * t502 - mrSges(3,2) * t501) + pkin(7) * (t524 * (m(3) * t486 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t502 - qJD(2) * t503 + t500 * t542 + t538) - t520 * (m(3) * t485 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t501 + qJD(2) * t504 - t500 * t543 + t370)); Ifges(3,3) * qJDD(2) + t527 + pkin(2) * t370 + Ifges(3,5) * t501 + Ifges(3,6) * t502 + mrSges(3,1) * t485 - mrSges(3,2) * t486 + (t520 * t493 - t524 * t494) * qJD(1); t527; t531; t556; t398;];
tauJ  = t1;
