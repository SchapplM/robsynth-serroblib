% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:20:03
% EndTime: 2019-05-06 18:20:13
% DurationCPUTime: 6.94s
% Computational Cost: add. (68440->324), mult. (148945->398), div. (0->0), fcn. (107786->10), ass. (0->127)
t552 = Ifges(6,1) + Ifges(7,1);
t543 = Ifges(6,4) - Ifges(7,5);
t550 = Ifges(7,4) + Ifges(6,5);
t551 = Ifges(6,2) + Ifges(7,3);
t549 = Ifges(6,6) - Ifges(7,6);
t548 = -Ifges(6,3) - Ifges(7,2);
t515 = sin(pkin(10));
t516 = cos(pkin(10));
t519 = sin(qJ(2));
t540 = qJD(1) * t519;
t500 = qJD(2) * t516 - t515 * t540;
t501 = qJD(2) * t515 + t516 * t540;
t518 = sin(qJ(4));
t521 = cos(qJ(4));
t473 = t500 * t521 - t501 * t518;
t474 = t500 * t518 + t501 * t521;
t517 = sin(qJ(5));
t545 = cos(qJ(5));
t454 = -t473 * t545 + t474 * t517;
t455 = t517 * t473 + t474 * t545;
t522 = cos(qJ(2));
t539 = qJD(1) * t522;
t511 = qJD(4) - t539;
t510 = qJD(5) + t511;
t547 = t551 * t454 - t543 * t455 - t549 * t510;
t546 = -t543 * t454 + t552 * t455 + t550 * t510;
t544 = -mrSges(6,3) - mrSges(7,2);
t525 = qJD(1) ^ 2;
t520 = sin(qJ(1));
t523 = cos(qJ(1));
t535 = g(1) * t520 - t523 * g(2);
t496 = -qJDD(1) * pkin(1) - pkin(7) * t525 - t535;
t538 = qJD(1) * qJD(2);
t536 = t522 * t538;
t505 = qJDD(1) * t519 + t536;
t512 = t519 * t538;
t506 = qJDD(1) * t522 - t512;
t459 = (-t505 - t536) * qJ(3) + (-t506 + t512) * pkin(2) + t496;
t530 = -g(1) * t523 - g(2) * t520;
t497 = -pkin(1) * t525 + qJDD(1) * pkin(7) + t530;
t477 = -g(3) * t519 + t522 * t497;
t503 = (-pkin(2) * t522 - qJ(3) * t519) * qJD(1);
t524 = qJD(2) ^ 2;
t462 = -pkin(2) * t524 + qJDD(2) * qJ(3) + t503 * t539 + t477;
t436 = -0.2e1 * qJD(3) * t501 + t516 * t459 - t462 * t515;
t482 = qJDD(2) * t515 + t505 * t516;
t420 = (-t500 * t539 - t482) * pkin(8) + (t500 * t501 - t506) * pkin(3) + t436;
t437 = 0.2e1 * qJD(3) * t500 + t515 * t459 + t516 * t462;
t481 = qJDD(2) * t516 - t505 * t515;
t483 = -pkin(3) * t539 - pkin(8) * t501;
t499 = t500 ^ 2;
t422 = -pkin(3) * t499 + pkin(8) * t481 + t483 * t539 + t437;
t402 = t521 * t420 - t422 * t518;
t448 = qJD(4) * t473 + t481 * t518 + t482 * t521;
t502 = qJDD(4) - t506;
t399 = (t473 * t511 - t448) * pkin(9) + (t473 * t474 + t502) * pkin(4) + t402;
t403 = t518 * t420 + t521 * t422;
t447 = -qJD(4) * t474 + t481 * t521 - t482 * t518;
t465 = pkin(4) * t511 - pkin(9) * t474;
t472 = t473 ^ 2;
t401 = -pkin(4) * t472 + pkin(9) * t447 - t465 * t511 + t403;
t395 = t517 * t399 + t545 * t401;
t414 = qJD(5) * t455 - t447 * t545 + t448 * t517;
t444 = mrSges(6,1) * t510 - mrSges(6,3) * t455;
t498 = qJDD(5) + t502;
t433 = pkin(5) * t454 - qJ(6) * t455;
t509 = t510 ^ 2;
t391 = -pkin(5) * t509 + qJ(6) * t498 + 0.2e1 * qJD(6) * t510 - t433 * t454 + t395;
t445 = -mrSges(7,1) * t510 + mrSges(7,2) * t455;
t537 = m(7) * t391 + t498 * mrSges(7,3) + t510 * t445;
t434 = mrSges(7,1) * t454 - mrSges(7,3) * t455;
t541 = -mrSges(6,1) * t454 - mrSges(6,2) * t455 - t434;
t380 = m(6) * t395 - mrSges(6,2) * t498 + t414 * t544 - t444 * t510 + t454 * t541 + t537;
t394 = t399 * t545 - t517 * t401;
t415 = -t454 * qJD(5) + t517 * t447 + t448 * t545;
t443 = -mrSges(6,2) * t510 - mrSges(6,3) * t454;
t392 = -t498 * pkin(5) - t509 * qJ(6) + t455 * t433 + qJDD(6) - t394;
t442 = -mrSges(7,2) * t454 + mrSges(7,3) * t510;
t531 = -m(7) * t392 + t498 * mrSges(7,1) + t510 * t442;
t382 = m(6) * t394 + mrSges(6,1) * t498 + t415 * t544 + t443 * t510 + t455 * t541 + t531;
t377 = t517 * t380 + t545 * t382;
t456 = -mrSges(5,1) * t473 + mrSges(5,2) * t474;
t463 = -mrSges(5,2) * t511 + mrSges(5,3) * t473;
t373 = m(5) * t402 + mrSges(5,1) * t502 - mrSges(5,3) * t448 - t456 * t474 + t463 * t511 + t377;
t464 = mrSges(5,1) * t511 - mrSges(5,3) * t474;
t532 = t545 * t380 - t382 * t517;
t374 = m(5) * t403 - mrSges(5,2) * t502 + mrSges(5,3) * t447 + t456 * t473 - t464 * t511 + t532;
t369 = t521 * t373 + t518 * t374;
t542 = t549 * t454 - t550 * t455 + t548 * t510;
t476 = -t522 * g(3) - t519 * t497;
t475 = -mrSges(4,1) * t500 + mrSges(4,2) * t501;
t479 = mrSges(4,2) * t539 + mrSges(4,3) * t500;
t367 = m(4) * t436 - mrSges(4,1) * t506 - mrSges(4,3) * t482 - t475 * t501 - t479 * t539 + t369;
t480 = -mrSges(4,1) * t539 - mrSges(4,3) * t501;
t533 = -t373 * t518 + t521 * t374;
t368 = m(4) * t437 + mrSges(4,2) * t506 + mrSges(4,3) * t481 + t475 * t500 + t480 * t539 + t533;
t534 = -t367 * t515 + t516 * t368;
t461 = -qJDD(2) * pkin(2) - qJ(3) * t524 + t503 * t540 + qJDD(3) - t476;
t438 = -pkin(3) * t481 - pkin(8) * t499 + t501 * t483 + t461;
t405 = -pkin(4) * t447 - pkin(9) * t472 + t474 * t465 + t438;
t397 = t405 + (t455 * t510 + t414) * pkin(5) + (t454 * t510 - t415) * qJ(6) - 0.2e1 * qJD(6) * t455;
t388 = m(7) * t397 + t414 * mrSges(7,1) - t415 * mrSges(7,3) + t454 * t442 - t455 * t445;
t363 = t367 * t516 + t368 * t515;
t529 = m(6) * t405 + t414 * mrSges(6,1) + t415 * mrSges(6,2) + t454 * t443 + t455 * t444 + t388;
t528 = m(5) * t438 - t447 * mrSges(5,1) + t448 * mrSges(5,2) - t473 * t463 + t474 * t464 + t529;
t387 = mrSges(7,2) * t415 + t434 * t455 - t531;
t527 = -mrSges(6,1) * t394 + mrSges(7,1) * t392 + mrSges(6,2) * t395 - mrSges(7,3) * t391 + pkin(5) * t387 - qJ(6) * t537 + t548 * t498 + t547 * t455 + (qJ(6) * t434 - t546) * t454 - t550 * t415 + (mrSges(7,2) * qJ(6) + t549) * t414;
t383 = m(4) * t461 - t481 * mrSges(4,1) + t482 * mrSges(4,2) - t500 * t479 + t501 * t480 + t528;
t450 = Ifges(5,4) * t474 + Ifges(5,2) * t473 + Ifges(5,6) * t511;
t451 = Ifges(5,1) * t474 + Ifges(5,4) * t473 + Ifges(5,5) * t511;
t526 = mrSges(5,1) * t402 - mrSges(5,2) * t403 + Ifges(5,5) * t448 + Ifges(5,6) * t447 + Ifges(5,3) * t502 + pkin(4) * t377 + t474 * t450 - t473 * t451 - t527;
t508 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t539;
t507 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t540;
t504 = (-mrSges(3,1) * t522 + mrSges(3,2) * t519) * qJD(1);
t495 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t519 + Ifges(3,4) * t522) * qJD(1);
t494 = Ifges(3,6) * qJD(2) + (t519 * Ifges(3,4) + Ifges(3,2) * t522) * qJD(1);
t469 = Ifges(4,1) * t501 + Ifges(4,4) * t500 - Ifges(4,5) * t539;
t468 = Ifges(4,4) * t501 + Ifges(4,2) * t500 - Ifges(4,6) * t539;
t467 = Ifges(4,5) * t501 + Ifges(4,6) * t500 - Ifges(4,3) * t539;
t449 = Ifges(5,5) * t474 + Ifges(5,6) * t473 + Ifges(5,3) * t511;
t376 = mrSges(6,2) * t405 + mrSges(7,2) * t392 - mrSges(6,3) * t394 - mrSges(7,3) * t397 - qJ(6) * t388 - t543 * t414 + t552 * t415 + t542 * t454 + t550 * t498 + t547 * t510;
t375 = -mrSges(6,1) * t405 - mrSges(7,1) * t397 + mrSges(7,2) * t391 + mrSges(6,3) * t395 - pkin(5) * t388 - t551 * t414 + t543 * t415 + t542 * t455 + t549 * t498 + t546 * t510;
t365 = mrSges(5,2) * t438 - mrSges(5,3) * t402 + Ifges(5,1) * t448 + Ifges(5,4) * t447 + Ifges(5,5) * t502 - pkin(9) * t377 - t517 * t375 + t376 * t545 + t473 * t449 - t511 * t450;
t364 = -mrSges(5,1) * t438 + mrSges(5,3) * t403 + Ifges(5,4) * t448 + Ifges(5,2) * t447 + Ifges(5,6) * t502 - pkin(4) * t529 + pkin(9) * t532 + t375 * t545 + t517 * t376 - t474 * t449 + t511 * t451;
t362 = mrSges(4,2) * t461 - mrSges(4,3) * t436 + Ifges(4,1) * t482 + Ifges(4,4) * t481 - Ifges(4,5) * t506 - pkin(8) * t369 - t364 * t518 + t365 * t521 + t467 * t500 + t468 * t539;
t361 = -mrSges(4,1) * t461 + mrSges(4,3) * t437 + Ifges(4,4) * t482 + Ifges(4,2) * t481 - Ifges(4,6) * t506 - pkin(3) * t528 + pkin(8) * t533 + t521 * t364 + t518 * t365 - t501 * t467 - t469 * t539;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t535 - mrSges(2,2) * t530 + t519 * (mrSges(3,2) * t496 - mrSges(3,3) * t476 + Ifges(3,1) * t505 + Ifges(3,4) * t506 + Ifges(3,5) * qJDD(2) - qJ(3) * t363 - qJD(2) * t494 - t515 * t361 + t516 * t362) + t522 * (-mrSges(3,1) * t496 - mrSges(4,1) * t436 + mrSges(4,2) * t437 + mrSges(3,3) * t477 + Ifges(3,4) * t505 - Ifges(4,5) * t482 + Ifges(3,6) * qJDD(2) - Ifges(4,6) * t481 - pkin(2) * t363 - pkin(3) * t369 + qJD(2) * t495 - t501 * t468 + t500 * t469 - t526 + (Ifges(3,2) + Ifges(4,3)) * t506) + pkin(1) * (-m(3) * t496 + mrSges(3,1) * t506 - mrSges(3,2) * t505 + (-t507 * t519 + t508 * t522) * qJD(1) - t363) + pkin(7) * (t522 * (m(3) * t477 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t506 - qJD(2) * t507 + t504 * t539 + t534) - t519 * (m(3) * t476 + qJDD(2) * mrSges(3,1) - t505 * mrSges(3,3) + qJD(2) * t508 - t504 * t540 - t383)); Ifges(3,5) * t505 + Ifges(3,6) * t506 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t476 - mrSges(3,2) * t477 + t515 * t362 + t516 * t361 - pkin(2) * t383 + qJ(3) * t534 + (t519 * t494 - t522 * t495) * qJD(1); t383; t526; -t527; t387;];
tauJ  = t1;
