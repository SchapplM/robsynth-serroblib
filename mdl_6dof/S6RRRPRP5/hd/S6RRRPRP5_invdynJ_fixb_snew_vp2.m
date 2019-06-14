% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP5
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
% Datum: 2019-05-07 07:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:52:00
% EndTime: 2019-05-07 07:52:12
% DurationCPUTime: 7.30s
% Computational Cost: add. (74819->324), mult. (156577->398), div. (0->0), fcn. (112532->10), ass. (0->127)
t554 = Ifges(6,1) + Ifges(7,1);
t543 = Ifges(6,4) - Ifges(7,5);
t552 = Ifges(7,4) + Ifges(6,5);
t553 = Ifges(6,2) + Ifges(7,3);
t551 = Ifges(6,6) - Ifges(7,6);
t550 = -Ifges(6,3) - Ifges(7,2);
t518 = sin(qJ(3));
t521 = cos(qJ(3));
t519 = sin(qJ(2));
t540 = qJD(1) * t519;
t501 = qJD(2) * t521 - t518 * t540;
t502 = qJD(2) * t518 + t521 * t540;
t515 = sin(pkin(10));
t516 = cos(pkin(10));
t478 = t501 * t516 - t502 * t515;
t479 = t501 * t515 + t502 * t516;
t517 = sin(qJ(5));
t545 = cos(qJ(5));
t454 = -t478 * t545 + t479 * t517;
t455 = t517 * t478 + t479 * t545;
t522 = cos(qJ(2));
t539 = qJD(1) * t522;
t511 = qJD(3) - t539;
t510 = qJD(5) + t511;
t549 = t454 * t553 - t455 * t543 - t510 * t551;
t548 = -t543 * t454 + t455 * t554 + t552 * t510;
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
t459 = (-t505 - t536) * pkin(8) + (-t506 + t512) * pkin(2) + t496;
t530 = -g(1) * t523 - g(2) * t520;
t497 = -pkin(1) * t525 + qJDD(1) * pkin(7) + t530;
t485 = -g(3) * t519 + t522 * t497;
t504 = (-pkin(2) * t522 - pkin(8) * t519) * qJD(1);
t524 = qJD(2) ^ 2;
t462 = -pkin(2) * t524 + qJDD(2) * pkin(8) + t504 * t539 + t485;
t437 = t521 * t459 - t518 * t462;
t476 = qJD(3) * t501 + qJDD(2) * t518 + t505 * t521;
t500 = qJDD(3) - t506;
t420 = (t501 * t511 - t476) * qJ(4) + (t501 * t502 + t500) * pkin(3) + t437;
t438 = t518 * t459 + t521 * t462;
t475 = -qJD(3) * t502 + qJDD(2) * t521 - t505 * t518;
t482 = pkin(3) * t511 - qJ(4) * t502;
t499 = t501 ^ 2;
t422 = -pkin(3) * t499 + qJ(4) * t475 - t482 * t511 + t438;
t402 = -0.2e1 * qJD(4) * t479 + t516 * t420 - t515 * t422;
t451 = t475 * t515 + t476 * t516;
t399 = (t478 * t511 - t451) * pkin(9) + (t478 * t479 + t500) * pkin(4) + t402;
t403 = 0.2e1 * qJD(4) * t478 + t515 * t420 + t516 * t422;
t450 = t475 * t516 - t476 * t515;
t465 = pkin(4) * t511 - pkin(9) * t479;
t477 = t478 ^ 2;
t401 = -pkin(4) * t477 + pkin(9) * t450 - t465 * t511 + t403;
t395 = t517 * t399 + t545 * t401;
t414 = qJD(5) * t455 - t450 * t545 + t451 * t517;
t443 = mrSges(6,1) * t510 - mrSges(6,3) * t455;
t498 = qJDD(5) + t500;
t433 = pkin(5) * t454 - qJ(6) * t455;
t509 = t510 ^ 2;
t391 = -pkin(5) * t509 + qJ(6) * t498 + 0.2e1 * qJD(6) * t510 - t433 * t454 + t395;
t444 = -mrSges(7,1) * t510 + mrSges(7,2) * t455;
t537 = m(7) * t391 + t498 * mrSges(7,3) + t510 * t444;
t434 = mrSges(7,1) * t454 - mrSges(7,3) * t455;
t541 = -mrSges(6,1) * t454 - mrSges(6,2) * t455 - t434;
t544 = -mrSges(6,3) - mrSges(7,2);
t380 = m(6) * t395 - mrSges(6,2) * t498 + t414 * t544 - t443 * t510 + t454 * t541 + t537;
t394 = t399 * t545 - t517 * t401;
t415 = -t454 * qJD(5) + t517 * t450 + t451 * t545;
t442 = -mrSges(6,2) * t510 - mrSges(6,3) * t454;
t392 = -t498 * pkin(5) - t509 * qJ(6) + t455 * t433 + qJDD(6) - t394;
t441 = -mrSges(7,2) * t454 + mrSges(7,3) * t510;
t531 = -m(7) * t392 + t498 * mrSges(7,1) + t510 * t441;
t382 = m(6) * t394 + mrSges(6,1) * t498 + t415 * t544 + t442 * t510 + t455 * t541 + t531;
t375 = t517 * t380 + t545 * t382;
t456 = -mrSges(5,1) * t478 + mrSges(5,2) * t479;
t463 = -mrSges(5,2) * t511 + mrSges(5,3) * t478;
t373 = m(5) * t402 + mrSges(5,1) * t500 - mrSges(5,3) * t451 - t456 * t479 + t463 * t511 + t375;
t464 = mrSges(5,1) * t511 - mrSges(5,3) * t479;
t532 = t545 * t380 - t382 * t517;
t374 = m(5) * t403 - mrSges(5,2) * t500 + mrSges(5,3) * t450 + t456 * t478 - t464 * t511 + t532;
t369 = t516 * t373 + t515 * t374;
t448 = Ifges(5,4) * t479 + Ifges(5,2) * t478 + Ifges(5,6) * t511;
t449 = Ifges(5,1) * t479 + Ifges(5,4) * t478 + Ifges(5,5) * t511;
t467 = Ifges(4,4) * t502 + Ifges(4,2) * t501 + Ifges(4,6) * t511;
t468 = Ifges(4,1) * t502 + Ifges(4,4) * t501 + Ifges(4,5) * t511;
t387 = mrSges(7,2) * t415 + t434 * t455 - t531;
t528 = -mrSges(6,1) * t394 + mrSges(7,1) * t392 + mrSges(6,2) * t395 - mrSges(7,3) * t391 + pkin(5) * t387 - qJ(6) * t537 + t550 * t498 + t549 * t455 + (qJ(6) * t434 - t548) * t454 - t552 * t415 + (mrSges(7,2) * qJ(6) + t551) * t414;
t547 = mrSges(4,1) * t437 + mrSges(5,1) * t402 - mrSges(4,2) * t438 - mrSges(5,2) * t403 + Ifges(4,5) * t476 + Ifges(5,5) * t451 + Ifges(4,6) * t475 + Ifges(5,6) * t450 + pkin(3) * t369 + pkin(4) * t375 + t479 * t448 - t478 * t449 + t502 * t467 - t501 * t468 + (Ifges(4,3) + Ifges(5,3)) * t500 - t528;
t542 = t551 * t454 - t455 * t552 + t550 * t510;
t484 = -t522 * g(3) - t519 * t497;
t480 = -mrSges(4,1) * t501 + mrSges(4,2) * t502;
t481 = -mrSges(4,2) * t511 + mrSges(4,3) * t501;
t367 = m(4) * t437 + mrSges(4,1) * t500 - mrSges(4,3) * t476 - t480 * t502 + t481 * t511 + t369;
t483 = mrSges(4,1) * t511 - mrSges(4,3) * t502;
t533 = -t373 * t515 + t516 * t374;
t368 = m(4) * t438 - mrSges(4,2) * t500 + mrSges(4,3) * t475 + t480 * t501 - t483 * t511 + t533;
t534 = -t367 * t518 + t521 * t368;
t461 = -qJDD(2) * pkin(2) - pkin(8) * t524 + t504 * t540 - t484;
t436 = -pkin(3) * t475 - qJ(4) * t499 + t502 * t482 + qJDD(4) + t461;
t405 = -pkin(4) * t450 - pkin(9) * t477 + t479 * t465 + t436;
t397 = t405 + (t454 * t510 - t415) * qJ(6) - 0.2e1 * qJD(6) * t455 + (t455 * t510 + t414) * pkin(5);
t388 = m(7) * t397 + t414 * mrSges(7,1) - t415 * mrSges(7,3) + t454 * t441 - t455 * t444;
t363 = t367 * t521 + t368 * t518;
t529 = m(6) * t405 + t414 * mrSges(6,1) + t415 * mrSges(6,2) + t454 * t442 + t455 * t443 + t388;
t383 = m(5) * t436 - t450 * mrSges(5,1) + t451 * mrSges(5,2) - t478 * t463 + t479 * t464 + t529;
t527 = -m(4) * t461 + t475 * mrSges(4,1) - t476 * mrSges(4,2) + t501 * t481 - t502 * t483 - t383;
t508 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t539;
t507 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t540;
t503 = (-mrSges(3,1) * t522 + mrSges(3,2) * t519) * qJD(1);
t495 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t519 + Ifges(3,4) * t522) * qJD(1);
t494 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t519 + Ifges(3,2) * t522) * qJD(1);
t466 = Ifges(4,5) * t502 + Ifges(4,6) * t501 + Ifges(4,3) * t511;
t447 = Ifges(5,5) * t479 + Ifges(5,6) * t478 + Ifges(5,3) * t511;
t377 = mrSges(6,2) * t405 + mrSges(7,2) * t392 - mrSges(6,3) * t394 - mrSges(7,3) * t397 - qJ(6) * t388 - t543 * t414 + t415 * t554 + t542 * t454 + t552 * t498 + t549 * t510;
t376 = -mrSges(6,1) * t405 - mrSges(7,1) * t397 + mrSges(7,2) * t391 + mrSges(6,3) * t395 - pkin(5) * t388 - t414 * t553 + t543 * t415 + t542 * t455 + t551 * t498 + t548 * t510;
t365 = mrSges(5,2) * t436 - mrSges(5,3) * t402 + Ifges(5,1) * t451 + Ifges(5,4) * t450 + Ifges(5,5) * t500 - pkin(9) * t375 - t517 * t376 + t377 * t545 + t478 * t447 - t511 * t448;
t364 = -mrSges(5,1) * t436 + mrSges(5,3) * t403 + Ifges(5,4) * t451 + Ifges(5,2) * t450 + Ifges(5,6) * t500 - pkin(4) * t529 + pkin(9) * t532 + t376 * t545 + t517 * t377 - t479 * t447 + t511 * t449;
t362 = mrSges(4,2) * t461 - mrSges(4,3) * t437 + Ifges(4,1) * t476 + Ifges(4,4) * t475 + Ifges(4,5) * t500 - qJ(4) * t369 - t364 * t515 + t365 * t516 + t466 * t501 - t467 * t511;
t361 = -mrSges(4,1) * t461 + mrSges(4,3) * t438 + Ifges(4,4) * t476 + Ifges(4,2) * t475 + Ifges(4,6) * t500 - pkin(3) * t383 + qJ(4) * t533 + t516 * t364 + t515 * t365 - t502 * t466 + t511 * t468;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t535 - mrSges(2,2) * t530 + t519 * (mrSges(3,2) * t496 - mrSges(3,3) * t484 + Ifges(3,1) * t505 + Ifges(3,4) * t506 + Ifges(3,5) * qJDD(2) - pkin(8) * t363 - qJD(2) * t494 - t518 * t361 + t521 * t362) + t522 * (-mrSges(3,1) * t496 + mrSges(3,3) * t485 + Ifges(3,4) * t505 + Ifges(3,2) * t506 + Ifges(3,6) * qJDD(2) - pkin(2) * t363 + qJD(2) * t495 - t547) + pkin(1) * (-m(3) * t496 + mrSges(3,1) * t506 - mrSges(3,2) * t505 + (-t507 * t519 + t508 * t522) * qJD(1) - t363) + pkin(7) * (t522 * (m(3) * t485 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t506 - qJD(2) * t507 + t503 * t539 + t534) - t519 * (m(3) * t484 + qJDD(2) * mrSges(3,1) - t505 * mrSges(3,3) + qJD(2) * t508 - t503 * t540 + t527)); Ifges(3,5) * t505 + Ifges(3,6) * t506 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t484 - mrSges(3,2) * t485 + t518 * t362 + t521 * t361 + pkin(2) * t527 + pkin(8) * t534 + (t519 * t494 - t522 * t495) * qJD(1); t547; t383; -t528; t387;];
tauJ  = t1;
