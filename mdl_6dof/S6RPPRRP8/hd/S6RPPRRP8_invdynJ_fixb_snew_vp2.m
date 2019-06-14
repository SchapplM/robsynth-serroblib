% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:06:30
% EndTime: 2019-05-05 15:06:33
% DurationCPUTime: 1.74s
% Computational Cost: add. (11778->234), mult. (26330->283), div. (0->0), fcn. (17835->8), ass. (0->105)
t511 = Ifges(6,1) + Ifges(7,1);
t501 = Ifges(6,4) - Ifges(7,5);
t500 = -Ifges(6,5) - Ifges(7,4);
t510 = Ifges(6,2) + Ifges(7,3);
t499 = Ifges(6,6) - Ifges(7,6);
t509 = -Ifges(6,3) - Ifges(7,2);
t465 = qJD(1) ^ 2;
t461 = sin(qJ(1));
t463 = cos(qJ(1));
t482 = t461 * g(1) - t463 * g(2);
t471 = -t465 * qJ(2) + qJDD(2) - t482;
t497 = -pkin(1) - qJ(3);
t508 = -(2 * qJD(1) * qJD(3)) + t497 * qJDD(1) + t471;
t458 = cos(pkin(9));
t507 = t458 ^ 2;
t477 = -t463 * g(1) - t461 * g(2);
t506 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t477;
t457 = sin(pkin(9));
t460 = sin(qJ(4));
t462 = cos(qJ(4));
t475 = t457 * t462 + t458 * t460;
t440 = t475 * qJD(1);
t474 = -t457 * t460 + t458 * t462;
t441 = t474 * qJD(1);
t487 = t441 * qJD(4);
t425 = -t475 * qJDD(1) - t487;
t488 = t440 * qJD(4);
t426 = t474 * qJDD(1) - t488;
t459 = sin(qJ(5));
t504 = cos(qJ(5));
t428 = -t504 * qJD(4) + t459 * t441;
t395 = -t428 * qJD(5) + t459 * qJDD(4) + t504 * t426;
t429 = t459 * qJD(4) + t504 * t441;
t400 = mrSges(7,1) * t428 - mrSges(7,3) * t429;
t422 = t457 * g(3) + t508 * t458;
t503 = pkin(3) * t465;
t405 = (-pkin(7) * qJDD(1) - t457 * t503) * t458 + t422;
t423 = -g(3) * t458 + t508 * t457;
t454 = t457 ^ 2;
t485 = qJDD(1) * t457;
t406 = -pkin(7) * t485 - t454 * t503 + t423;
t383 = t460 * t405 + t462 * t406;
t424 = pkin(4) * t440 - pkin(8) * t441;
t464 = qJD(4) ^ 2;
t379 = -pkin(4) * t464 + qJDD(4) * pkin(8) - t424 * t440 + t383;
t470 = qJDD(3) + t506;
t490 = -t454 - t507;
t412 = pkin(3) * t485 + (t490 * pkin(7) + t497) * t465 + t470;
t381 = (-t426 + t488) * pkin(8) + (-t425 + t487) * pkin(4) + t412;
t375 = -t459 * t379 + t504 * t381;
t399 = pkin(5) * t428 - qJ(6) * t429;
t421 = qJDD(5) - t425;
t438 = qJD(5) + t440;
t437 = t438 ^ 2;
t373 = -t421 * pkin(5) - t437 * qJ(6) + t429 * t399 + qJDD(6) - t375;
t407 = -mrSges(7,2) * t428 + mrSges(7,3) * t438;
t478 = -m(7) * t373 + t421 * mrSges(7,1) + t438 * t407;
t369 = t395 * mrSges(7,2) + t429 * t400 - t478;
t376 = t504 * t379 + t459 * t381;
t372 = -pkin(5) * t437 + qJ(6) * t421 + 0.2e1 * qJD(6) * t438 - t399 * t428 + t376;
t394 = t429 * qJD(5) - t504 * qJDD(4) + t459 * t426;
t410 = -mrSges(7,1) * t438 + mrSges(7,2) * t429;
t483 = m(7) * t372 + t421 * mrSges(7,3) + t438 * t410;
t492 = t501 * t428 - t511 * t429 + t500 * t438;
t493 = t510 * t428 - t501 * t429 - t499 * t438;
t505 = -t499 * t394 - t500 * t395 - t509 * t421 - t492 * t428 - t493 * t429 + mrSges(6,1) * t375 - mrSges(7,1) * t373 - mrSges(6,2) * t376 + mrSges(7,3) * t372 - pkin(5) * t369 + qJ(6) * (-t394 * mrSges(7,2) - t428 * t400 + t483);
t502 = -mrSges(6,3) - mrSges(7,2);
t496 = mrSges(4,2) * t458;
t419 = mrSges(5,1) * t440 + mrSges(5,2) * t441;
t434 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t441;
t409 = mrSges(6,1) * t438 - mrSges(6,3) * t429;
t491 = -mrSges(6,1) * t428 - mrSges(6,2) * t429 - t400;
t365 = m(6) * t376 - t421 * mrSges(6,2) + t502 * t394 - t438 * t409 + t491 * t428 + t483;
t408 = -mrSges(6,2) * t438 - mrSges(6,3) * t428;
t367 = m(6) * t375 + t421 * mrSges(6,1) + t502 * t395 + t438 * t408 + t491 * t429 + t478;
t479 = t504 * t365 - t367 * t459;
t358 = m(5) * t383 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t425 - qJD(4) * t434 - t419 * t440 + t479;
t382 = t462 * t405 - t460 * t406;
t433 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t440;
t378 = -qJDD(4) * pkin(4) - t464 * pkin(8) + t441 * t424 - t382;
t374 = -0.2e1 * qJD(6) * t429 + (t428 * t438 - t395) * qJ(6) + (t429 * t438 + t394) * pkin(5) + t378;
t370 = m(7) * t374 + mrSges(7,1) * t394 - t395 * mrSges(7,3) + t407 * t428 - t429 * t410;
t466 = -m(6) * t378 - t394 * mrSges(6,1) - mrSges(6,2) * t395 - t428 * t408 - t409 * t429 - t370;
t362 = m(5) * t382 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t426 + qJD(4) * t433 - t419 * t441 + t466;
t495 = t460 * t358 + t462 * t362;
t360 = t459 * t365 + t504 * t367;
t494 = t499 * t428 + t500 * t429 + t509 * t438;
t481 = t490 * mrSges(4,3);
t480 = t462 * t358 - t460 * t362;
t473 = -qJDD(1) * mrSges(4,3) - t465 * (mrSges(4,1) * t457 + t496);
t476 = (m(4) * t422 + t473 * t458 + t495) * t458 + (m(4) * t423 + t473 * t457 + t480) * t457;
t469 = m(5) * t412 - t425 * mrSges(5,1) + t426 * mrSges(5,2) + t440 * t433 + t441 * t434 + t360;
t432 = t497 * t465 + t470;
t468 = m(4) * t432 + mrSges(4,1) * t485 + qJDD(1) * t496 + t469;
t439 = -qJDD(1) * pkin(1) + t471;
t436 = t465 * pkin(1) - t506;
t415 = Ifges(5,1) * t441 - Ifges(5,4) * t440 + Ifges(5,5) * qJD(4);
t414 = Ifges(5,4) * t441 - Ifges(5,2) * t440 + Ifges(5,6) * qJD(4);
t413 = Ifges(5,5) * t441 - Ifges(5,6) * t440 + Ifges(5,3) * qJD(4);
t359 = mrSges(6,2) * t378 + mrSges(7,2) * t373 - mrSges(6,3) * t375 - mrSges(7,3) * t374 - qJ(6) * t370 - t501 * t394 + t511 * t395 - t500 * t421 + t494 * t428 + t493 * t438;
t355 = -mrSges(6,1) * t378 - mrSges(7,1) * t374 + mrSges(7,2) * t372 + mrSges(6,3) * t376 - pkin(5) * t370 - t510 * t394 + t501 * t395 + t499 * t421 + t494 * t429 - t492 * t438;
t352 = -mrSges(5,1) * t412 + mrSges(5,3) * t383 + Ifges(5,4) * t426 + Ifges(5,2) * t425 + Ifges(5,6) * qJDD(4) - pkin(4) * t360 + qJD(4) * t415 - t441 * t413 - t505;
t351 = m(3) * t439 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t465 + t476;
t350 = mrSges(5,2) * t412 - mrSges(5,3) * t382 + Ifges(5,1) * t426 + Ifges(5,4) * t425 + Ifges(5,5) * qJDD(4) - pkin(8) * t360 - qJD(4) * t414 - t459 * t355 + t504 * t359 - t440 * t413;
t1 = [mrSges(2,1) * t482 - mrSges(2,2) * t477 + mrSges(3,2) * t439 - mrSges(3,3) * t436 + t458 * (mrSges(4,2) * t432 - mrSges(4,3) * t422 - pkin(7) * t495 + t462 * t350 - t460 * t352) - t457 * (-mrSges(4,1) * t432 + mrSges(4,3) * t423 - pkin(3) * t469 + pkin(7) * t480 + t460 * t350 + t462 * t352) - qJ(3) * t476 - pkin(1) * t351 + qJ(2) * (-m(3) * t436 + (mrSges(3,2) + t481) * t465 + t468) + (Ifges(4,1) * t507 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t458 + Ifges(4,2) * t457) * t457) * qJDD(1); t351; t465 * t481 + t468; mrSges(5,1) * t382 - mrSges(5,2) * t383 + Ifges(5,5) * t426 + Ifges(5,6) * t425 + Ifges(5,3) * qJDD(4) + pkin(4) * t466 + pkin(8) * t479 + t504 * t355 + t459 * t359 + t441 * t414 + t440 * t415; t505; t369;];
tauJ  = t1;
