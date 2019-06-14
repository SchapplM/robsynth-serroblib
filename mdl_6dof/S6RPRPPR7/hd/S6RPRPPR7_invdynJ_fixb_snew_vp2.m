% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-05-05 17:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:15:23
% EndTime: 2019-05-05 17:15:25
% DurationCPUTime: 1.73s
% Computational Cost: add. (10502->275), mult. (23154->327), div. (0->0), fcn. (14283->8), ass. (0->114)
t521 = -2 * qJD(4);
t520 = Ifges(5,1) + Ifges(6,2);
t519 = -Ifges(6,1) - Ifges(5,3);
t513 = Ifges(5,4) + Ifges(6,6);
t512 = Ifges(5,5) - Ifges(6,4);
t518 = -Ifges(5,2) - Ifges(6,3);
t511 = Ifges(5,6) - Ifges(6,5);
t481 = qJD(1) ^ 2;
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t497 = t476 * g(1) - t479 * g(2);
t487 = -t481 * qJ(2) + qJDD(2) - t497;
t515 = -pkin(1) - pkin(7);
t444 = qJDD(1) * t515 + t487;
t475 = sin(qJ(3));
t478 = cos(qJ(3));
t431 = t475 * g(3) + t478 * t444;
t500 = qJD(1) * qJD(3);
t498 = t475 * t500;
t462 = qJDD(1) * t478 - t498;
t404 = (-t462 - t498) * qJ(4) + (-t475 * t478 * t481 + qJDD(3)) * pkin(3) + t431;
t432 = -g(3) * t478 + t475 * t444;
t461 = -qJDD(1) * t475 - t478 * t500;
t503 = qJD(1) * t478;
t464 = qJD(3) * pkin(3) - qJ(4) * t503;
t472 = t475 ^ 2;
t405 = -pkin(3) * t472 * t481 + qJ(4) * t461 - qJD(3) * t464 + t432;
t473 = cos(pkin(9));
t504 = qJD(1) * t475;
t510 = sin(pkin(9));
t454 = t473 * t503 - t510 * t504;
t389 = t473 * t404 - t510 * t405 + t454 * t521;
t494 = -t479 * g(1) - t476 * g(2);
t490 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t494;
t407 = -t461 * pkin(3) + qJDD(4) + t464 * t503 + (-qJ(4) * t472 + t515) * t481 + t490;
t430 = t461 * t510 + t473 * t462;
t440 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t454;
t517 = m(5) * t407 + t430 * mrSges(5,2) + t454 * t440;
t516 = -2 * qJD(5);
t514 = mrSges(5,1) - mrSges(6,2);
t453 = (t473 * t475 + t478 * t510) * qJD(1);
t422 = mrSges(5,1) * t453 + mrSges(5,2) * t454;
t421 = pkin(4) * t453 - qJ(5) * t454;
t480 = qJD(3) ^ 2;
t386 = -qJDD(3) * pkin(4) - t480 * qJ(5) + t454 * t421 + qJDD(5) - t389;
t502 = qJD(3) * t453;
t381 = (t453 * t454 - qJDD(3)) * pkin(8) + (t430 + t502) * pkin(5) + t386;
t429 = -t473 * t461 + t462 * t510;
t443 = pkin(5) * t454 - qJD(3) * pkin(8);
t452 = t453 ^ 2;
t482 = (-t430 + t502) * qJ(5) + t407 + (qJD(3) * pkin(4) + t516) * t454;
t384 = t482 - t454 * t443 - t452 * pkin(5) + (pkin(4) + pkin(8)) * t429;
t474 = sin(qJ(6));
t477 = cos(qJ(6));
t379 = t381 * t477 - t384 * t474;
t433 = -qJD(3) * t474 + t453 * t477;
t400 = qJD(6) * t433 + qJDD(3) * t477 + t429 * t474;
t434 = qJD(3) * t477 + t453 * t474;
t408 = -mrSges(7,1) * t433 + mrSges(7,2) * t434;
t450 = qJD(6) + t454;
t411 = -mrSges(7,2) * t450 + mrSges(7,3) * t433;
t428 = qJDD(6) + t430;
t376 = m(7) * t379 + mrSges(7,1) * t428 - mrSges(7,3) * t400 - t408 * t434 + t411 * t450;
t380 = t381 * t474 + t384 * t477;
t399 = -qJD(6) * t434 - qJDD(3) * t474 + t429 * t477;
t412 = mrSges(7,1) * t450 - mrSges(7,3) * t434;
t377 = m(7) * t380 - mrSges(7,2) * t428 + mrSges(7,3) * t399 + t408 * t433 - t412 * t450;
t369 = t477 * t376 + t474 * t377;
t423 = -mrSges(6,2) * t453 - mrSges(6,3) * t454;
t485 = -m(6) * t386 - t430 * mrSges(6,1) - t454 * t423 - t369;
t439 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t453;
t441 = mrSges(6,1) * t453 - qJD(3) * mrSges(6,3);
t505 = t439 - t441;
t365 = m(5) * t389 - t430 * mrSges(5,3) + qJD(3) * t505 + qJDD(3) * t514 - t454 * t422 + t485;
t448 = t453 * t521;
t509 = t510 * t404 + t473 * t405;
t390 = t448 + t509;
t488 = t480 * pkin(4) - qJDD(3) * qJ(5) - t509;
t385 = qJD(3) * t516 + ((2 * qJD(4)) + t421) * t453 + t488;
t442 = mrSges(6,1) * t454 + qJD(3) * mrSges(6,2);
t383 = -t429 * pkin(5) - t452 * pkin(8) - t453 * t421 + t448 + ((2 * qJD(5)) + t443) * qJD(3) - t488;
t486 = -m(7) * t383 + t399 * mrSges(7,1) - t400 * mrSges(7,2) + t433 * t411 - t434 * t412;
t483 = -m(6) * t385 + qJDD(3) * mrSges(6,3) + qJD(3) * t442 - t486;
t374 = m(5) * t390 - qJDD(3) * mrSges(5,2) - qJD(3) * t440 + (-t422 - t423) * t453 + (-mrSges(5,3) - mrSges(6,1)) * t429 + t483;
t363 = t473 * t365 + t510 * t374;
t508 = t519 * qJD(3) + t511 * t453 - t512 * t454;
t507 = t511 * qJD(3) + t518 * t453 + t513 * t454;
t506 = t512 * qJD(3) - t513 * t453 + t520 * t454;
t496 = -t474 * t376 + t477 * t377;
t493 = -t365 * t510 + t473 * t374;
t460 = (mrSges(4,1) * t475 + mrSges(4,2) * t478) * qJD(1);
t463 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t504;
t465 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t503;
t492 = (m(4) * t431 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t462 + qJD(3) * t463 - t460 * t503 + t363) * t478 + (m(4) * t432 - qJDD(3) * mrSges(4,2) + t461 * mrSges(4,3) - qJD(3) * t465 - t460 * t504 + t493) * t475;
t388 = t429 * pkin(4) + t482;
t489 = m(6) * t388 - t430 * mrSges(6,3) - t454 * t442 + t496;
t393 = Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t450;
t394 = Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t450;
t484 = mrSges(7,1) * t379 - mrSges(7,2) * t380 + Ifges(7,5) * t400 + Ifges(7,6) * t399 + Ifges(7,3) * t428 + t434 * t393 - t433 * t394;
t367 = -t429 * mrSges(6,2) - t453 * t441 + t489;
t457 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t478 - Ifges(4,4) * t475) * qJD(1);
t456 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t478 - Ifges(4,2) * t475) * qJD(1);
t451 = -qJDD(1) * pkin(1) + t487;
t445 = t481 * pkin(1) - t490;
t438 = t481 * t515 + t490;
t392 = Ifges(7,5) * t434 + Ifges(7,6) * t433 + Ifges(7,3) * t450;
t371 = mrSges(7,2) * t383 - mrSges(7,3) * t379 + Ifges(7,1) * t400 + Ifges(7,4) * t399 + Ifges(7,5) * t428 + t392 * t433 - t393 * t450;
t370 = -mrSges(7,1) * t383 + mrSges(7,3) * t380 + Ifges(7,4) * t400 + Ifges(7,2) * t399 + Ifges(7,6) * t428 - t392 * t434 + t394 * t450;
t368 = qJDD(3) * mrSges(6,2) + qJD(3) * t441 - t485;
t366 = t429 * t514 + t453 * t505 + t489 + t517;
t360 = mrSges(6,1) * t386 + mrSges(5,2) * t407 - mrSges(5,3) * t389 - mrSges(6,3) * t388 + pkin(5) * t369 - qJ(5) * t367 - t507 * qJD(3) + t512 * qJDD(3) - t513 * t429 + t520 * t430 + t508 * t453 + t484;
t359 = m(3) * t451 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t481 + t492;
t358 = -mrSges(5,1) * t407 - mrSges(6,1) * t385 + mrSges(6,2) * t388 + mrSges(5,3) * t390 - pkin(4) * t367 - pkin(5) * t486 - pkin(8) * t496 + t506 * qJD(3) + t511 * qJDD(3) - t477 * t370 - t474 * t371 + t518 * t429 + t513 * t430 + t508 * t454;
t1 = [mrSges(2,1) * t497 - mrSges(2,2) * t494 + mrSges(3,2) * t451 - mrSges(3,3) * t445 + t478 * (mrSges(4,2) * t438 - mrSges(4,3) * t431 + Ifges(4,1) * t462 + Ifges(4,4) * t461 + Ifges(4,5) * qJDD(3) - qJ(4) * t363 - qJD(3) * t456 - t358 * t510 + t473 * t360) - t475 * (-mrSges(4,1) * t438 + mrSges(4,3) * t432 + Ifges(4,4) * t462 + Ifges(4,2) * t461 + Ifges(4,6) * qJDD(3) - pkin(3) * t366 + qJ(4) * t493 + qJD(3) * t457 + t473 * t358 + t360 * t510) - pkin(7) * t492 - pkin(1) * t359 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t445 + m(4) * t438 - t461 * mrSges(4,1) + t429 * mrSges(5,1) + t481 * mrSges(3,2) + t462 * mrSges(4,2) + t453 * t439 + t367 + qJDD(1) * mrSges(3,3) + (t463 * t475 + t465 * t478) * qJD(1) + t517) * qJ(2); t359; -t474 * t370 + t477 * t371 + qJ(5) * t483 + Ifges(4,6) * t461 + Ifges(4,5) * t462 + mrSges(4,1) * t431 - mrSges(4,2) * t432 + mrSges(6,2) * t386 + mrSges(5,1) * t389 - mrSges(5,2) * t390 - mrSges(6,3) * t385 - pkin(8) * t369 - pkin(4) * t368 + pkin(3) * t363 + t507 * t454 + (-qJ(5) * t423 + t506) * t453 + t512 * t430 + (-mrSges(6,1) * qJ(5) - t511) * t429 + (t456 * t478 + t457 * t475) * qJD(1) + (Ifges(4,3) - t519) * qJDD(3); t366; t368; t484;];
tauJ  = t1;
