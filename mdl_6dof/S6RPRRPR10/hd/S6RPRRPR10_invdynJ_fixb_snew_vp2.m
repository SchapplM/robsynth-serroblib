% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-05-05 23:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:53:23
% EndTime: 2019-05-05 23:53:27
% DurationCPUTime: 2.08s
% Computational Cost: add. (13083->269), mult. (25008->320), div. (0->0), fcn. (15163->8), ass. (0->111)
t522 = Ifges(5,1) + Ifges(6,1);
t510 = Ifges(5,4) - Ifges(6,5);
t509 = Ifges(5,5) + Ifges(6,4);
t521 = -Ifges(5,2) - Ifges(6,3);
t508 = Ifges(5,6) - Ifges(6,6);
t519 = Ifges(5,3) + Ifges(6,2);
t474 = sin(qJ(4));
t478 = cos(qJ(3));
t501 = qJD(1) * t478;
t512 = cos(qJ(4));
t455 = -qJD(3) * t512 + t474 * t501;
t475 = sin(qJ(3));
t500 = qJD(1) * qJD(3);
t498 = t475 * t500;
t460 = qJDD(1) * t478 - t498;
t423 = -t455 * qJD(4) + t474 * qJDD(3) + t460 * t512;
t456 = t474 * qJD(3) + t501 * t512;
t429 = mrSges(6,1) * t455 - mrSges(6,3) * t456;
t481 = qJD(1) ^ 2;
t513 = -pkin(1) - pkin(7);
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t492 = -g(1) * t479 - g(2) * t476;
t516 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t492;
t440 = t481 * t513 - t516;
t497 = t478 * t500;
t459 = -qJDD(1) * t475 - t497;
t403 = (-t460 + t498) * pkin(8) + (-t459 + t497) * pkin(3) + t440;
t496 = g(1) * t476 - t479 * g(2);
t487 = -qJ(2) * t481 + qJDD(2) - t496;
t441 = qJDD(1) * t513 + t487;
t432 = -g(3) * t478 + t475 * t441;
t458 = (pkin(3) * t475 - pkin(8) * t478) * qJD(1);
t480 = qJD(3) ^ 2;
t502 = qJD(1) * t475;
t407 = -pkin(3) * t480 + qJDD(3) * pkin(8) - t458 * t502 + t432;
t387 = t403 * t512 - t474 * t407;
t428 = pkin(4) * t455 - qJ(5) * t456;
t454 = qJDD(4) - t459;
t466 = qJD(4) + t502;
t465 = t466 ^ 2;
t385 = -t454 * pkin(4) - t465 * qJ(5) + t456 * t428 + qJDD(5) - t387;
t507 = t455 * t466;
t379 = (-t423 - t507) * pkin(9) + (t455 * t456 - t454) * pkin(5) + t385;
t388 = t474 * t403 + t512 * t407;
t514 = 2 * qJD(5);
t384 = -pkin(4) * t465 + t454 * qJ(5) - t455 * t428 + t466 * t514 + t388;
t422 = qJD(4) * t456 - qJDD(3) * t512 + t460 * t474;
t439 = -pkin(5) * t466 - pkin(9) * t456;
t453 = t455 ^ 2;
t380 = -pkin(5) * t453 + pkin(9) * t422 + t439 * t466 + t384;
t473 = sin(qJ(6));
t477 = cos(qJ(6));
t377 = t379 * t477 - t380 * t473;
t424 = t455 * t477 - t456 * t473;
t395 = qJD(6) * t424 + t422 * t473 + t423 * t477;
t425 = t455 * t473 + t456 * t477;
t401 = -mrSges(7,1) * t424 + mrSges(7,2) * t425;
t464 = qJD(6) - t466;
t408 = -mrSges(7,2) * t464 + mrSges(7,3) * t424;
t452 = qJDD(6) - t454;
t374 = m(7) * t377 + mrSges(7,1) * t452 - mrSges(7,3) * t395 - t401 * t425 + t408 * t464;
t378 = t379 * t473 + t380 * t477;
t394 = -qJD(6) * t425 + t422 * t477 - t423 * t473;
t409 = mrSges(7,1) * t464 - mrSges(7,3) * t425;
t375 = m(7) * t378 - mrSges(7,2) * t452 + mrSges(7,3) * t394 + t401 * t424 - t409 * t464;
t368 = t374 * t477 + t375 * t473;
t436 = -mrSges(6,2) * t455 + mrSges(6,3) * t466;
t485 = -m(6) * t385 + t454 * mrSges(6,1) + t466 * t436 - t368;
t367 = mrSges(6,2) * t423 + t429 * t456 - t485;
t397 = Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t464;
t398 = Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t464;
t484 = -mrSges(7,1) * t377 + mrSges(7,2) * t378 - Ifges(7,5) * t395 - Ifges(7,6) * t394 - Ifges(7,3) * t452 - t425 * t397 + t424 * t398;
t435 = -mrSges(6,1) * t466 + mrSges(6,2) * t456;
t494 = -t374 * t473 + t477 * t375;
t488 = m(6) * t384 + t454 * mrSges(6,3) + t466 * t435 + t494;
t504 = -t510 * t455 + t522 * t456 + t509 * t466;
t505 = t521 * t455 + t510 * t456 + t508 * t466;
t520 = -t508 * t422 + t509 * t423 + t519 * t454 + mrSges(5,1) * t387 - mrSges(6,1) * t385 - mrSges(5,2) * t388 + mrSges(6,3) * t384 - pkin(4) * t367 - pkin(5) * t368 + qJ(5) * (-mrSges(6,2) * t422 - t429 * t455 + t488) + t484 + t505 * t456 + t504 * t455;
t431 = t475 * g(3) + t478 * t441;
t486 = qJDD(3) * pkin(3) + pkin(8) * t480 - t458 * t501 + t431;
t515 = (-t423 + t507) * qJ(5) - t486;
t511 = -mrSges(5,3) - mrSges(6,2);
t434 = mrSges(5,1) * t466 - mrSges(5,3) * t456;
t503 = -mrSges(5,1) * t455 - mrSges(5,2) * t456 - t429;
t363 = m(5) * t388 - mrSges(5,2) * t454 + t422 * t511 - t434 * t466 + t455 * t503 + t488;
t433 = -mrSges(5,2) * t466 - mrSges(5,3) * t455;
t365 = m(5) * t387 + mrSges(5,1) * t454 + t423 * t511 + t433 * t466 + t456 * t503 + t485;
t360 = t474 * t363 + t512 * t365;
t506 = t508 * t455 - t456 * t509 - t519 * t466;
t495 = t512 * t363 - t365 * t474;
t382 = -pkin(9) * t453 + (-pkin(4) - pkin(5)) * t422 + (-pkin(4) * t466 + t439 + t514) * t456 - t515;
t491 = -m(7) * t382 + t394 * mrSges(7,1) - t395 * mrSges(7,2) + t424 * t408 - t425 * t409;
t457 = (mrSges(4,1) * t475 + mrSges(4,2) * t478) * qJD(1);
t461 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t502;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t501;
t386 = -0.2e1 * qJD(5) * t456 + (t456 * t466 + t422) * pkin(4) + t515;
t372 = m(6) * t386 + t422 * mrSges(6,1) - t423 * mrSges(6,3) - t456 * t435 + t455 * t436 + t491;
t483 = m(5) * t486 - t422 * mrSges(5,1) - t423 * mrSges(5,2) - t455 * t433 - t456 * t434 - t372;
t490 = (m(4) * t432 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t459 - qJD(3) * t462 - t457 * t502 + t495) * t475 + (m(4) * t431 + qJDD(3) * mrSges(4,1) - t460 * mrSges(4,3) + qJD(3) * t461 - t457 * t501 + t483) * t478;
t448 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t478 - Ifges(4,4) * t475) * qJD(1);
t447 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t478 - Ifges(4,2) * t475) * qJD(1);
t443 = -qJDD(1) * pkin(1) + t487;
t442 = pkin(1) * t481 + t516;
t396 = Ifges(7,5) * t425 + Ifges(7,6) * t424 + Ifges(7,3) * t464;
t370 = mrSges(7,2) * t382 - mrSges(7,3) * t377 + Ifges(7,1) * t395 + Ifges(7,4) * t394 + Ifges(7,5) * t452 + t396 * t424 - t397 * t464;
t369 = -mrSges(7,1) * t382 + mrSges(7,3) * t378 + Ifges(7,4) * t395 + Ifges(7,2) * t394 + Ifges(7,6) * t452 - t396 * t425 + t398 * t464;
t358 = m(3) * t443 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t481 + t490;
t357 = -mrSges(5,2) * t486 + mrSges(6,2) * t385 - mrSges(5,3) * t387 - mrSges(6,3) * t386 - pkin(9) * t368 - qJ(5) * t372 - t369 * t473 + t370 * t477 - t510 * t422 + t522 * t423 + t509 * t454 + t506 * t455 - t505 * t466;
t356 = mrSges(5,1) * t486 - mrSges(6,1) * t386 + mrSges(6,2) * t384 + mrSges(5,3) * t388 - pkin(4) * t372 - pkin(5) * t491 - pkin(9) * t494 - t477 * t369 - t473 * t370 + t521 * t422 + t510 * t423 + t508 * t454 + t506 * t456 + t504 * t466;
t1 = [mrSges(2,1) * t496 - mrSges(2,2) * t492 + mrSges(3,2) * t443 - mrSges(3,3) * t442 + t478 * (mrSges(4,2) * t440 - mrSges(4,3) * t431 + Ifges(4,1) * t460 + Ifges(4,4) * t459 + Ifges(4,5) * qJDD(3) - pkin(8) * t360 - qJD(3) * t447 - t474 * t356 + t357 * t512) - t475 * (-mrSges(4,1) * t440 + mrSges(4,3) * t432 + Ifges(4,4) * t460 + Ifges(4,2) * t459 + Ifges(4,6) * qJDD(3) - pkin(3) * t360 + qJD(3) * t448 - t520) - pkin(7) * t490 - pkin(1) * t358 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t442 + m(4) * t440 - mrSges(4,1) * t459 + mrSges(3,2) * t481 + mrSges(4,2) * t460 + t360 + qJDD(1) * mrSges(3,3) + (t461 * t475 + t462 * t478) * qJD(1)) * qJ(2); t358; Ifges(4,5) * t460 + Ifges(4,6) * t459 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t431 - mrSges(4,2) * t432 + t474 * t357 + t512 * t356 + pkin(3) * t483 + pkin(8) * t495 + (t447 * t478 + t448 * t475) * qJD(1); t520; t367; -t484;];
tauJ  = t1;
