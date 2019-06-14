% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:26:59
% EndTime: 2019-05-06 04:27:06
% DurationCPUTime: 4.07s
% Computational Cost: add. (46452->291), mult. (92190->362), div. (0->0), fcn. (62432->10), ass. (0->118)
t474 = sin(qJ(1));
t479 = cos(qJ(1));
t491 = -t479 * g(1) - t474 * g(2);
t501 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t491;
t500 = (-pkin(1) - pkin(7));
t481 = qJD(1) ^ 2;
t439 = (t481 * t500) - t501;
t478 = cos(qJ(3));
t498 = qJD(1) * qJD(3);
t464 = t478 * t498;
t473 = sin(qJ(3));
t457 = -t473 * qJDD(1) - t464;
t496 = t473 * t498;
t458 = qJDD(1) * t478 - t496;
t414 = (-t458 + t496) * pkin(8) + (-t457 + t464) * pkin(3) + t439;
t495 = t474 * g(1) - t479 * g(2);
t487 = -t481 * qJ(2) + qJDD(2) - t495;
t440 = qJDD(1) * t500 + t487;
t434 = -g(3) * t478 + t473 * t440;
t456 = (pkin(3) * t473 - pkin(8) * t478) * qJD(1);
t466 = t473 * qJD(1);
t480 = qJD(3) ^ 2;
t417 = -pkin(3) * t480 + qJDD(3) * pkin(8) - t456 * t466 + t434;
t472 = sin(qJ(4));
t477 = cos(qJ(4));
t395 = t477 * t414 - t472 * t417;
t499 = qJD(1) * t478;
t453 = qJD(3) * t477 - t472 * t499;
t428 = qJD(4) * t453 + qJDD(3) * t472 + t458 * t477;
t452 = qJDD(4) - t457;
t454 = qJD(3) * t472 + t477 * t499;
t463 = t466 + qJD(4);
t385 = (t453 * t463 - t428) * pkin(9) + (t453 * t454 + t452) * pkin(4) + t395;
t396 = t472 * t414 + t477 * t417;
t427 = -qJD(4) * t454 + qJDD(3) * t477 - t458 * t472;
t438 = pkin(4) * t463 - pkin(9) * t454;
t451 = t453 ^ 2;
t387 = -pkin(4) * t451 + pkin(9) * t427 - t438 * t463 + t396;
t471 = sin(qJ(5));
t476 = cos(qJ(5));
t372 = t476 * t385 - t471 * t387;
t430 = t453 * t476 - t454 * t471;
t402 = qJD(5) * t430 + t427 * t471 + t428 * t476;
t431 = t453 * t471 + t454 * t476;
t450 = qJDD(5) + t452;
t462 = qJD(5) + t463;
t369 = (t430 * t462 - t402) * pkin(10) + (t430 * t431 + t450) * pkin(5) + t372;
t373 = t471 * t385 + t476 * t387;
t401 = -qJD(5) * t431 + t427 * t476 - t428 * t471;
t420 = pkin(5) * t462 - pkin(10) * t431;
t429 = t430 ^ 2;
t370 = -pkin(5) * t429 + pkin(10) * t401 - t420 * t462 + t373;
t470 = sin(qJ(6));
t475 = cos(qJ(6));
t367 = t369 * t475 - t370 * t470;
t409 = t430 * t475 - t431 * t470;
t381 = qJD(6) * t409 + t401 * t470 + t402 * t475;
t410 = t430 * t470 + t431 * t475;
t393 = -mrSges(7,1) * t409 + mrSges(7,2) * t410;
t459 = qJD(6) + t462;
t403 = -mrSges(7,2) * t459 + mrSges(7,3) * t409;
t446 = qJDD(6) + t450;
t364 = m(7) * t367 + mrSges(7,1) * t446 - mrSges(7,3) * t381 - t393 * t410 + t403 * t459;
t368 = t369 * t470 + t370 * t475;
t380 = -qJD(6) * t410 + t401 * t475 - t402 * t470;
t404 = mrSges(7,1) * t459 - mrSges(7,3) * t410;
t365 = m(7) * t368 - mrSges(7,2) * t446 + mrSges(7,3) * t380 + t393 * t409 - t404 * t459;
t357 = t475 * t364 + t470 * t365;
t411 = -mrSges(6,1) * t430 + mrSges(6,2) * t431;
t418 = -mrSges(6,2) * t462 + mrSges(6,3) * t430;
t354 = m(6) * t372 + mrSges(6,1) * t450 - mrSges(6,3) * t402 - t411 * t431 + t418 * t462 + t357;
t419 = mrSges(6,1) * t462 - mrSges(6,3) * t431;
t492 = -t364 * t470 + t475 * t365;
t355 = m(6) * t373 - mrSges(6,2) * t450 + mrSges(6,3) * t401 + t411 * t430 - t419 * t462 + t492;
t350 = t476 * t354 + t471 * t355;
t432 = -mrSges(5,1) * t453 + mrSges(5,2) * t454;
t435 = -mrSges(5,2) * t463 + mrSges(5,3) * t453;
t348 = m(5) * t395 + mrSges(5,1) * t452 - mrSges(5,3) * t428 - t432 * t454 + t435 * t463 + t350;
t436 = mrSges(5,1) * t463 - mrSges(5,3) * t454;
t493 = -t354 * t471 + t476 * t355;
t349 = m(5) * t396 - mrSges(5,2) * t452 + mrSges(5,3) * t427 + t432 * t453 - t436 * t463 + t493;
t342 = t477 * t348 + t472 * t349;
t494 = -t348 * t472 + t477 * t349;
t433 = g(3) * t473 + t440 * t478;
t416 = -qJDD(3) * pkin(3) - pkin(8) * t480 + t456 * t499 - t433;
t394 = -pkin(4) * t427 - pkin(9) * t451 + t454 * t438 + t416;
t375 = -pkin(5) * t401 - pkin(10) * t429 + t420 * t431 + t394;
t490 = m(7) * t375 - t380 * mrSges(7,1) + t381 * mrSges(7,2) - t409 * t403 + t410 * t404;
t455 = (mrSges(4,1) * t473 + mrSges(4,2) * t478) * qJD(1);
t460 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t466;
t461 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t485 = m(6) * t394 - t401 * mrSges(6,1) + t402 * mrSges(6,2) - t430 * t418 + t431 * t419 + t490;
t483 = -m(5) * t416 + t427 * mrSges(5,1) - t428 * mrSges(5,2) + t453 * t435 - t454 * t436 - t485;
t489 = (m(4) * t434 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t457 - qJD(3) * t461 - t455 * t466 + t494) * t473 + (m(4) * t433 + qJDD(3) * mrSges(4,1) - t458 * mrSges(4,3) + qJD(3) * t460 - t455 * t499 + t483) * t478;
t389 = Ifges(7,4) * t410 + Ifges(7,2) * t409 + Ifges(7,6) * t459;
t390 = Ifges(7,1) * t410 + Ifges(7,4) * t409 + Ifges(7,5) * t459;
t486 = -mrSges(7,1) * t367 + mrSges(7,2) * t368 - Ifges(7,5) * t381 - Ifges(7,6) * t380 - Ifges(7,3) * t446 - t410 * t389 + t409 * t390;
t406 = Ifges(6,4) * t431 + Ifges(6,2) * t430 + Ifges(6,6) * t462;
t407 = Ifges(6,1) * t431 + Ifges(6,4) * t430 + Ifges(6,5) * t462;
t484 = -mrSges(6,1) * t372 + mrSges(6,2) * t373 - Ifges(6,5) * t402 - Ifges(6,6) * t401 - Ifges(6,3) * t450 - pkin(5) * t357 - t431 * t406 + t430 * t407 + t486;
t422 = Ifges(5,4) * t454 + Ifges(5,2) * t453 + Ifges(5,6) * t463;
t423 = Ifges(5,1) * t454 + Ifges(5,4) * t453 + Ifges(5,5) * t463;
t482 = mrSges(5,1) * t395 - mrSges(5,2) * t396 + Ifges(5,5) * t428 + Ifges(5,6) * t427 + Ifges(5,3) * t452 + pkin(4) * t350 + t454 * t422 - t453 * t423 - t484;
t449 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t478 - Ifges(4,4) * t473) * qJD(1);
t448 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t478 - Ifges(4,2) * t473) * qJD(1);
t443 = -qJDD(1) * pkin(1) + t487;
t441 = t481 * pkin(1) + t501;
t421 = Ifges(5,5) * t454 + Ifges(5,6) * t453 + Ifges(5,3) * t463;
t405 = Ifges(6,5) * t431 + Ifges(6,6) * t430 + Ifges(6,3) * t462;
t388 = Ifges(7,5) * t410 + Ifges(7,6) * t409 + Ifges(7,3) * t459;
t359 = mrSges(7,2) * t375 - mrSges(7,3) * t367 + Ifges(7,1) * t381 + Ifges(7,4) * t380 + Ifges(7,5) * t446 + t388 * t409 - t389 * t459;
t358 = -mrSges(7,1) * t375 + mrSges(7,3) * t368 + Ifges(7,4) * t381 + Ifges(7,2) * t380 + Ifges(7,6) * t446 - t388 * t410 + t390 * t459;
t344 = mrSges(6,2) * t394 - mrSges(6,3) * t372 + Ifges(6,1) * t402 + Ifges(6,4) * t401 + Ifges(6,5) * t450 - pkin(10) * t357 - t358 * t470 + t359 * t475 + t405 * t430 - t406 * t462;
t343 = -mrSges(6,1) * t394 + mrSges(6,3) * t373 + Ifges(6,4) * t402 + Ifges(6,2) * t401 + Ifges(6,6) * t450 - pkin(5) * t490 + pkin(10) * t492 + t475 * t358 + t470 * t359 - t431 * t405 + t462 * t407;
t340 = m(3) * t443 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t481) + t489;
t339 = mrSges(5,2) * t416 - mrSges(5,3) * t395 + Ifges(5,1) * t428 + Ifges(5,4) * t427 + Ifges(5,5) * t452 - pkin(9) * t350 - t343 * t471 + t344 * t476 + t421 * t453 - t422 * t463;
t338 = -mrSges(5,1) * t416 + mrSges(5,3) * t396 + Ifges(5,4) * t428 + Ifges(5,2) * t427 + Ifges(5,6) * t452 - pkin(4) * t485 + pkin(9) * t493 + t476 * t343 + t471 * t344 - t454 * t421 + t463 * t423;
t1 = [mrSges(2,1) * t495 - mrSges(2,2) * t491 + mrSges(3,2) * t443 - mrSges(3,3) * t441 + t478 * (mrSges(4,2) * t439 - mrSges(4,3) * t433 + Ifges(4,1) * t458 + Ifges(4,4) * t457 + Ifges(4,5) * qJDD(3) - pkin(8) * t342 - qJD(3) * t448 - t338 * t472 + t339 * t477) - t473 * (-mrSges(4,1) * t439 + mrSges(4,3) * t434 + Ifges(4,4) * t458 + Ifges(4,2) * t457 + Ifges(4,6) * qJDD(3) - pkin(3) * t342 + qJD(3) * t449 - t482) - pkin(7) * t489 - pkin(1) * t340 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t441 + m(4) * t439 - t457 * mrSges(4,1) + t481 * mrSges(3,2) + t458 * mrSges(4,2) + t342 + qJDD(1) * mrSges(3,3) + (t460 * t473 + t461 * t478) * qJD(1)) * qJ(2); t340; Ifges(4,5) * t458 + Ifges(4,6) * t457 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t433 - mrSges(4,2) * t434 + t472 * t339 + t477 * t338 + pkin(3) * t483 + pkin(8) * t494 + (t448 * t478 + t449 * t473) * qJD(1); t482; -t484; -t486;];
tauJ  = t1;
