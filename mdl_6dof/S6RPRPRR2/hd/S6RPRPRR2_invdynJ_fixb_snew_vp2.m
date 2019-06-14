% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:22:58
% EndTime: 2019-05-05 18:23:01
% DurationCPUTime: 3.38s
% Computational Cost: add. (35027->293), mult. (76061->374), div. (0->0), fcn. (51757->12), ass. (0->119)
t476 = sin(qJ(1));
t480 = cos(qJ(1));
t495 = t476 * g(1) - g(2) * t480;
t455 = qJDD(1) * pkin(1) + t495;
t482 = qJD(1) ^ 2;
t489 = -g(1) * t480 - g(2) * t476;
t457 = -pkin(1) * t482 + t489;
t470 = sin(pkin(10));
t472 = cos(pkin(10));
t431 = t470 * t455 + t472 * t457;
t422 = -pkin(2) * t482 + qJDD(1) * pkin(7) + t431;
t468 = -g(3) + qJDD(2);
t475 = sin(qJ(3));
t479 = cos(qJ(3));
t412 = -t422 * t475 + t479 * t468;
t497 = qJD(1) * qJD(3);
t496 = t479 * t497;
t458 = qJDD(1) * t475 + t496;
t396 = (-t458 + t496) * qJ(4) + (t475 * t479 * t482 + qJDD(3)) * pkin(3) + t412;
t413 = t479 * t422 + t475 * t468;
t459 = qJDD(1) * t479 - t475 * t497;
t498 = t475 * qJD(1);
t460 = qJD(3) * pkin(3) - qJ(4) * t498;
t467 = t479 ^ 2;
t399 = -pkin(3) * t467 * t482 + qJ(4) * t459 - qJD(3) * t460 + t413;
t469 = sin(pkin(11));
t471 = cos(pkin(11));
t444 = (t479 * t469 + t475 * t471) * qJD(1);
t383 = -0.2e1 * qJD(4) * t444 + t396 * t471 - t469 * t399;
t500 = qJD(1) * t479;
t443 = -t469 * t498 + t471 * t500;
t384 = 0.2e1 * qJD(4) * t443 + t469 * t396 + t471 * t399;
t424 = -mrSges(5,1) * t443 + mrSges(5,2) * t444;
t432 = -t469 * t458 + t459 * t471;
t438 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t444;
t426 = -pkin(4) * t443 - pkin(8) * t444;
t481 = qJD(3) ^ 2;
t376 = -pkin(4) * t481 + qJDD(3) * pkin(8) + t426 * t443 + t384;
t430 = t455 * t472 - t470 * t457;
t488 = -qJDD(1) * pkin(2) - t430;
t400 = -pkin(3) * t459 + qJDD(4) + t460 * t498 + (-qJ(4) * t467 - pkin(7)) * t482 + t488;
t433 = t458 * t471 + t459 * t469;
t387 = (-qJD(3) * t443 - t433) * pkin(8) + (qJD(3) * t444 - t432) * pkin(4) + t400;
t474 = sin(qJ(5));
t478 = cos(qJ(5));
t371 = -t376 * t474 + t478 * t387;
t435 = qJD(3) * t478 - t444 * t474;
t407 = qJD(5) * t435 + qJDD(3) * t474 + t433 * t478;
t429 = qJDD(5) - t432;
t436 = qJD(3) * t474 + t444 * t478;
t442 = qJD(5) - t443;
t369 = (t435 * t442 - t407) * pkin(9) + (t435 * t436 + t429) * pkin(5) + t371;
t372 = t478 * t376 + t474 * t387;
t406 = -qJD(5) * t436 + qJDD(3) * t478 - t433 * t474;
t416 = pkin(5) * t442 - pkin(9) * t436;
t434 = t435 ^ 2;
t370 = -pkin(5) * t434 + pkin(9) * t406 - t416 * t442 + t372;
t473 = sin(qJ(6));
t477 = cos(qJ(6));
t367 = t369 * t477 - t370 * t473;
t408 = t435 * t477 - t436 * t473;
t381 = qJD(6) * t408 + t406 * t473 + t407 * t477;
t409 = t435 * t473 + t436 * t477;
t392 = -mrSges(7,1) * t408 + mrSges(7,2) * t409;
t439 = qJD(6) + t442;
t397 = -mrSges(7,2) * t439 + mrSges(7,3) * t408;
t427 = qJDD(6) + t429;
t364 = m(7) * t367 + mrSges(7,1) * t427 - t381 * mrSges(7,3) - t392 * t409 + t397 * t439;
t368 = t369 * t473 + t370 * t477;
t380 = -qJD(6) * t409 + t406 * t477 - t407 * t473;
t398 = mrSges(7,1) * t439 - mrSges(7,3) * t409;
t365 = m(7) * t368 - mrSges(7,2) * t427 + t380 * mrSges(7,3) + t392 * t408 - t398 * t439;
t356 = t477 * t364 + t473 * t365;
t410 = -mrSges(6,1) * t435 + mrSges(6,2) * t436;
t414 = -mrSges(6,2) * t442 + mrSges(6,3) * t435;
t354 = m(6) * t371 + mrSges(6,1) * t429 - mrSges(6,3) * t407 - t410 * t436 + t414 * t442 + t356;
t415 = mrSges(6,1) * t442 - mrSges(6,3) * t436;
t491 = -t473 * t364 + t477 * t365;
t355 = m(6) * t372 - mrSges(6,2) * t429 + mrSges(6,3) * t406 + t410 * t435 - t415 * t442 + t491;
t492 = -t474 * t354 + t478 * t355;
t348 = m(5) * t384 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t432 - qJD(3) * t438 + t424 * t443 + t492;
t437 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t443;
t375 = -qJDD(3) * pkin(4) - pkin(8) * t481 + t444 * t426 - t383;
t373 = -pkin(5) * t406 - pkin(9) * t434 + t416 * t436 + t375;
t487 = m(7) * t373 - t380 * mrSges(7,1) + t381 * mrSges(7,2) - t408 * t397 + t398 * t409;
t485 = -m(6) * t375 + t406 * mrSges(6,1) - mrSges(6,2) * t407 + t435 * t414 - t415 * t436 - t487;
t360 = m(5) * t383 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t433 + qJD(3) * t437 - t424 * t444 + t485;
t343 = t469 * t348 + t471 * t360;
t350 = t478 * t354 + t474 * t355;
t456 = (-mrSges(4,1) * t479 + mrSges(4,2) * t475) * qJD(1);
t462 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t500;
t341 = m(4) * t412 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t458 + qJD(3) * t462 - t456 * t498 + t343;
t461 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t498;
t493 = t471 * t348 - t360 * t469;
t342 = m(4) * t413 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t459 - qJD(3) * t461 + t456 * t500 + t493;
t494 = -t475 * t341 + t479 * t342;
t389 = Ifges(7,4) * t409 + Ifges(7,2) * t408 + Ifges(7,6) * t439;
t390 = Ifges(7,1) * t409 + Ifges(7,4) * t408 + Ifges(7,5) * t439;
t486 = -mrSges(7,1) * t367 + mrSges(7,2) * t368 - Ifges(7,5) * t381 - Ifges(7,6) * t380 - Ifges(7,3) * t427 - t409 * t389 + t408 * t390;
t349 = m(5) * t400 - t432 * mrSges(5,1) + mrSges(5,2) * t433 - t443 * t437 + t438 * t444 + t350;
t421 = -pkin(7) * t482 + t488;
t484 = -m(4) * t421 + t459 * mrSges(4,1) - mrSges(4,2) * t458 - t461 * t498 + t462 * t500 - t349;
t402 = Ifges(6,4) * t436 + Ifges(6,2) * t435 + Ifges(6,6) * t442;
t403 = Ifges(6,1) * t436 + Ifges(6,4) * t435 + Ifges(6,5) * t442;
t483 = mrSges(6,1) * t371 - mrSges(6,2) * t372 + Ifges(6,5) * t407 + Ifges(6,6) * t406 + Ifges(6,3) * t429 + pkin(5) * t356 + t436 * t402 - t435 * t403 - t486;
t450 = Ifges(4,5) * qJD(3) + (t475 * Ifges(4,1) + t479 * Ifges(4,4)) * qJD(1);
t449 = Ifges(4,6) * qJD(3) + (t475 * Ifges(4,4) + t479 * Ifges(4,2)) * qJD(1);
t420 = Ifges(5,1) * t444 + Ifges(5,4) * t443 + Ifges(5,5) * qJD(3);
t419 = Ifges(5,4) * t444 + Ifges(5,2) * t443 + Ifges(5,6) * qJD(3);
t418 = Ifges(5,5) * t444 + Ifges(5,6) * t443 + Ifges(5,3) * qJD(3);
t401 = Ifges(6,5) * t436 + Ifges(6,6) * t435 + Ifges(6,3) * t442;
t388 = Ifges(7,5) * t409 + Ifges(7,6) * t408 + Ifges(7,3) * t439;
t358 = mrSges(7,2) * t373 - mrSges(7,3) * t367 + Ifges(7,1) * t381 + Ifges(7,4) * t380 + Ifges(7,5) * t427 + t388 * t408 - t389 * t439;
t357 = -mrSges(7,1) * t373 + mrSges(7,3) * t368 + Ifges(7,4) * t381 + Ifges(7,2) * t380 + Ifges(7,6) * t427 - t388 * t409 + t390 * t439;
t345 = mrSges(6,2) * t375 - mrSges(6,3) * t371 + Ifges(6,1) * t407 + Ifges(6,4) * t406 + Ifges(6,5) * t429 - pkin(9) * t356 - t357 * t473 + t358 * t477 + t401 * t435 - t402 * t442;
t344 = -mrSges(6,1) * t375 + mrSges(6,3) * t372 + Ifges(6,4) * t407 + Ifges(6,2) * t406 + Ifges(6,6) * t429 - pkin(5) * t487 + pkin(9) * t491 + t477 * t357 + t473 * t358 - t436 * t401 + t442 * t403;
t339 = -mrSges(5,1) * t400 + mrSges(5,3) * t384 + Ifges(5,4) * t433 + Ifges(5,2) * t432 + Ifges(5,6) * qJDD(3) - pkin(4) * t350 + qJD(3) * t420 - t444 * t418 - t483;
t338 = mrSges(5,2) * t400 - mrSges(5,3) * t383 + Ifges(5,1) * t433 + Ifges(5,4) * t432 + Ifges(5,5) * qJDD(3) - pkin(8) * t350 - qJD(3) * t419 - t344 * t474 + t345 * t478 + t418 * t443;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t495 - mrSges(2,2) * t489 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t430 - mrSges(3,2) * t431 + t475 * (mrSges(4,2) * t421 - mrSges(4,3) * t412 + Ifges(4,1) * t458 + Ifges(4,4) * t459 + Ifges(4,5) * qJDD(3) - qJ(4) * t343 - qJD(3) * t449 + t338 * t471 - t339 * t469) + t479 * (-mrSges(4,1) * t421 + mrSges(4,3) * t413 + Ifges(4,4) * t458 + Ifges(4,2) * t459 + Ifges(4,6) * qJDD(3) - pkin(3) * t349 + qJ(4) * t493 + qJD(3) * t450 + t469 * t338 + t471 * t339) + pkin(2) * t484 + pkin(7) * t494 + pkin(1) * (t470 * (m(3) * t431 - mrSges(3,1) * t482 - qJDD(1) * mrSges(3,2) + t494) + t472 * (m(3) * t430 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t482 + t484)); m(3) * t468 + t341 * t479 + t342 * t475; Ifges(4,5) * t458 + Ifges(4,6) * t459 + mrSges(4,1) * t412 - mrSges(4,2) * t413 + Ifges(5,5) * t433 + Ifges(5,6) * t432 + t444 * t419 - t443 * t420 + mrSges(5,1) * t383 - mrSges(5,2) * t384 + t474 * t345 + t478 * t344 + pkin(4) * t485 + pkin(8) * t492 + pkin(3) * t343 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t475 * t449 - t479 * t450) * qJD(1); t349; t483; -t486;];
tauJ  = t1;
