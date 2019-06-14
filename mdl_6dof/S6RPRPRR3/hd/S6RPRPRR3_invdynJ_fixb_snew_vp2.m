% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR3
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
% Datum: 2019-05-05 18:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:33:39
% EndTime: 2019-05-05 18:33:43
% DurationCPUTime: 3.86s
% Computational Cost: add. (44889->292), mult. (93890->370), div. (0->0), fcn. (64024->12), ass. (0->119)
t462 = sin(qJ(1));
t466 = cos(qJ(1));
t479 = t462 * g(1) - g(2) * t466;
t441 = qJDD(1) * pkin(1) + t479;
t468 = qJD(1) ^ 2;
t474 = -g(1) * t466 - g(2) * t462;
t444 = -pkin(1) * t468 + t474;
t456 = sin(pkin(10));
t458 = cos(pkin(10));
t415 = t458 * t441 - t456 * t444;
t405 = -qJDD(1) * pkin(2) - t468 * pkin(7) - t415;
t461 = sin(qJ(3));
t465 = cos(qJ(3));
t481 = qJD(1) * qJD(3);
t480 = t465 * t481;
t445 = qJDD(1) * t461 + t480;
t452 = t461 * t481;
t446 = qJDD(1) * t465 - t452;
t391 = (-t445 - t480) * qJ(4) + (-t446 + t452) * pkin(3) + t405;
t416 = t456 * t441 + t458 * t444;
t406 = -pkin(2) * t468 + qJDD(1) * pkin(7) + t416;
t454 = -g(3) + qJDD(2);
t399 = t465 * t406 + t461 * t454;
t442 = (-pkin(3) * t465 - qJ(4) * t461) * qJD(1);
t467 = qJD(3) ^ 2;
t482 = qJD(1) * t465;
t397 = -pkin(3) * t467 + qJDD(3) * qJ(4) + t442 * t482 + t399;
t455 = sin(pkin(11));
t457 = cos(pkin(11));
t483 = qJD(1) * t461;
t438 = qJD(3) * t455 + t457 * t483;
t371 = -0.2e1 * qJD(4) * t438 + t457 * t391 - t455 * t397;
t422 = qJDD(3) * t455 + t445 * t457;
t437 = qJD(3) * t457 - t455 * t483;
t365 = (-t437 * t482 - t422) * pkin(8) + (t437 * t438 - t446) * pkin(4) + t371;
t372 = 0.2e1 * qJD(4) * t437 + t455 * t391 + t457 * t397;
t421 = qJDD(3) * t457 - t445 * t455;
t423 = -pkin(4) * t482 - pkin(8) * t438;
t436 = t437 ^ 2;
t370 = -pkin(4) * t436 + pkin(8) * t421 + t423 * t482 + t372;
t460 = sin(qJ(5));
t464 = cos(qJ(5));
t355 = t464 * t365 - t460 * t370;
t413 = t437 * t464 - t438 * t460;
t384 = qJD(5) * t413 + t421 * t460 + t422 * t464;
t414 = t437 * t460 + t438 * t464;
t440 = qJDD(5) - t446;
t450 = qJD(5) - t482;
t353 = (t413 * t450 - t384) * pkin(9) + (t413 * t414 + t440) * pkin(5) + t355;
t356 = t460 * t365 + t464 * t370;
t383 = -qJD(5) * t414 + t421 * t464 - t422 * t460;
t402 = pkin(5) * t450 - pkin(9) * t414;
t412 = t413 ^ 2;
t354 = -pkin(5) * t412 + pkin(9) * t383 - t402 * t450 + t356;
t459 = sin(qJ(6));
t463 = cos(qJ(6));
t351 = t353 * t463 - t354 * t459;
t393 = t413 * t463 - t414 * t459;
t367 = qJD(6) * t393 + t383 * t459 + t384 * t463;
t394 = t413 * t459 + t414 * t463;
t378 = -mrSges(7,1) * t393 + mrSges(7,2) * t394;
t449 = qJD(6) + t450;
t380 = -mrSges(7,2) * t449 + mrSges(7,3) * t393;
t434 = qJDD(6) + t440;
t346 = m(7) * t351 + mrSges(7,1) * t434 - mrSges(7,3) * t367 - t378 * t394 + t380 * t449;
t352 = t353 * t459 + t354 * t463;
t366 = -qJD(6) * t394 + t383 * t463 - t384 * t459;
t381 = mrSges(7,1) * t449 - mrSges(7,3) * t394;
t347 = m(7) * t352 - mrSges(7,2) * t434 + mrSges(7,3) * t366 + t378 * t393 - t381 * t449;
t340 = t463 * t346 + t459 * t347;
t396 = -mrSges(6,1) * t413 + mrSges(6,2) * t414;
t400 = -mrSges(6,2) * t450 + mrSges(6,3) * t413;
t338 = m(6) * t355 + mrSges(6,1) * t440 - mrSges(6,3) * t384 - t396 * t414 + t400 * t450 + t340;
t401 = mrSges(6,1) * t450 - mrSges(6,3) * t414;
t475 = -t346 * t459 + t463 * t347;
t339 = m(6) * t356 - mrSges(6,2) * t440 + mrSges(6,3) * t383 + t396 * t413 - t401 * t450 + t475;
t334 = t464 * t338 + t460 * t339;
t443 = (-mrSges(4,1) * t465 + mrSges(4,2) * t461) * qJD(1);
t447 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t483;
t417 = -mrSges(5,1) * t437 + mrSges(5,2) * t438;
t419 = mrSges(5,2) * t482 + mrSges(5,3) * t437;
t332 = m(5) * t371 - mrSges(5,1) * t446 - mrSges(5,3) * t422 - t417 * t438 - t419 * t482 + t334;
t420 = -mrSges(5,1) * t482 - mrSges(5,3) * t438;
t476 = -t338 * t460 + t464 * t339;
t333 = m(5) * t372 + mrSges(5,2) * t446 + mrSges(5,3) * t421 + t417 * t437 + t420 * t482 + t476;
t477 = -t332 * t455 + t457 * t333;
t327 = m(4) * t399 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t446 - qJD(3) * t447 + t443 * t482 + t477;
t398 = -t461 * t406 + t454 * t465;
t395 = -qJDD(3) * pkin(3) - qJ(4) * t467 + t442 * t483 + qJDD(4) - t398;
t379 = -pkin(4) * t421 - pkin(8) * t436 + t438 * t423 + t395;
t358 = -pkin(5) * t383 - pkin(9) * t412 + t402 * t414 + t379;
t473 = m(7) * t358 - t366 * mrSges(7,1) + t367 * mrSges(7,2) - t393 * t380 + t394 * t381;
t471 = m(6) * t379 - t383 * mrSges(6,1) + t384 * mrSges(6,2) - t413 * t400 + t414 * t401 + t473;
t349 = m(5) * t395 - t421 * mrSges(5,1) + t422 * mrSges(5,2) - t437 * t419 + t438 * t420 + t471;
t448 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t482;
t348 = m(4) * t398 + qJDD(3) * mrSges(4,1) - t445 * mrSges(4,3) + qJD(3) * t448 - t443 * t483 - t349;
t478 = t465 * t327 - t348 * t461;
t328 = t332 * t457 + t333 * t455;
t374 = Ifges(7,4) * t394 + Ifges(7,2) * t393 + Ifges(7,6) * t449;
t375 = Ifges(7,1) * t394 + Ifges(7,4) * t393 + Ifges(7,5) * t449;
t472 = -mrSges(7,1) * t351 + mrSges(7,2) * t352 - Ifges(7,5) * t367 - Ifges(7,6) * t366 - Ifges(7,3) * t434 - t394 * t374 + t393 * t375;
t470 = -m(4) * t405 + t446 * mrSges(4,1) - mrSges(4,2) * t445 - t447 * t483 + t448 * t482 - t328;
t388 = Ifges(6,4) * t414 + Ifges(6,2) * t413 + Ifges(6,6) * t450;
t389 = Ifges(6,1) * t414 + Ifges(6,4) * t413 + Ifges(6,5) * t450;
t469 = mrSges(6,1) * t355 - mrSges(6,2) * t356 + Ifges(6,5) * t384 + Ifges(6,6) * t383 + Ifges(6,3) * t440 + pkin(5) * t340 + t414 * t388 - t413 * t389 - t472;
t433 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t461 + Ifges(4,4) * t465) * qJD(1);
t432 = Ifges(4,6) * qJD(3) + (t461 * Ifges(4,4) + Ifges(4,2) * t465) * qJD(1);
t409 = Ifges(5,1) * t438 + Ifges(5,4) * t437 - Ifges(5,5) * t482;
t408 = Ifges(5,4) * t438 + Ifges(5,2) * t437 - Ifges(5,6) * t482;
t407 = Ifges(5,5) * t438 + Ifges(5,6) * t437 - Ifges(5,3) * t482;
t387 = Ifges(6,5) * t414 + Ifges(6,6) * t413 + Ifges(6,3) * t450;
t373 = Ifges(7,5) * t394 + Ifges(7,6) * t393 + Ifges(7,3) * t449;
t342 = mrSges(7,2) * t358 - mrSges(7,3) * t351 + Ifges(7,1) * t367 + Ifges(7,4) * t366 + Ifges(7,5) * t434 + t373 * t393 - t374 * t449;
t341 = -mrSges(7,1) * t358 + mrSges(7,3) * t352 + Ifges(7,4) * t367 + Ifges(7,2) * t366 + Ifges(7,6) * t434 - t373 * t394 + t375 * t449;
t330 = mrSges(6,2) * t379 - mrSges(6,3) * t355 + Ifges(6,1) * t384 + Ifges(6,4) * t383 + Ifges(6,5) * t440 - pkin(9) * t340 - t341 * t459 + t342 * t463 + t387 * t413 - t388 * t450;
t329 = -mrSges(6,1) * t379 + mrSges(6,3) * t356 + Ifges(6,4) * t384 + Ifges(6,2) * t383 + Ifges(6,6) * t440 - pkin(5) * t473 + pkin(9) * t475 + t463 * t341 + t459 * t342 - t414 * t387 + t450 * t389;
t325 = mrSges(5,2) * t395 - mrSges(5,3) * t371 + Ifges(5,1) * t422 + Ifges(5,4) * t421 - Ifges(5,5) * t446 - pkin(8) * t334 - t329 * t460 + t330 * t464 + t407 * t437 + t408 * t482;
t324 = -mrSges(5,1) * t395 + mrSges(5,3) * t372 + Ifges(5,4) * t422 + Ifges(5,2) * t421 - Ifges(5,6) * t446 - pkin(4) * t471 + pkin(8) * t476 + t464 * t329 + t460 * t330 - t438 * t407 - t409 * t482;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t479 - mrSges(2,2) * t474 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t415 - mrSges(3,2) * t416 + t461 * (mrSges(4,2) * t405 - mrSges(4,3) * t398 + Ifges(4,1) * t445 + Ifges(4,4) * t446 + Ifges(4,5) * qJDD(3) - qJ(4) * t328 - qJD(3) * t432 - t455 * t324 + t457 * t325) + t465 * (-mrSges(4,1) * t405 - mrSges(5,1) * t371 + mrSges(5,2) * t372 + mrSges(4,3) * t399 + Ifges(4,4) * t445 - Ifges(5,5) * t422 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t421 - pkin(3) * t328 - pkin(4) * t334 + qJD(3) * t433 - t438 * t408 + t437 * t409 - t469 + (Ifges(4,2) + Ifges(5,3)) * t446) + pkin(2) * t470 + pkin(7) * t478 + pkin(1) * (t456 * (m(3) * t416 - mrSges(3,1) * t468 - qJDD(1) * mrSges(3,2) + t478) + t458 * (m(3) * t415 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t468 + t470)); m(3) * t454 + t327 * t461 + t348 * t465; Ifges(4,5) * t445 + Ifges(4,6) * t446 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t398 - mrSges(4,2) * t399 + t455 * t325 + t457 * t324 - pkin(3) * t349 + qJ(4) * t477 + (t432 * t461 - t433 * t465) * qJD(1); t349; t469; -t472;];
tauJ  = t1;
