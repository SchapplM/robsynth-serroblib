% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:04:48
% EndTime: 2019-05-05 22:04:52
% DurationCPUTime: 4.25s
% Computational Cost: add. (48986->292), mult. (98726->370), div. (0->0), fcn. (67199->12), ass. (0->119)
t462 = sin(qJ(1));
t466 = cos(qJ(1));
t479 = t462 * g(1) - g(2) * t466;
t441 = qJDD(1) * pkin(1) + t479;
t468 = qJD(1) ^ 2;
t474 = -g(1) * t466 - g(2) * t462;
t443 = -pkin(1) * t468 + t474;
t456 = sin(pkin(10));
t458 = cos(pkin(10));
t420 = t458 * t441 - t456 * t443;
t405 = -qJDD(1) * pkin(2) - t468 * pkin(7) - t420;
t461 = sin(qJ(3));
t465 = cos(qJ(3));
t481 = qJD(1) * qJD(3);
t480 = t465 * t481;
t445 = qJDD(1) * t461 + t480;
t452 = t461 * t481;
t446 = qJDD(1) * t465 - t452;
t388 = (-t445 - t480) * pkin(8) + (-t446 + t452) * pkin(3) + t405;
t421 = t456 * t441 + t458 * t443;
t406 = -pkin(2) * t468 + qJDD(1) * pkin(7) + t421;
t454 = -g(3) + qJDD(2);
t399 = t465 * t406 + t461 * t454;
t444 = (-pkin(3) * t465 - pkin(8) * t461) * qJD(1);
t467 = qJD(3) ^ 2;
t482 = qJD(1) * t465;
t396 = -pkin(3) * t467 + qJDD(3) * pkin(8) + t444 * t482 + t399;
t460 = sin(qJ(4));
t464 = cos(qJ(4));
t377 = t464 * t388 - t460 * t396;
t483 = qJD(1) * t461;
t439 = qJD(3) * t464 - t460 * t483;
t416 = qJD(4) * t439 + qJDD(3) * t460 + t445 * t464;
t438 = qJDD(4) - t446;
t440 = qJD(3) * t460 + t464 * t483;
t450 = qJD(4) - t482;
t361 = (t439 * t450 - t416) * qJ(5) + (t439 * t440 + t438) * pkin(4) + t377;
t378 = t460 * t388 + t464 * t396;
t415 = -qJD(4) * t440 + qJDD(3) * t464 - t445 * t460;
t424 = pkin(4) * t450 - qJ(5) * t440;
t437 = t439 ^ 2;
t369 = -pkin(4) * t437 + qJ(5) * t415 - t424 * t450 + t378;
t455 = sin(pkin(11));
t457 = cos(pkin(11));
t419 = t439 * t455 + t440 * t457;
t355 = -0.2e1 * qJD(5) * t419 + t457 * t361 - t455 * t369;
t390 = t415 * t455 + t416 * t457;
t418 = t439 * t457 - t440 * t455;
t353 = (t418 * t450 - t390) * pkin(9) + (t418 * t419 + t438) * pkin(5) + t355;
t356 = 0.2e1 * qJD(5) * t418 + t455 * t361 + t457 * t369;
t389 = t415 * t457 - t416 * t455;
t402 = pkin(5) * t450 - pkin(9) * t419;
t417 = t418 ^ 2;
t354 = -pkin(5) * t417 + pkin(9) * t389 - t402 * t450 + t356;
t459 = sin(qJ(6));
t463 = cos(qJ(6));
t351 = t353 * t463 - t354 * t459;
t393 = t418 * t463 - t419 * t459;
t368 = qJD(6) * t393 + t389 * t459 + t390 * t463;
t394 = t418 * t459 + t419 * t463;
t379 = -mrSges(7,1) * t393 + mrSges(7,2) * t394;
t449 = qJD(6) + t450;
t380 = -mrSges(7,2) * t449 + mrSges(7,3) * t393;
t434 = qJDD(6) + t438;
t346 = m(7) * t351 + mrSges(7,1) * t434 - mrSges(7,3) * t368 - t379 * t394 + t380 * t449;
t352 = t353 * t459 + t354 * t463;
t367 = -qJD(6) * t394 + t389 * t463 - t390 * t459;
t381 = mrSges(7,1) * t449 - mrSges(7,3) * t394;
t347 = m(7) * t352 - mrSges(7,2) * t434 + mrSges(7,3) * t367 + t379 * t393 - t381 * t449;
t340 = t463 * t346 + t459 * t347;
t397 = -mrSges(6,1) * t418 + mrSges(6,2) * t419;
t400 = -mrSges(6,2) * t450 + mrSges(6,3) * t418;
t338 = m(6) * t355 + mrSges(6,1) * t438 - mrSges(6,3) * t390 - t397 * t419 + t400 * t450 + t340;
t401 = mrSges(6,1) * t450 - mrSges(6,3) * t419;
t475 = -t346 * t459 + t463 * t347;
t339 = m(6) * t356 - mrSges(6,2) * t438 + mrSges(6,3) * t389 + t397 * t418 - t401 * t450 + t475;
t334 = t457 * t338 + t455 * t339;
t386 = Ifges(6,4) * t419 + Ifges(6,2) * t418 + Ifges(6,6) * t450;
t387 = Ifges(6,1) * t419 + Ifges(6,4) * t418 + Ifges(6,5) * t450;
t408 = Ifges(5,4) * t440 + Ifges(5,2) * t439 + Ifges(5,6) * t450;
t409 = Ifges(5,1) * t440 + Ifges(5,4) * t439 + Ifges(5,5) * t450;
t372 = Ifges(7,4) * t394 + Ifges(7,2) * t393 + Ifges(7,6) * t449;
t373 = Ifges(7,1) * t394 + Ifges(7,4) * t393 + Ifges(7,5) * t449;
t472 = -mrSges(7,1) * t351 + mrSges(7,2) * t352 - Ifges(7,5) * t368 - Ifges(7,6) * t367 - Ifges(7,3) * t434 - t394 * t372 + t393 * t373;
t485 = mrSges(5,1) * t377 + mrSges(6,1) * t355 - mrSges(5,2) * t378 - mrSges(6,2) * t356 + Ifges(5,5) * t416 + Ifges(6,5) * t390 + Ifges(5,6) * t415 + Ifges(6,6) * t389 + pkin(4) * t334 + pkin(5) * t340 + t419 * t386 - t418 * t387 + t440 * t408 - t439 * t409 + (Ifges(5,3) + Ifges(6,3)) * t438 - t472;
t442 = (-mrSges(4,1) * t465 + mrSges(4,2) * t461) * qJD(1);
t447 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t483;
t422 = -mrSges(5,1) * t439 + mrSges(5,2) * t440;
t423 = -mrSges(5,2) * t450 + mrSges(5,3) * t439;
t332 = m(5) * t377 + mrSges(5,1) * t438 - mrSges(5,3) * t416 - t422 * t440 + t423 * t450 + t334;
t425 = mrSges(5,1) * t450 - mrSges(5,3) * t440;
t476 = -t338 * t455 + t457 * t339;
t333 = m(5) * t378 - mrSges(5,2) * t438 + mrSges(5,3) * t415 + t422 * t439 - t425 * t450 + t476;
t477 = -t332 * t460 + t464 * t333;
t327 = m(4) * t399 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t446 - qJD(3) * t447 + t442 * t482 + t477;
t398 = -t461 * t406 + t454 * t465;
t448 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t482;
t395 = -qJDD(3) * pkin(3) - pkin(8) * t467 + t444 * t483 - t398;
t374 = -pkin(4) * t415 - qJ(5) * t437 + t440 * t424 + qJDD(5) + t395;
t358 = -pkin(5) * t389 - pkin(9) * t417 + t402 * t419 + t374;
t473 = m(7) * t358 - t367 * mrSges(7,1) + t368 * mrSges(7,2) - t393 * t380 + t394 * t381;
t349 = m(6) * t374 - t389 * mrSges(6,1) + t390 * mrSges(6,2) - t418 * t400 + t419 * t401 + t473;
t470 = -m(5) * t395 + t415 * mrSges(5,1) - t416 * mrSges(5,2) + t439 * t423 - t440 * t425 - t349;
t348 = m(4) * t398 + qJDD(3) * mrSges(4,1) - t445 * mrSges(4,3) + qJD(3) * t448 - t442 * t483 + t470;
t478 = t465 * t327 - t348 * t461;
t328 = t332 * t464 + t333 * t460;
t471 = -m(4) * t405 + t446 * mrSges(4,1) - mrSges(4,2) * t445 - t447 * t483 + t448 * t482 - t328;
t433 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t461 + Ifges(4,4) * t465) * qJD(1);
t432 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t461 + Ifges(4,2) * t465) * qJD(1);
t407 = Ifges(5,5) * t440 + Ifges(5,6) * t439 + Ifges(5,3) * t450;
t385 = Ifges(6,5) * t419 + Ifges(6,6) * t418 + Ifges(6,3) * t450;
t371 = Ifges(7,5) * t394 + Ifges(7,6) * t393 + Ifges(7,3) * t449;
t342 = mrSges(7,2) * t358 - mrSges(7,3) * t351 + Ifges(7,1) * t368 + Ifges(7,4) * t367 + Ifges(7,5) * t434 + t371 * t393 - t372 * t449;
t341 = -mrSges(7,1) * t358 + mrSges(7,3) * t352 + Ifges(7,4) * t368 + Ifges(7,2) * t367 + Ifges(7,6) * t434 - t371 * t394 + t373 * t449;
t330 = mrSges(6,2) * t374 - mrSges(6,3) * t355 + Ifges(6,1) * t390 + Ifges(6,4) * t389 + Ifges(6,5) * t438 - pkin(9) * t340 - t341 * t459 + t342 * t463 + t385 * t418 - t386 * t450;
t329 = -mrSges(6,1) * t374 + mrSges(6,3) * t356 + Ifges(6,4) * t390 + Ifges(6,2) * t389 + Ifges(6,6) * t438 - pkin(5) * t473 + pkin(9) * t475 + t463 * t341 + t459 * t342 - t419 * t385 + t450 * t387;
t325 = mrSges(5,2) * t395 - mrSges(5,3) * t377 + Ifges(5,1) * t416 + Ifges(5,4) * t415 + Ifges(5,5) * t438 - qJ(5) * t334 - t329 * t455 + t330 * t457 + t407 * t439 - t408 * t450;
t324 = -mrSges(5,1) * t395 + mrSges(5,3) * t378 + Ifges(5,4) * t416 + Ifges(5,2) * t415 + Ifges(5,6) * t438 - pkin(4) * t349 + qJ(5) * t476 + t457 * t329 + t455 * t330 - t440 * t407 + t450 * t409;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t479 - mrSges(2,2) * t474 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t420 - mrSges(3,2) * t421 + t461 * (mrSges(4,2) * t405 - mrSges(4,3) * t398 + Ifges(4,1) * t445 + Ifges(4,4) * t446 + Ifges(4,5) * qJDD(3) - pkin(8) * t328 - qJD(3) * t432 - t460 * t324 + t464 * t325) + t465 * (-mrSges(4,1) * t405 + mrSges(4,3) * t399 + Ifges(4,4) * t445 + Ifges(4,2) * t446 + Ifges(4,6) * qJDD(3) - pkin(3) * t328 + qJD(3) * t433 - t485) + pkin(2) * t471 + pkin(7) * t478 + pkin(1) * (t456 * (m(3) * t421 - mrSges(3,1) * t468 - qJDD(1) * mrSges(3,2) + t478) + t458 * (m(3) * t420 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t468 + t471)); m(3) * t454 + t327 * t461 + t348 * t465; Ifges(4,5) * t445 + Ifges(4,6) * t446 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t398 - mrSges(4,2) * t399 + t460 * t325 + t464 * t324 + pkin(3) * t470 + pkin(8) * t477 + (t432 * t461 - t433 * t465) * qJD(1); t485; t349; -t472;];
tauJ  = t1;
