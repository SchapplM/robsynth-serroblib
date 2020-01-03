% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:22
% EndTime: 2019-12-31 19:38:25
% DurationCPUTime: 1.70s
% Computational Cost: add. (8609->259), mult. (19498->320), div. (0->0), fcn. (11518->8), ass. (0->102)
t475 = Ifges(3,1) + Ifges(4,1);
t470 = Ifges(3,4) - Ifges(4,5);
t469 = Ifges(3,5) + Ifges(4,4);
t474 = Ifges(3,2) + Ifges(4,3);
t468 = Ifges(3,6) - Ifges(4,6);
t473 = (Ifges(3,3) + Ifges(4,2));
t472 = 2 * qJD(3);
t471 = mrSges(3,3) + mrSges(4,2);
t445 = cos(qJ(2));
t448 = qJD(1) ^ 2;
t467 = t445 ^ 2 * t448;
t443 = sin(qJ(1));
t446 = cos(qJ(1));
t455 = -g(1) * t446 - g(2) * t443;
t408 = -pkin(1) * t448 + qJDD(1) * pkin(6) + t455;
t442 = sin(qJ(2));
t387 = -g(3) * t442 + t445 * t408;
t412 = (-t445 * pkin(2) - t442 * qJ(3)) * qJD(1);
t447 = qJD(2) ^ 2;
t461 = t445 * qJD(1);
t372 = -pkin(2) * t447 + qJDD(2) * qJ(3) + (qJD(2) * t472) + t412 * t461 + t387;
t460 = qJD(1) * qJD(2);
t458 = t442 * t460;
t416 = qJDD(1) * t445 - t458;
t462 = t442 * qJD(1);
t419 = -(qJD(2) * pkin(3)) - qJ(4) * t462;
t368 = -pkin(3) * t467 - qJ(4) * t416 + qJD(2) * t419 + t372;
t386 = -t445 * g(3) - t442 * t408;
t373 = -qJDD(2) * pkin(2) - qJ(3) * t447 + t412 * t462 + qJDD(3) - t386;
t459 = t445 * t460;
t415 = qJDD(1) * t442 + t459;
t369 = (-t415 + t459) * qJ(4) + (-t442 * t445 * t448 - qJDD(2)) * pkin(3) + t373;
t438 = sin(pkin(8));
t439 = cos(pkin(8));
t400 = (-t445 * t438 + t442 * t439) * qJD(1);
t350 = -0.2e1 * qJD(4) * t400 - t368 * t438 + t439 * t369;
t385 = t415 * t439 - t416 * t438;
t399 = (-t442 * t438 - t445 * t439) * qJD(1);
t348 = (-qJD(2) * t399 - t385) * pkin(7) + (t399 * t400 - qJDD(2)) * pkin(4) + t350;
t351 = 0.2e1 * qJD(4) * t399 + t439 * t368 + t438 * t369;
t384 = -t415 * t438 - t416 * t439;
t390 = -qJD(2) * pkin(4) - pkin(7) * t400;
t398 = t399 ^ 2;
t349 = -pkin(4) * t398 + pkin(7) * t384 + qJD(2) * t390 + t351;
t441 = sin(qJ(5));
t444 = cos(qJ(5));
t346 = t348 * t444 - t349 * t441;
t379 = t399 * t444 - t400 * t441;
t357 = qJD(5) * t379 + t384 * t441 + t385 * t444;
t380 = t399 * t441 + t400 * t444;
t367 = -mrSges(6,1) * t379 + mrSges(6,2) * t380;
t432 = -qJD(2) + qJD(5);
t374 = -mrSges(6,2) * t432 + mrSges(6,3) * t379;
t431 = -qJDD(2) + qJDD(5);
t343 = m(6) * t346 + mrSges(6,1) * t431 - mrSges(6,3) * t357 - t367 * t380 + t374 * t432;
t347 = t348 * t441 + t349 * t444;
t356 = -qJD(5) * t380 + t384 * t444 - t385 * t441;
t375 = mrSges(6,1) * t432 - mrSges(6,3) * t380;
t344 = m(6) * t347 - mrSges(6,2) * t431 + mrSges(6,3) * t356 + t367 * t379 - t375 * t432;
t336 = t444 * t343 + t441 * t344;
t466 = (t473 * qJD(2)) + (t442 * t469 + t445 * t468) * qJD(1);
t465 = -t468 * qJD(2) + (-t442 * t470 - t474 * t445) * qJD(1);
t464 = t469 * qJD(2) + (t442 * t475 + t445 * t470) * qJD(1);
t463 = t443 * g(1) - t446 * g(2);
t382 = -mrSges(5,1) * t399 + mrSges(5,2) * t400;
t388 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t399;
t334 = m(5) * t350 - qJDD(2) * mrSges(5,1) - mrSges(5,3) * t385 - qJD(2) * t388 - t382 * t400 + t336;
t389 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t400;
t456 = -t343 * t441 + t444 * t344;
t335 = m(5) * t351 + qJDD(2) * mrSges(5,2) + mrSges(5,3) * t384 + qJD(2) * t389 + t382 * t399 + t456;
t457 = -t438 * t334 + t439 * t335;
t407 = -qJDD(1) * pkin(1) - t448 * pkin(6) - t463;
t453 = -t416 * pkin(2) + t407 + (-t415 - t459) * qJ(3);
t358 = -pkin(2) * t458 + pkin(3) * t416 - qJ(4) * t467 + qJDD(4) - t453 + (t419 + t472) * t462;
t353 = -pkin(4) * t384 - pkin(7) * t398 + t390 * t400 + t358;
t454 = m(6) * t353 - t356 * mrSges(6,1) + t357 * mrSges(6,2) - t379 * t374 + t380 * t375;
t332 = t439 * t334 + t438 * t335;
t413 = (-t445 * mrSges(4,1) - t442 * mrSges(4,3)) * qJD(1);
t421 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t462;
t452 = m(4) * t372 + qJDD(2) * mrSges(4,3) + qJD(2) * t421 + t413 * t461 + t457;
t423 = mrSges(4,2) * t461 + qJD(2) * mrSges(4,3);
t451 = m(4) * t373 - qJDD(2) * mrSges(4,1) - qJD(2) * t423 + t332;
t345 = m(5) * t358 - t384 * mrSges(5,1) + mrSges(5,2) * t385 - t399 * t388 + t389 * t400 + t454;
t370 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t462 + t453;
t450 = m(4) * t370 - t345;
t360 = Ifges(6,4) * t380 + Ifges(6,2) * t379 + Ifges(6,6) * t432;
t361 = Ifges(6,1) * t380 + Ifges(6,4) * t379 + Ifges(6,5) * t432;
t449 = mrSges(6,1) * t346 - mrSges(6,2) * t347 + Ifges(6,5) * t357 + Ifges(6,6) * t356 + Ifges(6,3) * t431 + t380 * t360 - t379 * t361;
t422 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t461;
t420 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t462;
t414 = (-t445 * mrSges(3,1) + t442 * mrSges(3,2)) * qJD(1);
t378 = Ifges(5,1) * t400 + Ifges(5,4) * t399 - (Ifges(5,5) * qJD(2));
t377 = Ifges(5,4) * t400 + Ifges(5,2) * t399 - (Ifges(5,6) * qJD(2));
t376 = Ifges(5,5) * t400 + Ifges(5,6) * t399 - (Ifges(5,3) * qJD(2));
t359 = Ifges(6,5) * t380 + Ifges(6,6) * t379 + Ifges(6,3) * t432;
t339 = t450 + (-t421 * t442 - t423 * t445) * qJD(1) - mrSges(4,1) * t416 - mrSges(4,3) * t415;
t338 = mrSges(6,2) * t353 - mrSges(6,3) * t346 + Ifges(6,1) * t357 + Ifges(6,4) * t356 + Ifges(6,5) * t431 + t359 * t379 - t360 * t432;
t337 = -mrSges(6,1) * t353 + mrSges(6,3) * t347 + Ifges(6,4) * t357 + Ifges(6,2) * t356 + Ifges(6,6) * t431 - t359 * t380 + t361 * t432;
t331 = t415 * mrSges(4,2) + t413 * t462 + t451;
t330 = mrSges(5,2) * t358 - mrSges(5,3) * t350 + Ifges(5,1) * t385 + Ifges(5,4) * t384 - Ifges(5,5) * qJDD(2) - pkin(7) * t336 + qJD(2) * t377 - t337 * t441 + t338 * t444 + t376 * t399;
t329 = -mrSges(5,1) * t358 + mrSges(5,3) * t351 + Ifges(5,4) * t385 + Ifges(5,2) * t384 - Ifges(5,6) * qJDD(2) - pkin(4) * t454 + pkin(7) * t456 - qJD(2) * t378 + t444 * t337 + t441 * t338 - t400 * t376;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t463 - mrSges(2,2) * t455 + t442 * (mrSges(3,2) * t407 + mrSges(4,2) * t373 - mrSges(3,3) * t386 - mrSges(4,3) * t370 - qJ(3) * t339 - qJ(4) * t332 + t465 * qJD(2) + t469 * qJDD(2) - t438 * t329 + t439 * t330 + t475 * t415 + t470 * t416 + t466 * t461) + t445 * (-mrSges(3,1) * t407 - mrSges(4,1) * t370 + mrSges(4,2) * t372 + mrSges(3,3) * t387 - pkin(2) * t339 + pkin(3) * t345 - qJ(4) * t457 + t464 * qJD(2) + t468 * qJDD(2) - t439 * t329 - t438 * t330 + t470 * t415 + t474 * t416 - t466 * t462) + pkin(1) * (-t450 + (-mrSges(3,2) + mrSges(4,3)) * t415 + ((t422 + t423) * t445 + (-t420 + t421) * t442) * qJD(1) - m(3) * t407 + (mrSges(3,1) + mrSges(4,1)) * t416) + pkin(6) * (t445 * (m(3) * t387 - qJDD(2) * mrSges(3,2) - qJD(2) * t420 + t414 * t461 + t471 * t416 + t452) + (-m(3) * t386 - qJDD(2) * mrSges(3,1) - qJD(2) * t422 + t471 * t415 + (t413 + t414) * t462 + t451) * t442); qJ(3) * t452 + (qJ(3) * mrSges(4,2) + t468) * t416 + t469 * t415 + (Ifges(5,3) + t473) * qJDD(2) + t399 * t378 - Ifges(5,6) * t384 - Ifges(5,5) * t385 + mrSges(3,1) * t386 - mrSges(3,2) * t387 + mrSges(4,3) * t372 - mrSges(4,1) * t373 + mrSges(5,2) * t351 - mrSges(5,1) * t350 - pkin(4) * t336 - t449 - pkin(3) * t332 - pkin(2) * t331 - t400 * t377 + (-t465 * t442 - t464 * t445) * qJD(1); t331; t345; t449;];
tauJ = t1;
