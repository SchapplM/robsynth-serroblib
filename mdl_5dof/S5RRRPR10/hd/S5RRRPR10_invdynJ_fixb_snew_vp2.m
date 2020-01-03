% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:27:32
% EndTime: 2019-12-31 21:27:38
% DurationCPUTime: 4.88s
% Computational Cost: add. (51654->288), mult. (114347->382), div. (0->0), fcn. (89502->12), ass. (0->125)
t427 = sin(pkin(5));
t432 = sin(qJ(2));
t436 = cos(qJ(2));
t450 = qJD(1) * qJD(2);
t418 = (-qJDD(1) * t436 + t432 * t450) * t427;
t452 = qJD(1) * t427;
t416 = (-pkin(2) * t436 - pkin(8) * t432) * t452;
t429 = cos(pkin(5));
t423 = qJD(1) * t429 + qJD(2);
t421 = t423 ^ 2;
t422 = qJDD(1) * t429 + qJDD(2);
t451 = qJD(1) * t436;
t438 = qJD(1) ^ 2;
t433 = sin(qJ(1));
t437 = cos(qJ(1));
t444 = -g(1) * t437 - g(2) * t433;
t459 = pkin(7) * t427;
t414 = -pkin(1) * t438 + qJDD(1) * t459 + t444;
t447 = t433 * g(1) - g(2) * t437;
t413 = qJDD(1) * pkin(1) + t438 * t459 + t447;
t456 = t413 * t429;
t453 = t436 * t414 + t432 * t456;
t375 = -pkin(2) * t421 + pkin(8) * t422 + (-g(3) * t432 + t416 * t451) * t427 + t453;
t417 = (qJDD(1) * t432 + t436 * t450) * t427;
t458 = g(3) * t429;
t376 = pkin(2) * t418 - pkin(8) * t417 - t458 + (-t413 + (pkin(2) * t432 - pkin(8) * t436) * t423 * qJD(1)) * t427;
t431 = sin(qJ(3));
t435 = cos(qJ(3));
t352 = -t375 * t431 + t435 * t376;
t449 = t432 * t452;
t406 = t423 * t435 - t431 * t449;
t388 = qJD(3) * t406 + t417 * t435 + t422 * t431;
t407 = t423 * t431 + t435 * t449;
t410 = qJDD(3) + t418;
t448 = t427 * t451;
t420 = qJD(3) - t448;
t345 = (t406 * t420 - t388) * qJ(4) + (t406 * t407 + t410) * pkin(3) + t352;
t353 = t435 * t375 + t431 * t376;
t387 = -qJD(3) * t407 - t417 * t431 + t422 * t435;
t397 = pkin(3) * t420 - qJ(4) * t407;
t405 = t406 ^ 2;
t347 = -pkin(3) * t405 + qJ(4) * t387 - t397 * t420 + t353;
t426 = sin(pkin(10));
t428 = cos(pkin(10));
t393 = t406 * t428 - t407 * t426;
t460 = 2 * qJD(4);
t342 = t426 * t345 + t428 * t347 + t393 * t460;
t394 = t406 * t426 + t407 * t428;
t370 = -pkin(4) * t393 - pkin(9) * t394;
t419 = t420 ^ 2;
t340 = -pkin(4) * t419 + pkin(9) * t410 + t370 * t393 + t342;
t454 = t427 * t436;
t389 = -g(3) * t454 - t432 * t414 + t436 * t456;
t374 = -pkin(2) * t422 - pkin(8) * t421 + t416 * t449 - t389;
t348 = -pkin(3) * t387 - qJ(4) * t405 + t407 * t397 + qJDD(4) + t374;
t363 = t387 * t428 - t388 * t426;
t364 = t387 * t426 + t388 * t428;
t343 = (-t393 * t420 - t364) * pkin(9) + (t394 * t420 - t363) * pkin(4) + t348;
t430 = sin(qJ(5));
t434 = cos(qJ(5));
t337 = -t340 * t430 + t343 * t434;
t377 = -t394 * t430 + t420 * t434;
t351 = qJD(5) * t377 + t364 * t434 + t410 * t430;
t378 = t394 * t434 + t420 * t430;
t358 = -mrSges(6,1) * t377 + mrSges(6,2) * t378;
t392 = qJD(5) - t393;
t359 = -mrSges(6,2) * t392 + mrSges(6,3) * t377;
t362 = qJDD(5) - t363;
t334 = m(6) * t337 + mrSges(6,1) * t362 - mrSges(6,3) * t351 - t358 * t378 + t359 * t392;
t338 = t340 * t434 + t343 * t430;
t350 = -qJD(5) * t378 - t364 * t430 + t410 * t434;
t360 = mrSges(6,1) * t392 - mrSges(6,3) * t378;
t335 = m(6) * t338 - mrSges(6,2) * t362 + mrSges(6,3) * t350 + t358 * t377 - t360 * t392;
t326 = -t334 * t430 + t434 * t335;
t369 = -mrSges(5,1) * t393 + mrSges(5,2) * t394;
t380 = mrSges(5,1) * t420 - mrSges(5,3) * t394;
t323 = m(5) * t342 - mrSges(5,2) * t410 + mrSges(5,3) * t363 + t369 * t393 - t380 * t420 + t326;
t443 = -t345 * t428 + t347 * t426;
t339 = -pkin(4) * t410 - pkin(9) * t419 + (t460 + t370) * t394 + t443;
t336 = -m(6) * t339 + t350 * mrSges(6,1) - mrSges(6,2) * t351 + t377 * t359 - t360 * t378;
t341 = -0.2e1 * qJD(4) * t394 - t443;
t379 = -mrSges(5,2) * t420 + mrSges(5,3) * t393;
t330 = m(5) * t341 + mrSges(5,1) * t410 - mrSges(5,3) * t364 - t369 * t394 + t379 * t420 + t336;
t319 = t426 * t323 + t428 * t330;
t354 = Ifges(6,5) * t378 + Ifges(6,6) * t377 + Ifges(6,3) * t392;
t356 = Ifges(6,1) * t378 + Ifges(6,4) * t377 + Ifges(6,5) * t392;
t327 = -mrSges(6,1) * t339 + mrSges(6,3) * t338 + Ifges(6,4) * t351 + Ifges(6,2) * t350 + Ifges(6,6) * t362 - t354 * t378 + t356 * t392;
t355 = Ifges(6,4) * t378 + Ifges(6,2) * t377 + Ifges(6,6) * t392;
t328 = mrSges(6,2) * t339 - mrSges(6,3) * t337 + Ifges(6,1) * t351 + Ifges(6,4) * t350 + Ifges(6,5) * t362 + t354 * t377 - t355 * t392;
t366 = Ifges(5,4) * t394 + Ifges(5,2) * t393 + Ifges(5,6) * t420;
t367 = Ifges(5,1) * t394 + Ifges(5,4) * t393 + Ifges(5,5) * t420;
t382 = Ifges(4,4) * t407 + Ifges(4,2) * t406 + Ifges(4,6) * t420;
t383 = Ifges(4,1) * t407 + Ifges(4,4) * t406 + Ifges(4,5) * t420;
t461 = Ifges(4,5) * t388 + Ifges(4,6) * t387 + t407 * t382 - t406 * t383 + mrSges(4,1) * t352 - mrSges(4,2) * t353 + Ifges(5,5) * t364 + Ifges(5,6) * t363 + t394 * t366 - t393 * t367 + mrSges(5,1) * t341 - mrSges(5,2) * t342 + t430 * t328 + t434 * t327 + pkin(4) * t336 + pkin(9) * t326 + pkin(3) * t319 + (Ifges(4,3) + Ifges(5,3)) * t410;
t455 = t427 * t432;
t395 = -mrSges(4,1) * t406 + mrSges(4,2) * t407;
t396 = -mrSges(4,2) * t420 + mrSges(4,3) * t406;
t317 = m(4) * t352 + mrSges(4,1) * t410 - mrSges(4,3) * t388 - t395 * t407 + t396 * t420 + t319;
t398 = mrSges(4,1) * t420 - mrSges(4,3) * t407;
t445 = t428 * t323 - t330 * t426;
t318 = m(4) * t353 - mrSges(4,2) * t410 + mrSges(4,3) * t387 + t395 * t406 - t398 * t420 + t445;
t311 = t435 * t317 + t431 * t318;
t325 = t434 * t334 + t430 * t335;
t446 = -t317 * t431 + t435 * t318;
t324 = m(5) * t348 - t363 * mrSges(5,1) + mrSges(5,2) * t364 - t393 * t379 + t380 * t394 + t325;
t441 = mrSges(6,1) * t337 - mrSges(6,2) * t338 + Ifges(6,5) * t351 + Ifges(6,6) * t350 + Ifges(6,3) * t362 + t355 * t378 - t356 * t377;
t440 = -m(4) * t374 + t387 * mrSges(4,1) - mrSges(4,2) * t388 + t406 * t396 - t398 * t407 - t324;
t415 = (-mrSges(3,1) * t436 + mrSges(3,2) * t432) * t452;
t412 = -mrSges(3,2) * t423 + mrSges(3,3) * t448;
t411 = mrSges(3,1) * t423 - mrSges(3,3) * t449;
t402 = -t413 * t427 - t458;
t401 = Ifges(3,5) * t423 + (Ifges(3,1) * t432 + Ifges(3,4) * t436) * t452;
t400 = Ifges(3,6) * t423 + (Ifges(3,4) * t432 + Ifges(3,2) * t436) * t452;
t399 = Ifges(3,3) * t423 + (Ifges(3,5) * t432 + Ifges(3,6) * t436) * t452;
t390 = -g(3) * t455 + t453;
t381 = Ifges(4,5) * t407 + Ifges(4,6) * t406 + Ifges(4,3) * t420;
t365 = Ifges(5,5) * t394 + Ifges(5,6) * t393 + Ifges(5,3) * t420;
t320 = m(3) * t389 + mrSges(3,1) * t422 - mrSges(3,3) * t417 + t412 * t423 - t415 * t449 + t440;
t313 = -mrSges(5,1) * t348 + mrSges(5,3) * t342 + Ifges(5,4) * t364 + Ifges(5,2) * t363 + Ifges(5,6) * t410 - pkin(4) * t325 - t365 * t394 + t367 * t420 - t441;
t312 = mrSges(5,2) * t348 - mrSges(5,3) * t341 + Ifges(5,1) * t364 + Ifges(5,4) * t363 + Ifges(5,5) * t410 - pkin(9) * t325 - t327 * t430 + t328 * t434 + t365 * t393 - t366 * t420;
t310 = m(3) * t390 - mrSges(3,2) * t422 - mrSges(3,3) * t418 - t411 * t423 + t415 * t448 + t446;
t309 = mrSges(4,2) * t374 - mrSges(4,3) * t352 + Ifges(4,1) * t388 + Ifges(4,4) * t387 + Ifges(4,5) * t410 - qJ(4) * t319 + t312 * t428 - t313 * t426 + t381 * t406 - t382 * t420;
t308 = -mrSges(4,1) * t374 + mrSges(4,3) * t353 + Ifges(4,4) * t388 + Ifges(4,2) * t387 + Ifges(4,6) * t410 - pkin(3) * t324 + qJ(4) * t445 + t426 * t312 + t428 * t313 - t407 * t381 + t420 * t383;
t307 = Ifges(3,5) * t417 - Ifges(3,6) * t418 + Ifges(3,3) * t422 + mrSges(3,1) * t389 - mrSges(3,2) * t390 + t431 * t309 + t435 * t308 + pkin(2) * t440 + pkin(8) * t446 + (t400 * t432 - t401 * t436) * t452;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t447 - mrSges(2,2) * t444 + (mrSges(3,2) * t402 - mrSges(3,3) * t389 + Ifges(3,1) * t417 - Ifges(3,4) * t418 + Ifges(3,5) * t422 - pkin(8) * t311 - t308 * t431 + t309 * t435 + t399 * t448 - t400 * t423) * t455 + (-mrSges(3,1) * t402 + mrSges(3,3) * t390 + Ifges(3,4) * t417 - Ifges(3,2) * t418 + Ifges(3,6) * t422 - pkin(2) * t311 - t399 * t449 + t423 * t401 - t461) * t454 + t429 * t307 + pkin(1) * ((t310 * t432 + t320 * t436) * t429 + (-m(3) * t402 - t418 * mrSges(3,1) - t417 * mrSges(3,2) + (-t411 * t432 + t412 * t436) * t452 - t311) * t427) + (t310 * t436 - t432 * t320) * t459; t307; t461; t324; t441;];
tauJ = t1;
