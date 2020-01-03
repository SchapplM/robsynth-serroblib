% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:56
% EndTime: 2019-12-31 22:33:03
% DurationCPUTime: 5.00s
% Computational Cost: add. (56931->289), mult. (122181->380), div. (0->0), fcn. (96533->12), ass. (0->126)
t434 = sin(pkin(5));
t468 = t434 * pkin(7);
t435 = cos(pkin(5));
t467 = t435 * g(3);
t446 = qJD(1) ^ 2;
t440 = sin(qJ(1));
t445 = cos(qJ(1));
t457 = t440 * g(1) - t445 * g(2);
t419 = qJDD(1) * pkin(1) + t446 * t468 + t457;
t466 = t419 * t435;
t439 = sin(qJ(2));
t465 = t434 * t439;
t444 = cos(qJ(2));
t464 = t434 * t444;
t462 = qJD(1) * t434;
t422 = (-t444 * pkin(2) - t439 * pkin(8)) * t462;
t431 = t435 * qJD(1) + qJD(2);
t429 = t431 ^ 2;
t430 = t435 * qJDD(1) + qJDD(2);
t461 = qJD(1) * t444;
t453 = -t445 * g(1) - t440 * g(2);
t460 = qJDD(1) * t434;
t420 = -t446 * pkin(1) + pkin(7) * t460 + t453;
t463 = t444 * t420 + t439 * t466;
t380 = -t429 * pkin(2) + t430 * pkin(8) + (-g(3) * t439 + t422 * t461) * t434 + t463;
t423 = (qJD(2) * t461 + qJDD(1) * t439) * t434;
t459 = t439 * t462;
t424 = -qJD(2) * t459 + t444 * t460;
t381 = -t424 * pkin(2) - t423 * pkin(8) - t467 + (-t419 + (pkin(2) * t439 - pkin(8) * t444) * t431 * qJD(1)) * t434;
t438 = sin(qJ(3));
t443 = cos(qJ(3));
t354 = -t438 * t380 + t443 * t381;
t411 = t443 * t431 - t438 * t459;
t393 = t411 * qJD(3) + t443 * t423 + t438 * t430;
t412 = t438 * t431 + t443 * t459;
t416 = qJDD(3) - t424;
t458 = t434 * t461;
t427 = qJD(3) - t458;
t347 = (t411 * t427 - t393) * pkin(9) + (t411 * t412 + t416) * pkin(3) + t354;
t355 = t443 * t380 + t438 * t381;
t392 = -t412 * qJD(3) - t438 * t423 + t443 * t430;
t402 = t427 * pkin(3) - t412 * pkin(9);
t410 = t411 ^ 2;
t349 = -t410 * pkin(3) + t392 * pkin(9) - t427 * t402 + t355;
t437 = sin(qJ(4));
t442 = cos(qJ(4));
t345 = t437 * t347 + t442 * t349;
t398 = t437 * t411 + t442 * t412;
t364 = -t398 * qJD(4) + t442 * t392 - t437 * t393;
t397 = t442 * t411 - t437 * t412;
t374 = -t397 * mrSges(5,1) + t398 * mrSges(5,2);
t426 = qJD(4) + t427;
t385 = t426 * mrSges(5,1) - t398 * mrSges(5,3);
t415 = qJDD(4) + t416;
t375 = -t397 * pkin(4) - t398 * pkin(10);
t425 = t426 ^ 2;
t341 = -t425 * pkin(4) + t415 * pkin(10) + t397 * t375 + t345;
t394 = -g(3) * t464 - t439 * t420 + t444 * t466;
t379 = -t430 * pkin(2) - t429 * pkin(8) + t422 * t459 - t394;
t353 = -t392 * pkin(3) - t410 * pkin(9) + t412 * t402 + t379;
t365 = t397 * qJD(4) + t437 * t392 + t442 * t393;
t342 = (-t397 * t426 - t365) * pkin(10) + (t398 * t426 - t364) * pkin(4) + t353;
t436 = sin(qJ(5));
t441 = cos(qJ(5));
t338 = -t436 * t341 + t441 * t342;
t382 = -t436 * t398 + t441 * t426;
t352 = t382 * qJD(5) + t441 * t365 + t436 * t415;
t363 = qJDD(5) - t364;
t383 = t441 * t398 + t436 * t426;
t367 = -t382 * mrSges(6,1) + t383 * mrSges(6,2);
t396 = qJD(5) - t397;
t368 = -t396 * mrSges(6,2) + t382 * mrSges(6,3);
t334 = m(6) * t338 + t363 * mrSges(6,1) - t352 * mrSges(6,3) - t383 * t367 + t396 * t368;
t339 = t441 * t341 + t436 * t342;
t351 = -t383 * qJD(5) - t436 * t365 + t441 * t415;
t369 = t396 * mrSges(6,1) - t383 * mrSges(6,3);
t335 = m(6) * t339 - t363 * mrSges(6,2) + t351 * mrSges(6,3) + t382 * t367 - t396 * t369;
t454 = -t436 * t334 + t441 * t335;
t321 = m(5) * t345 - t415 * mrSges(5,2) + t364 * mrSges(5,3) + t397 * t374 - t426 * t385 + t454;
t344 = t442 * t347 - t437 * t349;
t384 = -t426 * mrSges(5,2) + t397 * mrSges(5,3);
t340 = -t415 * pkin(4) - t425 * pkin(10) + t398 * t375 - t344;
t452 = -m(6) * t340 + t351 * mrSges(6,1) - t352 * mrSges(6,2) + t382 * t368 - t383 * t369;
t330 = m(5) * t344 + t415 * mrSges(5,1) - t365 * mrSges(5,3) - t398 * t374 + t426 * t384 + t452;
t317 = t437 * t321 + t442 * t330;
t399 = -t411 * mrSges(4,1) + t412 * mrSges(4,2);
t400 = -t427 * mrSges(4,2) + t411 * mrSges(4,3);
t315 = m(4) * t354 + t416 * mrSges(4,1) - t393 * mrSges(4,3) - t412 * t399 + t427 * t400 + t317;
t401 = t427 * mrSges(4,1) - t412 * mrSges(4,3);
t455 = t442 * t321 - t437 * t330;
t316 = m(4) * t355 - t416 * mrSges(4,2) + t392 * mrSges(4,3) + t411 * t399 - t427 * t401 + t455;
t309 = t443 * t315 + t438 * t316;
t323 = t441 * t334 + t436 * t335;
t456 = -t438 * t315 + t443 * t316;
t451 = m(5) * t353 - t364 * mrSges(5,1) + t365 * mrSges(5,2) - t397 * t384 + t398 * t385 + t323;
t356 = Ifges(6,5) * t383 + Ifges(6,6) * t382 + Ifges(6,3) * t396;
t358 = Ifges(6,1) * t383 + Ifges(6,4) * t382 + Ifges(6,5) * t396;
t327 = -mrSges(6,1) * t340 + mrSges(6,3) * t339 + Ifges(6,4) * t352 + Ifges(6,2) * t351 + Ifges(6,6) * t363 - t383 * t356 + t396 * t358;
t357 = Ifges(6,4) * t383 + Ifges(6,2) * t382 + Ifges(6,6) * t396;
t328 = mrSges(6,2) * t340 - mrSges(6,3) * t338 + Ifges(6,1) * t352 + Ifges(6,4) * t351 + Ifges(6,5) * t363 + t382 * t356 - t396 * t357;
t371 = Ifges(5,4) * t398 + Ifges(5,2) * t397 + Ifges(5,6) * t426;
t372 = Ifges(5,1) * t398 + Ifges(5,4) * t397 + Ifges(5,5) * t426;
t450 = -mrSges(5,1) * t344 + mrSges(5,2) * t345 - Ifges(5,5) * t365 - Ifges(5,6) * t364 - Ifges(5,3) * t415 - pkin(4) * t452 - pkin(10) * t454 - t441 * t327 - t436 * t328 - t398 * t371 + t397 * t372;
t449 = mrSges(6,1) * t338 - mrSges(6,2) * t339 + Ifges(6,5) * t352 + Ifges(6,6) * t351 + Ifges(6,3) * t363 + t383 * t357 - t382 * t358;
t448 = -m(4) * t379 + t392 * mrSges(4,1) - t393 * mrSges(4,2) + t411 * t400 - t412 * t401 - t451;
t387 = Ifges(4,4) * t412 + Ifges(4,2) * t411 + Ifges(4,6) * t427;
t388 = Ifges(4,1) * t412 + Ifges(4,4) * t411 + Ifges(4,5) * t427;
t447 = mrSges(4,1) * t354 - mrSges(4,2) * t355 + Ifges(4,5) * t393 + Ifges(4,6) * t392 + Ifges(4,3) * t416 + pkin(3) * t317 + t412 * t387 - t411 * t388 - t450;
t421 = (-t444 * mrSges(3,1) + t439 * mrSges(3,2)) * t462;
t418 = -t431 * mrSges(3,2) + mrSges(3,3) * t458;
t417 = t431 * mrSges(3,1) - mrSges(3,3) * t459;
t406 = -t434 * t419 - t467;
t405 = Ifges(3,5) * t431 + (t439 * Ifges(3,1) + t444 * Ifges(3,4)) * t462;
t404 = Ifges(3,6) * t431 + (t439 * Ifges(3,4) + t444 * Ifges(3,2)) * t462;
t403 = Ifges(3,3) * t431 + (t439 * Ifges(3,5) + t444 * Ifges(3,6)) * t462;
t395 = -g(3) * t465 + t463;
t386 = Ifges(4,5) * t412 + Ifges(4,6) * t411 + Ifges(4,3) * t427;
t370 = Ifges(5,5) * t398 + Ifges(5,6) * t397 + Ifges(5,3) * t426;
t318 = m(3) * t394 + t430 * mrSges(3,1) - t423 * mrSges(3,3) + t431 * t418 - t421 * t459 + t448;
t311 = -mrSges(5,1) * t353 + mrSges(5,3) * t345 + Ifges(5,4) * t365 + Ifges(5,2) * t364 + Ifges(5,6) * t415 - pkin(4) * t323 - t398 * t370 + t426 * t372 - t449;
t310 = mrSges(5,2) * t353 - mrSges(5,3) * t344 + Ifges(5,1) * t365 + Ifges(5,4) * t364 + Ifges(5,5) * t415 - pkin(10) * t323 - t436 * t327 + t441 * t328 + t397 * t370 - t426 * t371;
t308 = m(3) * t395 - t430 * mrSges(3,2) + t424 * mrSges(3,3) - t431 * t417 + t421 * t458 + t456;
t307 = mrSges(4,2) * t379 - mrSges(4,3) * t354 + Ifges(4,1) * t393 + Ifges(4,4) * t392 + Ifges(4,5) * t416 - pkin(9) * t317 + t442 * t310 - t437 * t311 + t411 * t386 - t427 * t387;
t306 = -mrSges(4,1) * t379 + mrSges(4,3) * t355 + Ifges(4,4) * t393 + Ifges(4,2) * t392 + Ifges(4,6) * t416 - pkin(3) * t451 + pkin(9) * t455 + t437 * t310 + t442 * t311 - t412 * t386 + t427 * t388;
t305 = Ifges(3,5) * t423 + Ifges(3,6) * t424 + Ifges(3,3) * t430 + mrSges(3,1) * t394 - mrSges(3,2) * t395 + t438 * t307 + t443 * t306 + pkin(2) * t448 + pkin(8) * t456 + (t439 * t404 - t444 * t405) * t462;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t457 - mrSges(2,2) * t453 + (mrSges(3,2) * t406 - mrSges(3,3) * t394 + Ifges(3,1) * t423 + Ifges(3,4) * t424 + Ifges(3,5) * t430 - pkin(8) * t309 - t438 * t306 + t443 * t307 + t403 * t458 - t431 * t404) * t465 + (-mrSges(3,1) * t406 + mrSges(3,3) * t395 + Ifges(3,4) * t423 + Ifges(3,2) * t424 + Ifges(3,6) * t430 - pkin(2) * t309 - t403 * t459 + t431 * t405 - t447) * t464 + t435 * t305 + pkin(1) * ((t439 * t308 + t444 * t318) * t435 + (-m(3) * t406 + t424 * mrSges(3,1) - t423 * mrSges(3,2) + (-t417 * t439 + t418 * t444) * t462 - t309) * t434) + (t444 * t308 - t439 * t318) * t468; t305; t447; -t450; t449;];
tauJ = t1;
