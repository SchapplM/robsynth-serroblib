% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:39
% EndTime: 2019-12-31 22:08:46
% DurationCPUTime: 2.85s
% Computational Cost: add. (24947->266), mult. (53413->339), div. (0->0), fcn. (40719->10), ass. (0->117)
t456 = Ifges(5,1) + Ifges(6,1);
t449 = Ifges(5,4) + Ifges(6,4);
t448 = Ifges(5,5) + Ifges(6,5);
t455 = Ifges(5,2) + Ifges(6,2);
t447 = Ifges(5,6) + Ifges(6,6);
t454 = Ifges(5,3) + Ifges(6,3);
t412 = sin(pkin(5));
t416 = sin(qJ(2));
t420 = cos(qJ(2));
t435 = qJD(1) * qJD(2);
t404 = (-qJDD(1) * t420 + t416 * t435) * t412;
t413 = cos(pkin(5));
t409 = t413 * qJD(1) + qJD(2);
t415 = sin(qJ(3));
t419 = cos(qJ(3));
t437 = qJD(1) * t412;
t432 = t416 * t437;
t392 = t419 * t409 - t415 * t432;
t403 = (qJDD(1) * t416 + t420 * t435) * t412;
t408 = t413 * qJDD(1) + qJDD(2);
t375 = t392 * qJD(3) + t419 * t403 + t415 * t408;
t393 = t415 * t409 + t419 * t432;
t436 = qJD(1) * t420;
t431 = t412 * t436;
t406 = qJD(3) - t431;
t414 = sin(qJ(4));
t418 = cos(qJ(4));
t381 = -t414 * t393 + t418 * t406;
t396 = qJDD(3) + t404;
t342 = t381 * qJD(4) + t418 * t375 + t414 * t396;
t382 = t418 * t393 + t414 * t406;
t359 = -t381 * mrSges(6,1) + t382 * mrSges(6,2);
t402 = (-t420 * pkin(2) - t416 * pkin(8)) * t437;
t407 = t409 ^ 2;
t422 = qJD(1) ^ 2;
t417 = sin(qJ(1));
t421 = cos(qJ(1));
t427 = -t421 * g(1) - t417 * g(2);
t452 = t412 * pkin(7);
t400 = -t422 * pkin(1) + qJDD(1) * t452 + t427;
t430 = t417 * g(1) - t421 * g(2);
t399 = qJDD(1) * pkin(1) + t422 * t452 + t430;
t445 = t399 * t413;
t438 = t420 * t400 + t416 * t445;
t356 = -t407 * pkin(2) + t408 * pkin(8) + (-g(3) * t416 + t402 * t436) * t412 + t438;
t451 = t413 * g(3);
t357 = t404 * pkin(2) - t403 * pkin(8) - t451 + (-t399 + (pkin(2) * t416 - pkin(8) * t420) * t409 * qJD(1)) * t412;
t337 = t419 * t356 + t415 * t357;
t379 = -t392 * pkin(3) - t393 * pkin(9);
t405 = t406 ^ 2;
t332 = -t405 * pkin(3) + t396 * pkin(9) + t392 * t379 + t337;
t443 = t412 * t420;
t376 = -g(3) * t443 - t416 * t400 + t420 * t445;
t355 = -t408 * pkin(2) - t407 * pkin(8) + t402 * t432 - t376;
t374 = -t393 * qJD(3) - t415 * t403 + t419 * t408;
t335 = (-t392 * t406 - t375) * pkin(9) + (t393 * t406 - t374) * pkin(3) + t355;
t327 = -t414 * t332 + t418 * t335;
t372 = qJDD(4) - t374;
t391 = qJD(4) - t392;
t324 = -0.2e1 * qJD(5) * t382 + (t381 * t391 - t342) * qJ(5) + (t381 * t382 + t372) * pkin(4) + t327;
t362 = -t391 * mrSges(6,2) + t381 * mrSges(6,3);
t434 = m(6) * t324 + t372 * mrSges(6,1) + t391 * t362;
t321 = -t342 * mrSges(6,3) - t382 * t359 + t434;
t328 = t418 * t332 + t414 * t335;
t341 = -t382 * qJD(4) - t414 * t375 + t418 * t396;
t364 = t391 * pkin(4) - t382 * qJ(5);
t380 = t381 ^ 2;
t326 = -t380 * pkin(4) + t341 * qJ(5) + 0.2e1 * qJD(5) * t381 - t391 * t364 + t328;
t440 = t449 * t381 + t456 * t382 + t448 * t391;
t441 = -t455 * t381 - t449 * t382 - t447 * t391;
t453 = mrSges(5,1) * t327 + mrSges(6,1) * t324 - mrSges(5,2) * t328 - mrSges(6,2) * t326 + pkin(4) * t321 + t447 * t341 + t448 * t342 + t454 * t372 - t440 * t381 - t441 * t382;
t450 = -mrSges(5,2) - mrSges(6,2);
t444 = t412 * t416;
t360 = -t381 * mrSges(5,1) + t382 * mrSges(5,2);
t363 = -t391 * mrSges(5,2) + t381 * mrSges(5,3);
t315 = m(5) * t327 + t372 * mrSges(5,1) + t391 * t363 + (-t359 - t360) * t382 + (-mrSges(5,3) - mrSges(6,3)) * t342 + t434;
t433 = m(6) * t326 + t341 * mrSges(6,3) + t381 * t359;
t365 = t391 * mrSges(6,1) - t382 * mrSges(6,3);
t439 = -t391 * mrSges(5,1) + t382 * mrSges(5,3) - t365;
t317 = m(5) * t328 + t341 * mrSges(5,3) + t381 * t360 + t450 * t372 + t439 * t391 + t433;
t314 = -t414 * t315 + t418 * t317;
t378 = -t392 * mrSges(4,1) + t393 * mrSges(4,2);
t384 = t406 * mrSges(4,1) - t393 * mrSges(4,3);
t311 = m(4) * t337 - t396 * mrSges(4,2) + t374 * mrSges(4,3) + t392 * t378 - t406 * t384 + t314;
t336 = -t415 * t356 + t419 * t357;
t331 = -t396 * pkin(3) - t405 * pkin(9) + t393 * t379 - t336;
t329 = -t341 * pkin(4) - t380 * qJ(5) + t382 * t364 + qJDD(5) + t331;
t428 = -m(6) * t329 + t341 * mrSges(6,1) + t381 * t362;
t320 = -m(5) * t331 + t341 * mrSges(5,1) + t450 * t342 + t381 * t363 + t439 * t382 + t428;
t383 = -t406 * mrSges(4,2) + t392 * mrSges(4,3);
t319 = m(4) * t336 + t396 * mrSges(4,1) - t375 * mrSges(4,3) - t393 * t378 + t406 * t383 + t320;
t306 = t415 * t311 + t419 * t319;
t442 = -t447 * t381 - t448 * t382 - t454 * t391;
t429 = t419 * t311 - t415 * t319;
t313 = t418 * t315 + t414 * t317;
t424 = -m(4) * t355 + t374 * mrSges(4,1) - t375 * mrSges(4,2) + t392 * t383 - t393 * t384 - t313;
t322 = t342 * mrSges(6,2) + t382 * t365 - t428;
t307 = -mrSges(5,1) * t331 + mrSges(5,3) * t328 - mrSges(6,1) * t329 + mrSges(6,3) * t326 - pkin(4) * t322 + qJ(5) * t433 + (-qJ(5) * t365 + t440) * t391 + t442 * t382 + (-qJ(5) * mrSges(6,2) + t447) * t372 + t449 * t342 + t455 * t341;
t312 = mrSges(5,2) * t331 + mrSges(6,2) * t329 - mrSges(5,3) * t327 - mrSges(6,3) * t324 - qJ(5) * t321 + t449 * t341 + t456 * t342 + t448 * t372 - t442 * t381 + t441 * t391;
t369 = Ifges(4,4) * t393 + Ifges(4,2) * t392 + Ifges(4,6) * t406;
t370 = Ifges(4,1) * t393 + Ifges(4,4) * t392 + Ifges(4,5) * t406;
t423 = mrSges(4,1) * t336 - mrSges(4,2) * t337 + Ifges(4,5) * t375 + Ifges(4,6) * t374 + Ifges(4,3) * t396 + pkin(3) * t320 + pkin(9) * t314 + t418 * t307 + t414 * t312 + t393 * t369 - t392 * t370;
t401 = (-t420 * mrSges(3,1) + t416 * mrSges(3,2)) * t437;
t398 = -t409 * mrSges(3,2) + mrSges(3,3) * t431;
t397 = t409 * mrSges(3,1) - mrSges(3,3) * t432;
t388 = -t412 * t399 - t451;
t387 = Ifges(3,5) * t409 + (t416 * Ifges(3,1) + t420 * Ifges(3,4)) * t437;
t386 = Ifges(3,6) * t409 + (t416 * Ifges(3,4) + t420 * Ifges(3,2)) * t437;
t385 = Ifges(3,3) * t409 + (t416 * Ifges(3,5) + t420 * Ifges(3,6)) * t437;
t377 = -g(3) * t444 + t438;
t368 = Ifges(4,5) * t393 + Ifges(4,6) * t392 + Ifges(4,3) * t406;
t308 = m(3) * t376 + t408 * mrSges(3,1) - t403 * mrSges(3,3) + t409 * t398 - t401 * t432 + t424;
t305 = m(3) * t377 - t408 * mrSges(3,2) - t404 * mrSges(3,3) - t409 * t397 + t401 * t431 + t429;
t304 = -mrSges(4,1) * t355 + mrSges(4,3) * t337 + Ifges(4,4) * t375 + Ifges(4,2) * t374 + Ifges(4,6) * t396 - pkin(3) * t313 - t393 * t368 + t406 * t370 - t453;
t303 = mrSges(4,2) * t355 - mrSges(4,3) * t336 + Ifges(4,1) * t375 + Ifges(4,4) * t374 + Ifges(4,5) * t396 - pkin(9) * t313 - t414 * t307 + t418 * t312 + t392 * t368 - t406 * t369;
t302 = Ifges(3,5) * t403 - Ifges(3,6) * t404 + Ifges(3,3) * t408 + mrSges(3,1) * t376 - mrSges(3,2) * t377 + t415 * t303 + t419 * t304 + pkin(2) * t424 + pkin(8) * t429 + (t416 * t386 - t420 * t387) * t437;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t430 - mrSges(2,2) * t427 + (mrSges(3,2) * t388 - mrSges(3,3) * t376 + Ifges(3,1) * t403 - Ifges(3,4) * t404 + Ifges(3,5) * t408 - pkin(8) * t306 + t419 * t303 - t415 * t304 + t385 * t431 - t409 * t386) * t444 + (-mrSges(3,1) * t388 + mrSges(3,3) * t377 + Ifges(3,4) * t403 - Ifges(3,2) * t404 + Ifges(3,6) * t408 - pkin(2) * t306 - t385 * t432 + t409 * t387 - t423) * t443 + t413 * t302 + pkin(1) * ((t416 * t305 + t420 * t308) * t413 + (-m(3) * t388 - t404 * mrSges(3,1) - t403 * mrSges(3,2) + (-t397 * t416 + t398 * t420) * t437 - t306) * t412) + (t420 * t305 - t416 * t308) * t452; t302; t423; t453; t322;];
tauJ = t1;
