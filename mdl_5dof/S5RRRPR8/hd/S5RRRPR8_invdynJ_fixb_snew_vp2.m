% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:23
% EndTime: 2019-12-31 21:19:26
% DurationCPUTime: 1.94s
% Computational Cost: add. (10453->256), mult. (21748->311), div. (0->0), fcn. (14152->8), ass. (0->107)
t447 = Ifges(4,4) + Ifges(5,6);
t458 = -Ifges(4,2) - Ifges(5,3);
t455 = Ifges(4,6) - Ifges(5,5);
t457 = Ifges(4,1) + Ifges(5,2);
t456 = Ifges(4,5) - Ifges(5,4);
t454 = Ifges(4,3) + Ifges(5,1);
t418 = sin(qJ(3));
t422 = cos(qJ(2));
t440 = qJD(1) * t422;
t419 = sin(qJ(2));
t441 = qJD(1) * t419;
t450 = cos(qJ(3));
t396 = t418 * t441 - t450 * t440;
t397 = (t418 * t422 + t450 * t419) * qJD(1);
t415 = qJD(2) + qJD(3);
t453 = t458 * t396 + t447 * t397 + t455 * t415;
t439 = qJD(1) * qJD(2);
t403 = t419 * qJDD(1) + t422 * t439;
t404 = t422 * qJDD(1) - t419 * t439;
t363 = t397 * qJD(3) + t418 * t403 - t450 * t404;
t364 = -t396 * qJD(3) + t450 * t403 + t418 * t404;
t407 = qJD(2) * pkin(2) - pkin(7) * t441;
t416 = t422 ^ 2;
t424 = qJD(1) ^ 2;
t420 = sin(qJ(1));
t423 = cos(qJ(1));
t438 = t420 * g(1) - t423 * g(2);
t434 = -qJDD(1) * pkin(1) - t438;
t365 = -t404 * pkin(2) + t407 * t441 + (-pkin(7) * t416 - pkin(6)) * t424 + t434;
t384 = t415 * mrSges(4,1) - t397 * mrSges(4,3);
t446 = t396 * t415;
t451 = -2 * qJD(4);
t426 = (-t364 + t446) * qJ(4) + t365 + (t415 * pkin(3) + t451) * t397;
t328 = t363 * pkin(3) + t426;
t386 = t397 * mrSges(5,1) + t415 * mrSges(5,2);
t387 = t397 * pkin(4) - t415 * pkin(8);
t392 = t396 ^ 2;
t323 = -t392 * pkin(4) - t397 * t387 + (pkin(3) + pkin(8)) * t363 + t426;
t435 = -t423 * g(1) - t420 * g(2);
t399 = -t424 * pkin(1) + qJDD(1) * pkin(6) + t435;
t445 = t419 * t399;
t449 = pkin(2) * t424;
t351 = qJDD(2) * pkin(2) - t403 * pkin(7) - t445 + (pkin(7) * t439 + t419 * t449 - g(3)) * t422;
t382 = -t419 * g(3) + t422 * t399;
t352 = t404 * pkin(7) - qJD(2) * t407 - t416 * t449 + t382;
t333 = t450 * t351 - t418 * t352;
t375 = t396 * pkin(3) - t397 * qJ(4);
t413 = t415 ^ 2;
t414 = qJDD(2) + qJDD(3);
t331 = -t414 * pkin(3) - t413 * qJ(4) + t397 * t375 + qJDD(4) - t333;
t324 = (t396 * t397 - t414) * pkin(8) + (t364 + t446) * pkin(4) + t331;
t417 = sin(qJ(5));
t421 = cos(qJ(5));
t321 = -t417 * t323 + t421 * t324;
t379 = t421 * t396 - t417 * t415;
t340 = t379 * qJD(5) + t417 * t363 + t421 * t414;
t380 = t417 * t396 + t421 * t415;
t348 = -t379 * mrSges(6,1) + t380 * mrSges(6,2);
t362 = qJDD(5) + t364;
t391 = qJD(5) + t397;
t366 = -t391 * mrSges(6,2) + t379 * mrSges(6,3);
t318 = m(6) * t321 + t362 * mrSges(6,1) - t340 * mrSges(6,3) - t380 * t348 + t391 * t366;
t322 = t421 * t323 + t417 * t324;
t339 = -t380 * qJD(5) + t421 * t363 - t417 * t414;
t367 = t391 * mrSges(6,1) - t380 * mrSges(6,3);
t319 = m(6) * t322 - t362 * mrSges(6,2) + t339 * mrSges(6,3) + t379 * t348 - t391 * t367;
t436 = -t417 * t318 + t421 * t319;
t433 = -m(5) * t328 + t364 * mrSges(5,3) + t397 * t386 - t436;
t385 = t396 * mrSges(5,1) - t415 * mrSges(5,3);
t442 = -t415 * mrSges(4,2) - t396 * mrSges(4,3) - t385;
t448 = mrSges(4,1) - mrSges(5,2);
t452 = m(4) * t365 + t364 * mrSges(4,2) + t448 * t363 + t397 * t384 + t442 * t396 - t433;
t376 = t396 * mrSges(4,1) + t397 * mrSges(4,2);
t309 = t421 * t318 + t417 * t319;
t377 = -t396 * mrSges(5,2) - t397 * mrSges(5,3);
t431 = -m(5) * t331 - t364 * mrSges(5,1) - t397 * t377 - t309;
t305 = m(4) * t333 - t364 * mrSges(4,3) - t397 * t376 + t448 * t414 + t442 * t415 + t431;
t334 = t418 * t351 + t450 * t352;
t430 = -t413 * pkin(3) + t414 * qJ(4) - t396 * t375 + t334;
t329 = t415 * t451 - t430;
t326 = -t363 * pkin(4) - t392 * pkin(8) + ((2 * qJD(4)) + t387) * t415 + t430;
t432 = -m(6) * t326 + t339 * mrSges(6,1) - t340 * mrSges(6,2) + t379 * t366 - t380 * t367;
t428 = -m(5) * t329 + t414 * mrSges(5,3) + t415 * t386 - t432;
t315 = m(4) * t334 - t414 * mrSges(4,2) - t415 * t384 + (-t376 - t377) * t396 + (-mrSges(4,3) - mrSges(5,1)) * t363 + t428;
t303 = t450 * t305 + t418 * t315;
t444 = t455 * t396 - t456 * t397 - t454 * t415;
t443 = -t447 * t396 + t457 * t397 + t456 * t415;
t437 = -t418 * t305 + t450 * t315;
t342 = Ifges(6,4) * t380 + Ifges(6,2) * t379 + Ifges(6,6) * t391;
t343 = Ifges(6,1) * t380 + Ifges(6,4) * t379 + Ifges(6,5) * t391;
t429 = mrSges(6,1) * t321 - mrSges(6,2) * t322 + Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t362 + t380 * t342 - t379 * t343;
t308 = t414 * mrSges(5,2) + t415 * t385 - t431;
t341 = Ifges(6,5) * t380 + Ifges(6,6) * t379 + Ifges(6,3) * t391;
t311 = -mrSges(6,1) * t326 + mrSges(6,3) * t322 + Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t362 - t380 * t341 + t391 * t343;
t312 = mrSges(6,2) * t326 - mrSges(6,3) * t321 + Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t362 + t379 * t341 - t391 * t342;
t425 = -mrSges(4,2) * t334 - mrSges(5,3) * t329 - pkin(3) * t308 - pkin(8) * t309 - t417 * t311 + t421 * t312 + t443 * t396 + qJ(4) * (-t396 * t377 + t428) + mrSges(5,2) * t331 + mrSges(4,1) * t333 + t454 * t414 + t453 * t397 + t456 * t364 + (-qJ(4) * mrSges(5,1) - t455) * t363;
t406 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t440;
t405 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t441;
t402 = (-t422 * mrSges(3,1) + t419 * mrSges(3,2)) * qJD(1);
t398 = -t424 * pkin(6) + t434;
t395 = Ifges(3,5) * qJD(2) + (t419 * Ifges(3,1) + t422 * Ifges(3,4)) * qJD(1);
t394 = Ifges(3,6) * qJD(2) + (t419 * Ifges(3,4) + t422 * Ifges(3,2)) * qJD(1);
t381 = -t422 * g(3) - t445;
t306 = -t363 * mrSges(5,2) - t396 * t385 - t433;
t302 = mrSges(5,1) * t331 + mrSges(4,2) * t365 - mrSges(4,3) * t333 - mrSges(5,3) * t328 + pkin(4) * t309 - qJ(4) * t306 - t447 * t363 + t457 * t364 + t444 * t396 + t456 * t414 - t453 * t415 + t429;
t301 = -mrSges(4,1) * t365 - mrSges(5,1) * t329 + mrSges(5,2) * t328 + mrSges(4,3) * t334 - pkin(3) * t306 - pkin(4) * t432 - pkin(8) * t436 - t421 * t311 - t417 * t312 + t458 * t363 + t447 * t364 + t444 * t397 + t455 * t414 + t443 * t415;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t438 - mrSges(2,2) * t435 + t419 * (mrSges(3,2) * t398 - mrSges(3,3) * t381 + Ifges(3,1) * t403 + Ifges(3,4) * t404 + Ifges(3,5) * qJDD(2) - pkin(7) * t303 - qJD(2) * t394 - t418 * t301 + t450 * t302) + t422 * (-mrSges(3,1) * t398 + mrSges(3,3) * t382 + Ifges(3,4) * t403 + Ifges(3,2) * t404 + Ifges(3,6) * qJDD(2) - pkin(2) * t452 + pkin(7) * t437 + qJD(2) * t395 + t450 * t301 + t418 * t302) + pkin(1) * (-m(3) * t398 + t404 * mrSges(3,1) - t403 * mrSges(3,2) + (-t405 * t419 + t406 * t422) * qJD(1) - t452) + pkin(6) * (t422 * (m(3) * t382 - qJDD(2) * mrSges(3,2) + t404 * mrSges(3,3) - qJD(2) * t405 + t402 * t440 + t437) - t419 * (m(3) * t381 + qJDD(2) * mrSges(3,1) - t403 * mrSges(3,3) + qJD(2) * t406 - t402 * t441 + t303)); Ifges(3,3) * qJDD(2) + (t419 * t394 - t422 * t395) * qJD(1) + Ifges(3,5) * t403 + Ifges(3,6) * t404 + mrSges(3,1) * t381 - mrSges(3,2) * t382 + pkin(2) * t303 + t425; t425; t308; t429;];
tauJ = t1;
