% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:41
% EndTime: 2019-12-31 19:07:44
% DurationCPUTime: 2.37s
% Computational Cost: add. (21660->245), mult. (52216->315), div. (0->0), fcn. (39506->10), ass. (0->107)
t415 = qJD(1) ^ 2;
t438 = pkin(2) * t415;
t437 = pkin(6) * qJDD(1);
t410 = sin(qJ(1));
t414 = cos(qJ(1));
t426 = -t414 * g(1) - t410 * g(2);
t393 = -t415 * pkin(1) + qJDD(1) * qJ(2) + t426;
t405 = sin(pkin(9));
t406 = cos(pkin(9));
t432 = qJD(1) * qJD(2);
t430 = -t406 * g(3) - 0.2e1 * t405 * t432;
t370 = (t406 * t438 - t393 - t437) * t405 + t430;
t384 = -t405 * g(3) + (t393 + 0.2e1 * t432) * t406;
t403 = t406 ^ 2;
t371 = -t403 * t438 + t406 * t437 + t384;
t409 = sin(qJ(3));
t413 = cos(qJ(3));
t353 = t413 * t370 - t409 * t371;
t423 = t405 * t413 + t406 * t409;
t422 = -t405 * t409 + t406 * t413;
t391 = t422 * qJD(1);
t433 = t391 * qJD(3);
t382 = t423 * qJDD(1) + t433;
t392 = t423 * qJD(1);
t333 = (-t382 + t433) * pkin(7) + (t391 * t392 + qJDD(3)) * pkin(3) + t353;
t354 = t409 * t370 + t413 * t371;
t381 = -t392 * qJD(3) + t422 * qJDD(1);
t387 = qJD(3) * pkin(3) - t392 * pkin(7);
t390 = t391 ^ 2;
t338 = -t390 * pkin(3) + t381 * pkin(7) - qJD(3) * t387 + t354;
t408 = sin(qJ(4));
t412 = cos(qJ(4));
t331 = t408 * t333 + t412 * t338;
t377 = t408 * t391 + t412 * t392;
t350 = -t377 * qJD(4) + t412 * t381 - t408 * t382;
t376 = t412 * t391 - t408 * t392;
t361 = -t376 * mrSges(5,1) + t377 * mrSges(5,2);
t404 = qJD(3) + qJD(4);
t368 = t404 * mrSges(5,1) - t377 * mrSges(5,3);
t401 = qJDD(3) + qJDD(4);
t362 = -t376 * pkin(4) - t377 * pkin(8);
t400 = t404 ^ 2;
t327 = -t400 * pkin(4) + t401 * pkin(8) + t376 * t362 + t331;
t431 = t410 * g(1) - t414 * g(2);
t425 = qJDD(2) - t431;
t435 = -t405 ^ 2 - t403;
t380 = (-pkin(2) * t406 - pkin(1)) * qJDD(1) + (t435 * pkin(6) - qJ(2)) * t415 + t425;
t345 = -t381 * pkin(3) - t390 * pkin(7) + t392 * t387 + t380;
t351 = t376 * qJD(4) + t408 * t381 + t412 * t382;
t328 = (-t376 * t404 - t351) * pkin(8) + (t377 * t404 - t350) * pkin(4) + t345;
t407 = sin(qJ(5));
t411 = cos(qJ(5));
t324 = -t407 * t327 + t411 * t328;
t363 = -t407 * t377 + t411 * t404;
t336 = t363 * qJD(5) + t411 * t351 + t407 * t401;
t349 = qJDD(5) - t350;
t364 = t411 * t377 + t407 * t404;
t352 = -t363 * mrSges(6,1) + t364 * mrSges(6,2);
t372 = qJD(5) - t376;
t355 = -t372 * mrSges(6,2) + t363 * mrSges(6,3);
t321 = m(6) * t324 + t349 * mrSges(6,1) - t336 * mrSges(6,3) - t364 * t352 + t372 * t355;
t325 = t411 * t327 + t407 * t328;
t335 = -t364 * qJD(5) - t407 * t351 + t411 * t401;
t356 = t372 * mrSges(6,1) - t364 * mrSges(6,3);
t322 = m(6) * t325 - t349 * mrSges(6,2) + t335 * mrSges(6,3) + t363 * t352 - t372 * t356;
t427 = -t407 * t321 + t411 * t322;
t309 = m(5) * t331 - t401 * mrSges(5,2) + t350 * mrSges(5,3) + t376 * t361 - t404 * t368 + t427;
t330 = t412 * t333 - t408 * t338;
t367 = -t404 * mrSges(5,2) + t376 * mrSges(5,3);
t326 = -t401 * pkin(4) - t400 * pkin(8) + t377 * t362 - t330;
t420 = -m(6) * t326 + t335 * mrSges(6,1) - t336 * mrSges(6,2) + t363 * t355 - t364 * t356;
t317 = m(5) * t330 + t401 * mrSges(5,1) - t351 * mrSges(5,3) - t377 * t361 + t404 * t367 + t420;
t305 = t408 * t309 + t412 * t317;
t379 = -t391 * mrSges(4,1) + t392 * mrSges(4,2);
t385 = -qJD(3) * mrSges(4,2) + t391 * mrSges(4,3);
t303 = m(4) * t353 + qJDD(3) * mrSges(4,1) - t382 * mrSges(4,3) + qJD(3) * t385 - t392 * t379 + t305;
t386 = qJD(3) * mrSges(4,1) - t392 * mrSges(4,3);
t428 = t412 * t309 - t408 * t317;
t304 = m(4) * t354 - qJDD(3) * mrSges(4,2) + t381 * mrSges(4,3) - qJD(3) * t386 + t391 * t379 + t428;
t436 = t413 * t303 + t409 * t304;
t311 = t411 * t321 + t407 * t322;
t429 = -t409 * t303 + t413 * t304;
t424 = -t406 * mrSges(3,1) + t405 * mrSges(3,2);
t421 = mrSges(3,3) * qJDD(1) + t415 * t424;
t419 = m(5) * t345 - t350 * mrSges(5,1) + t351 * mrSges(5,2) - t376 * t367 + t377 * t368 + t311;
t339 = Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t372;
t341 = Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t372;
t314 = -mrSges(6,1) * t326 + mrSges(6,3) * t325 + Ifges(6,4) * t336 + Ifges(6,2) * t335 + Ifges(6,6) * t349 - t364 * t339 + t372 * t341;
t340 = Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t372;
t315 = mrSges(6,2) * t326 - mrSges(6,3) * t324 + Ifges(6,1) * t336 + Ifges(6,4) * t335 + Ifges(6,5) * t349 + t363 * t339 - t372 * t340;
t358 = Ifges(5,4) * t377 + Ifges(5,2) * t376 + Ifges(5,6) * t404;
t359 = Ifges(5,1) * t377 + Ifges(5,4) * t376 + Ifges(5,5) * t404;
t418 = mrSges(5,1) * t330 - mrSges(5,2) * t331 + Ifges(5,5) * t351 + Ifges(5,6) * t350 + Ifges(5,3) * t401 + pkin(4) * t420 + pkin(8) * t427 + t411 * t314 + t407 * t315 + t377 * t358 - t376 * t359;
t417 = mrSges(6,1) * t324 - mrSges(6,2) * t325 + Ifges(6,5) * t336 + Ifges(6,6) * t335 + Ifges(6,3) * t349 + t364 * t340 - t363 * t341;
t416 = m(4) * t380 - t381 * mrSges(4,1) + t382 * mrSges(4,2) - t391 * t385 + t392 * t386 + t419;
t389 = -qJDD(1) * pkin(1) - t415 * qJ(2) + t425;
t383 = -t405 * t393 + t430;
t375 = Ifges(4,1) * t392 + Ifges(4,4) * t391 + Ifges(4,5) * qJD(3);
t374 = Ifges(4,4) * t392 + Ifges(4,2) * t391 + Ifges(4,6) * qJD(3);
t373 = Ifges(4,5) * t392 + Ifges(4,6) * t391 + Ifges(4,3) * qJD(3);
t357 = Ifges(5,5) * t377 + Ifges(5,6) * t376 + Ifges(5,3) * t404;
t306 = t435 * t415 * mrSges(3,3) + m(3) * t389 + t424 * qJDD(1) + t416;
t299 = -mrSges(5,1) * t345 + mrSges(5,3) * t331 + Ifges(5,4) * t351 + Ifges(5,2) * t350 + Ifges(5,6) * t401 - pkin(4) * t311 - t377 * t357 + t404 * t359 - t417;
t298 = mrSges(5,2) * t345 - mrSges(5,3) * t330 + Ifges(5,1) * t351 + Ifges(5,4) * t350 + Ifges(5,5) * t401 - pkin(8) * t311 - t407 * t314 + t411 * t315 + t376 * t357 - t404 * t358;
t297 = mrSges(4,2) * t380 - mrSges(4,3) * t353 + Ifges(4,1) * t382 + Ifges(4,4) * t381 + Ifges(4,5) * qJDD(3) - pkin(7) * t305 - qJD(3) * t374 + t412 * t298 - t408 * t299 + t391 * t373;
t296 = -mrSges(4,1) * t380 + mrSges(4,3) * t354 + Ifges(4,4) * t382 + Ifges(4,2) * t381 + Ifges(4,6) * qJDD(3) - pkin(3) * t419 + pkin(7) * t428 + qJD(3) * t375 + t408 * t298 + t412 * t299 - t392 * t373;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t431 - mrSges(2,2) * t426 + t405 * (mrSges(3,2) * t389 - mrSges(3,3) * t383 + t413 * t297 - t409 * t296 - pkin(6) * t436 + (Ifges(3,1) * t405 + Ifges(3,4) * t406) * qJDD(1)) + t406 * (-mrSges(3,1) * t389 + mrSges(3,3) * t384 + t409 * t297 + t413 * t296 - pkin(2) * t416 + pkin(6) * t429 + (Ifges(3,4) * t405 + Ifges(3,2) * t406) * qJDD(1)) - pkin(1) * t306 + qJ(2) * ((m(3) * t384 + t421 * t406 + t429) * t406 + (-m(3) * t383 + t421 * t405 - t436) * t405); t306; mrSges(4,1) * t353 - mrSges(4,2) * t354 + Ifges(4,5) * t382 + Ifges(4,6) * t381 + Ifges(4,3) * qJDD(3) + pkin(3) * t305 + t392 * t374 - t391 * t375 + t418; t418; t417;];
tauJ = t1;
