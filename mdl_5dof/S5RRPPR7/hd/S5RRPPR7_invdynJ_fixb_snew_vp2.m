% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:35:12
% EndTime: 2019-12-31 19:35:14
% DurationCPUTime: 1.54s
% Computational Cost: add. (8005->258), mult. (18445->313), div. (0->0), fcn. (11713->8), ass. (0->107)
t455 = -2 * qJD(3);
t454 = Ifges(4,1) + Ifges(5,2);
t453 = -Ifges(5,1) - Ifges(4,3);
t448 = Ifges(4,4) + Ifges(5,6);
t447 = Ifges(4,5) - Ifges(5,4);
t452 = Ifges(4,2) + Ifges(5,3);
t446 = Ifges(4,6) - Ifges(5,5);
t413 = sin(qJ(2));
t416 = cos(qJ(2));
t434 = qJD(1) * qJD(2);
t402 = t413 * qJDD(1) + t416 * t434;
t419 = qJD(1) ^ 2;
t414 = sin(qJ(1));
t417 = cos(qJ(1));
t430 = -t417 * g(1) - t414 * g(2);
t399 = -t419 * pkin(1) + qJDD(1) * pkin(6) + t430;
t444 = t413 * t399;
t450 = pkin(2) * t419;
t351 = qJDD(2) * pkin(2) - t402 * qJ(3) - t444 + (qJ(3) * t434 + t413 * t450 - g(3)) * t416;
t379 = -t413 * g(3) + t416 * t399;
t403 = t416 * qJDD(1) - t413 * t434;
t438 = qJD(1) * t413;
t404 = qJD(2) * pkin(2) - qJ(3) * t438;
t410 = t416 ^ 2;
t352 = t403 * qJ(3) - qJD(2) * t404 - t410 * t450 + t379;
t411 = sin(pkin(8));
t445 = cos(pkin(8));
t393 = (t411 * t416 + t445 * t413) * qJD(1);
t336 = t445 * t351 - t411 * t352 + t393 * t455;
t433 = t414 * g(1) - t417 * g(2);
t428 = -qJDD(1) * pkin(1) - t433;
t356 = -t403 * pkin(2) + qJDD(3) + t404 * t438 + (-qJ(3) * t410 - pkin(6)) * t419 + t428;
t373 = t411 * t402 - t445 * t403;
t374 = t445 * t402 + t411 * t403;
t381 = qJD(2) * mrSges(4,1) - t393 * mrSges(4,3);
t437 = qJD(1) * t416;
t392 = t411 * t438 - t445 * t437;
t436 = qJD(2) * t392;
t451 = -2 * qJD(4);
t420 = (-t374 + t436) * qJ(4) + t356 + (qJD(2) * pkin(3) + t451) * t393;
t335 = t373 * pkin(3) + t420;
t383 = t393 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t366 = t392 * pkin(3) - t393 * qJ(4);
t418 = qJD(2) ^ 2;
t333 = -qJDD(2) * pkin(3) - t418 * qJ(4) + t393 * t366 + qJDD(4) - t336;
t328 = (t392 * t393 - qJDD(2)) * pkin(7) + (t374 + t436) * pkin(4) + t333;
t384 = t393 * pkin(4) - qJD(2) * pkin(7);
t391 = t392 ^ 2;
t331 = -t391 * pkin(4) - t393 * t384 + (pkin(3) + pkin(7)) * t373 + t420;
t412 = sin(qJ(5));
t415 = cos(qJ(5));
t326 = t415 * t328 - t412 * t331;
t375 = -t412 * qJD(2) + t415 * t392;
t347 = t375 * qJD(5) + t415 * qJDD(2) + t412 * t373;
t376 = t415 * qJD(2) + t412 * t392;
t353 = -t375 * mrSges(6,1) + t376 * mrSges(6,2);
t390 = qJD(5) + t393;
t357 = -t390 * mrSges(6,2) + t375 * mrSges(6,3);
t372 = qJDD(5) + t374;
t323 = m(6) * t326 + t372 * mrSges(6,1) - t347 * mrSges(6,3) - t376 * t353 + t390 * t357;
t327 = t412 * t328 + t415 * t331;
t346 = -t376 * qJD(5) - t412 * qJDD(2) + t415 * t373;
t358 = t390 * mrSges(6,1) - t376 * mrSges(6,3);
t324 = m(6) * t327 - t372 * mrSges(6,2) + t346 * mrSges(6,3) + t375 * t353 - t390 * t358;
t431 = -t412 * t323 + t415 * t324;
t427 = -m(5) * t335 + t374 * mrSges(5,3) + t393 * t383 - t431;
t382 = t392 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t439 = -qJD(2) * mrSges(4,2) - t392 * mrSges(4,3) - t382;
t449 = mrSges(4,1) - mrSges(5,2);
t313 = m(4) * t356 + t374 * mrSges(4,2) + t449 * t373 + t393 * t381 + t439 * t392 - t427;
t367 = t392 * mrSges(4,1) + t393 * mrSges(4,2);
t316 = t415 * t323 + t412 * t324;
t368 = -t392 * mrSges(5,2) - t393 * mrSges(5,3);
t424 = -m(5) * t333 - t374 * mrSges(5,1) - t393 * t368 - t316;
t312 = m(4) * t336 - t374 * mrSges(4,3) + t439 * qJD(2) + t449 * qJDD(2) - t393 * t367 + t424;
t387 = t392 * t455;
t443 = t411 * t351 + t445 * t352;
t337 = t387 + t443;
t426 = t418 * pkin(3) - qJDD(2) * qJ(4) - t443;
t332 = qJD(2) * t451 + ((2 * qJD(3)) + t366) * t392 + t426;
t330 = -t373 * pkin(4) - t391 * pkin(7) - t392 * t366 + t387 + ((2 * qJD(4)) + t384) * qJD(2) - t426;
t425 = -m(6) * t330 + t346 * mrSges(6,1) - t347 * mrSges(6,2) + t375 * t357 - t376 * t358;
t422 = -m(5) * t332 + qJDD(2) * mrSges(5,3) + qJD(2) * t383 - t425;
t321 = m(4) * t337 - qJDD(2) * mrSges(4,2) - qJD(2) * t381 + (-t367 - t368) * t392 + (-mrSges(4,3) - mrSges(5,1)) * t373 + t422;
t310 = t445 * t312 + t411 * t321;
t442 = t453 * qJD(2) + t446 * t392 - t447 * t393;
t441 = -t446 * qJD(2) + t452 * t392 - t448 * t393;
t440 = t447 * qJD(2) - t448 * t392 + t454 * t393;
t432 = -t411 * t312 + t445 * t321;
t340 = Ifges(6,4) * t376 + Ifges(6,2) * t375 + Ifges(6,6) * t390;
t341 = Ifges(6,1) * t376 + Ifges(6,4) * t375 + Ifges(6,5) * t390;
t423 = mrSges(6,1) * t326 - mrSges(6,2) * t327 + Ifges(6,5) * t347 + Ifges(6,6) * t346 + Ifges(6,3) * t372 + t376 * t340 - t375 * t341;
t406 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t437;
t405 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t438;
t401 = (-t416 * mrSges(3,1) + t413 * mrSges(3,2)) * qJD(1);
t398 = -t419 * pkin(6) + t428;
t396 = Ifges(3,5) * qJD(2) + (t413 * Ifges(3,1) + t416 * Ifges(3,4)) * qJD(1);
t395 = Ifges(3,6) * qJD(2) + (t413 * Ifges(3,4) + t416 * Ifges(3,2)) * qJD(1);
t378 = -t416 * g(3) - t444;
t339 = Ifges(6,5) * t376 + Ifges(6,6) * t375 + Ifges(6,3) * t390;
t318 = mrSges(6,2) * t330 - mrSges(6,3) * t326 + Ifges(6,1) * t347 + Ifges(6,4) * t346 + Ifges(6,5) * t372 + t375 * t339 - t390 * t340;
t317 = -mrSges(6,1) * t330 + mrSges(6,3) * t327 + Ifges(6,4) * t347 + Ifges(6,2) * t346 + Ifges(6,6) * t372 - t376 * t339 + t390 * t341;
t315 = qJDD(2) * mrSges(5,2) + qJD(2) * t382 - t424;
t314 = -t373 * mrSges(5,2) - t392 * t382 - t427;
t309 = mrSges(5,1) * t333 + mrSges(4,2) * t356 - mrSges(4,3) * t336 - mrSges(5,3) * t335 + pkin(4) * t316 - qJ(4) * t314 + t441 * qJD(2) + t447 * qJDD(2) - t448 * t373 + t454 * t374 + t442 * t392 + t423;
t308 = -mrSges(4,1) * t356 - mrSges(5,1) * t332 + mrSges(5,2) * t335 + mrSges(4,3) * t337 - pkin(3) * t314 - pkin(4) * t425 - pkin(7) * t431 + t440 * qJD(2) + t446 * qJDD(2) - t415 * t317 - t412 * t318 - t452 * t373 + t448 * t374 + t442 * t393;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t433 - mrSges(2,2) * t430 + t413 * (mrSges(3,2) * t398 - mrSges(3,3) * t378 + Ifges(3,1) * t402 + Ifges(3,4) * t403 + Ifges(3,5) * qJDD(2) - qJ(3) * t310 - qJD(2) * t395 - t411 * t308 + t445 * t309) + t416 * (-mrSges(3,1) * t398 + mrSges(3,3) * t379 + Ifges(3,4) * t402 + Ifges(3,2) * t403 + Ifges(3,6) * qJDD(2) - pkin(2) * t313 + qJ(3) * t432 + qJD(2) * t396 + t445 * t308 + t411 * t309) + pkin(1) * (-m(3) * t398 + t403 * mrSges(3,1) - t402 * mrSges(3,2) + (-t405 * t413 + t406 * t416) * qJD(1) - t313) + pkin(6) * (t416 * (m(3) * t379 - qJDD(2) * mrSges(3,2) + t403 * mrSges(3,3) - qJD(2) * t405 + t401 * t437 + t432) - t413 * (m(3) * t378 + qJDD(2) * mrSges(3,1) - t402 * mrSges(3,3) + qJD(2) * t406 - t401 * t438 + t310)); -pkin(3) * t315 - pkin(7) * t316 + pkin(2) * t310 + qJ(4) * t422 + Ifges(3,5) * t402 + Ifges(3,6) * t403 + mrSges(3,1) * t378 - mrSges(3,2) * t379 - t412 * t317 + t415 * t318 + mrSges(4,1) * t336 - mrSges(4,2) * t337 - mrSges(5,3) * t332 + mrSges(5,2) * t333 - t441 * t393 + (-qJ(4) * t368 + t440) * t392 + t447 * t374 + (-qJ(4) * mrSges(5,1) - t446) * t373 + (t413 * t395 - t416 * t396) * qJD(1) + (Ifges(3,3) - t453) * qJDD(2); t313; t315; t423;];
tauJ = t1;
