% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR10
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:06
% EndTime: 2019-12-31 19:10:08
% DurationCPUTime: 2.36s
% Computational Cost: add. (21359->250), mult. (50414->323), div. (0->0), fcn. (36914->10), ass. (0->110)
t418 = qJD(1) ^ 2;
t441 = pkin(2) * t418;
t412 = sin(qJ(1));
t416 = cos(qJ(1));
t428 = -t416 * g(1) - t412 * g(2);
t397 = -t418 * pkin(1) + qJDD(1) * qJ(2) + t428;
t407 = sin(pkin(9));
t408 = cos(pkin(9));
t435 = qJD(1) * qJD(2);
t432 = -t408 * g(3) - 0.2e1 * t407 * t435;
t370 = (-pkin(6) * qJDD(1) + t408 * t441 - t397) * t407 + t432;
t384 = -t407 * g(3) + (t397 + 0.2e1 * t435) * t408;
t406 = t408 ^ 2;
t434 = qJDD(1) * t408;
t371 = pkin(6) * t434 - t406 * t441 + t384;
t411 = sin(qJ(3));
t415 = cos(qJ(3));
t350 = t411 * t370 + t415 * t371;
t437 = qJD(1) * t408;
t438 = qJD(1) * t407;
t395 = -t411 * t438 + t415 * t437;
t425 = t407 * t415 + t408 * t411;
t396 = t425 * qJD(1);
t377 = -t395 * mrSges(4,1) + t396 * mrSges(4,2);
t392 = t396 * qJD(3);
t381 = -t411 * t407 * qJDD(1) + t415 * t434 - t392;
t389 = qJD(3) * mrSges(4,1) - t396 * mrSges(4,3);
t379 = -t395 * pkin(3) - t396 * pkin(7);
t417 = qJD(3) ^ 2;
t340 = -t417 * pkin(3) + qJDD(3) * pkin(7) + t395 * t379 + t350;
t433 = t412 * g(1) - t416 * g(2);
t427 = qJDD(2) - t433;
t439 = -t407 ^ 2 - t406;
t380 = (-pkin(2) * t408 - pkin(1)) * qJDD(1) + (t439 * pkin(6) - qJ(2)) * t418 + t427;
t436 = t395 * qJD(3);
t382 = t425 * qJDD(1) + t436;
t343 = (-t382 - t436) * pkin(7) + (-t381 + t392) * pkin(3) + t380;
t410 = sin(qJ(4));
t414 = cos(qJ(4));
t330 = -t410 * t340 + t414 * t343;
t386 = t414 * qJD(3) - t410 * t396;
t359 = t386 * qJD(4) + t410 * qJDD(3) + t414 * t382;
t378 = qJDD(4) - t381;
t387 = t410 * qJD(3) + t414 * t396;
t393 = qJD(4) - t395;
t327 = (t386 * t393 - t359) * pkin(8) + (t386 * t387 + t378) * pkin(4) + t330;
t331 = t414 * t340 + t410 * t343;
t358 = -t387 * qJD(4) + t414 * qJDD(3) - t410 * t382;
t369 = t393 * pkin(4) - t387 * pkin(8);
t385 = t386 ^ 2;
t328 = -t385 * pkin(4) + t358 * pkin(8) - t393 * t369 + t331;
t409 = sin(qJ(5));
t413 = cos(qJ(5));
t325 = t413 * t327 - t409 * t328;
t360 = t413 * t386 - t409 * t387;
t336 = t360 * qJD(5) + t409 * t358 + t413 * t359;
t361 = t409 * t386 + t413 * t387;
t348 = -t360 * mrSges(6,1) + t361 * mrSges(6,2);
t391 = qJD(5) + t393;
t351 = -t391 * mrSges(6,2) + t360 * mrSges(6,3);
t376 = qJDD(5) + t378;
t322 = m(6) * t325 + t376 * mrSges(6,1) - t336 * mrSges(6,3) - t361 * t348 + t391 * t351;
t326 = t409 * t327 + t413 * t328;
t335 = -t361 * qJD(5) + t413 * t358 - t409 * t359;
t352 = t391 * mrSges(6,1) - t361 * mrSges(6,3);
t323 = m(6) * t326 - t376 * mrSges(6,2) + t335 * mrSges(6,3) + t360 * t348 - t391 * t352;
t314 = t413 * t322 + t409 * t323;
t362 = -t386 * mrSges(5,1) + t387 * mrSges(5,2);
t365 = -t393 * mrSges(5,2) + t386 * mrSges(5,3);
t312 = m(5) * t330 + t378 * mrSges(5,1) - t359 * mrSges(5,3) - t387 * t362 + t393 * t365 + t314;
t366 = t393 * mrSges(5,1) - t387 * mrSges(5,3);
t429 = -t409 * t322 + t413 * t323;
t313 = m(5) * t331 - t378 * mrSges(5,2) + t358 * mrSges(5,3) + t386 * t362 - t393 * t366 + t429;
t430 = -t410 * t312 + t414 * t313;
t307 = m(4) * t350 - qJDD(3) * mrSges(4,2) + t381 * mrSges(4,3) - qJD(3) * t389 + t395 * t377 + t430;
t349 = t415 * t370 - t411 * t371;
t388 = -qJD(3) * mrSges(4,2) + t395 * mrSges(4,3);
t339 = -qJDD(3) * pkin(3) - t417 * pkin(7) + t396 * t379 - t349;
t329 = -t358 * pkin(4) - t385 * pkin(8) + t387 * t369 + t339;
t423 = m(6) * t329 - t335 * mrSges(6,1) + t336 * mrSges(6,2) - t360 * t351 + t361 * t352;
t420 = -m(5) * t339 + t358 * mrSges(5,1) - t359 * mrSges(5,2) + t386 * t365 - t387 * t366 - t423;
t318 = m(4) * t349 + qJDD(3) * mrSges(4,1) - t382 * mrSges(4,3) + qJD(3) * t388 - t396 * t377 + t420;
t440 = t411 * t307 + t415 * t318;
t308 = t414 * t312 + t410 * t313;
t431 = t415 * t307 - t411 * t318;
t426 = -t408 * mrSges(3,1) + t407 * mrSges(3,2);
t424 = mrSges(3,3) * qJDD(1) + t418 * t426;
t345 = Ifges(6,4) * t361 + Ifges(6,2) * t360 + Ifges(6,6) * t391;
t346 = Ifges(6,1) * t361 + Ifges(6,4) * t360 + Ifges(6,5) * t391;
t422 = -mrSges(6,1) * t325 + mrSges(6,2) * t326 - Ifges(6,5) * t336 - Ifges(6,6) * t335 - Ifges(6,3) * t376 - t361 * t345 + t360 * t346;
t421 = m(4) * t380 - t381 * mrSges(4,1) + t382 * mrSges(4,2) - t395 * t388 + t396 * t389 + t308;
t354 = Ifges(5,4) * t387 + Ifges(5,2) * t386 + Ifges(5,6) * t393;
t355 = Ifges(5,1) * t387 + Ifges(5,4) * t386 + Ifges(5,5) * t393;
t419 = mrSges(5,1) * t330 - mrSges(5,2) * t331 + Ifges(5,5) * t359 + Ifges(5,6) * t358 + Ifges(5,3) * t378 + pkin(4) * t314 + t387 * t354 - t386 * t355 - t422;
t399 = (Ifges(3,5) * t407 + Ifges(3,6) * t408) * qJD(1);
t394 = -qJDD(1) * pkin(1) - t418 * qJ(2) + t427;
t383 = -t407 * t397 + t432;
t374 = Ifges(4,1) * t396 + Ifges(4,4) * t395 + Ifges(4,5) * qJD(3);
t373 = Ifges(4,4) * t396 + Ifges(4,2) * t395 + Ifges(4,6) * qJD(3);
t372 = Ifges(4,5) * t396 + Ifges(4,6) * t395 + Ifges(4,3) * qJD(3);
t353 = Ifges(5,5) * t387 + Ifges(5,6) * t386 + Ifges(5,3) * t393;
t344 = Ifges(6,5) * t361 + Ifges(6,6) * t360 + Ifges(6,3) * t391;
t316 = mrSges(6,2) * t329 - mrSges(6,3) * t325 + Ifges(6,1) * t336 + Ifges(6,4) * t335 + Ifges(6,5) * t376 + t360 * t344 - t391 * t345;
t315 = -mrSges(6,1) * t329 + mrSges(6,3) * t326 + Ifges(6,4) * t336 + Ifges(6,2) * t335 + Ifges(6,6) * t376 - t361 * t344 + t391 * t346;
t304 = t439 * t418 * mrSges(3,3) + m(3) * t394 + t426 * qJDD(1) + t421;
t303 = mrSges(5,2) * t339 - mrSges(5,3) * t330 + Ifges(5,1) * t359 + Ifges(5,4) * t358 + Ifges(5,5) * t378 - pkin(8) * t314 - t409 * t315 + t413 * t316 + t386 * t353 - t393 * t354;
t302 = -mrSges(5,1) * t339 + mrSges(5,3) * t331 + Ifges(5,4) * t359 + Ifges(5,2) * t358 + Ifges(5,6) * t378 - pkin(4) * t423 + pkin(8) * t429 + t413 * t315 + t409 * t316 - t387 * t353 + t393 * t355;
t301 = -mrSges(4,1) * t380 + mrSges(4,3) * t350 + Ifges(4,4) * t382 + Ifges(4,2) * t381 + Ifges(4,6) * qJDD(3) - pkin(3) * t308 + qJD(3) * t374 - t396 * t372 - t419;
t300 = mrSges(4,2) * t380 - mrSges(4,3) * t349 + Ifges(4,1) * t382 + Ifges(4,4) * t381 + Ifges(4,5) * qJDD(3) - pkin(7) * t308 - qJD(3) * t373 - t410 * t302 + t414 * t303 + t395 * t372;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t433 - mrSges(2,2) * t428 + t407 * (t399 * t437 + mrSges(3,2) * t394 - mrSges(3,3) * t383 + t415 * t300 - t411 * t301 - pkin(6) * t440 + (Ifges(3,1) * t407 + Ifges(3,4) * t408) * qJDD(1)) + t408 * (-t399 * t438 - mrSges(3,1) * t394 + mrSges(3,3) * t384 + t411 * t300 + t415 * t301 - pkin(2) * t421 + pkin(6) * t431 + (Ifges(3,4) * t407 + Ifges(3,2) * t408) * qJDD(1)) - pkin(1) * t304 + qJ(2) * ((m(3) * t384 + t424 * t408 + t431) * t408 + (-m(3) * t383 + t424 * t407 - t440) * t407); t304; mrSges(4,1) * t349 - mrSges(4,2) * t350 + Ifges(4,5) * t382 + Ifges(4,6) * t381 + Ifges(4,3) * qJDD(3) + pkin(3) * t420 + pkin(7) * t430 + t414 * t302 + t410 * t303 + t396 * t373 - t395 * t374; t419; -t422;];
tauJ = t1;
