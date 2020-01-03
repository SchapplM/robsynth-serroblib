% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:25
% EndTime: 2019-12-31 21:56:28
% DurationCPUTime: 1.78s
% Computational Cost: add. (13428->252), mult. (27099->313), div. (0->0), fcn. (18198->8), ass. (0->102)
t441 = Ifges(5,1) + Ifges(6,1);
t434 = Ifges(5,4) - Ifges(6,5);
t433 = -Ifges(5,5) - Ifges(6,4);
t440 = Ifges(5,2) + Ifges(6,3);
t432 = Ifges(5,6) - Ifges(6,6);
t439 = -Ifges(5,3) - Ifges(6,2);
t404 = sin(qJ(3));
t405 = sin(qJ(2));
t407 = cos(qJ(3));
t408 = cos(qJ(2));
t384 = (t405 * t404 - t408 * t407) * qJD(1);
t423 = qJD(1) * qJD(2);
t390 = qJDD(1) * t405 + t408 * t423;
t391 = qJDD(1) * t408 - t405 * t423;
t361 = -qJD(3) * t384 + t390 * t407 + t391 * t404;
t385 = (t408 * t404 + t405 * t407) * qJD(1);
t401 = qJD(2) + qJD(3);
t403 = sin(qJ(4));
t437 = cos(qJ(4));
t372 = t385 * t403 - t401 * t437;
t400 = qJDD(2) + qJDD(3);
t334 = -t372 * qJD(4) + t361 * t437 + t403 * t400;
t373 = t385 * t437 + t403 * t401;
t349 = mrSges(6,1) * t372 - mrSges(6,3) * t373;
t360 = -qJD(3) * t385 - t390 * t404 + t391 * t407;
t425 = qJD(1) * t405;
t394 = qJD(2) * pkin(2) - pkin(7) * t425;
t402 = t408 ^ 2;
t410 = qJD(1) ^ 2;
t406 = sin(qJ(1));
t409 = cos(qJ(1));
t421 = g(1) * t406 - t409 * g(2);
t415 = -qJDD(1) * pkin(1) - t421;
t362 = -pkin(2) * t391 + t394 * t425 + (-pkin(7) * t402 - pkin(6)) * t410 + t415;
t324 = (t384 * t401 - t361) * pkin(8) + (t385 * t401 - t360) * pkin(3) + t362;
t417 = -g(1) * t409 - g(2) * t406;
t387 = -pkin(1) * t410 + qJDD(1) * pkin(6) + t417;
t430 = t387 * t405;
t436 = pkin(2) * t410;
t352 = qJDD(2) * pkin(2) - pkin(7) * t390 - t430 + (pkin(7) * t423 + t405 * t436 - g(3)) * t408;
t375 = -g(3) * t405 + t408 * t387;
t353 = pkin(7) * t391 - qJD(2) * t394 - t402 * t436 + t375;
t330 = t404 * t352 + t407 * t353;
t371 = pkin(3) * t384 - pkin(8) * t385;
t399 = t401 ^ 2;
t327 = -pkin(3) * t399 + pkin(8) * t400 - t371 * t384 + t330;
t321 = t324 * t437 - t403 * t327;
t348 = pkin(4) * t372 - qJ(5) * t373;
t359 = qJDD(4) - t360;
t380 = qJD(4) + t384;
t379 = t380 ^ 2;
t319 = -t359 * pkin(4) - t379 * qJ(5) + t373 * t348 + qJDD(5) - t321;
t363 = -mrSges(6,2) * t372 + mrSges(6,3) * t380;
t418 = -m(6) * t319 + t359 * mrSges(6,1) + t380 * t363;
t315 = mrSges(6,2) * t334 + t349 * t373 - t418;
t322 = t403 * t324 + t437 * t327;
t318 = -pkin(4) * t379 + qJ(5) * t359 + 0.2e1 * qJD(5) * t380 - t348 * t372 + t322;
t333 = qJD(4) * t373 + t361 * t403 - t400 * t437;
t366 = -mrSges(6,1) * t380 + mrSges(6,2) * t373;
t422 = m(6) * t318 + t359 * mrSges(6,3) + t380 * t366;
t427 = t434 * t372 - t441 * t373 + t433 * t380;
t428 = t440 * t372 - t434 * t373 - t432 * t380;
t438 = -t432 * t333 - t433 * t334 - t439 * t359 - t427 * t372 - t428 * t373 + mrSges(5,1) * t321 - mrSges(6,1) * t319 - mrSges(5,2) * t322 + mrSges(6,3) * t318 - pkin(4) * t315 + qJ(5) * (-mrSges(6,2) * t333 - t349 * t372 + t422);
t435 = -mrSges(5,3) - mrSges(6,2);
t370 = mrSges(4,1) * t384 + mrSges(4,2) * t385;
t377 = mrSges(4,1) * t401 - mrSges(4,3) * t385;
t365 = mrSges(5,1) * t380 - mrSges(5,3) * t373;
t426 = -mrSges(5,1) * t372 - mrSges(5,2) * t373 - t349;
t310 = m(5) * t322 - mrSges(5,2) * t359 + t435 * t333 - t365 * t380 + t426 * t372 + t422;
t364 = -mrSges(5,2) * t380 - mrSges(5,3) * t372;
t312 = m(5) * t321 + mrSges(5,1) * t359 + t435 * t334 + t364 * t380 + t426 * t373 + t418;
t419 = t437 * t310 - t312 * t403;
t299 = m(4) * t330 - mrSges(4,2) * t400 + mrSges(4,3) * t360 - t370 * t384 - t377 * t401 + t419;
t329 = t352 * t407 - t404 * t353;
t376 = -mrSges(4,2) * t401 - mrSges(4,3) * t384;
t326 = -pkin(3) * t400 - pkin(8) * t399 + t385 * t371 - t329;
t320 = -0.2e1 * qJD(5) * t373 + (t372 * t380 - t334) * qJ(5) + (t373 * t380 + t333) * pkin(4) + t326;
t316 = m(6) * t320 + mrSges(6,1) * t333 - t334 * mrSges(6,3) + t363 * t372 - t373 * t366;
t411 = -m(5) * t326 - t333 * mrSges(5,1) - mrSges(5,2) * t334 - t372 * t364 - t365 * t373 - t316;
t307 = m(4) * t329 + mrSges(4,1) * t400 - mrSges(4,3) * t361 - t370 * t385 + t376 * t401 + t411;
t296 = t404 * t299 + t407 * t307;
t305 = t403 * t310 + t437 * t312;
t429 = t432 * t372 + t433 * t373 + t439 * t380;
t424 = qJD(1) * t408;
t420 = t407 * t299 - t307 * t404;
t301 = -mrSges(5,1) * t326 - mrSges(6,1) * t320 + mrSges(6,2) * t318 + mrSges(5,3) * t322 - pkin(4) * t316 - t440 * t333 + t434 * t334 + t432 * t359 + t429 * t373 - t427 * t380;
t303 = mrSges(5,2) * t326 + mrSges(6,2) * t319 - mrSges(5,3) * t321 - mrSges(6,3) * t320 - qJ(5) * t316 - t434 * t333 + t441 * t334 - t433 * t359 + t429 * t372 + t428 * t380;
t368 = Ifges(4,4) * t385 - Ifges(4,2) * t384 + Ifges(4,6) * t401;
t369 = Ifges(4,1) * t385 - Ifges(4,4) * t384 + Ifges(4,5) * t401;
t414 = mrSges(4,1) * t329 - mrSges(4,2) * t330 + Ifges(4,5) * t361 + Ifges(4,6) * t360 + Ifges(4,3) * t400 + pkin(3) * t411 + pkin(8) * t419 + t437 * t301 + t403 * t303 + t385 * t368 + t369 * t384;
t413 = m(4) * t362 - mrSges(4,1) * t360 + mrSges(4,2) * t361 + t376 * t384 + t385 * t377 + t305;
t393 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t424;
t392 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t425;
t389 = (-t408 * mrSges(3,1) + t405 * mrSges(3,2)) * qJD(1);
t386 = -pkin(6) * t410 + t415;
t383 = Ifges(3,5) * qJD(2) + (t405 * Ifges(3,1) + t408 * Ifges(3,4)) * qJD(1);
t382 = Ifges(3,6) * qJD(2) + (t405 * Ifges(3,4) + t408 * Ifges(3,2)) * qJD(1);
t374 = -g(3) * t408 - t430;
t367 = Ifges(4,5) * t385 - Ifges(4,6) * t384 + Ifges(4,3) * t401;
t295 = -mrSges(4,1) * t362 + mrSges(4,3) * t330 + Ifges(4,4) * t361 + Ifges(4,2) * t360 + Ifges(4,6) * t400 - pkin(3) * t305 - t385 * t367 + t401 * t369 - t438;
t294 = mrSges(4,2) * t362 - mrSges(4,3) * t329 + Ifges(4,1) * t361 + Ifges(4,4) * t360 + Ifges(4,5) * t400 - pkin(8) * t305 - t403 * t301 + t303 * t437 - t384 * t367 - t401 * t368;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t421 - mrSges(2,2) * t417 + t405 * (mrSges(3,2) * t386 - mrSges(3,3) * t374 + Ifges(3,1) * t390 + Ifges(3,4) * t391 + Ifges(3,5) * qJDD(2) - pkin(7) * t296 - qJD(2) * t382 + t407 * t294 - t404 * t295) + t408 * (-mrSges(3,1) * t386 + mrSges(3,3) * t375 + Ifges(3,4) * t390 + Ifges(3,2) * t391 + Ifges(3,6) * qJDD(2) - pkin(2) * t413 + pkin(7) * t420 + qJD(2) * t383 + t404 * t294 + t407 * t295) + pkin(1) * (-m(3) * t386 + mrSges(3,1) * t391 - mrSges(3,2) * t390 + (-t392 * t405 + t393 * t408) * qJD(1) - t413) + pkin(6) * (t408 * (m(3) * t375 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t391 - qJD(2) * t392 + t389 * t424 + t420) - t405 * (m(3) * t374 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t390 + qJD(2) * t393 - t389 * t425 + t296)); t414 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t374 - mrSges(3,2) * t375 + Ifges(3,5) * t390 + Ifges(3,6) * t391 + pkin(2) * t296 + (t405 * t382 - t408 * t383) * qJD(1); t414; t438; t315;];
tauJ = t1;
