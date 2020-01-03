% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:57:02
% EndTime: 2019-12-31 19:57:05
% DurationCPUTime: 1.62s
% Computational Cost: add. (10454->253), mult. (23659->313), div. (0->0), fcn. (15605->8), ass. (0->101)
t442 = -2 * qJD(3);
t441 = Ifges(5,1) + Ifges(6,1);
t435 = Ifges(5,4) + Ifges(6,4);
t434 = Ifges(5,5) + Ifges(6,5);
t440 = Ifges(5,2) + Ifges(6,2);
t433 = Ifges(5,6) + Ifges(6,6);
t439 = Ifges(5,3) + Ifges(6,3);
t404 = sin(qJ(2));
t407 = cos(qJ(2));
t423 = qJD(1) * qJD(2);
t394 = qJDD(1) * t404 + t407 * t423;
t410 = qJD(1) ^ 2;
t405 = sin(qJ(1));
t408 = cos(qJ(1));
t415 = -g(1) * t408 - g(2) * t405;
t391 = -pkin(1) * t410 + qJDD(1) * pkin(6) + t415;
t431 = t404 * t391;
t437 = pkin(2) * t410;
t354 = qJDD(2) * pkin(2) - t394 * qJ(3) - t431 + (qJ(3) * t423 + t404 * t437 - g(3)) * t407;
t378 = -g(3) * t404 + t407 * t391;
t395 = qJDD(1) * t407 - t404 * t423;
t426 = qJD(1) * t404;
t396 = qJD(2) * pkin(2) - qJ(3) * t426;
t400 = t407 ^ 2;
t355 = qJ(3) * t395 - qJD(2) * t396 - t400 * t437 + t378;
t401 = sin(pkin(8));
t402 = cos(pkin(8));
t386 = (t401 * t407 + t402 * t404) * qJD(1);
t334 = t354 * t402 - t401 * t355 + t386 * t442;
t385 = (t401 * t404 - t402 * t407) * qJD(1);
t373 = t394 * t402 + t395 * t401;
t403 = sin(qJ(4));
t406 = cos(qJ(4));
t375 = qJD(2) * t406 - t386 * t403;
t349 = qJD(4) * t375 + qJDD(2) * t403 + t373 * t406;
t376 = qJD(2) * t403 + t386 * t406;
t356 = -mrSges(6,1) * t375 + mrSges(6,2) * t376;
t335 = t401 * t354 + t402 * t355 + t385 * t442;
t369 = pkin(3) * t385 - pkin(7) * t386;
t409 = qJD(2) ^ 2;
t330 = -pkin(3) * t409 + qJDD(2) * pkin(7) - t369 * t385 + t335;
t420 = t405 * g(1) - t408 * g(2);
t413 = -qJDD(1) * pkin(1) - t420;
t359 = -t395 * pkin(2) + qJDD(3) + t396 * t426 + (-qJ(3) * t400 - pkin(6)) * t410 + t413;
t372 = -t394 * t401 + t395 * t402;
t333 = (qJD(2) * t385 - t373) * pkin(7) + (qJD(2) * t386 - t372) * pkin(3) + t359;
t326 = -t403 * t330 + t406 * t333;
t371 = qJDD(4) - t372;
t384 = qJD(4) + t385;
t322 = -0.2e1 * qJD(5) * t376 + (t375 * t384 - t349) * qJ(5) + (t375 * t376 + t371) * pkin(4) + t326;
t360 = -mrSges(6,2) * t384 + mrSges(6,3) * t375;
t422 = m(6) * t322 + t371 * mrSges(6,1) + t384 * t360;
t319 = -t349 * mrSges(6,3) - t376 * t356 + t422;
t327 = t406 * t330 + t403 * t333;
t348 = -qJD(4) * t376 + qJDD(2) * t406 - t373 * t403;
t362 = pkin(4) * t384 - qJ(5) * t376;
t374 = t375 ^ 2;
t324 = -pkin(4) * t374 + qJ(5) * t348 + 0.2e1 * qJD(5) * t375 - t362 * t384 + t327;
t428 = t435 * t375 + t441 * t376 + t434 * t384;
t429 = -t440 * t375 - t435 * t376 - t433 * t384;
t438 = mrSges(5,1) * t326 + mrSges(6,1) * t322 - mrSges(5,2) * t327 - mrSges(6,2) * t324 + pkin(4) * t319 + t433 * t348 + t434 * t349 + t439 * t371 - t428 * t375 - t429 * t376;
t436 = -mrSges(5,2) - mrSges(6,2);
t368 = mrSges(4,1) * t385 + mrSges(4,2) * t386;
t380 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t386;
t357 = -mrSges(5,1) * t375 + mrSges(5,2) * t376;
t361 = -mrSges(5,2) * t384 + mrSges(5,3) * t375;
t313 = m(5) * t326 + t371 * mrSges(5,1) + t384 * t361 + (-t356 - t357) * t376 + (-mrSges(5,3) - mrSges(6,3)) * t349 + t422;
t421 = m(6) * t324 + t348 * mrSges(6,3) + t375 * t356;
t363 = mrSges(6,1) * t384 - mrSges(6,3) * t376;
t427 = -mrSges(5,1) * t384 + mrSges(5,3) * t376 - t363;
t316 = m(5) * t327 + t348 * mrSges(5,3) + t375 * t357 + t436 * t371 + t427 * t384 + t421;
t418 = -t313 * t403 + t406 * t316;
t308 = m(4) * t335 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t372 - qJD(2) * t380 - t368 * t385 + t418;
t379 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t385;
t329 = -qJDD(2) * pkin(3) - pkin(7) * t409 + t386 * t369 - t334;
t325 = -pkin(4) * t348 - qJ(5) * t374 + t362 * t376 + qJDD(5) + t329;
t416 = -m(6) * t325 + t348 * mrSges(6,1) + t375 * t360;
t411 = -m(5) * t329 + t348 * mrSges(5,1) + t436 * t349 + t375 * t361 + t427 * t376 + t416;
t318 = m(4) * t334 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t373 + qJD(2) * t379 - t368 * t386 + t411;
t304 = t401 * t308 + t402 * t318;
t311 = t406 * t313 + t403 * t316;
t430 = -t433 * t375 - t434 * t376 - t439 * t384;
t425 = qJD(1) * t407;
t419 = t402 * t308 - t318 * t401;
t310 = m(4) * t359 - mrSges(4,1) * t372 + mrSges(4,2) * t373 + t379 * t385 + t380 * t386 + t311;
t398 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t425;
t397 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t426;
t393 = (-mrSges(3,1) * t407 + mrSges(3,2) * t404) * qJD(1);
t390 = -t410 * pkin(6) + t413;
t389 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t404 + Ifges(3,4) * t407) * qJD(1);
t388 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t404 + Ifges(3,2) * t407) * qJD(1);
t377 = -t407 * g(3) - t431;
t367 = Ifges(4,1) * t386 - Ifges(4,4) * t385 + Ifges(4,5) * qJD(2);
t366 = Ifges(4,4) * t386 - Ifges(4,2) * t385 + Ifges(4,6) * qJD(2);
t365 = Ifges(4,5) * t386 - Ifges(4,6) * t385 + Ifges(4,3) * qJD(2);
t320 = t349 * mrSges(6,2) + t376 * t363 - t416;
t309 = mrSges(5,2) * t329 + mrSges(6,2) * t325 - mrSges(5,3) * t326 - mrSges(6,3) * t322 - qJ(5) * t319 + t435 * t348 + t441 * t349 + t434 * t371 - t430 * t375 + t429 * t384;
t305 = -mrSges(5,1) * t329 + mrSges(5,3) * t327 - mrSges(6,1) * t325 + mrSges(6,3) * t324 - pkin(4) * t320 + qJ(5) * t421 + (-qJ(5) * t363 + t428) * t384 + t430 * t376 + (-mrSges(6,2) * qJ(5) + t433) * t371 + t435 * t349 + t440 * t348;
t303 = -mrSges(4,1) * t359 + mrSges(4,3) * t335 + Ifges(4,4) * t373 + Ifges(4,2) * t372 + Ifges(4,6) * qJDD(2) - pkin(3) * t311 + qJD(2) * t367 - t386 * t365 - t438;
t302 = mrSges(4,2) * t359 - mrSges(4,3) * t334 + Ifges(4,1) * t373 + Ifges(4,4) * t372 + Ifges(4,5) * qJDD(2) - pkin(7) * t311 - qJD(2) * t366 - t305 * t403 + t309 * t406 - t365 * t385;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t420 - mrSges(2,2) * t415 + t404 * (mrSges(3,2) * t390 - mrSges(3,3) * t377 + Ifges(3,1) * t394 + Ifges(3,4) * t395 + Ifges(3,5) * qJDD(2) - qJ(3) * t304 - qJD(2) * t388 + t402 * t302 - t401 * t303) + t407 * (-mrSges(3,1) * t390 + mrSges(3,3) * t378 + Ifges(3,4) * t394 + Ifges(3,2) * t395 + Ifges(3,6) * qJDD(2) - pkin(2) * t310 + qJ(3) * t419 + qJD(2) * t389 + t401 * t302 + t402 * t303) + pkin(1) * (-m(3) * t390 + mrSges(3,1) * t395 - mrSges(3,2) * t394 + (-t397 * t404 + t398 * t407) * qJD(1) - t310) + pkin(6) * (t407 * (m(3) * t378 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t395 - qJD(2) * t397 + t393 * t425 + t419) - t404 * (m(3) * t377 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t394 + qJD(2) * t398 - t393 * t426 + t304)); Ifges(3,5) * t394 + Ifges(3,6) * t395 + mrSges(3,1) * t377 - mrSges(3,2) * t378 + Ifges(4,5) * t373 + Ifges(4,6) * t372 + t386 * t366 + t385 * t367 + mrSges(4,1) * t334 - mrSges(4,2) * t335 + t403 * t309 + t406 * t305 + pkin(3) * t411 + pkin(7) * t418 + pkin(2) * t304 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t404 * t388 - t407 * t389) * qJD(1); t310; t438; t320;];
tauJ = t1;
