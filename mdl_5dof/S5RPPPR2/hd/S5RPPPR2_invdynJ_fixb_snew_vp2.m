% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:41
% EndTime: 2020-01-03 11:22:43
% DurationCPUTime: 2.15s
% Computational Cost: add. (8382->231), mult. (23494->327), div. (0->0), fcn. (15854->10), ass. (0->111)
t404 = sin(pkin(7));
t407 = cos(pkin(7));
t412 = qJD(1) ^ 2;
t409 = sin(qJ(1));
t411 = cos(qJ(1));
t429 = -g(2) * t409 + t411 * g(3);
t450 = -pkin(1) * t412 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t429;
t363 = -g(1) * t404 + t450 * t407;
t422 = -t407 * pkin(2) - t404 * qJ(3);
t386 = t422 * qJD(1);
t439 = qJD(1) * t407;
t354 = t386 * t439 + t363;
t403 = sin(pkin(8));
t406 = cos(pkin(8));
t441 = -t411 * g(2) - t409 * g(3);
t416 = -qJ(2) * t412 + qJDD(2) - t441;
t440 = qJD(1) * t404;
t449 = (-pkin(1) + t422) * qJDD(1) + t416 - 0.2e1 * qJD(3) * t440;
t334 = -t403 * t354 + t449 * t406;
t362 = -t407 * g(1) - t450 * t404;
t402 = sin(pkin(9));
t405 = cos(pkin(9));
t444 = t404 * t406;
t415 = t402 * t444 + t407 * t405;
t375 = t415 * qJD(1);
t373 = t415 * qJDD(1);
t448 = 2 * qJD(4);
t447 = Ifges(4,4) * t406;
t401 = t407 ^ 2;
t446 = t401 * t412;
t445 = t403 * t404;
t443 = t407 * t412;
t335 = t406 * t354 + t449 * t403;
t379 = (t403 * pkin(3) - t406 * qJ(4)) * t440;
t434 = t403 * t440;
t436 = qJDD(1) * t407;
t332 = -pkin(3) * t446 - qJ(4) * t436 - t379 * t434 + t335;
t353 = t386 * t440 + qJDD(3) - t362;
t340 = ((-qJDD(1) * t406 - t403 * t443) * qJ(4) + (qJDD(1) * t403 - t406 * t443) * pkin(3)) * t404 + t353;
t328 = t405 * t332 + t402 * t340 - t375 * t448;
t433 = t406 * t440;
t376 = -t402 * t439 + t405 * t433;
t355 = mrSges(5,1) * t375 + mrSges(5,2) * t376;
t361 = mrSges(5,1) * t434 - mrSges(5,3) * t376;
t356 = pkin(4) * t375 - pkin(6) * t376;
t437 = qJDD(1) * t404;
t430 = t403 * t437;
t400 = t404 ^ 2;
t435 = t403 ^ 2 * t400 * t412;
t326 = -pkin(4) * t435 + pkin(6) * t430 - t356 * t375 + t328;
t331 = pkin(3) * t436 - qJ(4) * t446 + t379 * t433 + qJDD(4) - t334;
t374 = (-t407 * t402 + t405 * t444) * qJDD(1);
t329 = (t375 * t434 - t374) * pkin(6) + (t376 * t434 + t373) * pkin(4) + t331;
t408 = sin(qJ(5));
t410 = cos(qJ(5));
t323 = -t326 * t408 + t329 * t410;
t357 = -t376 * t408 + t410 * t434;
t358 = t376 * t410 + t408 * t434;
t342 = -mrSges(6,1) * t357 + mrSges(6,2) * t358;
t344 = qJD(5) * t357 + t374 * t410 + t408 * t430;
t372 = qJD(5) + t375;
t345 = -mrSges(6,2) * t372 + mrSges(6,3) * t357;
t371 = qJDD(5) + t373;
t321 = m(6) * t323 + mrSges(6,1) * t371 - mrSges(6,3) * t344 - t342 * t358 + t345 * t372;
t324 = t326 * t410 + t329 * t408;
t343 = -qJD(5) * t358 - t374 * t408 + t410 * t430;
t346 = mrSges(6,1) * t372 - mrSges(6,3) * t358;
t322 = m(6) * t324 - mrSges(6,2) * t371 + mrSges(6,3) * t343 + t342 * t357 - t346 * t372;
t427 = -t321 * t408 + t410 * t322;
t313 = m(5) * t328 - mrSges(5,3) * t373 - t355 * t375 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t361) * t445 + t427;
t421 = t332 * t402 - t340 * t405;
t327 = -0.2e1 * qJD(4) * t376 - t421;
t360 = -mrSges(5,2) * t434 - mrSges(5,3) * t375;
t325 = -pkin(4) * t430 - pkin(6) * t435 + (t448 + t356) * t376 + t421;
t414 = -m(6) * t325 + t343 * mrSges(6,1) - mrSges(6,2) * t344 + t357 * t345 - t346 * t358;
t319 = m(5) * t327 - mrSges(5,3) * t374 - t355 * t376 + (mrSges(5,1) * qJDD(1) + qJD(1) * t360) * t445 + t414;
t309 = t402 * t313 + t405 * t319;
t428 = t405 * t313 - t319 * t402;
t425 = m(4) * t353 + t309;
t424 = -t407 * mrSges(3,1) + t404 * mrSges(3,2);
t423 = t403 * mrSges(4,1) + t406 * mrSges(4,2);
t380 = t423 * t440;
t418 = -t407 * mrSges(4,1) - mrSges(4,3) * t444;
t384 = t418 * qJD(1);
t417 = t407 * mrSges(4,2) - mrSges(4,3) * t445;
t308 = m(4) * t335 + t417 * qJDD(1) + (-t380 * t445 + t384 * t407) * qJD(1) + t428;
t315 = t321 * t410 + t322 * t408;
t314 = m(5) * t331 + t373 * mrSges(5,1) + mrSges(5,2) * t374 + t375 * t360 + t361 * t376 + t315;
t383 = t417 * qJD(1);
t310 = m(4) * t334 + t418 * qJDD(1) + (-t380 * t444 - t383 * t407) * qJD(1) - t314;
t305 = t308 * t403 + t310 * t406;
t420 = t383 * t403 + t384 * t406;
t419 = -Ifges(4,5) * t406 + Ifges(4,6) * t403 + Ifges(3,4);
t337 = Ifges(6,4) * t358 + Ifges(6,2) * t357 + Ifges(6,6) * t372;
t338 = Ifges(6,1) * t358 + Ifges(6,4) * t357 + Ifges(6,5) * t372;
t413 = mrSges(6,1) * t323 - mrSges(6,2) * t324 + Ifges(6,5) * t344 + Ifges(6,6) * t343 + Ifges(6,3) * t371 + t358 * t337 - t357 * t338;
t388 = (Ifges(3,5) * t404 + Ifges(3,6) * t407) * qJD(1);
t387 = t424 * qJD(1);
t382 = -qJDD(1) * pkin(1) + t416;
t368 = (-Ifges(4,5) * t407 + (Ifges(4,1) * t406 - Ifges(4,4) * t403) * t404) * qJD(1);
t367 = (-Ifges(4,6) * t407 + (-Ifges(4,2) * t403 + t447) * t404) * qJD(1);
t349 = Ifges(5,1) * t376 - Ifges(5,4) * t375 + Ifges(5,5) * t434;
t348 = Ifges(5,4) * t376 - Ifges(5,2) * t375 + Ifges(5,6) * t434;
t347 = Ifges(5,5) * t376 - Ifges(5,6) * t375 + Ifges(5,3) * t434;
t336 = Ifges(6,5) * t358 + Ifges(6,6) * t357 + Ifges(6,3) * t372;
t317 = mrSges(6,2) * t325 - mrSges(6,3) * t323 + Ifges(6,1) * t344 + Ifges(6,4) * t343 + Ifges(6,5) * t371 + t336 * t357 - t337 * t372;
t316 = -mrSges(6,1) * t325 + mrSges(6,3) * t324 + Ifges(6,4) * t344 + Ifges(6,2) * t343 + Ifges(6,6) * t371 - t336 * t358 + t338 * t372;
t307 = -mrSges(5,1) * t331 + mrSges(5,3) * t328 + Ifges(5,4) * t374 - Ifges(5,2) * t373 - pkin(4) * t315 - t376 * t347 + (Ifges(5,6) * qJDD(1) + qJD(1) * t349) * t445 - t413;
t306 = mrSges(5,2) * t331 - mrSges(5,3) * t327 + Ifges(5,1) * t374 - Ifges(5,4) * t373 - pkin(6) * t315 - t316 * t408 + t317 * t410 - t347 * t375 + (Ifges(5,5) * qJDD(1) - qJD(1) * t348) * t445;
t304 = m(3) * t382 + t424 * qJDD(1) + (-t400 - t401) * t412 * mrSges(3,3) + t305;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t441 - mrSges(2,2) * t429 + t404 * (t388 * t439 + mrSges(3,2) * t382 - mrSges(3,3) * t362 + t406 * (mrSges(4,2) * t353 - mrSges(4,3) * t334 - qJ(4) * t309 + t306 * t405 - t307 * t402 + t367 * t439) - t403 * (-mrSges(4,1) * t353 - mrSges(5,1) * t327 + mrSges(5,2) * t328 + mrSges(4,3) * t335 - Ifges(5,5) * t374 + Ifges(5,6) * t373 - pkin(3) * t309 - pkin(4) * t414 - pkin(6) * t427 - t410 * t316 - t408 * t317 - t376 * t348 - t375 * t349 - t368 * t439) - qJ(3) * t305 + (t407 * t419 + (Ifges(4,1) * t406 ^ 2 + Ifges(3,1) + (-0.2e1 * t447 + (Ifges(4,2) + Ifges(5,3)) * t403) * t403) * t404) * qJDD(1)) + t407 * (-mrSges(3,1) * t382 + mrSges(3,3) * t363 - mrSges(4,1) * t334 + mrSges(4,2) * t335 - t402 * t306 - t405 * t307 + pkin(3) * t314 - qJ(4) * t428 - pkin(2) * t305 + (Ifges(3,2) + Ifges(4,3)) * t436 + (t419 * qJDD(1) + (-t367 * t406 - t368 * t403 - t388) * qJD(1)) * t404) - pkin(1) * t304 + qJ(2) * ((m(3) * t363 + t406 * t308 - t403 * t310 + (qJDD(1) * mrSges(3,3) + qJD(1) * t387) * t407) * t407 + (-m(3) * t362 + (mrSges(3,3) + t423) * t437 + (t387 + t420) * t440 + t425) * t404); t304; (t420 * qJD(1) + t423 * qJDD(1)) * t404 + t425; t314; t413;];
tauJ = t1;
