% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:09
% EndTime: 2019-12-31 20:55:12
% DurationCPUTime: 2.18s
% Computational Cost: add. (13866->252), mult. (31155->314), div. (0->0), fcn. (21194->8), ass. (0->101)
t453 = Ifges(5,1) + Ifges(6,1);
t443 = Ifges(5,4) - Ifges(6,5);
t451 = Ifges(6,4) + Ifges(5,5);
t452 = Ifges(5,2) + Ifges(6,3);
t449 = Ifges(5,6) - Ifges(6,6);
t450 = -Ifges(6,2) - Ifges(5,3);
t418 = sin(qJ(3));
t419 = sin(qJ(2));
t421 = cos(qJ(3));
t422 = cos(qJ(2));
t396 = (-t419 * t418 + t422 * t421) * qJD(1);
t397 = (t422 * t418 + t419 * t421) * qJD(1);
t417 = sin(pkin(8));
t442 = cos(pkin(8));
t382 = -t442 * t396 + t397 * t417;
t383 = t417 * t396 + t442 * t397;
t415 = qJD(2) + qJD(3);
t448 = t452 * t382 - t443 * t383 - t449 * t415;
t447 = -t443 * t382 + t453 * t383 + t451 * t415;
t446 = -2 * qJD(4);
t424 = qJD(1) ^ 2;
t445 = pkin(2) * t424;
t444 = -mrSges(5,3) - mrSges(6,2);
t420 = sin(qJ(1));
t423 = cos(qJ(1));
t430 = -g(1) * t423 - g(2) * t420;
t399 = -pkin(1) * t424 + qJDD(1) * pkin(6) + t430;
t441 = t399 * t419;
t436 = qJD(1) * qJD(2);
t402 = qJDD(1) * t419 + t422 * t436;
t363 = qJDD(2) * pkin(2) - pkin(7) * t402 - t441 + (pkin(7) * t436 + t419 * t445 - g(3)) * t422;
t386 = -g(3) * t419 + t422 * t399;
t403 = qJDD(1) * t422 - t419 * t436;
t438 = qJD(1) * t419;
t406 = qJD(2) * pkin(2) - pkin(7) * t438;
t416 = t422 ^ 2;
t364 = pkin(7) * t403 - qJD(2) * t406 - t416 * t445 + t386;
t338 = t421 * t363 - t364 * t418;
t371 = qJD(3) * t396 + t402 * t421 + t403 * t418;
t414 = qJDD(2) + qJDD(3);
t330 = (t396 * t415 - t371) * qJ(4) + (t396 * t397 + t414) * pkin(3) + t338;
t339 = t418 * t363 + t421 * t364;
t370 = -qJD(3) * t397 - t402 * t418 + t403 * t421;
t388 = pkin(3) * t415 - qJ(4) * t397;
t392 = t396 ^ 2;
t332 = -pkin(3) * t392 + qJ(4) * t370 - t388 * t415 + t339;
t328 = t417 * t330 + t442 * t332 + t382 * t446;
t346 = -t442 * t370 + t371 * t417;
t374 = mrSges(5,1) * t415 - mrSges(5,3) * t383;
t356 = pkin(4) * t382 - qJ(5) * t383;
t413 = t415 ^ 2;
t322 = -pkin(4) * t413 + qJ(5) * t414 + 0.2e1 * qJD(5) * t415 - t356 * t382 + t328;
t375 = -mrSges(6,1) * t415 + mrSges(6,2) * t383;
t435 = m(6) * t322 + t414 * mrSges(6,3) + t415 * t375;
t357 = mrSges(6,1) * t382 - mrSges(6,3) * t383;
t439 = -mrSges(5,1) * t382 - mrSges(5,2) * t383 - t357;
t312 = m(5) * t328 - mrSges(5,2) * t414 + t444 * t346 - t374 * t415 + t439 * t382 + t435;
t428 = t442 * t330 - t417 * t332;
t327 = t383 * t446 + t428;
t347 = t417 * t370 + t442 * t371;
t373 = -mrSges(5,2) * t415 - mrSges(5,3) * t382;
t323 = -t414 * pkin(4) - t413 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t356) * t383 - t428;
t376 = -mrSges(6,2) * t382 + mrSges(6,3) * t415;
t431 = -m(6) * t323 + t414 * mrSges(6,1) + t415 * t376;
t314 = m(5) * t327 + mrSges(5,1) * t414 + t444 * t347 + t373 * t415 + t439 * t383 + t431;
t309 = t417 * t312 + t442 * t314;
t384 = -mrSges(4,1) * t396 + mrSges(4,2) * t397;
t387 = -mrSges(4,2) * t415 + mrSges(4,3) * t396;
t304 = m(4) * t338 + mrSges(4,1) * t414 - mrSges(4,3) * t371 - t384 * t397 + t387 * t415 + t309;
t389 = mrSges(4,1) * t415 - mrSges(4,3) * t397;
t432 = t442 * t312 - t314 * t417;
t305 = m(4) * t339 - mrSges(4,2) * t414 + mrSges(4,3) * t370 + t384 * t396 - t389 * t415 + t432;
t300 = t421 * t304 + t418 * t305;
t440 = t449 * t382 - t451 * t383 + t450 * t415;
t437 = qJD(1) * t422;
t434 = g(1) * t420 - t423 * g(2);
t433 = -t304 * t418 + t421 * t305;
t429 = -qJDD(1) * pkin(1) - t434;
t372 = -pkin(2) * t403 + t406 * t438 + (-pkin(7) * t416 - pkin(6)) * t424 + t429;
t334 = -pkin(3) * t370 - qJ(4) * t392 + t397 * t388 + qJDD(4) + t372;
t325 = -0.2e1 * qJD(5) * t383 + (t382 * t415 - t347) * qJ(5) + (t383 * t415 + t346) * pkin(4) + t334;
t427 = -m(6) * t325 - t346 * mrSges(6,1) + t347 * mrSges(6,3) + t383 * t375 - t382 * t376;
t315 = m(5) * t334 + t346 * mrSges(5,1) + t347 * mrSges(5,2) + t382 * t373 + t383 * t374 - t427;
t426 = m(4) * t372 - t370 * mrSges(4,1) + t371 * mrSges(4,2) - t396 * t387 + t397 * t389 + t315;
t318 = mrSges(6,2) * t347 + t357 * t383 - t431;
t379 = Ifges(4,4) * t397 + Ifges(4,2) * t396 + Ifges(4,6) * t415;
t380 = Ifges(4,1) * t397 + Ifges(4,4) * t396 + Ifges(4,5) * t415;
t425 = mrSges(4,1) * t338 + mrSges(5,1) * t327 - mrSges(6,1) * t323 - mrSges(4,2) * t339 - mrSges(5,2) * t328 + mrSges(6,3) * t322 + Ifges(4,5) * t371 + Ifges(4,6) * t370 + pkin(3) * t309 - pkin(4) * t318 + qJ(5) * t435 + t397 * t379 - t396 * t380 - t448 * t383 + (-qJ(5) * t357 + t447) * t382 + t451 * t347 + (-qJ(5) * mrSges(6,2) - t449) * t346 + (Ifges(4,3) - t450) * t414;
t405 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t437;
t404 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t438;
t401 = (-t422 * mrSges(3,1) + t419 * mrSges(3,2)) * qJD(1);
t398 = -pkin(6) * t424 + t429;
t395 = Ifges(3,5) * qJD(2) + (t419 * Ifges(3,1) + t422 * Ifges(3,4)) * qJD(1);
t394 = Ifges(3,6) * qJD(2) + (t419 * Ifges(3,4) + t422 * Ifges(3,2)) * qJD(1);
t385 = -g(3) * t422 - t441;
t378 = Ifges(4,5) * t397 + Ifges(4,6) * t396 + Ifges(4,3) * t415;
t307 = mrSges(5,2) * t334 + mrSges(6,2) * t323 - mrSges(5,3) * t327 - mrSges(6,3) * t325 + qJ(5) * t427 - t443 * t346 + t453 * t347 + t440 * t382 + t451 * t414 + t448 * t415;
t306 = -mrSges(5,1) * t334 - mrSges(6,1) * t325 + mrSges(6,2) * t322 + mrSges(5,3) * t328 + pkin(4) * t427 - t452 * t346 + t443 * t347 + t440 * t383 + t449 * t414 + t447 * t415;
t299 = mrSges(4,2) * t372 - mrSges(4,3) * t338 + Ifges(4,1) * t371 + Ifges(4,4) * t370 + Ifges(4,5) * t414 - qJ(4) * t309 - t417 * t306 + t442 * t307 + t396 * t378 - t415 * t379;
t298 = -mrSges(4,1) * t372 + mrSges(4,3) * t339 + Ifges(4,4) * t371 + Ifges(4,2) * t370 + Ifges(4,6) * t414 - pkin(3) * t315 + qJ(4) * t432 + t442 * t306 + t417 * t307 - t397 * t378 + t415 * t380;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t434 - mrSges(2,2) * t430 + t419 * (mrSges(3,2) * t398 - mrSges(3,3) * t385 + Ifges(3,1) * t402 + Ifges(3,4) * t403 + Ifges(3,5) * qJDD(2) - pkin(7) * t300 - qJD(2) * t394 - t418 * t298 + t421 * t299) + t422 * (-mrSges(3,1) * t398 + mrSges(3,3) * t386 + Ifges(3,4) * t402 + Ifges(3,2) * t403 + Ifges(3,6) * qJDD(2) - pkin(2) * t426 + pkin(7) * t433 + qJD(2) * t395 + t421 * t298 + t418 * t299) + pkin(1) * ((-t419 * t404 + t422 * t405) * qJD(1) - t426 - t402 * mrSges(3,2) + t403 * mrSges(3,1) - m(3) * t398) + pkin(6) * (t422 * (m(3) * t386 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t403 - qJD(2) * t404 + t401 * t437 + t433) - t419 * (m(3) * t385 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t402 + qJD(2) * t405 - t401 * t438 + t300)); (t419 * t394 - t422 * t395) * qJD(1) + Ifges(3,3) * qJDD(2) + t425 + Ifges(3,5) * t402 + Ifges(3,6) * t403 + mrSges(3,1) * t385 - mrSges(3,2) * t386 + pkin(2) * t300; t425; t315; t318;];
tauJ = t1;
