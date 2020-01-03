% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:51:08
% EndTime: 2019-12-31 18:51:11
% DurationCPUTime: 1.41s
% Computational Cost: add. (8941->223), mult. (21038->272), div. (0->0), fcn. (14739->8), ass. (0->99)
t437 = Ifges(5,1) + Ifges(6,1);
t431 = Ifges(5,4) + Ifges(6,4);
t430 = Ifges(5,5) + Ifges(6,5);
t436 = Ifges(5,2) + Ifges(6,2);
t429 = Ifges(5,6) + Ifges(6,6);
t435 = Ifges(5,3) + Ifges(6,3);
t400 = qJD(1) ^ 2;
t391 = sin(pkin(8));
t392 = cos(pkin(8));
t394 = sin(qJ(3));
t397 = cos(qJ(3));
t405 = t391 * t394 - t392 * t397;
t381 = t405 * qJD(1);
t406 = t391 * t397 + t392 * t394;
t382 = t406 * qJD(1);
t418 = qJD(3) * t382;
t369 = -t405 * qJDD(1) - t418;
t419 = qJD(3) * t381;
t370 = t406 * qJDD(1) - t419;
t393 = sin(qJ(4));
t396 = cos(qJ(4));
t374 = qJD(3) * t396 - t382 * t393;
t346 = qJD(4) * t374 + qJDD(3) * t393 + t370 * t396;
t375 = qJD(3) * t393 + t382 * t396;
t348 = -mrSges(6,1) * t374 + mrSges(6,2) * t375;
t395 = sin(qJ(1));
t398 = cos(qJ(1));
t409 = -g(1) * t398 - g(2) * t395;
t383 = -pkin(1) * t400 + qJDD(1) * qJ(2) + t409;
t417 = qJD(1) * qJD(2);
t413 = -g(3) * t392 - 0.2e1 * t391 * t417;
t427 = pkin(6) * qJDD(1);
t433 = pkin(2) * t400;
t359 = (t392 * t433 - t383 - t427) * t391 + t413;
t372 = -g(3) * t391 + (t383 + 0.2e1 * t417) * t392;
t390 = t392 ^ 2;
t360 = -t390 * t433 + t392 * t427 + t372;
t332 = t394 * t359 + t397 * t360;
t367 = pkin(3) * t381 - pkin(7) * t382;
t399 = qJD(3) ^ 2;
t327 = -pkin(3) * t399 + qJDD(3) * pkin(7) - t367 * t381 + t332;
t414 = g(1) * t395 - t398 * g(2);
t408 = qJDD(2) - t414;
t421 = -t391 ^ 2 - t390;
t368 = (-pkin(2) * t392 - pkin(1)) * qJDD(1) + (t421 * pkin(6) - qJ(2)) * t400 + t408;
t330 = (-t370 + t419) * pkin(7) + (-t369 + t418) * pkin(3) + t368;
t323 = -t327 * t393 + t396 * t330;
t366 = qJDD(4) - t369;
t379 = qJD(4) + t381;
t319 = -0.2e1 * qJD(5) * t375 + (t374 * t379 - t346) * qJ(5) + (t374 * t375 + t366) * pkin(4) + t323;
t352 = -mrSges(6,2) * t379 + mrSges(6,3) * t374;
t416 = m(6) * t319 + t366 * mrSges(6,1) + t379 * t352;
t316 = -mrSges(6,3) * t346 - t348 * t375 + t416;
t324 = t396 * t327 + t393 * t330;
t345 = -qJD(4) * t375 + qJDD(3) * t396 - t370 * t393;
t354 = pkin(4) * t379 - qJ(5) * t375;
t373 = t374 ^ 2;
t321 = -pkin(4) * t373 + qJ(5) * t345 + 0.2e1 * qJD(5) * t374 - t354 * t379 + t324;
t423 = t431 * t374 + t437 * t375 + t430 * t379;
t424 = -t436 * t374 - t431 * t375 - t429 * t379;
t434 = mrSges(5,1) * t323 + mrSges(6,1) * t319 - mrSges(5,2) * t324 - mrSges(6,2) * t321 + pkin(4) * t316 + t429 * t345 + t430 * t346 + t435 * t366 - t423 * t374 - t424 * t375;
t432 = -mrSges(5,2) - mrSges(6,2);
t365 = mrSges(4,1) * t381 + mrSges(4,2) * t382;
t377 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t382;
t349 = -mrSges(5,1) * t374 + mrSges(5,2) * t375;
t353 = -mrSges(5,2) * t379 + mrSges(5,3) * t374;
t310 = m(5) * t323 + mrSges(5,1) * t366 + t353 * t379 + (-t348 - t349) * t375 + (-mrSges(5,3) - mrSges(6,3)) * t346 + t416;
t415 = m(6) * t321 + t345 * mrSges(6,3) + t374 * t348;
t355 = mrSges(6,1) * t379 - mrSges(6,3) * t375;
t422 = -mrSges(5,1) * t379 + mrSges(5,3) * t375 - t355;
t313 = m(5) * t324 + mrSges(5,3) * t345 + t349 * t374 + t432 * t366 + t422 * t379 + t415;
t411 = -t310 * t393 + t396 * t313;
t306 = m(4) * t332 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t369 - qJD(3) * t377 - t365 * t381 + t411;
t331 = t359 * t397 - t394 * t360;
t376 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t381;
t326 = -qJDD(3) * pkin(3) - pkin(7) * t399 + t382 * t367 - t331;
t322 = -pkin(4) * t345 - qJ(5) * t373 + t354 * t375 + qJDD(5) + t326;
t410 = -m(6) * t322 + t345 * mrSges(6,1) + t374 * t352;
t401 = -m(5) * t326 + t345 * mrSges(5,1) + t432 * t346 + t374 * t353 + t422 * t375 + t410;
t315 = m(4) * t331 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t370 + qJD(3) * t376 - t365 * t382 + t401;
t426 = t394 * t306 + t397 * t315;
t308 = t396 * t310 + t393 * t313;
t425 = -t429 * t374 - t430 * t375 - t435 * t379;
t412 = t397 * t306 - t394 * t315;
t407 = -mrSges(3,1) * t392 + mrSges(3,2) * t391;
t404 = mrSges(3,3) * qJDD(1) + t400 * t407;
t402 = m(4) * t368 - mrSges(4,1) * t369 + mrSges(4,2) * t370 + t376 * t381 + t377 * t382 + t308;
t380 = -qJDD(1) * pkin(1) - qJ(2) * t400 + t408;
t371 = -t383 * t391 + t413;
t363 = Ifges(4,1) * t382 - Ifges(4,4) * t381 + Ifges(4,5) * qJD(3);
t362 = Ifges(4,4) * t382 - Ifges(4,2) * t381 + Ifges(4,6) * qJD(3);
t361 = Ifges(4,5) * t382 - Ifges(4,6) * t381 + Ifges(4,3) * qJD(3);
t317 = mrSges(6,2) * t346 + t355 * t375 - t410;
t307 = mrSges(5,2) * t326 + mrSges(6,2) * t322 - mrSges(5,3) * t323 - mrSges(6,3) * t319 - qJ(5) * t316 + t431 * t345 + t437 * t346 + t430 * t366 - t425 * t374 + t424 * t379;
t303 = t421 * t400 * mrSges(3,3) + m(3) * t380 + t407 * qJDD(1) + t402;
t302 = -mrSges(5,1) * t326 + mrSges(5,3) * t324 - mrSges(6,1) * t322 + mrSges(6,3) * t321 - pkin(4) * t317 + qJ(5) * t415 + (-qJ(5) * t355 + t423) * t379 + t425 * t375 + (-mrSges(6,2) * qJ(5) + t429) * t366 + t431 * t346 + t436 * t345;
t301 = -mrSges(4,1) * t368 + mrSges(4,3) * t332 + Ifges(4,4) * t370 + Ifges(4,2) * t369 + Ifges(4,6) * qJDD(3) - pkin(3) * t308 + qJD(3) * t363 - t382 * t361 - t434;
t300 = mrSges(4,2) * t368 - mrSges(4,3) * t331 + Ifges(4,1) * t370 + Ifges(4,4) * t369 + Ifges(4,5) * qJDD(3) - pkin(7) * t308 - qJD(3) * t362 - t302 * t393 + t307 * t396 - t361 * t381;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t414 - mrSges(2,2) * t409 + t391 * (mrSges(3,2) * t380 - mrSges(3,3) * t371 + t397 * t300 - t394 * t301 - pkin(6) * t426 + (Ifges(3,1) * t391 + Ifges(3,4) * t392) * qJDD(1)) + t392 * (-mrSges(3,1) * t380 + mrSges(3,3) * t372 + t394 * t300 + t397 * t301 - pkin(2) * t402 + pkin(6) * t412 + (Ifges(3,4) * t391 + Ifges(3,2) * t392) * qJDD(1)) - pkin(1) * t303 + qJ(2) * ((m(3) * t372 + t404 * t392 + t412) * t392 + (-m(3) * t371 + t404 * t391 - t426) * t391); t303; mrSges(4,1) * t331 - mrSges(4,2) * t332 + Ifges(4,5) * t370 + Ifges(4,6) * t369 + Ifges(4,3) * qJDD(3) + pkin(3) * t401 + pkin(7) * t411 + t396 * t302 + t393 * t307 + t382 * t362 + t381 * t363; t434; t317;];
tauJ = t1;
