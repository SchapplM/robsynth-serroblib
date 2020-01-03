% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:40
% EndTime: 2019-12-31 18:53:42
% DurationCPUTime: 1.39s
% Computational Cost: add. (8787->221), mult. (20571->272), div. (0->0), fcn. (14353->8), ass. (0->98)
t435 = Ifges(5,1) + Ifges(6,1);
t428 = Ifges(5,4) - Ifges(6,5);
t427 = -Ifges(5,5) - Ifges(6,4);
t434 = Ifges(5,2) + Ifges(6,3);
t426 = Ifges(5,6) - Ifges(6,6);
t433 = -Ifges(5,3) - Ifges(6,2);
t398 = qJD(1) ^ 2;
t390 = sin(pkin(8));
t391 = cos(pkin(8));
t393 = sin(qJ(3));
t395 = cos(qJ(3));
t403 = t390 * t393 - t391 * t395;
t378 = t403 * qJD(1);
t404 = t390 * t395 + t391 * t393;
t379 = t404 * qJD(1);
t415 = t379 * qJD(3);
t366 = -qJDD(1) * t403 - t415;
t416 = t378 * qJD(3);
t367 = qJDD(1) * t404 - t416;
t392 = sin(qJ(4));
t431 = cos(qJ(4));
t370 = -qJD(3) * t431 + t379 * t392;
t341 = -qJD(4) * t370 + qJDD(3) * t392 + t367 * t431;
t371 = qJD(3) * t392 + t379 * t431;
t345 = mrSges(6,1) * t370 - mrSges(6,3) * t371;
t394 = sin(qJ(1));
t396 = cos(qJ(1));
t407 = -g(1) * t396 - g(2) * t394;
t380 = -pkin(1) * t398 + qJDD(1) * qJ(2) + t407;
t414 = qJD(1) * qJD(2);
t411 = -g(3) * t391 - 0.2e1 * t390 * t414;
t424 = pkin(6) * qJDD(1);
t430 = pkin(2) * t398;
t355 = (t391 * t430 - t380 - t424) * t390 + t411;
t369 = -g(3) * t390 + (t380 + 0.2e1 * t414) * t391;
t389 = t391 ^ 2;
t356 = -t389 * t430 + t391 * t424 + t369;
t329 = t393 * t355 + t395 * t356;
t364 = pkin(3) * t378 - pkin(7) * t379;
t397 = qJD(3) ^ 2;
t325 = -pkin(3) * t397 + qJDD(3) * pkin(7) - t364 * t378 + t329;
t412 = g(1) * t394 - t396 * g(2);
t406 = qJDD(2) - t412;
t418 = -t390 ^ 2 - t389;
t365 = (-pkin(2) * t391 - pkin(1)) * qJDD(1) + (pkin(6) * t418 - qJ(2)) * t398 + t406;
t327 = (-t367 + t416) * pkin(7) + (-t366 + t415) * pkin(3) + t365;
t321 = -t325 * t392 + t327 * t431;
t344 = pkin(4) * t370 - qJ(5) * t371;
t363 = qJDD(4) - t366;
t376 = qJD(4) + t378;
t375 = t376 ^ 2;
t319 = -pkin(4) * t363 - qJ(5) * t375 + t344 * t371 + qJDD(5) - t321;
t349 = -mrSges(6,2) * t370 + mrSges(6,3) * t376;
t408 = -m(6) * t319 + t363 * mrSges(6,1) + t376 * t349;
t316 = mrSges(6,2) * t341 + t345 * t371 - t408;
t322 = t325 * t431 + t327 * t392;
t318 = -pkin(4) * t375 + qJ(5) * t363 + 0.2e1 * qJD(5) * t376 - t344 * t370 + t322;
t340 = qJD(4) * t371 - qJDD(3) * t431 + t367 * t392;
t352 = -mrSges(6,1) * t376 + mrSges(6,2) * t371;
t413 = m(6) * t318 + t363 * mrSges(6,3) + t376 * t352;
t420 = t370 * t428 - t371 * t435 + t376 * t427;
t421 = t370 * t434 - t371 * t428 - t376 * t426;
t432 = -t426 * t340 - t427 * t341 - t433 * t363 - t420 * t370 - t421 * t371 + mrSges(5,1) * t321 - mrSges(6,1) * t319 - mrSges(5,2) * t322 + mrSges(6,3) * t318 - pkin(4) * t316 + qJ(5) * (-mrSges(6,2) * t340 - t345 * t370 + t413);
t429 = -mrSges(5,3) - mrSges(6,2);
t362 = mrSges(4,1) * t378 + mrSges(4,2) * t379;
t373 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t379;
t351 = mrSges(5,1) * t376 - mrSges(5,3) * t371;
t419 = -mrSges(5,1) * t370 - mrSges(5,2) * t371 - t345;
t311 = m(5) * t322 - mrSges(5,2) * t363 + t340 * t429 - t351 * t376 + t370 * t419 + t413;
t350 = -mrSges(5,2) * t376 - mrSges(5,3) * t370;
t313 = m(5) * t321 + mrSges(5,1) * t363 + t341 * t429 + t350 * t376 + t371 * t419 + t408;
t409 = t311 * t431 - t313 * t392;
t305 = m(4) * t329 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t366 - qJD(3) * t373 - t362 * t378 + t409;
t328 = t355 * t395 - t393 * t356;
t372 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t378;
t324 = -qJDD(3) * pkin(3) - pkin(7) * t397 + t379 * t364 - t328;
t320 = -0.2e1 * qJD(5) * t371 + (t370 * t376 - t341) * qJ(5) + (t371 * t376 + t340) * pkin(4) + t324;
t315 = m(6) * t320 + mrSges(6,1) * t340 - t341 * mrSges(6,3) + t349 * t370 - t371 * t352;
t399 = -m(5) * t324 - t340 * mrSges(5,1) - mrSges(5,2) * t341 - t370 * t350 - t351 * t371 - t315;
t308 = m(4) * t328 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t367 + qJD(3) * t372 - t362 * t379 + t399;
t423 = t305 * t393 + t308 * t395;
t306 = t311 * t392 + t313 * t431;
t422 = t370 * t426 + t371 * t427 + t376 * t433;
t410 = t305 * t395 - t393 * t308;
t405 = -mrSges(3,1) * t391 + mrSges(3,2) * t390;
t402 = mrSges(3,3) * qJDD(1) + t398 * t405;
t401 = m(4) * t365 - mrSges(4,1) * t366 + mrSges(4,2) * t367 + t372 * t378 + t373 * t379 + t306;
t377 = -qJDD(1) * pkin(1) - qJ(2) * t398 + t406;
t368 = -t380 * t390 + t411;
t359 = Ifges(4,1) * t379 - Ifges(4,4) * t378 + Ifges(4,5) * qJD(3);
t358 = Ifges(4,4) * t379 - Ifges(4,2) * t378 + Ifges(4,6) * qJD(3);
t357 = Ifges(4,5) * t379 - Ifges(4,6) * t378 + Ifges(4,3) * qJD(3);
t302 = mrSges(5,2) * t324 + mrSges(6,2) * t319 - mrSges(5,3) * t321 - mrSges(6,3) * t320 - qJ(5) * t315 - t428 * t340 + t341 * t435 - t427 * t363 + t422 * t370 + t421 * t376;
t301 = mrSges(3,3) * t398 * t418 + m(3) * t377 + qJDD(1) * t405 + t401;
t300 = -mrSges(5,1) * t324 - mrSges(6,1) * t320 + mrSges(6,2) * t318 + mrSges(5,3) * t322 - pkin(4) * t315 - t340 * t434 + t428 * t341 + t426 * t363 + t422 * t371 - t420 * t376;
t299 = -mrSges(4,1) * t365 + mrSges(4,3) * t329 + Ifges(4,4) * t367 + Ifges(4,2) * t366 + Ifges(4,6) * qJDD(3) - pkin(3) * t306 + qJD(3) * t359 - t379 * t357 - t432;
t298 = mrSges(4,2) * t365 - mrSges(4,3) * t328 + Ifges(4,1) * t367 + Ifges(4,4) * t366 + Ifges(4,5) * qJDD(3) - pkin(7) * t306 - qJD(3) * t358 - t300 * t392 + t302 * t431 - t357 * t378;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t412 - mrSges(2,2) * t407 + t390 * (mrSges(3,2) * t377 - mrSges(3,3) * t368 + t395 * t298 - t393 * t299 - pkin(6) * t423 + (Ifges(3,1) * t390 + Ifges(3,4) * t391) * qJDD(1)) + t391 * (-mrSges(3,1) * t377 + mrSges(3,3) * t369 + t393 * t298 + t395 * t299 - pkin(2) * t401 + pkin(6) * t410 + (Ifges(3,4) * t390 + Ifges(3,2) * t391) * qJDD(1)) - pkin(1) * t301 + qJ(2) * ((m(3) * t369 + t391 * t402 + t410) * t391 + (-m(3) * t368 + t390 * t402 - t423) * t390); t301; mrSges(4,1) * t328 - mrSges(4,2) * t329 + Ifges(4,5) * t367 + Ifges(4,6) * t366 + Ifges(4,3) * qJDD(3) + pkin(3) * t399 + pkin(7) * t409 + t300 * t431 + t392 * t302 + t379 * t358 + t378 * t359; t432; t316;];
tauJ = t1;
