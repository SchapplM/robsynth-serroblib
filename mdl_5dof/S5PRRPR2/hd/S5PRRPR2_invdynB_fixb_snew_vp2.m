% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:29
% EndTime: 2019-12-05 16:17:31
% DurationCPUTime: 2.32s
% Computational Cost: add. (28286->195), mult. (39158->261), div. (0->0), fcn. (25113->10), ass. (0->97)
t438 = 2 * qJD(4);
t400 = sin(pkin(8));
t402 = cos(pkin(8));
t384 = t400 * g(1) - t402 * g(2);
t386 = -t402 * g(1) - t400 * g(2);
t405 = sin(qJ(2));
t408 = cos(qJ(2));
t370 = t408 * t384 - t405 * t386;
t368 = qJDD(2) * pkin(2) + t370;
t371 = t405 * t384 + t408 * t386;
t409 = qJD(2) ^ 2;
t369 = -t409 * pkin(2) + t371;
t404 = sin(qJ(3));
t407 = cos(qJ(3));
t361 = t404 * t368 + t407 * t369;
t397 = (qJD(2) + qJD(3));
t395 = t397 ^ 2;
t396 = qJDD(2) + qJDD(3);
t359 = -t395 * pkin(3) + t396 * qJ(4) + t361;
t437 = (t397 * t438) + t359;
t399 = sin(pkin(9));
t436 = mrSges(5,2) * t399;
t434 = mrSges(5,3) * t396;
t433 = t399 * t397;
t403 = sin(qJ(5));
t432 = t399 * t403;
t406 = cos(qJ(5));
t431 = t399 * t406;
t401 = cos(pkin(9));
t430 = t401 * t396;
t429 = t401 * t397;
t398 = -g(3) + qJDD(1);
t428 = t401 * t398;
t355 = t399 * t398 + t437 * t401;
t378 = (-mrSges(5,1) * t401 + t436) * t397;
t415 = -pkin(4) * t401 - pkin(7) * t399;
t380 = t415 * t397;
t353 = t380 * t429 + t355;
t360 = t407 * t368 - t404 * t369;
t411 = -t395 * qJ(4) + qJDD(4) - t360;
t356 = (-pkin(3) + t415) * t396 + t411;
t350 = -t403 * t353 + t406 * t356;
t387 = qJD(5) - t429;
t424 = t397 * t432;
t373 = -t387 * mrSges(6,2) - mrSges(6,3) * t424;
t375 = (mrSges(6,1) * t403 + mrSges(6,2) * t406) * t433;
t425 = qJD(5) * t397;
t377 = (t396 * t406 - t403 * t425) * t399;
t385 = qJDD(5) - t430;
t423 = t397 * t431;
t348 = m(6) * t350 + t385 * mrSges(6,1) - t377 * mrSges(6,3) + t387 * t373 - t375 * t423;
t351 = t406 * t353 + t403 * t356;
t374 = t387 * mrSges(6,1) - mrSges(6,3) * t423;
t376 = (-t396 * t403 - t406 * t425) * t399;
t349 = m(6) * t351 - t385 * mrSges(6,2) + t376 * mrSges(6,3) - t387 * t374 - t375 * t424;
t417 = -t403 * t348 + t406 * t349;
t341 = m(5) * t355 + (t378 * t397 + t434) * t401 + t417;
t354 = -t437 * t399 + t428;
t352 = -t428 + (t359 + (t438 + t380) * t397) * t399;
t412 = -m(6) * t352 + t376 * mrSges(6,1) - t377 * mrSges(6,2);
t346 = m(5) * t354 + (-t434 + (-t373 * t403 - t374 * t406 - t378) * t397) * t399 + t412;
t418 = t401 * t341 - t399 * t346;
t335 = m(4) * t361 - t395 * mrSges(4,1) - t396 * mrSges(4,2) + t418;
t342 = t406 * t348 + t403 * t349;
t358 = -t396 * pkin(3) + t411;
t410 = -m(5) * t358 + mrSges(5,1) * t430 - t342 + (t399 ^ 2 + t401 ^ 2) * mrSges(5,3) * t395;
t338 = m(4) * t360 - t395 * mrSges(4,2) + (mrSges(4,1) - t436) * t396 + t410;
t330 = t404 * t335 + t407 * t338;
t328 = m(3) * t370 + qJDD(2) * mrSges(3,1) - t409 * mrSges(3,2) + t330;
t419 = t407 * t335 - t404 * t338;
t329 = m(3) * t371 - t409 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t419;
t323 = t408 * t328 + t405 * t329;
t321 = m(2) * t384 + t323;
t420 = -t405 * t328 + t408 * t329;
t322 = m(2) * t386 + t420;
t427 = t402 * t321 + t400 * t322;
t336 = t399 * t341 + t401 * t346;
t422 = m(4) * t398 + t336;
t421 = -t400 * t321 + t402 * t322;
t416 = m(3) * t398 + t422;
t414 = Ifges(5,1) * t399 + Ifges(5,4) * t401;
t413 = Ifges(5,5) * t399 + Ifges(5,6) * t401;
t379 = t413 * t397;
t364 = Ifges(6,5) * t387 + (Ifges(6,1) * t406 - Ifges(6,4) * t403) * t433;
t363 = Ifges(6,6) * t387 + (Ifges(6,4) * t406 - Ifges(6,2) * t403) * t433;
t362 = Ifges(6,3) * t387 + (Ifges(6,5) * t406 - Ifges(6,6) * t403) * t433;
t344 = mrSges(6,2) * t352 - mrSges(6,3) * t350 + Ifges(6,1) * t377 + Ifges(6,4) * t376 + Ifges(6,5) * t385 - t362 * t424 - t387 * t363;
t343 = -mrSges(6,1) * t352 + mrSges(6,3) * t351 + Ifges(6,4) * t377 + Ifges(6,2) * t376 + Ifges(6,6) * t385 - t362 * t423 + t387 * t364;
t332 = Ifges(5,2) * t430 - mrSges(5,1) * t358 - mrSges(6,1) * t350 + mrSges(6,2) * t351 + mrSges(5,3) * t355 - Ifges(6,5) * t377 - Ifges(6,6) * t376 - Ifges(6,3) * t385 - pkin(4) * t342 + (Ifges(5,4) * t396 + (-t363 * t406 - t364 * t403 - t379) * t397) * t399;
t331 = mrSges(5,2) * t358 - mrSges(5,3) * t354 - pkin(7) * t342 - t403 * t343 + t406 * t344 + t379 * t429 + t414 * t396;
t324 = t395 * Ifges(4,5) - mrSges(4,1) * t398 + mrSges(4,3) * t361 - mrSges(5,1) * t354 + mrSges(5,2) * t355 - t403 * t344 - t406 * t343 - pkin(4) * t412 - pkin(7) * t417 - pkin(3) * t336 + (Ifges(4,6) - t413) * t396 + (-pkin(4) * (-t373 * t432 - t374 * t431) + (-t399 * (Ifges(5,4) * t399 + Ifges(5,2) * t401) + t401 * t414) * t397) * t397;
t317 = mrSges(4,2) * t398 - mrSges(4,3) * t360 + Ifges(4,5) * t396 - t395 * Ifges(4,6) - qJ(4) * t336 + t401 * t331 - t399 * t332;
t316 = mrSges(3,2) * t398 - mrSges(3,3) * t370 + Ifges(3,5) * qJDD(2) - t409 * Ifges(3,6) - pkin(6) * t330 + t407 * t317 - t404 * t324;
t315 = -mrSges(3,1) * t398 + mrSges(3,3) * t371 + t409 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t422 + pkin(6) * t419 + t404 * t317 + t407 * t324;
t314 = mrSges(2,2) * t398 - mrSges(2,3) * t384 - pkin(5) * t323 - t405 * t315 + t408 * t316;
t313 = -mrSges(2,1) * t398 + mrSges(2,3) * t386 - pkin(1) * t416 + pkin(5) * t420 + t408 * t315 + t405 * t316;
t1 = [-m(1) * g(1) + t421; -m(1) * g(2) + t427; -m(1) * g(3) + m(2) * t398 + t416; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t427 - t400 * t313 + t402 * t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t421 + t402 * t313 + t400 * t314; pkin(1) * t323 + mrSges(2,1) * t384 - mrSges(2,2) * t386 + pkin(2) * t330 + mrSges(3,1) * t370 - mrSges(3,2) * t371 + t401 * t332 + pkin(3) * (-t396 * t436 + t410) + qJ(4) * t418 + t399 * t331 + mrSges(4,1) * t360 - mrSges(4,2) * t361 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t396 + Ifges(3,3) * qJDD(2);];
tauB = t1;
