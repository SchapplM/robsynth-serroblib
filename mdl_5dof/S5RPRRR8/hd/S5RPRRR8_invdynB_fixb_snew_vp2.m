% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:52
% EndTime: 2019-12-31 19:05:53
% DurationCPUTime: 1.55s
% Computational Cost: add. (21330->226), mult. (27593->276), div. (0->0), fcn. (13015->8), ass. (0->93)
t435 = -m(3) - m(4);
t434 = -pkin(1) - pkin(2);
t433 = -mrSges(2,1) - mrSges(3,1);
t432 = Ifges(3,4) + Ifges(2,5);
t431 = Ifges(2,6) - Ifges(3,6);
t400 = -qJD(1) + qJD(3);
t407 = sin(qJ(4));
t430 = t400 * t407;
t411 = cos(qJ(4));
t429 = t400 * t411;
t409 = sin(qJ(1));
t413 = cos(qJ(1));
t391 = -t413 * g(1) - t409 * g(2);
t414 = qJD(1) ^ 2;
t419 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t391;
t380 = -t414 * pkin(1) + t419;
t372 = t434 * t414 + t419;
t390 = t409 * g(1) - t413 * g(2);
t418 = -t414 * qJ(2) + qJDD(2) - t390;
t374 = t434 * qJDD(1) + t418;
t408 = sin(qJ(3));
t412 = cos(qJ(3));
t361 = t412 * t372 + t408 * t374;
t396 = t400 ^ 2;
t398 = -qJDD(1) + qJDD(3);
t359 = -t396 * pkin(3) + t398 * pkin(7) + t361;
t352 = t411 * g(3) - t407 * t359;
t427 = qJD(4) * t400;
t426 = t411 * t427;
t385 = t407 * t398 + t426;
t349 = (-t385 + t426) * pkin(8) + (t396 * t407 * t411 + qJDD(4)) * pkin(4) + t352;
t353 = t407 * g(3) + t411 * t359;
t386 = t411 * t398 - t407 * t427;
t389 = qJD(4) * pkin(4) - pkin(8) * t430;
t404 = t411 ^ 2;
t350 = -t404 * t396 * pkin(4) + t386 * pkin(8) - qJD(4) * t389 + t353;
t406 = sin(qJ(5));
t410 = cos(qJ(5));
t347 = t410 * t349 - t406 * t350;
t378 = (-t406 * t407 + t410 * t411) * t400;
t357 = t378 * qJD(5) + t410 * t385 + t406 * t386;
t379 = (t406 * t411 + t407 * t410) * t400;
t366 = -t378 * mrSges(6,1) + t379 * mrSges(6,2);
t399 = qJD(4) + qJD(5);
t367 = -t399 * mrSges(6,2) + t378 * mrSges(6,3);
t397 = qJDD(4) + qJDD(5);
t345 = m(6) * t347 + t397 * mrSges(6,1) - t357 * mrSges(6,3) - t379 * t366 + t399 * t367;
t348 = t406 * t349 + t410 * t350;
t356 = -t379 * qJD(5) - t406 * t385 + t410 * t386;
t368 = t399 * mrSges(6,1) - t379 * mrSges(6,3);
t346 = m(6) * t348 - t397 * mrSges(6,2) + t356 * mrSges(6,3) + t378 * t366 - t399 * t368;
t338 = t410 * t345 + t406 * t346;
t384 = (-mrSges(5,1) * t411 + mrSges(5,2) * t407) * t400;
t388 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t429;
t336 = m(5) * t352 + qJDD(4) * mrSges(5,1) - t385 * mrSges(5,3) + qJD(4) * t388 - t384 * t430 + t338;
t387 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t430;
t422 = -t406 * t345 + t410 * t346;
t337 = m(5) * t353 - qJDD(4) * mrSges(5,2) + t386 * mrSges(5,3) - qJD(4) * t387 + t384 * t429 + t422;
t423 = -t407 * t336 + t411 * t337;
t332 = m(4) * t361 - t396 * mrSges(4,1) - t398 * mrSges(4,2) + t423;
t360 = -t408 * t372 + t412 * t374;
t420 = -t398 * pkin(3) - t360;
t358 = -t396 * pkin(7) + t420;
t351 = t389 * t430 - t386 * pkin(4) + (-pkin(8) * t404 - pkin(7)) * t396 + t420;
t417 = m(6) * t351 - t356 * mrSges(6,1) + t357 * mrSges(6,2) - t378 * t367 + t379 * t368;
t415 = -m(5) * t358 + t386 * mrSges(5,1) - t385 * mrSges(5,2) - t387 * t430 + t388 * t429 - t417;
t341 = m(4) * t360 + t398 * mrSges(4,1) - t396 * mrSges(4,2) + t415;
t424 = t412 * t332 - t408 * t341;
t421 = m(3) * t380 + qJDD(1) * mrSges(3,3) + t424;
t327 = m(2) * t391 - qJDD(1) * mrSges(2,2) + t433 * t414 + t421;
t329 = t408 * t332 + t412 * t341;
t383 = -qJDD(1) * pkin(1) + t418;
t416 = -m(3) * t383 + qJDD(1) * mrSges(3,1) + t414 * mrSges(3,3) - t329;
t328 = m(2) * t390 + qJDD(1) * mrSges(2,1) - t414 * mrSges(2,2) + t416;
t428 = t409 * t327 + t413 * t328;
t425 = t413 * t327 - t409 * t328;
t334 = t411 * t336 + t407 * t337;
t377 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t407 + Ifges(5,4) * t411) * t400;
t376 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t407 + Ifges(5,2) * t411) * t400;
t375 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t407 + Ifges(5,6) * t411) * t400;
t364 = Ifges(6,1) * t379 + Ifges(6,4) * t378 + Ifges(6,5) * t399;
t363 = Ifges(6,4) * t379 + Ifges(6,2) * t378 + Ifges(6,6) * t399;
t362 = Ifges(6,5) * t379 + Ifges(6,6) * t378 + Ifges(6,3) * t399;
t340 = mrSges(6,2) * t351 - mrSges(6,3) * t347 + Ifges(6,1) * t357 + Ifges(6,4) * t356 + Ifges(6,5) * t397 + t378 * t362 - t399 * t363;
t339 = -mrSges(6,1) * t351 + mrSges(6,3) * t348 + Ifges(6,4) * t357 + Ifges(6,2) * t356 + Ifges(6,6) * t397 - t379 * t362 + t399 * t364;
t333 = t435 * g(3) - t334;
t330 = mrSges(5,2) * t358 - mrSges(5,3) * t352 + Ifges(5,1) * t385 + Ifges(5,4) * t386 + Ifges(5,5) * qJDD(4) - pkin(8) * t338 - qJD(4) * t376 - t406 * t339 + t410 * t340 + t375 * t429;
t323 = -mrSges(5,1) * t358 + mrSges(5,3) * t353 + Ifges(5,4) * t385 + Ifges(5,2) * t386 + Ifges(5,6) * qJDD(4) - pkin(4) * t417 + pkin(8) * t422 + qJD(4) * t377 + t410 * t339 + t406 * t340 - t375 * t430;
t322 = Ifges(4,6) * t398 + t396 * Ifges(4,5) - mrSges(4,1) * g(3) + mrSges(4,3) * t361 - Ifges(5,5) * t385 - Ifges(5,6) * t386 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t352 + mrSges(5,2) * t353 - Ifges(6,5) * t357 - Ifges(6,6) * t356 - Ifges(6,3) * t397 - t379 * t363 + t378 * t364 - mrSges(6,1) * t347 + mrSges(6,2) * t348 - pkin(4) * t338 - pkin(3) * t334 + (-t407 * t376 + t411 * t377) * t400;
t321 = mrSges(4,2) * g(3) - mrSges(4,3) * t360 + Ifges(4,5) * t398 - t396 * Ifges(4,6) - pkin(7) * t334 - t407 * t323 + t411 * t330;
t320 = mrSges(3,2) * t383 - mrSges(2,3) * t390 - pkin(6) * t329 - qJ(2) * t333 + t412 * t321 - t408 * t322 - t431 * t414 + t432 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t319 = mrSges(2,3) * t391 + mrSges(3,2) * t380 - t408 * t321 - t412 * t322 + pkin(2) * t334 - pkin(6) * t424 - pkin(1) * t333 + t432 * t414 + t431 * qJDD(1) + (pkin(2) * m(4) - t433) * g(3);
t1 = [-m(1) * g(1) + t425; -m(1) * g(2) + t428; (-m(1) - m(2) + t435) * g(3) - t334; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t428 - t409 * t319 + t413 * t320; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t425 + t413 * t319 + t409 * t320; qJ(2) * (-t414 * mrSges(3,1) + t421) + mrSges(2,1) * t390 - mrSges(2,2) * t391 + pkin(1) * t416 - pkin(2) * t329 - mrSges(3,1) * t383 + mrSges(3,3) * t380 - pkin(3) * t415 - pkin(7) * t423 - t407 * t330 - t411 * t323 - mrSges(4,1) * t360 + mrSges(4,2) * t361 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(4,3) * t398 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
