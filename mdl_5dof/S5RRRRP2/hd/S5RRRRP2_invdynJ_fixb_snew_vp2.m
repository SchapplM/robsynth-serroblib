% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP2
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
% m [6x1]
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:56
% EndTime: 2022-01-20 11:48:58
% DurationCPUTime: 1.21s
% Computational Cost: add. (9308->201), mult. (11771->247), div. (0->0), fcn. (6989->8), ass. (0->84)
t399 = Ifges(5,4) + Ifges(6,4);
t406 = Ifges(5,2) + Ifges(6,2);
t403 = Ifges(5,6) + Ifges(6,6);
t405 = Ifges(5,1) + Ifges(6,1);
t404 = Ifges(5,5) + Ifges(6,5);
t402 = Ifges(5,3) + Ifges(6,3);
t372 = qJD(1) + qJD(2);
t374 = sin(qJ(4));
t375 = sin(qJ(3));
t378 = cos(qJ(4));
t379 = cos(qJ(3));
t347 = (-t374 * t375 + t378 * t379) * t372;
t348 = (t374 * t379 + t375 * t378) * t372;
t371 = qJD(3) + qJD(4);
t401 = t406 * t347 + t399 * t348 + t403 * t371;
t368 = t372 ^ 2;
t400 = pkin(3) * t368;
t398 = t372 * t375;
t397 = t372 * t379;
t377 = sin(qJ(1));
t381 = cos(qJ(1));
t390 = t377 * g(1) - t381 * g(2);
t358 = qJDD(1) * pkin(1) + t390;
t387 = -t381 * g(1) - t377 * g(2);
t359 = -qJD(1) ^ 2 * pkin(1) + t387;
t376 = sin(qJ(2));
t380 = cos(qJ(2));
t336 = t376 * t358 + t380 * t359;
t370 = qJDD(1) + qJDD(2);
t333 = -t368 * pkin(2) + t370 * pkin(7) + t336;
t396 = t375 * t333;
t393 = qJD(3) * t372;
t353 = t375 * t370 + t379 * t393;
t302 = qJDD(3) * pkin(3) - t353 * pkin(8) - t396 + (pkin(8) * t393 + t375 * t400 - g(3)) * t379;
t319 = -t375 * g(3) + t379 * t333;
t354 = t379 * t370 - t375 * t393;
t362 = qJD(3) * pkin(3) - pkin(8) * t398;
t373 = t379 ^ 2;
t303 = t354 * pkin(8) - qJD(3) * t362 - t373 * t400 + t319;
t297 = t378 * t302 - t374 * t303;
t317 = t347 * qJD(4) + t378 * t353 + t374 * t354;
t330 = -t347 * mrSges(6,1) + t348 * mrSges(6,2);
t331 = -t347 * mrSges(5,1) + t348 * mrSges(5,2);
t339 = -t371 * mrSges(5,2) + t347 * mrSges(5,3);
t369 = qJDD(3) + qJDD(4);
t291 = -0.2e1 * qJD(5) * t348 + (t347 * t371 - t317) * qJ(5) + (t347 * t348 + t369) * pkin(4) + t297;
t338 = -t371 * mrSges(6,2) + t347 * mrSges(6,3);
t392 = m(6) * t291 + t369 * mrSges(6,1) + t371 * t338;
t282 = m(5) * t297 + t369 * mrSges(5,1) + t371 * t339 + (-t330 - t331) * t348 + (-mrSges(5,3) - mrSges(6,3)) * t317 + t392;
t298 = t374 * t302 + t378 * t303;
t316 = -t348 * qJD(4) - t374 * t353 + t378 * t354;
t341 = t371 * mrSges(6,1) - t348 * mrSges(6,3);
t342 = t371 * mrSges(5,1) - t348 * mrSges(5,3);
t340 = t371 * pkin(4) - t348 * qJ(5);
t343 = t347 ^ 2;
t293 = -t343 * pkin(4) + t316 * qJ(5) + 0.2e1 * qJD(5) * t347 - t371 * t340 + t298;
t391 = m(6) * t293 + t316 * mrSges(6,3) + t347 * t330;
t285 = m(5) * t298 + t316 * mrSges(5,3) + t347 * t331 + (-t341 - t342) * t371 + (-mrSges(5,2) - mrSges(6,2)) * t369 + t391;
t279 = t378 * t282 + t374 * t285;
t395 = -t403 * t347 - t404 * t348 - t402 * t371;
t394 = -t399 * t347 - t405 * t348 - t404 * t371;
t318 = -t379 * g(3) - t396;
t352 = (-mrSges(4,1) * t379 + mrSges(4,2) * t375) * t372;
t360 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t398;
t361 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t397;
t388 = -t374 * t282 + t378 * t285;
t389 = -t375 * (m(4) * t318 + qJDD(3) * mrSges(4,1) - t353 * mrSges(4,3) + qJD(3) * t361 - t352 * t398 + t279) + t379 * (m(4) * t319 - qJDD(3) * mrSges(4,2) + t354 * mrSges(4,3) - qJD(3) * t360 + t352 * t397 + t388);
t335 = t380 * t358 - t376 * t359;
t386 = -t370 * pkin(2) - t335;
t304 = -t354 * pkin(3) + t362 * t398 + (-pkin(8) * t373 - pkin(7)) * t368 + t386;
t295 = -t316 * pkin(4) - t343 * qJ(5) + t348 * t340 + qJDD(5) + t304;
t288 = m(6) * t295 - t316 * mrSges(6,1) + t317 * mrSges(6,2) - t347 * t338 + t348 * t341;
t275 = -mrSges(5,1) * t304 + mrSges(5,3) * t298 - mrSges(6,1) * t295 + mrSges(6,3) * t293 - pkin(4) * t288 + qJ(5) * t391 + (-qJ(5) * t341 - t394) * t371 + (-qJ(5) * mrSges(6,2) + t403) * t369 + t395 * t348 + t399 * t317 + t406 * t316;
t287 = -t317 * mrSges(6,3) - t348 * t330 + t392;
t276 = mrSges(5,2) * t304 + mrSges(6,2) * t295 - mrSges(5,3) * t297 - mrSges(6,3) * t291 - qJ(5) * t287 + t399 * t316 + t405 * t317 - t395 * t347 + t404 * t369 - t401 * t371;
t332 = -t368 * pkin(7) + t386;
t344 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t375 + Ifges(4,6) * t379) * t372;
t345 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t375 + Ifges(4,2) * t379) * t372;
t346 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t375 + Ifges(4,4) * t379) * t372;
t384 = m(5) * t304 - t316 * mrSges(5,1) + t317 * mrSges(5,2) - t347 * t339 + t348 * t342 + t288;
t382 = -m(4) * t332 + t354 * mrSges(4,1) - t353 * mrSges(4,2) - t360 * t398 + t361 * t397 - t384;
t385 = -mrSges(3,2) * t336 + t379 * (-mrSges(4,1) * t332 + mrSges(4,3) * t319 + Ifges(4,4) * t353 + Ifges(4,2) * t354 + Ifges(4,6) * qJDD(3) - pkin(3) * t384 + pkin(8) * t388 + qJD(3) * t346 + t378 * t275 + t374 * t276 - t344 * t398) + t375 * (mrSges(4,2) * t332 - mrSges(4,3) * t318 + Ifges(4,1) * t353 + Ifges(4,4) * t354 + Ifges(4,5) * qJDD(3) - pkin(8) * t279 - qJD(3) * t345 - t374 * t275 + t378 * t276 + t344 * t397) + pkin(7) * t389 + pkin(2) * t382 + mrSges(3,1) * t335 + Ifges(3,3) * t370;
t383 = mrSges(5,1) * t297 + mrSges(6,1) * t291 - mrSges(5,2) * t298 - mrSges(6,2) * t293 + pkin(4) * t287 + t403 * t316 + t404 * t317 + t394 * t347 + t401 * t348 + t402 * t369;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t390 - mrSges(2,2) * t387 + pkin(1) * (t376 * (m(3) * t336 - t368 * mrSges(3,1) - t370 * mrSges(3,2) + t389) + t380 * (m(3) * t335 + t370 * mrSges(3,1) - t368 * mrSges(3,2) + t382)) + t385; t385; t383 + Ifges(4,5) * t353 + Ifges(4,6) * t354 + mrSges(4,1) * t318 - mrSges(4,2) * t319 + pkin(3) * t279 + (t375 * t345 - t379 * t346) * t372 + Ifges(4,3) * qJDD(3); t383; t288;];
tauJ = t1;
