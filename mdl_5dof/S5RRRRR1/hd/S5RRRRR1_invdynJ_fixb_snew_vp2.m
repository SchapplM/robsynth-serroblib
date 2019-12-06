% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:50:10
% EndTime: 2019-12-05 18:50:13
% DurationCPUTime: 1.97s
% Computational Cost: add. (16240->244), mult. (35622->316), div. (0->0), fcn. (27118->10), ass. (0->98)
t396 = qJD(1) ^ 2;
t390 = sin(qJ(1));
t395 = cos(qJ(1));
t403 = -t395 * g(1) - t390 * g(2);
t372 = -t396 * pkin(1) + t403;
t389 = sin(qJ(2));
t394 = cos(qJ(2));
t363 = t394 * g(3) - t389 * t372;
t358 = (t389 * t394 * t396 + qJDD(2)) * pkin(2) + t363;
t364 = t389 * g(3) + t394 * t372;
t359 = (-t394 ^ 2 * t396 - qJD(2) ^ 2) * pkin(2) + t364;
t388 = sin(qJ(3));
t393 = cos(qJ(3));
t333 = t393 * t358 - t388 * t359;
t368 = (t389 * t388 - t394 * t393) * qJD(1);
t369 = (-t394 * t388 - t389 * t393) * qJD(1);
t383 = qJDD(2) + qJDD(3);
t318 = (t368 * t369 + t383) * pkin(3) + t333;
t334 = t388 * t358 + t393 * t359;
t384 = qJD(2) + qJD(3);
t325 = (-t368 ^ 2 - t384 ^ 2) * pkin(3) + t334;
t387 = sin(qJ(4));
t392 = cos(qJ(4));
t305 = t387 * t318 + t392 * t325;
t407 = qJD(1) * qJD(2);
t373 = -t389 * qJDD(1) - t394 * t407;
t406 = t389 * t407;
t374 = -t394 * qJDD(1) + t406;
t341 = -t369 * qJD(3) - t388 * t373 + t393 * t374;
t342 = t368 * qJD(3) + t393 * t373 + t388 * t374;
t351 = t387 * t368 + t392 * t369;
t314 = -t351 * qJD(4) + t392 * t341 - t387 * t342;
t350 = t392 * t368 - t387 * t369;
t330 = -t350 * mrSges(5,1) + t351 * mrSges(5,2);
t379 = qJD(4) + t384;
t344 = t379 * mrSges(5,1) - t351 * mrSges(5,3);
t378 = qJDD(4) + t383;
t315 = t350 * qJD(4) + t387 * t341 + t392 * t342;
t405 = t390 * g(1) - t395 * g(2);
t371 = qJDD(1) * pkin(1) + t405;
t356 = (-t374 - t406) * pkin(2) + t371;
t321 = t356 + (t369 * t384 - t341) * pkin(3);
t297 = (-t350 * t379 - t315) * pkin(6) + (t351 * t379 - t314) * pkin(4) + t321;
t331 = -t350 * pkin(4) - t351 * pkin(6);
t377 = t379 ^ 2;
t299 = -t377 * pkin(4) + t378 * pkin(6) + t350 * t331 + t305;
t386 = sin(qJ(5));
t391 = cos(qJ(5));
t295 = t391 * t297 - t386 * t299;
t336 = -t386 * t351 + t391 * t379;
t302 = t336 * qJD(5) + t391 * t315 + t386 * t378;
t312 = qJDD(5) - t314;
t337 = t391 * t351 + t386 * t379;
t319 = -t336 * mrSges(6,1) + t337 * mrSges(6,2);
t348 = qJD(5) - t350;
t322 = -t348 * mrSges(6,2) + t336 * mrSges(6,3);
t292 = m(6) * t295 + t312 * mrSges(6,1) - t302 * mrSges(6,3) - t337 * t319 + t348 * t322;
t296 = t386 * t297 + t391 * t299;
t301 = -t337 * qJD(5) - t386 * t315 + t391 * t378;
t323 = t348 * mrSges(6,1) - t337 * mrSges(6,3);
t293 = m(6) * t296 - t312 * mrSges(6,2) + t301 * mrSges(6,3) + t336 * t319 - t348 * t323;
t404 = -t386 * t292 + t391 * t293;
t280 = m(5) * t305 - t378 * mrSges(5,2) + t314 * mrSges(5,3) + t350 * t330 - t379 * t344 + t404;
t304 = t392 * t318 - t387 * t325;
t343 = -t379 * mrSges(5,2) + t350 * mrSges(5,3);
t298 = -t378 * pkin(4) - t377 * pkin(6) + t351 * t331 - t304;
t402 = -m(6) * t298 + t301 * mrSges(6,1) - t302 * mrSges(6,2) + t336 * t322 - t337 * t323;
t288 = m(5) * t304 + t378 * mrSges(5,1) - t315 * mrSges(5,3) - t351 * t330 + t379 * t343 + t402;
t410 = t387 * t280 + t392 * t288;
t282 = t391 * t292 + t386 * t293;
t409 = t389 * qJD(1);
t408 = t394 * qJD(1);
t401 = m(5) * t321 - t314 * mrSges(5,1) + t315 * mrSges(5,2) - t350 * t343 + t351 * t344 + t282;
t306 = Ifges(6,5) * t337 + Ifges(6,6) * t336 + Ifges(6,3) * t348;
t308 = Ifges(6,1) * t337 + Ifges(6,4) * t336 + Ifges(6,5) * t348;
t285 = -mrSges(6,1) * t298 + mrSges(6,3) * t296 + Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t312 - t337 * t306 + t348 * t308;
t307 = Ifges(6,4) * t337 + Ifges(6,2) * t336 + Ifges(6,6) * t348;
t286 = mrSges(6,2) * t298 - mrSges(6,3) * t295 + Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t312 + t336 * t306 - t348 * t307;
t327 = Ifges(5,4) * t351 + Ifges(5,2) * t350 + Ifges(5,6) * t379;
t328 = Ifges(5,1) * t351 + Ifges(5,4) * t350 + Ifges(5,5) * t379;
t400 = mrSges(5,1) * t304 - mrSges(5,2) * t305 + Ifges(5,5) * t315 + Ifges(5,6) * t314 + Ifges(5,3) * t378 + pkin(4) * t402 + pkin(6) * t404 + t391 * t285 + t386 * t286 + t351 * t327 - t350 * t328;
t399 = mrSges(6,1) * t295 - mrSges(6,2) * t296 + Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t312 + t337 * t307 - t336 * t308;
t360 = -t384 * mrSges(4,2) + t368 * mrSges(4,3);
t361 = t384 * mrSges(4,1) - t369 * mrSges(4,3);
t398 = m(4) * t356 - t341 * mrSges(4,1) + t342 * mrSges(4,2) - t368 * t360 + t369 * t361 + t401;
t346 = Ifges(4,4) * t369 + Ifges(4,2) * t368 + Ifges(4,6) * t384;
t347 = Ifges(4,1) * t369 + Ifges(4,4) * t368 + Ifges(4,5) * t384;
t397 = mrSges(4,1) * t333 - mrSges(4,2) * t334 + Ifges(4,5) * t342 + Ifges(4,6) * t341 + Ifges(4,3) * t383 + pkin(3) * t410 + t369 * t346 - t368 * t347 + t400;
t367 = Ifges(3,5) * qJD(2) + (-t389 * Ifges(3,1) - t394 * Ifges(3,4)) * qJD(1);
t366 = Ifges(3,6) * qJD(2) + (-t389 * Ifges(3,4) - t394 * Ifges(3,2)) * qJD(1);
t352 = -t368 * mrSges(4,1) + t369 * mrSges(4,2);
t345 = Ifges(4,5) * t369 + Ifges(4,6) * t368 + Ifges(4,3) * t384;
t326 = Ifges(5,5) * t351 + Ifges(5,6) * t350 + Ifges(5,3) * t379;
t277 = -mrSges(5,1) * t321 + mrSges(5,3) * t305 + Ifges(5,4) * t315 + Ifges(5,2) * t314 + Ifges(5,6) * t378 - pkin(4) * t282 - t351 * t326 + t379 * t328 - t399;
t276 = mrSges(5,2) * t321 - mrSges(5,3) * t304 + Ifges(5,1) * t315 + Ifges(5,4) * t314 + Ifges(5,5) * t378 - pkin(6) * t282 - t386 * t285 + t391 * t286 + t350 * t326 - t379 * t327;
t275 = mrSges(4,2) * t356 - mrSges(4,3) * t333 + Ifges(4,1) * t342 + Ifges(4,4) * t341 + Ifges(4,5) * t383 + t392 * t276 - t387 * t277 + t368 * t345 - t384 * t346;
t274 = -mrSges(4,1) * t356 + mrSges(4,3) * t334 + Ifges(4,4) * t342 + Ifges(4,2) * t341 + Ifges(4,6) * t383 - pkin(3) * t401 + t387 * t276 + t392 * t277 - t369 * t345 + t384 * t347;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t405 - mrSges(2,2) * t403 - t389 * (mrSges(3,2) * t371 - mrSges(3,3) * t363 + Ifges(3,1) * t373 + Ifges(3,4) * t374 + Ifges(3,5) * qJDD(2) - qJD(2) * t366 - t388 * t274 + t393 * t275) - t394 * (-mrSges(3,1) * t371 + mrSges(3,3) * t364 + Ifges(3,4) * t373 + Ifges(3,2) * t374 + Ifges(3,6) * qJDD(2) - pkin(2) * t398 + qJD(2) * t367 + t393 * t274 + t388 * t275) + pkin(1) * (t398 + m(3) * t371 + t373 * mrSges(3,2) - t374 * mrSges(3,1) - (qJD(2) * mrSges(3,1) + mrSges(3,3) * t409) * t409 + (-qJD(2) * mrSges(3,2) - mrSges(3,3) * t408) * t408); pkin(2) * (t388 * (m(4) * t334 - t383 * mrSges(4,2) + t341 * mrSges(4,3) + t392 * t280 - t387 * t288 + t368 * t352 - t384 * t361) + t393 * (m(4) * t333 + t383 * mrSges(4,1) - t342 * mrSges(4,3) - t369 * t352 + t384 * t360 + t410)) + Ifges(3,5) * t373 + Ifges(3,6) * t374 + mrSges(3,1) * t363 - mrSges(3,2) * t364 + Ifges(3,3) * qJDD(2) + (-t389 * t366 + t394 * t367) * qJD(1) + t397; t397; t400; t399;];
tauJ = t1;
