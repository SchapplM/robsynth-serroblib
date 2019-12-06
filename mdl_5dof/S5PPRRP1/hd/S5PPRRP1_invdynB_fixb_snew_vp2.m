% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:45
% EndTime: 2019-12-05 15:06:47
% DurationCPUTime: 1.25s
% Computational Cost: add. (9388->187), mult. (16474->229), div. (0->0), fcn. (9426->8), ass. (0->84)
t421 = Ifges(5,1) + Ifges(6,1);
t414 = Ifges(5,4) + Ifges(6,4);
t413 = Ifges(5,5) + Ifges(6,5);
t420 = Ifges(5,2) + Ifges(6,2);
t419 = Ifges(5,6) + Ifges(6,6);
t418 = Ifges(5,3) + Ifges(6,3);
t391 = qJD(3) ^ 2;
t384 = sin(pkin(7));
t386 = cos(pkin(7));
t374 = -t386 * g(1) - t384 * g(2);
t382 = -g(3) + qJDD(1);
t383 = sin(pkin(8));
t385 = cos(pkin(8));
t351 = -t383 * t374 + t385 * t382;
t352 = t385 * t374 + t383 * t382;
t388 = sin(qJ(3));
t390 = cos(qJ(3));
t346 = t390 * t351 - t388 * t352;
t393 = -qJDD(3) * pkin(3) - t346;
t344 = -t391 * pkin(6) + t393;
t387 = sin(qJ(4));
t389 = cos(qJ(4));
t404 = qJD(3) * qJD(4);
t400 = t389 * t404;
t369 = t387 * qJDD(3) + t400;
t370 = t389 * qJDD(3) - t387 * t404;
t405 = qJD(3) * t389;
t379 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t405;
t406 = qJD(3) * t387;
t375 = qJD(4) * pkin(4) - qJ(5) * t406;
t381 = t389 ^ 2;
t340 = t375 * t406 - t370 * pkin(4) + qJDD(5) + (-qJ(5) * t381 - pkin(6)) * t391 + t393;
t378 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t405;
t395 = m(6) * t340 - t370 * mrSges(6,1) - t378 * t405;
t415 = -mrSges(5,2) - mrSges(6,2);
t417 = -m(5) * t344 + t370 * mrSges(5,1) + t415 * t369 + t379 * t405 - t395;
t416 = pkin(4) * t391;
t347 = t388 * t351 + t390 * t352;
t345 = -t391 * pkin(3) + qJDD(3) * pkin(6) + t347;
t373 = t384 * g(1) - t386 * g(2);
t372 = qJDD(2) - t373;
t342 = t389 * t345 + t387 * t372;
t368 = (-mrSges(5,1) * t389 + mrSges(5,2) * t387) * qJD(3);
t403 = qJD(3) * qJD(5);
t339 = t370 * qJ(5) - qJD(4) * t375 - t381 * t416 + 0.2e1 * t389 * t403 + t342;
t367 = (-mrSges(6,1) * t389 + mrSges(6,2) * t387) * qJD(3);
t401 = m(6) * t339 + t370 * mrSges(6,3) + t367 * t405;
t376 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t406;
t407 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t406 - t376;
t334 = m(5) * t342 + t370 * mrSges(5,3) + t407 * qJD(4) + t415 * qJDD(4) + t368 * t405 + t401;
t332 = t389 * t334;
t361 = t389 * t372;
t341 = -t387 * t345 + t361;
t338 = qJDD(4) * pkin(4) + t361 + (-t369 + t400) * qJ(5) + (t389 * t416 - t345 - 0.2e1 * t403) * t387;
t402 = m(6) * t338 + qJDD(4) * mrSges(6,1) + qJD(4) * t378;
t333 = m(5) * t341 + qJDD(4) * mrSges(5,1) + qJD(4) * t379 + (-mrSges(5,3) - mrSges(6,3)) * t369 + (-t367 - t368) * t406 + t402;
t324 = m(4) * t347 - t391 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t387 * t333 + t332;
t396 = qJD(3) * t407;
t329 = m(4) * t346 + qJDD(3) * mrSges(4,1) - t391 * mrSges(4,2) + t387 * t396 + t417;
t319 = t388 * t324 + t390 * t329;
t317 = m(3) * t351 + t319;
t397 = t390 * t324 - t388 * t329;
t318 = m(3) * t352 + t397;
t398 = -t383 * t317 + t385 * t318;
t311 = m(2) * t374 + t398;
t327 = t389 * t333 + t387 * t334;
t394 = (-m(3) - m(4)) * t372 - t327;
t326 = m(2) * t373 + t394;
t411 = t384 * t311 + t386 * t326;
t312 = t385 * t317 + t383 * t318;
t410 = t418 * qJD(4) + (t413 * t387 + t419 * t389) * qJD(3);
t409 = -t419 * qJD(4) + (-t414 * t387 - t420 * t389) * qJD(3);
t408 = t413 * qJD(4) + (t421 * t387 + t414 * t389) * qJD(3);
t399 = t386 * t311 - t384 * t326;
t335 = -t369 * mrSges(6,3) - t367 * t406 + t402;
t321 = mrSges(5,2) * t344 + mrSges(6,2) * t340 - mrSges(5,3) * t341 - mrSges(6,3) * t338 - qJ(5) * t335 + t409 * qJD(4) + t413 * qJDD(4) + t421 * t369 + t414 * t370 + t410 * t405;
t320 = -mrSges(5,1) * t344 + mrSges(5,3) * t342 - mrSges(6,1) * t340 + mrSges(6,3) * t339 - pkin(4) * t395 + qJ(5) * t401 + t420 * t370 + (-pkin(4) * mrSges(6,2) + t414) * t369 + (-qJ(5) * mrSges(6,2) + t419) * qJDD(4) + (-qJ(5) * t376 + t408) * qJD(4) + (-pkin(4) * t376 - t410) * t406;
t313 = -mrSges(4,1) * t372 - mrSges(5,1) * t341 - mrSges(6,1) * t338 + mrSges(5,2) * t342 + mrSges(6,2) * t339 + mrSges(4,3) * t347 + t391 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t327 - pkin(4) * t335 - t419 * t370 - t413 * t369 - t418 * qJDD(4) + (t409 * t387 + t408 * t389) * qJD(3);
t308 = mrSges(4,2) * t372 - mrSges(4,3) * t346 + Ifges(4,5) * qJDD(3) - t391 * Ifges(4,6) - pkin(6) * t327 - t387 * t320 + t389 * t321;
t307 = mrSges(3,2) * t372 - mrSges(3,3) * t351 - pkin(5) * t319 + t390 * t308 - t388 * t313;
t306 = -mrSges(3,1) * t372 + mrSges(3,3) * t352 + t388 * t308 + t390 * t313 - pkin(2) * (m(4) * t372 + t327) + pkin(5) * t397;
t305 = -mrSges(2,1) * t382 - pkin(1) * t312 + mrSges(2,3) * t374 - pkin(2) * t319 - mrSges(3,1) * t351 + mrSges(3,2) * t352 - t389 * t320 - pkin(3) * t417 - pkin(6) * t332 - mrSges(4,1) * t346 + mrSges(4,2) * t347 - Ifges(4,3) * qJDD(3) + (-pkin(3) * t396 + pkin(6) * t333 - t321) * t387;
t304 = mrSges(2,2) * t382 - mrSges(2,3) * t373 - qJ(2) * t312 - t383 * t306 + t385 * t307;
t1 = [-m(1) * g(1) + t399; -m(1) * g(2) + t411; -m(1) * g(3) + m(2) * t382 + t312; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t411 + t386 * t304 - t384 * t305; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t399 + t384 * t304 + t386 * t305; -mrSges(1,1) * g(2) + mrSges(2,1) * t373 + mrSges(1,2) * g(1) - mrSges(2,2) * t374 + pkin(1) * t394 + qJ(2) * t398 + t385 * t306 + t383 * t307;];
tauB = t1;
