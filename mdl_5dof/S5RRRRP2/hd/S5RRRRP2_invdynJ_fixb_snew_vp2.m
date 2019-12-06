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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:47:45
% EndTime: 2019-12-05 18:47:47
% DurationCPUTime: 1.18s
% Computational Cost: add. (9308->201), mult. (11771->247), div. (0->0), fcn. (6989->8), ass. (0->84)
t405 = Ifges(5,4) + Ifges(6,4);
t412 = Ifges(5,2) + Ifges(6,2);
t409 = Ifges(5,6) + Ifges(6,6);
t411 = Ifges(5,1) + Ifges(6,1);
t410 = Ifges(5,5) + Ifges(6,5);
t408 = Ifges(5,3) + Ifges(6,3);
t378 = qJD(1) + qJD(2);
t380 = sin(qJ(4));
t381 = sin(qJ(3));
t384 = cos(qJ(4));
t385 = cos(qJ(3));
t351 = (-t380 * t381 + t384 * t385) * t378;
t352 = (t380 * t385 + t381 * t384) * t378;
t377 = qJD(3) + qJD(4);
t407 = t412 * t351 + t405 * t352 + t409 * t377;
t374 = t378 ^ 2;
t406 = pkin(3) * t374;
t404 = t378 * t381;
t403 = t378 * t385;
t383 = sin(qJ(1));
t387 = cos(qJ(1));
t399 = t387 * g(2) + t383 * g(3);
t362 = qJDD(1) * pkin(1) + t399;
t395 = t383 * g(2) - t387 * g(3);
t363 = -qJD(1) ^ 2 * pkin(1) + t395;
t382 = sin(qJ(2));
t386 = cos(qJ(2));
t340 = t382 * t362 + t386 * t363;
t376 = qJDD(1) + qJDD(2);
t337 = -t374 * pkin(2) + t376 * pkin(7) + t340;
t402 = t381 * t337;
t398 = qJD(3) * t378;
t357 = t381 * t376 + t385 * t398;
t306 = qJDD(3) * pkin(3) - t357 * pkin(8) - t402 + (pkin(8) * t398 + t381 * t406 - g(1)) * t385;
t323 = -t381 * g(1) + t385 * t337;
t358 = t385 * t376 - t381 * t398;
t366 = qJD(3) * pkin(3) - pkin(8) * t404;
t379 = t385 ^ 2;
t307 = t358 * pkin(8) - qJD(3) * t366 - t379 * t406 + t323;
t301 = t384 * t306 - t380 * t307;
t321 = t351 * qJD(4) + t384 * t357 + t380 * t358;
t334 = -t351 * mrSges(6,1) + t352 * mrSges(6,2);
t335 = -t351 * mrSges(5,1) + t352 * mrSges(5,2);
t343 = -t377 * mrSges(5,2) + t351 * mrSges(5,3);
t375 = qJDD(3) + qJDD(4);
t295 = -0.2e1 * qJD(5) * t352 + (t351 * t377 - t321) * qJ(5) + (t351 * t352 + t375) * pkin(4) + t301;
t342 = -t377 * mrSges(6,2) + t351 * mrSges(6,3);
t397 = m(6) * t295 + t375 * mrSges(6,1) + t377 * t342;
t286 = m(5) * t301 + t375 * mrSges(5,1) + t377 * t343 + (-t334 - t335) * t352 + (-mrSges(5,3) - mrSges(6,3)) * t321 + t397;
t302 = t380 * t306 + t384 * t307;
t320 = -t352 * qJD(4) - t380 * t357 + t384 * t358;
t345 = t377 * mrSges(6,1) - t352 * mrSges(6,3);
t346 = t377 * mrSges(5,1) - t352 * mrSges(5,3);
t344 = t377 * pkin(4) - t352 * qJ(5);
t347 = t351 ^ 2;
t297 = -t347 * pkin(4) + t320 * qJ(5) + 0.2e1 * qJD(5) * t351 - t377 * t344 + t302;
t396 = m(6) * t297 + t320 * mrSges(6,3) + t351 * t334;
t289 = m(5) * t302 + t320 * mrSges(5,3) + t351 * t335 + (-t345 - t346) * t377 + (-mrSges(5,2) - mrSges(6,2)) * t375 + t396;
t283 = t384 * t286 + t380 * t289;
t401 = -t409 * t351 - t410 * t352 - t408 * t377;
t400 = -t405 * t351 - t411 * t352 - t410 * t377;
t322 = -t385 * g(1) - t402;
t356 = (-mrSges(4,1) * t385 + mrSges(4,2) * t381) * t378;
t364 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t404;
t365 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t403;
t393 = -t380 * t286 + t384 * t289;
t394 = -t381 * (m(4) * t322 + qJDD(3) * mrSges(4,1) - t357 * mrSges(4,3) + qJD(3) * t365 - t356 * t404 + t283) + t385 * (m(4) * t323 - qJDD(3) * mrSges(4,2) + t358 * mrSges(4,3) - qJD(3) * t364 + t356 * t403 + t393);
t339 = t386 * t362 - t382 * t363;
t392 = -t376 * pkin(2) - t339;
t308 = -t358 * pkin(3) + t366 * t404 + (-pkin(8) * t379 - pkin(7)) * t374 + t392;
t299 = -t320 * pkin(4) - t347 * qJ(5) + t352 * t344 + qJDD(5) + t308;
t292 = m(6) * t299 - t320 * mrSges(6,1) + t321 * mrSges(6,2) - t351 * t342 + t352 * t345;
t279 = -mrSges(5,1) * t308 + mrSges(5,3) * t302 - mrSges(6,1) * t299 + mrSges(6,3) * t297 - pkin(4) * t292 + qJ(5) * t396 + (-qJ(5) * t345 - t400) * t377 + (-qJ(5) * mrSges(6,2) + t409) * t375 + t401 * t352 + t405 * t321 + t412 * t320;
t291 = -t321 * mrSges(6,3) - t352 * t334 + t397;
t280 = mrSges(5,2) * t308 + mrSges(6,2) * t299 - mrSges(5,3) * t301 - mrSges(6,3) * t295 - qJ(5) * t291 + t405 * t320 + t411 * t321 - t401 * t351 + t410 * t375 - t407 * t377;
t336 = -t374 * pkin(7) + t392;
t348 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t381 + Ifges(4,6) * t385) * t378;
t349 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t381 + Ifges(4,2) * t385) * t378;
t350 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t381 + Ifges(4,4) * t385) * t378;
t390 = m(5) * t308 - t320 * mrSges(5,1) + t321 * mrSges(5,2) - t351 * t343 + t352 * t346 + t292;
t388 = -m(4) * t336 + t358 * mrSges(4,1) - t357 * mrSges(4,2) - t364 * t404 + t365 * t403 - t390;
t391 = -mrSges(3,2) * t340 + t385 * (-mrSges(4,1) * t336 + mrSges(4,3) * t323 + Ifges(4,4) * t357 + Ifges(4,2) * t358 + Ifges(4,6) * qJDD(3) - pkin(3) * t390 + pkin(8) * t393 + qJD(3) * t350 + t384 * t279 + t380 * t280 - t348 * t404) + t381 * (mrSges(4,2) * t336 - mrSges(4,3) * t322 + Ifges(4,1) * t357 + Ifges(4,4) * t358 + Ifges(4,5) * qJDD(3) - pkin(8) * t283 - qJD(3) * t349 - t380 * t279 + t384 * t280 + t348 * t403) + pkin(7) * t394 + pkin(2) * t388 + mrSges(3,1) * t339 + Ifges(3,3) * t376;
t389 = mrSges(5,1) * t301 + mrSges(6,1) * t295 - mrSges(5,2) * t302 - mrSges(6,2) * t297 + pkin(4) * t291 + t409 * t320 + t410 * t321 + t400 * t351 + t407 * t352 + t408 * t375;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t399 - mrSges(2,2) * t395 + pkin(1) * (t382 * (m(3) * t340 - t374 * mrSges(3,1) - t376 * mrSges(3,2) + t394) + t386 * (m(3) * t339 + t376 * mrSges(3,1) - t374 * mrSges(3,2) + t388)) + t391; t391; t389 + Ifges(4,3) * qJDD(3) + (t381 * t349 - t385 * t350) * t378 + Ifges(4,5) * t357 + Ifges(4,6) * t358 + mrSges(4,1) * t322 - mrSges(4,2) * t323 + pkin(3) * t283; t389; t292;];
tauJ = t1;
