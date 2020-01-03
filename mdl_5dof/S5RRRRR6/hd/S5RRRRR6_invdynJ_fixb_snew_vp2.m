% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:05
% EndTime: 2020-01-03 12:15:06
% DurationCPUTime: 1.57s
% Computational Cost: add. (23359->223), mult. (29931->289), div. (0->0), fcn. (19448->10), ass. (0->95)
t372 = qJD(1) + qJD(2);
t368 = t372 ^ 2;
t400 = pkin(3) * t368;
t376 = sin(qJ(3));
t399 = t372 * t376;
t381 = cos(qJ(3));
t398 = t372 * t381;
t378 = sin(qJ(1));
t383 = cos(qJ(1));
t391 = -t383 * g(2) - t378 * g(3);
t357 = qJDD(1) * pkin(1) + t391;
t395 = -t378 * g(2) + t383 * g(3);
t358 = -qJD(1) ^ 2 * pkin(1) + t395;
t377 = sin(qJ(2));
t382 = cos(qJ(2));
t338 = t377 * t357 + t382 * t358;
t370 = qJDD(1) + qJDD(2);
t335 = -t368 * pkin(2) + t370 * pkin(7) + t338;
t397 = t376 * t335;
t396 = qJD(3) * t372;
t352 = t376 * t370 + t381 * t396;
t313 = qJDD(3) * pkin(3) - t352 * pkin(8) - t397 + (pkin(8) * t396 + t376 * t400 - g(1)) * t381;
t325 = -t376 * g(1) + t381 * t335;
t353 = t381 * t370 - t376 * t396;
t361 = qJD(3) * pkin(3) - pkin(8) * t399;
t373 = t381 ^ 2;
t314 = t353 * pkin(8) - qJD(3) * t361 - t373 * t400 + t325;
t375 = sin(qJ(4));
t380 = cos(qJ(4));
t295 = t380 * t313 - t375 * t314;
t346 = (-t375 * t376 + t380 * t381) * t372;
t321 = t346 * qJD(4) + t380 * t352 + t375 * t353;
t347 = (t375 * t381 + t376 * t380) * t372;
t369 = qJDD(3) + qJDD(4);
t371 = qJD(3) + qJD(4);
t290 = (t346 * t371 - t321) * pkin(9) + (t346 * t347 + t369) * pkin(4) + t295;
t296 = t375 * t313 + t380 * t314;
t320 = -t347 * qJD(4) - t375 * t352 + t380 * t353;
t341 = t371 * pkin(4) - t347 * pkin(9);
t342 = t346 ^ 2;
t291 = -t342 * pkin(4) + t320 * pkin(9) - t371 * t341 + t296;
t374 = sin(qJ(5));
t379 = cos(qJ(5));
t288 = t379 * t290 - t374 * t291;
t330 = t379 * t346 - t374 * t347;
t302 = t330 * qJD(5) + t374 * t320 + t379 * t321;
t331 = t374 * t346 + t379 * t347;
t309 = -t330 * mrSges(6,1) + t331 * mrSges(6,2);
t366 = qJD(5) + t371;
t322 = -t366 * mrSges(6,2) + t330 * mrSges(6,3);
t365 = qJDD(5) + t369;
t285 = m(6) * t288 + t365 * mrSges(6,1) - t302 * mrSges(6,3) - t331 * t309 + t366 * t322;
t289 = t374 * t290 + t379 * t291;
t301 = -t331 * qJD(5) + t379 * t320 - t374 * t321;
t323 = t366 * mrSges(6,1) - t331 * mrSges(6,3);
t286 = m(6) * t289 - t365 * mrSges(6,2) + t301 * mrSges(6,3) + t330 * t309 - t366 * t323;
t278 = t379 * t285 + t374 * t286;
t333 = -t346 * mrSges(5,1) + t347 * mrSges(5,2);
t339 = -t371 * mrSges(5,2) + t346 * mrSges(5,3);
t275 = m(5) * t295 + t369 * mrSges(5,1) - t321 * mrSges(5,3) - t347 * t333 + t371 * t339 + t278;
t340 = t371 * mrSges(5,1) - t347 * mrSges(5,3);
t392 = -t374 * t285 + t379 * t286;
t276 = m(5) * t296 - t369 * mrSges(5,2) + t320 * mrSges(5,3) + t346 * t333 - t371 * t340 + t392;
t271 = t380 * t275 + t375 * t276;
t324 = -t381 * g(1) - t397;
t351 = (-mrSges(4,1) * t381 + mrSges(4,2) * t376) * t372;
t359 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t399;
t360 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t398;
t393 = -t375 * t275 + t380 * t276;
t394 = -t376 * (m(4) * t324 + qJDD(3) * mrSges(4,1) - t352 * mrSges(4,3) + qJD(3) * t360 - t351 * t399 + t271) + t381 * (m(4) * t325 - qJDD(3) * mrSges(4,2) + t353 * mrSges(4,3) - qJD(3) * t359 + t351 * t398 + t393);
t337 = t382 * t357 - t377 * t358;
t389 = -t370 * pkin(2) - t337;
t315 = -t353 * pkin(3) + t361 * t399 + (-pkin(8) * t373 - pkin(7)) * t368 + t389;
t293 = -t320 * pkin(4) - t342 * pkin(9) + t347 * t341 + t315;
t390 = m(6) * t293 - t301 * mrSges(6,1) + t302 * mrSges(6,2) - t330 * t322 + t331 * t323;
t304 = Ifges(6,5) * t331 + Ifges(6,6) * t330 + Ifges(6,3) * t366;
t306 = Ifges(6,1) * t331 + Ifges(6,4) * t330 + Ifges(6,5) * t366;
t279 = -mrSges(6,1) * t293 + mrSges(6,3) * t289 + Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t365 - t331 * t304 + t366 * t306;
t305 = Ifges(6,4) * t331 + Ifges(6,2) * t330 + Ifges(6,6) * t366;
t280 = mrSges(6,2) * t293 - mrSges(6,3) * t288 + Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t365 + t330 * t304 - t366 * t305;
t326 = Ifges(5,5) * t347 + Ifges(5,6) * t346 + Ifges(5,3) * t371;
t328 = Ifges(5,1) * t347 + Ifges(5,4) * t346 + Ifges(5,5) * t371;
t267 = -mrSges(5,1) * t315 + mrSges(5,3) * t296 + Ifges(5,4) * t321 + Ifges(5,2) * t320 + Ifges(5,6) * t369 - pkin(4) * t390 + pkin(9) * t392 + t379 * t279 + t374 * t280 - t347 * t326 + t371 * t328;
t327 = Ifges(5,4) * t347 + Ifges(5,2) * t346 + Ifges(5,6) * t371;
t268 = mrSges(5,2) * t315 - mrSges(5,3) * t295 + Ifges(5,1) * t321 + Ifges(5,4) * t320 + Ifges(5,5) * t369 - pkin(9) * t278 - t374 * t279 + t379 * t280 + t346 * t326 - t371 * t327;
t334 = -t368 * pkin(7) + t389;
t343 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t376 + Ifges(4,6) * t381) * t372;
t344 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t376 + Ifges(4,2) * t381) * t372;
t345 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t376 + Ifges(4,4) * t381) * t372;
t386 = m(5) * t315 - t320 * mrSges(5,1) + t321 * mrSges(5,2) - t346 * t339 + t347 * t340 + t390;
t384 = -m(4) * t334 + t353 * mrSges(4,1) - t352 * mrSges(4,2) - t359 * t399 + t360 * t398 - t386;
t388 = -mrSges(3,2) * t338 + t381 * (-mrSges(4,1) * t334 + mrSges(4,3) * t325 + Ifges(4,4) * t352 + Ifges(4,2) * t353 + Ifges(4,6) * qJDD(3) - pkin(3) * t386 + pkin(8) * t393 + qJD(3) * t345 + t380 * t267 + t375 * t268 - t343 * t399) + t376 * (mrSges(4,2) * t334 - mrSges(4,3) * t324 + Ifges(4,1) * t352 + Ifges(4,4) * t353 + Ifges(4,5) * qJDD(3) - pkin(8) * t271 - qJD(3) * t344 - t375 * t267 + t380 * t268 + t343 * t398) + pkin(7) * t394 + pkin(2) * t384 + mrSges(3,1) * t337 + Ifges(3,3) * t370;
t387 = mrSges(6,1) * t288 - mrSges(6,2) * t289 + Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t365 + t331 * t305 - t330 * t306;
t385 = mrSges(5,1) * t295 - mrSges(5,2) * t296 + Ifges(5,5) * t321 + Ifges(5,6) * t320 + Ifges(5,3) * t369 + pkin(4) * t278 + t347 * t327 - t346 * t328 + t387;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t391 - mrSges(2,2) * t395 + pkin(1) * (t377 * (m(3) * t338 - t368 * mrSges(3,1) - t370 * mrSges(3,2) + t394) + t382 * (m(3) * t337 + t370 * mrSges(3,1) - t368 * mrSges(3,2) + t384)) + t388; t388; t385 + Ifges(4,3) * qJDD(3) + (t376 * t344 - t381 * t345) * t372 + Ifges(4,5) * t352 + Ifges(4,6) * t353 + pkin(3) * t271 + mrSges(4,1) * t324 - mrSges(4,2) * t325; t385; t387;];
tauJ = t1;
