% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:51
% EndTime: 2019-12-05 18:33:52
% DurationCPUTime: 1.38s
% Computational Cost: add. (17827->191), mult. (24620->246), div. (0->0), fcn. (16678->10), ass. (0->90)
t371 = qJD(1) + qJD(2);
t365 = t371 ^ 2;
t373 = cos(pkin(9));
t404 = pkin(3) * t373;
t372 = sin(pkin(9));
t403 = mrSges(4,2) * t372;
t369 = t373 ^ 2;
t401 = t365 * t369;
t367 = qJDD(1) + qJDD(2);
t400 = t367 * t373;
t377 = sin(qJ(1));
t381 = cos(qJ(1));
t398 = t381 * g(2) + t377 * g(3);
t353 = qJDD(1) * pkin(1) + t398;
t394 = t377 * g(2) - g(3) * t381;
t354 = -qJD(1) ^ 2 * pkin(1) + t394;
t376 = sin(qJ(2));
t380 = cos(qJ(2));
t341 = t376 * t353 + t380 * t354;
t338 = -pkin(2) * t365 + qJ(3) * t367 + t341;
t397 = qJD(3) * t371;
t395 = -t373 * g(1) - 0.2e1 * t372 * t397;
t319 = (-pkin(7) * t367 + t365 * t404 - t338) * t372 + t395;
t323 = -t372 * g(1) + (t338 + 0.2e1 * t397) * t373;
t320 = -pkin(3) * t401 + pkin(7) * t400 + t323;
t375 = sin(qJ(4));
t379 = cos(qJ(4));
t301 = t379 * t319 - t375 * t320;
t387 = t372 * t379 + t373 * t375;
t386 = -t372 * t375 + t373 * t379;
t346 = t386 * t371;
t396 = t346 * qJD(4);
t337 = t367 * t387 + t396;
t347 = t387 * t371;
t297 = (-t337 + t396) * pkin(8) + (t346 * t347 + qJDD(4)) * pkin(4) + t301;
t302 = t375 * t319 + t379 * t320;
t336 = -t347 * qJD(4) + t367 * t386;
t344 = qJD(4) * pkin(4) - pkin(8) * t347;
t345 = t346 ^ 2;
t298 = -pkin(4) * t345 + pkin(8) * t336 - qJD(4) * t344 + t302;
t374 = sin(qJ(5));
t378 = cos(qJ(5));
t295 = t297 * t378 - t298 * t374;
t329 = t346 * t378 - t347 * t374;
t309 = qJD(5) * t329 + t336 * t374 + t337 * t378;
t330 = t346 * t374 + t347 * t378;
t315 = -mrSges(6,1) * t329 + mrSges(6,2) * t330;
t370 = qJD(4) + qJD(5);
t324 = -mrSges(6,2) * t370 + mrSges(6,3) * t329;
t366 = qJDD(4) + qJDD(5);
t292 = m(6) * t295 + mrSges(6,1) * t366 - mrSges(6,3) * t309 - t315 * t330 + t324 * t370;
t296 = t297 * t374 + t298 * t378;
t308 = -qJD(5) * t330 + t336 * t378 - t337 * t374;
t325 = mrSges(6,1) * t370 - mrSges(6,3) * t330;
t293 = m(6) * t296 - mrSges(6,2) * t366 + mrSges(6,3) * t308 + t315 * t329 - t325 * t370;
t284 = t378 * t292 + t374 * t293;
t334 = -mrSges(5,1) * t346 + mrSges(5,2) * t347;
t342 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t346;
t282 = m(5) * t301 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t337 + qJD(4) * t342 - t334 * t347 + t284;
t343 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t347;
t391 = -t292 * t374 + t378 * t293;
t283 = m(5) * t302 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t336 - qJD(4) * t343 + t334 * t346 + t391;
t399 = t379 * t282 + t375 * t283;
t322 = -t372 * t338 + t395;
t388 = mrSges(4,3) * t367 + (-mrSges(4,1) * t373 + t403) * t365;
t392 = -t375 * t282 + t379 * t283;
t393 = -t372 * (m(4) * t322 - t372 * t388 + t399) + t373 * (m(4) * t323 + t373 * t388 + t392);
t340 = t380 * t353 - t376 * t354;
t390 = qJDD(3) - t340;
t368 = t372 ^ 2;
t321 = (-pkin(2) - t404) * t367 + (-qJ(3) + (-t368 - t369) * pkin(7)) * t365 + t390;
t300 = -t336 * pkin(4) - t345 * pkin(8) + t347 * t344 + t321;
t389 = m(6) * t300 - t308 * mrSges(6,1) + t309 * mrSges(6,2) - t329 * t324 + t330 * t325;
t310 = Ifges(6,5) * t330 + Ifges(6,6) * t329 + Ifges(6,3) * t370;
t312 = Ifges(6,1) * t330 + Ifges(6,4) * t329 + Ifges(6,5) * t370;
t285 = -mrSges(6,1) * t300 + mrSges(6,3) * t296 + Ifges(6,4) * t309 + Ifges(6,2) * t308 + Ifges(6,6) * t366 - t310 * t330 + t312 * t370;
t311 = Ifges(6,4) * t330 + Ifges(6,2) * t329 + Ifges(6,6) * t370;
t286 = mrSges(6,2) * t300 - mrSges(6,3) * t295 + Ifges(6,1) * t309 + Ifges(6,4) * t308 + Ifges(6,5) * t366 + t310 * t329 - t311 * t370;
t326 = Ifges(5,5) * t347 + Ifges(5,6) * t346 + Ifges(5,3) * qJD(4);
t328 = Ifges(5,1) * t347 + Ifges(5,4) * t346 + Ifges(5,5) * qJD(4);
t275 = -mrSges(5,1) * t321 + mrSges(5,3) * t302 + Ifges(5,4) * t337 + Ifges(5,2) * t336 + Ifges(5,6) * qJDD(4) - pkin(4) * t389 + pkin(8) * t391 + qJD(4) * t328 + t378 * t285 + t374 * t286 - t347 * t326;
t327 = Ifges(5,4) * t347 + Ifges(5,2) * t346 + Ifges(5,6) * qJD(4);
t276 = mrSges(5,2) * t321 - mrSges(5,3) * t301 + Ifges(5,1) * t337 + Ifges(5,4) * t336 + Ifges(5,5) * qJDD(4) - pkin(8) * t284 - qJD(4) * t327 - t285 * t374 + t286 * t378 + t326 * t346;
t335 = -t367 * pkin(2) - t365 * qJ(3) + t390;
t383 = m(5) * t321 - t336 * mrSges(5,1) + t337 * mrSges(5,2) - t346 * t342 + t347 * t343 + t389;
t382 = -m(4) * t335 + mrSges(4,1) * t400 - t383 + (t365 * t368 + t401) * mrSges(4,3);
t288 = t367 * t403 - t382;
t385 = -mrSges(3,2) * t341 + t373 * (-mrSges(4,1) * t335 + mrSges(4,3) * t323 + t375 * t276 + t379 * t275 - pkin(3) * t383 + pkin(7) * t392 + (Ifges(4,4) * t372 + Ifges(4,2) * t373) * t367) + t372 * (mrSges(4,2) * t335 - mrSges(4,3) * t322 + t379 * t276 - t375 * t275 - pkin(7) * t399 + (Ifges(4,1) * t372 + Ifges(4,4) * t373) * t367) + qJ(3) * t393 - pkin(2) * t288 + mrSges(3,1) * t340 + Ifges(3,3) * t367;
t384 = mrSges(6,1) * t295 - mrSges(6,2) * t296 + Ifges(6,5) * t309 + Ifges(6,6) * t308 + Ifges(6,3) * t366 + t330 * t311 - t312 * t329;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t398 - mrSges(2,2) * t394 + pkin(1) * (t376 * (m(3) * t341 - mrSges(3,1) * t365 - mrSges(3,2) * t367 + t393) + t380 * ((mrSges(3,1) - t403) * t367 - t365 * mrSges(3,2) + m(3) * t340 + t382)) + t385; t385; t288; mrSges(5,1) * t301 - mrSges(5,2) * t302 + Ifges(5,5) * t337 + Ifges(5,6) * t336 + Ifges(5,3) * qJDD(4) + pkin(4) * t284 + t327 * t347 - t328 * t346 + t384; t384;];
tauJ = t1;
