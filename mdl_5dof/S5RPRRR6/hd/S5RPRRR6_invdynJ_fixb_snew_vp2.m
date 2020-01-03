% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:12
% EndTime: 2019-12-31 19:01:13
% DurationCPUTime: 1.15s
% Computational Cost: add. (8453->219), mult. (16546->284), div. (0->0), fcn. (10514->10), ass. (0->93)
t361 = sin(qJ(4));
t362 = sin(qJ(3));
t365 = cos(qJ(4));
t366 = cos(qJ(3));
t336 = (t362 * t361 - t366 * t365) * qJD(1);
t363 = sin(qJ(1));
t367 = cos(qJ(1));
t380 = t363 * g(1) - g(2) * t367;
t341 = qJDD(1) * pkin(1) + t380;
t368 = qJD(1) ^ 2;
t376 = -g(1) * t367 - g(2) * t363;
t343 = -pkin(1) * t368 + t376;
t358 = sin(pkin(9));
t359 = cos(pkin(9));
t324 = t358 * t341 + t359 * t343;
t320 = -pkin(2) * t368 + qJDD(1) * pkin(6) + t324;
t357 = -g(3) + qJDD(2);
t308 = -t320 * t362 + t366 * t357;
t382 = qJD(1) * qJD(3);
t381 = t366 * t382;
t344 = qJDD(1) * t362 + t381;
t295 = (-t344 + t381) * pkin(7) + (t362 * t366 * t368 + qJDD(3)) * pkin(3) + t308;
t309 = t366 * t320 + t362 * t357;
t345 = qJDD(1) * t366 - t362 * t382;
t383 = t362 * qJD(1);
t348 = qJD(3) * pkin(3) - pkin(7) * t383;
t356 = t366 ^ 2;
t296 = -pkin(3) * t356 * t368 + pkin(7) * t345 - qJD(3) * t348 + t309;
t289 = t361 * t295 + t365 * t296;
t337 = (t366 * t361 + t362 * t365) * qJD(1);
t310 = -qJD(4) * t337 - t344 * t361 + t345 * t365;
t321 = mrSges(5,1) * t336 + mrSges(5,2) * t337;
t355 = qJD(3) + qJD(4);
t328 = mrSges(5,1) * t355 - mrSges(5,3) * t337;
t354 = qJDD(3) + qJDD(4);
t322 = pkin(4) * t336 - pkin(8) * t337;
t353 = t355 ^ 2;
t285 = -pkin(4) * t353 + pkin(8) * t354 - t322 * t336 + t289;
t323 = t341 * t359 - t358 * t343;
t374 = -qJDD(1) * pkin(2) - t323;
t301 = -pkin(3) * t345 + t348 * t383 + (-pkin(7) * t356 - pkin(6)) * t368 + t374;
t311 = -qJD(4) * t336 + t344 * t365 + t345 * t361;
t286 = (t336 * t355 - t311) * pkin(8) + (t337 * t355 - t310) * pkin(4) + t301;
t360 = sin(qJ(5));
t364 = cos(qJ(5));
t282 = -t285 * t360 + t286 * t364;
t325 = -t337 * t360 + t355 * t364;
t292 = qJD(5) * t325 + t311 * t364 + t354 * t360;
t326 = t337 * t364 + t355 * t360;
t302 = -mrSges(6,1) * t325 + mrSges(6,2) * t326;
t307 = qJDD(5) - t310;
t329 = qJD(5) + t336;
t312 = -mrSges(6,2) * t329 + mrSges(6,3) * t325;
t279 = m(6) * t282 + mrSges(6,1) * t307 - mrSges(6,3) * t292 - t302 * t326 + t312 * t329;
t283 = t285 * t364 + t286 * t360;
t291 = -qJD(5) * t326 - t311 * t360 + t354 * t364;
t313 = mrSges(6,1) * t329 - mrSges(6,3) * t326;
t280 = m(6) * t283 - mrSges(6,2) * t307 + mrSges(6,3) * t291 + t302 * t325 - t313 * t329;
t377 = -t279 * t360 + t364 * t280;
t267 = m(5) * t289 - mrSges(5,2) * t354 + mrSges(5,3) * t310 - t321 * t336 - t328 * t355 + t377;
t288 = t295 * t365 - t296 * t361;
t327 = -mrSges(5,2) * t355 - mrSges(5,3) * t336;
t284 = -pkin(4) * t354 - pkin(8) * t353 + t322 * t337 - t288;
t373 = -m(6) * t284 + t291 * mrSges(6,1) - mrSges(6,2) * t292 + t325 * t312 - t313 * t326;
t275 = m(5) * t288 + mrSges(5,1) * t354 - mrSges(5,3) * t311 - t321 * t337 + t327 * t355 + t373;
t264 = t361 * t267 + t365 * t275;
t269 = t364 * t279 + t360 * t280;
t384 = qJD(1) * t366;
t342 = (-t366 * mrSges(4,1) + t362 * mrSges(4,2)) * qJD(1);
t347 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t384;
t262 = m(4) * t308 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t344 + qJD(3) * t347 - t342 * t383 + t264;
t346 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t383;
t378 = t365 * t267 - t275 * t361;
t263 = m(4) * t309 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t345 - qJD(3) * t346 + t342 * t384 + t378;
t379 = -t362 * t262 + t366 * t263;
t372 = m(5) * t301 - t310 * mrSges(5,1) + mrSges(5,2) * t311 + t336 * t327 + t328 * t337 + t269;
t297 = Ifges(6,5) * t326 + Ifges(6,6) * t325 + Ifges(6,3) * t329;
t299 = Ifges(6,1) * t326 + Ifges(6,4) * t325 + Ifges(6,5) * t329;
t272 = -mrSges(6,1) * t284 + mrSges(6,3) * t283 + Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t307 - t297 * t326 + t299 * t329;
t298 = Ifges(6,4) * t326 + Ifges(6,2) * t325 + Ifges(6,6) * t329;
t273 = mrSges(6,2) * t284 - mrSges(6,3) * t282 + Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t307 + t297 * t325 - t298 * t329;
t316 = Ifges(5,4) * t337 - Ifges(5,2) * t336 + Ifges(5,6) * t355;
t317 = Ifges(5,1) * t337 - Ifges(5,4) * t336 + Ifges(5,5) * t355;
t371 = mrSges(5,1) * t288 - mrSges(5,2) * t289 + Ifges(5,5) * t311 + Ifges(5,6) * t310 + Ifges(5,3) * t354 + pkin(4) * t373 + pkin(8) * t377 + t364 * t272 + t360 * t273 + t337 * t316 + t317 * t336;
t370 = mrSges(6,1) * t282 - mrSges(6,2) * t283 + Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t307 + t298 * t326 - t299 * t325;
t319 = -pkin(6) * t368 + t374;
t369 = -m(4) * t319 + t345 * mrSges(4,1) - mrSges(4,2) * t344 - t346 * t383 + t347 * t384 - t372;
t335 = Ifges(4,5) * qJD(3) + (t362 * Ifges(4,1) + t366 * Ifges(4,4)) * qJD(1);
t334 = Ifges(4,6) * qJD(3) + (t362 * Ifges(4,4) + t366 * Ifges(4,2)) * qJD(1);
t315 = Ifges(5,5) * t337 - Ifges(5,6) * t336 + Ifges(5,3) * t355;
t260 = -mrSges(5,1) * t301 + mrSges(5,3) * t289 + Ifges(5,4) * t311 + Ifges(5,2) * t310 + Ifges(5,6) * t354 - pkin(4) * t269 - t315 * t337 + t317 * t355 - t370;
t259 = mrSges(5,2) * t301 - mrSges(5,3) * t288 + Ifges(5,1) * t311 + Ifges(5,4) * t310 + Ifges(5,5) * t354 - pkin(8) * t269 - t272 * t360 + t273 * t364 - t315 * t336 - t316 * t355;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t380 - mrSges(2,2) * t376 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t323 - mrSges(3,2) * t324 + t362 * (mrSges(4,2) * t319 - mrSges(4,3) * t308 + Ifges(4,1) * t344 + Ifges(4,4) * t345 + Ifges(4,5) * qJDD(3) - pkin(7) * t264 - qJD(3) * t334 + t259 * t365 - t260 * t361) + t366 * (-mrSges(4,1) * t319 + mrSges(4,3) * t309 + Ifges(4,4) * t344 + Ifges(4,2) * t345 + Ifges(4,6) * qJDD(3) - pkin(3) * t372 + pkin(7) * t378 + qJD(3) * t335 + t361 * t259 + t365 * t260) + pkin(2) * t369 + pkin(6) * t379 + pkin(1) * (t358 * (m(3) * t324 - mrSges(3,1) * t368 - qJDD(1) * mrSges(3,2) + t379) + t359 * (m(3) * t323 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t368 + t369)); m(3) * t357 + t262 * t366 + t263 * t362; t371 + mrSges(4,1) * t308 - mrSges(4,2) * t309 + Ifges(4,5) * t344 + Ifges(4,6) * t345 + pkin(3) * t264 + Ifges(4,3) * qJDD(3) + (t362 * t334 - t366 * t335) * qJD(1); t371; t370;];
tauJ = t1;
