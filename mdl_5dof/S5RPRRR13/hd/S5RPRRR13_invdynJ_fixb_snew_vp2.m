% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR13
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:46
% EndTime: 2019-12-31 19:14:47
% DurationCPUTime: 1.07s
% Computational Cost: add. (7961->216), mult. (15365->272), div. (0->0), fcn. (9330->8), ass. (0->89)
t355 = sin(qJ(1));
t359 = cos(qJ(1));
t369 = -t359 * g(1) - t355 * g(2);
t378 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t369;
t377 = (-pkin(1) - pkin(6));
t361 = qJD(1) ^ 2;
t324 = (t377 * t361) - t378;
t358 = cos(qJ(3));
t375 = qJD(1) * qJD(3);
t346 = t358 * t375;
t354 = sin(qJ(3));
t340 = -t354 * qJDD(1) - t346;
t373 = t354 * t375;
t341 = t358 * qJDD(1) - t373;
t302 = (-t341 + t373) * pkin(7) + (-t340 + t346) * pkin(3) + t324;
t372 = t355 * g(1) - t359 * g(2);
t366 = -t361 * qJ(2) + qJDD(2) - t372;
t325 = t377 * qJDD(1) + t366;
t319 = -t358 * g(3) + t354 * t325;
t339 = (t354 * pkin(3) - t358 * pkin(7)) * qJD(1);
t348 = t354 * qJD(1);
t360 = qJD(3) ^ 2;
t305 = -t360 * pkin(3) + qJDD(3) * pkin(7) - t339 * t348 + t319;
t353 = sin(qJ(4));
t357 = cos(qJ(4));
t287 = t357 * t302 - t353 * t305;
t376 = t358 * qJD(1);
t336 = t357 * qJD(3) - t353 * t376;
t314 = t336 * qJD(4) + t353 * qJDD(3) + t357 * t341;
t335 = qJDD(4) - t340;
t337 = t353 * qJD(3) + t357 * t376;
t345 = t348 + qJD(4);
t284 = (t336 * t345 - t314) * pkin(8) + (t336 * t337 + t335) * pkin(4) + t287;
t288 = t353 * t302 + t357 * t305;
t313 = -t337 * qJD(4) + t357 * qJDD(3) - t353 * t341;
t323 = t345 * pkin(4) - t337 * pkin(8);
t334 = t336 ^ 2;
t285 = -t334 * pkin(4) + t313 * pkin(8) - t345 * t323 + t288;
t352 = sin(qJ(5));
t356 = cos(qJ(5));
t282 = t356 * t284 - t352 * t285;
t315 = t356 * t336 - t352 * t337;
t294 = t315 * qJD(5) + t352 * t313 + t356 * t314;
t316 = t352 * t336 + t356 * t337;
t299 = -t315 * mrSges(6,1) + t316 * mrSges(6,2);
t344 = qJD(5) + t345;
t306 = -t344 * mrSges(6,2) + t315 * mrSges(6,3);
t333 = qJDD(5) + t335;
t279 = m(6) * t282 + t333 * mrSges(6,1) - t294 * mrSges(6,3) - t316 * t299 + t344 * t306;
t283 = t352 * t284 + t356 * t285;
t293 = -t316 * qJD(5) + t356 * t313 - t352 * t314;
t307 = t344 * mrSges(6,1) - t316 * mrSges(6,3);
t280 = m(6) * t283 - t333 * mrSges(6,2) + t293 * mrSges(6,3) + t315 * t299 - t344 * t307;
t272 = t356 * t279 + t352 * t280;
t317 = -t336 * mrSges(5,1) + t337 * mrSges(5,2);
t320 = -t345 * mrSges(5,2) + t336 * mrSges(5,3);
t270 = m(5) * t287 + t335 * mrSges(5,1) - t314 * mrSges(5,3) - t337 * t317 + t345 * t320 + t272;
t321 = t345 * mrSges(5,1) - t337 * mrSges(5,3);
t370 = -t352 * t279 + t356 * t280;
t271 = m(5) * t288 - t335 * mrSges(5,2) + t313 * mrSges(5,3) + t336 * t317 - t345 * t321 + t370;
t266 = t357 * t270 + t353 * t271;
t371 = -t353 * t270 + t357 * t271;
t318 = t354 * g(3) + t358 * t325;
t338 = (t354 * mrSges(4,1) + t358 * mrSges(4,2)) * qJD(1);
t342 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t348;
t343 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t376;
t304 = -qJDD(3) * pkin(3) - t360 * pkin(7) + t339 * t376 - t318;
t286 = -t313 * pkin(4) - t334 * pkin(8) + t337 * t323 + t304;
t365 = m(6) * t286 - t293 * mrSges(6,1) + t294 * mrSges(6,2) - t315 * t306 + t316 * t307;
t363 = -m(5) * t304 + t313 * mrSges(5,1) - t314 * mrSges(5,2) + t336 * t320 - t337 * t321 - t365;
t368 = t354 * (m(4) * t319 - qJDD(3) * mrSges(4,2) + t340 * mrSges(4,3) - qJD(3) * t343 - t338 * t348 + t371) + t358 * (m(4) * t318 + qJDD(3) * mrSges(4,1) - t341 * mrSges(4,3) + qJD(3) * t342 - t338 * t376 + t363);
t296 = Ifges(6,4) * t316 + Ifges(6,2) * t315 + Ifges(6,6) * t344;
t297 = Ifges(6,1) * t316 + Ifges(6,4) * t315 + Ifges(6,5) * t344;
t364 = -mrSges(6,1) * t282 + mrSges(6,2) * t283 - Ifges(6,5) * t294 - Ifges(6,6) * t293 - Ifges(6,3) * t333 - t316 * t296 + t315 * t297;
t309 = Ifges(5,4) * t337 + Ifges(5,2) * t336 + Ifges(5,6) * t345;
t310 = Ifges(5,1) * t337 + Ifges(5,4) * t336 + Ifges(5,5) * t345;
t362 = mrSges(5,1) * t287 - mrSges(5,2) * t288 + Ifges(5,5) * t314 + Ifges(5,6) * t313 + Ifges(5,3) * t335 + pkin(4) * t272 + t337 * t309 - t336 * t310 - t364;
t332 = (Ifges(4,5) * qJD(3)) + (t358 * Ifges(4,1) - t354 * Ifges(4,4)) * qJD(1);
t331 = (Ifges(4,6) * qJD(3)) + (t358 * Ifges(4,4) - t354 * Ifges(4,2)) * qJD(1);
t327 = -qJDD(1) * pkin(1) + t366;
t326 = t361 * pkin(1) + t378;
t308 = Ifges(5,5) * t337 + Ifges(5,6) * t336 + Ifges(5,3) * t345;
t295 = Ifges(6,5) * t316 + Ifges(6,6) * t315 + Ifges(6,3) * t344;
t274 = mrSges(6,2) * t286 - mrSges(6,3) * t282 + Ifges(6,1) * t294 + Ifges(6,4) * t293 + Ifges(6,5) * t333 + t315 * t295 - t344 * t296;
t273 = -mrSges(6,1) * t286 + mrSges(6,3) * t283 + Ifges(6,4) * t294 + Ifges(6,2) * t293 + Ifges(6,6) * t333 - t316 * t295 + t344 * t297;
t264 = mrSges(5,2) * t304 - mrSges(5,3) * t287 + Ifges(5,1) * t314 + Ifges(5,4) * t313 + Ifges(5,5) * t335 - pkin(8) * t272 - t352 * t273 + t356 * t274 + t336 * t308 - t345 * t309;
t263 = m(3) * t327 + qJDD(1) * mrSges(3,2) - (t361 * mrSges(3,3)) + t368;
t262 = -mrSges(5,1) * t304 + mrSges(5,3) * t288 + Ifges(5,4) * t314 + Ifges(5,2) * t313 + Ifges(5,6) * t335 - pkin(4) * t365 + pkin(8) * t370 + t356 * t273 + t352 * t274 - t337 * t308 + t345 * t310;
t1 = [mrSges(2,1) * t372 - mrSges(2,2) * t369 + mrSges(3,2) * t327 - mrSges(3,3) * t326 + t358 * (mrSges(4,2) * t324 - mrSges(4,3) * t318 + Ifges(4,1) * t341 + Ifges(4,4) * t340 + Ifges(4,5) * qJDD(3) - pkin(7) * t266 - qJD(3) * t331 - t353 * t262 + t357 * t264) - t354 * (-mrSges(4,1) * t324 + mrSges(4,3) * t319 + Ifges(4,4) * t341 + Ifges(4,2) * t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t266 + qJD(3) * t332 - t362) - pkin(6) * t368 - pkin(1) * t263 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t326 + m(4) * t324 - t340 * mrSges(4,1) + t361 * mrSges(3,2) + t341 * mrSges(4,2) + t266 + qJDD(1) * mrSges(3,3) + (t342 * t354 + t343 * t358) * qJD(1)) * qJ(2); t263; Ifges(4,5) * t341 + Ifges(4,6) * t340 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t318 - mrSges(4,2) * t319 + t353 * t264 + t357 * t262 + pkin(3) * t363 + pkin(7) * t371 + (t358 * t331 + t354 * t332) * qJD(1); t362; -t364;];
tauJ = t1;
