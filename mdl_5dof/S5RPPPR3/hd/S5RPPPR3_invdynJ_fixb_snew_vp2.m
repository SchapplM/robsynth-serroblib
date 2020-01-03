% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:56
% EndTime: 2019-12-31 17:43:57
% DurationCPUTime: 0.63s
% Computational Cost: add. (1782->146), mult. (3759->188), div. (0->0), fcn. (2109->8), ass. (0->73)
t336 = sin(pkin(8));
t333 = t336 ^ 2;
t344 = qJD(1) ^ 2;
t338 = cos(pkin(8));
t334 = t338 ^ 2;
t366 = t334 * t344;
t373 = -t333 * t344 - t366;
t372 = 0.2e1 * t338;
t341 = sin(qJ(1));
t343 = cos(qJ(1));
t356 = t341 * g(1) - g(2) * t343;
t320 = qJDD(1) * pkin(1) + t356;
t351 = -g(1) * t343 - g(2) * t341;
t321 = -pkin(1) * t344 + t351;
t337 = sin(pkin(7));
t339 = cos(pkin(7));
t304 = t337 * t320 + t339 * t321;
t295 = -pkin(2) * t344 + qJDD(1) * qJ(3) + t304;
t335 = -g(3) + qJDD(2);
t361 = (qJD(1) * qJD(3));
t287 = t335 * t338 + (-t295 - (2 * t361)) * t336;
t316 = (-t338 * mrSges(5,1) - t336 * mrSges(5,3)) * qJD(1);
t369 = mrSges(4,2) * t336;
t371 = (t316 + (-t338 * mrSges(4,1) + t369) * qJD(1)) * qJD(1) + (mrSges(5,2) + mrSges(4,3)) * qJDD(1);
t368 = qJ(3) * t344;
t365 = t336 * qJ(4);
t350 = -t338 * pkin(3) - t365;
t315 = t350 * qJD(1);
t362 = t336 * qJD(1);
t284 = t315 * t362 + qJDD(4) - t287;
t280 = (-pkin(4) * t338 * t344 - pkin(6) * qJDD(1)) * t336 + t284;
t288 = t338 * t295 + t336 * t335 + t361 * t372;
t285 = t338 * qJD(1) * t315 + t288;
t359 = t338 * qJDD(1);
t281 = -pkin(4) * t366 - pkin(6) * t359 + t285;
t340 = sin(qJ(5));
t342 = cos(qJ(5));
t278 = t280 * t342 - t281 * t340;
t347 = -t336 * t340 - t338 * t342;
t309 = t347 * qJD(1);
t348 = t336 * t342 - t338 * t340;
t310 = t348 * qJD(1);
t298 = -mrSges(6,1) * t309 + mrSges(6,2) * t310;
t302 = qJD(5) * t309 + t348 * qJDD(1);
t305 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t309;
t276 = m(6) * t278 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t302 + qJD(5) * t305 - t298 * t310;
t279 = t280 * t340 + t281 * t342;
t301 = -qJD(5) * t310 + t347 * qJDD(1);
t306 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t310;
t277 = m(6) * t279 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t301 - qJD(5) * t306 + t298 * t309;
t364 = t342 * t276 + t340 * t277;
t303 = t339 * t320 - t337 * t321;
t360 = qJDD(1) * t336;
t357 = -qJDD(3) + t303;
t352 = m(5) * t284 + t364;
t267 = m(4) * t287 - t371 * t336 - t352;
t354 = -t276 * t340 + t342 * t277;
t268 = m(4) * t288 + m(5) * t285 + t371 * t338 + t354;
t355 = -t267 * t336 + t338 * t268;
t346 = -0.2e1 * qJD(4) * t362 - t357;
t283 = (qJ(3) + (-t333 - t334) * pkin(6)) * t344 + (t365 + pkin(2) + (pkin(3) + pkin(4)) * t338) * qJDD(1) - t346;
t349 = -m(6) * t283 + t301 * mrSges(6,1) - t302 * mrSges(6,2) + t309 * t305 - t310 * t306;
t286 = -t368 + (-pkin(2) + t350) * qJDD(1) + t346;
t272 = m(5) * t286 - mrSges(5,1) * t359 + t373 * mrSges(5,2) - mrSges(5,3) * t360 + t349;
t291 = -qJDD(1) * pkin(2) - t357 - t368;
t345 = m(4) * t291 - mrSges(4,1) * t359 + t373 * mrSges(4,3) + t272;
t294 = Ifges(6,1) * t310 + Ifges(6,4) * t309 + Ifges(6,5) * qJD(5);
t293 = Ifges(6,4) * t310 + Ifges(6,2) * t309 + Ifges(6,6) * qJD(5);
t292 = Ifges(6,5) * t310 + Ifges(6,6) * t309 + Ifges(6,3) * qJD(5);
t271 = mrSges(4,2) * t360 + t345;
t270 = mrSges(6,2) * t283 - mrSges(6,3) * t278 + Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * qJDD(5) - qJD(5) * t293 + t292 * t309;
t269 = -mrSges(6,1) * t283 + mrSges(6,3) * t279 + Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * qJDD(5) + qJD(5) * t294 - t292 * t310;
t1 = [mrSges(2,1) * t356 - mrSges(2,2) * t351 + mrSges(3,1) * t303 - mrSges(3,2) * t304 + t336 * (mrSges(4,2) * t291 + mrSges(5,2) * t284 - mrSges(4,3) * t287 - mrSges(5,3) * t286 - pkin(6) * t364 - qJ(4) * t272 - t340 * t269 + t342 * t270) + t338 * (-mrSges(4,1) * t291 - mrSges(5,1) * t286 + mrSges(5,2) * t285 + mrSges(4,3) * t288 - pkin(3) * t272 - pkin(4) * t349 - pkin(6) * t354 - t342 * t269 - t340 * t270) - pkin(2) * t271 + qJ(3) * t355 + pkin(1) * (t337 * (m(3) * t304 - mrSges(3,1) * t344 + t355) + t339 * (m(3) * t303 - mrSges(3,2) * t344 - t345)) + (Ifges(2,3) + Ifges(3,3) + pkin(1) * (-t337 * mrSges(3,2) + t339 * (mrSges(3,1) - t369)) + (Ifges(4,2) + Ifges(5,3)) * t334 + ((Ifges(4,4) - Ifges(5,5)) * t372 + (Ifges(4,1) + Ifges(5,1)) * t336) * t336) * qJDD(1); m(3) * t335 + t267 * t338 + t268 * t336; t271; (qJDD(1) * mrSges(5,2) + qJD(1) * t316) * t336 + t352; mrSges(6,1) * t278 - mrSges(6,2) * t279 + Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * qJDD(5) + t293 * t310 - t294 * t309;];
tauJ = t1;
