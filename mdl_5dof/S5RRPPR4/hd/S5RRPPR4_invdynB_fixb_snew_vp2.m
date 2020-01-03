% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:45
% EndTime: 2019-12-31 19:27:46
% DurationCPUTime: 1.04s
% Computational Cost: add. (14278->182), mult. (16473->217), div. (0->0), fcn. (6962->8), ass. (0->77)
t335 = qJD(1) + qJD(2);
t366 = t335 ^ 2;
t365 = -m(3) - m(4);
t364 = -pkin(2) - pkin(3);
t363 = mrSges(3,1) + mrSges(4,1);
t362 = Ifges(4,4) + Ifges(3,5);
t361 = Ifges(3,6) - Ifges(4,6);
t340 = sin(qJ(5));
t360 = t335 * t340;
t343 = cos(qJ(5));
t359 = t335 * t343;
t342 = sin(qJ(1));
t345 = cos(qJ(1));
t324 = t342 * g(1) - t345 * g(2);
t320 = qJDD(1) * pkin(1) + t324;
t325 = -t345 * g(1) - t342 * g(2);
t346 = qJD(1) ^ 2;
t321 = -t346 * pkin(1) + t325;
t341 = sin(qJ(2));
t344 = cos(qJ(2));
t308 = t341 * t320 + t344 * t321;
t334 = qJDD(1) + qJDD(2);
t352 = t334 * qJ(3) + 0.2e1 * qJD(3) * t335 + t308;
t305 = -pkin(2) * t366 + t352;
t301 = t364 * t366 + t352;
t307 = t344 * t320 - t341 * t321;
t349 = -qJ(3) * t366 + qJDD(3) - t307;
t304 = t364 * t334 + t349;
t338 = sin(pkin(8));
t339 = cos(pkin(8));
t299 = t339 * t301 + t338 * t304;
t297 = -pkin(4) * t366 - t334 * pkin(7) + t299;
t337 = g(3) + qJDD(4);
t294 = -t340 * t297 + t343 * t337;
t314 = (mrSges(6,1) * t343 - mrSges(6,2) * t340) * t335;
t357 = qJD(5) * t335;
t315 = -t340 * t334 - t343 * t357;
t323 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t359;
t292 = m(6) * t294 + qJDD(5) * mrSges(6,1) - t315 * mrSges(6,3) + qJD(5) * t323 + t314 * t360;
t295 = t343 * t297 + t340 * t337;
t316 = -t343 * t334 + t340 * t357;
t322 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t360;
t293 = m(6) * t295 - qJDD(5) * mrSges(6,2) + t316 * mrSges(6,3) - qJD(5) * t322 - t314 * t359;
t353 = -t340 * t292 + t343 * t293;
t285 = m(5) * t299 - mrSges(5,1) * t366 + t334 * mrSges(5,2) + t353;
t298 = -t338 * t301 + t339 * t304;
t296 = t334 * pkin(4) - pkin(7) * t366 - t298;
t347 = -m(6) * t296 + t316 * mrSges(6,1) - t315 * mrSges(6,2) + t322 * t360 - t323 * t359;
t290 = m(5) * t298 - t334 * mrSges(5,1) - mrSges(5,2) * t366 + t347;
t354 = t339 * t285 - t338 * t290;
t351 = m(4) * t305 + t334 * mrSges(4,3) + t354;
t279 = m(3) * t308 - t334 * mrSges(3,2) - t363 * t366 + t351;
t283 = t338 * t285 + t339 * t290;
t306 = -t334 * pkin(2) + t349;
t348 = -m(4) * t306 + t334 * mrSges(4,1) + mrSges(4,3) * t366 - t283;
t281 = m(3) * t307 + t334 * mrSges(3,1) - mrSges(3,2) * t366 + t348;
t275 = t341 * t279 + t344 * t281;
t273 = m(2) * t324 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t275;
t355 = t344 * t279 - t341 * t281;
t274 = m(2) * t325 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t355;
t358 = t345 * t273 + t342 * t274;
t356 = -t342 * t273 + t345 * t274;
t287 = t343 * t292 + t340 * t293;
t350 = m(5) * t337 + t287;
t311 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t340 - Ifges(6,4) * t343) * t335;
t310 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t340 - Ifges(6,2) * t343) * t335;
t309 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t340 - Ifges(6,6) * t343) * t335;
t289 = mrSges(6,2) * t296 - mrSges(6,3) * t294 + Ifges(6,1) * t315 + Ifges(6,4) * t316 + Ifges(6,5) * qJDD(5) - qJD(5) * t310 - t309 * t359;
t288 = -mrSges(6,1) * t296 + mrSges(6,3) * t295 + Ifges(6,4) * t315 + Ifges(6,2) * t316 + Ifges(6,6) * qJDD(5) + qJD(5) * t311 + t309 * t360;
t286 = -m(4) * g(3) - t350;
t282 = -mrSges(5,1) * t337 - mrSges(6,1) * t294 + mrSges(6,2) * t295 + mrSges(5,3) * t299 + t366 * Ifges(5,5) - Ifges(6,5) * t315 - Ifges(5,6) * t334 - Ifges(6,6) * t316 - Ifges(6,3) * qJDD(5) - pkin(4) * t287 + (t310 * t340 - t311 * t343) * t335;
t276 = mrSges(5,2) * t337 - mrSges(5,3) * t298 - Ifges(5,5) * t334 - Ifges(5,6) * t366 - pkin(7) * t287 - t340 * t288 + t343 * t289;
t269 = mrSges(4,2) * t306 - mrSges(3,3) * t307 - qJ(3) * t286 - qJ(4) * t283 + t339 * t276 - t338 * t282 + t362 * t334 - t361 * t366 + (-mrSges(3,2) + mrSges(4,3)) * g(3);
t268 = mrSges(4,2) * t305 + mrSges(3,3) * t308 - pkin(2) * t286 + pkin(3) * t350 + t363 * g(3) - qJ(4) * t354 - t338 * t276 - t339 * t282 + t361 * t334 + t362 * t366;
t267 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(6) * t275 - t341 * t268 + t344 * t269;
t266 = Ifges(2,6) * qJDD(1) + t346 * Ifges(2,5) + mrSges(2,3) * t325 + t341 * t269 + t344 * t268 + pkin(1) * t350 + pkin(6) * t355 + (-pkin(1) * t365 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t356; -m(1) * g(2) + t358; (-m(1) - m(2) + t365) * g(3) - t350; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t358 - t342 * t266 + t345 * t267; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t356 + t345 * t266 + t342 * t267; pkin(1) * t275 - mrSges(2,2) * t325 + mrSges(2,1) * t324 + qJ(3) * (-mrSges(4,1) * t366 + t351) + pkin(2) * t348 + mrSges(3,1) * t307 - mrSges(3,2) * t308 - pkin(3) * t283 - mrSges(4,1) * t306 + mrSges(4,3) * t305 - t343 * t288 - pkin(4) * t347 - pkin(7) * t353 - mrSges(5,1) * t298 + mrSges(5,2) * t299 - t340 * t289 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (Ifges(5,3) + Ifges(3,3) + Ifges(4,2)) * t334;];
tauB = t1;
