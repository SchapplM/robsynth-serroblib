% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:19
% EndTime: 2019-12-31 16:50:21
% DurationCPUTime: 0.98s
% Computational Cost: add. (9375->196), mult. (17500->249), div. (0->0), fcn. (9613->8), ass. (0->81)
t329 = -g(3) + qJDD(2);
t336 = cos(qJ(3));
t353 = t336 * t329;
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t323 = t334 * g(1) - t337 * g(2);
t315 = qJDD(1) * pkin(1) + t323;
t324 = -t337 * g(1) - t334 * g(2);
t339 = qJD(1) ^ 2;
t317 = -t339 * pkin(1) + t324;
t330 = sin(pkin(7));
t331 = cos(pkin(7));
t300 = t330 * t315 + t331 * t317;
t291 = -t339 * pkin(2) + qJDD(1) * pkin(5) + t300;
t333 = sin(qJ(3));
t288 = t336 * t291 + t333 * t329;
t316 = (-mrSges(4,1) * t336 + mrSges(4,2) * t333) * qJD(1);
t349 = qJD(1) * qJD(3);
t347 = t333 * t349;
t320 = t336 * qJDD(1) - t347;
t351 = qJD(1) * t333;
t321 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t351;
t299 = t331 * t315 - t330 * t317;
t290 = -qJDD(1) * pkin(2) - t339 * pkin(5) - t299;
t346 = t336 * t349;
t319 = t333 * qJDD(1) + t346;
t284 = (-t319 - t346) * pkin(6) + (-t320 + t347) * pkin(3) + t290;
t318 = (-pkin(3) * t336 - pkin(6) * t333) * qJD(1);
t338 = qJD(3) ^ 2;
t350 = t336 * qJD(1);
t286 = -t338 * pkin(3) + qJDD(3) * pkin(6) + t318 * t350 + t288;
t332 = sin(qJ(4));
t335 = cos(qJ(4));
t282 = t335 * t284 - t332 * t286;
t313 = t335 * qJD(3) - t332 * t351;
t298 = t313 * qJD(4) + t332 * qJDD(3) + t335 * t319;
t314 = t332 * qJD(3) + t335 * t351;
t301 = -t313 * mrSges(5,1) + t314 * mrSges(5,2);
t325 = qJD(4) - t350;
t302 = -t325 * mrSges(5,2) + t313 * mrSges(5,3);
t312 = qJDD(4) - t320;
t280 = m(5) * t282 + t312 * mrSges(5,1) - t298 * mrSges(5,3) - t314 * t301 + t325 * t302;
t283 = t332 * t284 + t335 * t286;
t297 = -t314 * qJD(4) + t335 * qJDD(3) - t332 * t319;
t303 = t325 * mrSges(5,1) - t314 * mrSges(5,3);
t281 = m(5) * t283 - t312 * mrSges(5,2) + t297 * mrSges(5,3) + t313 * t301 - t325 * t303;
t342 = -t332 * t280 + t335 * t281;
t273 = m(4) * t288 - qJDD(3) * mrSges(4,2) + t320 * mrSges(4,3) - qJD(3) * t321 + t316 * t350 + t342;
t287 = -t333 * t291 + t353;
t322 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t350;
t285 = -qJDD(3) * pkin(3) - t338 * pkin(6) - t353 + (qJD(1) * t318 + t291) * t333;
t341 = -m(5) * t285 + t297 * mrSges(5,1) - t298 * mrSges(5,2) + t313 * t302 - t314 * t303;
t278 = m(4) * t287 + qJDD(3) * mrSges(4,1) - t319 * mrSges(4,3) + qJD(3) * t322 - t316 * t351 + t341;
t343 = t336 * t273 - t333 * t278;
t267 = m(3) * t300 - t339 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t343;
t274 = t335 * t280 + t332 * t281;
t340 = -m(4) * t290 + t320 * mrSges(4,1) - t319 * mrSges(4,2) - t321 * t351 + t322 * t350 - t274;
t270 = m(3) * t299 + qJDD(1) * mrSges(3,1) - t339 * mrSges(3,2) + t340;
t262 = t330 * t267 + t331 * t270;
t260 = m(2) * t323 + qJDD(1) * mrSges(2,1) - t339 * mrSges(2,2) + t262;
t344 = t331 * t267 - t330 * t270;
t261 = m(2) * t324 - t339 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t344;
t352 = t337 * t260 + t334 * t261;
t268 = t333 * t273 + t336 * t278;
t348 = m(3) * t329 + t268;
t345 = -t334 * t260 + t337 * t261;
t309 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t333 + Ifges(4,4) * t336) * qJD(1);
t308 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t333 + Ifges(4,2) * t336) * qJD(1);
t307 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t333 + Ifges(4,6) * t336) * qJD(1);
t294 = Ifges(5,1) * t314 + Ifges(5,4) * t313 + Ifges(5,5) * t325;
t293 = Ifges(5,4) * t314 + Ifges(5,2) * t313 + Ifges(5,6) * t325;
t292 = Ifges(5,5) * t314 + Ifges(5,6) * t313 + Ifges(5,3) * t325;
t276 = mrSges(5,2) * t285 - mrSges(5,3) * t282 + Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t312 + t313 * t292 - t325 * t293;
t275 = -mrSges(5,1) * t285 + mrSges(5,3) * t283 + Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t312 - t314 * t292 + t325 * t294;
t264 = -mrSges(4,1) * t290 - mrSges(5,1) * t282 + mrSges(5,2) * t283 + mrSges(4,3) * t288 + Ifges(4,4) * t319 - Ifges(5,5) * t298 + Ifges(4,2) * t320 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t297 - Ifges(5,3) * t312 - pkin(3) * t274 + qJD(3) * t309 - t314 * t293 + t313 * t294 - t307 * t351;
t263 = mrSges(4,2) * t290 - mrSges(4,3) * t287 + Ifges(4,1) * t319 + Ifges(4,4) * t320 + Ifges(4,5) * qJDD(3) - pkin(6) * t274 - qJD(3) * t308 - t332 * t275 + t335 * t276 + t307 * t350;
t256 = Ifges(3,6) * qJDD(1) + t339 * Ifges(3,5) - mrSges(3,1) * t329 + mrSges(3,3) * t300 - Ifges(4,5) * t319 - Ifges(4,6) * t320 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t287 + mrSges(4,2) * t288 - t332 * t276 - t335 * t275 - pkin(3) * t341 - pkin(6) * t342 - pkin(2) * t268 + (-t333 * t308 + t336 * t309) * qJD(1);
t255 = mrSges(3,2) * t329 - mrSges(3,3) * t299 + Ifges(3,5) * qJDD(1) - t339 * Ifges(3,6) - pkin(5) * t268 + t336 * t263 - t333 * t264;
t254 = -mrSges(2,2) * g(3) - mrSges(2,3) * t323 + Ifges(2,5) * qJDD(1) - t339 * Ifges(2,6) - qJ(2) * t262 + t331 * t255 - t330 * t256;
t253 = mrSges(2,1) * g(3) + mrSges(2,3) * t324 + t339 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t348 + qJ(2) * t344 + t330 * t255 + t331 * t256;
t1 = [-m(1) * g(1) + t345; -m(1) * g(2) + t352; (-m(1) - m(2)) * g(3) + t348; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t352 - t334 * t253 + t337 * t254; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t345 + t337 * t253 + t334 * t254; pkin(1) * t262 + mrSges(2,1) * t323 - mrSges(2,2) * t324 + t333 * t263 + t336 * t264 + pkin(2) * t340 + pkin(5) * t343 + mrSges(3,1) * t299 - mrSges(3,2) * t300 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
