% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:44
% EndTime: 2019-12-05 15:12:46
% DurationCPUTime: 1.41s
% Computational Cost: add. (18590->163), mult. (24222->212), div. (0->0), fcn. (16760->10), ass. (0->75)
t317 = qJD(3) + qJD(4);
t323 = sin(qJ(5));
t341 = t317 * t323;
t326 = cos(qJ(5));
t340 = t317 * t326;
t320 = sin(pkin(8));
t322 = cos(pkin(8));
t314 = -t322 * g(1) - t320 * g(2);
t318 = -g(3) + qJDD(1);
t319 = sin(pkin(9));
t321 = cos(pkin(9));
t299 = -t319 * t314 + t321 * t318;
t300 = t321 * t314 + t319 * t318;
t325 = sin(qJ(3));
t328 = cos(qJ(3));
t294 = t328 * t299 - t325 * t300;
t292 = qJDD(3) * pkin(3) + t294;
t295 = t325 * t299 + t328 * t300;
t329 = qJD(3) ^ 2;
t293 = -t329 * pkin(3) + t295;
t324 = sin(qJ(4));
t327 = cos(qJ(4));
t289 = t324 * t292 + t327 * t293;
t315 = t317 ^ 2;
t316 = qJDD(3) + qJDD(4);
t287 = -t315 * pkin(4) + t316 * pkin(7) + t289;
t313 = t320 * g(1) - t322 * g(2);
t312 = qJDD(2) - t313;
t284 = -t323 * t287 + t326 * t312;
t306 = (-mrSges(6,1) * t326 + mrSges(6,2) * t323) * t317;
t338 = qJD(5) * t317;
t307 = t323 * t316 + t326 * t338;
t311 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t340;
t282 = m(6) * t284 + qJDD(5) * mrSges(6,1) - t307 * mrSges(6,3) + qJD(5) * t311 - t306 * t341;
t285 = t326 * t287 + t323 * t312;
t308 = t326 * t316 - t323 * t338;
t310 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t341;
t283 = m(6) * t285 - qJDD(5) * mrSges(6,2) + t308 * mrSges(6,3) - qJD(5) * t310 + t306 * t340;
t332 = -t323 * t282 + t326 * t283;
t271 = m(5) * t289 - t315 * mrSges(5,1) - t316 * mrSges(5,2) + t332;
t288 = t327 * t292 - t324 * t293;
t286 = -t316 * pkin(4) - t315 * pkin(7) - t288;
t330 = -m(6) * t286 + t308 * mrSges(6,1) - t307 * mrSges(6,2) - t310 * t341 + t311 * t340;
t278 = m(5) * t288 + t316 * mrSges(5,1) - t315 * mrSges(5,2) + t330;
t268 = t324 * t271 + t327 * t278;
t266 = m(4) * t294 + qJDD(3) * mrSges(4,1) - t329 * mrSges(4,2) + t268;
t333 = t327 * t271 - t324 * t278;
t267 = m(4) * t295 - t329 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t333;
t260 = t328 * t266 + t325 * t267;
t258 = m(3) * t299 + t260;
t334 = -t325 * t266 + t328 * t267;
t259 = m(3) * t300 + t334;
t335 = -t319 * t258 + t321 * t259;
t251 = m(2) * t314 + t335;
t274 = t326 * t282 + t323 * t283;
t337 = m(5) * t312 + t274;
t331 = (-m(3) - m(4)) * t312 - t337;
t273 = m(2) * t313 + t331;
t339 = t320 * t251 + t322 * t273;
t252 = t321 * t258 + t319 * t259;
t336 = t322 * t251 - t320 * t273;
t303 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t323 + Ifges(6,4) * t326) * t317;
t302 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t323 + Ifges(6,2) * t326) * t317;
t301 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t323 + Ifges(6,6) * t326) * t317;
t276 = mrSges(6,2) * t286 - mrSges(6,3) * t284 + Ifges(6,1) * t307 + Ifges(6,4) * t308 + Ifges(6,5) * qJDD(5) - qJD(5) * t302 + t301 * t340;
t275 = -mrSges(6,1) * t286 + mrSges(6,3) * t285 + Ifges(6,4) * t307 + Ifges(6,2) * t308 + Ifges(6,6) * qJDD(5) + qJD(5) * t303 - t301 * t341;
t262 = -mrSges(5,1) * t312 - mrSges(6,1) * t284 + mrSges(6,2) * t285 + mrSges(5,3) * t289 + t315 * Ifges(5,5) - Ifges(6,5) * t307 + Ifges(5,6) * t316 - Ifges(6,6) * t308 - Ifges(6,3) * qJDD(5) - pkin(4) * t274 + (-t302 * t323 + t303 * t326) * t317;
t261 = mrSges(5,2) * t312 - mrSges(5,3) * t288 + Ifges(5,5) * t316 - t315 * Ifges(5,6) - pkin(7) * t274 - t323 * t275 + t326 * t276;
t254 = mrSges(4,2) * t312 - mrSges(4,3) * t294 + Ifges(4,5) * qJDD(3) - t329 * Ifges(4,6) - pkin(6) * t268 + t327 * t261 - t324 * t262;
t253 = -mrSges(4,1) * t312 + mrSges(4,3) * t295 + t329 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t337 + pkin(6) * t333 + t324 * t261 + t327 * t262;
t248 = -pkin(1) * t252 - pkin(2) * t260 - mrSges(3,1) * t299 + mrSges(3,2) * t300 - pkin(3) * t268 - mrSges(4,1) * t294 + mrSges(4,2) * t295 - pkin(7) * t332 - mrSges(5,1) * t288 + mrSges(5,2) * t289 - t323 * t276 - t326 * t275 - pkin(4) * t330 + mrSges(2,3) * t314 - Ifges(4,3) * qJDD(3) - Ifges(5,3) * t316 - mrSges(2,1) * t318;
t247 = mrSges(3,2) * t312 - mrSges(3,3) * t299 - pkin(5) * t260 - t325 * t253 + t328 * t254;
t246 = -mrSges(3,1) * t312 + mrSges(3,3) * t300 + t325 * t254 + t328 * t253 - pkin(2) * (m(4) * t312 + t337) + pkin(5) * t334;
t245 = mrSges(2,2) * t318 - mrSges(2,3) * t313 - qJ(2) * t252 - t319 * t246 + t321 * t247;
t1 = [-m(1) * g(1) + t336; -m(1) * g(2) + t339; -m(1) * g(3) + m(2) * t318 + t252; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t339 + t322 * t245 - t320 * t248; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t336 + t320 * t245 + t322 * t248; -mrSges(1,1) * g(2) + mrSges(2,1) * t313 + mrSges(1,2) * g(1) - mrSges(2,2) * t314 + pkin(1) * t331 + qJ(2) * t335 + t321 * t246 + t319 * t247;];
tauB = t1;
