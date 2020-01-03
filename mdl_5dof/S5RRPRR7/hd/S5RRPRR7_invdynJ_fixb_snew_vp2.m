% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:35
% EndTime: 2019-12-31 20:15:36
% DurationCPUTime: 0.59s
% Computational Cost: add. (5091->164), mult. (6175->211), div. (0->0), fcn. (3351->8), ass. (0->71)
t319 = sin(qJ(1));
t323 = cos(qJ(1));
t333 = t319 * g(1) - t323 * g(2);
t295 = qJDD(1) * pkin(1) + t333;
t330 = -t323 * g(1) - t319 * g(2);
t296 = -qJD(1) ^ 2 * pkin(1) + t330;
t318 = sin(qJ(2));
t322 = cos(qJ(2));
t278 = t318 * t295 + t322 * t296;
t312 = qJDD(1) + qJDD(2);
t314 = qJD(1) + qJD(2);
t331 = t312 * qJ(3) + 0.2e1 * qJD(3) * t314 + t278;
t340 = -pkin(2) - pkin(7);
t317 = sin(qJ(4));
t339 = t314 * t317;
t321 = cos(qJ(4));
t338 = t314 * t321;
t277 = t322 * t295 - t318 * t296;
t310 = t314 ^ 2;
t329 = -t310 * qJ(3) + qJDD(3) - t277;
t267 = t340 * t312 + t329;
t261 = t317 * g(3) + t321 * t267;
t336 = qJD(4) * t314;
t334 = t317 * t336;
t291 = t321 * t312 - t334;
t252 = (-t291 - t334) * pkin(8) + (-t310 * t317 * t321 + qJDD(4)) * pkin(4) + t261;
t262 = -t321 * g(3) + t317 * t267;
t290 = -t317 * t312 - t321 * t336;
t299 = qJD(4) * pkin(4) - pkin(8) * t338;
t315 = t317 ^ 2;
t254 = -t315 * t310 * pkin(4) + t290 * pkin(8) - qJD(4) * t299 + t262;
t316 = sin(qJ(5));
t320 = cos(qJ(5));
t249 = t320 * t252 - t316 * t254;
t284 = (-t316 * t321 - t317 * t320) * t314;
t260 = t284 * qJD(5) + t316 * t290 + t320 * t291;
t285 = (-t316 * t317 + t320 * t321) * t314;
t275 = -t284 * mrSges(6,1) + t285 * mrSges(6,2);
t313 = qJD(4) + qJD(5);
t279 = -t313 * mrSges(6,2) + t284 * mrSges(6,3);
t311 = qJDD(4) + qJDD(5);
t246 = m(6) * t249 + t311 * mrSges(6,1) - t260 * mrSges(6,3) - t285 * t275 + t313 * t279;
t250 = t316 * t252 + t320 * t254;
t259 = -t285 * qJD(5) + t320 * t290 - t316 * t291;
t280 = t313 * mrSges(6,1) - t285 * mrSges(6,3);
t247 = m(6) * t250 - t311 * mrSges(6,2) + t259 * mrSges(6,3) + t284 * t275 - t313 * t280;
t239 = t320 * t246 + t316 * t247;
t332 = -t316 * t246 + t320 * t247;
t289 = (mrSges(5,1) * t317 + mrSges(5,2) * t321) * t314;
t297 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t339;
t298 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t338;
t236 = t321 * (m(5) * t261 + qJDD(4) * mrSges(5,1) - t291 * mrSges(5,3) + qJD(4) * t297 - t289 * t338 + t239) + t317 * (m(5) * t262 - qJDD(4) * mrSges(5,2) + t290 * mrSges(5,3) - qJD(4) * t298 - t289 * t339 + t332);
t253 = t299 * t338 - t290 * pkin(4) + (-pkin(8) * t315 + t340) * t310 + t331;
t328 = m(6) * t253 - t259 * mrSges(6,1) + t260 * mrSges(6,2) - t284 * t279 + t285 * t280;
t274 = -t312 * pkin(2) + t329;
t327 = -m(4) * t274 + t310 * mrSges(4,3) - t236;
t269 = Ifges(6,4) * t285 + Ifges(6,2) * t284 + Ifges(6,6) * t313;
t270 = Ifges(6,1) * t285 + Ifges(6,4) * t284 + Ifges(6,5) * t313;
t326 = mrSges(6,1) * t249 - mrSges(6,2) * t250 + Ifges(6,5) * t260 + Ifges(6,6) * t259 + Ifges(6,3) * t311 + t285 * t269 - t284 * t270;
t235 = t312 * mrSges(4,2) - t327;
t268 = Ifges(6,5) * t285 + Ifges(6,6) * t284 + Ifges(6,3) * t313;
t240 = -mrSges(6,1) * t253 + mrSges(6,3) * t250 + Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t311 - t285 * t268 + t313 * t270;
t241 = mrSges(6,2) * t253 - mrSges(6,3) * t249 + Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t311 + t284 * t268 - t313 * t269;
t264 = t340 * t310 + t331;
t271 = t310 * pkin(2) - t331;
t281 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t321 - Ifges(5,6) * t317) * t314;
t282 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t321 - Ifges(5,2) * t317) * t314;
t283 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t321 - Ifges(5,4) * t317) * t314;
t324 = -m(4) * t271 + m(5) * t264 - t290 * mrSges(5,1) + t310 * mrSges(4,2) + t291 * mrSges(5,2) + t312 * mrSges(4,3) + t297 * t339 + t298 * t338 + t328;
t325 = -mrSges(3,2) * t278 - mrSges(4,3) * t271 - pkin(7) * t236 - t317 * (-mrSges(5,1) * t264 + mrSges(5,3) * t262 + Ifges(5,4) * t291 + Ifges(5,2) * t290 + Ifges(5,6) * qJDD(4) - pkin(4) * t328 + pkin(8) * t332 + qJD(4) * t283 + t320 * t240 + t316 * t241 - t281 * t338) + t321 * (mrSges(5,2) * t264 - mrSges(5,3) * t261 + Ifges(5,1) * t291 + Ifges(5,4) * t290 + Ifges(5,5) * qJDD(4) - pkin(8) * t239 - qJD(4) * t282 - t316 * t240 + t320 * t241 - t281 * t339) - pkin(2) * t235 + qJ(3) * t324 + mrSges(4,2) * t274 + mrSges(3,1) * t277 + (Ifges(3,3) + Ifges(4,1)) * t312;
t1 = [t325 + Ifges(2,3) * qJDD(1) + pkin(1) * (t318 * (m(3) * t278 - t310 * mrSges(3,1) - t312 * mrSges(3,2) + t324) + t322 * (m(3) * t277 - t310 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * t312 + t327)) + mrSges(2,1) * t333 - mrSges(2,2) * t330; t325; t235; mrSges(5,1) * t261 - mrSges(5,2) * t262 + Ifges(5,5) * t291 + Ifges(5,6) * t290 + Ifges(5,3) * qJDD(4) + pkin(4) * t239 + (t282 * t321 + t283 * t317) * t314 + t326; t326;];
tauJ = t1;
