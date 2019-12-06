% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:25
% EndTime: 2019-12-05 17:51:26
% DurationCPUTime: 0.66s
% Computational Cost: add. (4202->132), mult. (5880->179), div. (0->0), fcn. (3273->10), ass. (0->74)
t332 = 2 * qJD(4);
t305 = sin(qJ(1));
t308 = cos(qJ(1));
t322 = t308 * g(2) + t305 * g(3);
t282 = qJDD(1) * pkin(1) + t322;
t309 = qJD(1) ^ 2;
t317 = t305 * g(2) - t308 * g(3);
t283 = -t309 * pkin(1) + t317;
t300 = sin(pkin(8));
t302 = cos(pkin(8));
t269 = t302 * t282 - t300 * t283;
t264 = qJDD(1) * pkin(2) + t269;
t270 = t300 * t282 + t302 * t283;
t265 = -t309 * pkin(2) + t270;
t304 = sin(qJ(3));
t307 = cos(qJ(3));
t260 = t304 * t264 + t307 * t265;
t297 = (qJD(1) + qJD(3));
t295 = t297 ^ 2;
t296 = qJDD(1) + qJDD(3);
t257 = -t295 * pkin(3) + t296 * qJ(4) + t260;
t331 = (t297 * t332) + t257;
t299 = sin(pkin(9));
t330 = mrSges(5,2) * t299;
t328 = mrSges(5,3) * t296;
t327 = t297 * t299;
t301 = cos(pkin(9));
t326 = t301 * t296;
t325 = t301 * t297;
t298 = -g(1) + qJDD(2);
t324 = t301 * t298;
t253 = t299 * t298 + t331 * t301;
t315 = -pkin(4) * t301 - pkin(7) * t299;
t278 = t315 * t297;
t251 = t278 * t325 + t253;
t259 = t307 * t264 - t304 * t265;
t313 = -t295 * qJ(4) + qJDD(4) - t259;
t254 = (-pkin(3) + t315) * t296 + t313;
t303 = sin(qJ(5));
t306 = cos(qJ(5));
t248 = -t303 * t251 + t306 * t254;
t285 = qJD(5) - t325;
t319 = t303 * t327;
t271 = -t285 * mrSges(6,2) - mrSges(6,3) * t319;
t273 = (t303 * mrSges(6,1) + t306 * mrSges(6,2)) * t327;
t320 = qJD(5) * t297;
t275 = (t296 * t306 - t303 * t320) * t299;
t284 = qJDD(5) - t326;
t318 = t306 * t327;
t246 = m(6) * t248 + t284 * mrSges(6,1) - t275 * mrSges(6,3) + t285 * t271 - t273 * t318;
t249 = t306 * t251 + t303 * t254;
t272 = t285 * mrSges(6,1) - mrSges(6,3) * t318;
t274 = (-t296 * t303 - t306 * t320) * t299;
t247 = m(6) * t249 - t284 * mrSges(6,2) + t274 * mrSges(6,3) - t285 * t272 - t273 * t319;
t276 = (-mrSges(5,1) * t301 + t330) * t297;
t241 = m(5) * t253 - t303 * t246 + t306 * t247 + (t276 * t297 + t328) * t301;
t250 = -t324 + (t257 + (t332 + t278) * t297) * t299;
t252 = -t331 * t299 + t324;
t245 = m(5) * t252 - m(6) * t250 + t274 * mrSges(6,1) - t275 * mrSges(6,2) + (-t328 + (-t271 * t303 - t272 * t306 - t276) * t297) * t299;
t316 = t301 * t241 - t299 * t245;
t236 = m(4) * t260 - t295 * mrSges(4,1) - t296 * mrSges(4,2) + t316;
t244 = t306 * t246 + t303 * t247;
t256 = -t296 * pkin(3) + t313;
t311 = -m(5) * t256 + mrSges(5,1) * t326 - t244 + (t299 ^ 2 + t301 ^ 2) * mrSges(5,3) * t295;
t239 = m(4) * t259 - t295 * mrSges(4,2) + (mrSges(4,1) - t330) * t296 + t311;
t323 = t304 * t236 + t307 * t239;
t267 = Ifges(6,6) * t285 + (t306 * Ifges(6,4) - t303 * Ifges(6,2)) * t327;
t268 = Ifges(6,5) * t285 + (t306 * Ifges(6,1) - t303 * Ifges(6,4)) * t327;
t314 = t306 * t267 + t303 * t268;
t243 = t296 * t330 - t311;
t277 = (Ifges(5,5) * t299 + Ifges(5,6) * t301) * t297;
t310 = mrSges(6,1) * t248 - mrSges(6,2) * t249 + Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t284;
t312 = -mrSges(4,2) * t260 + t299 * (t277 * t325 + mrSges(5,2) * t256 - mrSges(5,3) * t252 + t306 * (mrSges(6,2) * t250 - mrSges(6,3) * t248 + Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t284 - t285 * t267) - t303 * (-mrSges(6,1) * t250 + mrSges(6,3) * t249 + Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t284 + t285 * t268) - pkin(7) * t244 + (Ifges(5,1) * t299 + Ifges(5,4) * t301) * t296) + t301 * (Ifges(5,2) * t326 - mrSges(5,1) * t256 + mrSges(5,3) * t253 - pkin(4) * t244 + (Ifges(5,4) * t296 + (-t277 - t314) * t297) * t299 - t310) + qJ(4) * t316 - pkin(3) * t243 + mrSges(4,1) * t259 + Ifges(4,3) * t296;
t1 = [mrSges(2,1) * t322 - mrSges(2,2) * t317 + pkin(1) * (t300 * (m(3) * t270 - t309 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t307 * t236 - t304 * t239) + t302 * (m(3) * t269 + qJDD(1) * mrSges(3,1) - t309 * mrSges(3,2) + t323)) + pkin(2) * t323 + mrSges(3,1) * t269 - mrSges(3,2) * t270 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t312; t299 * t241 + t301 * t245 + (m(3) + m(4)) * t298; t312; t243; t314 * t327 + t310;];
tauJ = t1;
