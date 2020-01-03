% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR6
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:35
% EndTime: 2019-12-31 16:52:36
% DurationCPUTime: 0.73s
% Computational Cost: add. (4318->171), mult. (10449->225), div. (0->0), fcn. (7242->8), ass. (0->78)
t307 = qJD(1) ^ 2;
t327 = pkin(2) * t307;
t326 = pkin(5) * qJDD(1);
t303 = sin(qJ(1));
t306 = cos(qJ(1));
t316 = -t306 * g(1) - t303 * g(2);
t288 = -t307 * pkin(1) + qJDD(1) * qJ(2) + t316;
t299 = sin(pkin(7));
t300 = cos(pkin(7));
t321 = qJD(1) * qJD(2);
t319 = -t300 * g(3) - 0.2e1 * t299 * t321;
t267 = (t300 * t327 - t288 - t326) * t299 + t319;
t279 = -t299 * g(3) + (t288 + 0.2e1 * t321) * t300;
t297 = t300 ^ 2;
t268 = -t297 * t327 + t300 * t326 + t279;
t302 = sin(qJ(3));
t305 = cos(qJ(3));
t255 = t305 * t267 - t302 * t268;
t313 = t299 * t305 + t300 * t302;
t312 = -t299 * t302 + t300 * t305;
t286 = t312 * qJD(1);
t322 = t286 * qJD(3);
t277 = t313 * qJDD(1) + t322;
t287 = t313 * qJD(1);
t245 = (-t277 + t322) * pkin(6) + (t286 * t287 + qJDD(3)) * pkin(3) + t255;
t256 = t302 * t267 + t305 * t268;
t276 = -t287 * qJD(3) + t312 * qJDD(1);
t282 = qJD(3) * pkin(3) - t287 * pkin(6);
t285 = t286 ^ 2;
t246 = -t285 * pkin(3) + t276 * pkin(6) - qJD(3) * t282 + t256;
t301 = sin(qJ(4));
t304 = cos(qJ(4));
t243 = t304 * t245 - t301 * t246;
t272 = t304 * t286 - t301 * t287;
t254 = t272 * qJD(4) + t301 * t276 + t304 * t277;
t273 = t301 * t286 + t304 * t287;
t261 = -t272 * mrSges(5,1) + t273 * mrSges(5,2);
t298 = qJD(3) + qJD(4);
t264 = -t298 * mrSges(5,2) + t272 * mrSges(5,3);
t295 = qJDD(3) + qJDD(4);
t240 = m(5) * t243 + t295 * mrSges(5,1) - t254 * mrSges(5,3) - t273 * t261 + t298 * t264;
t244 = t301 * t245 + t304 * t246;
t253 = -t273 * qJD(4) + t304 * t276 - t301 * t277;
t265 = t298 * mrSges(5,1) - t273 * mrSges(5,3);
t241 = m(5) * t244 - t295 * mrSges(5,2) + t253 * mrSges(5,3) + t272 * t261 - t298 * t265;
t233 = t304 * t240 + t301 * t241;
t274 = -t286 * mrSges(4,1) + t287 * mrSges(4,2);
t280 = -qJD(3) * mrSges(4,2) + t286 * mrSges(4,3);
t231 = m(4) * t255 + qJDD(3) * mrSges(4,1) - t277 * mrSges(4,3) + qJD(3) * t280 - t287 * t274 + t233;
t281 = qJD(3) * mrSges(4,1) - t287 * mrSges(4,3);
t317 = -t301 * t240 + t304 * t241;
t232 = m(4) * t256 - qJDD(3) * mrSges(4,2) + t276 * mrSges(4,3) - qJD(3) * t281 + t286 * t274 + t317;
t325 = t305 * t231 + t302 * t232;
t324 = -t299 ^ 2 - t297;
t320 = t303 * g(1) - t306 * g(2);
t318 = -t302 * t231 + t305 * t232;
t315 = qJDD(2) - t320;
t314 = -t300 * mrSges(3,1) + t299 * mrSges(3,2);
t311 = mrSges(3,3) * qJDD(1) + t307 * t314;
t275 = (-pkin(2) * t300 - pkin(1)) * qJDD(1) + (t324 * pkin(5) - qJ(2)) * t307 + t315;
t249 = -t276 * pkin(3) - t285 * pkin(6) + t287 * t282 + t275;
t310 = m(5) * t249 - t253 * mrSges(5,1) + t254 * mrSges(5,2) - t272 * t264 + t273 * t265;
t258 = Ifges(5,4) * t273 + Ifges(5,2) * t272 + Ifges(5,6) * t298;
t259 = Ifges(5,1) * t273 + Ifges(5,4) * t272 + Ifges(5,5) * t298;
t309 = mrSges(5,1) * t243 - mrSges(5,2) * t244 + Ifges(5,5) * t254 + Ifges(5,6) * t253 + Ifges(5,3) * t295 + t273 * t258 - t272 * t259;
t308 = m(4) * t275 - t276 * mrSges(4,1) + t277 * mrSges(4,2) - t286 * t280 + t287 * t281 + t310;
t284 = -qJDD(1) * pkin(1) - t307 * qJ(2) + t315;
t278 = -t299 * t288 + t319;
t271 = Ifges(4,1) * t287 + Ifges(4,4) * t286 + Ifges(4,5) * qJD(3);
t270 = Ifges(4,4) * t287 + Ifges(4,2) * t286 + Ifges(4,6) * qJD(3);
t269 = Ifges(4,5) * t287 + Ifges(4,6) * t286 + Ifges(4,3) * qJD(3);
t257 = Ifges(5,5) * t273 + Ifges(5,6) * t272 + Ifges(5,3) * t298;
t236 = t324 * t307 * mrSges(3,3) + m(3) * t284 + t314 * qJDD(1) + t308;
t235 = mrSges(5,2) * t249 - mrSges(5,3) * t243 + Ifges(5,1) * t254 + Ifges(5,4) * t253 + Ifges(5,5) * t295 + t272 * t257 - t298 * t258;
t234 = -mrSges(5,1) * t249 + mrSges(5,3) * t244 + Ifges(5,4) * t254 + Ifges(5,2) * t253 + Ifges(5,6) * t295 - t273 * t257 + t298 * t259;
t227 = mrSges(4,2) * t275 - mrSges(4,3) * t255 + Ifges(4,1) * t277 + Ifges(4,4) * t276 + Ifges(4,5) * qJDD(3) - pkin(6) * t233 - qJD(3) * t270 - t301 * t234 + t304 * t235 + t286 * t269;
t226 = -mrSges(4,1) * t275 + mrSges(4,3) * t256 + Ifges(4,4) * t277 + Ifges(4,2) * t276 + Ifges(4,6) * qJDD(3) - pkin(3) * t310 + pkin(6) * t317 + qJD(3) * t271 + t304 * t234 + t301 * t235 - t287 * t269;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t320 - mrSges(2,2) * t316 + t299 * (mrSges(3,2) * t284 - mrSges(3,3) * t278 + t305 * t227 - t302 * t226 - pkin(5) * t325 + (Ifges(3,1) * t299 + Ifges(3,4) * t300) * qJDD(1)) + t300 * (-mrSges(3,1) * t284 + mrSges(3,3) * t279 + t302 * t227 + t305 * t226 - pkin(2) * t308 + pkin(5) * t318 + (Ifges(3,4) * t299 + Ifges(3,2) * t300) * qJDD(1)) - pkin(1) * t236 + qJ(2) * ((m(3) * t279 + t311 * t300 + t318) * t300 + (-m(3) * t278 + t311 * t299 - t325) * t299); t236; mrSges(4,1) * t255 - mrSges(4,2) * t256 + Ifges(4,5) * t277 + Ifges(4,6) * t276 + Ifges(4,3) * qJDD(3) + pkin(3) * t233 + t287 * t270 - t286 * t271 + t309; t309;];
tauJ = t1;
