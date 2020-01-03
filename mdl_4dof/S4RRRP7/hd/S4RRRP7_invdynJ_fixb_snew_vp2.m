% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRP7
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:15
% EndTime: 2019-12-31 17:20:16
% DurationCPUTime: 0.73s
% Computational Cost: add. (2399->180), mult. (4668->223), div. (0->0), fcn. (2661->6), ass. (0->73)
t313 = Ifges(4,1) + Ifges(5,1);
t307 = Ifges(4,4) - Ifges(5,5);
t306 = -Ifges(4,5) - Ifges(5,4);
t312 = Ifges(4,2) + Ifges(5,3);
t305 = Ifges(4,6) - Ifges(5,6);
t311 = -Ifges(4,3) - Ifges(5,2);
t281 = sin(qJ(3));
t282 = sin(qJ(2));
t299 = qJD(1) * t282;
t309 = cos(qJ(3));
t267 = -t309 * qJD(2) + t281 * t299;
t284 = cos(qJ(2));
t297 = qJD(1) * qJD(2);
t294 = t284 * t297;
t271 = qJDD(1) * t282 + t294;
t244 = -t267 * qJD(3) + t281 * qJDD(2) + t309 * t271;
t268 = t281 * qJD(2) + t309 * t299;
t248 = mrSges(5,1) * t267 - mrSges(5,3) * t268;
t287 = qJD(1) ^ 2;
t283 = sin(qJ(1));
t285 = cos(qJ(1));
t293 = g(1) * t283 - t285 * g(2);
t262 = -qJDD(1) * pkin(1) - pkin(5) * t287 - t293;
t295 = t282 * t297;
t272 = qJDD(1) * t284 - t295;
t229 = (-t271 - t294) * pkin(6) + (-t272 + t295) * pkin(2) + t262;
t290 = -g(1) * t285 - t283 * g(2);
t263 = -pkin(1) * t287 + qJDD(1) * pkin(5) + t290;
t255 = -g(3) * t282 + t284 * t263;
t270 = (-t284 * pkin(2) - t282 * pkin(6)) * qJD(1);
t286 = qJD(2) ^ 2;
t298 = qJD(1) * t284;
t232 = -pkin(2) * t286 + qJDD(2) * pkin(6) + t270 * t298 + t255;
t226 = t309 * t229 - t281 * t232;
t247 = pkin(3) * t267 - qJ(4) * t268;
t266 = qJDD(3) - t272;
t276 = qJD(3) - t298;
t275 = t276 ^ 2;
t225 = -t266 * pkin(3) - t275 * qJ(4) + t268 * t247 + qJDD(4) - t226;
t253 = -mrSges(5,2) * t267 + mrSges(5,3) * t276;
t291 = -m(5) * t225 + t266 * mrSges(5,1) + t276 * t253;
t221 = mrSges(5,2) * t244 + t248 * t268 - t291;
t227 = t281 * t229 + t309 * t232;
t223 = -pkin(3) * t275 + qJ(4) * t266 + 0.2e1 * qJD(4) * t276 - t247 * t267 + t227;
t243 = qJD(3) * t268 - t309 * qJDD(2) + t271 * t281;
t252 = -mrSges(5,1) * t276 + mrSges(5,2) * t268;
t296 = m(5) * t223 + t266 * mrSges(5,3) + t276 * t252;
t301 = t307 * t267 - t313 * t268 + t306 * t276;
t302 = t312 * t267 - t307 * t268 - t305 * t276;
t310 = -t305 * t243 - t306 * t244 - t311 * t266 - t301 * t267 - t302 * t268 + mrSges(4,1) * t226 - mrSges(5,1) * t225 - mrSges(4,2) * t227 + mrSges(5,3) * t223 - pkin(3) * t221 + qJ(4) * (-mrSges(5,2) * t243 - t248 * t267 + t296);
t308 = -mrSges(4,3) - mrSges(5,2);
t303 = t305 * t267 + t306 * t268 + t311 * t276;
t300 = -mrSges(4,1) * t267 - mrSges(4,2) * t268 - t248;
t254 = -t284 * g(3) - t282 * t263;
t251 = mrSges(4,1) * t276 - mrSges(4,3) * t268;
t217 = m(4) * t227 - mrSges(4,2) * t266 + t308 * t243 - t251 * t276 + t300 * t267 + t296;
t250 = -mrSges(4,2) * t276 - mrSges(4,3) * t267;
t218 = m(4) * t226 + mrSges(4,1) * t266 + t308 * t244 + t250 * t276 + t300 * t268 + t291;
t292 = t309 * t217 - t218 * t281;
t215 = t281 * t217 + t309 * t218;
t231 = -qJDD(2) * pkin(2) - pkin(6) * t286 + t270 * t299 - t254;
t224 = -0.2e1 * qJD(4) * t268 + (t267 * t276 - t244) * qJ(4) + (t268 * t276 + t243) * pkin(3) + t231;
t219 = m(5) * t224 + mrSges(5,1) * t243 - t244 * mrSges(5,3) - t268 * t252 + t253 * t267;
t288 = -m(4) * t231 - t243 * mrSges(4,1) - mrSges(4,2) * t244 - t267 * t250 - t251 * t268 - t219;
t274 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t298;
t273 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t299;
t269 = (-t284 * mrSges(3,1) + t282 * mrSges(3,2)) * qJD(1);
t261 = Ifges(3,5) * qJD(2) + (t282 * Ifges(3,1) + t284 * Ifges(3,4)) * qJD(1);
t260 = Ifges(3,6) * qJD(2) + (t282 * Ifges(3,4) + t284 * Ifges(3,2)) * qJD(1);
t259 = Ifges(3,3) * qJD(2) + (t282 * Ifges(3,5) + t284 * Ifges(3,6)) * qJD(1);
t214 = mrSges(4,2) * t231 + mrSges(5,2) * t225 - mrSges(4,3) * t226 - mrSges(5,3) * t224 - qJ(4) * t219 - t307 * t243 + t313 * t244 - t306 * t266 + t303 * t267 + t302 * t276;
t213 = -mrSges(4,1) * t231 - mrSges(5,1) * t224 + mrSges(5,2) * t223 + mrSges(4,3) * t227 - pkin(3) * t219 - t312 * t243 + t307 * t244 + t305 * t266 + t303 * t268 - t301 * t276;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t293 - mrSges(2,2) * t290 + t282 * (mrSges(3,2) * t262 - mrSges(3,3) * t254 + Ifges(3,1) * t271 + Ifges(3,4) * t272 + Ifges(3,5) * qJDD(2) - pkin(6) * t215 - qJD(2) * t260 - t281 * t213 + t309 * t214 + t259 * t298) + t284 * (-mrSges(3,1) * t262 + mrSges(3,3) * t255 + Ifges(3,4) * t271 + Ifges(3,2) * t272 + Ifges(3,6) * qJDD(2) - pkin(2) * t215 + qJD(2) * t261 - t259 * t299 - t310) + pkin(1) * (-m(3) * t262 + t272 * mrSges(3,1) - t271 * mrSges(3,2) + (-t273 * t282 + t274 * t284) * qJD(1) - t215) + pkin(5) * (t284 * (m(3) * t255 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t272 - qJD(2) * t273 + t269 * t298 + t292) - t282 * (m(3) * t254 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t271 + qJD(2) * t274 - t269 * t299 + t288)); Ifges(3,5) * t271 + Ifges(3,6) * t272 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t254 - mrSges(3,2) * t255 + t281 * t214 + t309 * t213 + pkin(2) * t288 + pkin(6) * t292 + (t282 * t260 - t284 * t261) * qJD(1); t310; t221;];
tauJ = t1;
