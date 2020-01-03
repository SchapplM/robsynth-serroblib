% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR7
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:04
% EndTime: 2019-12-31 17:06:05
% DurationCPUTime: 0.80s
% Computational Cost: add. (4418->202), mult. (9959->266), div. (0->0), fcn. (6271->8), ass. (0->82)
t291 = sin(pkin(7));
t292 = cos(pkin(7));
t294 = sin(qJ(2));
t297 = cos(qJ(2));
t275 = (t294 * t291 - t297 * t292) * qJD(1);
t315 = 2 * qJD(3);
t300 = qJD(1) ^ 2;
t314 = pkin(2) * t300;
t295 = sin(qJ(1));
t298 = cos(qJ(1));
t306 = -t298 * g(1) - t295 * g(2);
t281 = -t300 * pkin(1) + qJDD(1) * pkin(5) + t306;
t313 = t294 * t281;
t310 = qJD(1) * qJD(2);
t284 = t294 * qJDD(1) + t297 * t310;
t252 = qJDD(2) * pkin(2) - t284 * qJ(3) - t313 + (qJ(3) * t310 + t294 * t314 - g(3)) * t297;
t269 = -t294 * g(3) + t297 * t281;
t285 = t297 * qJDD(1) - t294 * t310;
t312 = qJD(1) * t294;
t286 = qJD(2) * pkin(2) - qJ(3) * t312;
t290 = t297 ^ 2;
t253 = t285 * qJ(3) - qJD(2) * t286 - t290 * t314 + t269;
t242 = t291 * t252 + t292 * t253 - t275 * t315;
t276 = (t297 * t291 + t294 * t292) * qJD(1);
t261 = t275 * mrSges(4,1) + t276 * mrSges(4,2);
t264 = -t291 * t284 + t292 * t285;
t271 = qJD(2) * mrSges(4,1) - t276 * mrSges(4,3);
t262 = t275 * pkin(3) - t276 * pkin(6);
t299 = qJD(2) ^ 2;
t239 = -t299 * pkin(3) + qJDD(2) * pkin(6) - t275 * t262 + t242;
t309 = t295 * g(1) - t298 * g(2);
t303 = -qJDD(1) * pkin(1) - t309;
t255 = -t285 * pkin(2) + qJDD(3) + t286 * t312 + (-qJ(3) * t290 - pkin(5)) * t300 + t303;
t265 = t292 * t284 + t291 * t285;
t240 = (qJD(2) * t275 - t265) * pkin(6) + (qJD(2) * t276 - t264) * pkin(3) + t255;
t293 = sin(qJ(4));
t296 = cos(qJ(4));
t236 = -t293 * t239 + t296 * t240;
t266 = t296 * qJD(2) - t293 * t276;
t249 = t266 * qJD(4) + t293 * qJDD(2) + t296 * t265;
t267 = t293 * qJD(2) + t296 * t276;
t254 = -t266 * mrSges(5,1) + t267 * mrSges(5,2);
t274 = qJD(4) + t275;
t256 = -t274 * mrSges(5,2) + t266 * mrSges(5,3);
t263 = qJDD(4) - t264;
t234 = m(5) * t236 + t263 * mrSges(5,1) - t249 * mrSges(5,3) - t267 * t254 + t274 * t256;
t237 = t296 * t239 + t293 * t240;
t248 = -t267 * qJD(4) + t296 * qJDD(2) - t293 * t265;
t257 = t274 * mrSges(5,1) - t267 * mrSges(5,3);
t235 = m(5) * t237 - t263 * mrSges(5,2) + t248 * mrSges(5,3) + t266 * t254 - t274 * t257;
t307 = -t293 * t234 + t296 * t235;
t224 = m(4) * t242 - qJDD(2) * mrSges(4,2) + t264 * mrSges(4,3) - qJD(2) * t271 - t275 * t261 + t307;
t305 = -t292 * t252 + t291 * t253;
t241 = -0.2e1 * qJD(3) * t276 - t305;
t270 = -qJD(2) * mrSges(4,2) - t275 * mrSges(4,3);
t238 = -qJDD(2) * pkin(3) - t299 * pkin(6) + (t315 + t262) * t276 + t305;
t302 = -m(5) * t238 + t248 * mrSges(5,1) - t249 * mrSges(5,2) + t266 * t256 - t267 * t257;
t230 = m(4) * t241 + qJDD(2) * mrSges(4,1) - t265 * mrSges(4,3) + qJD(2) * t270 - t276 * t261 + t302;
t221 = t291 * t224 + t292 * t230;
t226 = t296 * t234 + t293 * t235;
t311 = qJD(1) * t297;
t308 = t292 * t224 - t291 * t230;
t225 = m(4) * t255 - t264 * mrSges(4,1) + t265 * mrSges(4,2) + t275 * t270 + t276 * t271 + t226;
t244 = Ifges(5,4) * t267 + Ifges(5,2) * t266 + Ifges(5,6) * t274;
t245 = Ifges(5,1) * t267 + Ifges(5,4) * t266 + Ifges(5,5) * t274;
t301 = mrSges(5,1) * t236 - mrSges(5,2) * t237 + Ifges(5,5) * t249 + Ifges(5,6) * t248 + Ifges(5,3) * t263 + t267 * t244 - t266 * t245;
t288 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t311;
t287 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t312;
t283 = (-t297 * mrSges(3,1) + t294 * mrSges(3,2)) * qJD(1);
t280 = -t300 * pkin(5) + t303;
t279 = Ifges(3,5) * qJD(2) + (t294 * Ifges(3,1) + t297 * Ifges(3,4)) * qJD(1);
t278 = Ifges(3,6) * qJD(2) + (t294 * Ifges(3,4) + t297 * Ifges(3,2)) * qJD(1);
t268 = -t297 * g(3) - t313;
t260 = Ifges(4,1) * t276 - Ifges(4,4) * t275 + Ifges(4,5) * qJD(2);
t259 = Ifges(4,4) * t276 - Ifges(4,2) * t275 + Ifges(4,6) * qJD(2);
t258 = Ifges(4,5) * t276 - Ifges(4,6) * t275 + Ifges(4,3) * qJD(2);
t243 = Ifges(5,5) * t267 + Ifges(5,6) * t266 + Ifges(5,3) * t274;
t228 = mrSges(5,2) * t238 - mrSges(5,3) * t236 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t263 + t266 * t243 - t274 * t244;
t227 = -mrSges(5,1) * t238 + mrSges(5,3) * t237 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t263 - t267 * t243 + t274 * t245;
t220 = -mrSges(4,1) * t255 + mrSges(4,3) * t242 + Ifges(4,4) * t265 + Ifges(4,2) * t264 + Ifges(4,6) * qJDD(2) - pkin(3) * t226 + qJD(2) * t260 - t276 * t258 - t301;
t219 = mrSges(4,2) * t255 - mrSges(4,3) * t241 + Ifges(4,1) * t265 + Ifges(4,4) * t264 + Ifges(4,5) * qJDD(2) - pkin(6) * t226 - qJD(2) * t259 - t293 * t227 + t296 * t228 - t275 * t258;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t309 - mrSges(2,2) * t306 + t294 * (mrSges(3,2) * t280 - mrSges(3,3) * t268 + Ifges(3,1) * t284 + Ifges(3,4) * t285 + Ifges(3,5) * qJDD(2) - qJ(3) * t221 - qJD(2) * t278 + t292 * t219 - t291 * t220) + t297 * (-mrSges(3,1) * t280 + mrSges(3,3) * t269 + Ifges(3,4) * t284 + Ifges(3,2) * t285 + Ifges(3,6) * qJDD(2) - pkin(2) * t225 + qJ(3) * t308 + qJD(2) * t279 + t291 * t219 + t292 * t220) + pkin(1) * (-m(3) * t280 + t285 * mrSges(3,1) - t284 * mrSges(3,2) + (-t287 * t294 + t288 * t297) * qJD(1) - t225) + pkin(5) * (t297 * (m(3) * t269 - qJDD(2) * mrSges(3,2) + t285 * mrSges(3,3) - qJD(2) * t287 + t283 * t311 + t308) - t294 * (m(3) * t268 + qJDD(2) * mrSges(3,1) - t284 * mrSges(3,3) + qJD(2) * t288 - t283 * t312 + t221)); Ifges(3,5) * t284 + Ifges(3,6) * t285 + mrSges(3,1) * t268 - mrSges(3,2) * t269 + Ifges(4,5) * t265 + Ifges(4,6) * t264 + t276 * t259 + t275 * t260 + mrSges(4,1) * t241 - mrSges(4,2) * t242 + t293 * t228 + t296 * t227 + pkin(3) * t302 + pkin(6) * t307 + pkin(2) * t221 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t294 * t278 - t297 * t279) * qJD(1); t225; t301;];
tauJ = t1;
