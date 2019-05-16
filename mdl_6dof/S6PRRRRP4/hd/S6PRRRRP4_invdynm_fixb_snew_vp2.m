% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:50:57
% EndTime: 2019-05-05 09:51:28
% DurationCPUTime: 15.12s
% Computational Cost: add. (288230->336), mult. (555432->416), div. (0->0), fcn. (394113->12), ass. (0->135)
t294 = sin(pkin(11));
t296 = cos(pkin(11));
t284 = g(1) * t294 - g(2) * t296;
t285 = -g(1) * t296 - g(2) * t294;
t293 = -g(3) + qJDD(1);
t304 = cos(qJ(2));
t297 = cos(pkin(6));
t301 = sin(qJ(2));
t328 = t297 * t301;
t295 = sin(pkin(6));
t329 = t295 * t301;
t241 = t284 * t328 + t304 * t285 + t293 * t329;
t306 = qJD(2) ^ 2;
t235 = -pkin(2) * t306 + qJDD(2) * pkin(8) + t241;
t258 = -t284 * t295 + t293 * t297;
t300 = sin(qJ(3));
t303 = cos(qJ(3));
t226 = t303 * t235 + t300 * t258;
t280 = (-pkin(3) * t303 - pkin(9) * t300) * qJD(2);
t305 = qJD(3) ^ 2;
t324 = qJD(2) * t303;
t211 = -pkin(3) * t305 + qJDD(3) * pkin(9) + t280 * t324 + t226;
t240 = -t301 * t285 + (t284 * t297 + t293 * t295) * t304;
t234 = -qJDD(2) * pkin(2) - t306 * pkin(8) - t240;
t323 = qJD(2) * qJD(3);
t321 = t303 * t323;
t281 = qJDD(2) * t300 + t321;
t292 = t300 * t323;
t282 = qJDD(2) * t303 - t292;
t215 = (-t281 - t321) * pkin(9) + (-t282 + t292) * pkin(3) + t234;
t299 = sin(qJ(4));
t302 = cos(qJ(4));
t193 = -t299 * t211 + t302 * t215;
t325 = qJD(2) * t300;
t277 = qJD(3) * t302 - t299 * t325;
t250 = qJD(4) * t277 + qJDD(3) * t299 + t281 * t302;
t274 = qJDD(4) - t282;
t278 = qJD(3) * t299 + t302 * t325;
t291 = qJD(4) - t324;
t190 = (t277 * t291 - t250) * pkin(10) + (t277 * t278 + t274) * pkin(4) + t193;
t194 = t302 * t211 + t299 * t215;
t249 = -qJD(4) * t278 + qJDD(3) * t302 - t281 * t299;
t257 = pkin(4) * t291 - pkin(10) * t278;
t273 = t277 ^ 2;
t192 = -pkin(4) * t273 + pkin(10) * t249 - t257 * t291 + t194;
t298 = sin(qJ(5));
t334 = cos(qJ(5));
t186 = t298 * t190 + t334 * t192;
t252 = t298 * t277 + t278 * t334;
t207 = t252 * qJD(5) - t249 * t334 + t298 * t250;
t290 = qJD(5) + t291;
t238 = mrSges(6,1) * t290 - mrSges(6,3) * t252;
t251 = -t277 * t334 + t298 * t278;
t270 = qJDD(5) + t274;
t227 = pkin(5) * t251 - qJ(6) * t252;
t289 = t290 ^ 2;
t181 = -pkin(5) * t289 + qJ(6) * t270 + 0.2e1 * qJD(6) * t290 - t227 * t251 + t186;
t239 = -mrSges(7,1) * t290 + mrSges(7,2) * t252;
t322 = m(7) * t181 + t270 * mrSges(7,3) + t290 * t239;
t228 = mrSges(7,1) * t251 - mrSges(7,3) * t252;
t326 = -mrSges(6,1) * t251 - mrSges(6,2) * t252 - t228;
t332 = -mrSges(6,3) - mrSges(7,2);
t170 = m(6) * t186 - t270 * mrSges(6,2) + t207 * t332 - t290 * t238 + t251 * t326 + t322;
t185 = t190 * t334 - t298 * t192;
t208 = -t251 * qJD(5) + t298 * t249 + t250 * t334;
t237 = -mrSges(6,2) * t290 - mrSges(6,3) * t251;
t183 = -t270 * pkin(5) - t289 * qJ(6) + t252 * t227 + qJDD(6) - t185;
t236 = -mrSges(7,2) * t251 + mrSges(7,3) * t290;
t318 = -m(7) * t183 + t270 * mrSges(7,1) + t290 * t236;
t172 = m(6) * t185 + t270 * mrSges(6,1) + t208 * t332 + t290 * t237 + t252 * t326 + t318;
t165 = t298 * t170 + t334 * t172;
t253 = -mrSges(5,1) * t277 + mrSges(5,2) * t278;
t255 = -mrSges(5,2) * t291 + mrSges(5,3) * t277;
t161 = m(5) * t193 + mrSges(5,1) * t274 - mrSges(5,3) * t250 - t253 * t278 + t255 * t291 + t165;
t256 = mrSges(5,1) * t291 - mrSges(5,3) * t278;
t319 = t334 * t170 - t172 * t298;
t162 = m(5) * t194 - mrSges(5,2) * t274 + mrSges(5,3) * t249 + t253 * t277 - t256 * t291 + t319;
t159 = -t161 * t299 + t302 * t162;
t279 = (-mrSges(4,1) * t303 + mrSges(4,2) * t300) * qJD(2);
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t325;
t157 = m(4) * t226 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t282 - qJD(3) * t286 + t279 * t324 + t159;
t225 = -t300 * t235 + t303 * t258;
t210 = -qJDD(3) * pkin(3) - t305 * pkin(9) + t280 * t325 - t225;
t195 = -t249 * pkin(4) - t273 * pkin(10) + t278 * t257 + t210;
t188 = -0.2e1 * qJD(6) * t252 + (t251 * t290 - t208) * qJ(6) + (t252 * t290 + t207) * pkin(5) + t195;
t178 = m(7) * t188 + t207 * mrSges(7,1) - t208 * mrSges(7,3) + t251 * t236 - t252 * t239;
t311 = m(6) * t195 + t207 * mrSges(6,1) + t208 * mrSges(6,2) + t251 * t237 + t252 * t238 + t178;
t173 = -m(5) * t210 + t249 * mrSges(5,1) - t250 * mrSges(5,2) + t277 * t255 - t278 * t256 - t311;
t287 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t324;
t167 = m(4) * t225 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,3) + qJD(3) * t287 - t279 * t325 + t173;
t152 = t300 * t157 + t303 * t167;
t220 = Ifges(7,1) * t252 + Ifges(7,4) * t290 + Ifges(7,5) * t251;
t221 = Ifges(6,1) * t252 - Ifges(6,4) * t251 + Ifges(6,5) * t290;
t317 = -mrSges(7,1) * t188 + mrSges(7,2) * t181;
t218 = Ifges(7,4) * t252 + Ifges(7,2) * t290 + Ifges(7,6) * t251;
t327 = -Ifges(6,5) * t252 + Ifges(6,6) * t251 - Ifges(6,3) * t290 - t218;
t163 = -mrSges(6,1) * t195 + mrSges(6,3) * t186 - pkin(5) * t178 + (t220 + t221) * t290 + (Ifges(6,6) - Ifges(7,6)) * t270 + t327 * t252 + (Ifges(6,4) - Ifges(7,5)) * t208 + (-Ifges(6,2) - Ifges(7,3)) * t207 + t317;
t219 = Ifges(6,4) * t252 - Ifges(6,2) * t251 + Ifges(6,6) * t290;
t216 = Ifges(7,5) * t252 + Ifges(7,6) * t290 + Ifges(7,3) * t251;
t314 = mrSges(7,2) * t183 - mrSges(7,3) * t188 + Ifges(7,1) * t208 + Ifges(7,4) * t270 + Ifges(7,5) * t207 + t290 * t216;
t164 = mrSges(6,2) * t195 - mrSges(6,3) * t185 + Ifges(6,1) * t208 - Ifges(6,4) * t207 + Ifges(6,5) * t270 - qJ(6) * t178 - t290 * t219 + t251 * t327 + t314;
t243 = Ifges(5,5) * t278 + Ifges(5,6) * t277 + Ifges(5,3) * t291;
t245 = Ifges(5,1) * t278 + Ifges(5,4) * t277 + Ifges(5,5) * t291;
t146 = -mrSges(5,1) * t210 + mrSges(5,3) * t194 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * t274 - pkin(4) * t311 + pkin(10) * t319 + t163 * t334 + t298 * t164 - t278 * t243 + t291 * t245;
t244 = Ifges(5,4) * t278 + Ifges(5,2) * t277 + Ifges(5,6) * t291;
t147 = mrSges(5,2) * t210 - mrSges(5,3) * t193 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * t274 - pkin(10) * t165 - t298 * t163 + t164 * t334 + t277 * t243 - t291 * t244;
t268 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t300 + Ifges(4,2) * t303) * qJD(2);
t269 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t300 + Ifges(4,4) * t303) * qJD(2);
t335 = mrSges(4,1) * t225 - mrSges(4,2) * t226 + Ifges(4,5) * t281 + Ifges(4,6) * t282 + Ifges(4,3) * qJDD(3) + pkin(3) * t173 + pkin(9) * t159 + t302 * t146 + t299 * t147 + (t268 * t300 - t269 * t303) * qJD(2);
t136 = -mrSges(3,1) * t258 + mrSges(3,3) * t241 + t306 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t152 - t335;
t320 = t303 * t157 - t167 * t300;
t150 = m(3) * t241 - mrSges(3,1) * t306 - qJDD(2) * mrSges(3,2) + t320;
t158 = t161 * t302 + t162 * t299;
t310 = -m(4) * t234 + t282 * mrSges(4,1) - mrSges(4,2) * t281 - t286 * t325 + t287 * t324 - t158;
t154 = m(3) * t240 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t306 + t310;
t144 = t304 * t150 - t154 * t301;
t336 = pkin(7) * t144 + t136 * t304;
t330 = t154 * t304;
t151 = m(3) * t258 + t152;
t141 = t150 * t328 - t151 * t295 + t297 * t330;
t267 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t300 + Ifges(4,6) * t303) * qJD(2);
t137 = mrSges(4,2) * t234 - mrSges(4,3) * t225 + Ifges(4,1) * t281 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - pkin(9) * t158 - qJD(3) * t268 - t146 * t299 + t147 * t302 + t267 * t324;
t312 = mrSges(7,1) * t183 - mrSges(7,3) * t181 - Ifges(7,4) * t208 - Ifges(7,2) * t270 - Ifges(7,6) * t207 + t252 * t216 - t251 * t220;
t309 = mrSges(6,2) * t186 - t251 * t221 - qJ(6) * (-t207 * mrSges(7,2) - t251 * t228 + t322) - pkin(5) * (-t208 * mrSges(7,2) - t252 * t228 + t318) - mrSges(6,1) * t185 + Ifges(6,6) * t207 - Ifges(6,5) * t208 - t252 * t219 - Ifges(6,3) * t270 + t312;
t307 = mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t250 + Ifges(5,6) * t249 + Ifges(5,3) * t274 + pkin(4) * t165 + t278 * t244 - t277 * t245 - t309;
t145 = -mrSges(4,1) * t234 + mrSges(4,3) * t226 + Ifges(4,4) * t281 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t158 + qJD(3) * t269 - t267 * t325 - t307;
t132 = mrSges(3,1) * t240 - mrSges(3,2) * t241 + Ifges(3,3) * qJDD(2) + pkin(2) * t310 + pkin(8) * t320 + t300 * t137 + t303 * t145;
t134 = mrSges(3,2) * t258 - mrSges(3,3) * t240 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t306 - pkin(8) * t152 + t137 * t303 - t145 * t300;
t313 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + pkin(1) * t141 + t297 * t132 + t134 * t329 + t295 * t336;
t142 = m(2) * t285 + t144;
t140 = t297 * t151 + (t150 * t301 + t330) * t295;
t138 = m(2) * t284 + t141;
t130 = mrSges(2,2) * t293 - mrSges(2,3) * t284 + t304 * t134 - t301 * t136 + (-t140 * t295 - t141 * t297) * pkin(7);
t129 = -mrSges(2,1) * t293 + mrSges(2,3) * t285 - pkin(1) * t140 - t295 * t132 + (t134 * t301 + t336) * t297;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t296 * t130 - t294 * t129 - qJ(1) * (t138 * t296 + t142 * t294), t130, t134, t137, t147, t164, -t218 * t251 + t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t294 * t130 + t296 * t129 + qJ(1) * (-t138 * t294 + t142 * t296), t129, t136, t145, t146, t163, -t312; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t313, t313, t132, t335, t307, -t309, Ifges(7,5) * t208 + Ifges(7,6) * t270 + Ifges(7,3) * t207 + t252 * t218 - t290 * t220 - t317;];
m_new  = t1;
