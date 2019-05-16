% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:13:32
% EndTime: 2019-05-05 18:13:58
% DurationCPUTime: 19.90s
% Computational Cost: add. (370009->341), mult. (815835->433), div. (0->0), fcn. (569717->12), ass. (0->137)
t305 = sin(qJ(1));
t309 = cos(qJ(1));
t284 = t305 * g(1) - g(2) * t309;
t275 = qJDD(1) * pkin(1) + t284;
t285 = -g(1) * t309 - g(2) * t305;
t310 = qJD(1) ^ 2;
t277 = -pkin(1) * t310 + t285;
t299 = sin(pkin(10));
t301 = cos(pkin(10));
t255 = t299 * t275 + t301 * t277;
t246 = -pkin(2) * t310 + qJDD(1) * pkin(7) + t255;
t297 = -g(3) + qJDD(2);
t304 = sin(qJ(3));
t308 = cos(qJ(3));
t234 = -t304 * t246 + t308 * t297;
t329 = qJD(1) * qJD(3);
t328 = t308 * t329;
t278 = qJDD(1) * t304 + t328;
t228 = (-t278 + t328) * qJ(4) + (t304 * t308 * t310 + qJDD(3)) * pkin(3) + t234;
t235 = t308 * t246 + t304 * t297;
t279 = qJDD(1) * t308 - t304 * t329;
t331 = qJD(1) * t304;
t281 = qJD(3) * pkin(3) - qJ(4) * t331;
t296 = t308 ^ 2;
t229 = -pkin(3) * t296 * t310 + qJ(4) * t279 - qJD(3) * t281 + t235;
t298 = sin(pkin(11));
t300 = cos(pkin(11));
t265 = (t298 * t308 + t300 * t304) * qJD(1);
t196 = -0.2e1 * qJD(4) * t265 + t300 * t228 - t298 * t229;
t257 = t278 * t300 + t279 * t298;
t264 = (-t298 * t304 + t300 * t308) * qJD(1);
t192 = (qJD(3) * t264 - t257) * pkin(8) + (t264 * t265 + qJDD(3)) * pkin(4) + t196;
t197 = 0.2e1 * qJD(4) * t264 + t298 * t228 + t300 * t229;
t256 = -t278 * t298 + t279 * t300;
t260 = qJD(3) * pkin(4) - pkin(8) * t265;
t263 = t264 ^ 2;
t194 = -pkin(4) * t263 + pkin(8) * t256 - qJD(3) * t260 + t197;
t303 = sin(qJ(5));
t307 = cos(qJ(5));
t189 = t303 * t192 + t307 * t194;
t244 = t264 * t303 + t265 * t307;
t213 = -qJD(5) * t244 + t256 * t307 - t257 * t303;
t243 = t264 * t307 - t265 * t303;
t226 = -mrSges(6,1) * t243 + mrSges(6,2) * t244;
t292 = qJD(3) + qJD(5);
t237 = mrSges(6,1) * t292 - mrSges(6,3) * t244;
t291 = qJDD(3) + qJDD(5);
t227 = -pkin(5) * t243 - pkin(9) * t244;
t290 = t292 ^ 2;
t186 = -pkin(5) * t290 + pkin(9) * t291 + t227 * t243 + t189;
t254 = t301 * t275 - t299 * t277;
t321 = -qJDD(1) * pkin(2) - t254;
t230 = -t279 * pkin(3) + qJDD(4) + t281 * t331 + (-qJ(4) * t296 - pkin(7)) * t310 + t321;
t202 = -t256 * pkin(4) - t263 * pkin(8) + t265 * t260 + t230;
t214 = qJD(5) * t243 + t256 * t303 + t257 * t307;
t190 = (-t243 * t292 - t214) * pkin(9) + (t244 * t292 - t213) * pkin(5) + t202;
t302 = sin(qJ(6));
t306 = cos(qJ(6));
t183 = -t186 * t302 + t190 * t306;
t232 = -t244 * t302 + t292 * t306;
t200 = qJD(6) * t232 + t214 * t306 + t291 * t302;
t212 = qJDD(6) - t213;
t233 = t244 * t306 + t292 * t302;
t215 = -mrSges(7,1) * t232 + mrSges(7,2) * t233;
t239 = qJD(6) - t243;
t216 = -mrSges(7,2) * t239 + mrSges(7,3) * t232;
t179 = m(7) * t183 + mrSges(7,1) * t212 - mrSges(7,3) * t200 - t215 * t233 + t216 * t239;
t184 = t186 * t306 + t190 * t302;
t199 = -qJD(6) * t233 - t214 * t302 + t291 * t306;
t217 = mrSges(7,1) * t239 - mrSges(7,3) * t233;
t180 = m(7) * t184 - mrSges(7,2) * t212 + mrSges(7,3) * t199 + t215 * t232 - t217 * t239;
t323 = -t179 * t302 + t306 * t180;
t166 = m(6) * t189 - mrSges(6,2) * t291 + mrSges(6,3) * t213 + t226 * t243 - t237 * t292 + t323;
t188 = t192 * t307 - t194 * t303;
t236 = -mrSges(6,2) * t292 + mrSges(6,3) * t243;
t185 = -pkin(5) * t291 - pkin(9) * t290 + t227 * t244 - t188;
t318 = -m(7) * t185 + t199 * mrSges(7,1) - mrSges(7,2) * t200 + t232 * t216 - t217 * t233;
t175 = m(6) * t188 + mrSges(6,1) * t291 - mrSges(6,3) * t214 - t226 * t244 + t236 * t292 + t318;
t160 = t303 * t166 + t307 * t175;
t249 = -mrSges(5,1) * t264 + mrSges(5,2) * t265;
t258 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t264;
t157 = m(5) * t196 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t257 + qJD(3) * t258 - t249 * t265 + t160;
t259 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t265;
t324 = t307 * t166 - t175 * t303;
t158 = m(5) * t197 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t256 - qJD(3) * t259 + t249 * t264 + t324;
t151 = t300 * t157 + t298 * t158;
t270 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t304 + Ifges(4,2) * t308) * qJD(1);
t271 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t304 + Ifges(4,4) * t308) * qJD(1);
t241 = Ifges(5,4) * t265 + Ifges(5,2) * t264 + Ifges(5,6) * qJD(3);
t242 = Ifges(5,1) * t265 + Ifges(5,4) * t264 + Ifges(5,5) * qJD(3);
t203 = Ifges(7,5) * t233 + Ifges(7,6) * t232 + Ifges(7,3) * t239;
t205 = Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t239;
t172 = -mrSges(7,1) * t185 + mrSges(7,3) * t184 + Ifges(7,4) * t200 + Ifges(7,2) * t199 + Ifges(7,6) * t212 - t203 * t233 + t205 * t239;
t204 = Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t239;
t173 = mrSges(7,2) * t185 - mrSges(7,3) * t183 + Ifges(7,1) * t200 + Ifges(7,4) * t199 + Ifges(7,5) * t212 + t203 * t232 - t204 * t239;
t219 = Ifges(6,4) * t244 + Ifges(6,2) * t243 + Ifges(6,6) * t292;
t220 = Ifges(6,1) * t244 + Ifges(6,4) * t243 + Ifges(6,5) * t292;
t316 = -mrSges(6,1) * t188 + mrSges(6,2) * t189 - Ifges(6,5) * t214 - Ifges(6,6) * t213 - Ifges(6,3) * t291 - pkin(5) * t318 - pkin(9) * t323 - t306 * t172 - t302 * t173 - t244 * t219 + t243 * t220;
t313 = -mrSges(5,1) * t196 + mrSges(5,2) * t197 - Ifges(5,5) * t257 - Ifges(5,6) * t256 - Ifges(5,3) * qJDD(3) - pkin(4) * t160 - t265 * t241 + t264 * t242 + t316;
t332 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t278 + Ifges(4,6) * t279 + Ifges(4,3) * qJDD(3) + pkin(3) * t151 + (t270 * t304 - t271 * t308) * qJD(1) - t313;
t276 = (-mrSges(4,1) * t308 + mrSges(4,2) * t304) * qJD(1);
t330 = qJD(1) * t308;
t283 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t330;
t149 = m(4) * t234 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t278 + qJD(3) * t283 - t276 * t331 + t151;
t282 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t331;
t325 = -t157 * t298 + t300 * t158;
t150 = m(4) * t235 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t279 - qJD(3) * t282 + t276 * t330 + t325;
t326 = -t149 * t304 + t308 * t150;
t142 = m(3) * t255 - mrSges(3,1) * t310 - qJDD(1) * mrSges(3,2) + t326;
t245 = -t310 * pkin(7) + t321;
t168 = t306 * t179 + t302 * t180;
t320 = m(6) * t202 - t213 * mrSges(6,1) + t214 * mrSges(6,2) - t243 * t236 + t244 * t237 + t168;
t315 = m(5) * t230 - t256 * mrSges(5,1) + mrSges(5,2) * t257 - t264 * t258 + t259 * t265 + t320;
t312 = -m(4) * t245 + t279 * mrSges(4,1) - mrSges(4,2) * t278 - t282 * t331 + t283 * t330 - t315;
t162 = m(3) * t254 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t310 + t312;
t138 = t299 * t142 + t301 * t162;
t144 = t308 * t149 + t304 * t150;
t327 = t301 * t142 - t162 * t299;
t218 = Ifges(6,5) * t244 + Ifges(6,6) * t243 + Ifges(6,3) * t292;
t152 = mrSges(6,2) * t202 - mrSges(6,3) * t188 + Ifges(6,1) * t214 + Ifges(6,4) * t213 + Ifges(6,5) * t291 - pkin(9) * t168 - t172 * t302 + t173 * t306 + t218 * t243 - t219 * t292;
t314 = mrSges(7,1) * t183 - mrSges(7,2) * t184 + Ifges(7,5) * t200 + Ifges(7,6) * t199 + Ifges(7,3) * t212 + t204 * t233 - t205 * t232;
t153 = -mrSges(6,1) * t202 + mrSges(6,3) * t189 + Ifges(6,4) * t214 + Ifges(6,2) * t213 + Ifges(6,6) * t291 - pkin(5) * t168 - t218 * t244 + t220 * t292 - t314;
t240 = Ifges(5,5) * t265 + Ifges(5,6) * t264 + Ifges(5,3) * qJD(3);
t139 = -mrSges(5,1) * t230 + mrSges(5,3) * t197 + Ifges(5,4) * t257 + Ifges(5,2) * t256 + Ifges(5,6) * qJDD(3) - pkin(4) * t320 + pkin(8) * t324 + qJD(3) * t242 + t303 * t152 + t307 * t153 - t265 * t240;
t145 = mrSges(5,2) * t230 - mrSges(5,3) * t196 + Ifges(5,1) * t257 + Ifges(5,4) * t256 + Ifges(5,5) * qJDD(3) - pkin(8) * t160 - qJD(3) * t241 + t152 * t307 - t153 * t303 + t240 * t264;
t269 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t304 + Ifges(4,6) * t308) * qJD(1);
t131 = -mrSges(4,1) * t245 + mrSges(4,3) * t235 + Ifges(4,4) * t278 + Ifges(4,2) * t279 + Ifges(4,6) * qJDD(3) - pkin(3) * t315 + qJ(4) * t325 + qJD(3) * t271 + t300 * t139 + t298 * t145 - t269 * t331;
t133 = mrSges(4,2) * t245 - mrSges(4,3) * t234 + Ifges(4,1) * t278 + Ifges(4,4) * t279 + Ifges(4,5) * qJDD(3) - qJ(4) * t151 - qJD(3) * t270 - t139 * t298 + t145 * t300 + t269 * t330;
t319 = mrSges(3,1) * t254 - mrSges(3,2) * t255 + Ifges(3,3) * qJDD(1) + pkin(2) * t312 + pkin(7) * t326 + t308 * t131 + t304 * t133;
t317 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + Ifges(2,3) * qJDD(1) + pkin(1) * t138 + t319;
t136 = m(2) * t285 - mrSges(2,1) * t310 - qJDD(1) * mrSges(2,2) + t327;
t135 = m(2) * t284 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t310 + t138;
t134 = -mrSges(3,1) * t297 + mrSges(3,3) * t255 + t310 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t144 - t332;
t129 = mrSges(3,2) * t297 - mrSges(3,3) * t254 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t310 - pkin(7) * t144 - t131 * t304 + t133 * t308;
t128 = -mrSges(2,2) * g(3) - mrSges(2,3) * t284 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t310 - qJ(2) * t138 + t129 * t301 - t134 * t299;
t127 = Ifges(2,6) * qJDD(1) + t310 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t285 + t299 * t129 + t301 * t134 - pkin(1) * (m(3) * t297 + t144) + qJ(2) * t327;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t128 - t305 * t127 - pkin(6) * (t135 * t309 + t136 * t305), t128, t129, t133, t145, t152, t173; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t128 + t309 * t127 + pkin(6) * (-t135 * t305 + t136 * t309), t127, t134, t131, t139, t153, t172; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t317, t317, t319, t332, -t313, -t316, t314;];
m_new  = t1;
