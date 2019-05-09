% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:17:43
% EndTime: 2019-05-06 01:18:05
% DurationCPUTime: 9.92s
% Computational Cost: add. (186494->337), mult. (359101->408), div. (0->0), fcn. (235583->10), ass. (0->128)
t292 = sin(qJ(1));
t295 = cos(qJ(1));
t276 = t292 * g(1) - g(2) * t295;
t267 = qJDD(1) * pkin(1) + t276;
t277 = -g(1) * t295 - g(2) * t292;
t297 = qJD(1) ^ 2;
t269 = -pkin(1) * t297 + t277;
t287 = sin(pkin(10));
t288 = cos(pkin(10));
t242 = t287 * t267 + t288 * t269;
t229 = -pkin(2) * t297 + qJDD(1) * pkin(7) + t242;
t286 = -g(3) + qJDD(2);
t291 = sin(qJ(3));
t294 = cos(qJ(3));
t219 = -t291 * t229 + t294 * t286;
t270 = (-pkin(3) * t294 - pkin(8) * t291) * qJD(1);
t296 = qJD(3) ^ 2;
t317 = qJD(1) * t291;
t212 = -qJDD(3) * pkin(3) - t296 * pkin(8) + t270 * t317 - t219;
t290 = sin(qJ(4));
t293 = cos(qJ(4));
t266 = qJD(3) * t290 + t293 * t317;
t315 = qJD(1) * qJD(3);
t313 = t294 * t315;
t271 = qJDD(1) * t291 + t313;
t236 = -qJD(4) * t266 + qJDD(3) * t293 - t271 * t290;
t316 = qJD(1) * t294;
t280 = qJD(4) - t316;
t246 = pkin(4) * t280 - pkin(9) * t266;
t265 = qJD(3) * t293 - t290 * t317;
t263 = t265 ^ 2;
t182 = -t236 * pkin(4) - t263 * pkin(9) + t266 * t246 + t212;
t237 = qJD(4) * t265 + qJDD(3) * t290 + t271 * t293;
t289 = sin(qJ(5));
t321 = cos(qJ(5));
t240 = t289 * t265 + t321 * t266;
t196 = t240 * qJD(5) - t321 * t236 + t289 * t237;
t239 = -t321 * t265 + t289 * t266;
t197 = -t239 * qJD(5) + t289 * t236 + t321 * t237;
t279 = qJD(5) + t280;
t174 = -0.2e1 * qJD(6) * t240 + (t239 * t279 - t197) * qJ(6) + (t240 * t279 + t196) * pkin(5) + t182;
t221 = -mrSges(7,2) * t239 + mrSges(7,3) * t279;
t224 = -mrSges(7,1) * t279 + mrSges(7,2) * t240;
t167 = m(7) * t174 + t196 * mrSges(7,1) - t197 * mrSges(7,3) + t239 * t221 - t240 * t224;
t241 = t288 * t267 - t287 * t269;
t228 = -qJDD(1) * pkin(2) - t297 * pkin(7) - t241;
t282 = t291 * t315;
t272 = qJDD(1) * t294 - t282;
t207 = (-t271 - t313) * pkin(8) + (-t272 + t282) * pkin(3) + t228;
t220 = t294 * t229 + t291 * t286;
t213 = -pkin(3) * t296 + qJDD(3) * pkin(8) + t270 * t316 + t220;
t183 = t293 * t207 - t290 * t213;
t264 = qJDD(4) - t272;
t179 = (t265 * t280 - t237) * pkin(9) + (t265 * t266 + t264) * pkin(4) + t183;
t184 = t290 * t207 + t293 * t213;
t181 = -pkin(4) * t263 + pkin(9) * t236 - t246 * t280 + t184;
t177 = t289 * t179 + t321 * t181;
t205 = Ifges(7,1) * t240 + Ifges(7,4) * t279 + Ifges(7,5) * t239;
t206 = Ifges(6,1) * t240 - Ifges(6,4) * t239 + Ifges(6,5) * t279;
t260 = qJDD(5) + t264;
t214 = pkin(5) * t239 - qJ(6) * t240;
t278 = t279 ^ 2;
t170 = -pkin(5) * t278 + qJ(6) * t260 + 0.2e1 * qJD(6) * t279 - t214 * t239 + t177;
t308 = -mrSges(7,1) * t174 + mrSges(7,2) * t170;
t203 = Ifges(7,4) * t240 + Ifges(7,2) * t279 + Ifges(7,6) * t239;
t319 = -Ifges(6,5) * t240 + Ifges(6,6) * t239 - Ifges(6,3) * t279 - t203;
t152 = -mrSges(6,1) * t182 + mrSges(6,3) * t177 - pkin(5) * t167 + (t205 + t206) * t279 + (Ifges(6,6) - Ifges(7,6)) * t260 + t319 * t240 + (Ifges(6,4) - Ifges(7,5)) * t197 + (-Ifges(6,2) - Ifges(7,3)) * t196 + t308;
t176 = t321 * t179 - t289 * t181;
t204 = Ifges(6,4) * t240 - Ifges(6,2) * t239 + Ifges(6,6) * t279;
t172 = -t260 * pkin(5) - t278 * qJ(6) + t240 * t214 + qJDD(6) - t176;
t201 = Ifges(7,5) * t240 + Ifges(7,6) * t279 + Ifges(7,3) * t239;
t306 = mrSges(7,2) * t172 - mrSges(7,3) * t174 + Ifges(7,1) * t197 + Ifges(7,4) * t260 + Ifges(7,5) * t196 + t279 * t201;
t153 = mrSges(6,2) * t182 - mrSges(6,3) * t176 + Ifges(6,1) * t197 - Ifges(6,4) * t196 + Ifges(6,5) * t260 - qJ(6) * t167 - t279 * t204 + t319 * t239 + t306;
t230 = Ifges(5,5) * t266 + Ifges(5,6) * t265 + Ifges(5,3) * t280;
t232 = Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * t280;
t222 = -mrSges(6,2) * t279 - mrSges(6,3) * t239;
t223 = mrSges(6,1) * t279 - mrSges(6,3) * t240;
t302 = m(6) * t182 + t196 * mrSges(6,1) + t197 * mrSges(6,2) + t239 * t222 + t240 * t223 + t167;
t314 = m(7) * t170 + t260 * mrSges(7,3) + t279 * t224;
t215 = mrSges(7,1) * t239 - mrSges(7,3) * t240;
t318 = -mrSges(6,1) * t239 - mrSges(6,2) * t240 - t215;
t320 = -mrSges(6,3) - mrSges(7,2);
t159 = m(6) * t177 - t260 * mrSges(6,2) + t320 * t196 - t279 * t223 + t318 * t239 + t314;
t309 = -m(7) * t172 + t260 * mrSges(7,1) + t279 * t221;
t161 = m(6) * t176 + t260 * mrSges(6,1) + t320 * t197 + t279 * t222 + t318 * t240 + t309;
t310 = t321 * t159 - t161 * t289;
t134 = -mrSges(5,1) * t212 + mrSges(5,3) * t184 + Ifges(5,4) * t237 + Ifges(5,2) * t236 + Ifges(5,6) * t264 - pkin(4) * t302 + pkin(9) * t310 + t321 * t152 + t289 * t153 - t266 * t230 + t280 * t232;
t154 = t289 * t159 + t321 * t161;
t231 = Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * t280;
t135 = mrSges(5,2) * t212 - mrSges(5,3) * t183 + Ifges(5,1) * t237 + Ifges(5,4) * t236 + Ifges(5,5) * t264 - pkin(9) * t154 - t289 * t152 + t321 * t153 + t265 * t230 - t280 * t231;
t243 = -mrSges(5,1) * t265 + mrSges(5,2) * t266;
t244 = -mrSges(5,2) * t280 + mrSges(5,3) * t265;
t150 = m(5) * t183 + mrSges(5,1) * t264 - mrSges(5,3) * t237 - t243 * t266 + t244 * t280 + t154;
t245 = mrSges(5,1) * t280 - mrSges(5,3) * t266;
t151 = m(5) * t184 - mrSges(5,2) * t264 + mrSges(5,3) * t236 + t243 * t265 - t245 * t280 + t310;
t148 = -t150 * t290 + t293 * t151;
t162 = -m(5) * t212 + t236 * mrSges(5,1) - t237 * mrSges(5,2) + t265 * t244 - t266 * t245 - t302;
t258 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t291 + Ifges(4,2) * t294) * qJD(1);
t259 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t291 + Ifges(4,4) * t294) * qJD(1);
t322 = mrSges(4,1) * t219 - mrSges(4,2) * t220 + Ifges(4,5) * t271 + Ifges(4,6) * t272 + Ifges(4,3) * qJDD(3) + pkin(3) * t162 + pkin(8) * t148 + t293 * t134 + t290 * t135 + (t258 * t291 - t259 * t294) * qJD(1);
t268 = (-mrSges(4,1) * t294 + mrSges(4,2) * t291) * qJD(1);
t274 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t317;
t146 = m(4) * t220 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t272 - qJD(3) * t274 + t268 * t316 + t148;
t275 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t316;
t156 = m(4) * t219 + qJDD(3) * mrSges(4,1) - t271 * mrSges(4,3) + qJD(3) * t275 - t268 * t317 + t162;
t311 = t294 * t146 - t156 * t291;
t138 = m(3) * t242 - mrSges(3,1) * t297 - qJDD(1) * mrSges(3,2) + t311;
t147 = t150 * t293 + t151 * t290;
t301 = -m(4) * t228 + t272 * mrSges(4,1) - mrSges(4,2) * t271 - t274 * t317 + t275 * t316 - t147;
t142 = m(3) * t241 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t297 + t301;
t131 = t287 * t138 + t288 * t142;
t140 = t291 * t146 + t294 * t156;
t312 = t288 * t138 - t142 * t287;
t257 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t291 + Ifges(4,6) * t294) * qJD(1);
t127 = mrSges(4,2) * t228 - mrSges(4,3) * t219 + Ifges(4,1) * t271 + Ifges(4,4) * t272 + Ifges(4,5) * qJDD(3) - pkin(8) * t147 - qJD(3) * t258 - t134 * t290 + t135 * t293 + t257 * t316;
t304 = mrSges(7,1) * t172 - mrSges(7,3) * t170 - Ifges(7,4) * t197 - Ifges(7,2) * t260 - Ifges(7,6) * t196 + t240 * t201 - t239 * t205;
t300 = mrSges(6,2) * t177 - t239 * t206 - qJ(6) * (-t196 * mrSges(7,2) - t239 * t215 + t314) - pkin(5) * (-t197 * mrSges(7,2) - t240 * t215 + t309) - mrSges(6,1) * t176 + Ifges(6,6) * t196 - Ifges(6,5) * t197 - t240 * t204 - Ifges(6,3) * t260 + t304;
t298 = mrSges(5,1) * t183 - mrSges(5,2) * t184 + Ifges(5,5) * t237 + Ifges(5,6) * t236 + Ifges(5,3) * t264 + pkin(4) * t154 + t266 * t231 - t265 * t232 - t300;
t133 = -mrSges(4,1) * t228 + mrSges(4,3) * t220 + Ifges(4,4) * t271 + Ifges(4,2) * t272 + Ifges(4,6) * qJDD(3) - pkin(3) * t147 + qJD(3) * t259 - t257 * t317 - t298;
t305 = mrSges(3,1) * t241 - mrSges(3,2) * t242 + Ifges(3,3) * qJDD(1) + pkin(2) * t301 + pkin(7) * t311 + t291 * t127 + t294 * t133;
t303 = mrSges(2,1) * t276 - mrSges(2,2) * t277 + Ifges(2,3) * qJDD(1) + pkin(1) * t131 + t305;
t129 = m(2) * t277 - mrSges(2,1) * t297 - qJDD(1) * mrSges(2,2) + t312;
t128 = m(2) * t276 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t297 + t131;
t125 = -mrSges(3,1) * t286 + mrSges(3,3) * t242 + t297 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t140 - t322;
t124 = mrSges(3,2) * t286 - mrSges(3,3) * t241 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t297 - pkin(7) * t140 + t127 * t294 - t133 * t291;
t123 = -mrSges(2,2) * g(3) - mrSges(2,3) * t276 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t297 - qJ(2) * t131 + t124 * t288 - t125 * t287;
t122 = Ifges(2,6) * qJDD(1) + t297 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t277 + t287 * t124 + t288 * t125 - pkin(1) * (m(3) * t286 + t140) + qJ(2) * t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t295 * t123 - t292 * t122 - pkin(6) * (t128 * t295 + t129 * t292), t123, t124, t127, t135, t153, -t203 * t239 + t306; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t292 * t123 + t295 * t122 + pkin(6) * (-t128 * t292 + t129 * t295), t122, t125, t133, t134, t152, -t304; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t303, t303, t305, t322, t298, -t300, Ifges(7,5) * t197 + Ifges(7,6) * t260 + Ifges(7,3) * t196 + t240 * t203 - t279 * t205 - t308;];
m_new  = t1;
