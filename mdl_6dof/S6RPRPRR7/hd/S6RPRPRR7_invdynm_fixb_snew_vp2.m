% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:12:25
% EndTime: 2019-05-05 19:12:46
% DurationCPUTime: 13.36s
% Computational Cost: add. (240035->341), mult. (534256->422), div. (0->0), fcn. (372577->10), ass. (0->136)
t302 = sin(qJ(1));
t306 = cos(qJ(1));
t280 = -t306 * g(1) - t302 * g(2);
t322 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t280;
t335 = -pkin(1) - pkin(7);
t334 = mrSges(2,1) - mrSges(3,2);
t333 = Ifges(2,5) - Ifges(3,4);
t332 = (-Ifges(2,6) + Ifges(3,5));
t279 = t302 * g(1) - t306 * g(2);
t307 = qJD(1) ^ 2;
t321 = -t307 * qJ(2) + qJDD(2) - t279;
t252 = qJDD(1) * t335 + t321;
t301 = sin(qJ(3));
t305 = cos(qJ(3));
t243 = t301 * g(3) + t305 * t252;
t329 = qJD(1) * qJD(3);
t327 = t301 * t329;
t274 = t305 * qJDD(1) - t327;
t219 = (-t274 - t327) * qJ(4) + (-t301 * t305 * t307 + qJDD(3)) * pkin(3) + t243;
t244 = -t305 * g(3) + t301 * t252;
t273 = -t301 * qJDD(1) - t305 * t329;
t330 = qJD(1) * t305;
t277 = qJD(3) * pkin(3) - qJ(4) * t330;
t294 = t301 ^ 2;
t220 = -t294 * t307 * pkin(3) + t273 * qJ(4) - qJD(3) * t277 + t244;
t297 = sin(pkin(10));
t298 = cos(pkin(10));
t262 = (-t297 * t301 + t298 * t305) * qJD(1);
t196 = -0.2e1 * qJD(4) * t262 + t298 * t219 - t297 * t220;
t241 = t297 * t273 + t298 * t274;
t261 = (-t297 * t305 - t298 * t301) * qJD(1);
t185 = (qJD(3) * t261 - t241) * pkin(8) + (t261 * t262 + qJDD(3)) * pkin(4) + t196;
t197 = 0.2e1 * qJD(4) * t261 + t297 * t219 + t298 * t220;
t240 = t298 * t273 - t297 * t274;
t251 = qJD(3) * pkin(4) - t262 * pkin(8);
t260 = t261 ^ 2;
t187 = -t260 * pkin(4) + t240 * pkin(8) - qJD(3) * t251 + t197;
t300 = sin(qJ(5));
t304 = cos(qJ(5));
t182 = t300 * t185 + t304 * t187;
t233 = t300 * t261 + t304 * t262;
t205 = -t233 * qJD(5) + t304 * t240 - t300 * t241;
t232 = t304 * t261 - t300 * t262;
t214 = -t232 * mrSges(6,1) + t233 * mrSges(6,2);
t287 = qJD(3) + qJD(5);
t227 = t287 * mrSges(6,1) - t233 * mrSges(6,3);
t286 = qJDD(3) + qJDD(5);
t215 = -t232 * pkin(5) - t233 * pkin(9);
t285 = t287 ^ 2;
t179 = -t285 * pkin(5) + t286 * pkin(9) + t232 * t215 + t182;
t222 = -t273 * pkin(3) + qJDD(4) + t277 * t330 + (-qJ(4) * t294 + t335) * t307 + t322;
t199 = -t240 * pkin(4) - t260 * pkin(8) + t262 * t251 + t222;
t206 = t232 * qJD(5) + t300 * t240 + t304 * t241;
t183 = (-t232 * t287 - t206) * pkin(9) + (t233 * t287 - t205) * pkin(5) + t199;
t299 = sin(qJ(6));
t303 = cos(qJ(6));
t176 = -t299 * t179 + t303 * t183;
t224 = -t299 * t233 + t303 * t287;
t190 = t224 * qJD(6) + t303 * t206 + t299 * t286;
t204 = qJDD(6) - t205;
t225 = t303 * t233 + t299 * t287;
t207 = -t224 * mrSges(7,1) + t225 * mrSges(7,2);
t228 = qJD(6) - t232;
t208 = -t228 * mrSges(7,2) + t224 * mrSges(7,3);
t172 = m(7) * t176 + t204 * mrSges(7,1) - t190 * mrSges(7,3) - t225 * t207 + t228 * t208;
t177 = t303 * t179 + t299 * t183;
t189 = -t225 * qJD(6) - t299 * t206 + t303 * t286;
t209 = t228 * mrSges(7,1) - t225 * mrSges(7,3);
t173 = m(7) * t177 - t204 * mrSges(7,2) + t189 * mrSges(7,3) + t224 * t207 - t228 * t209;
t324 = -t299 * t172 + t303 * t173;
t159 = m(6) * t182 - t286 * mrSges(6,2) + t205 * mrSges(6,3) + t232 * t214 - t287 * t227 + t324;
t181 = t304 * t185 - t300 * t187;
t226 = -t287 * mrSges(6,2) + t232 * mrSges(6,3);
t178 = -t286 * pkin(5) - t285 * pkin(9) + t233 * t215 - t181;
t319 = -m(7) * t178 + t189 * mrSges(7,1) - t190 * mrSges(7,2) + t224 * t208 - t225 * t209;
t168 = m(6) * t181 + t286 * mrSges(6,1) - t206 * mrSges(6,3) - t233 * t214 + t287 * t226 + t319;
t153 = t300 * t159 + t304 * t168;
t236 = -t261 * mrSges(5,1) + t262 * mrSges(5,2);
t249 = -qJD(3) * mrSges(5,2) + t261 * mrSges(5,3);
t150 = m(5) * t196 + qJDD(3) * mrSges(5,1) - t241 * mrSges(5,3) + qJD(3) * t249 - t262 * t236 + t153;
t250 = qJD(3) * mrSges(5,1) - t262 * mrSges(5,3);
t325 = t304 * t159 - t300 * t168;
t151 = m(5) * t197 - qJDD(3) * mrSges(5,2) + t240 * mrSges(5,3) - qJD(3) * t250 + t261 * t236 + t325;
t144 = t298 * t150 + t297 * t151;
t161 = t303 * t172 + t299 * t173;
t331 = qJD(1) * t301;
t272 = (mrSges(4,1) * t301 + mrSges(4,2) * t305) * qJD(1);
t276 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t331;
t141 = m(4) * t243 + qJDD(3) * mrSges(4,1) - t274 * mrSges(4,3) + qJD(3) * t276 - t272 * t330 + t144;
t278 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t326 = -t297 * t150 + t298 * t151;
t142 = m(4) * t244 - qJDD(3) * mrSges(4,2) + t273 * mrSges(4,3) - qJD(3) * t278 - t272 * t331 + t326;
t138 = -t301 * t141 + t305 * t142;
t137 = t305 * t141 + t301 * t142;
t259 = -qJDD(1) * pkin(1) + t321;
t320 = -m(3) * t259 + (t307 * mrSges(3,3)) - t137;
t318 = m(6) * t199 - t205 * mrSges(6,1) + t206 * mrSges(6,2) - t232 * t226 + t233 * t227 + t161;
t192 = Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t228;
t194 = Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t228;
t165 = -mrSges(7,1) * t178 + mrSges(7,3) * t177 + Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t204 - t225 * t192 + t228 * t194;
t193 = Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t228;
t166 = mrSges(7,2) * t178 - mrSges(7,3) * t176 + Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t204 + t224 * t192 - t228 * t193;
t210 = Ifges(6,5) * t233 + Ifges(6,6) * t232 + Ifges(6,3) * t287;
t211 = Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t287;
t145 = mrSges(6,2) * t199 - mrSges(6,3) * t181 + Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t286 - pkin(9) * t161 - t299 * t165 + t303 * t166 + t232 * t210 - t287 * t211;
t212 = Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t287;
t313 = mrSges(7,1) * t176 - mrSges(7,2) * t177 + Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t204 + t225 * t193 - t224 * t194;
t146 = -mrSges(6,1) * t199 + mrSges(6,3) * t182 + Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t286 - pkin(5) * t161 - t233 * t210 + t287 * t212 - t313;
t229 = Ifges(5,5) * t262 + Ifges(5,6) * t261 + (Ifges(5,3) * qJD(3));
t231 = Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * qJD(3);
t133 = -mrSges(5,1) * t222 + mrSges(5,3) * t197 + Ifges(5,4) * t241 + Ifges(5,2) * t240 + Ifges(5,6) * qJDD(3) - pkin(4) * t318 + pkin(8) * t325 + qJD(3) * t231 + t300 * t145 + t304 * t146 - t262 * t229;
t230 = Ifges(5,4) * t262 + Ifges(5,2) * t261 + Ifges(5,6) * qJD(3);
t139 = mrSges(5,2) * t222 - mrSges(5,3) * t196 + Ifges(5,1) * t241 + Ifges(5,4) * t240 + Ifges(5,5) * qJDD(3) - pkin(8) * t153 - qJD(3) * t230 + t304 * t145 - t300 * t146 + t261 * t229;
t248 = t307 * t335 + t322;
t263 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t305 - Ifges(4,6) * t301) * qJD(1);
t265 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t305 - Ifges(4,4) * t301) * qJD(1);
t312 = m(5) * t222 - t240 * mrSges(5,1) + t241 * mrSges(5,2) - t261 * t249 + t262 * t250 + t318;
t130 = -mrSges(4,1) * t248 + mrSges(4,3) * t244 + Ifges(4,4) * t274 + Ifges(4,2) * t273 + Ifges(4,6) * qJDD(3) - pkin(3) * t312 + qJ(4) * t326 + qJD(3) * t265 + t298 * t133 + t297 * t139 - t263 * t330;
t264 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t305 - Ifges(4,2) * t301) * qJD(1);
t132 = mrSges(4,2) * t248 - mrSges(4,3) * t243 + Ifges(4,1) * t274 + Ifges(4,4) * t273 + Ifges(4,5) * qJDD(3) - qJ(4) * t144 - qJD(3) * t264 - t297 * t133 + t298 * t139 - t263 * t331;
t255 = t307 * pkin(1) - t322;
t317 = mrSges(3,2) * t259 - mrSges(3,3) * t255 + Ifges(3,1) * qJDD(1) - pkin(7) * t137 - t301 * t130 + t305 * t132;
t156 = -m(4) * t248 + t273 * mrSges(4,1) - t274 * mrSges(4,2) - t276 * t331 - t278 * t330 - t312;
t316 = -mrSges(3,1) * t255 - pkin(2) * t156 - pkin(7) * t138 - t305 * t130 - t301 * t132;
t315 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t206 - Ifges(6,6) * t205 - Ifges(6,3) * t286 - pkin(5) * t319 - pkin(9) * t324 - t303 * t165 - t299 * t166 - t233 * t211 + t232 * t212;
t310 = -m(3) * t255 + t307 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t156;
t314 = -mrSges(2,2) * t280 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t320) + qJ(2) * t310 + mrSges(2,1) * t279 + Ifges(2,3) * qJDD(1) + t317;
t311 = -mrSges(5,1) * t196 + mrSges(5,2) * t197 - Ifges(5,5) * t241 - Ifges(5,6) * t240 - Ifges(5,3) * qJDD(3) - pkin(4) * t153 - t262 * t230 + t261 * t231 + t315;
t309 = mrSges(4,1) * t243 - mrSges(4,2) * t244 + Ifges(4,5) * t274 + Ifges(4,6) * t273 + Ifges(4,3) * qJDD(3) + pkin(3) * t144 + t264 * t330 + t265 * t331 - t311;
t308 = -mrSges(3,1) * t259 - pkin(2) * t137 - t309;
t154 = m(2) * t280 - t307 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t310;
t136 = -m(3) * g(3) + t138;
t134 = m(2) * t279 - t307 * mrSges(2,2) + qJDD(1) * t334 + t320;
t129 = -t308 + (t332 * t307) + t333 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t279 - qJ(2) * t136;
t128 = mrSges(2,3) * t280 - pkin(1) * t136 + g(3) * t334 - qJDD(1) * t332 + t307 * t333 + t316;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t306 * t129 - t302 * t128 - pkin(6) * (t306 * t134 + t302 * t154), t129, t317, t132, t139, t145, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t302 * t129 + t306 * t128 + pkin(6) * (-t302 * t134 + t306 * t154), t128, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t307 * Ifges(3,5)) + t308, t130, t133, t146, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t314, t314, mrSges(3,2) * g(3) + t307 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t316, t309, -t311, -t315, t313;];
m_new  = t1;
