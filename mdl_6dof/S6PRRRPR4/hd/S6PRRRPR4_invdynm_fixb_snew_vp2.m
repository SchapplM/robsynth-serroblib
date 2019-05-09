% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:41:43
% EndTime: 2019-05-05 07:42:31
% DurationCPUTime: 36.52s
% Computational Cost: add. (701212->339), mult. (1412491->436), div. (0->0), fcn. (1031503->14), ass. (0->144)
t299 = sin(pkin(11));
t302 = cos(pkin(11));
t289 = g(1) * t299 - g(2) * t302;
t290 = -g(1) * t302 - g(2) * t299;
t297 = -g(3) + qJDD(1);
t311 = cos(qJ(2));
t303 = cos(pkin(6));
t307 = sin(qJ(2));
t331 = t303 * t307;
t300 = sin(pkin(6));
t332 = t300 * t307;
t248 = t289 * t331 + t290 * t311 + t297 * t332;
t313 = qJD(2) ^ 2;
t243 = -pkin(2) * t313 + qJDD(2) * pkin(8) + t248;
t267 = -t289 * t300 + t297 * t303;
t306 = sin(qJ(3));
t310 = cos(qJ(3));
t238 = t243 * t310 + t267 * t306;
t285 = (-pkin(3) * t310 - pkin(9) * t306) * qJD(2);
t312 = qJD(3) ^ 2;
t329 = t310 * qJD(2);
t220 = -pkin(3) * t312 + qJDD(3) * pkin(9) + t285 * t329 + t238;
t247 = -t307 * t290 + (t289 * t303 + t297 * t300) * t311;
t242 = -qJDD(2) * pkin(2) - t313 * pkin(8) - t247;
t328 = qJD(2) * qJD(3);
t327 = t310 * t328;
t286 = qJDD(2) * t306 + t327;
t296 = t306 * t328;
t287 = qJDD(2) * t310 - t296;
t225 = (-t286 - t327) * pkin(9) + (-t287 + t296) * pkin(3) + t242;
t305 = sin(qJ(4));
t309 = cos(qJ(4));
t208 = -t305 * t220 + t225 * t309;
t330 = qJD(2) * t306;
t282 = qJD(3) * t309 - t305 * t330;
t258 = qJD(4) * t282 + qJDD(3) * t305 + t286 * t309;
t279 = qJDD(4) - t287;
t283 = qJD(3) * t305 + t309 * t330;
t295 = qJD(4) - t329;
t198 = (t282 * t295 - t258) * qJ(5) + (t282 * t283 + t279) * pkin(4) + t208;
t209 = t220 * t309 + t225 * t305;
t257 = -qJD(4) * t283 + qJDD(3) * t309 - t286 * t305;
t265 = pkin(4) * t295 - qJ(5) * t283;
t278 = t282 ^ 2;
t200 = -pkin(4) * t278 + qJ(5) * t257 - t265 * t295 + t209;
t298 = sin(pkin(12));
t301 = cos(pkin(12));
t261 = t282 * t298 + t283 * t301;
t192 = -0.2e1 * qJD(5) * t261 + t198 * t301 - t298 * t200;
t233 = t257 * t298 + t258 * t301;
t260 = t282 * t301 - t283 * t298;
t189 = (t260 * t295 - t233) * pkin(10) + (t260 * t261 + t279) * pkin(5) + t192;
t193 = 0.2e1 * qJD(5) * t260 + t198 * t298 + t200 * t301;
t232 = t257 * t301 - t258 * t298;
t246 = pkin(5) * t295 - pkin(10) * t261;
t259 = t260 ^ 2;
t190 = -pkin(5) * t259 + pkin(10) * t232 - t246 * t295 + t193;
t304 = sin(qJ(6));
t308 = cos(qJ(6));
t187 = t189 * t308 - t190 * t304;
t235 = t260 * t308 - t261 * t304;
t206 = qJD(6) * t235 + t232 * t304 + t233 * t308;
t236 = t260 * t304 + t261 * t308;
t216 = -mrSges(7,1) * t235 + mrSges(7,2) * t236;
t294 = qJD(6) + t295;
t223 = -mrSges(7,2) * t294 + mrSges(7,3) * t235;
t275 = qJDD(6) + t279;
t181 = m(7) * t187 + mrSges(7,1) * t275 - mrSges(7,3) * t206 - t216 * t236 + t223 * t294;
t188 = t189 * t304 + t190 * t308;
t205 = -qJD(6) * t236 + t232 * t308 - t233 * t304;
t224 = mrSges(7,1) * t294 - mrSges(7,3) * t236;
t182 = m(7) * t188 - mrSges(7,2) * t275 + mrSges(7,3) * t205 + t216 * t235 - t224 * t294;
t175 = t181 * t308 + t182 * t304;
t239 = -mrSges(6,1) * t260 + mrSges(6,2) * t261;
t244 = -mrSges(6,2) * t295 + mrSges(6,3) * t260;
t172 = m(6) * t192 + mrSges(6,1) * t279 - mrSges(6,3) * t233 - t239 * t261 + t244 * t295 + t175;
t245 = mrSges(6,1) * t295 - mrSges(6,3) * t261;
t324 = -t181 * t304 + t182 * t308;
t173 = m(6) * t193 - mrSges(6,2) * t279 + mrSges(6,3) * t232 + t239 * t260 - t245 * t295 + t324;
t168 = t172 * t301 + t173 * t298;
t262 = -mrSges(5,1) * t282 + mrSges(5,2) * t283;
t264 = -mrSges(5,2) * t295 + mrSges(5,3) * t282;
t166 = m(5) * t208 + mrSges(5,1) * t279 - mrSges(5,3) * t258 - t262 * t283 + t264 * t295 + t168;
t266 = mrSges(5,1) * t295 - mrSges(5,3) * t283;
t325 = -t172 * t298 + t173 * t301;
t167 = m(5) * t209 - mrSges(5,2) * t279 + mrSges(5,3) * t257 + t262 * t282 - t266 * t295 + t325;
t162 = -t166 * t305 + t167 * t309;
t284 = (-mrSges(4,1) * t310 + mrSges(4,2) * t306) * qJD(2);
t291 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t160 = m(4) * t238 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t287 - qJD(3) * t291 + t284 * t329 + t162;
t237 = -t243 * t306 + t267 * t310;
t219 = -qJDD(3) * pkin(3) - pkin(9) * t312 + t285 * t330 - t237;
t210 = -pkin(4) * t257 - qJ(5) * t278 + t265 * t283 + qJDD(5) + t219;
t195 = -pkin(5) * t232 - pkin(10) * t259 + t246 * t261 + t210;
t323 = m(7) * t195 - mrSges(7,1) * t205 + mrSges(7,2) * t206 - t223 * t235 + t224 * t236;
t318 = m(6) * t210 - mrSges(6,1) * t232 + mrSges(6,2) * t233 - t244 * t260 + t245 * t261 + t323;
t185 = -m(5) * t219 + mrSges(5,1) * t257 - mrSges(5,2) * t258 + t264 * t282 - t266 * t283 - t318;
t292 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t329;
t184 = m(4) * t237 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t286 + qJD(3) * t292 - t284 * t330 + t185;
t155 = t160 * t306 + t184 * t310;
t211 = Ifges(7,5) * t236 + Ifges(7,6) * t235 + Ifges(7,3) * t294;
t213 = Ifges(7,1) * t236 + Ifges(7,4) * t235 + Ifges(7,5) * t294;
t176 = -mrSges(7,1) * t195 + mrSges(7,3) * t188 + Ifges(7,4) * t206 + Ifges(7,2) * t205 + Ifges(7,6) * t275 - t211 * t236 + t213 * t294;
t212 = Ifges(7,4) * t236 + Ifges(7,2) * t235 + Ifges(7,6) * t294;
t177 = mrSges(7,2) * t195 - mrSges(7,3) * t187 + Ifges(7,1) * t206 + Ifges(7,4) * t205 + Ifges(7,5) * t275 + t211 * t235 - t212 * t294;
t229 = Ifges(6,5) * t261 + Ifges(6,6) * t260 + Ifges(6,3) * t295;
t231 = Ifges(6,1) * t261 + Ifges(6,4) * t260 + Ifges(6,5) * t295;
t163 = -mrSges(6,1) * t210 + mrSges(6,3) * t193 + Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t279 - pkin(5) * t323 + pkin(10) * t324 + t308 * t176 + t304 * t177 - t261 * t229 + t295 * t231;
t230 = Ifges(6,4) * t261 + Ifges(6,2) * t260 + Ifges(6,6) * t295;
t164 = mrSges(6,2) * t210 - mrSges(6,3) * t192 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t279 - pkin(10) * t175 - t176 * t304 + t177 * t308 + t229 * t260 - t230 * t295;
t249 = Ifges(5,5) * t283 + Ifges(5,6) * t282 + Ifges(5,3) * t295;
t251 = Ifges(5,1) * t283 + Ifges(5,4) * t282 + Ifges(5,5) * t295;
t149 = -mrSges(5,1) * t219 + mrSges(5,3) * t209 + Ifges(5,4) * t258 + Ifges(5,2) * t257 + Ifges(5,6) * t279 - pkin(4) * t318 + qJ(5) * t325 + t301 * t163 + t298 * t164 - t283 * t249 + t295 * t251;
t250 = Ifges(5,4) * t283 + Ifges(5,2) * t282 + Ifges(5,6) * t295;
t150 = mrSges(5,2) * t219 - mrSges(5,3) * t208 + Ifges(5,1) * t258 + Ifges(5,4) * t257 + Ifges(5,5) * t279 - qJ(5) * t168 - t163 * t298 + t164 * t301 + t249 * t282 - t250 * t295;
t272 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t306 + Ifges(4,2) * t310) * qJD(2);
t273 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t306 + Ifges(4,4) * t310) * qJD(2);
t336 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t286 + Ifges(4,6) * t287 + Ifges(4,3) * qJDD(3) + pkin(3) * t185 + pkin(9) * t162 + t309 * t149 + t305 * t150 + (t272 * t306 - t273 * t310) * qJD(2);
t139 = -mrSges(3,1) * t267 + mrSges(3,3) * t248 + t313 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t155 - t336;
t326 = t160 * t310 - t184 * t306;
t153 = m(3) * t248 - mrSges(3,1) * t313 - qJDD(2) * mrSges(3,2) + t326;
t161 = t166 * t309 + t167 * t305;
t317 = -m(4) * t242 + mrSges(4,1) * t287 - mrSges(4,2) * t286 - t291 * t330 + t292 * t329 - t161;
t157 = m(3) * t247 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t313 + t317;
t147 = t153 * t311 - t157 * t307;
t337 = pkin(7) * t147 + t139 * t311;
t333 = t157 * t311;
t154 = m(3) * t267 + t155;
t144 = t153 * t331 - t154 * t300 + t303 * t333;
t271 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t306 + Ifges(4,6) * t310) * qJD(2);
t140 = mrSges(4,2) * t242 - mrSges(4,3) * t237 + Ifges(4,1) * t286 + Ifges(4,4) * t287 + Ifges(4,5) * qJDD(3) - pkin(9) * t161 - qJD(3) * t272 - t149 * t305 + t150 * t309 + t271 * t329;
t319 = -mrSges(7,1) * t187 + mrSges(7,2) * t188 - Ifges(7,5) * t206 - Ifges(7,6) * t205 - Ifges(7,3) * t275 - t212 * t236 + t213 * t235;
t316 = -mrSges(6,1) * t192 + mrSges(6,2) * t193 - Ifges(6,5) * t233 - Ifges(6,6) * t232 - Ifges(6,3) * t279 - pkin(5) * t175 - t230 * t261 + t231 * t260 + t319;
t314 = mrSges(5,1) * t208 - mrSges(5,2) * t209 + Ifges(5,5) * t258 + Ifges(5,6) * t257 + Ifges(5,3) * t279 + pkin(4) * t168 + t250 * t283 - t251 * t282 - t316;
t148 = -mrSges(4,1) * t242 + mrSges(4,3) * t238 + Ifges(4,4) * t286 + Ifges(4,2) * t287 + Ifges(4,6) * qJDD(3) - pkin(3) * t161 + qJD(3) * t273 - t271 * t330 - t314;
t135 = mrSges(3,1) * t247 - mrSges(3,2) * t248 + Ifges(3,3) * qJDD(2) + pkin(2) * t317 + pkin(8) * t326 + t306 * t140 + t310 * t148;
t137 = mrSges(3,2) * t267 - mrSges(3,3) * t247 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t313 - pkin(8) * t155 + t140 * t310 - t148 * t306;
t320 = mrSges(2,1) * t289 - mrSges(2,2) * t290 + pkin(1) * t144 + t303 * t135 + t137 * t332 + t300 * t337;
t145 = m(2) * t290 + t147;
t143 = t303 * t154 + (t153 * t307 + t333) * t300;
t141 = m(2) * t289 + t144;
t133 = mrSges(2,2) * t297 - mrSges(2,3) * t289 + t311 * t137 - t307 * t139 + (-t143 * t300 - t144 * t303) * pkin(7);
t132 = -mrSges(2,1) * t297 + mrSges(2,3) * t290 - pkin(1) * t143 - t300 * t135 + (t137 * t307 + t337) * t303;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t302 * t133 - t299 * t132 - qJ(1) * (t141 * t302 + t145 * t299), t133, t137, t140, t150, t164, t177; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t299 * t133 + t302 * t132 + qJ(1) * (-t141 * t299 + t145 * t302), t132, t139, t148, t149, t163, t176; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t320, t320, t135, t336, t314, -t316, -t319;];
m_new  = t1;
