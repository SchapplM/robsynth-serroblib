% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP1
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
% Datum: 2019-05-06 01:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:07:39
% EndTime: 2019-05-06 01:07:57
% DurationCPUTime: 9.33s
% Computational Cost: add. (174209->338), mult. (339958->413), div. (0->0), fcn. (224592->10), ass. (0->128)
t297 = sin(qJ(4));
t298 = sin(qJ(3));
t300 = cos(qJ(4));
t301 = cos(qJ(3));
t265 = (t297 * t298 - t300 * t301) * qJD(1);
t299 = sin(qJ(1));
t302 = cos(qJ(1));
t279 = t299 * g(1) - t302 * g(2);
t270 = qJDD(1) * pkin(1) + t279;
t280 = -t302 * g(1) - t299 * g(2);
t303 = qJD(1) ^ 2;
t272 = -t303 * pkin(1) + t280;
t294 = sin(pkin(10));
t295 = cos(pkin(10));
t252 = t294 * t270 + t295 * t272;
t247 = -t303 * pkin(2) + qJDD(1) * pkin(7) + t252;
t293 = -g(3) + qJDD(2);
t233 = -t298 * t247 + t301 * t293;
t325 = qJD(1) * qJD(3);
t323 = t301 * t325;
t273 = t298 * qJDD(1) + t323;
t207 = (-t273 + t323) * pkin(8) + (t298 * t301 * t303 + qJDD(3)) * pkin(3) + t233;
t234 = t301 * t247 + t298 * t293;
t274 = t301 * qJDD(1) - t298 * t325;
t327 = qJD(1) * t298;
t278 = qJD(3) * pkin(3) - pkin(8) * t327;
t292 = t301 ^ 2;
t208 = -t292 * t303 * pkin(3) + t274 * pkin(8) - qJD(3) * t278 + t234;
t193 = t297 * t207 + t300 * t208;
t266 = (t297 * t301 + t298 * t300) * qJD(1);
t249 = t265 * pkin(4) - t266 * pkin(9);
t289 = qJD(3) + qJD(4);
t287 = t289 ^ 2;
t288 = qJDD(3) + qJDD(4);
t188 = -t287 * pkin(4) + t288 * pkin(9) - t265 * t249 + t193;
t251 = t295 * t270 - t294 * t272;
t314 = -qJDD(1) * pkin(2) - t251;
t217 = -t274 * pkin(3) + t278 * t327 + (-pkin(8) * t292 - pkin(7)) * t303 + t314;
t235 = -t266 * qJD(4) - t297 * t273 + t300 * t274;
t236 = -t265 * qJD(4) + t300 * t273 + t297 * t274;
t190 = (t265 * t289 - t236) * pkin(9) + (t266 * t289 - t235) * pkin(4) + t217;
t296 = sin(qJ(5));
t332 = cos(qJ(5));
t184 = -t296 * t188 + t190 * t332;
t185 = t188 * t332 + t296 * t190;
t254 = t266 * t332 + t296 * t289;
t202 = t254 * qJD(5) + t296 * t236 - t288 * t332;
t253 = t296 * t266 - t289 * t332;
t203 = -t253 * qJD(5) + t236 * t332 + t296 * t288;
t258 = qJD(5) + t265;
t209 = Ifges(7,5) * t254 + Ifges(7,6) * t258 + Ifges(7,3) * t253;
t212 = Ifges(6,4) * t254 - Ifges(6,2) * t253 + Ifges(6,6) * t258;
t214 = Ifges(6,1) * t254 - Ifges(6,4) * t253 + Ifges(6,5) * t258;
t221 = t253 * mrSges(7,1) - t254 * mrSges(7,3);
t232 = qJDD(5) - t235;
t220 = t253 * pkin(5) - t254 * qJ(6);
t257 = t258 ^ 2;
t180 = -t257 * pkin(5) + t232 * qJ(6) + 0.2e1 * qJD(6) * t258 - t253 * t220 + t185;
t182 = -t232 * pkin(5) - t257 * qJ(6) + t254 * t220 + qJDD(6) - t184;
t213 = Ifges(7,1) * t254 + Ifges(7,4) * t258 + Ifges(7,5) * t253;
t313 = mrSges(7,1) * t182 - mrSges(7,3) * t180 - Ifges(7,4) * t203 - Ifges(7,2) * t232 - Ifges(7,6) * t202 - t253 * t213;
t237 = -t253 * mrSges(7,2) + t258 * mrSges(7,3);
t318 = -m(7) * t182 + t232 * mrSges(7,1) + t258 * t237;
t240 = -t258 * mrSges(7,1) + t254 * mrSges(7,2);
t324 = m(7) * t180 + t232 * mrSges(7,3) + t258 * t240;
t334 = -(-t212 + t209) * t254 + mrSges(6,1) * t184 - mrSges(6,2) * t185 + Ifges(6,5) * t203 - Ifges(6,6) * t202 + Ifges(6,3) * t232 + pkin(5) * (-t203 * mrSges(7,2) - t254 * t221 + t318) + qJ(6) * (-t202 * mrSges(7,2) - t253 * t221 + t324) + t253 * t214 - t313;
t248 = t265 * mrSges(5,1) + t266 * mrSges(5,2);
t256 = t289 * mrSges(5,1) - t266 * mrSges(5,3);
t239 = t258 * mrSges(6,1) - t254 * mrSges(6,3);
t328 = -t253 * mrSges(6,1) - t254 * mrSges(6,2) - t221;
t331 = -mrSges(6,3) - mrSges(7,2);
t170 = m(6) * t185 - t232 * mrSges(6,2) + t202 * t331 - t258 * t239 + t253 * t328 + t324;
t238 = -t258 * mrSges(6,2) - t253 * mrSges(6,3);
t172 = m(6) * t184 + t232 * mrSges(6,1) + t203 * t331 + t258 * t238 + t254 * t328 + t318;
t319 = t170 * t332 - t296 * t172;
t158 = m(5) * t193 - t288 * mrSges(5,2) + t235 * mrSges(5,3) - t265 * t248 - t289 * t256 + t319;
t192 = t300 * t207 - t297 * t208;
t255 = -t289 * mrSges(5,2) - t265 * mrSges(5,3);
t187 = -t288 * pkin(4) - t287 * pkin(9) + t266 * t249 - t192;
t183 = -0.2e1 * qJD(6) * t254 + (t253 * t258 - t203) * qJ(6) + (t254 * t258 + t202) * pkin(5) + t187;
t177 = m(7) * t183 + t202 * mrSges(7,1) - t203 * mrSges(7,3) + t253 * t237 - t254 * t240;
t307 = -m(6) * t187 - t202 * mrSges(6,1) - t203 * mrSges(6,2) - t253 * t238 - t254 * t239 - t177;
t167 = m(5) * t192 + t288 * mrSges(5,1) - t236 * mrSges(5,3) - t266 * t248 + t289 * t255 + t307;
t152 = t297 * t158 + t300 * t167;
t263 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t298 + Ifges(4,2) * t301) * qJD(1);
t264 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t298 + Ifges(4,4) * t301) * qJD(1);
t317 = -mrSges(7,1) * t183 + mrSges(7,2) * t180;
t211 = Ifges(7,4) * t254 + Ifges(7,2) * t258 + Ifges(7,6) * t253;
t330 = -Ifges(6,5) * t254 + Ifges(6,6) * t253 - Ifges(6,3) * t258 - t211;
t160 = -mrSges(6,1) * t187 + mrSges(6,3) * t185 - pkin(5) * t177 + (t213 + t214) * t258 + t330 * t254 + (Ifges(6,6) - Ifges(7,6)) * t232 + (Ifges(6,4) - Ifges(7,5)) * t203 + (-Ifges(6,2) - Ifges(7,3)) * t202 + t317;
t312 = mrSges(7,2) * t182 - mrSges(7,3) * t183 + Ifges(7,1) * t203 + Ifges(7,4) * t232 + Ifges(7,5) * t202 + t258 * t209;
t162 = mrSges(6,2) * t187 - mrSges(6,3) * t184 + Ifges(6,1) * t203 - Ifges(6,4) * t202 + Ifges(6,5) * t232 - qJ(6) * t177 - t258 * t212 + t253 * t330 + t312;
t243 = Ifges(5,4) * t266 - Ifges(5,2) * t265 + Ifges(5,6) * t289;
t244 = Ifges(5,1) * t266 - Ifges(5,4) * t265 + Ifges(5,5) * t289;
t308 = -mrSges(5,1) * t192 + mrSges(5,2) * t193 - Ifges(5,5) * t236 - Ifges(5,6) * t235 - Ifges(5,3) * t288 - pkin(4) * t307 - pkin(9) * t319 - t160 * t332 - t296 * t162 - t266 * t243 - t265 * t244;
t333 = mrSges(4,1) * t233 - mrSges(4,2) * t234 + Ifges(4,5) * t273 + Ifges(4,6) * t274 + Ifges(4,3) * qJDD(3) + pkin(3) * t152 + (t298 * t263 - t301 * t264) * qJD(1) - t308;
t271 = (-mrSges(4,1) * t301 + mrSges(4,2) * t298) * qJD(1);
t326 = qJD(1) * t301;
t277 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t326;
t150 = m(4) * t233 + qJDD(3) * mrSges(4,1) - t273 * mrSges(4,3) + qJD(3) * t277 - t271 * t327 + t152;
t276 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t327;
t320 = t300 * t158 - t297 * t167;
t151 = m(4) * t234 - qJDD(3) * mrSges(4,2) + t274 * mrSges(4,3) - qJD(3) * t276 + t271 * t326 + t320;
t321 = -t298 * t150 + t301 * t151;
t142 = m(3) * t252 - t303 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t321;
t246 = -t303 * pkin(7) + t314;
t164 = t296 * t170 + t172 * t332;
t310 = m(5) * t217 - t235 * mrSges(5,1) + t236 * mrSges(5,2) + t265 * t255 + t266 * t256 + t164;
t306 = -m(4) * t246 + t274 * mrSges(4,1) - t273 * mrSges(4,2) - t276 * t327 + t277 * t326 - t310;
t154 = m(3) * t251 + qJDD(1) * mrSges(3,1) - t303 * mrSges(3,2) + t306;
t139 = t294 * t142 + t295 * t154;
t144 = t301 * t150 + t298 * t151;
t322 = t295 * t142 - t294 * t154;
t242 = Ifges(5,5) * t266 - Ifges(5,6) * t265 + Ifges(5,3) * t289;
t145 = mrSges(5,2) * t217 - mrSges(5,3) * t192 + Ifges(5,1) * t236 + Ifges(5,4) * t235 + Ifges(5,5) * t288 - pkin(9) * t164 - t296 * t160 + t162 * t332 - t265 * t242 - t289 * t243;
t146 = -mrSges(5,1) * t217 + mrSges(5,3) * t193 + Ifges(5,4) * t236 + Ifges(5,2) * t235 + Ifges(5,6) * t288 - pkin(4) * t164 - t266 * t242 + t289 * t244 - t334;
t262 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t298 + Ifges(4,6) * t301) * qJD(1);
t133 = -mrSges(4,1) * t246 + mrSges(4,3) * t234 + Ifges(4,4) * t273 + Ifges(4,2) * t274 + Ifges(4,6) * qJDD(3) - pkin(3) * t310 + pkin(8) * t320 + qJD(3) * t264 + t297 * t145 + t300 * t146 - t262 * t327;
t135 = mrSges(4,2) * t246 - mrSges(4,3) * t233 + Ifges(4,1) * t273 + Ifges(4,4) * t274 + Ifges(4,5) * qJDD(3) - pkin(8) * t152 - qJD(3) * t263 + t300 * t145 - t297 * t146 + t262 * t326;
t311 = mrSges(3,1) * t251 - mrSges(3,2) * t252 + Ifges(3,3) * qJDD(1) + pkin(2) * t306 + pkin(7) * t321 + t301 * t133 + t298 * t135;
t309 = mrSges(2,1) * t279 - mrSges(2,2) * t280 + Ifges(2,3) * qJDD(1) + pkin(1) * t139 + t311;
t137 = m(2) * t280 - t303 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t322;
t136 = m(2) * t279 + qJDD(1) * mrSges(2,1) - t303 * mrSges(2,2) + t139;
t131 = -mrSges(3,1) * t293 + mrSges(3,3) * t252 + t303 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t144 - t333;
t130 = mrSges(3,2) * t293 - mrSges(3,3) * t251 + Ifges(3,5) * qJDD(1) - t303 * Ifges(3,6) - pkin(7) * t144 - t298 * t133 + t301 * t135;
t129 = -mrSges(2,2) * g(3) - mrSges(2,3) * t279 + Ifges(2,5) * qJDD(1) - t303 * Ifges(2,6) - qJ(2) * t139 + t295 * t130 - t294 * t131;
t128 = Ifges(2,6) * qJDD(1) + t303 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t280 + t294 * t130 + t295 * t131 - pkin(1) * (m(3) * t293 + t144) + qJ(2) * t322;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t302 * t129 - t299 * t128 - pkin(6) * (t302 * t136 + t299 * t137), t129, t130, t135, t145, t162, -t253 * t211 + t312; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t299 * t129 + t302 * t128 + pkin(6) * (-t299 * t136 + t302 * t137), t128, t131, t133, t146, t160, -t254 * t209 - t313; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t309, t309, t311, t333, -t308, t334, Ifges(7,5) * t203 + Ifges(7,6) * t232 + Ifges(7,3) * t202 + t254 * t211 - t258 * t213 - t317;];
m_new  = t1;
