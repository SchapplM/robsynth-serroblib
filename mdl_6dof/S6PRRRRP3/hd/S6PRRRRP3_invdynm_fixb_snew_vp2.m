% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRP3
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
% Datum: 2019-05-05 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:42:41
% EndTime: 2019-05-05 09:43:15
% DurationCPUTime: 16.37s
% Computational Cost: add. (302959->334), mult. (588794->416), div. (0->0), fcn. (419806->12), ass. (0->135)
t307 = sin(qJ(4));
t311 = cos(qJ(4));
t308 = sin(qJ(3));
t336 = qJD(2) * t308;
t287 = t307 * qJD(3) + t311 * t336;
t312 = cos(qJ(3));
t334 = qJD(2) * qJD(3);
t331 = t312 * t334;
t290 = t308 * qJDD(2) + t331;
t258 = -t287 * qJD(4) + t311 * qJDD(3) - t307 * t290;
t286 = t311 * qJD(3) - t307 * t336;
t259 = t286 * qJD(4) + t307 * qJDD(3) + t311 * t290;
t306 = sin(qJ(5));
t310 = cos(qJ(5));
t261 = t310 * t286 - t306 * t287;
t218 = t261 * qJD(5) + t306 * t258 + t310 * t259;
t262 = t306 * t286 + t310 * t287;
t238 = -t261 * mrSges(7,1) + t262 * mrSges(7,2);
t302 = sin(pkin(11));
t304 = cos(pkin(11));
t293 = t302 * g(1) - t304 * g(2);
t294 = -t304 * g(1) - t302 * g(2);
t301 = -g(3) + qJDD(1);
t313 = cos(qJ(2));
t305 = cos(pkin(6));
t309 = sin(qJ(2));
t338 = t305 * t309;
t303 = sin(pkin(6));
t339 = t303 * t309;
t251 = t293 * t338 + t313 * t294 + t301 * t339;
t315 = qJD(2) ^ 2;
t244 = -t315 * pkin(2) + qJDD(2) * pkin(8) + t251;
t268 = -t303 * t293 + t305 * t301;
t237 = t312 * t244 + t308 * t268;
t289 = (-pkin(3) * t312 - pkin(9) * t308) * qJD(2);
t314 = qJD(3) ^ 2;
t335 = t312 * qJD(2);
t222 = -t314 * pkin(3) + qJDD(3) * pkin(9) + t289 * t335 + t237;
t250 = -t309 * t294 + (t293 * t305 + t301 * t303) * t313;
t243 = -qJDD(2) * pkin(2) - t315 * pkin(8) - t250;
t300 = t308 * t334;
t291 = t312 * qJDD(2) - t300;
t226 = (-t290 - t331) * pkin(9) + (-t291 + t300) * pkin(3) + t243;
t199 = -t307 * t222 + t311 * t226;
t283 = qJDD(4) - t291;
t299 = qJD(4) - t335;
t196 = (t286 * t299 - t259) * pkin(10) + (t286 * t287 + t283) * pkin(4) + t199;
t200 = t311 * t222 + t307 * t226;
t267 = t299 * pkin(4) - t287 * pkin(10);
t282 = t286 ^ 2;
t198 = -t282 * pkin(4) + t258 * pkin(10) - t299 * t267 + t200;
t189 = t310 * t196 - t306 * t198;
t279 = qJDD(5) + t283;
t298 = qJD(5) + t299;
t184 = -0.2e1 * qJD(6) * t262 + (t261 * t298 - t218) * qJ(6) + (t261 * t262 + t279) * pkin(5) + t189;
t245 = -t298 * mrSges(7,2) + t261 * mrSges(7,3);
t333 = m(7) * t184 + t279 * mrSges(7,1) + t298 * t245;
t181 = -t218 * mrSges(7,3) - t262 * t238 + t333;
t190 = t306 * t196 + t310 * t198;
t217 = -t262 * qJD(5) + t310 * t258 - t306 * t259;
t230 = Ifges(6,4) * t262 + Ifges(6,2) * t261 + Ifges(6,6) * t298;
t231 = Ifges(7,1) * t262 + Ifges(7,4) * t261 + Ifges(7,5) * t298;
t232 = Ifges(6,1) * t262 + Ifges(6,4) * t261 + Ifges(6,5) * t298;
t247 = t298 * pkin(5) - t262 * qJ(6);
t260 = t261 ^ 2;
t187 = -t260 * pkin(5) + t217 * qJ(6) + 0.2e1 * qJD(6) * t261 - t298 * t247 + t190;
t229 = Ifges(7,4) * t262 + Ifges(7,2) * t261 + Ifges(7,6) * t298;
t323 = -mrSges(7,1) * t184 + mrSges(7,2) * t187 - Ifges(7,5) * t218 - Ifges(7,6) * t217 - Ifges(7,3) * t279 - t262 * t229;
t346 = mrSges(6,1) * t189 - mrSges(6,2) * t190 + Ifges(6,5) * t218 + Ifges(6,6) * t217 + Ifges(6,3) * t279 + pkin(5) * t181 + t262 * t230 - t323 + (-t232 - t231) * t261;
t239 = -t261 * mrSges(6,1) + t262 * mrSges(6,2);
t246 = -t298 * mrSges(6,2) + t261 * mrSges(6,3);
t173 = m(6) * t189 + t279 * mrSges(6,1) + t298 * t246 + (-t238 - t239) * t262 + (-mrSges(6,3) - mrSges(7,3)) * t218 + t333;
t248 = t298 * mrSges(7,1) - t262 * mrSges(7,3);
t249 = t298 * mrSges(6,1) - t262 * mrSges(6,3);
t332 = m(7) * t187 + t217 * mrSges(7,3) + t261 * t238;
t176 = m(6) * t190 + t217 * mrSges(6,3) + t261 * t239 + (-t248 - t249) * t298 + (-mrSges(6,2) - mrSges(7,2)) * t279 + t332;
t171 = t310 * t173 + t306 * t176;
t253 = Ifges(5,4) * t287 + Ifges(5,2) * t286 + Ifges(5,6) * t299;
t254 = Ifges(5,1) * t287 + Ifges(5,4) * t286 + Ifges(5,5) * t299;
t345 = mrSges(5,1) * t199 - mrSges(5,2) * t200 + Ifges(5,5) * t259 + Ifges(5,6) * t258 + Ifges(5,3) * t283 + pkin(4) * t171 + t287 * t253 - t286 * t254 + t346;
t263 = -t286 * mrSges(5,1) + t287 * mrSges(5,2);
t265 = -t299 * mrSges(5,2) + t286 * mrSges(5,3);
t168 = m(5) * t199 + t283 * mrSges(5,1) - t259 * mrSges(5,3) - t287 * t263 + t299 * t265 + t171;
t266 = t299 * mrSges(5,1) - t287 * mrSges(5,3);
t328 = -t306 * t173 + t310 * t176;
t169 = m(5) * t200 - t283 * mrSges(5,2) + t258 * mrSges(5,3) + t286 * t263 - t299 * t266 + t328;
t165 = -t307 * t168 + t311 * t169;
t288 = (-mrSges(4,1) * t312 + mrSges(4,2) * t308) * qJD(2);
t295 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t336;
t163 = m(4) * t237 - qJDD(3) * mrSges(4,2) + t291 * mrSges(4,3) - qJD(3) * t295 + t288 * t335 + t165;
t236 = -t308 * t244 + t312 * t268;
t221 = -qJDD(3) * pkin(3) - t314 * pkin(9) + t289 * t336 - t236;
t201 = -t258 * pkin(4) - t282 * pkin(10) + t287 * t267 + t221;
t193 = -t217 * pkin(5) - t260 * qJ(6) + t262 * t247 + qJDD(6) + t201;
t327 = m(7) * t193 - t217 * mrSges(7,1) + t218 * mrSges(7,2) - t261 * t245 + t262 * t248;
t320 = m(6) * t201 - t217 * mrSges(6,1) + t218 * mrSges(6,2) - t261 * t246 + t262 * t249 + t327;
t179 = -m(5) * t221 + t258 * mrSges(5,1) - t259 * mrSges(5,2) + t286 * t265 - t287 * t266 - t320;
t296 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t335;
t178 = m(4) * t236 + qJDD(3) * mrSges(4,1) - t290 * mrSges(4,3) + qJD(3) * t296 - t288 * t336 + t179;
t158 = t308 * t163 + t312 * t178;
t227 = Ifges(7,5) * t262 + Ifges(7,6) * t261 + Ifges(7,3) * t298;
t228 = Ifges(6,5) * t262 + Ifges(6,6) * t261 + Ifges(6,3) * t298;
t324 = -mrSges(7,1) * t193 + mrSges(7,3) * t187 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t279 + t298 * t231;
t166 = Ifges(6,4) * t218 + Ifges(6,2) * t217 + Ifges(6,6) * t279 + t298 * t232 - mrSges(6,1) * t201 + mrSges(6,3) * t190 - pkin(5) * t327 + qJ(6) * (-t279 * mrSges(7,2) - t298 * t248 + t332) + (-t228 - t227) * t262 + t324;
t322 = mrSges(7,2) * t193 - mrSges(7,3) * t184 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t279 + t261 * t227;
t170 = mrSges(6,2) * t201 - mrSges(6,3) * t189 + Ifges(6,1) * t218 + Ifges(6,4) * t217 + Ifges(6,5) * t279 - qJ(6) * t181 + t261 * t228 + (-t229 - t230) * t298 + t322;
t252 = Ifges(5,5) * t287 + Ifges(5,6) * t286 + Ifges(5,3) * t299;
t152 = -mrSges(5,1) * t221 + mrSges(5,3) * t200 + Ifges(5,4) * t259 + Ifges(5,2) * t258 + Ifges(5,6) * t283 - pkin(4) * t320 + pkin(10) * t328 + t310 * t166 + t306 * t170 - t287 * t252 + t299 * t254;
t153 = mrSges(5,2) * t221 - mrSges(5,3) * t199 + Ifges(5,1) * t259 + Ifges(5,4) * t258 + Ifges(5,5) * t283 - pkin(10) * t171 - t306 * t166 + t310 * t170 + t286 * t252 - t299 * t253;
t277 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t308 + Ifges(4,2) * t312) * qJD(2);
t278 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t308 + Ifges(4,4) * t312) * qJD(2);
t343 = mrSges(4,1) * t236 - mrSges(4,2) * t237 + Ifges(4,5) * t290 + Ifges(4,6) * t291 + Ifges(4,3) * qJDD(3) + pkin(3) * t179 + pkin(9) * t165 + t311 * t152 + t307 * t153 + (t308 * t277 - t312 * t278) * qJD(2);
t142 = -mrSges(3,1) * t268 + mrSges(3,3) * t251 + t315 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t158 - t343;
t329 = t312 * t163 - t308 * t178;
t156 = m(3) * t251 - t315 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t329;
t164 = t311 * t168 + t307 * t169;
t319 = -m(4) * t243 + t291 * mrSges(4,1) - t290 * mrSges(4,2) - t295 * t336 + t296 * t335 - t164;
t160 = m(3) * t250 + qJDD(2) * mrSges(3,1) - t315 * mrSges(3,2) + t319;
t150 = t313 * t156 - t309 * t160;
t344 = pkin(7) * t150 + t142 * t313;
t340 = t160 * t313;
t157 = m(3) * t268 + t158;
t147 = t156 * t338 - t303 * t157 + t305 * t340;
t276 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t308 + Ifges(4,6) * t312) * qJD(2);
t143 = mrSges(4,2) * t243 - mrSges(4,3) * t236 + Ifges(4,1) * t290 + Ifges(4,4) * t291 + Ifges(4,5) * qJDD(3) - pkin(9) * t164 - qJD(3) * t277 - t307 * t152 + t311 * t153 + t276 * t335;
t151 = -mrSges(4,1) * t243 + mrSges(4,3) * t237 + Ifges(4,4) * t290 + Ifges(4,2) * t291 + Ifges(4,6) * qJDD(3) - pkin(3) * t164 + qJD(3) * t278 - t276 * t336 - t345;
t138 = mrSges(3,1) * t250 - mrSges(3,2) * t251 + Ifges(3,3) * qJDD(2) + pkin(2) * t319 + pkin(8) * t329 + t308 * t143 + t312 * t151;
t140 = mrSges(3,2) * t268 - mrSges(3,3) * t250 + Ifges(3,5) * qJDD(2) - t315 * Ifges(3,6) - pkin(8) * t158 + t312 * t143 - t308 * t151;
t321 = mrSges(2,1) * t293 - mrSges(2,2) * t294 + pkin(1) * t147 + t305 * t138 + t140 * t339 + t303 * t344;
t148 = m(2) * t294 + t150;
t146 = t305 * t157 + (t156 * t309 + t340) * t303;
t144 = m(2) * t293 + t147;
t136 = mrSges(2,2) * t301 - mrSges(2,3) * t293 + t313 * t140 - t309 * t142 + (-t146 * t303 - t147 * t305) * pkin(7);
t135 = -mrSges(2,1) * t301 + mrSges(2,3) * t294 - pkin(1) * t146 - t303 * t138 + (t140 * t309 + t344) * t305;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t136 - t302 * t135 - qJ(1) * (t304 * t144 + t302 * t148), t136, t140, t143, t153, t170, -t298 * t229 + t322; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t302 * t136 + t304 * t135 + qJ(1) * (-t302 * t144 + t304 * t148), t135, t142, t151, t152, t166, -t262 * t227 + t324; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t321, t321, t138, t343, t345, t346, -t261 * t231 - t323;];
m_new  = t1;
