% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRP2
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
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:34:44
% EndTime: 2019-05-05 09:35:12
% DurationCPUTime: 14.69s
% Computational Cost: add. (275220->337), mult. (537649->421), div. (0->0), fcn. (382771->12), ass. (0->135)
t301 = sin(pkin(11));
t303 = cos(pkin(11));
t287 = g(1) * t301 - g(2) * t303;
t288 = -g(1) * t303 - g(2) * t301;
t300 = -g(3) + qJDD(1);
t311 = cos(qJ(2));
t304 = cos(pkin(6));
t308 = sin(qJ(2));
t339 = t304 * t308;
t302 = sin(pkin(6));
t340 = t302 * t308;
t259 = t287 * t339 + t288 * t311 + t300 * t340;
t312 = qJD(2) ^ 2;
t254 = -pkin(2) * t312 + qJDD(2) * pkin(8) + t259;
t269 = -t287 * t302 + t300 * t304;
t307 = sin(qJ(3));
t310 = cos(qJ(3));
t229 = -t307 * t254 + t269 * t310;
t333 = qJD(2) * qJD(3);
t331 = t310 * t333;
t284 = qJDD(2) * t307 + t331;
t208 = (-t284 + t331) * pkin(9) + (t307 * t310 * t312 + qJDD(3)) * pkin(3) + t229;
t230 = t254 * t310 + t269 * t307;
t285 = qJDD(2) * t310 - t307 * t333;
t335 = qJD(2) * t307;
t292 = qJD(3) * pkin(3) - pkin(9) * t335;
t299 = t310 ^ 2;
t209 = -pkin(3) * t299 * t312 + pkin(9) * t285 - qJD(3) * t292 + t230;
t306 = sin(qJ(4));
t309 = cos(qJ(4));
t204 = t208 * t306 + t209 * t309;
t277 = (t306 * t310 + t307 * t309) * qJD(2);
t247 = -qJD(4) * t277 - t284 * t306 + t285 * t309;
t276 = (t306 * t307 - t309 * t310) * qJD(2);
t261 = mrSges(5,1) * t276 + mrSges(5,2) * t277;
t298 = qJD(3) + qJD(4);
t268 = mrSges(5,1) * t298 - mrSges(5,3) * t277;
t297 = qJDD(3) + qJDD(4);
t262 = pkin(4) * t276 - pkin(10) * t277;
t296 = t298 ^ 2;
t199 = -pkin(4) * t296 + pkin(10) * t297 - t262 * t276 + t204;
t258 = -t308 * t288 + (t287 * t304 + t300 * t302) * t311;
t319 = -qJDD(2) * pkin(2) - t258;
t220 = -t285 * pkin(3) + t292 * t335 + (-pkin(9) * t299 - pkin(8)) * t312 + t319;
t248 = -qJD(4) * t276 + t284 * t309 + t285 * t306;
t201 = (t276 * t298 - t248) * pkin(10) + (t277 * t298 - t247) * pkin(4) + t220;
t305 = sin(qJ(5));
t345 = cos(qJ(5));
t196 = t199 * t345 + t201 * t305;
t264 = t277 * t345 + t298 * t305;
t218 = qJD(5) * t264 + t248 * t305 - t297 * t345;
t245 = qJDD(5) - t247;
t271 = qJD(5) + t276;
t251 = mrSges(6,1) * t271 - mrSges(6,3) * t264;
t263 = t277 * t305 - t298 * t345;
t233 = pkin(5) * t263 - qJ(6) * t264;
t270 = t271 ^ 2;
t191 = -pkin(5) * t270 + qJ(6) * t245 + 0.2e1 * qJD(6) * t271 - t233 * t263 + t196;
t252 = -mrSges(7,1) * t271 + mrSges(7,2) * t264;
t332 = m(7) * t191 + mrSges(7,3) * t245 + t252 * t271;
t234 = mrSges(7,1) * t263 - mrSges(7,3) * t264;
t336 = -mrSges(6,1) * t263 - mrSges(6,2) * t264 - t234;
t343 = -mrSges(6,3) - mrSges(7,2);
t181 = m(6) * t196 - t245 * mrSges(6,2) + t218 * t343 - t271 * t251 + t263 * t336 + t332;
t195 = -t199 * t305 + t201 * t345;
t219 = -qJD(5) * t263 + t248 * t345 + t297 * t305;
t250 = -mrSges(6,2) * t271 - mrSges(6,3) * t263;
t193 = -pkin(5) * t245 - qJ(6) * t270 + t233 * t264 + qJDD(6) - t195;
t249 = -mrSges(7,2) * t263 + mrSges(7,3) * t271;
t327 = -m(7) * t193 + mrSges(7,1) * t245 + t249 * t271;
t183 = m(6) * t195 + t245 * mrSges(6,1) + t219 * t343 + t271 * t250 + t264 * t336 + t327;
t328 = t181 * t345 - t183 * t305;
t169 = m(5) * t204 - mrSges(5,2) * t297 + mrSges(5,3) * t247 - t261 * t276 - t268 * t298 + t328;
t203 = t309 * t208 - t209 * t306;
t267 = -mrSges(5,2) * t298 - mrSges(5,3) * t276;
t198 = -t297 * pkin(4) - t296 * pkin(10) + t262 * t277 - t203;
t194 = -0.2e1 * qJD(6) * t264 + (t263 * t271 - t219) * qJ(6) + (t264 * t271 + t218) * pkin(5) + t198;
t188 = m(7) * t194 + mrSges(7,1) * t218 - mrSges(7,3) * t219 + t249 * t263 - t252 * t264;
t316 = -m(6) * t198 - mrSges(6,1) * t218 - mrSges(6,2) * t219 - t250 * t263 - t251 * t264 - t188;
t178 = m(5) * t203 + mrSges(5,1) * t297 - mrSges(5,3) * t248 - t261 * t277 + t267 * t298 + t316;
t164 = t169 * t306 + t178 * t309;
t283 = (-mrSges(4,1) * t310 + mrSges(4,2) * t307) * qJD(2);
t334 = qJD(2) * t310;
t290 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t334;
t162 = m(4) * t229 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t284 + qJD(3) * t290 - t283 * t335 + t164;
t289 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t335;
t329 = t169 * t309 - t178 * t306;
t163 = m(4) * t230 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t285 - qJD(3) * t289 + t283 * t334 + t329;
t156 = t162 * t310 + t163 * t307;
t274 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t307 + Ifges(4,2) * t310) * qJD(2);
t275 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t307 + Ifges(4,4) * t310) * qJD(2);
t225 = Ifges(7,1) * t264 + Ifges(7,4) * t271 + Ifges(7,5) * t263;
t226 = Ifges(6,1) * t264 - Ifges(6,4) * t263 + Ifges(6,5) * t271;
t326 = -mrSges(7,1) * t194 + mrSges(7,2) * t191;
t223 = Ifges(7,4) * t264 + Ifges(7,2) * t271 + Ifges(7,6) * t263;
t338 = -Ifges(6,5) * t264 + Ifges(6,6) * t263 - Ifges(6,3) * t271 - t223;
t171 = -mrSges(6,1) * t198 + mrSges(6,3) * t196 - pkin(5) * t188 + (t225 + t226) * t271 + t338 * t264 + (Ifges(6,6) - Ifges(7,6)) * t245 + (Ifges(6,4) - Ifges(7,5)) * t219 + (-Ifges(6,2) - Ifges(7,3)) * t218 + t326;
t224 = Ifges(6,4) * t264 - Ifges(6,2) * t263 + Ifges(6,6) * t271;
t221 = Ifges(7,5) * t264 + Ifges(7,6) * t271 + Ifges(7,3) * t263;
t321 = mrSges(7,2) * t193 - mrSges(7,3) * t194 + Ifges(7,1) * t219 + Ifges(7,4) * t245 + Ifges(7,5) * t218 + t221 * t271;
t173 = mrSges(6,2) * t198 - mrSges(6,3) * t195 + Ifges(6,1) * t219 - Ifges(6,4) * t218 + Ifges(6,5) * t245 - qJ(6) * t188 - t271 * t224 + t263 * t338 + t321;
t256 = Ifges(5,4) * t277 - Ifges(5,2) * t276 + Ifges(5,6) * t298;
t257 = Ifges(5,1) * t277 - Ifges(5,4) * t276 + Ifges(5,5) * t298;
t317 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t248 - Ifges(5,6) * t247 - Ifges(5,3) * t297 - pkin(4) * t316 - pkin(10) * t328 - t171 * t345 - t173 * t305 - t256 * t277 - t276 * t257;
t346 = mrSges(4,1) * t229 - mrSges(4,2) * t230 + Ifges(4,5) * t284 + Ifges(4,6) * t285 + Ifges(4,3) * qJDD(3) + pkin(3) * t164 + (t274 * t307 - t275 * t310) * qJD(2) - t317;
t142 = -mrSges(3,1) * t269 + mrSges(3,3) * t259 + t312 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t156 - t346;
t330 = -t162 * t307 + t163 * t310;
t154 = m(3) * t259 - mrSges(3,1) * t312 - qJDD(2) * mrSges(3,2) + t330;
t253 = -t312 * pkin(8) + t319;
t175 = t181 * t305 + t183 * t345;
t318 = m(5) * t220 - mrSges(5,1) * t247 + mrSges(5,2) * t248 + t267 * t276 + t268 * t277 + t175;
t315 = -m(4) * t253 + mrSges(4,1) * t285 - mrSges(4,2) * t284 - t289 * t335 + t290 * t334 - t318;
t166 = m(3) * t258 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t312 + t315;
t151 = t154 * t311 - t166 * t308;
t348 = pkin(7) * t151 + t142 * t311;
t322 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t219 - Ifges(7,2) * t245 - Ifges(7,6) * t218 - t225 * t263;
t347 = -(-t224 + t221) * t264 + mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t219 - Ifges(6,6) * t218 + Ifges(6,3) * t245 + pkin(5) * (-t219 * mrSges(7,2) - t264 * t234 + t327) + qJ(6) * (-t218 * mrSges(7,2) - t263 * t234 + t332) + t263 * t226 - t322;
t341 = t166 * t311;
t155 = m(3) * t269 + t156;
t147 = t154 * t339 - t155 * t302 + t304 * t341;
t255 = Ifges(5,5) * t277 - Ifges(5,6) * t276 + Ifges(5,3) * t298;
t157 = mrSges(5,2) * t220 - mrSges(5,3) * t203 + Ifges(5,1) * t248 + Ifges(5,4) * t247 + Ifges(5,5) * t297 - pkin(10) * t175 - t171 * t305 + t173 * t345 - t255 * t276 - t256 * t298;
t158 = -mrSges(5,1) * t220 + mrSges(5,3) * t204 + Ifges(5,4) * t248 + Ifges(5,2) * t247 + Ifges(5,6) * t297 - pkin(4) * t175 - t277 * t255 + t298 * t257 - t347;
t273 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t307 + Ifges(4,6) * t310) * qJD(2);
t143 = -mrSges(4,1) * t253 + mrSges(4,3) * t230 + Ifges(4,4) * t284 + Ifges(4,2) * t285 + Ifges(4,6) * qJDD(3) - pkin(3) * t318 + pkin(9) * t329 + qJD(3) * t275 + t306 * t157 + t309 * t158 - t273 * t335;
t148 = mrSges(4,2) * t253 - mrSges(4,3) * t229 + Ifges(4,1) * t284 + Ifges(4,4) * t285 + Ifges(4,5) * qJDD(3) - pkin(9) * t164 - qJD(3) * t274 + t157 * t309 - t158 * t306 + t273 * t334;
t138 = mrSges(3,1) * t258 - mrSges(3,2) * t259 + Ifges(3,3) * qJDD(2) + pkin(2) * t315 + pkin(8) * t330 + t310 * t143 + t307 * t148;
t140 = mrSges(3,2) * t269 - mrSges(3,3) * t258 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t312 - pkin(8) * t156 - t143 * t307 + t148 * t310;
t320 = mrSges(2,1) * t287 - mrSges(2,2) * t288 + pkin(1) * t147 + t304 * t138 + t140 * t340 + t302 * t348;
t149 = m(2) * t288 + t151;
t146 = t304 * t155 + (t154 * t308 + t341) * t302;
t144 = m(2) * t287 + t147;
t136 = mrSges(2,2) * t300 - mrSges(2,3) * t287 + t311 * t140 - t308 * t142 + (-t146 * t302 - t147 * t304) * pkin(7);
t135 = -mrSges(2,1) * t300 + mrSges(2,3) * t288 - pkin(1) * t146 - t302 * t138 + (t140 * t308 + t348) * t304;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t303 * t136 - t301 * t135 - qJ(1) * (t144 * t303 + t149 * t301), t136, t140, t148, t157, t173, -t223 * t263 + t321; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t136 + t303 * t135 + qJ(1) * (-t144 * t301 + t149 * t303), t135, t142, t143, t158, t171, -t264 * t221 - t322; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t320, t320, t138, t346, -t317, t347, Ifges(7,5) * t219 + Ifges(7,6) * t245 + Ifges(7,3) * t218 + t264 * t223 - t271 * t225 - t326;];
m_new  = t1;
