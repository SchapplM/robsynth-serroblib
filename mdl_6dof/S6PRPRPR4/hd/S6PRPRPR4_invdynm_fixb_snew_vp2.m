% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:45:47
% EndTime: 2019-05-04 22:46:20
% DurationCPUTime: 27.35s
% Computational Cost: add. (471259->321), mult. (1083349->414), div. (0->0), fcn. (823240->14), ass. (0->148)
t305 = sin(pkin(10));
t309 = cos(pkin(10));
t289 = g(1) * t305 - g(2) * t309;
t290 = -g(1) * t309 - g(2) * t305;
t302 = -g(3) + qJDD(1);
t315 = cos(qJ(2));
t310 = cos(pkin(6));
t313 = sin(qJ(2));
t347 = t310 * t313;
t306 = sin(pkin(6));
t348 = t306 * t313;
t257 = t289 * t347 + t315 * t290 + t302 * t348;
t317 = qJD(2) ^ 2;
t247 = -pkin(2) * t317 + qJDD(2) * qJ(3) + t257;
t304 = sin(pkin(11));
t277 = -t289 * t306 + t302 * t310;
t308 = cos(pkin(11));
t342 = qJD(2) * qJD(3);
t346 = t308 * t277 - 0.2e1 * t304 * t342;
t354 = pkin(3) * t308;
t228 = (-pkin(8) * qJDD(2) + t317 * t354 - t247) * t304 + t346;
t231 = t304 * t277 + (t247 + 0.2e1 * t342) * t308;
t340 = qJDD(2) * t308;
t300 = t308 ^ 2;
t349 = t300 * t317;
t229 = -pkin(3) * t349 + pkin(8) * t340 + t231;
t312 = sin(qJ(4));
t355 = cos(qJ(4));
t211 = t312 * t228 + t355 * t229;
t339 = t308 * t355;
t345 = qJD(2) * t304;
t279 = -qJD(2) * t339 + t312 * t345;
t328 = t355 * t304 + t308 * t312;
t280 = t328 * qJD(2);
t261 = mrSges(5,1) * t279 + mrSges(5,2) * t280;
t341 = qJDD(2) * t304;
t344 = qJD(4) * t280;
t266 = -qJDD(2) * t339 + t312 * t341 + t344;
t276 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t280;
t260 = pkin(4) * t279 - qJ(5) * t280;
t316 = qJD(4) ^ 2;
t207 = -pkin(4) * t316 + qJDD(4) * qJ(5) - t260 * t279 + t211;
t299 = t304 ^ 2;
t256 = -t313 * t290 + (t289 * t310 + t302 * t306) * t315;
t325 = qJDD(3) - t256;
t238 = (-pkin(2) - t354) * qJDD(2) + (-qJ(3) + (-t299 - t300) * pkin(8)) * t317 + t325;
t343 = t279 * qJD(4);
t267 = t328 * qJDD(2) - t343;
t214 = (-t267 + t343) * qJ(5) + (t266 + t344) * pkin(4) + t238;
t303 = sin(pkin(12));
t307 = cos(pkin(12));
t272 = qJD(4) * t303 + t280 * t307;
t202 = -0.2e1 * qJD(5) * t272 - t303 * t207 + t307 * t214;
t252 = qJDD(4) * t303 + t267 * t307;
t271 = qJD(4) * t307 - t280 * t303;
t200 = (t271 * t279 - t252) * pkin(9) + (t271 * t272 + t266) * pkin(5) + t202;
t203 = 0.2e1 * qJD(5) * t271 + t307 * t207 + t303 * t214;
t250 = pkin(5) * t279 - pkin(9) * t272;
t251 = qJDD(4) * t307 - t267 * t303;
t270 = t271 ^ 2;
t201 = -pkin(5) * t270 + pkin(9) * t251 - t250 * t279 + t203;
t311 = sin(qJ(6));
t314 = cos(qJ(6));
t198 = t200 * t314 - t201 * t311;
t239 = t271 * t314 - t272 * t311;
t219 = qJD(6) * t239 + t251 * t311 + t252 * t314;
t240 = t271 * t311 + t272 * t314;
t224 = -mrSges(7,1) * t239 + mrSges(7,2) * t240;
t278 = qJD(6) + t279;
t232 = -mrSges(7,2) * t278 + mrSges(7,3) * t239;
t265 = qJDD(6) + t266;
t193 = m(7) * t198 + mrSges(7,1) * t265 - mrSges(7,3) * t219 - t224 * t240 + t232 * t278;
t199 = t200 * t311 + t201 * t314;
t218 = -qJD(6) * t240 + t251 * t314 - t252 * t311;
t233 = mrSges(7,1) * t278 - mrSges(7,3) * t240;
t194 = m(7) * t199 - mrSges(7,2) * t265 + mrSges(7,3) * t218 + t224 * t239 - t233 * t278;
t185 = t314 * t193 + t311 * t194;
t242 = -mrSges(6,1) * t271 + mrSges(6,2) * t272;
t248 = -mrSges(6,2) * t279 + mrSges(6,3) * t271;
t183 = m(6) * t202 + mrSges(6,1) * t266 - mrSges(6,3) * t252 - t242 * t272 + t248 * t279 + t185;
t249 = mrSges(6,1) * t279 - mrSges(6,3) * t272;
t335 = -t193 * t311 + t314 * t194;
t184 = m(6) * t203 - mrSges(6,2) * t266 + mrSges(6,3) * t251 + t242 * t271 - t249 * t279 + t335;
t336 = -t183 * t303 + t307 * t184;
t176 = m(5) * t211 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t266 - qJD(4) * t276 - t261 * t279 + t336;
t210 = t355 * t228 - t312 * t229;
t275 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t279;
t206 = -qJDD(4) * pkin(4) - t316 * qJ(5) + t280 * t260 + qJDD(5) - t210;
t204 = -t251 * pkin(5) - t270 * pkin(9) + t272 * t250 + t206;
t326 = m(7) * t204 - t218 * mrSges(7,1) + mrSges(7,2) * t219 - t239 * t232 + t233 * t240;
t320 = -m(6) * t206 + t251 * mrSges(6,1) - mrSges(6,2) * t252 + t271 * t248 - t249 * t272 - t326;
t189 = m(5) * t210 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t267 + qJD(4) * t275 - t261 * t280 + t320;
t167 = t312 * t176 + t355 * t189;
t230 = -t247 * t304 + t346;
t352 = mrSges(4,2) * t304;
t329 = mrSges(4,3) * qJDD(2) + t317 * (-mrSges(4,1) * t308 + t352);
t165 = m(4) * t230 - t329 * t304 + t167;
t337 = t355 * t176 - t312 * t189;
t166 = m(4) * t231 + t329 * t308 + t337;
t160 = t308 * t165 + t304 * t166;
t332 = Ifges(4,5) * t304 + Ifges(4,6) * t308;
t220 = Ifges(7,5) * t240 + Ifges(7,6) * t239 + Ifges(7,3) * t278;
t222 = Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t278;
t186 = -mrSges(7,1) * t204 + mrSges(7,3) * t199 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t265 - t220 * t240 + t222 * t278;
t221 = Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t278;
t187 = mrSges(7,2) * t204 - mrSges(7,3) * t198 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t265 + t220 * t239 - t221 * t278;
t234 = Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t279;
t236 = Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t279;
t169 = -mrSges(6,1) * t206 + mrSges(6,3) * t203 + Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t266 - pkin(5) * t326 + pkin(9) * t335 + t314 * t186 + t311 * t187 - t272 * t234 + t279 * t236;
t235 = Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t279;
t171 = mrSges(6,2) * t206 - mrSges(6,3) * t202 + Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t266 - pkin(9) * t185 - t186 * t311 + t187 * t314 + t234 * t271 - t235 * t279;
t254 = Ifges(5,4) * t280 - Ifges(5,2) * t279 + Ifges(5,6) * qJD(4);
t255 = Ifges(5,1) * t280 - Ifges(5,4) * t279 + Ifges(5,5) * qJD(4);
t322 = -mrSges(5,1) * t210 + mrSges(5,2) * t211 - Ifges(5,5) * t267 + Ifges(5,6) * t266 - Ifges(5,3) * qJDD(4) - pkin(4) * t320 - qJ(5) * t336 - t307 * t169 - t303 * t171 - t280 * t254 - t279 * t255;
t333 = Ifges(4,4) * t304 + Ifges(4,2) * t308;
t334 = Ifges(4,1) * t304 + Ifges(4,4) * t308;
t357 = qJD(2) * t308;
t356 = -mrSges(4,1) * t230 + mrSges(4,2) * t231 - pkin(3) * t167 - (t333 * t345 - t334 * t357) * qJD(2) + t322;
t146 = (Ifges(3,6) - t332) * qJDD(2) + t317 * Ifges(3,5) - mrSges(3,1) * t277 + mrSges(3,3) * t257 - pkin(2) * t160 + t356;
t338 = -t165 * t304 + t308 * t166;
t158 = m(3) * t257 - mrSges(3,1) * t317 - qJDD(2) * mrSges(3,2) + t338;
t245 = -qJDD(2) * pkin(2) - t317 * qJ(3) + t325;
t178 = t307 * t183 + t303 * t184;
t323 = m(5) * t238 + t266 * mrSges(5,1) + t267 * mrSges(5,2) + t279 * t275 + t280 * t276 + t178;
t321 = -m(4) * t245 + mrSges(4,1) * t340 - t323 + (t299 * t317 + t349) * mrSges(4,3);
t173 = (mrSges(3,1) - t352) * qJDD(2) + t321 - t317 * mrSges(3,2) + m(3) * t256;
t154 = t315 * t158 - t173 * t313;
t358 = pkin(7) * t154 + t146 * t315;
t350 = t173 * t315;
t159 = m(3) * t277 + t160;
t151 = t158 * t347 - t159 * t306 + t310 * t350;
t253 = Ifges(5,5) * t280 - Ifges(5,6) * t279 + Ifges(5,3) * qJD(4);
t155 = mrSges(5,2) * t238 - mrSges(5,3) * t210 + Ifges(5,1) * t267 - Ifges(5,4) * t266 + Ifges(5,5) * qJDD(4) - qJ(5) * t178 - qJD(4) * t254 - t169 * t303 + t171 * t307 - t253 * t279;
t324 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t219 - Ifges(7,6) * t218 - Ifges(7,3) * t265 - t240 * t221 + t239 * t222;
t318 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t252 - Ifges(6,6) * t251 - pkin(5) * t185 - t272 * t235 + t271 * t236 + t324;
t161 = Ifges(5,6) * qJDD(4) + (-Ifges(5,2) - Ifges(6,3)) * t266 + t318 - t280 * t253 + Ifges(5,4) * t267 + qJD(4) * t255 - mrSges(5,1) * t238 + mrSges(5,3) * t211 - pkin(4) * t178;
t285 = t332 * qJD(2);
t144 = -mrSges(4,1) * t245 + mrSges(4,3) * t231 - pkin(3) * t323 + pkin(8) * t337 + t333 * qJDD(2) + t312 * t155 + t355 * t161 - t285 * t345;
t147 = mrSges(4,2) * t245 - mrSges(4,3) * t230 - pkin(8) * t167 + t334 * qJDD(2) + t355 * t155 - t312 * t161 + t285 * t357;
t141 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t256 - mrSges(3,2) * t257 + t304 * t147 + t308 * t144 + pkin(2) * (-mrSges(4,2) * t341 + t321) + qJ(3) * t338;
t143 = mrSges(3,2) * t277 - mrSges(3,3) * t256 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t317 - qJ(3) * t160 - t144 * t304 + t147 * t308;
t327 = mrSges(2,1) * t289 - mrSges(2,2) * t290 + pkin(1) * t151 + t310 * t141 + t143 * t348 + t358 * t306;
t152 = m(2) * t290 + t154;
t150 = t310 * t159 + (t158 * t313 + t350) * t306;
t148 = m(2) * t289 + t151;
t139 = mrSges(2,2) * t302 - mrSges(2,3) * t289 + t315 * t143 - t313 * t146 + (-t150 * t306 - t151 * t310) * pkin(7);
t138 = -mrSges(2,1) * t302 + mrSges(2,3) * t290 - pkin(1) * t150 - t306 * t141 + (t143 * t313 + t358) * t310;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t139 - t305 * t138 - qJ(1) * (t148 * t309 + t152 * t305), t139, t143, t147, t155, t171, t187; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t139 + t309 * t138 + qJ(1) * (-t148 * t305 + t152 * t309), t138, t146, t144, t161, t169, t186; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t327, t327, t141, t332 * qJDD(2) - t356, -t322, Ifges(6,3) * t266 - t318, -t324;];
m_new  = t1;
