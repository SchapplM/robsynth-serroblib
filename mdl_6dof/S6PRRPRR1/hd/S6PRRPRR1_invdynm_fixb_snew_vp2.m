% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:14:28
% EndTime: 2019-05-05 04:15:12
% DurationCPUTime: 32.98s
% Computational Cost: add. (590586->340), mult. (1296652->441), div. (0->0), fcn. (965933->14), ass. (0->144)
t306 = sin(pkin(11));
t309 = cos(pkin(11));
t292 = t306 * g(1) - t309 * g(2);
t293 = -t309 * g(1) - t306 * g(2);
t304 = -g(3) + qJDD(1);
t318 = cos(qJ(2));
t310 = cos(pkin(6));
t314 = sin(qJ(2));
t340 = t310 * t314;
t307 = sin(pkin(6));
t341 = t307 * t314;
t259 = t292 * t340 + t318 * t293 + t304 * t341;
t319 = qJD(2) ^ 2;
t251 = -t319 * pkin(2) + qJDD(2) * pkin(8) + t259;
t273 = -t307 * t292 + t310 * t304;
t313 = sin(qJ(3));
t317 = cos(qJ(3));
t242 = -t313 * t251 + t317 * t273;
t337 = qJD(2) * qJD(3);
t336 = t317 * t337;
t289 = t313 * qJDD(2) + t336;
t229 = (-t289 + t336) * qJ(4) + (t313 * t317 * t319 + qJDD(3)) * pkin(3) + t242;
t243 = t317 * t251 + t313 * t273;
t290 = t317 * qJDD(2) - t313 * t337;
t339 = qJD(2) * t313;
t294 = qJD(3) * pkin(3) - qJ(4) * t339;
t303 = t317 ^ 2;
t231 = -t303 * t319 * pkin(3) + t290 * qJ(4) - qJD(3) * t294 + t243;
t305 = sin(pkin(12));
t308 = cos(pkin(12));
t278 = (t305 * t317 + t308 * t313) * qJD(2);
t207 = -0.2e1 * qJD(4) * t278 + t308 * t229 - t305 * t231;
t267 = t308 * t289 + t305 * t290;
t277 = (-t305 * t313 + t308 * t317) * qJD(2);
t203 = (qJD(3) * t277 - t267) * pkin(9) + (t277 * t278 + qJDD(3)) * pkin(4) + t207;
t208 = 0.2e1 * qJD(4) * t277 + t305 * t229 + t308 * t231;
t266 = -t305 * t289 + t308 * t290;
t272 = qJD(3) * pkin(4) - t278 * pkin(9);
t276 = t277 ^ 2;
t205 = -t276 * pkin(4) + t266 * pkin(9) - qJD(3) * t272 + t208;
t312 = sin(qJ(5));
t316 = cos(qJ(5));
t200 = t312 * t203 + t316 * t205;
t257 = t312 * t277 + t316 * t278;
t224 = -t257 * qJD(5) + t316 * t266 - t312 * t267;
t256 = t316 * t277 - t312 * t278;
t239 = -t256 * mrSges(6,1) + t257 * mrSges(6,2);
t301 = qJD(3) + qJD(5);
t249 = t301 * mrSges(6,1) - t257 * mrSges(6,3);
t300 = qJDD(3) + qJDD(5);
t240 = -t256 * pkin(5) - t257 * pkin(10);
t299 = t301 ^ 2;
t197 = -t299 * pkin(5) + t300 * pkin(10) + t256 * t240 + t200;
t258 = -t314 * t293 + (t292 * t310 + t304 * t307) * t318;
t326 = -qJDD(2) * pkin(2) - t258;
t241 = -t290 * pkin(3) + qJDD(4) + t294 * t339 + (-qJ(4) * t303 - pkin(8)) * t319 + t326;
t213 = -t266 * pkin(4) - t276 * pkin(9) + t278 * t272 + t241;
t225 = t256 * qJD(5) + t312 * t266 + t316 * t267;
t201 = (-t256 * t301 - t225) * pkin(10) + (t257 * t301 - t224) * pkin(5) + t213;
t311 = sin(qJ(6));
t315 = cos(qJ(6));
t194 = -t311 * t197 + t315 * t201;
t245 = -t311 * t257 + t315 * t301;
t211 = t245 * qJD(6) + t315 * t225 + t311 * t300;
t223 = qJDD(6) - t224;
t246 = t315 * t257 + t311 * t301;
t230 = -t245 * mrSges(7,1) + t246 * mrSges(7,2);
t252 = qJD(6) - t256;
t232 = -t252 * mrSges(7,2) + t245 * mrSges(7,3);
t190 = m(7) * t194 + t223 * mrSges(7,1) - t211 * mrSges(7,3) - t246 * t230 + t252 * t232;
t195 = t315 * t197 + t311 * t201;
t210 = -t246 * qJD(6) - t311 * t225 + t315 * t300;
t233 = t252 * mrSges(7,1) - t246 * mrSges(7,3);
t191 = m(7) * t195 - t223 * mrSges(7,2) + t210 * mrSges(7,3) + t245 * t230 - t252 * t233;
t332 = -t311 * t190 + t315 * t191;
t177 = m(6) * t200 - t300 * mrSges(6,2) + t224 * mrSges(6,3) + t256 * t239 - t301 * t249 + t332;
t199 = t316 * t203 - t312 * t205;
t248 = -t301 * mrSges(6,2) + t256 * mrSges(6,3);
t196 = -t300 * pkin(5) - t299 * pkin(10) + t257 * t240 - t199;
t327 = -m(7) * t196 + t210 * mrSges(7,1) - t211 * mrSges(7,2) + t245 * t232 - t246 * t233;
t186 = m(6) * t199 + t300 * mrSges(6,1) - t225 * mrSges(6,3) - t257 * t239 + t301 * t248 + t327;
t172 = t312 * t177 + t316 * t186;
t262 = -t277 * mrSges(5,1) + t278 * mrSges(5,2);
t270 = -qJD(3) * mrSges(5,2) + t277 * mrSges(5,3);
t169 = m(5) * t207 + qJDD(3) * mrSges(5,1) - t267 * mrSges(5,3) + qJD(3) * t270 - t278 * t262 + t172;
t271 = qJD(3) * mrSges(5,1) - t278 * mrSges(5,3);
t333 = t316 * t177 - t312 * t186;
t170 = m(5) * t208 - qJDD(3) * mrSges(5,2) + t266 * mrSges(5,3) - qJD(3) * t271 + t277 * t262 + t333;
t163 = t308 * t169 + t305 * t170;
t288 = (-mrSges(4,1) * t317 + mrSges(4,2) * t313) * qJD(2);
t338 = qJD(2) * t317;
t296 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t338;
t161 = m(4) * t242 + qJDD(3) * mrSges(4,1) - t289 * mrSges(4,3) + qJD(3) * t296 - t288 * t339 + t163;
t295 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t334 = -t305 * t169 + t308 * t170;
t162 = m(4) * t243 - qJDD(3) * mrSges(4,2) + t290 * mrSges(4,3) - qJD(3) * t295 + t288 * t338 + t334;
t156 = t317 * t161 + t313 * t162;
t281 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t313 + Ifges(4,2) * t317) * qJD(2);
t282 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t313 + Ifges(4,4) * t317) * qJD(2);
t254 = Ifges(5,4) * t278 + Ifges(5,2) * t277 + Ifges(5,6) * qJD(3);
t255 = Ifges(5,1) * t278 + Ifges(5,4) * t277 + Ifges(5,5) * qJD(3);
t214 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t252;
t216 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t252;
t183 = -mrSges(7,1) * t196 + mrSges(7,3) * t195 + Ifges(7,4) * t211 + Ifges(7,2) * t210 + Ifges(7,6) * t223 - t246 * t214 + t252 * t216;
t215 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t252;
t184 = mrSges(7,2) * t196 - mrSges(7,3) * t194 + Ifges(7,1) * t211 + Ifges(7,4) * t210 + Ifges(7,5) * t223 + t245 * t214 - t252 * t215;
t235 = Ifges(6,4) * t257 + Ifges(6,2) * t256 + Ifges(6,6) * t301;
t236 = Ifges(6,1) * t257 + Ifges(6,4) * t256 + Ifges(6,5) * t301;
t325 = -mrSges(6,1) * t199 + mrSges(6,2) * t200 - Ifges(6,5) * t225 - Ifges(6,6) * t224 - Ifges(6,3) * t300 - pkin(5) * t327 - pkin(10) * t332 - t315 * t183 - t311 * t184 - t257 * t235 + t256 * t236;
t322 = -mrSges(5,1) * t207 + mrSges(5,2) * t208 - Ifges(5,5) * t267 - Ifges(5,6) * t266 - Ifges(5,3) * qJDD(3) - pkin(4) * t172 - t278 * t254 + t277 * t255 + t325;
t345 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t289 + Ifges(4,6) * t290 + Ifges(4,3) * qJDD(3) + pkin(3) * t163 + (t313 * t281 - t317 * t282) * qJD(2) - t322;
t143 = -mrSges(3,1) * t273 + mrSges(3,3) * t259 + t319 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t156 - t345;
t335 = -t313 * t161 + t317 * t162;
t154 = m(3) * t259 - t319 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t335;
t250 = -t319 * pkin(8) + t326;
t179 = t315 * t190 + t311 * t191;
t329 = m(6) * t213 - t224 * mrSges(6,1) + t225 * mrSges(6,2) - t256 * t248 + t257 * t249 + t179;
t324 = m(5) * t241 - t266 * mrSges(5,1) + t267 * mrSges(5,2) - t277 * t270 + t278 * t271 + t329;
t321 = -m(4) * t250 + t290 * mrSges(4,1) - t289 * mrSges(4,2) - t295 * t339 + t296 * t338 - t324;
t174 = m(3) * t258 + qJDD(2) * mrSges(3,1) - t319 * mrSges(3,2) + t321;
t150 = t318 * t154 - t314 * t174;
t346 = pkin(7) * t150 + t143 * t318;
t342 = t174 * t318;
t155 = m(3) * t273 + t156;
t147 = t154 * t340 - t307 * t155 + t310 * t342;
t234 = Ifges(6,5) * t257 + Ifges(6,6) * t256 + Ifges(6,3) * t301;
t164 = mrSges(6,2) * t213 - mrSges(6,3) * t199 + Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t300 - pkin(10) * t179 - t311 * t183 + t315 * t184 + t256 * t234 - t301 * t235;
t323 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t211 + Ifges(7,6) * t210 + Ifges(7,3) * t223 + t246 * t215 - t245 * t216;
t165 = -mrSges(6,1) * t213 + mrSges(6,3) * t200 + Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t300 - pkin(5) * t179 - t257 * t234 + t301 * t236 - t323;
t253 = Ifges(5,5) * t278 + Ifges(5,6) * t277 + Ifges(5,3) * qJD(3);
t151 = -mrSges(5,1) * t241 + mrSges(5,3) * t208 + Ifges(5,4) * t267 + Ifges(5,2) * t266 + Ifges(5,6) * qJDD(3) - pkin(4) * t329 + pkin(9) * t333 + qJD(3) * t255 + t312 * t164 + t316 * t165 - t278 * t253;
t157 = mrSges(5,2) * t241 - mrSges(5,3) * t207 + Ifges(5,1) * t267 + Ifges(5,4) * t266 + Ifges(5,5) * qJDD(3) - pkin(9) * t172 - qJD(3) * t254 + t316 * t164 - t312 * t165 + t277 * t253;
t280 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t313 + Ifges(4,6) * t317) * qJD(2);
t140 = -mrSges(4,1) * t250 + mrSges(4,3) * t243 + Ifges(4,4) * t289 + Ifges(4,2) * t290 + Ifges(4,6) * qJDD(3) - pkin(3) * t324 + qJ(4) * t334 + qJD(3) * t282 + t308 * t151 + t305 * t157 - t280 * t339;
t141 = mrSges(4,2) * t250 - mrSges(4,3) * t242 + Ifges(4,1) * t289 + Ifges(4,4) * t290 + Ifges(4,5) * qJDD(3) - qJ(4) * t163 - qJD(3) * t281 - t305 * t151 + t308 * t157 + t280 * t338;
t137 = mrSges(3,1) * t258 - mrSges(3,2) * t259 + Ifges(3,3) * qJDD(2) + pkin(2) * t321 + pkin(8) * t335 + t317 * t140 + t313 * t141;
t139 = mrSges(3,2) * t273 - mrSges(3,3) * t258 + Ifges(3,5) * qJDD(2) - t319 * Ifges(3,6) - pkin(8) * t156 - t313 * t140 + t317 * t141;
t328 = mrSges(2,1) * t292 - mrSges(2,2) * t293 + pkin(1) * t147 + t310 * t137 + t139 * t341 + t307 * t346;
t148 = m(2) * t293 + t150;
t146 = t310 * t155 + (t154 * t314 + t342) * t307;
t144 = m(2) * t292 + t147;
t135 = mrSges(2,2) * t304 - mrSges(2,3) * t292 + t318 * t139 - t314 * t143 + (-t146 * t307 - t147 * t310) * pkin(7);
t134 = -mrSges(2,1) * t304 + mrSges(2,3) * t293 - pkin(1) * t146 - t307 * t137 + (t139 * t314 + t346) * t310;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t135 - t306 * t134 - qJ(1) * (t309 * t144 + t306 * t148), t135, t139, t141, t157, t164, t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t306 * t135 + t309 * t134 + qJ(1) * (-t306 * t144 + t309 * t148), t134, t143, t140, t151, t165, t183; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t328, t328, t137, t345, -t322, -t325, t323;];
m_new  = t1;
