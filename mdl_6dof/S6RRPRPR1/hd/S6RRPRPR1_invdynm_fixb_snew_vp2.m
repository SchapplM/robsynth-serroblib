% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 12:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:47:34
% EndTime: 2019-05-06 12:48:30
% DurationCPUTime: 42.83s
% Computational Cost: add. (734968->386), mult. (1726755->490), div. (0->0), fcn. (1297475->12), ass. (0->151)
t338 = sin(qJ(2));
t341 = cos(qJ(2));
t362 = qJD(1) * qJD(2);
t316 = qJDD(1) * t338 + t341 * t362;
t339 = sin(qJ(1));
t342 = cos(qJ(1));
t323 = -g(1) * t342 - g(2) * t339;
t343 = qJD(1) ^ 2;
t311 = -pkin(1) * t343 + qJDD(1) * pkin(7) + t323;
t365 = t338 * t311;
t366 = pkin(2) * t343;
t270 = qJDD(2) * pkin(2) - t316 * qJ(3) - t365 + (qJ(3) * t362 + t338 * t366 - g(3)) * t341;
t296 = -g(3) * t338 + t341 * t311;
t317 = qJDD(1) * t341 - t338 * t362;
t364 = qJD(1) * t338;
t319 = qJD(2) * pkin(2) - qJ(3) * t364;
t331 = t341 ^ 2;
t271 = qJ(3) * t317 - qJD(2) * t319 - t331 * t366 + t296;
t333 = sin(pkin(10));
t335 = cos(pkin(10));
t306 = (t333 * t341 + t335 * t338) * qJD(1);
t238 = -0.2e1 * qJD(3) * t306 + t335 * t270 - t333 * t271;
t294 = t316 * t335 + t317 * t333;
t305 = (-t333 * t338 + t335 * t341) * qJD(1);
t224 = (qJD(2) * t305 - t294) * pkin(8) + (t305 * t306 + qJDD(2)) * pkin(3) + t238;
t239 = 0.2e1 * qJD(3) * t305 + t333 * t270 + t335 * t271;
t293 = -t316 * t333 + t317 * t335;
t299 = qJD(2) * pkin(3) - pkin(8) * t306;
t304 = t305 ^ 2;
t228 = -pkin(3) * t304 + pkin(8) * t293 - qJD(2) * t299 + t239;
t337 = sin(qJ(4));
t367 = cos(qJ(4));
t212 = t337 * t224 + t367 * t228;
t285 = t337 * t305 + t367 * t306;
t251 = qJD(4) * t285 - t367 * t293 + t294 * t337;
t284 = -t367 * t305 + t306 * t337;
t265 = mrSges(5,1) * t284 + mrSges(5,2) * t285;
t328 = qJD(2) + qJD(4);
t279 = mrSges(5,1) * t328 - mrSges(5,3) * t285;
t327 = qJDD(2) + qJDD(4);
t264 = pkin(4) * t284 - qJ(5) * t285;
t326 = t328 ^ 2;
t206 = -pkin(4) * t326 + qJ(5) * t327 - t264 * t284 + t212;
t322 = t339 * g(1) - t342 * g(2);
t355 = -qJDD(1) * pkin(1) - t322;
t277 = -t317 * pkin(2) + qJDD(3) + t319 * t364 + (-qJ(3) * t331 - pkin(7)) * t343 + t355;
t236 = -t293 * pkin(3) - t304 * pkin(8) + t306 * t299 + t277;
t252 = -t284 * qJD(4) + t337 * t293 + t367 * t294;
t209 = (t284 * t328 - t252) * qJ(5) + (t285 * t328 + t251) * pkin(4) + t236;
t332 = sin(pkin(11));
t334 = cos(pkin(11));
t276 = t285 * t334 + t328 * t332;
t201 = -0.2e1 * qJD(5) * t276 - t332 * t206 + t334 * t209;
t243 = t252 * t334 + t327 * t332;
t275 = -t285 * t332 + t328 * t334;
t199 = (t275 * t284 - t243) * pkin(9) + (t275 * t276 + t251) * pkin(5) + t201;
t202 = 0.2e1 * qJD(5) * t275 + t334 * t206 + t332 * t209;
t242 = -t252 * t332 + t327 * t334;
t258 = pkin(5) * t284 - pkin(9) * t276;
t274 = t275 ^ 2;
t200 = -pkin(5) * t274 + pkin(9) * t242 - t258 * t284 + t202;
t336 = sin(qJ(6));
t340 = cos(qJ(6));
t197 = t199 * t340 - t200 * t336;
t253 = t275 * t340 - t276 * t336;
t218 = qJD(6) * t253 + t242 * t336 + t243 * t340;
t254 = t275 * t336 + t276 * t340;
t225 = -mrSges(7,1) * t253 + mrSges(7,2) * t254;
t280 = qJD(6) + t284;
t229 = -mrSges(7,2) * t280 + mrSges(7,3) * t253;
t250 = qJDD(6) + t251;
t192 = m(7) * t197 + mrSges(7,1) * t250 - mrSges(7,3) * t218 - t225 * t254 + t229 * t280;
t198 = t199 * t336 + t200 * t340;
t217 = -qJD(6) * t254 + t242 * t340 - t243 * t336;
t230 = mrSges(7,1) * t280 - mrSges(7,3) * t254;
t193 = m(7) * t198 - mrSges(7,2) * t250 + mrSges(7,3) * t217 + t225 * t253 - t230 * t280;
t184 = t340 * t192 + t336 * t193;
t255 = -mrSges(6,1) * t275 + mrSges(6,2) * t276;
t256 = -mrSges(6,2) * t284 + mrSges(6,3) * t275;
t182 = m(6) * t201 + mrSges(6,1) * t251 - mrSges(6,3) * t243 - t255 * t276 + t256 * t284 + t184;
t257 = mrSges(6,1) * t284 - mrSges(6,3) * t276;
t357 = -t192 * t336 + t340 * t193;
t183 = m(6) * t202 - mrSges(6,2) * t251 + mrSges(6,3) * t242 + t255 * t275 - t257 * t284 + t357;
t358 = -t182 * t332 + t334 * t183;
t175 = m(5) * t212 - mrSges(5,2) * t327 - mrSges(5,3) * t251 - t265 * t284 - t279 * t328 + t358;
t211 = t367 * t224 - t337 * t228;
t278 = -mrSges(5,2) * t328 - mrSges(5,3) * t284;
t205 = -t327 * pkin(4) - t326 * qJ(5) + t285 * t264 + qJDD(5) - t211;
t203 = -t242 * pkin(5) - t274 * pkin(9) + t276 * t258 + t205;
t352 = m(7) * t203 - t217 * mrSges(7,1) + mrSges(7,2) * t218 - t253 * t229 + t230 * t254;
t348 = -m(6) * t205 + t242 * mrSges(6,1) - mrSges(6,2) * t243 + t275 * t256 - t257 * t276 - t352;
t188 = m(5) * t211 + mrSges(5,1) * t327 - mrSges(5,3) * t252 - t265 * t285 + t278 * t328 + t348;
t166 = t337 * t175 + t367 * t188;
t288 = -mrSges(4,1) * t305 + mrSges(4,2) * t306;
t297 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t305;
t163 = m(4) * t238 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t294 + qJD(2) * t297 - t288 * t306 + t166;
t298 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t306;
t359 = t367 * t175 - t188 * t337;
t164 = m(4) * t239 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t293 - qJD(2) * t298 + t288 * t305 + t359;
t158 = t335 * t163 + t333 * t164;
t295 = -t341 * g(3) - t365;
t308 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t341) * qJD(1);
t309 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t341) * qJD(1);
t282 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * qJD(2);
t283 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * qJD(2);
t219 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t280;
t221 = Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t280;
t185 = -mrSges(7,1) * t203 + mrSges(7,3) * t198 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t250 - t219 * t254 + t221 * t280;
t220 = Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t280;
t186 = mrSges(7,2) * t203 - mrSges(7,3) * t197 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t250 + t219 * t253 - t220 * t280;
t231 = Ifges(6,5) * t276 + Ifges(6,6) * t275 + Ifges(6,3) * t284;
t233 = Ifges(6,1) * t276 + Ifges(6,4) * t275 + Ifges(6,5) * t284;
t168 = -mrSges(6,1) * t205 + mrSges(6,3) * t202 + Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t251 - pkin(5) * t352 + pkin(9) * t357 + t340 * t185 + t336 * t186 - t276 * t231 + t284 * t233;
t232 = Ifges(6,4) * t276 + Ifges(6,2) * t275 + Ifges(6,6) * t284;
t170 = mrSges(6,2) * t205 - mrSges(6,3) * t201 + Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t251 - pkin(9) * t184 - t185 * t336 + t186 * t340 + t231 * t275 - t232 * t284;
t260 = Ifges(5,4) * t285 - Ifges(5,2) * t284 + Ifges(5,6) * t328;
t261 = Ifges(5,1) * t285 - Ifges(5,4) * t284 + Ifges(5,5) * t328;
t350 = -mrSges(5,1) * t211 + mrSges(5,2) * t212 - Ifges(5,5) * t252 + Ifges(5,6) * t251 - Ifges(5,3) * t327 - pkin(4) * t348 - qJ(5) * t358 - t334 * t168 - t332 * t170 - t285 * t260 - t284 * t261;
t347 = -mrSges(4,1) * t238 + mrSges(4,2) * t239 - Ifges(4,5) * t294 - Ifges(4,6) * t293 - Ifges(4,3) * qJDD(2) - pkin(3) * t166 - t306 * t282 + t305 * t283 + t350;
t368 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * t158 + (t308 * t338 - t309 * t341) * qJD(1) - t347;
t177 = t334 * t182 + t332 * t183;
t363 = qJD(1) * t341;
t315 = (-mrSges(3,1) * t341 + mrSges(3,2) * t338) * qJD(1);
t321 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t156 = m(3) * t295 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t316 + qJD(2) * t321 - t315 * t364 + t158;
t320 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t360 = -t163 * t333 + t335 * t164;
t157 = m(3) * t296 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t317 - qJD(2) * t320 + t315 * t363 + t360;
t361 = -t156 * t338 + t341 * t157;
t354 = m(5) * t236 + t251 * mrSges(5,1) + t252 * mrSges(5,2) + t284 * t278 + t285 * t279 + t177;
t259 = Ifges(5,5) * t285 - Ifges(5,6) * t284 + Ifges(5,3) * t328;
t154 = mrSges(5,2) * t236 - mrSges(5,3) * t211 + Ifges(5,1) * t252 - Ifges(5,4) * t251 + Ifges(5,5) * t327 - qJ(5) * t177 - t168 * t332 + t170 * t334 - t259 * t284 - t260 * t328;
t351 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t218 - Ifges(7,6) * t217 - Ifges(7,3) * t250 - t254 * t220 + t253 * t221;
t345 = -mrSges(6,1) * t201 + mrSges(6,2) * t202 - Ifges(6,5) * t243 - Ifges(6,6) * t242 - pkin(5) * t184 - t276 * t232 + t275 * t233 + t351;
t159 = (-Ifges(5,2) - Ifges(6,3)) * t251 + Ifges(5,6) * t327 + t328 * t261 - t285 * t259 + Ifges(5,4) * t252 - mrSges(5,1) * t236 + mrSges(5,3) * t212 - pkin(4) * t177 + t345;
t281 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * qJD(2);
t149 = -mrSges(4,1) * t277 + mrSges(4,3) * t239 + Ifges(4,4) * t294 + Ifges(4,2) * t293 + Ifges(4,6) * qJDD(2) - pkin(3) * t354 + pkin(8) * t359 + qJD(2) * t283 + t337 * t154 + t367 * t159 - t306 * t281;
t150 = mrSges(4,2) * t277 - mrSges(4,3) * t238 + Ifges(4,1) * t294 + Ifges(4,4) * t293 + Ifges(4,5) * qJDD(2) - pkin(8) * t166 - qJD(2) * t282 + t367 * t154 - t337 * t159 + t305 * t281;
t307 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t341) * qJD(1);
t310 = -t343 * pkin(7) + t355;
t349 = m(4) * t277 - t293 * mrSges(4,1) + mrSges(4,2) * t294 - t305 * t297 + t298 * t306 + t354;
t145 = -mrSges(3,1) * t310 + mrSges(3,3) * t296 + Ifges(3,4) * t316 + Ifges(3,2) * t317 + Ifges(3,6) * qJDD(2) - pkin(2) * t349 + qJ(3) * t360 + qJD(2) * t309 + t335 * t149 + t333 * t150 - t307 * t364;
t147 = mrSges(3,2) * t310 - mrSges(3,3) * t295 + Ifges(3,1) * t316 + Ifges(3,4) * t317 + Ifges(3,5) * qJDD(2) - qJ(3) * t158 - qJD(2) * t308 - t149 * t333 + t150 * t335 + t307 * t363;
t346 = -m(3) * t310 + t317 * mrSges(3,1) - mrSges(3,2) * t316 - t320 * t364 + t321 * t363 - t349;
t353 = mrSges(2,1) * t322 - mrSges(2,2) * t323 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t361 + t341 * t145 + t338 * t147;
t171 = m(2) * t322 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t343 + t346;
t153 = t156 * t341 + t157 * t338;
t151 = m(2) * t323 - mrSges(2,1) * t343 - qJDD(1) * mrSges(2,2) + t361;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t323 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t153 - t368;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t343 - pkin(7) * t153 - t145 * t338 + t147 * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t143 - t339 * t148 - pkin(6) * (t151 * t339 + t171 * t342), t143, t147, t150, t154, t170, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t143 + t342 * t148 + pkin(6) * (t151 * t342 - t171 * t339), t148, t145, t149, t159, t168, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, t368, -t347, -t350, Ifges(6,3) * t251 - t345, -t351;];
m_new  = t1;
