% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:22:39
% EndTime: 2019-05-06 17:23:16
% DurationCPUTime: 20.14s
% Computational Cost: add. (335091->382), mult. (775738->470), div. (0->0), fcn. (570747->10), ass. (0->142)
t330 = sin(qJ(2));
t333 = cos(qJ(2));
t356 = qJD(1) * qJD(2);
t308 = t330 * qJDD(1) + t333 * t356;
t331 = sin(qJ(1));
t334 = cos(qJ(1));
t315 = -t334 * g(1) - t331 * g(2);
t335 = qJD(1) ^ 2;
t303 = -t335 * pkin(1) + qJDD(1) * pkin(7) + t315;
t362 = t330 * t303;
t364 = pkin(2) * t335;
t266 = qJDD(2) * pkin(2) - t308 * qJ(3) - t362 + (qJ(3) * t356 + t330 * t364 - g(3)) * t333;
t289 = -t330 * g(3) + t333 * t303;
t309 = t333 * qJDD(1) - t330 * t356;
t358 = qJD(1) * t330;
t311 = qJD(2) * pkin(2) - qJ(3) * t358;
t325 = t333 ^ 2;
t267 = t309 * qJ(3) - qJD(2) * t311 - t325 * t364 + t289;
t326 = sin(pkin(10));
t327 = cos(pkin(10));
t298 = (t326 * t333 + t327 * t330) * qJD(1);
t233 = -0.2e1 * qJD(3) * t298 + t327 * t266 - t326 * t267;
t287 = t327 * t308 + t326 * t309;
t297 = (-t326 * t330 + t327 * t333) * qJD(1);
t206 = (qJD(2) * t297 - t287) * pkin(8) + (t297 * t298 + qJDD(2)) * pkin(3) + t233;
t234 = 0.2e1 * qJD(3) * t297 + t326 * t266 + t327 * t267;
t286 = -t326 * t308 + t327 * t309;
t292 = qJD(2) * pkin(3) - t298 * pkin(8);
t296 = t297 ^ 2;
t209 = -t296 * pkin(3) + t286 * pkin(8) - qJD(2) * t292 + t234;
t329 = sin(qJ(4));
t332 = cos(qJ(4));
t204 = t329 * t206 + t332 * t209;
t278 = t332 * t297 - t329 * t298;
t279 = t329 * t297 + t332 * t298;
t261 = -t278 * pkin(4) - t279 * pkin(9);
t322 = qJD(2) + qJD(4);
t320 = t322 ^ 2;
t321 = qJDD(2) + qJDD(4);
t199 = -t320 * pkin(4) + t321 * pkin(9) + t278 * t261 + t204;
t314 = t331 * g(1) - t334 * g(2);
t347 = -qJDD(1) * pkin(1) - t314;
t268 = -t309 * pkin(2) + qJDD(3) + t311 * t358 + (-qJ(3) * t325 - pkin(7)) * t335 + t347;
t229 = -t286 * pkin(3) - t296 * pkin(8) + t298 * t292 + t268;
t246 = -t279 * qJD(4) + t332 * t286 - t329 * t287;
t247 = t278 * qJD(4) + t329 * t286 + t332 * t287;
t201 = (-t278 * t322 - t247) * pkin(9) + (t279 * t322 - t246) * pkin(4) + t229;
t328 = sin(qJ(5));
t365 = cos(qJ(5));
t195 = -t328 * t199 + t201 * t365;
t196 = t199 * t365 + t328 * t201;
t270 = t279 * t365 + t328 * t322;
t218 = t270 * qJD(5) + t328 * t247 - t321 * t365;
t269 = t328 * t279 - t322 * t365;
t219 = -t269 * qJD(5) + t247 * t365 + t328 * t321;
t274 = qJD(5) - t278;
t220 = Ifges(7,5) * t270 + Ifges(7,6) * t274 + Ifges(7,3) * t269;
t223 = Ifges(6,4) * t270 - Ifges(6,2) * t269 + Ifges(6,6) * t274;
t225 = Ifges(6,1) * t270 - Ifges(6,4) * t269 + Ifges(6,5) * t274;
t245 = qJDD(5) - t246;
t249 = t269 * mrSges(7,1) - t270 * mrSges(7,3);
t248 = t269 * pkin(5) - t270 * qJ(6);
t273 = t274 ^ 2;
t191 = -t273 * pkin(5) + t245 * qJ(6) + 0.2e1 * qJD(6) * t274 - t269 * t248 + t196;
t193 = -t245 * pkin(5) - t273 * qJ(6) + t270 * t248 + qJDD(6) - t195;
t224 = Ifges(7,1) * t270 + Ifges(7,4) * t274 + Ifges(7,5) * t269;
t345 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t219 - Ifges(7,2) * t245 - Ifges(7,6) * t218 - t269 * t224;
t251 = -t269 * mrSges(7,2) + t274 * mrSges(7,3);
t350 = -m(7) * t193 + t245 * mrSges(7,1) + t274 * t251;
t254 = -t274 * mrSges(7,1) + t270 * mrSges(7,2);
t355 = m(7) * t191 + t245 * mrSges(7,3) + t274 * t254;
t367 = -(-t223 + t220) * t270 + mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t219 - Ifges(6,6) * t218 + Ifges(6,3) * t245 + pkin(5) * (-t219 * mrSges(7,2) - t270 * t249 + t350) + qJ(6) * (-t218 * mrSges(7,2) - t269 * t249 + t355) + t269 * t225 - t345;
t260 = -t278 * mrSges(5,1) + t279 * mrSges(5,2);
t272 = t322 * mrSges(5,1) - t279 * mrSges(5,3);
t253 = t274 * mrSges(6,1) - t270 * mrSges(6,3);
t359 = -t269 * mrSges(6,1) - t270 * mrSges(6,2) - t249;
t363 = -mrSges(6,3) - mrSges(7,2);
t181 = m(6) * t196 - t245 * mrSges(6,2) + t218 * t363 - t274 * t253 + t269 * t359 + t355;
t252 = -t274 * mrSges(6,2) - t269 * mrSges(6,3);
t183 = m(6) * t195 + t245 * mrSges(6,1) + t219 * t363 + t274 * t252 + t270 * t359 + t350;
t351 = t181 * t365 - t328 * t183;
t169 = m(5) * t204 - t321 * mrSges(5,2) + t246 * mrSges(5,3) + t278 * t260 - t322 * t272 + t351;
t203 = t332 * t206 - t329 * t209;
t271 = -t322 * mrSges(5,2) + t278 * mrSges(5,3);
t198 = -t321 * pkin(4) - t320 * pkin(9) + t279 * t261 - t203;
t194 = -0.2e1 * qJD(6) * t270 + (t269 * t274 - t219) * qJ(6) + (t270 * t274 + t218) * pkin(5) + t198;
t188 = m(7) * t194 + t218 * mrSges(7,1) - t219 * mrSges(7,3) + t269 * t251 - t270 * t254;
t340 = -m(6) * t198 - t218 * mrSges(6,1) - t219 * mrSges(6,2) - t269 * t252 - t270 * t253 - t188;
t178 = m(5) * t203 + t321 * mrSges(5,1) - t247 * mrSges(5,3) - t279 * t260 + t322 * t271 + t340;
t164 = t329 * t169 + t332 * t178;
t282 = -t297 * mrSges(4,1) + t298 * mrSges(4,2);
t290 = -qJD(2) * mrSges(4,2) + t297 * mrSges(4,3);
t161 = m(4) * t233 + qJDD(2) * mrSges(4,1) - t287 * mrSges(4,3) + qJD(2) * t290 - t298 * t282 + t164;
t291 = qJD(2) * mrSges(4,1) - t298 * mrSges(4,3);
t352 = t332 * t169 - t329 * t178;
t162 = m(4) * t234 - qJDD(2) * mrSges(4,2) + t286 * mrSges(4,3) - qJD(2) * t291 + t297 * t282 + t352;
t155 = t327 * t161 + t326 * t162;
t288 = -t333 * g(3) - t362;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t330 + Ifges(3,2) * t333) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t330 + Ifges(3,4) * t333) * qJD(1);
t276 = Ifges(4,4) * t298 + Ifges(4,2) * t297 + Ifges(4,6) * qJD(2);
t277 = Ifges(4,1) * t298 + Ifges(4,4) * t297 + Ifges(4,5) * qJD(2);
t349 = -mrSges(7,1) * t194 + mrSges(7,2) * t191;
t222 = Ifges(7,4) * t270 + Ifges(7,2) * t274 + Ifges(7,6) * t269;
t361 = -Ifges(6,5) * t270 + Ifges(6,6) * t269 - Ifges(6,3) * t274 - t222;
t171 = -mrSges(6,1) * t198 + mrSges(6,3) * t196 - pkin(5) * t188 + (t224 + t225) * t274 + t361 * t270 + (Ifges(6,6) - Ifges(7,6)) * t245 + (Ifges(6,4) - Ifges(7,5)) * t219 + (-Ifges(6,2) - Ifges(7,3)) * t218 + t349;
t344 = mrSges(7,2) * t193 - mrSges(7,3) * t194 + Ifges(7,1) * t219 + Ifges(7,4) * t245 + Ifges(7,5) * t218 + t274 * t220;
t173 = mrSges(6,2) * t198 - mrSges(6,3) * t195 + Ifges(6,1) * t219 - Ifges(6,4) * t218 + Ifges(6,5) * t245 - qJ(6) * t188 - t274 * t223 + t269 * t361 + t344;
t256 = Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t322;
t257 = Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t322;
t342 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t247 - Ifges(5,6) * t246 - Ifges(5,3) * t321 - pkin(4) * t340 - pkin(9) * t351 - t171 * t365 - t328 * t173 - t279 * t256 + t278 * t257;
t339 = -mrSges(4,1) * t233 + mrSges(4,2) * t234 - Ifges(4,5) * t287 - Ifges(4,6) * t286 - Ifges(4,3) * qJDD(2) - pkin(3) * t164 - t298 * t276 + t297 * t277 + t342;
t366 = mrSges(3,1) * t288 - mrSges(3,2) * t289 + Ifges(3,5) * t308 + Ifges(3,6) * t309 + Ifges(3,3) * qJDD(2) + pkin(2) * t155 + (t330 * t300 - t333 * t301) * qJD(1) - t339;
t175 = t328 * t181 + t183 * t365;
t357 = qJD(1) * t333;
t307 = (-mrSges(3,1) * t333 + mrSges(3,2) * t330) * qJD(1);
t313 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t357;
t153 = m(3) * t288 + qJDD(2) * mrSges(3,1) - t308 * mrSges(3,3) + qJD(2) * t313 - t307 * t358 + t155;
t312 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t358;
t353 = -t326 * t161 + t327 * t162;
t154 = m(3) * t289 - qJDD(2) * mrSges(3,2) + t309 * mrSges(3,3) - qJD(2) * t312 + t307 * t357 + t353;
t354 = -t330 * t153 + t333 * t154;
t346 = m(5) * t229 - t246 * mrSges(5,1) + t247 * mrSges(5,2) - t278 * t271 + t279 * t272 + t175;
t255 = Ifges(5,5) * t279 + Ifges(5,6) * t278 + Ifges(5,3) * t322;
t156 = mrSges(5,2) * t229 - mrSges(5,3) * t203 + Ifges(5,1) * t247 + Ifges(5,4) * t246 + Ifges(5,5) * t321 - pkin(9) * t175 - t328 * t171 + t173 * t365 + t278 * t255 - t322 * t256;
t157 = -mrSges(5,1) * t229 + mrSges(5,3) * t204 + Ifges(5,4) * t247 + Ifges(5,2) * t246 + Ifges(5,6) * t321 - pkin(4) * t175 - t279 * t255 + t322 * t257 - t367;
t275 = Ifges(4,5) * t298 + Ifges(4,6) * t297 + Ifges(4,3) * qJD(2);
t147 = -mrSges(4,1) * t268 + mrSges(4,3) * t234 + Ifges(4,4) * t287 + Ifges(4,2) * t286 + Ifges(4,6) * qJDD(2) - pkin(3) * t346 + pkin(8) * t352 + qJD(2) * t277 + t329 * t156 + t332 * t157 - t298 * t275;
t148 = mrSges(4,2) * t268 - mrSges(4,3) * t233 + Ifges(4,1) * t287 + Ifges(4,4) * t286 + Ifges(4,5) * qJDD(2) - pkin(8) * t164 - qJD(2) * t276 + t332 * t156 - t329 * t157 + t297 * t275;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t330 + Ifges(3,6) * t333) * qJD(1);
t302 = -t335 * pkin(7) + t347;
t341 = m(4) * t268 - t286 * mrSges(4,1) + t287 * mrSges(4,2) - t297 * t290 + t298 * t291 + t346;
t143 = -mrSges(3,1) * t302 + mrSges(3,3) * t289 + Ifges(3,4) * t308 + Ifges(3,2) * t309 + Ifges(3,6) * qJDD(2) - pkin(2) * t341 + qJ(3) * t353 + qJD(2) * t301 + t327 * t147 + t326 * t148 - t299 * t358;
t145 = mrSges(3,2) * t302 - mrSges(3,3) * t288 + Ifges(3,1) * t308 + Ifges(3,4) * t309 + Ifges(3,5) * qJDD(2) - qJ(3) * t155 - qJD(2) * t300 - t326 * t147 + t327 * t148 + t299 * t357;
t338 = -m(3) * t302 + t309 * mrSges(3,1) - t308 * mrSges(3,2) - t312 * t358 + t313 * t357 - t341;
t343 = mrSges(2,1) * t314 - mrSges(2,2) * t315 + Ifges(2,3) * qJDD(1) + pkin(1) * t338 + pkin(7) * t354 + t333 * t143 + t330 * t145;
t165 = m(2) * t314 + qJDD(1) * mrSges(2,1) - t335 * mrSges(2,2) + t338;
t151 = t333 * t153 + t330 * t154;
t149 = m(2) * t315 - t335 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t354;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t315 + t335 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t366;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t314 + Ifges(2,5) * qJDD(1) - t335 * Ifges(2,6) - pkin(7) * t151 - t330 * t143 + t333 * t145;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t334 * t141 - t331 * t146 - pkin(6) * (t331 * t149 + t334 * t165), t141, t145, t148, t156, t173, -t269 * t222 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t331 * t141 + t334 * t146 + pkin(6) * (t334 * t149 - t331 * t165), t146, t143, t147, t157, t171, -t270 * t220 - t345; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t366, -t339, -t342, t367, Ifges(7,5) * t219 + Ifges(7,6) * t245 + Ifges(7,3) * t218 + t270 * t222 - t274 * t224 - t349;];
m_new  = t1;
