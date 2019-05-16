% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:38:05
% EndTime: 2019-05-08 05:39:52
% DurationCPUTime: 40.32s
% Computational Cost: add. (736570->395), mult. (1565264->497), div. (0->0), fcn. (1262278->12), ass. (0->156)
t341 = sin(qJ(2));
t345 = cos(qJ(2));
t336 = sin(pkin(6));
t367 = qJD(1) * t336;
t319 = (-pkin(2) * t345 - pkin(9) * t341) * t367;
t337 = cos(pkin(6));
t332 = qJD(1) * t337 + qJD(2);
t330 = t332 ^ 2;
t331 = qJDD(1) * t337 + qJDD(2);
t366 = qJD(1) * t345;
t342 = sin(qJ(1));
t346 = cos(qJ(1));
t327 = t342 * g(1) - g(2) * t346;
t347 = qJD(1) ^ 2;
t378 = pkin(8) * t336;
t316 = qJDD(1) * pkin(1) + t347 * t378 + t327;
t328 = -g(1) * t346 - g(2) * t342;
t365 = qJDD(1) * t336;
t317 = -pkin(1) * t347 + pkin(8) * t365 + t328;
t373 = t337 * t341;
t368 = t316 * t373 + t345 * t317;
t275 = -t330 * pkin(2) + t331 * pkin(9) + (-g(3) * t341 + t319 * t366) * t336 + t368;
t320 = (qJD(2) * t366 + qJDD(1) * t341) * t336;
t363 = t341 * t367;
t321 = -qJD(2) * t363 + t345 * t365;
t377 = t337 * g(3);
t276 = -t321 * pkin(2) - t320 * pkin(9) - t377 + (-t316 + (pkin(2) * t341 - pkin(9) * t345) * t332 * qJD(1)) * t336;
t340 = sin(qJ(3));
t344 = cos(qJ(3));
t234 = -t340 * t275 + t344 * t276;
t308 = t332 * t344 - t340 * t363;
t288 = qJD(3) * t308 + t320 * t344 + t331 * t340;
t309 = t332 * t340 + t344 * t363;
t313 = qJDD(3) - t321;
t362 = t336 * t366;
t326 = qJD(3) - t362;
t219 = (t308 * t326 - t288) * pkin(10) + (t308 * t309 + t313) * pkin(3) + t234;
t235 = t344 * t275 + t340 * t276;
t287 = -qJD(3) * t309 - t320 * t340 + t331 * t344;
t298 = pkin(3) * t326 - pkin(10) * t309;
t307 = t308 ^ 2;
t222 = -pkin(3) * t307 + pkin(10) * t287 - t298 * t326 + t235;
t339 = sin(qJ(4));
t343 = cos(qJ(4));
t217 = t339 * t219 + t343 * t222;
t293 = t308 * t343 - t309 * t339;
t294 = t308 * t339 + t309 * t343;
t270 = -pkin(4) * t293 - pkin(11) * t294;
t312 = qJDD(4) + t313;
t324 = qJD(4) + t326;
t323 = t324 ^ 2;
t212 = -pkin(4) * t323 + pkin(11) * t312 + t270 * t293 + t217;
t372 = t337 * t345;
t374 = t336 * t345;
t289 = -g(3) * t374 + t316 * t372 - t341 * t317;
t274 = -t331 * pkin(2) - t330 * pkin(9) + t319 * t363 - t289;
t233 = -t287 * pkin(3) - t307 * pkin(10) + t309 * t298 + t274;
t254 = -qJD(4) * t294 + t287 * t343 - t288 * t339;
t255 = qJD(4) * t293 + t287 * t339 + t288 * t343;
t214 = (-t293 * t324 - t255) * pkin(11) + (t294 * t324 - t254) * pkin(4) + t233;
t338 = sin(qJ(5));
t379 = cos(qJ(5));
t208 = -t338 * t212 + t214 * t379;
t209 = t212 * t379 + t338 * t214;
t278 = t294 * t379 + t338 * t324;
t229 = qJD(5) * t278 + t255 * t338 - t312 * t379;
t277 = t294 * t338 - t324 * t379;
t230 = -t277 * qJD(5) + t255 * t379 + t338 * t312;
t292 = qJD(5) - t293;
t236 = Ifges(7,5) * t278 + Ifges(7,6) * t292 + Ifges(7,3) * t277;
t239 = Ifges(6,4) * t278 - Ifges(6,2) * t277 + Ifges(6,6) * t292;
t241 = Ifges(6,1) * t278 - Ifges(6,4) * t277 + Ifges(6,5) * t292;
t253 = qJDD(5) - t254;
t259 = mrSges(7,1) * t277 - mrSges(7,3) * t278;
t258 = pkin(5) * t277 - qJ(6) * t278;
t291 = t292 ^ 2;
t204 = -pkin(5) * t291 + qJ(6) * t253 + 0.2e1 * qJD(6) * t292 - t258 * t277 + t209;
t206 = -t253 * pkin(5) - t291 * qJ(6) + t278 * t258 + qJDD(6) - t208;
t240 = Ifges(7,1) * t278 + Ifges(7,4) * t292 + Ifges(7,5) * t277;
t356 = mrSges(7,1) * t206 - mrSges(7,3) * t204 - Ifges(7,4) * t230 - Ifges(7,2) * t253 - Ifges(7,6) * t229 - t277 * t240;
t261 = -mrSges(7,2) * t277 + mrSges(7,3) * t292;
t358 = -m(7) * t206 + t253 * mrSges(7,1) + t292 * t261;
t264 = -mrSges(7,1) * t292 + mrSges(7,2) * t278;
t364 = m(7) * t204 + t253 * mrSges(7,3) + t292 * t264;
t380 = -(-t239 + t236) * t278 + mrSges(6,1) * t208 - mrSges(6,2) * t209 + Ifges(6,5) * t230 - Ifges(6,6) * t229 + Ifges(6,3) * t253 + pkin(5) * (-t230 * mrSges(7,2) - t278 * t259 + t358) + qJ(6) * (-t229 * mrSges(7,2) - t277 * t259 + t364) + t277 * t241 - t356;
t376 = -mrSges(6,3) - mrSges(7,2);
t375 = t336 * t341;
t269 = -mrSges(5,1) * t293 + mrSges(5,2) * t294;
t280 = mrSges(5,1) * t324 - mrSges(5,3) * t294;
t263 = mrSges(6,1) * t292 - mrSges(6,3) * t278;
t369 = -mrSges(6,1) * t277 - mrSges(6,2) * t278 - t259;
t194 = m(6) * t209 - t253 * mrSges(6,2) + t229 * t376 - t292 * t263 + t277 * t369 + t364;
t262 = -mrSges(6,2) * t292 - mrSges(6,3) * t277;
t196 = m(6) * t208 + t253 * mrSges(6,1) + t230 * t376 + t292 * t262 + t278 * t369 + t358;
t359 = t194 * t379 - t196 * t338;
t182 = m(5) * t217 - mrSges(5,2) * t312 + mrSges(5,3) * t254 + t269 * t293 - t280 * t324 + t359;
t216 = t343 * t219 - t339 * t222;
t279 = -mrSges(5,2) * t324 + mrSges(5,3) * t293;
t211 = -t312 * pkin(4) - t323 * pkin(11) + t294 * t270 - t216;
t207 = -0.2e1 * qJD(6) * t278 + (t277 * t292 - t230) * qJ(6) + (t278 * t292 + t229) * pkin(5) + t211;
t201 = m(7) * t207 + mrSges(7,1) * t229 - t230 * mrSges(7,3) + t261 * t277 - t278 * t264;
t351 = -m(6) * t211 - t229 * mrSges(6,1) - mrSges(6,2) * t230 - t277 * t262 - t263 * t278 - t201;
t191 = m(5) * t216 + mrSges(5,1) * t312 - mrSges(5,3) * t255 - t269 * t294 + t279 * t324 + t351;
t177 = t339 * t182 + t343 * t191;
t295 = -mrSges(4,1) * t308 + mrSges(4,2) * t309;
t296 = -mrSges(4,2) * t326 + mrSges(4,3) * t308;
t175 = m(4) * t234 + mrSges(4,1) * t313 - mrSges(4,3) * t288 - t295 * t309 + t296 * t326 + t177;
t297 = mrSges(4,1) * t326 - mrSges(4,3) * t309;
t360 = t343 * t182 - t191 * t339;
t176 = m(4) * t235 - mrSges(4,2) * t313 + mrSges(4,3) * t287 + t295 * t308 - t297 * t326 + t360;
t169 = t344 * t175 + t340 * t176;
t188 = t338 * t194 + t196 * t379;
t238 = Ifges(7,4) * t278 + Ifges(7,2) * t292 + Ifges(7,6) * t277;
t371 = -Ifges(6,5) * t278 + Ifges(6,6) * t277 - Ifges(6,3) * t292 - t238;
t290 = -g(3) * t375 + t368;
t314 = mrSges(3,1) * t332 - mrSges(3,3) * t363;
t318 = (-mrSges(3,1) * t345 + mrSges(3,2) * t341) * t367;
t361 = -t175 * t340 + t344 * t176;
t167 = m(3) * t290 - mrSges(3,2) * t331 + mrSges(3,3) * t321 - t314 * t332 + t318 * t362 + t361;
t315 = -mrSges(3,2) * t332 + mrSges(3,3) * t362;
t353 = m(5) * t233 - t254 * mrSges(5,1) + mrSges(5,2) * t255 - t293 * t279 + t280 * t294 + t188;
t350 = -m(4) * t274 + t287 * mrSges(4,1) - mrSges(4,2) * t288 + t308 * t296 - t297 * t309 - t353;
t179 = m(3) * t289 + mrSges(3,1) * t331 - mrSges(3,3) * t320 + t315 * t332 - t318 * t363 + t350;
t164 = t345 * t167 - t179 * t341;
t302 = -t336 * t316 - t377;
t168 = m(3) * t302 - t321 * mrSges(3,1) + t320 * mrSges(3,2) + (t314 * t341 - t315 * t345) * t367 + t169;
t160 = t167 * t373 - t168 * t336 + t179 * t372;
t357 = -mrSges(7,1) * t207 + mrSges(7,2) * t204;
t355 = mrSges(7,2) * t206 - mrSges(7,3) * t207 + Ifges(7,1) * t230 + Ifges(7,4) * t253 + Ifges(7,5) * t229 + t292 * t236;
t184 = -mrSges(6,1) * t211 + mrSges(6,3) * t209 - pkin(5) * t201 + (t240 + t241) * t292 + t371 * t278 + (Ifges(6,6) - Ifges(7,6)) * t253 + (Ifges(6,4) - Ifges(7,5)) * t230 + (-Ifges(6,2) - Ifges(7,3)) * t229 + t357;
t186 = mrSges(6,2) * t211 - mrSges(6,3) * t208 + Ifges(6,1) * t230 - Ifges(6,4) * t229 + Ifges(6,5) * t253 - qJ(6) * t201 - t292 * t239 + t277 * t371 + t355;
t265 = Ifges(5,5) * t294 + Ifges(5,6) * t293 + Ifges(5,3) * t324;
t266 = Ifges(5,4) * t294 + Ifges(5,2) * t293 + Ifges(5,6) * t324;
t170 = mrSges(5,2) * t233 - mrSges(5,3) * t216 + Ifges(5,1) * t255 + Ifges(5,4) * t254 + Ifges(5,5) * t312 - pkin(11) * t188 - t338 * t184 + t186 * t379 + t293 * t265 - t324 * t266;
t267 = Ifges(5,1) * t294 + Ifges(5,4) * t293 + Ifges(5,5) * t324;
t171 = -mrSges(5,1) * t233 + mrSges(5,3) * t217 + Ifges(5,4) * t255 + Ifges(5,2) * t254 + Ifges(5,6) * t312 - pkin(4) * t188 - t294 * t265 + t324 * t267 - t380;
t281 = Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * t326;
t283 = Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * t326;
t156 = -mrSges(4,1) * t274 + mrSges(4,3) * t235 + Ifges(4,4) * t288 + Ifges(4,2) * t287 + Ifges(4,6) * t313 - pkin(3) * t353 + pkin(10) * t360 + t339 * t170 + t343 * t171 - t309 * t281 + t326 * t283;
t282 = Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * t326;
t161 = mrSges(4,2) * t274 - mrSges(4,3) * t234 + Ifges(4,1) * t288 + Ifges(4,4) * t287 + Ifges(4,5) * t313 - pkin(10) * t177 + t170 * t343 - t171 * t339 + t281 * t308 - t282 * t326;
t300 = Ifges(3,6) * t332 + (Ifges(3,4) * t341 + Ifges(3,2) * t345) * t367;
t301 = Ifges(3,5) * t332 + (Ifges(3,1) * t341 + Ifges(3,4) * t345) * t367;
t151 = Ifges(3,5) * t320 + Ifges(3,6) * t321 + Ifges(3,3) * t331 + mrSges(3,1) * t289 - mrSges(3,2) * t290 + t340 * t161 + t344 * t156 + pkin(2) * t350 + pkin(9) * t361 + (t300 * t341 - t301 * t345) * t367;
t299 = Ifges(3,3) * t332 + (Ifges(3,5) * t341 + Ifges(3,6) * t345) * t367;
t153 = mrSges(3,2) * t302 - mrSges(3,3) * t289 + Ifges(3,1) * t320 + Ifges(3,4) * t321 + Ifges(3,5) * t331 - pkin(9) * t169 - t156 * t340 + t161 * t344 + t299 * t362 - t300 * t332;
t352 = -mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t255 - Ifges(5,6) * t254 - Ifges(5,3) * t312 - pkin(4) * t351 - pkin(11) * t359 - t184 * t379 - t338 * t186 - t294 * t266 + t293 * t267;
t348 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t288 + Ifges(4,6) * t287 + Ifges(4,3) * t313 + pkin(3) * t177 + t309 * t282 - t308 * t283 - t352;
t155 = -mrSges(3,1) * t302 + mrSges(3,3) * t290 + Ifges(3,4) * t320 + Ifges(3,2) * t321 + Ifges(3,6) * t331 - pkin(2) * t169 - t299 * t363 + t332 * t301 - t348;
t354 = mrSges(2,1) * t327 - mrSges(2,2) * t328 + Ifges(2,3) * qJDD(1) + pkin(1) * t160 + t337 * t151 + t153 * t375 + t155 * t374 + t164 * t378;
t162 = m(2) * t328 - mrSges(2,1) * t347 - qJDD(1) * mrSges(2,2) + t164;
t159 = t337 * t168 + (t167 * t341 + t179 * t345) * t336;
t157 = m(2) * t327 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t347 + t160;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - t347 * Ifges(2,6) + t345 * t153 - t341 * t155 + (-t159 * t336 - t160 * t337) * pkin(8);
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t159 - t336 * t151 + (pkin(8) * t164 + t153 * t341 + t155 * t345) * t337;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t149 - t342 * t148 - pkin(7) * (t157 * t346 + t162 * t342), t149, t153, t161, t170, t186, -t238 * t277 + t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t149 + t346 * t148 + pkin(7) * (-t157 * t342 + t162 * t346), t148, t155, t156, t171, t184, -t278 * t236 - t356; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t354, t354, t151, t348, -t352, t380, Ifges(7,5) * t230 + Ifges(7,6) * t253 + Ifges(7,3) * t229 + t278 * t238 - t292 * t240 - t357;];
m_new  = t1;
