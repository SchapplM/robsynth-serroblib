% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 21:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:46:33
% EndTime: 2019-05-07 21:47:19
% DurationCPUTime: 16.08s
% Computational Cost: add. (298613->383), mult. (599517->462), div. (0->0), fcn. (427791->10), ass. (0->146)
t338 = sin(qJ(3));
t342 = cos(qJ(3));
t339 = sin(qJ(2));
t368 = qJD(1) * t339;
t316 = qJD(2) * t338 + t342 * t368;
t343 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t365 = t343 * t366;
t319 = qJDD(1) * t339 + t365;
t282 = -qJD(3) * t316 + qJDD(2) * t342 - t319 * t338;
t315 = qJD(2) * t342 - t338 * t368;
t283 = qJD(3) * t315 + qJDD(2) * t338 + t319 * t342;
t337 = sin(qJ(4));
t373 = cos(qJ(4));
t285 = -t373 * t315 + t337 * t316;
t240 = -t285 * qJD(4) + t337 * t282 + t373 * t283;
t367 = qJD(1) * t343;
t330 = qJD(3) - t367;
t292 = pkin(3) * t330 - pkin(9) * t316;
t313 = t315 ^ 2;
t340 = sin(qJ(1));
t344 = cos(qJ(1));
t327 = -g(1) * t344 - g(2) * t340;
t346 = qJD(1) ^ 2;
t309 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t327;
t290 = -t343 * g(3) - t339 * t309;
t318 = (-pkin(2) * t343 - pkin(8) * t339) * qJD(1);
t345 = qJD(2) ^ 2;
t358 = qJDD(2) * pkin(2) + t345 * pkin(8) - t318 * t368 + t290;
t354 = t282 * pkin(3) + t313 * pkin(9) - t316 * t292 + t358;
t329 = qJD(4) + t330;
t371 = t285 * t329;
t376 = (-t240 + t371) * qJ(5) - t354;
t286 = t337 * t315 + t373 * t316;
t239 = t286 * qJD(4) - t373 * t282 + t337 * t283;
t274 = -pkin(5) * t329 - pkin(10) * t286;
t284 = t285 ^ 2;
t374 = 2 * qJD(5);
t197 = -t284 * pkin(10) + (-pkin(4) - pkin(5)) * t239 + (-pkin(4) * t329 + t274 + t374) * t286 - t376;
t336 = sin(qJ(6));
t341 = cos(qJ(6));
t258 = t285 * t336 + t286 * t341;
t212 = -qJD(6) * t258 + t239 * t341 - t240 * t336;
t257 = t285 * t341 - t286 * t336;
t213 = qJD(6) * t257 + t239 * t336 + t240 * t341;
t323 = qJD(6) - t329;
t246 = -mrSges(7,2) * t323 + mrSges(7,3) * t257;
t247 = mrSges(7,1) * t323 - mrSges(7,3) * t258;
t189 = -m(7) * t197 + t212 * mrSges(7,1) - t213 * mrSges(7,2) + t257 * t246 - t258 * t247;
t204 = -0.2e1 * qJD(5) * t286 + (t286 * t329 + t239) * pkin(4) + t376;
t270 = -mrSges(6,2) * t285 + mrSges(6,3) * t329;
t273 = -mrSges(6,1) * t329 + mrSges(6,2) * t286;
t185 = m(6) * t204 + t239 * mrSges(6,1) - t240 * mrSges(6,3) + t285 * t270 - t286 * t273 + t189;
t326 = t340 * g(1) - t344 * g(2);
t308 = -qJDD(1) * pkin(1) - t346 * pkin(7) - t326;
t331 = t339 * t366;
t320 = t343 * qJDD(1) - t331;
t264 = (-t319 - t365) * pkin(8) + (-t320 + t331) * pkin(2) + t308;
t291 = -g(3) * t339 + t343 * t309;
t269 = -pkin(2) * t345 + qJDD(2) * pkin(8) + t318 * t367 + t291;
t241 = t342 * t264 - t338 * t269;
t314 = qJDD(3) - t320;
t217 = (t315 * t330 - t283) * pkin(9) + (t315 * t316 + t314) * pkin(3) + t241;
t242 = t338 * t264 + t342 * t269;
t220 = -pkin(3) * t313 + pkin(9) * t282 - t292 * t330 + t242;
t207 = t337 * t217 + t373 * t220;
t252 = Ifges(6,1) * t286 + Ifges(6,4) * t329 + Ifges(6,5) * t285;
t253 = Ifges(5,1) * t286 - Ifges(5,4) * t285 + Ifges(5,5) * t329;
t310 = qJDD(4) + t314;
t206 = t373 * t217 - t337 * t220;
t259 = pkin(4) * t285 - qJ(5) * t286;
t328 = t329 ^ 2;
t202 = -t310 * pkin(4) - t328 * qJ(5) + t286 * t259 + qJDD(5) - t206;
t194 = (-t240 - t371) * pkin(10) + (t285 * t286 - t310) * pkin(5) + t202;
t200 = -pkin(4) * t328 + t310 * qJ(5) - t285 * t259 + t329 * t374 + t207;
t195 = -pkin(5) * t284 + pkin(10) * t239 + t274 * t329 + t200;
t192 = t194 * t341 - t195 * t336;
t226 = -mrSges(7,1) * t257 + mrSges(7,2) * t258;
t303 = qJDD(6) - t310;
t187 = m(7) * t192 + mrSges(7,1) * t303 - mrSges(7,3) * t213 - t226 * t258 + t246 * t323;
t193 = t194 * t336 + t195 * t341;
t188 = m(7) * t193 - mrSges(7,2) * t303 + mrSges(7,3) * t212 + t226 * t257 - t247 * t323;
t179 = -t336 * t187 + t341 * t188;
t221 = Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t323;
t223 = Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t323;
t182 = -mrSges(7,1) * t197 + mrSges(7,3) * t193 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t303 - t221 * t258 + t223 * t323;
t222 = Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t323;
t183 = mrSges(7,2) * t197 - mrSges(7,3) * t192 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t303 + t221 * t257 - t222 * t323;
t353 = -mrSges(6,1) * t204 + mrSges(6,2) * t200 - pkin(5) * t189 - pkin(10) * t179 - t341 * t182 - t336 * t183;
t250 = Ifges(6,4) * t286 + Ifges(6,2) * t329 + Ifges(6,6) * t285;
t370 = -Ifges(5,5) * t286 + Ifges(5,6) * t285 - Ifges(5,3) * t329 - t250;
t163 = mrSges(5,1) * t354 + mrSges(5,3) * t207 - pkin(4) * t185 + (t253 + t252) * t329 + (Ifges(5,6) - Ifges(6,6)) * t310 + t370 * t286 + (Ifges(5,4) - Ifges(6,5)) * t240 + (-Ifges(5,2) - Ifges(6,3)) * t239 + t353;
t251 = Ifges(5,4) * t286 - Ifges(5,2) * t285 + Ifges(5,6) * t329;
t178 = t341 * t187 + t336 * t188;
t248 = Ifges(6,5) * t286 + Ifges(6,6) * t329 + Ifges(6,3) * t285;
t355 = mrSges(6,2) * t202 - mrSges(6,3) * t204 + Ifges(6,1) * t240 + Ifges(6,4) * t310 + Ifges(6,5) * t239 - pkin(10) * t178 - t336 * t182 + t341 * t183 + t329 * t248;
t164 = -mrSges(5,2) * t354 - mrSges(5,3) * t206 + Ifges(5,1) * t240 - Ifges(5,4) * t239 + Ifges(5,5) * t310 - qJ(5) * t185 - t329 * t251 + t370 * t285 + t355;
t276 = Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t330;
t278 = Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t330;
t271 = -mrSges(5,2) * t329 - mrSges(5,3) * t285;
t272 = mrSges(5,1) * t329 - mrSges(5,3) * t286;
t351 = -m(5) * t354 + t239 * mrSges(5,1) + t240 * mrSges(5,2) + t285 * t271 + t286 * t272 + t185;
t360 = m(6) * t200 + t310 * mrSges(6,3) + t329 * t273 + t179;
t260 = mrSges(6,1) * t285 - mrSges(6,3) * t286;
t369 = -mrSges(5,1) * t285 - mrSges(5,2) * t286 - t260;
t372 = -mrSges(5,3) - mrSges(6,2);
t171 = m(5) * t207 - t310 * mrSges(5,2) + t372 * t239 - t329 * t272 + t369 * t285 + t360;
t356 = -m(6) * t202 + t310 * mrSges(6,1) + t329 * t270 - t178;
t173 = m(5) * t206 + t310 * mrSges(5,1) + t372 * t240 + t329 * t271 + t369 * t286 + t356;
t363 = t373 * t171 - t173 * t337;
t152 = mrSges(4,1) * t358 + mrSges(4,3) * t242 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * t314 - pkin(3) * t351 + pkin(9) * t363 + t373 * t163 + t337 * t164 - t316 * t276 + t330 * t278;
t168 = t337 * t171 + t373 * t173;
t277 = Ifges(4,4) * t316 + Ifges(4,2) * t315 + Ifges(4,6) * t330;
t153 = -mrSges(4,2) * t358 - mrSges(4,3) * t241 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * t314 - pkin(9) * t168 - t337 * t163 + t373 * t164 + t315 * t276 - t330 * t277;
t287 = -mrSges(4,1) * t315 + mrSges(4,2) * t316;
t288 = -mrSges(4,2) * t330 + mrSges(4,3) * t315;
t166 = m(4) * t241 + mrSges(4,1) * t314 - mrSges(4,3) * t283 - t287 * t316 + t288 * t330 + t168;
t289 = mrSges(4,1) * t330 - mrSges(4,3) * t316;
t167 = m(4) * t242 - mrSges(4,2) * t314 + mrSges(4,3) * t282 + t287 * t315 - t289 * t330 + t363;
t162 = -t166 * t338 + t342 * t167;
t184 = m(4) * t358 + t282 * mrSges(4,1) - t283 * mrSges(4,2) + t315 * t288 - t316 * t289 - t351;
t306 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t343) * qJD(1);
t307 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t343) * qJD(1);
t375 = mrSges(3,1) * t290 - mrSges(3,2) * t291 + Ifges(3,5) * t319 + Ifges(3,6) * t320 + Ifges(3,3) * qJDD(2) + pkin(2) * t184 + pkin(8) * t162 + t342 * t152 + t338 * t153 + (t306 * t339 - t307 * t343) * qJD(1);
t317 = (-mrSges(3,1) * t343 + mrSges(3,2) * t339) * qJD(1);
t324 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t368;
t160 = m(3) * t291 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t320 - qJD(2) * t324 + t317 * t367 + t162;
t325 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t367;
t180 = m(3) * t290 + qJDD(2) * mrSges(3,1) - t319 * mrSges(3,3) + qJD(2) * t325 - t317 * t368 + t184;
t364 = t343 * t160 - t180 * t339;
t161 = t166 * t342 + t167 * t338;
t359 = mrSges(7,1) * t192 - mrSges(7,2) * t193 + Ifges(7,5) * t213 + Ifges(7,6) * t212 + Ifges(7,3) * t303 + t258 * t222 - t257 * t223;
t305 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t343) * qJD(1);
t149 = mrSges(3,2) * t308 - mrSges(3,3) * t290 + Ifges(3,1) * t319 + Ifges(3,4) * t320 + Ifges(3,5) * qJDD(2) - pkin(8) * t161 - qJD(2) * t306 - t152 * t338 + t153 * t342 + t305 * t367;
t350 = mrSges(6,1) * t202 - mrSges(6,3) * t200 - Ifges(6,4) * t240 - Ifges(6,2) * t310 - Ifges(6,6) * t239 + pkin(5) * t178 + t286 * t248 - t285 * t252 + t359;
t348 = mrSges(5,2) * t207 - t285 * t253 - qJ(5) * (-t239 * mrSges(6,2) - t285 * t260 + t360) - pkin(4) * (-t240 * mrSges(6,2) - t286 * t260 + t356) - mrSges(5,1) * t206 + Ifges(5,6) * t239 - Ifges(5,5) * t240 - t286 * t251 - Ifges(5,3) * t310 + t350;
t347 = mrSges(4,1) * t241 - mrSges(4,2) * t242 + Ifges(4,5) * t283 + Ifges(4,6) * t282 + Ifges(4,3) * t314 + pkin(3) * t168 + t316 * t277 - t315 * t278 - t348;
t151 = -mrSges(3,1) * t308 + mrSges(3,3) * t291 + Ifges(3,4) * t319 + Ifges(3,2) * t320 + Ifges(3,6) * qJDD(2) - pkin(2) * t161 + qJD(2) * t307 - t305 * t368 - t347;
t352 = -m(3) * t308 + t320 * mrSges(3,1) - mrSges(3,2) * t319 - t324 * t368 + t325 * t367 - t161;
t357 = mrSges(2,1) * t326 - mrSges(2,2) * t327 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t364 + t339 * t149 + t343 * t151;
t157 = m(2) * t326 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t346 + t352;
t156 = t160 * t339 + t180 * t343;
t154 = m(2) * t327 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t364;
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t327 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t156 - t375;
t146 = -mrSges(2,2) * g(3) - mrSges(2,3) * t326 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t156 + t149 * t343 - t151 * t339;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t344 * t146 - t340 * t147 - pkin(6) * (t154 * t340 + t157 * t344), t146, t149, t153, t164, -t250 * t285 + t355, t183; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t146 + t344 * t147 + pkin(6) * (t154 * t344 - t157 * t340), t147, t151, t152, t163, -t350, t182; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t357, t357, t375, t347, -t348, Ifges(6,5) * t240 + Ifges(6,6) * t310 + Ifges(6,3) * t239 + t286 * t250 - t329 * t252 - t353, t359;];
m_new  = t1;
