% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:40:52
% EndTime: 2019-05-07 04:41:39
% DurationCPUTime: 15.98s
% Computational Cost: add. (276420->382), mult. (574986->463), div. (0->0), fcn. (404439->10), ass. (0->145)
t380 = -2 * qJD(4);
t340 = sin(qJ(1));
t344 = cos(qJ(1));
t327 = t340 * g(1) - t344 * g(2);
t346 = qJD(1) ^ 2;
t304 = -qJDD(1) * pkin(1) - t346 * pkin(7) - t327;
t339 = sin(qJ(2));
t343 = cos(qJ(2));
t368 = qJD(1) * qJD(2);
t366 = t343 * t368;
t321 = t339 * qJDD(1) + t366;
t367 = t339 * t368;
t322 = t343 * qJDD(1) - t367;
t265 = (-t321 - t366) * pkin(8) + (-t322 + t367) * pkin(2) + t304;
t328 = -t344 * g(1) - t340 * g(2);
t305 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t328;
t295 = -t339 * g(3) + t343 * t305;
t320 = (-pkin(2) * t343 - pkin(8) * t339) * qJD(1);
t345 = qJD(2) ^ 2;
t369 = t343 * qJD(1);
t270 = -t345 * pkin(2) + qJDD(2) * pkin(8) + t320 * t369 + t295;
t338 = sin(qJ(3));
t342 = cos(qJ(3));
t231 = t342 * t265 - t338 * t270;
t371 = qJD(1) * t339;
t317 = t342 * qJD(2) - t338 * t371;
t286 = t317 * qJD(3) + t338 * qJDD(2) + t342 * t321;
t316 = qJDD(3) - t322;
t318 = t338 * qJD(2) + t342 * t371;
t331 = qJD(3) - t369;
t218 = (t317 * t331 - t286) * qJ(4) + (t317 * t318 + t316) * pkin(3) + t231;
t232 = t338 * t265 + t342 * t270;
t285 = -t318 * qJD(3) + t342 * qJDD(2) - t338 * t321;
t292 = t331 * pkin(3) - t318 * qJ(4);
t315 = t317 ^ 2;
t221 = -t315 * pkin(3) + t285 * qJ(4) - t331 * t292 + t232;
t336 = sin(pkin(10));
t375 = cos(pkin(10));
t289 = t336 * t317 + t318 * t375;
t205 = t218 * t375 - t336 * t221 + t289 * t380;
t254 = t336 * t285 + t286 * t375;
t294 = -t343 * g(3) - t339 * t305;
t358 = qJDD(2) * pkin(2) + t345 * pkin(8) - t320 * t371 + t294;
t353 = t285 * pkin(3) + t315 * qJ(4) - t318 * t292 - qJDD(4) + t358;
t288 = -t317 * t375 + t336 * t318;
t374 = t288 * t331;
t379 = (-t254 + t374) * qJ(5) - t353;
t253 = -t285 * t375 + t336 * t286;
t275 = -t331 * pkin(5) - t289 * pkin(9);
t287 = t288 ^ 2;
t377 = 2 * qJD(5);
t198 = -t287 * pkin(9) + (-pkin(4) - pkin(5)) * t253 + (-pkin(4) * t331 + t275 + t377) * t289 - t379;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t259 = t337 * t288 + t341 * t289;
t213 = -t259 * qJD(6) + t341 * t253 - t337 * t254;
t258 = t341 * t288 - t337 * t289;
t214 = t258 * qJD(6) + t337 * t253 + t341 * t254;
t329 = qJD(6) - t331;
t236 = -t329 * mrSges(7,2) + t258 * mrSges(7,3);
t237 = t329 * mrSges(7,1) - t259 * mrSges(7,3);
t194 = -m(7) * t198 + t213 * mrSges(7,1) - t214 * mrSges(7,2) + t258 * t236 - t259 * t237;
t208 = -0.2e1 * qJD(5) * t289 + (t289 * t331 + t253) * pkin(4) + t379;
t273 = -t331 * mrSges(6,1) + t289 * mrSges(6,2);
t274 = -t288 * mrSges(6,2) + t331 * mrSges(6,3);
t186 = m(6) * t208 + t253 * mrSges(6,1) - t254 * mrSges(6,3) - t289 * t273 + t288 * t274 + t194;
t206 = t336 * t218 + t375 * t221 + t288 * t380;
t251 = Ifges(6,1) * t289 + Ifges(6,4) * t331 + Ifges(6,5) * t288;
t252 = Ifges(5,1) * t289 - Ifges(5,4) * t288 + Ifges(5,5) * t331;
t260 = t288 * pkin(4) - t289 * qJ(5);
t330 = t331 ^ 2;
t203 = -t316 * pkin(4) - t330 * qJ(5) + t289 * t260 + qJDD(5) - t205;
t195 = (-t254 - t374) * pkin(9) + (t288 * t289 - t316) * pkin(5) + t203;
t201 = -t330 * pkin(4) + t316 * qJ(5) - t288 * t260 + t331 * t377 + t206;
t196 = -t287 * pkin(5) + t253 * pkin(9) + t331 * t275 + t201;
t192 = t341 * t195 - t337 * t196;
t227 = -t258 * mrSges(7,1) + t259 * mrSges(7,2);
t312 = qJDD(6) - t316;
t188 = m(7) * t192 + t312 * mrSges(7,1) - t214 * mrSges(7,3) - t259 * t227 + t329 * t236;
t193 = t337 * t195 + t341 * t196;
t189 = m(7) * t193 - t312 * mrSges(7,2) + t213 * mrSges(7,3) + t258 * t227 - t329 * t237;
t180 = -t337 * t188 + t341 * t189;
t222 = Ifges(7,5) * t259 + Ifges(7,6) * t258 + Ifges(7,3) * t329;
t224 = Ifges(7,1) * t259 + Ifges(7,4) * t258 + Ifges(7,5) * t329;
t182 = -mrSges(7,1) * t198 + mrSges(7,3) * t193 + Ifges(7,4) * t214 + Ifges(7,2) * t213 + Ifges(7,6) * t312 - t259 * t222 + t329 * t224;
t223 = Ifges(7,4) * t259 + Ifges(7,2) * t258 + Ifges(7,6) * t329;
t183 = mrSges(7,2) * t198 - mrSges(7,3) * t192 + Ifges(7,1) * t214 + Ifges(7,4) * t213 + Ifges(7,5) * t312 + t258 * t222 - t329 * t223;
t354 = -mrSges(6,1) * t208 + mrSges(6,2) * t201 - pkin(5) * t194 - pkin(9) * t180 - t341 * t182 - t337 * t183;
t249 = Ifges(6,4) * t289 + Ifges(6,2) * t331 + Ifges(6,6) * t288;
t373 = -Ifges(5,5) * t289 + Ifges(5,6) * t288 - Ifges(5,3) * t331 - t249;
t164 = mrSges(5,1) * t353 + mrSges(5,3) * t206 - pkin(4) * t186 + (t252 + t251) * t331 + (Ifges(5,6) - Ifges(6,6)) * t316 + t373 * t289 + (Ifges(5,4) - Ifges(6,5)) * t254 + (-Ifges(5,2) - Ifges(6,3)) * t253 + t354;
t250 = Ifges(5,4) * t289 - Ifges(5,2) * t288 + Ifges(5,6) * t331;
t179 = t341 * t188 + t337 * t189;
t247 = Ifges(6,5) * t289 + Ifges(6,6) * t331 + Ifges(6,3) * t288;
t355 = mrSges(6,2) * t203 - mrSges(6,3) * t208 + Ifges(6,1) * t254 + Ifges(6,4) * t316 + Ifges(6,5) * t253 - pkin(9) * t179 - t337 * t182 + t341 * t183 + t331 * t247;
t165 = -mrSges(5,2) * t353 - mrSges(5,3) * t205 + Ifges(5,1) * t254 - Ifges(5,4) * t253 + Ifges(5,5) * t316 - qJ(5) * t186 - t331 * t250 + t288 * t373 + t355;
t277 = Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t331;
t279 = Ifges(4,1) * t318 + Ifges(4,4) * t317 + Ifges(4,5) * t331;
t271 = -t331 * mrSges(5,2) - t288 * mrSges(5,3);
t272 = t331 * mrSges(5,1) - t289 * mrSges(5,3);
t351 = -m(5) * t353 + t253 * mrSges(5,1) + t254 * mrSges(5,2) + t288 * t271 + t289 * t272 + t186;
t360 = m(6) * t201 + t316 * mrSges(6,3) + t331 * t273 + t180;
t261 = t288 * mrSges(6,1) - t289 * mrSges(6,3);
t372 = -t288 * mrSges(5,1) - t289 * mrSges(5,2) - t261;
t376 = -mrSges(5,3) - mrSges(6,2);
t172 = m(5) * t206 - t316 * mrSges(5,2) + t253 * t376 - t331 * t272 + t288 * t372 + t360;
t356 = -m(6) * t203 + t316 * mrSges(6,1) + t331 * t274 - t179;
t174 = m(5) * t205 + t316 * mrSges(5,1) + t254 * t376 + t331 * t271 + t289 * t372 + t356;
t364 = t375 * t172 - t336 * t174;
t153 = mrSges(4,1) * t358 + mrSges(4,3) * t232 + Ifges(4,4) * t286 + Ifges(4,2) * t285 + Ifges(4,6) * t316 - pkin(3) * t351 + qJ(4) * t364 + t164 * t375 + t336 * t165 - t318 * t277 + t331 * t279;
t169 = t336 * t172 + t375 * t174;
t278 = Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t331;
t154 = -mrSges(4,2) * t358 - mrSges(4,3) * t231 + Ifges(4,1) * t286 + Ifges(4,4) * t285 + Ifges(4,5) * t316 - qJ(4) * t169 - t336 * t164 + t165 * t375 + t317 * t277 - t331 * t278;
t290 = -t317 * mrSges(4,1) + t318 * mrSges(4,2);
t291 = -t331 * mrSges(4,2) + t317 * mrSges(4,3);
t167 = m(4) * t231 + t316 * mrSges(4,1) - t286 * mrSges(4,3) - t318 * t290 + t331 * t291 + t169;
t293 = t331 * mrSges(4,1) - t318 * mrSges(4,3);
t168 = m(4) * t232 - t316 * mrSges(4,2) + t285 * mrSges(4,3) + t317 * t290 - t331 * t293 + t364;
t163 = -t338 * t167 + t342 * t168;
t185 = m(4) * t358 + t285 * mrSges(4,1) - t286 * mrSges(4,2) + t317 * t291 - t318 * t293 - t351;
t302 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t343) * qJD(1);
t303 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t343) * qJD(1);
t378 = mrSges(3,1) * t294 - mrSges(3,2) * t295 + Ifges(3,5) * t321 + Ifges(3,6) * t322 + Ifges(3,3) * qJDD(2) + pkin(2) * t185 + pkin(8) * t163 + t342 * t153 + t338 * t154 + (t339 * t302 - t343 * t303) * qJD(1);
t319 = (-mrSges(3,1) * t343 + mrSges(3,2) * t339) * qJD(1);
t324 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t371;
t161 = m(3) * t295 - qJDD(2) * mrSges(3,2) + t322 * mrSges(3,3) - qJD(2) * t324 + t319 * t369 + t163;
t325 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t369;
t184 = m(3) * t294 + qJDD(2) * mrSges(3,1) - t321 * mrSges(3,3) + qJD(2) * t325 - t319 * t371 + t185;
t365 = t343 * t161 - t339 * t184;
t162 = t342 * t167 + t338 * t168;
t359 = mrSges(7,1) * t192 - mrSges(7,2) * t193 + Ifges(7,5) * t214 + Ifges(7,6) * t213 + Ifges(7,3) * t312 + t259 * t223 - t258 * t224;
t301 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t343) * qJD(1);
t150 = mrSges(3,2) * t304 - mrSges(3,3) * t294 + Ifges(3,1) * t321 + Ifges(3,4) * t322 + Ifges(3,5) * qJDD(2) - pkin(8) * t162 - qJD(2) * t302 - t338 * t153 + t342 * t154 + t301 * t369;
t350 = mrSges(6,1) * t203 - mrSges(6,3) * t201 - Ifges(6,4) * t254 - Ifges(6,2) * t316 - Ifges(6,6) * t253 + pkin(5) * t179 + t289 * t247 - t288 * t251 + t359;
t348 = mrSges(5,2) * t206 - t288 * t252 - qJ(5) * (-t253 * mrSges(6,2) - t288 * t261 + t360) - pkin(4) * (-t254 * mrSges(6,2) - t289 * t261 + t356) - mrSges(5,1) * t205 - t289 * t250 + Ifges(5,6) * t253 - Ifges(5,5) * t254 - Ifges(5,3) * t316 + t350;
t347 = mrSges(4,1) * t231 - mrSges(4,2) * t232 + Ifges(4,5) * t286 + Ifges(4,6) * t285 + Ifges(4,3) * t316 + pkin(3) * t169 + t318 * t278 - t317 * t279 - t348;
t152 = -mrSges(3,1) * t304 + mrSges(3,3) * t295 + Ifges(3,4) * t321 + Ifges(3,2) * t322 + Ifges(3,6) * qJDD(2) - pkin(2) * t162 + qJD(2) * t303 - t301 * t371 - t347;
t352 = -m(3) * t304 + t322 * mrSges(3,1) - t321 * mrSges(3,2) - t324 * t371 + t325 * t369 - t162;
t357 = mrSges(2,1) * t327 - mrSges(2,2) * t328 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t365 + t339 * t150 + t343 * t152;
t158 = m(2) * t327 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t352;
t157 = t339 * t161 + t343 * t184;
t155 = m(2) * t328 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t365;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t157 - t378;
t147 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t157 + t343 * t150 - t339 * t152;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t344 * t147 - t340 * t148 - pkin(6) * (t340 * t155 + t344 * t158), t147, t150, t154, t165, -t288 * t249 + t355, t183; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t147 + t344 * t148 + pkin(6) * (t344 * t155 - t340 * t158), t148, t152, t153, t164, -t350, t182; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t357, t357, t378, t347, -t348, Ifges(6,5) * t254 + Ifges(6,6) * t316 + Ifges(6,3) * t253 + t289 * t249 - t331 * t251 - t354, t359;];
m_new  = t1;
