% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP5
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
% Datum: 2019-05-06 01:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:28:21
% EndTime: 2019-05-06 01:28:56
% DurationCPUTime: 19.23s
% Computational Cost: add. (310993->361), mult. (746124->437), div. (0->0), fcn. (581411->10), ass. (0->143)
t330 = qJD(1) ^ 2;
t326 = sin(qJ(1));
t329 = cos(qJ(1));
t304 = -t329 * g(1) - t326 * g(2);
t297 = -t330 * pkin(1) + qJDD(1) * qJ(2) + t304;
t321 = sin(pkin(10));
t322 = cos(pkin(10));
t359 = qJD(1) * qJD(2);
t356 = -t322 * g(3) - 0.2e1 * t321 * t359;
t368 = pkin(2) * t322;
t269 = (-pkin(7) * qJDD(1) + t330 * t368 - t297) * t321 + t356;
t288 = -t321 * g(3) + (t297 + 0.2e1 * t359) * t322;
t358 = qJDD(1) * t322;
t317 = t322 ^ 2;
t365 = t317 * t330;
t270 = -pkin(2) * t365 + pkin(7) * t358 + t288;
t325 = sin(qJ(3));
t328 = cos(qJ(3));
t248 = t328 * t269 - t325 * t270;
t344 = t321 * t328 + t322 * t325;
t343 = -t321 * t325 + t322 * t328;
t295 = t343 * qJD(1);
t360 = t295 * qJD(3);
t286 = qJDD(1) * t344 + t360;
t296 = t344 * qJD(1);
t207 = (-t286 + t360) * pkin(8) + (t295 * t296 + qJDD(3)) * pkin(3) + t248;
t249 = t325 * t269 + t328 * t270;
t285 = -t296 * qJD(3) + qJDD(1) * t343;
t291 = qJD(3) * pkin(3) - t296 * pkin(8);
t294 = t295 ^ 2;
t218 = -t294 * pkin(3) + t285 * pkin(8) - qJD(3) * t291 + t249;
t324 = sin(qJ(4));
t327 = cos(qJ(4));
t203 = t324 * t207 + t327 * t218;
t276 = t327 * t295 - t324 * t296;
t277 = t324 * t295 + t327 * t296;
t260 = -t276 * pkin(4) - t277 * pkin(9);
t318 = qJD(3) + qJD(4);
t314 = t318 ^ 2;
t315 = qJDD(3) + qJDD(4);
t198 = -t314 * pkin(4) + t315 * pkin(9) + t276 * t260 + t203;
t316 = t321 ^ 2;
t303 = t326 * g(1) - t329 * g(2);
t350 = qJDD(2) - t303;
t284 = (-pkin(1) - t368) * qJDD(1) + (-qJ(2) + (-t316 - t317) * pkin(7)) * t330 + t350;
t231 = -t285 * pkin(3) - t294 * pkin(8) + t296 * t291 + t284;
t242 = -t277 * qJD(4) + t327 * t285 - t324 * t286;
t243 = t276 * qJD(4) + t324 * t285 + t327 * t286;
t200 = (-t276 * t318 - t243) * pkin(9) + (t277 * t318 - t242) * pkin(4) + t231;
t323 = sin(qJ(5));
t369 = cos(qJ(5));
t194 = -t323 * t198 + t369 * t200;
t195 = t369 * t198 + t323 * t200;
t263 = t369 * t277 + t323 * t318;
t214 = t263 * qJD(5) + t323 * t243 - t369 * t315;
t262 = t323 * t277 - t369 * t318;
t215 = -t262 * qJD(5) + t369 * t243 + t323 * t315;
t272 = qJD(5) - t276;
t219 = Ifges(7,5) * t263 + Ifges(7,6) * t272 + Ifges(7,3) * t262;
t222 = Ifges(6,4) * t263 - Ifges(6,2) * t262 + Ifges(6,6) * t272;
t224 = Ifges(6,1) * t263 - Ifges(6,4) * t262 + Ifges(6,5) * t272;
t241 = qJDD(5) - t242;
t245 = t262 * mrSges(7,1) - t263 * mrSges(7,3);
t244 = t262 * pkin(5) - t263 * qJ(6);
t271 = t272 ^ 2;
t190 = -t271 * pkin(5) + t241 * qJ(6) + 0.2e1 * qJD(6) * t272 - t262 * t244 + t195;
t192 = -t241 * pkin(5) - t271 * qJ(6) + t263 * t244 + qJDD(6) - t194;
t223 = Ifges(7,1) * t263 + Ifges(7,4) * t272 + Ifges(7,5) * t262;
t340 = mrSges(7,1) * t192 - mrSges(7,3) * t190 - Ifges(7,4) * t215 - Ifges(7,2) * t241 - Ifges(7,6) * t214 - t262 * t223;
t250 = -t262 * mrSges(7,2) + t272 * mrSges(7,3);
t351 = -m(7) * t192 + t241 * mrSges(7,1) + t272 * t250;
t253 = -t272 * mrSges(7,1) + t263 * mrSges(7,2);
t357 = m(7) * t190 + t241 * mrSges(7,3) + t272 * t253;
t371 = -(-t222 + t219) * t263 + mrSges(6,1) * t194 - mrSges(6,2) * t195 + Ifges(6,5) * t215 - Ifges(6,6) * t214 + Ifges(6,3) * t241 + pkin(5) * (-t215 * mrSges(7,2) - t263 * t245 + t351) + qJ(6) * (-t214 * mrSges(7,2) - t262 * t245 + t357) + t262 * t224 - t340;
t259 = -t276 * mrSges(5,1) + t277 * mrSges(5,2);
t267 = t318 * mrSges(5,1) - t277 * mrSges(5,3);
t252 = t272 * mrSges(6,1) - t263 * mrSges(6,3);
t362 = -t262 * mrSges(6,1) - t263 * mrSges(6,2) - t245;
t367 = -mrSges(6,3) - mrSges(7,2);
t180 = m(6) * t195 - t241 * mrSges(6,2) + t367 * t214 - t272 * t252 + t362 * t262 + t357;
t251 = -t272 * mrSges(6,2) - t262 * mrSges(6,3);
t182 = m(6) * t194 + t241 * mrSges(6,1) + t367 * t215 + t272 * t251 + t362 * t263 + t351;
t352 = t369 * t180 - t323 * t182;
t168 = m(5) * t203 - t315 * mrSges(5,2) + t242 * mrSges(5,3) + t276 * t259 - t318 * t267 + t352;
t202 = t327 * t207 - t324 * t218;
t266 = -t318 * mrSges(5,2) + t276 * mrSges(5,3);
t197 = -t315 * pkin(4) - t314 * pkin(9) + t277 * t260 - t202;
t193 = -0.2e1 * qJD(6) * t263 + (t262 * t272 - t215) * qJ(6) + (t263 * t272 + t214) * pkin(5) + t197;
t187 = m(7) * t193 + t214 * mrSges(7,1) - t215 * mrSges(7,3) + t262 * t250 - t263 * t253;
t335 = -m(6) * t197 - t214 * mrSges(6,1) - t215 * mrSges(6,2) - t262 * t251 - t263 * t252 - t187;
t177 = m(5) * t202 + t315 * mrSges(5,1) - t243 * mrSges(5,3) - t277 * t259 + t318 * t266 + t335;
t163 = t324 * t168 + t327 * t177;
t280 = -t295 * mrSges(4,1) + t296 * mrSges(4,2);
t289 = -qJD(3) * mrSges(4,2) + t295 * mrSges(4,3);
t160 = m(4) * t248 + qJDD(3) * mrSges(4,1) - t286 * mrSges(4,3) + qJD(3) * t289 - t296 * t280 + t163;
t290 = qJD(3) * mrSges(4,1) - t296 * mrSges(4,3);
t353 = t327 * t168 - t324 * t177;
t161 = m(4) * t249 - qJDD(3) * mrSges(4,2) + t285 * mrSges(4,3) - qJD(3) * t290 + t295 * t280 + t353;
t154 = t328 * t160 + t325 * t161;
t287 = -t321 * t297 + t356;
t274 = Ifges(4,4) * t296 + Ifges(4,2) * t295 + Ifges(4,6) * qJD(3);
t275 = Ifges(4,1) * t296 + Ifges(4,4) * t295 + Ifges(4,5) * qJD(3);
t349 = -mrSges(7,1) * t193 + mrSges(7,2) * t190;
t221 = Ifges(7,4) * t263 + Ifges(7,2) * t272 + Ifges(7,6) * t262;
t364 = -Ifges(6,5) * t263 + Ifges(6,6) * t262 - Ifges(6,3) * t272 - t221;
t170 = -mrSges(6,1) * t197 + mrSges(6,3) * t195 - pkin(5) * t187 + (t223 + t224) * t272 + t364 * t263 + (Ifges(6,6) - Ifges(7,6)) * t241 + (Ifges(6,4) - Ifges(7,5)) * t215 + (-Ifges(6,2) - Ifges(7,3)) * t214 + t349;
t339 = mrSges(7,2) * t192 - mrSges(7,3) * t193 + Ifges(7,1) * t215 + Ifges(7,4) * t241 + Ifges(7,5) * t214 + t272 * t219;
t172 = mrSges(6,2) * t197 - mrSges(6,3) * t194 + Ifges(6,1) * t215 - Ifges(6,4) * t214 + Ifges(6,5) * t241 - qJ(6) * t187 - t272 * t222 + t262 * t364 + t339;
t255 = Ifges(5,4) * t277 + Ifges(5,2) * t276 + Ifges(5,6) * t318;
t256 = Ifges(5,1) * t277 + Ifges(5,4) * t276 + Ifges(5,5) * t318;
t337 = -mrSges(5,1) * t202 + mrSges(5,2) * t203 - Ifges(5,5) * t243 - Ifges(5,6) * t242 - Ifges(5,3) * t315 - pkin(4) * t335 - pkin(9) * t352 - t369 * t170 - t323 * t172 - t277 * t255 + t276 * t256;
t333 = -mrSges(4,1) * t248 + mrSges(4,2) * t249 - Ifges(4,5) * t286 - Ifges(4,6) * t285 - Ifges(4,3) * qJDD(3) - pkin(3) * t163 - t296 * t274 + t295 * t275 + t337;
t347 = Ifges(3,4) * t321 + Ifges(3,2) * t322;
t348 = Ifges(3,1) * t321 + Ifges(3,4) * t322;
t370 = -mrSges(3,1) * t287 + mrSges(3,2) * t288 - pkin(2) * t154 - (t321 * t347 - t322 * t348) * t330 + t333;
t366 = mrSges(3,2) * t321;
t174 = t323 * t180 + t369 * t182;
t346 = Ifges(3,5) * t321 + Ifges(3,6) * t322;
t361 = t330 * t346;
t342 = mrSges(3,3) * qJDD(1) + t330 * (-mrSges(3,1) * t322 + t366);
t152 = m(3) * t287 - t321 * t342 + t154;
t354 = -t325 * t160 + t328 * t161;
t153 = m(3) * t288 + t322 * t342 + t354;
t355 = -t321 * t152 + t322 * t153;
t341 = m(5) * t231 - t242 * mrSges(5,1) + t243 * mrSges(5,2) - t276 * t266 + t277 * t267 + t174;
t254 = Ifges(5,5) * t277 + Ifges(5,6) * t276 + Ifges(5,3) * t318;
t155 = mrSges(5,2) * t231 - mrSges(5,3) * t202 + Ifges(5,1) * t243 + Ifges(5,4) * t242 + Ifges(5,5) * t315 - pkin(9) * t174 - t323 * t170 + t369 * t172 + t276 * t254 - t318 * t255;
t156 = -mrSges(5,1) * t231 + mrSges(5,3) * t203 + Ifges(5,4) * t243 + Ifges(5,2) * t242 + Ifges(5,6) * t315 - pkin(4) * t174 - t277 * t254 + t318 * t256 - t371;
t273 = Ifges(4,5) * t296 + Ifges(4,6) * t295 + Ifges(4,3) * qJD(3);
t146 = -mrSges(4,1) * t284 + mrSges(4,3) * t249 + Ifges(4,4) * t286 + Ifges(4,2) * t285 + Ifges(4,6) * qJDD(3) - pkin(3) * t341 + pkin(8) * t353 + qJD(3) * t275 + t324 * t155 + t327 * t156 - t296 * t273;
t147 = mrSges(4,2) * t284 - mrSges(4,3) * t248 + Ifges(4,1) * t286 + Ifges(4,4) * t285 + Ifges(4,5) * qJDD(3) - pkin(8) * t163 - qJD(3) * t274 + t327 * t155 - t324 * t156 + t295 * t273;
t293 = -qJDD(1) * pkin(1) - t330 * qJ(2) + t350;
t336 = m(4) * t284 - t285 * mrSges(4,1) + t286 * mrSges(4,2) - t295 * t289 + t296 * t290 + t341;
t142 = -mrSges(3,1) * t293 + mrSges(3,3) * t288 - pkin(2) * t336 + pkin(7) * t354 + qJDD(1) * t347 + t328 * t146 + t325 * t147 - t321 * t361;
t144 = mrSges(3,2) * t293 - mrSges(3,3) * t287 - pkin(7) * t154 + qJDD(1) * t348 - t325 * t146 + t328 * t147 + t322 * t361;
t334 = -m(3) * t293 + mrSges(3,1) * t358 - t336 + (t316 * t330 + t365) * mrSges(3,3);
t338 = -mrSges(2,2) * t304 + qJ(2) * t355 + t322 * t142 + t321 * t144 + pkin(1) * (-qJDD(1) * t366 + t334) + mrSges(2,1) * t303 + Ifges(2,3) * qJDD(1);
t164 = m(2) * t303 - t330 * mrSges(2,2) + t334 + (mrSges(2,1) - t366) * qJDD(1);
t150 = t322 * t152 + t321 * t153;
t148 = m(2) * t304 - t330 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t355;
t145 = mrSges(2,3) * t304 + mrSges(2,1) * g(3) + t330 * Ifges(2,5) - pkin(1) * t150 + (Ifges(2,6) - t346) * qJDD(1) + t370;
t140 = -mrSges(2,2) * g(3) - mrSges(2,3) * t303 + Ifges(2,5) * qJDD(1) - t330 * Ifges(2,6) - qJ(2) * t150 - t321 * t142 + t322 * t144;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t329 * t140 - t326 * t145 - pkin(6) * (t326 * t148 + t329 * t164), t140, t144, t147, t155, t172, -t262 * t221 + t339; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t326 * t140 + t329 * t145 + pkin(6) * (t329 * t148 - t326 * t164), t145, t142, t146, t156, t170, -t263 * t219 - t340; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t338, t338, qJDD(1) * t346 - t370, -t333, -t337, t371, Ifges(7,5) * t215 + Ifges(7,6) * t241 + Ifges(7,3) * t214 + t263 * t221 - t272 * t223 - t349;];
m_new  = t1;
