% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR3
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
% Datum: 2019-05-07 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:54:47
% EndTime: 2019-05-07 19:55:23
% DurationCPUTime: 18.99s
% Computational Cost: add. (312685->387), mult. (690451->467), div. (0->0), fcn. (509998->10), ass. (0->148)
t337 = sin(qJ(2));
t341 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t314 = qJDD(1) * t337 + t341 * t365;
t338 = sin(qJ(1));
t342 = cos(qJ(1));
t321 = -g(1) * t342 - g(2) * t338;
t343 = qJD(1) ^ 2;
t309 = -pkin(1) * t343 + qJDD(1) * pkin(7) + t321;
t370 = t337 * t309;
t373 = pkin(2) * t343;
t266 = qJDD(2) * pkin(2) - t314 * pkin(8) - t370 + (pkin(8) * t365 + t337 * t373 - g(3)) * t341;
t296 = -g(3) * t337 + t341 * t309;
t315 = qJDD(1) * t341 - t337 * t365;
t367 = qJD(1) * t337;
t319 = qJD(2) * pkin(2) - pkin(8) * t367;
t333 = t341 ^ 2;
t267 = pkin(8) * t315 - qJD(2) * t319 - t333 * t373 + t296;
t336 = sin(qJ(3));
t340 = cos(qJ(3));
t244 = t340 * t266 - t336 * t267;
t306 = (-t336 * t337 + t340 * t341) * qJD(1);
t276 = qJD(3) * t306 + t314 * t340 + t315 * t336;
t307 = (t336 * t341 + t337 * t340) * qJD(1);
t330 = qJDD(2) + qJDD(3);
t331 = qJD(2) + qJD(3);
t210 = (t306 * t331 - t276) * pkin(9) + (t306 * t307 + t330) * pkin(3) + t244;
t245 = t336 * t266 + t340 * t267;
t275 = -qJD(3) * t307 - t314 * t336 + t315 * t340;
t299 = pkin(3) * t331 - pkin(9) * t307;
t302 = t306 ^ 2;
t216 = -pkin(3) * t302 + pkin(9) * t275 - t299 * t331 + t245;
t335 = sin(qJ(4));
t374 = cos(qJ(4));
t207 = t210 * t374 - t335 * t216;
t208 = t335 * t210 + t374 * t216;
t293 = t335 * t306 + t307 * t374;
t237 = qJD(4) * t293 - t275 * t374 + t276 * t335;
t292 = -t306 * t374 + t307 * t335;
t238 = -t292 * qJD(4) + t335 * t275 + t276 * t374;
t328 = qJD(4) + t331;
t252 = Ifges(5,4) * t293 - Ifges(5,2) * t292 + Ifges(5,6) * t328;
t261 = -mrSges(6,2) * t292 - mrSges(6,3) * t293;
t278 = mrSges(6,1) * t292 - mrSges(6,3) * t328;
t327 = qJDD(4) + t330;
t259 = pkin(4) * t292 - qJ(5) * t293;
t326 = t328 ^ 2;
t204 = -t327 * pkin(4) - t326 * qJ(5) + t293 * t259 + qJDD(5) - t207;
t371 = t292 * t328;
t197 = (t292 * t293 - t327) * pkin(10) + (t238 + t371) * pkin(5) + t204;
t282 = pkin(5) * t293 - pkin(10) * t328;
t291 = t292 ^ 2;
t320 = t338 * g(1) - t342 * g(2);
t359 = -qJDD(1) * pkin(1) - t320;
t277 = -t315 * pkin(2) + t319 * t367 + (-pkin(8) * t333 - pkin(7)) * t343 + t359;
t221 = -t275 * pkin(3) - t302 * pkin(9) + t307 * t299 + t277;
t375 = -2 * qJD(5);
t347 = (-t238 + t371) * qJ(5) + t221 + (pkin(4) * t328 + t375) * t293;
t198 = t347 - pkin(5) * t291 - t282 * t293 + (pkin(4) + pkin(10)) * t237;
t334 = sin(qJ(6));
t339 = cos(qJ(6));
t194 = t197 * t339 - t198 * t334;
t270 = t292 * t339 - t328 * t334;
t218 = qJD(6) * t270 + t237 * t334 + t327 * t339;
t235 = qJDD(6) + t238;
t271 = t292 * t334 + t328 * t339;
t246 = -mrSges(7,1) * t270 + mrSges(7,2) * t271;
t288 = qJD(6) + t293;
t247 = -mrSges(7,2) * t288 + mrSges(7,3) * t270;
t191 = m(7) * t194 + mrSges(7,1) * t235 - mrSges(7,3) * t218 - t246 * t271 + t247 * t288;
t195 = t197 * t334 + t198 * t339;
t217 = -qJD(6) * t271 + t237 * t339 - t327 * t334;
t248 = mrSges(7,1) * t288 - mrSges(7,3) * t271;
t192 = m(7) * t195 - mrSges(7,2) * t235 + mrSges(7,3) * t217 + t246 * t270 - t248 * t288;
t179 = t339 * t191 + t334 * t192;
t356 = -t326 * pkin(4) + t327 * qJ(5) - t292 * t259 + t208;
t200 = -t237 * pkin(5) - t291 * pkin(10) + ((2 * qJD(5)) + t282) * t328 + t356;
t223 = Ifges(7,5) * t271 + Ifges(7,6) * t270 + Ifges(7,3) * t288;
t225 = Ifges(7,1) * t271 + Ifges(7,4) * t270 + Ifges(7,5) * t288;
t182 = -mrSges(7,1) * t200 + mrSges(7,3) * t195 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t235 - t223 * t271 + t225 * t288;
t224 = Ifges(7,4) * t271 + Ifges(7,2) * t270 + Ifges(7,6) * t288;
t183 = mrSges(7,2) * t200 - mrSges(7,3) * t194 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t235 + t223 * t270 - t224 * t288;
t201 = t328 * t375 - t356;
t249 = Ifges(6,5) * t328 - Ifges(6,6) * t293 + Ifges(6,3) * t292;
t352 = -mrSges(6,2) * t204 + mrSges(6,3) * t201 - Ifges(6,1) * t327 + Ifges(6,4) * t238 - Ifges(6,5) * t237 + pkin(10) * t179 + t334 * t182 - t339 * t183 + t293 * t249;
t196 = -m(7) * t200 + t217 * mrSges(7,1) - t218 * mrSges(7,2) + t270 * t247 - t271 * t248;
t279 = mrSges(6,1) * t293 + mrSges(6,2) * t328;
t353 = -m(6) * t201 + t327 * mrSges(6,3) + t328 * t279 - t196;
t357 = -m(6) * t204 - t238 * mrSges(6,1) - t293 * t261 - t179;
t251 = Ifges(6,4) * t328 - Ifges(6,2) * t293 + Ifges(6,6) * t292;
t368 = Ifges(5,1) * t293 - Ifges(5,4) * t292 + Ifges(5,5) * t328 - t251;
t379 = -mrSges(5,2) * t208 + pkin(4) * (-t327 * mrSges(6,2) - t328 * t278 + t357) + qJ(5) * (-t237 * mrSges(6,1) - t292 * t261 + t353) + mrSges(5,1) * t207 - Ifges(5,6) * t237 + Ifges(5,5) * t238 + t293 * t252 + Ifges(5,3) * t327 - t352 + t368 * t292;
t260 = mrSges(5,1) * t292 + mrSges(5,2) * t293;
t280 = -mrSges(5,2) * t328 - mrSges(5,3) * t292;
t175 = m(5) * t207 - t238 * mrSges(5,3) - t293 * t260 + (-t278 + t280) * t328 + (mrSges(5,1) - mrSges(6,2)) * t327 + t357;
t281 = mrSges(5,1) * t328 - mrSges(5,3) * t293;
t186 = m(5) * t208 - t327 * mrSges(5,2) - t328 * t281 + (-t260 - t261) * t292 + (-mrSges(5,3) - mrSges(6,1)) * t237 + t353;
t171 = t374 * t175 + t335 * t186;
t286 = Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t331;
t287 = Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t331;
t378 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t276 + Ifges(4,6) * t275 + Ifges(4,3) * t330 + pkin(3) * t171 + t307 * t286 - t306 * t287 + t379;
t294 = -mrSges(4,1) * t306 + mrSges(4,2) * t307;
t297 = -mrSges(4,2) * t331 + mrSges(4,3) * t306;
t168 = m(4) * t244 + mrSges(4,1) * t330 - mrSges(4,3) * t276 - t294 * t307 + t297 * t331 + t171;
t298 = mrSges(4,1) * t331 - mrSges(4,3) * t307;
t361 = -t175 * t335 + t374 * t186;
t169 = m(4) * t245 - mrSges(4,2) * t330 + mrSges(4,3) * t275 + t294 * t306 - t298 * t331 + t361;
t163 = t340 * t168 + t336 * t169;
t295 = -t341 * g(3) - t370;
t304 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t337 + Ifges(3,2) * t341) * qJD(1);
t305 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t337 + Ifges(3,4) * t341) * qJD(1);
t377 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t163 + (t304 * t337 - t305 * t341) * qJD(1) + t378;
t372 = Ifges(5,4) + Ifges(6,6);
t180 = -t334 * t191 + t339 * t192;
t253 = Ifges(6,1) * t328 - Ifges(6,4) * t293 + Ifges(6,5) * t292;
t369 = -Ifges(5,5) * t293 + Ifges(5,6) * t292 - Ifges(5,3) * t328 - t253;
t366 = qJD(1) * t341;
t313 = (-mrSges(3,1) * t341 + mrSges(3,2) * t337) * qJD(1);
t318 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t161 = m(3) * t295 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t314 + qJD(2) * t318 - t313 * t367 + t163;
t317 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t362 = -t168 * t336 + t340 * t169;
t162 = m(3) * t296 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t317 + t313 * t366 + t362;
t363 = -t161 * t337 + t341 * t162;
t205 = t237 * pkin(4) + t347;
t176 = m(6) * t205 - t237 * mrSges(6,2) - t238 * mrSges(6,3) - t292 * t278 - t293 * t279 + t180;
t351 = -mrSges(6,1) * t201 + mrSges(6,2) * t205 - pkin(5) * t196 - pkin(10) * t180 - t339 * t182 - t334 * t183;
t159 = -mrSges(5,1) * t221 + mrSges(5,3) * t208 - pkin(4) * t176 + t368 * t328 + (Ifges(5,6) - Ifges(6,5)) * t327 + t369 * t293 + t372 * t238 + (-Ifges(5,2) - Ifges(6,3)) * t237 + t351;
t355 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t218 + Ifges(7,6) * t217 + Ifges(7,3) * t235 + t271 * t224 - t270 * t225;
t350 = mrSges(6,1) * t204 - mrSges(6,3) * t205 + pkin(5) * t179 + t355;
t164 = t350 + (-t252 + t249) * t328 + (Ifges(5,5) - Ifges(6,4)) * t327 + t369 * t292 + (Ifges(5,1) + Ifges(6,2)) * t238 - t372 * t237 + mrSges(5,2) * t221 - mrSges(5,3) * t207 - qJ(5) * t176;
t285 = Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t331;
t354 = m(5) * t221 + t237 * mrSges(5,1) + t238 * mrSges(5,2) + t292 * t280 + t293 * t281 + t176;
t154 = -mrSges(4,1) * t277 + mrSges(4,3) * t245 + Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * t330 - pkin(3) * t354 + pkin(9) * t361 + t159 * t374 + t335 * t164 - t307 * t285 + t331 * t287;
t155 = mrSges(4,2) * t277 - mrSges(4,3) * t244 + Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * t330 - pkin(9) * t171 - t335 * t159 + t164 * t374 + t306 * t285 - t331 * t286;
t303 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t337 + Ifges(3,6) * t341) * qJD(1);
t308 = -t343 * pkin(7) + t359;
t349 = m(4) * t277 - t275 * mrSges(4,1) + t276 * mrSges(4,2) - t306 * t297 + t307 * t298 + t354;
t150 = -mrSges(3,1) * t308 + mrSges(3,3) * t296 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t349 + pkin(8) * t362 + qJD(2) * t305 + t340 * t154 + t336 * t155 - t303 * t367;
t152 = mrSges(3,2) * t308 - mrSges(3,3) * t295 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t163 - qJD(2) * t304 - t154 * t336 + t155 * t340 + t303 * t366;
t346 = -m(3) * t308 + t315 * mrSges(3,1) - t314 * mrSges(3,2) - t317 * t367 + t318 * t366 - t349;
t358 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t363 + t341 * t150 + t337 * t152;
t172 = m(2) * t320 + qJDD(1) * mrSges(2,1) - t343 * mrSges(2,2) + t346;
t158 = t161 * t341 + t162 * t337;
t156 = m(2) * t321 - mrSges(2,1) * t343 - qJDD(1) * mrSges(2,2) + t363;
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t158 - t377;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t343 - pkin(7) * t158 - t150 * t337 + t152 * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t148 - t338 * t153 - pkin(6) * (t156 * t338 + t172 * t342), t148, t152, t155, t164, -t292 * t251 - t352, t183; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t338 * t148 + t342 * t153 + pkin(6) * (t156 * t342 - t172 * t338), t153, t150, t154, t159, Ifges(6,4) * t327 - Ifges(6,2) * t238 + Ifges(6,6) * t237 - t328 * t249 + t292 * t253 - t350, t182; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t358, t358, t377, t378, t379, Ifges(6,5) * t327 - Ifges(6,6) * t238 + Ifges(6,3) * t237 + t328 * t251 + t293 * t253 - t351, t355;];
m_new  = t1;
