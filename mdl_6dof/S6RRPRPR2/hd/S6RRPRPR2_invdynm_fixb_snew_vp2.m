% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 13:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:01:27
% EndTime: 2019-05-06 13:02:09
% DurationCPUTime: 17.46s
% Computational Cost: add. (269890->386), mult. (634633->469), div. (0->0), fcn. (463022->10), ass. (0->146)
t338 = sin(qJ(2));
t341 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t316 = t338 * qJDD(1) + t341 * t365;
t339 = sin(qJ(1));
t342 = cos(qJ(1));
t323 = -t342 * g(1) - t339 * g(2);
t343 = qJD(1) ^ 2;
t311 = -t343 * pkin(1) + qJDD(1) * pkin(7) + t323;
t370 = t338 * t311;
t373 = pkin(2) * t343;
t266 = qJDD(2) * pkin(2) - t316 * qJ(3) - t370 + (qJ(3) * t365 + t338 * t373 - g(3)) * t341;
t296 = -t338 * g(3) + t341 * t311;
t317 = t341 * qJDD(1) - t338 * t365;
t367 = qJD(1) * t338;
t319 = qJD(2) * pkin(2) - qJ(3) * t367;
t333 = t341 ^ 2;
t267 = t317 * qJ(3) - qJD(2) * t319 - t333 * t373 + t296;
t334 = sin(pkin(10));
t335 = cos(pkin(10));
t306 = (t334 * t341 + t335 * t338) * qJD(1);
t227 = -0.2e1 * qJD(3) * t306 + t335 * t266 - t334 * t267;
t294 = t335 * t316 + t334 * t317;
t305 = (-t334 * t338 + t335 * t341) * qJD(1);
t210 = (qJD(2) * t305 - t294) * pkin(8) + (t305 * t306 + qJDD(2)) * pkin(3) + t227;
t228 = 0.2e1 * qJD(3) * t305 + t334 * t266 + t335 * t267;
t293 = -t334 * t316 + t335 * t317;
t299 = qJD(2) * pkin(3) - t306 * pkin(8);
t304 = t305 ^ 2;
t213 = -t304 * pkin(3) + t293 * pkin(8) - qJD(2) * t299 + t228;
t337 = sin(qJ(4));
t374 = cos(qJ(4));
t207 = t210 * t374 - t337 * t213;
t208 = t337 * t210 + t374 * t213;
t285 = t337 * t305 + t306 * t374;
t244 = t285 * qJD(4) - t293 * t374 + t337 * t294;
t284 = -t305 * t374 + t337 * t306;
t245 = -t284 * qJD(4) + t337 * t293 + t294 * t374;
t330 = qJD(2) + qJD(4);
t252 = Ifges(5,4) * t285 - Ifges(5,2) * t284 + Ifges(5,6) * t330;
t261 = -t284 * mrSges(6,2) - t285 * mrSges(6,3);
t274 = t284 * mrSges(6,1) - t330 * mrSges(6,3);
t329 = qJDD(2) + qJDD(4);
t259 = t284 * pkin(4) - t285 * qJ(5);
t328 = t330 ^ 2;
t203 = -t329 * pkin(4) - t328 * qJ(5) + t285 * t259 + qJDD(5) - t207;
t371 = t284 * t330;
t197 = (t284 * t285 - t329) * pkin(9) + (t245 + t371) * pkin(5) + t203;
t276 = t285 * pkin(5) - t330 * pkin(9);
t280 = t284 ^ 2;
t322 = t339 * g(1) - t342 * g(2);
t359 = -qJDD(1) * pkin(1) - t322;
t269 = -t317 * pkin(2) + qJDD(3) + t319 * t367 + (-qJ(3) * t333 - pkin(7)) * t343 + t359;
t225 = -t293 * pkin(3) - t304 * pkin(8) + t306 * t299 + t269;
t375 = -2 * qJD(5);
t347 = (-t245 + t371) * qJ(5) + t225 + (t330 * pkin(4) + t375) * t285;
t200 = -t285 * t276 - t280 * pkin(5) + (pkin(4) + pkin(9)) * t244 + t347;
t336 = sin(qJ(6));
t340 = cos(qJ(6));
t194 = t340 * t197 - t336 * t200;
t270 = t340 * t284 - t336 * t330;
t219 = t270 * qJD(6) + t336 * t244 + t340 * t329;
t243 = qJDD(6) + t245;
t271 = t336 * t284 + t340 * t330;
t246 = -mrSges(7,1) * t270 + mrSges(7,2) * t271;
t279 = qJD(6) + t285;
t247 = -t279 * mrSges(7,2) + t270 * mrSges(7,3);
t191 = m(7) * t194 + t243 * mrSges(7,1) - t219 * mrSges(7,3) - t271 * t246 + t279 * t247;
t195 = t336 * t197 + t340 * t200;
t218 = -t271 * qJD(6) + t340 * t244 - t336 * t329;
t248 = t279 * mrSges(7,1) - t271 * mrSges(7,3);
t192 = m(7) * t195 - t243 * mrSges(7,2) + t218 * mrSges(7,3) + t270 * t246 - t279 * t248;
t179 = t340 * t191 + t336 * t192;
t356 = -t328 * pkin(4) + t329 * qJ(5) - t284 * t259 + t208;
t199 = -t244 * pkin(5) - t280 * pkin(9) + ((2 * qJD(5)) + t276) * t330 + t356;
t220 = Ifges(7,5) * t271 + Ifges(7,6) * t270 + Ifges(7,3) * t279;
t222 = Ifges(7,1) * t271 + Ifges(7,4) * t270 + Ifges(7,5) * t279;
t182 = -mrSges(7,1) * t199 + mrSges(7,3) * t195 + Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t243 - t271 * t220 + t279 * t222;
t221 = Ifges(7,4) * t271 + Ifges(7,2) * t270 + Ifges(7,6) * t279;
t183 = mrSges(7,2) * t199 - mrSges(7,3) * t194 + Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t243 + t270 * t220 - t279 * t221;
t201 = t330 * t375 - t356;
t249 = Ifges(6,5) * t330 - Ifges(6,6) * t285 + Ifges(6,3) * t284;
t352 = -mrSges(6,2) * t203 + mrSges(6,3) * t201 - Ifges(6,1) * t329 + Ifges(6,4) * t245 - Ifges(6,5) * t244 + pkin(9) * t179 + t336 * t182 - t340 * t183 + t285 * t249;
t196 = -m(7) * t199 + t218 * mrSges(7,1) - t219 * mrSges(7,2) + t270 * t247 - t271 * t248;
t275 = t285 * mrSges(6,1) + t330 * mrSges(6,2);
t353 = -m(6) * t201 + t329 * mrSges(6,3) + t330 * t275 - t196;
t357 = -m(6) * t203 - t245 * mrSges(6,1) - t285 * t261 - t179;
t251 = Ifges(6,4) * t330 - Ifges(6,2) * t285 + Ifges(6,6) * t284;
t368 = Ifges(5,1) * t285 - Ifges(5,4) * t284 + Ifges(5,5) * t330 - t251;
t379 = -mrSges(5,2) * t208 + pkin(4) * (-t329 * mrSges(6,2) - t330 * t274 + t357) + qJ(5) * (-t244 * mrSges(6,1) - t284 * t261 + t353) + mrSges(5,1) * t207 + t285 * t252 - Ifges(5,6) * t244 + Ifges(5,5) * t245 + Ifges(5,3) * t329 - t352 + t368 * t284;
t260 = t284 * mrSges(5,1) + t285 * mrSges(5,2);
t272 = -t330 * mrSges(5,2) - t284 * mrSges(5,3);
t175 = m(5) * t207 - t245 * mrSges(5,3) - t285 * t260 + (t272 - t274) * t330 + (mrSges(5,1) - mrSges(6,2)) * t329 + t357;
t273 = t330 * mrSges(5,1) - t285 * mrSges(5,3);
t186 = m(5) * t208 - t329 * mrSges(5,2) - t330 * t273 + (-t260 - t261) * t284 + (-mrSges(5,3) - mrSges(6,1)) * t244 + t353;
t171 = t374 * t175 + t337 * t186;
t282 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * qJD(2);
t283 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * qJD(2);
t378 = mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t294 + Ifges(4,6) * t293 + Ifges(4,3) * qJDD(2) + pkin(3) * t171 + t306 * t282 - t305 * t283 + t379;
t288 = -t305 * mrSges(4,1) + t306 * mrSges(4,2);
t297 = -qJD(2) * mrSges(4,2) + t305 * mrSges(4,3);
t168 = m(4) * t227 + qJDD(2) * mrSges(4,1) - t294 * mrSges(4,3) + qJD(2) * t297 - t306 * t288 + t171;
t298 = qJD(2) * mrSges(4,1) - t306 * mrSges(4,3);
t361 = -t337 * t175 + t374 * t186;
t169 = m(4) * t228 - qJDD(2) * mrSges(4,2) + t293 * mrSges(4,3) - qJD(2) * t298 + t305 * t288 + t361;
t163 = t335 * t168 + t334 * t169;
t295 = -t341 * g(3) - t370;
t308 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t341) * qJD(1);
t309 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t341) * qJD(1);
t377 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * t163 + (t338 * t308 - t341 * t309) * qJD(1) + t378;
t372 = Ifges(5,4) + Ifges(6,6);
t180 = -t336 * t191 + t340 * t192;
t253 = Ifges(6,1) * t330 - Ifges(6,4) * t285 + Ifges(6,5) * t284;
t369 = -Ifges(5,5) * t285 + Ifges(5,6) * t284 - Ifges(5,3) * t330 - t253;
t366 = qJD(1) * t341;
t315 = (-mrSges(3,1) * t341 + mrSges(3,2) * t338) * qJD(1);
t321 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t161 = m(3) * t295 + qJDD(2) * mrSges(3,1) - t316 * mrSges(3,3) + qJD(2) * t321 - t315 * t367 + t163;
t320 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t362 = -t334 * t168 + t335 * t169;
t162 = m(3) * t296 - qJDD(2) * mrSges(3,2) + t317 * mrSges(3,3) - qJD(2) * t320 + t315 * t366 + t362;
t363 = -t338 * t161 + t341 * t162;
t205 = t244 * pkin(4) + t347;
t176 = m(6) * t205 - t244 * mrSges(6,2) - t245 * mrSges(6,3) - t284 * t274 - t285 * t275 + t180;
t351 = -mrSges(6,1) * t201 + mrSges(6,2) * t205 - pkin(5) * t196 - pkin(9) * t180 - t340 * t182 - t336 * t183;
t159 = -mrSges(5,1) * t225 + mrSges(5,3) * t208 - pkin(4) * t176 + t368 * t330 + (Ifges(5,6) - Ifges(6,5)) * t329 + t369 * t285 + t372 * t245 + (-Ifges(5,2) - Ifges(6,3)) * t244 + t351;
t355 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t219 + Ifges(7,6) * t218 + Ifges(7,3) * t243 + t271 * t221 - t270 * t222;
t350 = mrSges(6,1) * t203 - mrSges(6,3) * t205 + pkin(5) * t179 + t355;
t164 = (-t252 + t249) * t330 + (Ifges(5,5) - Ifges(6,4)) * t329 + t369 * t284 + (Ifges(5,1) + Ifges(6,2)) * t245 - t372 * t244 + mrSges(5,2) * t225 - mrSges(5,3) * t207 - qJ(5) * t176 + t350;
t281 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * qJD(2);
t354 = m(5) * t225 + t244 * mrSges(5,1) + t245 * mrSges(5,2) + t284 * t272 + t285 * t273 + t176;
t154 = -mrSges(4,1) * t269 + mrSges(4,3) * t228 + Ifges(4,4) * t294 + Ifges(4,2) * t293 + Ifges(4,6) * qJDD(2) - pkin(3) * t354 + pkin(8) * t361 + qJD(2) * t283 + t159 * t374 + t337 * t164 - t306 * t281;
t155 = mrSges(4,2) * t269 - mrSges(4,3) * t227 + Ifges(4,1) * t294 + Ifges(4,4) * t293 + Ifges(4,5) * qJDD(2) - pkin(8) * t171 - qJD(2) * t282 - t337 * t159 + t164 * t374 + t305 * t281;
t307 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t341) * qJD(1);
t310 = -t343 * pkin(7) + t359;
t349 = m(4) * t269 - t293 * mrSges(4,1) + t294 * mrSges(4,2) - t305 * t297 + t306 * t298 + t354;
t150 = -mrSges(3,1) * t310 + mrSges(3,3) * t296 + Ifges(3,4) * t316 + Ifges(3,2) * t317 + Ifges(3,6) * qJDD(2) - pkin(2) * t349 + qJ(3) * t362 + qJD(2) * t309 + t335 * t154 + t334 * t155 - t307 * t367;
t152 = mrSges(3,2) * t310 - mrSges(3,3) * t295 + Ifges(3,1) * t316 + Ifges(3,4) * t317 + Ifges(3,5) * qJDD(2) - qJ(3) * t163 - qJD(2) * t308 - t334 * t154 + t335 * t155 + t307 * t366;
t346 = -m(3) * t310 + t317 * mrSges(3,1) - t316 * mrSges(3,2) - t320 * t367 + t321 * t366 - t349;
t358 = mrSges(2,1) * t322 - mrSges(2,2) * t323 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t363 + t341 * t150 + t338 * t152;
t172 = m(2) * t322 + qJDD(1) * mrSges(2,1) - t343 * mrSges(2,2) + t346;
t158 = t341 * t161 + t338 * t162;
t156 = m(2) * t323 - t343 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t363;
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t323 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t158 - t377;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - t343 * Ifges(2,6) - pkin(7) * t158 - t338 * t150 + t341 * t152;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t148 - t339 * t153 - pkin(6) * (t339 * t156 + t342 * t172), t148, t152, t155, t164, -t284 * t251 - t352, t183; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t148 + t342 * t153 + pkin(6) * (t342 * t156 - t339 * t172), t153, t150, t154, t159, Ifges(6,4) * t329 - Ifges(6,2) * t245 + Ifges(6,6) * t244 - t330 * t249 + t284 * t253 - t350, t182; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t358, t358, t377, t378, t379, Ifges(6,5) * t329 - Ifges(6,6) * t245 + Ifges(6,3) * t244 + t330 * t251 + t285 * t253 - t351, t355;];
m_new  = t1;
