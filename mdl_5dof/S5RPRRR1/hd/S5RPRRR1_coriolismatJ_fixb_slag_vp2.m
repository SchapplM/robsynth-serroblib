% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:40
% EndTime: 2019-12-05 18:08:51
% DurationCPUTime: 3.81s
% Computational Cost: add. (3103->420), mult. (9286->642), div. (0->0), fcn. (8162->6), ass. (0->246)
t248 = sin(qJ(4));
t242 = t248 ^ 2;
t251 = cos(qJ(4));
t245 = t251 ^ 2;
t416 = t245 + t242;
t298 = t416 * mrSges(5,3);
t252 = cos(qJ(3));
t327 = t248 * t252;
t307 = t327 / 0.2e1;
t423 = Ifges(6,3) * t307;
t247 = sin(qJ(5));
t250 = cos(qJ(5));
t291 = mrSges(6,1) * t250 - mrSges(6,2) * t247;
t345 = qJ(2) * t291;
t249 = sin(qJ(3));
t321 = t252 * t247;
t181 = t249 * t250 - t251 * t321;
t322 = t251 * t252;
t182 = t247 * t249 + t250 * t322;
t324 = t250 * t251;
t180 = t249 * t324 - t321;
t158 = t180 * qJ(2);
t337 = t158 * t250;
t326 = t249 * t251;
t280 = t247 * t326 + t250 * t252;
t156 = t280 * qJ(2);
t340 = t156 * t247;
t395 = -t250 / 0.4e1;
t401 = -t247 / 0.4e1;
t238 = Ifges(6,4) * t250;
t212 = Ifges(6,1) * t247 + t238;
t403 = -t212 / 0.4e1;
t382 = Ifges(6,4) * t247;
t208 = Ifges(6,2) * t250 + t382;
t404 = -t208 / 0.4e1;
t43 = Ifges(6,4) * t182 + Ifges(6,2) * t181 + Ifges(6,6) * t327;
t45 = Ifges(6,1) * t182 + Ifges(6,4) * t181 + Ifges(6,5) * t327;
t422 = (t340 / 0.2e1 + t337 / 0.2e1) * mrSges(6,3) - Ifges(5,3) * t249 / 0.2e1 + t181 * t404 + t182 * t403 + t45 * t401 + t43 * t395;
t389 = t252 / 0.2e1;
t306 = t326 / 0.2e1;
t421 = qJ(2) * mrSges(5,2);
t358 = t250 * mrSges(6,2);
t363 = t247 * mrSges(6,1);
t290 = t358 + t363;
t420 = t248 * t290;
t199 = -t252 * mrSges(5,1) - mrSges(5,3) * t326;
t57 = mrSges(6,1) * t280 + mrSges(6,2) * t180;
t346 = t199 - t57;
t419 = t346 * qJ(2);
t418 = -Ifges(6,2) * t247 + t238;
t239 = Ifges(5,4) * t251;
t417 = -Ifges(5,2) * t248 + t239;
t214 = Ifges(5,1) * t248 + t239;
t414 = -Ifges(6,6) * t181 / 0.2e1 - Ifges(6,5) * t182 / 0.2e1;
t241 = t247 ^ 2;
t244 = t250 ^ 2;
t320 = t241 + t244;
t413 = mrSges(6,3) * t320 - mrSges(5,2);
t373 = Ifges(6,6) * t247;
t379 = Ifges(6,5) * t250;
t207 = -t373 + t379;
t370 = Ifges(6,3) * t251;
t140 = t248 * t207 - t370;
t213 = Ifges(6,1) * t250 - t382;
t378 = Ifges(6,5) * t251;
t146 = t248 * t213 - t378;
t372 = Ifges(6,6) * t251;
t142 = t248 * t418 - t372;
t412 = m(6) / 0.2e1;
t411 = -t280 / 0.2e1;
t410 = t180 / 0.2e1;
t409 = t181 / 0.2e1;
t408 = t182 / 0.2e1;
t407 = t420 / 0.2e1;
t331 = t247 * t248;
t353 = t251 * mrSges(6,2);
t284 = mrSges(6,3) * t331 - t353;
t406 = -t284 / 0.2e1;
t328 = t248 * t250;
t355 = t251 * mrSges(6,1);
t285 = mrSges(6,3) * t328 + t355;
t405 = -t285 / 0.2e1;
t402 = -t247 / 0.2e1;
t400 = t247 / 0.2e1;
t399 = t248 / 0.2e1;
t398 = t249 / 0.2e1;
t397 = t249 / 0.4e1;
t396 = -t250 / 0.2e1;
t394 = t250 / 0.2e1;
t393 = t250 / 0.4e1;
t392 = -t251 / 0.2e1;
t391 = -t251 / 0.4e1;
t390 = t251 / 0.2e1;
t282 = -t358 / 0.2e1 - t363 / 0.2e1;
t276 = t282 * t251;
t279 = -t284 * t396 - t285 * t400;
t22 = (t276 - t279) * t251 + (t282 * t248 + t407) * t248;
t25 = (t244 / 0.2e1 + t241 / 0.2e1) * t242 * mrSges(6,3) + (-t284 * t400 - t285 * t394 + t291 * t390) * t248;
t388 = t22 * qJD(4) - t25 * qJD(5);
t387 = m(6) * qJD(2);
t384 = Ifges(5,4) * t248;
t383 = Ifges(6,4) * t180;
t371 = Ifges(6,3) * t248;
t369 = qJ(2) * mrSges(5,1);
t157 = t181 * qJ(2);
t159 = t182 * qJ(2);
t354 = t251 * mrSges(5,2);
t184 = t249 * (t248 * mrSges(5,1) + t354);
t329 = t248 * t249;
t352 = t252 * mrSges(5,2);
t197 = -mrSges(5,3) * t329 + t352;
t253 = qJ(2) ^ 2;
t350 = t252 * Ifges(5,6);
t144 = t249 * t417 - t350;
t366 = t159 * mrSges(6,2);
t367 = t157 * mrSges(6,1);
t40 = Ifges(6,5) * t180 - Ifges(6,6) * t280 + Ifges(6,3) * t329;
t264 = t367 - t144 / 0.2e1 + t40 / 0.2e1 - t366 + t350 / 0.2e1;
t281 = t423 - t417 * t252 / 0.2e1 - Ifges(5,6) * t249 - t414;
t364 = t182 * mrSges(6,2);
t365 = t181 * mrSges(6,1);
t292 = t364 - t365;
t215 = Ifges(5,1) * t251 - t384;
t294 = Ifges(5,5) * t249 + t215 * t389;
t351 = t252 * Ifges(5,5);
t148 = t249 * t215 - t351;
t295 = t148 / 0.2e1 - t351 / 0.2e1;
t42 = -Ifges(6,2) * t280 + Ifges(6,6) * t329 + t383;
t162 = Ifges(6,4) * t280;
t44 = Ifges(6,1) * t180 + Ifges(6,5) * t329 - t162;
t97 = -mrSges(6,2) * t329 - mrSges(6,3) * t280;
t98 = mrSges(6,1) * t329 - t180 * mrSges(6,3);
t1 = t156 * t98 - t158 * t97 + t43 * t411 + t45 * t410 + t42 * t409 + t44 * t408 + (-t157 * t182 + t159 * t181) * mrSges(6,3) + (qJ(2) * t184 + Ifges(4,4) * t252 + t295 * t251 + (qJ(2) * t292 + t264) * t248) * t252 + (-Ifges(4,4) * t249 + (-qJ(2) * t197 + t294) * t251 + (Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + m(5) * (0.1e1 - t416) * t253) * t252 + (t281 + t419) * t248) * t249 + m(6) * (-t242 * t253 * t252 * t249 + t156 * t157 - t158 * t159);
t368 = t1 * qJD(1);
t362 = t247 * t42;
t361 = t247 * t98;
t360 = t249 * mrSges(5,1);
t359 = t249 * mrSges(5,2);
t357 = t250 * t44;
t356 = t250 * t97;
t210 = Ifges(5,2) * t251 + t384;
t187 = t249 * t210;
t189 = t249 * t214;
t246 = t252 ^ 2;
t335 = t159 * t250;
t339 = t157 * t247;
t286 = -t335 + t339;
t336 = t159 * t247;
t338 = t157 * t250;
t343 = qJ(2) * t252;
t82 = t249 * t140;
t83 = t249 * t142;
t84 = t249 * t146;
t4 = -t83 * t411 - t84 * t410 - t322 * t419 + (m(6) * t246 * t251 * t253 + (m(6) * t286 - t197 - t356 + t361) * t343) * t248 + (-t298 * t343 + (-t189 / 0.2e1 + qJ(2) * t360 + t264) * t251 + (t362 / 0.2e1 - t357 / 0.2e1 - t82 / 0.2e1 + t187 / 0.2e1 + (-t252 * t420 - t359) * qJ(2) + (t336 + t338) * mrSges(6,3) - t295) * t248) * t249;
t349 = t4 * qJD(1);
t56 = mrSges(6,1) * t180 - mrSges(6,2) * t280;
t58 = -Ifges(6,5) * t280 - Ifges(6,6) * t180;
t59 = -Ifges(6,2) * t180 - t162;
t60 = -Ifges(6,1) * t280 - t383;
t7 = t157 * t97 - t159 * t98 + (t56 * t343 + t58 * t398) * t248 + (-t42 / 0.2e1 + t60 / 0.2e1 - t159 * mrSges(6,3)) * t180 - (t59 / 0.2e1 + t44 / 0.2e1 - t157 * mrSges(6,3)) * t280;
t348 = t7 * qJD(1);
t347 = -mrSges(5,1) - t291;
t344 = qJ(2) * t246;
t258 = t292 * t392 + (-mrSges(5,3) * t322 + t360) * t390;
t265 = m(6) * (qJ(2) * t326 - t337 - t340);
t205 = -t251 * mrSges(5,1) + t248 * mrSges(5,2);
t266 = t181 * t405 + t182 * t406 + t205 * t398;
t12 = -t252 * mrSges(4,2) - t249 * mrSges(4,1) + (-t265 / 0.2e1 + t359 / 0.2e1 + (t181 * t396 + t182 * t402) * mrSges(6,3) + (t407 + (mrSges(5,3) / 0.2e1 - t282) * t248) * t252) * t248 - t258 + t266 + t389 * t298;
t342 = qJD(1) * t12;
t230 = t242 * t344;
t243 = t249 ^ 2;
t234 = t243 * qJ(2);
t323 = t251 * t197;
t13 = m(3) * qJ(2) + t181 * t98 + t182 * t97 + t249 * t184 + mrSges(3,3) + (t246 + t243) * mrSges(4,3) + (-t346 * t248 + t323) * t252 + m(6) * (t157 * t181 + t159 * t182 + t230) + m(5) * (t245 * t344 + t230 + t234) + m(4) * (t234 + t344);
t341 = t13 * qJD(1);
t334 = t181 * t247;
t333 = t182 * t250;
t332 = t247 * t142;
t330 = t248 * t420;
t325 = t250 * t146;
t318 = qJD(4) * t248;
t317 = qJD(4) * t251;
t206 = Ifges(6,5) * t247 + Ifges(6,6) * t250;
t316 = t206 / 0.2e1 - Ifges(5,6);
t315 = -0.1e1 + t320;
t314 = t44 / 0.4e1 + t59 / 0.4e1;
t313 = t60 / 0.4e1 - t42 / 0.4e1;
t311 = t247 * t392;
t310 = -t329 / 0.2e1;
t309 = t329 / 0.2e1;
t308 = t328 / 0.2e1;
t305 = t324 / 0.2e1;
t304 = t207 * t391;
t303 = t140 / 0.2e1 - t210 / 0.2e1;
t188 = t248 * t212;
t302 = -t142 / 0.4e1 - t188 / 0.4e1;
t186 = t248 * t208;
t301 = -t186 / 0.4e1 + t146 / 0.4e1;
t300 = t403 - t418 / 0.4e1;
t299 = t213 / 0.4e1 + t404;
t296 = (-mrSges(5,1) / 0.2e1 - t291 / 0.2e1) * t248;
t293 = t215 / 0.4e1 - t210 / 0.4e1 + t140 / 0.4e1;
t255 = (t97 * t402 + t98 * t396 + (t180 * t396 - t280 * t400) * mrSges(6,3)) * t248 + t56 * t392;
t283 = t365 / 0.2e1 - t364 / 0.2e1;
t15 = t255 - t283;
t287 = -t15 * qJD(1) + t25 * qJD(3);
t277 = qJ(2) * (-t315 * t242 - t245);
t141 = t251 * t207 + t371;
t143 = Ifges(6,6) * t248 + t251 * t418;
t147 = Ifges(6,5) * t248 + t213 * t251;
t17 = (t325 / 0.2e1 - t332 / 0.2e1 - t141 / 0.2e1 + t214 / 0.2e1 + t417 / 0.2e1) * t251 + (t147 * t394 + t143 * t402 + t215 / 0.2e1 + t303) * t248;
t256 = (-t336 / 0.2e1 - t338 / 0.2e1) * mrSges(6,3) + t148 / 0.4e1 - t187 / 0.4e1 + t82 / 0.4e1 - t362 / 0.4e1 + t357 / 0.4e1;
t257 = -t83 * t401 - t84 * t393 - t189 / 0.4e1 - t144 / 0.4e1 + t40 / 0.4e1 - t366 / 0.2e1 + t367 / 0.2e1;
t263 = t141 / 0.4e1 - t417 / 0.4e1 - t214 / 0.4e1 + t332 / 0.4e1 - t325 / 0.4e1;
t237 = Ifges(5,5) * t251;
t267 = -t280 * t143 / 0.4e1 + t180 * t147 / 0.4e1 - t252 * t237 / 0.4e1;
t3 = ((qJ(2) * t407 - Ifges(5,5) / 0.2e1) * t252 + t293 * t249 + t256) * t251 + ((0.3e1 / 0.4e1 * Ifges(5,6) - t206 / 0.4e1 + ((t284 / 0.2e1 + t353 / 0.2e1) * t250 + (t355 / 0.2e1 + t405) * t247) * qJ(2)) * t252 + (-t345 / 0.2e1 + t263) * t249 + t257) * t248 + t267 + t422;
t275 = t3 * qJD(1) + t22 * qJD(2) + t17 * qJD(3);
t185 = t248 * t206;
t21 = -t185 * t390 + ((t188 / 0.2e1 + t142 / 0.2e1) * t250 + (t146 / 0.2e1 - t186 / 0.2e1) * t247) * t248;
t254 = t302 * t180 - t301 * t280 + t157 * t406 + t159 * t285 / 0.2e1 + t58 * t391;
t261 = -t156 * mrSges(6,1) / 0.2e1 - t158 * mrSges(6,2) / 0.2e1 + t414;
t262 = (t339 / 0.2e1 - t335 / 0.2e1) * mrSges(6,3) - t185 * t397;
t269 = t399 * t345;
t6 = (t313 * t250 - t314 * t247 + (-Ifges(6,3) / 0.2e1 + t269) * t252 + t262) * t248 + t254 + t261;
t274 = -qJD(1) * t6 + qJD(2) * t25 + qJD(3) * t21;
t270 = t299 * t248 + t301;
t271 = t300 * t248 + t302;
t19 = t304 - t371 / 0.2e1 + (-t378 / 0.2e1 + t270) * t250 + (t372 / 0.2e1 + t271) * t247;
t27 = (t212 / 0.2e1 + t418 / 0.2e1) * t250 + (t213 / 0.2e1 - t208 / 0.2e1) * t247;
t259 = t299 * t180 + t314 * t250 + t280 * t300;
t9 = (-t370 / 0.2e1 + (t207 / 0.4e1 + t379 / 0.2e1) * t248) * t249 + (Ifges(6,6) * t310 + t313) * t247 + t259;
t273 = -qJD(1) * t9 - qJD(3) * t19 - qJD(4) * t27;
t272 = (t333 / 0.2e1 - t334 / 0.2e1) * mrSges(6,3);
t10 = (-t57 / 0.2e1 + t199 / 0.2e1) * t248 + t272 + (t296 - m(6) * t277 / 0.2e1) * t252 + (-t352 / 0.2e1 + t286 * t412 - t356 / 0.2e1 + t361 / 0.2e1 - t197 / 0.2e1) * t251;
t102 = t315 * t251 * t248;
t268 = -t10 * qJD(1) + t22 * qJD(3) + t102 * t387;
t77 = t290 * t392 + t276;
t18 = t304 + Ifges(6,5) * t305 + Ifges(6,6) * t311 + t371 / 0.2e1 + t270 * t250 + t271 * t247;
t16 = t255 + t283;
t14 = (-mrSges(6,2) * t327 + t181 * mrSges(6,3)) * t308 - (mrSges(6,1) * t327 - t182 * mrSges(6,3)) * t331 / 0.2e1 + (t330 / 0.2e1 + (t245 / 0.2e1 + t242 / 0.2e1) * mrSges(5,3)) * t252 + t258 + t266 + (-mrSges(5,3) * t327 + t265 - t359) * t399;
t11 = (-t251 * t286 + t252 * t277) * t412 + t97 * t305 + t284 * t249 * t308 + t98 * t311 + t247 * t285 * t310 + t57 * t399 + t420 * t306 - t248 * t199 / 0.2e1 + t323 / 0.2e1 + t272 + (-t354 / 0.2e1 + t296) * t252;
t8 = Ifges(6,3) * t306 + t309 * t373 + t310 * t379 + t313 * t247 + (t290 * t343 / 0.2e1 + t207 * t397) * t248 + t259 + t290 * qJ(2) * t307;
t5 = t423 + (t252 * t269 + t60 * t393 + t42 * t395 + t262 + (t44 + t59) * t401) * t248 + t254 - t261;
t2 = t306 * t421 + t309 * t369 - t310 * t345 + ((Ifges(5,6) / 0.4e1 + t279 * qJ(2)) * t252 + t257) * t248 + (t343 * t420 + t256) * t251 + ((t421 / 0.2e1 + t293) * t251 + (t369 / 0.2e1 + t263) * t248) * t249 + Ifges(5,5) * t322 / 0.2e1 + t267 + (-Ifges(5,6) / 0.2e1 + t206 / 0.4e1) * t327 - t422;
t20 = [qJD(2) * t13 + qJD(3) * t1 + qJD(4) * t4 + qJD(5) * t7, t341 + t14 * qJD(3) + t11 * qJD(4) + t16 * qJD(5) + (-t322 + t333 - t334) * t248 * t387, t14 * qJD(2) + t2 * qJD(4) + t5 * qJD(5) + t368 + (Ifges(4,5) * t252 - Ifges(4,6) * t249 + t142 * t409 + t146 * t408 - t156 * t285 + t158 * t284 + (t214 * t389 - t281) * t251 + (t303 * t252 + t45 * t394 + t43 * t402 + t294) * t248 + ((-mrSges(4,1) + t205) * t252 + (mrSges(4,2) - t298 - t330) * t249) * qJ(2)) * qJD(3), t349 + t11 * qJD(2) + t2 * qJD(3) + (-t83 * t394 - t84 * t400) * qJD(4) + t8 * qJD(5) + (t249 * t316 + t343 * t347) * t317 + ((t208 * t400 + t212 * t396 - Ifges(5,5)) * t249 - t413 * t343) * t318, t348 + t16 * qJD(2) + t5 * qJD(3) + t8 * qJD(4) + (-mrSges(6,1) * t159 - mrSges(6,2) * t157 + t58) * qJD(5); -qJD(3) * t12 - qJD(4) * t10 + qJD(5) * t15 - t341, m(6) * t102 * qJD(4), -t342 + t388, t77 * qJD(5) + t317 * t413 + t347 * t318 + t268, -qJD(5) * t248 * t291 + t77 * qJD(4) - t287; qJD(2) * t12 + qJD(4) * t3 + qJD(5) * t6 - t368, t342 + t388, qJD(4) * t17 - qJD(5) * t21, (t237 + (t212 * t390 + t143 / 0.2e1) * t250 + t316 * t248 + (t208 * t392 + t147 / 0.2e1) * t247) * qJD(4) + t18 * qJD(5) + t275, qJD(4) * t18 - qJD(5) * t185 - t274; qJD(2) * t10 - qJD(3) * t3 + qJD(5) * t9 - t349, -t268, qJD(5) * t19 - t275, t27 * qJD(5), qJD(5) * t207 - t273; -qJD(2) * t15 - qJD(3) * t6 - qJD(4) * t9 - t348, t287, -qJD(4) * t19 + t274, t273, 0;];
Cq = t20;
