% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:31:00
% EndTime: 2019-03-09 03:31:07
% DurationCPUTime: 3.76s
% Computational Cost: add. (5241->460), mult. (9777->589), div. (0->0), fcn. (7491->4), ass. (0->230)
t398 = Ifges(7,4) + Ifges(6,5);
t245 = sin(qJ(5));
t243 = t245 ^ 2;
t247 = cos(qJ(5));
t244 = t247 ^ 2;
t320 = t243 + t244;
t400 = mrSges(7,2) + mrSges(6,3);
t406 = t400 * t320;
t235 = t247 * mrSges(7,3);
t346 = t245 * mrSges(7,1);
t179 = -t235 + t346;
t344 = t247 * mrSges(6,2);
t347 = t245 * mrSges(6,1);
t180 = t344 + t347;
t405 = t179 + t180;
t358 = mrSges(6,1) + mrSges(7,1);
t404 = (t244 / 0.2e1 + t243 / 0.2e1) * t400;
t248 = cos(qJ(3));
t232 = t248 * qJ(4);
t246 = sin(qJ(3));
t321 = -t246 * pkin(3) + t232;
t174 = qJ(2) - t321;
t181 = -t246 * mrSges(5,2) - t248 * mrSges(5,3);
t403 = m(5) * t174 + t181;
t330 = t246 * t247;
t402 = t358 * t330;
t339 = qJ(6) * t245;
t360 = pkin(5) * t247;
t182 = t339 + t360;
t401 = -t182 / 0.2e1;
t370 = -t246 / 0.2e1;
t369 = t246 / 0.2e1;
t366 = t248 / 0.2e1;
t338 = qJ(6) * t247;
t285 = pkin(5) * t245 - t338;
t273 = m(7) * t285;
t399 = -Ifges(4,4) - Ifges(5,6);
t397 = Ifges(7,2) + Ifges(6,3);
t396 = t246 * (m(6) / 0.4e1 + m(7) / 0.4e1);
t332 = t245 * t246;
t160 = t248 * mrSges(6,1) - mrSges(6,3) * t332;
t161 = -t248 * mrSges(7,1) + mrSges(7,2) * t332;
t395 = -t160 + t161;
t162 = -t248 * mrSges(6,2) + mrSges(6,3) * t330;
t343 = t248 * mrSges(7,3);
t163 = mrSges(7,2) * t330 + t343;
t394 = t162 + t163;
t393 = Ifges(7,6) * t332 + t398 * t330;
t392 = -t273 - t405;
t239 = Ifges(7,5) * t247;
t185 = Ifges(7,3) * t245 + t239;
t188 = -Ifges(7,1) * t245 + t239;
t231 = m(7) * qJ(6) + mrSges(7,3);
t142 = pkin(8) * t246 + t174;
t241 = t248 * pkin(4);
t250 = -pkin(1) - pkin(7);
t176 = t248 * t250 - t241;
t59 = -t142 * t245 - t176 * t247;
t328 = t247 * t142;
t333 = t245 * t176;
t60 = t328 - t333;
t281 = t59 * t245 - t60 * t247;
t326 = t248 * qJ(6);
t52 = t60 + t326;
t359 = pkin(5) * t248;
t53 = -t59 - t359;
t283 = t53 * t245 + t52 * t247;
t390 = -m(6) * t281 + m(7) * t283 + t247 * t394 + t403;
t389 = t246 ^ 2;
t388 = t248 ^ 2;
t387 = 2 * qJD(3);
t386 = -m(6) / 0.2e1;
t385 = m(6) / 0.2e1;
t384 = -m(7) / 0.2e1;
t383 = m(7) / 0.2e1;
t382 = mrSges(7,1) / 0.2e1;
t381 = Ifges(7,6) / 0.2e1;
t322 = pkin(3) * t248 + qJ(4) * t246;
t158 = pkin(8) * t248 + t322;
t227 = t246 * t250;
t175 = -pkin(4) * t246 + t227;
t61 = -t158 * t245 + t175 * t247;
t56 = pkin(5) * t246 - t61;
t380 = t56 / 0.2e1;
t379 = t60 / 0.2e1;
t68 = t227 + (-pkin(4) - t182) * t246;
t378 = t68 / 0.2e1;
t377 = qJ(6) / 0.2e1;
t144 = t285 * t246;
t376 = t144 / 0.2e1;
t145 = t179 * t246;
t375 = t145 / 0.2e1;
t146 = t180 * t246;
t374 = t146 / 0.2e1;
t373 = t175 / 0.2e1;
t372 = -t245 / 0.2e1;
t368 = -t247 / 0.2e1;
t367 = t247 / 0.2e1;
t364 = t248 * t273;
t173 = qJ(4) + t285;
t363 = m(7) * t173;
t362 = m(7) * t182;
t361 = m(7) * t245;
t357 = mrSges(6,2) - mrSges(7,3);
t356 = t52 - t60;
t355 = t53 + t59;
t354 = m(7) * qJD(6);
t352 = Ifges(6,4) * t245;
t351 = Ifges(6,4) * t247;
t350 = Ifges(7,5) * t245;
t349 = Ifges(6,6) * t248;
t345 = t246 * mrSges(7,1);
t177 = t247 * mrSges(7,1) + t245 * mrSges(7,3);
t143 = t246 * t177;
t148 = t246 * t185;
t219 = Ifges(6,4) * t330;
t149 = -Ifges(6,2) * t332 + t219;
t216 = Ifges(7,5) * t332;
t150 = Ifges(7,1) * t330 + t216;
t191 = Ifges(6,1) * t247 - t352;
t151 = t246 * t191;
t111 = Ifges(7,4) * t248 - t246 * t188;
t113 = Ifges(6,1) * t332 + Ifges(6,5) * t248 + t219;
t267 = t53 * mrSges(7,2) + t111 / 0.2e1 + t113 / 0.2e1 - t59 * mrSges(6,3);
t107 = Ifges(7,6) * t248 - Ifges(7,3) * t330 + t216;
t287 = Ifges(6,2) * t247 + t352;
t109 = t246 * t287 + t349;
t268 = t107 / 0.2e1 - t109 / 0.2e1 - t52 * mrSges(7,2) - t60 * mrSges(6,3);
t4 = t68 * t145 + t175 * t146 + ((t149 / 0.2e1 - t148 / 0.2e1 + t267) * t247 + (-t349 / 0.2e1 + t150 / 0.2e1 + t151 / 0.2e1 + t268) * t245) * t246 + (m(7) * t68 - t143) * t144 + (m(7) * t53 + t395) * t60 + (m(7) * t52 + t394) * t59 + t393 * t366;
t342 = t4 * qJD(1);
t331 = t245 * t248;
t116 = t160 * t331;
t261 = t246 * t404;
t334 = t245 * t161;
t264 = -t334 / 0.2e1 + t394 * t368;
t251 = (m(7) * t376 + t374 + t375) * t246 + ((-t281 - t283) * t383 + t261 + t264) * t248 + t116 / 0.2e1;
t254 = -t235 / 0.2e1 + t273 / 0.2e1 + t347 / 0.2e1 + t346 / 0.2e1 + t344 / 0.2e1;
t8 = t251 + t254;
t341 = t8 * qJD(1);
t311 = mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t312 = t382 + mrSges(6,1) / 0.2e1;
t258 = (t245 * t312 - t247 * t311) * t248 + t364 / 0.2e1;
t265 = m(7) * (t245 * t355 + t247 * t356);
t296 = -t163 / 0.2e1 - t162 / 0.2e1;
t297 = -t161 / 0.2e1 + t160 / 0.2e1;
t9 = t296 * t247 + t297 * t245 - t265 / 0.2e1 + t261 + t258;
t340 = t9 * qJD(1);
t13 = -t116 + (t334 + t390) * t248;
t337 = qJD(1) * t13;
t14 = t246 * mrSges(4,1) + t248 * mrSges(4,2) + mrSges(3,3) + t395 * t245 + (m(4) + m(3)) * qJ(2) + t390;
t336 = qJD(1) * t14;
t23 = m(7) * (t248 * t52 - t332 * t68) + t143 * t332 + t248 * t163;
t335 = qJD(1) * t23;
t327 = t247 * t248;
t76 = (-t389 / 0.2e1 - t388 / 0.2e1 - 0.1e1 / 0.2e1) * t361;
t325 = t76 * qJD(1);
t62 = t247 * t158 + t245 * t175;
t249 = -pkin(3) - pkin(8);
t324 = t320 * t246 * t249;
t319 = qJD(3) * t247;
t318 = qJD(3) * t248;
t317 = qJD(5) * t245;
t315 = m(7) * t331;
t314 = t245 * t354;
t310 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t309 = t381 - Ifges(6,6) / 0.2e1;
t308 = -t59 / 0.2e1 - t53 / 0.2e1;
t307 = t379 - t52 / 0.2e1;
t304 = -t332 / 0.2e1;
t303 = t332 / 0.2e1;
t286 = -Ifges(7,3) * t247 + t350;
t299 = Ifges(7,6) * t370 + t286 * t366 + Ifges(6,6) * t369 - t287 * t248 / 0.2e1;
t289 = Ifges(6,1) * t245 + t351;
t298 = t398 * t370 + (-t188 + t289) * t366;
t187 = -Ifges(6,2) * t245 + t351;
t295 = -t185 / 0.2e1 + t187 / 0.2e1;
t189 = Ifges(7,1) * t247 + t350;
t294 = t189 / 0.2e1 + t191 / 0.2e1;
t292 = m(5) * t250 - mrSges(5,1);
t291 = mrSges(6,2) * t304 + mrSges(7,3) * t303 + t402 / 0.2e1;
t290 = 0.2e1 * t320 * t396;
t178 = t247 * mrSges(6,1) - t245 * mrSges(6,2);
t54 = -qJ(6) * t246 + t62;
t69 = -t241 + (-t182 + t250) * t248;
t3 = -t69 * t143 + t61 * t160 + t56 * t161 + t62 * t162 + t54 * t163 + t403 * t322 + m(7) * (t52 * t54 + t53 * t56 + t68 * t69) + m(6) * (t175 * t176 + t59 * t61 + t60 * t62) + (-t59 * mrSges(6,1) + t53 * mrSges(7,1) - qJ(2) * mrSges(4,2) + t60 * mrSges(6,2) + t174 * mrSges(5,3) - t52 * mrSges(7,3) + (-t176 * mrSges(6,1) - t299) * t247 + (t176 * mrSges(6,2) + t298) * t245 + (-t245 * t310 + t247 * t309 - t399) * t246) * t246 + (qJ(2) * mrSges(4,1) - t174 * mrSges(5,2) + t399 * t248 + (-t175 * mrSges(6,1) - t68 * mrSges(7,1) - t248 * t309 - t268) * t247 + (t175 * mrSges(6,2) - t68 * mrSges(7,3) + t248 * t310 + t267) * t245 + (-Ifges(4,1) - Ifges(5,2) + Ifges(5,3) + Ifges(4,2) - t397) * t246) * t248;
t280 = t245 * t62 + t247 * t61;
t282 = t245 * t54 - t247 * t56;
t6 = (t143 / 0.2e1 + (-t282 + t68) * t384 + (t175 - t280) * t386) * t248 + (t161 * t367 + t178 * t366 + t160 * t368 + (t245 * t52 - t247 * t53 + t69) * t384 + (t245 * t60 + t247 * t59 + t176) * t386 + t394 * t372) * t246;
t284 = t3 * qJD(1) - t6 * qJD(2);
t37 = 0.4e1 * (0.1e1 - t320) * t248 * t396;
t279 = -t6 * qJD(1) + t37 * qJD(2);
t257 = (-t247 * t68 + (-t173 * t246 + t248 * t249) * t245) * t383 - t143 * t368;
t276 = m(7) * t380 + t345 / 0.2e1;
t21 = (mrSges(7,2) * t248 + t179 * t369) * t245 - t257 + t276;
t64 = (t179 + t363) * t247;
t278 = qJD(1) * t21 + qJD(3) * t64;
t29 = t343 + 0.2e1 * (t326 / 0.2e1 + t328 / 0.4e1 - t333 / 0.4e1 - t60 / 0.4e1) * m(7);
t277 = qJD(1) * t29 + qJD(5) * t231;
t275 = mrSges(7,2) * t331 + t345;
t274 = mrSges(7,2) * t327 - t246 * mrSges(7,3);
t272 = -t289 / 0.4e1 + t188 / 0.4e1 - t187 / 0.4e1 + t185 / 0.4e1;
t271 = t191 / 0.4e1 + t189 / 0.4e1 - t287 / 0.4e1 + t286 / 0.4e1;
t237 = Ifges(7,6) * t247;
t252 = (t144 * t173 + t182 * t68) * t383 + qJ(4) * t374 + t179 * t376 + t173 * t375 + t178 * t373 + t143 * t401 + t248 * t237 / 0.4e1 + t177 * t378;
t253 = (-pkin(5) * t56 + qJ(6) * t54) * t384 - t54 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t380 - t61 * mrSges(6,1) / 0.2e1 + t62 * mrSges(6,2) / 0.2e1;
t255 = (t356 * t383 - t296) * t249 + t107 / 0.4e1 - t109 / 0.4e1 + t150 / 0.4e1 + t151 / 0.4e1;
t256 = (t355 * t383 - t297) * t249 - t111 / 0.4e1 - t113 / 0.4e1 + t148 / 0.4e1 - t149 / 0.4e1;
t262 = t249 * t404;
t1 = t252 + (pkin(5) * t382 + mrSges(7,3) * t377 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - t262) * t246 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t381) * t248 + (-t326 / 0.2e1 + t307) * mrSges(7,2) + t271 * t246 + t255) * t247 + ((-0.3e1 / 0.4e1 * Ifges(7,4) - 0.3e1 / 0.4e1 * Ifges(6,5)) * t248 + (t359 / 0.2e1 + t308) * mrSges(7,2) + t272 * t246 + t256) * t245 + t253;
t19 = qJ(4) * t178 + t182 * t179 + (t177 + t362) * t173 + (t188 / 0.2e1 - t289 / 0.2e1 - t295) * t247 + (t287 / 0.2e1 - t286 / 0.2e1 - t294) * t245;
t25 = (-t178 / 0.2e1 - t177 / 0.2e1 + (t360 / 0.2e1 + t339 / 0.2e1 + t401) * m(7)) * t246 + t291;
t269 = t1 * qJD(1) - t25 * qJD(2) + t19 * qJD(3);
t259 = t280 * t386 + t282 * t384;
t263 = m(6) * t373 + m(7) * t378 + mrSges(6,2) * t303 + mrSges(7,3) * t304 - t402 / 0.2e1;
t12 = (t245 * t311 + t247 * t312) * t246 + t259 + t263;
t39 = (t386 + t384) * t246 + t290;
t44 = t363 - t235 + t344 + mrSges(5,3) + t358 * t245 + (m(6) + m(5)) * qJ(4);
t266 = qJD(1) * t12 - qJD(2) * t39 + qJD(3) * t44;
t164 = (m(7) * t249 - mrSges(7,2)) * t245;
t75 = (0.1e1 - t388 - t389) * t361 / 0.2e1;
t38 = m(5) * t246 + t290 + (m(6) + m(7)) * t369;
t27 = (t60 + 0.2e1 * t326) * t383 + m(7) * t379 + t163;
t26 = t246 * t362 + t291 + (t177 + t178) * t369;
t22 = t179 * t304 + t257 + t276;
t11 = (-t246 * mrSges(6,1) - mrSges(6,3) * t331) * t367 + t275 * t368 + t292 * t246 - t259 + t263 + (t246 * mrSges(6,2) + mrSges(6,3) * t327 + t274) * t245 / 0.2e1;
t10 = t265 / 0.2e1 + t160 * t372 + t258 - t264 + t370 * t406;
t7 = t251 - t254;
t5 = t6 * qJD(3);
t2 = (t272 * t245 + t271 * t247 - t262) * t246 + t252 - pkin(5) * t275 / 0.2e1 + t274 * t377 + ((-Ifges(7,4) / 0.4e1 - Ifges(6,5) / 0.4e1) * t248 + t308 * mrSges(7,2) + t256) * t245 + (-t349 / 0.4e1 + t307 * mrSges(7,2) + t255) * t247 - t253 + t397 * t370 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t327 + t398 * t331 / 0.2e1;
t15 = [qJD(2) * t14 + qJD(3) * t3 - qJD(4) * t13 + qJD(5) * t4 + qJD(6) * t23, qJD(5) * t7 + qJD(6) * t75 + t336 - t5, t11 * qJD(4) + t2 * qJD(5) + t22 * qJD(6) + (t56 * mrSges(7,2) - t61 * mrSges(6,3) + t298) * t319 + ((t173 * t69 + t249 * t282) * t383 + (qJ(4) * t176 + t249 * t280) * t385) * t387 + (Ifges(5,5) - Ifges(4,6) + (-mrSges(4,2) + mrSges(5,3)) * t250 + (-t173 * mrSges(7,1) + t295) * t247 + (-t173 * mrSges(7,3) + t294) * t245 + (-t178 + t292) * qJ(4)) * t318 + t284 + (t176 * t180 + t69 * t179 + (-t54 * mrSges(7,2) - t62 * mrSges(6,3) + t299) * t245 + (pkin(3) * mrSges(5,1) + Ifges(5,4) - Ifges(4,5) + (-m(5) * pkin(3) - mrSges(4,1) + mrSges(5,2)) * t250 + (-t249 * t358 - t310) * t247 + (t249 * t357 - t309) * t245) * t246) * qJD(3), qJD(3) * t11 + qJD(5) * t10 - t337, t7 * qJD(2) + t2 * qJD(3) + t10 * qJD(4) + t27 * qJD(6) + t342 + ((-mrSges(7,2) * t182 - Ifges(6,6) * t245) * t246 + (-m(7) * pkin(5) - t358) * t60 + (-mrSges(6,2) + t231) * t59 + t393) * qJD(5), qJD(2) * t75 + qJD(3) * t22 + qJD(5) * t27 + t335; qJD(5) * t8 + qJD(6) * t76 - t336 - t5, t37 * qJD(3), -t181 * qJD(3) + t38 * qJD(4) + t26 * qJD(5) + (-mrSges(4,2) + t405) * t318 + ((t173 * t248 + t324) * t383 + (t232 + t324) * t385 + m(5) * t321 / 0.2e1) * t387 + ((-mrSges(4,1) - t406) * qJD(3) - t247 * t354) * t246 + t279, t38 * qJD(3), t341 + t26 * qJD(3) + qJD(5) * t364 + (t357 * qJD(5) * t247 + (qJD(5) * t358 - t354) * t245) * t248, t325 + (-t246 * t319 - t248 * t317) * m(7); qJD(4) * t12 + qJD(5) * t1 - qJD(6) * t21 - t284, -qJD(4) * t39 - qJD(5) * t25 - t279, qJD(4) * t44 + qJD(5) * t19 - qJD(6) * t64, t266, t164 * qJD(6) + (pkin(5) * mrSges(7,2) - t398) * t317 + t269 + (-mrSges(7,2) * t338 - Ifges(6,6) * t247 + t249 * t392 + t237) * qJD(5), qJD(5) * t164 - t278; -qJD(3) * t12 - qJD(5) * t9 + t248 * t314 + t337, t39 * qJD(3), -t266, 0, qJD(5) * t392 + t314 - t340 (qJD(1) * t248 + qJD(5)) * t361; -qJD(2) * t8 - qJD(3) * t1 + qJD(4) * t9 + qJD(6) * t29 - t342, qJD(3) * t25 - t341, -t269, t340, t231 * qJD(6), t277; -qJD(2) * t76 + qJD(3) * t21 - qJD(4) * t315 - qJD(5) * t29 - t335, -t325, t278, -qJD(1) * t315, -t277, 0;];
Cq  = t15;
