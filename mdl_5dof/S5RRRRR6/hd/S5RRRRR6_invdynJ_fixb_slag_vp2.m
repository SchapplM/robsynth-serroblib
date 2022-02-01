% Calculate vector of inverse dynamics joint torques for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:40
% EndTime: 2022-01-20 12:07:55
% DurationCPUTime: 5.28s
% Computational Cost: add. (8737->478), mult. (13204->640), div. (0->0), fcn. (8699->16), ass. (0->245)
t280 = sin(qJ(4));
t281 = sin(qJ(3));
t285 = cos(qJ(4));
t286 = cos(qJ(3));
t207 = -t280 * t281 + t285 * t286;
t289 = -pkin(8) - pkin(7);
t314 = qJD(3) * t289;
t215 = t281 * t314;
t216 = t286 * t314;
t228 = t289 * t281;
t268 = t286 * pkin(8);
t229 = pkin(7) * t286 + t268;
t287 = cos(qJ(2));
t348 = pkin(1) * qJD(1);
t317 = t287 * t348;
t321 = qJD(4) * t285;
t322 = qJD(4) * t280;
t380 = -t207 * t317 + t215 * t285 + t216 * t280 + t228 * t321 - t229 * t322;
t160 = t228 * t280 + t229 * t285;
t208 = t280 * t286 + t281 * t285;
t379 = -qJD(4) * t160 + t208 * t317 - t215 * t280 + t216 * t285;
t277 = qJ(3) + qJ(4);
t267 = qJ(5) + t277;
t250 = sin(t267);
t251 = cos(t267);
t398 = t251 * mrSges(6,1) - t250 * mrSges(6,2);
t263 = sin(t277);
t265 = cos(t277);
t297 = -mrSges(5,1) * t265 + mrSges(5,2) * t263 - t398;
t227 = -t286 * mrSges(4,1) + t281 * mrSges(4,2);
t390 = -mrSges(3,1) + t227;
t397 = t390 + t297;
t298 = t207 * qJD(4);
t141 = qJD(3) * t207 + t298;
t362 = pkin(9) * t141;
t396 = -t362 + t379;
t299 = t208 * qJD(4);
t142 = -qJD(3) * t208 - t299;
t139 = t142 * pkin(9);
t395 = -t139 - t380;
t367 = t281 / 0.2e1;
t273 = qJD(1) + qJD(2);
t180 = t207 * t273;
t181 = t208 * t273;
t279 = sin(qJ(5));
t284 = cos(qJ(5));
t310 = t180 * t284 - t181 * t279;
t113 = Ifges(6,4) * t310;
t117 = t180 * t279 + t181 * t284;
t272 = qJD(3) + qJD(4);
t262 = qJD(5) + t272;
t54 = Ifges(6,1) * t117 + Ifges(6,5) * t262 + t113;
t394 = t113 + t54;
t253 = pkin(3) * t285 + pkin(4);
t319 = qJD(5) * t284;
t320 = qJD(5) * t279;
t328 = t280 * t284;
t282 = sin(qJ(2));
t318 = t282 * t348;
t220 = pkin(7) * t273 + t318;
t312 = pkin(8) * t273 + t220;
t166 = t312 * t286;
t155 = t285 * t166;
t165 = t312 * t281;
t104 = t165 * t280 - t155;
t361 = pkin(9) * t180;
t79 = t104 - t361;
t153 = t280 * t166;
t105 = -t165 * t285 - t153;
t174 = t181 * pkin(9);
t80 = -t174 + t105;
t393 = t279 * t80 - t284 * t79 - t253 * t320 + (-t280 * t319 + (-t279 * t285 - t328) * qJD(4)) * pkin(3);
t329 = t279 * t280;
t392 = -t279 * t79 - t284 * t80 + t253 * t319 + (-t280 * t320 + (t284 * t285 - t329) * qJD(4)) * pkin(3);
t331 = t273 * t281;
t218 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t331;
t330 = t273 * t286;
t219 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t330;
t391 = (mrSges(3,2) * t273 + t218 * t281 - t219 * t286) * t287;
t278 = qJ(1) + qJ(2);
t264 = sin(t278);
t266 = cos(t278);
t377 = g(1) * t266 + g(2) * t264;
t389 = mrSges(3,2) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3);
t221 = -pkin(2) * t273 - t317;
t355 = mrSges(4,2) * t286;
t388 = t221 * (mrSges(4,1) * t281 + t355) + qJD(3) * (Ifges(4,5) * t286 - Ifges(4,6) * t281) / 0.2e1;
t349 = Ifges(6,4) * t117;
t53 = Ifges(6,2) * t310 + Ifges(6,6) * t262 + t349;
t387 = t53 / 0.2e1;
t373 = -t310 / 0.2e1;
t386 = mrSges(6,2) * t310;
t385 = Ifges(6,1) * t310;
t384 = Ifges(6,5) * t310;
t346 = qJD(3) * pkin(3);
t156 = -t165 + t346;
t101 = t156 * t280 + t155;
t77 = t101 + t361;
t343 = t279 * t77;
t100 = t156 * t285 - t153;
t76 = t100 - t174;
t72 = pkin(4) * t272 + t76;
t35 = t284 * t72 - t343;
t383 = t310 * t35;
t159 = t228 * t285 - t229 * t280;
t360 = pkin(9) * t208;
t125 = t159 - t360;
t200 = t207 * pkin(9);
t126 = t200 + t160;
t67 = t125 * t284 - t126 * t279;
t382 = t67 * qJD(5) + t279 * t396 - t284 * t395;
t68 = t125 * t279 + t126 * t284;
t381 = -t68 * qJD(5) + t279 * t395 + t284 * t396;
t252 = pkin(1) * t282 + pkin(7);
t356 = -pkin(8) - t252;
t201 = t356 * t281;
t332 = t252 * t286;
t202 = t268 + t332;
t135 = t201 * t280 + t202 * t285;
t376 = t390 * t273;
t374 = t221 * t282 + (t281 ^ 2 + t286 ^ 2) * t220 * t287;
t372 = -t117 / 0.2e1;
t371 = t117 / 0.2e1;
t369 = t181 / 0.2e1;
t368 = -t262 / 0.2e1;
t341 = t284 * t77;
t36 = t279 * t72 + t341;
t365 = mrSges(6,3) * t36;
t364 = pkin(1) * t287;
t363 = pkin(3) * t281;
t357 = g(3) * t286;
t353 = mrSges(5,3) * t180;
t352 = Ifges(4,4) * t281;
t351 = Ifges(4,4) * t286;
t350 = Ifges(5,4) * t181;
t347 = pkin(1) * qJD(2);
t345 = t181 * mrSges(5,3);
t339 = t286 * Ifges(4,2);
t338 = pkin(1) * qJDD(1);
t315 = qJD(1) * t347;
t205 = t282 * t338 + t287 * t315;
t271 = qJDD(1) + qJDD(2);
t184 = pkin(7) * t271 + t205;
t324 = qJD(3) * t281;
t127 = t184 * t286 - t220 * t324;
t337 = t127 * t286;
t323 = qJD(3) * t286;
t128 = -t184 * t281 - t220 * t323;
t336 = t128 * t281;
t325 = pkin(2) * t266 + pkin(7) * t264;
t270 = qJDD(3) + qJDD(4);
t316 = t287 * t347;
t257 = pkin(3) * t324;
t254 = pkin(3) * t286 + pkin(2);
t129 = -pkin(4) * t142 + t257;
t311 = qJD(3) * t356;
t134 = t201 * t285 - t202 * t280;
t249 = pkin(4) * t265;
t214 = t249 + t254;
t274 = pkin(9) - t289;
t309 = t214 * t266 + t264 * t274;
t308 = t254 * t266 - t264 * t289;
t204 = -t282 * t315 + t287 * t338;
t307 = mrSges(6,1) * t250 + mrSges(6,2) * t251;
t305 = t339 + t352;
t108 = t134 - t360;
t109 = t200 + t135;
t50 = t108 * t284 - t109 * t279;
t51 = t108 * t279 + t109 * t284;
t137 = t207 * t284 - t208 * t279;
t138 = t207 * t279 + t208 * t284;
t303 = t218 * t286 + t219 * t281;
t177 = -pkin(4) * t207 - t254;
t261 = qJDD(5) + t270;
t194 = t271 * t281 + t273 * t323;
t96 = qJDD(3) * pkin(3) - pkin(8) * t194 + t128;
t193 = t271 * t286 - t273 * t324;
t99 = pkin(8) * t193 + t127;
t29 = -qJD(4) * t101 - t280 * t99 + t285 * t96;
t93 = t193 * t280 + t194 * t285 + t273 * t298;
t15 = pkin(4) * t270 - pkin(9) * t93 + t29;
t28 = t156 * t321 - t166 * t322 + t280 * t96 + t285 * t99;
t94 = t193 * t285 - t194 * t280 - t273 * t299;
t19 = pkin(9) * t94 + t28;
t3 = qJD(5) * t35 + t15 * t279 + t19 * t284;
t32 = qJD(5) * t310 + t279 * t94 + t284 * t93;
t33 = -qJD(5) * t117 - t279 * t93 + t284 * t94;
t4 = -qJD(5) * t36 + t15 * t284 - t19 * t279;
t302 = mrSges(6,1) * t4 - t3 * mrSges(6,2) + Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * t261;
t300 = t281 * (Ifges(4,1) * t286 - t352);
t161 = t281 * t311 + t286 * t316;
t162 = -t281 * t316 + t286 * t311;
t60 = t161 * t285 + t162 * t280 + t201 * t321 - t202 * t322;
t183 = -pkin(2) * t271 - t204;
t182 = -t254 * t273 - t317;
t130 = -pkin(3) * t193 + t183;
t296 = mrSges(5,1) * t263 + mrSges(5,2) * t265 + t307;
t61 = -qJD(4) * t135 - t161 * t280 + t162 * t285;
t295 = -t303 * qJD(3) + t286 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t193) - t281 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t194);
t294 = t389 * t264 + t266 * t397;
t111 = Ifges(5,2) * t180 + Ifges(5,6) * t272 + t350;
t173 = Ifges(5,4) * t180;
t112 = Ifges(5,1) * t181 + Ifges(5,5) * t272 + t173;
t124 = -pkin(4) * t180 + t182;
t292 = mrSges(6,3) * t383 + t29 * mrSges(5,1) - t28 * mrSges(5,2) - t182 * (mrSges(5,1) * t181 + mrSges(5,2) * t180) - t124 * t386 - t272 * (Ifges(5,5) * t180 - Ifges(5,6) * t181) / 0.2e1 + Ifges(5,3) * t270 + t302 + t385 * t372 + t384 * t368 + t100 * t353 + t111 * t369 - t181 * (Ifges(5,1) * t180 - t350) / 0.2e1 + Ifges(5,6) * t94 + Ifges(5,5) * t93 + t394 * t373 - (-Ifges(5,2) * t181 + t112 + t173) * t180 / 0.2e1 + (-t124 * mrSges(6,1) - Ifges(6,4) * t372 - Ifges(6,2) * t373 - Ifges(6,6) * t368 + t365 + t387) * t117;
t291 = (m(4) * pkin(2) + m(5) * t254 + m(6) * t214 - t397) * t264 + (-m(4) * pkin(7) + m(5) * t289 - m(6) * t274 + t389) * t266;
t178 = Ifges(4,6) * qJD(3) + t273 * t305;
t233 = Ifges(4,4) * t330;
t179 = Ifges(4,1) * t331 + Ifges(4,5) * qJD(3) + t233;
t55 = qJD(5) * t137 + t141 * t284 + t142 * t279;
t56 = -qJD(5) * t138 - t141 * t279 + t142 * t284;
t57 = -pkin(4) * t94 + t130;
t290 = t310 * (Ifges(6,4) * t55 + Ifges(6,2) * t56) / 0.2e1 + (mrSges(5,2) * t130 - mrSges(5,3) * t29 + Ifges(5,1) * t93 + Ifges(5,4) * t94 + Ifges(5,5) * t270) * t208 + (-mrSges(5,1) * t130 + mrSges(5,3) * t28 + Ifges(5,4) * t93 + Ifges(5,2) * t94 + Ifges(5,6) * t270) * t207 + (mrSges(6,2) * t57 - mrSges(6,3) * t4 + Ifges(6,1) * t32 + Ifges(6,4) * t33 + Ifges(6,5) * t261) * t138 + (-mrSges(6,1) * t57 + mrSges(6,3) * t3 + Ifges(6,4) * t32 + Ifges(6,2) * t33 + Ifges(6,6) * t261) * t137 + t286 * (Ifges(4,4) * t194 + Ifges(4,2) * t193) / 0.2e1 + (t179 + t273 * (-Ifges(4,2) * t281 + t351)) * t323 / 0.2e1 + (Ifges(4,1) * t194 + Ifges(4,4) * t193) * t367 + (-t100 * t141 + t101 * t142) * mrSges(5,3) + (0.2e1 * Ifges(4,5) * t367 + Ifges(4,6) * t286) * qJDD(3) + (t300 * t273 / 0.2e1 + t388) * qJD(3) + t56 * t387 + t272 * (Ifges(5,5) * t141 + Ifges(5,6) * t142) / 0.2e1 + Ifges(3,3) * t271 + t262 * (Ifges(6,5) * t55 + Ifges(6,6) * t56) / 0.2e1 + t183 * t227 + t204 * mrSges(3,1) - t205 * mrSges(3,2) + t182 * (-mrSges(5,1) * t142 + mrSges(5,2) * t141) + t180 * (Ifges(5,4) * t141 + Ifges(5,2) * t142) / 0.2e1 + (Ifges(5,1) * t141 + Ifges(5,4) * t142) * t369 + (Ifges(6,1) * t55 + Ifges(6,4) * t56) * t371 + t56 * t365 + mrSges(4,3) * t337 - t35 * t55 * mrSges(6,3) + t193 * t305 / 0.2e1 - t178 * t324 / 0.2e1 + t194 * (t281 * Ifges(4,1) + t351) / 0.2e1 + t55 * t54 / 0.2e1 + t124 * (-mrSges(6,1) * t56 + mrSges(6,2) * t55) + t141 * t112 / 0.2e1 + t142 * t111 / 0.2e1;
t288 = cos(qJ(1));
t283 = sin(qJ(1));
t269 = t288 * pkin(1);
t258 = t282 * t347;
t255 = -pkin(2) - t364;
t226 = -t254 - t364;
t217 = t258 + t257;
t191 = pkin(3) * t328 + t253 * t279;
t190 = -pkin(3) * t329 + t253 * t284;
t164 = t177 - t364;
t152 = mrSges(5,1) * t272 - t345;
t151 = -mrSges(5,2) * t272 + t353;
t143 = pkin(3) * t331 + pkin(4) * t181;
t131 = -mrSges(4,1) * t193 + mrSges(4,2) * t194;
t122 = -mrSges(5,1) * t180 + mrSges(5,2) * t181;
t120 = t129 + t258;
t103 = mrSges(6,1) * t262 - mrSges(6,3) * t117;
t102 = -mrSges(6,2) * t262 + mrSges(6,3) * t310;
t85 = -mrSges(5,2) * t270 + mrSges(5,3) * t94;
t84 = mrSges(5,1) * t270 - mrSges(5,3) * t93;
t59 = -mrSges(6,1) * t310 + mrSges(6,2) * t117;
t49 = t61 - t362;
t48 = t139 + t60;
t47 = -mrSges(5,1) * t94 + mrSges(5,2) * t93;
t38 = t284 * t76 - t343;
t37 = -t279 * t76 - t341;
t26 = -mrSges(6,2) * t261 + mrSges(6,3) * t33;
t25 = mrSges(6,1) * t261 - mrSges(6,3) * t32;
t9 = -mrSges(6,1) * t33 + mrSges(6,2) * t32;
t8 = -qJD(5) * t51 - t279 * t48 + t284 * t49;
t7 = qJD(5) * t50 + t279 * t49 + t284 * t48;
t1 = [t290 + m(6) * (t120 * t124 + t164 * t57 + t3 * t51 + t35 * t8 + t36 * t7 + t4 * t50) + m(5) * (t100 * t61 + t101 * t60 + t130 * t226 + t134 * t29 + t135 * t28 + t182 * t217) + m(4) * (t127 * t332 + t183 * t255 - t252 * t336) + ((mrSges(3,1) * t287 - mrSges(3,2) * t282) * t271 + (-g(2) * t288 + t204 * t287 + t205 * t282) * m(3) + (m(3) + m(4) + m(5) + m(6)) * t283 * g(1) + (m(4) * t374 + t282 * t376 - t391) * qJD(2)) * pkin(1) + Ifges(2,3) * qJDD(1) + t255 * t131 + t217 * t122 + t226 * t47 - mrSges(4,3) * t336 + t295 * t252 + (-mrSges(2,1) * t288 + mrSges(2,2) * t283 - m(5) * (t269 + t308) - m(4) * (t269 + t325) - m(6) * (t269 + t309) + t294) * g(2) + (mrSges(2,1) * t283 + mrSges(2,2) * t288 + t291) * g(1) + t50 * t25 + t51 * t26 + t7 * t102 + t8 * t103 + t120 * t59 + t134 * t84 + t135 * t85 + t60 * t151 + t61 * t152 + t164 * t9; t290 + t379 * t152 + t380 * t151 + t381 * t103 + (-mrSges(4,3) * t128 + t122 * t346) * t281 + t382 * t102 - t254 * t47 + t177 * t9 + (t391 + (-t122 - t376 - t59) * t282) * t348 + t294 * g(2) + t291 * g(1) + t295 * pkin(7) + t67 * t25 + t68 * t26 + t129 * t59 - pkin(2) * t131 + t159 * t84 + t160 * t85 + (-t309 * g(2) + t177 * t57 + t3 * t68 + t4 * t67 + t382 * t36 + t381 * t35 + (t129 - t318) * t124) * m(6) + (-t308 * g(2) - t130 * t254 + t159 * t29 + t160 * t28 + (t257 - t318) * t182 + t380 * t101 + t379 * t100) * m(5) + (-pkin(2) * t183 + (-t336 + t337) * pkin(7) - t374 * t348 - t325 * g(2)) * m(4); -m(5) * (t100 * t104 + t101 * t105) + (t227 + t297) * g(3) + t303 * t220 + t393 * t103 + t392 * t102 + (t178 * t367 + (-t300 / 0.2e1 + t339 * t367) * t273 + (-m(5) * t182 - t122) * t363 - (t179 + t233) * t286 / 0.2e1 - t388) * t273 + (t280 * t85 + t285 * t84 + (t151 * t285 - t152 * t280) * qJD(4) + (-t100 * t322 + t101 * t321 + t28 * t280 + t285 * t29 - t357) * m(5)) * pkin(3) + Ifges(4,5) * t194 + t190 * t25 + t191 * t26 + Ifges(4,6) * t193 + t101 * t345 + t292 + Ifges(4,3) * qJDD(3) - t127 * mrSges(4,2) + t128 * mrSges(4,1) - t143 * t59 - t105 * t151 - t104 * t152 + t377 * (-m(6) * (-pkin(4) * t263 - t363) + t355 + (m(5) * pkin(3) + mrSges(4,1)) * t281 + t296) + (-t357 * pkin(3) - t249 * g(3) - t124 * t143 + t190 * t4 + t191 * t3 + t393 * t35 + t392 * t36) * m(6); t297 * g(3) + (t152 + t345) * t101 - m(6) * (t35 * t37 + t36 * t38) + (-t181 * t59 + t284 * t25 + t279 * t26 + (t102 * t284 - t103 * t279) * qJD(5) + (-g(3) * t265 - t124 * t181 + t263 * t377 + t279 * t3 + t284 * t4 + t36 * t319 - t35 * t320) * m(6)) * pkin(4) + t292 - t38 * t102 - t37 * t103 - t100 * t151 + t377 * t296; -t124 * (mrSges(6,1) * t117 + t386) + (-t349 + t385) * t372 + t53 * t371 + (-Ifges(6,6) * t117 + t384) * t368 - t35 * t102 + t36 * t103 - g(3) * t398 + (t117 * t36 + t383) * mrSges(6,3) + t302 + (-Ifges(6,2) * t117 + t394) * t373 + t377 * t307;];
tau = t1;
