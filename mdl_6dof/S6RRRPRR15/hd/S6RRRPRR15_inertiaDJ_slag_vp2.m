% Calculate time derivative of joint inertia matrix for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR15_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:39
% EndTime: 2019-03-09 20:25:03
% DurationCPUTime: 10.35s
% Computational Cost: add. (15998->883), mult. (47587->1250), div. (0->0), fcn. (47241->12), ass. (0->362)
t289 = sin(pkin(6));
t439 = 0.2e1 * t289;
t292 = sin(qJ(6));
t296 = cos(qJ(6));
t297 = cos(qJ(5));
t355 = qJD(6) * t297;
t293 = sin(qJ(5));
t362 = qJD(5) * t293;
t308 = t292 * t362 - t296 * t355;
t288 = sin(pkin(7));
t294 = sin(qJ(3));
t379 = t288 * t294;
t275 = pkin(10) * t379;
t290 = cos(pkin(7));
t298 = cos(qJ(3));
t347 = -pkin(2) * t298 - pkin(3);
t166 = pkin(4) * t379 + t275 + (-pkin(11) + t347) * t290;
t336 = -qJ(4) * t294 - pkin(2);
t426 = pkin(3) + pkin(11);
t185 = (-t298 * t426 + t336) * t288;
t437 = t293 * t166 + t297 * t185;
t438 = qJD(5) * t437;
t405 = t292 / 0.2e1;
t403 = t296 / 0.2e1;
t370 = t292 ^ 2 + t296 ^ 2;
t291 = cos(pkin(6));
t295 = sin(qJ(2));
t371 = t295 * t298;
t299 = cos(qJ(2));
t372 = t294 * t299;
t313 = t290 * t372 + t371;
t168 = t289 * t313 + t291 * t379;
t376 = t289 * t299;
t224 = -t288 * t376 + t290 * t291;
t401 = pkin(1) * t291;
t280 = t295 * t401;
t231 = pkin(9) * t376 + t280;
t351 = t290 * t376;
t314 = t288 * t291 + t351;
t162 = pkin(10) * t314 + t231;
t281 = t299 * t401;
t331 = t289 * (-pkin(10) * t290 - pkin(9));
t317 = t295 * t331;
t175 = pkin(2) * t291 + t281 + t317;
t399 = pkin(10) * t288;
t198 = (-pkin(2) * t299 - t295 * t399 - pkin(1)) * t289;
t87 = -t294 * t162 + t298 * (t175 * t290 + t198 * t288);
t52 = pkin(4) * t168 - t224 * t426 - t87;
t377 = t289 * t295;
t350 = t294 * t377;
t378 = t288 * t298;
t167 = -t291 * t378 - t298 * t351 + t350;
t121 = -t175 * t288 + t290 * t198;
t312 = -qJ(4) * t168 + t121;
t54 = t167 * t426 + t312;
t397 = t293 * t52 + t297 * t54;
t436 = qJD(5) * t397;
t367 = qJD(3) * t294;
t344 = t288 * t367;
t271 = pkin(3) * t344;
t365 = qJD(4) * t294;
t159 = t271 + (-t365 + (pkin(11) * t294 - qJ(4) * t298) * qJD(3)) * t288;
t375 = t290 * t294;
t279 = pkin(2) * t375;
t425 = pkin(4) + pkin(10);
t187 = (t378 * t425 + t279) * qJD(3);
t67 = -t159 * t293 + t187 * t297 - t438;
t369 = qJD(2) * t289;
t345 = t299 * t369;
t346 = t295 * t369;
t128 = -t346 * t375 - qJD(3) * t350 + (qJD(3) * t314 + t345) * t298;
t178 = (t299 * t331 - t280) * qJD(2);
t153 = t175 * t375;
t273 = qJD(2) * t281;
t177 = qJD(2) * t317 + t273;
t366 = qJD(3) * t298;
t333 = -qJD(3) * t153 - t162 * t366 - t294 * t177 - t198 * t344;
t374 = t290 * t298;
t307 = -t178 * t374 - t333;
t199 = (pkin(2) * t295 - t299 * t399) * t369;
t380 = t199 * t298;
t31 = pkin(4) * t128 + (-t346 * t426 - t380) * t288 + t307;
t127 = t291 * t344 + (t313 * qJD(3) + (t290 * t371 + t372) * qJD(2)) * t289;
t129 = -t178 * t288 + t290 * t199;
t305 = -qJ(4) * t128 - qJD(4) * t168 + t129;
t33 = t127 * t426 + t305;
t6 = -t293 * t33 + t297 * t31 - t436;
t435 = 2 * m(4);
t434 = 2 * m(5);
t433 = 2 * m(6);
t432 = 2 * m(7);
t431 = 2 * mrSges(5,1);
t430 = -2 * mrSges(3,3);
t429 = -2 * mrSges(4,3);
t4 = -pkin(5) * t128 - t6;
t428 = m(7) * t4;
t131 = t167 * t293 + t224 * t297;
t332 = t288 * t346;
t77 = qJD(5) * t131 - t127 * t297 + t293 * t332;
t319 = t167 * t297 - t224 * t293;
t78 = qJD(5) * t319 + t127 * t293 + t297 * t332;
t24 = Ifges(6,1) * t78 - Ifges(6,4) * t77 + Ifges(6,5) * t128;
t427 = t24 / 0.2e1;
t343 = t288 * t366;
t63 = -pkin(5) * t343 - t67;
t424 = m(7) * t63;
t225 = t290 * t293 + t297 * t378;
t179 = -qJD(5) * t225 + t293 * t344;
t352 = t293 * t378;
t360 = qJD(5) * t297;
t180 = -qJD(5) * t352 + t290 * t360 - t297 * t344;
t113 = Ifges(6,1) * t179 - Ifges(6,4) * t180 + Ifges(6,5) * t343;
t423 = t113 / 0.2e1;
t390 = Ifges(7,4) * t292;
t262 = Ifges(7,2) * t296 + t390;
t389 = Ifges(7,4) * t296;
t327 = -Ifges(7,2) * t292 + t389;
t146 = -t262 * t355 + (Ifges(7,6) * t297 - t293 * t327) * qJD(5);
t422 = t146 / 0.2e1;
t264 = Ifges(7,1) * t292 + t389;
t328 = Ifges(7,1) * t296 - t390;
t147 = -t264 * t355 + (Ifges(7,5) * t297 - t293 * t328) * qJD(5);
t421 = t147 / 0.2e1;
t226 = t290 * t297 - t352;
t181 = -t226 * t292 + t296 * t379;
t420 = t181 / 0.2e1;
t182 = t226 * t296 + t292 * t379;
t419 = t182 / 0.2e1;
t208 = Ifges(7,6) * t293 + t297 * t327;
t418 = t208 / 0.2e1;
t209 = Ifges(7,5) * t293 + t297 * t328;
t417 = t209 / 0.2e1;
t356 = qJD(6) * t296;
t283 = Ifges(7,5) * t356;
t357 = qJD(6) * t292;
t416 = -Ifges(7,6) * t357 / 0.2e1 + t283 / 0.2e1;
t415 = (-Ifges(6,5) * t293 - Ifges(6,6) * t297) * qJD(5) / 0.2e1;
t243 = t327 * qJD(6);
t414 = t243 / 0.2e1;
t245 = t328 * qJD(6);
t413 = t245 / 0.2e1;
t391 = Ifges(6,4) * t297;
t246 = (-Ifges(6,1) * t293 - t391) * qJD(5);
t412 = t246 / 0.2e1;
t411 = Ifges(7,5) * t405 + Ifges(7,6) * t403;
t410 = Ifges(6,5) * t297 / 0.2e1 - Ifges(6,6) * t293 / 0.2e1;
t409 = t262 / 0.2e1;
t408 = t264 / 0.2e1;
t392 = Ifges(6,4) * t293;
t265 = Ifges(6,1) * t297 - t392;
t407 = t265 / 0.2e1;
t406 = -t292 / 0.2e1;
t404 = -t296 / 0.2e1;
t402 = m(7) * t293;
t400 = pkin(5) * t293;
t398 = pkin(12) * t297;
t73 = Ifges(5,1) * t332 - Ifges(5,4) * t128 + Ifges(5,5) * t127;
t74 = Ifges(4,5) * t128 - Ifges(4,6) * t127 + Ifges(4,3) * t332;
t396 = t73 + t74;
t395 = mrSges(7,3) * t297;
t394 = Ifges(4,4) * t294;
t393 = Ifges(4,4) * t298;
t388 = Ifges(5,6) * t294;
t387 = Ifges(5,6) * t298;
t386 = Ifges(7,6) * t292;
t385 = t129 * mrSges(4,1);
t384 = t129 * mrSges(4,2);
t222 = -pkin(9) * t346 + t273;
t383 = t222 * mrSges(3,2);
t223 = t231 * qJD(2);
t382 = t223 * mrSges(3,1);
t20 = -t293 * t54 + t297 * t52;
t18 = -pkin(5) * t168 - t20;
t381 = qJD(5) * t18;
t373 = t293 * t426;
t230 = pkin(10) * t378 + t279;
t368 = qJD(3) * t288;
t114 = t166 * t297 - t185 * t293;
t105 = -pkin(5) * t379 - t114;
t364 = qJD(5) * t105;
t363 = qJD(5) * t292;
t361 = qJD(5) * t296;
t254 = qJ(4) - t398 + t400;
t196 = t296 * t254 + t292 * t373;
t359 = qJD(6) * t196;
t197 = t292 * t254 - t296 * t373;
t358 = qJD(6) * t197;
t90 = t131 * t296 + t168 * t292;
t27 = -qJD(6) * t90 + t128 * t296 - t292 * t78;
t89 = -t131 * t292 + t168 * t296;
t28 = qJD(6) * t89 + t128 * t292 + t296 * t78;
t7 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t77;
t23 = Ifges(6,4) * t78 - Ifges(6,2) * t77 + Ifges(6,6) * t128;
t354 = t7 / 0.2e1 - t23 / 0.2e1;
t22 = Ifges(6,5) * t78 - Ifges(6,6) * t77 + Ifges(6,3) * t128;
t36 = Ifges(7,5) * t90 + Ifges(7,6) * t89 - Ifges(7,3) * t319;
t60 = Ifges(6,4) * t131 + Ifges(6,2) * t319 + Ifges(6,6) * t168;
t353 = t36 / 0.2e1 - t60 / 0.2e1;
t112 = Ifges(6,4) * t179 - Ifges(6,2) * t180 + Ifges(6,6) * t343;
t103 = qJD(6) * t181 + t179 * t296 + t292 * t343;
t104 = -qJD(6) * t182 - t179 * t292 + t296 * t343;
t47 = Ifges(7,5) * t103 + Ifges(7,6) * t104 + Ifges(7,3) * t180;
t349 = t47 / 0.2e1 - t112 / 0.2e1;
t141 = Ifges(6,4) * t226 - Ifges(6,2) * t225 + Ifges(6,6) * t379;
t98 = Ifges(7,5) * t182 + Ifges(7,6) * t181 + Ifges(7,3) * t225;
t348 = t98 / 0.2e1 - t141 / 0.2e1;
t88 = t298 * t162 + t198 * t379 + t153;
t111 = Ifges(6,5) * t179 - Ifges(6,6) * t180 + Ifges(6,3) * t343;
t204 = -t290 * qJ(4) - t230;
t342 = t290 * t366;
t340 = t426 * t360;
t309 = t292 * t355 + t293 * t361;
t145 = -Ifges(7,5) * t309 + Ifges(7,6) * t308 + Ifges(7,3) * t360;
t244 = (-Ifges(6,2) * t297 - t392) * qJD(5);
t338 = t145 / 0.2e1 - t244 / 0.2e1;
t207 = Ifges(7,3) * t293 + (Ifges(7,5) * t296 - t386) * t297;
t263 = -Ifges(6,2) * t293 + t391;
t337 = t207 / 0.2e1 - t263 / 0.2e1;
t108 = t128 * mrSges(5,1) + mrSges(5,2) * t332;
t233 = qJD(4) + (pkin(5) * t297 + pkin(12) * t293) * qJD(5);
t132 = t292 * t233 - t296 * t340 + t359;
t335 = t132 - t359;
t133 = t296 * t233 + t292 * t340 - t358;
t334 = -t133 - t358;
t82 = -t224 * qJ(4) - t88;
t184 = pkin(4) * t378 - t204;
t19 = pkin(12) * t168 + t397;
t56 = -pkin(4) * t167 - t82;
t35 = -pkin(5) * t319 - pkin(12) * t131 + t56;
t10 = -t19 * t292 + t296 * t35;
t42 = -t162 * t367 + t175 * t342 + t298 * t177 + t178 * t375 + t198 * t343 + t199 * t379;
t39 = -qJ(4) * t332 - t224 * qJD(4) - t42;
t29 = -pkin(4) * t127 - t39;
t13 = pkin(5) * t77 - pkin(12) * t78 + t29;
t5 = t293 * t31 + t297 * t33 + t52 * t360 - t362 * t54;
t3 = pkin(12) * t128 + t5;
t1 = qJD(6) * t10 + t13 * t292 + t296 * t3;
t11 = t19 * t296 + t292 * t35;
t2 = -qJD(6) * t11 + t13 * t296 - t292 * t3;
t330 = t1 * t296 - t2 * t292;
t259 = t293 * mrSges(6,1) + t297 * mrSges(6,2);
t258 = -mrSges(7,1) * t296 + mrSges(7,2) * t292;
t329 = mrSges(7,1) * t292 + mrSges(7,2) * t296;
t14 = -mrSges(7,2) * t77 + mrSges(7,3) * t27;
t15 = mrSges(7,1) * t77 - mrSges(7,3) * t28;
t326 = t296 * t14 - t292 * t15;
t66 = t297 * t159 + t166 * t360 - t185 * t362 + t293 * t187;
t62 = pkin(12) * t343 + t66;
t106 = pkin(12) * t379 + t437;
t118 = pkin(5) * t225 - pkin(12) * t226 + t184;
t64 = -t106 * t292 + t118 * t296;
t272 = pkin(2) * t342;
t284 = t290 * qJD(4);
t165 = -t344 * t425 + t272 + t284;
t91 = pkin(5) * t180 - pkin(12) * t179 + t165;
t16 = qJD(6) * t64 + t292 * t91 + t296 * t62;
t65 = t106 * t296 + t118 * t292;
t17 = -qJD(6) * t65 - t292 * t62 + t296 * t91;
t325 = t16 * t296 - t17 * t292;
t85 = mrSges(7,1) * t180 - mrSges(7,3) * t103;
t86 = -mrSges(7,2) * t180 + mrSges(7,3) * t104;
t324 = -t292 * t85 + t296 * t86;
t220 = -pkin(10) * t344 + t272;
t37 = Ifges(7,4) * t90 + Ifges(7,2) * t89 - Ifges(7,6) * t319;
t38 = Ifges(7,1) * t90 + Ifges(7,4) * t89 - Ifges(7,5) * t319;
t316 = t37 * t406 + t38 * t403;
t100 = Ifges(7,1) * t182 + Ifges(7,4) * t181 + Ifges(7,5) * t225;
t99 = Ifges(7,4) * t182 + Ifges(7,2) * t181 + Ifges(7,6) * t225;
t315 = t100 * t403 + t406 * t99;
t12 = -mrSges(7,1) * t27 + mrSges(7,2) * t28;
t45 = mrSges(6,1) * t128 - mrSges(6,3) * t78;
t311 = t45 + m(6) * (t6 + t436) - t12;
t143 = mrSges(6,1) * t343 - mrSges(6,3) * t179;
t55 = -mrSges(7,1) * t104 + mrSges(7,2) * t103;
t310 = t143 + m(6) * (t67 + t438) - t55;
t44 = -mrSges(6,2) * t128 - mrSges(6,3) * t77;
t46 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t93 = mrSges(6,1) * t168 - mrSges(6,3) * t131;
t304 = t44 + m(6) * (-qJD(5) * t20 + t5) + (t46 - t93) * qJD(5);
t119 = -mrSges(7,1) * t181 + mrSges(7,2) * t182;
t144 = -mrSges(6,2) * t343 - mrSges(6,3) * t180;
t191 = mrSges(6,1) * t379 - mrSges(6,3) * t226;
t303 = t144 + m(6) * (-qJD(5) * t114 + t66) + (t119 - t191) * qJD(5);
t302 = -t10 * t356 - t11 * t357 + t330;
t301 = -t356 * t64 - t357 * t65 + t325;
t270 = Ifges(3,5) * t345;
t269 = Ifges(4,5) * t343;
t268 = Ifges(5,5) * t344;
t248 = mrSges(7,1) * t293 - t296 * t395;
t247 = -mrSges(7,2) * t293 - t292 * t395;
t240 = (mrSges(6,1) * t297 - mrSges(6,2) * t293) * qJD(5);
t239 = t329 * qJD(6);
t238 = mrSges(5,1) * t379 + mrSges(5,2) * t290;
t237 = -mrSges(5,1) * t378 - mrSges(5,3) * t290;
t236 = -mrSges(4,2) * t290 + mrSges(4,3) * t378;
t235 = mrSges(4,1) * t290 - mrSges(4,3) * t379;
t232 = t329 * t297;
t229 = -pkin(9) * t377 + t281;
t228 = pkin(2) * t374 - t275;
t227 = (mrSges(5,2) * t298 - mrSges(5,3) * t294) * t288;
t221 = t230 * qJD(3);
t219 = (Ifges(4,1) * t298 - t394) * t368;
t218 = (-Ifges(4,2) * t294 + t393) * t368;
t217 = -Ifges(5,4) * t343 + t268;
t216 = -Ifges(4,6) * t344 + t269;
t215 = (-Ifges(5,2) * t298 + t388) * t368;
t214 = (Ifges(5,3) * t294 - t387) * t368;
t213 = (mrSges(4,1) * t294 + mrSges(4,2) * t298) * t368;
t212 = (-mrSges(5,2) * t294 - mrSges(5,3) * t298) * t368;
t206 = t290 * t347 + t275;
t205 = (-pkin(3) * t298 + t336) * t288;
t203 = Ifges(5,4) * t290 + (-Ifges(5,2) * t294 - t387) * t288;
t202 = Ifges(5,5) * t290 + (-Ifges(5,3) * t298 - t388) * t288;
t201 = Ifges(4,5) * t290 + (Ifges(4,1) * t294 + t393) * t288;
t200 = Ifges(4,6) * t290 + (Ifges(4,2) * t298 + t394) * t288;
t195 = -t220 - t284;
t194 = -mrSges(7,2) * t360 + mrSges(7,3) * t308;
t193 = mrSges(7,1) * t360 + mrSges(7,3) * t309;
t190 = -mrSges(6,2) * t379 - mrSges(6,3) * t225;
t189 = t271 + (-qJ(4) * t366 - t365) * t288;
t160 = -mrSges(7,1) * t308 - mrSges(7,2) * t309;
t152 = mrSges(6,1) * t225 + mrSges(6,2) * t226;
t142 = Ifges(6,1) * t226 - Ifges(6,4) * t225 + Ifges(6,5) * t379;
t140 = Ifges(6,5) * t226 - Ifges(6,6) * t225 + Ifges(6,3) * t379;
t139 = mrSges(7,1) * t225 - mrSges(7,3) * t182;
t138 = -mrSges(7,2) * t225 + mrSges(7,3) * t181;
t137 = mrSges(5,1) * t168 + mrSges(5,2) * t224;
t136 = mrSges(5,1) * t167 - mrSges(5,3) * t224;
t135 = mrSges(4,1) * t224 - mrSges(4,3) * t168;
t134 = -mrSges(4,2) * t224 - mrSges(4,3) * t167;
t117 = mrSges(6,1) * t180 + mrSges(6,2) * t179;
t116 = -mrSges(5,2) * t167 - mrSges(5,3) * t168;
t110 = mrSges(4,1) * t332 - mrSges(4,3) * t128;
t109 = -mrSges(4,2) * t332 - mrSges(4,3) * t127;
t107 = mrSges(5,1) * t127 - mrSges(5,3) * t332;
t97 = Ifges(4,1) * t168 - Ifges(4,4) * t167 + Ifges(4,5) * t224;
t96 = Ifges(4,4) * t168 - Ifges(4,2) * t167 + Ifges(4,6) * t224;
t95 = Ifges(5,4) * t224 - Ifges(5,2) * t168 + Ifges(5,6) * t167;
t94 = Ifges(5,5) * t224 - Ifges(5,6) * t168 + Ifges(5,3) * t167;
t92 = -mrSges(6,2) * t168 + mrSges(6,3) * t319;
t84 = -mrSges(6,1) * t319 + mrSges(6,2) * t131;
t83 = -pkin(3) * t224 - t87;
t81 = -mrSges(5,2) * t127 - mrSges(5,3) * t128;
t80 = mrSges(4,1) * t127 + mrSges(4,2) * t128;
t79 = pkin(3) * t167 + t312;
t76 = Ifges(4,1) * t128 - Ifges(4,4) * t127 + Ifges(4,5) * t332;
t75 = Ifges(4,4) * t128 - Ifges(4,2) * t127 + Ifges(4,6) * t332;
t72 = Ifges(5,4) * t332 - Ifges(5,2) * t128 + Ifges(5,6) * t127;
t71 = Ifges(5,5) * t332 - Ifges(5,6) * t128 + Ifges(5,3) * t127;
t61 = Ifges(6,1) * t131 + Ifges(6,4) * t319 + Ifges(6,5) * t168;
t59 = Ifges(6,5) * t131 + Ifges(6,6) * t319 + Ifges(6,3) * t168;
t58 = -mrSges(7,1) * t319 - mrSges(7,3) * t90;
t57 = mrSges(7,2) * t319 + mrSges(7,3) * t89;
t49 = Ifges(7,1) * t103 + Ifges(7,4) * t104 + Ifges(7,5) * t180;
t48 = Ifges(7,4) * t103 + Ifges(7,2) * t104 + Ifges(7,6) * t180;
t43 = (t178 * t290 + t199 * t288) * t298 + t333;
t41 = pkin(3) * t127 + t305;
t40 = (-pkin(3) * t346 - t380) * t288 + t307;
t34 = mrSges(6,1) * t77 + mrSges(6,2) * t78;
t9 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t77;
t8 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t77;
t21 = [-(t7 - t23) * t319 + (t59 + t97 - t95) * t128 + (t94 - t96) * t127 + (t36 - t60) * t77 + (0.2e1 * (t222 * t299 + t223 * t295) * mrSges(3,3) + ((t229 * t430 + Ifges(3,5) * t291 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t299) * t439) * t299 + (t231 * t430 - 0.2e1 * Ifges(3,6) * t291 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t295 + (Ifges(3,1) - Ifges(3,2)) * t299) * t439 + ((Ifges(5,1) + Ifges(4,3)) * t224 + (-Ifges(5,4) + Ifges(4,5)) * t168 + (Ifges(5,5) - Ifges(4,6)) * t167) * t288) * t295) * qJD(2)) * t289 + (t270 - 0.2e1 * t382 - 0.2e1 * t383) * t291 + 0.2e1 * m(3) * (t222 * t231 - t223 * t229) + (t20 * t6 + t29 * t56 + t397 * t5) * t433 + 0.2e1 * t397 * t44 + (t22 - t72 + t76 + 0.2e1 * t384) * t168 + (t71 - t75 + 0.2e1 * t385) * t167 + t396 * t224 + (t1 * t11 + t10 * t2 + t18 * t4) * t432 + (t39 * t82 + t40 * t83 + t41 * t79) * t434 + (t121 * t129 + t42 * t88 + t43 * t87) * t435 + 0.2e1 * t11 * t14 + 0.2e1 * t10 * t15 + 0.2e1 * t18 * t12 + t27 * t37 + t28 * t38 + 0.2e1 * t20 * t45 + 0.2e1 * t4 * t46 + 0.2e1 * t56 * t34 + 0.2e1 * t1 * t57 + 0.2e1 * t2 * t58 + t78 * t61 + 0.2e1 * t79 * t81 + 0.2e1 * t29 * t84 + t89 * t8 + t90 * t9 + 0.2e1 * t5 * t92 + 0.2e1 * t6 * t93 + 0.2e1 * t82 * t107 + 0.2e1 * t83 * t108 + 0.2e1 * t88 * t109 + 0.2e1 * t87 * t110 + 0.2e1 * t41 * t116 + 0.2e1 * t121 * t80 + t131 * t24 + 0.2e1 * t42 * t134 + 0.2e1 * t43 * t135 + 0.2e1 * t39 * t136 + 0.2e1 * t40 * t137; -t349 * t319 + m(6) * (t114 * t6 + t165 * t56 + t184 * t29 + t20 * t67 + t397 * t66 + t437 * t5) + t437 * t44 + t348 * t77 + t353 * t180 - t383 + (t137 - t135) * t221 + m(4) * (t220 * t88 - t221 * t87 + t228 * t43 + t230 * t42) + (t111 / 0.2e1 - t215 / 0.2e1 + t219 / 0.2e1) * t168 + (t73 / 0.2e1 + t74 / 0.2e1) * t290 + (t140 / 0.2e1 + t201 / 0.2e1 - t203 / 0.2e1) * t128 - Ifges(3,6) * t346 + (t214 / 0.2e1 - t218 / 0.2e1) * t167 + m(7) * (t1 * t65 + t10 * t17 + t105 * t4 + t11 * t16 + t18 * t63 + t2 * t64) + m(5) * (t189 * t79 + t195 * t82 + t204 * t39 + t205 * t41 + t206 * t40 + t221 * t83) + t354 * t225 + t397 * t144 + (-t200 / 0.2e1 + t202 / 0.2e1) * t127 - t382 + ((-m(4) * t129 - t80) * pkin(2) + (-t71 / 0.2e1 + t75 / 0.2e1 - t385) * t298 + (t22 / 0.2e1 - t72 / 0.2e1 + t76 / 0.2e1 + t384) * t294 + (((-Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1) * t298 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t294) * t288 + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t290) * t346 + ((-t88 * mrSges(4,3) + t82 * mrSges(5,1) + t94 / 0.2e1 - t96 / 0.2e1) * t294 + (t83 * mrSges(5,1) - t87 * mrSges(4,3) - t95 / 0.2e1 + t97 / 0.2e1 + t59 / 0.2e1) * t298) * qJD(3)) * t288 + t270 + t9 * t419 + t8 * t420 + t131 * t423 + t226 * t427 + (t216 / 0.2e1 + t217 / 0.2e1) * t224 + t18 * t55 + t16 * t57 + t17 * t58 + t63 * t46 + t64 * t15 + t65 * t14 + t10 * t85 + t11 * t86 + t89 * t48 / 0.2e1 + t90 * t49 / 0.2e1 + t66 * t92 + t67 * t93 + t27 * t99 / 0.2e1 + t28 * t100 / 0.2e1 + t103 * t38 / 0.2e1 + t104 * t37 / 0.2e1 + t105 * t12 + t114 * t45 + t56 * t117 + t4 * t119 + t1 * t138 + t2 * t139 + t78 * t142 / 0.2e1 + t20 * t143 + t29 * t152 + t165 * t84 + t179 * t61 / 0.2e1 + t184 * t34 + t189 * t116 + t5 * t190 + t6 * t191 + t195 * t136 + t204 * t107 + t205 * t81 + t206 * t108 + t79 * t212 + t121 * t213 + t220 * t134 + t41 * t227 + t228 * t110 + t230 * t109 + t43 * t235 + t42 * t236 + t39 * t237 + t40 * t238; (t114 * t67 + t165 * t184 + t437 * t66) * t433 + 0.2e1 * t437 * t144 + (t47 - t112) * t225 + (t98 - t141) * t180 + 0.2e1 * (-t235 + t238) * t221 + (t216 + t217) * t290 + (t105 * t63 + t16 * t65 + t17 * t64) * t432 + (t189 * t205 + t195 * t204 + t206 * t221) * t434 + (t220 * t230 - t221 * t228) * t435 + (-0.2e1 * pkin(2) * t213 + (-t214 + t218) * t298 + (t111 - t215 + t219) * t294 + ((t204 * t431 + t230 * t429 - t200 + t202) * t294 + (t206 * t431 + t228 * t429 + t140 + t201 - t203) * t298) * qJD(3)) * t288 + 0.2e1 * t64 * t85 + 0.2e1 * t65 * t86 + t103 * t100 + t104 * t99 + 0.2e1 * t105 * t55 + 0.2e1 * t63 * t119 + 0.2e1 * t16 * t138 + 0.2e1 * t17 * t139 + 0.2e1 * t114 * t143 + 0.2e1 * t165 * t152 + t179 * t142 + t181 * t48 + t182 * t49 + 0.2e1 * t184 * t117 + 0.2e1 * t66 * t190 + 0.2e1 * t67 * t191 + 0.2e1 * t205 * t212 + t226 * t113 + 0.2e1 * t189 * t227 + 0.2e1 * t220 * t236 + 0.2e1 * t195 * t237; -t338 * t319 + (-t5 * mrSges(6,3) + (t20 * mrSges(6,3) - t61 / 0.2e1 - t316) * qJD(5) - (m(7) * t381 + t304) * t426 + t354) * t293 + (-t6 * mrSges(6,3) + t8 * t406 + t9 * t403 + t427 + (t37 * t404 + t38 * t406) * qJD(6) + (-mrSges(6,3) * t397 + t353) * qJD(5) - (qJD(5) * t92 + t311 - t428) * t426) * t297 + t396 + (t84 - t136) * qJD(4) + (t34 - t107) * qJ(4) + m(7) * (t197 * t1 + t133 * t10 + t132 * t11 + t196 * t2) + m(6) * (qJ(4) * t29 + qJD(4) * t56) + m(5) * (-pkin(3) * t40 - qJ(4) * t39 - qJD(4) * t82) + t337 * t77 + t78 * t407 + t128 * t410 + t131 * t412 + t168 * t415 + t28 * t417 + t27 * t418 + t90 * t421 + t89 * t422 - t39 * mrSges(5,3) + t40 * mrSges(5,2) - t42 * mrSges(4,2) + t43 * mrSges(4,1) - pkin(3) * t108 + t132 * t57 + t133 * t58 + t18 * t160 + t10 * t193 + t11 * t194 + t196 * t15 + t197 * t14 + t4 * t232 + t56 * t240 + t1 * t247 + t2 * t248 + t29 * t259; (-t67 * mrSges(6,3) + t48 * t406 + t49 * t403 + t423 + (t100 * t406 + t404 * t99) * qJD(6) + (-mrSges(6,3) * t437 + t348) * qJD(5) - (qJD(5) * t190 + t310 - t424) * t426) * t297 + (-t66 * mrSges(6,3) + (t114 * mrSges(6,3) - t142 / 0.2e1 - t315) * qJD(5) - (m(7) * t364 + t303) * t426 + t349) * t293 + (mrSges(5,2) - mrSges(4,1)) * t221 + (t152 - t237) * qJD(4) + m(6) * (qJ(4) * t165 + qJD(4) * t184) + m(7) * (t132 * t65 + t133 * t64 + t197 * t16 + t196 * t17) + t337 * t180 + t338 * t225 + m(5) * (-pkin(3) * t221 - qJ(4) * t195 - qJD(4) * t204) + t268 + t269 + (t294 * t415 + ((-mrSges(5,1) * qJ(4) - Ifges(4,6)) * t294 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t410) * t298) * qJD(3)) * t288 + t179 * t407 + t226 * t412 + t103 * t417 + t104 * t418 + t147 * t419 + t146 * t420 + qJ(4) * t117 + t132 * t138 + t133 * t139 + t105 * t160 + t64 * t193 + t65 * t194 - t195 * mrSges(5,3) + t196 * t85 + t197 * t86 - t220 * mrSges(4,2) + t63 * t232 + t184 * t240 + t16 * t247 + t17 * t248 + t165 * t259; 0.2e1 * t133 * t248 + 0.2e1 * t196 * t193 + 0.2e1 * t132 * t247 + 0.2e1 * t197 * t194 + (t197 * t132 + t196 * t133) * t432 + 0.2e1 * qJ(4) * t240 + 0.2e1 * (mrSges(5,3) + t259 + (m(5) + m(6)) * qJ(4)) * qJD(4) + (t145 - t244 + (t208 * t292 - t209 * t296 - 0.2e1 * t232 * t426 - t265) * qJD(5)) * t293 + (-t292 * t146 + t296 * t147 + 0.2e1 * t426 * t160 + t246 + (-t208 * t296 - t209 * t292) * qJD(6) + (-0.2e1 * t402 * t426 ^ 2 + t207 - t263) * qJD(5)) * t297; m(5) * t40 + ((-t292 * t58 + t296 * t57 + t92) * qJD(5) + m(7) * (-t10 * t363 + t11 * t361 - t4) + t311) * t297 + ((-t292 * t57 - t296 * t58) * qJD(6) + m(7) * (t302 + t381) + t304 + t326) * t293 + t108; m(5) * t221 + mrSges(5,1) * t343 + ((t138 * t296 - t139 * t292 + t190) * qJD(5) + m(7) * (t361 * t65 - t363 * t64 - t63) + t310) * t297 + ((-t138 * t292 - t139 * t296) * qJD(6) + m(7) * (t301 + t364) + t303 + t324) * t293; (-t160 + (m(7) * (-t196 * t292 + t197 * t296) + t296 * t247 - t292 * t248) * qJD(5)) * t297 + (m(7) * (t132 * t296 - t133 * t292 - t196 * t356 - t197 * t357 + 0.2e1 * t340) - t247 * t357 + t296 * t194 - t248 * t356 - t292 * t193 + qJD(5) * t232) * t293; 0.2e1 * (-0.1e1 + t370) * t360 * t402; -t5 * mrSges(6,2) + t6 * mrSges(6,1) + t18 * t239 - t319 * t416 + t89 * t414 + t90 * t413 + t4 * t258 + t77 * t411 + t27 * t409 + t28 * t408 + t9 * t405 + t8 * t403 + t316 * qJD(6) + (-t12 - t428) * pkin(5) + ((-t10 * t296 - t11 * t292) * qJD(6) + t330) * mrSges(7,3) + (m(7) * t302 - t356 * t58 - t357 * t57 + t326) * pkin(12) + t22; -t66 * mrSges(6,2) + t67 * mrSges(6,1) + t105 * t239 + t225 * t416 + t181 * t414 + t182 * t413 + t63 * t258 + t180 * t411 + t104 * t409 + t103 * t408 + t49 * t405 + t48 * t403 + t315 * qJD(6) + (-t55 - t424) * pkin(5) + ((-t292 * t65 - t296 * t64) * qJD(6) + t325) * mrSges(7,3) + (m(7) * t301 - t138 * t357 - t139 * t356 + t324) * pkin(12) + t111; -pkin(5) * t160 + (t416 + (-Ifges(6,5) - (-m(7) * pkin(5) - mrSges(6,1) + t258) * t426) * qJD(5)) * t293 + (-t264 * t362 / 0.2e1 + qJD(6) * t417 + t422 + t335 * mrSges(7,3) + (m(7) * t335 - qJD(6) * t248 + t194) * pkin(12)) * t296 + (t362 * t409 - qJD(6) * t208 / 0.2e1 + t421 + t334 * mrSges(7,3) + (m(7) * t334 - qJD(6) * t247 - t193) * pkin(12)) * t292 + (t243 * t406 + t245 * t403 + t426 * t239 + (t262 * t404 + t264 * t406) * qJD(6) + (mrSges(6,2) * t426 - Ifges(6,6) + t411) * qJD(5)) * t297; -t297 * t239 + (t293 * t258 + m(7) * (t370 * t398 - t400) + t370 * t395 - t259) * qJD(5); -0.2e1 * pkin(5) * t239 + t243 * t296 + t245 * t292 + (-t262 * t292 + t264 * t296) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t47; mrSges(7,1) * t133 - mrSges(7,2) * t132 + t145; (t293 * t357 - t296 * t360) * mrSges(7,2) + (-t292 * t360 - t293 * t356) * mrSges(7,1); t283 + (pkin(12) * t258 - t386) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
