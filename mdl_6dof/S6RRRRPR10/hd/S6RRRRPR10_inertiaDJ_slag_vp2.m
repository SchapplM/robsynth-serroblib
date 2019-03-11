% Calculate time derivative of joint inertia matrix for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:01:02
% EndTime: 2019-03-09 23:01:24
% DurationCPUTime: 9.20s
% Computational Cost: add. (11226->652), mult. (28791->927), div. (0->0), fcn. (27973->10), ass. (0->298)
t435 = Ifges(5,1) + Ifges(7,3);
t434 = Ifges(6,1) + Ifges(5,3);
t423 = Ifges(5,5) - Ifges(6,4);
t422 = -Ifges(5,6) + Ifges(6,5);
t433 = -mrSges(5,1) + mrSges(6,2);
t246 = sin(qJ(4));
t383 = cos(qJ(4));
t384 = cos(qJ(3));
t296 = t383 * t384;
t382 = sin(qJ(3));
t206 = t246 * t382 - t296;
t245 = sin(qJ(6));
t344 = qJD(6) * t245;
t309 = -t344 / 0.2e1;
t207 = t246 * t384 + t382 * t383;
t409 = qJD(3) + qJD(4);
t166 = t409 * t207;
t248 = cos(qJ(6));
t352 = t248 * t166;
t432 = t206 * t309 + t352 / 0.2e1;
t359 = t206 * t248;
t316 = t359 / 0.2e1;
t353 = t245 * t166;
t431 = qJD(6) * t316 + t353 / 0.2e1;
t343 = qJD(6) * t248;
t276 = t206 * t343 + t353;
t275 = t206 * t344 - t352;
t339 = t384 * pkin(9);
t222 = pkin(10) * t384 + t339;
t267 = qJD(3) * t222;
t336 = t382 * pkin(9);
t270 = -pkin(10) * t382 - t336;
t268 = t246 * t270;
t311 = qJD(4) * t383;
t119 = t222 * t311 + t383 * t267 + t268 * t409;
t310 = t382 * qJD(3);
t165 = (qJD(4) * t382 + t310) * t246 - t409 * t296;
t372 = t165 * mrSges(6,1);
t430 = m(6) * t119 - t372;
t349 = t245 ^ 2 + t248 ^ 2;
t429 = (-mrSges(5,2) * t383 + (-t349 * mrSges(7,3) + t433) * t246) * pkin(3) * qJD(4);
t378 = mrSges(5,3) + mrSges(6,1);
t428 = 0.2e1 * t378;
t244 = sin(pkin(6));
t247 = sin(qJ(2));
t347 = qJD(2) * t247;
t323 = t244 * t347;
t364 = cos(pkin(6));
t289 = t364 * t382;
t355 = t244 * t247;
t198 = t355 * t384 + t289;
t261 = t198 * qJD(3);
t249 = cos(qJ(2));
t346 = qJD(2) * t249;
t294 = t382 * t346;
t256 = t244 * t294 + t261;
t290 = t364 * t384;
t264 = t355 * t382 - t290;
t295 = t384 * t346;
t257 = -qJD(3) * t264 + t244 * t295;
t258 = t383 * t264;
t345 = qJD(4) * t246;
t75 = qJD(4) * t258 + t198 * t345 + t246 * t256 - t257 * t383;
t262 = t246 * t264;
t76 = -qJD(4) * t262 + t198 * t311 + t246 * t257 + t256 * t383;
t146 = t198 * t246 + t258;
t354 = t244 * t249;
t278 = -t146 * t245 + t248 * t354;
t40 = qJD(6) * t278 - t245 * t323 + t248 * t76;
t124 = t146 * t248 + t245 * t354;
t41 = qJD(6) * t124 + t245 * t76 + t248 * t323;
t8 = Ifges(7,5) * t41 + Ifges(7,6) * t40 - Ifges(7,3) * t75;
t427 = -Ifges(5,1) * t75 - Ifges(5,4) * t76 + Ifges(5,5) * t323 + t8;
t239 = -pkin(3) * t384 - pkin(2);
t271 = -t207 * qJ(5) + t239;
t394 = pkin(4) + pkin(11);
t133 = t206 * t394 + t271;
t208 = t383 * t270;
t180 = t222 * t246 - t208;
t144 = pkin(5) * t207 + t180;
t81 = t133 * t248 + t144 * t245;
t361 = qJD(6) * t81;
t241 = pkin(3) * t310;
t274 = qJ(5) * t165 - qJD(5) * t207 + t241;
t68 = t166 * t394 + t274;
t87 = -t165 * pkin(5) + t119;
t23 = -t245 * t68 + t248 * t87 - t361;
t80 = -t133 * t245 + t144 * t248;
t22 = qJD(6) * t80 + t245 * t87 + t248 * t68;
t365 = t245 * t22;
t411 = -t23 * t248 - t365;
t426 = m(7) * ((-t245 * t80 + t248 * t81) * qJD(6) - t411);
t147 = t198 * t383 - t262;
t324 = pkin(1) * t364;
t234 = t247 * t324;
t204 = pkin(8) * t354 + t234;
t191 = pkin(9) * t364 + t204;
t192 = (-pkin(2) * t249 - pkin(9) * t247 - pkin(1)) * t244;
t135 = -t191 * t382 + t384 * t192;
t100 = -pkin(3) * t354 - t198 * pkin(10) + t135;
t136 = t384 * t191 + t382 * t192;
t114 = -pkin(10) * t264 + t136;
t52 = t100 * t383 - t246 * t114;
t46 = pkin(4) * t354 - t52;
t34 = t147 * pkin(5) + pkin(11) * t354 + t46;
t203 = -pkin(8) * t355 + t249 * t324;
t190 = -pkin(2) * t364 - t203;
t148 = pkin(3) * t264 + t190;
t251 = -t147 * qJ(5) + t148;
t42 = t146 * t394 + t251;
t19 = -t245 * t42 + t248 * t34;
t20 = t245 * t34 + t248 * t42;
t322 = t244 * t346;
t24 = pkin(3) * t256 + t76 * pkin(4) + pkin(8) * t322 + t75 * qJ(5) + qJD(2) * t234 - t147 * qJD(5);
t13 = t76 * pkin(11) + t24;
t362 = qJD(6) * t20;
t348 = qJD(2) * t244;
t195 = (pkin(2) * t247 - pkin(9) * t249) * t348;
t196 = t203 * qJD(2);
t288 = t384 * t195 - t196 * t382;
t312 = qJD(3) * t384;
t58 = pkin(3) * t323 - pkin(10) * t257 - t191 * t312 - t192 * t310 + t288;
t84 = -t191 * t310 + t192 * t312 + t382 * t195 + t384 * t196;
t61 = -pkin(10) * t256 + t84;
t277 = t100 * t345 + t114 * t311 + t246 * t61 - t383 * t58;
t4 = -t75 * pkin(5) - t323 * t394 + t277;
t2 = -t13 * t245 + t248 * t4 - t362;
t1 = qJD(6) * t19 + t13 * t248 + t245 * t4;
t379 = t245 * t1;
t410 = -t2 * t248 - t379;
t425 = m(7) * ((-t19 * t245 + t20 * t248) * qJD(6) - t410);
t424 = mrSges(5,2) - mrSges(6,3);
t421 = -Ifges(5,4) * t146 - Ifges(5,5) * t354 - Ifges(7,5) * t278 + Ifges(7,6) * t124 + t147 * t435;
t420 = t196 * mrSges(3,2);
t216 = mrSges(7,1) * t245 + mrSges(7,2) * t248;
t419 = mrSges(6,3) + t216;
t54 = t276 * Ifges(7,5) - t275 * Ifges(7,6) - Ifges(7,3) * t165;
t418 = -Ifges(5,1) * t165 - Ifges(5,4) * t166 + t54;
t284 = Ifges(7,5) * t245 + Ifges(7,6) * t248;
t417 = t435 * t207 + (-Ifges(5,4) + t284) * t206;
t360 = t206 * t245;
t156 = mrSges(7,1) * t207 - mrSges(7,3) * t360;
t157 = -mrSges(7,2) * t207 + mrSges(7,3) * t359;
t416 = -t156 * t344 + t157 * t343;
t82 = -mrSges(7,2) * t147 + mrSges(7,3) * t124;
t83 = mrSges(7,1) * t147 + mrSges(7,3) * t278;
t415 = t82 * t343 - t83 * t344;
t85 = -qJD(3) * t136 + t288;
t414 = -t382 * t85 + t384 * t84;
t413 = t323 * t434 + t422 * t76 - t423 * t75;
t412 = -t165 * t423 + t166 * t422;
t273 = -t100 * t311 + t114 * t345 - t246 * t58 - t383 * t61;
t11 = -t244 * (qJ(5) * t347 - qJD(5) * t249) + t273;
t64 = mrSges(6,1) * t76 - mrSges(6,3) * t323;
t408 = -m(6) * t11 - t64;
t14 = -pkin(4) * t323 + t277;
t65 = -t75 * mrSges(6,1) + mrSges(6,2) * t323;
t407 = m(6) * t14 + t65;
t129 = mrSges(6,1) * t146 + mrSges(6,3) * t354;
t53 = t246 * t100 + t383 * t114;
t45 = qJ(5) * t354 - t53;
t406 = m(6) * t45 + t129;
t130 = mrSges(6,1) * t147 - mrSges(6,2) * t354;
t132 = -mrSges(5,1) * t354 - mrSges(5,3) * t147;
t405 = m(6) * t46 + t130 - t132;
t403 = 2 * m(5);
t402 = 0.2e1 * m(6);
t401 = 0.2e1 * m(7);
t287 = mrSges(7,1) * t248 - mrSges(7,2) * t245;
t209 = t287 * qJD(6);
t399 = 0.2e1 * t209;
t398 = m(5) * pkin(3);
t397 = t40 / 0.2e1;
t396 = t41 / 0.2e1;
t393 = t124 / 0.2e1;
t392 = -t278 / 0.2e1;
t388 = -t284 * qJD(6) / 0.2e1;
t387 = -Ifges(7,5) * t248 / 0.2e1 + Ifges(7,6) * t245 / 0.2e1;
t386 = -t245 / 0.2e1;
t385 = t248 / 0.2e1;
t381 = pkin(3) * t246;
t377 = Ifges(3,4) * t247;
t376 = Ifges(3,4) * t249;
t375 = Ifges(7,4) * t245;
t374 = Ifges(7,4) * t248;
t371 = t165 * mrSges(5,3);
t370 = t166 * mrSges(6,1);
t369 = t166 * mrSges(5,3);
t368 = t206 * mrSges(6,1);
t367 = t206 * mrSges(5,3);
t363 = mrSges(4,3) * qJD(3);
t298 = pkin(3) * t311;
t228 = t298 + qJD(5);
t236 = qJ(5) + t381;
t358 = t228 * t236;
t338 = t383 * pkin(3);
t238 = -t338 - pkin(4);
t235 = -pkin(11) + t238;
t357 = t235 * t245;
t356 = t235 * t248;
t337 = t382 * pkin(3);
t335 = pkin(3) * t345;
t333 = Ifges(4,4) * t384;
t332 = Ifges(4,4) * t382;
t329 = mrSges(7,3) * t344;
t327 = mrSges(7,3) * t343;
t325 = Ifges(4,5) * t257 - Ifges(4,6) * t256 + Ifges(4,3) * t323;
t317 = t360 / 0.2e1;
t308 = -t343 / 0.2e1;
t302 = t245 * t335;
t301 = t248 * t335;
t300 = pkin(9) * t312;
t299 = pkin(9) * t310;
t293 = t323 / 0.2e1;
t286 = Ifges(7,1) * t245 + t374;
t285 = Ifges(7,2) * t248 + t375;
t118 = -t409 * t208 + t222 * t345 + t246 * t267;
t181 = t222 * t383 + t268;
t283 = -t181 * t118 + t180 * t119;
t282 = qJ(5) * t228 + qJD(5) * t236;
t281 = t349 * t335;
t280 = Ifges(3,5) * t322 - Ifges(3,6) * t323;
t279 = Ifges(4,5) * t312 - Ifges(4,6) * t310;
t269 = mrSges(4,1) * t382 + mrSges(4,2) * t384;
t266 = -t247 * t312 - t294;
t265 = -t247 * t310 + t295;
t212 = t285 * qJD(6);
t214 = t286 * qJD(6);
t218 = -Ifges(7,2) * t245 + t374;
t220 = Ifges(7,1) * t248 - t375;
t263 = -t248 * t214 + (-t218 * t248 - t220 * t245) * qJD(6) + t245 * t212;
t25 = mrSges(7,2) * t75 + mrSges(7,3) * t40;
t26 = -mrSges(7,1) * t75 - mrSges(7,3) * t41;
t255 = t245 * t25 + t248 * t26 + t425;
t92 = mrSges(7,2) * t165 - mrSges(7,3) * t275;
t93 = -mrSges(7,1) * t165 - mrSges(7,3) * t276;
t254 = t245 * t92 + t248 * t93 + t426;
t127 = Ifges(7,6) * t207 + t206 * t285;
t128 = Ifges(7,5) * t207 + t206 * t286;
t145 = -t206 * pkin(5) + t181;
t55 = Ifges(7,4) * t276 - Ifges(7,2) * t275 - Ifges(7,6) * t165;
t56 = Ifges(7,1) * t276 - Ifges(7,4) * t275 - Ifges(7,5) * t165;
t86 = -pkin(5) * t166 - t118;
t253 = t433 * t119 + t127 * t308 + t128 * t309 + t145 * t209 + t165 * t387 + t207 * t388 - t212 * t316 - t214 * t317 + t86 * t216 + t432 * t218 + t431 * t220 + t80 * t329 + t56 * t385 + t55 * t386 + t412;
t10 = Ifges(7,1) * t41 + Ifges(7,4) * t40 - Ifges(7,5) * t75;
t35 = -pkin(5) * t146 - t45;
t48 = -Ifges(7,4) * t278 + Ifges(7,2) * t124 + Ifges(7,6) * t147;
t49 = -Ifges(7,1) * t278 + Ifges(7,4) * t124 + Ifges(7,5) * t147;
t5 = -pkin(5) * t76 - t11;
t9 = Ifges(7,4) * t41 + Ifges(7,2) * t40 - Ifges(7,6) * t75;
t252 = -t277 * mrSges(5,1) + t273 * mrSges(5,2) + t14 * mrSges(6,2) - t11 * mrSges(6,3) + t10 * t385 + t19 * t329 + t35 * t209 + t5 * t216 + t218 * t397 + t220 * t396 + t49 * t309 + t9 * t386 - t212 * t393 - t214 * t392 + t147 * t388 + t48 * t308 + t75 * t387 + t413;
t221 = Ifges(4,1) * t382 + t333;
t219 = Ifges(4,2) * t384 + t332;
t215 = (Ifges(4,1) * t384 - t332) * qJD(3);
t213 = (-Ifges(4,2) * t382 + t333) * qJD(3);
t210 = t269 * qJD(3);
t197 = t204 * qJD(2);
t179 = -mrSges(4,1) * t354 - mrSges(4,3) * t198;
t178 = mrSges(4,2) * t354 - mrSges(4,3) * t264;
t173 = Ifges(5,4) * t207 - Ifges(5,2) * t206;
t172 = -Ifges(6,2) * t207 + Ifges(6,6) * t206;
t171 = -Ifges(6,6) * t207 + Ifges(6,3) * t206;
t170 = -mrSges(6,2) * t206 - mrSges(6,3) * t207;
t169 = mrSges(5,1) * t206 + mrSges(5,2) * t207;
t150 = t206 * pkin(4) + t271;
t149 = t287 * t206;
t143 = -t290 * t363 + (mrSges(4,1) * t347 - mrSges(4,3) * t265) * t244;
t142 = -t289 * t363 + (-mrSges(4,2) * t347 + mrSges(4,3) * t266) * t244;
t140 = Ifges(4,1) * t198 - Ifges(4,4) * t264 - Ifges(4,5) * t354;
t139 = Ifges(4,4) * t198 - Ifges(4,2) * t264 - Ifges(4,6) * t354;
t134 = pkin(3) * t261 + (t234 + (t337 + pkin(8)) * t354) * qJD(2);
t131 = mrSges(5,2) * t354 - mrSges(5,3) * t146;
t115 = t269 * t322 + (mrSges(4,1) * t198 - mrSges(4,2) * t264) * qJD(3);
t112 = -Ifges(5,4) * t165 - Ifges(5,2) * t166;
t111 = Ifges(6,2) * t165 + Ifges(6,6) * t166;
t110 = Ifges(6,6) * t165 + Ifges(6,3) * t166;
t109 = mrSges(5,1) * t166 - mrSges(5,2) * t165;
t108 = -mrSges(6,2) * t166 + mrSges(6,3) * t165;
t102 = (Ifges(4,1) * t290 - Ifges(4,4) * t289) * qJD(3) + (Ifges(4,1) * t265 + Ifges(4,4) * t266 + Ifges(4,5) * t347) * t244;
t101 = (Ifges(4,4) * t290 - Ifges(4,2) * t289) * qJD(3) + (Ifges(4,4) * t265 + Ifges(4,2) * t266 + Ifges(4,6) * t347) * t244;
t96 = -mrSges(6,2) * t146 - mrSges(6,3) * t147;
t95 = mrSges(5,1) * t146 + mrSges(5,2) * t147;
t94 = pkin(4) * t166 + t274;
t90 = Ifges(5,4) * t147 - Ifges(5,2) * t146 - Ifges(5,6) * t354;
t89 = -Ifges(6,4) * t354 - Ifges(6,2) * t147 + Ifges(6,6) * t146;
t88 = -Ifges(6,5) * t354 - Ifges(6,6) * t147 + Ifges(6,3) * t146;
t77 = mrSges(7,1) * t275 + mrSges(7,2) * t276;
t69 = -mrSges(7,1) * t124 - mrSges(7,2) * t278;
t67 = -mrSges(5,2) * t323 - mrSges(5,3) * t76;
t66 = mrSges(5,1) * t323 + mrSges(5,3) * t75;
t62 = t146 * pkin(4) + t251;
t32 = mrSges(5,1) * t76 - mrSges(5,2) * t75;
t31 = -mrSges(6,2) * t76 + mrSges(6,3) * t75;
t29 = -Ifges(5,4) * t75 - Ifges(5,2) * t76 + Ifges(5,6) * t323;
t28 = Ifges(6,4) * t323 + Ifges(6,2) * t75 + Ifges(6,6) * t76;
t27 = Ifges(6,5) * t323 + Ifges(6,6) * t75 + Ifges(6,3) * t76;
t21 = -mrSges(7,1) * t40 + mrSges(7,2) * t41;
t3 = [-0.2e1 * (t203 * t249 + t204 * t247) * mrSges(3,3) * t348 - t264 * t101 + 0.2e1 * t197 * (mrSges(4,1) * t264 + t198 * mrSges(4,2)) - 0.2e1 * t273 * t131 + (t134 * t148 - t273 * t53 - t277 * t52) * t403 - 0.2e1 * t277 * t132 - t278 * t10 + (t1 * t20 + t19 * t2 + t35 * t5) * t401 + (t11 * t45 + t14 * t46 + t24 * t62) * t402 - t256 * t139 + t257 * t140 + (-t28 + t427) * t147 + (t89 - t421) * t75 - 0.2e1 * t364 * t420 + (0.2e1 * mrSges(3,3) * t196 - t325 - t413) * t354 + (Ifges(3,5) * t364 + (Ifges(3,1) * t247 + t376) * t244) * t322 - 0.2e1 * t197 * (mrSges(3,1) * t364 - mrSges(3,3) * t355) + t364 * t280 + (t27 - t29) * t146 + ((-Ifges(3,2) * t247 + t376) * t346 + (Ifges(3,1) * t249 - t377) * t347 - 0.2e1 * pkin(1) * (mrSges(3,1) * t247 + mrSges(3,2) * t249) * qJD(2)) * t244 ^ 2 + 0.2e1 * t20 * t25 + 0.2e1 * t19 * t26 + (-Ifges(3,6) * t364 - (Ifges(3,2) * t249 + t377) * t244 + Ifges(4,5) * t198 - Ifges(4,6) * t264 + t423 * t147 + t422 * t146 + (-Ifges(4,3) - t434) * t354) * t323 + 0.2e1 * t35 * t21 + t40 * t48 + t41 * t49 + 0.2e1 * t62 * t31 + 0.2e1 * t45 * t64 + 0.2e1 * t46 * t65 + 0.2e1 * t52 * t66 + 0.2e1 * t53 * t67 + 0.2e1 * t5 * t69 + 0.2e1 * t1 * t82 + 0.2e1 * t2 * t83 + 0.2e1 * t24 * t96 + t124 * t9 + 0.2e1 * t11 * t129 + 0.2e1 * t14 * t130 + 0.2e1 * t134 * t95 + (t88 - t90) * t76 + 0.2e1 * t136 * t142 + 0.2e1 * t135 * t143 + 0.2e1 * t148 * t32 + 0.2e1 * t84 * t178 + 0.2e1 * t85 * t179 + 0.2e1 * t190 * t115 + t198 * t102 + 0.2e1 * m(4) * (t135 * t85 + t136 * t84 + t190 * t197) + 0.2e1 * m(3) * (t196 * t204 - t197 * t203); (-t421 / 0.2e1 + t89 / 0.2e1) * t165 - t420 + (t418 / 0.2e1 - t111 / 0.2e1) * t147 - t264 * t213 / 0.2e1 + t95 * t241 + (t110 / 0.2e1 - t112 / 0.2e1) * t146 + t280 + t431 * t49 + t432 * t48 + (-m(5) * t273 + t408 + t67) * t181 + t273 * t367 + (t427 / 0.2e1 + t277 * mrSges(5,3) + t14 * mrSges(6,1) - t28 / 0.2e1 + t423 * t293) * t207 + (m(5) * t277 + t407 - t66) * t180 + t56 * t392 + t55 * t393 + t128 * t396 + t127 * t397 + t11 * t368 + t45 * t370 + t52 * t371 - t256 * t219 / 0.2e1 + t257 * t221 / 0.2e1 + m(5) * (t239 * t134 + t148 * t241) + t9 * t316 + t10 * t317 + (-m(4) * pkin(2) - mrSges(4,1) * t384 + mrSges(4,2) * t382 - mrSges(3,1)) * t197 + (t88 / 0.2e1 - t90 / 0.2e1) * t166 + m(4) * ((-t135 * t384 - t136 * t382) * qJD(3) + t414) * pkin(9) + (t27 / 0.2e1 - t29 / 0.2e1 + t422 * t293) * t206 - (t279 + t412) * t354 / 0.2e1 + (-t135 * t312 - t136 * t310 + t414) * mrSges(4,3) + (-m(5) * t53 - t131 + t406) * t118 + (-m(5) * t52 + t405) * t119 + (Ifges(4,5) * t382 + Ifges(4,6) * t384) * t293 + t384 * t101 / 0.2e1 + t382 * t102 / 0.2e1 - t46 * t372 - t53 * t369 - t143 * t336 + m(7) * (t1 * t81 + t145 * t5 + t19 * t23 + t2 * t80 + t20 * t22 + t35 * t86) + t140 * t312 / 0.2e1 - t139 * t310 / 0.2e1 - t179 * t300 - t178 * t299 + m(6) * (t150 * t24 + t62 * t94) + (t171 / 0.2e1 - t173 / 0.2e1) * t76 + t142 * t339 + t35 * t77 + t80 * t26 + t81 * t25 + t22 * t82 + t23 * t83 + t86 * t69 + t20 * t92 + t19 * t93 + t94 * t96 + t62 * t108 - pkin(2) * t115 + t145 * t21 + t148 * t109 - t5 * t149 + t150 * t31 + t2 * t156 + t1 * t157 + t134 * t169 + t24 * t170 + t190 * t210 + t198 * t215 / 0.2e1 + t239 * t32 + (-t417 / 0.2e1 + t172 / 0.2e1) * t75; t382 * t215 + t384 * t213 - 0.2e1 * pkin(2) * t210 + 0.2e1 * t150 * t108 + 0.2e1 * t239 * t109 + 0.2e1 * t145 * t77 - 0.2e1 * t86 * t149 + 0.2e1 * t23 * t156 + 0.2e1 * t22 * t157 + 0.2e1 * t94 * t170 + 0.2e1 * t80 * t93 + 0.2e1 * t81 * t92 + (0.2e1 * t169 * t337 - t219 * t382 + t221 * t384) * qJD(3) + (t239 * t241 + t283) * t403 + (t145 * t86 + t22 * t81 + t23 * t80) * t401 + (t150 * t94 + t283) * t402 + (t119 * t428 - t111 + t418) * t207 + (-t180 * t428 + t172 - t417) * t165 + (t127 * t248 + t128 * t245 - t181 * t428 + t171 - t173) * t166 + (t245 * t56 + t248 * t55 + t110 - t112 + t118 * t428 + (-t127 * t245 + t128 * t248) * qJD(6)) * t206; t131 * t298 + t83 * t301 + t82 * t302 + t325 + (-t383 * t277 - t273 * t246 + (-t246 * t52 + t383 * t53) * qJD(4)) * t398 + t67 * t381 + t26 * t356 + t25 * t357 - t20 * t327 + t252 + t66 * t338 - t84 * mrSges(4,2) + t85 * mrSges(4,1) + t407 * t238 + (m(7) * (t19 * t248 + t20 * t245) + t405) * t335 + (m(7) * t5 + t21 + t408) * t236 + (t415 + t425) * t235 + (m(7) * t35 - t406 + t69) * t228 + t410 * mrSges(7,3); mrSges(4,2) * t299 + t279 + t156 * t301 + t157 * t302 + (-t383 * t119 + (t180 * t246 + t181 * t383) * qJD(4)) * t398 + t93 * t356 + t92 * t357 + t338 * t371 + t253 - t369 * t381 - t298 * t367 - t81 * t327 - mrSges(4,1) * t300 + (m(7) * t86 - t370 + t77) * t236 + (t416 + t426) * t235 + (m(6) * t181 + m(7) * t145 - t149 - t368) * t228 + (-m(6) * t236 - t246 * t398 + t424) * t118 + t411 * mrSges(7,3) + t430 * t238 + (m(6) * t180 + m(7) * (t245 * t81 + t248 * t80) + t378 * t207) * t335; t236 * t399 + 0.2e1 * t419 * t228 + 0.2e1 * t429 + (t235 * t281 + t358) * t401 + (t238 * t335 + t358) * t402 + t263; -(t255 + t415) * t394 + (-t379 + (-t2 - t362) * t248) * mrSges(7,3) + (t21 - t64) * qJ(5) + (t69 - t129) * qJD(5) + m(7) * (qJ(5) * t5 + qJD(5) * t35) + m(6) * (-pkin(4) * t14 - qJ(5) * t11 - qJD(5) * t45) + t252 - pkin(4) * t65; m(7) * (qJ(5) * t86 + qJD(5) * t145) + m(6) * (-pkin(4) * t119 - qJ(5) * t118 + qJD(5) * t181) - (t254 + t416) * t394 + (pkin(4) * t165 - qJ(5) * t166 - qJD(5) * t206) * mrSges(6,1) + t253 + t424 * t118 + (-t365 + (-t23 - t361) * t248) * mrSges(7,3) + qJ(5) * t77 - qJD(5) * t149; (qJ(5) + t236) * t209 + t429 + m(7) * (-t281 * t394 + t282) + m(6) * (-pkin(4) * t335 + t282) + t263 + t419 * (qJD(5) + t228); qJ(5) * t399 + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t419) * qJD(5) + t263; (-t245 * t83 + t248 * t82) * qJD(6) + t255 + t407; (-t156 * t245 + t157 * t248) * qJD(6) + t254 + t430; (m(7) * t349 + m(6)) * t335; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t8; mrSges(7,1) * t23 - mrSges(7,2) * t22 + t54; t287 * t335 + (-t216 * t235 - t284) * qJD(6); ((mrSges(7,2) * t394 - Ifges(7,6)) * t248 + (mrSges(7,1) * t394 - Ifges(7,5)) * t245) * qJD(6); -t216 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
