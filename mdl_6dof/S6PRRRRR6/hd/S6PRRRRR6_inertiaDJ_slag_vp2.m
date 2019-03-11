% Calculate time derivative of joint inertia matrix for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_inertiaDJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:14:15
% EndTime: 2019-03-09 01:14:38
% DurationCPUTime: 9.65s
% Computational Cost: add. (16343->848), mult. (50877->1258), div. (0->0), fcn. (52267->16), ass. (0->349)
t267 = sin(pkin(7));
t269 = cos(pkin(8));
t270 = cos(pkin(7));
t275 = sin(qJ(3));
t279 = cos(qJ(4));
t347 = t275 * t279;
t274 = sin(qJ(4));
t280 = cos(qJ(3));
t349 = t274 * t280;
t291 = t269 * t349 + t347;
t266 = sin(pkin(8));
t339 = qJD(4) * t274;
t319 = t266 * t339;
t122 = t270 * t319 + (t291 * qJD(4) + (t269 * t347 + t349) * qJD(3)) * t267;
t272 = sin(qJ(6));
t277 = cos(qJ(6));
t338 = qJD(4) * t279;
t318 = t266 * t338;
t344 = t279 * t280;
t350 = t274 * t275;
t414 = t269 * t344 - t350;
t121 = t270 * t318 + (t414 * qJD(4) + (-t269 * t350 + t344) * qJD(3)) * t267;
t357 = t266 * t274;
t162 = t267 * t291 + t270 * t357;
t353 = t267 * t280;
t212 = -t266 * t353 + t269 * t270;
t273 = sin(qJ(5));
t278 = cos(qJ(5));
t125 = t162 * t273 - t278 * t212;
t341 = qJD(3) * t267;
t322 = t275 * t341;
t310 = t266 * t322;
t69 = -qJD(5) * t125 + t121 * t278 + t273 * t310;
t126 = t162 * t278 + t212 * t273;
t356 = t266 * t279;
t159 = -t414 * t267 - t270 * t356;
t89 = -t126 * t272 + t159 * t277;
t27 = qJD(6) * t89 + t122 * t272 + t277 * t69;
t90 = t126 * t277 + t159 * t272;
t28 = -qJD(6) * t90 + t122 * t277 - t272 * t69;
t12 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t379 = pkin(2) * t270;
t261 = t280 * t379;
t308 = t267 * (-pkin(11) * t269 - pkin(10));
t296 = t275 * t308;
t168 = pkin(3) * t270 + t261 + t296;
t377 = pkin(11) * t266;
t187 = (-pkin(3) * t280 - t275 * t377 - pkin(2)) * t267;
t115 = -t168 * t266 + t269 * t187;
t73 = pkin(4) * t159 - pkin(12) * t162 + t115;
t260 = t275 * t379;
t221 = pkin(10) * t353 + t260;
t156 = (t266 * t270 + t269 * t353) * pkin(11) + t221;
t352 = t269 * t274;
t84 = t279 * t156 + t168 * t352 + t187 * t357;
t77 = pkin(12) * t212 + t84;
t373 = t273 * t73 + t278 * t77;
t412 = qJD(5) * t373;
t170 = (t280 * t308 - t260) * qJD(3);
t190 = (pkin(3) * t275 - t280 * t377) * t341;
t254 = qJD(3) * t261;
t358 = t168 * t269;
t422 = qJD(3) * t296 + qJD(4) * t358 + t254;
t45 = -t156 * t339 + t170 * t352 + t187 * t318 + t190 * t357 + t279 * t422;
t43 = pkin(12) * t310 + t45;
t123 = -t170 * t266 + t269 * t190;
t54 = pkin(4) * t122 - pkin(12) * t121 + t123;
t11 = -t273 * t43 + t278 * t54 - t412;
t9 = -pkin(5) * t122 - t11;
t328 = -m(7) * t9 - t12;
t47 = mrSges(6,1) * t122 - mrSges(6,3) * t69;
t424 = -t328 - t47;
t215 = -t278 * t269 + t273 * t357;
t171 = -qJD(5) * t215 + t278 * t318;
t141 = mrSges(6,1) * t319 - mrSges(6,3) * t171;
t216 = t269 * t273 + t278 * t357;
t173 = -t216 * t272 - t277 * t356;
t104 = qJD(6) * t173 + t171 * t277 + t272 * t319;
t292 = -t216 * t277 + t272 * t356;
t105 = qJD(6) * t292 - t171 * t272 + t277 * t319;
t57 = -mrSges(7,1) * t105 + mrSges(7,2) * t104;
t340 = qJD(4) * t266;
t207 = (pkin(4) * t274 - pkin(12) * t279) * t340;
t256 = pkin(11) * t357;
t351 = t269 * t279;
t218 = pkin(3) * t351 - t256;
t208 = t218 * qJD(4);
t220 = pkin(3) * t352 + pkin(11) * t356;
t197 = pkin(12) * t269 + t220;
t198 = (-pkin(4) * t279 - pkin(12) * t274 - pkin(3)) * t266;
t416 = t278 * t197 + t273 * t198;
t420 = qJD(5) * t416;
t92 = t207 * t278 - t208 * t273 - t420;
t86 = -pkin(5) * t319 - t92;
t324 = -m(7) * t86 - t57;
t423 = -t141 - t324;
t421 = -m(6) * pkin(4) - mrSges(6,1) * t278 + mrSges(6,2) * t273;
t334 = qJD(6) * t273;
t336 = qJD(5) * t278;
t287 = -t272 * t334 + t277 * t336;
t382 = t272 / 0.2e1;
t380 = t277 / 0.2e1;
t49 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t94 = mrSges(6,1) * t159 - mrSges(6,3) * t126;
t419 = -t94 + t49;
t131 = mrSges(5,1) * t212 - mrSges(5,3) * t162;
t80 = mrSges(6,1) * t125 + mrSges(6,2) * t126;
t418 = -t131 + t80;
t114 = -mrSges(7,1) * t173 - mrSges(7,2) * t292;
t179 = -mrSges(6,1) * t356 - mrSges(6,3) * t216;
t417 = -t179 + t114;
t150 = mrSges(6,1) * t215 + mrSges(6,2) * t216;
t224 = mrSges(5,1) * t269 - mrSges(5,3) * t357;
t415 = -t224 + t150;
t413 = -mrSges(5,1) + t421;
t83 = -t274 * t156 + t279 * (t187 * t266 + t358);
t242 = -mrSges(7,1) * t277 + mrSges(7,2) * t272;
t411 = -m(7) * pkin(5) - mrSges(6,1) + t242;
t410 = 2 * m(5);
t409 = 0.2e1 * m(6);
t408 = 0.2e1 * m(7);
t407 = 0.2e1 * pkin(12);
t406 = -2 * mrSges(4,3);
t405 = -0.2e1 * mrSges(5,3);
t404 = t267 ^ 2;
t70 = qJD(5) * t126 + t121 * t273 - t278 * t310;
t22 = Ifges(6,1) * t69 - Ifges(6,4) * t70 + Ifges(6,5) * t122;
t402 = t22 / 0.2e1;
t62 = Ifges(6,1) * t126 - Ifges(6,4) * t125 + Ifges(6,5) * t159;
t401 = t62 / 0.2e1;
t172 = qJD(5) * t216 + t273 * t318;
t110 = Ifges(6,1) * t171 - Ifges(6,4) * t172 + Ifges(6,5) * t319;
t400 = t110 / 0.2e1;
t140 = Ifges(6,1) * t216 - Ifges(6,4) * t215 - Ifges(6,5) * t356;
t399 = t140 / 0.2e1;
t366 = Ifges(7,4) * t272;
t246 = Ifges(7,2) * t277 + t366;
t365 = Ifges(7,4) * t277;
t304 = -Ifges(7,2) * t272 + t365;
t144 = -t246 * t334 + (Ifges(7,6) * t273 + t278 * t304) * qJD(5);
t398 = t144 / 0.2e1;
t248 = Ifges(7,1) * t272 + t365;
t305 = Ifges(7,1) * t277 - t366;
t145 = -t248 * t334 + (Ifges(7,5) * t273 + t278 * t305) * qJD(5);
t397 = t145 / 0.2e1;
t396 = t173 / 0.2e1;
t395 = -t292 / 0.2e1;
t200 = -Ifges(7,6) * t278 + t273 * t304;
t394 = t200 / 0.2e1;
t201 = -Ifges(7,5) * t278 + t273 * t305;
t393 = t201 / 0.2e1;
t333 = qJD(6) * t277;
t264 = Ifges(7,5) * t333;
t335 = qJD(6) * t272;
t230 = -Ifges(7,6) * t335 + t264;
t392 = t230 / 0.2e1;
t232 = t304 * qJD(6);
t391 = t232 / 0.2e1;
t234 = t305 * qJD(6);
t390 = t234 / 0.2e1;
t368 = Ifges(6,4) * t273;
t235 = (Ifges(6,1) * t278 - t368) * qJD(5);
t389 = t235 / 0.2e1;
t388 = Ifges(7,5) * t382 + Ifges(7,6) * t380;
t387 = Ifges(6,5) * t273 / 0.2e1 + Ifges(6,6) * t278 / 0.2e1;
t386 = t246 / 0.2e1;
t385 = t248 / 0.2e1;
t367 = Ifges(6,4) * t278;
t249 = Ifges(6,1) * t273 + t367;
t384 = t249 / 0.2e1;
t383 = -t272 / 0.2e1;
t381 = -t277 / 0.2e1;
t378 = pkin(3) * t266;
t376 = pkin(12) * t278;
t268 = sin(pkin(6));
t276 = sin(qJ(2));
t345 = t276 * t280;
t281 = cos(qJ(2));
t346 = t275 * t281;
t288 = t270 * t346 + t345;
t285 = t288 * qJD(3);
t289 = -t270 * t345 - t346;
t271 = cos(pkin(6));
t309 = t271 * t322;
t354 = t267 * t276;
t282 = t266 * t309 + (t266 * t285 + (-t266 * t289 + t269 * t354) * qJD(2)) * t268;
t343 = t280 * t281;
t348 = t275 * t276;
t290 = t270 * t343 - t348;
t321 = t280 * t341;
t120 = t271 * t321 + (t290 * qJD(3) + (-t270 * t348 + t343) * qJD(2)) * t268;
t160 = t268 * t290 + t271 * t353;
t213 = -t267 * t268 * t281 + t271 * t270;
t300 = t160 * t269 + t213 * t266;
t355 = t267 * t275;
t161 = t268 * t288 + t271 * t355;
t359 = t161 * t274;
t40 = -t309 * t352 + t120 * t279 + (-t269 * t285 + (t266 * t354 + t269 * t289) * qJD(2)) * t274 * t268 + (t279 * t300 - t359) * qJD(4);
t124 = -t160 * t266 + t213 * t269;
t88 = t161 * t279 + t274 * t300;
t56 = t124 * t273 + t278 * t88;
t15 = qJD(5) * t56 + t273 * t40 - t278 * t282;
t55 = -t124 * t278 + t273 * t88;
t375 = t15 * t55;
t283 = -t309 + (qJD(2) * t289 - t285) * t268;
t342 = qJD(2) * t268;
t323 = t276 * t342;
t39 = -t267 * t323 * t356 + qJD(4) * t88 + t120 * t274 - t283 * t351;
t87 = -t160 * t351 - t213 * t356 + t359;
t374 = t87 * t39;
t372 = m(7) * qJD(5);
t371 = mrSges(7,3) * t273;
t370 = Ifges(5,4) * t274;
t369 = Ifges(5,4) * t279;
t364 = Ifges(7,6) * t272;
t363 = t15 * t273;
t337 = qJD(5) * t273;
t16 = t124 * t336 + t273 * t282 + t40 * t278 - t337 * t88;
t362 = t16 * t278;
t209 = t220 * qJD(4);
t361 = t209 * t87;
t332 = qJD(6) * t278;
t5 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t70;
t21 = Ifges(6,4) * t69 - Ifges(6,2) * t70 + Ifges(6,6) * t122;
t330 = t5 / 0.2e1 - t21 / 0.2e1;
t20 = Ifges(6,5) * t69 - Ifges(6,6) * t70 + Ifges(6,3) * t122;
t36 = Ifges(7,5) * t90 + Ifges(7,6) * t89 + Ifges(7,3) * t125;
t61 = Ifges(6,4) * t126 - Ifges(6,2) * t125 + Ifges(6,6) * t159;
t329 = -t61 / 0.2e1 + t36 / 0.2e1;
t139 = Ifges(6,4) * t216 - Ifges(6,2) * t215 - Ifges(6,6) * t356;
t98 = -Ifges(7,5) * t292 + Ifges(7,6) * t173 + Ifges(7,3) * t215;
t326 = -t139 / 0.2e1 + t98 / 0.2e1;
t109 = Ifges(6,4) * t171 - Ifges(6,2) * t172 + Ifges(6,6) * t319;
t50 = Ifges(7,5) * t104 + Ifges(7,6) * t105 + Ifges(7,3) * t172;
t325 = t50 / 0.2e1 - t109 / 0.2e1;
t66 = Ifges(5,5) * t121 - Ifges(5,6) * t122 + Ifges(5,3) * t310;
t108 = Ifges(6,5) * t171 - Ifges(6,6) * t172 + Ifges(6,3) * t319;
t286 = t272 * t336 + t273 * t333;
t143 = Ifges(7,5) * t287 - Ifges(7,6) * t286 + Ifges(7,3) * t337;
t233 = (-Ifges(6,2) * t273 + t367) * qJD(5);
t315 = -t233 / 0.2e1 + t143 / 0.2e1;
t199 = -Ifges(7,3) * t278 + (Ifges(7,5) * t277 - t364) * t273;
t247 = Ifges(6,2) * t278 + t368;
t314 = -t247 / 0.2e1 + t199 / 0.2e1;
t238 = (pkin(5) * t273 - pkin(13) * t278) * qJD(5);
t241 = -pkin(5) * t278 - pkin(13) * t273 - pkin(4);
t134 = t241 * t333 + t238 * t272 + (-t272 * t332 - t277 * t337) * pkin(12);
t191 = t241 * t277 - t272 * t376;
t313 = -qJD(6) * t191 + t134;
t135 = -t241 * t335 + t238 * t277 + (t272 * t337 - t277 * t332) * pkin(12);
t192 = t241 * t272 + t277 * t376;
t312 = -qJD(6) * t192 - t135;
t311 = -t156 * t338 - t187 * t319 - t274 * t422;
t30 = pkin(13) * t159 + t373;
t76 = -pkin(4) * t212 - t83;
t41 = pkin(5) * t125 - pkin(13) * t126 + t76;
t13 = -t272 * t30 + t277 * t41;
t44 = -t170 * t351 + (-pkin(4) * t322 - t190 * t279) * t266 - t311;
t17 = pkin(5) * t70 - pkin(13) * t69 + t44;
t10 = t273 * t54 + t278 * t43 + t73 * t336 - t337 * t77;
t8 = pkin(13) * t122 + t10;
t1 = qJD(6) * t13 + t17 * t272 + t277 * t8;
t14 = t272 * t41 + t277 * t30;
t2 = -qJD(6) * t14 + t17 * t277 - t272 * t8;
t307 = t1 * t277 - t2 * t272;
t306 = mrSges(7,1) * t272 + mrSges(7,2) * t277;
t101 = pkin(5) * t172 - pkin(13) * t171 + t209;
t196 = t256 + (-pkin(3) * t279 - pkin(4)) * t269;
t127 = pkin(5) * t215 - pkin(13) * t216 + t196;
t129 = -pkin(13) * t356 + t416;
t78 = t127 * t277 - t129 * t272;
t91 = -t197 * t337 + t198 * t336 + t273 * t207 + t278 * t208;
t85 = pkin(13) * t319 + t91;
t23 = qJD(6) * t78 + t101 * t272 + t277 * t85;
t79 = t127 * t272 + t129 * t277;
t24 = -qJD(6) * t79 + t101 * t277 - t272 * t85;
t303 = t23 * t277 - t24 * t272;
t33 = t272 * t87 + t277 * t56;
t32 = -t272 * t56 + t277 * t87;
t34 = -t273 * t77 + t278 * t73;
t136 = -t197 * t273 + t198 * t278;
t37 = Ifges(7,4) * t90 + Ifges(7,2) * t89 + Ifges(7,6) * t125;
t38 = Ifges(7,1) * t90 + Ifges(7,4) * t89 + Ifges(7,5) * t125;
t295 = t37 * t383 + t38 * t380;
t294 = t336 * t55 + t363;
t100 = -Ifges(7,1) * t292 + Ifges(7,4) * t173 + Ifges(7,5) * t215;
t99 = -Ifges(7,4) * t292 + Ifges(7,2) * t173 + Ifges(7,6) * t215;
t293 = t100 * t380 + t383 * t99;
t265 = Ifges(6,5) * t336;
t253 = Ifges(4,5) * t321;
t252 = Ifges(5,5) * t318;
t237 = -mrSges(7,1) * t278 - t277 * t371;
t236 = mrSges(7,2) * t278 - t272 * t371;
t231 = -Ifges(6,6) * t337 + t265;
t229 = (mrSges(6,1) * t273 + mrSges(6,2) * t278) * qJD(5);
t228 = t306 * qJD(6);
t227 = -mrSges(4,2) * t270 + mrSges(4,3) * t353;
t226 = mrSges(4,1) * t270 - mrSges(4,3) * t355;
t225 = -mrSges(5,2) * t269 + mrSges(5,3) * t356;
t222 = t306 * t273;
t219 = -pkin(10) * t355 + t261;
t217 = (-mrSges(5,1) * t279 + mrSges(5,2) * t274) * t266;
t211 = t221 * qJD(3);
t210 = -pkin(10) * t322 + t254;
t206 = (Ifges(5,1) * t279 - t370) * t340;
t205 = (-Ifges(5,2) * t274 + t369) * t340;
t204 = -Ifges(5,6) * t319 + t252;
t203 = (mrSges(4,1) * t275 + mrSges(4,2) * t280) * t341;
t202 = (mrSges(5,1) * t274 + mrSges(5,2) * t279) * t340;
t194 = Ifges(5,5) * t269 + (Ifges(5,1) * t274 + t369) * t266;
t193 = Ifges(5,6) * t269 + (Ifges(5,2) * t279 + t370) * t266;
t183 = -mrSges(7,2) * t337 - mrSges(7,3) * t286;
t182 = mrSges(7,1) * t337 - mrSges(7,3) * t287;
t178 = mrSges(6,2) * t356 - mrSges(6,3) * t215;
t155 = mrSges(7,1) * t286 + mrSges(7,2) * t287;
t142 = -mrSges(6,2) * t319 - mrSges(6,3) * t172;
t138 = Ifges(6,5) * t216 - Ifges(6,6) * t215 - Ifges(6,3) * t356;
t133 = mrSges(7,1) * t215 + mrSges(7,3) * t292;
t132 = -mrSges(7,2) * t215 + mrSges(7,3) * t173;
t130 = -mrSges(5,2) * t212 - mrSges(5,3) * t159;
t128 = pkin(5) * t356 - t136;
t113 = mrSges(6,1) * t172 + mrSges(6,2) * t171;
t111 = mrSges(5,1) * t159 + mrSges(5,2) * t162;
t107 = -mrSges(5,2) * t310 - mrSges(5,3) * t122;
t106 = mrSges(5,1) * t310 - mrSges(5,3) * t121;
t96 = Ifges(5,1) * t162 - Ifges(5,4) * t159 + Ifges(5,5) * t212;
t95 = Ifges(5,4) * t162 - Ifges(5,2) * t159 + Ifges(5,6) * t212;
t93 = -mrSges(6,2) * t159 - mrSges(6,3) * t125;
t82 = -mrSges(7,2) * t172 + mrSges(7,3) * t105;
t81 = mrSges(7,1) * t172 - mrSges(7,3) * t104;
t75 = mrSges(5,1) * t122 + mrSges(5,2) * t121;
t68 = Ifges(5,1) * t121 - Ifges(5,4) * t122 + Ifges(5,5) * t310;
t67 = Ifges(5,4) * t121 - Ifges(5,2) * t122 + Ifges(5,6) * t310;
t60 = Ifges(6,5) * t126 - Ifges(6,6) * t125 + Ifges(6,3) * t159;
t59 = mrSges(7,1) * t125 - mrSges(7,3) * t90;
t58 = -mrSges(7,2) * t125 + mrSges(7,3) * t89;
t52 = Ifges(7,1) * t104 + Ifges(7,4) * t105 + Ifges(7,5) * t172;
t51 = Ifges(7,4) * t104 + Ifges(7,2) * t105 + Ifges(7,6) * t172;
t48 = -mrSges(6,2) * t122 - mrSges(6,3) * t70;
t46 = (t170 * t269 + t190 * t266) * t279 + t311;
t31 = mrSges(6,1) * t70 + mrSges(6,2) * t69;
t29 = -pkin(5) * t159 - t34;
t19 = -mrSges(7,2) * t70 + mrSges(7,3) * t28;
t18 = mrSges(7,1) * t70 - mrSges(7,3) * t27;
t7 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t70;
t6 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t70;
t4 = qJD(6) * t32 + t16 * t277 + t272 * t39;
t3 = -qJD(6) * t33 - t16 * t272 + t277 * t39;
t25 = [0.2e1 * m(7) * (t3 * t32 + t33 * t4 + t375) + 0.2e1 * m(6) * (t16 * t56 + t374 + t375) + 0.2e1 * m(5) * (t124 * t282 + t88 * t40 + t374) + 0.2e1 * m(4) * (-t160 * t309 + t161 * t120 + (-t160 * t285 + (t160 * t289 + t213 * t354) * qJD(2)) * t268); -t281 * mrSges(3,2) * t342 + t283 * t226 + m(5) * (t115 * t282 + t123 * t124 + t84 * t40 + t45 * t88) + t282 * t111 + m(6) * (t10 * t56 + t16 * t373) + m(7) * (t1 * t33 + t13 * t3 + t14 * t4 + t2 * t32) + t120 * t227 + t213 * t203 + t40 * t130 + t124 * t75 + t88 * t107 + t16 * t93 + t3 * t59 + t56 * t48 + t4 * t58 + m(4) * (-t219 * t309 + t221 * t120 - t211 * t160 + t210 * t161 + (-t219 * t285 + (-pkin(2) * t276 * t404 + t219 * t289) * qJD(2)) * t268) + t32 * t18 + t33 * t19 + (-m(5) * t46 + m(6) * t44 - t106 + t31) * t87 + (-m(6) * t11 + t424) * t55 + (-m(5) * t83 + m(6) * t76 + t418) * t39 + (t404 * (-mrSges(4,1) * t280 + mrSges(4,2) * t275) - mrSges(3,1)) * t323 + (-m(6) * t34 + m(7) * t29 + t419) * t15 + (-t160 * t321 - t161 * t322) * mrSges(4,3); 0.2e1 * t373 * t48 + (t10 * t373 + t11 * t34 + t44 * t76) * t409 + 0.2e1 * m(4) * (t210 * t221 - t211 * t219) + (-0.2e1 * pkin(2) * t203 + ((0.2e1 * Ifges(4,4) * t353 + Ifges(4,5) * t270 + t219 * t406) * t280 + (-0.2e1 * Ifges(4,6) * t270 + t266 * (Ifges(5,5) * t162 - Ifges(5,6) * t159 + Ifges(5,3) * t212) + t221 * t406 + 0.2e1 * (-Ifges(4,4) * t275 + (Ifges(4,1) - Ifges(4,2)) * t280) * t267) * t275) * qJD(3)) * t267 + (t20 - t67) * t159 + (t36 - t61) * t70 - 0.2e1 * t211 * t226 + 0.2e1 * t210 * t227 + t212 * t66 + t162 * t68 + 0.2e1 * t45 * t130 + 0.2e1 * t46 * t131 + t126 * t22 + t121 * t96 + 0.2e1 * t123 * t111 + 0.2e1 * t115 * t75 + (t5 - t21) * t125 + 0.2e1 * t83 * t106 + 0.2e1 * t84 * t107 + 0.2e1 * t10 * t93 + 0.2e1 * t11 * t94 + t89 * t6 + t90 * t7 + 0.2e1 * t44 * t80 + 0.2e1 * t76 * t31 + t69 * t62 + 0.2e1 * t2 * t59 + (t60 - t95) * t122 + 0.2e1 * t1 * t58 + 0.2e1 * t9 * t49 + 0.2e1 * t34 * t47 + t28 * t37 + t27 * t38 + (t1 * t14 + t13 * t2 + t29 * t9) * t408 + (t115 * t123 + t45 * t84 + t46 * t83) * t410 + t270 * t253 + 0.2e1 * t13 * t18 + 0.2e1 * t14 * t19 + 0.2e1 * t29 * t12; m(7) * (t23 * t33 + t24 * t32 + t3 * t78 + t4 * t79) + t16 * t178 + t56 * t142 + t3 * t133 + t32 * t81 + t4 * t132 + t33 * t82 + t87 * t113 + m(6) * (t16 * t416 + t56 * t91 + t361) + m(5) * (t208 * t88 + t220 * t40 - t282 * t378 + t361) + t282 * t217 + t124 * t202 + t40 * t225 + t283 * mrSges(4,1) - t120 * mrSges(4,2) + (-m(6) * t92 + t423) * t55 + (-m(5) * t218 + m(6) * t196 + t415) * t39 + (-m(6) * t136 + m(7) * t128 + t417) * t15 + (t318 * t87 - t319 * t88) * mrSges(5,3); t373 * t142 + ((-t83 * mrSges(5,3) + t96 / 0.2e1) * t279 + (-t84 * mrSges(5,3) - t95 / 0.2e1 + t60 / 0.2e1) * t274) * t340 + (t274 * t68 / 0.2e1 - pkin(3) * t75 + (t67 / 0.2e1 - t20 / 0.2e1) * t279) * t266 + t253 + (t266 * (Ifges(5,3) * t269 + (Ifges(5,5) * t274 + Ifges(5,6) * t279) * t266) / 0.2e1 - Ifges(4,6)) * t322 + t418 * t209 + m(7) * (t1 * t79 + t128 * t9 + t13 * t24 + t14 * t23 + t2 * t78 + t29 * t86) + m(6) * (t10 * t416 + t11 * t136 + t196 * t44 + t209 * t76 + t34 * t92 + t373 * t91) + t416 * t48 + (-t205 / 0.2e1 + t108 / 0.2e1) * t159 + (-t193 / 0.2e1 + t138 / 0.2e1) * t122 + t269 * t66 / 0.2e1 + t123 * t217 + t218 * t106 + t220 * t107 + t46 * t224 + t45 * t225 + t212 * t204 / 0.2e1 + t115 * t202 + t162 * t206 / 0.2e1 + t208 * t130 - t210 * mrSges(4,2) - t211 * mrSges(4,1) + t121 * t194 / 0.2e1 + t196 * t31 + t10 * t178 + t11 * t179 + t34 * t141 + t44 * t150 + t1 * t132 + t2 * t133 + t136 * t47 + t128 * t12 + t76 * t113 + t9 * t114 + t104 * t38 / 0.2e1 + t105 * t37 / 0.2e1 + t91 * t93 + t92 * t94 + t28 * t99 / 0.2e1 + t27 * t100 / 0.2e1 + t89 * t51 / 0.2e1 + t90 * t52 / 0.2e1 + t79 * t19 + t13 * t81 + t14 * t82 + t86 * t49 + t78 * t18 + t24 * t59 + t29 * t57 + t23 * t58 + t7 * t395 + t6 * t396 + t69 * t399 + t126 * t400 + t171 * t401 + t216 * t402 + t325 * t125 + t326 * t70 + t329 * t172 + t330 * t215 + m(5) * (-t123 * t378 + t208 * t84 - t209 * t83 + t218 * t46 + t220 * t45); (t50 - t109) * t215 + (-t139 + t98) * t172 + 0.2e1 * t415 * t209 + 0.2e1 * t416 * t142 + (t136 * t92 + t196 * t209 + t416 * t91) * t409 - t292 * t52 + t269 * t204 + 0.2e1 * t208 * t225 + t216 * t110 + 0.2e1 * t196 * t113 + t173 * t51 + 0.2e1 * t91 * t178 + 0.2e1 * t92 * t179 + t171 * t140 + 0.2e1 * t136 * t141 + 0.2e1 * t23 * t132 + 0.2e1 * t24 * t133 + 0.2e1 * t128 * t57 + 0.2e1 * t86 * t114 + t104 * t100 + t105 * t99 + 0.2e1 * t78 * t81 + 0.2e1 * t79 * t82 + (t128 * t86 + t23 * t79 + t24 * t78) * t408 + (t208 * t220 - t209 * t218) * t410 + (-0.2e1 * pkin(3) * t202 + t206 * t274 + (-t108 + t205) * t279 + ((t218 * t405 + t194) * t279 + (t220 * t405 + t138 - t193) * t274) * qJD(4)) * t266; -t40 * mrSges(5,2) + t15 * t222 + t55 * t155 + t32 * t182 + t33 * t183 + t87 * t229 + t4 * t236 + t3 * t237 + m(7) * (t134 * t33 + t135 * t32 + t191 * t3 + t192 * t4) + (m(7) * t294 / 0.2e1 + m(6) * (-t337 * t56 + t294 + t362) / 0.2e1) * t407 + (t363 + t362 + (-t273 * t56 + t278 * t55) * qJD(5)) * mrSges(6,3) + t413 * t39; t66 + (t10 * mrSges(6,3) + (-t34 * mrSges(6,3) + t295 + t401) * qJD(5) + (t48 + t419 * qJD(5) + t29 * t372 + m(6) * (-qJD(5) * t34 + t10)) * pkin(12) - t330) * t278 + t1 * t236 + t2 * t237 + t9 * t222 + t76 * t229 + t159 * t231 / 0.2e1 + t14 * t183 + t191 * t18 + t192 * t19 + t13 * t182 + t29 * t155 + t134 * t58 + t135 * t59 - t45 * mrSges(5,2) + t46 * mrSges(5,1) + t314 * t70 + t315 * t125 + t69 * t384 + t122 * t387 + t126 * t389 + t27 * t393 + t28 * t394 + t90 * t397 + t89 * t398 + t421 * t44 + (t7 * t380 + t6 * t383 - t11 * mrSges(6,3) + t402 + (t37 * t381 + t38 * t383) * qJD(6) + (-mrSges(6,3) * t373 + t329) * qJD(5) + (-qJD(5) * t93 + m(6) * (-t11 - t412) + t424) * pkin(12)) * t273 + m(7) * (t1 * t192 + t13 * t135 + t134 * t14 + t191 * t2) - pkin(4) * t31; t252 + (t91 * mrSges(6,3) + (-t136 * mrSges(6,3) + t293 + t399) * qJD(5) + (t142 + t417 * qJD(5) + t128 * t372 + m(6) * (-qJD(5) * t136 + t91)) * pkin(12) - t325) * t278 + t413 * t209 + m(7) * (t134 * t79 + t135 * t78 + t191 * t24 + t192 * t23) + t23 * t236 + t24 * t237 + t86 * t222 + t196 * t229 - t208 * mrSges(5,2) + t79 * t183 + t191 * t81 + t192 * t82 + t78 * t182 + t128 * t155 + t134 * t132 + t135 * t133 - pkin(4) * t113 + t314 * t172 + t315 * t215 + t171 * t384 + t216 * t389 + t104 * t393 + t105 * t394 + t145 * t395 + t144 * t396 + (t52 * t380 + t51 * t383 - t92 * mrSges(6,3) + t400 + (t100 * t383 + t381 * t99) * qJD(6) + (-mrSges(6,3) * t416 + t326) * qJD(5) + (-qJD(5) * t178 + m(6) * (-t92 - t420) + t423) * pkin(12)) * t273 + (-t279 * t231 / 0.2e1 + (-Ifges(5,6) + t387) * t339) * t266; 0.2e1 * t135 * t237 + 0.2e1 * t191 * t182 + 0.2e1 * t134 * t236 + 0.2e1 * t192 * t183 + (t134 * t192 + t135 * t191) * t408 - 0.2e1 * pkin(4) * t229 + (-t143 + t233 + (-t200 * t272 + t201 * t277 + t222 * t407 + t249) * qJD(5)) * t278 + (t155 * t407 - t272 * t144 + t277 * t145 + t235 + (-t200 * t277 - t201 * t272) * qJD(6) + (pkin(12) ^ 2 * t278 * t408 + t199 - t247) * qJD(5)) * t273; -t16 * mrSges(6,2) + t55 * t228 + (m(7) * pkin(13) + mrSges(7,3)) * (-t3 * t272 + t4 * t277 + (-t272 * t33 - t277 * t32) * qJD(6)) + t411 * t15; t6 * t380 + t7 * t382 + t90 * t390 + t9 * t242 + t70 * t388 + t28 * t386 + t27 * t385 + t29 * t228 + t125 * t392 + t89 * t391 - t10 * mrSges(6,2) + t11 * mrSges(6,1) + t295 * qJD(6) + t328 * pkin(5) + ((-t13 * t277 - t14 * t272) * qJD(6) + t307) * mrSges(7,3) + (m(7) * (-t13 * t333 - t14 * t335 + t307) + t277 * t19 - t272 * t18 - t58 * t335 - t59 * t333) * pkin(13) + t20; t51 * t380 + t52 * t382 - t292 * t390 + t86 * t242 + t172 * t388 + t105 * t386 + t104 * t385 + t128 * t228 + t215 * t392 + t173 * t391 + t92 * mrSges(6,1) - t91 * mrSges(6,2) + t293 * qJD(6) + t324 * pkin(5) + ((-t272 * t79 - t277 * t78) * qJD(6) + t303) * mrSges(7,3) + (m(7) * (-t333 * t78 - t335 * t79 + t303) + t277 * t82 - t272 * t81 - t132 * t335 - t133 * t333) * pkin(13) + t108; -pkin(5) * t155 + t265 + (-t230 / 0.2e1 + t411 * qJD(5) * pkin(12)) * t278 + (qJD(6) * t393 + t398 + t336 * t385 + t313 * mrSges(7,3) + (m(7) * t313 - qJD(6) * t237 + t183) * pkin(13)) * t277 + (-qJD(6) * t200 / 0.2e1 + t397 - t246 * t336 / 0.2e1 + t312 * mrSges(7,3) + (m(7) * t312 - qJD(6) * t236 - t182) * pkin(13)) * t272 + (t234 * t380 + t232 * t383 + pkin(12) * t228 + (t246 * t381 + t248 * t383) * qJD(6) + (pkin(12) * mrSges(6,2) - Ifges(6,6) + t388) * qJD(5)) * t273; -0.2e1 * pkin(5) * t228 + t232 * t277 + t234 * t272 + (-t246 * t272 + t248 * t277) * qJD(6); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t24 - mrSges(7,2) * t23 + t50; mrSges(7,1) * t135 - mrSges(7,2) * t134 + t143; t264 + (pkin(13) * t242 - t364) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
