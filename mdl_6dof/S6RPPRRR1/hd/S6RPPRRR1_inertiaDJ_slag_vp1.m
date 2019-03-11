% Calculate time derivative of joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:57
% EndTime: 2019-03-09 02:18:11
% DurationCPUTime: 8.14s
% Computational Cost: add. (33759->632), mult. (24296->889), div. (0->0), fcn. (22731->12), ass. (0->331)
t259 = pkin(11) + qJ(4);
t255 = qJ(5) + t259;
t246 = sin(t255);
t260 = qJD(4) + qJD(5);
t267 = cos(qJ(6));
t265 = sin(qJ(6));
t406 = Icges(7,4) * t267;
t309 = -Icges(7,2) * t265 + t406;
t247 = cos(t255);
t382 = t247 * t260;
t407 = Icges(7,4) * t265;
t135 = t309 * t382 + (Icges(7,6) * t260 + (-Icges(7,2) * t267 - t407) * qJD(6)) * t246;
t195 = -Icges(7,6) * t247 + t246 * t309;
t313 = Icges(7,1) * t267 - t407;
t196 = -Icges(7,5) * t247 + t246 * t313;
t442 = -t135 * t265 + (-t195 * t267 - t196 * t265) * qJD(6);
t306 = Icges(7,5) * t267 - Icges(7,6) * t265;
t134 = t306 * t382 + (Icges(7,3) * t260 + (-Icges(7,5) * t265 - Icges(7,6) * t267) * qJD(6)) * t246;
t375 = t260 * t265;
t441 = -t195 * t375 - t134;
t261 = qJ(1) + pkin(10);
t254 = cos(t261);
t217 = Icges(6,5) * t246 + Icges(6,6) * t247;
t285 = t217 * t260;
t252 = sin(t261);
t307 = Icges(6,5) * t247 - Icges(6,6) * t246;
t175 = -Icges(6,3) * t254 + t252 * t307;
t431 = qJD(1) * t175;
t117 = -t254 * t285 - t431;
t176 = Icges(6,3) * t252 + t254 * t307;
t362 = qJD(1) * t176;
t118 = -t252 * t285 + t362;
t408 = Icges(6,4) * t247;
t310 = -Icges(6,2) * t246 + t408;
t177 = -Icges(6,6) * t254 + t252 * t310;
t409 = Icges(6,4) * t246;
t218 = Icges(6,2) * t247 + t409;
t388 = t218 * t260;
t119 = -qJD(1) * t177 - t254 * t388;
t314 = Icges(6,1) * t247 - t409;
t179 = -Icges(6,5) * t254 + t252 * t314;
t219 = Icges(6,1) * t246 + t408;
t387 = t219 * t260;
t121 = -qJD(1) * t179 - t254 * t387;
t178 = Icges(6,6) * t252 + t254 * t310;
t180 = Icges(6,5) * t252 + t254 * t314;
t298 = t178 * t246 - t180 * t247;
t122 = qJD(1) * t180 - t252 * t387;
t334 = t177 * t260 - t122;
t120 = qJD(1) * t178 - t252 * t388;
t336 = t179 * t260 + t120;
t384 = t246 * t260;
t299 = t177 * t246 - t179 * t247;
t435 = t254 * t299;
t85 = -t175 * t254 - t252 * t299;
t436 = t252 * t298;
t86 = -t176 * t254 - t436;
t333 = -qJD(6) * t247 + qJD(1);
t275 = t246 * t375 + t267 * t333;
t332 = qJD(1) * t247 - qJD(6);
t377 = t254 * t265;
t126 = t252 * t275 - t332 * t377;
t374 = t260 * t267;
t274 = -t246 * t374 + t265 * t333;
t376 = t254 * t267;
t127 = t252 * t274 + t332 * t376;
t380 = t252 * t265;
t206 = -t247 * t380 - t376;
t379 = t252 * t267;
t207 = t247 * t379 - t377;
t386 = t246 * t252;
t139 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t386;
t141 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t386;
t143 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t386;
t359 = qJD(1) * t254;
t381 = t252 * t260;
t278 = t246 * t359 + t247 * t381;
t71 = Icges(7,5) * t127 + Icges(7,6) * t126 + Icges(7,3) * t278;
t73 = Icges(7,4) * t127 + Icges(7,2) * t126 + Icges(7,6) * t278;
t75 = Icges(7,1) * t127 + Icges(7,4) * t126 + Icges(7,5) * t278;
t15 = t126 * t141 + t127 * t143 + t139 * t278 + t206 * t73 + t207 * t75 + t386 * t71;
t208 = -t247 * t377 + t379;
t209 = t247 * t376 + t380;
t385 = t246 * t254;
t140 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t385;
t142 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t385;
t144 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t385;
t124 = t254 * t275 + t332 * t380;
t125 = t254 * t274 - t332 * t379;
t360 = qJD(1) * t252;
t343 = t246 * t360;
t378 = t254 * t260;
t349 = t247 * t378;
t277 = -t343 + t349;
t70 = Icges(7,5) * t125 + Icges(7,6) * t124 + Icges(7,3) * t277;
t72 = Icges(7,4) * t125 + Icges(7,2) * t124 + Icges(7,6) * t277;
t74 = Icges(7,1) * t125 + Icges(7,4) * t124 + Icges(7,5) * t277;
t16 = t126 * t142 + t127 * t144 + t140 * t278 + t206 * t72 + t207 * t74 + t386 * t70;
t52 = t139 * t386 + t141 * t206 + t143 * t207;
t53 = t140 * t386 + t142 * t206 + t144 * t207;
t321 = t252 * t52 + t254 * t53;
t9 = qJD(1) * t321 - t15 * t254 + t16 * t252;
t439 = -(t254 * t118 + (t86 + t435) * qJD(1)) * t254 - (t85 * qJD(1) + (-t119 * t246 + t121 * t247 - t178 * t382 - t180 * t384 + t362) * t252 + (-t117 + t334 * t247 + t336 * t246 + (-t175 - t298) * qJD(1)) * t254) * t252 - t9;
t429 = 2 * m(5);
t251 = sin(t259);
t415 = rSges(5,2) * t251;
t253 = cos(t259);
t417 = rSges(5,1) * t253;
t327 = -t415 + t417;
t192 = -rSges(5,3) * t254 + t252 * t327;
t355 = t254 * t415;
t245 = t252 * rSges(5,3);
t365 = t254 * t417 + t245;
t193 = -t355 + t365;
t228 = rSges(5,1) * t251 + rSges(5,2) * t253;
t288 = qJD(4) * t228;
t342 = t251 * t360;
t272 = rSges(5,2) * t342 + rSges(5,3) * t359 - t254 * t288;
t437 = t252 * t288;
t63 = (qJD(1) * t192 + t272) * t254 + (-t437 + (-t193 - t355 + t245) * qJD(1)) * t252;
t438 = t429 * t63;
t410 = Icges(5,4) * t253;
t312 = -Icges(5,2) * t251 + t410;
t189 = Icges(5,6) * t252 + t254 * t312;
t411 = Icges(5,4) * t251;
t316 = Icges(5,1) * t253 - t411;
t191 = Icges(5,5) * t252 + t254 * t316;
t296 = t189 * t251 - t191 * t253;
t434 = t296 * t254;
t324 = -rSges(7,1) * t207 - rSges(7,2) * t206;
t146 = rSges(7,3) * t386 - t324;
t147 = t209 * rSges(7,1) + t208 * rSges(7,2) + rSges(7,3) * t385;
t433 = -t252 * t146 - t254 * t147;
t264 = -pkin(7) - qJ(3);
t263 = cos(pkin(11));
t248 = t263 * pkin(3) + pkin(2);
t291 = -t248 - t327;
t421 = sin(qJ(1)) * pkin(1);
t159 = -t421 + (rSges(5,3) - t264) * t254 + t291 * t252;
t257 = cos(qJ(1)) * pkin(1);
t160 = -t252 * t264 + t257 + (t248 - t415) * t254 + t365;
t432 = t159 * t254 + t160 * t252;
t308 = Icges(5,5) * t253 - Icges(5,6) * t251;
t186 = -Icges(5,3) * t254 + t252 * t308;
t188 = -Icges(5,6) * t254 + t252 * t312;
t190 = -Icges(5,5) * t254 + t252 * t316;
t295 = t218 * t246 - t219 * t247;
t430 = qJD(1) * t295 + t307 * t260;
t258 = -pkin(8) + t264;
t363 = t258 - t264;
t234 = pkin(4) * t253 + t248;
t366 = t234 - t248;
t163 = t252 * t366 + t254 * t363;
t293 = rSges(4,1) * t263 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t412 = rSges(4,3) + qJ(3);
t166 = t252 * t412 + t254 * t293 + t257;
t428 = 2 * m(6);
t427 = 2 * m(7);
t249 = t252 ^ 2;
t250 = t254 ^ 2;
t426 = t252 / 0.2e1;
t425 = -t254 / 0.2e1;
t424 = -rSges(7,3) - pkin(9);
t423 = m(5) * t228;
t220 = rSges(6,1) * t246 + rSges(6,2) * t247;
t422 = m(6) * t220;
t420 = pkin(4) * t251;
t419 = pkin(5) * t246;
t418 = pkin(5) * t247;
t416 = rSges(6,1) * t247;
t304 = -t141 * t265 + t143 * t267;
t20 = (t260 * t304 - t71) * t247 + (t139 * t260 - t265 * t73 + t267 * t75 + (-t141 * t267 - t143 * t265) * qJD(6)) * t246;
t414 = t20 * t254;
t303 = -t142 * t265 + t144 * t267;
t21 = (t260 * t303 - t70) * t247 + (t140 * t260 - t265 * t72 + t267 * t74 + (-t142 * t267 - t144 * t265) * qJD(6)) * t246;
t413 = t21 * t252;
t244 = t252 * rSges(6,3);
t338 = -t220 - t420;
t174 = t338 * t254;
t396 = t174 * t254;
t395 = t188 * t253;
t394 = t189 * t253;
t393 = t190 * t251;
t392 = t191 * t251;
t326 = -rSges(6,2) * t246 + t416;
t205 = t326 * t260;
t389 = t205 * t252;
t383 = t247 * t254;
t323 = rSges(7,1) * t267 - rSges(7,2) * t265;
t145 = t323 * t382 + (rSges(7,3) * t260 + (-rSges(7,1) * t265 - rSges(7,2) * t267) * qJD(6)) * t246;
t330 = pkin(9) * t246 + t418;
t373 = -t330 * t260 - t145;
t199 = pkin(5) * t383 + pkin(9) * t385;
t372 = -t147 - t199;
t212 = t254 * t234;
t164 = -t248 * t254 - t252 * t363 + t212;
t371 = t252 * t163 + t254 * t164;
t181 = -rSges(6,3) * t254 + t252 * t326;
t182 = rSges(6,1) * t383 - rSges(6,2) * t385 + t244;
t109 = t252 * t181 + t254 * t182;
t197 = -rSges(7,3) * t247 + t246 * t323;
t172 = t197 * t360;
t222 = -pkin(9) * t247 + t419;
t370 = t222 * t360 + t172;
t369 = -t197 - t222;
t368 = rSges(6,2) * t343 + rSges(6,3) * t359;
t358 = qJD(4) * t251;
t354 = pkin(4) * t358;
t367 = -t252 * t354 - t258 * t360;
t364 = t249 + t250;
t187 = Icges(5,3) * t252 + t254 * t308;
t361 = qJD(1) * t187;
t357 = qJD(4) * t253;
t353 = pkin(4) * t357;
t62 = -t140 * t247 + t246 * t303;
t194 = -Icges(7,3) * t247 + t246 * t306;
t84 = t194 * t385 + t195 * t208 + t196 * t209;
t352 = t62 / 0.2e1 + t84 / 0.2e1;
t61 = -t139 * t247 + t246 * t304;
t83 = t194 * t386 + t195 * t206 + t196 * t207;
t351 = t83 / 0.2e1 + t61 / 0.2e1;
t279 = -t246 * t378 - t247 * t360;
t289 = t220 * t260;
t348 = t252 * (-t252 * t289 + (t254 * t326 + t244) * qJD(1)) + t254 * (rSges(6,1) * t279 - rSges(6,2) * t349 + t368) + t181 * t359;
t136 = t313 * t382 + (Icges(7,5) * t260 + (-Icges(7,1) * t265 - t406) * qJD(6)) * t246;
t347 = t246 * t267 * t136 + t247 * t196 * t374 + t194 * t384;
t239 = t264 * t360;
t346 = t252 * (t359 * t366 + t239 + t367) + t254 * (-qJD(1) * t163 - t254 * t354) + t163 * t359;
t345 = t125 * rSges(7,1) + t124 * rSges(7,2) + rSges(7,3) * t349;
t243 = qJD(3) * t254;
t344 = t243 - t367;
t340 = t360 / 0.2e1;
t339 = t359 / 0.2e1;
t157 = t369 * t254;
t337 = t180 * t260 + t119;
t335 = -t178 * t260 + t121;
t198 = t330 * t252;
t64 = t252 * t198 + t254 * t199 - t433;
t331 = t369 - t420;
t329 = -t252 * t258 + t212 + t257;
t325 = t127 * rSges(7,1) + t126 * rSges(7,2);
t322 = -t254 * t258 - t421;
t38 = t252 * t53 - t254 * t52;
t54 = t139 * t385 + t141 * t208 + t143 * t209;
t55 = t140 * t385 + t142 * t208 + t144 * t209;
t39 = t252 * t55 - t254 * t54;
t320 = t252 * t54 + t254 * t55;
t319 = t62 * t252 - t61 * t254;
t318 = t252 * t61 + t62 * t254;
t13 = t124 * t141 + t125 * t143 + t139 * t277 + t208 * t73 + t209 * t75 + t385 * t71;
t14 = t124 * t142 + t125 * t144 + t140 * t277 + t208 * t72 + t209 * t74 + t385 * t70;
t8 = qJD(1) * t320 - t13 * t254 + t14 * t252;
t87 = t175 * t252 - t435;
t88 = t176 * t252 - t254 * t298;
t317 = t39 * t359 + t38 * t360 + (-t87 * t359 - t85 * t360) * t254 + (t8 + (t88 * qJD(1) + (t120 * t246 - t122 * t247 + t177 * t382 + t179 * t384 - t431) * t254) * t254 + t86 * t360 + t88 * t359 + ((t87 + t436) * qJD(1) + (-t118 + t335 * t247 - t337 * t246 + (t176 - t299) * qJD(1)) * t254 + t252 * t117) * t252) * t252;
t315 = Icges(5,1) * t251 + t410;
t311 = Icges(5,2) * t253 + t411;
t290 = -t234 - t326;
t137 = -t421 + (rSges(6,3) - t258) * t254 + t290 * t252;
t138 = t182 + t329;
t305 = t137 * t254 + t138 * t252;
t302 = t146 * t254 - t147 * t252;
t297 = t188 * t251 - t190 * t253;
t294 = -t353 + t373;
t223 = pkin(9) * t349;
t76 = -rSges(7,3) * t343 + t345;
t77 = rSges(7,3) * t278 + t325;
t292 = (t146 + t198) * t359 + (pkin(5) * t279 - pkin(9) * t343 + t223 + t76) * t254 + (t77 + t278 * pkin(9) + (-t246 * t381 + t247 * t359) * pkin(5)) * t252;
t131 = t331 * t254;
t284 = qJD(4) * t315;
t283 = qJD(4) * t311;
t282 = qJD(4) * (-Icges(5,5) * t251 - Icges(5,6) * t253);
t281 = t246 * t424 - t234 - t418;
t276 = t254 * t439 + t317;
t25 = t246 * t321 - t247 * t83;
t26 = t246 * t320 - t247 * t84;
t30 = t124 * t195 + t125 * t196 + t134 * t385 + t135 * t208 + t136 * t209 + t194 * t277;
t3 = (t260 * t320 - t30) * t247 + (-qJD(1) * t39 + t13 * t252 + t14 * t254 + t260 * t84) * t246;
t31 = t126 * t195 + t127 * t196 + t134 * t386 + t135 * t206 + t136 * t207 + t194 * t278;
t4 = (t260 * t321 - t31) * t247 + (-qJD(1) * t38 + t15 * t252 + t16 * t254 + t260 * t83) * t246;
t273 = t3 * t426 + t4 * t425 - t247 * (qJD(1) * t318 + t413 - t414) / 0.2e1 + t25 * t340 - t39 * t343 / 0.2e1 + t319 * t384 / 0.2e1 + t9 * t386 / 0.2e1 + t8 * t385 / 0.2e1 + (t252 * t38 + t254 * t39) * t382 / 0.2e1 + (t246 * t38 + t26) * t339;
t271 = t252 * t281 + t322;
t165 = -t252 * t293 + t254 * t412 - t421;
t203 = t310 * t260;
t204 = t314 * t260;
t270 = qJD(1) * t217 + (t204 - t388) * t247 + (-t203 - t387) * t246;
t269 = -t414 / 0.2e1 + t413 / 0.2e1 + (t246 * t335 + t247 * t337 + t252 * t430 + t270 * t254 + t30) * t426 + (-t246 * t334 + t247 * t336 + t270 * t252 - t254 * t430 + t31) * t425 + (t177 * t247 + t179 * t246 - t217 * t254 - t252 * t295 + t61 + t83) * t340 + (t178 * t247 + t180 * t246 + t217 * t252 - t254 * t295 + t62 + t84) * t339;
t242 = qJD(3) * t252;
t231 = pkin(4) * t342;
t216 = t327 * qJD(4);
t173 = t338 * t252;
t156 = t369 * t252;
t155 = -qJD(1) * t166 + t243;
t154 = qJD(1) * t165 + t242;
t149 = t252 * t282 + t361;
t148 = -qJD(1) * t186 + t254 * t282;
t130 = t331 * t252;
t108 = -t220 * t359 - t389 + (-t251 * t359 - t252 * t357) * pkin(4);
t107 = t220 * t360 + t231 + (-t205 - t353) * t254;
t102 = t239 + t243 + t437 + (t254 * t291 - t245 - t257) * qJD(1);
t101 = t242 + (-t421 - t254 * t264 + (-t248 - t417) * t252) * qJD(1) + t272;
t100 = -t147 * t247 - t197 * t385;
t99 = t146 * t247 + t197 * t386;
t98 = -t194 * t247 + (-t195 * t265 + t196 * t267) * t246;
t97 = t187 * t252 - t434;
t96 = t186 * t252 - t254 * t297;
t95 = -t187 * t254 - t252 * t296;
t94 = -t186 * t254 - t252 * t297;
t93 = t329 - t372;
t92 = t271 + t324;
t91 = t220 * t381 + (t254 * t290 - t244 - t257) * qJD(1) + t344;
t90 = t242 + (-t289 - t354) * t254 + ((-t234 - t416) * t252 + t322) * qJD(1) + t368;
t89 = t98 * t384;
t82 = t302 * t246;
t81 = qJD(1) * t157 + t252 * t373;
t80 = t254 * t373 + t370;
t67 = qJD(1) * t131 + t252 * t294;
t66 = t254 * t294 + t231 + t370;
t65 = t109 + t371;
t58 = -t182 * t360 + t348;
t47 = (t247 * t424 + t419) * t381 + (t254 * t281 - t257) * qJD(1) - t325 + t344;
t46 = t223 + t242 + (-pkin(5) * t384 - t354) * t254 + t271 * qJD(1) + t345;
t45 = t64 + t371;
t44 = (t197 * t381 + t77) * t247 + (t145 * t252 - t146 * t260 + t197 * t359) * t246;
t43 = (-t197 * t378 - t76) * t247 + (-t145 * t254 + t147 * t260 + t172) * t246;
t42 = t246 * t442 + t441 * t247 + t347;
t40 = (-t164 - t182) * t360 + t346 + t348;
t27 = t302 * t382 + (qJD(1) * t433 - t252 * t76 + t254 * t77) * t246;
t22 = t360 * t372 + t292;
t17 = (-t164 + t372) * t360 + t292 + t346;
t1 = [t347 + (t46 * t93 + t47 * t92) * t427 + (t137 * t91 + t138 * t90) * t428 + (t101 * t160 + t102 * t159) * t429 + 0.2e1 * m(4) * (t154 * t166 + t155 * t165) + t219 * t382 - t218 * t384 + (-t311 + t316) * t358 + (t312 + t315) * t357 + (t203 + t441) * t247 + (t204 + t442) * t246; 0; 0; m(7) * (t252 * t47 - t254 * t46 + (t252 * t93 + t254 * t92) * qJD(1)) + m(6) * (qJD(1) * t305 + t252 * t91 - t254 * t90) + m(5) * (qJD(1) * t432 - t101 * t254 + t102 * t252) + m(4) * (-t154 * t254 + t155 * t252 + (t165 * t254 + t166 * t252) * qJD(1)); 0; 0; t269 + m(5) * ((-t101 * t252 - t102 * t254) * t228 - t432 * t216) + (-qJD(4) * t297 + (qJD(1) * t189 - t252 * t283) * t253 + (qJD(1) * t191 - t252 * t284) * t251) * t425 + (-qJD(4) * t296 + (-qJD(1) * t188 - t254 * t283) * t253 + (-qJD(1) * t190 - t254 * t284) * t251) * t426 + ((-t160 * t423 + t394 / 0.2e1 + t392 / 0.2e1) * t254 + (t159 * t423 + t395 / 0.2e1 + t393 / 0.2e1) * t252) * qJD(1) + (t250 / 0.2e1 + t249 / 0.2e1) * t308 * qJD(4) + m(6) * (t107 * t137 + t108 * t138 + t173 * t90 + t174 * t91) + m(7) * (t130 * t46 + t131 * t47 + t66 * t92 + t67 * t93); m(5) * t63 + m(6) * t40 + m(7) * t17; m(6) * (t107 * t252 - t108 * t254 + (t173 * t252 + t396) * qJD(1)) + m(7) * (t252 * t66 - t254 * t67 + (t130 * t252 + t131 * t254) * qJD(1)); (t130 * t67 + t131 * t66 + t45 * t17) * t427 + (t107 * t174 + t108 * t173 + t65 * t40) * t428 + t364 * t228 * t216 * t429 + t317 + (t192 * t438 + t97 * t359 + t95 * t360 + (t96 * qJD(1) + (t296 * qJD(1) + t148) * t252) * t252) * t252 + (-t96 * t359 - t94 * t360 + t193 * t438 + (-t95 * qJD(1) + (-qJD(1) * t297 - t149) * t254) * t254 + ((t188 * t357 + t190 * t358 + t148 - (t393 + t395) * qJD(4)) * t254 + (-t149 + (-t392 - t394) * qJD(4) + t189 * t357 + t191 * t358 - t361) * t252 + (t97 - t94 + t434 + (t187 - t297) * t252) * qJD(1)) * t252 + t439) * t254; t269 + m(7) * (t156 * t46 + t157 * t47 + t80 * t92 + t81 * t93) - m(6) * t305 * t205 + (-t252 * t90 - t254 * t91 + (t137 * t252 - t138 * t254) * qJD(1)) * t422; m(6) * t58 + m(7) * t22; m(7) * (t252 * t80 - t254 * t81 + (t156 * t252 + t157 * t254) * qJD(1)); m(7) * (t130 * t81 + t131 * t80 + t156 * t67 + t157 * t66 + t64 * t17 + t22 * t45) + m(6) * (t109 * t40 - t173 * t389 - t205 * t396 + t58 * t65) + (-t107 * t254 - t108 * t252 + (-t173 * t254 + t174 * t252) * qJD(1)) * t422 + t276; (t205 * t220 * t364 + t109 * t58) * t428 + (t156 * t81 + t157 * t80 + t64 * t22) * t427 + t276; t89 + m(7) * (t100 * t46 + t43 * t93 + t44 * t92 + t47 * t99) + (-t42 + (t252 * t351 + t254 * t352) * t260) * t247 + ((t21 / 0.2e1 + t30 / 0.2e1) * t254 + (t31 / 0.2e1 + t20 / 0.2e1) * t252 + (-t252 * t352 + t254 * t351) * qJD(1)) * t246; m(7) * t27; m(7) * (t252 * t44 - t254 * t43 + (t100 * t252 + t254 * t99) * qJD(1)); t273 + m(7) * (t100 * t67 + t130 * t43 + t131 * t44 + t17 * t82 + t27 * t45 + t66 * t99); t273 + m(7) * (t100 * t81 + t156 * t43 + t157 * t44 + t22 * t82 + t27 * t64 + t80 * t99); (t100 * t43 + t27 * t82 + t44 * t99) * t427 + (t42 * t247 - t89 + (-t247 * t318 + t252 * t25 + t254 * t26) * t260) * t247 + (t254 * t3 + t252 * t4 + t318 * t384 + (-t20 * t252 - t21 * t254 - t260 * t98) * t247 + (t247 * t319 + t254 * t25 - t252 * t26) * qJD(1)) * t246;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
