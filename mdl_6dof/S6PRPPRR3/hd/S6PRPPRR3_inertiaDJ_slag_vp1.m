% Calculate time derivative of joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:30
% EndTime: 2019-03-08 19:21:51
% DurationCPUTime: 16.27s
% Computational Cost: add. (49574->842), mult. (141575->1211), div. (0->0), fcn. (167637->12), ass. (0->338)
t425 = Icges(3,4) - Icges(4,5);
t424 = Icges(3,1) + Icges(4,1);
t423 = Icges(4,4) + Icges(3,5);
t422 = Icges(3,2) + Icges(4,3);
t421 = Icges(3,6) - Icges(4,6);
t319 = sin(qJ(2));
t420 = t425 * t319;
t322 = cos(qJ(2));
t419 = t425 * t322;
t313 = sin(pkin(11));
t378 = cos(pkin(11));
t340 = t378 * t319;
t315 = sin(pkin(6));
t353 = qJD(2) * t322;
t342 = t315 * t353;
t355 = qJD(2) * t315;
t276 = t313 * t342 - t340 * t355;
t286 = (t313 * t319 + t322 * t378) * t315;
t277 = qJD(2) * t286;
t398 = -Icges(5,5) * t277 + Icges(5,6) * t276 + (-t319 * t421 + t322 * t423) * t355;
t314 = sin(pkin(10));
t316 = cos(pkin(10));
t317 = cos(pkin(6));
t369 = t317 * t319;
t301 = t314 * t322 + t316 * t369;
t368 = t317 * t322;
t327 = -t314 * t319 + t316 * t368;
t371 = t315 * t316;
t418 = t301 * t425 + t327 * t422 - t371 * t421;
t302 = t314 * t368 + t316 * t319;
t303 = -t314 * t369 + t316 * t322;
t372 = t314 * t315;
t417 = t302 * t422 - t303 * t425 - t372 * t421;
t416 = -t301 * t424 - t425 * t327 + t423 * t371;
t415 = -t425 * t302 + t303 * t424 + t423 * t372;
t288 = t327 * qJD(2);
t289 = t301 * qJD(2);
t414 = -t288 * t425 + t289 * t422;
t290 = t302 * qJD(2);
t354 = qJD(2) * t319;
t341 = t317 * t354;
t291 = -t314 * t341 + t316 * t353;
t413 = t290 * t425 + t291 * t422;
t412 = -t288 * t424 + t425 * t289;
t411 = -t290 * t424 - t425 * t291;
t410 = 2 * m(3);
t409 = 0.2e1 * t315;
t219 = t288 * t313 - t289 * t378;
t221 = -t290 * t313 - t291 * t378;
t254 = t301 * t313 + t327 * t378;
t256 = -t302 * t378 + t303 * t313;
t255 = t301 * t378 - t313 * t327;
t321 = cos(qJ(5));
t381 = sin(qJ(5));
t344 = t315 * t381;
t208 = t255 * t321 + t316 * t344;
t318 = sin(qJ(6));
t320 = cos(qJ(6));
t164 = -t208 * t318 + t254 * t320;
t165 = t208 * t320 + t254 * t318;
t370 = t315 * t321;
t323 = -t255 * t381 + t316 * t370;
t115 = Icges(7,5) * t165 + Icges(7,6) * t164 - Icges(7,3) * t323;
t117 = Icges(7,4) * t165 + Icges(7,2) * t164 - Icges(7,6) * t323;
t119 = Icges(7,1) * t165 + Icges(7,4) * t164 - Icges(7,5) * t323;
t220 = t288 * t378 + t289 * t313;
t159 = qJD(5) * t323 + t220 * t321;
t121 = -qJD(6) * t165 - t159 * t318 + t219 * t320;
t122 = qJD(6) * t164 + t159 * t320 + t219 * t318;
t158 = qJD(5) * t208 + t220 * t381;
t73 = Icges(7,5) * t122 + Icges(7,6) * t121 + Icges(7,3) * t158;
t75 = Icges(7,4) * t122 + Icges(7,2) * t121 + Icges(7,6) * t158;
t77 = Icges(7,1) * t122 + Icges(7,4) * t121 + Icges(7,5) * t158;
t16 = t115 * t158 + t117 * t121 + t119 * t122 + t164 * t75 + t165 * t77 - t323 * t73;
t257 = t302 * t313 + t303 * t378;
t334 = t314 * t344;
t210 = t257 * t321 - t334;
t166 = -t210 * t318 + t256 * t320;
t167 = t210 * t320 + t256 * t318;
t324 = -t257 * t381 - t314 * t370;
t116 = Icges(7,5) * t167 + Icges(7,6) * t166 - Icges(7,3) * t324;
t118 = Icges(7,4) * t167 + Icges(7,2) * t166 - Icges(7,6) * t324;
t120 = Icges(7,1) * t167 + Icges(7,4) * t166 - Icges(7,5) * t324;
t222 = -t290 * t378 + t291 * t313;
t161 = qJD(5) * t324 + t222 * t321;
t123 = -qJD(6) * t167 - t161 * t318 + t221 * t320;
t124 = qJD(6) * t166 + t161 * t320 + t221 * t318;
t351 = qJD(5) * t321;
t160 = -qJD(5) * t334 + t222 * t381 + t257 * t351;
t74 = Icges(7,5) * t124 + Icges(7,6) * t123 + Icges(7,3) * t160;
t76 = Icges(7,4) * t124 + Icges(7,2) * t123 + Icges(7,6) * t160;
t78 = Icges(7,1) * t124 + Icges(7,4) * t123 + Icges(7,5) * t160;
t17 = t116 * t158 + t118 * t121 + t120 * t122 + t164 * t76 + t165 * t78 - t323 * t74;
t287 = (-t313 * t322 + t340) * t315;
t343 = t317 * t381;
t268 = t287 * t321 - t343;
t205 = -t268 * t318 + t286 * t320;
t206 = t268 * t320 + t286 * t318;
t267 = t287 * t381 + t317 * t321;
t143 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t267;
t144 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t267;
t145 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t267;
t204 = -qJD(5) * t267 + t277 * t321;
t147 = -qJD(6) * t206 - t204 * t318 + t276 * t320;
t148 = qJD(6) * t205 + t204 * t320 + t276 * t318;
t203 = -qJD(5) * t343 + t277 * t381 + t287 * t351;
t97 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t203;
t98 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t203;
t99 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t203;
t32 = t121 * t144 + t122 * t145 + t143 * t158 + t164 * t98 + t165 * t99 - t323 * t97;
t52 = -t115 * t323 + t117 * t164 + t119 * t165;
t53 = -t116 * t323 + t118 * t164 + t120 * t165;
t70 = -t143 * t323 + t144 * t164 + t145 * t165;
t3 = t16 * t254 + t17 * t256 + t219 * t52 + t221 * t53 + t276 * t70 + t286 * t32;
t107 = Icges(6,5) * t159 - Icges(6,6) * t158 + Icges(6,3) * t219;
t109 = Icges(6,4) * t159 - Icges(6,2) * t158 + Icges(6,6) * t219;
t111 = Icges(6,1) * t159 - Icges(6,4) * t158 + Icges(6,5) * t219;
t135 = Icges(6,5) * t208 + Icges(6,6) * t323 + Icges(6,3) * t254;
t137 = Icges(6,4) * t208 + Icges(6,2) * t323 + Icges(6,6) * t254;
t139 = Icges(6,1) * t208 + Icges(6,4) * t323 + Icges(6,5) * t254;
t39 = t107 * t254 + t109 * t323 + t111 * t208 + t135 * t219 - t137 * t158 + t139 * t159;
t108 = Icges(6,5) * t161 - Icges(6,6) * t160 + Icges(6,3) * t221;
t110 = Icges(6,4) * t161 - Icges(6,2) * t160 + Icges(6,6) * t221;
t112 = Icges(6,1) * t161 - Icges(6,4) * t160 + Icges(6,5) * t221;
t136 = Icges(6,5) * t210 + Icges(6,6) * t324 + Icges(6,3) * t256;
t138 = Icges(6,4) * t210 + Icges(6,2) * t324 + Icges(6,6) * t256;
t140 = Icges(6,1) * t210 + Icges(6,4) * t324 + Icges(6,5) * t256;
t40 = t108 * t254 + t110 * t323 + t112 * t208 + t136 * t219 - t138 * t158 + t140 * t159;
t149 = Icges(6,5) * t204 - Icges(6,6) * t203 + Icges(6,3) * t276;
t150 = Icges(6,4) * t204 - Icges(6,2) * t203 + Icges(6,6) * t276;
t151 = Icges(6,1) * t204 - Icges(6,4) * t203 + Icges(6,5) * t276;
t191 = Icges(6,5) * t268 - Icges(6,6) * t267 + Icges(6,3) * t286;
t192 = Icges(6,4) * t268 - Icges(6,2) * t267 + Icges(6,6) * t286;
t193 = Icges(6,1) * t268 - Icges(6,4) * t267 + Icges(6,5) * t286;
t50 = t149 * t254 + t150 * t323 + t151 * t208 - t158 * t192 + t159 * t193 + t191 * t219;
t82 = t135 * t254 + t137 * t323 + t139 * t208;
t83 = t136 * t254 + t138 * t323 + t140 * t208;
t94 = t191 * t254 + t192 * t323 + t193 * t208;
t408 = t219 * t82 + t221 * t83 + t254 * t39 + t256 * t40 + t276 * t94 + t286 * t50 + t3;
t18 = t115 * t160 + t117 * t123 + t119 * t124 + t166 * t75 + t167 * t77 - t324 * t73;
t19 = t116 * t160 + t118 * t123 + t120 * t124 + t166 * t76 + t167 * t78 - t324 * t74;
t33 = t123 * t144 + t124 * t145 + t143 * t160 + t166 * t98 + t167 * t99 - t324 * t97;
t54 = -t115 * t324 + t117 * t166 + t119 * t167;
t55 = -t116 * t324 + t118 * t166 + t120 * t167;
t71 = -t143 * t324 + t144 * t166 + t145 * t167;
t4 = t18 * t254 + t19 * t256 + t219 * t54 + t221 * t55 + t276 * t71 + t286 * t33;
t41 = t107 * t256 + t109 * t324 + t111 * t210 + t135 * t221 - t137 * t160 + t139 * t161;
t42 = t108 * t256 + t110 * t324 + t112 * t210 + t136 * t221 - t138 * t160 + t140 * t161;
t51 = t149 * t256 + t150 * t324 + t151 * t210 - t160 * t192 + t161 * t193 + t191 * t221;
t84 = t135 * t256 + t137 * t324 + t139 * t210;
t85 = t136 * t256 + t138 * t324 + t140 * t210;
t95 = t191 * t256 + t192 * t324 + t193 * t210;
t407 = t219 * t84 + t221 * t85 + t254 * t41 + t256 * t42 + t276 * t95 + t286 * t51 + t4;
t250 = rSges(3,1) * t288 - rSges(3,2) * t289;
t253 = -rSges(3,1) * t290 - rSges(3,2) * t291;
t406 = m(3) * (t250 * t314 + t253 * t316) * t315;
t405 = -t421 * t317 + (-t322 * t422 - t420) * t315;
t404 = t423 * t317 + (t319 * t424 + t419) * t315;
t402 = (t319 * t422 - t419) * t355;
t401 = (t322 * t424 - t420) * t355;
t400 = -Icges(5,5) * t222 + Icges(5,6) * t221 - t290 * t423 - t291 * t421;
t399 = Icges(5,5) * t220 - Icges(5,6) * t219 - t288 * t423 + t289 * t421;
t283 = t317 * rSges(3,3) + (rSges(3,1) * t319 + rSges(3,2) * t322) * t315;
t299 = (rSges(3,1) * t322 - rSges(3,2) * t319) * t355;
t397 = t283 * t299 * t410;
t396 = 2 * m(6);
t395 = 2 * m(7);
t20 = t115 * t203 + t117 * t147 + t119 * t148 + t205 * t75 + t206 * t77 + t267 * t73;
t21 = t116 * t203 + t118 * t147 + t120 * t148 + t205 * t76 + t206 * t78 + t267 * t74;
t38 = t143 * t203 + t144 * t147 + t145 * t148 + t205 * t98 + t206 * t99 + t267 * t97;
t58 = t115 * t267 + t117 * t205 + t119 * t206;
t59 = t116 * t267 + t118 * t205 + t120 * t206;
t86 = t143 * t267 + t144 * t205 + t145 * t206;
t5 = t158 * t58 + t160 * t59 - t20 * t323 + t203 * t86 - t21 * t324 + t267 * t38;
t393 = t5 / 0.2e1;
t392 = t158 / 0.2e1;
t391 = t160 / 0.2e1;
t390 = t203 / 0.2e1;
t389 = -t323 / 0.2e1;
t388 = -t324 / 0.2e1;
t387 = t219 / 0.2e1;
t386 = t221 / 0.2e1;
t385 = t267 / 0.2e1;
t384 = t276 / 0.2e1;
t383 = t314 / 0.2e1;
t382 = -t316 / 0.2e1;
t79 = rSges(7,1) * t122 + rSges(7,2) * t121 + rSges(7,3) * t158;
t380 = pkin(5) * t159 + pkin(9) * t158 + t79;
t80 = rSges(7,1) * t124 + rSges(7,2) * t123 + rSges(7,3) * t160;
t379 = pkin(5) * t161 + pkin(9) * t160 + t80;
t373 = (t191 * t286 - t192 * t267 + t193 * t268) * t276;
t100 = rSges(7,1) * t148 + rSges(7,2) * t147 + rSges(7,3) * t203;
t367 = pkin(5) * t204 + pkin(9) * t203 + t100;
t127 = rSges(7,1) * t165 + rSges(7,2) * t164 - rSges(7,3) * t323;
t366 = pkin(5) * t208 - pkin(9) * t323 + t127;
t128 = rSges(7,1) * t167 + rSges(7,2) * t166 - rSges(7,3) * t324;
t365 = pkin(5) * t210 - pkin(9) * t324 + t128;
t146 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t267;
t364 = pkin(5) * t268 + pkin(9) * t267 + t146;
t201 = pkin(2) * t288 + qJ(3) * t289 - qJD(3) * t327;
t202 = -pkin(2) * t290 + qJ(3) * t291 + qJD(3) * t302;
t363 = t201 * t372 + t202 * t371;
t198 = t317 * t202;
t352 = qJD(4) * t315;
t270 = -pkin(3) * t290 - t314 * t352;
t362 = t317 * t270 + t198;
t269 = pkin(3) * t288 + t316 * t352;
t361 = -t201 - t269;
t259 = pkin(2) * t301 - qJ(3) * t327;
t260 = pkin(2) * t303 + qJ(3) * t302;
t360 = t259 * t372 + t260 * t371;
t258 = t317 * t260;
t272 = pkin(3) * t303 - qJ(4) * t372;
t359 = t317 * t272 + t258;
t271 = pkin(3) * t301 + qJ(4) * t371;
t358 = -t259 - t271;
t273 = (-qJD(3) * t322 + (pkin(2) * t322 + qJ(3) * t319) * qJD(2)) * t315;
t357 = -pkin(3) * t342 + qJD(4) * t317 - t273;
t304 = (pkin(2) * t319 - qJ(3) * t322) * t315;
t356 = -pkin(3) * t315 * t319 + qJ(4) * t317 - t304;
t179 = pkin(4) * t222 + pkin(8) * t221;
t350 = t317 * t179 + t362;
t178 = pkin(4) * t220 + pkin(8) * t219;
t349 = -t178 + t361;
t190 = pkin(4) * t257 + pkin(8) * t256;
t348 = t317 * t190 + t359;
t189 = pkin(4) * t255 + pkin(8) * t254;
t347 = -t189 + t358;
t346 = -pkin(4) * t277 - pkin(8) * t276 + t357;
t345 = -pkin(4) * t287 - pkin(8) * t286 + t356;
t339 = (-t273 - (rSges(4,1) * t322 + rSges(4,3) * t319) * t355) * t315;
t338 = (-t317 * rSges(4,2) - (rSges(4,1) * t319 - rSges(4,3) * t322) * t315 - t304) * t315;
t337 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t336 = t269 * t372 + t270 * t371 + t363;
t335 = t271 * t372 + t272 * t371 + t360;
t333 = (-rSges(5,1) * t277 + rSges(5,2) * t276 + t357) * t315;
t332 = (-rSges(5,1) * t287 + rSges(5,2) * t286 + rSges(5,3) * t317 + t356) * t315;
t152 = rSges(6,1) * t204 - rSges(6,2) * t203 + rSges(6,3) * t276;
t331 = (-t152 + t346) * t315;
t194 = rSges(6,1) * t268 - rSges(6,2) * t267 + rSges(6,3) * t286;
t330 = (-t194 + t345) * t315;
t329 = t178 * t372 + t179 * t371 + t336;
t328 = t189 * t372 + t190 * t371 + t335;
t326 = (t346 - t367) * t315;
t325 = (t345 - t364) * t315;
t252 = -rSges(4,1) * t290 + rSges(4,3) * t291;
t249 = rSges(4,1) * t288 + rSges(4,3) * t289;
t234 = rSges(3,1) * t303 - rSges(3,2) * t302 + rSges(3,3) * t372;
t233 = rSges(4,1) * t303 + rSges(4,2) * t372 + rSges(4,3) * t302;
t232 = rSges(3,1) * t301 + rSges(3,2) * t327 - rSges(3,3) * t371;
t231 = rSges(4,1) * t301 - rSges(4,2) * t371 - rSges(4,3) * t327;
t216 = Icges(5,1) * t287 - Icges(5,4) * t286 - Icges(5,5) * t317;
t215 = Icges(5,4) * t287 - Icges(5,2) * t286 - Icges(5,6) * t317;
t213 = Icges(5,1) * t277 - Icges(5,4) * t276;
t212 = Icges(5,4) * t277 - Icges(5,2) * t276;
t185 = rSges(5,1) * t257 - rSges(5,2) * t256 - rSges(5,3) * t372;
t184 = rSges(5,1) * t255 - rSges(5,2) * t254 + rSges(5,3) * t371;
t183 = Icges(5,1) * t257 - Icges(5,4) * t256 - Icges(5,5) * t372;
t182 = Icges(5,1) * t255 - Icges(5,4) * t254 + Icges(5,5) * t371;
t181 = Icges(5,4) * t257 - Icges(5,2) * t256 - Icges(5,6) * t372;
t180 = Icges(5,4) * t255 - Icges(5,2) * t254 + Icges(5,6) * t371;
t177 = rSges(5,1) * t222 - rSges(5,2) * t221;
t176 = rSges(5,1) * t220 - rSges(5,2) * t219;
t175 = Icges(5,1) * t222 - Icges(5,4) * t221;
t174 = Icges(5,1) * t220 - Icges(5,4) * t219;
t173 = Icges(5,4) * t222 - Icges(5,2) * t221;
t172 = Icges(5,4) * t220 - Icges(5,2) * t219;
t154 = (-t231 - t259) * t317 + t316 * t338;
t153 = t233 * t317 + t314 * t338 + t258;
t142 = rSges(6,1) * t210 + rSges(6,2) * t324 + rSges(6,3) * t256;
t141 = rSges(6,1) * t208 + rSges(6,2) * t323 + rSges(6,3) * t254;
t134 = (-t201 - t249) * t317 + t316 * t339;
t133 = t252 * t317 + t314 * t339 + t198;
t132 = (t231 * t314 + t233 * t316) * t315 + t360;
t129 = (t249 * t314 + t252 * t316) * t315 + t363;
t126 = (-t184 + t358) * t317 + t316 * t332;
t125 = t185 * t317 + t314 * t332 + t359;
t114 = rSges(6,1) * t161 - rSges(6,2) * t160 + rSges(6,3) * t221;
t113 = rSges(6,1) * t159 - rSges(6,2) * t158 + rSges(6,3) * t219;
t106 = t142 * t286 - t194 * t256;
t105 = -t141 * t286 + t194 * t254;
t104 = (-t176 + t361) * t317 + t316 * t333;
t103 = t177 * t317 + t314 * t333 + t362;
t101 = (t184 * t314 + t185 * t316) * t315 + t335;
t96 = t141 * t256 - t142 * t254;
t93 = (t176 * t314 + t177 * t316) * t315 + t336;
t92 = t128 * t267 + t146 * t324;
t91 = -t127 * t267 - t146 * t323;
t90 = (-t141 + t347) * t317 + t316 * t330;
t89 = t142 * t317 + t314 * t330 + t348;
t88 = t136 * t286 - t138 * t267 + t140 * t268;
t87 = t135 * t286 - t137 * t267 + t139 * t268;
t81 = -t127 * t324 + t128 * t323;
t72 = (t141 * t314 + t142 * t316) * t315 + t328;
t69 = -t256 * t364 + t286 * t365;
t68 = t254 * t364 - t286 * t366;
t67 = (-t113 + t349) * t317 + t316 * t331;
t66 = t114 * t317 + t314 * t331 + t350;
t65 = t114 * t286 + t142 * t276 - t152 * t256 - t194 * t221;
t64 = -t113 * t286 - t141 * t276 + t152 * t254 + t194 * t219;
t63 = (t347 - t366) * t317 + t316 * t325;
t62 = t314 * t325 + t317 * t365 + t348;
t61 = t149 * t286 - t150 * t267 + t151 * t268 + t191 * t276 - t192 * t203 + t193 * t204;
t60 = -t254 * t365 + t256 * t366;
t57 = (t113 * t314 + t114 * t316) * t315 + t329;
t56 = t113 * t256 - t114 * t254 + t141 * t221 - t142 * t219;
t49 = (t314 * t366 + t316 * t365) * t315 + t328;
t48 = t108 * t286 - t110 * t267 + t112 * t268 + t136 * t276 - t138 * t203 + t140 * t204;
t47 = t107 * t286 - t109 * t267 + t111 * t268 + t135 * t276 - t137 * t203 + t139 * t204;
t46 = t100 * t324 + t128 * t203 - t146 * t160 + t267 * t80;
t45 = -t100 * t323 - t127 * t203 + t146 * t158 - t267 * t79;
t44 = (t349 - t380) * t317 + t316 * t326;
t43 = t314 * t326 + t317 * t379 + t350;
t37 = t127 * t160 - t128 * t158 + t323 * t80 - t324 * t79;
t36 = (t314 * t380 + t316 * t379) * t315 + t329;
t35 = -t221 * t364 - t256 * t367 + t276 * t365 + t286 * t379;
t34 = t219 * t364 + t254 * t367 - t276 * t366 - t286 * t380;
t31 = t317 * t86 + (t314 * t59 - t316 * t58) * t315;
t30 = t254 * t58 + t256 * t59 + t286 * t86;
t29 = t267 * t86 - t323 * t58 - t324 * t59;
t28 = -t219 * t365 + t221 * t366 - t254 * t379 + t256 * t380;
t27 = t317 * t71 + (t314 * t55 - t316 * t54) * t315;
t26 = t317 * t70 + (t314 * t53 - t316 * t52) * t315;
t25 = t254 * t54 + t256 * t55 + t286 * t71;
t24 = t254 * t52 + t256 * t53 + t286 * t70;
t23 = t267 * t71 - t323 * t54 - t324 * t55;
t22 = t267 * t70 - t323 * t52 - t324 * t53;
t15 = t317 * t61 + (t314 * t48 - t316 * t47) * t315;
t14 = t317 * t51 + (t314 * t42 - t316 * t41) * t315;
t13 = t317 * t50 + (t314 * t40 - t316 * t39) * t315;
t12 = t219 * t87 + t221 * t88 + t254 * t47 + t256 * t48 + t286 * t61 + t373;
t9 = t317 * t38 + (-t20 * t316 + t21 * t314) * t315;
t8 = t317 * t33 + (-t18 * t316 + t19 * t314) * t315;
t7 = t317 * t32 + (-t16 * t316 + t17 * t314) * t315;
t6 = t20 * t254 + t21 * t256 + t219 * t58 + t221 * t59 + t276 * t86 + t286 * t38;
t2 = t158 * t54 + t160 * t55 - t18 * t323 - t19 * t324 + t203 * t71 + t267 * t33;
t1 = t158 * t52 - t16 * t323 + t160 * t53 - t17 * t324 + t203 * t70 + t267 * t32;
t10 = [0; m(4) * t129 + m(5) * t93 + m(6) * t57 + m(7) * t36 + t406; (t36 * t49 + t43 * t62 + t44 * t63) * t395 + (t57 * t72 + t66 * t89 + t67 * t90) * t396 + 0.2e1 * m(5) * (t101 * t93 + t103 * t125 + t104 * t126) + 0.2e1 * m(4) * (t129 * t132 + t133 * t153 + t134 * t154) + (t232 * t314 + t234 * t316) * t406 * t409 + (t9 + t15 + ((t172 * t286 - t174 * t287 + t180 * t276 - t182 * t277 + (t414 * t322 + t412 * t319 + (t319 * t418 + t322 * t416) * qJD(2)) * t315) * t316 + (-t173 * t286 + t175 * t287 - t181 * t276 + t183 * t277 + (-t413 * t322 + t411 * t319 + (t319 * t417 + t322 * t415) * qJD(2)) * t315) * t314) * t315 + (-t212 * t286 + t213 * t287 - t215 * t276 + t216 * t277 + (t232 * t250 + t234 * t253) * t410 + (-t402 * t322 + t401 * t319 + (t319 * t405 + t322 * t404) * qJD(2) + t399 * t316 + t400 * t314) * t315 + t398 * t317) * t317) * t317 + (t8 + t14 + (-t173 * t256 + t175 * t257 - t181 * t221 + t183 * t222 - t290 * t415 + t291 * t417 + t302 * t413 + t303 * t411 + t372 * t400 + t397) * t372 + (-t212 * t256 + t213 * t257 - t215 * t221 + t216 * t222 + (-t234 * t299 - t253 * t283) * t410 + t398 * t372 + t401 * t303 + t402 * t302 + t405 * t291 - t404 * t290) * t317) * t372 + (-t7 - t13 + (-t172 * t254 + t174 * t255 - t180 * t219 + t182 * t220 - t288 * t416 - t289 * t418 - t301 * t412 - t327 * t414 + t371 * t399 + t397) * t371 + (t212 * t254 - t213 * t255 + t215 * t219 - t216 * t220 + (t232 * t299 + t250 * t283) * t410 + t398 * t371 - t401 * t301 + t402 * t327 - t405 * t289 - t404 * t288) * t317 + (t173 * t254 - t175 * t255 + t181 * t219 - t183 * t220 + t172 * t256 - t174 * t257 + t180 * t221 - t182 * t222 + t400 * t371 + t399 * t372 + t412 * t303 - t414 * t302 - t411 * t301 + t413 * t327 + t418 * t291 - t416 * t290 - t417 * t289 - t415 * t288) * t372) * t371; (m(4) + m(5) + m(6) + m(7)) * t315 * t354; m(7) * (t289 * t62 + t291 * t63 - t327 * t43 + t302 * t44 + (-t322 * t36 + t354 * t49) * t315) + m(6) * (t289 * t89 + t291 * t90 - t327 * t66 + t302 * t67 + (-t322 * t57 + t354 * t72) * t315) + m(5) * (-t327 * t103 + t302 * t104 + t289 * t125 + t291 * t126 + (t101 * t354 - t322 * t93) * t315) + m(4) * (-t327 * t133 + t302 * t134 + t289 * t153 + t291 * t154 + (-t129 * t322 + t132 * t354) * t315); 0.4e1 * (m(4) / 0.2e1 + t337) * (-t315 ^ 2 * t319 * t353 - t289 * t327 + t302 * t291); 0; m(7) * (-t317 * t36 + (-t314 * t44 + t316 * t43) * t315) + m(6) * (-t317 * t57 + (-t314 * t67 + t316 * t66) * t315) + m(5) * (-t317 * t93 + (t103 * t316 - t104 * t314) * t315); t337 * (t289 * t316 - t291 * t314 - t341) * t409; 0; m(6) * t56 + m(7) * t28; t31 * t384 + t26 * t387 + t27 * t386 + (t9 / 0.2e1 + t15 / 0.2e1) * t286 + (t8 / 0.2e1 + t14 / 0.2e1) * t256 + (t7 / 0.2e1 + t13 / 0.2e1) * t254 + m(7) * (t28 * t49 + t34 * t63 + t35 * t62 + t36 * t60 + t43 * t69 + t44 * t68) + m(6) * (t105 * t67 + t106 * t66 + t56 * t72 + t57 * t96 + t64 * t90 + t65 * t89) + (t6 / 0.2e1 + t373 / 0.2e1 + t12 / 0.2e1 + t95 * t386 + t94 * t387) * t317 + ((t314 * t88 - t316 * t87) * t384 + (t314 * t85 - t316 * t84) * t386 + (t314 * t83 - t316 * t82) * t387 + t407 * t383 + t408 * t382) * t315; m(6) * (t105 * t291 + t106 * t289 - t65 * t327 + t64 * t302 + (-t322 * t56 + t354 * t96) * t315) + m(7) * (t69 * t289 + t68 * t291 - t35 * t327 + t34 * t302 + (-t28 * t322 + t354 * t60) * t315); m(6) * (-t317 * t56 + (-t314 * t64 + t316 * t65) * t315) + m(7) * (-t28 * t317 + (-t314 * t34 + t316 * t35) * t315); t276 * t30 + (t12 + t6 + t373) * t286 + (t276 * t88 + t407) * t256 + (t276 * t87 + t408) * t254 + (t254 * t84 + t256 * t85 + t286 * t95 + t25) * t221 + (t254 * t82 + t256 * t83 + t286 * t94 + t24) * t219 + (t28 * t60 + t34 * t68 + t35 * t69) * t395 + (t105 * t64 + t106 * t65 + t56 * t96) * t396; m(7) * t37; m(7) * (t36 * t81 + t37 * t49 + t43 * t92 + t44 * t91 + t45 * t63 + t46 * t62) + t31 * t390 + t9 * t385 + t26 * t392 + t7 * t389 + t317 * t393 + t27 * t391 + t8 * t388 + (t1 * t382 + t2 * t383) * t315; m(7) * (t92 * t289 + t91 * t291 - t46 * t327 + t45 * t302 + (-t322 * t37 + t354 * t81) * t315); m(7) * (-t317 * t37 + (-t314 * t45 + t316 * t46) * t315); t29 * t384 + t286 * t393 + t22 * t387 + t254 * t1 / 0.2e1 + t24 * t392 + t3 * t389 + t25 * t391 + t4 * t388 + m(7) * (t28 * t81 + t34 * t91 + t35 * t92 + t37 * t60 + t45 * t68 + t46 * t69) + t23 * t386 + t256 * t2 / 0.2e1 + t30 * t390 + t6 * t385; t160 * t23 - t324 * t2 + t158 * t22 - t323 * t1 + t203 * t29 + t267 * t5 + (t37 * t81 + t45 * t91 + t46 * t92) * t395;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
