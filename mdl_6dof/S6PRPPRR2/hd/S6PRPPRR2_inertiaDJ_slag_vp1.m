% Calculate time derivative of joint inertia matrix for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:50
% EndTime: 2019-03-08 19:18:08
% DurationCPUTime: 16.29s
% Computational Cost: add. (50306->834), mult. (143570->1200), div. (0->0), fcn. (169942->12), ass. (0->341)
t430 = Icges(4,1) + Icges(5,2);
t429 = Icges(4,4) + Icges(5,6);
t427 = Icges(5,4) - Icges(4,5);
t426 = Icges(5,5) - Icges(4,6);
t428 = -Icges(4,2) - Icges(5,3);
t320 = sin(qJ(2));
t322 = cos(qJ(2));
t379 = sin(pkin(11));
t380 = cos(pkin(11));
t310 = -t320 * t379 + t322 * t380;
t305 = t310 * qJD(2);
t315 = sin(pkin(10));
t317 = cos(pkin(10));
t327 = t320 * t380 + t322 * t379;
t306 = t327 * qJD(2);
t318 = cos(pkin(6));
t323 = t318 * t306;
t242 = -t315 * t305 - t317 * t323;
t326 = t318 * t310;
t324 = qJD(2) * t326;
t243 = -t315 * t306 + t317 * t324;
t419 = t242 * t428 - t243 * t429;
t244 = -t305 * t317 + t315 * t323;
t245 = -t306 * t317 - t315 * t324;
t418 = t244 * t428 - t245 * t429;
t417 = -t429 * t242 - t243 * t430;
t416 = t429 * t244 + t245 * t430;
t265 = -t315 * t327 + t317 * t326;
t325 = t318 * t327;
t266 = t315 * t310 + t317 * t325;
t316 = sin(pkin(6));
t374 = t316 * t317;
t415 = t265 * t428 - t266 * t429 - t374 * t426;
t267 = -t315 * t326 - t317 * t327;
t268 = t310 * t317 - t315 * t325;
t375 = t315 * t316;
t414 = t267 * t428 - t268 * t429 + t375 * t426;
t413 = -t429 * t265 - t266 * t430 - t427 * t374;
t412 = t429 * t267 + t268 * t430 - t427 * t375;
t291 = t327 * t316;
t285 = qJD(2) * t291;
t290 = t310 * t316;
t286 = t290 * qJD(2);
t411 = t285 * t428 + t286 * t429;
t410 = t429 * t285 - t286 * t430;
t409 = -t290 * t428 + t291 * t429 - t318 * t426;
t408 = -t429 * t290 - t291 * t430 + t427 * t318;
t425 = 2 * m(3);
t424 = 2 * m(4);
t342 = m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1;
t423 = 0.2e1 * t342;
t388 = cos(qJ(5));
t348 = t316 * t388;
t387 = sin(qJ(5));
t218 = -t267 * t387 + t315 * t348;
t319 = sin(qJ(6));
t321 = cos(qJ(6));
t166 = -t218 * t319 + t268 * t321;
t167 = t218 * t321 + t268 * t319;
t347 = t316 * t387;
t217 = t267 * t388 + t315 * t347;
t118 = Icges(7,5) * t167 + Icges(7,6) * t166 + Icges(7,3) * t217;
t120 = Icges(7,4) * t167 + Icges(7,2) * t166 + Icges(7,6) * t217;
t122 = Icges(7,1) * t167 + Icges(7,4) * t166 + Icges(7,5) * t217;
t162 = -qJD(5) * t217 - t244 * t387;
t124 = -qJD(6) * t167 - t162 * t319 + t245 * t321;
t125 = qJD(6) * t166 + t162 * t321 + t245 * t319;
t163 = qJD(5) * t218 + t244 * t388;
t73 = Icges(7,5) * t125 + Icges(7,6) * t124 + Icges(7,3) * t163;
t75 = Icges(7,4) * t125 + Icges(7,2) * t124 + Icges(7,6) * t163;
t77 = Icges(7,1) * t125 + Icges(7,4) * t124 + Icges(7,5) * t163;
t16 = t118 * t163 + t120 * t124 + t122 * t125 + t166 * t75 + t167 * t77 + t217 * t73;
t407 = t265 * t387 + t317 * t348;
t168 = t266 * t321 + t319 * t407;
t169 = t266 * t319 - t321 * t407;
t219 = -t265 * t388 + t317 * t347;
t119 = Icges(7,5) * t169 + Icges(7,6) * t168 - Icges(7,3) * t219;
t121 = Icges(7,4) * t169 + Icges(7,2) * t168 - Icges(7,6) * t219;
t123 = Icges(7,1) * t169 + Icges(7,4) * t168 - Icges(7,5) * t219;
t164 = qJD(5) * t219 - t242 * t387;
t126 = -qJD(6) * t169 - t164 * t319 + t243 * t321;
t127 = qJD(6) * t168 + t164 * t321 + t243 * t319;
t165 = qJD(5) * t407 - t242 * t388;
t74 = Icges(7,5) * t127 + Icges(7,6) * t126 - Icges(7,3) * t165;
t76 = Icges(7,4) * t127 + Icges(7,2) * t126 - Icges(7,6) * t165;
t78 = Icges(7,1) * t127 + Icges(7,4) * t126 - Icges(7,5) * t165;
t17 = t119 * t163 + t121 * t124 + t123 * t125 + t166 * t76 + t167 * t78 + t217 * t74;
t272 = t290 * t388 + t318 * t387;
t210 = -qJD(5) * t272 + t285 * t387;
t273 = -t290 * t387 + t318 * t388;
t215 = t273 * t321 + t291 * t319;
t149 = -qJD(6) * t215 - t210 * t319 + t286 * t321;
t214 = -t273 * t319 + t291 * t321;
t150 = qJD(6) * t214 + t210 * t321 + t286 * t319;
t211 = qJD(5) * t273 - t285 * t388;
t102 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t211;
t103 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t211;
t104 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t211;
t143 = Icges(7,5) * t215 + Icges(7,6) * t214 + Icges(7,3) * t272;
t144 = Icges(7,4) * t215 + Icges(7,2) * t214 + Icges(7,6) * t272;
t145 = Icges(7,1) * t215 + Icges(7,4) * t214 + Icges(7,5) * t272;
t32 = t102 * t217 + t103 * t166 + t104 * t167 + t124 * t144 + t125 * t145 + t143 * t163;
t52 = t118 * t217 + t120 * t166 + t122 * t167;
t53 = t119 * t217 + t121 * t166 + t123 * t167;
t71 = t143 * t217 + t144 * t166 + t145 * t167;
t3 = t16 * t268 + t17 * t266 + t243 * t53 + t245 * t52 + t286 * t71 + t291 * t32;
t110 = Icges(6,5) * t162 - Icges(6,6) * t163 + Icges(6,3) * t245;
t112 = Icges(6,4) * t162 - Icges(6,2) * t163 + Icges(6,6) * t245;
t114 = Icges(6,1) * t162 - Icges(6,4) * t163 + Icges(6,5) * t245;
t135 = Icges(6,5) * t218 - Icges(6,6) * t217 + Icges(6,3) * t268;
t137 = Icges(6,4) * t218 - Icges(6,2) * t217 + Icges(6,6) * t268;
t139 = Icges(6,1) * t218 - Icges(6,4) * t217 + Icges(6,5) * t268;
t41 = t110 * t268 - t112 * t217 + t114 * t218 + t135 * t245 - t137 * t163 + t139 * t162;
t111 = Icges(6,5) * t164 + Icges(6,6) * t165 + Icges(6,3) * t243;
t113 = Icges(6,4) * t164 + Icges(6,2) * t165 + Icges(6,6) * t243;
t115 = Icges(6,1) * t164 + Icges(6,4) * t165 + Icges(6,5) * t243;
t136 = -Icges(6,5) * t407 + Icges(6,6) * t219 + Icges(6,3) * t266;
t138 = -Icges(6,4) * t407 + Icges(6,2) * t219 + Icges(6,6) * t266;
t140 = -Icges(6,1) * t407 + Icges(6,4) * t219 + Icges(6,5) * t266;
t42 = t111 * t268 - t113 * t217 + t115 * t218 + t136 * t245 - t138 * t163 + t140 * t162;
t151 = Icges(6,5) * t210 - Icges(6,6) * t211 + Icges(6,3) * t286;
t152 = Icges(6,4) * t210 - Icges(6,2) * t211 + Icges(6,6) * t286;
t153 = Icges(6,1) * t210 - Icges(6,4) * t211 + Icges(6,5) * t286;
t204 = Icges(6,5) * t273 - Icges(6,6) * t272 + Icges(6,3) * t291;
t205 = Icges(6,4) * t273 - Icges(6,2) * t272 + Icges(6,6) * t291;
t206 = Icges(6,1) * t273 - Icges(6,4) * t272 + Icges(6,5) * t291;
t50 = t151 * t268 - t152 * t217 + t153 * t218 + t162 * t206 - t163 * t205 + t204 * t245;
t82 = t135 * t268 - t137 * t217 + t139 * t218;
t83 = t136 * t268 - t138 * t217 + t140 * t218;
t97 = t204 * t268 - t205 * t217 + t206 * t218;
t422 = t243 * t83 + t245 * t82 + t266 * t42 + t268 * t41 + t286 * t97 + t291 * t50 + t3;
t18 = -t118 * t165 + t120 * t126 + t122 * t127 + t168 * t75 + t169 * t77 - t219 * t73;
t19 = -t119 * t165 + t121 * t126 + t123 * t127 + t168 * t76 + t169 * t78 - t219 * t74;
t33 = -t102 * t219 + t103 * t168 + t104 * t169 + t126 * t144 + t127 * t145 - t143 * t165;
t54 = -t118 * t219 + t120 * t168 + t122 * t169;
t55 = -t119 * t219 + t121 * t168 + t123 * t169;
t72 = -t143 * t219 + t144 * t168 + t145 * t169;
t4 = t18 * t268 + t19 * t266 + t243 * t55 + t245 * t54 + t286 * t72 + t291 * t33;
t43 = t110 * t266 + t112 * t219 - t114 * t407 + t135 * t243 + t137 * t165 + t139 * t164;
t44 = t111 * t266 + t113 * t219 - t115 * t407 + t136 * t243 + t138 * t165 + t140 * t164;
t51 = t151 * t266 + t152 * t219 - t153 * t407 + t164 * t206 + t165 * t205 + t204 * t243;
t84 = t135 * t266 + t137 * t219 - t139 * t407;
t85 = t136 * t266 + t138 * t219 - t140 * t407;
t98 = t204 * t266 + t205 * t219 - t206 * t407;
t421 = t243 * t85 + t245 * t84 + t266 * t44 + t268 * t43 + t286 * t98 + t291 * t51 + t4;
t372 = t318 * t322;
t301 = -t315 * t320 + t317 * t372;
t292 = t301 * qJD(2);
t373 = t318 * t320;
t302 = t315 * t322 + t317 * t373;
t293 = t302 * qJD(2);
t263 = rSges(3,1) * t292 - rSges(3,2) * t293;
t303 = -t315 * t372 - t317 * t320;
t294 = t303 * qJD(2);
t329 = t315 * t373 - t317 * t322;
t295 = t329 * qJD(2);
t264 = rSges(3,1) * t294 + rSges(3,2) * t295;
t420 = m(3) * (t263 * t315 + t264 * t317) * t316;
t406 = Icges(3,5) * t294 + Icges(3,6) * t295 - t426 * t244 - t245 * t427;
t357 = qJD(2) * t316;
t405 = (Icges(3,5) * t322 - Icges(3,6) * t320) * t357 - t427 * t286 + t426 * t285;
t404 = -Icges(3,5) * t292 + Icges(3,6) * t293 + t426 * t242 + t243 * t427;
t289 = t318 * rSges(3,3) + (rSges(3,1) * t320 + rSges(3,2) * t322) * t316;
t299 = (rSges(3,1) * t322 - rSges(3,2) * t320) * t357;
t403 = t289 * t299 * t425;
t402 = 0.2e1 * m(6);
t401 = 0.2e1 * m(7);
t20 = t118 * t211 + t120 * t149 + t122 * t150 + t214 * t75 + t215 * t77 + t272 * t73;
t21 = t119 * t211 + t121 * t149 + t123 * t150 + t214 * t76 + t215 * t78 + t272 * t74;
t38 = t102 * t272 + t103 * t214 + t104 * t215 + t143 * t211 + t144 * t149 + t145 * t150;
t58 = t118 * t272 + t120 * t214 + t122 * t215;
t59 = t119 * t272 + t121 * t214 + t123 * t215;
t86 = t143 * t272 + t144 * t214 + t145 * t215;
t5 = t163 * t58 - t165 * t59 + t20 * t217 - t21 * t219 + t211 * t86 + t272 * t38;
t400 = t5 / 0.2e1;
t399 = t163 / 0.2e1;
t398 = -t165 / 0.2e1;
t397 = t211 / 0.2e1;
t396 = t217 / 0.2e1;
t395 = -t219 / 0.2e1;
t394 = t243 / 0.2e1;
t393 = t245 / 0.2e1;
t392 = t272 / 0.2e1;
t391 = t286 / 0.2e1;
t390 = t315 / 0.2e1;
t389 = -t317 / 0.2e1;
t386 = pkin(8) * t243;
t385 = pkin(8) * t245;
t384 = pkin(2) * t322;
t383 = pkin(2) * qJD(2);
t79 = rSges(7,1) * t125 + rSges(7,2) * t124 + rSges(7,3) * t163;
t382 = pkin(5) * t162 + pkin(9) * t163 + t79;
t80 = rSges(7,1) * t127 + rSges(7,2) * t126 - rSges(7,3) * t165;
t381 = pkin(5) * t164 - pkin(9) * t165 + t80;
t378 = Icges(3,4) * t320;
t377 = Icges(3,4) * t322;
t376 = (t204 * t291 - t205 * t272 + t206 * t273) * t286;
t105 = rSges(7,1) * t150 + rSges(7,2) * t149 + rSges(7,3) * t211;
t371 = pkin(5) * t210 + pkin(9) * t211 + t105;
t128 = rSges(7,1) * t167 + rSges(7,2) * t166 + rSges(7,3) * t217;
t370 = pkin(5) * t218 + pkin(9) * t217 + t128;
t129 = rSges(7,1) * t169 + rSges(7,2) * t168 - rSges(7,3) * t219;
t369 = -pkin(5) * t407 - pkin(9) * t219 + t129;
t146 = rSges(7,1) * t215 + rSges(7,2) * t214 + rSges(7,3) * t272;
t368 = pkin(5) * t273 + pkin(9) * t272 + t146;
t157 = pkin(3) * t245 - qJ(4) * t244 - qJD(4) * t267;
t354 = t322 * t383;
t308 = -t316 * qJD(3) + t318 * t354;
t355 = t320 * t383;
t277 = -t308 * t315 - t317 * t355;
t271 = t318 * t277;
t367 = t318 * t157 + t271;
t156 = pkin(3) * t243 - qJ(4) * t242 - qJD(4) * t265;
t276 = t308 * t317 - t315 * t355;
t366 = -t156 - t276;
t203 = pkin(3) * t268 - qJ(4) * t267;
t346 = pkin(2) * t373 - qJ(3) * t316;
t255 = -t315 * t346 + t317 * t384;
t241 = t318 * t255;
t365 = t318 * t203 + t241;
t202 = pkin(3) * t266 - qJ(4) * t265;
t254 = t315 * t384 + t317 * t346;
t364 = -t202 - t254;
t307 = qJD(3) * t318 + t316 * t354;
t363 = -pkin(3) * t286 - qJ(4) * t285 + qJD(4) * t290 - t307;
t356 = pkin(8) * t286 * t316;
t362 = -t317 * t356 - t318 * t386;
t361 = t254 * t375 + t255 * t374;
t311 = pkin(2) * t316 * t320 + qJ(3) * t318;
t360 = -rSges(4,1) * t291 - rSges(4,2) * t290 - rSges(4,3) * t318 - t311;
t359 = -pkin(3) * t291 + qJ(4) * t290 - t311;
t358 = t276 * t375 + t277 * t374;
t353 = -t315 * t356 + t367;
t222 = pkin(4) * t375 + pkin(8) * t268;
t352 = t318 * t222 + t365;
t223 = -pkin(4) * t374 + pkin(8) * t266;
t351 = -t223 + t364;
t350 = -pkin(4) * t318 - pkin(8) * t291 + t359;
t343 = (-rSges(4,1) * t286 + rSges(4,2) * t285 - t307) * t316;
t341 = t156 * t375 + t157 * t374 + t358;
t340 = t202 * t375 + t203 * t374 + t361;
t338 = t382 + t385;
t154 = rSges(6,1) * t210 - rSges(6,2) * t211 + rSges(6,3) * t286;
t335 = (-t154 + t363) * t316;
t334 = (rSges(5,2) * t286 - rSges(5,3) * t285 + t363) * t316;
t333 = (-rSges(5,1) * t318 + rSges(5,2) * t291 + rSges(5,3) * t290 + t359) * t316;
t332 = (t363 - t371) * t316;
t207 = rSges(6,1) * t273 - rSges(6,2) * t272 + rSges(6,3) * t291;
t331 = (-t207 + t350) * t316;
t330 = t222 * t374 + t223 * t375 + t340;
t328 = (t350 - t368) * t316;
t298 = (Icges(3,1) * t322 - t378) * t357;
t297 = (-Icges(3,2) * t320 + t377) * t357;
t288 = Icges(3,5) * t318 + (Icges(3,1) * t320 + t377) * t316;
t287 = Icges(3,6) * t318 + (Icges(3,2) * t322 + t378) * t316;
t262 = Icges(3,1) * t294 + Icges(3,4) * t295;
t261 = Icges(3,1) * t292 - Icges(3,4) * t293;
t260 = Icges(3,4) * t294 + Icges(3,2) * t295;
t259 = Icges(3,4) * t292 - Icges(3,2) * t293;
t252 = -rSges(3,1) * t329 + rSges(3,2) * t303 + rSges(3,3) * t375;
t251 = rSges(3,1) * t302 + rSges(3,2) * t301 - rSges(3,3) * t374;
t250 = -Icges(3,1) * t329 + Icges(3,4) * t303 + Icges(3,5) * t375;
t249 = Icges(3,1) * t302 + Icges(3,4) * t301 - Icges(3,5) * t374;
t248 = -Icges(3,4) * t329 + Icges(3,2) * t303 + Icges(3,6) * t375;
t247 = Icges(3,4) * t302 + Icges(3,2) * t301 - Icges(3,6) * t374;
t198 = rSges(4,1) * t268 + rSges(4,2) * t267 + rSges(4,3) * t375;
t197 = rSges(4,1) * t266 + rSges(4,2) * t265 - rSges(4,3) * t374;
t196 = -rSges(5,1) * t374 - rSges(5,2) * t266 - rSges(5,3) * t265;
t195 = rSges(5,1) * t375 - rSges(5,2) * t268 - rSges(5,3) * t267;
t186 = -rSges(5,2) * t245 - rSges(5,3) * t244;
t185 = rSges(4,1) * t245 + rSges(4,2) * t244;
t184 = -rSges(5,2) * t243 - rSges(5,3) * t242;
t183 = rSges(4,1) * t243 + rSges(4,2) * t242;
t142 = -rSges(6,1) * t407 + rSges(6,2) * t219 + rSges(6,3) * t266;
t141 = rSges(6,1) * t218 - rSges(6,2) * t217 + rSges(6,3) * t268;
t134 = (-t183 - t276) * t318 + t317 * t343;
t133 = t185 * t318 + t315 * t343 + t271;
t117 = rSges(6,1) * t164 + rSges(6,2) * t165 + rSges(6,3) * t243;
t116 = rSges(6,1) * t162 - rSges(6,2) * t163 + rSges(6,3) * t245;
t109 = (t183 * t315 + t185 * t317) * t316 + t358;
t108 = t141 * t291 - t207 * t268;
t107 = -t142 * t291 + t207 * t266;
t101 = (-t196 + t364) * t318 + t317 * t333;
t100 = t195 * t318 + t315 * t333 + t365;
t99 = -t141 * t266 + t142 * t268;
t96 = (-t184 + t366) * t318 + t317 * t334;
t95 = t186 * t318 + t315 * t334 + t367;
t94 = (t195 * t317 + t196 * t315) * t316 + t340;
t93 = -t129 * t272 - t146 * t219;
t92 = t128 * t272 - t146 * t217;
t91 = (-t142 + t351) * t318 + t317 * t331;
t90 = t141 * t318 + t315 * t331 + t352;
t89 = t136 * t291 - t138 * t272 + t140 * t273;
t88 = t135 * t291 - t137 * t272 + t139 * t273;
t87 = (t184 * t315 + t186 * t317) * t316 + t341;
t81 = t128 * t219 + t129 * t217;
t70 = -t268 * t368 + t291 * t370;
t69 = t266 * t368 - t291 * t369;
t68 = (t141 * t317 + t142 * t315) * t316 + t330;
t67 = (-t117 + t366) * t318 + t317 * t335 + t362;
t66 = (t116 + t385) * t318 + t315 * t335 + t353;
t65 = t116 * t291 + t141 * t286 - t154 * t268 - t207 * t245;
t64 = -t117 * t291 - t142 * t286 + t154 * t266 + t207 * t243;
t63 = (t351 - t369) * t318 + t317 * t328;
t62 = t315 * t328 + t318 * t370 + t352;
t61 = t151 * t291 - t152 * t272 + t153 * t273 + t204 * t286 - t205 * t211 + t206 * t210;
t60 = -t266 * t370 + t268 * t369;
t57 = (t116 * t317 + t117 * t315 + (t243 * t315 + t245 * t317) * pkin(8)) * t316 + t341;
t56 = -t116 * t266 + t117 * t268 - t141 * t243 + t142 * t245;
t49 = (t315 * t369 + t317 * t370) * t316 + t330;
t48 = t111 * t291 - t113 * t272 + t115 * t273 + t136 * t286 - t138 * t211 + t140 * t210;
t47 = t110 * t291 - t112 * t272 + t114 * t273 + t135 * t286 - t137 * t211 + t139 * t210;
t46 = -t105 * t219 - t129 * t211 - t146 * t165 - t272 * t80;
t45 = -t105 * t217 + t128 * t211 - t146 * t163 + t272 * t79;
t40 = (t366 - t381) * t318 + t317 * t332 + t362;
t39 = t315 * t332 + t318 * t338 + t353;
t37 = t128 * t165 + t129 * t163 + t217 * t80 + t219 * t79;
t36 = -t245 * t368 - t268 * t371 + t286 * t370 + t291 * t382;
t35 = t243 * t368 + t266 * t371 - t286 * t369 - t291 * t381;
t34 = (t338 * t317 + (t381 + t386) * t315) * t316 + t341;
t31 = t318 * t86 + (t315 * t58 - t317 * t59) * t316;
t30 = t266 * t59 + t268 * t58 + t291 * t86;
t29 = t217 * t58 - t219 * t59 + t272 * t86;
t28 = -t243 * t370 + t245 * t369 - t266 * t382 + t268 * t381;
t27 = t318 * t72 + (t315 * t54 - t317 * t55) * t316;
t26 = t318 * t71 + (t315 * t52 - t317 * t53) * t316;
t25 = t266 * t55 + t268 * t54 + t291 * t72;
t24 = t266 * t53 + t268 * t52 + t291 * t71;
t23 = t217 * t54 - t219 * t55 + t272 * t72;
t22 = t217 * t52 - t219 * t53 + t272 * t71;
t15 = t318 * t61 + (t315 * t47 - t317 * t48) * t316;
t14 = t318 * t51 + (t315 * t43 - t317 * t44) * t316;
t13 = t318 * t50 + (t315 * t41 - t317 * t42) * t316;
t12 = t243 * t89 + t245 * t88 + t266 * t48 + t268 * t47 + t291 * t61 + t376;
t9 = t318 * t38 + (t20 * t315 - t21 * t317) * t316;
t8 = t318 * t33 + (t18 * t315 - t19 * t317) * t316;
t7 = t318 * t32 + (t16 * t315 - t17 * t317) * t316;
t6 = t20 * t268 + t21 * t266 + t243 * t59 + t245 * t58 + t286 * t86 + t291 * t38;
t2 = t163 * t54 - t165 * t55 + t18 * t217 - t19 * t219 + t211 * t72 + t272 * t33;
t1 = t16 * t217 + t163 * t52 - t165 * t53 - t17 * t219 + t211 * t71 + t272 * t32;
t10 = [0; m(4) * t109 + m(5) * t87 + m(6) * t57 + m(7) * t34 + t420; (t34 * t49 + t39 * t62 + t40 * t63) * t401 + (t57 * t68 + t66 * t90 + t67 * t91) * t402 + 0.2e1 * m(5) * (t100 * t95 + t101 * t96 + t87 * t94) + (t241 * t133 + t361 * t109 + ((t198 * t109 + t134 * t360) * t317 + (t197 * t109 + t133 * t360) * t315) * t316) * t424 + 0.2e1 * (t251 * t315 + t252 * t317) * t316 * t420 + (t9 + t15 + ((-t197 - t254) * t134 + t198 * t133) * t424 + ((t251 * t263 + t252 * t264) * t425 + t405 * t318 - t410 * t291 + t411 * t290 - t408 * t286 - t409 * t285) * t318 + ((t297 * t322 + t298 * t320 + (-t287 * t320 + t288 * t322) * qJD(2)) * t318 + (t404 * t318 + (-t259 * t322 - t261 * t320 - (-t247 * t320 + t249 * t322) * qJD(2)) * t316 + t417 * t291 + t419 * t290 + t413 * t286 - t415 * t285) * t317 + (t406 * t318 + (t260 * t322 + t262 * t320 + (-t248 * t320 + t250 * t322) * qJD(2)) * t316 + t416 * t291 - t418 * t290 + t412 * t286 + t414 * t285) * t315) * t316) * t318 + (t7 + t13 + (-t244 * t414 + t245 * t412 + t248 * t295 + t250 * t294 + t260 * t303 - t262 * t329 - t267 * t418 + t268 * t416 + t375 * t406 + t403) * t375 + (t287 * t295 + t288 * t294 + t297 * t303 - t298 * t329 + (-t252 * t299 - t264 * t289) * t425 + t405 * t375 - t410 * t268 + t411 * t267 - t408 * t245 + t409 * t244) * t318) * t375 + (-t8 - t14 + (-t415 * t242 - t413 * t243 - t247 * t293 + t249 * t292 + t259 * t301 + t261 * t302 - t265 * t419 - t417 * t266 + t374 * t404 + t403) * t374 + (t287 * t293 - t288 * t292 - t297 * t301 - t298 * t302 + (t251 * t299 + t263 * t289) * t425 + t405 * t374 + t410 * t266 - t411 * t265 + t408 * t243 - t409 * t242) * t318 + (t406 * t374 + t404 * t375 + t248 * t293 - t250 * t292 - t260 * t301 - t262 * t302 - t247 * t295 - t249 * t294 - t259 * t303 + t261 * t329 + t417 * t268 + t419 * t267 - t416 * t266 + t418 * t265 + t413 * t245 + t415 * t244 - t412 * t243 + t414 * t242) * t375) * t374; 0; m(7) * (t318 * t34 + (t315 * t40 - t317 * t39) * t316) + m(6) * (t318 * t57 + (t315 * t67 - t317 * t66) * t316) + m(5) * (t318 * t87 + (t315 * t96 - t317 * t95) * t316) + m(4) * (t109 * t318 + (-t133 * t317 + t134 * t315) * t316); 0; t285 * t423; m(7) * (-t242 * t62 - t244 * t63 - t265 * t39 - t267 * t40 + t285 * t49 - t290 * t34) + m(6) * (-t242 * t90 - t244 * t91 - t265 * t66 - t267 * t67 + t285 * t68 - t290 * t57) + m(5) * (-t100 * t242 - t101 * t244 - t265 * t95 - t267 * t96 + t285 * t94 - t290 * t87); (t285 * t318 + (t242 * t317 - t244 * t315) * t316) * t423; 0.4e1 * t342 * (t242 * t265 + t244 * t267 - t285 * t290); m(6) * t56 + m(7) * t28; t27 * t394 + t26 * t393 + t31 * t391 + (t9 / 0.2e1 + t15 / 0.2e1) * t291 + (t7 / 0.2e1 + t13 / 0.2e1) * t268 + (t8 / 0.2e1 + t14 / 0.2e1) * t266 + m(7) * (t28 * t49 + t34 * t60 + t35 * t63 + t36 * t62 + t39 * t70 + t40 * t69) + m(6) * (t107 * t67 + t108 * t66 + t56 * t68 + t57 * t99 + t64 * t91 + t65 * t90) + (t6 / 0.2e1 + t98 * t394 + t97 * t393 + t376 / 0.2e1 + t12 / 0.2e1) * t318 + ((t315 * t84 - t317 * t85) * t394 + (t315 * t82 - t317 * t83) * t393 + (t315 * t88 - t317 * t89) * t391 + t422 * t390 + t421 * t389) * t316; m(6) * (t318 * t56 + (t315 * t64 - t317 * t65) * t316) + m(7) * (t28 * t318 + (t315 * t35 - t317 * t36) * t316); m(6) * (-t107 * t244 - t108 * t242 - t265 * t65 - t267 * t64 + t285 * t99 - t290 * t56) + m(7) * (-t242 * t70 - t244 * t69 - t265 * t36 - t267 * t35 - t28 * t290 + t285 * t60); t286 * t30 + (t12 + t6 + t376) * t291 + (t286 * t88 + t422) * t268 + (t286 * t89 + t421) * t266 + (t266 * t83 + t268 * t82 + t291 * t97 + t24) * t245 + (t266 * t85 + t268 * t84 + t291 * t98 + t25) * t243 + (t28 * t60 + t35 * t69 + t36 * t70) * t401 + (t107 * t64 + t108 * t65 + t56 * t99) * t402; m(7) * t37; m(7) * (t34 * t81 + t37 * t49 + t39 * t92 + t40 * t93 + t45 * t62 + t46 * t63) + t27 * t398 + t8 * t395 + t26 * t399 + t7 * t396 + t31 * t397 + t9 * t392 + t318 * t400 + (t1 * t390 + t2 * t389) * t316; m(7) * (t318 * t37 + (t315 * t46 - t317 * t45) * t316); m(7) * (-t242 * t92 - t244 * t93 - t265 * t45 - t267 * t46 + t285 * t81 - t290 * t37); t23 * t394 + t266 * t2 / 0.2e1 + t25 * t398 + t4 * t395 + t24 * t399 + t3 * t396 + t29 * t391 + t291 * t400 + t22 * t393 + t268 * t1 / 0.2e1 + m(7) * (t28 * t81 + t35 * t93 + t36 * t92 + t37 * t60 + t45 * t70 + t46 * t69) + t30 * t397 + t6 * t392; (t37 * t81 + t45 * t92 + t46 * t93) * t401 + t163 * t22 + t217 * t1 - t165 * t23 - t219 * t2 + t211 * t29 + t272 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
