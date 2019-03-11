% Calculate time derivative of joint inertia matrix for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:12
% EndTime: 2019-03-09 02:10:30
% DurationCPUTime: 11.96s
% Computational Cost: add. (10417->749), mult. (28224->1058), div. (0->0), fcn. (27018->6), ass. (0->338)
t247 = sin(qJ(4));
t250 = cos(qJ(4));
t246 = sin(qJ(5));
t249 = cos(qJ(5));
t288 = Icges(6,5) * t249 - Icges(6,6) * t246;
t173 = Icges(6,3) * t247 + t250 * t288;
t290 = Icges(7,4) * t249 + Icges(7,6) * t246;
t176 = Icges(7,2) * t247 + t250 * t290;
t437 = t176 + t173;
t344 = qJD(4) * t247;
t436 = qJD(1) * t250 / 0.2e1;
t251 = cos(qJ(1));
t360 = t251 * t249;
t248 = sin(qJ(1));
t365 = t248 * t246;
t196 = t247 * t365 - t360;
t364 = t248 * t249;
t197 = t246 * t251 + t247 * t364;
t431 = rSges(7,3) + qJ(6);
t433 = rSges(7,1) + pkin(5);
t435 = t431 * t196 + t197 * t433;
t311 = qJD(1) * t247 + qJD(5);
t274 = t311 * t251;
t312 = qJD(5) * t247 + qJD(1);
t342 = qJD(4) * t250;
t324 = t248 * t342;
t119 = t312 * t364 + (t274 + t324) * t246;
t275 = t312 * t246;
t120 = t249 * t274 + (t249 * t342 - t275) * t248;
t434 = -t196 * qJD(6) - t431 * t119 - t120 * t433;
t386 = Icges(6,4) * t249;
t291 = -Icges(6,2) * t246 + t386;
t177 = Icges(6,6) * t247 + t250 * t291;
t387 = Icges(6,4) * t246;
t294 = Icges(6,1) * t249 - t387;
t181 = Icges(6,5) * t247 + t250 * t294;
t383 = Icges(7,5) * t249;
t287 = Icges(7,3) * t246 + t383;
t172 = Icges(7,6) * t247 + t250 * t287;
t384 = Icges(7,5) * t246;
t293 = Icges(7,1) * t249 + t384;
t180 = Icges(7,4) * t247 + t250 * t293;
t279 = t172 * t246 + t180 * t249;
t432 = (-t177 * t246 + t181 * t249 + t279) * t250 + t437 * t247;
t340 = qJD(5) * t250;
t140 = (Icges(7,3) * t249 - t384) * t340 + (Icges(7,6) * t250 - t247 * t287) * qJD(4);
t141 = (-Icges(6,5) * t246 - Icges(6,6) * t249) * t340 + (Icges(6,3) * t250 - t247 * t288) * qJD(4);
t144 = (-Icges(7,4) * t246 + Icges(7,6) * t249) * t340 + (Icges(7,2) * t250 - t247 * t290) * qJD(4);
t148 = (-Icges(7,1) * t246 + t383) * t340 + (Icges(7,4) * t250 - t247 * t293) * qJD(4);
t149 = (-Icges(6,1) * t246 - t386) * t340 + (Icges(6,5) * t250 - t247 * t294) * qJD(4);
t321 = t249 * t340;
t322 = t246 * t340;
t325 = t249 * t344;
t326 = t246 * t344;
t368 = t246 * t250;
t430 = t140 * t368 + t172 * t321 + t177 * t326 - t180 * t322 - t181 * t325 + (t148 + t149) * t249 * t250 + t437 * t342 + (t141 + t144) * t247;
t341 = qJD(4) * t251;
t319 = t247 * t341;
t346 = qJD(1) * t248;
t256 = t250 * t346 + t319;
t367 = t247 * t251;
t198 = t246 * t367 + t364;
t334 = t247 * t360;
t199 = t334 - t365;
t361 = t250 * t251;
t128 = Icges(7,4) * t199 - Icges(7,2) * t361 + Icges(7,6) * t198;
t124 = Icges(7,5) * t199 - Icges(7,6) * t361 + Icges(7,3) * t198;
t132 = Icges(7,1) * t199 - Icges(7,4) * t361 + Icges(7,5) * t198;
t285 = t124 * t246 + t132 * t249;
t60 = t128 * t247 + t250 * t285;
t126 = Icges(6,5) * t199 - Icges(6,6) * t198 - Icges(6,3) * t361;
t130 = Icges(6,4) * t199 - Icges(6,2) * t198 - Icges(6,6) * t361;
t134 = Icges(6,1) * t199 - Icges(6,4) * t198 - Icges(6,5) * t361;
t283 = t130 * t246 - t134 * t249;
t62 = t126 * t247 - t250 * t283;
t396 = t60 + t62;
t363 = t248 * t250;
t127 = Icges(7,4) * t197 - Icges(7,2) * t363 + Icges(7,6) * t196;
t123 = Icges(7,5) * t197 - Icges(7,6) * t363 + Icges(7,3) * t196;
t131 = Icges(7,1) * t197 - Icges(7,4) * t363 + Icges(7,5) * t196;
t286 = t123 * t246 + t131 * t249;
t59 = t127 * t247 + t250 * t286;
t125 = Icges(6,5) * t197 - Icges(6,6) * t196 - Icges(6,3) * t363;
t129 = Icges(6,4) * t197 - Icges(6,2) * t196 - Icges(6,6) * t363;
t133 = Icges(6,1) * t197 - Icges(6,4) * t196 - Icges(6,5) * t363;
t284 = t129 * t246 - t133 * t249;
t61 = t125 * t247 - t250 * t284;
t397 = t59 + t61;
t429 = -t248 * t397 - t251 * t396;
t428 = t248 * t396 - t251 * t397;
t343 = qJD(4) * t248;
t320 = t247 * t343;
t345 = qJD(1) * t251;
t327 = t250 * t345;
t257 = t320 - t327;
t67 = Icges(7,5) * t120 + Icges(7,6) * t257 + Icges(7,3) * t119;
t71 = Icges(7,4) * t120 + Icges(7,2) * t257 + Icges(7,6) * t119;
t75 = Icges(7,1) * t120 + Icges(7,4) * t257 + Icges(7,5) * t119;
t19 = (-qJD(4) * t286 + t71) * t247 + (qJD(4) * t127 + t246 * t67 + t249 * t75 + (t123 * t249 - t131 * t246) * qJD(5)) * t250;
t69 = Icges(6,5) * t120 - Icges(6,6) * t119 + Icges(6,3) * t257;
t73 = Icges(6,4) * t120 - Icges(6,2) * t119 + Icges(6,6) * t257;
t77 = Icges(6,1) * t120 - Icges(6,4) * t119 + Icges(6,5) * t257;
t21 = (qJD(4) * t284 + t69) * t247 + (qJD(4) * t125 - t246 * t73 + t249 * t77 + (-t129 * t249 - t133 * t246) * qJD(5)) * t250;
t427 = t19 + t21;
t323 = t250 * t341;
t117 = -qJD(5) * t334 - t246 * t323 - t249 * t345 + t311 * t365;
t118 = -t251 * t275 + (-t248 * t311 + t323) * t249;
t66 = Icges(7,5) * t118 + Icges(7,6) * t256 - Icges(7,3) * t117;
t70 = Icges(7,4) * t118 + Icges(7,2) * t256 - Icges(7,6) * t117;
t74 = Icges(7,1) * t118 + Icges(7,4) * t256 - Icges(7,5) * t117;
t20 = (-qJD(4) * t285 + t70) * t247 + (qJD(4) * t128 + t246 * t66 + t249 * t74 + (t124 * t249 - t132 * t246) * qJD(5)) * t250;
t68 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t256;
t72 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t256;
t76 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t256;
t22 = (qJD(4) * t283 + t68) * t247 + (qJD(4) * t126 - t246 * t72 + t249 * t76 + (-t130 * t249 - t134 * t246) * qJD(5)) * t250;
t426 = -t20 - t22;
t51 = -t125 * t363 - t129 * t196 + t133 * t197;
t52 = -t126 * t363 - t130 * t196 + t134 * t197;
t300 = t248 * t51 + t251 * t52;
t49 = t123 * t196 - t127 * t363 + t131 * t197;
t50 = t124 * t196 - t128 * t363 + t132 * t197;
t302 = t248 * t49 + t251 * t50;
t86 = t172 * t196 - t176 * t363 + t180 * t197;
t87 = -t173 * t363 - t177 * t196 + t181 * t197;
t425 = (-t300 - t302) * t250 + (t86 + t87) * t247;
t55 = -t125 * t361 - t198 * t129 + t199 * t133;
t56 = -t126 * t361 - t198 * t130 + t199 * t134;
t296 = t248 * t55 + t251 * t56;
t53 = t198 * t123 - t127 * t361 + t199 * t131;
t54 = t198 * t124 - t128 * t361 + t199 * t132;
t298 = t248 * t53 + t251 * t54;
t88 = t198 * t172 - t176 * t361 + t199 * t180;
t89 = -t173 * t361 - t198 * t177 + t199 * t181;
t398 = (-t296 - t298) * t250 + (t88 + t89) * t247;
t304 = rSges(7,1) * t249 + rSges(7,3) * t246;
t392 = rSges(7,2) * t250;
t359 = (-pkin(5) * t344 + qJ(6) * t340) * t249 + (-qJ(6) * t344 + (-pkin(5) * qJD(5) + qJD(6)) * t250) * t246 + (-rSges(7,1) * t246 + rSges(7,3) * t249) * t340 + (-t247 * t304 + t392) * qJD(4);
t424 = t248 * t359;
t423 = t431 * t198 + t199 * t433;
t244 = t248 ^ 2;
t245 = t251 ^ 2;
t422 = t244 + t245;
t222 = rSges(5,1) * t323;
t308 = rSges(5,1) * t247 + rSges(5,2) * t250;
t400 = -pkin(1) - qJ(3);
t264 = -t308 + t400;
t408 = -rSges(5,3) - pkin(7);
t253 = t248 * t264 + t251 * t408;
t350 = qJ(2) * t345 + qJD(2) * t248;
t330 = qJD(3) * t251 + t350;
t103 = -rSges(5,2) * t319 + qJD(1) * t253 + t222 + t330;
t393 = rSges(5,2) * t247;
t218 = rSges(5,1) * t250 - t393;
t239 = qJD(2) * t251;
t349 = pkin(7) * t346 + t239;
t104 = (-t218 * qJD(4) - qJD(3)) * t248 + ((rSges(5,3) - qJ(2)) * t248 + t264 * t251) * qJD(1) + t349;
t421 = t248 * t103 + t104 * t251;
t242 = t251 * qJ(2);
t156 = t242 + t253;
t240 = t248 * qJ(2);
t348 = t251 * pkin(1) + t240;
t329 = t251 * qJ(3) + t348;
t351 = rSges(5,1) * t367 + rSges(5,2) * t361;
t157 = t248 * t408 + t329 + t351;
t420 = t156 * t251 + t157 * t248;
t289 = Icges(5,5) * t247 + Icges(5,6) * t250;
t174 = Icges(5,3) * t251 + t248 * t289;
t419 = qJD(1) * t174;
t389 = Icges(5,4) * t247;
t292 = Icges(5,2) * t250 + t389;
t178 = Icges(5,6) * t251 + t248 * t292;
t388 = Icges(5,4) * t250;
t295 = Icges(5,1) * t247 + t388;
t182 = Icges(5,5) * t251 + t248 * t295;
t356 = -rSges(7,2) * t361 + t423;
t358 = -rSges(7,2) * t363 + t435;
t418 = t248 * t356 - t251 * t358;
t417 = 2 * m(5);
t416 = 2 * m(6);
t415 = 2 * m(7);
t414 = t247 / 0.2e1;
t410 = rSges(3,2) - pkin(1);
t409 = -rSges(7,2) - pkin(8);
t407 = -rSges(6,3) - pkin(8);
t406 = m(5) * t218;
t405 = pkin(4) * t247;
t404 = pkin(4) * t250;
t403 = pkin(7) * t251;
t395 = t256 * rSges(7,2) + t198 * qJD(6) - t431 * t117 + t118 * t433;
t394 = -rSges(7,2) * t257 + t434;
t391 = rSges(6,3) * t250;
t305 = rSges(6,1) * t249 - rSges(6,2) * t246;
t153 = (-rSges(6,1) * t246 - rSges(6,2) * t249) * t340 + (-t247 * t305 + t391) * qJD(4);
t379 = t153 * t251;
t376 = t178 * t247;
t375 = t178 * t250;
t179 = -Icges(5,6) * t248 + t251 * t292;
t374 = t179 * t247;
t373 = t179 * t250;
t372 = t182 * t247;
t371 = t182 * t250;
t183 = -Icges(5,5) * t248 + t251 * t295;
t370 = t183 * t247;
t369 = t183 * t250;
t306 = -t197 * rSges(6,1) + t196 * rSges(6,2);
t137 = -rSges(6,3) * t363 - t306;
t232 = pkin(8) * t363;
t337 = t248 * t405;
t203 = -t232 + t337;
t357 = -t137 - t203;
t355 = rSges(7,2) * t247 + (pkin(5) * t249 + qJ(6) * t246 + t304) * t250;
t353 = t199 * rSges(6,1) - t198 * rSges(6,2);
t212 = (pkin(8) * t250 - t405) * qJD(4);
t219 = pkin(8) * t247 + t404;
t352 = t248 * t212 + t219 * t345;
t175 = -Icges(5,3) * t248 + t251 * t289;
t347 = qJD(1) * t175;
t338 = -rSges(4,3) + t400;
t301 = t52 * t248 - t251 * t51;
t303 = t50 * t248 - t251 * t49;
t336 = -t303 / 0.2e1 - t301 / 0.2e1;
t297 = t56 * t248 - t251 * t55;
t299 = t54 * t248 - t251 * t53;
t335 = -t297 / 0.2e1 - t299 / 0.2e1;
t333 = -t203 - t358;
t332 = pkin(4) * t323 + pkin(8) * t256;
t228 = pkin(8) * t327;
t331 = t228 + t349;
t145 = (-Icges(6,2) * t249 - t387) * t340 + (Icges(6,6) * t250 - t247 * t291) * qJD(4);
t318 = t432 * t342 + (-t279 * t344 + (-t145 * t246 + (-t177 * t249 - t181 * t246) * qJD(5)) * t250 + t430) * t247;
t211 = t308 * qJD(4);
t317 = t211 * t422;
t316 = t359 * t251;
t315 = t355 * t248;
t314 = t355 * t251;
t313 = qJD(1) * t355;
t81 = t118 * rSges(6,1) + t117 * rSges(6,2) + rSges(6,3) * t256;
t310 = t400 - t405;
t309 = t232 + t242 - t403;
t307 = t120 * rSges(6,1) - t119 * rSges(6,2);
t139 = -rSges(6,3) * t361 + t353;
t282 = t137 * t251 - t139 * t248;
t278 = t372 + t375;
t277 = t370 + t373;
t33 = t119 * t172 + t120 * t180 + t196 * t140 - t144 * t363 + t197 * t148 + t176 * t257;
t34 = -t119 * t177 + t120 * t181 - t141 * t363 - t196 * t145 + t197 * t149 + t173 * t257;
t273 = -t19 / 0.2e1 - t21 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1;
t31 = -t117 * t172 + t118 * t180 + t198 * t140 - t144 * t361 + t199 * t148 + t176 * t256;
t32 = t117 * t177 + t118 * t181 - t141 * t361 - t198 * t145 + t199 * t149 + t173 * t256;
t272 = -t20 / 0.2e1 - t22 / 0.2e1 - t31 / 0.2e1 - t32 / 0.2e1;
t271 = -t86 / 0.2e1 - t87 / 0.2e1 - t59 / 0.2e1 - t61 / 0.2e1;
t270 = t88 / 0.2e1 + t89 / 0.2e1 + t60 / 0.2e1 + t62 / 0.2e1;
t269 = rSges(3,3) * t251 + t248 * t410;
t186 = rSges(6,3) * t247 + t250 * t305;
t268 = -t248 * t153 - t186 * t345;
t233 = pkin(4) * t367;
t267 = -t248 * pkin(7) + t233 + t329;
t266 = t310 + t392;
t265 = t310 + t391;
t263 = t278 * t251;
t262 = t277 * t248;
t260 = qJD(4) * (-Icges(5,2) * t247 + t388);
t259 = qJD(4) * (Icges(5,5) * t250 - Icges(5,6) * t247);
t258 = rSges(4,2) * t251 + t248 * t338;
t252 = (t248 * t310 - t403) * qJD(1) + t330 + t332;
t210 = t251 * t219;
t209 = t248 * t219;
t204 = -pkin(8) * t361 + t233;
t201 = t251 * t212;
t190 = -rSges(3,2) * t251 + t248 * rSges(3,3) + t348;
t189 = t242 + t269;
t188 = t204 * t346;
t187 = -t248 * rSges(5,3) + t351;
t184 = rSges(5,3) * t251 + t248 * t308;
t171 = t248 * rSges(4,2) + rSges(4,3) * t251 + t329;
t170 = t242 + t258;
t165 = t239 + (t410 * t251 + (-rSges(3,3) - qJ(2)) * t248) * qJD(1);
t164 = qJD(1) * t269 + t350;
t163 = pkin(8) * t320 - t228 + (t247 * t345 + t324) * pkin(4);
t162 = -qJD(1) * t337 + t332;
t159 = t186 * t251 + t210;
t158 = t186 * t248 + t209;
t155 = -qJD(3) * t248 + t239 + ((-rSges(4,2) - qJ(2)) * t248 + t338 * t251) * qJD(1);
t154 = qJD(1) * t258 + t330;
t143 = t248 * t259 + t347;
t142 = t251 * t259 - t419;
t108 = t210 + t314;
t107 = t209 + t315;
t106 = t247 * t139 + t186 * t361;
t105 = -t137 * t247 - t186 * t363;
t102 = t361 * t407 + t267 + t353;
t101 = t248 * t265 + t306 + t309;
t100 = -t248 * t175 + t251 * t277;
t99 = -t248 * t174 + t263;
t96 = t175 * t251 + t262;
t95 = t174 * t251 + t248 * t278;
t92 = t282 * t250;
t91 = -t268 + t352;
t90 = t379 + t201 + (-t186 - t219) * t346;
t85 = t361 * t409 + t267 + t423;
t84 = t248 * t266 + t309 - t435;
t83 = rSges(6,3) * t257 + t307;
t65 = (-t139 - t204) * t251 + t357 * t248;
t64 = t247 * t356 + t250 * t314;
t63 = -t247 * t358 - t355 * t363;
t58 = qJD(1) * t314 + t352 + t424;
t57 = t201 + t316 + (-t219 - t355) * t346;
t48 = t418 * t250;
t47 = (-qJD(3) + (t247 * t407 - t404) * qJD(4)) * t248 + (t251 * t265 - t240) * qJD(1) - t307 + t331;
t46 = t252 + t81;
t45 = (-t204 - t356) * t251 + t333 * t248;
t44 = (t186 * t343 - t83) * t247 + (-qJD(4) * t137 + t268) * t250;
t43 = (-t186 * t341 + t81) * t247 + (qJD(4) * t139 - t186 * t346 + t379) * t250;
t40 = (-qJD(3) + (t247 * t409 - t404) * qJD(4)) * t248 + (t251 * t266 - t240) * qJD(1) + t331 + t434;
t39 = t252 + t395;
t30 = t282 * t344 + (t248 * t81 - t251 * t83 + (t248 * t137 + t139 * t251) * qJD(1)) * t250;
t29 = t188 + (qJD(1) * t139 - t163 - t83) * t248 + (qJD(1) * t357 - t162 - t81) * t251;
t24 = (qJD(4) * t315 + t394) * t247 + (-qJD(4) * t358 - t251 * t313 - t424) * t250;
t23 = (-t341 * t355 + t395) * t247 + (qJD(4) * t356 - t248 * t313 + t316) * t250;
t18 = -t119 * t130 + t120 * t134 + t126 * t257 - t196 * t72 + t197 * t76 - t363 * t68;
t17 = -t119 * t129 + t120 * t133 + t125 * t257 - t196 * t73 + t197 * t77 - t363 * t69;
t16 = t119 * t124 + t120 * t132 + t128 * t257 + t196 * t66 + t197 * t74 - t363 * t70;
t15 = t119 * t123 + t120 * t131 + t127 * t257 + t196 * t67 + t197 * t75 - t363 * t71;
t14 = t117 * t130 + t118 * t134 + t126 * t256 - t198 * t72 + t199 * t76 - t361 * t68;
t13 = t117 * t129 + t118 * t133 + t125 * t256 - t198 * t73 + t199 * t77 - t361 * t69;
t12 = -t117 * t124 + t118 * t132 + t128 * t256 + t198 * t66 + t199 * t74 - t361 * t70;
t11 = -t117 * t123 + t118 * t131 + t127 * t256 + t198 * t67 + t199 * t75 - t361 * t71;
t10 = t188 + (qJD(1) * t356 - t163 + t394) * t248 + (qJD(1) * t333 - t162 - t395) * t251;
t9 = -t418 * t344 + (t394 * t251 + t395 * t248 + (t248 * t358 + t251 * t356) * qJD(1)) * t250;
t8 = -qJD(1) * t300 + t17 * t251 - t18 * t248;
t7 = -qJD(1) * t302 + t15 * t251 - t16 * t248;
t6 = -qJD(1) * t296 + t13 * t251 - t14 * t248;
t5 = -qJD(1) * t298 + t11 * t251 - t12 * t248;
t4 = (qJD(4) * t300 + t34) * t247 + (qJD(1) * t301 + qJD(4) * t87 - t17 * t248 - t18 * t251) * t250;
t3 = (qJD(4) * t302 + t33) * t247 + (qJD(1) * t303 + qJD(4) * t86 - t15 * t248 - t16 * t251) * t250;
t2 = (qJD(4) * t296 + t32) * t247 + (qJD(1) * t297 + qJD(4) * t89 - t13 * t248 - t14 * t251) * t250;
t1 = (qJD(4) * t298 + t31) * t247 + (qJD(1) * t299 + qJD(4) * t88 - t11 * t248 - t12 * t251) * t250;
t25 = [-t250 * t260 - t295 * t342 - t181 * t322 - t177 * t321 - t180 * t325 - t172 * t326 + (t39 * t85 + t40 * t84) * t415 + (t101 * t47 + t102 * t46) * t416 + (t103 * t157 + t104 * t156) * t417 + 0.2e1 * m(4) * (t154 * t171 + t155 * t170) + 0.2e1 * m(3) * (t164 * t190 + t165 * t189) - t145 * t368 + t430 + (-Icges(5,1) * t250 + t292 + t389) * t344; m(7) * (t248 * t40 - t251 * t39 + (t248 * t85 + t251 * t84) * qJD(1)) + m(6) * (t248 * t47 - t251 * t46 + (t101 * t251 + t102 * t248) * qJD(1)) + m(5) * (qJD(1) * t420 - t103 * t251 + t248 * t104) + m(3) * (-t164 * t251 + t248 * t165 + (t189 * t251 + t190 * t248) * qJD(1)) + m(4) * (-t154 * t251 + t248 * t155 + (t170 * t251 + t171 * t248) * qJD(1)); 0; m(7) * (t248 * t39 + t251 * t40 + (-t248 * t84 + t251 * t85) * qJD(1)) + m(6) * (t248 * t46 + t251 * t47 + (-t101 * t248 + t102 * t251) * qJD(1)) + m(5) * ((-t156 * t248 + t157 * t251) * qJD(1) + t421) + m(4) * (t248 * t154 + t155 * t251 + (-t170 * t248 + t171 * t251) * qJD(1)); 0; 0; (-(qJD(1) * t179 + t248 * t260) * t247 / 0.2e1 + t183 * t436 + (-t375 / 0.2e1 - t372 / 0.2e1) * qJD(4) - t273) * t251 + ((-qJD(1) * t178 + t251 * t260) * t414 + t182 * t436 + (t373 / 0.2e1 + t370 / 0.2e1) * qJD(4) + t272) * t248 + m(5) * (-t211 * t420 + t218 * t421) + m(7) * (t107 * t39 + t108 * t40 + t57 * t84 + t58 * t85) + m(6) * (t101 * t90 + t102 * t91 + t158 * t46 + t159 * t47) - (t244 / 0.2e1 + t245 / 0.2e1) * t289 * qJD(4) + ((t157 * t406 + t374 / 0.2e1 - t369 / 0.2e1 - t270) * t251 + (-t156 * t406 + t376 / 0.2e1 - t371 / 0.2e1 + t271) * t248) * qJD(1); m(6) * (t90 * t248 - t251 * t91 + (t158 * t248 + t159 * t251) * qJD(1)) + m(7) * (t57 * t248 - t251 * t58 + (t107 * t248 + t108 * t251) * qJD(1)); m(6) * (t91 * t248 + t251 * t90 + (t158 * t251 - t159 * t248) * qJD(1)) + m(7) * (t58 * t248 + t251 * t57 + (t107 * t251 - t108 * t248) * qJD(1)) - m(5) * t317; (t45 * t10 + t107 * t58 + t108 * t57) * t415 + t251 * t8 + t251 * t7 - t248 * t6 - t248 * t5 + (t158 * t91 + t159 * t90 + t65 * t29) * t416 - t248 * ((t248 * t142 + (-t99 + t262) * qJD(1)) * t248 + (-t100 * qJD(1) + (-t178 * t344 + t182 * t342 - t419) * t251 + (-t143 + (-t369 + t374) * qJD(4) + (t175 - t278) * qJD(1)) * t248) * t251) + (-t218 * t317 + (-t251 * t222 + (-t218 * t244 + t245 * t393) * qJD(4) + (rSges(5,3) * t422 - t251 * t184 + t248 * t187) * qJD(1)) * (-t248 * t184 - t187 * t251)) * t417 + t251 * ((t251 * t143 + (-t96 + t263) * qJD(1)) * t251 + (-t95 * qJD(1) + (t179 * t344 - t183 * t342 + t347) * t248 + (-t142 + (t371 - t376) * qJD(4) + (-t174 - t277) * qJD(1)) * t251) * t248) + (t96 * t248 - t251 * t95 + t301 + t303) * t346 + (t100 * t248 - t251 * t99 + t297 + t299) * t345; m(7) * (t23 * t85 + t24 * t84 + t64 * t39 + t63 * t40) + m(6) * (t101 * t44 + t102 * t43 + t105 * t47 + t106 * t46) + (-t248 * t271 + t251 * t270) * t344 + (t272 * t251 + t273 * t248 + (t248 * t270 + t251 * t271) * qJD(1)) * t250 + t318; m(6) * (t44 * t248 - t251 * t43 + (t105 * t251 + t106 * t248) * qJD(1)) + m(7) * (-t23 * t251 + t24 * t248 + (t248 * t64 + t251 * t63) * qJD(1)); m(6) * (t43 * t248 + t251 * t44 + (-t105 * t248 + t106 * t251) * qJD(1)) + m(7) * (t23 * t248 + t24 * t251 + (-t248 * t63 + t251 * t64) * qJD(1)); m(7) * (t48 * t10 + t107 * t23 + t108 * t24 + t9 * t45 + t63 * t57 + t64 * t58) + m(6) * (t105 * t90 + t106 * t91 + t158 * t43 + t159 * t44 - t29 * t92 + t30 * t65) + (t4 / 0.2e1 + t3 / 0.2e1 + t335 * t344) * t251 + (-t2 / 0.2e1 - t1 / 0.2e1 + t336 * t344) * t248 + ((t248 * t335 - t251 * t336) * qJD(1) - (t7 + t8) * t248 / 0.2e1 - (t5 + t6) * t251 / 0.2e1 - t428 * qJD(4) / 0.2e1) * t250 + (qJD(1) * t429 + t426 * t248 + t427 * t251) * t414 - (t425 * t248 + t398 * t251) * qJD(1) / 0.2e1; (t23 * t64 + t24 * t63 + t48 * t9) * t415 + (t105 * t44 + t106 * t43 - t30 * t92) * t416 + (((t247 * t396 + t398) * t251 + (t247 * t397 + t425) * t248) * qJD(4) + t318) * t247 + ((-t1 - t2) * t251 + (-t3 - t4) * t248 + (-t427 * t248 + t426 * t251) * t247 + (t247 * t432 + t250 * t429) * qJD(4) + (t247 * t428 + t398 * t248 - t251 * t425) * qJD(1)) * t250; m(7) * (-t117 * t84 + t119 * t85 + t196 * t39 + t198 * t40); m(7) * (-t117 * t248 - t119 * t251 + (t196 * t248 + t198 * t251) * qJD(1)); m(7) * (-t117 * t251 + t119 * t248 + (t196 * t251 - t198 * t248) * qJD(1)); m(7) * (t45 * t321 + t107 * t119 - t108 * t117 + t196 * t58 + t198 * t57 + (t10 * t250 - t344 * t45) * t246); m(7) * (t48 * t321 - t117 * t63 + t119 * t64 + t196 * t23 + t198 * t24 + (t250 * t9 - t344 * t48) * t246); (-t117 * t198 + t119 * t196 + (t321 - t326) * t368) * t415;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
