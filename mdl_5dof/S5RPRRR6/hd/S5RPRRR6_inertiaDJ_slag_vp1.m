% Calculate time derivative of joint inertia matrix for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:01:12
% DurationCPUTime: 7.55s
% Computational Cost: add. (25338->586), mult. (23270->836), div. (0->0), fcn. (21919->10), ass. (0->324)
t249 = qJ(3) + qJ(4);
t244 = sin(t249);
t247 = qJD(3) + qJD(4);
t253 = cos(qJ(5));
t250 = sin(qJ(5));
t391 = Icges(6,4) * t253;
t295 = -Icges(6,2) * t250 + t391;
t245 = cos(t249);
t366 = t245 * t247;
t392 = Icges(6,4) * t250;
t141 = t295 * t366 + (Icges(6,6) * t247 + (-Icges(6,2) * t253 - t392) * qJD(5)) * t244;
t191 = -Icges(6,6) * t245 + t244 * t295;
t299 = Icges(6,1) * t253 - t392;
t192 = -Icges(6,5) * t245 + t244 * t299;
t432 = -t250 * t141 + (-t191 * t253 - t192 * t250) * qJD(5);
t292 = Icges(6,5) * t253 - Icges(6,6) * t250;
t140 = t292 * t366 + (Icges(6,3) * t247 + (-Icges(6,5) * t250 - Icges(6,6) * t253) * qJD(5)) * t244;
t362 = t247 * t250;
t431 = -t191 * t362 - t140;
t394 = Icges(5,4) * t244;
t300 = Icges(5,1) * t245 - t394;
t211 = Icges(5,2) * t245 + t394;
t363 = t247 * t211;
t430 = t300 * t247 - t363;
t248 = qJ(1) + pkin(9);
t243 = cos(t248);
t210 = Icges(5,5) * t244 + Icges(5,6) * t245;
t273 = t210 * t247;
t242 = sin(t248);
t293 = Icges(5,5) * t245 - Icges(5,6) * t244;
t169 = -Icges(5,3) * t243 + t242 * t293;
t419 = qJD(1) * t169;
t122 = -t243 * t273 - t419;
t170 = Icges(5,3) * t242 + t243 * t293;
t350 = qJD(1) * t170;
t123 = -t242 * t273 + t350;
t393 = Icges(5,4) * t245;
t296 = -Icges(5,2) * t244 + t393;
t171 = -Icges(5,6) * t243 + t242 * t296;
t124 = -qJD(1) * t171 - t243 * t363;
t173 = -Icges(5,5) * t243 + t242 * t300;
t212 = Icges(5,1) * t244 + t393;
t377 = t212 * t247;
t126 = -qJD(1) * t173 - t243 * t377;
t172 = Icges(5,6) * t242 + t243 * t296;
t174 = Icges(5,5) * t242 + t243 * t300;
t285 = t172 * t244 - t174 * t245;
t127 = qJD(1) * t174 - t242 * t377;
t320 = t171 * t247 - t127;
t125 = qJD(1) * t172 - t242 * t363;
t322 = t173 * t247 + t125;
t367 = t244 * t247;
t286 = t171 * t244 - t173 * t245;
t423 = t243 * t286;
t85 = -t169 * t243 - t242 * t286;
t425 = t242 * t285;
t86 = -t170 * t243 - t425;
t319 = -qJD(5) * t245 + qJD(1);
t263 = t244 * t362 + t253 * t319;
t345 = qJD(1) * t245;
t318 = -qJD(5) + t345;
t370 = t243 * t250;
t118 = t242 * t263 - t318 * t370;
t361 = t247 * t253;
t262 = -t244 * t361 + t250 * t319;
t369 = t243 * t253;
t119 = t242 * t262 + t318 * t369;
t365 = t245 * t250;
t198 = -t242 * t365 - t369;
t364 = t245 * t253;
t199 = t242 * t364 - t370;
t376 = t242 * t244;
t130 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t376;
t132 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t376;
t134 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t376;
t346 = qJD(1) * t244;
t267 = t242 * t366 + t243 * t346;
t68 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t267;
t70 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t267;
t72 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t267;
t15 = t118 * t132 + t119 * t134 + t130 * t267 + t198 * t70 + t199 * t72 + t376 * t68;
t373 = t242 * t253;
t200 = -t243 * t365 + t373;
t374 = t242 * t250;
t201 = t243 * t364 + t374;
t372 = t243 * t244;
t131 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t372;
t133 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t372;
t135 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t372;
t116 = t243 * t263 + t318 * t374;
t117 = t243 * t262 - t318 * t373;
t329 = t242 * t346;
t334 = t243 * t366;
t265 = -t329 + t334;
t67 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t265;
t69 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t265;
t71 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t265;
t16 = t118 * t133 + t119 * t135 + t131 * t267 + t198 * t69 + t199 * t71 + t376 * t67;
t50 = t130 * t376 + t132 * t198 + t134 * t199;
t51 = t131 * t376 + t133 * t198 + t135 * t199;
t307 = t242 * t50 + t243 * t51;
t9 = qJD(1) * t307 - t15 * t243 + t16 * t242;
t428 = -(t123 * t243 + (t86 + t423) * qJD(1)) * t243 - (t85 * qJD(1) + (-t124 * t244 + t126 * t245 - t172 * t366 - t174 * t367 + t350) * t242 + (-t122 + t320 * t245 + t322 * t244 + (-t169 - t285) * qJD(1)) * t243) * t242 - t9;
t417 = 2 * m(4);
t251 = sin(qJ(3));
t400 = rSges(4,2) * t251;
t254 = cos(qJ(3));
t402 = rSges(4,1) * t254;
t313 = -t400 + t402;
t399 = rSges(4,3) * t243;
t188 = t242 * t313 - t399;
t340 = t243 * t400;
t237 = t242 * rSges(4,3);
t352 = t243 * t402 + t237;
t189 = -t340 + t352;
t231 = rSges(4,1) * t251 + rSges(4,2) * t254;
t276 = qJD(3) * t231;
t344 = qJD(1) * t251;
t328 = t242 * t344;
t347 = qJD(1) * t243;
t260 = rSges(4,2) * t328 + rSges(4,3) * t347 - t243 * t276;
t426 = t242 * t276;
t64 = (qJD(1) * t188 + t260) * t243 + (-t426 + (-t189 - t340 + t237) * qJD(1)) * t242;
t427 = t417 * t64;
t239 = pkin(3) * t254 + pkin(2);
t404 = pkin(2) - t239;
t424 = t242 * t404;
t395 = Icges(4,4) * t254;
t298 = -Icges(4,2) * t251 + t395;
t185 = Icges(4,6) * t242 + t243 * t298;
t396 = Icges(4,4) * t251;
t302 = Icges(4,1) * t254 - t396;
t187 = Icges(4,5) * t242 + t243 * t302;
t283 = t185 * t251 - t187 * t254;
t422 = t283 * t243;
t310 = -t199 * rSges(6,1) - t198 * rSges(6,2);
t138 = rSges(6,3) * t376 - t310;
t139 = t201 * rSges(6,1) + t200 * rSges(6,2) + rSges(6,3) * t372;
t421 = -t242 * t138 - t243 * t139;
t238 = t243 * pkin(6);
t405 = sin(qJ(1)) * pkin(1);
t420 = t238 - t405;
t294 = Icges(4,5) * t254 - Icges(4,6) * t251;
t182 = -Icges(4,3) * t243 + t242 * t294;
t184 = -Icges(4,6) * t243 + t242 * t298;
t186 = -Icges(4,5) * t243 + t242 * t302;
t282 = t211 * t244 - t212 * t245;
t418 = qJD(1) * t282 + t293 * t247;
t416 = 2 * m(5);
t415 = 2 * m(6);
t240 = t242 ^ 2;
t241 = t243 ^ 2;
t414 = t242 / 0.2e1;
t413 = -t243 / 0.2e1;
t412 = -rSges(6,3) - pkin(8);
t411 = m(4) * t231;
t213 = rSges(5,1) * t244 + rSges(5,2) * t245;
t410 = m(5) * t213;
t409 = pkin(3) * t251;
t408 = pkin(4) * t244;
t407 = pkin(4) * t245;
t406 = t242 * pkin(6);
t246 = cos(qJ(1)) * pkin(1);
t256 = -pkin(7) - pkin(6);
t403 = -pkin(6) - t256;
t401 = rSges(5,1) * t245;
t291 = -t132 * t250 + t134 * t253;
t19 = (t247 * t291 - t68) * t245 + (t130 * t247 - t250 * t70 + t253 * t72 + (-t132 * t253 - t134 * t250) * qJD(5)) * t244;
t398 = t19 * t243;
t290 = -t133 * t250 + t135 * t253;
t20 = (t247 * t290 - t67) * t245 + (t131 * t247 - t250 * t69 + t253 * t71 + (-t133 * t253 - t135 * t250) * qJD(5)) * t244;
t397 = t20 * t242;
t236 = t242 * rSges(5,3);
t384 = t184 * t254;
t383 = t185 * t254;
t382 = t186 * t251;
t381 = t187 * t251;
t312 = -rSges(5,2) * t244 + t401;
t205 = t312 * t247;
t378 = t205 * t242;
t375 = t242 * t247;
t371 = t243 * t245;
t368 = t243 * t256;
t197 = pkin(4) * t371 + pkin(8) * t372;
t359 = -t139 - t197;
t309 = rSges(6,1) * t253 - rSges(6,2) * t250;
t143 = t309 * t366 + (rSges(6,3) * t247 + (-rSges(6,1) * t250 - rSges(6,2) * t253) * qJD(5)) * t244;
t316 = pkin(8) * t244 + t407;
t358 = -t316 * t247 - t143;
t166 = t238 + t368 - t424;
t220 = t243 * t239;
t167 = -t243 * pkin(2) + t242 * t403 + t220;
t357 = t242 * t166 + t243 * t167;
t175 = -t243 * rSges(5,3) + t242 * t312;
t176 = rSges(5,1) * t371 - rSges(5,2) * t372 + t236;
t107 = t242 * t175 + t243 * t176;
t193 = -rSges(6,3) * t245 + t244 * t309;
t348 = qJD(1) * t242;
t168 = t193 * t348;
t214 = -pkin(8) * t245 + t408;
t356 = t214 * t348 + t168;
t355 = -t193 - t214;
t354 = rSges(5,2) * t329 + rSges(5,3) * t347;
t343 = qJD(3) * t251;
t339 = pkin(3) * t343;
t353 = t242 * t339 + t256 * t348;
t351 = t240 + t241;
t183 = Icges(4,3) * t242 + t243 * t294;
t349 = qJD(1) * t183;
t342 = qJD(3) * t254;
t338 = pkin(3) * t342;
t61 = -t130 * t245 + t244 * t291;
t190 = -Icges(6,3) * t245 + t244 * t292;
t83 = t190 * t376 + t191 * t198 + t192 * t199;
t337 = t61 / 0.2e1 + t83 / 0.2e1;
t62 = -t131 * t245 + t244 * t290;
t84 = t190 * t372 + t191 * t200 + t192 * t201;
t336 = t62 / 0.2e1 + t84 / 0.2e1;
t266 = -t242 * t345 - t243 * t367;
t277 = t213 * t247;
t333 = t242 * (-t242 * t277 + (t243 * t312 + t236) * qJD(1)) + t243 * (rSges(5,1) * t266 - rSges(5,2) * t334 + t354) + t175 * t347;
t332 = t117 * rSges(6,1) + t116 * rSges(6,2) + rSges(6,3) * t334;
t142 = t299 * t366 + (Icges(6,5) * t247 + (-Icges(6,1) * t250 - t391) * qJD(5)) * t244;
t331 = t244 * t253 * t142 + t245 * t192 * t361 + t190 * t367;
t330 = t242 * ((-t243 * t404 - t406) * qJD(1) - t353) + t243 * (-t243 * t339 + (t243 * t403 + t424) * qJD(1)) + t166 * t347;
t326 = t348 / 0.2e1;
t325 = t347 / 0.2e1;
t324 = -t213 - t409;
t155 = t355 * t243;
t323 = t174 * t247 + t124;
t321 = -t172 * t247 + t126;
t196 = t316 * t242;
t63 = t242 * t196 + t243 * t197 - t421;
t317 = t355 - t409;
t315 = -t242 * t256 + t220 + t246;
t311 = rSges(6,1) * t119 + rSges(6,2) * t118;
t308 = -t368 - t405;
t38 = t242 * t51 - t243 * t50;
t52 = t130 * t372 + t132 * t200 + t134 * t201;
t53 = t131 * t372 + t133 * t200 + t135 * t201;
t39 = t242 * t53 - t243 * t52;
t306 = t242 * t52 + t243 * t53;
t305 = t242 * t62 - t243 * t61;
t304 = t242 * t61 + t243 * t62;
t13 = t116 * t132 + t117 * t134 + t130 * t265 + t200 * t70 + t201 * t72 + t372 * t68;
t14 = t116 * t133 + t117 * t135 + t131 * t265 + t200 * t69 + t201 * t71 + t372 * t67;
t8 = qJD(1) * t306 - t13 * t243 + t14 * t242;
t87 = t169 * t242 - t423;
t88 = t170 * t242 - t243 * t285;
t303 = t39 * t347 + t38 * t348 + (-t87 * t347 - t85 * t348) * t243 + (t8 + (t88 * qJD(1) + (t125 * t244 - t127 * t245 + t171 * t366 + t173 * t367 - t419) * t243) * t243 + t86 * t348 + t88 * t347 + ((t87 + t425) * qJD(1) + (-t123 + t321 * t245 - t323 * t244 + (t170 - t286) * qJD(1)) * t243 + t122 * t242) * t242) * t242;
t301 = Icges(4,1) * t251 + t395;
t297 = Icges(4,2) * t254 + t396;
t289 = t138 * t243 - t139 * t242;
t284 = t184 * t251 - t186 * t254;
t281 = -t338 + t358;
t280 = -pkin(2) - t313;
t209 = pkin(8) * t334;
t75 = -rSges(6,3) * t329 + t332;
t76 = rSges(6,3) * t267 + t311;
t279 = t242 * t76 + t243 * t75 + t242 * (t267 * pkin(8) + (-t242 * t367 + t243 * t345) * pkin(4)) + t243 * (pkin(4) * t266 - pkin(8) * t329 + t209) + (t138 + t196) * t347;
t137 = t317 * t243;
t278 = -t239 - t312;
t272 = qJD(3) * t301;
t271 = qJD(3) * t297;
t270 = qJD(3) * (-Icges(4,5) * t251 - Icges(4,6) * t254);
t269 = t244 * t412 - t239 - t407;
t264 = t428 * t243 + t303;
t26 = t244 * t307 - t245 * t83;
t27 = t244 * t306 - t245 * t84;
t36 = t116 * t191 + t117 * t192 + t140 * t372 + t141 * t200 + t142 * t201 + t190 * t265;
t3 = (t247 * t306 - t36) * t245 + (-qJD(1) * t39 + t13 * t242 + t14 * t243 + t247 * t84) * t244;
t37 = t118 * t191 + t119 * t192 + t140 * t376 + t141 * t198 + t142 * t199 + t190 * t267;
t4 = (t247 * t307 - t37) * t245 + (-qJD(1) * t38 + t15 * t242 + t16 * t243 + t247 * t83) * t244;
t261 = t3 * t414 + t9 * t376 / 0.2e1 + t4 * t413 - t245 * (qJD(1) * t304 + t397 - t398) / 0.2e1 + t26 * t326 - t39 * t329 / 0.2e1 + t305 * t367 / 0.2e1 + t8 * t372 / 0.2e1 + (t242 * t38 + t243 * t39) * t366 / 0.2e1 + (t244 * t38 + t27) * t325;
t259 = t242 * t269 + t308;
t203 = t296 * t247;
t258 = qJD(1) * t210 + t430 * t245 + (-t203 - t377) * t244;
t257 = -t398 / 0.2e1 + t397 / 0.2e1 + (t418 * t242 + t258 * t243 + t244 * t321 + t245 * t323 + t36) * t414 + (t258 * t242 - t243 * t418 - t244 * t320 + t245 * t322 + t37) * t413 + (t171 * t245 + t173 * t244 - t210 * t243 - t242 * t282 + t61 + t83) * t326 + (t172 * t245 + t174 * t244 + t210 * t242 - t243 * t282 + t62 + t84) * t325;
t227 = pkin(3) * t328;
t219 = t313 * qJD(3);
t179 = t324 * t243;
t178 = t324 * t242;
t157 = t406 + t246 + (pkin(2) - t400) * t243 + t352;
t156 = t242 * t280 + t399 + t420;
t154 = t355 * t242;
t149 = t242 * t270 + t349;
t148 = -qJD(1) * t182 + t243 * t270;
t147 = t176 + t315;
t146 = -t405 + (rSges(5,3) - t256) * t243 + t278 * t242;
t136 = t317 * t242;
t111 = -t213 * t347 - t378 + (-t242 * t342 - t243 * t344) * pkin(3);
t110 = t213 * t348 + t227 + (-t205 - t338) * t243;
t106 = t426 + (-t246 + (-rSges(4,3) - pkin(6)) * t242 + t280 * t243) * qJD(1);
t105 = ((-pkin(2) - t402) * t242 + t420) * qJD(1) + t260;
t100 = -t190 * t245 + (-t191 * t250 + t192 * t253) * t244;
t99 = -t139 * t245 - t193 * t372;
t98 = t138 * t245 + t193 * t376;
t97 = t213 * t375 + (t243 * t278 - t236 - t246) * qJD(1) + t353;
t96 = (-t277 - t339) * t243 + ((-t239 - t401) * t242 + t308) * qJD(1) + t354;
t95 = t183 * t242 - t422;
t94 = t182 * t242 - t284 * t243;
t93 = -t183 * t243 - t242 * t283;
t92 = -t182 * t243 - t242 * t284;
t91 = t315 - t359;
t90 = t259 + t310;
t89 = t100 * t367;
t82 = t289 * t244;
t81 = qJD(1) * t155 + t242 * t358;
t80 = t243 * t358 + t356;
t77 = t107 + t357;
t74 = qJD(1) * t137 + t242 * t281;
t73 = t243 * t281 + t227 + t356;
t58 = -t176 * t348 + t333;
t49 = (t245 * t412 + t408) * t375 + (t243 * t269 - t246) * qJD(1) - t311 + t353;
t48 = t209 + (-pkin(4) * t367 - t339) * t243 + t259 * qJD(1) + t332;
t45 = t63 + t357;
t44 = t432 * t244 + t431 * t245 + t331;
t43 = (t193 * t375 + t76) * t245 + (-t138 * t247 + t143 * t242 + t193 * t347) * t244;
t42 = (-t193 * t243 * t247 - t75) * t245 + (t139 * t247 - t143 * t243 + t168) * t244;
t41 = (-t167 - t176) * t348 + t330 + t333;
t25 = t289 * t366 + (qJD(1) * t421 - t242 * t75 + t243 * t76) * t244;
t22 = t348 * t359 + t279;
t21 = (-t167 + t359) * t348 + t279 + t330;
t1 = [(t48 * t91 + t49 * t90) * t415 + t212 * t366 + (t146 * t97 + t147 * t96) * t416 + (t105 * t157 + t106 * t156) * t417 + t331 + (t302 - t297) * t343 + (t301 + t298) * t342 + (t203 + t431) * t245 + (t430 + t432) * t244; 0; 0; (t241 / 0.2e1 + t240 / 0.2e1) * t294 * qJD(3) + m(6) * (t136 * t48 + t137 * t49 + t73 * t90 + t74 * t91) + m(5) * (t110 * t146 + t111 * t147 + t178 * t96 + t179 * t97) + t257 + (-qJD(3) * t284 + (qJD(1) * t185 - t242 * t271) * t254 + (qJD(1) * t187 - t242 * t272) * t251) * t413 + (-qJD(3) * t283 + (-qJD(1) * t184 - t243 * t271) * t254 + (-qJD(1) * t186 - t243 * t272) * t251) * t414 + ((t383 / 0.2e1 + t381 / 0.2e1 - t157 * t411) * t243 + (t384 / 0.2e1 + t382 / 0.2e1 + t156 * t411) * t242) * qJD(1) + m(4) * ((-t105 * t242 - t106 * t243) * t231 + (-t156 * t243 - t157 * t242) * t219); m(4) * t64 + m(5) * t41 + m(6) * t21; (t136 * t74 + t137 * t73 + t21 * t45) * t415 + (t110 * t179 + t111 * t178 + t41 * t77) * t416 + t351 * t231 * t219 * t417 + t303 + (t188 * t427 + t95 * t347 + t93 * t348 + (t94 * qJD(1) + (t283 * qJD(1) + t148) * t242) * t242) * t242 + (t189 * t427 - t92 * t348 - t94 * t347 + (-t93 * qJD(1) + (-qJD(1) * t284 - t149) * t243) * t243 + ((t148 - (t382 + t384) * qJD(3) + t184 * t342 + t186 * t343) * t243 + (t185 * t342 + t187 * t343 - t349 - t149 + (-t381 - t383) * qJD(3)) * t242 + (-t92 + t95 + t422 + (t183 - t284) * t242) * qJD(1)) * t242 + t428) * t243; m(6) * (t154 * t48 + t155 * t49 + t80 * t90 + t81 * t91) + t257 + (-t242 * t96 - t243 * t97 + (t146 * t242 - t147 * t243) * qJD(1)) * t410 + m(5) * (-t146 * t243 - t147 * t242) * t205; m(5) * t58 + m(6) * t22; m(6) * (t136 * t81 + t137 * t80 + t154 * t74 + t155 * t73 + t21 * t63 + t22 * t45) + m(5) * (-t179 * t205 * t243 + t107 * t41 - t178 * t378 + t58 * t77) + (-t110 * t243 - t111 * t242 + (-t178 * t243 + t179 * t242) * qJD(1)) * t410 + t264; (t205 * t213 * t351 + t107 * t58) * t416 + (t154 * t81 + t155 * t80 + t22 * t63) * t415 + t264; m(6) * (t42 * t91 + t43 * t90 + t48 * t99 + t49 * t98) + t89 + (-t44 + (t242 * t337 + t243 * t336) * t247) * t245 + ((t20 / 0.2e1 + t36 / 0.2e1) * t243 + (t19 / 0.2e1 + t37 / 0.2e1) * t242 + (-t242 * t336 + t243 * t337) * qJD(1)) * t244; m(6) * t25; t261 + m(6) * (t136 * t42 + t137 * t43 + t21 * t82 + t25 * t45 + t73 * t98 + t74 * t99); t261 + m(6) * (t154 * t42 + t155 * t43 + t22 * t82 + t25 * t63 + t80 * t98 + t81 * t99); (t25 * t82 + t42 * t99 + t43 * t98) * t415 + (t44 * t245 - t89 + (t242 * t26 + t243 * t27 - t245 * t304) * t247) * t245 + (t243 * t3 + t242 * t4 + t304 * t367 + (-t100 * t247 - t19 * t242 - t20 * t243) * t245 + (-t242 * t27 + t243 * t26 + t245 * t305) * qJD(1)) * t244;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
