% Calculate time derivative of joint inertia matrix for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP11_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:28
% EndTime: 2019-12-31 18:53:46
% DurationCPUTime: 9.92s
% Computational Cost: add. (17846->680), mult. (26636->979), div. (0->0), fcn. (25738->8), ass. (0->330)
t237 = pkin(8) + qJ(3);
t232 = sin(t237);
t233 = cos(t237);
t241 = sin(qJ(4));
t243 = cos(qJ(4));
t277 = Icges(5,5) * t243 - Icges(5,6) * t241;
t171 = -Icges(5,3) * t233 + t232 * t277;
t279 = Icges(6,4) * t243 + Icges(6,6) * t241;
t172 = -Icges(6,2) * t233 + t232 * t279;
t424 = t171 + t172;
t244 = cos(qJ(1));
t348 = t243 * t244;
t242 = sin(qJ(1));
t350 = t242 * t241;
t198 = t233 * t350 + t348;
t349 = t242 * t243;
t351 = t241 * t244;
t199 = t233 * t349 - t351;
t418 = rSges(6,3) + qJ(5);
t420 = rSges(6,1) + pkin(4);
t423 = -t418 * t198 - t420 * t199;
t333 = qJD(3) * t242;
t313 = t232 * t333;
t328 = qJD(4) * t243;
t329 = qJD(4) * t241;
t337 = qJD(1) * t244;
t338 = qJD(1) * t242;
t126 = -t241 * t313 - t244 * t329 - t243 * t338 + (t241 * t337 + t242 * t328) * t233;
t265 = (-qJD(4) * t233 + qJD(1)) * t241;
t299 = qJD(1) * t233 - qJD(4);
t332 = qJD(3) * t243;
t127 = t299 * t348 + (-t232 * t332 + t265) * t242;
t422 = t198 * qJD(5) + t418 * t126 + t420 * t127;
t421 = -qJD(1) * t232 / 0.2e1;
t372 = Icges(6,5) * t243;
t276 = Icges(6,3) * t241 + t372;
t170 = -Icges(6,6) * t233 + t232 * t276;
t375 = Icges(5,4) * t243;
t280 = -Icges(5,2) * t241 + t375;
t173 = -Icges(5,6) * t233 + t232 * t280;
t373 = Icges(6,5) * t241;
t283 = Icges(6,1) * t243 + t373;
t174 = -Icges(6,4) * t233 + t232 * t283;
t376 = Icges(5,4) * t241;
t284 = Icges(5,1) * t243 - t376;
t175 = -Icges(5,5) * t233 + t232 * t284;
t419 = t424 * t233 + ((-t174 - t175) * t243 + (-t170 + t173) * t241) * t232;
t334 = qJD(3) * t241;
t312 = t233 * t334;
t408 = t232 * t328 + t312;
t200 = t233 * t351 - t349;
t323 = t233 * t348;
t201 = t323 + t350;
t355 = t232 * t244;
t141 = Icges(6,4) * t201 + Icges(6,2) * t355 + Icges(6,6) * t200;
t137 = Icges(6,5) * t201 + Icges(6,6) * t355 + Icges(6,3) * t200;
t145 = Icges(6,1) * t201 + Icges(6,4) * t355 + Icges(6,5) * t200;
t272 = t137 * t241 + t145 * t243;
t60 = -t141 * t233 + t232 * t272;
t139 = Icges(5,5) * t201 - Icges(5,6) * t200 + Icges(5,3) * t355;
t143 = Icges(5,4) * t201 - Icges(5,2) * t200 + Icges(5,6) * t355;
t147 = Icges(5,1) * t201 - Icges(5,4) * t200 + Icges(5,5) * t355;
t270 = -t143 * t241 + t147 * t243;
t62 = -t139 * t233 + t232 * t270;
t385 = t60 + t62;
t357 = t232 * t242;
t140 = Icges(6,4) * t199 + Icges(6,2) * t357 + Icges(6,6) * t198;
t136 = Icges(6,5) * t199 + Icges(6,6) * t357 + Icges(6,3) * t198;
t144 = Icges(6,1) * t199 + Icges(6,4) * t357 + Icges(6,5) * t198;
t273 = t136 * t241 + t144 * t243;
t59 = -t140 * t233 + t232 * t273;
t138 = Icges(5,5) * t199 - Icges(5,6) * t198 + Icges(5,3) * t357;
t142 = Icges(5,4) * t199 - Icges(5,2) * t198 + Icges(5,6) * t357;
t146 = Icges(5,1) * t199 - Icges(5,4) * t198 + Icges(5,5) * t357;
t271 = -t142 * t241 + t146 * t243;
t61 = -t138 * t233 + t232 * t271;
t386 = t59 + t61;
t417 = t242 * t386 + t244 * t385;
t416 = t242 * t385 - t244 * t386;
t307 = t233 * t333;
t250 = t232 * t337 + t307;
t67 = Icges(6,5) * t127 + Icges(6,6) * t250 + Icges(6,3) * t126;
t71 = Icges(6,4) * t127 + Icges(6,2) * t250 + Icges(6,6) * t126;
t75 = Icges(6,1) * t127 + Icges(6,4) * t250 + Icges(6,5) * t126;
t19 = (qJD(3) * t273 - t71) * t233 + (qJD(3) * t140 + t241 * t67 + t243 * t75 + (t136 * t243 - t144 * t241) * qJD(4)) * t232;
t69 = Icges(5,5) * t127 - Icges(5,6) * t126 + Icges(5,3) * t250;
t73 = Icges(5,4) * t127 - Icges(5,2) * t126 + Icges(5,6) * t250;
t77 = Icges(5,1) * t127 - Icges(5,4) * t126 + Icges(5,5) * t250;
t21 = (qJD(3) * t271 - t69) * t233 + (qJD(3) * t138 - t241 * t73 + t243 * t77 + (-t142 * t243 - t146 * t241) * qJD(4)) * t232;
t415 = -t19 - t21;
t331 = qJD(3) * t244;
t310 = t232 * t331;
t124 = qJD(1) * t198 - qJD(4) * t323 + t241 * t310 - t242 * t329;
t125 = t244 * t265 + (-t242 * t299 - t310) * t243;
t306 = t233 * t331;
t314 = t232 * t338;
t249 = t306 - t314;
t66 = Icges(6,5) * t125 + Icges(6,6) * t249 - Icges(6,3) * t124;
t70 = Icges(6,4) * t125 + Icges(6,2) * t249 - Icges(6,6) * t124;
t74 = Icges(6,1) * t125 + Icges(6,4) * t249 - Icges(6,5) * t124;
t20 = (qJD(3) * t272 - t70) * t233 + (qJD(3) * t141 + t241 * t66 + t243 * t74 + (t137 * t243 - t145 * t241) * qJD(4)) * t232;
t68 = Icges(5,5) * t125 + Icges(5,6) * t124 + Icges(5,3) * t249;
t72 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t249;
t76 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t249;
t22 = (qJD(3) * t270 - t68) * t233 + (qJD(3) * t139 - t241 * t72 + t243 * t76 + (-t143 * t243 - t147 * t241) * qJD(4)) * t232;
t414 = t20 + t22;
t53 = t138 * t357 - t142 * t198 + t146 * t199;
t54 = t139 * t357 - t143 * t198 + t147 * t199;
t289 = t242 * t53 + t244 * t54;
t51 = t136 * t198 + t140 * t357 + t144 * t199;
t52 = t137 * t198 + t141 * t357 + t145 * t199;
t290 = t242 * t51 + t244 * t52;
t84 = t170 * t198 + t172 * t357 + t174 * t199;
t85 = t171 * t357 - t173 * t198 + t175 * t199;
t389 = (-t84 - t85) * t233 + (t289 + t290) * t232;
t57 = t138 * t355 - t200 * t142 + t201 * t146;
t58 = t139 * t355 - t200 * t143 + t201 * t147;
t287 = t242 * t57 + t244 * t58;
t55 = t200 * t136 + t140 * t355 + t201 * t144;
t56 = t200 * t137 + t141 * t355 + t201 * t145;
t288 = t242 * t55 + t244 * t56;
t86 = t200 * t170 + t172 * t355 + t201 * t174;
t87 = t171 * t355 - t200 * t173 + t201 * t175;
t413 = (-t86 - t87) * t233 + (t287 + t288) * t232;
t412 = rSges(6,2) * t306 + t200 * qJD(5) - t418 * t124 + t125 * t420;
t330 = qJD(4) * t232;
t118 = (Icges(6,3) * t243 - t373) * t330 + (Icges(6,6) * t232 + t233 * t276) * qJD(3);
t120 = (-Icges(6,4) * t241 + Icges(6,6) * t243) * t330 + (Icges(6,2) * t232 + t233 * t279) * qJD(3);
t122 = (-Icges(6,1) * t241 + t372) * t330 + (Icges(6,4) * t232 + t233 * t283) * qJD(3);
t123 = (-Icges(5,1) * t241 - t375) * t330 + (Icges(5,5) * t232 + t233 * t284) * qJD(3);
t309 = t232 * t329;
t311 = t233 * t332;
t336 = qJD(3) * t232;
t358 = t232 * t241;
t411 = t175 * t311 + t118 * t358 - t233 * t120 + (-t309 + t311) * t174 + t408 * t170 + (t123 + t122) * t232 * t243 + t424 * t336;
t378 = Icges(4,4) * t232;
t286 = Icges(4,1) * t233 - t378;
t183 = Icges(4,5) * t242 + t244 * t286;
t359 = t183 * t233;
t377 = Icges(4,4) * t233;
t282 = -Icges(4,2) * t232 + t377;
t181 = Icges(4,6) * t242 + t244 * t282;
t364 = t181 * t232;
t266 = -t359 + t364;
t410 = t242 * t266;
t182 = -Icges(4,5) * t244 + t242 * t286;
t361 = t182 * t233;
t180 = -Icges(4,6) * t244 + t242 * t282;
t366 = t180 * t232;
t267 = -t361 + t366;
t409 = t244 * t267;
t236 = t242 * rSges(4,3);
t407 = -rSges(4,2) * t355 + t236;
t240 = -pkin(6) - qJ(2);
t239 = cos(pkin(8));
t229 = pkin(2) * t239 + pkin(1);
t381 = rSges(4,1) * t233;
t296 = -rSges(4,2) * t232 + t381;
t258 = -t229 - t296;
t158 = (rSges(4,3) - t240) * t244 + t258 * t242;
t354 = t233 * t244;
t185 = rSges(4,1) * t354 + t407;
t301 = t244 * t229 - t242 * t240;
t159 = t185 + t301;
t406 = t158 * t244 + t159 * t242;
t278 = Icges(4,5) * t233 - Icges(4,6) * t232;
t178 = -Icges(4,3) * t244 + t242 * t278;
t264 = rSges(3,1) * t239 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t380 = rSges(3,3) + qJ(2);
t169 = t242 * t380 + t244 * t264;
t344 = rSges(6,2) * t355 + t418 * t200 + t420 * t201;
t345 = rSges(6,2) * t357 - t423;
t405 = -t242 * t345 - t244 * t344;
t404 = 2 * m(4);
t403 = 2 * m(5);
t402 = 2 * m(6);
t401 = t242 ^ 2;
t400 = t244 ^ 2;
t399 = -t233 / 0.2e1;
t395 = -rSges(6,2) - pkin(7);
t394 = -rSges(5,3) - pkin(7);
t210 = rSges(4,1) * t232 + rSges(4,2) * t233;
t393 = m(4) * t210;
t392 = pkin(3) * t233;
t119 = (-Icges(5,5) * t241 - Icges(5,6) * t243) * t330 + (Icges(5,3) * t232 + t233 * t277) * qJD(3);
t121 = (-Icges(5,2) * t243 - t376) * t330 + (Icges(5,6) * t232 + t233 * t280) * qJD(3);
t352 = t241 * t121;
t387 = (-t173 * t334 - t119) * t233 + (-t352 + (-t173 * t243 - t175 * t241) * qJD(4)) * t232 + t411;
t384 = -rSges(6,2) * t314 + t412;
t383 = rSges(6,2) * t250 + t422;
t382 = t419 * t336;
t294 = -t199 * rSges(5,1) + t198 * rSges(5,2);
t149 = rSges(5,3) * t357 - t294;
t369 = t149 * t244;
t365 = t180 * t233;
t363 = t181 * t233;
t362 = t182 * t232;
t360 = t183 * t232;
t353 = t240 * t244;
t291 = pkin(4) * t243 + qJ(5) * t241;
t292 = rSges(6,1) * t243 + rSges(6,3) * t241;
t335 = qJD(3) * t233;
t347 = t291 * t335 + (qJD(5) * t241 + (-pkin(4) * t241 + qJ(5) * t243) * qJD(4)) * t232 + (-rSges(6,1) * t241 + rSges(6,3) * t243) * t330 + (rSges(6,2) * t232 + t233 * t292) * qJD(3);
t293 = rSges(5,1) * t243 - rSges(5,2) * t241;
t129 = (-rSges(5,1) * t241 - rSges(5,2) * t243) * t330 + (rSges(5,3) * t232 + t233 * t293) * qJD(3);
t298 = pkin(7) * t232 + t392;
t206 = t298 * qJD(3);
t346 = -t129 - t206;
t343 = -rSges(6,2) * t233 + (t291 + t292) * t232;
t177 = -rSges(5,3) * t233 + t232 * t293;
t211 = pkin(3) * t232 - pkin(7) * t233;
t342 = -t177 - t211;
t195 = t298 * t242;
t224 = pkin(3) * t354;
t196 = pkin(7) * t355 + t224;
t341 = t242 * t195 + t244 * t196;
t235 = qJD(2) * t244;
t340 = t240 * t338 + t235;
t179 = Icges(4,3) * t242 + t244 * t278;
t339 = qJD(1) * t179;
t35 = t52 * t242 - t244 * t51;
t36 = t54 * t242 - t244 * t53;
t325 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t56 * t242 - t244 * t55;
t38 = t58 * t242 - t244 * t57;
t324 = t37 / 0.2e1 + t38 / 0.2e1;
t321 = t125 * rSges(5,1) + t124 * rSges(5,2) + rSges(5,3) * t306;
t319 = -t206 - t347;
t217 = pkin(3) * t313;
t218 = pkin(7) * t306;
t318 = t242 * (pkin(7) * t250 + qJD(1) * t224 - t217) + t244 * (-pkin(7) * t314 + t218 + (-t233 * t338 - t310) * pkin(3)) + t195 * t337;
t317 = -t211 - t343;
t151 = t201 * rSges(5,1) - t200 * rSges(5,2) + rSges(5,3) * t355;
t316 = t217 + t340;
t315 = t177 * t338;
t305 = -t229 - t392;
t304 = t242 * t343;
t303 = t244 * t343;
t302 = t345 * t244;
t153 = t342 * t244;
t300 = qJD(1) * t343;
t108 = t317 * t244;
t297 = t242 * t300;
t295 = t127 * rSges(5,1) - t126 * rSges(5,2);
t281 = Icges(4,2) * t233 + t378;
t269 = -t151 * t242 + t369;
t268 = -t242 * t149 - t151 * t244;
t263 = t301 + t196;
t31 = t200 * t118 + t120 * t355 + t201 * t122 - t124 * t170 + t125 * t174 + t172 * t249;
t32 = t119 * t355 - t200 * t121 + t201 * t123 + t124 * t173 + t125 * t175 + t171 * t249;
t262 = t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 + t22 / 0.2e1;
t33 = t198 * t118 + t120 * t357 + t199 * t122 + t126 * t170 + t127 * t174 + t172 * t250;
t34 = t119 * t357 - t198 * t121 + t199 * t123 - t126 * t173 + t127 * t175 + t171 * t250;
t261 = t21 / 0.2e1 + t19 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1;
t260 = t61 / 0.2e1 + t59 / 0.2e1 + t84 / 0.2e1 + t85 / 0.2e1;
t259 = t86 / 0.2e1 + t62 / 0.2e1 + t60 / 0.2e1 + t87 / 0.2e1;
t234 = qJD(2) * t242;
t257 = -pkin(3) * t310 + t218 + t234;
t256 = qJD(3) * t210;
t254 = qJD(3) * t281;
t253 = qJD(3) * (-Icges(4,5) * t232 - Icges(4,6) * t233);
t252 = t232 * t395 + t305;
t251 = t232 * t394 + t305;
t248 = -t242 * t344 + t302;
t247 = rSges(4,2) * t314 + rSges(4,3) * t337 - t244 * t256;
t246 = t242 * t252 - t353;
t245 = t242 * t251 - t353;
t168 = -t242 * t264 + t244 * t380;
t205 = t296 * qJD(3);
t197 = t211 * t338;
t184 = -rSges(4,3) * t244 + t242 * t296;
t157 = -qJD(1) * t169 + t235;
t156 = qJD(1) * t168 + t234;
t152 = t342 * t242;
t131 = t242 * t253 + t339;
t130 = -qJD(1) * t178 + t244 * t253;
t107 = t317 * t242;
t106 = t210 * t333 + (t244 * t258 - t236) * qJD(1) + t340;
t105 = t234 + (-t353 + (-t229 - t381) * t242) * qJD(1) + t247;
t104 = t263 + t151;
t103 = t245 + t294;
t102 = -t233 * t151 - t177 * t355;
t101 = t149 * t233 + t177 * t357;
t100 = t242 * t179 - t244 * t266;
t99 = t242 * t178 - t409;
t98 = -t179 * t244 - t410;
t97 = -t178 * t244 - t242 * t267;
t92 = t269 * t232;
t91 = t263 + t344;
t90 = t246 + t423;
t89 = qJD(1) * t153 + t242 * t346;
t88 = t244 * t346 + t197 + t315;
t81 = rSges(5,3) * t250 + t295;
t79 = -rSges(5,3) * t314 + t321;
t65 = -t268 + t341;
t64 = -t232 * t303 - t233 * t344;
t63 = t232 * t304 + t233 * t345;
t50 = t251 * t337 + t307 * t394 - t295 + t316;
t49 = qJD(1) * t245 + t257 + t321;
t48 = qJD(1) * t108 + t242 * t319;
t47 = t244 * t319 + t197 + t297;
t46 = t248 * t232;
t45 = t341 - t405;
t44 = (t177 * t333 + t81) * t233 + (-qJD(3) * t149 + t242 * t129 + t177 * t337) * t232;
t43 = (-t177 * t331 - t79) * t233 + (qJD(3) * t151 - t129 * t244 + t315) * t232;
t42 = t252 * t337 + t307 * t395 + t316 - t422;
t41 = qJD(1) * t246 + t257 + t412;
t30 = t269 * t335 + (qJD(1) * t268 - t242 * t79 + t244 * t81) * t232;
t29 = t242 * t81 + t244 * t79 + (t369 + (-t151 - t196) * t242) * qJD(1) + t318;
t24 = (qJD(3) * t304 + t383) * t233 + (-qJD(3) * t345 + t242 * t347 + t244 * t300) * t232;
t23 = (-qJD(3) * t303 - t384) * t233 + (qJD(3) * t344 - t244 * t347 + t297) * t232;
t18 = -t126 * t143 + t127 * t147 + t139 * t250 - t198 * t72 + t199 * t76 + t357 * t68;
t17 = -t126 * t142 + t127 * t146 + t138 * t250 - t198 * t73 + t199 * t77 + t357 * t69;
t16 = t126 * t137 + t127 * t145 + t141 * t250 + t198 * t66 + t199 * t74 + t357 * t70;
t15 = t126 * t136 + t127 * t144 + t140 * t250 + t198 * t67 + t199 * t75 + t357 * t71;
t14 = t124 * t143 + t125 * t147 + t139 * t249 - t200 * t72 + t201 * t76 + t355 * t68;
t13 = t124 * t142 + t125 * t146 + t138 * t249 - t200 * t73 + t201 * t77 + t355 * t69;
t12 = -t124 * t137 + t125 * t145 + t141 * t249 + t200 * t66 + t201 * t74 + t355 * t70;
t11 = -t124 * t136 + t125 * t144 + t140 * t249 + t200 * t67 + t201 * t75 + t355 * t71;
t10 = t384 * t244 + t383 * t242 + (t302 + (-t196 - t344) * t242) * qJD(1) + t318;
t9 = t248 * t335 + (qJD(1) * t405 - t384 * t242 + t383 * t244) * t232;
t8 = qJD(1) * t289 - t17 * t244 + t18 * t242;
t7 = qJD(1) * t290 - t15 * t244 + t16 * t242;
t6 = qJD(1) * t287 - t13 * t244 + t14 * t242;
t5 = qJD(1) * t288 - t11 * t244 + t12 * t242;
t4 = (qJD(3) * t289 - t34) * t233 + (-qJD(1) * t36 + qJD(3) * t85 + t17 * t242 + t18 * t244) * t232;
t3 = (qJD(3) * t290 - t33) * t233 + (-qJD(1) * t35 + qJD(3) * t84 + t15 * t242 + t16 * t244) * t232;
t2 = (qJD(3) * t287 - t32) * t233 + (-qJD(1) * t38 + qJD(3) * t87 + t13 * t242 + t14 * t244) * t232;
t1 = (qJD(3) * t288 - t31) * t233 + (-qJD(1) * t37 + qJD(3) * t86 + t11 * t242 + t12 * t244) * t232;
t25 = [-t175 * t309 - t233 * t119 + 0.2e1 * m(3) * (t156 * t169 + t157 * t168) + (t41 * t91 + t42 * t90) * t402 + (t103 * t50 + t104 * t49) * t403 + (t105 * t159 + t106 * t158) * t404 - t232 * t352 + (-t281 + t286) * t336 + (Icges(4,1) * t232 + t282 + t377) * t335 - t408 * t173 + t411; m(6) * (t242 * t42 - t244 * t41 + (t242 * t91 + t244 * t90) * qJD(1)) + m(5) * (t242 * t50 - t244 * t49 + (t103 * t244 + t104 * t242) * qJD(1)) + m(4) * (qJD(1) * t406 - t105 * t244 + t242 * t106) + m(3) * (-t156 * t244 + t242 * t157 + (t168 * t244 + t169 * t242) * qJD(1)); 0; ((qJD(1) * t181 - t242 * t254) * t399 + t183 * t421 + (t366 / 0.2e1 - t361 / 0.2e1) * qJD(3) - t261) * t244 + ((-qJD(1) * t180 - t244 * t254) * t233 / 0.2e1 + t182 * t421 + (-t364 / 0.2e1 + t359 / 0.2e1) * qJD(3) + t262) * t242 + m(4) * ((-t105 * t242 - t106 * t244) * t210 - t406 * t205) + m(6) * (t107 * t41 + t108 * t42 + t47 * t90 + t48 * t91) + m(5) * (t103 * t88 + t104 * t89 + t152 * t49 + t153 * t50) + (t401 / 0.2e1 + t400 / 0.2e1) * t278 * qJD(3) + ((-t159 * t393 + t363 / 0.2e1 + t360 / 0.2e1 + t259) * t244 + (t158 * t393 + t365 / 0.2e1 + t362 / 0.2e1 + t260) * t242) * qJD(1); m(5) * (t88 * t242 - t244 * t89 + (t152 * t242 + t153 * t244) * qJD(1)) + m(6) * (t47 * t242 - t244 * t48 + (t107 * t242 + t108 * t244) * qJD(1)); (t152 * t89 + t153 * t88 + t29 * t65) * t403 + t242 * t5 + t242 * t6 + (t10 * t45 + t107 * t48 + t108 * t47) * t402 - t244 * t7 - t244 * t8 + ((t242 * t184 + t185 * t244) * ((qJD(1) * t184 + t247) * t244 + (-t242 * t256 + (-t185 + t407) * qJD(1)) * t242) + (t400 + t401) * t210 * t205) * t404 + t242 * ((t242 * t130 + (t99 + t410) * qJD(1)) * t242 + (t100 * qJD(1) + (t180 * t335 + t182 * t336) * t244 + (-t131 + (-t360 - t363) * qJD(3) + (t179 - t267) * qJD(1)) * t242) * t244) - t244 * ((t244 * t131 + (t98 + t409) * qJD(1)) * t244 + (t97 * qJD(1) + (-t181 * t335 - t183 * t336 + t339) * t242 + (-t130 + (t362 + t365) * qJD(3) - t266 * qJD(1)) * t244) * t242) + (t98 * t242 - t244 * t97 + t35 + t36) * t338 + (t100 * t242 - t244 * t99 + t37 + t38) * t337; m(6) * (t23 * t91 + t24 * t90 + t41 * t64 + t42 * t63) + m(5) * (t101 * t50 + t102 * t49 + t103 * t44 + t104 * t43) + ((t242 * t260 + t244 * t259) * qJD(3) - t387) * t233 + (t262 * t244 + t261 * t242 + (-t242 * t259 + t244 * t260) * qJD(1)) * t232 - t382; m(5) * (t44 * t242 - t244 * t43 + (t101 * t244 + t102 * t242) * qJD(1)) + m(6) * (-t23 * t244 + t24 * t242 + (t242 * t64 + t244 * t63) * qJD(1)); m(6) * (t10 * t46 + t107 * t23 + t108 * t24 + t45 * t9 + t47 * t63 + t48 * t64) + m(5) * (t101 * t88 + t102 * t89 + t152 * t43 + t153 * t44 + t29 * t92 + t30 * t65) + (-t4 / 0.2e1 - t3 / 0.2e1 + t324 * t335) * t244 + (t2 / 0.2e1 + t1 / 0.2e1 + t325 * t335) * t242 + ((-t242 * t324 + t244 * t325) * qJD(1) + (t7 + t8) * t242 / 0.2e1 + (t5 + t6) * t244 / 0.2e1 + t416 * qJD(3) / 0.2e1) * t232 + (qJD(1) * t417 + t414 * t242 + t415 * t244) * t399 + (t389 * t242 + t413 * t244) * qJD(1) / 0.2e1; (t23 * t64 + t24 * t63 + t46 * t9) * t402 + (t101 * t44 + t102 * t43 + t30 * t92) * t403 + (t387 * t233 + ((-t233 * t385 + t413) * t244 + (-t233 * t386 + t389) * t242) * qJD(3) + t382) * t233 + ((t1 + t2) * t244 + (t3 + t4) * t242 + (t415 * t242 - t414 * t244) * t233 + (t417 * t232 + t419 * t233) * qJD(3) + (t233 * t416 - t242 * t413 + t389 * t244) * qJD(1)) * t232; m(6) * (-t124 * t90 + t126 * t91 + t198 * t41 + t200 * t42); m(6) * (-t124 * t242 - t126 * t244 + (t198 * t242 + t200 * t244) * qJD(1)); m(6) * (t45 * t312 + t107 * t126 - t108 * t124 + t198 * t48 + t200 * t47 + (t10 * t241 + t328 * t45) * t232); m(6) * (t46 * t312 - t124 * t63 + t126 * t64 + t198 * t23 + t200 * t24 + (t241 * t9 + t328 * t46) * t232); (-t124 * t200 + t126 * t198 + t358 * t408) * t402;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;
