% Calculate time derivative of joint inertia matrix for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP13_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP13_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP13_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:50
% DurationCPUTime: 10.26s
% Computational Cost: add. (9771->695), mult. (26648->989), div. (0->0), fcn. (25720->6), ass. (0->333)
t243 = sin(qJ(3));
t246 = cos(qJ(3));
t242 = sin(qJ(4));
t245 = cos(qJ(4));
t278 = Icges(5,5) * t245 - Icges(5,6) * t242;
t169 = Icges(5,3) * t243 + t246 * t278;
t281 = Icges(6,4) * t245 + Icges(6,6) * t242;
t172 = Icges(6,2) * t243 + t246 * t281;
t438 = t172 + t169;
t247 = cos(qJ(1));
t355 = t247 * t245;
t330 = t243 * t355;
t342 = qJD(1) * t247;
t244 = sin(qJ(1));
t301 = qJD(1) * t243 + qJD(4);
t338 = qJD(3) * t247;
t311 = t246 * t338;
t417 = t301 * t244 - t311;
t117 = -qJD(4) * t330 + t242 * t417 - t245 * t342;
t302 = qJD(4) * t243 + qJD(1);
t364 = t242 * t247;
t118 = t245 * t417 + t302 * t364;
t359 = t244 * t245;
t362 = t243 * t247;
t194 = t242 * t362 + t359;
t385 = rSges(6,3) + qJ(5);
t433 = rSges(6,1) + pkin(4);
t437 = t194 * qJD(5) - t385 * t117 - t433 * t118;
t360 = t244 * t242;
t195 = -t330 + t360;
t436 = -t385 * t194 + t433 * t195;
t192 = t243 * t360 - t355;
t193 = t243 * t359 + t364;
t358 = t244 * t246;
t125 = Icges(5,5) * t193 - Icges(5,6) * t192 - Icges(5,3) * t358;
t129 = Icges(5,4) * t193 - Icges(5,2) * t192 - Icges(5,6) * t358;
t133 = Icges(5,1) * t193 - Icges(5,4) * t192 - Icges(5,5) * t358;
t51 = -t125 * t358 - t129 * t192 + t133 * t193;
t356 = t246 * t247;
t126 = Icges(5,5) * t195 + Icges(5,6) * t194 + Icges(5,3) * t356;
t130 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t356;
t134 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t356;
t52 = -t126 * t358 - t130 * t192 + t134 * t193;
t292 = t244 * t51 - t247 * t52;
t123 = Icges(6,5) * t193 - Icges(6,6) * t358 + Icges(6,3) * t192;
t127 = Icges(6,4) * t193 - Icges(6,2) * t358 + Icges(6,6) * t192;
t131 = Icges(6,1) * t193 - Icges(6,4) * t358 + Icges(6,5) * t192;
t49 = t123 * t192 - t127 * t358 + t131 * t193;
t124 = Icges(6,5) * t195 + Icges(6,6) * t356 - Icges(6,3) * t194;
t128 = Icges(6,4) * t195 + Icges(6,2) * t356 - Icges(6,6) * t194;
t132 = Icges(6,1) * t195 + Icges(6,4) * t356 - Icges(6,5) * t194;
t50 = t124 * t192 - t128 * t358 + t132 * t193;
t293 = t244 * t49 - t247 * t50;
t378 = Icges(6,5) * t245;
t277 = Icges(6,3) * t242 + t378;
t168 = Icges(6,6) * t243 + t246 * t277;
t379 = Icges(6,5) * t242;
t285 = Icges(6,1) * t245 + t379;
t176 = Icges(6,4) * t243 + t246 * t285;
t86 = t168 * t192 - t172 * t358 + t176 * t193;
t381 = Icges(5,4) * t245;
t282 = -Icges(5,2) * t242 + t381;
t173 = Icges(5,6) * t243 + t246 * t282;
t382 = Icges(5,4) * t242;
t286 = Icges(5,1) * t245 - t382;
t177 = Icges(5,5) * t243 + t246 * t286;
t87 = -t169 * t358 - t173 * t192 + t177 * t193;
t432 = (-t292 - t293) * t246 + (t86 + t87) * t243;
t55 = t125 * t356 + t194 * t129 + t195 * t133;
t56 = t126 * t356 + t194 * t130 + t195 * t134;
t290 = t244 * t55 - t247 * t56;
t53 = -t194 * t123 + t127 * t356 + t195 * t131;
t54 = -t194 * t124 + t128 * t356 + t195 * t132;
t291 = t244 * t53 - t247 * t54;
t88 = -t194 * t168 + t172 * t356 + t195 * t176;
t89 = t169 * t356 + t194 * t173 + t195 * t177;
t431 = (-t290 - t291) * t246 + (t88 + t89) * t243;
t269 = t168 * t242 + t176 * t245;
t430 = (-t173 * t242 + t177 * t245 + t269) * t246 + t438 * t243;
t337 = qJD(4) * t246;
t140 = (Icges(6,3) * t245 - t379) * t337 + (Icges(6,6) * t246 - t243 * t277) * qJD(3);
t141 = (-Icges(5,5) * t242 - Icges(5,6) * t245) * t337 + (Icges(5,3) * t246 - t243 * t278) * qJD(3);
t144 = (-Icges(6,4) * t242 + Icges(6,6) * t245) * t337 + (Icges(6,2) * t246 - t243 * t281) * qJD(3);
t148 = (-Icges(6,1) * t242 + t378) * t337 + (Icges(6,4) * t246 - t243 * t285) * qJD(3);
t149 = (-Icges(5,1) * t242 - t381) * t337 + (Icges(5,5) * t246 - t243 * t286) * qJD(3);
t309 = t245 * t337;
t310 = t242 * t337;
t341 = qJD(3) * t243;
t314 = t245 * t341;
t316 = t242 * t341;
t339 = qJD(3) * t246;
t365 = t242 * t246;
t429 = t140 * t365 + t168 * t309 + t173 * t316 - t176 * t310 - t177 * t314 + (t148 + t149) * t245 * t246 + t438 * t339 + (t141 + t144) * t243;
t312 = t244 * t339;
t428 = t243 * t342 + t312;
t313 = t243 * t338;
t343 = qJD(1) * t244;
t427 = t246 * t343 + t313;
t275 = t124 * t242 + t132 * t245;
t60 = t128 * t243 + t246 * t275;
t273 = t130 * t242 - t134 * t245;
t62 = t126 * t243 - t246 * t273;
t392 = t60 + t62;
t276 = t123 * t242 + t131 * t245;
t59 = t127 * t243 + t246 * t276;
t274 = t129 * t242 - t133 * t245;
t61 = t125 * t243 - t246 * t274;
t393 = t59 + t61;
t426 = -t244 * t393 + t247 * t392;
t425 = t244 * t392 + t247 * t393;
t265 = t301 * t247;
t119 = t302 * t359 + (t265 + t312) * t242;
t120 = t245 * t265 + (-t242 * t302 + t245 * t339) * t244;
t340 = qJD(3) * t244;
t315 = t243 * t340;
t317 = t246 * t342;
t254 = t315 - t317;
t67 = Icges(6,5) * t120 + Icges(6,6) * t254 + Icges(6,3) * t119;
t71 = Icges(6,4) * t120 + Icges(6,2) * t254 + Icges(6,6) * t119;
t75 = Icges(6,1) * t120 + Icges(6,4) * t254 + Icges(6,5) * t119;
t19 = (-qJD(3) * t276 + t71) * t243 + (qJD(3) * t127 + t242 * t67 + t245 * t75 + (t123 * t245 - t131 * t242) * qJD(4)) * t246;
t69 = Icges(5,5) * t120 - Icges(5,6) * t119 + Icges(5,3) * t254;
t73 = Icges(5,4) * t120 - Icges(5,2) * t119 + Icges(5,6) * t254;
t77 = Icges(5,1) * t120 - Icges(5,4) * t119 + Icges(5,5) * t254;
t21 = (qJD(3) * t274 + t69) * t243 + (qJD(3) * t125 - t242 * t73 + t245 * t77 + (-t129 * t245 - t133 * t242) * qJD(4)) * t246;
t424 = t19 + t21;
t66 = Icges(6,5) * t118 - Icges(6,6) * t427 + Icges(6,3) * t117;
t70 = Icges(6,4) * t118 - Icges(6,2) * t427 + Icges(6,6) * t117;
t74 = Icges(6,1) * t118 - Icges(6,4) * t427 + Icges(6,5) * t117;
t20 = (-qJD(3) * t275 + t70) * t243 + (qJD(3) * t128 + t242 * t66 + t245 * t74 + (t124 * t245 - t132 * t242) * qJD(4)) * t246;
t68 = Icges(5,5) * t118 - Icges(5,6) * t117 - Icges(5,3) * t427;
t72 = Icges(5,4) * t118 - Icges(5,2) * t117 - Icges(5,6) * t427;
t76 = Icges(5,1) * t118 - Icges(5,4) * t117 - Icges(5,5) * t427;
t22 = (qJD(3) * t273 + t68) * t243 + (qJD(3) * t126 - t242 * t72 + t245 * t76 + (-t130 * t245 - t134 * t242) * qJD(4)) * t246;
t423 = t20 + t22;
t422 = -rSges(6,2) * t315 - t192 * qJD(5) - t385 * t119 - t433 * t120;
t383 = Icges(4,4) * t246;
t287 = Icges(4,1) * t243 + t383;
t178 = Icges(4,5) * t247 + t244 * t287;
t369 = t178 * t243;
t384 = Icges(4,4) * t243;
t283 = Icges(4,2) * t246 + t384;
t174 = Icges(4,6) * t247 + t244 * t283;
t372 = t174 * t246;
t268 = t369 + t372;
t421 = t247 * t268;
t298 = rSges(4,1) * t243 + rSges(4,2) * t246;
t258 = t247 * t298;
t351 = rSges(6,2) * t356 + t436;
t420 = t247 * t351;
t419 = t385 * t192 + t433 * t193;
t230 = pkin(3) * t362;
t236 = t247 * qJ(2);
t418 = t230 + t236;
t326 = t428 * rSges(4,1) + rSges(4,2) * t317;
t408 = -pkin(1) - pkin(6);
t335 = -rSges(4,3) + t408;
t346 = qJ(2) * t342 + qJD(2) * t244;
t105 = (-rSges(4,2) * t341 + qJD(1) * t335) * t244 + t326 + t346;
t389 = rSges(4,2) * t243;
t213 = rSges(4,1) * t246 - t389;
t234 = qJD(2) * t247;
t106 = t234 + t213 * t338 + (t335 * t247 + (-qJ(2) - t298) * t244) * qJD(1);
t416 = -t105 * t247 + t244 * t106;
t279 = Icges(4,5) * t243 + Icges(4,6) * t246;
t415 = -Icges(4,3) * t244 + t247 * t279;
t414 = -Icges(4,6) * t244 + t247 * t283;
t413 = -Icges(4,5) * t244 + t247 * t287;
t353 = -rSges(6,2) * t358 + t419;
t412 = t244 * t351 + t247 * t353;
t411 = 2 * m(4);
t410 = 2 * m(5);
t409 = 2 * m(6);
t240 = t244 ^ 2;
t241 = t247 ^ 2;
t407 = -t243 / 0.2e1;
t404 = t246 / 0.2e1;
t401 = rSges(3,2) - pkin(1);
t400 = m(4) * t213;
t399 = pkin(3) * t243;
t391 = -rSges(6,2) * t427 - t437;
t390 = rSges(6,2) * t317 + t422;
t388 = rSges(6,2) * t246;
t387 = rSges(5,3) * t246;
t386 = t244 * rSges(4,3);
t237 = t247 * rSges(4,3);
t296 = -t195 * rSges(5,1) - t194 * rSges(5,2);
t139 = rSges(5,3) * t356 - t296;
t374 = t139 * t247;
t373 = t174 * t243;
t371 = t414 * t243;
t370 = t414 * t246;
t368 = t178 * t246;
t367 = t413 * t243;
t366 = t413 * t246;
t363 = t243 * t244;
t294 = rSges(6,1) * t245 + rSges(6,3) * t242;
t354 = (-pkin(4) * t341 + qJ(5) * t337) * t245 + (-qJ(5) * t341 + (-pkin(4) * qJD(4) + qJD(5)) * t246) * t242 + (-rSges(6,1) * t242 + rSges(6,3) * t245) * t337 + (-t243 * t294 + t388) * qJD(3);
t348 = t193 * rSges(5,1) - t192 * rSges(5,2);
t137 = -rSges(5,3) * t358 + t348;
t229 = pkin(3) * t363;
t198 = -pkin(7) * t358 + t229;
t352 = -t137 - t198;
t350 = rSges(6,2) * t243 + (pkin(4) * t245 + qJ(5) * t242 + t294) * t246;
t207 = (pkin(7) * t246 - t399) * qJD(3);
t214 = pkin(3) * t246 + pkin(7) * t243;
t347 = t244 * t207 + t214 * t342;
t345 = t247 * pkin(1) + t244 * qJ(2);
t170 = Icges(4,3) * t247 + t244 * t279;
t344 = qJD(1) * t170;
t35 = t50 * t244 + t247 * t49;
t36 = t52 * t244 + t247 * t51;
t334 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t54 * t244 + t247 * t53;
t38 = t56 * t244 + t247 * t55;
t333 = -t37 / 0.2e1 - t38 / 0.2e1;
t332 = t408 * t244;
t331 = t408 * t247;
t328 = t120 * rSges(5,1) - t119 * rSges(5,2) + rSges(5,3) * t315;
t327 = -t198 - t353;
t325 = t428 * pkin(3) + pkin(7) * t315;
t324 = pkin(3) * t311 + t427 * pkin(7);
t180 = rSges(4,1) * t363 + rSges(4,2) * t358 + t237;
t323 = t247 * pkin(6) + t345;
t322 = (-rSges(6,2) - pkin(7)) * t246;
t321 = (-rSges(5,3) - pkin(7)) * t246;
t295 = rSges(5,1) * t245 - rSges(5,2) * t242;
t182 = rSges(5,3) * t243 + t246 * t295;
t319 = t182 * t343;
t308 = -qJ(2) - t399;
t145 = (-Icges(5,2) * t245 - t382) * t337 + (Icges(5,6) * t246 - t243 * t282) * qJD(3);
t307 = t430 * t339 + (-t269 * t341 + (-t242 * t145 + (-t173 * t245 - t177 * t242) * qJD(4)) * t246 + t429) * t243;
t206 = t298 * qJD(3);
t306 = t206 * (t240 + t241);
t305 = t350 * t244;
t304 = qJD(1) * t350;
t303 = qJD(3) * t350;
t300 = t234 + t324;
t299 = t229 + t323;
t297 = t118 * rSges(5,1) - t117 * rSges(5,2);
t289 = t325 + t346;
t288 = Icges(4,1) * t246 - t384;
t284 = -Icges(4,2) * t243 + t383;
t280 = Icges(4,5) * t246 - Icges(4,6) * t243;
t272 = t137 * t247 + t139 * t244;
t267 = -t367 - t370;
t33 = t119 * t168 + t120 * t176 + t192 * t140 - t144 * t358 + t193 * t148 + t172 * t254;
t34 = -t119 * t173 + t120 * t177 - t141 * t358 - t192 * t145 + t193 * t149 + t169 * t254;
t264 = -t19 / 0.2e1 - t21 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1;
t31 = t117 * t168 + t118 * t176 - t194 * t140 + t144 * t356 + t195 * t148 - t172 * t427;
t32 = -t117 * t173 + t118 * t177 + t141 * t356 + t194 * t145 + t195 * t149 - t169 * t427;
t263 = t20 / 0.2e1 + t22 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t262 = t59 / 0.2e1 + t61 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t261 = -t60 / 0.2e1 - t62 / 0.2e1 - t89 / 0.2e1 - t88 / 0.2e1;
t260 = rSges(3,3) * t247 + t401 * t244;
t153 = (-rSges(5,1) * t242 - rSges(5,2) * t245) * t337 + (-t243 * t295 + t387) * qJD(3);
t259 = t244 * t153 + t182 * t342;
t257 = t267 * t244;
t256 = qJD(3) * t288;
t255 = qJD(3) * t284;
t252 = t247 * t322 + t332;
t251 = t247 * t321 + t332;
t249 = t244 * t354 + t247 * t304;
t205 = t244 * t214;
t200 = t214 * t343;
t199 = pkin(7) * t356 - t230;
t186 = t247 * t199;
t185 = -rSges(3,2) * t247 + t244 * rSges(3,3) + t345;
t184 = t236 + t260;
t183 = t386 - t258;
t163 = t234 + (t401 * t247 + (-rSges(3,3) - qJ(2)) * t244) * qJD(1);
t162 = qJD(1) * t260 + t346;
t161 = t323 + t180;
t160 = t244 * t335 + t236 + t258;
t159 = -pkin(7) * t317 + t325;
t156 = (-t182 - t214) * t247;
t155 = t182 * t244 + t205;
t154 = t247 * (qJD(1) * t229 - t324);
t143 = qJD(1) * t415 + t280 * t340;
t142 = -t280 * t338 + t344;
t108 = (-t214 - t350) * t247;
t107 = t205 + t305;
t104 = t244 * t321 + t299 + t348;
t103 = t251 + t296 + t418;
t102 = -t243 * t139 + t182 * t356;
t101 = t137 * t243 + t182 * t358;
t100 = -t244 * t415 - t247 * t267;
t99 = t244 * t170 - t421;
t96 = -t247 * t415 + t257;
t95 = t170 * t247 + t244 * t268;
t92 = t272 * t246;
t91 = t259 + t347;
t90 = t319 + t200 + (-t153 - t207) * t247;
t85 = t244 * t322 + t299 + t419;
t84 = t252 + t418 - t436;
t83 = -rSges(5,3) * t317 + t328;
t81 = -rSges(5,3) * t427 + t297;
t65 = t244 * t352 + t186 + t374;
t64 = -t243 * t351 + t350 * t356;
t63 = t243 * t353 + t246 * t305;
t58 = t249 + t347;
t57 = t200 + t244 * t304 + (-t207 - t354) * t247;
t48 = rSges(5,3) * t313 + (t331 + (t308 + t387) * t244) * qJD(1) - t297 + t300;
t47 = qJD(1) * t251 + t289 + t328;
t46 = t412 * t246;
t45 = t244 * t327 + t186 + t420;
t44 = (-t182 * t340 + t83) * t243 + (qJD(3) * t137 + t259) * t246;
t43 = (-t182 * t338 - t81) * t243 + (-qJD(3) * t139 + t153 * t247 - t319) * t246;
t40 = rSges(6,2) * t313 + (t331 + (t308 + t388) * t244) * qJD(1) + t300 + t437;
t39 = qJD(1) * t252 + t289 - t422;
t30 = t272 * t341 + (-t244 * t81 - t247 * t83 + (t244 * t137 - t374) * qJD(1)) * t246;
t29 = t247 * t81 + t154 + (-t159 - t83) * t244 + (t352 * t247 + (-t139 - t199) * t244) * qJD(1);
t24 = (-t244 * t303 - t390) * t243 + (qJD(3) * t353 + t249) * t246;
t23 = (-t247 * t303 - t391) * t243 + (-qJD(3) * t351 + t247 * t354 - t343 * t350) * t246;
t18 = -t119 * t130 + t120 * t134 + t126 * t254 - t192 * t72 + t193 * t76 - t358 * t68;
t17 = -t119 * t129 + t120 * t133 + t125 * t254 - t192 * t73 + t193 * t77 - t358 * t69;
t16 = t119 * t124 + t120 * t132 + t128 * t254 + t192 * t66 + t193 * t74 - t358 * t70;
t15 = t119 * t123 + t120 * t131 + t127 * t254 + t192 * t67 + t193 * t75 - t358 * t71;
t14 = -t117 * t130 + t118 * t134 - t126 * t427 + t194 * t72 + t195 * t76 + t356 * t68;
t13 = -t117 * t129 + t118 * t133 - t125 * t427 + t194 * t73 + t195 * t77 + t356 * t69;
t12 = t117 * t124 + t118 * t132 - t128 * t427 - t194 * t66 + t195 * t74 + t356 * t70;
t11 = t117 * t123 + t118 * t131 - t127 * t427 - t194 * t67 + t195 * t75 + t356 * t71;
t10 = t154 + t391 * t247 + (-t159 + t390) * t244 + (t327 * t247 + (-t199 - t351) * t244) * qJD(1);
t9 = t412 * t341 + (t390 * t247 - t391 * t244 + (t244 * t353 - t420) * qJD(1)) * t246;
t8 = -qJD(1) * t292 + t17 * t247 + t18 * t244;
t7 = -qJD(1) * t293 + t15 * t247 + t16 * t244;
t6 = -qJD(1) * t290 + t13 * t247 + t14 * t244;
t5 = -qJD(1) * t291 + t11 * t247 + t12 * t244;
t4 = (qJD(3) * t292 + t34) * t243 + (-qJD(1) * t36 + qJD(3) * t87 - t17 * t244 + t18 * t247) * t246;
t3 = (qJD(3) * t293 + t33) * t243 + (-qJD(1) * t35 + qJD(3) * t86 - t15 * t244 + t16 * t247) * t246;
t2 = (qJD(3) * t290 + t32) * t243 + (-qJD(1) * t38 + qJD(3) * t89 - t13 * t244 + t14 * t247) * t246;
t1 = (qJD(3) * t291 + t31) * t243 + (-qJD(1) * t37 + qJD(3) * t88 - t11 * t244 + t12 * t247) * t246;
t25 = [0.2e1 * m(3) * (t162 * t185 + t163 * t184) + (t105 * t161 + t106 * t160) * t411 + (t103 * t48 + t104 * t47) * t410 + (t39 * t85 + t40 * t84) * t409 - t246 * t255 - t243 * t256 + t283 * t341 - t287 * t339 - t177 * t310 - t173 * t309 - t176 * t314 - t168 * t316 - t145 * t365 + t429; m(6) * (t244 * t40 - t247 * t39 + (t244 * t85 + t247 * t84) * qJD(1)) + m(5) * (t244 * t48 - t247 * t47 + (t103 * t247 + t104 * t244) * qJD(1)) + m(4) * ((t160 * t247 + t161 * t244) * qJD(1) + t416) + m(3) * (-t162 * t247 + t244 * t163 + (t184 * t247 + t185 * t244) * qJD(1)); 0; ((qJD(1) * t414 + t244 * t255) * t407 + (qJD(1) * t413 + t244 * t256) * t404 + (-t372 / 0.2e1 - t369 / 0.2e1) * qJD(3) - t264) * t247 + ((qJD(1) * t174 - t284 * t338) * t407 + (qJD(1) * t178 - t288 * t338) * t404 + (t370 / 0.2e1 + t367 / 0.2e1) * qJD(3) + t263) * t244 + m(5) * (t103 * t91 + t104 * t90 + t155 * t48 + t156 * t47) + m(6) * (t107 * t40 + t108 * t39 + t57 * t85 + t58 * t84) + m(4) * (t416 * t213 - (t160 * t244 - t161 * t247) * t206) - (t240 / 0.2e1 + t241 / 0.2e1) * t279 * qJD(3) + ((t373 / 0.2e1 - t368 / 0.2e1 + t161 * t400 - t262) * t244 + (t160 * t400 + t371 / 0.2e1 - t366 / 0.2e1 - t261) * t247) * qJD(1); m(5) * (t91 * t244 - t247 * t90 + (t155 * t247 + t156 * t244) * qJD(1)) + m(6) * (t58 * t244 - t247 * t57 + (t107 * t247 + t108 * t244) * qJD(1)) - m(4) * t306; (t45 * t10 + t107 * t58 + t108 * t57) * t409 + (t155 * t91 + t156 * t90 + t65 * t29) * t410 + t244 * t6 + t247 * t7 + t244 * t5 + t247 * t8 + ((-t244 * t180 + t183 * t247) * (-t244 * t326 + (-t213 * t241 + t240 * t389) * qJD(3) + ((-t180 + t237) * t247 + (-t183 + t258 + t386) * t244) * qJD(1)) - t213 * t306) * t411 + t247 * ((t247 * t143 + (t96 + t421) * qJD(1)) * t247 + (-t95 * qJD(1) + (-t339 * t413 + t341 * t414) * t244 + (t142 + (t368 - t373) * qJD(3) + (-t170 + t267) * qJD(1)) * t247) * t244) + t244 * ((t244 * t142 + (-t99 + t257) * qJD(1)) * t244 + (t100 * qJD(1) + (t174 * t341 - t178 * t339 + t344) * t247 + (t143 + (t366 - t371) * qJD(3) + t268 * qJD(1)) * t244) * t247) + (-t96 * t244 - t247 * t95 - t35 - t36) * t343 + (t100 * t244 + t247 * t99 + t37 + t38) * t342; m(5) * (t101 * t47 + t102 * t48 + t103 * t43 + t104 * t44) + m(6) * (t23 * t84 + t24 * t85 + t63 * t39 + t64 * t40) + (t244 * t262 + t247 * t261) * t341 + (t263 * t247 + t264 * t244 + (t244 * t261 - t247 * t262) * qJD(1)) * t246 + t307; m(5) * (t43 * t244 - t247 * t44 + (t101 * t244 + t102 * t247) * qJD(1)) + m(6) * (t23 * t244 - t24 * t247 + (t244 * t63 + t247 * t64) * qJD(1)); m(5) * (t101 * t90 + t102 * t91 + t155 * t43 + t156 * t44 - t29 * t92 + t30 * t65) + m(6) * (-t46 * t10 + t107 * t23 + t108 * t24 + t9 * t45 + t63 * t57 + t64 * t58) + (t4 / 0.2e1 + t3 / 0.2e1 + t333 * t341 + t431 * qJD(1) / 0.2e1) * t247 + (t2 / 0.2e1 + t1 / 0.2e1 + t334 * t341 - t432 * qJD(1) / 0.2e1) * t244 + ((t244 * t333 - t247 * t334) * qJD(1) - (t7 + t8) * t244 / 0.2e1 + (t6 + t5) * t247 / 0.2e1 + t425 * qJD(3) / 0.2e1) * t246 + (t426 * qJD(1) + t423 * t244 + t424 * t247) * t243 / 0.2e1; (t23 * t64 + t24 * t63 - t46 * t9) * t409 + (t101 * t44 + t102 * t43 - t30 * t92) * t410 + (((-t243 * t392 - t431) * t247 + (t243 * t393 + t432) * t244) * qJD(3) + t307) * t243 + ((t1 + t2) * t247 + (-t3 - t4) * t244 + (-t244 * t424 + t423 * t247) * t243 + (t430 * t243 + t426 * t246) * qJD(3) + (-t425 * t243 - t244 * t431 - t247 * t432) * qJD(1)) * t246; m(6) * (t117 * t85 + t119 * t84 + t192 * t40 - t194 * t39); m(6) * (-t117 * t247 + t119 * t244 + (t192 * t247 - t194 * t244) * qJD(1)); m(6) * (t45 * t309 + t107 * t119 + t108 * t117 + t192 * t58 - t194 * t57 + (t10 * t246 - t341 * t45) * t242); m(6) * (-t46 * t309 + t117 * t63 + t119 * t64 + t192 * t23 - t194 * t24 + (t246 * t9 + t341 * t46) * t242); (-t117 * t194 + t119 * t192 + (t309 - t316) * t365) * t409;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;
