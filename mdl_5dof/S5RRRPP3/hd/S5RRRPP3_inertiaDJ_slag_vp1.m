% Calculate time derivative of joint inertia matrix for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:22
% DurationCPUTime: 8.26s
% Computational Cost: add. (6639->452), mult. (7464->613), div. (0->0), fcn. (5582->6), ass. (0->246)
t218 = sin(qJ(3));
t220 = cos(qJ(3));
t343 = Icges(6,6) * t220;
t347 = Icges(4,4) * t220;
t394 = -t343 + t347 + (Icges(4,1) + Icges(6,3)) * t218;
t344 = Icges(6,6) * t218;
t346 = Icges(5,6) * t218;
t393 = t344 - t346 + (-Icges(6,2) - Icges(5,3)) * t220;
t265 = Icges(6,3) * t220 + t344;
t348 = Icges(4,4) * t218;
t279 = Icges(4,1) * t220 - t348;
t392 = t265 + t279;
t345 = Icges(5,6) * t220;
t268 = -Icges(5,3) * t218 + t345;
t269 = Icges(6,2) * t218 + t343;
t391 = t268 - t269;
t317 = qJD(3) * t218;
t316 = qJD(3) * t220;
t174 = Icges(4,2) * t220 + t348;
t271 = Icges(5,2) * t218 + t345;
t382 = (t271 + t394) * t220 + (-t174 + t393) * t218;
t360 = rSges(6,1) + pkin(4);
t217 = qJ(1) + qJ(2);
t211 = sin(t217);
t389 = t211 * t360;
t212 = cos(t217);
t214 = qJD(1) + qJD(2);
t338 = t212 * t214;
t388 = rSges(6,3) + qJ(5);
t302 = t212 * t317;
t334 = t214 * t220;
t387 = t211 * t334 + t302;
t335 = t214 * t218;
t386 = t211 * t316 + t212 * t335;
t171 = Icges(4,5) * t218 + Icges(4,6) * t220;
t275 = Icges(6,4) * t220 - Icges(6,5) * t218;
t276 = Icges(5,4) * t218 + Icges(5,5) * t220;
t384 = t171 - t275 - t276;
t273 = Icges(4,5) * t220 - Icges(4,6) * t218;
t274 = Icges(6,4) * t218 + Icges(6,5) * t220;
t277 = Icges(5,4) * t220 - Icges(5,5) * t218;
t381 = t382 * t214 + (-t273 - t274 + t277) * qJD(3);
t322 = t211 ^ 2 + t212 ^ 2;
t380 = -0.1e1 + t322;
t244 = t268 * t212;
t105 = Icges(5,5) * t211 - t244;
t272 = Icges(5,2) * t220 - t346;
t246 = t272 * t212;
t109 = Icges(5,4) * t211 - t246;
t257 = t105 * t218 - t109 * t220;
t379 = t211 * t257;
t243 = t265 * t212;
t103 = Icges(6,5) * t211 + t243;
t245 = t269 * t212;
t107 = Icges(6,4) * t211 + t245;
t262 = t103 * t220 + t107 * t218;
t378 = t211 * t262;
t278 = -Icges(4,2) * t218 + t347;
t250 = t278 * t212;
t100 = Icges(4,6) * t211 + t250;
t251 = t279 * t212;
t102 = Icges(4,5) * t211 + t251;
t263 = t100 * t218 - t102 * t220;
t377 = t211 * t263;
t106 = -Icges(5,5) * t212 - t211 * t268;
t110 = -Icges(5,4) * t212 - t211 * t272;
t255 = t106 * t218 - t110 * t220;
t376 = t212 * t255;
t104 = -Icges(6,5) * t212 + t211 * t265;
t108 = -Icges(6,4) * t212 + t211 * t269;
t260 = t104 * t220 + t108 * t218;
t375 = t212 * t260;
t101 = -Icges(4,5) * t212 + t211 * t279;
t99 = -Icges(4,6) * t212 + t211 * t278;
t281 = t101 * t220 - t218 * t99;
t374 = t212 * t281;
t373 = t391 * qJD(3);
t372 = t360 * t212;
t368 = 2 * m(3);
t367 = 2 * m(4);
t366 = 2 * m(5);
t365 = 2 * m(6);
t364 = m(5) / 0.2e1;
t363 = m(6) / 0.2e1;
t355 = rSges(4,1) * t220;
t157 = (-rSges(4,2) * t218 + t355) * qJD(3);
t359 = m(4) * t157;
t184 = rSges(4,1) * t218 + rSges(4,2) * t220;
t358 = m(4) * t184;
t219 = sin(qJ(1));
t357 = pkin(1) * t219;
t356 = pkin(3) * t220;
t354 = rSges(5,2) * t218;
t353 = pkin(1) * qJD(1);
t199 = t211 * rSges(5,1);
t197 = t211 * rSges(4,3);
t182 = pkin(3) * t218 - qJ(4) * t220;
t297 = -rSges(6,2) * t220 + t218 * t388;
t293 = -t182 - t297;
t88 = t293 * t212;
t352 = t214 * t88;
t351 = -rSges(6,2) - qJ(4);
t350 = -rSges(5,3) - qJ(4);
t290 = qJ(4) * t218 + t356;
t125 = t290 * t211;
t336 = t212 * t220;
t337 = t212 * t218;
t126 = pkin(3) * t336 + qJ(4) * t337;
t331 = t211 * t125 + t212 * t126;
t291 = rSges(6,2) * t218 + rSges(6,3) * t220;
t339 = t211 * t220;
t332 = qJ(5) * t339 + t211 * t291 - t372;
t333 = rSges(6,2) * t337 + t388 * t336 + t389;
t19 = t211 * t332 + t212 * t333 + t331;
t342 = qJD(3) * t19;
t341 = t211 * t214;
t340 = t211 * t218;
t130 = qJD(3) * t290 - qJD(4) * t220;
t330 = -t130 - (-rSges(5,2) * t220 + rSges(5,3) * t218) * qJD(3);
t310 = t211 * t335;
t329 = rSges(4,2) * t310 + rSges(4,3) * t338;
t304 = t211 * t317;
t328 = t388 * t304;
t301 = t212 * t316;
t315 = qJD(4) * t218;
t327 = qJ(4) * t301 + t212 * t315;
t292 = rSges(5,3) * t220 + t354;
t326 = -t182 + t292;
t325 = rSges(4,2) * t340 + t212 * rSges(4,3);
t324 = t212 * rSges(5,1) + rSges(5,2) * t339;
t323 = t212 * pkin(2) + t211 * pkin(7);
t319 = qJD(3) * t211;
t318 = qJD(3) * t212;
t314 = -pkin(3) - t388;
t166 = pkin(3) * t304;
t307 = t212 * t334;
t313 = t211 * (pkin(3) * t307 + qJ(4) * t386 + t211 * t315 - t166) + t212 * (-pkin(3) * t387 - qJ(4) * t310 + t327) + t125 * t338;
t312 = t219 * t353;
t221 = cos(qJ(1));
t311 = t221 * t353;
t306 = -rSges(4,1) * t304 - rSges(4,2) * t386;
t180 = pkin(7) * t338;
t305 = t180 + t327;
t298 = -pkin(2) - t355;
t141 = t212 * rSges(3,1) - rSges(3,2) * t211;
t93 = t326 * t212;
t296 = rSges(5,1) * t338 + t387 * rSges(5,2) + rSges(5,3) * t301;
t295 = rSges(6,2) * t301 + qJD(5) * t336 + t360 * t338;
t294 = t323 + t126;
t124 = -rSges(3,1) * t338 + rSges(3,2) * t341;
t140 = -rSges(3,1) * t211 - rSges(3,2) * t212;
t206 = t212 * pkin(7);
t237 = t218 * t351 + t220 * t314 - pkin(2);
t227 = t237 * t211;
t47 = t206 + t227 + t372;
t45 = t47 - t357;
t213 = t221 * pkin(1);
t48 = t294 + t333;
t46 = t213 + t48;
t289 = t211 * t46 + t212 * t45;
t288 = t211 * t48 + t212 * t47;
t280 = t101 * t218 + t220 * t99;
t264 = t100 * t220 + t102 * t218;
t261 = t103 * t218 - t107 * t220;
t259 = t104 * t218 - t108 * t220;
t258 = t105 * t220 + t109 * t218;
t256 = -t106 * t220 - t110 * t218;
t116 = rSges(4,1) * t336 - rSges(4,2) * t337 + t197;
t118 = -rSges(5,2) * t336 + rSges(5,3) * t337 + t199;
t123 = t140 * t214;
t249 = t277 * t212;
t248 = t274 * t212;
t247 = t273 * t212;
t242 = t218 * t350 - pkin(2) - t356;
t241 = (t278 + t394) * t316 + (t392 + t393) * t317;
t90 = t116 + t323;
t89 = t211 * t298 + t206 + t325;
t239 = -qJ(5) * t316 - t291 * qJD(3) - qJD(5) * t218 - t130;
t238 = t242 * t211;
t84 = t118 + t294;
t236 = Icges(5,1) * t214 + qJD(3) * t276;
t235 = Icges(6,1) * t214 + qJD(3) * t275;
t228 = Icges(4,3) * t214 - qJD(3) * t171;
t83 = t206 + t238 + t324;
t44 = (t298 * t212 + (-rSges(4,3) - pkin(7)) * t211) * t214 - t306;
t145 = t272 * qJD(3);
t223 = (-qJD(3) * t174 + t145) * t218 + (qJD(3) * t271 + t373) * t220 + t241;
t222 = (t211 * t384 + t212 * t382 + t261 + t264) * t338 / 0.2e1 + (-t381 * t211 + (t257 + t262 - t263) * qJD(3) + ((-t278 - t391) * t220 + (-t272 - t392) * t218) * t341) * t211 / 0.2e1 - (t381 * t212 + (t255 + t260 + t281) * qJD(3)) * t212 / 0.2e1 + (t211 * t382 - t212 * t384 + t256 + t259 + t280) * t341 / 0.2e1 - (t258 + (t244 - t245 + t250) * t220 + (t243 + t246 + t251) * t218) * t338 / 0.2e1;
t43 = -rSges(4,1) * t387 - rSges(4,2) * t301 - pkin(2) * t341 + t180 + t329;
t22 = -pkin(3) * t302 + t214 * t238 + t296 + t305;
t11 = t214 * t227 + t302 * t314 + t295 + t305;
t154 = rSges(5,2) * t307;
t23 = t154 + t166 + t242 * t338 + (-t315 + (-rSges(5,1) - pkin(7)) * t214 + (t220 * t350 - t354) * qJD(3)) * t211;
t12 = t166 + (-t315 + (qJD(3) * t351 - qJD(5)) * t220) * t211 + ((-pkin(7) - t360) * t211 + t237 * t212) * t214 + t328;
t129 = t141 + t213;
t128 = t140 - t357;
t127 = t182 * t341;
t120 = rSges(5,3) * t340 - t324;
t115 = rSges(4,1) * t339 - t325;
t114 = -Icges(5,1) * t212 - t211 * t277;
t113 = Icges(5,1) * t211 - t249;
t112 = -Icges(6,1) * t212 + t211 * t274;
t111 = Icges(6,1) * t211 + t248;
t98 = Icges(4,3) * t211 + t247;
t97 = -Icges(4,3) * t212 + t211 * t273;
t96 = t124 - t311;
t95 = t123 - t312;
t92 = t326 * t211;
t91 = t380 * t218 * t316;
t87 = t293 * t211;
t86 = t213 + t90;
t85 = t89 - t357;
t82 = t212 * t236 + t277 * t341;
t81 = t211 * t236 - t214 * t249;
t80 = t212 * t235 - t274 * t341;
t79 = t211 * t235 + t214 * t248;
t66 = t211 * t228 + t214 * t247;
t65 = t212 * t228 - t273 * t341;
t64 = t213 + t84;
t63 = t83 - t357;
t42 = t211 * t330 + t214 * t93;
t41 = t212 * t330 - t292 * t341 + t127;
t40 = t44 - t311;
t39 = t43 - t312;
t38 = -t114 * t212 + t211 * t255;
t37 = -t113 * t212 + t379;
t36 = -t112 * t212 + t211 * t260;
t35 = -t111 * t212 + t378;
t34 = t114 * t211 + t376;
t33 = t113 * t211 + t212 * t257;
t32 = t112 * t211 + t375;
t31 = t111 * t211 + t212 * t262;
t30 = t211 * t98 - t212 * t263;
t29 = t211 * t97 + t374;
t28 = -t212 * t98 - t377;
t27 = t211 * t281 - t212 * t97;
t26 = t118 * t212 + t120 * t211 + t331;
t25 = t211 * t239 + t352;
t24 = t212 * t239 + t297 * t341 + t127;
t21 = t23 - t311;
t20 = t22 - t312;
t10 = t12 - t311;
t9 = t11 - t312;
t2 = (t120 * t214 + t296) * t212 + (-t154 + t292 * t319 + (-t118 - t126 + t199) * t214) * t211 + t313;
t1 = (t214 * t332 - t302 * t388 + t295) * t212 + ((rSges(6,2) * qJD(3) + qJD(5)) * t339 + (-t126 - t333 + t389) * t214 - t328) * t211 + t313;
t3 = [(t10 * t45 + t46 * t9) * t365 - t174 * t317 + t271 * t316 + t218 * t145 + (t20 * t64 + t21 * t63) * t366 + (t39 * t86 + t40 * t85) * t367 + (t128 * t96 + t129 * t95) * t368 + t241 + t373 * t220; m(6) * (t10 * t47 + t11 * t46 + t12 * t45 + t48 * t9) + m(5) * (t20 * t84 + t21 * t83 + t22 * t64 + t23 * t63) + m(4) * (t39 * t90 + t40 * t89 + t43 * t86 + t44 * t85) + m(3) * (t123 * t129 + t124 * t128 + t140 * t96 + t141 * t95) + t223; (t22 * t84 + t23 * t83) * t366 + (t11 * t48 + t12 * t47) * t365 + (t43 * t90 + t44 * t89) * t367 + (t123 * t141 + t124 * t140) * t368 + t223; t222 + m(5) * (t20 * t92 + t21 * t93 + t41 * t63 + t42 * t64) + m(6) * (t10 * t88 + t24 * t45 + t25 * t46 + t87 * t9) + (-t211 * t86 - t212 * t85) * t359 + ((-t214 * t86 - t40) * t212 + (t214 * t85 - t39) * t211) * t358; t222 + m(5) * (t22 * t92 + t23 * t93 + t41 * t83 + t42 * t84) + m(6) * (t11 * t87 + t12 * t88 + t24 * t47 + t25 * t48) + (-t211 * t90 - t212 * t89) * t359 + ((-t214 * t90 - t44) * t212 + (t214 * t89 - t43) * t211) * t358; (t1 * t19 + t24 * t88 + t25 * t87) * t365 + (t2 * t26 + t41 * t93 + t42 * t92) * t366 + t211 * ((t211 * t80 + (t32 - t378) * t214) * t211 + (t31 * t214 + (t104 * t317 - t108 * t316) * t212 + (-t261 * qJD(3) + t214 * t260 - t79) * t211) * t212) + t211 * ((t211 * t82 + (t34 - t379) * t214) * t211 + (t33 * t214 + (-t106 * t316 - t110 * t317) * t212 + (t258 * qJD(3) + t214 * t255 - t81) * t211) * t212) + ((t115 * t211 + t116 * t212) * (((-t116 + t197) * t214 + t306) * t211 + (t214 * t115 - t184 * t318 + t329) * t212) + t322 * t184 * t157) * t367 - t212 * ((t212 * t81 + (t37 - t376) * t214) * t212 + (t38 * t214 + (t105 * t316 + t109 * t317) * t211 + (t256 * qJD(3) + t214 * t257 - t82) * t212) * t211) + t211 * ((t211 * t65 + (t29 + t377) * t214) * t211 + (t30 * t214 + (t101 * t317 + t316 * t99) * t212 + (-t264 * qJD(3) + t214 * t281 - t66) * t211) * t212) - t212 * ((t212 * t79 + (t35 - t375) * t214) * t212 + (t36 * t214 + (-t103 * t317 + t107 * t316) * t211 + (t259 * qJD(3) + t214 * t262 - t80) * t212) * t211) - t212 * ((t212 * t66 + (t28 - t374) * t214) * t212 + (t27 * t214 + (-t100 * t316 - t102 * t317) * t211 + (t280 * qJD(3) - t214 * t263 - t65) * t212) * t211) + ((-t27 - t36 - t38) * t212 + (t28 + t35 + t37) * t211) * t341 + ((-t29 - t32 - t34) * t212 + (t30 + t31 + t33) * t211) * t338; 0.2e1 * (t289 * t363 + (t211 * t64 + t212 * t63) * t364) * t316 + 0.2e1 * ((t10 * t212 + t211 * t9 + t338 * t46 - t341 * t45) * t363 + (t20 * t211 + t21 * t212 + t338 * t64 - t341 * t63) * t364) * t218; 0.2e1 * ((t211 * t84 + t212 * t83) * t364 + t288 * t363) * t316 + 0.2e1 * ((t211 * t22 + t212 * t23 + t338 * t84 - t341 * t83) * t364 + (t11 * t211 + t12 * t212 + t338 * t48 - t341 * t47) * t363) * t218; 0.2e1 * ((t318 * t88 + t319 * t87 - t1) * t363 + (t318 * t93 + t319 * t92 - t2) * t364) * t220 + 0.2e1 * ((t211 * t25 + t212 * t24 + t338 * t87 - t341 * t88 + t342) * t363 + (qJD(3) * t26 + t211 * t42 + t212 * t41 + t338 * t92 - t341 * t93) * t364) * t218; 0.4e1 * (t364 + t363) * t91; m(6) * (-t289 * t317 + ((t214 * t46 + t10) * t212 + (-t214 * t45 + t9) * t211) * t220); m(6) * (-t288 * t317 + ((t214 * t48 + t12) * t212 + (-t214 * t47 + t11) * t211) * t220); m(6) * ((t1 + (-t211 * t87 - t212 * t88) * qJD(3)) * t218 + (t342 + (t214 * t87 + t24) * t212 + (t25 - t352) * t211) * t220); m(6) * t380 * (-t218 ^ 2 + t220 ^ 2) * qJD(3); -0.2e1 * m(6) * t91;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
