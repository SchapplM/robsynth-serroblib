% Calculate joint inertia matrix for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR13_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR13_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:09
% EndTime: 2019-03-09 23:51:26
% DurationCPUTime: 8.05s
% Computational Cost: add. (32296->699), mult. (83343->955), div. (0->0), fcn. (108924->12), ass. (0->317)
t295 = cos(pkin(6));
t300 = sin(qJ(1));
t302 = cos(qJ(2));
t360 = t300 * t302;
t299 = sin(qJ(2));
t303 = cos(qJ(1));
t361 = t299 * t303;
t281 = t295 * t361 + t360;
t298 = sin(qJ(3));
t294 = sin(pkin(6));
t371 = cos(qJ(3));
t336 = t294 * t371;
t260 = t281 * t298 + t303 * t336;
t359 = t302 * t303;
t362 = t299 * t300;
t283 = -t295 * t362 + t359;
t262 = t283 * t298 - t300 * t336;
t366 = t294 * t299;
t278 = -t295 * t371 + t298 * t366;
t279 = t295 * t298 + t299 * t336;
t297 = sin(qJ(4));
t364 = t294 * t302;
t370 = cos(qJ(4));
t258 = t279 * t297 + t364 * t370;
t259 = t279 * t370 - t297 * t364;
t296 = sin(qJ(6));
t301 = cos(qJ(6));
t210 = t258 * t301 - t259 * t296;
t211 = t258 * t296 + t259 * t301;
t132 = Icges(7,5) * t211 + Icges(7,6) * t210 - Icges(7,3) * t278;
t133 = Icges(7,4) * t211 + Icges(7,2) * t210 - Icges(7,6) * t278;
t134 = Icges(7,1) * t211 + Icges(7,4) * t210 - Icges(7,5) * t278;
t59 = -t278 * t132 + t210 * t133 + t211 * t134;
t368 = t59 * t278;
t363 = t294 * t303;
t261 = t281 * t371 - t298 * t363;
t280 = -t295 * t359 + t362;
t228 = t261 * t297 - t280 * t370;
t229 = t261 * t370 + t280 * t297;
t168 = t228 * t301 - t229 * t296;
t169 = t228 * t296 + t229 * t301;
t111 = Icges(7,5) * t169 + Icges(7,6) * t168 - Icges(7,3) * t260;
t113 = Icges(7,4) * t169 + Icges(7,2) * t168 - Icges(7,6) * t260;
t115 = Icges(7,1) * t169 + Icges(7,4) * t168 - Icges(7,5) * t260;
t38 = -t111 * t278 + t113 * t210 + t115 * t211;
t365 = t294 * t300;
t263 = t283 * t371 + t298 * t365;
t282 = t295 * t360 + t361;
t230 = t263 * t297 - t282 * t370;
t231 = t263 * t370 + t282 * t297;
t170 = t230 * t301 - t231 * t296;
t171 = t230 * t296 + t231 * t301;
t112 = Icges(7,5) * t171 + Icges(7,6) * t170 - Icges(7,3) * t262;
t114 = Icges(7,4) * t171 + Icges(7,2) * t170 - Icges(7,6) * t262;
t116 = Icges(7,1) * t171 + Icges(7,4) * t170 - Icges(7,5) * t262;
t39 = -t112 * t278 + t114 * t210 + t116 * t211;
t9 = t38 * t260 + t39 * t262 + t368;
t375 = -t9 / 0.2e1;
t381 = t9 / 0.2e1;
t34 = -t111 * t260 + t113 * t168 + t115 * t169;
t35 = -t112 * t260 + t114 * t168 + t116 * t169;
t53 = -t132 * t260 + t133 * t168 + t134 * t169;
t1 = t260 * t34 + t262 * t35 + t278 * t53;
t377 = -t1 / 0.2e1;
t380 = t1 / 0.2e1;
t36 = -t111 * t262 + t113 * t170 + t115 * t171;
t37 = -t112 * t262 + t114 * t170 + t116 * t171;
t54 = -t132 * t262 + t133 * t170 + t134 * t171;
t3 = t260 * t36 + t262 * t37 + t278 * t54;
t376 = -t3 / 0.2e1;
t118 = t171 * rSges(7,1) + t170 * rSges(7,2) - t262 * rSges(7,3);
t379 = t231 * pkin(5) + t118;
t187 = Icges(6,5) * t259 + Icges(6,6) * t278 + Icges(6,3) * t258;
t189 = Icges(6,4) * t259 + Icges(6,2) * t278 + Icges(6,6) * t258;
t191 = Icges(6,1) * t259 + Icges(6,4) * t278 + Icges(6,5) * t258;
t101 = t258 * t187 + t278 * t189 + t259 * t191;
t188 = Icges(5,5) * t259 - Icges(5,6) * t258 + Icges(5,3) * t278;
t190 = Icges(5,4) * t259 - Icges(5,2) * t258 + Icges(5,6) * t278;
t192 = Icges(5,1) * t259 - Icges(5,4) * t258 + Icges(5,5) * t278;
t102 = t278 * t188 - t258 * t190 + t259 * t192;
t233 = Icges(4,5) * t279 - Icges(4,6) * t278 - Icges(4,3) * t364;
t234 = Icges(4,4) * t279 - Icges(4,2) * t278 - Icges(4,6) * t364;
t235 = Icges(4,1) * t279 - Icges(4,4) * t278 - Icges(4,5) * t364;
t130 = -t233 * t364 - t278 * t234 + t279 * t235;
t378 = -t101 - t102 - t130 - t59;
t374 = -t260 / 0.2e1;
t373 = -t262 / 0.2e1;
t372 = -t278 / 0.2e1;
t369 = pkin(10) * t262;
t237 = Icges(3,5) * t281 - Icges(3,6) * t280 - Icges(3,3) * t363;
t367 = t237 * t303;
t319 = -rSges(7,1) * t169 - rSges(7,2) * t168;
t117 = -rSges(7,3) * t260 - t319;
t358 = pkin(5) * t229 - pkin(11) * t260 + t117;
t357 = -pkin(11) * t262 + t379;
t135 = rSges(7,1) * t211 + rSges(7,2) * t210 - rSges(7,3) * t278;
t356 = pkin(5) * t259 - pkin(11) * t278 + t135;
t154 = t231 * rSges(6,1) + t262 * rSges(6,2) + t230 * rSges(6,3);
t177 = t231 * pkin(4) + qJ(5) * t230;
t355 = -t154 - t177;
t155 = t231 * rSges(5,1) - t230 * rSges(5,2) + t262 * rSges(5,3);
t257 = t263 * pkin(3);
t215 = t257 + t369;
t354 = -t155 - t215;
t220 = t228 * qJ(5);
t176 = pkin(4) * t229 + t220;
t214 = pkin(3) * t261 + t260 * pkin(10);
t204 = t282 * t214;
t353 = t282 * t176 + t204;
t193 = rSges(6,1) * t259 + rSges(6,2) * t278 + rSges(6,3) * t258;
t213 = pkin(4) * t259 + qJ(5) * t258;
t352 = -t193 - t213;
t194 = rSges(5,1) * t259 - rSges(5,2) * t258 + rSges(5,3) * t278;
t248 = pkin(3) * t279 + pkin(10) * t278;
t351 = -t194 - t248;
t350 = t214 * t364 + t280 * t248;
t250 = t283 * pkin(2) + pkin(9) * t282;
t247 = t295 * t250;
t349 = t295 * t215 + t247;
t249 = pkin(2) * t281 + t280 * pkin(9);
t348 = -t214 - t249;
t347 = t249 * t365 + t250 * t363;
t346 = t303 * pkin(1) + pkin(8) * t365;
t345 = -t38 / 0.2e1 - t53 / 0.2e1;
t344 = -t39 / 0.2e1 - t54 / 0.2e1;
t343 = -t177 - t357;
t342 = -t213 - t356;
t341 = -t215 + t355;
t340 = t295 * t177 + t349;
t339 = -t176 + t348;
t338 = -t248 + t352;
t202 = t263 * rSges(4,1) - t262 * rSges(4,2) + t282 * rSges(4,3);
t267 = Icges(3,3) * t295 + (Icges(3,5) * t299 + Icges(3,6) * t302) * t294;
t268 = Icges(3,6) * t295 + (Icges(3,4) * t299 + Icges(3,2) * t302) * t294;
t269 = Icges(3,5) * t295 + (Icges(3,1) * t299 + Icges(3,4) * t302) * t294;
t337 = t295 * t267 + t268 * t364 + t269 * t366;
t244 = t283 * rSges(3,1) - t282 * rSges(3,2) + rSges(3,3) * t365;
t335 = -t300 * pkin(1) + pkin(8) * t363;
t236 = rSges(4,1) * t279 - rSges(4,2) * t278 - rSges(4,3) * t364;
t284 = (pkin(2) * t299 - pkin(9) * t302) * t294;
t334 = t294 * (-t236 - t284);
t333 = -t215 + t343;
t332 = -t248 + t342;
t331 = t176 * t364 + t280 * t213 + t350;
t330 = t214 * t365 + t215 * t363 + t347;
t140 = Icges(6,5) * t229 + Icges(6,6) * t260 + Icges(6,3) * t228;
t144 = Icges(6,4) * t229 + Icges(6,2) * t260 + Icges(6,6) * t228;
t148 = Icges(6,1) * t229 + Icges(6,4) * t260 + Icges(6,5) * t228;
t68 = t140 * t228 + t144 * t260 + t148 * t229;
t141 = Icges(6,5) * t231 + Icges(6,6) * t262 + Icges(6,3) * t230;
t145 = Icges(6,4) * t231 + Icges(6,2) * t262 + Icges(6,6) * t230;
t149 = Icges(6,1) * t231 + Icges(6,4) * t262 + Icges(6,5) * t230;
t69 = t141 * t228 + t145 * t260 + t149 * t229;
t88 = t187 * t228 + t189 * t260 + t191 * t229;
t13 = t260 * t68 + t262 * t69 + t278 * t88;
t142 = Icges(5,5) * t229 - Icges(5,6) * t228 + Icges(5,3) * t260;
t146 = Icges(5,4) * t229 - Icges(5,2) * t228 + Icges(5,6) * t260;
t150 = Icges(5,1) * t229 - Icges(5,4) * t228 + Icges(5,5) * t260;
t70 = t142 * t260 - t146 * t228 + t150 * t229;
t143 = Icges(5,5) * t231 - Icges(5,6) * t230 + Icges(5,3) * t262;
t147 = Icges(5,4) * t231 - Icges(5,2) * t230 + Icges(5,6) * t262;
t151 = Icges(5,1) * t231 - Icges(5,4) * t230 + Icges(5,5) * t262;
t71 = t143 * t260 - t147 * t228 + t151 * t229;
t89 = t188 * t260 - t190 * t228 + t192 * t229;
t14 = t260 * t70 + t262 * t71 + t278 * t89;
t329 = t380 + t14 / 0.2e1 + t13 / 0.2e1;
t72 = t140 * t230 + t144 * t262 + t148 * t231;
t73 = t141 * t230 + t145 * t262 + t149 * t231;
t90 = t187 * t230 + t189 * t262 + t191 * t231;
t15 = t260 * t72 + t262 * t73 + t278 * t90;
t74 = t142 * t262 - t146 * t230 + t150 * t231;
t75 = t143 * t262 - t147 * t230 + t151 * t231;
t91 = t188 * t262 - t190 * t230 + t192 * t231;
t16 = t260 * t74 + t262 * t75 + t278 * t91;
t328 = t3 / 0.2e1 + t16 / 0.2e1 + t15 / 0.2e1;
t17 = t280 * t68 + t282 * t69 - t364 * t88;
t18 = t280 * t70 + t282 * t71 - t364 * t89;
t5 = t280 * t34 + t282 * t35 - t364 * t53;
t327 = t5 / 0.2e1 + t18 / 0.2e1 + t17 / 0.2e1;
t19 = t280 * t72 + t282 * t73 - t364 * t90;
t20 = t280 * t74 + t282 * t75 - t364 * t91;
t6 = t280 * t36 + t282 * t37 - t364 * t54;
t326 = t6 / 0.2e1 + t20 / 0.2e1 + t19 / 0.2e1;
t21 = t88 * t295 + (t300 * t69 - t303 * t68) * t294;
t22 = t89 * t295 + (t300 * t71 - t303 * t70) * t294;
t7 = t53 * t295 + (t300 * t35 - t303 * t34) * t294;
t325 = t7 / 0.2e1 + t21 / 0.2e1 + t22 / 0.2e1;
t23 = t90 * t295 + (t300 * t73 - t303 * t72) * t294;
t24 = t91 * t295 + (t300 * t75 - t303 * t74) * t294;
t8 = t54 * t295 + (t300 * t37 - t303 * t36) * t294;
t324 = t8 / 0.2e1 + t23 / 0.2e1 + t24 / 0.2e1;
t77 = t140 * t258 + t144 * t278 + t148 * t259;
t78 = t141 * t258 + t145 * t278 + t149 * t259;
t97 = t101 * t278;
t25 = t77 * t260 + t78 * t262 + t97;
t79 = t142 * t278 - t146 * t258 + t150 * t259;
t80 = t143 * t278 - t147 * t258 + t151 * t259;
t98 = t102 * t278;
t26 = t79 * t260 + t80 * t262 + t98;
t323 = t381 + t25 / 0.2e1 + t26 / 0.2e1;
t11 = t38 * t280 + t39 * t282 - t364 * t59;
t27 = -t101 * t364 + t77 * t280 + t78 * t282;
t28 = -t102 * t364 + t79 * t280 + t80 * t282;
t322 = t11 / 0.2e1 + t28 / 0.2e1 + t27 / 0.2e1;
t58 = t59 * t295;
t12 = t58 + (t39 * t300 - t38 * t303) * t294;
t99 = t101 * t295;
t29 = t99 + (t78 * t300 - t77 * t303) * t294;
t100 = t102 * t295;
t30 = t100 + (t80 * t300 - t79 * t303) * t294;
t321 = t12 / 0.2e1 + t29 / 0.2e1 + t30 / 0.2e1;
t320 = t294 * (-t284 + t351);
t318 = -rSges(6,2) * t260 - rSges(6,3) * t228;
t316 = t250 + t346;
t315 = t294 * (-t284 + t338);
t314 = t176 * t365 + t177 * t363 + t330;
t313 = t257 + t316;
t312 = t294 * (-t284 + t332);
t311 = -t249 + t335;
t201 = rSges(4,1) * t261 - rSges(4,2) * t260 + rSges(4,3) * t280;
t153 = rSges(5,1) * t229 - rSges(5,2) * t228 + rSges(5,3) * t260;
t243 = t281 * rSges(3,1) - t280 * rSges(3,2) - rSges(3,3) * t363;
t310 = t77 / 0.2e1 + t89 / 0.2e1 + t88 / 0.2e1 + t79 / 0.2e1 - t345;
t309 = t80 / 0.2e1 + t78 / 0.2e1 + t91 / 0.2e1 + t90 / 0.2e1 - t344;
t308 = -t214 + t311;
t307 = t177 + t313;
t306 = -t220 + t308;
t195 = Icges(4,5) * t261 - Icges(4,6) * t260 + Icges(4,3) * t280;
t197 = Icges(4,4) * t261 - Icges(4,2) * t260 + Icges(4,6) * t280;
t199 = Icges(4,1) * t261 - Icges(4,4) * t260 + Icges(4,5) * t280;
t107 = -t195 * t364 - t197 * t278 + t199 * t279;
t124 = t233 * t280 - t234 * t260 + t235 * t261;
t305 = t124 / 0.2e1 + t107 / 0.2e1 + t310;
t196 = Icges(4,5) * t263 - Icges(4,6) * t262 + Icges(4,3) * t282;
t198 = Icges(4,4) * t263 - Icges(4,2) * t262 + Icges(4,6) * t282;
t200 = Icges(4,1) * t263 - Icges(4,4) * t262 + Icges(4,5) * t282;
t108 = -t196 * t364 - t198 * t278 + t200 * t279;
t125 = t233 * t282 - t234 * t262 + t235 * t263;
t304 = t125 / 0.2e1 + t108 / 0.2e1 + t309;
t286 = rSges(2,1) * t303 - t300 * rSges(2,2);
t285 = -t300 * rSges(2,1) - rSges(2,2) * t303;
t270 = rSges(3,3) * t295 + (rSges(3,1) * t299 + rSges(3,2) * t302) * t294;
t242 = Icges(3,1) * t283 - Icges(3,4) * t282 + Icges(3,5) * t365;
t241 = Icges(3,1) * t281 - Icges(3,4) * t280 - Icges(3,5) * t363;
t240 = Icges(3,4) * t283 - Icges(3,2) * t282 + Icges(3,6) * t365;
t239 = Icges(3,4) * t281 - Icges(3,2) * t280 - Icges(3,6) * t363;
t238 = Icges(3,5) * t283 - Icges(3,6) * t282 + Icges(3,3) * t365;
t219 = t244 + t346;
t218 = -t243 + t335;
t206 = -t295 * t243 - t270 * t363;
t205 = t244 * t295 - t270 * t365;
t186 = t337 * t295;
t183 = t260 * t213;
t180 = (t243 * t300 + t244 * t303) * t294;
t179 = t267 * t365 - t268 * t282 + t269 * t283;
t178 = -t267 * t363 - t280 * t268 + t281 * t269;
t160 = t278 * t177;
t159 = t316 + t202;
t158 = -t201 + t311;
t156 = t262 * t176;
t152 = rSges(6,1) * t229 - t318;
t139 = -t202 * t364 - t236 * t282;
t138 = t201 * t364 + t236 * t280;
t137 = t238 * t295 + (t240 * t302 + t242 * t299) * t294;
t136 = t237 * t295 + (t239 * t302 + t241 * t299) * t294;
t129 = t130 * t295;
t128 = t201 * t282 - t202 * t280;
t127 = (-t201 - t249) * t295 + t303 * t334;
t126 = t202 * t295 + t300 * t334 + t247;
t121 = t155 + t313 + t369;
t120 = -t153 + t308;
t119 = (t201 * t300 + t202 * t303) * t294 + t347;
t110 = t155 * t278 - t194 * t262;
t109 = -t153 * t278 + t194 * t260;
t106 = t196 * t282 - t198 * t262 + t200 * t263;
t105 = t195 * t282 - t197 * t262 + t199 * t263;
t104 = t196 * t280 - t198 * t260 + t200 * t261;
t103 = t195 * t280 - t197 * t260 + t199 * t261;
t96 = t153 * t262 - t155 * t260;
t95 = t154 + t307 + t369;
t94 = (-rSges(6,1) - pkin(4)) * t229 + t306 + t318;
t93 = t282 * t351 + t354 * t364;
t92 = t153 * t364 + t194 * t280 + t350;
t87 = (-t153 + t348) * t295 + t303 * t320;
t86 = t155 * t295 + t300 * t320 + t349;
t85 = -t118 * t278 + t135 * t262;
t84 = t117 * t278 - t135 * t260;
t83 = t153 * t282 + t280 * t354 + t204;
t82 = t154 * t278 + t262 * t352 + t160;
t81 = t193 * t260 + t183 + (-t152 - t176) * t278;
t76 = (t153 * t300 + t155 * t303) * t294 + t330;
t67 = (pkin(10) - pkin(11)) * t262 + t307 + t379;
t66 = (rSges(7,3) + pkin(11)) * t260 + (-pkin(4) - pkin(5)) * t229 + t306 + t319;
t65 = t282 * t338 + t341 * t364;
t64 = t152 * t364 + t193 * t280 + t331;
t63 = -t117 * t262 + t118 * t260;
t62 = (-t152 + t339) * t295 + t303 * t315;
t61 = t154 * t295 + t300 * t315 + t340;
t60 = t152 * t262 + t260 * t355 + t156;
t55 = t152 * t282 + t280 * t341 + t353;
t52 = (t152 * t300 + t154 * t303) * t294 + t314;
t51 = t129 + (-t107 * t303 + t108 * t300) * t294;
t50 = t107 * t280 + t108 * t282 - t130 * t364;
t49 = t262 * t342 + t278 * t357 + t160;
t48 = t183 + t356 * t260 + (-t176 - t358) * t278;
t47 = t125 * t295 + (-t105 * t303 + t106 * t300) * t294;
t46 = t124 * t295 + (-t103 * t303 + t104 * t300) * t294;
t45 = t105 * t280 + t106 * t282 - t125 * t364;
t44 = t103 * t280 + t104 * t282 - t124 * t364;
t43 = t282 * t332 + t333 * t364;
t42 = t280 * t356 + t358 * t364 + t331;
t41 = (t339 - t358) * t295 + t303 * t312;
t40 = t295 * t357 + t300 * t312 + t340;
t33 = t260 * t343 + t262 * t358 + t156;
t32 = t280 * t333 + t282 * t358 + t353;
t31 = (t300 * t358 + t303 * t357) * t294 + t314;
t2 = [m(7) * (t66 ^ 2 + t67 ^ 2) + m(6) * (t94 ^ 2 + t95 ^ 2) + m(5) * (t120 ^ 2 + t121 ^ 2) + m(4) * (t158 ^ 2 + t159 ^ 2) + m(3) * (t218 ^ 2 + t219 ^ 2) + m(2) * (t285 ^ 2 + t286 ^ 2) + Icges(2,3) + t337 - t378; t58 + t186 + t99 + t100 + t129 + m(4) * (t126 * t159 + t127 * t158) + m(3) * (t205 * t219 + t206 * t218) + m(7) * (t40 * t67 + t41 * t66) + m(6) * (t61 * t95 + t62 * t94) + m(5) * (t120 * t87 + t121 * t86) + ((-t136 / 0.2e1 - t178 / 0.2e1 - t305) * t303 + (t137 / 0.2e1 + t179 / 0.2e1 + t304) * t300) * t294; (t12 + t29 + t30 + t51 + t186) * t295 + m(7) * (t31 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t52 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t76 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t119 ^ 2 + t126 ^ 2 + t127 ^ 2) + m(3) * (t180 ^ 2 + t205 ^ 2 + t206 ^ 2) + ((-t7 - t22 - t21 - t46 + (-t280 * t239 + t281 * t241 - t294 * t367) * t363) * t303 + (t8 + t24 + t23 + t47 + ((-t240 * t282 + t242 * t283 + (t238 * t300 - t367) * t294) * t300 + (t238 * t363 + t239 * t282 + t280 * t240 - t241 * t283 - t281 * t242) * t303) * t294) * t300 + ((-t136 - t178) * t303 + (t137 + t179) * t300) * t295) * t294; t378 * t364 + m(7) * (t42 * t66 + t43 * t67) + m(6) * (t64 * t94 + t65 * t95) + m(5) * (t120 * t92 + t121 * t93) + m(4) * (t138 * t158 + t139 * t159) + t304 * t282 + t305 * t280; (t50 / 0.2e1 + t322) * t295 + (t47 / 0.2e1 + t324) * t282 + (t46 / 0.2e1 + t325) * t280 + m(4) * (t119 * t128 + t126 * t139 + t127 * t138) + m(7) * (t31 * t32 + t40 * t43 + t41 * t42) + m(6) * (t52 * t55 + t61 * t65 + t62 * t64) + m(5) * (t76 * t83 + t86 * t93 + t87 * t92) + ((-t44 / 0.2e1 - t327) * t303 + (-t51 / 0.2e1 - t321) * t302 + (t45 / 0.2e1 + t326) * t300) * t294; (-t11 - t27 - t28 - t50) * t364 + (t6 + t20 + t19 + t45) * t282 + (t5 + t18 + t17 + t44) * t280 + m(7) * (t32 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t55 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t83 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(4) * (t128 ^ 2 + t138 ^ 2 + t139 ^ 2); t368 + t98 + t97 + m(7) * (t48 * t66 + t49 * t67) + m(6) * (t81 * t94 + t82 * t95) + m(5) * (t109 * t120 + t110 * t121) + t309 * t262 + t310 * t260; t323 * t295 + t321 * t278 + t324 * t262 + t325 * t260 + m(7) * (t31 * t33 + t40 * t49 + t41 * t48) + m(6) * (t60 * t52 + t61 * t82 + t62 * t81) + m(5) * (t109 * t87 + t110 * t86 + t76 * t96) + (t300 * t328 - t303 * t329) * t294; -t323 * t364 + t328 * t282 + t329 * t280 + t322 * t278 + t326 * t262 + t327 * t260 + m(7) * (t32 * t33 + t42 * t48 + t43 * t49) + m(6) * (t60 * t55 + t64 * t81 + t65 * t82) + m(5) * (t109 * t92 + t110 * t93 + t83 * t96); (t9 + t26 + t25) * t278 + (t3 + t16 + t15) * t262 + (t1 + t14 + t13) * t260 + m(7) * (t33 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t60 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t96 ^ 2); m(7) * (t228 * t67 + t230 * t66) + m(6) * (t228 * t95 + t230 * t94); m(7) * (t228 * t40 + t230 * t41 + t258 * t31) + m(6) * (t228 * t61 + t230 * t62 + t258 * t52); m(7) * (t228 * t43 + t230 * t42 + t258 * t32) + m(6) * (t228 * t65 + t230 * t64 + t258 * t55); m(7) * (t228 * t49 + t230 * t48 + t258 * t33) + m(6) * (t228 * t82 + t230 * t81 + t258 * t60); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t228 ^ 2 + t230 ^ 2 + t258 ^ 2); m(7) * (t66 * t84 + t67 * t85) - t368 + t344 * t262 + t345 * t260; t8 * t373 + t7 * t374 + m(7) * (t63 * t31 + t40 * t85 + t41 * t84) + t295 * t375 + t12 * t372 + (t300 * t376 + t303 * t380) * t294; t364 * t381 + m(7) * (t63 * t32 + t42 * t84 + t43 * t85) + t280 * t377 + t5 * t374 + t6 * t373 + t11 * t372 + t282 * t376; m(7) * (t63 * t33 + t48 * t84 + t49 * t85) + 0.2e1 * t375 * t278 + 0.2e1 * t376 * t262 + 0.2e1 * t377 * t260; m(7) * (t228 * t85 + t230 * t84 + t258 * t63); t260 * t1 + t278 * t9 + t262 * t3 + m(7) * (t63 ^ 2 + t84 ^ 2 + t85 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
