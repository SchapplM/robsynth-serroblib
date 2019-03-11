% Calculate joint inertia matrix for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:11
% EndTime: 2019-03-09 21:11:24
% DurationCPUTime: 6.26s
% Computational Cost: add. (11511->510), mult. (16238->717), div. (0->0), fcn. (17839->8), ass. (0->252)
t274 = qJ(3) + qJ(4);
t267 = sin(t274);
t268 = cos(t274);
t280 = cos(qJ(1));
t277 = sin(qJ(1));
t279 = cos(qJ(2));
t335 = t277 * t279;
t230 = t267 * t335 + t268 * t280;
t333 = t280 * t267;
t231 = t268 * t335 - t333;
t276 = sin(qJ(2));
t338 = t276 * t277;
t134 = Icges(7,5) * t338 + Icges(7,6) * t230 + Icges(7,3) * t231;
t140 = Icges(6,4) * t338 - Icges(6,2) * t231 + Icges(6,6) * t230;
t150 = Icges(5,1) * t231 - Icges(5,4) * t230 + Icges(5,5) * t338;
t385 = t134 - t140 + t150;
t232 = -t277 * t268 + t279 * t333;
t334 = t279 * t280;
t233 = t267 * t277 + t268 * t334;
t337 = t276 * t280;
t135 = Icges(7,5) * t337 + Icges(7,6) * t232 + Icges(7,3) * t233;
t141 = Icges(6,4) * t337 - Icges(6,2) * t233 + Icges(6,6) * t232;
t151 = Icges(5,1) * t233 - Icges(5,4) * t232 + Icges(5,5) * t337;
t384 = t135 - t141 + t151;
t136 = Icges(6,5) * t338 - Icges(6,6) * t231 + Icges(6,3) * t230;
t138 = Icges(7,4) * t338 + Icges(7,2) * t230 + Icges(7,6) * t231;
t148 = Icges(5,4) * t231 - Icges(5,2) * t230 + Icges(5,6) * t338;
t383 = t136 + t138 - t148;
t137 = Icges(6,5) * t337 - Icges(6,6) * t233 + Icges(6,3) * t232;
t139 = Icges(7,4) * t337 + Icges(7,2) * t232 + Icges(7,6) * t233;
t149 = Icges(5,4) * t233 - Icges(5,2) * t232 + Icges(5,6) * t337;
t382 = t137 + t139 - t149;
t142 = Icges(7,1) * t338 + Icges(7,4) * t230 + Icges(7,5) * t231;
t144 = Icges(6,1) * t338 - Icges(6,4) * t231 + Icges(6,5) * t230;
t146 = Icges(5,5) * t231 - Icges(5,6) * t230 + Icges(5,3) * t338;
t381 = t142 + t144 + t146;
t143 = Icges(7,1) * t337 + Icges(7,4) * t232 + Icges(7,5) * t233;
t145 = Icges(6,1) * t337 - Icges(6,4) * t233 + Icges(6,5) * t232;
t147 = Icges(5,5) * t233 - Icges(5,6) * t232 + Icges(5,3) * t337;
t380 = t143 + t145 + t147;
t368 = rSges(7,1) + pkin(5);
t367 = rSges(7,3) + qJ(6);
t200 = -Icges(5,6) * t279 + (Icges(5,4) * t268 - Icges(5,2) * t267) * t276;
t205 = -Icges(6,4) * t279 + (-Icges(6,2) * t268 + Icges(6,6) * t267) * t276;
t341 = t268 * t276;
t342 = t267 * t276;
t201 = -Icges(5,5) * t279 + (Icges(5,1) * t268 - Icges(5,4) * t267) * t276;
t202 = -Icges(7,5) * t279 + (Icges(7,6) * t267 + Icges(7,3) * t268) * t276;
t203 = -Icges(6,5) * t279 + (-Icges(6,6) * t268 + Icges(6,3) * t267) * t276;
t204 = -Icges(7,4) * t279 + (Icges(7,2) * t267 + Icges(7,6) * t268) * t276;
t370 = (t203 + t204) * t342 + (t201 + t202) * t341;
t199 = -Icges(5,3) * t279 + (Icges(5,5) * t268 - Icges(5,6) * t267) * t276;
t206 = -Icges(7,1) * t279 + (Icges(7,4) * t267 + Icges(7,5) * t268) * t276;
t207 = -Icges(6,1) * t279 + (-Icges(6,4) * t268 + Icges(6,5) * t267) * t276;
t372 = -t199 - t206 - t207;
t358 = -t200 * t342 - t205 * t341 + t279 * t372 + t370;
t379 = t358 * t279;
t378 = t383 * t230 + t385 * t231 + t381 * t338;
t377 = t382 * t230 + t384 * t231 + t380 * t338;
t376 = t383 * t232 + t385 * t233 + t381 * t337;
t375 = t382 * t232 + t384 * t233 + t380 * t337;
t100 = t202 * t231 + t204 * t230 + t206 * t338;
t101 = t203 * t230 - t205 * t231 + t207 * t338;
t104 = t199 * t338 - t200 * t230 + t201 * t231;
t374 = -t100 - t101 - t104;
t102 = t202 * t233 + t204 * t232 + t206 * t337;
t103 = t203 * t232 - t205 * t233 + t207 * t337;
t105 = t199 * t337 - t200 * t232 + t201 * t233;
t373 = -t102 - t103 - t105;
t371 = Icges(3,5) * t276;
t369 = t371 / 0.2e1;
t348 = rSges(7,2) * t230;
t366 = t367 * t231 + t338 * t368 + t348;
t365 = (rSges(7,2) * t267 + rSges(7,3) * t268) * t276 + qJ(6) * t341 - t368 * t279;
t364 = t374 * t279 + (t277 * t378 + t280 * t377) * t276;
t363 = t373 * t279 + (t277 * t376 + t280 * t375) * t276;
t362 = t277 * t377 - t280 * t378;
t361 = t277 * t375 - t280 * t376;
t73 = -t279 * t146 + (-t148 * t267 + t150 * t268) * t276;
t75 = -t279 * t142 + (t134 * t268 + t138 * t267) * t276;
t77 = -t279 * t144 + (t136 * t267 - t140 * t268) * t276;
t360 = -t73 - t75 - t77;
t74 = -t279 * t147 + (-t149 * t267 + t151 * t268) * t276;
t76 = -t279 * t143 + (t135 * t268 + t139 * t267) * t276;
t78 = -t279 * t145 + (t137 * t267 - t141 * t268) * t276;
t359 = t74 + t76 + t78;
t357 = t277 ^ 2;
t356 = t280 ^ 2;
t355 = t277 / 0.2e1;
t354 = -t279 / 0.2e1;
t353 = -t280 / 0.2e1;
t352 = t280 / 0.2e1;
t351 = pkin(2) * t279;
t350 = pkin(8) * t276;
t278 = cos(qJ(3));
t266 = pkin(3) * t278 + pkin(2);
t349 = -pkin(2) + t266;
t347 = rSges(6,3) * t230;
t346 = t280 * rSges(3,3);
t344 = Icges(3,4) * t279;
t275 = sin(qJ(3));
t223 = -Icges(4,6) * t279 + (Icges(4,4) * t278 - Icges(4,2) * t275) * t276;
t343 = t223 * t275;
t340 = t275 * t277;
t339 = t275 * t280;
t281 = -pkin(9) - pkin(8);
t336 = t276 * t281;
t153 = rSges(6,1) * t338 - rSges(6,2) * t231 + t347;
t211 = t230 * qJ(5);
t168 = pkin(4) * t231 + t211;
t158 = t168 * t337;
t332 = t153 * t337 + t158;
t331 = t232 * rSges(7,2) + t367 * t233 + t368 * t337;
t155 = rSges(6,1) * t337 - t233 * rSges(6,2) + t232 * rSges(6,3);
t169 = t233 * pkin(4) + t232 * qJ(5);
t330 = -t155 - t169;
t157 = t233 * rSges(5,1) - t232 * rSges(5,2) + rSges(5,3) * t337;
t288 = pkin(3) * t340 + t266 * t334 - t280 * t336;
t320 = pkin(2) * t334 + pkin(8) * t337;
t179 = t288 - t320;
t329 = -t157 - t179;
t234 = (pkin(4) * t268 + qJ(5) * t267) * t276;
t328 = t279 * t168 + t234 * t338;
t321 = -pkin(3) * t339 - t277 * t336;
t178 = (t279 * t349 - t350) * t277 + t321;
t198 = (pkin(8) + t281) * t279 + t349 * t276;
t327 = t279 * t178 + t198 * t338;
t296 = -rSges(5,1) * t231 + rSges(5,2) * t230;
t156 = rSges(5,3) * t338 - t296;
t208 = -t279 * rSges(5,3) + (rSges(5,1) * t268 - rSges(5,2) * t267) * t276;
t119 = t279 * t156 + t208 * t338;
t325 = -t198 - t208;
t210 = -t279 * rSges(6,1) + (-rSges(6,2) * t268 + rSges(6,3) * t267) * t276;
t324 = -t210 - t234;
t229 = -t279 * rSges(4,3) + (rSges(4,1) * t278 - rSges(4,2) * t275) * t276;
t252 = t276 * pkin(2) - t279 * pkin(8);
t323 = -t229 - t252;
t322 = t357 * (t350 + t351) + t280 * t320;
t319 = t280 * pkin(1) + t277 * pkin(7);
t318 = t379 + (t360 * t277 - t359 * t280) * t276;
t220 = -Icges(4,3) * t279 + (Icges(4,5) * t278 - Icges(4,6) * t275) * t276;
t226 = -Icges(4,5) * t279 + (Icges(4,1) * t278 - Icges(4,4) * t275) * t276;
t241 = -t275 * t334 + t277 * t278;
t242 = t278 * t334 + t340;
t110 = t220 * t337 + t223 * t241 + t226 * t242;
t171 = Icges(4,5) * t242 + Icges(4,6) * t241 + Icges(4,3) * t337;
t173 = Icges(4,4) * t242 + Icges(4,2) * t241 + Icges(4,6) * t337;
t175 = Icges(4,1) * t242 + Icges(4,4) * t241 + Icges(4,5) * t337;
t90 = -t279 * t171 + (-t173 * t275 + t175 * t278) * t276;
t317 = t110 / 0.2e1 + t90 / 0.2e1;
t239 = -t275 * t335 - t278 * t280;
t240 = t278 * t335 - t339;
t109 = t220 * t338 + t223 * t239 + t226 * t240;
t170 = Icges(4,5) * t240 + Icges(4,6) * t239 + Icges(4,3) * t338;
t172 = Icges(4,4) * t240 + Icges(4,2) * t239 + Icges(4,6) * t338;
t174 = Icges(4,1) * t240 + Icges(4,4) * t239 + Icges(4,5) * t338;
t89 = -t279 * t170 + (-t172 * t275 + t174 * t278) * t276;
t316 = t89 / 0.2e1 + t109 / 0.2e1;
t314 = t366 * t337 + t158;
t313 = -t169 - t331;
t312 = -t179 + t330;
t311 = -t252 + t325;
t310 = -t198 + t324;
t309 = -t234 - t365;
t177 = t242 * rSges(4,1) + t241 * rSges(4,2) + rSges(4,3) * t337;
t271 = t280 * pkin(7);
t308 = t271 - t321;
t307 = t338 / 0.2e1;
t306 = t337 / 0.2e1;
t305 = -t266 * t279 - pkin(1);
t304 = -t179 + t313;
t303 = t277 * t178 + t280 * t179 + t322;
t87 = t279 * t153 + t210 * t338 + t328;
t302 = -t198 + t309;
t301 = -t252 + t310;
t300 = -t211 + t308;
t299 = t363 * t337 + t364 * t338;
t298 = rSges(3,1) * t279 - rSges(3,2) * t276;
t297 = -rSges(4,1) * t240 - rSges(4,2) * t239;
t295 = -t252 + t302;
t293 = -Icges(3,2) * t276 + t344;
t292 = Icges(3,5) * t279 - Icges(3,6) * t276;
t289 = rSges(3,1) * t334 - rSges(3,2) * t337 + t277 * rSges(3,3);
t287 = t277 * t168 + t280 * t169 + t303;
t53 = t366 * t279 + t365 * t338 + t328;
t286 = t288 + t319;
t285 = t279 * t318 + t299;
t284 = (-t360 - t374) * t307 + (t359 - t373) * t306;
t283 = t169 + t286;
t282 = t363 * t355 + (t359 * t277 + t360 * t280) * t354 + t364 * t353 + t362 * t307 + t361 * t306;
t273 = t276 ^ 2;
t251 = rSges(2,1) * t280 - rSges(2,2) * t277;
t250 = -rSges(2,1) * t277 - rSges(2,2) * t280;
t249 = rSges(3,1) * t276 + rSges(3,2) * t279;
t246 = Icges(3,6) * t279 + t371;
t222 = Icges(3,3) * t277 + t280 * t292;
t221 = -Icges(3,3) * t280 + t277 * t292;
t196 = t276 * t278 * t226;
t195 = t289 + t319;
t194 = t346 + t271 + (-pkin(1) - t298) * t277;
t181 = t323 * t280;
t180 = t323 * t277;
t176 = rSges(4,3) * t338 - t297;
t164 = t280 * t289 + (t277 * t298 - t346) * t277;
t160 = t178 * t337;
t130 = t156 * t337;
t127 = t177 + t319 + t320;
t126 = t271 + (-t351 - pkin(1) + (-rSges(4,3) - pkin(8)) * t276) * t277 + t297;
t125 = t311 * t280;
t124 = t311 * t277;
t123 = -t177 * t279 - t229 * t337;
t122 = t176 * t279 + t229 * t338;
t121 = -t279 * t220 - t276 * t343 + t196;
t120 = -t279 * t157 - t208 * t337;
t118 = t286 + t157;
t117 = (-rSges(5,3) * t276 + t305) * t277 + t296 + t308;
t116 = t301 * t280;
t115 = t301 * t277;
t111 = (t176 * t280 - t177 * t277) * t276;
t108 = t295 * t280;
t107 = t295 * t277;
t106 = -t157 * t338 + t130;
t93 = t176 * t277 + t177 * t280 + t322;
t92 = t283 + t155;
t91 = -t347 + (rSges(6,2) - pkin(4)) * t231 + (-rSges(6,1) * t276 + t305) * t277 + t300;
t88 = t279 * t330 + t324 * t337;
t86 = t283 + t331;
t85 = -t348 + (-pkin(4) - t367) * t231 + (-t368 * t276 + t305) * t277 + t300;
t84 = t279 * t329 + t325 * t337;
t83 = t119 + t327;
t82 = t171 * t337 + t173 * t241 + t175 * t242;
t81 = t170 * t337 + t172 * t241 + t174 * t242;
t80 = t171 * t338 + t173 * t239 + t175 * t240;
t79 = t170 * t338 + t172 * t239 + t174 * t240;
t54 = t279 * t313 + t309 * t337;
t52 = t329 * t338 + t130 + t160;
t51 = t330 * t338 + t332;
t50 = t279 * t312 + t310 * t337;
t49 = t87 + t327;
t48 = t156 * t277 + t157 * t280 + t303;
t47 = t313 * t338 + t314;
t46 = t279 * t304 + t302 * t337;
t45 = t53 + t327;
t44 = t312 * t338 + t160 + t332;
t43 = t153 * t277 + t155 * t280 + t287;
t42 = t277 * t82 - t280 * t81;
t41 = t277 * t80 - t280 * t79;
t37 = t304 * t338 + t160 + t314;
t36 = t366 * t277 + t331 * t280 + t287;
t23 = -t110 * t279 + (t277 * t81 + t280 * t82) * t276;
t22 = -t109 * t279 + (t277 * t79 + t280 * t80) * t276;
t1 = [Icges(2,3) + t196 + (Icges(3,1) * t276 - t200 * t267 - t205 * t268 - t343 + t344) * t276 + (Icges(3,4) * t276 + Icges(3,2) * t279 - t220 + t372) * t279 + m(7) * (t85 ^ 2 + t86 ^ 2) + m(6) * (t91 ^ 2 + t92 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t194 ^ 2 + t195 ^ 2) + m(2) * (t250 ^ 2 + t251 ^ 2) + t370; (-t77 / 0.2e1 - t75 / 0.2e1 - t73 / 0.2e1 - t100 / 0.2e1 - t101 / 0.2e1 - t104 / 0.2e1 + (-Icges(3,6) * t280 + t277 * t293) * t354 + t280 * t369 + t246 * t352 - t316) * t280 + (t78 / 0.2e1 + t76 / 0.2e1 + t74 / 0.2e1 + t102 / 0.2e1 + t103 / 0.2e1 + t105 / 0.2e1 + t279 * (Icges(3,6) * t277 + t280 * t293) / 0.2e1 + t277 * t369 + t246 * t355 + t317) * t277 + m(6) * (t115 * t92 + t116 * t91) + m(7) * (t107 * t86 + t108 * t85) + m(5) * (t117 * t125 + t118 * t124) + m(4) * (t126 * t181 + t127 * t180) + m(3) * (-t194 * t280 - t195 * t277) * t249; m(7) * (t107 ^ 2 + t108 ^ 2 + t36 ^ 2) + m(6) * (t115 ^ 2 + t116 ^ 2 + t43 ^ 2) + m(5) * (t124 ^ 2 + t125 ^ 2 + t48 ^ 2) + m(4) * (t180 ^ 2 + t181 ^ 2 + t93 ^ 2) + m(3) * (t164 ^ 2 + (t356 + t357) * t249 ^ 2) + (-t356 * t221 - t362 - t41) * t280 + (t357 * t222 + t42 + (-t277 * t221 + t280 * t222) * t280 + t361) * t277; (t277 * t316 + t280 * t317) * t276 + m(7) * (t45 * t85 + t46 * t86) + m(6) * (t49 * t91 + t50 * t92) + m(5) * (t117 * t83 + t118 * t84) + m(4) * (t122 * t126 + t123 * t127) + (-t121 - t358) * t279 + t284; m(7) * (t107 * t46 + t108 * t45 + t36 * t37) + m(6) * (t115 * t50 + t116 * t49 + t43 * t44) + m(5) * (t124 * t84 + t125 * t83 + t48 * t52) + m(4) * (t111 * t93 + t122 * t181 + t123 * t180) + (t352 * t42 + t355 * t41) * t276 + t282 + t23 * t355 + t22 * t353 + (t90 * t277 - t89 * t280) * t354; (t277 * t22 + t280 * t23) * t276 + (t121 * t279 + (-t277 * t89 - t280 * t90) * t276 + t318) * t279 + m(5) * (t52 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t111 ^ 2 + t122 ^ 2 + t123 ^ 2) + m(7) * (t37 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t44 ^ 2 + t49 ^ 2 + t50 ^ 2) + t299; m(7) * (t53 * t85 + t54 * t86) + m(6) * (t87 * t91 + t88 * t92) + m(5) * (t117 * t119 + t118 * t120) - t379 + t284; t282 + m(7) * (t107 * t54 + t108 * t53 + t36 * t47) + m(6) * (t115 * t88 + t116 * t87 + t43 * t51) + m(5) * (t106 * t48 + t119 * t125 + t120 * t124); m(7) * (t37 * t47 + t45 * t53 + t46 * t54) + m(6) * (t44 * t51 + t49 * t87 + t50 * t88) + m(5) * (t106 * t52 + t119 * t83 + t120 * t84) + t285; m(7) * (t47 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t51 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t106 ^ 2 + t119 ^ 2 + t120 ^ 2) + t285; m(7) * (t230 * t86 + t232 * t85) + m(6) * (t230 * t92 + t232 * t91); m(7) * (t107 * t230 + t108 * t232 + t342 * t36) + m(6) * (t115 * t230 + t116 * t232 + t342 * t43); m(7) * (t230 * t46 + t232 * t45 + t342 * t37) + m(6) * (t230 * t50 + t232 * t49 + t342 * t44); m(7) * (t230 * t54 + t232 * t53 + t342 * t47) + m(6) * (t230 * t88 + t232 * t87 + t342 * t51); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t267 ^ 2 * t273 + t230 ^ 2 + t232 ^ 2); m(7) * (t231 * t86 + t233 * t85); m(7) * (t107 * t231 + t108 * t233 + t341 * t36); m(7) * (t231 * t46 + t233 * t45 + t341 * t37); m(7) * (t231 * t54 + t233 * t53 + t341 * t47); m(7) * (t267 * t268 * t273 + t230 * t231 + t232 * t233); m(7) * (t268 ^ 2 * t273 + t231 ^ 2 + t233 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
