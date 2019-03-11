% Calculate joint inertia matrix for
% S6RRRRPP5
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:18
% EndTime: 2019-03-09 21:05:32
% DurationCPUTime: 6.24s
% Computational Cost: add. (11470->510), mult. (16210->718), div. (0->0), fcn. (17799->8), ass. (0->255)
t276 = qJ(3) + qJ(4);
t267 = sin(t276);
t268 = cos(t276);
t282 = cos(qJ(1));
t279 = sin(qJ(1));
t281 = cos(qJ(2));
t339 = t279 * t281;
t233 = t267 * t339 + t268 * t282;
t234 = -t267 * t282 + t268 * t339;
t278 = sin(qJ(2));
t342 = t278 * t279;
t134 = Icges(7,5) * t234 + Icges(7,6) * t233 - Icges(7,3) * t342;
t138 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t342;
t142 = Icges(6,4) * t234 + Icges(6,2) * t342 + Icges(6,6) * t233;
t388 = -t134 + t138 + t142;
t338 = t281 * t282;
t235 = t267 * t338 - t279 * t268;
t236 = t267 * t279 + t268 * t338;
t341 = t278 * t282;
t135 = Icges(7,5) * t236 + Icges(7,6) * t235 - Icges(7,3) * t341;
t139 = Icges(5,5) * t236 - Icges(5,6) * t235 + Icges(5,3) * t341;
t143 = Icges(6,4) * t236 + Icges(6,2) * t341 + Icges(6,6) * t235;
t387 = -t135 + t139 + t143;
t136 = Icges(6,5) * t234 + Icges(6,6) * t342 + Icges(6,3) * t233;
t140 = Icges(7,4) * t234 + Icges(7,2) * t233 - Icges(7,6) * t342;
t144 = Icges(5,4) * t234 - Icges(5,2) * t233 + Icges(5,6) * t342;
t386 = t136 + t140 - t144;
t137 = Icges(6,5) * t236 + Icges(6,6) * t341 + Icges(6,3) * t235;
t141 = Icges(7,4) * t236 + Icges(7,2) * t235 - Icges(7,6) * t341;
t145 = Icges(5,4) * t236 - Icges(5,2) * t235 + Icges(5,6) * t341;
t385 = t137 + t141 - t145;
t146 = Icges(7,1) * t234 + Icges(7,4) * t233 - Icges(7,5) * t342;
t148 = Icges(6,1) * t234 + Icges(6,4) * t342 + Icges(6,5) * t233;
t150 = Icges(5,1) * t234 - Icges(5,4) * t233 + Icges(5,5) * t342;
t384 = t146 + t148 + t150;
t147 = Icges(7,1) * t236 + Icges(7,4) * t235 - Icges(7,5) * t341;
t149 = Icges(6,1) * t236 + Icges(6,4) * t341 + Icges(6,5) * t235;
t151 = Icges(5,1) * t236 - Icges(5,4) * t235 + Icges(5,5) * t341;
t383 = t147 + t149 + t151;
t375 = rSges(7,1) + pkin(5);
t350 = rSges(7,3) + qJ(6);
t203 = -Icges(5,3) * t281 + (Icges(5,5) * t268 - Icges(5,6) * t267) * t278;
t205 = -Icges(6,2) * t281 + (Icges(6,4) * t268 + Icges(6,6) * t267) * t278;
t382 = t203 + t205;
t381 = t386 * t233 + t384 * t234 + t342 * t388;
t380 = t385 * t233 + t234 * t383 + t342 * t387;
t379 = t386 * t235 + t384 * t236 + t341 * t388;
t378 = t235 * t385 + t236 * t383 + t341 * t387;
t201 = Icges(7,3) * t281 + (Icges(7,5) * t268 + Icges(7,6) * t267) * t278;
t204 = Icges(7,6) * t281 + (Icges(7,4) * t268 + Icges(7,2) * t267) * t278;
t207 = Icges(7,5) * t281 + (Icges(7,1) * t268 + Icges(7,4) * t267) * t278;
t100 = -t201 * t342 + t204 * t233 + t207 * t234;
t202 = -Icges(6,6) * t281 + (Icges(6,5) * t268 + Icges(6,3) * t267) * t278;
t208 = -Icges(6,4) * t281 + (Icges(6,1) * t268 + Icges(6,5) * t267) * t278;
t101 = t202 * t233 + t205 * t342 + t208 * t234;
t206 = -Icges(5,6) * t281 + (Icges(5,4) * t268 - Icges(5,2) * t267) * t278;
t209 = -Icges(5,5) * t281 + (Icges(5,1) * t268 - Icges(5,4) * t267) * t278;
t102 = t203 * t342 - t206 * t233 + t209 * t234;
t377 = -t100 - t101 - t102;
t103 = -t201 * t341 + t204 * t235 + t207 * t236;
t104 = t202 * t235 + t205 * t341 + t208 * t236;
t105 = t203 * t341 - t206 * t235 + t209 * t236;
t376 = -t103 - t104 - t105;
t374 = Icges(3,5) * t278;
t345 = t268 * t278;
t346 = t267 * t278;
t373 = t281 * t201 + (t202 + t204) * t346 + (t207 + t208 + t209) * t345;
t372 = t374 / 0.2e1;
t316 = t206 * t346 + t281 * t382 - t373;
t371 = t316 * t281;
t353 = rSges(7,2) * t233;
t370 = t234 * t375 - t350 * t342 + t353;
t369 = (rSges(7,1) * t268 + rSges(7,2) * t267) * t278 + pkin(5) * t345 + t350 * t281;
t368 = t377 * t281 + (t279 * t381 + t282 * t380) * t278;
t367 = t376 * t281 + (t279 * t379 + t282 * t378) * t278;
t366 = t279 * t380 - t282 * t381;
t365 = t279 * t378 - t282 * t379;
t73 = t281 * t134 + (t140 * t267 + t146 * t268) * t278;
t75 = -t281 * t142 + (t136 * t267 + t148 * t268) * t278;
t77 = -t281 * t138 + (-t144 * t267 + t150 * t268) * t278;
t364 = -t73 - t75 - t77;
t74 = t281 * t135 + (t141 * t267 + t147 * t268) * t278;
t76 = -t281 * t143 + (t137 * t267 + t149 * t268) * t278;
t78 = -t281 * t139 + (-t145 * t267 + t151 * t268) * t278;
t363 = t74 + t76 + t78;
t362 = t235 * rSges(7,2) + t236 * t375;
t274 = t279 ^ 2;
t275 = t282 ^ 2;
t361 = t279 / 0.2e1;
t360 = -t281 / 0.2e1;
t359 = -t282 / 0.2e1;
t358 = t282 / 0.2e1;
t357 = m(7) * t278;
t356 = pkin(2) * t281;
t355 = pkin(8) * t278;
t280 = cos(qJ(3));
t266 = pkin(3) * t280 + pkin(2);
t354 = -pkin(2) + t266;
t352 = t233 * rSges(6,3);
t351 = t282 * rSges(3,3);
t348 = Icges(3,4) * t281;
t277 = sin(qJ(3));
t226 = -Icges(4,6) * t281 + (Icges(4,4) * t280 - Icges(4,2) * t277) * t278;
t347 = t226 * t277;
t344 = t277 * t279;
t343 = t277 * t282;
t283 = -pkin(9) - pkin(8);
t340 = t278 * t283;
t153 = t234 * rSges(6,1) + rSges(6,2) * t342 + t352;
t213 = t233 * qJ(5);
t168 = t234 * pkin(4) + t213;
t158 = t168 * t341;
t337 = t153 * t341 + t158;
t336 = -t341 * t350 + t362;
t156 = t236 * rSges(6,1) + rSges(6,2) * t341 + t235 * rSges(6,3);
t169 = t236 * pkin(4) + t235 * qJ(5);
t335 = -t156 - t169;
t157 = t236 * rSges(5,1) - t235 * rSges(5,2) + rSges(5,3) * t341;
t319 = t282 * t340;
t325 = pkin(3) * t344 + t266 * t338;
t289 = -t319 + t325;
t323 = pkin(2) * t338 + pkin(8) * t341;
t179 = t289 - t323;
t334 = -t157 - t179;
t237 = (pkin(4) * t268 + qJ(5) * t267) * t278;
t333 = t281 * t168 + t237 * t342;
t324 = -pkin(3) * t343 - t279 * t340;
t178 = (t281 * t354 - t355) * t279 + t324;
t200 = (pkin(8) + t283) * t281 + t354 * t278;
t332 = t281 * t178 + t200 * t342;
t297 = -t234 * rSges(5,1) + t233 * rSges(5,2);
t154 = rSges(5,3) * t342 - t297;
t212 = -t281 * rSges(5,3) + (rSges(5,1) * t268 - rSges(5,2) * t267) * t278;
t119 = t281 * t154 + t212 * t342;
t330 = -t200 - t212;
t211 = -t281 * rSges(6,2) + (rSges(6,1) * t268 + rSges(6,3) * t267) * t278;
t329 = -t211 - t237;
t232 = -t281 * rSges(4,3) + (rSges(4,1) * t280 - rSges(4,2) * t277) * t278;
t255 = t278 * pkin(2) - t281 * pkin(8);
t327 = -t232 - t255;
t326 = t274 * (t355 + t356) + t282 * t323;
t322 = t282 * pkin(1) + t279 * pkin(7);
t321 = t274 + t275;
t320 = -t371 + (t279 * t364 - t282 * t363) * t278;
t223 = -Icges(4,3) * t281 + (Icges(4,5) * t280 - Icges(4,6) * t277) * t278;
t229 = -Icges(4,5) * t281 + (Icges(4,1) * t280 - Icges(4,4) * t277) * t278;
t242 = -t277 * t339 - t280 * t282;
t243 = t280 * t339 - t343;
t109 = t223 * t342 + t226 * t242 + t229 * t243;
t170 = Icges(4,5) * t243 + Icges(4,6) * t242 + Icges(4,3) * t342;
t172 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t342;
t174 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t342;
t89 = -t281 * t170 + (-t172 * t277 + t174 * t280) * t278;
t318 = t109 / 0.2e1 + t89 / 0.2e1;
t244 = -t277 * t338 + t279 * t280;
t245 = t280 * t338 + t344;
t110 = t223 * t341 + t226 * t244 + t229 * t245;
t171 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t341;
t173 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t341;
t175 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t341;
t90 = -t281 * t171 + (-t173 * t277 + t175 * t280) * t278;
t317 = t110 / 0.2e1 + t90 / 0.2e1;
t315 = t341 * t370 + t158;
t314 = -t169 - t336;
t313 = -t179 + t335;
t312 = -t200 + t329;
t311 = -t255 + t330;
t310 = -t237 - t369;
t177 = t245 * rSges(4,1) + t244 * rSges(4,2) + rSges(4,3) * t341;
t271 = t282 * pkin(7);
t309 = t271 - t324;
t308 = t342 / 0.2e1;
t307 = t341 / 0.2e1;
t306 = -t266 * t281 - pkin(1);
t305 = -t179 + t314;
t304 = t279 * t178 + t282 * t179 + t326;
t87 = t281 * t153 + t211 * t342 + t333;
t303 = -t200 + t310;
t302 = -t255 + t312;
t301 = -t213 + t309;
t300 = t341 * t367 + t342 * t368;
t299 = rSges(3,1) * t281 - rSges(3,2) * t278;
t298 = -t243 * rSges(4,1) - t242 * rSges(4,2);
t296 = -t255 + t303;
t294 = -Icges(3,2) * t278 + t348;
t293 = Icges(3,5) * t281 - Icges(3,6) * t278;
t290 = rSges(3,1) * t338 - rSges(3,2) * t341 + t279 * rSges(3,3);
t288 = t279 * t168 + t282 * t169 + t304;
t53 = t281 * t370 + t342 * t369 + t333;
t287 = t169 + t322 + t325;
t286 = t281 * t320 + t300;
t285 = (-t364 - t377) * t308 + (t363 - t376) * t307;
t284 = t367 * t361 + (t279 * t363 + t282 * t364) * t360 + t368 * t359 + t366 * t308 + t365 * t307;
t273 = t278 ^ 2;
t254 = rSges(2,1) * t282 - rSges(2,2) * t279;
t253 = -rSges(2,1) * t279 - rSges(2,2) * t282;
t252 = rSges(3,1) * t278 + rSges(3,2) * t281;
t249 = Icges(3,6) * t281 + t374;
t225 = Icges(3,3) * t279 + t282 * t293;
t224 = -Icges(3,3) * t282 + t279 * t293;
t198 = t278 * t280 * t229;
t196 = t290 + t322;
t195 = t351 + t271 + (-pkin(1) - t299) * t279;
t181 = t327 * t282;
t180 = t327 * t279;
t176 = rSges(4,3) * t342 - t298;
t164 = t282 * t290 + (t279 * t299 - t351) * t279;
t160 = t178 * t341;
t130 = t154 * t341;
t127 = t177 + t322 + t323;
t126 = t271 + (-t356 - pkin(1) + (-rSges(4,3) - pkin(8)) * t278) * t279 + t298;
t125 = t311 * t282;
t124 = t311 * t279;
t123 = -t177 * t281 - t232 * t341;
t122 = t176 * t281 + t232 * t342;
t121 = -t281 * t223 - t278 * t347 + t198;
t120 = -t281 * t157 - t212 * t341;
t118 = t289 + t157 + t322;
t117 = (-rSges(5,3) * t278 + t306) * t279 + t297 + t309;
t116 = t302 * t282;
t115 = t302 * t279;
t111 = (t176 * t282 - t177 * t279) * t278;
t108 = t296 * t282;
t107 = t296 * t279;
t106 = -t157 * t342 + t130;
t93 = t176 * t279 + t177 * t282 + t326;
t92 = t156 + t287 - t319;
t91 = -t352 + (-rSges(6,1) - pkin(4)) * t234 + (-rSges(6,2) * t278 + t306) * t279 + t301;
t88 = t281 * t335 + t329 * t341;
t86 = (-t283 - t350) * t341 + t287 + t362;
t85 = -t353 + (-pkin(4) - t375) * t234 + (t278 * t350 + t306) * t279 + t301;
t84 = t281 * t334 + t330 * t341;
t83 = t119 + t332;
t82 = t171 * t341 + t173 * t244 + t175 * t245;
t81 = t170 * t341 + t172 * t244 + t174 * t245;
t80 = t171 * t342 + t173 * t242 + t175 * t243;
t79 = t170 * t342 + t172 * t242 + t174 * t243;
t54 = t281 * t314 + t310 * t341;
t52 = t334 * t342 + t130 + t160;
t51 = t335 * t342 + t337;
t50 = t281 * t313 + t312 * t341;
t49 = t87 + t332;
t48 = t154 * t279 + t157 * t282 + t304;
t47 = t314 * t342 + t315;
t46 = t281 * t305 + t303 * t341;
t45 = t53 + t332;
t44 = t313 * t342 + t160 + t337;
t43 = t153 * t279 + t156 * t282 + t288;
t42 = t279 * t82 - t282 * t81;
t41 = t279 * t80 - t282 * t79;
t37 = t305 * t342 + t160 + t315;
t36 = t279 * t370 + t336 * t282 + t288;
t23 = -t110 * t281 + (t279 * t81 + t282 * t82) * t278;
t22 = -t109 * t281 + (t279 * t79 + t282 * t80) * t278;
t1 = [Icges(2,3) + t198 + (Icges(3,1) * t278 - t206 * t267 - t347 + t348) * t278 + (Icges(3,4) * t278 + Icges(3,2) * t281 - t223 - t382) * t281 + m(6) * (t91 ^ 2 + t92 ^ 2) + m(7) * (t85 ^ 2 + t86 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t195 ^ 2 + t196 ^ 2) + m(2) * (t253 ^ 2 + t254 ^ 2) + t373; (-t77 / 0.2e1 - t75 / 0.2e1 - t73 / 0.2e1 - t100 / 0.2e1 - t101 / 0.2e1 - t102 / 0.2e1 + (-Icges(3,6) * t282 + t279 * t294) * t360 + t282 * t372 + t249 * t358 - t318) * t282 + (t78 / 0.2e1 + t76 / 0.2e1 + t74 / 0.2e1 + t103 / 0.2e1 + t104 / 0.2e1 + t105 / 0.2e1 + t281 * (Icges(3,6) * t279 + t282 * t294) / 0.2e1 + t279 * t372 + t249 * t361 + t317) * t279 + m(5) * (t117 * t125 + t118 * t124) + m(6) * (t115 * t92 + t116 * t91) + m(7) * (t107 * t86 + t108 * t85) + m(4) * (t126 * t181 + t127 * t180) + m(3) * (-t195 * t282 - t196 * t279) * t252; m(6) * (t115 ^ 2 + t116 ^ 2 + t43 ^ 2) + m(7) * (t107 ^ 2 + t108 ^ 2 + t36 ^ 2) + m(5) * (t124 ^ 2 + t125 ^ 2 + t48 ^ 2) + m(4) * (t180 ^ 2 + t181 ^ 2 + t93 ^ 2) + m(3) * (t252 ^ 2 * t321 + t164 ^ 2) + (-t275 * t224 - t366 - t41) * t282 + (t274 * t225 + t42 + (-t279 * t224 + t282 * t225) * t282 + t365) * t279; (t279 * t318 + t282 * t317) * t278 + m(6) * (t49 * t91 + t50 * t92) + m(7) * (t45 * t85 + t46 * t86) + m(5) * (t117 * t83 + t118 * t84) + m(4) * (t122 * t126 + t123 * t127) + (-t121 + t316) * t281 + t285; (t358 * t42 + t361 * t41) * t278 + t284 + m(6) * (t115 * t50 + t116 * t49 + t43 * t44) + m(7) * (t107 * t46 + t108 * t45 + t36 * t37) + m(5) * (t124 * t84 + t125 * t83 + t48 * t52) + m(4) * (t111 * t93 + t122 * t181 + t123 * t180) + t23 * t361 + t22 * t359 + (t90 * t279 - t89 * t282) * t360; (t279 * t22 + t282 * t23) * t278 + (t121 * t281 + (-t279 * t89 - t282 * t90) * t278 + t320) * t281 + m(5) * (t52 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t111 ^ 2 + t122 ^ 2 + t123 ^ 2) + m(6) * (t44 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(7) * (t37 ^ 2 + t45 ^ 2 + t46 ^ 2) + t300; t371 + m(6) * (t87 * t91 + t88 * t92) + m(7) * (t53 * t85 + t54 * t86) + m(5) * (t117 * t119 + t118 * t120) + t285; t284 + m(6) * (t115 * t88 + t116 * t87 + t43 * t51) + m(7) * (t107 * t54 + t108 * t53 + t36 * t47) + m(5) * (t106 * t48 + t119 * t125 + t120 * t124); m(6) * (t44 * t51 + t49 * t87 + t50 * t88) + m(7) * (t37 * t47 + t45 * t53 + t46 * t54) + m(5) * (t106 * t52 + t119 * t83 + t120 * t84) + t286; m(7) * (t47 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t51 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t106 ^ 2 + t119 ^ 2 + t120 ^ 2) + t286; m(6) * (t233 * t92 + t235 * t91) + m(7) * (t233 * t86 + t235 * t85); m(6) * (t115 * t233 + t116 * t235 + t346 * t43) + m(7) * (t107 * t233 + t108 * t235 + t346 * t36); m(6) * (t233 * t50 + t235 * t49 + t346 * t44) + m(7) * (t233 * t46 + t235 * t45 + t346 * t37); m(7) * (t233 * t54 + t235 * t53 + t346 * t47) + m(6) * (t233 * t88 + t235 * t87 + t346 * t51); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t267 ^ 2 * t273 + t233 ^ 2 + t235 ^ 2); (-t279 * t86 - t282 * t85) * t357; m(7) * (t281 * t36 + (-t107 * t279 - t108 * t282) * t278); m(7) * (t281 * t37 + (-t279 * t46 - t282 * t45) * t278); m(7) * (t281 * t47 + (-t279 * t54 - t282 * t53) * t278); (-t233 * t279 - t235 * t282 + t267 * t281) * t357; m(7) * (t273 * t321 + t281 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
