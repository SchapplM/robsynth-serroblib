% Calculate joint inertia matrix for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:12:04
% EndTime: 2019-03-09 22:12:16
% DurationCPUTime: 5.82s
% Computational Cost: add. (14329->518), mult. (19933->728), div. (0->0), fcn. (23402->10), ass. (0->255)
t258 = qJ(2) + qJ(3);
t252 = cos(t258);
t264 = cos(qJ(4));
t266 = cos(qJ(1));
t326 = t264 * t266;
t260 = sin(qJ(4));
t262 = sin(qJ(1));
t329 = t260 * t262;
t224 = t252 * t329 + t326;
t327 = t262 * t264;
t328 = t260 * t266;
t225 = t252 * t327 - t328;
t251 = sin(t258);
t333 = t251 * t262;
t142 = Icges(6,5) * t225 + Icges(6,6) * t333 + Icges(6,3) * t224;
t148 = Icges(5,4) * t225 - Icges(5,2) * t224 + Icges(5,6) * t333;
t381 = t142 - t148;
t226 = t252 * t328 - t327;
t227 = t252 * t326 + t329;
t331 = t251 * t266;
t143 = Icges(6,5) * t227 + Icges(6,6) * t331 + Icges(6,3) * t226;
t149 = Icges(5,4) * t227 - Icges(5,2) * t226 + Icges(5,6) * t331;
t380 = t143 - t149;
t144 = Icges(5,5) * t225 - Icges(5,6) * t224 + Icges(5,3) * t333;
t146 = Icges(6,4) * t225 + Icges(6,2) * t333 + Icges(6,6) * t224;
t379 = t144 + t146;
t145 = Icges(5,5) * t227 - Icges(5,6) * t226 + Icges(5,3) * t331;
t147 = Icges(6,4) * t227 + Icges(6,2) * t331 + Icges(6,6) * t226;
t378 = t145 + t147;
t150 = Icges(6,1) * t225 + Icges(6,4) * t333 + Icges(6,5) * t224;
t152 = Icges(5,1) * t225 - Icges(5,4) * t224 + Icges(5,5) * t333;
t377 = t150 + t152;
t151 = Icges(6,1) * t227 + Icges(6,4) * t331 + Icges(6,5) * t226;
t153 = Icges(5,1) * t227 - Icges(5,4) * t226 + Icges(5,5) * t331;
t376 = t151 + t153;
t375 = t224 * t381 + t377 * t225 + t379 * t333;
t374 = t224 * t380 + t225 * t376 + t333 * t378;
t373 = t226 * t381 + t377 * t227 + t379 * t331;
t372 = t226 * t380 + t227 * t376 + t331 * t378;
t190 = -Icges(6,6) * t252 + (Icges(6,5) * t264 + Icges(6,3) * t260) * t251;
t192 = -Icges(6,2) * t252 + (Icges(6,4) * t264 + Icges(6,6) * t260) * t251;
t194 = -Icges(6,4) * t252 + (Icges(6,1) * t264 + Icges(6,5) * t260) * t251;
t94 = t190 * t224 + t192 * t333 + t194 * t225;
t191 = -Icges(5,3) * t252 + (Icges(5,5) * t264 - Icges(5,6) * t260) * t251;
t193 = -Icges(5,6) * t252 + (Icges(5,4) * t264 - Icges(5,2) * t260) * t251;
t195 = -Icges(5,5) * t252 + (Icges(5,1) * t264 - Icges(5,4) * t260) * t251;
t95 = t191 * t333 - t193 * t224 + t195 * t225;
t371 = -t94 - t95;
t96 = t190 * t226 + t192 * t331 + t194 * t227;
t97 = t191 * t331 - t193 * t226 + t195 * t227;
t370 = -t96 - t97;
t369 = t371 * t252 + (t262 * t375 + t374 * t266) * t251;
t368 = t370 * t252 + (t262 * t373 + t266 * t372) * t251;
t259 = sin(qJ(6));
t263 = cos(qJ(6));
t168 = t224 * t263 - t225 * t259;
t169 = t224 * t259 + t225 * t263;
t287 = -t169 * rSges(7,1) - t168 * rSges(7,2);
t114 = -rSges(7,3) * t333 - t287;
t323 = t225 * pkin(5) - pkin(10) * t333 + t114;
t170 = t226 * t263 - t227 * t259;
t171 = t226 * t259 + t227 * t263;
t319 = t171 * rSges(7,1) + t170 * rSges(7,2);
t115 = -rSges(7,3) * t331 + t319;
t221 = t227 * pkin(5);
t367 = -pkin(10) * t331 + t115 + t221;
t366 = -t191 - t192;
t108 = Icges(7,5) * t169 + Icges(7,6) * t168 - Icges(7,3) * t333;
t110 = Icges(7,4) * t169 + Icges(7,2) * t168 - Icges(7,6) * t333;
t112 = Icges(7,1) * t169 + Icges(7,4) * t168 - Icges(7,5) * t333;
t36 = -t108 * t331 + t110 * t170 + t112 * t171;
t109 = Icges(7,5) * t171 + Icges(7,6) * t170 - Icges(7,3) * t331;
t111 = Icges(7,4) * t171 + Icges(7,2) * t170 - Icges(7,6) * t331;
t113 = Icges(7,1) * t171 + Icges(7,4) * t170 - Icges(7,5) * t331;
t37 = -t109 * t331 + t111 * t170 + t113 * t171;
t14 = t262 * t37 - t266 * t36;
t365 = t262 * t372 - t266 * t373 + t14;
t34 = -t108 * t333 + t110 * t168 + t112 * t169;
t35 = -t109 * t333 + t111 * t168 + t113 * t169;
t13 = t262 * t35 - t266 * t34;
t364 = t374 * t262 - t266 * t375 + t13;
t334 = t251 * t260;
t332 = t251 * t264;
t362 = t190 * t334 + (t194 + t195) * t332;
t363 = (-t193 * t334 + t366 * t252 + t362) * t252;
t281 = Icges(4,5) * t252 - Icges(4,6) * t251;
t200 = -Icges(4,3) * t266 + t281 * t262;
t201 = Icges(4,3) * t262 + t281 * t266;
t257 = t266 ^ 2;
t335 = Icges(4,4) * t252;
t283 = -Icges(4,2) * t251 + t335;
t203 = Icges(4,6) * t262 + t283 * t266;
t336 = Icges(4,4) * t251;
t285 = Icges(4,1) * t252 - t336;
t205 = Icges(4,5) * t262 + t285 * t266;
t279 = -t203 * t251 + t205 * t252;
t202 = -Icges(4,6) * t266 + t283 * t262;
t204 = -Icges(4,5) * t266 + t285 * t262;
t280 = t202 * t251 - t204 * t252;
t361 = -t257 * t200 - (t279 * t262 + (-t201 + t280) * t266) * t262 - t364;
t206 = (-t259 * t264 + t260 * t263) * t251;
t207 = (t259 * t260 + t263 * t264) * t251;
t131 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t252;
t132 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t252;
t133 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t252;
t56 = -t131 * t333 + t132 * t168 + t133 * t169;
t5 = -t56 * t252 + (t262 * t34 + t266 * t35) * t251;
t57 = -t131 * t331 + t132 * t170 + t133 * t171;
t6 = -t57 * t252 + (t262 * t36 + t266 * t37) * t251;
t47 = t108 * t252 + t110 * t206 + t112 * t207;
t48 = t109 * t252 + t111 * t206 + t113 * t207;
t307 = t252 * t131 + t206 * t132 + t207 * t133;
t59 = t307 * t252;
t7 = -t59 + (t47 * t262 + t48 * t266) * t251;
t360 = (t262 * t5 + t266 * t6) * t251 - t252 * t7;
t343 = t48 * t262;
t344 = t47 * t266;
t359 = t252 * (t343 - t344) / 0.2e1 + t266 * t5 / 0.2e1 - t262 * t6 / 0.2e1;
t256 = t262 ^ 2;
t357 = t262 / 0.2e1;
t356 = -t266 / 0.2e1;
t355 = -rSges(7,3) - pkin(10);
t261 = sin(qJ(2));
t354 = pkin(2) * t261;
t353 = pkin(3) * t252;
t265 = cos(qJ(2));
t349 = rSges(3,1) * t265;
t348 = rSges(3,2) * t261;
t347 = t224 * rSges(6,3);
t345 = t266 * rSges(3,3);
t79 = -t252 * t146 + (t142 * t260 + t150 * t264) * t251;
t342 = t79 * t266;
t80 = -t252 * t147 + (t143 * t260 + t151 * t264) * t251;
t341 = t80 * t262;
t81 = -t252 * t144 + (-t148 * t260 + t152 * t264) * t251;
t340 = t81 * t266;
t82 = -t252 * t145 + (-t149 * t260 + t153 * t264) * t251;
t339 = t82 * t262;
t338 = Icges(3,4) * t261;
t337 = Icges(3,4) * t265;
t330 = t252 * t266;
t267 = -pkin(8) - pkin(7);
t325 = t266 * t267;
t134 = rSges(7,1) * t207 + rSges(7,2) * t206 + rSges(7,3) * t252;
t322 = pkin(5) * t332 + pkin(10) * t252 + t134;
t156 = t227 * rSges(6,1) + rSges(6,2) * t331 + t226 * rSges(6,3);
t173 = t227 * pkin(4) + t226 * qJ(5);
t321 = -t156 - t173;
t216 = t224 * qJ(5);
t172 = t225 * pkin(4) + t216;
t223 = (pkin(4) * t264 + qJ(5) * t260) * t251;
t320 = t252 * t172 + t223 * t333;
t250 = pkin(2) * t265 + pkin(1);
t242 = t266 * t250;
t255 = t266 * pkin(7);
t317 = t262 * (t325 + t255 + (-pkin(1) + t250) * t262) + t266 * (-t266 * pkin(1) + t242 + (-pkin(7) - t267) * t262);
t274 = rSges(4,1) * t330 - rSges(4,2) * t331 + t262 * rSges(4,3);
t289 = rSges(4,1) * t252 - rSges(4,2) * t251;
t135 = t262 * (-t266 * rSges(4,3) + t289 * t262) + t266 * t274;
t196 = -t252 * rSges(6,2) + (rSges(6,1) * t264 + rSges(6,3) * t260) * t251;
t316 = -t196 - t223;
t197 = -t252 * rSges(5,3) + (rSges(5,1) * t264 - rSges(5,2) * t260) * t251;
t234 = pkin(3) * t251 - pkin(9) * t252;
t315 = -t197 - t234;
t313 = pkin(3) * t330 + pkin(9) * t331;
t314 = t256 * (pkin(9) * t251 + t353) + t266 * t313;
t312 = t262 * rSges(3,3) + t266 * t349;
t311 = t256 + t257;
t310 = -t48 / 0.2e1 - t57 / 0.2e1;
t309 = -t56 / 0.2e1 - t47 / 0.2e1;
t308 = -t173 - t367;
t306 = -t223 - t322;
t305 = -t234 + t316;
t157 = t227 * rSges(5,1) - t226 * rSges(5,2) + rSges(5,3) * t331;
t233 = rSges(4,1) * t251 + rSges(4,2) * t252;
t302 = -t233 - t354;
t301 = -t234 - t354;
t300 = -t250 - t353;
t299 = (t256 * t201 + (t280 * t266 + (-t200 + t279) * t262) * t266 + t365) * t262;
t298 = -t216 - t325;
t297 = -t262 * t267 + t242;
t296 = -t14 * t331 / 0.2e1 - t13 * t333 / 0.2e1 + t359;
t295 = -t234 + t306;
t294 = t262 * t172 + t266 * t173 + t314;
t288 = -t225 * rSges(5,1) + t224 * rSges(5,2);
t155 = rSges(5,3) * t333 - t288;
t87 = t262 * t155 + t266 * t157 + t314;
t293 = -t197 + t301;
t292 = -t223 + t301;
t290 = -t348 + t349;
t286 = Icges(3,1) * t265 - t338;
t284 = -Icges(3,2) * t261 + t337;
t282 = Icges(3,5) * t265 - Icges(3,6) * t261;
t231 = Icges(4,2) * t252 + t336;
t232 = Icges(4,1) * t251 + t335;
t276 = -t231 * t251 + t232 * t252;
t275 = -t196 + t292;
t273 = t297 + t313;
t154 = t225 * rSges(6,1) + rSges(6,2) * t333 + t347;
t58 = t262 * t154 + t266 * t156 + t294;
t272 = t292 - t322;
t31 = t323 * t262 + t367 * t266 + t294;
t271 = t266 * t361 + t299;
t270 = t173 + t273;
t269 = -t359 - (t341 - t342 + t339 - t340) * t252 / 0.2e1 + t368 * t357 + t369 * t356 + t364 * t333 / 0.2e1 + t365 * t331 / 0.2e1;
t230 = Icges(4,5) * t251 + Icges(4,6) * t252;
t268 = -t344 / 0.2e1 + t343 / 0.2e1 - t342 / 0.2e1 + t341 / 0.2e1 - t340 / 0.2e1 + t339 / 0.2e1 + (t203 * t252 + t205 * t251 + t262 * t230 + t266 * t276 - t370 + t57) * t357 + (t202 * t252 + t204 * t251 - t266 * t230 + t262 * t276 - t371 + t56) * t356;
t241 = rSges(2,1) * t266 - rSges(2,2) * t262;
t240 = -rSges(2,1) * t262 - rSges(2,2) * t266;
t239 = rSges(3,1) * t261 + rSges(3,2) * t265;
t211 = Icges(3,3) * t262 + t282 * t266;
t210 = -Icges(3,3) * t266 + t282 * t262;
t199 = t302 * t266;
t198 = t302 * t262;
t184 = t262 * pkin(7) + (pkin(1) - t348) * t266 + t312;
t183 = t345 + t255 + (-pkin(1) - t290) * t262;
t175 = t274 + t297;
t174 = (rSges(4,3) - t267) * t266 + (-t250 - t289) * t262;
t161 = t315 * t266;
t160 = t315 * t262;
t159 = t266 * (-t266 * t348 + t312) + (t290 * t262 - t345) * t262;
t158 = t172 * t331;
t137 = t293 * t266;
t136 = t293 * t262;
t125 = t305 * t266;
t124 = t305 * t262;
t123 = t275 * t266;
t122 = t275 * t262;
t119 = t273 + t157;
t118 = -t325 + ((-rSges(5,3) - pkin(9)) * t251 + t300) * t262 + t288;
t117 = -t157 * t252 - t197 * t331;
t116 = t155 * t252 + t197 * t333;
t103 = t135 + t317;
t102 = (t155 * t266 - t157 * t262) * t251;
t101 = t295 * t266;
t100 = t295 * t262;
t99 = t270 + t156;
t98 = -t347 + (-rSges(6,1) - pkin(4)) * t225 + ((-rSges(6,2) - pkin(9)) * t251 + t300) * t262 + t298;
t93 = t272 * t266;
t92 = t272 * t262;
t86 = t321 * t252 + t316 * t331;
t85 = t154 * t252 + t196 * t333 + t320;
t84 = t252 * t115 + t134 * t331;
t83 = -t252 * t114 - t134 * t333;
t64 = t355 * t331 + t221 + t270 + t319;
t63 = (-pkin(4) - pkin(5)) * t225 + ((-pkin(9) - t355) * t251 + t300) * t262 + t287 + t298;
t62 = t87 + t317;
t61 = t158 + (t154 * t266 + t321 * t262) * t251;
t60 = (-t114 * t266 + t115 * t262) * t251;
t51 = t58 + t317;
t50 = t252 * t308 + t306 * t331;
t49 = t252 * t323 + t322 * t333 + t320;
t42 = t158 + (t262 * t308 + t266 * t323) * t251;
t26 = t31 + t317;
t1 = [t265 * (Icges(3,2) * t265 + t338) + t261 * (Icges(3,1) * t261 + t337) + Icges(2,3) + (-t193 * t260 + t232) * t251 + (t231 + t366) * t252 + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t98 ^ 2 + t99 ^ 2) + m(5) * (t118 ^ 2 + t119 ^ 2) + m(4) * (t174 ^ 2 + t175 ^ 2) + m(3) * (t183 ^ 2 + t184 ^ 2) + m(2) * (t240 ^ 2 + t241 ^ 2) + t307 + t362; m(3) * (-t183 * t266 - t184 * t262) * t239 + (t256 / 0.2e1 + t257 / 0.2e1) * (Icges(3,5) * t261 + Icges(3,6) * t265) + t268 + m(7) * (t63 * t93 + t64 * t92) + m(6) * (t122 * t99 + t123 * t98) + m(5) * (t118 * t137 + t119 * t136) + m(4) * (t174 * t199 + t175 * t198) + (t265 * (Icges(3,6) * t262 + t284 * t266) + t261 * (Icges(3,5) * t262 + t286 * t266)) * t357 + (t265 * (-Icges(3,6) * t266 + t284 * t262) + t261 * (-Icges(3,5) * t266 + t286 * t262)) * t356; m(7) * (t26 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(6) * (t122 ^ 2 + t123 ^ 2 + t51 ^ 2) + m(5) * (t136 ^ 2 + t137 ^ 2 + t62 ^ 2) + m(4) * (t103 ^ 2 + t198 ^ 2 + t199 ^ 2) + t262 * t256 * t211 + m(3) * (t239 ^ 2 * t311 + t159 ^ 2) + t299 + (-t257 * t210 + (-t262 * t210 + t266 * t211) * t262 + t361) * t266; t268 + m(4) * (-t174 * t266 - t175 * t262) * t233 + m(7) * (t100 * t64 + t101 * t63) + m(6) * (t124 * t99 + t125 * t98) + m(5) * (t118 * t161 + t119 * t160); m(7) * (t100 * t92 + t101 * t93 + t31 * t26) + m(6) * (t122 * t124 + t123 * t125 + t58 * t51) + m(5) * (t136 * t160 + t137 * t161 + t62 * t87) + m(4) * (t135 * t103 + (-t198 * t262 - t199 * t266) * t233) + t271; m(7) * (t100 ^ 2 + t101 ^ 2 + t31 ^ 2) + m(6) * (t124 ^ 2 + t125 ^ 2 + t58 ^ 2) + m(5) * (t160 ^ 2 + t161 ^ 2 + t87 ^ 2) + m(4) * (t233 ^ 2 * t311 + t135 ^ 2) + t271; -t59 - t363 + m(7) * (t49 * t63 + t50 * t64) + m(6) * (t85 * t98 + t86 * t99) + m(5) * (t116 * t118 + t117 * t119) + ((t82 / 0.2e1 + t80 / 0.2e1 + t97 / 0.2e1 + t96 / 0.2e1 - t310) * t266 + (t79 / 0.2e1 + t95 / 0.2e1 + t94 / 0.2e1 + t81 / 0.2e1 - t309) * t262) * t251; t269 + m(7) * (t42 * t26 + t49 * t93 + t50 * t92) + m(6) * (t122 * t86 + t123 * t85 + t61 * t51) + m(5) * (t102 * t62 + t116 * t137 + t117 * t136); t269 + m(5) * (t102 * t87 + t116 * t161 + t117 * t160) + m(7) * (t100 * t50 + t101 * t49 + t42 * t31) + m(6) * (t124 * t86 + t125 * t85 + t61 * t58); m(7) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t61 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(5) * (t102 ^ 2 + t116 ^ 2 + t117 ^ 2) + (-t7 + t363) * t252 + ((t6 + (-t80 - t82) * t252 + t368) * t266 + (t5 + (-t79 - t81) * t252 + t369) * t262) * t251; m(7) * (t224 * t64 + t226 * t63) + m(6) * (t224 * t99 + t226 * t98); m(7) * (t224 * t92 + t226 * t93 + t26 * t334) + m(6) * (t122 * t224 + t123 * t226 + t334 * t51); m(7) * (t100 * t224 + t101 * t226 + t31 * t334) + m(6) * (t124 * t224 + t125 * t226 + t334 * t58); m(7) * (t224 * t50 + t226 * t49 + t334 * t42) + m(6) * (t224 * t86 + t226 * t85 + t334 * t61); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t251 ^ 2 * t260 ^ 2 + t224 ^ 2 + t226 ^ 2); m(7) * (t63 * t83 + t64 * t84) + t59 + (t262 * t309 + t266 * t310) * t251; m(7) * (t60 * t26 + t83 * t93 + t84 * t92) + t296; m(7) * (t100 * t84 + t101 * t83 + t60 * t31) + t296; m(7) * (t60 * t42 + t49 * t83 + t50 * t84) - t360; m(7) * (t224 * t84 + t226 * t83 + t334 * t60); m(7) * (t60 ^ 2 + t83 ^ 2 + t84 ^ 2) + t360;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
