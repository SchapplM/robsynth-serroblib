% Calculate joint inertia matrix for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:39
% EndTime: 2019-03-09 18:20:56
% DurationCPUTime: 6.95s
% Computational Cost: add. (10864->464), mult. (11487->662), div. (0->0), fcn. (12251->10), ass. (0->234)
t354 = Icges(4,4) + Icges(5,6);
t353 = Icges(4,1) + Icges(5,2);
t352 = -Icges(4,2) - Icges(5,3);
t231 = qJ(2) + qJ(3);
t221 = cos(t231);
t351 = t354 * t221;
t219 = sin(t231);
t350 = t354 * t219;
t349 = Icges(5,4) - Icges(4,5);
t348 = Icges(5,5) - Icges(4,6);
t347 = t352 * t219 + t351;
t346 = t353 * t221 - t350;
t345 = Icges(5,1) + Icges(4,3);
t234 = sin(qJ(1));
t237 = cos(qJ(1));
t344 = t347 * t234 + t348 * t237;
t343 = -t348 * t234 + t347 * t237;
t342 = t346 * t234 + t349 * t237;
t341 = -t349 * t234 + t346 * t237;
t340 = t348 * t219 - t349 * t221;
t339 = t352 * t221 - t350;
t338 = t353 * t219 + t351;
t337 = t234 * t345 + t340 * t237;
t336 = -t340 * t234 + t237 * t345;
t335 = t344 * t219 - t221 * t342;
t334 = t343 * t219 - t221 * t341;
t230 = qJ(5) + qJ(6);
t220 = cos(t230);
t302 = t220 * t237;
t218 = sin(t230);
t305 = t218 * t234;
t169 = t219 * t302 - t305;
t303 = t220 * t234;
t304 = t219 * t237;
t170 = t218 * t304 + t303;
t300 = t221 * t237;
t104 = t170 * rSges(7,1) + t169 * rSges(7,2) + rSges(7,3) * t300;
t235 = cos(qJ(5));
t215 = pkin(5) * t235 + pkin(4);
t238 = -pkin(10) - pkin(9);
t232 = sin(qJ(5));
t298 = t232 * t237;
t285 = t219 * t298;
t246 = pkin(5) * t285 + t234 * t215 - t238 * t300;
t288 = t234 * pkin(4) + pkin(9) * t300;
t122 = t246 - t288;
t333 = t104 + t122;
t171 = t218 * t237 + t219 * t303;
t172 = t219 * t305 - t302;
t270 = -t172 * rSges(7,1) - t171 * rSges(7,2);
t301 = t221 * t234;
t105 = rSges(7,3) * t301 - t270;
t227 = t237 * pkin(4);
t320 = -pkin(9) - t238;
t321 = pkin(5) * t232;
t123 = -t237 * t215 + t227 + (t219 * t321 + t221 * t320) * t234;
t332 = t105 + t123;
t331 = -t349 * t219 - t348 * t221;
t330 = t339 * t219 + t338 * t221;
t100 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t300;
t102 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t300;
t98 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t300;
t30 = t100 * t171 + t102 * t172 + t301 * t98;
t101 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t301;
t103 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t301;
t99 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t301;
t31 = t101 * t171 + t103 * t172 + t301 * t99;
t16 = t234 * t30 - t237 * t31;
t296 = t235 * t237;
t299 = t232 * t234;
t181 = t219 * t296 - t299;
t297 = t234 * t235;
t182 = t285 + t297;
t111 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t300;
t113 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t300;
t115 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t300;
t183 = t219 * t297 + t298;
t184 = t219 * t299 - t296;
t42 = t111 * t301 + t113 * t183 + t115 * t184;
t112 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t301;
t114 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t301;
t116 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t301;
t43 = t112 * t301 + t114 * t183 + t116 * t184;
t22 = t234 * t42 - t237 * t43;
t229 = t237 ^ 2;
t329 = -t16 - t22 + t336 * t229 + (t334 * t234 + (-t335 + t337) * t237) * t234;
t228 = t234 ^ 2;
t328 = m(5) / 0.2e1;
t327 = m(6) / 0.2e1;
t326 = m(7) / 0.2e1;
t325 = t219 / 0.2e1;
t324 = t234 / 0.2e1;
t323 = -t237 / 0.2e1;
t233 = sin(qJ(2));
t322 = pkin(2) * t233;
t236 = cos(qJ(2));
t319 = rSges(3,1) * t236;
t318 = rSges(3,2) * t233;
t317 = t237 * rSges(3,3);
t36 = t219 * t98 + (-t100 * t220 - t102 * t218) * t221;
t316 = t36 * t234;
t37 = t219 * t99 + (-t101 * t220 - t103 * t218) * t221;
t315 = t37 * t237;
t51 = t111 * t219 + (-t113 * t235 - t115 * t232) * t221;
t314 = t51 * t234;
t52 = t112 * t219 + (-t114 * t235 - t116 * t232) * t221;
t313 = t52 * t237;
t312 = Icges(3,4) * t233;
t311 = Icges(3,4) * t236;
t306 = qJ(4) * t219;
t239 = -pkin(8) - pkin(7);
t295 = t237 * t239;
t136 = rSges(7,3) * t219 + (-rSges(7,1) * t218 - rSges(7,2) * t220) * t221;
t163 = t219 * t320 - t221 * t321;
t294 = -t136 - t163;
t216 = pkin(2) * t236 + pkin(1);
t207 = t237 * t216;
t226 = t237 * pkin(7);
t293 = t234 * (t295 + t226 + (-pkin(1) + t216) * t234) + t237 * (-t237 * pkin(1) + t207 + (-pkin(7) - t239) * t234);
t249 = rSges(4,1) * t300 - rSges(4,2) * t304 + t234 * rSges(4,3);
t272 = rSges(4,1) * t221 - rSges(4,2) * t219;
t108 = t234 * (-t237 * rSges(4,3) + t234 * t272) + t237 * t249;
t290 = pkin(3) * t300 + qJ(4) * t304;
t292 = t228 * (pkin(3) * t221 + t306) + t237 * t290;
t194 = pkin(3) * t219 - qJ(4) * t221;
t291 = rSges(5,2) * t219 + rSges(5,3) * t221 - t194;
t289 = t234 * rSges(3,3) + t237 * t319;
t287 = t228 + t229;
t28 = t100 * t169 + t102 * t170 + t300 * t98;
t29 = t101 * t169 + t103 * t170 + t300 * t99;
t133 = Icges(7,3) * t219 + (-Icges(7,5) * t218 - Icges(7,6) * t220) * t221;
t134 = Icges(7,6) * t219 + (-Icges(7,4) * t218 - Icges(7,2) * t220) * t221;
t135 = Icges(7,5) * t219 + (-Icges(7,1) * t218 - Icges(7,4) * t220) * t221;
t58 = t133 * t300 + t134 * t169 + t135 * t170;
t5 = t219 * t58 + (t234 * t29 + t237 * t28) * t221;
t59 = t133 * t301 + t134 * t171 + t135 * t172;
t6 = t219 * t59 + (t234 * t31 + t237 * t30) * t221;
t129 = t219 * t133;
t260 = -t220 * t134 - t218 * t135;
t67 = (t221 * t260 + t129) * t219;
t286 = t5 * t300 + t6 * t301 + t219 * (t67 + (t234 * t37 + t237 * t36) * t221);
t119 = t182 * rSges(6,1) + t181 * rSges(6,2) + rSges(6,3) * t300;
t284 = t301 / 0.2e1;
t283 = t300 / 0.2e1;
t196 = rSges(4,1) * t219 + rSges(4,2) * t221;
t282 = -t196 - t322;
t281 = -pkin(9) * t219 - t194;
t15 = t234 * t28 - t237 * t29;
t40 = t111 * t300 + t113 * t181 + t115 * t182;
t41 = t112 * t300 + t114 * t181 + t116 * t182;
t21 = t234 * t40 - t237 * t41;
t280 = (t15 + t21 + t337 * t228 + ((-t334 + t336) * t234 + t335 * t237) * t237) * t234;
t279 = -t234 * t239 + t207;
t248 = t234 * rSges(5,1) - rSges(5,2) * t300 + rSges(5,3) * t304;
t73 = t234 * (-t237 * rSges(5,1) + (-rSges(5,2) * t221 + rSges(5,3) * t219) * t234) + t237 * t248 + t292;
t278 = t234 * (pkin(9) * t301 - t227) + t237 * t288 + t292;
t277 = t15 * t283 + t16 * t284 + t6 * t323 + t5 * t324 + (-t315 + t316) * t325;
t276 = t67 + (t37 + t59) * t284 + (t36 + t58) * t283;
t275 = t291 - t322;
t148 = rSges(6,3) * t219 + (-rSges(6,1) * t232 - rSges(6,2) * t235) * t221;
t274 = -t148 + t281;
t273 = -t318 + t319;
t271 = -t184 * rSges(6,1) - t183 * rSges(6,2);
t269 = Icges(3,1) * t236 - t312;
t267 = -Icges(3,2) * t233 + t311;
t264 = Icges(3,5) * t236 - Icges(3,6) * t233;
t146 = Icges(6,6) * t219 + (-Icges(6,4) * t232 - Icges(6,2) * t235) * t221;
t147 = Icges(6,5) * t219 + (-Icges(6,1) * t232 - Icges(6,4) * t235) * t221;
t259 = -t235 * t146 - t232 * t147;
t250 = t281 + t294;
t247 = t281 - t322;
t245 = t279 + t290;
t120 = rSges(6,3) * t301 - t271;
t46 = t237 * t119 + t234 * t120 + t278;
t244 = -t148 + t247;
t243 = t247 + t294;
t25 = t332 * t234 + t237 * t333 + t278;
t145 = Icges(6,3) * t219 + (-Icges(6,5) * t232 - Icges(6,6) * t235) * t221;
t63 = t145 * t300 + t146 * t181 + t147 * t182;
t10 = t219 * t63 + (t234 * t41 + t237 * t40) * t221;
t64 = t145 * t301 + t146 * t183 + t147 * t184;
t11 = t219 * t64 + (t234 * t43 + t237 * t42) * t221;
t242 = t10 * t324 + t11 * t323 + t21 * t283 + t22 * t284 + (-t313 + t314) * t325 + t277;
t241 = t329 * t237 + t280;
t240 = -t315 / 0.2e1 + t316 / 0.2e1 - t313 / 0.2e1 + t314 / 0.2e1 + (t341 * t219 + t221 * t343 + t331 * t234 + t330 * t237 + t58 + t63) * t324 + (t342 * t219 + t221 * t344 + t330 * t234 - t331 * t237 + t59 + t64) * t323;
t204 = rSges(2,1) * t237 - rSges(2,2) * t234;
t203 = -rSges(2,1) * t234 - rSges(2,2) * t237;
t202 = rSges(3,1) * t233 + rSges(3,2) * t236;
t174 = Icges(3,3) * t234 + t237 * t264;
t173 = -Icges(3,3) * t237 + t234 * t264;
t150 = t282 * t237;
t149 = t282 * t234;
t138 = t234 * pkin(7) + (pkin(1) - t318) * t237 + t289;
t137 = t317 + t226 + (-pkin(1) - t273) * t234;
t132 = t219 * t145;
t131 = t291 * t237;
t130 = t291 * t234;
t128 = t136 * t301;
t127 = t249 + t279;
t126 = (rSges(4,3) - t239) * t237 + (-t216 - t272) * t234;
t125 = t275 * t237;
t124 = t275 * t234;
t121 = t237 * (-t237 * t318 + t289) + (t234 * t273 - t317) * t234;
t107 = t274 * t237;
t106 = t274 * t234;
t97 = t245 + t248;
t96 = (rSges(5,1) - t239) * t237 + (-t216 + (rSges(5,2) - pkin(3)) * t221 + (-rSges(5,3) - qJ(4)) * t219) * t234;
t93 = t219 * t104;
t92 = t244 * t237;
t91 = t244 * t234;
t86 = t105 * t300;
t81 = t250 * t237;
t80 = t250 * t234;
t79 = t119 * t219 - t148 * t300;
t78 = -t120 * t219 + t148 * t301;
t77 = t243 * t237;
t76 = t243 * t234;
t75 = t245 + t119 + t288;
t74 = -t295 + t227 + (-t306 - t216 + (-rSges(6,3) - pkin(3) - pkin(9)) * t221) * t234 + t271;
t72 = -t136 * t300 + t93;
t71 = -t105 * t219 + t128;
t70 = (t221 * t259 + t132) * t219;
t69 = t108 + t293;
t68 = (-t119 * t234 + t120 * t237) * t221;
t66 = t245 + t246 + t104;
t65 = (t215 - t239) * t237 + (-t216 + (-qJ(4) - t321) * t219 + (-rSges(7,3) - pkin(3) + t238) * t221) * t234 + t270;
t60 = -t104 * t301 + t86;
t53 = t73 + t293;
t48 = t122 * t219 + t294 * t300 + t93;
t47 = t163 * t301 - t219 * t332 + t128;
t27 = t86 + (t123 * t237 - t234 * t333) * t221;
t26 = t46 + t293;
t23 = t25 + t293;
t1 = [t236 * (Icges(3,2) * t236 + t312) + t233 * (Icges(3,1) * t233 + t311) + Icges(2,3) + t129 + t132 + t338 * t219 + (t259 + t260 - t339) * t221 + m(7) * (t65 ^ 2 + t66 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t137 ^ 2 + t138 ^ 2) + m(2) * (t203 ^ 2 + t204 ^ 2); ((-Icges(3,6) * t237 + t234 * t267) * t236 + (-Icges(3,5) * t237 + t234 * t269) * t233) * t323 + ((Icges(3,6) * t234 + t237 * t267) * t236 + (Icges(3,5) * t234 + t237 * t269) * t233) * t324 + t240 + m(3) * (-t137 * t237 - t138 * t234) * t202 + m(7) * (t65 * t77 + t66 * t76) + m(6) * (t74 * t92 + t75 * t91) + m(5) * (t124 * t97 + t125 * t96) + m(4) * (t126 * t150 + t127 * t149) + (t228 / 0.2e1 + t229 / 0.2e1) * (Icges(3,5) * t233 + Icges(3,6) * t236); m(7) * (t23 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(6) * (t26 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t124 ^ 2 + t125 ^ 2 + t53 ^ 2) + m(4) * (t149 ^ 2 + t150 ^ 2 + t69 ^ 2) + m(3) * (t202 ^ 2 * t287 + t121 ^ 2) + t234 * t228 * t174 + t280 + (-t173 * t229 + (-t234 * t173 + t174 * t237) * t234 + t329) * t237; m(4) * (-t126 * t237 - t127 * t234) * t196 + m(7) * (t65 * t81 + t66 * t80) + m(6) * (t106 * t75 + t107 * t74) + m(5) * (t130 * t97 + t131 * t96) + t240; m(7) * (t25 * t23 + t76 * t80 + t77 * t81) + m(6) * (t106 * t91 + t107 * t92 + t46 * t26) + m(5) * (t124 * t130 + t125 * t131 + t53 * t73) + m(4) * (t108 * t69 + (-t149 * t234 - t150 * t237) * t196) + t241; m(7) * (t25 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(6) * (t106 ^ 2 + t107 ^ 2 + t46 ^ 2) + m(5) * (t130 ^ 2 + t131 ^ 2 + t73 ^ 2) + m(4) * (t196 ^ 2 * t287 + t108 ^ 2) + t241; 0.2e1 * ((t234 * t66 + t237 * t65) * t326 + (t234 * t75 + t237 * t74) * t327 + (t234 * t97 + t237 * t96) * t328) * t219; m(7) * (-t221 * t23 + (t234 * t76 + t237 * t77) * t219) + m(6) * (-t221 * t26 + (t234 * t91 + t237 * t92) * t219) + m(5) * (-t221 * t53 + (t124 * t234 + t125 * t237) * t219); m(7) * (-t221 * t25 + (t234 * t80 + t237 * t81) * t219) + m(6) * (-t221 * t46 + (t106 * t234 + t107 * t237) * t219) + m(5) * (-t221 * t73 + (t130 * t234 + t131 * t237) * t219); 0.2e1 * (t328 + t327 + t326) * (t219 ^ 2 * t287 + t221 ^ 2); t70 + m(7) * (t47 * t65 + t48 * t66) + m(6) * (t74 * t78 + t75 * t79) + ((t63 / 0.2e1 + t51 / 0.2e1) * t237 + (t64 / 0.2e1 + t52 / 0.2e1) * t234) * t221 + t276; m(7) * (t27 * t23 + t47 * t77 + t48 * t76) + m(6) * (t68 * t26 + t78 * t92 + t79 * t91) + t242; m(7) * (t27 * t25 + t47 * t81 + t48 * t80) + m(6) * (t106 * t79 + t107 * t78 + t68 * t46) + t242; m(6) * (-t221 * t68 + (t234 * t79 + t237 * t78) * t219) + m(7) * (-t221 * t27 + (t234 * t48 + t237 * t47) * t219); t219 * t70 + m(7) * (t27 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t68 ^ 2 + t78 ^ 2 + t79 ^ 2) + (t237 * t10 + t234 * t11 + t219 * (t234 * t52 + t237 * t51)) * t221 + t286; m(7) * (t71 * t65 + t66 * t72) + t276; m(7) * (t60 * t23 + t71 * t77 + t72 * t76) + t277; m(7) * (t60 * t25 + t71 * t81 + t72 * t80) + t277; m(7) * (-t221 * t60 + (t234 * t72 + t237 * t71) * t219); m(7) * (t60 * t27 + t71 * t47 + t48 * t72) + t286; m(7) * (t60 ^ 2 + t71 ^ 2 + t72 ^ 2) + t286;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
