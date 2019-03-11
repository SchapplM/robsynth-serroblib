% Calculate joint inertia matrix for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:10
% EndTime: 2019-03-09 18:25:24
% DurationCPUTime: 5.87s
% Computational Cost: add. (18910->591), mult. (17473->817), div. (0->0), fcn. (18790->12), ass. (0->289)
t288 = sin(qJ(2));
t386 = Icges(3,5) * t288;
t385 = t386 / 0.2e1;
t282 = qJ(3) + pkin(11);
t271 = sin(t282);
t272 = cos(t282);
t291 = cos(qJ(2));
t200 = -Icges(5,3) * t291 + (Icges(5,5) * t272 - Icges(5,6) * t271) * t288;
t287 = sin(qJ(3));
t290 = cos(qJ(3));
t218 = -Icges(4,3) * t291 + (Icges(4,5) * t290 - Icges(4,6) * t287) * t288;
t384 = -t200 - t218;
t201 = -Icges(5,6) * t291 + (Icges(5,4) * t272 - Icges(5,2) * t271) * t288;
t221 = -Icges(4,6) * t291 + (Icges(4,4) * t290 - Icges(4,2) * t287) * t288;
t383 = -t201 * t271 - t221 * t287;
t224 = -Icges(4,5) * t291 + (Icges(4,1) * t290 - Icges(4,4) * t287) * t288;
t292 = cos(qJ(1));
t289 = sin(qJ(1));
t350 = t289 * t291;
t235 = -t287 * t350 - t290 * t292;
t354 = t287 * t292;
t236 = t290 * t350 - t354;
t353 = t288 * t289;
t105 = t218 * t353 + t221 * t235 + t224 * t236;
t214 = -t271 * t350 - t272 * t292;
t215 = -t271 * t292 + t272 * t350;
t148 = Icges(5,5) * t215 + Icges(5,6) * t214 + Icges(5,3) * t353;
t150 = Icges(5,4) * t215 + Icges(5,2) * t214 + Icges(5,6) * t353;
t152 = Icges(5,1) * t215 + Icges(5,4) * t214 + Icges(5,5) * t353;
t58 = t148 * t353 + t150 * t214 + t152 * t215;
t349 = t291 * t292;
t216 = -t271 * t349 + t289 * t272;
t217 = t289 * t271 + t272 * t349;
t352 = t288 * t292;
t149 = Icges(5,5) * t217 + Icges(5,6) * t216 + Icges(5,3) * t352;
t151 = Icges(5,4) * t217 + Icges(5,2) * t216 + Icges(5,6) * t352;
t153 = Icges(5,1) * t217 + Icges(5,4) * t216 + Icges(5,5) * t352;
t59 = t149 * t353 + t151 * t214 + t153 * t215;
t162 = Icges(4,5) * t236 + Icges(4,6) * t235 + Icges(4,3) * t353;
t164 = Icges(4,4) * t236 + Icges(4,2) * t235 + Icges(4,6) * t353;
t166 = Icges(4,1) * t236 + Icges(4,4) * t235 + Icges(4,5) * t353;
t70 = t162 * t353 + t164 * t235 + t166 * t236;
t237 = -t287 * t349 + t289 * t290;
t351 = t289 * t287;
t238 = t290 * t349 + t351;
t163 = Icges(4,5) * t238 + Icges(4,6) * t237 + Icges(4,3) * t352;
t165 = Icges(4,4) * t238 + Icges(4,2) * t237 + Icges(4,6) * t352;
t167 = Icges(4,1) * t238 + Icges(4,4) * t237 + Icges(4,5) * t352;
t71 = t163 * t353 + t165 * t235 + t167 * t236;
t202 = -Icges(5,5) * t291 + (Icges(5,1) * t272 - Icges(5,4) * t271) * t288;
t91 = t200 * t353 + t201 * t214 + t202 * t215;
t382 = (-t91 - t105) * t291 + ((t59 + t71) * t292 + (t58 + t70) * t289) * t288;
t106 = t218 * t352 + t237 * t221 + t238 * t224;
t60 = t148 * t352 + t216 * t150 + t217 * t152;
t61 = t149 * t352 + t216 * t151 + t217 * t153;
t72 = t162 * t352 + t237 * t164 + t238 * t166;
t73 = t163 * t352 + t237 * t165 + t238 * t167;
t92 = t200 * t352 + t216 * t201 + t217 * t202;
t381 = (-t92 - t106) * t291 + ((t61 + t73) * t292 + (t60 + t72) * t289) * t288;
t66 = -t148 * t291 + (-t150 * t271 + t152 * t272) * t288;
t80 = -t162 * t291 + (-t164 * t287 + t166 * t290) * t288;
t380 = -t66 - t80;
t67 = -t149 * t291 + (-t151 * t271 + t153 * t272) * t288;
t81 = -t163 * t291 + (-t165 * t287 + t167 * t290) * t288;
t379 = t67 + t81;
t273 = qJ(5) + t282;
t269 = qJ(6) + t273;
t254 = sin(t269);
t255 = cos(t269);
t197 = -t254 * t349 + t289 * t255;
t198 = t289 * t254 + t255 * t349;
t136 = t198 * rSges(7,1) + t197 * rSges(7,2) + rSges(7,3) * t352;
t270 = t290 * pkin(3) + pkin(2);
t241 = pkin(4) * t272 + t270;
t268 = cos(t273);
t213 = pkin(5) * t268 + t241;
t243 = t287 * pkin(3) + pkin(4) * t271;
t267 = sin(t273);
t228 = pkin(5) * t267 + t243;
t378 = t213 * t349 + t289 * t228 + t136;
t377 = (t202 * t272 + t224 * t290) * t288;
t284 = t289 ^ 2;
t285 = t292 ^ 2;
t376 = m(5) / 0.2e1;
t375 = m(6) / 0.2e1;
t374 = m(7) / 0.2e1;
t195 = -t254 * t350 - t255 * t292;
t196 = -t254 * t292 + t255 * t350;
t129 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t353;
t131 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t353;
t133 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t353;
t44 = t129 * t353 + t131 * t195 + t133 * t196;
t130 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t352;
t132 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t352;
t134 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t352;
t45 = t130 * t353 + t132 * t195 + t134 * t196;
t181 = -Icges(7,3) * t291 + (Icges(7,5) * t255 - Icges(7,6) * t254) * t288;
t182 = -Icges(7,6) * t291 + (Icges(7,4) * t255 - Icges(7,2) * t254) * t288;
t183 = -Icges(7,5) * t291 + (Icges(7,1) * t255 - Icges(7,4) * t254) * t288;
t78 = t181 * t353 + t182 * t195 + t183 * t196;
t5 = -t78 * t291 + (t289 * t44 + t292 * t45) * t288;
t46 = t129 * t352 + t197 * t131 + t198 * t133;
t47 = t130 * t352 + t197 * t132 + t198 * t134;
t79 = t181 * t352 + t197 * t182 + t198 * t183;
t6 = -t79 * t291 + (t289 * t46 + t292 * t47) * t288;
t373 = t6 * t352 + t5 * t353;
t372 = t289 / 0.2e1;
t371 = -t291 / 0.2e1;
t370 = -t292 / 0.2e1;
t369 = pkin(2) * t291;
t368 = pkin(8) * t288;
t367 = -pkin(2) + t270;
t286 = -qJ(4) - pkin(8);
t189 = -Icges(6,5) * t291 + (Icges(6,1) * t268 - Icges(6,4) * t267) * t288;
t174 = t288 * t268 * t189;
t187 = -Icges(6,3) * t291 + (Icges(6,5) * t268 - Icges(6,6) * t267) * t288;
t188 = -Icges(6,6) * t291 + (Icges(6,4) * t268 - Icges(6,2) * t267) * t288;
t358 = t188 * t267;
t100 = -t291 * t187 - t288 * t358 + t174;
t172 = t288 * t255 * t183;
t359 = t182 * t254;
t93 = -t291 * t181 - t288 * t359 + t172;
t364 = t93 * t291;
t56 = -t129 * t291 + (-t131 * t254 + t133 * t255) * t288;
t57 = -t130 * t291 + (-t132 * t254 + t134 * t255) * t288;
t15 = -t364 + (t289 * t56 + t292 * t57) * t288;
t207 = -t267 * t350 - t268 * t292;
t208 = -t267 * t292 + t268 * t350;
t139 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t353;
t141 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t353;
t143 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t353;
t64 = -t139 * t291 + (-t141 * t267 + t143 * t268) * t288;
t209 = -t267 * t349 + t289 * t268;
t210 = t289 * t267 + t268 * t349;
t140 = Icges(6,5) * t210 + Icges(6,6) * t209 + Icges(6,3) * t352;
t142 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t352;
t144 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t352;
t65 = -t140 * t291 + (-t142 * t267 + t144 * t268) * t288;
t366 = -t15 + t100 * t291 - (t289 * t64 + t292 * t65) * t288;
t365 = t292 * rSges(3,3);
t281 = -pkin(9) + t286;
t274 = -pkin(10) + t281;
t334 = t292 * t243 + t281 * t353;
t338 = t213 - t241;
t355 = t228 * t292;
t101 = -t355 + (-t274 * t288 + t291 * t338) * t289 + t334;
t307 = -t196 * rSges(7,1) - t195 * rSges(7,2);
t135 = rSges(7,3) * t353 - t307;
t124 = t135 * t352;
t363 = t101 * t352 + t124;
t362 = -t93 - t100;
t360 = Icges(3,4) * t291;
t329 = t274 - t281;
t336 = t241 * t349 + t289 * t243;
t348 = -t329 * t352 - t336 + t378;
t331 = pkin(3) * t354 + t286 * t353;
t333 = t241 - t270;
t120 = t333 * t350 + t331 - t334;
t168 = (t291 * t367 - t368) * t289 - t331;
t156 = t168 * t352;
t347 = t120 * t352 + t156;
t346 = t383 * t288 + t291 * t384 + t377;
t327 = t281 - t286;
t332 = -pkin(3) * t351 - t270 * t349;
t121 = -t327 * t352 + t332 + t336;
t299 = -t286 * t352 - t332;
t330 = pkin(2) * t349 + pkin(8) * t352;
t169 = t299 - t330;
t345 = -t121 - t169;
t147 = t288 * t338 + t291 * t329;
t184 = -rSges(7,3) * t291 + (rSges(7,1) * t255 - rSges(7,2) * t254) * t288;
t344 = -t147 - t184;
t155 = t217 * rSges(5,1) + t216 * rSges(5,2) + rSges(5,3) * t352;
t343 = -t155 - t169;
t199 = (pkin(8) + t286) * t291 + t367 * t288;
t342 = t291 * t168 + t199 * t353;
t173 = t288 * t333 + t291 * t327;
t341 = -t173 - t199;
t103 = t291 * t135 + t184 * t353;
t308 = -rSges(6,1) * t208 - rSges(6,2) * t207;
t145 = rSges(6,3) * t353 - t308;
t192 = -rSges(6,3) * t291 + (rSges(6,1) * t268 - rSges(6,2) * t267) * t288;
t109 = t291 * t145 + t192 * t353;
t205 = -rSges(5,3) * t291 + (rSges(5,1) * t272 - rSges(5,2) * t271) * t288;
t339 = -t199 - t205;
t227 = -rSges(4,3) * t291 + (rSges(4,1) * t290 - rSges(4,2) * t287) * t288;
t252 = pkin(2) * t288 - pkin(8) * t291;
t337 = -t227 - t252;
t335 = t284 * (t368 + t369) + t292 * t330;
t328 = t292 * pkin(1) + t289 * pkin(7);
t326 = t284 + t285;
t146 = t210 * rSges(6,1) + t209 * rSges(6,2) + rSges(6,3) * t352;
t325 = -t146 + t345;
t324 = -t192 + t341;
t323 = -t252 + t339;
t171 = t238 * rSges(4,1) + t237 * rSges(4,2) + rSges(4,3) * t352;
t49 = t139 * t353 + t141 * t207 + t143 * t208;
t50 = t140 * t353 + t142 * t207 + t144 * t208;
t85 = t187 * t353 + t188 * t207 + t189 * t208;
t13 = -t85 * t291 + (t289 * t49 + t292 * t50) * t288;
t51 = t139 * t352 + t209 * t141 + t210 * t143;
t52 = t140 * t352 + t209 * t142 + t210 * t144;
t86 = t187 * t352 + t209 * t188 + t210 * t189;
t14 = -t86 * t291 + (t289 * t51 + t292 * t52) * t288;
t322 = t13 * t353 + t14 * t352 + t373;
t321 = t353 / 0.2e1;
t320 = t352 / 0.2e1;
t319 = (t56 + t78) * t321 + (t57 + t79) * t320;
t42 = t291 * t101 + t147 * t353 + t103;
t318 = -t291 * t15 + t373;
t317 = t345 - t348;
t316 = t291 * t120 + t173 * t353 + t342;
t315 = t341 + t344;
t314 = t289 * t168 + t292 * t169 + t335;
t313 = -t252 + t324;
t21 = t45 * t289 - t292 * t44;
t22 = t47 * t289 - t292 * t46;
t312 = t21 * t321 + t22 * t320 + t5 * t370 + t6 * t372 + (t57 * t289 - t56 * t292) * t371;
t311 = rSges(3,1) * t291 - rSges(3,2) * t288;
t310 = -rSges(4,1) * t236 - rSges(4,2) * t235;
t309 = -rSges(5,1) * t215 - rSges(5,2) * t214;
t306 = -t252 + t315;
t304 = -Icges(3,2) * t288 + t360;
t303 = Icges(3,5) * t291 - Icges(3,6) * t288;
t300 = rSges(3,1) * t349 - rSges(3,2) * t352 + t289 * rSges(3,3);
t298 = t105 / 0.2e1 + t91 / 0.2e1 + t80 / 0.2e1 + t66 / 0.2e1;
t297 = t106 / 0.2e1 + t92 / 0.2e1 + t81 / 0.2e1 + t67 / 0.2e1;
t296 = t289 * t120 + t292 * t121 + t314;
t295 = t291 * t366 + t322;
t294 = t319 + (t64 + t85) * t321 + (t65 + t86) * t320;
t27 = t50 * t289 - t292 * t49;
t28 = t52 * t289 - t292 * t51;
t293 = t13 * t370 + t14 * t372 + t27 * t321 + t28 * t320 + t312 + (t65 * t289 - t64 * t292) * t371;
t279 = t292 * pkin(7);
t251 = rSges(2,1) * t292 - t289 * rSges(2,2);
t250 = -t289 * rSges(2,1) - rSges(2,2) * t292;
t249 = rSges(3,1) * t288 + rSges(3,2) * t291;
t245 = Icges(3,6) * t291 + t386;
t220 = Icges(3,3) * t289 + t292 * t303;
t219 = -Icges(3,3) * t292 + t289 * t303;
t186 = t300 + t328;
t185 = t365 + t279 + (-pkin(1) - t311) * t289;
t176 = t337 * t292;
t175 = t337 * t289;
t170 = rSges(4,3) * t353 - t310;
t157 = t292 * t300 + (t289 * t311 - t365) * t289;
t154 = rSges(5,3) * t353 - t309;
t128 = t145 * t352;
t127 = t171 + t328 + t330;
t126 = t279 + (-t369 - pkin(1) + (-rSges(4,3) - pkin(8)) * t288) * t289 + t310;
t123 = t323 * t292;
t122 = t323 * t289;
t119 = -t291 * t171 - t227 * t352;
t118 = t170 * t291 + t227 * t353;
t112 = t299 + t155 + t328;
t111 = t279 + (-rSges(5,3) * t288 - t270 * t291 - pkin(1)) * t289 + t309 + t331;
t110 = -t291 * t146 - t192 * t352;
t107 = (t170 * t292 - t171 * t289) * t288;
t104 = -t291 * t136 - t184 * t352;
t98 = -t281 * t352 + t146 + t328 + t336;
t97 = t279 + (-rSges(6,3) * t288 - t241 * t291 - pkin(1)) * t289 + t308 + t334;
t95 = t313 * t292;
t94 = t313 * t289;
t90 = -t146 * t353 + t128;
t89 = t289 * t170 + t171 * t292 + t335;
t88 = -t274 * t352 + t328 + t378;
t87 = t355 + t279 + (-t213 * t291 - pkin(1) + (-rSges(7,3) + t274) * t288) * t289 + t307;
t84 = -t136 * t353 + t124;
t75 = t291 * t343 + t339 * t352;
t74 = t154 * t291 + t205 * t353 + t342;
t69 = t306 * t292;
t68 = t306 * t289;
t55 = t156 + (t154 * t292 + t289 * t343) * t288;
t48 = t289 * t154 + t155 * t292 + t314;
t43 = -t291 * t348 + t344 * t352;
t41 = t291 * t325 + t324 * t352;
t40 = t316 + t109;
t39 = t73 * t289 - t292 * t72;
t38 = t71 * t289 - t292 * t70;
t37 = -t348 * t353 + t363;
t36 = t325 * t353 + t128 + t347;
t35 = t289 * t145 + t146 * t292 + t296;
t33 = t61 * t289 - t292 * t60;
t32 = t59 * t289 - t292 * t58;
t24 = t291 * t317 + t315 * t352;
t23 = t316 + t42;
t8 = t317 * t353 + t347 + t363;
t7 = t348 * t292 + (t101 + t135) * t289 + t296;
t1 = [Icges(2,3) + t172 + t174 + (Icges(3,4) * t288 + Icges(3,2) * t291 - t181 - t187 + t384) * t291 + (Icges(3,1) * t288 - t358 - t359 + t360 + t383) * t288 + m(7) * (t87 ^ 2 + t88 ^ 2) + m(6) * (t97 ^ 2 + t98 ^ 2) + m(5) * (t111 ^ 2 + t112 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t185 ^ 2 + t186 ^ 2) + m(2) * (t250 ^ 2 + t251 ^ 2) + t377; (-t64 / 0.2e1 - t56 / 0.2e1 - t78 / 0.2e1 - t85 / 0.2e1 + t289 * t304 * t371 - t298 + (-Icges(3,6) * t371 + t385 + t245 / 0.2e1) * t292) * t292 + (t65 / 0.2e1 + t57 / 0.2e1 + t79 / 0.2e1 + t86 / 0.2e1 + (Icges(3,6) * t289 + t292 * t304) * t291 / 0.2e1 + t289 * t385 + t245 * t372 + t297) * t289 + m(6) * (t94 * t98 + t95 * t97) + m(7) * (t68 * t88 + t69 * t87) + m(5) * (t111 * t123 + t112 * t122) + m(4) * (t126 * t176 + t127 * t175) + m(3) * (-t185 * t292 - t186 * t289) * t249; m(7) * (t68 ^ 2 + t69 ^ 2 + t7 ^ 2) + m(6) * (t35 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(5) * (t122 ^ 2 + t123 ^ 2 + t48 ^ 2) + m(4) * (t175 ^ 2 + t176 ^ 2 + t89 ^ 2) + m(3) * (t249 ^ 2 * t326 + t157 ^ 2) + (-t285 * t219 - t21 - t27 - t32 - t38) * t292 + (t284 * t220 + t22 + t28 + t33 + t39 + (-t289 * t219 + t292 * t220) * t292) * t289; (-t346 + t362) * t291 + m(7) * (t23 * t87 + t24 * t88) + m(6) * (t40 * t97 + t41 * t98) + m(5) * (t111 * t74 + t112 * t75) + m(4) * (t118 * t126 + t119 * t127) + (t289 * t298 + t292 * t297) * t288 + t294; ((t33 / 0.2e1 + t39 / 0.2e1) * t292 + (t38 / 0.2e1 + t32 / 0.2e1) * t289) * t288 + m(7) * (t23 * t69 + t24 * t68 + t7 * t8) + m(6) * (t35 * t36 + t40 * t95 + t41 * t94) + m(5) * (t122 * t75 + t123 * t74 + t48 * t55) + m(4) * (t107 * t89 + t118 * t176 + t119 * t175) + t293 + t381 * t372 + (t289 * t379 + t292 * t380) * t371 + t382 * t370; (t346 * t291 + t366) * t291 + m(5) * (t55 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t107 ^ 2 + t118 ^ 2 + t119 ^ 2) + m(7) * (t23 ^ 2 + t24 ^ 2 + t8 ^ 2) + m(6) * (t36 ^ 2 + t40 ^ 2 + t41 ^ 2) + ((-t291 * t379 + t381) * t292 + (t291 * t380 + t382) * t289) * t288 + t322; 0.2e1 * ((t289 * t88 + t292 * t87) * t374 + (t289 * t98 + t292 * t97) * t375 + (t111 * t292 + t112 * t289) * t376) * t288; m(7) * (-t291 * t7 + (t289 * t68 + t292 * t69) * t288) + m(6) * (-t291 * t35 + (t289 * t94 + t292 * t95) * t288) + m(5) * (-t291 * t48 + (t122 * t289 + t123 * t292) * t288); m(7) * (-t291 * t8 + (t23 * t292 + t24 * t289) * t288) + m(6) * (-t291 * t36 + (t289 * t41 + t292 * t40) * t288) + m(5) * (-t291 * t55 + (t289 * t75 + t292 * t74) * t288); 0.2e1 * (t376 + t375 + t374) * (t288 ^ 2 * t326 + t291 ^ 2); t362 * t291 + m(7) * (t42 * t87 + t43 * t88) + m(6) * (t109 * t97 + t110 * t98) + t294; m(7) * (t37 * t7 + t42 * t69 + t43 * t68) + m(6) * (t109 * t95 + t110 * t94 + t35 * t90) + t293; m(7) * (t23 * t42 + t24 * t43 + t37 * t8) + m(6) * (t109 * t40 + t110 * t41 + t36 * t90) + t295; m(6) * (-t90 * t291 + (t109 * t292 + t110 * t289) * t288) + m(7) * (-t37 * t291 + (t289 * t43 + t292 * t42) * t288); m(7) * (t37 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2 + t90 ^ 2) + t295; m(7) * (t103 * t87 + t104 * t88) - t364 + t319; m(7) * (t103 * t69 + t104 * t68 + t7 * t84) + t312; m(7) * (t103 * t23 + t104 * t24 + t8 * t84) + t318; m(7) * (-t84 * t291 + (t103 * t292 + t104 * t289) * t288); m(7) * (t103 * t42 + t104 * t43 + t37 * t84) + t318; m(7) * (t103 ^ 2 + t104 ^ 2 + t84 ^ 2) + t318;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
