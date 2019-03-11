% Calculate joint inertia matrix for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:39:14
% EndTime: 2019-03-10 03:39:23
% DurationCPUTime: 5.26s
% Computational Cost: add. (23841->569), mult. (19960->778), div. (0->0), fcn. (21738->12), ass. (0->289)
t289 = qJ(4) + qJ(5);
t281 = qJ(6) + t289;
t273 = sin(t281);
t274 = cos(t281);
t293 = sin(qJ(1));
t290 = qJ(2) + qJ(3);
t280 = cos(t290);
t296 = cos(qJ(1));
t365 = t280 * t296;
t222 = -t273 * t365 + t274 * t293;
t223 = t273 * t293 + t274 * t365;
t278 = sin(t290);
t368 = t278 * t296;
t148 = t223 * rSges(7,1) + t222 * rSges(7,2) + rSges(7,3) * t368;
t294 = cos(qJ(4));
t275 = t294 * pkin(4) + pkin(3);
t279 = cos(t289);
t253 = pkin(5) * t279 + t275;
t277 = sin(t289);
t291 = sin(qJ(4));
t254 = pkin(4) * t291 + pkin(5) * t277;
t398 = t253 * t365 + t293 * t254 + t148;
t297 = -pkin(10) - pkin(9);
t286 = -pkin(11) + t297;
t343 = t286 - t297;
t364 = t291 * t293;
t348 = -pkin(4) * t364 - t275 * t365;
t397 = -t343 * t368 + t348 + t398;
t314 = Icges(4,5) * t280 - Icges(4,6) * t278;
t214 = -Icges(4,3) * t296 + t293 * t314;
t215 = Icges(4,3) * t293 + t296 * t314;
t366 = t280 * t293;
t220 = -t273 * t366 - t274 * t296;
t221 = -t273 * t296 + t274 * t366;
t369 = t278 * t293;
t140 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t369;
t142 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t369;
t144 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t369;
t48 = t140 * t369 + t142 * t220 + t144 * t221;
t141 = Icges(7,5) * t223 + Icges(7,6) * t222 + Icges(7,3) * t368;
t143 = Icges(7,4) * t223 + Icges(7,2) * t222 + Icges(7,6) * t368;
t145 = Icges(7,1) * t223 + Icges(7,4) * t222 + Icges(7,5) * t368;
t49 = t141 * t369 + t143 * t220 + t145 * t221;
t23 = t293 * t49 - t296 * t48;
t228 = -t277 * t366 - t279 * t296;
t229 = -t277 * t296 + t279 * t366;
t152 = Icges(6,5) * t229 + Icges(6,6) * t228 + Icges(6,3) * t369;
t154 = Icges(6,4) * t229 + Icges(6,2) * t228 + Icges(6,6) * t369;
t156 = Icges(6,1) * t229 + Icges(6,4) * t228 + Icges(6,5) * t369;
t52 = t152 * t369 + t154 * t228 + t156 * t229;
t230 = -t277 * t365 + t279 * t293;
t231 = t277 * t293 + t279 * t365;
t153 = Icges(6,5) * t231 + Icges(6,6) * t230 + Icges(6,3) * t368;
t155 = Icges(6,4) * t231 + Icges(6,2) * t230 + Icges(6,6) * t368;
t157 = Icges(6,1) * t231 + Icges(6,4) * t230 + Icges(6,5) * t368;
t53 = t153 * t369 + t155 * t228 + t157 * t229;
t28 = t293 * t53 - t296 * t52;
t288 = t296 ^ 2;
t374 = Icges(4,4) * t280;
t316 = -Icges(4,2) * t278 + t374;
t217 = Icges(4,6) * t293 + t296 * t316;
t375 = Icges(4,4) * t278;
t318 = Icges(4,1) * t280 - t375;
t219 = Icges(4,5) * t293 + t296 * t318;
t312 = -t217 * t278 + t219 * t280;
t216 = -Icges(4,6) * t296 + t293 * t316;
t218 = -Icges(4,5) * t296 + t293 * t318;
t313 = t216 * t278 - t218 * t280;
t361 = t294 * t296;
t241 = -t280 * t364 - t361;
t362 = t293 * t294;
t363 = t291 * t296;
t242 = t280 * t362 - t363;
t171 = Icges(5,5) * t242 + Icges(5,6) * t241 + Icges(5,3) * t369;
t173 = Icges(5,4) * t242 + Icges(5,2) * t241 + Icges(5,6) * t369;
t175 = Icges(5,1) * t242 + Icges(5,4) * t241 + Icges(5,5) * t369;
t72 = t171 * t369 + t173 * t241 + t175 * t242;
t243 = -t280 * t363 + t362;
t244 = t280 * t361 + t364;
t172 = Icges(5,5) * t244 + Icges(5,6) * t243 + Icges(5,3) * t368;
t174 = Icges(5,4) * t244 + Icges(5,2) * t243 + Icges(5,6) * t368;
t176 = Icges(5,1) * t244 + Icges(5,4) * t243 + Icges(5,5) * t368;
t73 = t172 * t369 + t174 * t241 + t176 * t242;
t37 = t293 * t73 - t296 * t72;
t396 = -t23 - t28 - t37 - t288 * t214 - (t312 * t293 + (-t215 + t313) * t296) * t293;
t287 = t293 ^ 2;
t192 = -Icges(7,3) * t280 + (Icges(7,5) * t274 - Icges(7,6) * t273) * t278;
t193 = -Icges(7,6) * t280 + (Icges(7,4) * t274 - Icges(7,2) * t273) * t278;
t194 = -Icges(7,5) * t280 + (Icges(7,1) * t274 - Icges(7,4) * t273) * t278;
t85 = t192 * t369 + t193 * t220 + t194 * t221;
t5 = -t85 * t280 + (t293 * t48 + t296 * t49) * t278;
t50 = t140 * t368 + t142 * t222 + t144 * t223;
t51 = t141 * t368 + t143 * t222 + t145 * t223;
t86 = t192 * t368 + t193 * t222 + t194 * t223;
t6 = -t86 * t280 + (t293 * t50 + t296 * t51) * t278;
t395 = t6 * t368 + t5 * t369;
t394 = -t280 / 0.2e1;
t393 = t293 / 0.2e1;
t392 = -t296 / 0.2e1;
t292 = sin(qJ(2));
t391 = pkin(2) * t292;
t390 = pkin(3) * t280;
t389 = pkin(9) * t278;
t388 = -pkin(3) + t275;
t198 = -Icges(6,5) * t280 + (Icges(6,1) * t279 - Icges(6,4) * t277) * t278;
t185 = t278 * t279 * t198;
t196 = -Icges(6,3) * t280 + (Icges(6,5) * t279 - Icges(6,6) * t277) * t278;
t197 = -Icges(6,6) * t280 + (Icges(6,4) * t279 - Icges(6,2) * t277) * t278;
t371 = t197 * t277;
t107 = -t280 * t196 - t278 * t371 + t185;
t183 = t278 * t274 * t194;
t372 = t193 * t273;
t102 = -t280 * t192 - t278 * t372 + t183;
t373 = t102 * t280;
t60 = -t280 * t140 + (-t142 * t273 + t144 * t274) * t278;
t61 = -t280 * t141 + (-t143 * t273 + t145 * t274) * t278;
t13 = -t373 + (t293 * t60 + t296 * t61) * t278;
t67 = -t280 * t152 + (-t154 * t277 + t156 * t279) * t278;
t68 = -t280 * t153 + (-t155 * t277 + t157 * t279) * t278;
t387 = -t13 + t107 * t280 - (t293 * t67 + t296 * t68) * t278;
t295 = cos(qJ(2));
t386 = rSges(3,1) * t295;
t385 = rSges(3,2) * t292;
t384 = t296 * rSges(3,3);
t383 = t60 * t296;
t382 = t61 * t293;
t381 = t67 * t296;
t380 = t68 * t293;
t79 = -t280 * t171 + (-t173 * t291 + t175 * t294) * t278;
t379 = t79 * t296;
t80 = -t280 * t172 + (-t174 * t291 + t176 * t294) * t278;
t378 = t80 * t293;
t377 = Icges(3,4) * t292;
t376 = Icges(3,4) * t295;
t207 = -Icges(5,6) * t280 + (Icges(5,4) * t294 - Icges(5,2) * t291) * t278;
t370 = t207 * t291;
t367 = t278 * t297;
t298 = -pkin(8) - pkin(7);
t360 = t296 * t298;
t359 = -t102 - t107;
t346 = pkin(4) * t363 + t293 * t367;
t347 = t253 - t275;
t129 = -t296 * t254 + (-t278 * t286 + t280 * t347) * t293 + t346;
t320 = -t221 * rSges(7,1) - t220 * rSges(7,2);
t147 = rSges(7,3) * t369 - t320;
t131 = t147 * t368;
t358 = t129 * t368 + t131;
t159 = t231 * rSges(6,1) + t230 * rSges(6,2) + rSges(6,3) * t368;
t306 = -t296 * t367 - t348;
t345 = pkin(3) * t365 + pkin(9) * t368;
t170 = t306 - t345;
t356 = -t159 - t170;
t169 = (t280 * t388 - t389) * t293 - t346;
t191 = (pkin(9) + t297) * t280 + t388 * t278;
t355 = t280 * t169 + t191 * t369;
t180 = t278 * t347 + t280 * t343;
t195 = -t280 * rSges(7,3) + (rSges(7,1) * t274 - rSges(7,2) * t273) * t278;
t354 = -t180 - t195;
t110 = t280 * t147 + t195 * t369;
t321 = -t229 * rSges(6,1) - t228 * rSges(6,2);
t158 = rSges(6,3) * t369 - t321;
t199 = -t280 * rSges(6,3) + (rSges(6,1) * t279 - rSges(6,2) * t277) * t278;
t115 = t280 * t158 + t199 * t369;
t353 = -t191 - t199;
t276 = pkin(2) * t295 + pkin(1);
t263 = t296 * t276;
t285 = t296 * pkin(7);
t352 = t293 * (t360 + t285 + (-pkin(1) + t276) * t293) + t296 * (-t296 * pkin(1) + t263 + (-pkin(7) - t298) * t293);
t307 = rSges(4,1) * t365 - rSges(4,2) * t368 + t293 * rSges(4,3);
t323 = rSges(4,1) * t280 - rSges(4,2) * t278;
t161 = t293 * (-t296 * rSges(4,3) + t293 * t323) + t296 * t307;
t211 = -t280 * rSges(5,3) + (rSges(5,1) * t294 - rSges(5,2) * t291) * t278;
t252 = t278 * pkin(3) - t280 * pkin(9);
t351 = -t211 - t252;
t350 = t287 * (t389 + t390) + t296 * t345;
t344 = t293 * rSges(3,3) + t296 * t386;
t342 = t287 + t288;
t94 = t196 * t369 + t197 * t228 + t198 * t229;
t11 = -t94 * t280 + (t293 * t52 + t296 * t53) * t278;
t54 = t152 * t368 + t154 * t230 + t156 * t231;
t55 = t153 * t368 + t155 * t230 + t157 * t231;
t95 = t196 * t368 + t197 * t230 + t198 * t231;
t12 = -t95 * t280 + (t293 * t54 + t296 * t55) * t278;
t341 = t11 * t369 + t12 * t368 + t395;
t340 = -t170 - t397;
t339 = -t191 + t354;
t338 = -t252 + t353;
t178 = t244 * rSges(5,1) + t243 * rSges(5,2) + rSges(5,3) * t368;
t337 = t369 / 0.2e1;
t336 = t368 / 0.2e1;
t251 = rSges(4,1) * t278 + rSges(4,2) * t280;
t335 = -t251 - t391;
t334 = -t252 - t391;
t24 = t293 * t51 - t296 * t50;
t29 = t293 * t55 - t296 * t54;
t74 = t171 * t368 + t173 * t243 + t175 * t244;
t75 = t172 * t368 + t174 * t243 + t176 * t244;
t38 = t293 * t75 - t296 * t74;
t333 = (t287 * t215 + t24 + t29 + t38 + (t313 * t296 + (-t214 + t312) * t293) * t296) * t293;
t332 = (t60 + t85) * t337 + (t61 + t86) * t336;
t331 = -t293 * t298 + t263;
t330 = -t280 * t13 + t395;
t46 = t280 * t129 + t180 * t369 + t110;
t329 = t293 * t169 + t296 * t170 + t350;
t328 = -t252 + t339;
t322 = -t242 * rSges(5,1) - t241 * rSges(5,2);
t177 = rSges(5,3) * t369 - t322;
t89 = t293 * t177 + t296 * t178 + t350;
t327 = t23 * t337 + t24 * t336 + t5 * t392 + t6 * t393 + (t382 - t383) * t394;
t326 = -t191 + t334;
t325 = -t211 + t334;
t324 = -t385 + t386;
t319 = Icges(3,1) * t295 - t377;
t317 = -Icges(3,2) * t292 + t376;
t315 = Icges(3,5) * t295 - Icges(3,6) * t292;
t248 = Icges(4,2) * t280 + t375;
t249 = Icges(4,1) * t278 + t374;
t309 = -t248 * t278 + t249 * t280;
t308 = -t199 + t326;
t44 = t293 * t158 + t296 * t159 + t329;
t305 = t326 + t354;
t304 = t280 * t387 + t341;
t303 = t332 + (t67 + t94) * t337 + (t68 + t95) * t336;
t30 = t329 + t397 * t296 + (t129 + t147) * t293;
t302 = t11 * t392 + t12 * t393 + t28 * t337 + t29 * t336 + t327 + (t380 - t381) * t394;
t301 = t296 * t396 + t333;
t206 = -Icges(5,3) * t280 + (Icges(5,5) * t294 - Icges(5,6) * t291) * t278;
t208 = -Icges(5,5) * t280 + (Icges(5,1) * t294 - Icges(5,4) * t291) * t278;
t103 = t206 * t369 + t207 * t241 + t208 * t242;
t17 = -t103 * t280 + (t293 * t72 + t296 * t73) * t278;
t104 = t206 * t368 + t207 * t243 + t208 * t244;
t18 = -t104 * t280 + (t293 * t74 + t296 * t75) * t278;
t300 = t17 * t392 + t18 * t393 + t38 * t336 + t37 * t337 + t302 + (t378 - t379) * t394;
t247 = Icges(4,5) * t278 + Icges(4,6) * t280;
t299 = -t383 / 0.2e1 + t382 / 0.2e1 - t381 / 0.2e1 + t380 / 0.2e1 - t379 / 0.2e1 + t378 / 0.2e1 + (t217 * t280 + t219 * t278 + t293 * t247 + t296 * t309 + t104 + t86 + t95) * t393 + (t216 * t280 + t218 * t278 - t296 * t247 + t293 * t309 + t103 + t85 + t94) * t392;
t262 = rSges(2,1) * t296 - rSges(2,2) * t293;
t261 = -rSges(2,1) * t293 - rSges(2,2) * t296;
t260 = rSges(3,1) * t292 + rSges(3,2) * t295;
t234 = Icges(3,3) * t293 + t296 * t315;
t233 = -Icges(3,3) * t296 + t293 * t315;
t213 = t335 * t296;
t212 = t335 * t293;
t201 = t293 * pkin(7) + (pkin(1) - t385) * t296 + t344;
t200 = t384 + t285 + (-pkin(1) - t324) * t293;
t190 = t278 * t294 * t208;
t188 = t307 + t331;
t187 = (rSges(4,3) - t298) * t296 + (-t276 - t323) * t293;
t182 = t351 * t296;
t181 = t351 * t293;
t179 = t296 * (-t296 * t385 + t344) + (t293 * t324 - t384) * t293;
t166 = t325 * t296;
t165 = t325 * t293;
t151 = t169 * t368;
t137 = t158 * t368;
t128 = t331 + t178 + t345;
t127 = -t360 + (-t390 - t276 + (-rSges(5,3) - pkin(9)) * t278) * t293 + t322;
t123 = t338 * t296;
t122 = t338 * t293;
t120 = -t178 * t280 - t211 * t368;
t119 = t177 * t280 + t211 * t369;
t118 = t308 * t296;
t117 = t308 * t293;
t116 = -t280 * t159 - t199 * t368;
t114 = -t280 * t206 - t278 * t370 + t190;
t113 = t306 + t331 + t159;
t112 = -t360 + (-rSges(6,3) * t278 - t275 * t280 - t276) * t293 + t321 + t346;
t111 = -t280 * t148 - t195 * t368;
t109 = t161 + t352;
t108 = (t177 * t296 - t178 * t293) * t278;
t106 = -t286 * t368 + t331 + t398;
t105 = (t254 - t298) * t296 + (-t253 * t280 - t276 + (-rSges(7,3) + t286) * t278) * t293 + t320;
t99 = -t159 * t369 + t137;
t98 = t328 * t296;
t97 = t328 * t293;
t96 = -t148 * t369 + t131;
t93 = t305 * t296;
t92 = t305 * t293;
t70 = t280 * t356 + t353 * t368;
t69 = t115 + t355;
t64 = t89 + t352;
t47 = -t280 * t397 + t354 * t368;
t45 = t356 * t369 + t137 + t151;
t43 = -t369 * t397 + t358;
t42 = t44 + t352;
t41 = t280 * t340 + t339 * t368;
t40 = t46 + t355;
t32 = t340 * t369 + t151 + t358;
t21 = t30 + t352;
t1 = [t295 * (Icges(3,2) * t295 + t377) + t292 * (Icges(3,1) * t292 + t376) + Icges(2,3) + t183 + t185 + t190 + (-t192 - t196 - t206 + t248) * t280 + (t249 - t370 - t371 - t372) * t278 + m(7) * (t105 ^ 2 + t106 ^ 2) + m(6) * (t112 ^ 2 + t113 ^ 2) + m(5) * (t127 ^ 2 + t128 ^ 2) + m(4) * (t187 ^ 2 + t188 ^ 2) + m(3) * (t200 ^ 2 + t201 ^ 2) + m(2) * (t261 ^ 2 + t262 ^ 2); m(7) * (t105 * t93 + t106 * t92) + m(6) * (t112 * t118 + t113 * t117) + m(5) * (t127 * t166 + t128 * t165) + m(4) * (t187 * t213 + t188 * t212) + t299 + (t288 / 0.2e1 + t287 / 0.2e1) * (Icges(3,5) * t292 + Icges(3,6) * t295) + m(3) * (-t200 * t296 - t201 * t293) * t260 + (t295 * (Icges(3,6) * t293 + t296 * t317) + t292 * (Icges(3,5) * t293 + t296 * t319)) * t393 + (t295 * (-Icges(3,6) * t296 + t293 * t317) + t292 * (-Icges(3,5) * t296 + t293 * t319)) * t392; m(7) * (t21 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(6) * (t117 ^ 2 + t118 ^ 2 + t42 ^ 2) + m(5) * (t165 ^ 2 + t166 ^ 2 + t64 ^ 2) + m(4) * (t109 ^ 2 + t212 ^ 2 + t213 ^ 2) + m(3) * (t260 ^ 2 * t342 + t179 ^ 2) + t293 * t287 * t234 + t333 + (-t288 * t233 + (-t293 * t233 + t296 * t234) * t293 + t396) * t296; m(7) * (t105 * t98 + t106 * t97) + m(6) * (t112 * t123 + t113 * t122) + m(5) * (t127 * t182 + t128 * t181) + m(4) * (-t187 * t296 - t188 * t293) * t251 + t299; m(7) * (t30 * t21 + t92 * t97 + t93 * t98) + m(6) * (t117 * t122 + t118 * t123 + t44 * t42) + m(5) * (t165 * t181 + t166 * t182 + t64 * t89) + m(4) * (t161 * t109 + (-t212 * t293 - t213 * t296) * t251) + t301; m(7) * (t30 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(6) * (t122 ^ 2 + t123 ^ 2 + t44 ^ 2) + m(5) * (t181 ^ 2 + t182 ^ 2 + t89 ^ 2) + m(4) * (t251 ^ 2 * t342 + t161 ^ 2) + t301; (-t114 + t359) * t280 + m(7) * (t105 * t40 + t106 * t41) + m(6) * (t112 * t69 + t113 * t70) + m(5) * (t119 * t127 + t120 * t128) + ((t104 / 0.2e1 + t80 / 0.2e1) * t296 + (t103 / 0.2e1 + t79 / 0.2e1) * t293) * t278 + t303; m(7) * (t32 * t21 + t40 * t93 + t41 * t92) + m(6) * (t117 * t70 + t118 * t69 + t45 * t42) + m(5) * (t108 * t64 + t119 * t166 + t120 * t165) + t300; m(7) * (t32 * t30 + t40 * t98 + t41 * t97) + m(6) * (t122 * t70 + t123 * t69 + t45 * t44) + m(5) * (t108 * t89 + t119 * t182 + t120 * t181) + t300; (t114 * t280 + t387) * t280 + (t293 * t17 + t296 * t18 - t280 * (t293 * t79 + t296 * t80)) * t278 + m(7) * (t32 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t45 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t108 ^ 2 + t119 ^ 2 + t120 ^ 2) + t341; t359 * t280 + m(7) * (t105 * t46 + t106 * t47) + m(6) * (t112 * t115 + t113 * t116) + t303; m(7) * (t43 * t21 + t46 * t93 + t47 * t92) + m(6) * (t115 * t118 + t116 * t117 + t42 * t99) + t302; m(7) * (t43 * t30 + t46 * t98 + t47 * t97) + m(6) * (t115 * t123 + t116 * t122 + t44 * t99) + t302; m(7) * (t32 * t43 + t40 * t46 + t41 * t47) + m(6) * (t115 * t69 + t116 * t70 + t45 * t99) + t304; m(7) * (t43 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t115 ^ 2 + t116 ^ 2 + t99 ^ 2) + t304; m(7) * (t105 * t110 + t106 * t111) - t373 + t332; m(7) * (t110 * t93 + t111 * t92 + t21 * t96) + t327; m(7) * (t110 * t98 + t111 * t97 + t30 * t96) + t327; m(7) * (t110 * t40 + t111 * t41 + t32 * t96) + t330; m(7) * (t110 * t46 + t111 * t47 + t43 * t96) + t330; m(7) * (t110 ^ 2 + t111 ^ 2 + t96 ^ 2) + t330;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
