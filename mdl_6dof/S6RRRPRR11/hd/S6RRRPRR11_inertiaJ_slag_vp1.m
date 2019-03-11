% Calculate joint inertia matrix for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:55
% EndTime: 2019-03-09 19:25:10
% DurationCPUTime: 7.54s
% Computational Cost: add. (27960->658), mult. (72150->905), div. (0->0), fcn. (94051->12), ass. (0->307)
t285 = cos(pkin(6));
t290 = sin(qJ(1));
t292 = cos(qJ(2));
t349 = t290 * t292;
t289 = sin(qJ(2));
t293 = cos(qJ(1));
t350 = t289 * t293;
t271 = t285 * t350 + t349;
t288 = sin(qJ(3));
t284 = sin(pkin(6));
t362 = cos(qJ(3));
t316 = t284 * t362;
t247 = t271 * t288 + t293 * t316;
t352 = t284 * t293;
t248 = t271 * t362 - t288 * t352;
t287 = sin(qJ(5));
t361 = cos(qJ(5));
t191 = -t247 * t361 + t248 * t287;
t348 = t292 * t293;
t351 = t289 * t290;
t273 = -t285 * t351 + t348;
t249 = t273 * t288 - t290 * t316;
t354 = t284 * t290;
t250 = t273 * t362 + t288 * t354;
t193 = -t249 * t361 + t250 * t287;
t355 = t284 * t289;
t268 = -t285 * t362 + t288 * t355;
t269 = t285 * t288 + t289 * t316;
t232 = -t268 * t361 + t269 * t287;
t192 = t247 * t287 + t248 * t361;
t270 = -t285 * t348 + t351;
t286 = sin(qJ(6));
t291 = cos(qJ(6));
t157 = -t192 * t286 - t270 * t291;
t158 = t192 * t291 - t270 * t286;
t101 = Icges(7,5) * t158 + Icges(7,6) * t157 + Icges(7,3) * t191;
t103 = Icges(7,4) * t158 + Icges(7,2) * t157 + Icges(7,6) * t191;
t105 = Icges(7,1) * t158 + Icges(7,4) * t157 + Icges(7,5) * t191;
t24 = t101 * t191 + t103 * t157 + t105 * t158;
t194 = t249 * t287 + t250 * t361;
t272 = t285 * t349 + t350;
t159 = -t194 * t286 - t272 * t291;
t160 = t194 * t291 - t272 * t286;
t102 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t193;
t104 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t193;
t106 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t193;
t25 = t102 * t191 + t104 * t157 + t106 * t158;
t233 = t268 * t287 + t269 * t361;
t353 = t284 * t292;
t208 = -t233 * t286 + t291 * t353;
t209 = t233 * t291 + t286 * t353;
t133 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t232;
t134 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t232;
t135 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t232;
t47 = t133 * t191 + t134 * t157 + t135 * t158;
t1 = t191 * t24 + t193 * t25 + t232 * t47;
t373 = -t1 / 0.2e1;
t28 = t101 * t232 + t103 * t208 + t105 * t209;
t29 = t102 * t232 + t104 * t208 + t106 * t209;
t61 = t232 * t133 + t208 * t134 + t209 * t135;
t52 = t61 * t232;
t9 = t28 * t191 + t29 * t193 + t52;
t372 = t9 / 0.2e1;
t125 = Icges(6,5) * t192 - Icges(6,6) * t191 - Icges(6,3) * t270;
t127 = Icges(6,4) * t192 - Icges(6,2) * t191 - Icges(6,6) * t270;
t129 = Icges(6,1) * t192 - Icges(6,4) * t191 - Icges(6,5) * t270;
t57 = -t125 * t270 - t127 * t191 + t129 * t192;
t126 = Icges(6,5) * t194 - Icges(6,6) * t193 - Icges(6,3) * t272;
t128 = Icges(6,4) * t194 - Icges(6,2) * t193 - Icges(6,6) * t272;
t130 = Icges(6,1) * t194 - Icges(6,4) * t193 - Icges(6,5) * t272;
t58 = -t126 * t270 - t128 * t191 + t130 * t192;
t151 = Icges(6,5) * t233 - Icges(6,6) * t232 + Icges(6,3) * t353;
t152 = Icges(6,4) * t233 - Icges(6,2) * t232 + Icges(6,6) * t353;
t153 = Icges(6,1) * t233 - Icges(6,4) * t232 + Icges(6,5) * t353;
t72 = -t151 * t270 - t152 * t191 + t153 * t192;
t14 = -t270 * t57 - t272 * t58 + t353 * t72;
t4 = -t24 * t270 - t25 * t272 + t353 * t47;
t371 = -t4 - t14;
t59 = -t125 * t272 - t127 * t193 + t129 * t194;
t60 = -t126 * t272 - t128 * t193 + t130 * t194;
t73 = -t151 * t272 - t152 * t193 + t153 * t194;
t16 = -t270 * t59 - t272 * t60 + t353 * t73;
t26 = t101 * t193 + t103 * t159 + t105 * t160;
t27 = t102 * t193 + t104 * t159 + t106 * t160;
t48 = t133 * t193 + t134 * t159 + t135 * t160;
t6 = -t26 * t270 - t27 * t272 + t353 * t48;
t370 = -t6 - t16;
t55 = t61 * t353;
t10 = t28 * t270 + t29 * t272 - t55;
t62 = t125 * t353 - t127 * t232 + t129 * t233;
t63 = t126 * t353 - t128 * t232 + t130 * t233;
t77 = t151 * t353 - t232 * t152 + t233 * t153;
t75 = t77 * t353;
t19 = t62 * t270 + t63 * t272 - t75;
t369 = -t10 - t19;
t330 = -t4 / 0.2e1 - t14 / 0.2e1;
t329 = t6 / 0.2e1 + t16 / 0.2e1;
t326 = t10 / 0.2e1 + t19 / 0.2e1;
t363 = t232 / 0.2e1;
t364 = t193 / 0.2e1;
t365 = t191 / 0.2e1;
t2 = t191 * t26 + t193 * t27 + t232 * t48;
t366 = t2 / 0.2e1;
t368 = -t10 * t363 + t270 * t373 - t272 * t366 + t353 * t372 + t6 * t364 + t4 * t365;
t308 = -rSges(7,1) * t158 - rSges(7,2) * t157;
t107 = rSges(7,3) * t191 - t308;
t360 = pkin(5) * t192;
t347 = pkin(11) * t191 + t107 + t360;
t367 = t272 * t347;
t136 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t232;
t345 = pkin(5) * t233 + pkin(11) * t232 + t136;
t49 = -t270 * t345 - t347 * t353;
t359 = pkin(9) * t272;
t358 = t55 + t75;
t131 = rSges(6,1) * t192 - rSges(6,2) * t191 - rSges(6,3) * t270;
t357 = t131 * t272;
t221 = Icges(3,5) * t271 - Icges(3,6) * t270 - Icges(3,3) * t352;
t356 = t221 * t293;
t108 = t160 * rSges(7,1) + t159 * rSges(7,2) + t193 * rSges(7,3);
t346 = t194 * pkin(5) + pkin(11) * t193 + t108;
t177 = t250 * rSges(5,1) + t272 * rSges(5,2) + t249 * rSges(5,3);
t197 = t250 * pkin(3) + qJ(4) * t249;
t344 = -t177 - t197;
t240 = t247 * qJ(4);
t196 = pkin(3) * t248 + t240;
t179 = t272 * t196;
t265 = t270 * pkin(10);
t211 = pkin(4) * t248 - t265;
t343 = t272 * t211 + t179;
t235 = pkin(3) * t269 + qJ(4) * t268;
t342 = t196 * t353 + t270 * t235;
t267 = t273 * pkin(2);
t237 = t267 + t359;
t234 = t285 * t237;
t341 = t285 * t197 + t234;
t236 = pkin(2) * t271 + t270 * pkin(9);
t340 = -t196 - t236;
t245 = t250 * pkin(4);
t212 = -pkin(10) * t272 + t245;
t339 = -t197 - t212;
t213 = Icges(5,5) * t269 - Icges(5,6) * t353 + Icges(5,3) * t268;
t217 = Icges(5,1) * t269 - Icges(5,4) * t353 + Icges(5,5) * t268;
t338 = t268 * t213 + t269 * t217;
t216 = Icges(4,4) * t269 - Icges(4,2) * t268 - Icges(4,6) * t353;
t218 = Icges(4,1) * t269 - Icges(4,4) * t268 - Icges(4,5) * t353;
t337 = -t268 * t216 + t269 * t218;
t219 = rSges(5,1) * t269 - rSges(5,2) * t353 + rSges(5,3) * t268;
t336 = -t219 - t235;
t335 = t236 * t354 + t237 * t352;
t253 = pkin(4) * t269 + pkin(10) * t353;
t334 = -t235 - t253;
t333 = t293 * pkin(1) + pkin(8) * t354;
t17 = t72 * t285 + (t290 * t58 - t293 * t57) * t284;
t7 = t47 * t285 + (-t24 * t293 + t25 * t290) * t284;
t328 = -t7 / 0.2e1 - t17 / 0.2e1;
t18 = t73 * t285 + (t290 * t60 - t293 * t59) * t284;
t8 = t48 * t285 + (-t26 * t293 + t27 * t290) * t284;
t327 = -t8 / 0.2e1 - t18 / 0.2e1;
t56 = t61 * t285;
t12 = t56 + (-t28 * t293 + t29 * t290) * t284;
t76 = t77 * t285;
t21 = t76 + (t63 * t290 - t62 * t293) * t284;
t324 = t12 / 0.2e1 + t21 / 0.2e1;
t323 = t47 / 0.2e1 + t28 / 0.2e1;
t322 = t48 / 0.2e1 + t29 / 0.2e1;
t132 = t194 * rSges(6,1) - t193 * rSges(6,2) - t272 * rSges(6,3);
t321 = -t132 + t339;
t154 = rSges(6,1) * t233 - rSges(6,2) * t232 + rSges(6,3) * t353;
t320 = -t154 + t334;
t319 = t285 * t212 + t341;
t318 = -t211 + t340;
t178 = t250 * rSges(4,1) - t249 * rSges(4,2) + t272 * rSges(4,3);
t255 = Icges(3,3) * t285 + (Icges(3,5) * t289 + Icges(3,6) * t292) * t284;
t256 = Icges(3,6) * t285 + (Icges(3,4) * t289 + Icges(3,2) * t292) * t284;
t257 = Icges(3,5) * t285 + (Icges(3,1) * t289 + Icges(3,4) * t292) * t284;
t317 = t285 * t255 + t256 * t353 + t257 * t355;
t228 = t273 * rSges(3,1) - t272 * rSges(3,2) + rSges(3,3) * t354;
t315 = -t290 * pkin(1) + pkin(8) * t352;
t220 = rSges(4,1) * t269 - rSges(4,2) * t268 - rSges(4,3) * t353;
t274 = (pkin(2) * t289 - pkin(9) * t292) * t284;
t314 = t284 * (-t220 - t274);
t313 = t339 - t346;
t312 = t334 - t345;
t311 = t196 * t354 + t197 * t352 + t335;
t310 = t211 * t353 + t270 * t253 + t342;
t309 = t284 * (-t274 + t336);
t307 = -rSges(5,2) * t270 - rSges(5,3) * t247;
t306 = t284 * (-t274 + t320);
t305 = -t62 / 0.2e1 - t72 / 0.2e1 - t323;
t304 = -t63 / 0.2e1 - t73 / 0.2e1 - t322;
t303 = t211 * t354 + t212 * t352 + t311;
t98 = -t131 * t353 - t154 * t270;
t302 = t284 * (-t274 + t312);
t301 = -t236 + t315;
t300 = t197 + t267 + t333;
t176 = rSges(4,1) * t248 - rSges(4,2) * t247 + rSges(4,3) * t270;
t299 = -t240 + t301;
t227 = t271 * rSges(3,1) - t270 * rSges(3,2) - rSges(3,3) * t352;
t297 = t245 + (pkin(9) - pkin(10)) * t272 + t300;
t215 = Icges(5,4) * t269 - Icges(5,2) * t353 + Icges(5,6) * t268;
t109 = t213 * t247 + t215 * t270 + t217 * t248;
t214 = Icges(4,5) * t269 - Icges(4,6) * t268 - Icges(4,3) * t353;
t110 = t214 * t270 - t216 * t247 + t218 * t248;
t163 = Icges(5,5) * t248 + Icges(5,6) * t270 + Icges(5,3) * t247;
t167 = Icges(5,4) * t248 + Icges(5,2) * t270 + Icges(5,6) * t247;
t171 = Icges(5,1) * t248 + Icges(5,4) * t270 + Icges(5,5) * t247;
t92 = t163 * t268 - t167 * t353 + t171 * t269;
t165 = Icges(4,5) * t248 - Icges(4,6) * t247 + Icges(4,3) * t270;
t169 = Icges(4,4) * t248 - Icges(4,2) * t247 + Icges(4,6) * t270;
t173 = Icges(4,1) * t248 - Icges(4,4) * t247 + Icges(4,5) * t270;
t94 = -t165 * t353 - t169 * t268 + t173 * t269;
t296 = t94 / 0.2e1 + t92 / 0.2e1 + t110 / 0.2e1 + t109 / 0.2e1 - t305;
t111 = t213 * t249 + t215 * t272 + t217 * t250;
t112 = t214 * t272 - t216 * t249 + t218 * t250;
t164 = Icges(5,5) * t250 + Icges(5,6) * t272 + Icges(5,3) * t249;
t168 = Icges(5,4) * t250 + Icges(5,2) * t272 + Icges(5,6) * t249;
t172 = Icges(5,1) * t250 + Icges(5,4) * t272 + Icges(5,5) * t249;
t93 = t164 * t268 - t168 * t353 + t172 * t269;
t166 = Icges(4,5) * t250 - Icges(4,6) * t249 + Icges(4,3) * t272;
t170 = Icges(4,4) * t250 - Icges(4,2) * t249 + Icges(4,6) * t272;
t174 = Icges(4,1) * t250 - Icges(4,4) * t249 + Icges(4,5) * t272;
t95 = -t166 * t353 - t170 * t268 + t174 * t269;
t295 = t95 / 0.2e1 + t93 / 0.2e1 + t112 / 0.2e1 + t111 / 0.2e1 - t304;
t294 = t265 + (-pkin(3) - pkin(4)) * t248 + t299;
t276 = rSges(2,1) * t293 - t290 * rSges(2,2);
t275 = -t290 * rSges(2,1) - rSges(2,2) * t293;
t258 = rSges(3,3) * t285 + (rSges(3,1) * t289 + rSges(3,2) * t292) * t284;
t226 = Icges(3,1) * t273 - Icges(3,4) * t272 + Icges(3,5) * t354;
t225 = Icges(3,1) * t271 - Icges(3,4) * t270 - Icges(3,5) * t352;
t224 = Icges(3,4) * t273 - Icges(3,2) * t272 + Icges(3,6) * t354;
t223 = Icges(3,4) * t271 - Icges(3,2) * t270 - Icges(3,6) * t352;
t222 = Icges(3,5) * t273 - Icges(3,6) * t272 + Icges(3,3) * t354;
t206 = t228 + t333;
t205 = -t227 + t315;
t182 = -t285 * t227 - t258 * t352;
t181 = t228 * t285 - t258 * t354;
t175 = rSges(5,1) * t248 - t307;
t161 = t317 * t285;
t149 = (t227 * t290 + t228 * t293) * t284;
t148 = t255 * t354 - t256 * t272 + t257 * t273;
t147 = -t255 * t352 - t270 * t256 + t271 * t257;
t146 = t237 + t178 + t333;
t145 = -t176 + t301;
t140 = -t178 * t353 - t220 * t272;
t139 = t176 * t353 + t220 * t270;
t138 = t222 * t285 + (t224 * t292 + t226 * t289) * t284;
t137 = t221 * t285 + (t223 * t292 + t225 * t289) * t284;
t124 = -t214 * t353 + t337;
t123 = -t215 * t353 + t338;
t122 = t124 * t285;
t121 = t123 * t285;
t120 = t177 + t300 + t359;
t119 = (-rSges(5,1) - pkin(3)) * t248 + t299 + t307;
t117 = t176 * t272 - t178 * t270;
t116 = (-t176 - t236) * t285 + t293 * t314;
t115 = t178 * t285 + t290 * t314 + t234;
t100 = (t176 * t290 + t178 * t293) * t284 + t335;
t99 = t132 * t353 + t154 * t272;
t97 = t272 * t336 + t344 * t353;
t96 = t175 * t353 + t219 * t270 + t342;
t91 = (-t175 + t340) * t285 + t293 * t309;
t90 = t177 * t285 + t290 * t309 + t341;
t89 = t297 + t132;
t88 = -t131 + t294;
t87 = t166 * t272 - t170 * t249 + t174 * t250;
t86 = t165 * t272 - t169 * t249 + t173 * t250;
t85 = t164 * t249 + t168 * t272 + t172 * t250;
t84 = t163 * t249 + t167 * t272 + t171 * t250;
t83 = t166 * t270 - t170 * t247 + t174 * t248;
t82 = t165 * t270 - t169 * t247 + t173 * t248;
t81 = t164 * t247 + t168 * t270 + t172 * t248;
t80 = t163 * t247 + t167 * t270 + t171 * t248;
t79 = t132 * t270 - t357;
t78 = t175 * t272 + t270 * t344 + t179;
t74 = (t175 * t290 + t177 * t293) * t284 + t311;
t71 = t272 * t320 + t321 * t353;
t70 = -t98 + t310;
t69 = t108 * t232 - t136 * t193;
t68 = -t107 * t232 + t136 * t191;
t67 = (-t131 + t318) * t285 + t293 * t306;
t66 = t132 * t285 + t290 * t306 + t319;
t65 = t297 + t346;
t64 = -t360 + (-rSges(7,3) - pkin(11)) * t191 + t294 + t308;
t54 = t270 * t321 + t343 + t357;
t53 = t107 * t193 - t108 * t191;
t51 = (t131 * t290 + t132 * t293) * t284 + t303;
t50 = t272 * t345 + t346 * t353;
t46 = t122 + (t95 * t290 - t94 * t293) * t284;
t45 = t121 + (t93 * t290 - t92 * t293) * t284;
t44 = -t124 * t353 + t94 * t270 + t95 * t272;
t43 = -t123 * t353 + t92 * t270 + t93 * t272;
t42 = t270 * t346 - t367;
t41 = t112 * t285 + (t290 * t87 - t293 * t86) * t284;
t40 = t111 * t285 + (t290 * t85 - t293 * t84) * t284;
t39 = t110 * t285 + (t290 * t83 - t293 * t82) * t284;
t38 = t109 * t285 + (t290 * t81 - t293 * t80) * t284;
t37 = t272 * t312 + t313 * t353;
t36 = t310 - t49;
t35 = -t112 * t353 + t270 * t86 + t272 * t87;
t34 = -t111 * t353 + t270 * t84 + t272 * t85;
t33 = -t110 * t353 + t270 * t82 + t272 * t83;
t32 = -t109 * t353 + t270 * t80 + t272 * t81;
t31 = (t318 - t347) * t285 + t293 * t302;
t30 = t285 * t346 + t290 * t302 + t319;
t23 = t270 * t313 + t343 + t367;
t22 = (t290 * t347 + t293 * t346) * t284 + t303;
t3 = [(-t214 - t215) * t353 + t317 + m(7) * (t64 ^ 2 + t65 ^ 2) + m(6) * (t88 ^ 2 + t89 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2) + m(4) * (t145 ^ 2 + t146 ^ 2) + m(3) * (t205 ^ 2 + t206 ^ 2) + m(2) * (t275 ^ 2 + t276 ^ 2) + Icges(2,3) + t77 + t61 + t337 + t338; t122 + t161 + t121 + t56 + t76 + m(6) * (t66 * t89 + t67 * t88) + m(7) * (t30 * t65 + t31 * t64) + m(5) * (t119 * t91 + t120 * t90) + m(4) * (t115 * t146 + t116 * t145) + m(3) * (t181 * t206 + t182 * t205) + ((-t137 / 0.2e1 - t147 / 0.2e1 - t296) * t293 + (t138 / 0.2e1 + t148 / 0.2e1 + t295) * t290) * t284; (t12 + t21 + t46 + t45 + t161) * t285 + m(7) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t51 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t74 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(4) * (t100 ^ 2 + t115 ^ 2 + t116 ^ 2) + m(3) * (t149 ^ 2 + t181 ^ 2 + t182 ^ 2) + ((-t7 - t17 - t39 - t38 + (-t270 * t223 + t271 * t225 - t284 * t356) * t352) * t293 + (t8 + t18 + t41 + t40 + ((-t224 * t272 + t226 * t273 + (t222 * t290 - t356) * t284) * t290 + (t222 * t352 + t223 * t272 + t270 * t224 - t225 * t273 - t271 * t226) * t293) * t284) * t290 + ((-t137 - t147) * t293 + (t138 + t148) * t290) * t285) * t284; (-t123 - t124) * t353 + m(7) * (t36 * t64 + t37 * t65) + m(6) * (t70 * t88 + t71 * t89) + m(5) * (t119 * t96 + t120 * t97) + m(4) * (t139 * t145 + t140 * t146) + t295 * t272 + t296 * t270 - t358; (t44 / 0.2e1 + t43 / 0.2e1 + t326) * t285 + (t41 / 0.2e1 + t40 / 0.2e1 - t327) * t272 + (t38 / 0.2e1 + t39 / 0.2e1 - t328) * t270 + m(5) * (t74 * t78 + t90 * t97 + t91 * t96) + m(6) * (t51 * t54 + t66 * t71 + t67 * t70) + m(7) * (t22 * t23 + t30 * t37 + t31 * t36) + m(4) * (t100 * t117 + t115 * t140 + t116 * t139) + ((-t32 / 0.2e1 - t33 / 0.2e1 - t330) * t293 + (-t45 / 0.2e1 - t46 / 0.2e1 - t324) * t292 + (t34 / 0.2e1 + t35 / 0.2e1 - t329) * t290) * t284; (-t43 - t44 + t369) * t353 + (t35 + t34 + t370) * t272 + (t33 + t32 + t371) * t270 + m(7) * (t23 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t54 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t78 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(4) * (t117 ^ 2 + t139 ^ 2 + t140 ^ 2); m(7) * (t247 * t65 + t249 * t64) + m(6) * (t247 * t89 + t249 * t88) + m(5) * (t119 * t249 + t120 * t247); m(7) * (t22 * t268 + t247 * t30 + t249 * t31) + m(6) * (t247 * t66 + t249 * t67 + t268 * t51) + m(5) * (t247 * t90 + t249 * t91 + t268 * t74); m(7) * (t23 * t268 + t247 * t37 + t249 * t36) + m(6) * (t247 * t71 + t249 * t70 + t268 * t54) + m(5) * (t247 * t97 + t249 * t96 + t268 * t78); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t247 ^ 2 + t249 ^ 2 + t268 ^ 2); m(7) * (t49 * t64 + t50 * t65) + m(6) * (t88 * t98 + t89 * t99) + t304 * t272 + t305 * t270 + t358; -t326 * t285 + t327 * t272 + t328 * t270 + m(7) * (t22 * t42 + t30 * t50 + t31 * t49) + m(6) * (t51 * t79 + t66 * t99 + t67 * t98) + (t290 * t329 + t292 * t324 + t293 * t330) * t284; m(7) * (t23 * t42 + t36 * t49 + t37 * t50) + m(6) * (t54 * t79 + t70 * t98 + t71 * t99) + 0.2e1 * t326 * t353 + 0.2e1 * t329 * t272 - 0.2e1 * t330 * t270; m(6) * (t247 * t99 + t249 * t98 + t268 * t79) + m(7) * (t247 * t50 + t249 * t49 + t268 * t42); t369 * t353 + t370 * t272 + t371 * t270 + m(7) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t79 ^ 2 + t98 ^ 2 + t99 ^ 2); m(7) * (t64 * t68 + t65 * t69) + t52 + t322 * t193 + t323 * t191; t12 * t363 + t285 * t372 + t7 * t365 + t8 * t364 + m(7) * (t22 * t53 + t30 * t69 + t31 * t68) + (t290 * t366 + t293 * t373) * t284; m(7) * (t23 * t53 + t36 * t68 + t37 * t69) - t368; m(7) * (t247 * t69 + t249 * t68 + t268 * t53); m(7) * (t42 * t53 + t49 * t68 + t50 * t69) + t368; t193 * t2 + t191 * t1 + t232 * t9 + m(7) * (t53 ^ 2 + t68 ^ 2 + t69 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
