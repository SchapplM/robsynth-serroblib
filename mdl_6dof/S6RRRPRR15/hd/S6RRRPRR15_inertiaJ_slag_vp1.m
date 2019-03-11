% Calculate joint inertia matrix for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR15_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR15_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR15_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:39
% EndTime: 2019-03-09 20:24:58
% DurationCPUTime: 8.02s
% Computational Cost: add. (43747->698), mult. (121108->957), div. (0->0), fcn. (158660->14), ass. (0->322)
t299 = cos(pkin(6));
t303 = sin(qJ(1));
t305 = cos(qJ(2));
t370 = t303 * t305;
t302 = sin(qJ(2));
t306 = cos(qJ(1));
t372 = t302 * t306;
t284 = -t299 * t370 - t372;
t298 = sin(pkin(6));
t379 = cos(pkin(7));
t336 = t298 * t379;
t378 = sin(pkin(7));
t266 = -t284 * t378 + t303 * t336;
t369 = t305 * t306;
t371 = t303 * t302;
t282 = t299 * t369 - t371;
t265 = -t282 * t378 - t306 * t336;
t373 = t298 * t306;
t283 = t299 * t372 + t370;
t384 = cos(qJ(3));
t326 = t384 * t378;
t320 = t298 * t326;
t328 = t379 * t384;
t382 = sin(qJ(3));
t243 = -t282 * t328 + t283 * t382 + t306 * t320;
t301 = sin(qJ(5));
t383 = cos(qJ(5));
t208 = -t243 * t383 + t265 * t301;
t285 = -t299 * t371 + t369;
t245 = -t284 * t328 + t285 * t382 - t303 * t320;
t210 = -t245 * t383 + t266 * t301;
t374 = t298 * t305;
t376 = t298 * t302;
t262 = -t299 * t326 - t328 * t374 + t376 * t382;
t281 = t299 * t379 - t374 * t378;
t237 = -t262 * t383 + t281 * t301;
t209 = t243 * t301 + t265 * t383;
t325 = t378 * t382;
t319 = t298 * t325;
t327 = t379 * t382;
t244 = t282 * t327 + t283 * t384 - t306 * t319;
t300 = sin(qJ(6));
t304 = cos(qJ(6));
t152 = -t209 * t300 + t244 * t304;
t153 = t209 * t304 + t244 * t300;
t103 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t208;
t105 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t208;
t107 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t208;
t28 = t103 * t208 + t105 * t152 + t107 * t153;
t211 = t245 * t301 + t266 * t383;
t246 = t284 * t327 + t285 * t384 + t303 * t319;
t154 = -t211 * t300 + t246 * t304;
t155 = t211 * t304 + t246 * t300;
t104 = Icges(7,5) * t155 + Icges(7,6) * t154 + Icges(7,3) * t210;
t106 = Icges(7,4) * t155 + Icges(7,2) * t154 + Icges(7,6) * t210;
t108 = Icges(7,1) * t155 + Icges(7,4) * t154 + Icges(7,5) * t210;
t29 = t104 * t208 + t106 * t152 + t108 * t153;
t238 = t262 * t301 + t281 * t383;
t263 = t299 * t325 + (t302 * t384 + t305 * t327) * t298;
t199 = -t238 * t300 + t263 * t304;
t200 = t238 * t304 + t263 * t300;
t131 = Icges(7,5) * t200 + Icges(7,6) * t199 + Icges(7,3) * t237;
t132 = Icges(7,4) * t200 + Icges(7,2) * t199 + Icges(7,6) * t237;
t133 = Icges(7,1) * t200 + Icges(7,4) * t199 + Icges(7,5) * t237;
t47 = t131 * t208 + t132 * t152 + t133 * t153;
t1 = t208 * t28 + t210 * t29 + t237 * t47;
t390 = t1 / 0.2e1;
t30 = t103 * t210 + t105 * t154 + t107 * t155;
t31 = t104 * t210 + t106 * t154 + t108 * t155;
t48 = t131 * t210 + t132 * t154 + t133 * t155;
t2 = t208 * t30 + t210 * t31 + t237 * t48;
t389 = t2 / 0.2e1;
t38 = t103 * t237 + t105 * t199 + t107 * t200;
t39 = t104 * t237 + t106 * t199 + t108 * t200;
t57 = t237 * t131 + t199 * t132 + t200 * t133;
t52 = t57 * t237;
t9 = t38 * t208 + t39 * t210 + t52;
t388 = t9 / 0.2e1;
t387 = t208 / 0.2e1;
t386 = t210 / 0.2e1;
t385 = t237 / 0.2e1;
t381 = pkin(3) * t244;
t380 = pkin(5) * t209;
t251 = Icges(3,5) * t283 + Icges(3,6) * t282 - Icges(3,3) * t373;
t377 = t251 * t306;
t375 = t298 * t303;
t322 = -rSges(7,1) * t153 - rSges(7,2) * t152;
t109 = rSges(7,3) * t208 - t322;
t368 = pkin(12) * t208 + t109 + t380;
t110 = t155 * rSges(7,1) + t154 * rSges(7,2) + t210 * rSges(7,3);
t367 = t211 * pkin(5) + t210 * pkin(12) + t110;
t134 = rSges(7,1) * t200 + rSges(7,2) * t199 + rSges(7,3) * t237;
t366 = pkin(5) * t238 + pkin(12) * t237 + t134;
t324 = -rSges(5,1) * t265 - rSges(5,3) * t243;
t176 = -rSges(5,2) * t244 - t324;
t229 = t243 * qJ(4);
t190 = t229 + t381;
t365 = -t176 - t190;
t178 = t266 * t190;
t212 = pkin(4) * t265 + t244 * pkin(11);
t364 = t266 * t212 + t178;
t191 = t246 * pkin(3) + t245 * qJ(4);
t180 = t281 * t191;
t213 = t266 * pkin(4) + t246 * pkin(11);
t363 = t281 * t213 + t180;
t250 = t285 * pkin(2) + pkin(10) * t266;
t247 = t299 * t250;
t362 = t299 * t191 + t247;
t361 = -t190 - t212;
t360 = -t191 - t213;
t226 = pkin(3) * t263 + qJ(4) * t262;
t198 = t265 * t226;
t248 = pkin(4) * t281 + pkin(11) * t263;
t359 = t265 * t248 + t198;
t221 = rSges(5,1) * t281 - rSges(5,2) * t263 + rSges(5,3) * t262;
t358 = -t221 - t226;
t357 = -t226 - t248;
t249 = pkin(2) * t283 + pkin(10) * t265;
t356 = t249 * t375 + t250 * t373;
t355 = t306 * pkin(1) + pkin(9) * t375;
t135 = Icges(6,5) * t209 - Icges(6,6) * t208 + Icges(6,3) * t244;
t137 = Icges(6,4) * t209 - Icges(6,2) * t208 + Icges(6,6) * t244;
t139 = Icges(6,1) * t209 - Icges(6,4) * t208 + Icges(6,5) * t244;
t59 = t135 * t244 - t137 * t208 + t139 * t209;
t136 = Icges(6,5) * t211 - Icges(6,6) * t210 + Icges(6,3) * t246;
t138 = Icges(6,4) * t211 - Icges(6,2) * t210 + Icges(6,6) * t246;
t140 = Icges(6,1) * t211 - Icges(6,4) * t210 + Icges(6,5) * t246;
t60 = t136 * t244 - t138 * t208 + t140 * t209;
t156 = Icges(6,5) * t238 - Icges(6,6) * t237 + Icges(6,3) * t263;
t157 = Icges(6,4) * t238 - Icges(6,2) * t237 + Icges(6,6) * t263;
t158 = Icges(6,1) * t238 - Icges(6,4) * t237 + Icges(6,5) * t263;
t73 = t156 * t244 - t157 * t208 + t158 * t209;
t13 = t244 * t59 + t246 * t60 + t263 * t73;
t3 = t244 * t28 + t246 * t29 + t263 * t47;
t354 = t3 / 0.2e1 + t13 / 0.2e1;
t61 = t135 * t246 - t137 * t210 + t139 * t211;
t62 = t136 * t246 - t138 * t210 + t140 * t211;
t74 = t156 * t246 - t157 * t210 + t158 * t211;
t14 = t244 * t61 + t246 * t62 + t263 * t74;
t4 = t244 * t30 + t246 * t31 + t263 * t48;
t353 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t265 * t59 + t266 * t60 + t281 * t73;
t5 = t265 * t28 + t266 * t29 + t281 * t47;
t352 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t265 * t61 + t266 * t62 + t281 * t74;
t6 = t265 * t30 + t266 * t31 + t281 * t48;
t351 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t73 * t299 + (t303 * t60 - t306 * t59) * t298;
t7 = t47 * t299 + (-t28 * t306 + t29 * t303) * t298;
t350 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t74 * t299 + (t303 * t62 - t306 * t61) * t298;
t8 = t48 * t299 + (-t30 * t306 + t303 * t31) * t298;
t349 = t8 / 0.2e1 + t18 / 0.2e1;
t54 = t57 * t263;
t10 = t38 * t244 + t39 * t246 + t54;
t63 = t135 * t263 - t137 * t237 + t139 * t238;
t64 = t136 * t263 - t138 * t237 + t140 * t238;
t81 = t263 * t156 - t237 * t157 + t238 * t158;
t77 = t81 * t263;
t19 = t63 * t244 + t64 * t246 + t77;
t348 = t10 / 0.2e1 + t19 / 0.2e1;
t55 = t57 * t281;
t11 = t38 * t265 + t39 * t266 + t55;
t78 = t81 * t281;
t20 = t63 * t265 + t64 * t266 + t78;
t347 = t11 / 0.2e1 + t20 / 0.2e1;
t56 = t57 * t299;
t12 = t56 + (t39 * t303 - t38 * t306) * t298;
t80 = t81 * t299;
t21 = t80 + (t64 * t303 - t63 * t306) * t298;
t346 = t12 / 0.2e1 + t21 / 0.2e1;
t345 = t39 / 0.2e1 + t48 / 0.2e1;
t344 = t47 / 0.2e1 + t38 / 0.2e1;
t323 = -rSges(6,1) * t209 + rSges(6,2) * t208;
t141 = rSges(6,3) * t244 - t323;
t343 = -t141 + t361;
t161 = rSges(6,1) * t238 - rSges(6,2) * t237 + rSges(6,3) * t263;
t342 = -t161 + t357;
t215 = Icges(4,5) * t263 - Icges(4,6) * t262 + Icges(4,3) * t281;
t217 = Icges(4,4) * t263 - Icges(4,2) * t262 + Icges(4,6) * t281;
t219 = Icges(4,1) * t263 - Icges(4,4) * t262 + Icges(4,5) * t281;
t122 = t281 * t215 - t262 * t217 + t263 * t219;
t341 = t299 * t213 + t362;
t214 = Icges(5,5) * t281 - Icges(5,6) * t263 + Icges(5,3) * t262;
t216 = Icges(5,4) * t281 - Icges(5,2) * t263 + Icges(5,6) * t262;
t218 = Icges(5,1) * t281 - Icges(5,4) * t263 + Icges(5,5) * t262;
t123 = t262 * t214 - t263 * t216 + t281 * t218;
t142 = t211 * rSges(6,1) - t210 * rSges(6,2) + t246 * rSges(6,3);
t175 = t246 * rSges(4,1) - t245 * rSges(4,2) + t266 * rSges(4,3);
t177 = t266 * rSges(5,1) - t246 * rSges(5,2) + t245 * rSges(5,3);
t273 = Icges(3,3) * t299 + (Icges(3,5) * t302 + Icges(3,6) * t305) * t298;
t274 = Icges(3,6) * t299 + (Icges(3,4) * t302 + Icges(3,2) * t305) * t298;
t275 = Icges(3,5) * t299 + (Icges(3,1) * t302 + Icges(3,4) * t305) * t298;
t340 = t299 * t273 + t274 * t374 + t275 * t376;
t258 = t285 * rSges(3,1) + t284 * rSges(3,2) + rSges(3,3) * t375;
t339 = -t303 * pkin(1) + pkin(9) * t373;
t220 = rSges(4,1) * t263 - rSges(4,2) * t262 + rSges(4,3) * t281;
t269 = pkin(2) * t376 + pkin(10) * t281;
t335 = t298 * (-t220 - t269);
t334 = t361 - t368;
t333 = t357 - t366;
t332 = t190 * t375 + t191 * t373 + t356;
t329 = t298 * (-t269 + t358);
t321 = t298 * (-t269 + t342);
t318 = t64 / 0.2e1 + t74 / 0.2e1 + t345;
t317 = t73 / 0.2e1 + t63 / 0.2e1 + t344;
t316 = t212 * t375 + t213 * t373 + t332;
t315 = t298 * (-t269 + t333);
t174 = rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t265;
t314 = -t249 + t339;
t257 = t283 * rSges(3,1) + t282 * rSges(3,2) - rSges(3,3) * t373;
t313 = -t229 + t314;
t312 = t250 + t355;
t112 = t215 * t266 - t217 * t245 + t219 * t246;
t114 = t214 * t245 - t216 * t246 + t218 * t266;
t165 = Icges(4,5) * t246 - Icges(4,6) * t245 + Icges(4,3) * t266;
t169 = Icges(4,4) * t246 - Icges(4,2) * t245 + Icges(4,6) * t266;
t173 = Icges(4,1) * t246 - Icges(4,4) * t245 + Icges(4,5) * t266;
t91 = t165 * t281 - t169 * t262 + t173 * t263;
t163 = Icges(5,5) * t266 - Icges(5,6) * t246 + Icges(5,3) * t245;
t167 = Icges(5,4) * t266 - Icges(5,2) * t246 + Icges(5,6) * t245;
t171 = Icges(5,1) * t266 - Icges(5,4) * t246 + Icges(5,5) * t245;
t93 = t163 * t262 - t167 * t263 + t171 * t281;
t311 = t114 / 0.2e1 + t93 / 0.2e1 + t91 / 0.2e1 + t112 / 0.2e1 + t318;
t111 = t215 * t265 - t217 * t243 + t219 * t244;
t113 = t214 * t243 - t216 * t244 + t218 * t265;
t164 = Icges(4,5) * t244 - Icges(4,6) * t243 + Icges(4,3) * t265;
t168 = Icges(4,4) * t244 - Icges(4,2) * t243 + Icges(4,6) * t265;
t172 = Icges(4,1) * t244 - Icges(4,4) * t243 + Icges(4,5) * t265;
t90 = t164 * t281 - t168 * t262 + t172 * t263;
t162 = Icges(5,5) * t265 - Icges(5,6) * t244 + Icges(5,3) * t243;
t166 = Icges(5,4) * t265 - Icges(5,2) * t244 + Icges(5,6) * t243;
t170 = Icges(5,1) * t265 - Icges(5,4) * t244 + Icges(5,5) * t243;
t92 = t162 * t262 - t166 * t263 + t170 * t281;
t310 = t111 / 0.2e1 + t113 / 0.2e1 + t90 / 0.2e1 + t92 / 0.2e1 + t317;
t309 = -t212 + t313;
t308 = t191 + t312;
t307 = t213 + t308;
t292 = rSges(2,1) * t306 - t303 * rSges(2,2);
t291 = -t303 * rSges(2,1) - rSges(2,2) * t306;
t276 = t299 * rSges(3,3) + (rSges(3,1) * t302 + rSges(3,2) * t305) * t298;
t256 = Icges(3,1) * t285 + Icges(3,4) * t284 + Icges(3,5) * t375;
t255 = Icges(3,1) * t283 + Icges(3,4) * t282 - Icges(3,5) * t373;
t254 = Icges(3,4) * t285 + Icges(3,2) * t284 + Icges(3,6) * t375;
t253 = Icges(3,4) * t283 + Icges(3,2) * t282 - Icges(3,6) * t373;
t252 = Icges(3,5) * t285 + Icges(3,6) * t284 + Icges(3,3) * t375;
t242 = t258 + t355;
t241 = -t257 + t339;
t225 = -t299 * t257 - t276 * t373;
t224 = t258 * t299 - t276 * t375;
t222 = t340 * t299;
t197 = (t257 * t303 + t258 * t306) * t298;
t196 = t273 * t375 + t274 * t284 + t275 * t285;
t195 = -t273 * t373 + t282 * t274 + t283 * t275;
t160 = t299 * t252 + (t254 * t305 + t256 * t302) * t298;
t159 = t299 * t251 + (t253 * t305 + t255 * t302) * t298;
t146 = t312 + t175;
t145 = -t174 + t314;
t130 = t175 * t281 - t220 * t266;
t129 = -t174 * t281 + t220 * t265;
t127 = t308 + t177;
t126 = (rSges(5,2) - pkin(3)) * t244 + t313 + t324;
t125 = (-t174 - t249) * t299 + t306 * t335;
t124 = t299 * t175 + t303 * t335 + t247;
t121 = t123 * t299;
t120 = t122 * t299;
t117 = t174 * t266 - t175 * t265;
t116 = t123 * t281;
t115 = t122 * t281;
t102 = (t174 * t303 + t175 * t306) * t298 + t356;
t101 = t142 * t263 - t161 * t246;
t100 = -t141 * t263 + t161 * t244;
t99 = t307 + t142;
t98 = (-rSges(6,3) - pkin(3)) * t244 + t309 + t323;
t97 = t281 * t177 + t266 * t358 + t180;
t96 = t265 * t221 + t281 * t365 + t198;
t95 = (-t249 + t365) * t299 + t306 * t329;
t94 = t299 * t177 + t303 * t329 + t362;
t89 = t163 * t245 - t167 * t246 + t171 * t266;
t88 = t162 * t245 - t166 * t246 + t170 * t266;
t87 = t163 * t243 - t167 * t244 + t171 * t265;
t86 = t162 * t243 - t166 * t244 + t170 * t265;
t85 = t165 * t266 - t169 * t245 + t173 * t246;
t84 = t164 * t266 - t168 * t245 + t172 * t246;
t83 = t165 * t265 - t169 * t243 + t173 * t244;
t82 = t164 * t265 - t168 * t243 + t172 * t244;
t79 = t141 * t246 - t142 * t244;
t76 = t266 * t176 + t178 + (-t177 - t191) * t265;
t75 = (t176 * t303 + t177 * t306) * t298 + t332;
t72 = t110 * t237 - t134 * t210;
t71 = -t109 * t237 + t134 * t208;
t70 = (-t249 + t343) * t299 + t306 * t321;
t69 = t299 * t142 + t303 * t321 + t341;
t68 = t307 + t367;
t67 = -t381 - t380 + (-rSges(7,3) - pkin(12)) * t208 + t309 + t322;
t66 = t281 * t142 + t266 * t342 + t363;
t65 = t265 * t161 + t281 * t343 + t359;
t58 = t109 * t210 - t110 * t208;
t53 = (t141 * t303 + t142 * t306) * t298 + t316;
t51 = t266 * t141 + (-t142 + t360) * t265 + t364;
t50 = -t246 * t366 + t263 * t367;
t49 = t244 * t366 - t263 * t368;
t46 = t121 + (t93 * t303 - t92 * t306) * t298;
t45 = t120 + (t91 * t303 - t90 * t306) * t298;
t44 = -t244 * t367 + t246 * t368;
t43 = (-t249 + t334) * t299 + t306 * t315;
t42 = t299 * t367 + t303 * t315 + t341;
t41 = t92 * t265 + t93 * t266 + t116;
t40 = t90 * t265 + t91 * t266 + t115;
t37 = t266 * t333 + t281 * t367 + t363;
t36 = t265 * t366 + t281 * t334 + t359;
t35 = t114 * t299 + (t303 * t89 - t306 * t88) * t298;
t34 = t113 * t299 + (t303 * t87 - t306 * t86) * t298;
t33 = t112 * t299 + (t303 * t85 - t306 * t84) * t298;
t32 = t111 * t299 + (t303 * t83 - t306 * t82) * t298;
t27 = t114 * t281 + t265 * t88 + t266 * t89;
t26 = t113 * t281 + t265 * t86 + t266 * t87;
t25 = t112 * t281 + t265 * t84 + t266 * t85;
t24 = t111 * t281 + t265 * t82 + t266 * t83;
t23 = (t303 * t368 + t306 * t367) * t298 + t316;
t22 = t368 * t266 + (t360 - t367) * t265 + t364;
t118 = [m(7) * (t67 ^ 2 + t68 ^ 2) + m(6) * (t98 ^ 2 + t99 ^ 2) + m(5) * (t126 ^ 2 + t127 ^ 2) + m(4) * (t145 ^ 2 + t146 ^ 2) + m(3) * (t241 ^ 2 + t242 ^ 2) + m(2) * (t291 ^ 2 + t292 ^ 2) + t340 + Icges(2,3) + t123 + t122 + t81 + t57; t222 + t56 + t120 + t121 + t80 + m(3) * (t224 * t242 + t225 * t241) + m(7) * (t42 * t68 + t43 * t67) + m(6) * (t69 * t99 + t70 * t98) + m(5) * (t126 * t95 + t127 * t94) + m(4) * (t124 * t146 + t125 * t145) + ((-t159 / 0.2e1 - t195 / 0.2e1 - t310) * t306 + (t160 / 0.2e1 + t196 / 0.2e1 + t311) * t303) * t298; (t12 + t21 + t45 + t46 + t222) * t299 + m(7) * (t23 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t53 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t75 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(4) * (t102 ^ 2 + t124 ^ 2 + t125 ^ 2) + m(3) * (t197 ^ 2 + t224 ^ 2 + t225 ^ 2) + ((-t7 - t17 - t34 - t32 + (t282 * t253 + t283 * t255 - t298 * t377) * t373) * t306 + (t8 + t18 + t35 + t33 + ((t254 * t284 + t256 * t285 + (t252 * t303 - t377) * t298) * t303 + (t252 * t373 - t253 * t284 - t282 * t254 - t255 * t285 - t283 * t256) * t306) * t298) * t303 + ((-t159 - t195) * t306 + (t160 + t196) * t303) * t299) * t298; t55 + t78 + t116 + t115 + m(7) * (t36 * t67 + t37 * t68) + m(6) * (t65 * t98 + t66 * t99) + m(5) * (t126 * t96 + t127 * t97) + m(4) * (t129 * t145 + t130 * t146) + t311 * t266 + t310 * t265; (t40 / 0.2e1 + t41 / 0.2e1 + t347) * t299 + (t46 / 0.2e1 + t45 / 0.2e1 + t346) * t281 + (t35 / 0.2e1 + t33 / 0.2e1 + t349) * t266 + (t32 / 0.2e1 + t34 / 0.2e1 + t350) * t265 + m(7) * (t22 * t23 + t36 * t43 + t37 * t42) + m(6) * (t51 * t53 + t65 * t70 + t66 * t69) + m(5) * (t75 * t76 + t94 * t97 + t95 * t96) + m(4) * (t102 * t117 + t124 * t130 + t125 * t129) + ((-t24 / 0.2e1 - t26 / 0.2e1 - t352) * t306 + (t25 / 0.2e1 + t27 / 0.2e1 + t351) * t303) * t298; (t11 + t20 + t40 + t41) * t281 + (t6 + t16 + t25 + t27) * t266 + (t5 + t15 + t24 + t26) * t265 + m(7) * (t22 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t51 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t76 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(4) * (t117 ^ 2 + t129 ^ 2 + t130 ^ 2); m(7) * (t243 * t68 + t245 * t67) + m(6) * (t243 * t99 + t245 * t98) + m(5) * (t126 * t245 + t127 * t243); m(7) * (t23 * t262 + t243 * t42 + t245 * t43) + m(6) * (t243 * t69 + t245 * t70 + t262 * t53) + m(5) * (t243 * t94 + t245 * t95 + t262 * t75); m(7) * (t22 * t262 + t243 * t37 + t245 * t36) + m(6) * (t243 * t66 + t245 * t65 + t262 * t51) + m(5) * (t243 * t97 + t245 * t96 + t262 * t76); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t243 ^ 2 + t245 ^ 2 + t262 ^ 2); t54 + t77 + m(7) * (t49 * t67 + t50 * t68) + m(6) * (t100 * t98 + t101 * t99) + t318 * t246 + t317 * t244; t348 * t299 + t346 * t263 + t349 * t246 + t350 * t244 + m(7) * (t23 * t44 + t42 * t50 + t43 * t49) + m(6) * (t100 * t70 + t101 * t69 + t53 * t79) + (t303 * t353 - t306 * t354) * t298; t348 * t281 + t353 * t266 + t354 * t265 + t347 * t263 + t351 * t246 + t352 * t244 + m(7) * (t22 * t44 + t36 * t49 + t37 * t50) + m(6) * (t100 * t65 + t101 * t66 + t51 * t79); m(6) * (t100 * t245 + t101 * t243 + t262 * t79) + m(7) * (t243 * t50 + t245 * t49 + t262 * t44); (t10 + t19) * t263 + (t4 + t14) * t246 + (t3 + t13) * t244 + m(7) * (t44 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t100 ^ 2 + t101 ^ 2 + t79 ^ 2); t52 + m(7) * (t67 * t71 + t68 * t72) + t345 * t210 + t344 * t208; t8 * t386 + t12 * t385 + m(7) * (t58 * t23 + t42 * t72 + t43 * t71) + t299 * t388 + t7 * t387 + (t303 * t389 - t306 * t1 / 0.2e1) * t298; m(7) * (t58 * t22 + t36 * t71 + t37 * t72) + t266 * t389 + t5 * t387 + t11 * t385 + t265 * t390 + t6 * t386 + t281 * t388; m(7) * (t243 * t72 + t245 * t71 + t262 * t58); t244 * t390 + t4 * t386 + m(7) * (t58 * t44 + t49 * t71 + t50 * t72) + t3 * t387 + t263 * t388 + t10 * t385 + t246 * t389; t210 * t2 + t208 * t1 + t237 * t9 + m(7) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t118(1) t118(2) t118(4) t118(7) t118(11) t118(16); t118(2) t118(3) t118(5) t118(8) t118(12) t118(17); t118(4) t118(5) t118(6) t118(9) t118(13) t118(18); t118(7) t118(8) t118(9) t118(10) t118(14) t118(19); t118(11) t118(12) t118(13) t118(14) t118(15) t118(20); t118(16) t118(17) t118(18) t118(19) t118(20) t118(21);];
Mq  = res;
