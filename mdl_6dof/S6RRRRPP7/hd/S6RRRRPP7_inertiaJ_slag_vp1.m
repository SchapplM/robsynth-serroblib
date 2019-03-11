% Calculate joint inertia matrix for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:11
% EndTime: 2019-03-09 21:17:24
% DurationCPUTime: 6.88s
% Computational Cost: add. (27535->676), mult. (58220->916), div. (0->0), fcn. (74398->12), ass. (0->307)
t291 = cos(pkin(6));
t296 = sin(qJ(1));
t298 = cos(qJ(2));
t355 = t296 * t298;
t295 = sin(qJ(2));
t299 = cos(qJ(1));
t356 = t295 * t299;
t275 = t291 * t356 + t355;
t294 = sin(qJ(3));
t290 = sin(pkin(6));
t358 = t290 * t299;
t366 = cos(qJ(3));
t254 = t275 * t366 - t294 * t358;
t354 = t298 * t299;
t357 = t295 * t296;
t274 = -t291 * t354 + t357;
t340 = qJ(4) + pkin(11);
t288 = sin(t340);
t327 = cos(t340);
t209 = t254 * t288 - t274 * t327;
t210 = t254 * t327 + t274 * t288;
t368 = rSges(7,3) + qJ(6);
t369 = rSges(7,1) + pkin(5);
t370 = -t368 * t209 - t369 * t210;
t361 = t290 * t295;
t272 = -t291 * t366 + t294 * t361;
t329 = t290 * t366;
t273 = t291 * t294 + t295 * t329;
t359 = t290 * t298;
t222 = Icges(4,5) * t273 - Icges(4,6) * t272 - Icges(4,3) * t359;
t223 = Icges(4,4) * t273 - Icges(4,2) * t272 - Icges(4,6) * t359;
t224 = Icges(4,1) * t273 - Icges(4,4) * t272 - Icges(4,5) * t359;
t114 = -t222 * t359 - t272 * t223 + t273 * t224;
t241 = t273 * t288 + t327 * t359;
t242 = t273 * t327 - t288 * t359;
t169 = Icges(7,5) * t242 + Icges(7,6) * t272 + Icges(7,3) * t241;
t171 = Icges(7,4) * t242 + Icges(7,2) * t272 + Icges(7,6) * t241;
t173 = Icges(7,1) * t242 + Icges(7,4) * t272 + Icges(7,5) * t241;
t88 = t241 * t169 + t272 * t171 + t242 * t173;
t170 = Icges(6,5) * t242 - Icges(6,6) * t241 + Icges(6,3) * t272;
t172 = Icges(6,4) * t242 - Icges(6,2) * t241 + Icges(6,6) * t272;
t174 = Icges(6,1) * t242 - Icges(6,4) * t241 + Icges(6,5) * t272;
t89 = t272 * t170 - t241 * t172 + t242 * t174;
t293 = sin(qJ(4));
t297 = cos(qJ(4));
t251 = -t273 * t293 - t297 * t359;
t338 = t293 * t359;
t252 = t273 * t297 - t338;
t179 = Icges(5,5) * t252 + Icges(5,6) * t251 + Icges(5,3) * t272;
t180 = Icges(5,4) * t252 + Icges(5,2) * t251 + Icges(5,6) * t272;
t181 = Icges(5,1) * t252 + Icges(5,4) * t251 + Icges(5,5) * t272;
t93 = t272 * t179 + t251 * t180 + t252 * t181;
t367 = -t114 - t88 - t89 - t93;
t287 = pkin(4) * t297 + pkin(3);
t365 = -pkin(3) + t287;
t226 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t358;
t364 = t226 * t299;
t363 = t274 * t293;
t276 = t291 * t355 + t356;
t362 = t276 * t293;
t360 = t290 * t296;
t253 = t275 * t294 + t299 * t329;
t249 = t253 * pkin(10);
t292 = -qJ(5) - pkin(10);
t339 = pkin(4) * t363;
t121 = -t253 * t292 + t254 * t365 - t249 + t339;
t200 = pkin(3) * t254 + t249;
t191 = t276 * t200;
t353 = t276 * t121 + t191;
t277 = -t291 * t357 + t354;
t255 = t277 * t294 - t296 * t329;
t256 = t277 * t366 + t294 * t360;
t201 = t256 * pkin(3) + pkin(10) * t255;
t331 = pkin(4) * t362 - t255 * t292 + t256 * t287;
t122 = -t201 + t331;
t211 = t256 * t288 - t276 * t327;
t212 = t256 * t327 + t276 * t288;
t140 = t212 * rSges(6,1) - t211 * rSges(6,2) + t255 * rSges(6,3);
t352 = -t122 - t140;
t351 = rSges(7,2) * t253 - t370;
t350 = t255 * rSges(7,2) + t368 * t211 + t212 * t369;
t220 = -t256 * t293 + t276 * t297;
t221 = t256 * t297 + t362;
t150 = t221 * rSges(5,1) + t220 * rSges(5,2) + t255 * rSges(5,3);
t349 = -t150 - t201;
t348 = rSges(7,2) * t272 + t368 * t241 + t242 * t369;
t176 = rSges(6,1) * t242 - rSges(6,2) * t241 + rSges(6,3) * t272;
t177 = -pkin(4) * t338 + t365 * t273 + (-pkin(10) - t292) * t272;
t347 = -t176 - t177;
t182 = rSges(5,1) * t252 + rSges(5,2) * t251 + rSges(5,3) * t272;
t237 = pkin(3) * t273 + pkin(10) * t272;
t346 = -t182 - t237;
t345 = t200 * t359 + t274 * t237;
t239 = t277 * pkin(2) + pkin(9) * t276;
t236 = t291 * t239;
t344 = t291 * t201 + t236;
t238 = pkin(2) * t275 + t274 * pkin(9);
t343 = -t200 - t238;
t342 = t238 * t360 + t239 * t358;
t341 = t299 * pkin(1) + pkin(8) * t360;
t337 = t291 * t122 + t344;
t336 = -t121 + t343;
t335 = -t122 - t350;
t334 = -t201 + t352;
t333 = -t177 - t348;
t332 = -t237 + t347;
t190 = t256 * rSges(4,1) - t255 * rSges(4,2) + t276 * rSges(4,3);
t262 = Icges(3,3) * t291 + (Icges(3,5) * t295 + Icges(3,6) * t298) * t290;
t263 = Icges(3,6) * t291 + (Icges(3,4) * t295 + Icges(3,2) * t298) * t290;
t264 = Icges(3,5) * t291 + (Icges(3,1) * t295 + Icges(3,4) * t298) * t290;
t330 = t291 * t262 + t263 * t359 + t264 * t361;
t233 = t277 * rSges(3,1) - t276 * rSges(3,2) + rSges(3,3) * t360;
t328 = -t296 * pkin(1) + pkin(8) * t358;
t225 = rSges(4,1) * t273 - rSges(4,2) * t272 - rSges(4,3) * t359;
t278 = (pkin(2) * t295 - pkin(9) * t298) * t290;
t326 = t290 * (-t225 - t278);
t218 = -t254 * t293 + t274 * t297;
t219 = t254 * t297 + t363;
t143 = Icges(5,5) * t219 + Icges(5,6) * t218 + Icges(5,3) * t253;
t145 = Icges(5,4) * t219 + Icges(5,2) * t218 + Icges(5,6) * t253;
t147 = Icges(5,1) * t219 + Icges(5,4) * t218 + Icges(5,5) * t253;
t62 = t143 * t255 + t145 * t220 + t147 * t221;
t144 = Icges(5,5) * t221 + Icges(5,6) * t220 + Icges(5,3) * t255;
t146 = Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t255;
t148 = Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t255;
t63 = t144 * t255 + t146 * t220 + t148 * t221;
t81 = t179 * t255 + t180 * t220 + t181 * t221;
t14 = t253 * t62 + t255 * t63 + t272 * t81;
t123 = Icges(7,5) * t210 + Icges(7,6) * t253 + Icges(7,3) * t209;
t127 = Icges(7,4) * t210 + Icges(7,2) * t253 + Icges(7,6) * t209;
t131 = Icges(7,1) * t210 + Icges(7,4) * t253 + Icges(7,5) * t209;
t54 = t123 * t211 + t127 * t255 + t131 * t212;
t124 = Icges(7,5) * t212 + Icges(7,6) * t255 + Icges(7,3) * t211;
t128 = Icges(7,4) * t212 + Icges(7,2) * t255 + Icges(7,6) * t211;
t132 = Icges(7,1) * t212 + Icges(7,4) * t255 + Icges(7,5) * t211;
t55 = t124 * t211 + t128 * t255 + t132 * t212;
t74 = t169 * t211 + t171 * t255 + t173 * t212;
t3 = t253 * t54 + t255 * t55 + t272 * t74;
t125 = Icges(6,5) * t210 - Icges(6,6) * t209 + Icges(6,3) * t253;
t129 = Icges(6,4) * t210 - Icges(6,2) * t209 + Icges(6,6) * t253;
t133 = Icges(6,1) * t210 - Icges(6,4) * t209 + Icges(6,5) * t253;
t56 = t125 * t255 - t129 * t211 + t133 * t212;
t126 = Icges(6,5) * t212 - Icges(6,6) * t211 + Icges(6,3) * t255;
t130 = Icges(6,4) * t212 - Icges(6,2) * t211 + Icges(6,6) * t255;
t134 = Icges(6,1) * t212 - Icges(6,4) * t211 + Icges(6,5) * t255;
t57 = t126 * t255 - t130 * t211 + t134 * t212;
t75 = t170 * t255 - t172 * t211 + t174 * t212;
t4 = t253 * t56 + t255 * t57 + t272 * t75;
t325 = t14 / 0.2e1 + t3 / 0.2e1 + t4 / 0.2e1;
t60 = t143 * t253 + t145 * t218 + t147 * t219;
t61 = t144 * t253 + t146 * t218 + t148 * t219;
t80 = t179 * t253 + t180 * t218 + t181 * t219;
t15 = t274 * t60 + t276 * t61 - t359 * t80;
t50 = t123 * t209 + t127 * t253 + t131 * t210;
t51 = t124 * t209 + t128 * t253 + t132 * t210;
t72 = t169 * t209 + t171 * t253 + t173 * t210;
t5 = t274 * t50 + t276 * t51 - t359 * t72;
t52 = t125 * t253 - t129 * t209 + t133 * t210;
t53 = t126 * t253 - t130 * t209 + t134 * t210;
t73 = t170 * t253 - t172 * t209 + t174 * t210;
t6 = t274 * t52 + t276 * t53 - t359 * t73;
t324 = t15 / 0.2e1 + t6 / 0.2e1 + t5 / 0.2e1;
t1 = t253 * t50 + t255 * t51 + t272 * t72;
t13 = t253 * t60 + t255 * t61 + t272 * t80;
t2 = t253 * t52 + t255 * t53 + t272 * t73;
t323 = -t2 / 0.2e1 - t1 / 0.2e1 - t13 / 0.2e1;
t16 = t274 * t62 + t276 * t63 - t359 * t81;
t7 = t274 * t54 + t276 * t55 - t359 * t74;
t8 = t274 * t56 + t276 * t57 - t359 * t75;
t322 = t7 / 0.2e1 + t16 / 0.2e1 + t8 / 0.2e1;
t321 = t121 * t359 + t274 * t177 + t345;
t320 = -t201 + t335;
t319 = -t237 + t333;
t318 = t200 * t360 + t201 * t358 + t342;
t10 = t73 * t291 + (t296 * t53 - t299 * t52) * t290;
t17 = t80 * t291 + (t296 * t61 - t299 * t60) * t290;
t9 = t72 * t291 + (t296 * t51 - t299 * t50) * t290;
t317 = t17 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1;
t11 = t74 * t291 + (t296 * t55 - t299 * t54) * t290;
t12 = t75 * t291 + (t296 * t57 - t299 * t56) * t290;
t18 = t81 * t291 + (t296 * t63 - t299 * t62) * t290;
t316 = t18 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1;
t64 = t123 * t241 + t127 * t272 + t131 * t242;
t65 = t124 * t241 + t128 * t272 + t132 * t242;
t86 = t88 * t291;
t23 = t86 + (t65 * t296 - t64 * t299) * t290;
t66 = t125 * t272 - t129 * t241 + t133 * t242;
t67 = t126 * t272 - t130 * t241 + t134 * t242;
t87 = t89 * t291;
t24 = t87 + (t67 * t296 - t66 * t299) * t290;
t69 = t143 * t272 + t145 * t251 + t147 * t252;
t70 = t144 * t272 + t146 * t251 + t148 * t252;
t92 = t93 * t291;
t27 = t92 + (t70 * t296 - t69 * t299) * t290;
t315 = t23 / 0.2e1 + t27 / 0.2e1 + t24 / 0.2e1;
t82 = t88 * t272;
t19 = t64 * t253 + t65 * t255 + t82;
t83 = t89 * t272;
t20 = t66 * t253 + t67 * t255 + t83;
t91 = t93 * t272;
t25 = t69 * t253 + t70 * t255 + t91;
t314 = t25 / 0.2e1 + t20 / 0.2e1 + t19 / 0.2e1;
t21 = t64 * t274 + t65 * t276 - t359 * t88;
t22 = t66 * t274 + t67 * t276 - t359 * t89;
t26 = t69 * t274 + t70 * t276 - t359 * t93;
t313 = t26 / 0.2e1 + t22 / 0.2e1 + t21 / 0.2e1;
t312 = t290 * (-t278 + t346);
t311 = -rSges(6,1) * t210 + rSges(6,2) * t209;
t310 = t239 + t341;
t309 = t290 * (-t278 + t332);
t308 = t121 * t360 + t122 * t358 + t318;
t307 = t290 * (-t278 + t319);
t306 = -t238 + t328;
t189 = rSges(4,1) * t254 - rSges(4,2) * t253 + rSges(4,3) * t274;
t149 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t253;
t232 = t275 * rSges(3,1) - t274 * rSges(3,2) - rSges(3,3) * t358;
t305 = t310 + t331;
t304 = t69 / 0.2e1 + t66 / 0.2e1 + t64 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1 + t80 / 0.2e1;
t303 = t70 / 0.2e1 + t67 / 0.2e1 + t65 / 0.2e1 + t75 / 0.2e1 + t74 / 0.2e1 + t81 / 0.2e1;
t183 = Icges(4,5) * t254 - Icges(4,6) * t253 + Icges(4,3) * t274;
t185 = Icges(4,4) * t254 - Icges(4,2) * t253 + Icges(4,6) * t274;
t187 = Icges(4,1) * t254 - Icges(4,4) * t253 + Icges(4,5) * t274;
t100 = -t183 * t359 - t185 * t272 + t187 * t273;
t107 = t222 * t274 - t223 * t253 + t224 * t254;
t302 = t100 / 0.2e1 + t107 / 0.2e1 + t304;
t184 = Icges(4,5) * t256 - Icges(4,6) * t255 + Icges(4,3) * t276;
t186 = Icges(4,4) * t256 - Icges(4,2) * t255 + Icges(4,6) * t276;
t188 = Icges(4,1) * t256 - Icges(4,4) * t255 + Icges(4,5) * t276;
t101 = -t184 * t359 - t186 * t272 + t188 * t273;
t108 = t222 * t276 - t223 * t255 + t224 * t256;
t301 = t101 / 0.2e1 + t108 / 0.2e1 + t303;
t300 = -t254 * t287 + t306 - t339;
t280 = rSges(2,1) * t299 - t296 * rSges(2,2);
t279 = -t296 * rSges(2,1) - rSges(2,2) * t299;
t265 = rSges(3,3) * t291 + (rSges(3,1) * t295 + rSges(3,2) * t298) * t290;
t231 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t360;
t230 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t358;
t229 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t360;
t228 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t358;
t227 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t360;
t214 = t233 + t341;
t213 = -t232 + t328;
t194 = -t291 * t232 - t265 * t358;
t193 = t233 * t291 - t265 * t360;
t178 = t330 * t291;
t167 = (t232 * t296 + t233 * t299) * t290;
t166 = t262 * t360 - t263 * t276 + t264 * t277;
t165 = -t262 * t358 - t274 * t263 + t275 * t264;
t157 = t253 * t177;
t152 = t310 + t190;
t151 = -t189 + t306;
t142 = -t190 * t359 - t225 * t276;
t141 = t189 * t359 + t225 * t274;
t138 = rSges(6,3) * t253 - t311;
t136 = t227 * t291 + (t229 * t298 + t231 * t295) * t290;
t135 = t226 * t291 + (t228 * t298 + t230 * t295) * t290;
t115 = t272 * t122;
t113 = t114 * t291;
t112 = t255 * t121;
t111 = t189 * t276 - t190 * t274;
t110 = (-t189 - t238) * t291 + t299 * t326;
t109 = t190 * t291 + t296 * t326 + t236;
t106 = t310 - t349;
t105 = -t149 - t200 + t306;
t104 = (t189 * t296 + t190 * t299) * t290 + t342;
t103 = t150 * t272 - t182 * t255;
t102 = -t149 * t272 + t182 * t253;
t99 = t305 + t140;
t98 = (-rSges(6,3) + t292) * t253 + t300 + t311;
t97 = t184 * t276 - t186 * t255 + t188 * t256;
t96 = t183 * t276 - t185 * t255 + t187 * t256;
t95 = t184 * t274 - t186 * t253 + t188 * t254;
t94 = t183 * t274 - t185 * t253 + t187 * t254;
t90 = t149 * t255 - t150 * t253;
t85 = t276 * t346 + t349 * t359;
t84 = t149 * t359 + t182 * t274 + t345;
t79 = (-t149 + t343) * t291 + t299 * t312;
t78 = t150 * t291 + t296 * t312 + t344;
t77 = t305 + t350;
t76 = (-rSges(7,2) + t292) * t253 + t300 + t370;
t71 = t149 * t276 + t274 * t349 + t191;
t68 = (t149 * t296 + t150 * t299) * t290 + t318;
t59 = t140 * t272 + t255 * t347 + t115;
t58 = t176 * t253 + t157 + (-t121 - t138) * t272;
t49 = t276 * t332 + t334 * t359;
t48 = t138 * t359 + t176 * t274 + t321;
t47 = (-t138 + t336) * t291 + t299 * t309;
t46 = t140 * t291 + t296 * t309 + t337;
t45 = t138 * t255 + t253 * t352 + t112;
t44 = t113 + (-t100 * t299 + t101 * t296) * t290;
t43 = t100 * t274 + t101 * t276 - t114 * t359;
t42 = t255 * t333 + t272 * t350 + t115;
t41 = t157 + t348 * t253 + (-t121 - t351) * t272;
t40 = t138 * t276 + t274 * t334 + t353;
t39 = (t138 * t296 + t140 * t299) * t290 + t308;
t38 = t276 * t319 + t320 * t359;
t37 = t274 * t348 + t351 * t359 + t321;
t36 = t108 * t291 + (t296 * t97 - t299 * t96) * t290;
t35 = t107 * t291 + (t296 * t95 - t299 * t94) * t290;
t34 = (t336 - t351) * t291 + t299 * t307;
t33 = t291 * t350 + t296 * t307 + t337;
t32 = -t108 * t359 + t274 * t96 + t276 * t97;
t31 = -t107 * t359 + t274 * t94 + t276 * t95;
t30 = t253 * t335 + t255 * t351 + t112;
t29 = t274 * t320 + t276 * t351 + t353;
t28 = (t296 * t351 + t299 * t350) * t290 + t308;
t116 = [m(6) * (t98 ^ 2 + t99 ^ 2) + m(7) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t105 ^ 2 + t106 ^ 2) + m(4) * (t151 ^ 2 + t152 ^ 2) + m(3) * (t213 ^ 2 + t214 ^ 2) + m(2) * (t279 ^ 2 + t280 ^ 2) + Icges(2,3) + t330 - t367; t92 + t87 + t178 + t86 + t113 + m(7) * (t33 * t77 + t34 * t76) + m(6) * (t46 * t99 + t47 * t98) + m(5) * (t105 * t79 + t106 * t78) + m(4) * (t109 * t152 + t110 * t151) + m(3) * (t193 * t214 + t194 * t213) + ((-t135 / 0.2e1 - t165 / 0.2e1 - t302) * t299 + (t136 / 0.2e1 + t166 / 0.2e1 + t301) * t296) * t290; (t27 + t23 + t24 + t44 + t178) * t291 + m(6) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(7) * (t28 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(5) * (t68 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t104 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(3) * (t167 ^ 2 + t193 ^ 2 + t194 ^ 2) + ((-t10 - t9 - t17 - t35 + (-t274 * t228 + t275 * t230 - t290 * t364) * t358) * t299 + (t12 + t11 + t18 + t36 + ((-t229 * t276 + t231 * t277 + (t227 * t296 - t364) * t290) * t296 + (t227 * t358 + t228 * t276 + t274 * t229 - t230 * t277 - t275 * t231) * t299) * t290) * t296 + ((-t135 - t165) * t299 + (t136 + t166) * t296) * t291) * t290; t367 * t359 + m(6) * (t48 * t98 + t49 * t99) + m(7) * (t37 * t76 + t38 * t77) + m(5) * (t105 * t84 + t106 * t85) + m(4) * (t141 * t151 + t142 * t152) + t301 * t276 + t302 * t274; (t43 / 0.2e1 + t313) * t291 + (t36 / 0.2e1 + t316) * t276 + (t35 / 0.2e1 + t317) * t274 + m(7) * (t28 * t29 + t33 * t38 + t34 * t37) + m(6) * (t39 * t40 + t46 * t49 + t47 * t48) + m(5) * (t68 * t71 + t78 * t85 + t79 * t84) + m(4) * (t104 * t111 + t109 * t142 + t110 * t141) + ((-t31 / 0.2e1 - t324) * t299 + (-t44 / 0.2e1 - t315) * t298 + (t32 / 0.2e1 + t322) * t296) * t290; (-t21 - t22 - t26 - t43) * t359 + (t7 + t16 + t8 + t32) * t276 + (t15 + t6 + t5 + t31) * t274 + m(7) * (t29 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t71 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(4) * (t111 ^ 2 + t141 ^ 2 + t142 ^ 2); t91 + t82 + t83 + m(6) * (t58 * t98 + t59 * t99) + m(7) * (t41 * t76 + t42 * t77) + m(5) * (t102 * t105 + t103 * t106) + t303 * t255 + t304 * t253; t314 * t291 + t315 * t272 + t316 * t255 + t317 * t253 + m(6) * (t39 * t45 + t46 * t59 + t47 * t58) + m(7) * (t28 * t30 + t33 * t42 + t34 * t41) + m(5) * (t102 * t79 + t103 * t78 + t68 * t90) + (t296 * t325 + t299 * t323) * t290; -t314 * t359 + t325 * t276 - t323 * t274 + t313 * t272 + t322 * t255 + t324 * t253 + m(7) * (t29 * t30 + t37 * t41 + t38 * t42) + m(6) * (t40 * t45 + t48 * t58 + t49 * t59) + m(5) * (t102 * t84 + t103 * t85 + t71 * t90); (t25 + t20 + t19) * t272 + (t3 + t14 + t4) * t255 + (t13 + t2 + t1) * t253 + m(7) * (t30 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t45 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2 + t90 ^ 2); m(6) * (t253 * t99 + t255 * t98) + m(7) * (t253 * t77 + t255 * t76); m(6) * (t253 * t46 + t255 * t47 + t272 * t39) + m(7) * (t253 * t33 + t255 * t34 + t272 * t28); m(7) * (t253 * t38 + t255 * t37 + t272 * t29) + m(6) * (t253 * t49 + t255 * t48 + t272 * t40); m(7) * (t253 * t42 + t255 * t41 + t272 * t30) + m(6) * (t253 * t59 + t255 * t58 + t272 * t45); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t253 ^ 2 + t255 ^ 2 + t272 ^ 2); m(7) * (t209 * t77 + t211 * t76); m(7) * (t209 * t33 + t211 * t34 + t241 * t28); m(7) * (t209 * t38 + t211 * t37 + t241 * t29); m(7) * (t209 * t42 + t211 * t41 + t241 * t30); m(7) * (t209 * t253 + t211 * t255 + t241 * t272); m(7) * (t209 ^ 2 + t211 ^ 2 + t241 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t116(1) t116(2) t116(4) t116(7) t116(11) t116(16); t116(2) t116(3) t116(5) t116(8) t116(12) t116(17); t116(4) t116(5) t116(6) t116(9) t116(13) t116(18); t116(7) t116(8) t116(9) t116(10) t116(14) t116(19); t116(11) t116(12) t116(13) t116(14) t116(15) t116(20); t116(16) t116(17) t116(18) t116(19) t116(20) t116(21);];
Mq  = res;
