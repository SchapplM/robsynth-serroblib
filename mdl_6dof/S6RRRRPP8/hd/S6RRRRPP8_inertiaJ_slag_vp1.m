% Calculate joint inertia matrix for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:22
% EndTime: 2019-03-09 21:30:39
% DurationCPUTime: 6.85s
% Computational Cost: add. (22023->664), mult. (56584->901), div. (0->0), fcn. (72628->10), ass. (0->295)
t359 = rSges(7,1) + pkin(5);
t285 = cos(pkin(6));
t289 = sin(qJ(1));
t290 = cos(qJ(2));
t344 = t289 * t290;
t288 = sin(qJ(2));
t291 = cos(qJ(1));
t345 = t288 * t291;
t271 = t285 * t345 + t344;
t287 = sin(qJ(3));
t284 = sin(pkin(6));
t347 = t284 * t291;
t355 = cos(qJ(3));
t249 = t271 * t355 - t287 * t347;
t343 = t290 * t291;
t346 = t288 * t289;
t270 = -t285 * t343 + t346;
t286 = sin(qJ(4));
t354 = cos(qJ(4));
t215 = t249 * t286 - t270 * t354;
t322 = t284 * t355;
t247 = t271 * t287 + t291 * t322;
t358 = rSges(7,3) + qJ(6);
t360 = -rSges(7,2) * t215 + t358 * t247;
t273 = -t285 * t346 + t343;
t349 = t284 * t289;
t252 = t273 * t355 + t287 * t349;
t272 = t285 * t344 + t345;
t217 = t252 * t286 - t272 * t354;
t218 = t252 * t354 + t272 * t286;
t250 = t273 * t287 - t289 * t322;
t357 = t217 * rSges(7,2) - t250 * rSges(7,3) + t359 * t218;
t350 = t284 * t288;
t267 = -t285 * t355 + t287 * t350;
t269 = t285 * t287 + t288 * t322;
t348 = t284 * t290;
t220 = Icges(4,5) * t269 - Icges(4,6) * t267 - Icges(4,3) * t348;
t221 = Icges(4,4) * t269 - Icges(4,2) * t267 - Icges(4,6) * t348;
t222 = Icges(4,1) * t269 - Icges(4,4) * t267 - Icges(4,5) * t348;
t113 = -t220 * t348 - t267 * t221 + t269 * t222;
t245 = t269 * t286 + t348 * t354;
t246 = t269 * t354 - t286 * t348;
t170 = Icges(7,5) * t246 + Icges(7,6) * t245 - Icges(7,3) * t267;
t173 = Icges(7,4) * t246 + Icges(7,2) * t245 - Icges(7,6) * t267;
t176 = Icges(7,1) * t246 + Icges(7,4) * t245 - Icges(7,5) * t267;
t93 = -t267 * t170 + t245 * t173 + t246 * t176;
t171 = Icges(6,5) * t246 + Icges(6,6) * t267 + Icges(6,3) * t245;
t174 = Icges(6,4) * t246 + Icges(6,2) * t267 + Icges(6,6) * t245;
t177 = Icges(6,1) * t246 + Icges(6,4) * t267 + Icges(6,5) * t245;
t94 = t245 * t171 + t267 * t174 + t246 * t177;
t172 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t267;
t175 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t267;
t178 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t267;
t95 = t267 * t172 - t245 * t175 + t246 * t178;
t356 = -t113 - t93 - t94 - t95;
t353 = pkin(10) * t250;
t224 = Icges(3,5) * t271 - Icges(3,6) * t270 - Icges(3,3) * t347;
t351 = t224 * t291;
t216 = t249 * t354 + t270 * t286;
t342 = t359 * t216 - t360;
t341 = -qJ(6) * t250 + t357;
t140 = t218 * rSges(6,1) + t250 * rSges(6,2) + t217 * rSges(6,3);
t159 = t218 * pkin(4) + qJ(5) * t217;
t340 = -t140 - t159;
t141 = t218 * rSges(5,1) - t217 * rSges(5,2) + t250 * rSges(5,3);
t244 = t252 * pkin(3);
t200 = t244 + t353;
t339 = -t141 - t200;
t205 = t215 * qJ(5);
t158 = pkin(4) * t216 + t205;
t199 = pkin(3) * t249 + t247 * pkin(10);
t191 = t272 * t199;
t338 = t272 * t158 + t191;
t337 = rSges(7,2) * t245 + t359 * t246 - t358 * t267;
t180 = rSges(6,1) * t246 + rSges(6,2) * t267 + rSges(6,3) * t245;
t198 = pkin(4) * t246 + qJ(5) * t245;
t336 = -t180 - t198;
t181 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t267;
t235 = pkin(3) * t269 + pkin(10) * t267;
t335 = -t181 - t235;
t334 = t199 * t348 + t270 * t235;
t237 = t273 * pkin(2) + pkin(9) * t272;
t234 = t285 * t237;
t333 = t285 * t200 + t234;
t236 = pkin(2) * t271 + t270 * pkin(9);
t332 = -t199 - t236;
t331 = t236 * t349 + t237 * t347;
t330 = t291 * pkin(1) + pkin(8) * t349;
t329 = -t159 - t341;
t328 = -t200 + t340;
t327 = t285 * t159 + t333;
t326 = -t158 + t332;
t325 = -t198 - t337;
t324 = -t235 + t336;
t189 = t252 * rSges(4,1) - t250 * rSges(4,2) + t272 * rSges(4,3);
t256 = Icges(3,3) * t285 + (Icges(3,5) * t288 + Icges(3,6) * t290) * t284;
t257 = Icges(3,6) * t285 + (Icges(3,4) * t288 + Icges(3,2) * t290) * t284;
t258 = Icges(3,5) * t285 + (Icges(3,1) * t288 + Icges(3,4) * t290) * t284;
t323 = t285 * t256 + t257 * t348 + t258 * t350;
t231 = t273 * rSges(3,1) - t272 * rSges(3,2) + rSges(3,3) * t349;
t321 = -t289 * pkin(1) + pkin(8) * t347;
t223 = rSges(4,1) * t269 - rSges(4,2) * t267 - rSges(4,3) * t348;
t274 = (pkin(2) * t288 - pkin(9) * t290) * t284;
t320 = t284 * (-t223 - t274);
t118 = Icges(7,5) * t216 + Icges(7,6) * t215 - Icges(7,3) * t247;
t124 = Icges(7,4) * t216 + Icges(7,2) * t215 - Icges(7,6) * t247;
t130 = Icges(7,1) * t216 + Icges(7,4) * t215 - Icges(7,5) * t247;
t50 = -t118 * t247 + t124 * t215 + t130 * t216;
t119 = Icges(7,5) * t218 + Icges(7,6) * t217 - Icges(7,3) * t250;
t125 = Icges(7,4) * t218 + Icges(7,2) * t217 - Icges(7,6) * t250;
t131 = Icges(7,1) * t218 + Icges(7,4) * t217 - Icges(7,5) * t250;
t51 = -t119 * t247 + t125 * t215 + t131 * t216;
t76 = -t170 * t247 + t173 * t215 + t176 * t216;
t1 = t247 * t50 + t250 * t51 + t267 * t76;
t120 = Icges(6,5) * t216 + Icges(6,6) * t247 + Icges(6,3) * t215;
t126 = Icges(6,4) * t216 + Icges(6,2) * t247 + Icges(6,6) * t215;
t132 = Icges(6,1) * t216 + Icges(6,4) * t247 + Icges(6,5) * t215;
t52 = t120 * t215 + t126 * t247 + t132 * t216;
t121 = Icges(6,5) * t218 + Icges(6,6) * t250 + Icges(6,3) * t217;
t127 = Icges(6,4) * t218 + Icges(6,2) * t250 + Icges(6,6) * t217;
t133 = Icges(6,1) * t218 + Icges(6,4) * t250 + Icges(6,5) * t217;
t53 = t121 * t215 + t127 * t247 + t133 * t216;
t77 = t171 * t215 + t174 * t247 + t177 * t216;
t2 = t247 * t52 + t250 * t53 + t267 * t77;
t122 = Icges(5,5) * t216 - Icges(5,6) * t215 + Icges(5,3) * t247;
t128 = Icges(5,4) * t216 - Icges(5,2) * t215 + Icges(5,6) * t247;
t134 = Icges(5,1) * t216 - Icges(5,4) * t215 + Icges(5,5) * t247;
t54 = t122 * t247 - t128 * t215 + t134 * t216;
t123 = Icges(5,5) * t218 - Icges(5,6) * t217 + Icges(5,3) * t250;
t129 = Icges(5,4) * t218 - Icges(5,2) * t217 + Icges(5,6) * t250;
t135 = Icges(5,1) * t218 - Icges(5,4) * t217 + Icges(5,5) * t250;
t55 = t123 * t247 - t129 * t215 + t135 * t216;
t78 = t172 * t247 - t175 * t215 + t178 * t216;
t3 = t247 * t54 + t250 * t55 + t267 * t78;
t319 = -t3 / 0.2e1 - t2 / 0.2e1 - t1 / 0.2e1;
t56 = -t118 * t250 + t124 * t217 + t130 * t218;
t57 = -t119 * t250 + t125 * t217 + t131 * t218;
t79 = -t170 * t250 + t173 * t217 + t176 * t218;
t4 = t247 * t56 + t250 * t57 + t267 * t79;
t58 = t120 * t217 + t126 * t250 + t132 * t218;
t59 = t121 * t217 + t127 * t250 + t133 * t218;
t80 = t171 * t217 + t174 * t250 + t177 * t218;
t5 = t247 * t58 + t250 * t59 + t267 * t80;
t60 = t122 * t250 - t128 * t217 + t134 * t218;
t61 = t123 * t250 - t129 * t217 + t135 * t218;
t81 = t172 * t250 - t175 * t217 + t178 * t218;
t6 = t247 * t60 + t250 * t61 + t267 * t81;
t318 = t4 / 0.2e1 + t6 / 0.2e1 + t5 / 0.2e1;
t7 = t270 * t50 + t272 * t51 - t348 * t76;
t8 = t270 * t52 + t272 * t53 - t348 * t77;
t9 = t270 * t54 + t272 * t55 - t348 * t78;
t317 = t9 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1;
t316 = -t200 + t329;
t315 = t158 * t348 + t270 * t198 + t334;
t314 = -t235 + t325;
t313 = t199 * t349 + t200 * t347 + t331;
t10 = t270 * t56 + t272 * t57 - t348 * t79;
t11 = t270 * t58 + t272 * t59 - t348 * t80;
t12 = t270 * t60 + t272 * t61 - t348 * t81;
t312 = t12 / 0.2e1 + t11 / 0.2e1 + t10 / 0.2e1;
t13 = t76 * t285 + (t289 * t51 - t291 * t50) * t284;
t14 = t77 * t285 + (t289 * t53 - t291 * t52) * t284;
t15 = t78 * t285 + (t289 * t55 - t291 * t54) * t284;
t311 = t13 / 0.2e1 + t15 / 0.2e1 + t14 / 0.2e1;
t16 = t79 * t285 + (t289 * t57 - t291 * t56) * t284;
t17 = t80 * t285 + (t289 * t59 - t291 * t58) * t284;
t18 = t81 * t285 + (t289 * t61 - t291 * t60) * t284;
t310 = t18 / 0.2e1 + t17 / 0.2e1 + t16 / 0.2e1;
t63 = -t118 * t267 + t124 * t245 + t130 * t246;
t64 = -t119 * t267 + t125 * t245 + t131 * t246;
t87 = t93 * t267;
t19 = t63 * t247 + t64 * t250 + t87;
t65 = t120 * t245 + t126 * t267 + t132 * t246;
t66 = t121 * t245 + t127 * t267 + t133 * t246;
t88 = t94 * t267;
t20 = t65 * t247 + t66 * t250 + t88;
t67 = t122 * t267 - t128 * t245 + t134 * t246;
t68 = t123 * t267 - t129 * t245 + t135 * t246;
t89 = t95 * t267;
t21 = t67 * t247 + t68 * t250 + t89;
t309 = -t20 / 0.2e1 - t19 / 0.2e1 - t21 / 0.2e1;
t22 = t63 * t270 + t64 * t272 - t348 * t93;
t23 = t65 * t270 + t66 * t272 - t348 * t94;
t24 = t67 * t270 + t68 * t272 - t348 * t95;
t308 = t22 / 0.2e1 + t24 / 0.2e1 + t23 / 0.2e1;
t90 = t93 * t285;
t25 = t90 + (t64 * t289 - t63 * t291) * t284;
t91 = t94 * t285;
t26 = t91 + (t66 * t289 - t65 * t291) * t284;
t92 = t95 * t285;
t27 = t92 + (t68 * t289 - t67 * t291) * t284;
t307 = t25 / 0.2e1 + t27 / 0.2e1 + t26 / 0.2e1;
t306 = t284 * (-t274 + t335);
t305 = -rSges(6,2) * t247 - rSges(6,3) * t215;
t304 = t237 + t330;
t303 = t284 * (-t274 + t324);
t302 = t158 * t349 + t159 * t347 + t313;
t301 = t244 + t304;
t300 = t284 * (-t274 + t314);
t299 = -t236 + t321;
t188 = rSges(4,1) * t249 - rSges(4,2) * t247 + rSges(4,3) * t270;
t138 = rSges(5,1) * t216 - rSges(5,2) * t215 + rSges(5,3) * t247;
t230 = t271 * rSges(3,1) - t270 * rSges(3,2) - rSges(3,3) * t347;
t298 = t78 / 0.2e1 + t67 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1 + t77 / 0.2e1 + t76 / 0.2e1;
t297 = t81 / 0.2e1 + t80 / 0.2e1 + t64 / 0.2e1 + t68 / 0.2e1 + t66 / 0.2e1 + t79 / 0.2e1;
t296 = -t199 + t299;
t295 = t159 + t301;
t294 = -t205 + t296;
t182 = Icges(4,5) * t249 - Icges(4,6) * t247 + Icges(4,3) * t270;
t184 = Icges(4,4) * t249 - Icges(4,2) * t247 + Icges(4,6) * t270;
t186 = Icges(4,1) * t249 - Icges(4,4) * t247 + Icges(4,5) * t270;
t100 = -t182 * t348 - t184 * t267 + t186 * t269;
t107 = t220 * t270 - t221 * t247 + t222 * t249;
t293 = t100 / 0.2e1 + t107 / 0.2e1 + t298;
t183 = Icges(4,5) * t252 - Icges(4,6) * t250 + Icges(4,3) * t272;
t185 = Icges(4,4) * t252 - Icges(4,2) * t250 + Icges(4,6) * t272;
t187 = Icges(4,1) * t252 - Icges(4,4) * t250 + Icges(4,5) * t272;
t101 = -t183 * t348 - t185 * t267 + t187 * t269;
t108 = t220 * t272 - t221 * t250 + t222 * t252;
t292 = t108 / 0.2e1 + t101 / 0.2e1 + t297;
t276 = rSges(2,1) * t291 - t289 * rSges(2,2);
t275 = -t289 * rSges(2,1) - rSges(2,2) * t291;
t259 = rSges(3,3) * t285 + (rSges(3,1) * t288 + rSges(3,2) * t290) * t284;
t229 = Icges(3,1) * t273 - Icges(3,4) * t272 + Icges(3,5) * t349;
t228 = Icges(3,1) * t271 - Icges(3,4) * t270 - Icges(3,5) * t347;
t227 = Icges(3,4) * t273 - Icges(3,2) * t272 + Icges(3,6) * t349;
t226 = Icges(3,4) * t271 - Icges(3,2) * t270 - Icges(3,6) * t347;
t225 = Icges(3,5) * t273 - Icges(3,6) * t272 + Icges(3,3) * t349;
t204 = t231 + t330;
t203 = -t230 + t321;
t193 = -t285 * t230 - t259 * t347;
t192 = t231 * t285 - t259 * t349;
t169 = t323 * t285;
t166 = t247 * t198;
t162 = (t230 * t289 + t231 * t291) * t284;
t161 = t256 * t349 - t257 * t272 + t258 * t273;
t160 = -t256 * t347 - t270 * t257 + t271 * t258;
t146 = t267 * t159;
t145 = t304 + t189;
t144 = -t188 + t299;
t142 = t250 * t158;
t137 = rSges(6,1) * t216 - t305;
t117 = -t189 * t348 - t223 * t272;
t116 = t188 * t348 + t223 * t270;
t115 = t225 * t285 + (t227 * t290 + t229 * t288) * t284;
t114 = t224 * t285 + (t226 * t290 + t228 * t288) * t284;
t112 = t113 * t285;
t111 = t188 * t272 - t189 * t270;
t110 = (-t188 - t236) * t285 + t291 * t320;
t109 = t189 * t285 + t289 * t320 + t234;
t106 = t141 + t301 + t353;
t105 = -t138 + t296;
t104 = (t188 * t289 + t189 * t291) * t284 + t331;
t103 = t141 * t267 - t181 * t250;
t102 = -t138 * t267 + t181 * t247;
t99 = t183 * t272 - t185 * t250 + t187 * t252;
t98 = t182 * t272 - t184 * t250 + t186 * t252;
t97 = t183 * t270 - t185 * t247 + t187 * t249;
t96 = t182 * t270 - t184 * t247 + t186 * t249;
t86 = t138 * t250 - t141 * t247;
t85 = t140 + t295 + t353;
t84 = (-rSges(6,1) - pkin(4)) * t216 + t294 + t305;
t83 = t272 * t335 + t339 * t348;
t82 = t138 * t348 + t181 * t270 + t334;
t75 = (-t138 + t332) * t285 + t291 * t306;
t74 = t141 * t285 + t289 * t306 + t333;
t73 = (pkin(10) - qJ(6)) * t250 + t295 + t357;
t72 = (-pkin(4) - t359) * t216 + t294 + t360;
t71 = t138 * t272 + t270 * t339 + t191;
t70 = t140 * t267 + t250 * t336 + t146;
t69 = t180 * t247 + t166 + (-t137 - t158) * t267;
t62 = (t138 * t289 + t141 * t291) * t284 + t313;
t49 = t272 * t324 + t328 * t348;
t48 = t137 * t348 + t180 * t270 + t315;
t47 = (-t137 + t326) * t285 + t291 * t303;
t46 = t140 * t285 + t289 * t303 + t327;
t45 = t137 * t250 + t247 * t340 + t142;
t44 = t250 * t325 + t267 * t341 + t146;
t43 = t166 + t337 * t247 + (-t158 - t342) * t267;
t42 = t272 * t314 + t316 * t348;
t41 = t270 * t337 + t342 * t348 + t315;
t40 = t137 * t272 + t270 * t328 + t338;
t39 = (t326 - t342) * t285 + t291 * t300;
t38 = t285 * t341 + t289 * t300 + t327;
t37 = (t137 * t289 + t140 * t291) * t284 + t302;
t36 = t112 + (-t100 * t291 + t101 * t289) * t284;
t35 = t100 * t270 + t101 * t272 - t113 * t348;
t34 = t247 * t329 + t250 * t342 + t142;
t33 = t108 * t285 + (t289 * t99 - t291 * t98) * t284;
t32 = t107 * t285 + (t289 * t97 - t291 * t96) * t284;
t31 = -t108 * t348 + t270 * t98 + t272 * t99;
t30 = -t107 * t348 + t270 * t96 + t272 * t97;
t29 = t270 * t316 + t272 * t342 + t338;
t28 = (t289 * t342 + t291 * t341) * t284 + t302;
t136 = [m(7) * (t72 ^ 2 + t73 ^ 2) + m(6) * (t84 ^ 2 + t85 ^ 2) + m(5) * (t105 ^ 2 + t106 ^ 2) + m(4) * (t144 ^ 2 + t145 ^ 2) + m(3) * (t203 ^ 2 + t204 ^ 2) + m(2) * (t275 ^ 2 + t276 ^ 2) + Icges(2,3) + t323 - t356; t90 + t169 + t112 + t91 + t92 + m(7) * (t38 * t73 + t39 * t72) + m(6) * (t46 * t85 + t47 * t84) + m(5) * (t105 * t75 + t106 * t74) + m(4) * (t109 * t145 + t110 * t144) + m(3) * (t192 * t204 + t193 * t203) + ((-t114 / 0.2e1 - t160 / 0.2e1 - t293) * t291 + (t115 / 0.2e1 + t161 / 0.2e1 + t292) * t289) * t284; (t25 + t27 + t26 + t36 + t169) * t285 + m(7) * (t28 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t37 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t62 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t104 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(3) * (t162 ^ 2 + t192 ^ 2 + t193 ^ 2) + ((-t15 - t14 - t13 - t32 + (-t270 * t226 + t271 * t228 - t351 * t284) * t347) * t291 + (t17 + t16 + t18 + t33 + ((-t227 * t272 + t229 * t273 + (t225 * t289 - t351) * t284) * t289 + (t225 * t347 + t226 * t272 + t270 * t227 - t228 * t273 - t271 * t229) * t291) * t284) * t289 + ((-t114 - t160) * t291 + (t115 + t161) * t289) * t285) * t284; t356 * t348 + m(7) * (t41 * t72 + t42 * t73) + m(6) * (t48 * t84 + t49 * t85) + m(5) * (t105 * t82 + t106 * t83) + m(4) * (t116 * t144 + t117 * t145) + t292 * t272 + t293 * t270; (t35 / 0.2e1 + t308) * t285 + (t33 / 0.2e1 + t310) * t272 + (t32 / 0.2e1 + t311) * t270 + m(7) * (t28 * t29 + t38 * t42 + t39 * t41) + m(6) * (t37 * t40 + t46 * t49 + t47 * t48) + m(5) * (t62 * t71 + t74 * t83 + t75 * t82) + m(4) * (t104 * t111 + t109 * t117 + t110 * t116) + ((-t30 / 0.2e1 - t317) * t291 + (-t36 / 0.2e1 - t307) * t290 + (t31 / 0.2e1 + t312) * t289) * t284; (-t22 - t23 - t24 - t35) * t348 + (t10 + t12 + t11 + t31) * t272 + (t7 + t9 + t8 + t30) * t270 + m(7) * (t29 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t71 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(4) * (t111 ^ 2 + t116 ^ 2 + t117 ^ 2); t89 + t88 + t87 + m(7) * (t43 * t72 + t44 * t73) + m(6) * (t69 * t84 + t70 * t85) + m(5) * (t102 * t105 + t103 * t106) + t297 * t250 + t298 * t247; -t309 * t285 + t307 * t267 + t310 * t250 + t311 * t247 + m(7) * (t28 * t34 + t38 * t44 + t39 * t43) + m(6) * (t45 * t37 + t46 * t70 + t69 * t47) + m(5) * (t102 * t75 + t103 * t74 + t62 * t86) + (t289 * t318 + t291 * t319) * t284; t309 * t348 + t318 * t272 - t319 * t270 + t308 * t267 + t312 * t250 + t317 * t247 + m(7) * (t29 * t34 + t41 * t43 + t42 * t44) + m(6) * (t45 * t40 + t69 * t48 + t49 * t70) + m(5) * (t102 * t82 + t103 * t83 + t71 * t86); (t20 + t21 + t19) * t267 + (t5 + t4 + t6) * t250 + (t2 + t3 + t1) * t247 + m(6) * (t45 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2 + t86 ^ 2); m(7) * (t215 * t73 + t217 * t72) + m(6) * (t215 * t85 + t217 * t84); m(7) * (t215 * t38 + t217 * t39 + t245 * t28) + m(6) * (t215 * t46 + t217 * t47 + t245 * t37); m(7) * (t215 * t42 + t217 * t41 + t245 * t29) + m(6) * (t215 * t49 + t217 * t48 + t245 * t40); m(6) * (t215 * t70 + t217 * t69 + t245 * t45) + m(7) * (t215 * t44 + t217 * t43 + t245 * t34); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t215 ^ 2 + t217 ^ 2 + t245 ^ 2); m(7) * (-t247 * t73 - t250 * t72); m(7) * (-t247 * t38 - t250 * t39 - t267 * t28); m(7) * (-t247 * t42 - t250 * t41 - t267 * t29); m(7) * (-t247 * t44 - t250 * t43 - t267 * t34); m(7) * (-t215 * t247 - t217 * t250 - t245 * t267); m(7) * (t247 ^ 2 + t250 ^ 2 + t267 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t136(1) t136(2) t136(4) t136(7) t136(11) t136(16); t136(2) t136(3) t136(5) t136(8) t136(12) t136(17); t136(4) t136(5) t136(6) t136(9) t136(13) t136(18); t136(7) t136(8) t136(9) t136(10) t136(14) t136(19); t136(11) t136(12) t136(13) t136(14) t136(15) t136(20); t136(16) t136(17) t136(18) t136(19) t136(20) t136(21);];
Mq  = res;
