% Calculate joint inertia matrix for
% S6RRRRPP9
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:55
% EndTime: 2019-03-09 21:42:11
% DurationCPUTime: 6.58s
% Computational Cost: add. (22052->661), mult. (56659->899), div. (0->0), fcn. (72732->10), ass. (0->292)
t283 = cos(pkin(6));
t286 = sin(qJ(1));
t287 = cos(qJ(2));
t341 = t286 * t287;
t285 = sin(qJ(2));
t288 = cos(qJ(1));
t342 = t285 * t288;
t268 = t283 * t342 + t341;
t284 = sin(qJ(3));
t282 = sin(pkin(6));
t344 = t282 * t288;
t352 = cos(qJ(3));
t248 = t268 * t352 - t284 * t344;
t340 = t287 * t288;
t343 = t285 * t286;
t267 = -t283 * t340 + t343;
t350 = sin(qJ(4));
t351 = cos(qJ(4));
t214 = t248 * t350 - t267 * t351;
t319 = t282 * t352;
t247 = t268 * t284 + t288 * t319;
t355 = rSges(7,1) + pkin(5);
t356 = -rSges(7,2) * t214 - t355 * t247;
t354 = rSges(7,3) + qJ(6);
t347 = t282 * t285;
t265 = -t283 * t352 + t284 * t347;
t266 = t283 * t284 + t285 * t319;
t345 = t282 * t287;
t219 = Icges(4,5) * t266 - Icges(4,6) * t265 - Icges(4,3) * t345;
t220 = Icges(4,4) * t266 - Icges(4,2) * t265 - Icges(4,6) * t345;
t221 = Icges(4,1) * t266 - Icges(4,4) * t265 - Icges(4,5) * t345;
t113 = -t219 * t345 - t265 * t220 + t266 * t221;
t245 = t266 * t350 + t345 * t351;
t246 = t266 * t351 - t345 * t350;
t170 = Icges(7,5) * t265 + Icges(7,6) * t245 + Icges(7,3) * t246;
t173 = Icges(7,4) * t265 + Icges(7,2) * t245 + Icges(7,6) * t246;
t176 = Icges(7,1) * t265 + Icges(7,4) * t245 + Icges(7,5) * t246;
t93 = t246 * t170 + t245 * t173 + t265 * t176;
t171 = Icges(6,5) * t265 - Icges(6,6) * t246 + Icges(6,3) * t245;
t174 = Icges(6,4) * t265 - Icges(6,2) * t246 + Icges(6,6) * t245;
t177 = Icges(6,1) * t265 - Icges(6,4) * t246 + Icges(6,5) * t245;
t94 = t245 * t171 - t246 * t174 + t265 * t177;
t172 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t265;
t175 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t265;
t178 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t265;
t95 = t265 * t172 - t245 * t175 + t246 * t178;
t353 = -t113 - t93 - t94 - t95;
t223 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t344;
t348 = t223 * t288;
t346 = t282 * t286;
t215 = t248 * t351 + t267 * t350;
t339 = t354 * t215 - t356;
t270 = -t283 * t343 + t340;
t250 = t270 * t352 + t284 * t346;
t269 = t283 * t341 + t342;
t216 = t250 * t350 - t269 * t351;
t217 = t250 * t351 + t269 * t350;
t249 = t270 * t284 - t286 * t319;
t338 = t216 * rSges(7,2) + t354 * t217 + t249 * t355;
t139 = t249 * rSges(6,1) - t217 * rSges(6,2) + t216 * rSges(6,3);
t159 = t217 * pkin(4) + qJ(5) * t216;
t337 = -t139 - t159;
t141 = t217 * rSges(5,1) - t216 * rSges(5,2) + t249 * rSges(5,3);
t200 = t250 * pkin(3) + pkin(10) * t249;
t336 = -t141 - t200;
t205 = t214 * qJ(5);
t158 = pkin(4) * t215 + t205;
t199 = pkin(3) * t248 + t247 * pkin(10);
t191 = t269 * t199;
t335 = t269 * t158 + t191;
t334 = rSges(7,2) * t245 + t354 * t246 + t265 * t355;
t180 = rSges(6,1) * t265 - rSges(6,2) * t246 + rSges(6,3) * t245;
t198 = pkin(4) * t246 + qJ(5) * t245;
t333 = -t180 - t198;
t181 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t265;
t234 = pkin(3) * t266 + pkin(10) * t265;
t332 = -t181 - t234;
t331 = t199 * t345 + t267 * t234;
t236 = t270 * pkin(2) + pkin(9) * t269;
t233 = t283 * t236;
t330 = t283 * t200 + t233;
t235 = pkin(2) * t268 + t267 * pkin(9);
t329 = -t199 - t235;
t328 = t235 * t346 + t236 * t344;
t327 = t288 * pkin(1) + pkin(8) * t346;
t326 = -t159 - t338;
t325 = -t200 + t337;
t324 = t283 * t159 + t330;
t323 = -t158 + t329;
t322 = -t198 - t334;
t321 = -t234 + t333;
t189 = t250 * rSges(4,1) - t249 * rSges(4,2) + t269 * rSges(4,3);
t254 = Icges(3,3) * t283 + (Icges(3,5) * t285 + Icges(3,6) * t287) * t282;
t255 = Icges(3,6) * t283 + (Icges(3,4) * t285 + Icges(3,2) * t287) * t282;
t256 = Icges(3,5) * t283 + (Icges(3,1) * t285 + Icges(3,4) * t287) * t282;
t320 = t283 * t254 + t255 * t345 + t256 * t347;
t230 = t270 * rSges(3,1) - t269 * rSges(3,2) + rSges(3,3) * t346;
t318 = -t286 * pkin(1) + pkin(8) * t344;
t222 = rSges(4,1) * t266 - rSges(4,2) * t265 - rSges(4,3) * t345;
t271 = (pkin(2) * t285 - pkin(9) * t287) * t282;
t317 = t282 * (-t222 - t271);
t118 = Icges(7,5) * t247 + Icges(7,6) * t214 + Icges(7,3) * t215;
t124 = Icges(7,4) * t247 + Icges(7,2) * t214 + Icges(7,6) * t215;
t130 = Icges(7,1) * t247 + Icges(7,4) * t214 + Icges(7,5) * t215;
t50 = t118 * t215 + t124 * t214 + t130 * t247;
t119 = Icges(7,5) * t249 + Icges(7,6) * t216 + Icges(7,3) * t217;
t125 = Icges(7,4) * t249 + Icges(7,2) * t216 + Icges(7,6) * t217;
t131 = Icges(7,1) * t249 + Icges(7,4) * t216 + Icges(7,5) * t217;
t51 = t119 * t215 + t125 * t214 + t131 * t247;
t76 = t170 * t215 + t173 * t214 + t176 * t247;
t1 = t247 * t50 + t249 * t51 + t265 * t76;
t120 = Icges(6,5) * t247 - Icges(6,6) * t215 + Icges(6,3) * t214;
t126 = Icges(6,4) * t247 - Icges(6,2) * t215 + Icges(6,6) * t214;
t132 = Icges(6,1) * t247 - Icges(6,4) * t215 + Icges(6,5) * t214;
t52 = t120 * t214 - t126 * t215 + t132 * t247;
t121 = Icges(6,5) * t249 - Icges(6,6) * t217 + Icges(6,3) * t216;
t127 = Icges(6,4) * t249 - Icges(6,2) * t217 + Icges(6,6) * t216;
t133 = Icges(6,1) * t249 - Icges(6,4) * t217 + Icges(6,5) * t216;
t53 = t121 * t214 - t127 * t215 + t133 * t247;
t77 = t171 * t214 - t174 * t215 + t177 * t247;
t2 = t247 * t52 + t249 * t53 + t265 * t77;
t122 = Icges(5,5) * t215 - Icges(5,6) * t214 + Icges(5,3) * t247;
t128 = Icges(5,4) * t215 - Icges(5,2) * t214 + Icges(5,6) * t247;
t134 = Icges(5,1) * t215 - Icges(5,4) * t214 + Icges(5,5) * t247;
t58 = t122 * t247 - t128 * t214 + t134 * t215;
t123 = Icges(5,5) * t217 - Icges(5,6) * t216 + Icges(5,3) * t249;
t129 = Icges(5,4) * t217 - Icges(5,2) * t216 + Icges(5,6) * t249;
t135 = Icges(5,1) * t217 - Icges(5,4) * t216 + Icges(5,5) * t249;
t59 = t123 * t247 - t129 * t214 + t135 * t215;
t80 = t172 * t247 - t175 * t214 + t178 * t215;
t5 = t247 * t58 + t249 * t59 + t265 * t80;
t316 = t2 / 0.2e1 + t1 / 0.2e1 + t5 / 0.2e1;
t54 = t118 * t217 + t124 * t216 + t130 * t249;
t55 = t119 * t217 + t125 * t216 + t131 * t249;
t78 = t170 * t217 + t173 * t216 + t176 * t249;
t3 = t247 * t54 + t249 * t55 + t265 * t78;
t56 = t120 * t216 - t126 * t217 + t132 * t249;
t57 = t121 * t216 - t127 * t217 + t133 * t249;
t79 = t171 * t216 - t174 * t217 + t177 * t249;
t4 = t247 * t56 + t249 * t57 + t265 * t79;
t60 = t122 * t249 - t128 * t216 + t134 * t217;
t61 = t123 * t249 - t129 * t216 + t135 * t217;
t81 = t172 * t249 - t175 * t216 + t178 * t217;
t6 = t247 * t60 + t249 * t61 + t265 * t81;
t315 = t3 / 0.2e1 + t4 / 0.2e1 + t6 / 0.2e1;
t11 = t267 * t58 + t269 * t59 - t345 * t80;
t7 = t267 * t50 + t269 * t51 - t345 * t76;
t8 = t267 * t52 + t269 * t53 - t345 * t77;
t314 = t7 / 0.2e1 + t11 / 0.2e1 + t8 / 0.2e1;
t313 = -t200 + t326;
t312 = t158 * t345 + t267 * t198 + t331;
t311 = -t234 + t322;
t310 = t199 * t346 + t200 * t344 + t328;
t10 = t267 * t56 + t269 * t57 - t345 * t79;
t12 = t267 * t60 + t269 * t61 - t345 * t81;
t9 = t267 * t54 + t269 * t55 - t345 * t78;
t309 = t10 / 0.2e1 + t9 / 0.2e1 + t12 / 0.2e1;
t13 = t76 * t283 + (t286 * t51 - t288 * t50) * t282;
t14 = t77 * t283 + (t286 * t53 - t288 * t52) * t282;
t17 = t80 * t283 + (t286 * t59 - t288 * t58) * t282;
t308 = t14 / 0.2e1 + t13 / 0.2e1 + t17 / 0.2e1;
t15 = t78 * t283 + (t286 * t55 - t288 * t54) * t282;
t16 = t79 * t283 + (t286 * t57 - t288 * t56) * t282;
t18 = t81 * t283 + (t286 * t61 - t288 * t60) * t282;
t307 = t18 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1;
t63 = t118 * t246 + t124 * t245 + t130 * t265;
t64 = t119 * t246 + t125 * t245 + t131 * t265;
t87 = t93 * t265;
t19 = t63 * t247 + t64 * t249 + t87;
t65 = t120 * t245 - t126 * t246 + t132 * t265;
t66 = t121 * t245 - t127 * t246 + t133 * t265;
t88 = t94 * t265;
t20 = t65 * t247 + t66 * t249 + t88;
t67 = t122 * t265 - t128 * t245 + t134 * t246;
t68 = t123 * t265 - t129 * t245 + t135 * t246;
t89 = t95 * t265;
t21 = t67 * t247 + t68 * t249 + t89;
t306 = t19 / 0.2e1 + t21 / 0.2e1 + t20 / 0.2e1;
t22 = t63 * t267 + t64 * t269 - t345 * t93;
t23 = t65 * t267 + t66 * t269 - t345 * t94;
t24 = t67 * t267 + t68 * t269 - t345 * t95;
t305 = t22 / 0.2e1 + t23 / 0.2e1 + t24 / 0.2e1;
t90 = t93 * t283;
t25 = t90 + (t64 * t286 - t63 * t288) * t282;
t91 = t94 * t283;
t26 = t91 + (t66 * t286 - t65 * t288) * t282;
t92 = t95 * t283;
t27 = t92 + (t68 * t286 - t67 * t288) * t282;
t304 = t25 / 0.2e1 + t26 / 0.2e1 + t27 / 0.2e1;
t303 = t282 * (-t271 + t332);
t302 = -rSges(6,1) * t247 - rSges(6,3) * t214;
t301 = t236 + t327;
t300 = t282 * (-t271 + t321);
t299 = t158 * t346 + t159 * t344 + t310;
t298 = t282 * (-t271 + t311);
t297 = -t235 + t318;
t188 = rSges(4,1) * t248 - rSges(4,2) * t247 + rSges(4,3) * t267;
t140 = rSges(5,1) * t215 - rSges(5,2) * t214 + rSges(5,3) * t247;
t229 = t268 * rSges(3,1) - t267 * rSges(3,2) - rSges(3,3) * t344;
t296 = t67 / 0.2e1 + t76 / 0.2e1 + t80 / 0.2e1 + t77 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1;
t295 = t81 / 0.2e1 + t78 / 0.2e1 + t64 / 0.2e1 + t66 / 0.2e1 + t68 / 0.2e1 + t79 / 0.2e1;
t294 = t200 + t301;
t293 = -t199 + t297;
t292 = -t205 + t293;
t182 = Icges(4,5) * t248 - Icges(4,6) * t247 + Icges(4,3) * t267;
t184 = Icges(4,4) * t248 - Icges(4,2) * t247 + Icges(4,6) * t267;
t186 = Icges(4,1) * t248 - Icges(4,4) * t247 + Icges(4,5) * t267;
t100 = -t182 * t345 - t184 * t265 + t186 * t266;
t107 = t219 * t267 - t220 * t247 + t221 * t248;
t291 = t107 / 0.2e1 + t100 / 0.2e1 + t296;
t183 = Icges(4,5) * t250 - Icges(4,6) * t249 + Icges(4,3) * t269;
t185 = Icges(4,4) * t250 - Icges(4,2) * t249 + Icges(4,6) * t269;
t187 = Icges(4,1) * t250 - Icges(4,4) * t249 + Icges(4,5) * t269;
t101 = -t183 * t345 - t185 * t265 + t187 * t266;
t108 = t219 * t269 - t220 * t249 + t221 * t250;
t290 = t108 / 0.2e1 + t101 / 0.2e1 + t295;
t289 = t159 + t294;
t273 = rSges(2,1) * t288 - t286 * rSges(2,2);
t272 = -t286 * rSges(2,1) - rSges(2,2) * t288;
t257 = rSges(3,3) * t283 + (rSges(3,1) * t285 + rSges(3,2) * t287) * t282;
t228 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t346;
t227 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t344;
t226 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t346;
t225 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t344;
t224 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t346;
t204 = t230 + t327;
t203 = -t229 + t318;
t193 = -t283 * t229 - t257 * t344;
t192 = t230 * t283 - t257 * t346;
t169 = t320 * t283;
t166 = t247 * t198;
t162 = (t229 * t286 + t230 * t288) * t282;
t161 = t254 * t346 - t255 * t269 + t256 * t270;
t160 = -t254 * t344 - t267 * t255 + t268 * t256;
t146 = t265 * t159;
t145 = t301 + t189;
t144 = -t188 + t297;
t142 = t249 * t158;
t137 = -rSges(6,2) * t215 - t302;
t117 = -t189 * t345 - t222 * t269;
t116 = t188 * t345 + t222 * t267;
t115 = t224 * t283 + (t226 * t287 + t228 * t285) * t282;
t114 = t223 * t283 + (t225 * t287 + t227 * t285) * t282;
t112 = t113 * t283;
t111 = t188 * t269 - t189 * t267;
t110 = (-t188 - t235) * t283 + t288 * t317;
t109 = t189 * t283 + t286 * t317 + t233;
t106 = t294 + t141;
t105 = -t140 + t293;
t104 = (t188 * t286 + t189 * t288) * t282 + t328;
t103 = t141 * t265 - t181 * t249;
t102 = -t140 * t265 + t181 * t247;
t99 = t183 * t269 - t185 * t249 + t187 * t250;
t98 = t182 * t269 - t184 * t249 + t186 * t250;
t97 = t183 * t267 - t185 * t247 + t187 * t248;
t96 = t182 * t267 - t184 * t247 + t186 * t248;
t86 = t140 * t249 - t141 * t247;
t85 = t289 + t139;
t84 = (rSges(6,2) - pkin(4)) * t215 + t292 + t302;
t83 = t269 * t332 + t336 * t345;
t82 = t140 * t345 + t181 * t267 + t331;
t75 = (-t140 + t329) * t283 + t288 * t303;
t74 = t141 * t283 + t286 * t303 + t330;
t73 = t289 + t338;
t72 = (-pkin(4) - t354) * t215 + t292 + t356;
t71 = t140 * t269 + t267 * t336 + t191;
t70 = t139 * t265 + t249 * t333 + t146;
t69 = t180 * t247 + t166 + (-t137 - t158) * t265;
t62 = (t140 * t286 + t141 * t288) * t282 + t310;
t49 = t269 * t321 + t325 * t345;
t48 = t137 * t345 + t180 * t267 + t312;
t47 = (-t137 + t323) * t283 + t288 * t300;
t46 = t139 * t283 + t286 * t300 + t324;
t45 = t137 * t249 + t247 * t337 + t142;
t44 = t249 * t322 + t265 * t338 + t146;
t43 = t166 + t334 * t247 + (-t158 - t339) * t265;
t42 = t269 * t311 + t313 * t345;
t41 = t267 * t334 + t339 * t345 + t312;
t40 = t137 * t269 + t267 * t325 + t335;
t39 = (t323 - t339) * t283 + t288 * t298;
t38 = t283 * t338 + t286 * t298 + t324;
t37 = (t137 * t286 + t139 * t288) * t282 + t299;
t36 = t112 + (-t100 * t288 + t101 * t286) * t282;
t35 = t100 * t267 + t101 * t269 - t113 * t345;
t34 = t247 * t326 + t249 * t339 + t142;
t33 = t108 * t283 + (t286 * t99 - t288 * t98) * t282;
t32 = t107 * t283 + (t286 * t97 - t288 * t96) * t282;
t31 = -t108 * t345 + t267 * t98 + t269 * t99;
t30 = -t107 * t345 + t267 * t96 + t269 * t97;
t29 = t267 * t313 + t269 * t339 + t335;
t28 = (t286 * t339 + t288 * t338) * t282 + t299;
t136 = [m(6) * (t84 ^ 2 + t85 ^ 2) + m(7) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t105 ^ 2 + t106 ^ 2) + m(4) * (t144 ^ 2 + t145 ^ 2) + m(3) * (t203 ^ 2 + t204 ^ 2) + m(2) * (t272 ^ 2 + t273 ^ 2) + Icges(2,3) + t320 - t353; t112 + t169 + t90 + t92 + t91 + m(3) * (t192 * t204 + t193 * t203) + m(4) * (t109 * t145 + t110 * t144) + m(5) * (t105 * t75 + t106 * t74) + m(6) * (t46 * t85 + t47 * t84) + m(7) * (t38 * t73 + t39 * t72) + ((-t114 / 0.2e1 - t160 / 0.2e1 - t291) * t288 + (t115 / 0.2e1 + t161 / 0.2e1 + t290) * t286) * t282; (t25 + t26 + t27 + t36 + t169) * t283 + m(7) * (t28 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t37 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t62 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t104 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(3) * (t162 ^ 2 + t192 ^ 2 + t193 ^ 2) + ((-t13 - t17 - t14 - t32 + (-t267 * t225 + t268 * t227 - t282 * t348) * t344) * t288 + (t18 + t16 + t15 + t33 + ((-t226 * t269 + t228 * t270 + (t224 * t286 - t348) * t282) * t286 + (t224 * t344 + t225 * t269 + t267 * t226 - t227 * t270 - t268 * t228) * t288) * t282) * t286 + ((-t114 - t160) * t288 + (t115 + t161) * t286) * t283) * t282; t353 * t345 + m(6) * (t48 * t84 + t49 * t85) + m(7) * (t41 * t72 + t42 * t73) + m(5) * (t105 * t82 + t106 * t83) + m(4) * (t116 * t144 + t117 * t145) + t290 * t269 + t291 * t267; (t35 / 0.2e1 + t305) * t283 + (t33 / 0.2e1 + t307) * t269 + (t32 / 0.2e1 + t308) * t267 + m(4) * (t104 * t111 + t109 * t117 + t110 * t116) + m(5) * (t62 * t71 + t74 * t83 + t75 * t82) + m(7) * (t28 * t29 + t38 * t42 + t39 * t41) + m(6) * (t37 * t40 + t46 * t49 + t47 * t48) + ((-t30 / 0.2e1 - t314) * t288 + (-t36 / 0.2e1 - t304) * t287 + (t31 / 0.2e1 + t309) * t286) * t282; (-t22 - t23 - t24 - t35) * t345 + (t9 + t10 + t12 + t31) * t269 + (t8 + t7 + t11 + t30) * t267 + m(7) * (t29 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t71 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(4) * (t111 ^ 2 + t116 ^ 2 + t117 ^ 2); t89 + t88 + t87 + m(6) * (t69 * t84 + t70 * t85) + m(7) * (t43 * t72 + t44 * t73) + m(5) * (t102 * t105 + t103 * t106) + t295 * t249 + t296 * t247; t306 * t283 + t304 * t265 + t307 * t249 + t308 * t247 + m(7) * (t28 * t34 + t38 * t44 + t39 * t43) + m(6) * (t45 * t37 + t46 * t70 + t69 * t47) + m(5) * (t102 * t75 + t103 * t74 + t62 * t86) + (t286 * t315 - t288 * t316) * t282; -t306 * t345 + t315 * t269 + t316 * t267 + t305 * t265 + t309 * t249 + t314 * t247 + m(7) * (t29 * t34 + t41 * t43 + t42 * t44) + m(6) * (t45 * t40 + t69 * t48 + t49 * t70) + m(5) * (t102 * t82 + t103 * t83 + t71 * t86); (t19 + t20 + t21) * t265 + (t3 + t4 + t6) * t249 + (t2 + t5 + t1) * t247 + m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t45 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2 + t86 ^ 2); m(6) * (t214 * t85 + t216 * t84) + m(7) * (t214 * t73 + t216 * t72); m(7) * (t214 * t38 + t216 * t39 + t245 * t28) + m(6) * (t214 * t46 + t216 * t47 + t245 * t37); m(7) * (t214 * t42 + t216 * t41 + t245 * t29) + m(6) * (t214 * t49 + t216 * t48 + t245 * t40); m(7) * (t214 * t44 + t216 * t43 + t245 * t34) + m(6) * (t214 * t70 + t216 * t69 + t245 * t45); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t214 ^ 2 + t216 ^ 2 + t245 ^ 2); m(7) * (t215 * t73 + t217 * t72); m(7) * (t215 * t38 + t217 * t39 + t246 * t28); m(7) * (t215 * t42 + t217 * t41 + t246 * t29); m(7) * (t215 * t44 + t217 * t43 + t246 * t34); m(7) * (t214 * t215 + t216 * t217 + t245 * t246); m(7) * (t215 ^ 2 + t217 ^ 2 + t246 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t136(1) t136(2) t136(4) t136(7) t136(11) t136(16); t136(2) t136(3) t136(5) t136(8) t136(12) t136(17); t136(4) t136(5) t136(6) t136(9) t136(13) t136(18); t136(7) t136(8) t136(9) t136(10) t136(14) t136(19); t136(11) t136(12) t136(13) t136(14) t136(15) t136(20); t136(16) t136(17) t136(18) t136(19) t136(20) t136(21);];
Mq  = res;
