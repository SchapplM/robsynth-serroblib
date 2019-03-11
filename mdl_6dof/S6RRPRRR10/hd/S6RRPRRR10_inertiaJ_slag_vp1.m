% Calculate joint inertia matrix for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:23
% EndTime: 2019-03-09 14:18:39
% DurationCPUTime: 7.03s
% Computational Cost: add. (33383->668), mult. (49831->917), div. (0->0), fcn. (63278->14), ass. (0->303)
t298 = sin(pkin(12));
t300 = cos(pkin(12));
t301 = cos(pkin(6));
t299 = sin(pkin(6));
t304 = sin(qJ(2));
t363 = t299 * t304;
t274 = -t298 * t363 + t300 * t301;
t364 = t298 * t301;
t275 = t300 * t363 + t364;
t307 = cos(qJ(2));
t361 = t299 * t307;
t215 = Icges(4,4) * t275 + Icges(4,2) * t274 - Icges(4,6) * t361;
t216 = Icges(4,1) * t275 + Icges(4,4) * t274 - Icges(4,5) * t361;
t260 = Icges(3,3) * t301 + (Icges(3,5) * t304 + Icges(3,6) * t307) * t299;
t261 = Icges(3,6) * t301 + (Icges(3,4) * t304 + Icges(3,2) * t307) * t299;
t262 = Icges(3,5) * t301 + (Icges(3,1) * t304 + Icges(3,4) * t307) * t299;
t380 = t274 * t215 + t275 * t216 + t301 * t260 + t261 * t361 + t262 * t363;
t342 = pkin(12) + qJ(4);
t293 = sin(t342);
t326 = cos(t342);
t317 = t299 * t326;
t264 = t301 * t293 + t304 * t317;
t297 = qJ(5) + qJ(6);
t294 = sin(t297);
t295 = cos(t297);
t233 = -t264 * t294 - t295 * t361;
t234 = t264 * t295 - t294 * t361;
t263 = t293 * t363 - t301 * t326;
t151 = Icges(7,5) * t234 + Icges(7,6) * t233 + Icges(7,3) * t263;
t152 = Icges(7,4) * t234 + Icges(7,2) * t233 + Icges(7,6) * t263;
t153 = Icges(7,1) * t234 + Icges(7,4) * t233 + Icges(7,5) * t263;
t78 = t263 * t151 + t233 * t152 + t234 * t153;
t303 = sin(qJ(5));
t306 = cos(qJ(5));
t244 = -t264 * t303 - t306 * t361;
t336 = t303 * t361;
t245 = t264 * t306 - t336;
t156 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t263;
t157 = Icges(6,4) * t245 + Icges(6,2) * t244 + Icges(6,6) * t263;
t158 = Icges(6,1) * t245 + Icges(6,4) * t244 + Icges(6,5) * t263;
t82 = t263 * t156 + t244 * t157 + t245 * t158;
t379 = -t78 - t82;
t214 = Icges(4,5) * t275 + Icges(4,6) * t274 - Icges(4,3) * t361;
t378 = (-t214 * t361 + t380) * t301;
t305 = sin(qJ(1));
t357 = t305 * t307;
t308 = cos(qJ(1));
t358 = t304 * t308;
t277 = t301 * t358 + t357;
t246 = t277 * t293 + t308 * t317;
t377 = t246 / 0.2e1;
t356 = t307 * t308;
t359 = t304 * t305;
t279 = -t301 * t359 + t356;
t248 = t279 * t293 - t305 * t317;
t376 = t248 / 0.2e1;
t375 = t263 / 0.2e1;
t276 = -t301 * t356 + t359;
t374 = t276 / 0.2e1;
t278 = t301 * t357 + t358;
t373 = t278 / 0.2e1;
t372 = t301 / 0.2e1;
t291 = pkin(3) * t300 + pkin(2);
t371 = -pkin(2) + t291;
t292 = pkin(5) * t306 + pkin(4);
t370 = -pkin(4) + t292;
t360 = t299 * t308;
t247 = t277 * t326 - t293 * t360;
t198 = -t247 * t294 + t276 * t295;
t199 = t247 * t295 + t276 * t294;
t119 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t246;
t121 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t246;
t123 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t246;
t62 = t119 * t263 + t121 * t233 + t123 * t234;
t369 = t62 * t246;
t362 = t299 * t305;
t249 = t279 * t326 + t293 * t362;
t200 = -t249 * t294 + t278 * t295;
t201 = t249 * t295 + t278 * t294;
t120 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t248;
t122 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t248;
t124 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t248;
t63 = t120 * t263 + t122 * t233 + t124 * t234;
t368 = t63 * t248;
t302 = -pkin(9) - qJ(3);
t367 = t276 * t302;
t366 = t276 * t303;
t365 = t278 * t303;
t242 = t246 * pkin(10);
t309 = -pkin(11) - pkin(10);
t115 = pkin(5) * t366 - t246 * t309 + t247 * t370 - t242;
t318 = -t199 * rSges(7,1) - t198 * rSges(7,2);
t125 = t246 * rSges(7,3) - t318;
t355 = t115 + t125;
t192 = t249 * pkin(4) + pkin(10) * t248;
t331 = pkin(5) * t365 - t248 * t309 + t249 * t292;
t116 = -t192 + t331;
t126 = t201 * rSges(7,1) + t200 * rSges(7,2) + t248 * rSges(7,3);
t354 = t116 + t126;
t206 = -t249 * t303 + t278 * t306;
t207 = t249 * t306 + t365;
t134 = t207 * rSges(6,1) + t206 * rSges(6,2) + t248 * rSges(6,3);
t353 = -t134 - t192;
t148 = -pkin(5) * t336 + t370 * t264 + (-pkin(10) - t309) * t263;
t154 = rSges(7,1) * t234 + rSges(7,2) * t233 + rSges(7,3) * t263;
t352 = t148 + t154;
t237 = t279 * pkin(2) + qJ(3) * t278;
t338 = t298 * t362;
t329 = pkin(3) * t338 - t278 * t302 + t279 * t291;
t172 = -t237 + t329;
t232 = t301 * t237;
t351 = t301 * t172 + t232;
t268 = t276 * qJ(3);
t337 = t298 * t360;
t283 = pkin(3) * t337;
t171 = t277 * t371 - t268 - t283 - t367;
t236 = pkin(2) * t277 + t268;
t350 = -t171 - t236;
t252 = -t277 * t298 - t300 * t360;
t253 = t277 * t300 - t337;
t181 = rSges(4,1) * t253 + rSges(4,2) * t252 + rSges(4,3) * t276;
t349 = -t181 - t236;
t191 = t247 * pkin(4) + t242;
t218 = t264 * pkin(4) + t263 * pkin(10);
t348 = t191 * t361 + t276 * t218;
t211 = Icges(5,4) * t264 - Icges(5,2) * t263 - Icges(5,6) * t361;
t212 = Icges(5,1) * t264 - Icges(5,4) * t263 - Icges(5,5) * t361;
t347 = -t263 * t211 + t264 * t212;
t280 = (pkin(2) * t304 - qJ(3) * t307) * t299;
t345 = -pkin(3) * t364 - ((qJ(3) + t302) * t307 + t371 * t304) * t299 - t280;
t344 = t236 * t362 + t237 * t360;
t343 = t308 * pkin(1) + pkin(8) * t362;
t74 = t78 * t263;
t26 = t368 + t74 + t369;
t50 = t119 * t246 + t121 * t198 + t123 * t199;
t51 = t120 * t246 + t122 * t198 + t124 * t199;
t69 = t151 * t246 + t152 * t198 + t153 * t199;
t7 = t246 * t50 + t248 * t51 + t263 * t69;
t52 = t119 * t248 + t121 * t200 + t123 * t201;
t53 = t120 * t248 + t122 * t200 + t124 * t201;
t70 = t151 * t248 + t152 * t200 + t153 * t201;
t8 = t246 * t52 + t248 * t53 + t263 * t70;
t341 = t246 * t7 + t248 * t8 + t263 * t26;
t204 = -t247 * t303 + t276 * t306;
t205 = t247 * t306 + t366;
t127 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t246;
t129 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t246;
t131 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t246;
t64 = t127 * t263 + t129 * t244 + t131 * t245;
t72 = t156 * t246 + t157 * t204 + t158 * t205;
t340 = t64 / 0.2e1 + t72 / 0.2e1;
t128 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t248;
t130 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t248;
t132 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t248;
t65 = t128 * t263 + t130 * t244 + t132 * t245;
t73 = t156 * t248 + t157 * t206 + t158 * t207;
t339 = t65 / 0.2e1 + t73 / 0.2e1;
t335 = -t192 - t354;
t334 = t301 * t192 + t351;
t333 = -t191 + t350;
t332 = -t218 + t345;
t170 = t249 * rSges(5,1) - t248 * rSges(5,2) + t278 * rSges(5,3);
t254 = -t279 * t298 + t300 * t362;
t255 = t279 * t300 + t338;
t182 = t255 * rSges(4,1) + t254 * rSges(4,2) + t278 * rSges(4,3);
t227 = t279 * rSges(3,1) - t278 * rSges(3,2) + rSges(3,3) * t362;
t328 = -t361 / 0.2e1;
t327 = -t305 * pkin(1) + pkin(8) * t360;
t325 = t299 * (-rSges(4,1) * t275 - rSges(4,2) * t274 + rSges(4,3) * t361 - t280);
t324 = t171 * t362 + t172 * t360 + t344;
t323 = t369 / 0.2e1 + t368 / 0.2e1 + t69 * t377 + t70 * t376 + t74;
t213 = rSges(5,1) * t264 - rSges(5,2) * t263 - rSges(5,3) * t361;
t322 = t299 * (-t213 + t345);
t13 = t276 * t50 + t278 * t51 - t361 * t69;
t14 = t276 * t52 + t278 * t53 - t361 * t70;
t28 = t62 * t276 + t63 * t278 - t361 * t78;
t321 = t13 * t377 + t14 * t376 + t26 * t328 + t28 * t375 + t8 * t373 + t7 * t374;
t15 = t301 * t69 + (t305 * t51 - t308 * t50) * t299;
t16 = t301 * t70 + (t305 * t53 - t308 * t52) * t299;
t75 = t78 * t301;
t30 = t75 + (t63 * t305 - t62 * t308) * t299;
t320 = t15 * t377 + t16 * t376 + t26 * t372 + t30 * t375 + t8 * t362 / 0.2e1 - t7 * t360 / 0.2e1;
t319 = -rSges(5,1) * t247 + rSges(5,2) * t246;
t316 = t329 + t343;
t159 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t263;
t315 = t299 * (-t159 + t332);
t314 = t191 * t362 + t192 * t360 + t324;
t313 = t299 * (t332 - t352);
t312 = -t277 * t291 + t283 + t327;
t133 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t246;
t226 = rSges(3,1) * t277 - rSges(3,2) * t276 - rSges(3,3) * t360;
t210 = Icges(5,5) * t264 - Icges(5,6) * t263 - Icges(5,3) * t361;
t103 = t210 * t276 - t211 * t246 + t212 * t247;
t163 = Icges(5,5) * t247 - Icges(5,6) * t246 + Icges(5,3) * t276;
t165 = Icges(5,4) * t247 - Icges(5,2) * t246 + Icges(5,6) * t276;
t167 = Icges(5,1) * t247 - Icges(5,4) * t246 + Icges(5,5) * t276;
t94 = -t163 * t361 - t165 * t263 + t167 * t264;
t311 = t62 / 0.2e1 + t69 / 0.2e1 + t103 / 0.2e1 + t94 / 0.2e1 + t340;
t104 = t210 * t278 - t211 * t248 + t212 * t249;
t164 = Icges(5,5) * t249 - Icges(5,6) * t248 + Icges(5,3) * t278;
t166 = Icges(5,4) * t249 - Icges(5,2) * t248 + Icges(5,6) * t278;
t168 = Icges(5,1) * t249 - Icges(5,4) * t248 + Icges(5,5) * t278;
t95 = -t164 * t361 - t166 * t263 + t168 * t264;
t310 = t63 / 0.2e1 + t70 / 0.2e1 + t95 / 0.2e1 + t104 / 0.2e1 + t339;
t285 = rSges(2,1) * t308 - rSges(2,2) * t305;
t284 = -rSges(2,1) * t305 - rSges(2,2) * t308;
t265 = rSges(3,3) * t301 + (rSges(3,1) * t304 + rSges(3,2) * t307) * t299;
t225 = Icges(3,1) * t279 - Icges(3,4) * t278 + Icges(3,5) * t362;
t224 = Icges(3,1) * t277 - Icges(3,4) * t276 - Icges(3,5) * t360;
t223 = Icges(3,4) * t279 - Icges(3,2) * t278 + Icges(3,6) * t362;
t222 = Icges(3,4) * t277 - Icges(3,2) * t276 - Icges(3,6) * t360;
t221 = Icges(3,5) * t279 - Icges(3,6) * t278 + Icges(3,3) * t362;
t220 = Icges(3,5) * t277 - Icges(3,6) * t276 - Icges(3,3) * t360;
t209 = t227 + t343;
t208 = -t226 + t327;
t190 = -t226 * t301 - t265 * t360;
t189 = t227 * t301 - t265 * t362;
t180 = Icges(4,1) * t255 + Icges(4,4) * t254 + Icges(4,5) * t278;
t179 = Icges(4,1) * t253 + Icges(4,4) * t252 + Icges(4,5) * t276;
t178 = Icges(4,4) * t255 + Icges(4,2) * t254 + Icges(4,6) * t278;
t177 = Icges(4,4) * t253 + Icges(4,2) * t252 + Icges(4,6) * t276;
t176 = Icges(4,5) * t255 + Icges(4,6) * t254 + Icges(4,3) * t278;
t175 = Icges(4,5) * t253 + Icges(4,6) * t252 + Icges(4,3) * t276;
t173 = t278 * t191;
t169 = rSges(5,3) * t276 - t319;
t155 = (t226 * t305 + t227 * t308) * t299;
t150 = t260 * t362 - t261 * t278 + t262 * t279;
t149 = -t260 * t360 - t261 * t276 + t262 * t277;
t145 = t237 + t182 + t343;
t144 = t327 + t349;
t141 = t246 * t154;
t138 = t221 * t301 + (t223 * t307 + t225 * t304) * t299;
t137 = t220 * t301 + (t222 * t307 + t224 * t304) * t299;
t136 = t316 + t170;
t135 = (-rSges(5,3) + t302) * t276 + t312 + t319;
t118 = -t170 * t361 - t213 * t278;
t117 = t169 * t361 + t213 * t276;
t114 = t263 * t126;
t112 = t248 * t125;
t111 = t301 * t349 + t308 * t325;
t110 = t182 * t301 + t305 * t325 + t232;
t109 = -t210 * t361 + t347;
t108 = t169 * t278 - t170 * t276;
t107 = t109 * t301;
t106 = t214 * t278 + t215 * t254 + t216 * t255;
t105 = t214 * t276 + t215 * t252 + t216 * t253;
t102 = (t181 * t305 + t182 * t308) * t299 + t344;
t101 = -t176 * t361 + t178 * t274 + t180 * t275;
t100 = -t175 * t361 + t177 * t274 + t179 * t275;
t99 = t316 - t353;
t98 = -t133 - t191 + t312 + t367;
t97 = t134 * t263 - t159 * t248;
t96 = -t133 * t263 + t159 * t246;
t93 = -t154 * t248 + t114;
t92 = -t125 * t263 + t141;
t91 = t316 + t331 + t126;
t90 = -t247 * t292 + (-pkin(5) * t303 + t302) * t276 + (-rSges(7,3) + t309) * t246 + t312 + t318;
t89 = t164 * t278 - t166 * t248 + t168 * t249;
t88 = t163 * t278 - t165 * t248 + t167 * t249;
t87 = t164 * t276 - t166 * t246 + t168 * t247;
t86 = t163 * t276 - t165 * t246 + t167 * t247;
t85 = (-t169 + t350) * t301 + t308 * t322;
t84 = t170 * t301 + t305 * t322 + t351;
t83 = t133 * t248 - t134 * t246;
t81 = t82 * t301;
t80 = -t126 * t246 + t112;
t79 = t82 * t263;
t77 = t353 * t361 + (-t159 - t218) * t278;
t76 = t133 * t361 + t159 * t276 + t348;
t71 = (t169 * t305 + t170 * t308) * t299 + t324;
t66 = t133 * t278 + t276 * t353 + t173;
t61 = (-t133 + t333) * t301 + t308 * t315;
t60 = t134 * t301 + t305 * t315 + t334;
t59 = t128 * t248 + t130 * t206 + t132 * t207;
t58 = t127 * t248 + t129 * t206 + t131 * t207;
t57 = t128 * t246 + t130 * t204 + t132 * t205;
t56 = t127 * t246 + t129 * t204 + t131 * t205;
t49 = t116 * t263 - t248 * t352 + t114;
t48 = t148 * t246 - t263 * t355 + t141;
t47 = t335 * t361 + (-t218 - t352) * t278;
t46 = t276 * t352 + t355 * t361 + t348;
t45 = (t133 * t305 + t134 * t308) * t299 + t314;
t44 = t115 * t248 - t246 * t354 + t112;
t43 = t107 + (t95 * t305 - t94 * t308) * t299;
t42 = -t109 * t361 + t94 * t276 + t95 * t278;
t41 = (t333 - t355) * t301 + t308 * t313;
t40 = t301 * t354 + t305 * t313 + t334;
t39 = t276 * t335 + t278 * t355 + t173;
t38 = t104 * t301 + (t305 * t89 - t308 * t88) * t299;
t37 = t103 * t301 + (t305 * t87 - t308 * t86) * t299;
t36 = -t104 * t361 + t276 * t88 + t278 * t89;
t35 = -t103 * t361 + t276 * t86 + t278 * t87;
t34 = (t305 * t355 + t308 * t354) * t299 + t314;
t33 = t81 + (t65 * t305 - t64 * t308) * t299;
t32 = t64 * t276 + t65 * t278 - t361 * t82;
t31 = t64 * t246 + t65 * t248 + t79;
t23 = t301 * t73 + (t305 * t59 - t308 * t58) * t299;
t22 = t301 * t72 + (t305 * t57 - t308 * t56) * t299;
t20 = t276 * t58 + t278 * t59 - t361 * t73;
t19 = t276 * t56 + t278 * t57 - t361 * t72;
t18 = t246 * t58 + t248 * t59 + t263 * t73;
t17 = t246 * t56 + t248 * t57 + t263 * t72;
t1 = [m(7) * (t90 ^ 2 + t91 ^ 2) + m(6) * (t98 ^ 2 + t99 ^ 2) + m(5) * (t135 ^ 2 + t136 ^ 2) + m(4) * (t144 ^ 2 + t145 ^ 2) + m(3) * (t208 ^ 2 + t209 ^ 2) + m(2) * (t284 ^ 2 + t285 ^ 2) + (-t210 - t214) * t361 + Icges(2,3) + t347 - t379 + t380; t107 + t75 + t81 + m(7) * (t40 * t91 + t41 * t90) + m(6) * (t60 * t99 + t61 * t98) + m(5) * (t135 * t85 + t136 * t84) + m(4) * (t110 * t145 + t111 * t144) + m(3) * (t189 * t209 + t190 * t208) + ((-t137 / 0.2e1 - t100 / 0.2e1 - t105 / 0.2e1 - t149 / 0.2e1 - t311) * t308 + (t138 / 0.2e1 + t101 / 0.2e1 + t106 / 0.2e1 + t150 / 0.2e1 + t310) * t305) * t299 + t378; m(7) * (t34 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t45 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t71 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(4) * (t102 ^ 2 + t110 ^ 2 + t111 ^ 2) + m(3) * (t155 ^ 2 + t189 ^ 2 + t190 ^ 2) + (t16 + t23 + t38 + ((t176 * t278 + t178 * t254 + t180 * t255) * t305 - (t175 * t278 + t177 * t254 + t179 * t255) * t308) * t299 + (t221 * t362 - t223 * t278 + t225 * t279) * t362) * t362 + (-t15 - t22 - t37 - ((t176 * t276 + t178 * t252 + t180 * t253) * t305 - (t175 * t276 + t177 * t252 + t179 * t253) * t308) * t299 + (-t220 * t360 - t222 * t276 + t224 * t277) * t360 + (-t220 * t362 + t221 * t360 + t222 * t278 + t276 * t223 - t224 * t279 - t277 * t225) * t362) * t360 + (t30 + t33 + t43 + (t150 + t106) * t362 + (-t149 - t105) * t360 + ((-t100 - t137) * t308 + (t101 + t138) * t305) * t299 + t378) * t301; m(7) * (t276 * t91 + t278 * t90) + m(6) * (t276 * t99 + t278 * t98) + m(5) * (t135 * t278 + t136 * t276) + m(4) * (t144 * t278 + t145 * t276); m(7) * (t276 * t40 + t278 * t41 - t34 * t361) + m(6) * (t276 * t60 + t278 * t61 - t361 * t45) + m(5) * (t276 * t84 + t278 * t85 - t361 * t71) + m(4) * (-t102 * t361 + t110 * t276 + t111 * t278); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t299 ^ 2 * t307 ^ 2 + t276 ^ 2 + t278 ^ 2); (-t109 + t379) * t361 + m(7) * (t46 * t90 + t47 * t91) + m(6) * (t76 * t98 + t77 * t99) + m(5) * (t117 * t135 + t118 * t136) + t310 * t278 + t311 * t276; (t28 / 0.2e1 + t32 / 0.2e1 + t42 / 0.2e1) * t301 + (t16 / 0.2e1 + t23 / 0.2e1 + t38 / 0.2e1) * t278 + (t15 / 0.2e1 + t22 / 0.2e1 + t37 / 0.2e1) * t276 + m(7) * (t34 * t39 + t40 * t47 + t41 * t46) + m(6) * (t45 * t66 + t60 * t77 + t61 * t76) + m(5) * (t108 * t71 + t117 * t85 + t118 * t84) + ((-t13 / 0.2e1 - t19 / 0.2e1 - t35 / 0.2e1) * t308 + (-t30 / 0.2e1 - t33 / 0.2e1 - t43 / 0.2e1) * t307 + (t14 / 0.2e1 + t20 / 0.2e1 + t36 / 0.2e1) * t305) * t299; m(5) * (-t108 * t361 + t117 * t278 + t118 * t276) + m(6) * (t276 * t77 + t278 * t76 - t361 * t66) + m(7) * (t276 * t47 + t278 * t46 - t361 * t39); (-t28 - t32 - t42) * t361 + (t14 + t20 + t36) * t278 + (t13 + t19 + t35) * t276 + m(7) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t66 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t108 ^ 2 + t117 ^ 2 + t118 ^ 2); t79 + t339 * t248 + t340 * t246 + m(7) * (t48 * t90 + t49 * t91) + m(6) * (t96 * t98 + t97 * t99) + t323; t31 * t372 + t23 * t376 + t22 * t377 + t33 * t375 + (t305 * t18 / 0.2e1 - t308 * t17 / 0.2e1) * t299 + m(7) * (t34 * t44 + t40 * t49 + t41 * t48) + m(6) * (t45 * t83 + t60 * t97 + t61 * t96) + t320; m(6) * (t276 * t97 + t278 * t96 - t361 * t83) + m(7) * (t276 * t49 + t278 * t48 - t361 * t44); t31 * t328 + t18 * t373 + t20 * t376 + t32 * t375 + t17 * t374 + t19 * t377 + m(7) * (t39 * t44 + t46 * t48 + t47 * t49) + m(6) * (t66 * t83 + t76 * t96 + t77 * t97) + t321; t246 * t17 + t248 * t18 + t263 * t31 + m(7) * (t44 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t83 ^ 2 + t96 ^ 2 + t97 ^ 2) + t341; m(7) * (t90 * t92 + t91 * t93) + t323; m(7) * (t34 * t80 + t40 * t93 + t41 * t92) + t320; m(7) * (t276 * t93 + t278 * t92 - t361 * t80); m(7) * (t39 * t80 + t46 * t92 + t47 * t93) + t321; m(7) * (t44 * t80 + t48 * t92 + t49 * t93) + t341; m(7) * (t80 ^ 2 + t92 ^ 2 + t93 ^ 2) + t341;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
