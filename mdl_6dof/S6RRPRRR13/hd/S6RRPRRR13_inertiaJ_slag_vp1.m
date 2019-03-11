% Calculate joint inertia matrix for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR13_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:45:11
% EndTime: 2019-03-09 14:45:25
% DurationCPUTime: 6.43s
% Computational Cost: add. (22175->651), mult. (48341->889), div. (0->0), fcn. (61444->12), ass. (0->297)
t292 = sin(pkin(6));
t293 = cos(pkin(6));
t296 = sin(qJ(2));
t299 = cos(qJ(2));
t251 = Icges(3,3) * t293 + (Icges(3,5) * t296 + Icges(3,6) * t299) * t292;
t252 = Icges(3,6) * t293 + (Icges(3,4) * t296 + Icges(3,2) * t299) * t292;
t253 = Icges(3,5) * t293 + (Icges(3,1) * t296 + Icges(3,4) * t299) * t292;
t254 = Icges(4,5) * t293 + (-Icges(4,6) * t296 - Icges(4,3) * t299) * t292;
t255 = Icges(4,4) * t293 + (-Icges(4,2) * t296 - Icges(4,6) * t299) * t292;
t256 = Icges(4,1) * t293 + (-Icges(4,4) * t296 - Icges(4,5) * t299) * t292;
t351 = t292 * t299;
t353 = t292 * t296;
t372 = (-t299 * t254 - t296 * t255) * t292 + t252 * t351 + t253 * t353 + (t256 + t251) * t293;
t297 = sin(qJ(1));
t347 = t297 * t299;
t300 = cos(qJ(1));
t348 = t296 * t300;
t273 = t293 * t347 + t348;
t346 = t299 * t300;
t349 = t296 * t297;
t274 = -t293 * t349 + t346;
t352 = t292 * t297;
t203 = Icges(4,5) * t352 - Icges(4,6) * t274 + Icges(4,3) * t273;
t212 = Icges(3,4) * t274 - Icges(3,2) * t273 + Icges(3,6) * t352;
t371 = t203 - t212;
t271 = -t293 * t346 + t349;
t272 = t293 * t348 + t347;
t350 = t292 * t300;
t204 = -Icges(4,5) * t350 - Icges(4,6) * t272 + Icges(4,3) * t271;
t211 = Icges(3,4) * t272 - Icges(3,2) * t271 - Icges(3,6) * t350;
t370 = t204 - t211;
t205 = Icges(4,4) * t352 - Icges(4,2) * t274 + Icges(4,6) * t273;
t214 = Icges(3,1) * t274 - Icges(3,4) * t273 + Icges(3,5) * t352;
t369 = t205 - t214;
t206 = -Icges(4,4) * t350 - Icges(4,2) * t272 + Icges(4,6) * t271;
t213 = Icges(3,1) * t272 - Icges(3,4) * t271 - Icges(3,5) * t350;
t368 = t206 - t213;
t367 = t292 ^ 2;
t295 = sin(qJ(4));
t360 = cos(qJ(4));
t241 = -t273 * t360 + t295 * t352;
t366 = t241 / 0.2e1;
t243 = t271 * t360 + t295 * t350;
t365 = -t243 / 0.2e1;
t323 = t292 * t360;
t269 = t293 * t295 + t299 * t323;
t364 = t269 / 0.2e1;
t363 = t272 / 0.2e1;
t362 = t274 / 0.2e1;
t361 = t293 / 0.2e1;
t359 = t272 * pkin(2);
t298 = cos(qJ(5));
t287 = pkin(5) * t298 + pkin(4);
t358 = -pkin(4) + t287;
t242 = t273 * t295 + t297 * t323;
t291 = qJ(5) + qJ(6);
t288 = sin(t291);
t289 = cos(t291);
t185 = -t242 * t288 + t274 * t289;
t186 = t242 * t289 + t274 * t288;
t119 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t241;
t121 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t241;
t123 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t241;
t270 = t293 * t360 - t295 * t351;
t227 = -t270 * t288 + t289 * t353;
t228 = t270 * t289 + t288 * t353;
t60 = t119 * t269 + t121 * t227 + t123 * t228;
t357 = t60 * t241;
t244 = t271 * t295 - t300 * t323;
t187 = -t244 * t288 + t272 * t289;
t188 = t244 * t289 + t272 * t288;
t120 = Icges(7,5) * t188 + Icges(7,6) * t187 - Icges(7,3) * t243;
t122 = Icges(7,4) * t188 + Icges(7,2) * t187 - Icges(7,6) * t243;
t124 = Icges(7,1) * t188 + Icges(7,4) * t187 - Icges(7,5) * t243;
t61 = t120 * t269 + t122 * t227 + t124 * t228;
t356 = t61 * t243;
t294 = sin(qJ(5));
t355 = t272 * t294;
t354 = t274 * t294;
t179 = t242 * pkin(4) + pkin(10) * t241;
t301 = -pkin(11) - pkin(10);
t325 = pkin(5) * t354 - t241 * t301 + t242 * t287;
t113 = -t179 + t325;
t127 = t186 * rSges(7,1) + t185 * rSges(7,2) + t241 * rSges(7,3);
t345 = t113 + t127;
t237 = t243 * pkin(10);
t114 = pkin(5) * t355 + t243 * t301 + t244 * t358 + t237;
t312 = -t188 * rSges(7,1) - t187 * rSges(7,2);
t128 = -t243 * rSges(7,3) - t312;
t344 = t114 + t128;
t194 = -t242 * t294 + t274 * t298;
t195 = t242 * t298 + t354;
t135 = t195 * rSges(6,1) + t194 * rSges(6,2) + t241 * rSges(6,3);
t343 = -t135 - t179;
t196 = -t244 * t294 + t272 * t298;
t197 = t244 * t298 + t355;
t136 = rSges(6,1) * t197 + rSges(6,2) * t196 - rSges(6,3) * t243;
t180 = t244 * pkin(4) - t237;
t342 = -t136 - t180;
t156 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t269;
t329 = t294 * t353;
t157 = pkin(5) * t329 + t358 * t270 + (-pkin(10) - t301) * t269;
t341 = t156 + t157;
t340 = t372 * t293;
t207 = Icges(4,1) * t352 - Icges(4,4) * t274 + Icges(4,5) * t273;
t210 = Icges(3,5) * t274 - Icges(3,6) * t273 + Icges(3,3) * t352;
t339 = t207 + t210;
t208 = -Icges(4,1) * t350 - Icges(4,4) * t272 + Icges(4,5) * t271;
t209 = Icges(3,5) * t272 - Icges(3,6) * t271 - Icges(3,3) * t350;
t338 = -t208 - t209;
t261 = t271 * qJ(3);
t224 = t261 + t359;
t225 = t274 * pkin(2) + qJ(3) * t273;
t337 = t224 * t352 + t225 * t350;
t222 = t293 * t225;
t247 = pkin(3) * t352 + pkin(9) * t274;
t336 = t293 * t247 + t222;
t248 = -pkin(3) * t350 + t272 * pkin(9);
t335 = -t224 - t248;
t275 = (pkin(2) * t296 - qJ(3) * t299) * t292;
t334 = -pkin(3) * t293 - pkin(9) * t353 - t275;
t333 = t300 * pkin(1) + pkin(8) * t352;
t153 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t269;
t154 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t269;
t155 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t269;
t79 = t269 * t153 + t227 * t154 + t228 * t155;
t73 = t79 * t269;
t26 = -t356 + t73 + t357;
t47 = t119 * t241 + t121 * t185 + t123 * t186;
t48 = t120 * t241 + t122 * t185 + t124 * t186;
t69 = t153 * t241 + t154 * t185 + t155 * t186;
t7 = t241 * t47 - t243 * t48 + t269 * t69;
t49 = -t119 * t243 + t121 * t187 + t123 * t188;
t50 = -t120 * t243 + t122 * t187 + t124 * t188;
t70 = -t153 * t243 + t154 * t187 + t155 * t188;
t8 = t241 * t49 - t243 * t50 + t269 * t70;
t332 = t241 * t7 - t243 * t8 + t269 * t26;
t129 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t241;
t131 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t241;
t133 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t241;
t238 = -t270 * t294 + t298 * t353;
t239 = t270 * t298 + t329;
t62 = t129 * t269 + t131 * t238 + t133 * t239;
t160 = Icges(6,5) * t239 + Icges(6,6) * t238 + Icges(6,3) * t269;
t161 = Icges(6,4) * t239 + Icges(6,2) * t238 + Icges(6,6) * t269;
t162 = Icges(6,1) * t239 + Icges(6,4) * t238 + Icges(6,5) * t269;
t71 = t160 * t241 + t161 * t194 + t162 * t195;
t331 = t62 / 0.2e1 + t71 / 0.2e1;
t130 = Icges(6,5) * t197 + Icges(6,6) * t196 - Icges(6,3) * t243;
t132 = Icges(6,4) * t197 + Icges(6,2) * t196 - Icges(6,6) * t243;
t134 = Icges(6,1) * t197 + Icges(6,4) * t196 - Icges(6,5) * t243;
t63 = t130 * t269 + t132 * t238 + t134 * t239;
t72 = -t160 * t243 + t161 * t196 + t162 * t197;
t330 = -t72 / 0.2e1 - t63 / 0.2e1;
t84 = t269 * t160 + t238 * t161 + t239 * t162;
t328 = t293 * t179 + t336;
t327 = -t180 + t335;
t199 = Icges(5,5) * t270 - Icges(5,6) * t269 + Icges(5,3) * t353;
t200 = Icges(5,4) * t270 - Icges(5,2) * t269 + Icges(5,6) * t353;
t201 = Icges(5,1) * t270 - Icges(5,4) * t269 + Icges(5,5) * t353;
t109 = t199 * t353 - t269 * t200 + t270 * t201;
t223 = t270 * pkin(4) + t269 * pkin(10);
t326 = -t223 + t334;
t170 = t242 * rSges(5,1) - t241 * rSges(5,2) + t274 * rSges(5,3);
t218 = t274 * rSges(3,1) - t273 * rSges(3,2) + rSges(3,3) * t352;
t215 = rSges(4,1) * t352 - t274 * rSges(4,2) + t273 * rSges(4,3);
t322 = t353 / 0.2e1;
t321 = -t297 * pkin(1) + pkin(8) * t350;
t320 = t292 * (-rSges(4,1) * t293 - (-rSges(4,2) * t296 - rSges(4,3) * t299) * t292 - t275);
t319 = t247 * t350 + t248 * t352 + t337;
t318 = t357 / 0.2e1 - t356 / 0.2e1 + t69 * t366 + t70 * t365 + t73;
t317 = -t261 + t321;
t202 = rSges(5,1) * t270 - rSges(5,2) * t269 + rSges(5,3) * t353;
t316 = t292 * (-t202 + t334);
t13 = t272 * t48 + t274 * t47 + t353 * t69;
t14 = t272 * t50 + t274 * t49 + t353 * t70;
t76 = t79 * t353;
t28 = t61 * t272 + t60 * t274 + t76;
t315 = t13 * t366 + t14 * t365 + t26 * t322 + t28 * t364 + t7 * t362 + t8 * t363;
t15 = t293 * t69 + (t297 * t47 - t300 * t48) * t292;
t16 = t293 * t70 + (t297 * t49 - t300 * t50) * t292;
t77 = t79 * t293;
t30 = t77 + (t60 * t297 - t61 * t300) * t292;
t314 = t15 * t366 + t16 * t365 + t26 * t361 + t30 * t364 + t7 * t352 / 0.2e1 - t8 * t350 / 0.2e1;
t313 = -rSges(5,1) * t244 - rSges(5,2) * t243;
t163 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t269;
t311 = t292 * (-t163 + t326);
t310 = t225 + t333;
t309 = rSges(4,1) * t350 - rSges(4,3) * t271;
t308 = t179 * t350 + t180 * t352 + t319;
t307 = t317 - t248;
t306 = t292 * (t326 - t341);
t217 = rSges(3,1) * t272 - rSges(3,2) * t271 - rSges(3,3) * t350;
t103 = t199 * t272 + t200 * t243 + t201 * t244;
t165 = Icges(5,5) * t244 + Icges(5,6) * t243 + Icges(5,3) * t272;
t167 = Icges(5,4) * t244 + Icges(5,2) * t243 + Icges(5,6) * t272;
t169 = Icges(5,1) * t244 + Icges(5,4) * t243 + Icges(5,5) * t272;
t95 = t165 * t353 - t167 * t269 + t169 * t270;
t304 = t61 / 0.2e1 + t70 / 0.2e1 + t95 / 0.2e1 + t103 / 0.2e1 - t330;
t102 = t199 * t274 - t200 * t241 + t201 * t242;
t164 = Icges(5,5) * t242 - Icges(5,6) * t241 + Icges(5,3) * t274;
t166 = Icges(5,4) * t242 - Icges(5,2) * t241 + Icges(5,6) * t274;
t168 = Icges(5,1) * t242 - Icges(5,4) * t241 + Icges(5,5) * t274;
t94 = t164 * t353 - t166 * t269 + t168 * t270;
t303 = t69 / 0.2e1 + t60 / 0.2e1 + t102 / 0.2e1 + t94 / 0.2e1 + t331;
t302 = t247 + t310;
t278 = rSges(2,1) * t300 - rSges(2,2) * t297;
t277 = -rSges(2,1) * t297 - rSges(2,2) * t300;
t257 = rSges(3,3) * t293 + (rSges(3,1) * t296 + rSges(3,2) * t299) * t292;
t216 = -rSges(4,2) * t272 - t309;
t193 = t272 * t223;
t190 = t218 + t333;
t189 = -t217 + t321;
t176 = t179 * t353;
t174 = -t217 * t293 - t257 * t350;
t173 = t218 * t293 - t257 * t352;
t172 = t274 * t180;
t171 = rSges(5,3) * t272 - t313;
t152 = t310 + t215;
t151 = (rSges(4,2) - pkin(2)) * t272 + t309 + t317;
t149 = (t217 * t297 + t218 * t300) * t292;
t148 = t251 * t352 - t252 * t273 + t253 * t274;
t147 = -t251 * t350 - t252 * t271 + t253 * t272;
t146 = t254 * t271 - t255 * t272 - t256 * t350;
t145 = t254 * t273 - t255 * t274 + t256 * t352;
t141 = t243 * t156;
t138 = (-t216 - t224) * t293 + t300 * t320;
t137 = t215 * t293 + t297 * t320 + t222;
t126 = t170 * t353 - t202 * t274;
t125 = -t171 * t353 + t202 * t272;
t118 = t208 * t293 + (-t204 * t299 - t206 * t296) * t292;
t117 = t207 * t293 + (-t203 * t299 - t205 * t296) * t292;
t116 = t210 * t293 + (t212 * t299 + t214 * t296) * t292;
t115 = t209 * t293 + (t211 * t299 + t213 * t296) * t292;
t112 = t302 + t170;
t111 = (-rSges(5,3) - pkin(2)) * t272 + t307 + t313;
t110 = t269 * t127;
t108 = t241 * t128;
t107 = t109 * t293;
t106 = t109 * t353;
t105 = (t215 * t300 + t216 * t297) * t292 + t337;
t104 = -t170 * t272 + t171 * t274;
t101 = (-t171 + t335) * t293 + t300 * t316;
t100 = t170 * t293 + t297 * t316 + t336;
t99 = -t136 * t269 - t163 * t243;
t98 = t135 * t269 - t163 * t241;
t97 = t302 - t343;
t96 = t307 + t342 - t359;
t93 = -t128 * t269 - t141;
t92 = -t156 * t241 + t110;
t91 = t302 + t325 + t127;
t90 = -t244 * t287 + (-pkin(5) * t294 - pkin(2)) * t272 + (rSges(7,3) - t301) * t243 + t307 + t312;
t89 = t165 * t272 + t167 * t243 + t169 * t244;
t88 = t164 * t272 + t166 * t243 + t168 * t244;
t87 = t165 * t274 - t167 * t241 + t169 * t242;
t86 = t164 * t274 - t166 * t241 + t168 * t242;
t85 = (t170 * t300 + t171 * t297) * t292 + t319;
t83 = t84 * t293;
t82 = t84 * t353;
t81 = t84 * t269;
t80 = t135 * t243 + t136 * t241;
t78 = t127 * t243 + t108;
t75 = t135 * t353 + t176 + (-t163 - t223) * t274;
t74 = t163 * t272 + t342 * t353 + t193;
t66 = (-t136 + t327) * t293 + t300 * t311;
t65 = t135 * t293 + t297 * t311 + t328;
t64 = t136 * t274 + t272 * t343 + t172;
t59 = -t130 * t243 + t132 * t196 + t134 * t197;
t58 = -t129 * t243 + t131 * t196 + t133 * t197;
t57 = t130 * t241 + t132 * t194 + t134 * t195;
t56 = t129 * t241 + t131 * t194 + t133 * t195;
t53 = -t157 * t243 - t269 * t344 - t141;
t52 = t113 * t269 - t241 * t341 + t110;
t51 = (t135 * t300 + t136 * t297) * t292 + t308;
t46 = t176 + t345 * t353 + (-t223 - t341) * t274;
t45 = t193 + t341 * t272 + (-t180 - t344) * t353;
t44 = (t327 - t344) * t293 + t300 * t306;
t43 = t293 * t345 + t297 * t306 + t328;
t42 = t114 * t241 + t243 * t345 + t108;
t41 = t107 + (t94 * t297 - t95 * t300) * t292;
t40 = t95 * t272 + t94 * t274 + t106;
t39 = t172 + t344 * t274 + (-t179 - t345) * t272;
t38 = t103 * t293 + (t297 * t88 - t300 * t89) * t292;
t37 = t102 * t293 + (t297 * t86 - t300 * t87) * t292;
t36 = t103 * t353 + t272 * t89 + t274 * t88;
t35 = t102 * t353 + t272 * t87 + t274 * t86;
t34 = (t297 * t344 + t300 * t345) * t292 + t308;
t33 = t83 + (t62 * t297 - t63 * t300) * t292;
t32 = t63 * t272 + t62 * t274 + t82;
t31 = t62 * t241 - t63 * t243 + t81;
t22 = t293 * t72 + (t297 * t58 - t300 * t59) * t292;
t21 = t293 * t71 + (t297 * t56 - t300 * t57) * t292;
t20 = t272 * t59 + t274 * t58 + t353 * t72;
t19 = t272 * t57 + t274 * t56 + t353 * t71;
t18 = t241 * t58 - t243 * t59 + t269 * t72;
t17 = t241 * t56 - t243 * t57 + t269 * t71;
t1 = [Icges(2,3) + m(7) * (t90 ^ 2 + t91 ^ 2) + m(6) * (t96 ^ 2 + t97 ^ 2) + m(5) * (t111 ^ 2 + t112 ^ 2) + m(4) * (t151 ^ 2 + t152 ^ 2) + m(3) * (t189 ^ 2 + t190 ^ 2) + m(2) * (t277 ^ 2 + t278 ^ 2) + t109 + t84 + t79 + t372; t77 + t83 + t107 + m(7) * (t43 * t91 + t44 * t90) + m(6) * (t65 * t97 + t66 * t96) + m(5) * (t100 * t112 + t101 * t111) + m(4) * (t137 * t152 + t138 * t151) + m(3) * (t173 * t190 + t174 * t189) + ((-t118 / 0.2e1 - t115 / 0.2e1 - t146 / 0.2e1 - t147 / 0.2e1 - t304) * t300 + (t117 / 0.2e1 + t116 / 0.2e1 + t145 / 0.2e1 + t148 / 0.2e1 + t303) * t297) * t292 + t340; (t30 + t33 + t41 + t340) * t293 + m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t51 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t85 ^ 2) + m(4) * (t105 ^ 2 + t137 ^ 2 + t138 ^ 2) + m(3) * (t149 ^ 2 + t173 ^ 2 + t174 ^ 2) + ((-t16 - t22 - t38 + ((t271 * t370 - t272 * t368) * t292 + t338 * t367 * t300) * t300 + (-t115 - t118 - t146 - t147) * t293) * t300 + (t15 + t21 + t37 + ((t273 * t371 - t369 * t274) * t292 + t339 * t367 * t297) * t297 + (t145 + t148 + t117 + t116) * t293 + ((t297 * t338 + t300 * t339) * t292 + t368 * t274 - t370 * t273 + t369 * t272 - t371 * t271) * t350) * t297) * t292; m(7) * (t271 * t91 + t273 * t90) + m(6) * (t271 * t97 + t273 * t96) + m(5) * (t111 * t273 + t112 * t271) + m(4) * (t151 * t273 + t152 * t271); m(7) * (t271 * t43 + t273 * t44 - t34 * t351) + m(6) * (t271 * t65 + t273 * t66 - t351 * t51) + m(5) * (t100 * t271 + t101 * t273 - t351 * t85) + m(4) * (-t105 * t351 + t137 * t271 + t138 * t273); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t299 ^ 2 * t367 + t271 ^ 2 + t273 ^ 2); t76 + t82 + t106 + m(7) * (t45 * t90 + t46 * t91) + m(6) * (t74 * t96 + t75 * t97) + m(5) * (t111 * t125 + t112 * t126) + t303 * t274 + t304 * t272; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1) * t293 + (t15 / 0.2e1 + t21 / 0.2e1 + t37 / 0.2e1) * t274 + (t16 / 0.2e1 + t22 / 0.2e1 + t38 / 0.2e1) * t272 + m(7) * (t34 * t39 + t43 * t46 + t44 * t45) + m(6) * (t64 * t51 + t65 * t75 + t66 * t74) + m(5) * (t100 * t126 + t101 * t125 + t104 * t85) + ((-t14 / 0.2e1 - t20 / 0.2e1 - t36 / 0.2e1) * t300 + (t13 / 0.2e1 + t19 / 0.2e1 + t35 / 0.2e1) * t297 + (t30 / 0.2e1 + t33 / 0.2e1 + t41 / 0.2e1) * t296) * t292; m(5) * (-t104 * t351 + t125 * t273 + t126 * t271) + m(6) * (t271 * t75 + t273 * t74 - t351 * t64) + m(7) * (t271 * t46 + t273 * t45 - t351 * t39); (t28 + t32 + t40) * t353 + (t13 + t19 + t35) * t274 + (t14 + t20 + t36) * t272 + m(7) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t64 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(5) * (t104 ^ 2 + t125 ^ 2 + t126 ^ 2); t81 + t330 * t243 + t331 * t241 + m(7) * (t52 * t91 + t53 * t90) + m(6) * (t96 * t99 + t97 * t98) + t318; t21 * t366 + t22 * t365 + t33 * t364 + t31 * t361 + (-t300 * t18 / 0.2e1 + t297 * t17 / 0.2e1) * t292 + m(7) * (t34 * t42 + t43 * t52 + t44 * t53) + m(6) * (t51 * t80 + t65 * t98 + t66 * t99) + t314; m(6) * (t271 * t98 + t273 * t99 - t351 * t80) + m(7) * (t271 * t52 + t273 * t53 - t351 * t42); t31 * t322 + t32 * t364 + t17 * t362 + t20 * t365 + t19 * t366 + t18 * t363 + m(7) * (t39 * t42 + t45 * t53 + t46 * t52) + m(6) * (t64 * t80 + t74 * t99 + t75 * t98) + t315; t241 * t17 - t243 * t18 + t269 * t31 + m(7) * (t42 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t80 ^ 2 + t98 ^ 2 + t99 ^ 2) + t332; m(7) * (t90 * t93 + t91 * t92) + t318; m(7) * (t34 * t78 + t43 * t92 + t44 * t93) + t314; m(7) * (t271 * t92 + t273 * t93 - t351 * t78); m(7) * (t39 * t78 + t45 * t93 + t46 * t92) + t315; m(7) * (t42 * t78 + t52 * t92 + t53 * t93) + t332; m(7) * (t78 ^ 2 + t92 ^ 2 + t93 ^ 2) + t332;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
