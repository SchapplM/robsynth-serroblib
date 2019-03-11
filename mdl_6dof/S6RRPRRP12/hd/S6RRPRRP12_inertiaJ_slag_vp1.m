% Calculate joint inertia matrix for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:26
% EndTime: 2019-03-09 12:50:38
% DurationCPUTime: 5.25s
% Computational Cost: add. (7957->479), mult. (11887->673), div. (0->0), fcn. (12871->8), ass. (0->228)
t234 = qJ(4) + qJ(5);
t221 = sin(t234);
t222 = cos(t234);
t237 = sin(qJ(1));
t236 = sin(qJ(2));
t240 = cos(qJ(1));
t296 = t236 * t240;
t177 = t221 * t237 - t222 * t296;
t178 = t221 * t296 + t222 * t237;
t239 = cos(qJ(2));
t292 = t239 * t240;
t101 = Icges(7,4) * t178 + Icges(7,2) * t292 + Icges(7,6) * t177;
t99 = Icges(6,5) * t178 - Icges(6,6) * t177 + Icges(6,3) * t292;
t347 = t101 + t99;
t103 = Icges(6,4) * t178 - Icges(6,2) * t177 + Icges(6,6) * t292;
t97 = Icges(7,5) * t178 + Icges(7,6) * t292 + Icges(7,3) * t177;
t346 = t103 - t97;
t297 = t236 * t237;
t179 = t221 * t240 + t222 * t297;
t181 = t221 * t297 - t222 * t240;
t294 = t237 * t239;
t104 = Icges(6,4) * t181 + Icges(6,2) * t179 + Icges(6,6) * t294;
t98 = Icges(7,5) * t181 + Icges(7,6) * t294 - Icges(7,3) * t179;
t345 = t104 - t98;
t100 = Icges(6,5) * t181 + Icges(6,6) * t179 + Icges(6,3) * t294;
t102 = Icges(7,4) * t181 + Icges(7,2) * t294 - Icges(7,6) * t179;
t344 = t100 + t102;
t105 = Icges(7,1) * t178 + Icges(7,4) * t292 + Icges(7,5) * t177;
t107 = Icges(6,1) * t178 - Icges(6,4) * t177 + Icges(6,5) * t292;
t343 = t105 + t107;
t106 = Icges(7,1) * t181 + Icges(7,4) * t294 - Icges(7,5) * t179;
t108 = Icges(6,1) * t181 + Icges(6,4) * t179 + Icges(6,5) * t294;
t342 = t106 + t108;
t333 = rSges(7,3) + qJ(6);
t340 = rSges(7,1) + pkin(5);
t341 = t179 * t333 - t181 * t340;
t339 = -t346 * t177 + t343 * t178 + t292 * t347;
t338 = -t177 * t345 + t178 * t342 + t292 * t344;
t337 = t346 * t179 + t343 * t181 + t294 * t347;
t336 = t179 * t345 + t181 * t342 + t294 * t344;
t146 = Icges(7,6) * t236 + (-Icges(7,5) * t221 + Icges(7,3) * t222) * t239;
t148 = Icges(7,2) * t236 + (-Icges(7,4) * t221 + Icges(7,6) * t222) * t239;
t150 = Icges(7,4) * t236 + (-Icges(7,1) * t221 + Icges(7,5) * t222) * t239;
t69 = t146 * t177 + t148 * t292 + t150 * t178;
t147 = Icges(6,3) * t236 + (-Icges(6,5) * t221 - Icges(6,6) * t222) * t239;
t149 = Icges(6,6) * t236 + (-Icges(6,4) * t221 - Icges(6,2) * t222) * t239;
t151 = Icges(6,5) * t236 + (-Icges(6,1) * t221 - Icges(6,4) * t222) * t239;
t70 = t147 * t292 - t149 * t177 + t151 * t178;
t335 = t70 + t69;
t71 = -t146 * t179 + t148 * t294 + t150 * t181;
t72 = t147 * t294 + t149 * t179 + t151 * t181;
t334 = t71 + t72;
t298 = t222 * t239;
t332 = t146 * t298 + (t147 + t148) * t236;
t331 = -t222 * t149 + (-t150 - t151) * t221;
t330 = Icges(4,1) + Icges(3,3);
t329 = (-Icges(4,4) + Icges(3,5)) * t239 + (Icges(4,5) - Icges(3,6)) * t236;
t328 = -t237 / 0.2e1;
t311 = t237 / 0.2e1;
t310 = -t240 / 0.2e1;
t327 = t240 / 0.2e1;
t326 = (t237 * t338 + t240 * t339) * t239 + t335 * t236;
t325 = (t237 * t336 + t240 * t337) * t239 + t334 * t236;
t324 = t237 * t339 - t240 * t338;
t323 = t237 * t337 - t240 * t336;
t48 = t101 * t236 + (-t105 * t221 + t222 * t97) * t239;
t50 = t236 * t99 + (-t103 * t222 - t107 * t221) * t239;
t322 = t48 + t50;
t49 = t102 * t236 + (-t106 * t221 + t222 * t98) * t239;
t51 = t100 * t236 + (-t104 * t222 - t108 * t221) * t239;
t321 = t49 + t51;
t320 = (t239 * t331 + t332) * t236;
t290 = rSges(7,2) * t292 + t177 * t333 + t178 * t340;
t289 = rSges(7,2) * t294 - t341;
t319 = rSges(7,2) * t236 + (-t221 * t340 + t222 * t333) * t239;
t312 = t236 / 0.2e1;
t318 = t237 * t330 + t240 * t329;
t317 = -t237 * t329 + t240 * t330;
t231 = t237 ^ 2;
t233 = t240 ^ 2;
t316 = m(4) / 0.2e1;
t315 = m(5) / 0.2e1;
t314 = m(6) / 0.2e1;
t313 = m(7) / 0.2e1;
t235 = sin(qJ(4));
t309 = pkin(4) * t235;
t308 = t240 * rSges(4,1);
t307 = t240 * rSges(3,3);
t306 = t289 * t292;
t305 = t290 * t236;
t304 = Icges(3,4) * t236;
t303 = Icges(3,4) * t239;
t302 = Icges(4,6) * t236;
t301 = Icges(4,6) * t239;
t300 = qJ(3) * t236;
t238 = cos(qJ(4));
t295 = t237 * t238;
t293 = t238 * t240;
t241 = -pkin(9) - pkin(8);
t291 = t239 * t241;
t287 = t319 * t294;
t282 = pkin(2) * t292 + qJ(3) * t296;
t285 = t231 * (pkin(2) * t239 + t300) + t240 * t282;
t203 = pkin(2) * t236 - qJ(3) * t239;
t284 = rSges(4,2) * t236 + rSges(4,3) * t239 - t203;
t220 = pkin(4) * t238 + pkin(3);
t283 = -t220 * t240 - t237 * t291;
t281 = pkin(1) * t240 + pkin(7) * t237;
t280 = pkin(3) * t237 + pkin(8) * t292;
t279 = t231 + t233;
t189 = -t235 * t237 + t236 * t293;
t276 = t235 * t296;
t190 = t276 + t295;
t121 = Icges(5,5) * t190 + Icges(5,6) * t189 + Icges(5,3) * t292;
t123 = Icges(5,4) * t190 + Icges(5,2) * t189 + Icges(5,6) * t292;
t125 = Icges(5,1) * t190 + Icges(5,4) * t189 + Icges(5,5) * t292;
t61 = t121 * t236 + (-t123 * t238 - t125 * t235) * t239;
t161 = Icges(5,3) * t236 + (-Icges(5,5) * t235 - Icges(5,6) * t238) * t239;
t164 = Icges(5,6) * t236 + (-Icges(5,4) * t235 - Icges(5,2) * t238) * t239;
t167 = Icges(5,5) * t236 + (-Icges(5,1) * t235 - Icges(5,4) * t238) * t239;
t76 = t161 * t292 + t164 * t189 + t167 * t190;
t278 = t76 / 0.2e1 + t61 / 0.2e1;
t191 = t235 * t240 + t236 * t295;
t192 = t235 * t297 - t293;
t122 = Icges(5,5) * t192 + Icges(5,6) * t191 + Icges(5,3) * t294;
t124 = Icges(5,4) * t192 + Icges(5,2) * t191 + Icges(5,6) * t294;
t126 = Icges(5,1) * t192 + Icges(5,4) * t191 + Icges(5,5) * t294;
t62 = t122 * t236 + (-t124 * t238 - t126 * t235) * t239;
t77 = t161 * t294 + t164 * t191 + t167 * t192;
t277 = t77 / 0.2e1 + t62 / 0.2e1;
t110 = rSges(6,1) * t178 - rSges(6,2) * t177 + rSges(6,3) * t292;
t130 = rSges(5,1) * t190 + rSges(5,2) * t189 + rSges(5,3) * t292;
t227 = t240 * pkin(7);
t275 = t227 - t283;
t274 = t294 / 0.2e1;
t273 = t292 / 0.2e1;
t272 = -Icges(4,4) * t236 / 0.2e1 + Icges(3,5) * t312 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t239;
t271 = -pkin(8) * t236 - t203;
t228 = t240 * pkin(3);
t270 = t237 * (pkin(8) * t294 - t228) + t240 * t280 + t285;
t269 = t281 + t282;
t176 = rSges(5,3) * t236 + (-rSges(5,1) * t235 - rSges(5,2) * t238) * t239;
t268 = -t176 + t271;
t185 = -t239 * t309 + (-pkin(8) - t241) * t236;
t267 = -t185 + t271;
t266 = t325 * t294 + t326 * t292 + ((t237 * t321 + t240 * t322) * t239 + t320) * t236;
t265 = rSges(3,1) * t239 - rSges(3,2) * t236;
t264 = -rSges(5,1) * t192 - rSges(5,2) * t191;
t263 = -rSges(6,1) * t181 - rSges(6,2) * t179;
t262 = Icges(3,1) * t239 - t304;
t261 = -Icges(3,2) * t236 + t303;
t258 = -Icges(4,2) * t239 + t302;
t257 = Icges(4,3) * t236 - t301;
t256 = -t238 * t164 - t235 * t167;
t154 = rSges(6,3) * t236 + (-rSges(6,1) * t221 - rSges(6,2) * t222) * t239;
t251 = -t154 + t267;
t250 = rSges(3,1) * t292 - rSges(3,2) * t296 + rSges(3,3) * t237;
t249 = rSges(4,1) * t237 - rSges(4,2) * t292 + rSges(4,3) * t296;
t248 = pkin(4) * t276 + t220 * t237 - t240 * t291;
t134 = t248 - t280;
t135 = t228 + (-pkin(8) * t239 + t236 * t309) * t237 + t283;
t247 = t134 * t240 + t135 * t237 + t270;
t246 = (-qJ(3) - t309) * t236 - pkin(1);
t245 = t267 - t319;
t244 = (t237 * t322 - t240 * t321) * t312 + t326 * t311 + t325 * t310 + t323 * t274 + t324 * t273;
t243 = (t321 + t334) * t274 + (t322 + t335) * t273 + t320;
t242 = t248 + t269;
t232 = t239 ^ 2;
t207 = rSges(2,1) * t240 - rSges(2,2) * t237;
t206 = -rSges(2,1) * t237 - rSges(2,2) * t240;
t205 = rSges(3,1) * t236 + rSges(3,2) * t239;
t155 = t236 * t161;
t152 = t185 * t294;
t144 = t284 * t240;
t143 = t284 * t237;
t140 = t250 + t281;
t139 = t307 + t227 + (-pkin(1) - t265) * t237;
t138 = t154 * t294;
t133 = t249 + t269;
t132 = t308 + t227 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t239 + (-rSges(4,3) - qJ(3)) * t236) * t237;
t131 = rSges(5,3) * t294 - t264;
t127 = t236 * t134;
t118 = t268 * t240;
t117 = t268 * t237;
t116 = t135 * t292;
t115 = t240 * t250 + (t237 * t265 - t307) * t237;
t112 = rSges(6,3) * t294 - t263;
t96 = t236 * t110;
t94 = t112 * t292;
t92 = t251 * t240;
t91 = t251 * t237;
t90 = t130 * t236 - t176 * t292;
t89 = -t131 * t236 + t176 * t294;
t88 = t269 + t130 + t280;
t87 = t227 + t228 + (-t300 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(8)) * t239) * t237 + t264;
t86 = (t239 * t256 + t155) * t236;
t85 = t240 * t249 + (-t308 + (-rSges(4,2) * t239 + rSges(4,3) * t236) * t237) * t237 + t285;
t84 = -t154 * t292 + t96;
t83 = -t112 * t236 + t138;
t82 = t245 * t240;
t81 = t245 * t237;
t78 = (-t130 * t237 + t131 * t240) * t239;
t75 = t242 + t110;
t74 = ((-rSges(6,3) - pkin(2)) * t239 + t246) * t237 + t263 + t275;
t73 = -t110 * t294 + t94;
t64 = t242 + t290;
t63 = ((-rSges(7,2) - pkin(2)) * t239 + t246) * t237 + t275 + t341;
t60 = t127 + t96 + (-t154 - t185) * t292;
t59 = t138 + t152 + (-t112 - t135) * t236;
t58 = -t292 * t319 + t305;
t57 = -t236 * t289 + t287;
t56 = t130 * t240 + t131 * t237 + t270;
t55 = t122 * t294 + t124 * t191 + t126 * t192;
t54 = t121 * t294 + t123 * t191 + t125 * t192;
t53 = t122 * t292 + t124 * t189 + t126 * t190;
t52 = t121 * t292 + t123 * t189 + t125 * t190;
t35 = t116 + t94 + (-t110 - t134) * t294;
t34 = -t290 * t294 + t306;
t33 = t127 + (-t185 - t319) * t292 + t305;
t32 = t152 + (-t135 - t289) * t236 + t287;
t31 = t110 * t240 + t112 * t237 + t247;
t30 = t116 + (-t134 - t290) * t294 + t306;
t29 = t237 * t54 - t240 * t55;
t28 = t237 * t52 - t240 * t53;
t27 = t237 * t289 + t240 * t290 + t247;
t16 = t236 * t77 + (t237 * t55 + t240 * t54) * t239;
t15 = t236 * t76 + (t237 * t53 + t240 * t52) * t239;
t1 = [Icges(2,3) + t155 + m(6) * (t74 ^ 2 + t75 ^ 2) + m(7) * (t63 ^ 2 + t64 ^ 2) + m(5) * (t87 ^ 2 + t88 ^ 2) + m(4) * (t132 ^ 2 + t133 ^ 2) + m(3) * (t139 ^ 2 + t140 ^ 2) + m(2) * (t206 ^ 2 + t207 ^ 2) + (t302 + t304 + t256 + (Icges(4,3) + Icges(3,2)) * t239 + t331) * t239 + (t301 + t303 + (Icges(3,1) + Icges(4,2)) * t236) * t236 + t332; (-t71 / 0.2e1 - t72 / 0.2e1 - t51 / 0.2e1 - t49 / 0.2e1 + t272 * t240 + (Icges(4,5) * t310 + Icges(3,6) * t327 + t257 * t311 + t261 * t328) * t239 + (Icges(4,4) * t310 + Icges(3,5) * t327 + t258 * t311 + t262 * t328) * t236 - t277) * t240 + (t50 / 0.2e1 + t48 / 0.2e1 + t69 / 0.2e1 + t70 / 0.2e1 + t272 * t237 + (Icges(4,5) * t328 + Icges(3,6) * t311 + t257 * t310 + t261 * t327) * t239 + (Icges(4,4) * t328 + Icges(3,5) * t311 + t258 * t310 + t262 * t327) * t236 + t278) * t237 + m(7) * (t63 * t82 + t64 * t81) + m(6) * (t74 * t92 + t75 * t91) + m(5) * (t117 * t88 + t118 * t87) + m(4) * (t132 * t144 + t133 * t143) + m(3) * (-t139 * t240 - t140 * t237) * t205; m(7) * (t27 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(6) * (t31 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2 + t56 ^ 2) + m(4) * (t143 ^ 2 + t144 ^ 2 + t85 ^ 2) + m(3) * (t205 ^ 2 * t279 + t115 ^ 2) + (t233 * t317 - t29 - t323) * t240 + (t28 + t318 * t231 + (t237 * t317 + t240 * t318) * t240 + t324) * t237; 0.2e1 * ((t237 * t75 + t240 * t74) * t314 + (t237 * t64 + t240 * t63) * t313 + (t237 * t88 + t240 * t87) * t315 + (t132 * t240 + t133 * t237) * t316) * t236; m(7) * (-t239 * t27 + (t237 * t81 + t240 * t82) * t236) + m(6) * (-t239 * t31 + (t237 * t91 + t240 * t92) * t236) + m(5) * (-t239 * t56 + (t117 * t237 + t118 * t240) * t236) + m(4) * (-t239 * t85 + (t143 * t237 + t144 * t240) * t236); 0.2e1 * (t316 + t315 + t314 + t313) * (t236 ^ 2 * t279 + t232); t243 + t86 + (t237 * t277 + t240 * t278) * t239 + m(6) * (t59 * t74 + t60 * t75) + m(7) * (t32 * t63 + t33 * t64) + m(5) * (t87 * t89 + t88 * t90); t15 * t311 + m(7) * (t27 * t30 + t32 * t82 + t33 * t81) + m(6) * (t31 * t35 + t59 * t92 + t60 * t91) + m(5) * (t117 * t90 + t118 * t89 + t56 * t78) + (t28 * t327 + t29 * t311) * t239 + t16 * t310 + (t61 * t237 - t62 * t240) * t312 + t244; m(5) * (-t239 * t78 + (t237 * t90 + t240 * t89) * t236) + m(6) * (-t239 * t35 + (t237 * t60 + t240 * t59) * t236) + m(7) * (-t239 * t30 + (t237 * t33 + t240 * t32) * t236); t236 * t86 + (t240 * t15 + t237 * t16 + t236 * (t62 * t237 + t61 * t240)) * t239 + m(7) * (t30 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t78 ^ 2 + t89 ^ 2 + t90 ^ 2) + t266; m(6) * (t74 * t83 + t75 * t84) + m(7) * (t57 * t63 + t58 * t64) + t243; m(7) * (t27 * t34 + t57 * t82 + t58 * t81) + m(6) * (t31 * t73 + t83 * t92 + t84 * t91) + t244; m(6) * (-t239 * t73 + (t237 * t84 + t240 * t83) * t236) + m(7) * (-t239 * t34 + (t237 * t58 + t240 * t57) * t236); m(7) * (t30 * t34 + t32 * t57 + t33 * t58) + m(6) * (t35 * t73 + t59 * t83 + t60 * t84) + t266; m(7) * (t34 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t73 ^ 2 + t83 ^ 2 + t84 ^ 2) + t266; m(7) * (t177 * t63 - t179 * t64); m(7) * (t177 * t82 - t179 * t81 + t27 * t298); m(7) * (-t222 * t232 + (t177 * t240 - t179 * t237) * t236); m(7) * (t177 * t32 - t179 * t33 + t298 * t30); m(7) * (t177 * t57 - t179 * t58 + t298 * t34); m(7) * (t222 ^ 2 * t232 + t177 ^ 2 + t179 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
