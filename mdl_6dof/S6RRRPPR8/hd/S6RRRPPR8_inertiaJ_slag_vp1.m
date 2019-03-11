% Calculate joint inertia matrix for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:43
% EndTime: 2019-03-09 16:03:53
% DurationCPUTime: 4.56s
% Computational Cost: add. (13233->627), mult. (33806->868), div. (0->0), fcn. (42332->10), ass. (0->280)
t303 = m(6) / 0.2e1 + m(7) / 0.2e1;
t338 = 0.2e1 * t303;
t267 = cos(pkin(6));
t269 = sin(qJ(2));
t266 = sin(pkin(6));
t329 = sin(qJ(3));
t294 = t266 * t329;
t330 = cos(qJ(3));
t247 = -t267 * t330 + t269 * t294;
t268 = sin(qJ(6));
t271 = cos(qJ(6));
t272 = cos(qJ(2));
t324 = t266 * t272;
t225 = -t247 * t268 + t271 * t324;
t226 = t247 * t271 + t268 * t324;
t295 = t266 * t330;
t248 = t267 * t329 + t269 * t295;
t122 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t248;
t123 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t248;
t124 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t248;
t55 = t248 * t122 + t225 * t123 + t226 * t124;
t185 = Icges(6,5) * t247 - Icges(6,6) * t248 + Icges(6,3) * t324;
t188 = Icges(6,4) * t247 - Icges(6,2) * t248 + Icges(6,6) * t324;
t191 = Icges(6,1) * t247 - Icges(6,4) * t248 + Icges(6,5) * t324;
t97 = t185 * t324 - t248 * t188 + t247 * t191;
t337 = -t55 - t97;
t273 = cos(qJ(1));
t323 = t266 * t273;
t336 = t266 ^ 2;
t270 = sin(qJ(1));
t320 = t270 * t272;
t321 = t269 * t273;
t251 = t267 * t321 + t320;
t228 = t251 * t330 - t273 * t294;
t319 = t272 * t273;
t322 = t269 * t270;
t254 = -t267 * t322 + t319;
t230 = t254 * t330 + t270 * t294;
t227 = t251 * t329 + t273 * t295;
t249 = -t267 * t319 + t322;
t178 = -t227 * t268 - t249 * t271;
t179 = t227 * t271 - t249 * t268;
t104 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t228;
t106 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t228;
t108 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t228;
t229 = t254 * t329 - t270 * t295;
t252 = t267 * t320 + t321;
t180 = -t229 * t268 - t252 * t271;
t181 = t229 * t271 - t252 * t268;
t36 = t104 * t230 + t106 * t180 + t108 * t181;
t105 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t230;
t107 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t230;
t109 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t230;
t37 = t105 * t230 + t107 * t180 + t109 * t181;
t46 = t122 * t230 + t123 * t180 + t124 * t181;
t2 = t228 * t36 + t230 * t37 + t248 * t46;
t335 = t2 / 0.2e1;
t334 = -pkin(3) - pkin(4);
t333 = t228 / 0.2e1;
t332 = t230 / 0.2e1;
t331 = t248 / 0.2e1;
t328 = pkin(5) * t227;
t327 = pkin(9) * t252;
t197 = Icges(3,5) * t251 - Icges(3,6) * t249 - Icges(3,3) * t323;
t326 = t197 * t273;
t325 = t266 * t270;
t285 = -rSges(7,1) * t179 - rSges(7,2) * t178;
t110 = rSges(7,3) * t228 - t285;
t318 = pkin(10) * t228 + t110 + t328;
t111 = t181 * rSges(7,1) + t180 * rSges(7,2) + t230 * rSges(7,3);
t317 = t229 * pkin(5) + pkin(10) * t230 + t111;
t125 = rSges(7,1) * t226 + rSges(7,2) * t225 + rSges(7,3) * t248;
t316 = pkin(5) * t247 + pkin(10) * t248 + t125;
t148 = t230 * rSges(5,1) + t252 * rSges(5,2) + t229 * rSges(5,3);
t160 = t230 * pkin(3) + qJ(4) * t229;
t315 = -t148 - t160;
t214 = t227 * qJ(4);
t159 = pkin(3) * t228 + t214;
t150 = t252 * t159;
t239 = t249 * qJ(5);
t183 = pkin(4) * t228 - t239;
t314 = t252 * t183 + t150;
t210 = pkin(3) * t248 + qJ(4) * t247;
t313 = t159 * t324 + t249 * t210;
t246 = t254 * pkin(2);
t213 = t246 + t327;
t209 = t267 * t213;
t312 = t267 * t160 + t209;
t212 = pkin(2) * t251 + t249 * pkin(9);
t311 = -t159 - t212;
t223 = t230 * pkin(4);
t184 = -qJ(5) * t252 + t223;
t310 = -t160 - t184;
t186 = Icges(5,5) * t248 - Icges(5,6) * t324 + Icges(5,3) * t247;
t192 = Icges(5,1) * t248 - Icges(5,4) * t324 + Icges(5,5) * t247;
t309 = t247 * t186 + t248 * t192;
t190 = Icges(4,4) * t248 - Icges(4,2) * t247 - Icges(4,6) * t324;
t193 = Icges(4,1) * t248 - Icges(4,4) * t247 - Icges(4,5) * t324;
t308 = -t247 * t190 + t248 * t193;
t195 = rSges(5,1) * t248 - rSges(5,2) * t324 + rSges(5,3) * t247;
t307 = -t195 - t210;
t306 = t212 * t325 + t213 * t323;
t233 = pkin(4) * t248 + qJ(5) * t324;
t305 = -t210 - t233;
t304 = t273 * pkin(1) + pkin(8) * t325;
t39 = t104 * t248 + t106 * t225 + t108 * t226;
t45 = t122 * t228 + t123 * t178 + t124 * t179;
t302 = t45 / 0.2e1 + t39 / 0.2e1;
t40 = t105 * t248 + t107 * t225 + t109 * t226;
t301 = t46 / 0.2e1 + t40 / 0.2e1;
t147 = t229 * rSges(6,1) - t230 * rSges(6,2) - t252 * rSges(6,3);
t300 = -t147 + t310;
t299 = t267 * t184 + t312;
t298 = -t183 + t311;
t194 = rSges(6,1) * t247 - rSges(6,2) * t248 + rSges(6,3) * t324;
t297 = -t194 + t305;
t149 = t230 * rSges(4,1) - t229 * rSges(4,2) + t252 * rSges(4,3);
t235 = Icges(3,3) * t267 + (Icges(3,5) * t269 + Icges(3,6) * t272) * t266;
t236 = Icges(3,6) * t267 + (Icges(3,4) * t269 + Icges(3,2) * t272) * t266;
t237 = Icges(3,5) * t267 + (Icges(3,1) * t269 + Icges(3,4) * t272) * t266;
t296 = t266 * t269 * t237 + t267 * t235 + t236 * t324;
t204 = t254 * rSges(3,1) - t252 * rSges(3,2) + rSges(3,3) * t325;
t293 = -t270 * pkin(1) + pkin(8) * t323;
t196 = rSges(4,1) * t248 - rSges(4,2) * t247 - rSges(4,3) * t324;
t255 = (pkin(2) * t269 - pkin(9) * t272) * t266;
t292 = t266 * (-t196 - t255);
t291 = t310 - t317;
t290 = t305 - t316;
t289 = t159 * t325 + t160 * t323 + t306;
t288 = t183 * t324 + t249 * t233 + t313;
t287 = t266 * (-t255 + t307);
t286 = -rSges(6,1) * t227 + rSges(6,3) * t249;
t284 = -rSges(5,2) * t249 - rSges(5,3) * t227;
t283 = t266 * (-t255 + t297);
t282 = t183 * t325 + t184 * t323 + t289;
t281 = t266 * (-t255 + t290);
t280 = -t212 + t293;
t279 = t160 + t246 + t304;
t146 = rSges(4,1) * t228 - rSges(4,2) * t227 + rSges(4,3) * t249;
t278 = -t214 + t280;
t203 = t251 * rSges(3,1) - t249 * rSges(3,2) - rSges(3,3) * t323;
t277 = t239 + t278;
t276 = t223 + (pkin(9) - qJ(5)) * t252 + t279;
t126 = Icges(6,5) * t227 - Icges(6,6) * t228 - Icges(6,3) * t249;
t132 = Icges(6,4) * t227 - Icges(6,2) * t228 - Icges(6,6) * t249;
t138 = Icges(6,1) * t227 - Icges(6,4) * t228 - Icges(6,5) * t249;
t70 = t126 * t324 - t132 * t248 + t138 * t247;
t128 = Icges(5,5) * t228 + Icges(5,6) * t249 + Icges(5,3) * t227;
t134 = Icges(5,4) * t228 + Icges(5,2) * t249 + Icges(5,6) * t227;
t140 = Icges(5,1) * t228 + Icges(5,4) * t249 + Icges(5,5) * t227;
t72 = t128 * t247 - t134 * t324 + t140 * t248;
t130 = Icges(4,5) * t228 - Icges(4,6) * t227 + Icges(4,3) * t249;
t136 = Icges(4,4) * t228 - Icges(4,2) * t227 + Icges(4,6) * t249;
t142 = Icges(4,1) * t228 - Icges(4,4) * t227 + Icges(4,5) * t249;
t74 = -t130 * t324 - t136 * t247 + t142 * t248;
t83 = -t185 * t249 - t188 * t228 + t191 * t227;
t189 = Icges(5,4) * t248 - Icges(5,2) * t324 + Icges(5,6) * t247;
t84 = t186 * t227 + t189 * t249 + t192 * t228;
t187 = Icges(4,5) * t248 - Icges(4,6) * t247 - Icges(4,3) * t324;
t85 = t187 * t249 - t190 * t227 + t193 * t228;
t275 = t74 / 0.2e1 + t72 / 0.2e1 + t70 / 0.2e1 + t85 / 0.2e1 + t84 / 0.2e1 + t83 / 0.2e1 + t302;
t127 = Icges(6,5) * t229 - Icges(6,6) * t230 - Icges(6,3) * t252;
t133 = Icges(6,4) * t229 - Icges(6,2) * t230 - Icges(6,6) * t252;
t139 = Icges(6,1) * t229 - Icges(6,4) * t230 - Icges(6,5) * t252;
t71 = t127 * t324 - t133 * t248 + t139 * t247;
t129 = Icges(5,5) * t230 + Icges(5,6) * t252 + Icges(5,3) * t229;
t135 = Icges(5,4) * t230 + Icges(5,2) * t252 + Icges(5,6) * t229;
t141 = Icges(5,1) * t230 + Icges(5,4) * t252 + Icges(5,5) * t229;
t73 = t129 * t247 - t135 * t324 + t141 * t248;
t131 = Icges(4,5) * t230 - Icges(4,6) * t229 + Icges(4,3) * t252;
t137 = Icges(4,4) * t230 - Icges(4,2) * t229 + Icges(4,6) * t252;
t143 = Icges(4,1) * t230 - Icges(4,4) * t229 + Icges(4,5) * t252;
t75 = -t131 * t324 - t137 * t247 + t143 * t248;
t86 = -t185 * t252 - t188 * t230 + t191 * t229;
t87 = t186 * t229 + t189 * t252 + t192 * t230;
t88 = t187 * t252 - t190 * t229 + t193 * t230;
t274 = t73 / 0.2e1 + t71 / 0.2e1 + t88 / 0.2e1 + t87 / 0.2e1 + t86 / 0.2e1 + t75 / 0.2e1 + t301;
t257 = rSges(2,1) * t273 - t270 * rSges(2,2);
t256 = -t270 * rSges(2,1) - rSges(2,2) * t273;
t238 = t267 * rSges(3,3) + (rSges(3,1) * t269 + rSges(3,2) * t272) * t266;
t202 = Icges(3,1) * t254 - Icges(3,4) * t252 + Icges(3,5) * t325;
t201 = Icges(3,1) * t251 - Icges(3,4) * t249 - Icges(3,5) * t323;
t200 = Icges(3,4) * t254 - Icges(3,2) * t252 + Icges(3,6) * t325;
t199 = Icges(3,4) * t251 - Icges(3,2) * t249 - Icges(3,6) * t323;
t198 = Icges(3,5) * t254 - Icges(3,6) * t252 + Icges(3,3) * t325;
t173 = t204 + t304;
t172 = -t203 + t293;
t153 = -t267 * t203 - t238 * t323;
t152 = t204 * t267 - t238 * t325;
t145 = rSges(5,1) * t228 - t284;
t144 = -rSges(6,2) * t228 - t286;
t121 = t296 * t267;
t119 = (t203 * t270 + t204 * t273) * t266;
t118 = t235 * t325 - t236 * t252 + t237 * t254;
t117 = -t235 * t323 - t249 * t236 + t251 * t237;
t113 = t213 + t149 + t304;
t112 = -t146 + t280;
t103 = -t149 * t324 - t196 * t252;
t102 = t146 * t324 + t196 * t249;
t101 = t267 * t198 + (t200 * t272 + t202 * t269) * t266;
t100 = t267 * t197 + (t199 * t272 + t201 * t269) * t266;
t99 = -t187 * t324 + t308;
t98 = -t189 * t324 + t309;
t96 = t99 * t267;
t95 = t98 * t267;
t94 = t97 * t267;
t93 = t148 + t279 + t327;
t92 = (-rSges(5,1) - pkin(3)) * t228 + t278 + t284;
t91 = t146 * t252 - t149 * t249;
t90 = (-t146 - t212) * t267 + t273 * t292;
t89 = t149 * t267 + t270 * t292 + t209;
t82 = t276 + t147;
t81 = (rSges(6,2) + t334) * t228 + t277 + t286;
t80 = (t146 * t270 + t149 * t273) * t266 + t306;
t79 = t111 * t248 - t125 * t230;
t78 = -t110 * t248 + t125 * t228;
t77 = t252 * t307 + t315 * t324;
t76 = t145 * t324 + t195 * t249 + t313;
t69 = (-t145 + t311) * t267 + t273 * t287;
t68 = t148 * t267 + t270 * t287 + t312;
t67 = t131 * t252 - t137 * t229 + t143 * t230;
t66 = t130 * t252 - t136 * t229 + t142 * t230;
t65 = t129 * t229 + t135 * t252 + t141 * t230;
t64 = t128 * t229 + t134 * t252 + t140 * t230;
t63 = -t127 * t252 - t133 * t230 + t139 * t229;
t62 = -t126 * t252 - t132 * t230 + t138 * t229;
t61 = t131 * t249 - t137 * t227 + t143 * t228;
t60 = t130 * t249 - t136 * t227 + t142 * t228;
t59 = t129 * t227 + t135 * t249 + t141 * t228;
t58 = t128 * t227 + t134 * t249 + t140 * t228;
t57 = -t127 * t249 - t133 * t228 + t139 * t227;
t56 = -t126 * t249 - t132 * t228 + t138 * t227;
t54 = t55 * t267;
t53 = t55 * t248;
t52 = t110 * t230 - t111 * t228;
t51 = t145 * t252 + t249 * t315 + t150;
t50 = t276 + t317;
t49 = -t328 + (-rSges(7,3) - pkin(10) + t334) * t228 + t277 + t285;
t48 = t252 * t297 + t300 * t324;
t47 = t144 * t324 + t194 * t249 + t288;
t44 = (-t144 + t298) * t267 + t273 * t283;
t43 = t147 * t267 + t270 * t283 + t299;
t42 = (t145 * t270 + t148 * t273) * t266 + t289;
t41 = t144 * t252 + t249 * t300 + t314;
t38 = (t144 * t270 + t147 * t273) * t266 + t282;
t35 = t105 * t228 + t107 * t178 + t109 * t179;
t34 = t104 * t228 + t106 * t178 + t108 * t179;
t33 = t252 * t290 + t291 * t324;
t32 = t249 * t316 + t318 * t324 + t288;
t31 = (t298 - t318) * t267 + t273 * t281;
t30 = t267 * t317 + t270 * t281 + t299;
t29 = t96 + (t75 * t270 - t74 * t273) * t266;
t28 = t95 + (t73 * t270 - t72 * t273) * t266;
t27 = t94 + (t71 * t270 - t70 * t273) * t266;
t26 = t249 * t291 + t252 * t318 + t314;
t25 = (t270 * t318 + t273 * t317) * t266 + t282;
t24 = t74 * t249 + t75 * t252 - t324 * t99;
t23 = t72 * t249 + t73 * t252 - t324 * t98;
t22 = t70 * t249 + t71 * t252 - t324 * t97;
t21 = t88 * t267 + (t270 * t67 - t273 * t66) * t266;
t20 = t87 * t267 + (t270 * t65 - t273 * t64) * t266;
t19 = t86 * t267 + (t270 * t63 - t273 * t62) * t266;
t18 = t85 * t267 + (t270 * t61 - t273 * t60) * t266;
t17 = t84 * t267 + (t270 * t59 - t273 * t58) * t266;
t16 = t83 * t267 + (t270 * t57 - t273 * t56) * t266;
t15 = t249 * t66 + t252 * t67 - t324 * t88;
t14 = t249 * t64 + t252 * t65 - t324 * t87;
t13 = t249 * t62 + t252 * t63 - t324 * t86;
t12 = t249 * t60 + t252 * t61 - t324 * t85;
t11 = t249 * t58 + t252 * t59 - t324 * t84;
t10 = t249 * t56 + t252 * t57 - t324 * t83;
t9 = t54 + (t40 * t270 - t39 * t273) * t266;
t8 = t39 * t249 + t40 * t252 - t324 * t55;
t7 = t39 * t228 + t40 * t230 + t53;
t6 = t46 * t267 + (t270 * t37 - t273 * t36) * t266;
t5 = t45 * t267 + (t270 * t35 - t273 * t34) * t266;
t4 = t249 * t36 + t252 * t37 - t324 * t46;
t3 = t249 * t34 + t252 * t35 - t324 * t45;
t1 = t228 * t34 + t230 * t35 + t248 * t45;
t114 = [t296 + (-t187 - t189) * t324 + m(7) * (t49 ^ 2 + t50 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(6) * (t81 ^ 2 + t82 ^ 2) + m(4) * (t112 ^ 2 + t113 ^ 2) + m(3) * (t172 ^ 2 + t173 ^ 2) + m(2) * (t256 ^ 2 + t257 ^ 2) + Icges(2,3) + t308 + t309 - t337; t94 + t54 + t121 + t96 + t95 + m(7) * (t30 * t50 + t31 * t49) + m(6) * (t43 * t82 + t44 * t81) + m(5) * (t68 * t93 + t69 * t92) + m(4) * (t112 * t90 + t113 * t89) + m(3) * (t152 * t173 + t153 * t172) + ((-t100 / 0.2e1 - t117 / 0.2e1 - t275) * t273 + (t101 / 0.2e1 + t118 / 0.2e1 + t274) * t270) * t266; (t9 + t28 + t27 + t29 + t121) * t267 + m(7) * (t25 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t42 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t80 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(3) * (t119 ^ 2 + t152 ^ 2 + t153 ^ 2) + ((-t5 - t17 - t16 - t18 + ((-t249 * t199 + t251 * t201) * t266 - t336 * t326) * t273) * t273 + (t6 + t20 + t19 + t21 + ((-t200 * t252 + t202 * t254 + (t198 * t270 - t326) * t266) * t270 + (t198 * t323 + t199 * t252 + t249 * t200 - t201 * t254 - t251 * t202) * t273) * t266) * t270 + ((-t100 - t117) * t273 + (t101 + t118) * t270) * t267) * t266; (-t98 - t99 + t337) * t324 + m(7) * (t32 * t49 + t33 * t50) + m(5) * (t76 * t92 + t77 * t93) + m(6) * (t47 * t81 + t48 * t82) + m(4) * (t102 * t112 + t103 * t113) + t274 * t252 + t275 * t249; (t22 / 0.2e1 + t23 / 0.2e1 + t24 / 0.2e1 + t8 / 0.2e1) * t267 + (t21 / 0.2e1 + t19 / 0.2e1 + t20 / 0.2e1 + t6 / 0.2e1) * t252 + (t18 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1 + t5 / 0.2e1) * t249 + m(7) * (t25 * t26 + t30 * t33 + t31 * t32) + m(6) * (t38 * t41 + t43 * t48 + t44 * t47) + m(5) * (t51 * t42 + t68 * t77 + t69 * t76) + m(4) * (t102 * t90 + t103 * t89 + t80 * t91) + ((-t3 / 0.2e1 - t10 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1) * t273 + (-t9 / 0.2e1 - t27 / 0.2e1 - t28 / 0.2e1 - t29 / 0.2e1) * t272 + (t4 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1) * t270) * t266; (-t22 - t23 - t24 - t8) * t324 + (t4 + t14 + t13 + t15) * t252 + (t3 + t11 + t12 + t10) * t249 + m(7) * (t26 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t41 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(5) * (t51 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2 + t91 ^ 2); m(7) * (t227 * t50 + t229 * t49) + m(5) * (t227 * t93 + t229 * t92) + m(6) * (t227 * t82 + t229 * t81); m(7) * (t227 * t30 + t229 * t31 + t247 * t25) + m(6) * (t227 * t43 + t229 * t44 + t247 * t38) + m(5) * (t227 * t68 + t229 * t69 + t247 * t42); m(7) * (t227 * t33 + t229 * t32 + t247 * t26) + m(6) * (t227 * t48 + t229 * t47 + t247 * t41) + m(5) * (t227 * t77 + t229 * t76 + t247 * t51); 0.2e1 * (m(5) / 0.2e1 + t303) * (t227 ^ 2 + t229 ^ 2 + t247 ^ 2); m(7) * (-t249 * t50 - t252 * t49) + m(6) * (-t249 * t82 - t252 * t81); m(7) * (-t249 * t30 + t25 * t324 - t252 * t31) + m(6) * (-t249 * t43 - t252 * t44 + t324 * t38); m(7) * (-t249 * t33 - t252 * t32 + t26 * t324) + m(6) * (-t249 * t48 - t252 * t47 + t324 * t41); (-t227 * t249 - t229 * t252 + t247 * t324) * t338; (t272 ^ 2 * t336 + t249 ^ 2 + t252 ^ 2) * t338; m(7) * (t49 * t78 + t50 * t79) + t53 + t301 * t230 + t302 * t228; m(7) * (t52 * t25 + t30 * t79 + t31 * t78) + t6 * t332 + t267 * t7 / 0.2e1 + t5 * t333 + t9 * t331 + (-t273 * t1 / 0.2e1 + t270 * t335) * t266; -t7 * t324 / 0.2e1 + m(7) * (t52 * t26 + t32 * t78 + t33 * t79) + t4 * t332 + t249 * t1 / 0.2e1 + t3 * t333 + t252 * t335 + t8 * t331; m(7) * (t227 * t79 + t229 * t78 + t247 * t52); m(7) * (-t249 * t79 - t252 * t78 + t324 * t52); t230 * t2 + t228 * t1 + t248 * t7 + m(7) * (t52 ^ 2 + t78 ^ 2 + t79 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t114(1) t114(2) t114(4) t114(7) t114(11) t114(16); t114(2) t114(3) t114(5) t114(8) t114(12) t114(17); t114(4) t114(5) t114(6) t114(9) t114(13) t114(18); t114(7) t114(8) t114(9) t114(10) t114(14) t114(19); t114(11) t114(12) t114(13) t114(14) t114(15) t114(20); t114(16) t114(17) t114(18) t114(19) t114(20) t114(21);];
Mq  = res;
