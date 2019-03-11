% Calculate joint inertia matrix for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:22
% EndTime: 2019-03-09 07:38:40
% DurationCPUTime: 7.71s
% Computational Cost: add. (54203->613), mult. (141063->840), div. (0->0), fcn. (186512->16), ass. (0->288)
t282 = sin(qJ(1));
t278 = cos(pkin(6));
t340 = cos(pkin(13));
t308 = t282 * t340;
t276 = sin(pkin(13));
t284 = cos(qJ(1));
t331 = t284 * t276;
t296 = t278 * t308 + t331;
t277 = sin(pkin(6));
t341 = cos(pkin(7));
t310 = t277 * t341;
t339 = sin(pkin(7));
t251 = t282 * t310 + t296 * t339;
t307 = t284 * t340;
t332 = t282 * t276;
t297 = -t278 * t307 + t332;
t250 = -t284 * t310 + t297 * t339;
t360 = m(3) / 0.2e1;
t359 = m(4) / 0.2e1;
t358 = m(5) / 0.2e1;
t357 = m(6) / 0.2e1;
t356 = m(7) / 0.2e1;
t259 = t278 * t331 + t308;
t281 = sin(qJ(3));
t292 = t297 * t341;
t309 = t277 * t339;
t346 = cos(qJ(3));
t242 = t259 * t346 + (-t284 * t309 - t292) * t281;
t280 = sin(qJ(4));
t345 = cos(qJ(4));
t222 = t242 * t280 - t250 * t345;
t355 = t222 / 0.2e1;
t260 = -t278 * t332 + t307;
t290 = t296 * t341;
t244 = t260 * t346 + (t282 * t309 - t290) * t281;
t224 = t244 * t280 - t251 * t345;
t354 = t224 / 0.2e1;
t301 = t341 * t340;
t249 = t278 * t339 * t281 + (t346 * t276 + t281 * t301) * t277;
t258 = t278 * t341 - t340 * t309;
t239 = t249 * t280 - t258 * t345;
t353 = t239 / 0.2e1;
t303 = t346 * t339;
t298 = t277 * t303;
t241 = t259 * t281 + t284 * t298 + t346 * t292;
t352 = t241 / 0.2e1;
t243 = t260 * t281 - t282 * t298 + t346 * t290;
t351 = t243 / 0.2e1;
t335 = t276 * t277;
t248 = -t277 * t346 * t301 - t278 * t303 + t281 * t335;
t350 = t248 / 0.2e1;
t349 = t250 / 0.2e1;
t348 = t251 / 0.2e1;
t347 = t258 / 0.2e1;
t283 = cos(qJ(5));
t270 = pkin(5) * t283 + pkin(4);
t344 = -pkin(4) + t270;
t223 = t242 * t345 + t250 * t280;
t275 = qJ(5) + qJ(6);
t271 = sin(t275);
t272 = cos(t275);
t178 = -t223 * t271 + t241 * t272;
t179 = t223 * t272 + t241 * t271;
t118 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t222;
t120 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t222;
t122 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t222;
t240 = t249 * t345 + t258 * t280;
t209 = -t240 * t271 + t248 * t272;
t210 = t240 * t272 + t248 * t271;
t60 = t118 * t239 + t120 * t209 + t122 * t210;
t343 = t60 * t222;
t225 = t244 * t345 + t251 * t280;
t180 = -t225 * t271 + t243 * t272;
t181 = t225 * t272 + t243 * t271;
t119 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t224;
t121 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t224;
t123 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t224;
t61 = t119 * t239 + t121 * t209 + t123 * t210;
t342 = t61 * t224;
t279 = sin(qJ(5));
t338 = t241 * t279;
t337 = t243 * t279;
t336 = t248 * t279;
t334 = t277 * t282;
t333 = t277 * t284;
t220 = t222 * pkin(11);
t285 = -pkin(12) - pkin(11);
t319 = pkin(5) * t338;
t114 = -t222 * t285 + t344 * t223 - t220 + t319;
t299 = -t179 * rSges(7,1) - t178 * rSges(7,2);
t124 = t222 * rSges(7,3) - t299;
t330 = t114 + t124;
t173 = t225 * pkin(4) + t224 * pkin(11);
t312 = pkin(5) * t337 - t224 * t285 + t225 * t270;
t115 = -t173 + t312;
t125 = t181 * rSges(7,1) + t180 * rSges(7,2) + t224 * rSges(7,3);
t329 = t115 + t125;
t184 = -t223 * t279 + t241 * t283;
t185 = t223 * t283 + t338;
t132 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t222;
t172 = t223 * pkin(4) + t220;
t328 = -t132 - t172;
t186 = -t225 * t279 + t243 * t283;
t187 = t225 * t283 + t337;
t133 = t187 * rSges(6,1) + t186 * rSges(6,2) + t224 * rSges(6,3);
t327 = -t133 - t173;
t147 = pkin(5) * t336 + t344 * t240 + (-pkin(11) - t285) * t239;
t151 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t239;
t326 = t147 + t151;
t214 = -t240 * t279 + t248 * t283;
t215 = t240 * t283 + t336;
t155 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t239;
t205 = t240 * pkin(4) + t239 * pkin(11);
t325 = -t155 - t205;
t206 = t242 * pkin(3) + t241 * pkin(10);
t201 = t251 * t206;
t324 = t251 * t172 + t201;
t207 = t244 * pkin(3) + t243 * pkin(10);
t202 = t258 * t207;
t323 = t258 * t173 + t202;
t230 = pkin(3) * t249 + pkin(10) * t248;
t212 = t250 * t230;
t322 = t250 * t205 + t212;
t321 = t284 * pkin(1) + qJ(2) * t334;
t148 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t239;
t149 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t239;
t150 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t239;
t80 = t239 * t148 + t209 * t149 + t210 * t150;
t75 = t80 * t239;
t26 = t342 + t75 + t343;
t49 = t118 * t222 + t120 * t178 + t122 * t179;
t50 = t119 * t222 + t121 * t178 + t123 * t179;
t69 = t148 * t222 + t149 * t178 + t150 * t179;
t7 = t222 * t49 + t224 * t50 + t239 * t69;
t51 = t118 * t224 + t120 * t180 + t122 * t181;
t52 = t119 * t224 + t121 * t180 + t123 * t181;
t70 = t148 * t224 + t149 * t180 + t150 * t181;
t8 = t222 * t51 + t224 * t52 + t239 * t70;
t320 = t222 * t7 + t224 * t8 + t239 * t26;
t126 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t222;
t128 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t222;
t130 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t222;
t62 = t126 * t239 + t128 * t214 + t130 * t215;
t152 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t239;
t153 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t239;
t154 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t239;
t71 = t152 * t222 + t153 * t184 + t154 * t185;
t318 = t62 / 0.2e1 + t71 / 0.2e1;
t127 = Icges(6,5) * t187 + Icges(6,6) * t186 + Icges(6,3) * t224;
t129 = Icges(6,4) * t187 + Icges(6,2) * t186 + Icges(6,6) * t224;
t131 = Icges(6,1) * t187 + Icges(6,4) * t186 + Icges(6,5) * t224;
t63 = t127 * t239 + t129 * t214 + t131 * t215;
t72 = t152 * t224 + t153 * t186 + t154 * t187;
t317 = t63 / 0.2e1 + t72 / 0.2e1;
t316 = -t172 - t330;
t315 = -t173 - t329;
t83 = t239 * t152 + t214 * t153 + t215 * t154;
t314 = -t205 - t326;
t188 = Icges(5,5) * t240 - Icges(5,6) * t239 + Icges(5,3) * t248;
t189 = Icges(5,4) * t240 - Icges(5,2) * t239 + Icges(5,6) * t248;
t190 = Icges(5,1) * t240 - Icges(5,4) * t239 + Icges(5,5) * t248;
t107 = t248 * t188 - t239 * t189 + t240 * t190;
t226 = Icges(4,5) * t249 - Icges(4,6) * t248 + Icges(4,3) * t258;
t227 = Icges(4,4) * t249 - Icges(4,2) * t248 + Icges(4,6) * t258;
t228 = Icges(4,1) * t249 - Icges(4,4) * t248 + Icges(4,5) * t258;
t313 = t258 * t226 - t248 * t227 + t249 * t228;
t163 = t225 * rSges(5,1) - t224 * rSges(5,2) + t243 * rSges(5,3);
t200 = t244 * rSges(4,1) - t243 * rSges(4,2) + t251 * rSges(4,3);
t311 = -t282 * pkin(1) + qJ(2) * t333;
t306 = t343 / 0.2e1 + t342 / 0.2e1 + t69 * t355 + t70 * t354 + t75;
t15 = t241 * t49 + t243 * t50 + t248 * t69;
t16 = t241 * t51 + t243 * t52 + t248 * t70;
t76 = t80 * t248;
t29 = t60 * t241 + t61 * t243 + t76;
t302 = t15 * t355 + t16 * t354 + t26 * t350 + t29 * t353 + t8 * t351 + t7 * t352;
t17 = t250 * t49 + t251 * t50 + t258 * t69;
t18 = t250 * t51 + t251 * t52 + t258 * t70;
t78 = t80 * t258;
t31 = t60 * t250 + t61 * t251 + t78;
t300 = t17 * t355 + t18 * t354 + t26 * t347 + t31 * t353 + t8 * t348 + t7 * t349;
t199 = rSges(4,1) * t242 - rSges(4,2) * t241 + rSges(4,3) * t250;
t162 = rSges(5,1) * t223 - rSges(5,2) * t222 + rSges(5,3) * t241;
t295 = -t259 * pkin(2) - pkin(9) * t250 + t311;
t100 = t188 * t241 - t189 * t222 + t190 * t223;
t156 = Icges(5,5) * t223 - Icges(5,6) * t222 + Icges(5,3) * t241;
t158 = Icges(5,4) * t223 - Icges(5,2) * t222 + Icges(5,6) * t241;
t160 = Icges(5,1) * t223 - Icges(5,4) * t222 + Icges(5,5) * t241;
t89 = t156 * t248 - t158 * t239 + t160 * t240;
t294 = t60 / 0.2e1 + t69 / 0.2e1 + t100 / 0.2e1 + t89 / 0.2e1 + t318;
t101 = t188 * t243 - t189 * t224 + t190 * t225;
t157 = Icges(5,5) * t225 - Icges(5,6) * t224 + Icges(5,3) * t243;
t159 = Icges(5,4) * t225 - Icges(5,2) * t224 + Icges(5,6) * t243;
t161 = Icges(5,1) * t225 - Icges(5,4) * t224 + Icges(5,5) * t243;
t90 = t157 * t248 - t159 * t239 + t161 * t240;
t293 = t61 / 0.2e1 + t70 / 0.2e1 + t90 / 0.2e1 + t101 / 0.2e1 + t317;
t288 = -t206 + t295;
t287 = t260 * pkin(2) + pkin(9) * t251 + t321;
t286 = t207 + t287;
t266 = rSges(2,1) * t284 - rSges(2,2) * t282;
t265 = -rSges(2,1) * t282 - rSges(2,2) * t284;
t238 = t260 * rSges(3,1) - rSges(3,2) * t296 + rSges(3,3) * t334 + t321;
t237 = -t259 * rSges(3,1) + rSges(3,2) * t297 + rSges(3,3) * t333 + t311;
t229 = rSges(4,1) * t249 - rSges(4,2) * t248 + rSges(4,3) * t258;
t198 = Icges(4,1) * t244 - Icges(4,4) * t243 + Icges(4,5) * t251;
t197 = Icges(4,1) * t242 - Icges(4,4) * t241 + Icges(4,5) * t250;
t196 = Icges(4,4) * t244 - Icges(4,2) * t243 + Icges(4,6) * t251;
t195 = Icges(4,4) * t242 - Icges(4,2) * t241 + Icges(4,6) * t250;
t194 = Icges(4,5) * t244 - Icges(4,6) * t243 + Icges(4,3) * t251;
t193 = Icges(4,5) * t242 - Icges(4,6) * t241 + Icges(4,3) * t250;
t191 = rSges(5,1) * t240 - rSges(5,2) * t239 + rSges(5,3) * t248;
t175 = t241 * t205;
t170 = t287 + t200;
t169 = -t199 + t295;
t165 = t248 * t173;
t164 = t243 * t172;
t146 = t200 * t258 - t229 * t251;
t145 = -t199 * t258 + t229 * t250;
t140 = t222 * t151;
t139 = t199 * t251 - t200 * t250;
t136 = t313 * t258;
t135 = t226 * t251 - t227 * t243 + t228 * t244;
t134 = t226 * t250 - t227 * t241 + t228 * t242;
t117 = t286 + t163;
t116 = -t162 + t288;
t113 = t239 * t125;
t112 = t163 * t248 - t191 * t243;
t111 = -t162 * t248 + t191 * t241;
t110 = t224 * t124;
t109 = t194 * t258 - t196 * t248 + t198 * t249;
t108 = t193 * t258 - t195 * t248 + t197 * t249;
t106 = t162 * t243 - t163 * t241;
t105 = t107 * t258;
t104 = t107 * t248;
t103 = t258 * t163 + t202 + (-t191 - t230) * t251;
t102 = t250 * t191 + t212 + (-t162 - t206) * t258;
t99 = t286 - t327;
t98 = t288 + t328;
t97 = t133 * t239 - t155 * t224;
t96 = -t132 * t239 + t155 * t222;
t95 = -t151 * t224 + t113;
t94 = -t124 * t239 + t140;
t93 = t251 * t162 + t201 + (-t163 - t207) * t250;
t92 = t286 + t312 + t125;
t91 = -t319 - t223 * t270 + (-rSges(7,3) + t285) * t222 + t288 + t299;
t88 = t157 * t243 - t159 * t224 + t161 * t225;
t87 = t156 * t243 - t158 * t224 + t160 * t225;
t86 = t157 * t241 - t159 * t222 + t161 * t223;
t85 = t156 * t241 - t158 * t222 + t160 * t223;
t84 = t132 * t224 - t133 * t222;
t82 = -t125 * t222 + t110;
t81 = t83 * t258;
t79 = t83 * t248;
t77 = t83 * t239;
t74 = t248 * t133 + t243 * t325 + t165;
t73 = t241 * t155 + t248 * t328 + t175;
t66 = t258 * t133 + (-t230 + t325) * t251 + t323;
t65 = t250 * t155 + (-t206 + t328) * t258 + t322;
t64 = t243 * t132 + t241 * t327 + t164;
t57 = t251 * t132 + (-t207 + t327) * t250 + t324;
t56 = t127 * t224 + t129 * t186 + t131 * t187;
t55 = t126 * t224 + t128 * t186 + t130 * t187;
t54 = t127 * t222 + t129 * t184 + t131 * t185;
t53 = t126 * t222 + t128 * t184 + t130 * t185;
t48 = t239 * t115 - t224 * t326 + t113;
t47 = t222 * t147 - t239 * t330 + t140;
t46 = t243 * t314 + t248 * t329 + t165;
t45 = t241 * t326 + t248 * t316 + t175;
t44 = t329 * t258 + (-t230 + t314) * t251 + t323;
t43 = t326 * t250 + (-t206 + t316) * t258 + t322;
t42 = t224 * t114 - t222 * t329 + t110;
t41 = t241 * t315 + t243 * t330 + t164;
t40 = t89 * t250 + t90 * t251 + t105;
t39 = t89 * t241 + t90 * t243 + t104;
t38 = t330 * t251 + (-t207 + t315) * t250 + t324;
t37 = t101 * t258 + t250 * t87 + t251 * t88;
t36 = t100 * t258 + t250 * t85 + t251 * t86;
t35 = t101 * t248 + t241 * t87 + t243 * t88;
t34 = t100 * t248 + t241 * t85 + t243 * t86;
t33 = t62 * t250 + t63 * t251 + t81;
t32 = t62 * t241 + t63 * t243 + t79;
t28 = t62 * t222 + t63 * t224 + t77;
t22 = t250 * t55 + t251 * t56 + t258 * t72;
t21 = t250 * t53 + t251 * t54 + t258 * t71;
t20 = t241 * t55 + t243 * t56 + t248 * t72;
t19 = t241 * t53 + t243 * t54 + t248 * t71;
t14 = t222 * t55 + t224 * t56 + t239 * t72;
t13 = t222 * t53 + t224 * t54 + t239 * t71;
t1 = [t277 * t340 * (Icges(3,6) * t278 + (Icges(3,4) * t276 + Icges(3,2) * t340) * t277) + t313 + t278 * (Icges(3,3) * t278 + (Icges(3,5) * t276 + Icges(3,6) * t340) * t277) + m(7) * (t91 ^ 2 + t92 ^ 2) + m(6) * (t98 ^ 2 + t99 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2) + m(4) * (t169 ^ 2 + t170 ^ 2) + m(3) * (t237 ^ 2 + t238 ^ 2) + m(2) * (t265 ^ 2 + t266 ^ 2) + Icges(2,3) + (Icges(3,5) * t278 + (Icges(3,1) * t276 + Icges(3,4) * t340) * t277) * t335 + t107 + t83 + t80; 0.2e1 * ((t282 * t91 - t284 * t92) * t356 + (t282 * t98 - t284 * t99) * t357 + (t116 * t282 - t117 * t284) * t358 + (t169 * t282 - t170 * t284) * t359 + (t237 * t282 - t238 * t284) * t360) * t277; 0.2e1 * (t360 + t359 + t358 + t357 + t356) * (t278 ^ 2 + (t282 ^ 2 + t284 ^ 2) * t277 ^ 2); t78 + t81 + t105 + t136 + m(7) * (t43 * t91 + t44 * t92) + m(6) * (t65 * t98 + t66 * t99) + m(5) * (t102 * t116 + t103 * t117) + m(4) * (t145 * t169 + t146 * t170) + (t109 / 0.2e1 + t135 / 0.2e1 + t293) * t251 + (t108 / 0.2e1 + t134 / 0.2e1 + t294) * t250; m(4) * (t139 * t278 + (t145 * t282 - t146 * t284) * t277) + m(5) * (t93 * t278 + (t102 * t282 - t103 * t284) * t277) + m(6) * (t57 * t278 + (t282 * t65 - t284 * t66) * t277) + m(7) * (t38 * t278 + (t282 * t43 - t284 * t44) * t277); (t136 + t31 + t33 + t40) * t258 + (t18 + t22 + t37 + (t194 * t251 - t196 * t243 + t198 * t244) * t251 + (t109 + t135) * t258) * t251 + (t17 + t21 + t36 + (t193 * t250 - t195 * t241 + t197 * t242) * t250 + (t108 + t134) * t258 + (t193 * t251 + t194 * t250 - t195 * t243 - t196 * t241 + t197 * t244 + t198 * t242) * t251) * t250 + m(7) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t57 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2 + t93 ^ 2) + m(4) * (t139 ^ 2 + t145 ^ 2 + t146 ^ 2); t76 + t79 + t104 + m(7) * (t45 * t91 + t46 * t92) + m(6) * (t73 * t98 + t74 * t99) + m(5) * (t111 * t116 + t112 * t117) + t293 * t243 + t294 * t241; m(5) * (t106 * t278 + (t111 * t282 - t112 * t284) * t277) + m(6) * (t64 * t278 + (t282 * t73 - t284 * t74) * t277) + m(7) * (t41 * t278 + (t282 * t45 - t284 * t46) * t277); (t29 / 0.2e1 + t32 / 0.2e1 + t39 / 0.2e1) * t258 + (t16 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1) * t251 + (t15 / 0.2e1 + t19 / 0.2e1 + t34 / 0.2e1) * t250 + (t31 / 0.2e1 + t33 / 0.2e1 + t40 / 0.2e1) * t248 + (t18 / 0.2e1 + t22 / 0.2e1 + t37 / 0.2e1) * t243 + (t17 / 0.2e1 + t21 / 0.2e1 + t36 / 0.2e1) * t241 + m(7) * (t38 * t41 + t43 * t45 + t44 * t46) + m(6) * (t64 * t57 + t65 * t73 + t66 * t74) + m(5) * (t102 * t111 + t103 * t112 + t106 * t93); (t29 + t32 + t39) * t248 + (t16 + t20 + t35) * t243 + (t15 + t19 + t34) * t241 + m(7) * (t41 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t64 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t106 ^ 2 + t111 ^ 2 + t112 ^ 2); t77 + t317 * t224 + t318 * t222 + m(7) * (t47 * t91 + t48 * t92) + m(6) * (t96 * t98 + t97 * t99) + t306; m(6) * (t84 * t278 + (t282 * t96 - t284 * t97) * t277) + m(7) * (t42 * t278 + (t282 * t47 - t284 * t48) * t277); t13 * t349 + t33 * t353 + t21 * t355 + t22 * t354 + t14 * t348 + t28 * t347 + m(7) * (t38 * t42 + t43 * t47 + t44 * t48) + m(6) * (t57 * t84 + t65 * t96 + t66 * t97) + t300; t32 * t353 + t28 * t350 + t13 * t352 + t19 * t355 + t20 * t354 + t14 * t351 + m(7) * (t41 * t42 + t45 * t47 + t46 * t48) + m(6) * (t64 * t84 + t73 * t96 + t74 * t97) + t302; t222 * t13 + t224 * t14 + t239 * t28 + m(7) * (t42 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t84 ^ 2 + t96 ^ 2 + t97 ^ 2) + t320; m(7) * (t91 * t94 + t92 * t95) + t306; m(7) * (t82 * t278 + (t282 * t94 - t284 * t95) * t277); m(7) * (t38 * t82 + t43 * t94 + t44 * t95) + t300; m(7) * (t41 * t82 + t45 * t94 + t46 * t95) + t302; m(7) * (t42 * t82 + t47 * t94 + t48 * t95) + t320; m(7) * (t82 ^ 2 + t94 ^ 2 + t95 ^ 2) + t320;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
