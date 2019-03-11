% Calculate joint inertia matrix for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:37
% EndTime: 2019-03-09 17:41:50
% DurationCPUTime: 6.14s
% Computational Cost: add. (17865->651), mult. (45158->886), div. (0->0), fcn. (57178->10), ass. (0->287)
t352 = -qJ(6) - pkin(10) - rSges(7,3);
t273 = cos(pkin(6));
t276 = sin(qJ(2));
t272 = sin(pkin(6));
t347 = sin(qJ(3));
t302 = t272 * t347;
t348 = cos(qJ(3));
t254 = -t273 * t348 + t276 * t302;
t275 = sin(qJ(5));
t278 = cos(qJ(5));
t279 = cos(qJ(2));
t340 = t272 * t279;
t232 = t254 * t278 + t275 * t340;
t342 = t254 * t275;
t233 = -t278 * t340 + t342;
t303 = t272 * t348;
t255 = t273 * t347 + t276 * t303;
t141 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t255;
t143 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t255;
t145 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t255;
t71 = t255 * t141 + t232 * t143 + t233 * t145;
t142 = Icges(6,5) * t233 + Icges(6,6) * t232 + Icges(6,3) * t255;
t144 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t255;
t146 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t255;
t72 = t255 * t142 + t232 * t144 + t233 * t146;
t351 = -t71 - t72;
t280 = cos(qJ(1));
t339 = t272 * t280;
t335 = t279 * t280;
t277 = sin(qJ(1));
t338 = t276 * t277;
t259 = -t273 * t338 + t335;
t236 = t259 * t347 - t277 * t303;
t336 = t277 * t279;
t337 = t276 * t280;
t258 = t273 * t336 + t337;
t194 = t236 * t278 - t258 * t275;
t343 = t236 * t275;
t195 = t258 * t278 + t343;
t237 = t259 * t348 + t277 * t302;
t270 = pkin(5) * t278 + pkin(4);
t350 = t195 * rSges(7,1) + t194 * rSges(7,2) + pkin(5) * t343 - t237 * t352 + t258 * t270;
t257 = t273 * t337 + t336;
t234 = t257 * t347 + t280 * t303;
t256 = -t273 * t335 + t338;
t192 = t234 * t278 - t256 * t275;
t344 = t234 * t275;
t193 = t256 * t278 + t344;
t349 = -t193 * rSges(7,1) - t192 * rSges(7,2) - pkin(5) * t344;
t346 = pkin(4) - t270;
t206 = Icges(3,5) * t257 - Icges(3,6) * t256 - Icges(3,3) * t339;
t345 = t206 * t280;
t341 = t272 * t277;
t235 = t257 * t348 - t280 * t302;
t230 = t235 * pkin(10);
t334 = -t352 * t235 - t256 * t346 - t230 - t349;
t197 = t258 * pkin(4) + pkin(10) * t237;
t333 = -t197 + t350;
t332 = rSges(7,1) * t233 + rSges(7,2) * t232 + pkin(5) * t342 + t346 * t340 + (-pkin(10) - t352) * t255;
t162 = t258 * rSges(5,1) - t237 * rSges(5,2) + t236 * rSges(5,3);
t176 = t237 * pkin(3) + qJ(4) * t236;
t331 = -t162 - t176;
t223 = t234 * qJ(4);
t175 = pkin(3) * t235 + t223;
t165 = t258 * t175;
t196 = t256 * pkin(4) + t230;
t330 = t258 * t196 + t165;
t218 = pkin(3) * t255 + qJ(4) * t254;
t329 = t175 * t340 + t256 * t218;
t220 = t259 * pkin(2) + pkin(9) * t258;
t217 = t273 * t220;
t328 = t273 * t176 + t217;
t219 = t257 * pkin(2) + t256 * pkin(9);
t327 = -t175 - t219;
t326 = -t176 - t197;
t198 = -Icges(5,5) * t340 - Icges(5,6) * t255 + Icges(5,3) * t254;
t199 = -Icges(5,4) * t340 - Icges(5,2) * t255 + Icges(5,6) * t254;
t325 = t254 * t198 - t255 * t199;
t202 = Icges(4,4) * t255 - Icges(4,2) * t254 - Icges(4,6) * t340;
t203 = Icges(4,1) * t255 - Icges(4,4) * t254 - Icges(4,5) * t340;
t324 = -t254 * t202 + t255 * t203;
t204 = -rSges(5,1) * t340 - rSges(5,2) * t255 + rSges(5,3) * t254;
t323 = -t204 - t218;
t322 = t219 * t341 + t220 * t339;
t240 = -pkin(4) * t340 + t255 * pkin(10);
t321 = -t218 - t240;
t320 = t280 * pkin(1) + pkin(8) * t341;
t111 = Icges(7,5) * t193 + Icges(7,6) * t192 + Icges(7,3) * t235;
t115 = Icges(7,4) * t193 + Icges(7,2) * t192 + Icges(7,6) * t235;
t119 = Icges(7,1) * t193 + Icges(7,4) * t192 + Icges(7,5) * t235;
t40 = t111 * t235 + t115 * t192 + t119 * t193;
t112 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t237;
t116 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t237;
t120 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t237;
t41 = t112 * t235 + t116 * t192 + t120 * t193;
t59 = t141 * t235 + t143 * t192 + t145 * t193;
t1 = t235 * t40 + t237 * t41 + t255 * t59;
t113 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t235;
t117 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t235;
t121 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t235;
t42 = t113 * t235 + t117 * t192 + t121 * t193;
t114 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t237;
t118 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t237;
t122 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t237;
t43 = t114 * t235 + t118 * t192 + t122 * t193;
t60 = t142 * t235 + t144 * t192 + t146 * t193;
t2 = t235 * t42 + t237 * t43 + t255 * t60;
t319 = t2 / 0.2e1 + t1 / 0.2e1;
t44 = t111 * t237 + t115 * t194 + t119 * t195;
t45 = t112 * t237 + t116 * t194 + t120 * t195;
t61 = t141 * t237 + t143 * t194 + t145 * t195;
t3 = t235 * t44 + t237 * t45 + t255 * t61;
t46 = t113 * t237 + t117 * t194 + t121 * t195;
t47 = t114 * t237 + t118 * t194 + t122 * t195;
t62 = t142 * t237 + t144 * t194 + t146 * t195;
t4 = t235 * t46 + t237 * t47 + t255 * t62;
t318 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t256 * t40 + t258 * t41 - t340 * t59;
t6 = t256 * t42 + t258 * t43 - t340 * t60;
t317 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t256 * t44 + t258 * t45 - t340 * t61;
t8 = t256 * t46 + t258 * t47 - t340 * t62;
t316 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t60 * t273 + (t277 * t43 - t280 * t42) * t272;
t9 = t59 * t273 + (t277 * t41 - t280 * t40) * t272;
t314 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t61 * t273 + (t277 * t45 - t280 * t44) * t272;
t12 = t62 * t273 + (t277 * t47 - t280 * t46) * t272;
t313 = t12 / 0.2e1 + t11 / 0.2e1;
t50 = t111 * t255 + t115 * t232 + t119 * t233;
t51 = t112 * t255 + t116 * t232 + t120 * t233;
t65 = t71 * t255;
t13 = t50 * t235 + t51 * t237 + t65;
t52 = t113 * t255 + t117 * t232 + t121 * t233;
t53 = t114 * t255 + t118 * t232 + t122 * t233;
t66 = t72 * t255;
t14 = t52 * t235 + t53 * t237 + t66;
t312 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t50 * t256 + t51 * t258 - t340 * t71;
t16 = t52 * t256 + t53 * t258 - t340 * t72;
t311 = t16 / 0.2e1 + t15 / 0.2e1;
t69 = t71 * t273;
t17 = t69 + (t51 * t277 - t50 * t280) * t272;
t70 = t72 * t273;
t18 = t70 + (t53 * t277 - t52 * t280) * t272;
t310 = t18 / 0.2e1 + t17 / 0.2e1;
t126 = t195 * rSges(6,1) + t194 * rSges(6,2) + t237 * rSges(6,3);
t309 = -t126 + t326;
t148 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t255;
t308 = -t148 + t321;
t307 = t273 * t197 + t328;
t306 = -t196 + t327;
t164 = t237 * rSges(4,1) - t236 * rSges(4,2) + t258 * rSges(4,3);
t243 = Icges(3,3) * t273 + (Icges(3,5) * t276 + Icges(3,6) * t279) * t272;
t244 = Icges(3,6) * t273 + (Icges(3,4) * t276 + Icges(3,2) * t279) * t272;
t245 = Icges(3,5) * t273 + (Icges(3,1) * t276 + Icges(3,4) * t279) * t272;
t304 = t272 * t276 * t245 + t273 * t243 + t244 * t340;
t213 = t259 * rSges(3,1) - t258 * rSges(3,2) + rSges(3,3) * t341;
t301 = -t277 * pkin(1) + pkin(8) * t339;
t205 = rSges(4,1) * t255 - rSges(4,2) * t254 - rSges(4,3) * t340;
t260 = (pkin(2) * t276 - pkin(9) * t279) * t272;
t300 = t272 * (-t205 - t260);
t299 = t326 - t333;
t298 = t321 - t332;
t297 = t175 * t341 + t176 * t339 + t322;
t296 = t196 * t340 + t256 * t240 + t329;
t295 = t272 * (-t260 + t323);
t294 = -t256 * rSges(5,1) - t234 * rSges(5,3);
t293 = -t193 * rSges(6,1) - t192 * rSges(6,2);
t291 = t220 + t320;
t290 = t272 * (-t260 + t308);
t289 = t52 / 0.2e1 + t50 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1;
t288 = t62 / 0.2e1 + t53 / 0.2e1 + t51 / 0.2e1 + t61 / 0.2e1;
t287 = t196 * t341 + t197 * t339 + t297;
t286 = t272 * (-t260 + t298);
t285 = -t219 + t301;
t163 = rSges(4,1) * t235 - rSges(4,2) * t234 + rSges(4,3) * t256;
t284 = -t223 + t285;
t212 = t257 * rSges(3,1) - t256 * rSges(3,2) - rSges(3,3) * t339;
t283 = t176 + t291;
t149 = Icges(5,5) * t256 - Icges(5,6) * t235 + Icges(5,3) * t234;
t153 = Icges(5,4) * t256 - Icges(5,2) * t235 + Icges(5,6) * t234;
t157 = Icges(5,1) * t256 - Icges(5,4) * t235 + Icges(5,5) * t234;
t85 = t149 * t254 - t153 * t255 - t157 * t340;
t151 = Icges(4,5) * t235 - Icges(4,6) * t234 + Icges(4,3) * t256;
t155 = Icges(4,4) * t235 - Icges(4,2) * t234 + Icges(4,6) * t256;
t159 = Icges(4,1) * t235 - Icges(4,4) * t234 + Icges(4,5) * t256;
t87 = -t151 * t340 - t155 * t254 + t159 * t255;
t200 = -Icges(5,1) * t340 - Icges(5,4) * t255 + Icges(5,5) * t254;
t94 = t198 * t234 - t199 * t235 + t200 * t256;
t201 = Icges(4,5) * t255 - Icges(4,6) * t254 - Icges(4,3) * t340;
t96 = t201 * t256 - t202 * t234 + t203 * t235;
t282 = t85 / 0.2e1 + t96 / 0.2e1 + t94 / 0.2e1 + t87 / 0.2e1 + t289;
t150 = Icges(5,5) * t258 - Icges(5,6) * t237 + Icges(5,3) * t236;
t154 = Icges(5,4) * t258 - Icges(5,2) * t237 + Icges(5,6) * t236;
t158 = Icges(5,1) * t258 - Icges(5,4) * t237 + Icges(5,5) * t236;
t86 = t150 * t254 - t154 * t255 - t158 * t340;
t152 = Icges(4,5) * t237 - Icges(4,6) * t236 + Icges(4,3) * t258;
t156 = Icges(4,4) * t237 - Icges(4,2) * t236 + Icges(4,6) * t258;
t160 = Icges(4,1) * t237 - Icges(4,4) * t236 + Icges(4,5) * t258;
t88 = -t152 * t340 - t156 * t254 + t160 * t255;
t95 = t198 * t236 - t199 * t237 + t200 * t258;
t97 = t201 * t258 - t202 * t236 + t203 * t237;
t281 = t86 / 0.2e1 + t97 / 0.2e1 + t95 / 0.2e1 + t88 / 0.2e1 + t288;
t262 = rSges(2,1) * t280 - t277 * rSges(2,2);
t261 = -t277 * rSges(2,1) - rSges(2,2) * t280;
t246 = t273 * rSges(3,3) + (rSges(3,1) * t276 + rSges(3,2) * t279) * t272;
t211 = Icges(3,1) * t259 - Icges(3,4) * t258 + Icges(3,5) * t341;
t210 = Icges(3,1) * t257 - Icges(3,4) * t256 - Icges(3,5) * t339;
t209 = Icges(3,4) * t259 - Icges(3,2) * t258 + Icges(3,6) * t341;
t208 = Icges(3,4) * t257 - Icges(3,2) * t256 - Icges(3,6) * t339;
t207 = Icges(3,5) * t259 - Icges(3,6) * t258 + Icges(3,3) * t341;
t185 = t213 + t320;
t184 = -t212 + t301;
t169 = -t273 * t212 - t246 * t339;
t168 = t213 * t273 - t246 * t341;
t161 = -rSges(5,2) * t235 - t294;
t140 = t304 * t273;
t137 = (t212 * t277 + t213 * t280) * t272;
t136 = t243 * t341 - t244 * t258 + t245 * t259;
t135 = -t243 * t339 - t256 * t244 + t257 * t245;
t130 = t291 + t164;
t129 = -t163 + t285;
t124 = rSges(6,3) * t235 - t293;
t110 = -t164 * t340 - t205 * t258;
t109 = t163 * t340 + t205 * t256;
t108 = t273 * t207 + (t209 * t279 + t211 * t276) * t272;
t107 = t273 * t206 + (t208 * t279 + t210 * t276) * t272;
t106 = -t201 * t340 + t324;
t105 = -t200 * t340 + t325;
t104 = t106 * t273;
t103 = t105 * t273;
t102 = t283 + t162;
t101 = (rSges(5,2) - pkin(3)) * t235 + t284 + t294;
t100 = t163 * t258 - t164 * t256;
t99 = (-t163 - t219) * t273 + t280 * t300;
t98 = t273 * t164 + t277 * t300 + t217;
t93 = (t163 * t277 + t164 * t280) * t272 + t322;
t92 = t126 * t255 - t148 * t237;
t91 = -t124 * t255 + t148 * t235;
t90 = t258 * t323 + t331 * t340;
t89 = t161 * t340 + t204 * t256 + t329;
t84 = t197 + t283 + t126;
t83 = (-rSges(6,3) - pkin(3)) * t235 + t284 + t293 - t196;
t82 = (-t161 + t327) * t273 + t280 * t295;
t81 = t273 * t162 + t277 * t295 + t328;
t80 = t152 * t258 - t156 * t236 + t160 * t237;
t79 = t151 * t258 - t155 * t236 + t159 * t237;
t78 = t152 * t256 - t156 * t234 + t160 * t235;
t77 = t151 * t256 - t155 * t234 + t159 * t235;
t76 = t150 * t236 - t154 * t237 + t158 * t258;
t75 = t149 * t236 - t153 * t237 + t157 * t258;
t74 = t150 * t234 - t154 * t235 + t158 * t256;
t73 = t149 * t234 - t153 * t235 + t157 * t256;
t68 = t283 + t350;
t67 = -t256 * t270 + (-pkin(3) + t352) * t235 + t284 + t349;
t64 = t124 * t237 - t126 * t235;
t63 = t258 * t161 + t256 * t331 + t165;
t58 = (t161 * t277 + t162 * t280) * t272 + t297;
t57 = t258 * t308 + t309 * t340;
t56 = t124 * t340 + t148 * t256 + t296;
t55 = (-t124 + t306) * t273 + t280 * t290;
t54 = t273 * t126 + t277 * t290 + t307;
t49 = -t237 * t332 + t255 * t333;
t48 = t235 * t332 - t255 * t334;
t39 = t258 * t124 + t256 * t309 + t330;
t38 = (t124 * t277 + t126 * t280) * t272 + t287;
t37 = -t235 * t333 + t237 * t334;
t36 = t258 * t298 + t299 * t340;
t35 = t256 * t332 + t334 * t340 + t296;
t34 = (t306 - t334) * t273 + t280 * t286;
t33 = t273 * t333 + t277 * t286 + t307;
t32 = t104 + (t88 * t277 - t87 * t280) * t272;
t31 = t103 + (t86 * t277 - t85 * t280) * t272;
t30 = -t106 * t340 + t87 * t256 + t88 * t258;
t29 = -t105 * t340 + t85 * t256 + t86 * t258;
t28 = t97 * t273 + (t277 * t80 - t280 * t79) * t272;
t27 = t96 * t273 + (t277 * t78 - t280 * t77) * t272;
t26 = t95 * t273 + (t277 * t76 - t280 * t75) * t272;
t25 = t94 * t273 + (t277 * t74 - t280 * t73) * t272;
t24 = t256 * t79 + t258 * t80 - t340 * t97;
t23 = t256 * t77 + t258 * t78 - t340 * t96;
t22 = t256 * t75 + t258 * t76 - t340 * t95;
t21 = t256 * t73 + t258 * t74 - t340 * t94;
t20 = t256 * t299 + t258 * t334 + t330;
t19 = (t277 * t334 + t280 * t333) * t272 + t287;
t123 = [(-t200 - t201) * t340 + m(7) * (t67 ^ 2 + t68 ^ 2) + m(6) * (t83 ^ 2 + t84 ^ 2) + m(5) * (t101 ^ 2 + t102 ^ 2) + m(4) * (t129 ^ 2 + t130 ^ 2) + m(3) * (t184 ^ 2 + t185 ^ 2) + m(2) * (t261 ^ 2 + t262 ^ 2) + t304 + Icges(2,3) + t324 + t325 - t351; t140 + t69 + t70 + t104 + t103 + m(7) * (t33 * t68 + t34 * t67) + m(6) * (t54 * t84 + t55 * t83) + m(5) * (t101 * t82 + t102 * t81) + m(4) * (t129 * t99 + t130 * t98) + m(3) * (t168 * t185 + t169 * t184) + ((-t107 / 0.2e1 - t135 / 0.2e1 - t282) * t280 + (t136 / 0.2e1 + t108 / 0.2e1 + t281) * t277) * t272; (t18 + t17 + t32 + t31 + t140) * t273 + m(7) * (t19 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t38 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t58 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t93 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(3) * (t137 ^ 2 + t168 ^ 2 + t169 ^ 2) + ((-t10 - t9 - t27 - t25 + (-t256 * t208 + t257 * t210 - t272 * t345) * t339) * t280 + (t11 + t12 + t28 + t26 + ((-t209 * t258 + t211 * t259 + (t207 * t277 - t345) * t272) * t277 + (t207 * t339 + t208 * t258 + t256 * t209 - t210 * t259 - t257 * t211) * t280) * t272) * t277 + ((-t107 - t135) * t280 + (t108 + t136) * t277) * t273) * t272; (-t105 - t106 + t351) * t340 + m(7) * (t35 * t67 + t36 * t68) + m(6) * (t56 * t83 + t57 * t84) + m(5) * (t101 * t89 + t102 * t90) + m(4) * (t109 * t129 + t110 * t130) + t281 * t258 + t282 * t256; (t29 / 0.2e1 + t30 / 0.2e1 + t311) * t273 + (t26 / 0.2e1 + t28 / 0.2e1 + t313) * t258 + (t25 / 0.2e1 + t27 / 0.2e1 + t314) * t256 + m(7) * (t19 * t20 + t33 * t36 + t34 * t35) + m(6) * (t38 * t39 + t54 * t57 + t55 * t56) + m(5) * (t63 * t58 + t81 * t90 + t82 * t89) + m(4) * (t100 * t93 + t109 * t99 + t110 * t98) + ((-t21 / 0.2e1 - t23 / 0.2e1 - t317) * t280 + (-t31 / 0.2e1 - t32 / 0.2e1 - t310) * t279 + (t22 / 0.2e1 + t24 / 0.2e1 + t316) * t277) * t272; (-t15 - t16 - t29 - t30) * t340 + (t22 + t7 + t8 + t24) * t258 + (t21 + t6 + t5 + t23) * t256 + m(4) * (t100 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(7) * (t20 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t39 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t63 ^ 2 + t89 ^ 2 + t90 ^ 2); m(7) * (t234 * t68 + t236 * t67) + m(6) * (t234 * t84 + t236 * t83) + m(5) * (t101 * t236 + t102 * t234); m(7) * (t19 * t254 + t234 * t33 + t236 * t34) + m(6) * (t234 * t54 + t236 * t55 + t254 * t38) + m(5) * (t234 * t81 + t236 * t82 + t254 * t58); m(7) * (t20 * t254 + t234 * t36 + t236 * t35) + m(6) * (t234 * t57 + t236 * t56 + t254 * t39) + m(5) * (t234 * t90 + t236 * t89 + t254 * t63); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t234 ^ 2 + t236 ^ 2 + t254 ^ 2); t66 + t65 + m(7) * (t48 * t67 + t49 * t68) + m(6) * (t83 * t91 + t84 * t92) + t288 * t237 + t289 * t235; t312 * t273 + t310 * t255 + t313 * t237 + t314 * t235 + m(7) * (t19 * t37 + t33 * t49 + t34 * t48) + m(6) * (t64 * t38 + t54 * t92 + t55 * t91) + (t277 * t318 - t280 * t319) * t272; -t312 * t340 + t318 * t258 + t319 * t256 + t311 * t255 + t316 * t237 + t317 * t235 + m(7) * (t20 * t37 + t35 * t48 + t36 * t49) + m(6) * (t64 * t39 + t56 * t91 + t57 * t92); m(6) * (t234 * t92 + t236 * t91 + t254 * t64) + m(7) * (t234 * t49 + t236 * t48 + t254 * t37); (t13 + t14) * t255 + (t4 + t3) * t237 + (t1 + t2) * t235 + m(7) * (t37 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t64 ^ 2 + t91 ^ 2 + t92 ^ 2); m(7) * (t235 * t68 + t237 * t67); m(7) * (t19 * t255 + t235 * t33 + t237 * t34); m(7) * (t20 * t255 + t235 * t36 + t237 * t35); m(7) * (t234 * t235 + t236 * t237 + t254 * t255); m(7) * (t235 * t49 + t237 * t48 + t255 * t37); m(7) * (t235 ^ 2 + t237 ^ 2 + t255 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t123(1) t123(2) t123(4) t123(7) t123(11) t123(16); t123(2) t123(3) t123(5) t123(8) t123(12) t123(17); t123(4) t123(5) t123(6) t123(9) t123(13) t123(18); t123(7) t123(8) t123(9) t123(10) t123(14) t123(19); t123(11) t123(12) t123(13) t123(14) t123(15) t123(20); t123(16) t123(17) t123(18) t123(19) t123(20) t123(21);];
Mq  = res;
