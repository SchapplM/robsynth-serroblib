% Calculate joint inertia matrix for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:19
% EndTime: 2019-03-09 11:17:30
% DurationCPUTime: 5.24s
% Computational Cost: add. (16311->609), mult. (29354->841), div. (0->0), fcn. (36326->12), ass. (0->274)
t261 = sin(pkin(6));
t262 = cos(pkin(6));
t266 = sin(qJ(2));
t270 = cos(qJ(2));
t216 = Icges(3,3) * t262 + (Icges(3,5) * t266 + Icges(3,6) * t270) * t261;
t217 = Icges(3,6) * t262 + (Icges(3,4) * t266 + Icges(3,2) * t270) * t261;
t218 = Icges(3,5) * t262 + (Icges(3,1) * t266 + Icges(3,4) * t270) * t261;
t219 = Icges(4,5) * t262 + (-Icges(4,6) * t266 - Icges(4,3) * t270) * t261;
t220 = Icges(4,4) * t262 + (-Icges(4,2) * t266 - Icges(4,6) * t270) * t261;
t221 = Icges(4,1) * t262 + (-Icges(4,4) * t266 - Icges(4,5) * t270) * t261;
t318 = t261 * t270;
t320 = t261 * t266;
t337 = (-t219 * t270 - t220 * t266) * t261 + t217 * t318 + t218 * t320 + (t221 + t216) * t262;
t300 = m(6) / 0.2e1 + m(7) / 0.2e1;
t336 = 0.2e1 * t300;
t271 = cos(qJ(1));
t312 = t270 * t271;
t267 = sin(qJ(1));
t315 = t266 * t267;
t239 = -t262 * t312 + t315;
t313 = t267 * t270;
t314 = t266 * t271;
t240 = t262 * t314 + t313;
t317 = t261 * t271;
t167 = -Icges(4,5) * t317 - Icges(4,6) * t240 + Icges(4,3) * t239;
t174 = Icges(3,4) * t240 - Icges(3,2) * t239 - Icges(3,6) * t317;
t335 = t167 - t174;
t169 = -Icges(4,4) * t317 - Icges(4,2) * t240 + Icges(4,6) * t239;
t176 = Icges(3,1) * t240 - Icges(3,4) * t239 - Icges(3,5) * t317;
t334 = t169 - t176;
t241 = t262 * t313 + t314;
t242 = -t262 * t315 + t312;
t319 = t261 * t267;
t166 = Icges(4,5) * t319 - Icges(4,6) * t242 + Icges(4,3) * t241;
t175 = Icges(3,4) * t242 - Icges(3,2) * t241 + Icges(3,6) * t319;
t333 = -t175 + t166;
t168 = Icges(4,4) * t319 - Icges(4,2) * t242 + Icges(4,6) * t241;
t177 = Icges(3,1) * t242 - Icges(3,4) * t241 + Icges(3,5) * t319;
t332 = t177 - t168;
t260 = t261 ^ 2;
t301 = qJ(4) + pkin(11);
t258 = sin(t301);
t291 = cos(t301);
t196 = -t241 * t291 + t258 * t319;
t198 = t239 * t291 + t258 * t317;
t282 = t261 * t291;
t222 = t262 * t258 + t270 * t282;
t197 = t241 * t258 + t267 * t282;
t264 = sin(qJ(6));
t268 = cos(qJ(6));
t150 = -t197 * t264 + t242 * t268;
t151 = t197 * t268 + t242 * t264;
t80 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t196;
t82 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t196;
t84 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t196;
t26 = t150 * t82 + t151 * t84 + t196 * t80;
t199 = t239 * t258 - t271 * t282;
t152 = -t199 * t264 + t240 * t268;
t153 = t199 * t268 + t240 * t264;
t81 = Icges(7,5) * t153 + Icges(7,6) * t152 - Icges(7,3) * t198;
t83 = Icges(7,4) * t153 + Icges(7,2) * t152 - Icges(7,6) * t198;
t85 = Icges(7,1) * t153 + Icges(7,4) * t152 - Icges(7,5) * t198;
t27 = t150 * t83 + t151 * t85 + t196 * t81;
t223 = -t258 * t318 + t262 * t291;
t194 = -t223 * t264 + t268 * t320;
t195 = t223 * t268 + t264 * t320;
t106 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t222;
t107 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t222;
t108 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t222;
t35 = t106 * t196 + t107 * t150 + t108 * t151;
t1 = t196 * t26 - t198 * t27 + t222 * t35;
t331 = t1 / 0.2e1;
t30 = t194 * t82 + t195 * t84 + t222 * t80;
t31 = t194 * t83 + t195 * t85 + t222 * t81;
t43 = t222 * t106 + t194 * t107 + t195 * t108;
t38 = t43 * t222;
t7 = t30 * t196 - t31 * t198 + t38;
t330 = t7 / 0.2e1;
t329 = t196 / 0.2e1;
t328 = -t198 / 0.2e1;
t327 = t222 / 0.2e1;
t326 = t199 * pkin(5);
t263 = -qJ(5) - pkin(9);
t325 = -pkin(2) + t263;
t86 = t151 * rSges(7,1) + t150 * rSges(7,2) + t196 * rSges(7,3);
t324 = t197 * pkin(5) + pkin(10) * t196 + t86;
t283 = -t153 * rSges(7,1) - t152 * rSges(7,2);
t87 = -t198 * rSges(7,3) - t283;
t323 = -t198 * pkin(10) + t326 + t87;
t265 = sin(qJ(4));
t322 = t239 * t265;
t321 = t241 * t265;
t316 = t265 * t270;
t109 = rSges(7,1) * t195 + rSges(7,2) * t194 + rSges(7,3) * t222;
t311 = pkin(5) * t223 + pkin(10) * t222 + t109;
t310 = t337 * t262;
t171 = -Icges(4,1) * t317 - Icges(4,4) * t240 + Icges(4,5) * t239;
t172 = Icges(3,5) * t240 - Icges(3,6) * t239 - Icges(3,3) * t317;
t309 = -t172 - t171;
t170 = Icges(4,1) * t319 - Icges(4,4) * t242 + Icges(4,5) * t241;
t173 = Icges(3,5) * t242 - Icges(3,6) * t241 + Icges(3,3) * t319;
t308 = t173 + t170;
t228 = t239 * qJ(3);
t188 = t240 * pkin(2) + t228;
t189 = t242 * pkin(2) + qJ(3) * t241;
t307 = t188 * t319 + t189 * t317;
t187 = t262 * t189;
t211 = pkin(3) * t319 + pkin(9) * t242;
t306 = t262 * t211 + t187;
t302 = pkin(3) * t317 - t240 * pkin(9);
t305 = -t188 + t302;
t243 = (pkin(2) * t266 - qJ(3) * t270) * t261;
t304 = -t262 * pkin(3) - pkin(9) * t320 - t243;
t303 = t271 * pkin(1) + pkin(8) * t319;
t299 = t35 / 0.2e1 + t30 / 0.2e1;
t36 = -t106 * t198 + t107 * t152 + t108 * t153;
t298 = -t36 / 0.2e1 - t31 / 0.2e1;
t269 = cos(qJ(4));
t257 = pkin(4) * t269 + pkin(3);
t293 = pkin(4) * t321 - t242 * t263 + t257 * t319;
t135 = -t211 + t293;
t297 = t262 * t135 + t306;
t288 = -pkin(4) * t322 + t257 * t317;
t136 = -t240 * t263 - t288 + t302;
t296 = -t136 + t305;
t156 = Icges(6,5) * t223 - Icges(6,6) * t222 + Icges(6,3) * t320;
t157 = Icges(6,4) * t223 - Icges(6,2) * t222 + Icges(6,6) * t320;
t158 = Icges(6,1) * t223 - Icges(6,4) * t222 + Icges(6,5) * t320;
t70 = t156 * t320 - t222 * t157 + t223 * t158;
t237 = -t262 * t265 - t269 * t318;
t238 = -t261 * t316 + t262 * t269;
t162 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t320;
t163 = Icges(5,4) * t238 + Icges(5,2) * t237 + Icges(5,6) * t320;
t164 = Icges(5,1) * t238 + Icges(5,4) * t237 + Icges(5,5) * t320;
t79 = t162 * t320 + t237 * t163 + t238 * t164;
t184 = (-pkin(3) + t257) * t262 + (-pkin(4) * t316 + (-pkin(9) - t263) * t266) * t261;
t295 = -t184 + t304;
t119 = t197 * rSges(6,1) - t196 * rSges(6,2) + t242 * rSges(6,3);
t205 = t241 * t269 - t265 * t319;
t206 = t269 * t319 + t321;
t133 = t206 * rSges(5,1) + t205 * rSges(5,2) + t242 * rSges(5,3);
t181 = t242 * rSges(3,1) - t241 * rSges(3,2) + rSges(3,3) * t319;
t178 = rSges(4,1) * t319 - t242 * rSges(4,2) + t241 * rSges(4,3);
t292 = -t267 * pkin(1) + pkin(8) * t317;
t290 = t261 * (-t262 * rSges(4,1) - (-rSges(4,2) * t266 - rSges(4,3) * t270) * t261 - t243);
t289 = t211 * t317 - t302 * t319 + t307;
t287 = -t228 + t292;
t165 = rSges(5,1) * t238 + rSges(5,2) * t237 + rSges(5,3) * t320;
t286 = t261 * (-t165 + t304);
t207 = t239 * t269 + t265 * t317;
t208 = -t269 * t317 + t322;
t285 = -t208 * rSges(5,1) - t207 * rSges(5,2);
t284 = -t199 * rSges(6,1) - t198 * rSges(6,2);
t159 = rSges(6,1) * t223 - rSges(6,2) * t222 + rSges(6,3) * t320;
t281 = t261 * (-t159 + t295);
t280 = t189 + t303;
t279 = rSges(4,1) * t317 - t239 * rSges(4,3);
t278 = t135 * t317 + t136 * t319 + t289;
t277 = t261 * (t295 - t311);
t180 = t240 * rSges(3,1) - t239 * rSges(3,2) - rSges(3,3) * t317;
t113 = Icges(6,5) * t197 - Icges(6,6) * t196 + Icges(6,3) * t242;
t115 = Icges(6,4) * t197 - Icges(6,2) * t196 + Icges(6,6) * t242;
t117 = Icges(6,1) * t197 - Icges(6,4) * t196 + Icges(6,5) * t242;
t52 = t113 * t320 - t115 * t222 + t117 * t223;
t127 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t242;
t129 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t242;
t131 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t242;
t62 = t127 * t320 + t129 * t237 + t131 * t238;
t64 = t156 * t242 - t157 * t196 + t158 * t197;
t71 = t162 * t242 + t163 * t205 + t164 * t206;
t275 = t71 / 0.2e1 + t64 / 0.2e1 + t62 / 0.2e1 + t52 / 0.2e1 + t299;
t114 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t240;
t116 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t240;
t118 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t240;
t53 = t114 * t320 - t116 * t222 + t118 * t223;
t128 = Icges(5,5) * t208 + Icges(5,6) * t207 + Icges(5,3) * t240;
t130 = Icges(5,4) * t208 + Icges(5,2) * t207 + Icges(5,6) * t240;
t132 = Icges(5,1) * t208 + Icges(5,4) * t207 + Icges(5,5) * t240;
t63 = t128 * t320 + t130 * t237 + t132 * t238;
t65 = t156 * t240 + t157 * t198 + t158 * t199;
t72 = t162 * t240 + t163 * t207 + t164 * t208;
t274 = t72 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1 + t53 / 0.2e1 - t298;
t273 = t287 + t288;
t272 = t280 + t293;
t248 = rSges(2,1) * t271 - t267 * rSges(2,2);
t247 = -t267 * rSges(2,1) - rSges(2,2) * t271;
t224 = t262 * rSges(3,3) + (rSges(3,1) * t266 + rSges(3,2) * t270) * t261;
t179 = -t240 * rSges(4,2) - t279;
t155 = t181 + t303;
t154 = -t180 + t292;
t149 = t240 * t184;
t140 = -t262 * t180 - t224 * t317;
t139 = t181 * t262 - t224 * t319;
t134 = rSges(5,3) * t240 - t285;
t124 = t135 * t320;
t120 = rSges(6,3) * t240 - t284;
t112 = t280 + t178;
t111 = (rSges(4,2) - pkin(2)) * t240 + t279 + t287;
t110 = t242 * t136;
t105 = (t180 * t267 + t181 * t271) * t261;
t104 = t216 * t319 - t217 * t241 + t218 * t242;
t103 = -t216 * t317 - t239 * t217 + t240 * t218;
t102 = t239 * t219 - t240 * t220 - t221 * t317;
t101 = t219 * t241 - t220 * t242 + t221 * t319;
t97 = (-t179 - t188) * t262 + t271 * t290;
t96 = t262 * t178 + t267 * t290 + t187;
t95 = t133 * t320 - t165 * t242;
t94 = -t134 * t320 + t165 * t240;
t93 = t262 * t171 + (-t167 * t270 - t169 * t266) * t261;
t92 = t262 * t170 + (-t166 * t270 - t168 * t266) * t261;
t91 = t262 * t173 + (t175 * t270 + t177 * t266) * t261;
t90 = t262 * t172 + (t174 * t270 + t176 * t266) * t261;
t89 = t211 + t280 + t133;
t88 = (-rSges(5,3) - pkin(2)) * t240 + t285 + t287 + t302;
t78 = t79 * t262;
t77 = t79 * t320;
t76 = (t178 * t271 + t179 * t267) * t261 + t307;
t75 = t272 + t119;
t74 = (-rSges(6,3) + t325) * t240 + t273 + t284;
t73 = -t133 * t240 + t134 * t242;
t69 = t70 * t262;
t68 = t70 * t320;
t67 = (-t134 + t305) * t262 + t271 * t286;
t66 = t262 * t133 + t267 * t286 + t306;
t61 = -t109 * t198 - t222 * t87;
t60 = -t109 * t196 + t222 * t86;
t59 = t119 * t320 + t124 + (-t159 - t184) * t242;
t58 = t240 * t159 + t149 + (-t120 - t136) * t320;
t57 = t128 * t240 + t130 * t207 + t132 * t208;
t56 = t127 * t240 + t129 * t207 + t131 * t208;
t55 = t128 * t242 + t130 * t205 + t132 * t206;
t54 = t127 * t242 + t129 * t205 + t131 * t206;
t51 = (t133 * t271 + t134 * t267) * t261 + t289;
t50 = t272 + t324;
t49 = -t326 + t325 * t240 + (rSges(7,3) + pkin(10)) * t198 + t273 + t283;
t48 = t114 * t240 + t116 * t198 + t118 * t199;
t47 = t113 * t240 + t115 * t198 + t117 * t199;
t46 = t114 * t242 - t116 * t196 + t118 * t197;
t45 = t113 * t242 - t115 * t196 + t117 * t197;
t44 = t196 * t87 + t198 * t86;
t42 = t43 * t262;
t41 = t43 * t320;
t40 = (-t120 + t296) * t262 + t271 * t281;
t39 = t262 * t119 + t267 * t281 + t297;
t37 = t242 * t120 + t110 + (-t119 - t135) * t240;
t34 = (t119 * t271 + t120 * t267) * t261 + t278;
t33 = t124 + t324 * t320 + (-t184 - t311) * t242;
t32 = t149 + t311 * t240 + (-t136 - t323) * t320;
t29 = t152 * t83 + t153 * t85 - t198 * t81;
t28 = t152 * t82 + t153 * t84 - t198 * t80;
t25 = (t296 - t323) * t262 + t271 * t277;
t24 = t262 * t324 + t267 * t277 + t297;
t23 = t110 + t323 * t242 + (-t135 - t324) * t240;
t22 = t78 + (t62 * t267 - t63 * t271) * t261;
t21 = t63 * t240 + t62 * t242 + t77;
t20 = (t267 * t323 + t271 * t324) * t261 + t278;
t19 = t72 * t262 + (t267 * t56 - t271 * t57) * t261;
t18 = t71 * t262 + (t267 * t54 - t271 * t55) * t261;
t17 = t69 + (t52 * t267 - t53 * t271) * t261;
t16 = t240 * t57 + t242 * t56 + t320 * t72;
t15 = t240 * t55 + t242 * t54 + t320 * t71;
t14 = t53 * t240 + t52 * t242 + t68;
t13 = t65 * t262 + (t267 * t47 - t271 * t48) * t261;
t12 = t64 * t262 + (t267 * t45 - t271 * t46) * t261;
t11 = t240 * t48 + t242 * t47 + t320 * t65;
t10 = t240 * t46 + t242 * t45 + t320 * t64;
t9 = t42 + (t30 * t267 - t31 * t271) * t261;
t8 = t31 * t240 + t30 * t242 + t41;
t6 = t36 * t262 + (t267 * t28 - t271 * t29) * t261;
t5 = t35 * t262 + (t26 * t267 - t27 * t271) * t261;
t4 = t240 * t29 + t242 * t28 + t320 * t36;
t3 = t240 * t27 + t242 * t26 + t320 * t35;
t2 = t196 * t28 - t198 * t29 + t222 * t36;
t98 = [Icges(2,3) + m(7) * (t49 ^ 2 + t50 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t88 ^ 2 + t89 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2) + m(3) * (t154 ^ 2 + t155 ^ 2) + m(2) * (t247 ^ 2 + t248 ^ 2) + t79 + t70 + t43 + t337; t42 + t69 + t78 + m(7) * (t24 * t50 + t25 * t49) + m(6) * (t39 * t75 + t40 * t74) + m(5) * (t66 * t89 + t67 * t88) + m(4) * (t111 * t97 + t112 * t96) + m(3) * (t139 * t155 + t140 * t154) + ((-t90 / 0.2e1 - t93 / 0.2e1 - t102 / 0.2e1 - t103 / 0.2e1 - t274) * t271 + (t91 / 0.2e1 + t92 / 0.2e1 + t101 / 0.2e1 + t104 / 0.2e1 + t275) * t267) * t261 + t310; (t9 + t22 + t17 + t310) * t262 + m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t34 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t51 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(4) * (t76 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(3) * (t105 ^ 2 + t139 ^ 2 + t140 ^ 2) + ((-t13 - t19 - t6 + ((t335 * t239 - t334 * t240) * t261 + t309 * t260 * t271) * t271 + (-t102 - t103 - t90 - t93) * t262) * t271 + (t5 + t18 + t12 + ((t333 * t241 + t332 * t242) * t261 + t308 * t260 * t267) * t267 + (t104 + t101 + t91 + t92) * t262 + ((t267 * t309 + t271 * t308) * t261 + t334 * t242 - t335 * t241 - t332 * t240 - t333 * t239) * t317) * t267) * t261; m(7) * (t239 * t50 + t241 * t49) + m(6) * (t239 * t75 + t241 * t74) + m(5) * (t239 * t89 + t241 * t88) + m(4) * (t111 * t241 + t112 * t239); m(7) * (-t20 * t318 + t239 * t24 + t241 * t25) + m(6) * (t239 * t39 + t241 * t40 - t318 * t34) + m(5) * (t239 * t66 + t241 * t67 - t318 * t51) + m(4) * (t239 * t96 + t241 * t97 - t318 * t76); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t300) * (t260 * t270 ^ 2 + t239 ^ 2 + t241 ^ 2); t41 + t77 + t68 + m(7) * (t32 * t49 + t33 * t50) + m(6) * (t58 * t74 + t59 * t75) + m(5) * (t88 * t94 + t89 * t95) + t275 * t242 + t274 * t240; (t8 / 0.2e1 + t21 / 0.2e1 + t14 / 0.2e1) * t262 + (t5 / 0.2e1 + t18 / 0.2e1 + t12 / 0.2e1) * t242 + (t6 / 0.2e1 + t19 / 0.2e1 + t13 / 0.2e1) * t240 + m(7) * (t20 * t23 + t24 * t33 + t25 * t32) + m(6) * (t34 * t37 + t39 * t59 + t40 * t58) + m(5) * (t51 * t73 + t66 * t95 + t67 * t94) + ((-t4 / 0.2e1 - t16 / 0.2e1 - t11 / 0.2e1) * t271 + (t3 / 0.2e1 + t15 / 0.2e1 + t10 / 0.2e1) * t267 + (t9 / 0.2e1 + t22 / 0.2e1 + t17 / 0.2e1) * t266) * t261; m(5) * (t239 * t95 + t241 * t94 - t318 * t73) + m(6) * (t239 * t59 + t241 * t58 - t318 * t37) + m(7) * (-t23 * t318 + t239 * t33 + t241 * t32); (t14 + t21 + t8) * t320 + (t3 + t10 + t15) * t242 + (t4 + t16 + t11) * t240 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t37 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t73 ^ 2 + t94 ^ 2 + t95 ^ 2); m(7) * (t240 * t50 + t242 * t49) + m(6) * (t240 * t75 + t242 * t74); m(7) * (t20 * t320 + t24 * t240 + t242 * t25) + m(6) * (t240 * t39 + t242 * t40 + t320 * t34); (-t260 * t266 * t270 + t239 * t240 + t241 * t242) * t336; m(7) * (t23 * t320 + t240 * t33 + t242 * t32) + m(6) * (t240 * t59 + t242 * t58 + t320 * t37); (t260 * t266 ^ 2 + t240 ^ 2 + t242 ^ 2) * t336; m(7) * (t49 * t61 + t50 * t60) + t38 + t298 * t198 + t299 * t196; t262 * t330 + m(7) * (t20 * t44 + t24 * t60 + t25 * t61) + t9 * t327 + t6 * t328 + t5 * t329 + (t267 * t331 - t271 * t2 / 0.2e1) * t261; m(7) * (t239 * t60 + t241 * t61 - t318 * t44); t320 * t330 + t240 * t2 / 0.2e1 + t4 * t328 + m(7) * (t23 * t44 + t32 * t61 + t33 * t60) + t242 * t331 + t8 * t327 + t3 * t329; m(7) * (t240 * t60 + t242 * t61 + t320 * t44); t196 * t1 - t198 * t2 + t222 * t7 + m(7) * (t44 ^ 2 + t60 ^ 2 + t61 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t98(1) t98(2) t98(4) t98(7) t98(11) t98(16); t98(2) t98(3) t98(5) t98(8) t98(12) t98(17); t98(4) t98(5) t98(6) t98(9) t98(13) t98(18); t98(7) t98(8) t98(9) t98(10) t98(14) t98(19); t98(11) t98(12) t98(13) t98(14) t98(15) t98(20); t98(16) t98(17) t98(18) t98(19) t98(20) t98(21);];
Mq  = res;
