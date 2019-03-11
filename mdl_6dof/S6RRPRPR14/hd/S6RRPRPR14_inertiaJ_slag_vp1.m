% Calculate joint inertia matrix for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR14_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR14_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:37
% EndTime: 2019-03-09 11:32:48
% DurationCPUTime: 5.04s
% Computational Cost: add. (11356->596), mult. (28988->823), div. (0->0), fcn. (36149->10), ass. (0->260)
t253 = sin(pkin(6));
t254 = cos(pkin(6));
t256 = sin(qJ(2));
t259 = cos(qJ(2));
t213 = Icges(3,3) * t254 + (Icges(3,5) * t256 + Icges(3,6) * t259) * t253;
t214 = Icges(3,6) * t254 + (Icges(3,4) * t256 + Icges(3,2) * t259) * t253;
t215 = Icges(3,5) * t254 + (Icges(3,1) * t256 + Icges(3,4) * t259) * t253;
t216 = Icges(4,5) * t254 + (-Icges(4,6) * t256 - Icges(4,3) * t259) * t253;
t217 = Icges(4,4) * t254 + (-Icges(4,2) * t256 - Icges(4,6) * t259) * t253;
t218 = Icges(4,1) * t254 + (-Icges(4,4) * t256 - Icges(4,5) * t259) * t253;
t303 = t253 * t259;
t305 = t253 * t256;
t322 = (-t259 * t216 - t256 * t217) * t253 + t214 * t303 + t215 * t305 + (t218 + t213) * t254;
t288 = m(6) / 0.2e1 + m(7) / 0.2e1;
t321 = 0.2e1 * t288;
t260 = cos(qJ(1));
t298 = t259 * t260;
t257 = sin(qJ(1));
t301 = t256 * t257;
t235 = -t254 * t298 + t301;
t299 = t257 * t259;
t300 = t256 * t260;
t236 = t254 * t300 + t299;
t302 = t253 * t260;
t167 = -Icges(4,5) * t302 - Icges(4,6) * t236 + Icges(4,3) * t235;
t174 = Icges(3,4) * t236 - Icges(3,2) * t235 - Icges(3,6) * t302;
t320 = t167 - t174;
t169 = -Icges(4,4) * t302 - Icges(4,2) * t236 + Icges(4,6) * t235;
t176 = Icges(3,1) * t236 - Icges(3,4) * t235 - Icges(3,5) * t302;
t319 = t169 - t176;
t237 = t254 * t299 + t300;
t238 = -t254 * t301 + t298;
t304 = t253 * t257;
t166 = Icges(4,5) * t304 - Icges(4,6) * t238 + Icges(4,3) * t237;
t175 = Icges(3,4) * t238 - Icges(3,2) * t237 + Icges(3,6) * t304;
t318 = -t175 + t166;
t168 = Icges(4,4) * t304 - Icges(4,2) * t238 + Icges(4,6) * t237;
t177 = Icges(3,1) * t238 - Icges(3,4) * t237 + Icges(3,5) * t304;
t317 = t177 - t168;
t316 = t253 ^ 2;
t310 = cos(qJ(4));
t281 = t253 * t310;
t309 = sin(qJ(4));
t202 = t237 * t309 + t257 * t281;
t205 = -t235 * t309 + t260 * t281;
t280 = t253 * t309;
t234 = t254 * t310 - t259 * t280;
t201 = -t237 * t310 + t257 * t280;
t255 = sin(qJ(6));
t258 = cos(qJ(6));
t152 = t201 * t258 - t238 * t255;
t153 = t201 * t255 + t238 * t258;
t203 = t235 * t310 + t260 * t280;
t150 = -t203 * t258 - t236 * t255;
t151 = -t203 * t255 + t236 * t258;
t88 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t205;
t90 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t205;
t92 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t205;
t26 = t152 * t90 + t153 * t92 + t202 * t88;
t89 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t202;
t91 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t202;
t93 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t202;
t27 = t152 * t91 + t153 * t93 + t202 * t89;
t233 = t254 * t309 + t259 * t281;
t198 = t233 * t258 - t255 * t305;
t199 = t233 * t255 + t258 * t305;
t111 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t234;
t112 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t234;
t113 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t234;
t36 = t111 * t202 + t112 * t152 + t113 * t153;
t2 = t202 * t27 - t205 * t26 + t234 * t36;
t315 = t2 / 0.2e1;
t30 = t198 * t90 + t199 * t92 + t234 * t88;
t31 = t198 * t91 + t199 * t93 + t234 * t89;
t44 = t234 * t111 + t198 * t112 + t199 * t113;
t41 = t44 * t234;
t7 = t31 * t202 - t30 * t205 + t41;
t314 = t7 / 0.2e1;
t313 = t202 / 0.2e1;
t312 = -t205 / 0.2e1;
t311 = t234 / 0.2e1;
t308 = rSges(6,3) * t203;
t273 = -rSges(7,1) * t151 - rSges(7,2) * t150;
t94 = -rSges(7,3) * t205 - t273;
t307 = pkin(5) * t236 - pkin(10) * t205 + t94;
t95 = t153 * rSges(7,1) + t152 * rSges(7,2) + t202 * rSges(7,3);
t306 = t238 * pkin(5) + pkin(10) * t202 + t95;
t297 = t322 * t254;
t114 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t234;
t296 = pkin(5) * t305 + pkin(10) * t234 + t114;
t171 = -Icges(4,1) * t302 - Icges(4,4) * t236 + Icges(4,5) * t235;
t172 = Icges(3,5) * t236 - Icges(3,6) * t235 - Icges(3,3) * t302;
t295 = -t172 - t171;
t170 = Icges(4,1) * t304 - Icges(4,4) * t238 + Icges(4,5) * t237;
t173 = Icges(3,5) * t238 - Icges(3,6) * t237 + Icges(3,3) * t304;
t294 = t173 + t170;
t223 = t235 * qJ(3);
t187 = t236 * pkin(2) + t223;
t188 = t238 * pkin(2) + qJ(3) * t237;
t293 = t187 * t304 + t188 * t302;
t185 = t254 * t188;
t209 = pkin(3) * t304 + pkin(9) * t238;
t292 = t254 * t209 + t185;
t210 = -pkin(3) * t302 + t236 * pkin(9);
t291 = -t187 - t210;
t239 = (pkin(2) * t256 - qJ(3) * t259) * t253;
t290 = -pkin(3) * t254 - pkin(9) * t305 - t239;
t289 = t260 * pkin(1) + pkin(8) * t304;
t35 = -t111 * t205 + t112 * t150 + t113 * t151;
t287 = -t30 / 0.2e1 - t35 / 0.2e1;
t286 = t31 / 0.2e1 + t36 / 0.2e1;
t139 = t202 * pkin(4) + qJ(5) * t201;
t285 = t254 * t139 + t292;
t191 = t203 * qJ(5);
t140 = -pkin(4) * t205 - t191;
t284 = -t140 + t291;
t161 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t305;
t162 = Icges(5,4) * t234 - Icges(5,2) * t233 + Icges(5,6) * t305;
t163 = Icges(5,1) * t234 - Icges(5,4) * t233 + Icges(5,5) * t305;
t79 = t161 * t305 - t233 * t162 + t234 * t163;
t158 = Icges(6,5) * t305 - Icges(6,6) * t234 + Icges(6,3) * t233;
t159 = Icges(6,4) * t305 - Icges(6,2) * t234 + Icges(6,6) * t233;
t160 = Icges(6,1) * t305 - Icges(6,4) * t234 + Icges(6,5) * t233;
t78 = t233 * t158 - t234 * t159 + t160 * t305;
t186 = pkin(4) * t234 + qJ(5) * t233;
t283 = -t186 + t290;
t127 = t202 * rSges(5,1) - t201 * rSges(5,2) + t238 * rSges(5,3);
t130 = t238 * rSges(6,1) - t202 * rSges(6,2) + t201 * rSges(6,3);
t181 = t238 * rSges(3,1) - t237 * rSges(3,2) + rSges(3,3) * t304;
t178 = rSges(4,1) * t304 - t238 * rSges(4,2) + t237 * rSges(4,3);
t279 = -t257 * pkin(1) + pkin(8) * t302;
t278 = t253 * (-rSges(4,1) * t254 - (-rSges(4,2) * t256 - rSges(4,3) * t259) * t253 - t239);
t277 = t209 * t302 + t210 * t304 + t293;
t276 = -t223 + t279;
t165 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t305;
t275 = t253 * (-t165 + t290);
t274 = rSges(5,1) * t205 - rSges(5,2) * t203;
t164 = rSges(6,1) * t305 - rSges(6,2) * t234 + rSges(6,3) * t233;
t272 = t253 * (-t164 + t283);
t271 = t188 + t289;
t270 = rSges(4,1) * t302 - t235 * rSges(4,3);
t269 = t139 * t302 + t140 * t304 + t277;
t268 = t276 - t210;
t267 = t253 * (t283 - t296);
t266 = t191 + t268;
t180 = t236 * rSges(3,1) - t235 * rSges(3,2) - rSges(3,3) * t302;
t115 = Icges(6,5) * t236 + Icges(6,6) * t205 - Icges(6,3) * t203;
t119 = Icges(6,4) * t236 + Icges(6,2) * t205 - Icges(6,6) * t203;
t123 = Icges(6,1) * t236 + Icges(6,4) * t205 - Icges(6,5) * t203;
t56 = t115 * t233 - t119 * t234 + t123 * t305;
t118 = -Icges(5,5) * t205 + Icges(5,6) * t203 + Icges(5,3) * t236;
t122 = -Icges(5,4) * t205 + Icges(5,2) * t203 + Icges(5,6) * t236;
t126 = -Icges(5,1) * t205 + Icges(5,4) * t203 + Icges(5,5) * t236;
t59 = t118 * t305 - t122 * t233 + t126 * t234;
t69 = t161 * t236 + t162 * t203 - t163 * t205;
t70 = -t158 * t203 + t159 * t205 + t160 * t236;
t264 = t59 / 0.2e1 + t56 / 0.2e1 + t69 / 0.2e1 + t70 / 0.2e1 - t287;
t116 = Icges(6,5) * t238 - Icges(6,6) * t202 + Icges(6,3) * t201;
t120 = Icges(6,4) * t238 - Icges(6,2) * t202 + Icges(6,6) * t201;
t124 = Icges(6,1) * t238 - Icges(6,4) * t202 + Icges(6,5) * t201;
t57 = t116 * t233 - t120 * t234 + t124 * t305;
t117 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t238;
t121 = Icges(5,4) * t202 - Icges(5,2) * t201 + Icges(5,6) * t238;
t125 = Icges(5,1) * t202 - Icges(5,4) * t201 + Icges(5,5) * t238;
t58 = t117 * t305 - t121 * t233 + t125 * t234;
t68 = t161 * t238 - t162 * t201 + t163 * t202;
t71 = t158 * t201 - t159 * t202 + t160 * t238;
t263 = t68 / 0.2e1 + t71 / 0.2e1 + t58 / 0.2e1 + t57 / 0.2e1 + t286;
t262 = t209 + t271;
t261 = t139 + t262;
t242 = rSges(2,1) * t260 - t257 * rSges(2,2);
t241 = -t257 * rSges(2,1) - rSges(2,2) * t260;
t219 = rSges(3,3) * t254 + (rSges(3,1) * t256 + rSges(3,2) * t259) * t253;
t179 = -t236 * rSges(4,2) - t270;
t149 = t236 * t186;
t146 = t181 + t289;
t145 = -t180 + t279;
t136 = t139 * t305;
t134 = -t254 * t180 - t219 * t302;
t133 = t181 * t254 - t219 * t304;
t131 = t238 * t140;
t129 = rSges(6,1) * t236 + rSges(6,2) * t205 - t308;
t128 = rSges(5,3) * t236 - t274;
t108 = t271 + t178;
t107 = (rSges(4,2) - pkin(2)) * t236 + t270 + t276;
t105 = (t180 * t257 + t181 * t260) * t253;
t104 = t213 * t304 - t214 * t237 + t215 * t238;
t103 = -t213 * t302 - t235 * t214 + t236 * t215;
t102 = t235 * t216 - t236 * t217 - t218 * t302;
t101 = t216 * t237 - t217 * t238 + t218 * t304;
t97 = (-t179 - t187) * t254 + t260 * t278;
t96 = t178 * t254 + t257 * t278 + t185;
t87 = t127 * t305 - t165 * t238;
t86 = -t128 * t305 + t165 * t236;
t85 = t171 * t254 + (-t167 * t259 - t169 * t256) * t253;
t84 = t170 * t254 + (-t166 * t259 - t168 * t256) * t253;
t83 = t173 * t254 + (t175 * t259 + t177 * t256) * t253;
t82 = t172 * t254 + (t174 * t259 + t176 * t256) * t253;
t81 = t262 + t127;
t80 = (-rSges(5,3) - pkin(2)) * t236 + t268 + t274;
t77 = t79 * t254;
t76 = t78 * t254;
t75 = t79 * t305;
t74 = t78 * t305;
t73 = (t178 * t260 + t179 * t257) * t253 + t293;
t72 = -t127 * t236 + t128 * t238;
t67 = t261 + t130;
t66 = t308 + (-rSges(6,1) - pkin(2)) * t236 + (-rSges(6,2) + pkin(4)) * t205 + t266;
t65 = (-t128 + t291) * t254 + t260 * t275;
t64 = t127 * t254 + t257 * t275 + t292;
t63 = -t114 * t205 - t234 * t94;
t62 = -t114 * t202 + t234 * t95;
t61 = t130 * t305 + t136 + (-t164 - t186) * t238;
t60 = t164 * t236 + t149 + (-t129 - t140) * t305;
t55 = t116 * t201 - t120 * t202 + t124 * t238;
t54 = t115 * t201 - t119 * t202 + t123 * t238;
t53 = -t116 * t203 + t120 * t205 + t124 * t236;
t52 = -t115 * t203 + t119 * t205 + t123 * t236;
t51 = t118 * t236 + t122 * t203 - t126 * t205;
t50 = t117 * t236 + t121 * t203 - t125 * t205;
t49 = t118 * t238 - t122 * t201 + t126 * t202;
t48 = t117 * t238 - t121 * t201 + t125 * t202;
t47 = t261 + t306;
t46 = (-pkin(2) - pkin(5)) * t236 + (rSges(7,3) + pkin(4) + pkin(10)) * t205 + t266 + t273;
t45 = (t127 * t260 + t128 * t257) * t253 + t277;
t43 = t44 * t254;
t42 = t44 * t305;
t40 = t202 * t94 + t205 * t95;
t39 = t129 * t238 + t131 + (-t130 - t139) * t236;
t38 = (-t129 + t284) * t254 + t260 * t272;
t37 = t130 * t254 + t257 * t272 + t285;
t34 = (t129 * t257 + t130 * t260) * t253 + t269;
t33 = t136 + t306 * t305 + (-t186 - t296) * t238;
t32 = t149 + t296 * t236 + (-t140 - t307) * t305;
t29 = (t284 - t307) * t254 + t260 * t267;
t28 = t254 * t306 + t257 * t267 + t285;
t25 = t150 * t91 + t151 * t93 - t205 * t89;
t24 = t150 * t90 + t151 * t92 - t205 * t88;
t23 = t131 + t307 * t238 + (-t139 - t306) * t236;
t22 = (t257 * t307 + t260 * t306) * t253 + t269;
t21 = t77 + (t58 * t257 - t59 * t260) * t253;
t20 = t76 + (t57 * t257 - t56 * t260) * t253;
t19 = t59 * t236 + t58 * t238 + t75;
t18 = t56 * t236 + t57 * t238 + t74;
t17 = t71 * t254 + (t257 * t55 - t260 * t54) * t253;
t16 = t70 * t254 + (t257 * t53 - t260 * t52) * t253;
t15 = t69 * t254 + (t257 * t50 - t260 * t51) * t253;
t14 = t68 * t254 + (t257 * t48 - t260 * t49) * t253;
t13 = t236 * t54 + t238 * t55 + t305 * t71;
t12 = t236 * t52 + t238 * t53 + t305 * t70;
t11 = t236 * t51 + t238 * t50 + t305 * t69;
t10 = t236 * t49 + t238 * t48 + t305 * t68;
t9 = t43 + (t31 * t257 - t30 * t260) * t253;
t8 = t30 * t236 + t31 * t238 + t42;
t6 = t36 * t254 + (t257 * t27 - t26 * t260) * t253;
t5 = t35 * t254 + (-t24 * t260 + t25 * t257) * t253;
t4 = t236 * t26 + t238 * t27 + t305 * t36;
t3 = t236 * t24 + t238 * t25 + t305 * t35;
t1 = t202 * t25 - t205 * t24 + t234 * t35;
t98 = [Icges(2,3) + m(7) * (t46 ^ 2 + t47 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(6) * (t66 ^ 2 + t67 ^ 2) + m(4) * (t107 ^ 2 + t108 ^ 2) + m(3) * (t145 ^ 2 + t146 ^ 2) + m(2) * (t241 ^ 2 + t242 ^ 2) + t78 + t79 + t44 + t322; t43 + t77 + t76 + m(7) * (t28 * t47 + t29 * t46) + m(6) * (t37 * t67 + t38 * t66) + m(5) * (t64 * t81 + t65 * t80) + m(4) * (t107 * t97 + t108 * t96) + m(3) * (t133 * t146 + t134 * t145) + ((-t82 / 0.2e1 - t85 / 0.2e1 - t102 / 0.2e1 - t103 / 0.2e1 - t264) * t260 + (t83 / 0.2e1 + t84 / 0.2e1 + t101 / 0.2e1 + t104 / 0.2e1 + t263) * t257) * t253 + t297; (t9 + t20 + t21 + t297) * t254 + m(7) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t34 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t45 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t73 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(3) * (t105 ^ 2 + t133 ^ 2 + t134 ^ 2) + ((-t15 - t16 - t5 + ((t235 * t320 - t236 * t319) * t253 + t295 * t316 * t260) * t260 + (-t102 - t103 - t82 - t85) * t254) * t260 + (t6 + t17 + t14 + ((t237 * t318 + t238 * t317) * t253 + t294 * t316 * t257) * t257 + (t104 + t101 + t83 + t84) * t254 + ((t257 * t295 + t260 * t294) * t253 + t319 * t238 - t320 * t237 - t317 * t236 - t318 * t235) * t302) * t257) * t253; m(7) * (t235 * t47 + t237 * t46) + m(5) * (t235 * t81 + t237 * t80) + m(6) * (t235 * t67 + t237 * t66) + m(4) * (t107 * t237 + t108 * t235); m(7) * (-t22 * t303 + t235 * t28 + t237 * t29) + m(6) * (t235 * t37 + t237 * t38 - t303 * t34) + m(5) * (t235 * t64 + t237 * t65 - t303 * t45) + m(4) * (t235 * t96 + t237 * t97 - t303 * t73); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t288) * (t259 ^ 2 * t316 + t235 ^ 2 + t237 ^ 2); t42 + t75 + t74 + m(7) * (t32 * t46 + t33 * t47) + m(5) * (t80 * t86 + t81 * t87) + m(6) * (t60 * t66 + t61 * t67) + t263 * t238 + t264 * t236; (t8 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1) * t254 + (t6 / 0.2e1 + t17 / 0.2e1 + t14 / 0.2e1) * t238 + (t5 / 0.2e1 + t16 / 0.2e1 + t15 / 0.2e1) * t236 + m(7) * (t22 * t23 + t28 * t33 + t29 * t32) + m(6) * (t34 * t39 + t37 * t61 + t38 * t60) + m(5) * (t45 * t72 + t64 * t87 + t65 * t86) + ((-t3 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1) * t260 + (t4 / 0.2e1 + t10 / 0.2e1 + t13 / 0.2e1) * t257 + (t9 / 0.2e1 + t21 / 0.2e1 + t20 / 0.2e1) * t256) * t253; m(5) * (t235 * t87 + t237 * t86 - t303 * t72) + m(6) * (t235 * t61 + t237 * t60 - t303 * t39) + m(7) * (-t23 * t303 + t235 * t33 + t237 * t32); (t18 + t19 + t8) * t305 + (t4 + t10 + t13) * t238 + (t3 + t11 + t12) * t236 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t39 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t72 ^ 2 + t86 ^ 2 + t87 ^ 2); m(7) * (t201 * t46 - t203 * t47) + m(6) * (t201 * t66 - t203 * t67); m(7) * (t201 * t29 - t203 * t28 + t22 * t233) + m(6) * (t201 * t38 - t203 * t37 + t233 * t34); (t201 * t237 - t203 * t235 - t233 * t303) * t321; m(7) * (t201 * t32 - t203 * t33 + t23 * t233) + m(6) * (t201 * t60 - t203 * t61 + t233 * t39); (t201 ^ 2 + t203 ^ 2 + t233 ^ 2) * t321; m(7) * (t46 * t63 + t47 * t62) + t41 + t287 * t205 + t286 * t202; m(7) * (t22 * t40 + t28 * t62 + t29 * t63) + t9 * t311 + t5 * t312 + t6 * t313 + t254 * t314 + (t257 * t315 - t260 * t1 / 0.2e1) * t253; m(7) * (t235 * t62 + t237 * t63 - t303 * t40); t305 * t314 + m(7) * (t23 * t40 + t32 * t63 + t33 * t62) + t4 * t313 + t3 * t312 + t238 * t315 + t236 * t1 / 0.2e1 + t8 * t311; m(7) * (t201 * t63 - t203 * t62 + t233 * t40); t202 * t2 - t205 * t1 + t234 * t7 + m(7) * (t40 ^ 2 + t62 ^ 2 + t63 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t98(1) t98(2) t98(4) t98(7) t98(11) t98(16); t98(2) t98(3) t98(5) t98(8) t98(12) t98(17); t98(4) t98(5) t98(6) t98(9) t98(13) t98(18); t98(7) t98(8) t98(9) t98(10) t98(14) t98(19); t98(11) t98(12) t98(13) t98(14) t98(15) t98(20); t98(16) t98(17) t98(18) t98(19) t98(20) t98(21);];
Mq  = res;
