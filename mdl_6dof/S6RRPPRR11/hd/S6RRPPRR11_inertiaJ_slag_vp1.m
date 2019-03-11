% Calculate joint inertia matrix for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:15
% EndTime: 2019-03-09 09:39:26
% DurationCPUTime: 5.54s
% Computational Cost: add. (14766->570), mult. (25786->799), div. (0->0), fcn. (32009->12), ass. (0->252)
t243 = sin(pkin(11));
t245 = cos(pkin(11));
t246 = cos(pkin(6));
t244 = sin(pkin(6));
t252 = cos(qJ(2));
t296 = t244 * t252;
t219 = -t243 * t246 - t245 * t296;
t220 = -t243 * t296 + t245 * t246;
t249 = sin(qJ(2));
t298 = t244 * t249;
t143 = Icges(5,5) * t220 + Icges(5,6) * t219 + Icges(5,3) * t298;
t144 = Icges(5,4) * t220 + Icges(5,2) * t219 + Icges(5,6) * t298;
t145 = Icges(5,1) * t220 + Icges(5,4) * t219 + Icges(5,5) * t298;
t198 = Icges(3,3) * t246 + (Icges(3,5) * t249 + Icges(3,6) * t252) * t244;
t199 = Icges(3,6) * t246 + (Icges(3,4) * t249 + Icges(3,2) * t252) * t244;
t200 = Icges(3,5) * t246 + (Icges(3,1) * t249 + Icges(3,4) * t252) * t244;
t201 = Icges(4,5) * t246 + (-Icges(4,6) * t249 - Icges(4,3) * t252) * t244;
t202 = Icges(4,4) * t246 + (-Icges(4,2) * t249 - Icges(4,6) * t252) * t244;
t203 = Icges(4,1) * t246 + (-Icges(4,4) * t249 - Icges(4,5) * t252) * t244;
t319 = (-t252 * t201 - t249 * t202) * t244 + t199 * t296 + t219 * t144 + t220 * t145 + (t200 + t143) * t298 + (t203 + t198) * t246;
t253 = cos(qJ(1));
t291 = t252 * t253;
t250 = sin(qJ(1));
t294 = t249 * t250;
t221 = -t246 * t291 + t294;
t292 = t250 * t252;
t293 = t249 * t253;
t222 = t246 * t293 + t292;
t295 = t244 * t253;
t149 = -Icges(4,5) * t295 - Icges(4,6) * t222 + Icges(4,3) * t221;
t156 = Icges(3,4) * t222 - Icges(3,2) * t221 - Icges(3,6) * t295;
t318 = t149 - t156;
t223 = t246 * t292 + t293;
t224 = -t246 * t294 + t291;
t297 = t244 * t250;
t148 = Icges(4,5) * t297 - Icges(4,6) * t224 + Icges(4,3) * t223;
t157 = Icges(3,4) * t224 - Icges(3,2) * t223 + Icges(3,6) * t297;
t317 = -t157 + t148;
t151 = -Icges(4,4) * t295 - Icges(4,2) * t222 + Icges(4,6) * t221;
t158 = Icges(3,1) * t222 - Icges(3,4) * t221 - Icges(3,5) * t295;
t316 = -t158 + t151;
t150 = Icges(4,4) * t297 - Icges(4,2) * t224 + Icges(4,6) * t223;
t159 = Icges(3,1) * t224 - Icges(3,4) * t223 + Icges(3,5) * t297;
t315 = -t159 + t150;
t272 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t314 = 0.2e1 * t272;
t153 = -Icges(4,1) * t295 - Icges(4,4) * t222 + Icges(4,5) * t221;
t154 = Icges(3,5) * t222 - Icges(3,6) * t221 - Icges(3,3) * t295;
t313 = -t154 - t153;
t152 = Icges(4,1) * t297 - Icges(4,4) * t224 + Icges(4,5) * t223;
t155 = Icges(3,5) * t224 - Icges(3,6) * t223 + Icges(3,3) * t297;
t312 = t155 + t152;
t311 = t319 * t246;
t284 = pkin(11) + qJ(5);
t240 = sin(t284);
t274 = cos(t284);
t178 = -t223 * t274 + t240 * t297;
t180 = t221 * t274 + t240 * t295;
t264 = t244 * t274;
t179 = t223 * t240 + t250 * t264;
t248 = sin(qJ(6));
t251 = cos(qJ(6));
t132 = -t179 * t248 + t224 * t251;
t133 = t179 * t251 + t224 * t248;
t67 = Icges(7,5) * t133 + Icges(7,6) * t132 + Icges(7,3) * t178;
t69 = Icges(7,4) * t133 + Icges(7,2) * t132 + Icges(7,6) * t178;
t71 = Icges(7,1) * t133 + Icges(7,4) * t132 + Icges(7,5) * t178;
t19 = t132 * t69 + t133 * t71 + t178 * t67;
t181 = t221 * t240 - t253 * t264;
t134 = -t181 * t248 + t222 * t251;
t135 = t181 * t251 + t222 * t248;
t68 = Icges(7,5) * t135 + Icges(7,6) * t134 - Icges(7,3) * t180;
t70 = Icges(7,4) * t135 + Icges(7,2) * t134 - Icges(7,6) * t180;
t72 = Icges(7,1) * t135 + Icges(7,4) * t134 - Icges(7,5) * t180;
t20 = t132 * t70 + t133 * t72 + t178 * t68;
t204 = t246 * t240 + t252 * t264;
t205 = -t240 * t296 + t246 * t274;
t176 = -t205 * t248 + t251 * t298;
t177 = t205 * t251 + t248 * t298;
t91 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t204;
t92 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t204;
t93 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t204;
t27 = t132 * t92 + t133 * t93 + t178 * t91;
t1 = t178 * t19 - t180 * t20 + t204 * t27;
t310 = t1 / 0.2e1;
t23 = t176 * t69 + t177 * t71 + t204 * t67;
t24 = t176 * t70 + t177 * t72 + t204 * t68;
t36 = t176 * t92 + t177 * t93 + t204 * t91;
t31 = t36 * t204;
t7 = t23 * t178 - t24 * t180 + t31;
t309 = t7 / 0.2e1;
t308 = t178 / 0.2e1;
t307 = -t180 / 0.2e1;
t306 = t204 / 0.2e1;
t305 = pkin(5) * t181;
t247 = -pkin(9) - qJ(4);
t304 = -pkin(2) + t247;
t73 = t133 * rSges(7,1) + t132 * rSges(7,2) + t178 * rSges(7,3);
t303 = t179 * pkin(5) + pkin(10) * t178 + t73;
t265 = -rSges(7,1) * t135 - rSges(7,2) * t134;
t74 = -rSges(7,3) * t180 - t265;
t302 = -pkin(10) * t180 + t305 + t74;
t94 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t204;
t301 = pkin(5) * t205 + pkin(10) * t204 + t94;
t300 = t221 * t243;
t299 = t223 * t243;
t210 = t221 * qJ(3);
t170 = t222 * pkin(2) + t210;
t171 = t224 * pkin(2) + qJ(3) * t223;
t290 = t170 * t297 + t171 * t295;
t169 = t246 * t171;
t193 = pkin(3) * t297 + qJ(4) * t224;
t289 = t246 * t193 + t169;
t285 = pkin(3) * t295 - t222 * qJ(4);
t288 = -t170 + t285;
t225 = (pkin(2) * t249 - qJ(3) * t252) * t244;
t287 = -pkin(3) * t246 - qJ(4) * t298 - t225;
t286 = t253 * pkin(1) + pkin(8) * t297;
t283 = t27 / 0.2e1 + t23 / 0.2e1;
t28 = t134 * t92 + t135 * t93 - t180 * t91;
t282 = -t28 / 0.2e1 - t24 / 0.2e1;
t239 = pkin(4) * t245 + pkin(3);
t276 = pkin(4) * t299 - t224 * t247 + t239 * t297;
t118 = -t193 + t276;
t281 = t246 * t118 + t289;
t270 = -pkin(4) * t300 + t239 * t295;
t119 = -t222 * t247 - t270 + t285;
t280 = -t119 + t288;
t138 = Icges(6,5) * t205 - Icges(6,6) * t204 + Icges(6,3) * t298;
t139 = Icges(6,4) * t205 - Icges(6,2) * t204 + Icges(6,6) * t298;
t140 = Icges(6,1) * t205 - Icges(6,4) * t204 + Icges(6,5) * t298;
t60 = t138 * t298 - t204 * t139 + t205 * t140;
t278 = -(-pkin(3) + t239) * t246 - (-pkin(4) * t243 * t252 + (-qJ(4) - t247) * t249) * t244 + t287;
t103 = t179 * rSges(6,1) - t178 * rSges(6,2) + t224 * rSges(6,3);
t186 = t223 * t245 - t243 * t297;
t187 = t245 * t297 + t299;
t116 = t187 * rSges(5,1) + t186 * rSges(5,2) + t224 * rSges(5,3);
t163 = t224 * rSges(3,1) - t223 * rSges(3,2) + rSges(3,3) * t297;
t160 = rSges(4,1) * t297 - t224 * rSges(4,2) + t223 * rSges(4,3);
t275 = -t250 * pkin(1) + pkin(8) * t295;
t273 = t244 * (-rSges(4,1) * t246 - (-rSges(4,2) * t249 - rSges(4,3) * t252) * t244 - t225);
t271 = t193 * t295 - t285 * t297 + t290;
t269 = -t210 + t275;
t268 = t244 * (-rSges(5,1) * t220 - rSges(5,2) * t219 - rSges(5,3) * t298 + t287);
t188 = t221 * t245 + t243 * t295;
t189 = -t245 * t295 + t300;
t267 = -rSges(5,1) * t189 - rSges(5,2) * t188;
t266 = -rSges(6,1) * t181 - rSges(6,2) * t180;
t141 = rSges(6,1) * t205 - rSges(6,2) * t204 + rSges(6,3) * t298;
t263 = t244 * (-t141 + t278);
t262 = t171 + t286;
t101 = Icges(6,1) * t179 - Icges(6,4) * t178 + Icges(6,5) * t224;
t97 = Icges(6,5) * t179 - Icges(6,6) * t178 + Icges(6,3) * t224;
t99 = Icges(6,4) * t179 - Icges(6,2) * t178 + Icges(6,6) * t224;
t45 = t101 * t205 - t204 * t99 + t298 * t97;
t53 = t138 * t224 - t139 * t178 + t140 * t179;
t261 = t53 / 0.2e1 + t45 / 0.2e1 + t283;
t100 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t222;
t102 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t222;
t98 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t222;
t46 = -t100 * t204 + t102 * t205 + t298 * t98;
t54 = t138 * t222 + t139 * t180 + t140 * t181;
t260 = t54 / 0.2e1 + t46 / 0.2e1 - t282;
t259 = rSges(4,1) * t295 - t221 * rSges(4,3);
t258 = t118 * t295 + t119 * t297 + t271;
t257 = t244 * (t278 - t301);
t162 = t222 * rSges(3,1) - t221 * rSges(3,2) - rSges(3,3) * t295;
t255 = t269 + t270;
t254 = t262 + t276;
t242 = t244 ^ 2;
t230 = rSges(2,1) * t253 - t250 * rSges(2,2);
t229 = -t250 * rSges(2,1) - rSges(2,2) * t253;
t206 = rSges(3,3) * t246 + (rSges(3,1) * t249 + rSges(3,2) * t252) * t244;
t161 = -t222 * rSges(4,2) - t259;
t137 = t163 + t286;
t136 = -t162 + t275;
t123 = -t246 * t162 - t206 * t295;
t122 = t163 * t246 - t206 * t297;
t117 = rSges(5,3) * t222 - t267;
t115 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t222;
t114 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t224;
t113 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t222;
t112 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t224;
t111 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t222;
t110 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t224;
t104 = rSges(6,3) * t222 - t266;
t96 = t262 + t160;
t95 = (rSges(4,2) - pkin(2)) * t222 + t259 + t269;
t90 = (t162 * t250 + t163 * t253) * t244;
t89 = t198 * t297 - t199 * t223 + t200 * t224;
t88 = -t198 * t295 - t221 * t199 + t222 * t200;
t87 = t221 * t201 - t222 * t202 - t203 * t295;
t86 = t201 * t223 - t202 * t224 + t203 * t297;
t82 = (-t161 - t170) * t246 + t253 * t273;
t81 = t160 * t246 + t250 * t273 + t169;
t80 = t153 * t246 + (-t149 * t252 - t151 * t249) * t244;
t79 = t152 * t246 + (-t148 * t252 - t150 * t249) * t244;
t78 = t155 * t246 + (t157 * t252 + t159 * t249) * t244;
t77 = t154 * t246 + (t156 * t252 + t158 * t249) * t244;
t76 = t193 + t262 + t116;
t75 = (-rSges(5,3) - pkin(2)) * t222 + t267 + t269 + t285;
t66 = t103 * t298 - t141 * t224;
t65 = -t104 * t298 + t141 * t222;
t63 = (t160 * t253 + t161 * t250) * t244 + t290;
t62 = t254 + t103;
t61 = (-rSges(6,3) + t304) * t222 + t255 + t266;
t59 = -t103 * t222 + t104 * t224;
t58 = t60 * t246;
t57 = t143 * t222 + t144 * t188 + t145 * t189;
t56 = t143 * t224 + t144 * t186 + t145 * t187;
t55 = t60 * t298;
t52 = (-t117 + t288) * t246 + t253 * t268;
t51 = t116 * t246 + t250 * t268 + t289;
t50 = t111 * t298 + t113 * t219 + t115 * t220;
t49 = t110 * t298 + t112 * t219 + t114 * t220;
t48 = -t180 * t94 - t204 * t74;
t47 = -t178 * t94 + t204 * t73;
t44 = t254 + t303;
t43 = -t305 + t304 * t222 + (rSges(7,3) + pkin(10)) * t180 + t255 + t265;
t42 = (t116 * t253 + t117 * t250) * t244 + t271;
t41 = t100 * t180 + t102 * t181 + t222 * t98;
t40 = t101 * t181 + t180 * t99 + t222 * t97;
t39 = -t100 * t178 + t102 * t179 + t224 * t98;
t38 = t101 * t179 - t178 * t99 + t224 * t97;
t37 = t178 * t74 + t180 * t73;
t35 = t36 * t246;
t34 = t36 * t298;
t33 = (-t104 + t280) * t246 + t253 * t263;
t32 = t103 * t246 + t250 * t263 + t281;
t30 = -t224 * t301 + t298 * t303;
t29 = t222 * t301 - t298 * t302;
t26 = (t103 * t253 + t104 * t250) * t244 + t258;
t25 = -t222 * t303 + t224 * t302;
t22 = t134 * t70 + t135 * t72 - t180 * t68;
t21 = t134 * t69 + t135 * t71 - t180 * t67;
t18 = (t280 - t302) * t246 + t253 * t257;
t17 = t246 * t303 + t250 * t257 + t281;
t16 = (t250 * t302 + t253 * t303) * t244 + t258;
t15 = t58 + (t45 * t250 - t46 * t253) * t244;
t14 = t46 * t222 + t45 * t224 + t55;
t13 = t54 * t246 + (t250 * t40 - t253 * t41) * t244;
t12 = t53 * t246 + (t250 * t38 - t253 * t39) * t244;
t11 = t222 * t41 + t224 * t40 + t298 * t54;
t10 = t222 * t39 + t224 * t38 + t298 * t53;
t9 = t35 + (t23 * t250 - t24 * t253) * t244;
t8 = t24 * t222 + t23 * t224 + t34;
t6 = t28 * t246 + (t21 * t250 - t22 * t253) * t244;
t5 = t27 * t246 + (t19 * t250 - t20 * t253) * t244;
t4 = t21 * t224 + t22 * t222 + t28 * t298;
t3 = t19 * t224 + t20 * t222 + t27 * t298;
t2 = t178 * t21 - t180 * t22 + t204 * t28;
t64 = [Icges(2,3) + m(7) * (t43 ^ 2 + t44 ^ 2) + m(6) * (t61 ^ 2 + t62 ^ 2) + m(5) * (t75 ^ 2 + t76 ^ 2) + m(4) * (t95 ^ 2 + t96 ^ 2) + m(3) * (t136 ^ 2 + t137 ^ 2) + m(2) * (t229 ^ 2 + t230 ^ 2) + t60 + t36 + t319; t58 + t35 + m(7) * (t17 * t44 + t18 * t43) + m(6) * (t32 * t62 + t33 * t61) + m(5) * (t51 * t76 + t52 * t75) + m(4) * (t81 * t96 + t82 * t95) + m(3) * (t122 * t137 + t123 * t136) + ((-t50 / 0.2e1 - t77 / 0.2e1 - t80 / 0.2e1 - t57 / 0.2e1 - t87 / 0.2e1 - t88 / 0.2e1 - t260) * t253 + (t56 / 0.2e1 + t86 / 0.2e1 + t89 / 0.2e1 + t49 / 0.2e1 + t78 / 0.2e1 + t79 / 0.2e1 + t261) * t250) * t244 + t311; m(7) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(6) * (t26 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t42 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(4) * (t63 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(3) * (t122 ^ 2 + t123 ^ 2 + t90 ^ 2) + (t9 + t15 + ((-t50 - t77 - t80) * t253 + (t49 + t78 + t79) * t250) * t244 + t311) * t246 + (t5 + t12 + ((t110 * t224 + t112 * t186 + t114 * t187) * t250 - (t111 * t224 + t113 * t186 + t115 * t187) * t253) * t244 + (t223 * t317 - t224 * t315 + t312 * t297) * t297 + (t56 + t89 + t86) * t246) * t297 + (-t6 - t13 - ((t110 * t222 + t112 * t188 + t114 * t189) * t250 - (t111 * t222 + t113 * t188 + t115 * t189) * t253) * t244 + (t221 * t318 - t316 * t222 + t313 * t295) * t295 + (-t57 - t88 - t87) * t246 + (-t317 * t221 + t315 * t222 - t223 * t318 + t316 * t224 + t312 * t295 + t313 * t297) * t297) * t295; m(7) * (t221 * t44 + t223 * t43) + m(6) * (t221 * t62 + t223 * t61) + m(5) * (t221 * t76 + t223 * t75) + m(4) * (t221 * t96 + t223 * t95); m(7) * (-t16 * t296 + t17 * t221 + t18 * t223) + m(6) * (t221 * t32 + t223 * t33 - t26 * t296) + m(5) * (t221 * t51 + t223 * t52 - t296 * t42) + m(4) * (t221 * t81 + t223 * t82 - t296 * t63); 0.2e1 * (m(4) / 0.2e1 + t272) * (t242 * t252 ^ 2 + t221 ^ 2 + t223 ^ 2); m(7) * (t222 * t44 + t224 * t43) + m(6) * (t222 * t62 + t224 * t61) + m(5) * (t222 * t76 + t224 * t75); m(7) * (t16 * t298 + t17 * t222 + t18 * t224) + m(6) * (t222 * t32 + t224 * t33 + t26 * t298) + m(5) * (t222 * t51 + t224 * t52 + t298 * t42); (-t242 * t249 * t252 + t221 * t222 + t223 * t224) * t314; (t242 * t249 ^ 2 + t222 ^ 2 + t224 ^ 2) * t314; t34 + t55 + m(7) * (t29 * t43 + t30 * t44) + m(6) * (t61 * t65 + t62 * t66) + t261 * t224 + t260 * t222; (t8 / 0.2e1 + t14 / 0.2e1) * t246 + (t5 / 0.2e1 + t12 / 0.2e1) * t224 + (t6 / 0.2e1 + t13 / 0.2e1) * t222 + m(7) * (t16 * t25 + t17 * t30 + t18 * t29) + m(6) * (t26 * t59 + t32 * t66 + t33 * t65) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t253 + (t3 / 0.2e1 + t10 / 0.2e1) * t250 + (t9 / 0.2e1 + t15 / 0.2e1) * t249) * t244; m(6) * (t221 * t66 + t223 * t65 - t296 * t59) + m(7) * (t221 * t30 + t223 * t29 - t25 * t296); m(6) * (t222 * t66 + t224 * t65 + t298 * t59) + m(7) * (t222 * t30 + t224 * t29 + t25 * t298); (t14 + t8) * t298 + (t3 + t10) * t224 + (t4 + t11) * t222 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t59 ^ 2 + t65 ^ 2 + t66 ^ 2); t31 + m(7) * (t43 * t48 + t44 * t47) + t282 * t180 + t283 * t178; m(7) * (t16 * t37 + t17 * t47 + t18 * t48) + t6 * t307 + t9 * t306 + t246 * t309 + t5 * t308 + (t250 * t310 - t253 * t2 / 0.2e1) * t244; m(7) * (t221 * t47 + t223 * t48 - t296 * t37); m(7) * (t222 * t47 + t224 * t48 + t298 * t37); t8 * t306 + m(7) * (t25 * t37 + t29 * t48 + t30 * t47) + t224 * t310 + t222 * t2 / 0.2e1 + t4 * t307 + t3 * t308 + t298 * t309; t178 * t1 - t180 * t2 + t204 * t7 + m(7) * (t37 ^ 2 + t47 ^ 2 + t48 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t64(1) t64(2) t64(4) t64(7) t64(11) t64(16); t64(2) t64(3) t64(5) t64(8) t64(12) t64(17); t64(4) t64(5) t64(6) t64(9) t64(13) t64(18); t64(7) t64(8) t64(9) t64(10) t64(14) t64(19); t64(11) t64(12) t64(13) t64(14) t64(15) t64(20); t64(16) t64(17) t64(18) t64(19) t64(20) t64(21);];
Mq  = res;
