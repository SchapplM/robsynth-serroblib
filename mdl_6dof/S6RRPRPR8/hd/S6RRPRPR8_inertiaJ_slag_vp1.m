% Calculate joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:35
% EndTime: 2019-03-09 10:49:44
% DurationCPUTime: 3.84s
% Computational Cost: add. (10658->517), mult. (14475->750), div. (0->0), fcn. (16846->10), ass. (0->235)
t216 = sin(qJ(1));
t290 = -t216 / 0.2e1;
t289 = t216 / 0.2e1;
t219 = cos(qJ(1));
t286 = t219 / 0.2e1;
t207 = pkin(10) + qJ(4);
t201 = sin(t207);
t202 = cos(t207);
t215 = sin(qJ(2));
t218 = cos(qJ(2));
t140 = -Icges(5,3) * t218 + (Icges(5,5) * t202 - Icges(5,6) * t201) * t215;
t141 = -Icges(6,2) * t218 + (Icges(6,4) * t202 + Icges(6,6) * t201) * t215;
t303 = -t140 - t141;
t265 = t216 * t218;
t160 = t201 * t265 + t202 * t219;
t161 = -t201 * t219 + t202 * t265;
t268 = t215 * t216;
t90 = Icges(6,5) * t161 + Icges(6,6) * t268 + Icges(6,3) * t160;
t94 = Icges(6,4) * t161 + Icges(6,2) * t268 + Icges(6,6) * t160;
t98 = Icges(6,1) * t161 + Icges(6,4) * t268 + Icges(6,5) * t160;
t36 = t160 * t90 + t161 * t98 + t268 * t94;
t264 = t218 * t219;
t162 = t201 * t264 - t216 * t202;
t163 = t216 * t201 + t202 * t264;
t267 = t215 * t219;
t91 = Icges(6,5) * t163 + Icges(6,6) * t267 + Icges(6,3) * t162;
t95 = Icges(6,4) * t163 + Icges(6,2) * t267 + Icges(6,6) * t162;
t99 = Icges(6,1) * t163 + Icges(6,4) * t267 + Icges(6,5) * t162;
t37 = t160 * t91 + t161 * t99 + t268 * t95;
t100 = Icges(5,1) * t161 - Icges(5,4) * t160 + Icges(5,5) * t268;
t92 = Icges(5,5) * t161 - Icges(5,6) * t160 + Icges(5,3) * t268;
t96 = Icges(5,4) * t161 - Icges(5,2) * t160 + Icges(5,6) * t268;
t38 = t100 * t161 - t160 * t96 + t268 * t92;
t101 = Icges(5,1) * t163 - Icges(5,4) * t162 + Icges(5,5) * t267;
t93 = Icges(5,5) * t163 - Icges(5,6) * t162 + Icges(5,3) * t267;
t97 = Icges(5,4) * t163 - Icges(5,2) * t162 + Icges(5,6) * t267;
t39 = t101 * t161 - t160 * t97 + t268 * t93;
t139 = -Icges(6,6) * t218 + (Icges(6,5) * t202 + Icges(6,3) * t201) * t215;
t143 = -Icges(6,4) * t218 + (Icges(6,1) * t202 + Icges(6,5) * t201) * t215;
t57 = t139 * t160 + t141 * t268 + t143 * t161;
t142 = -Icges(5,6) * t218 + (Icges(5,4) * t202 - Icges(5,2) * t201) * t215;
t144 = -Icges(5,5) * t218 + (Icges(5,1) * t202 - Icges(5,4) * t201) * t215;
t58 = t140 * t268 - t142 * t160 + t144 * t161;
t302 = (-t58 - t57) * t218 + ((t37 + t39) * t219 + (t36 + t38) * t216) * t215;
t40 = t162 * t90 + t163 * t98 + t267 * t94;
t41 = t162 * t91 + t163 * t99 + t267 * t95;
t42 = t163 * t100 - t162 * t96 + t267 * t92;
t43 = t163 * t101 - t162 * t97 + t267 * t93;
t59 = t162 * t139 + t141 * t267 + t163 * t143;
t60 = t140 * t267 - t162 * t142 + t163 * t144;
t301 = (-t59 - t60) * t218 + ((t41 + t43) * t219 + (t40 + t42) * t216) * t215;
t44 = -t218 * t94 + (t201 * t90 + t202 * t98) * t215;
t46 = -t218 * t92 + (t100 * t202 - t201 * t96) * t215;
t300 = -t44 - t46;
t45 = -t218 * t95 + (t201 * t91 + t202 * t99) * t215;
t47 = -t218 * t93 + (t101 * t202 - t201 * t97) * t215;
t299 = t45 + t47;
t271 = t201 * t215;
t270 = t202 * t215;
t297 = t139 * t271 + (t143 + t144) * t270;
t298 = (-t142 * t271 + t218 * t303 + t297) * t218;
t214 = sin(qJ(6));
t217 = cos(qJ(6));
t113 = t160 * t217 - t161 * t214;
t114 = t160 * t214 + t161 * t217;
t62 = Icges(7,5) * t114 + Icges(7,6) * t113 - Icges(7,3) * t268;
t64 = Icges(7,4) * t114 + Icges(7,2) * t113 - Icges(7,6) * t268;
t66 = Icges(7,1) * t114 + Icges(7,4) * t113 - Icges(7,5) * t268;
t12 = t113 * t64 + t114 * t66 - t268 * t62;
t115 = t162 * t217 - t163 * t214;
t116 = t162 * t214 + t163 * t217;
t63 = Icges(7,5) * t116 + Icges(7,6) * t115 - Icges(7,3) * t267;
t65 = Icges(7,4) * t116 + Icges(7,2) * t115 - Icges(7,6) * t267;
t67 = Icges(7,1) * t116 + Icges(7,4) * t115 - Icges(7,5) * t267;
t13 = t113 * t65 + t114 * t67 - t268 * t63;
t147 = (t201 * t217 - t202 * t214) * t215;
t148 = (t201 * t214 + t202 * t217) * t215;
t86 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t218;
t87 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t218;
t88 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t218;
t28 = t113 * t87 + t114 * t88 - t268 * t86;
t1 = -t28 * t218 + (t12 * t216 + t13 * t219) * t215;
t14 = t115 * t64 + t116 * t66 - t267 * t62;
t15 = t115 * t65 + t116 * t67 - t267 * t63;
t29 = t115 * t87 + t116 * t88 - t267 * t86;
t2 = -t29 * t218 + (t14 * t216 + t15 * t219) * t215;
t21 = t147 * t64 + t148 * t66 + t218 * t62;
t22 = t147 * t65 + t148 * t67 + t218 * t63;
t250 = t147 * t87 + t148 * t88 + t218 * t86;
t32 = t250 * t218;
t3 = -t32 + (t21 * t216 + t22 * t219) * t215;
t296 = (t216 * t1 + t219 * t2) * t215 - t218 * t3;
t209 = t216 ^ 2;
t210 = t219 ^ 2;
t295 = 0.2e1 * t215;
t294 = m(4) / 0.2e1;
t293 = m(5) / 0.2e1;
t292 = m(6) / 0.2e1;
t291 = m(7) / 0.2e1;
t287 = -t219 / 0.2e1;
t285 = rSges(7,3) + pkin(9);
t187 = rSges(3,1) * t215 + rSges(3,2) * t218;
t284 = m(3) * t187;
t283 = pkin(2) * t218;
t212 = cos(pkin(10));
t200 = pkin(3) * t212 + pkin(2);
t281 = -pkin(2) + t200;
t279 = t160 * rSges(6,3);
t278 = t219 * rSges(3,3);
t232 = -t114 * rSges(7,1) - t113 * rSges(7,2);
t68 = -rSges(7,3) * t268 - t232;
t277 = t161 * pkin(5) - pkin(9) * t268 + t68;
t154 = t163 * pkin(5);
t262 = t116 * rSges(7,1) + t115 * rSges(7,2);
t69 = -rSges(7,3) * t267 + t262;
t276 = -pkin(9) * t267 + t154 + t69;
t89 = rSges(7,1) * t148 + rSges(7,2) * t147 + rSges(7,3) * t218;
t275 = pkin(5) * t270 + pkin(9) * t218 + t89;
t274 = Icges(3,4) * t215;
t273 = Icges(3,4) * t218;
t272 = qJ(3) * t215;
t211 = sin(pkin(10));
t269 = t211 * t219;
t266 = t216 * t211;
t104 = t163 * rSges(6,1) + rSges(6,2) * t267 + t162 * rSges(6,3);
t121 = t163 * pkin(4) + t162 * qJ(5);
t263 = -t104 - t121;
t149 = t160 * qJ(5);
t120 = t161 * pkin(4) + t149;
t170 = (pkin(4) * t202 + qJ(5) * t201) * t215;
t261 = t218 * t120 + t170 * t268;
t186 = t215 * pkin(2) - t218 * qJ(3);
t213 = -pkin(8) - qJ(3);
t259 = -(qJ(3) + t213) * t218 - t281 * t215 - t186;
t258 = t218 * rSges(4,3) - (rSges(4,1) * t212 - rSges(4,2) * t211) * t215 - t186;
t254 = pkin(2) * t264 + qJ(3) * t267;
t257 = t209 * (t272 + t283) + t219 * t254;
t256 = pkin(3) * t266 + t200 * t264;
t255 = -pkin(3) * t269 - t213 * t268;
t253 = t219 * pkin(1) + t216 * pkin(7);
t252 = t209 + t210;
t251 = t292 + t291;
t249 = -t22 / 0.2e1 - t29 / 0.2e1;
t248 = -t28 / 0.2e1 - t21 / 0.2e1;
t247 = -t121 - t276;
t246 = t213 * t267;
t146 = -t218 * rSges(5,3) + (rSges(5,1) * t202 - rSges(5,2) * t201) * t215;
t245 = -t146 + t259;
t244 = -t170 + t259;
t105 = t163 * rSges(5,1) - t162 * rSges(5,2) + rSges(5,3) * t267;
t178 = -t211 * t264 + t216 * t212;
t179 = t212 * t264 + t266;
t243 = t179 * rSges(4,1) + t178 * rSges(4,2) + rSges(4,3) * t267;
t205 = t219 * pkin(7);
t242 = t205 - t255;
t241 = -t200 * t218 - pkin(1);
t225 = -t246 + t256;
t240 = t216 * ((t218 * t281 - t272) * t216 + t255) + t219 * (t225 - t254) + t257;
t145 = -t218 * rSges(6,2) + (rSges(6,1) * t202 + rSges(6,3) * t201) * t215;
t239 = -t145 + t244;
t238 = -t149 + t242;
t236 = t244 - t275;
t235 = rSges(3,1) * t218 - rSges(3,2) * t215;
t176 = -t211 * t265 - t212 * t219;
t177 = t212 * t265 - t269;
t234 = -t177 * rSges(4,1) - t176 * rSges(4,2);
t233 = -t161 * rSges(5,1) + t160 * rSges(5,2);
t231 = Icges(3,1) * t218 - t274;
t230 = -Icges(3,2) * t215 + t273;
t229 = Icges(3,5) * t218 - Icges(3,6) * t215;
t226 = rSges(3,1) * t264 - rSges(3,2) * t267 + t216 * rSges(3,3);
t224 = t216 * t120 + t219 * t121 + t240;
t223 = t121 + t253 + t256;
t222 = t2 * t290 + t218 * (-t21 * t219 + t22 * t216) / 0.2e1 + t1 * t286;
t221 = t47 / 0.2e1 + t45 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1 - t249;
t220 = t44 / 0.2e1 + t58 / 0.2e1 + t57 / 0.2e1 + t46 / 0.2e1 - t248;
t208 = t215 ^ 2;
t189 = rSges(2,1) * t219 - t216 * rSges(2,2);
t188 = -t216 * rSges(2,1) - rSges(2,2) * t219;
t183 = Icges(3,5) * t215 + Icges(3,6) * t218;
t165 = Icges(3,3) * t216 + t219 * t229;
t164 = -Icges(3,3) * t219 + t216 * t229;
t158 = -Icges(4,5) * t218 + (Icges(4,1) * t212 - Icges(4,4) * t211) * t215;
t157 = -Icges(4,6) * t218 + (Icges(4,4) * t212 - Icges(4,2) * t211) * t215;
t136 = t226 + t253;
t135 = t278 + t205 + (-pkin(1) - t235) * t216;
t129 = t258 * t219;
t128 = t258 * t216;
t127 = Icges(4,1) * t179 + Icges(4,4) * t178 + Icges(4,5) * t267;
t126 = Icges(4,1) * t177 + Icges(4,4) * t176 + Icges(4,5) * t268;
t125 = Icges(4,4) * t179 + Icges(4,2) * t178 + Icges(4,6) * t267;
t124 = Icges(4,4) * t177 + Icges(4,2) * t176 + Icges(4,6) * t268;
t123 = Icges(4,5) * t179 + Icges(4,6) * t178 + Icges(4,3) * t267;
t122 = Icges(4,5) * t177 + Icges(4,6) * t176 + Icges(4,3) * t268;
t119 = t219 * t226 + (t216 * t235 - t278) * t216;
t107 = t120 * t267;
t103 = rSges(5,3) * t268 - t233;
t102 = t161 * rSges(6,1) + rSges(6,2) * t268 + t279;
t83 = t243 + t253 + t254;
t82 = t205 + (-t283 - pkin(1) + (-rSges(4,3) - qJ(3)) * t215) * t216 + t234;
t81 = t245 * t219;
t80 = t245 * t216;
t77 = -t218 * t105 - t146 * t267;
t76 = t103 * t218 + t146 * t268;
t75 = t225 + t105 + t253;
t74 = (-rSges(5,3) * t215 + t241) * t216 + t233 + t242;
t73 = t239 * t219;
t72 = t239 * t216;
t61 = (t103 * t219 - t105 * t216) * t215;
t56 = t216 * (rSges(4,3) * t268 - t234) + t219 * t243 + t257;
t55 = t104 + t223 - t246;
t54 = -t279 + (-rSges(6,1) - pkin(4)) * t161 + (-rSges(6,2) * t215 + t241) * t216 + t238;
t53 = t236 * t219;
t52 = t236 * t216;
t51 = t263 * t218 + (-t145 - t170) * t267;
t50 = t102 * t218 + t145 * t268 + t261;
t49 = t218 * t69 + t267 * t89;
t48 = -t218 * t68 - t268 * t89;
t35 = t154 + (-t213 - t285) * t267 + t223 + t262;
t34 = (-pkin(4) - pkin(5)) * t161 + (t215 * t285 + t241) * t216 + t232 + t238;
t33 = t107 + (t102 * t219 + t216 * t263) * t215;
t31 = (t216 * t69 - t219 * t68) * t215;
t30 = t216 * t103 + t105 * t219 + t240;
t25 = t247 * t218 + (-t170 - t275) * t267;
t24 = t218 * t277 + t268 * t275 + t261;
t23 = t216 * t102 + t104 * t219 + t224;
t20 = t107 + (t216 * t247 + t219 * t277) * t215;
t19 = t43 * t216 - t219 * t42;
t18 = t41 * t216 - t219 * t40;
t17 = t39 * t216 - t219 * t38;
t16 = t37 * t216 - t219 * t36;
t11 = t216 * t277 + t219 * t276 + t224;
t5 = -t14 * t219 + t15 * t216;
t4 = -t12 * t219 + t13 * t216;
t6 = [Icges(2,3) + (t274 - (Icges(4,5) * t212 - Icges(4,6) * t211) * t215 + (Icges(3,2) + Icges(4,3)) * t218 + t303) * t218 + (Icges(3,1) * t215 - t142 * t201 - t157 * t211 + t158 * t212 + t273) * t215 + m(7) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t74 ^ 2 + t75 ^ 2) + m(6) * (t54 ^ 2 + t55 ^ 2) + m(4) * (t82 ^ 2 + t83 ^ 2) + m(3) * (t135 ^ 2 + t136 ^ 2) + m(2) * (t188 ^ 2 + t189 ^ 2) + t250 + t297; m(7) * (t34 * t53 + t35 * t52) + m(6) * (t54 * t73 + t55 * t72) + m(5) * (t74 * t81 + t75 * t80) + m(4) * (t128 * t83 + t129 * t82) + (-t176 * t157 / 0.2e1 - t177 * t158 / 0.2e1 - t135 * t284 + t183 * t286 + (Icges(3,6) * t286 + t230 * t290 + t122 / 0.2e1) * t218 - t220) * t219 + (t178 * t157 / 0.2e1 + t179 * t158 / 0.2e1 - t136 * t284 + t183 * t289 + (Icges(3,6) * t289 + t230 * t286 - t123 / 0.2e1) * t218 + t221) * t216 + ((Icges(3,5) * t216 - t125 * t211 + t127 * t212 + t219 * t231) * t289 + (-Icges(3,5) * t219 - t124 * t211 + t126 * t212 + t216 * t231) * t287) * t215; m(7) * (t11 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t23 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t30 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t128 ^ 2 + t129 ^ 2 + t56 ^ 2) + m(3) * (t187 ^ 2 * t252 + t119 ^ 2) + (-t210 * t164 - t16 - t17 - t4 + (t122 * t268 + t176 * t124 + t177 * t126) * t219) * t219 + (t5 + t19 + t18 + t209 * t165 + (t123 * t267 + t178 * t125 + t179 * t127) * t216 + (-t122 * t267 - t123 * t268 - t178 * t124 - t125 * t176 - t179 * t126 - t127 * t177 - t216 * t164 + t219 * t165) * t219) * t216; ((t216 * t35 + t219 * t34) * t291 + (t216 * t75 + t219 * t74) * t293 + (t216 * t55 + t219 * t54) * t292 + (t216 * t83 + t219 * t82) * t294) * t295; m(7) * (-t218 * t11 + (t216 * t52 + t219 * t53) * t215) + m(6) * (-t218 * t23 + (t216 * t72 + t219 * t73) * t215) + m(5) * (-t218 * t30 + (t216 * t80 + t219 * t81) * t215) + m(4) * (-t218 * t56 + (t128 * t216 + t129 * t219) * t215); 0.2e1 * (t294 + t293 + t251) * (t208 * t252 + t218 ^ 2); -t32 - t298 + m(7) * (t24 * t34 + t25 * t35) + m(5) * (t74 * t76 + t75 * t77) + m(6) * (t50 * t54 + t51 * t55) + (t216 * t220 + t219 * t221) * t215; m(7) * (t11 * t20 + t24 * t53 + t25 * t52) + m(6) * (t23 * t33 + t50 * t73 + t51 * t72) + m(5) * (t30 * t61 + t76 * t81 + t77 * t80) + ((t5 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t219 + (t4 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1) * t216) * t215 - t222 + t301 * t289 - (t216 * t299 + t219 * t300) * t218 / 0.2e1 + t302 * t287; m(5) * (-t61 * t218 + (t216 * t77 + t219 * t76) * t215) + m(6) * (-t33 * t218 + (t216 * t51 + t219 * t50) * t215) + m(7) * (-t20 * t218 + (t216 * t25 + t219 * t24) * t215); m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t33 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t61 ^ 2 + t76 ^ 2 + t77 ^ 2) + (-t3 + t298) * t218 + ((-t299 * t218 + t2 + t301) * t219 + (t300 * t218 + t1 + t302) * t216) * t215; m(7) * (t160 * t35 + t162 * t34) + m(6) * (t160 * t55 + t162 * t54); m(7) * (t11 * t271 + t160 * t52 + t162 * t53) + m(6) * (t160 * t72 + t162 * t73 + t23 * t271); t251 * (t160 * t216 + t162 * t219 - t201 * t218) * t295; m(7) * (t160 * t25 + t162 * t24 + t20 * t271) + m(6) * (t160 * t51 + t162 * t50 + t271 * t33); 0.2e1 * t251 * (t201 ^ 2 * t208 + t160 ^ 2 + t162 ^ 2); m(7) * (t34 * t48 + t35 * t49) + t32 + (t216 * t248 + t219 * t249) * t215; m(7) * (t11 * t31 + t48 * t53 + t49 * t52) + (t287 * t5 + t290 * t4) * t215 + t222; m(7) * (-t31 * t218 + (t216 * t49 + t219 * t48) * t215); m(7) * (t20 * t31 + t24 * t48 + t25 * t49) - t296; m(7) * (t160 * t49 + t162 * t48 + t271 * t31); m(7) * (t31 ^ 2 + t48 ^ 2 + t49 ^ 2) + t296;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
