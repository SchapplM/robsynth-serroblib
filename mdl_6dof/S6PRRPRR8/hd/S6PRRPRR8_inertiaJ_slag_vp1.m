% Calculate joint inertia matrix for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:36:06
% EndTime: 2019-03-08 22:36:21
% DurationCPUTime: 6.54s
% Computational Cost: add. (37392->600), mult. (103760->865), div. (0->0), fcn. (136172->14), ass. (0->270)
t294 = m(5) + m(6) + m(7);
t228 = sin(pkin(12));
t230 = cos(pkin(12));
t235 = sin(qJ(2));
t231 = cos(pkin(6));
t237 = cos(qJ(2));
t277 = t231 * t237;
t220 = -t228 * t235 + t230 * t277;
t278 = t231 * t235;
t221 = t228 * t237 + t230 * t278;
t234 = sin(qJ(3));
t229 = sin(pkin(6));
t282 = sin(pkin(7));
t287 = cos(qJ(3));
t242 = t287 * t282;
t240 = t229 * t242;
t283 = cos(pkin(7));
t243 = t283 * t287;
t189 = -t220 * t243 + t221 * t234 + t230 * t240;
t253 = t229 * t283;
t210 = -t220 * t282 - t230 * t253;
t233 = sin(qJ(5));
t286 = cos(qJ(5));
t167 = -t189 * t286 + t210 * t233;
t222 = -t228 * t277 - t230 * t235;
t223 = -t228 * t278 + t230 * t237;
t191 = -t222 * t243 + t223 * t234 - t228 * t240;
t211 = -t222 * t282 + t228 * t253;
t169 = -t191 * t286 + t211 * t233;
t279 = t229 * t235;
t208 = -t229 * t237 * t243 - t231 * t242 + t234 * t279;
t252 = t229 * t282;
t219 = t231 * t283 - t237 * t252;
t193 = -t208 * t286 + t219 * t233;
t168 = t189 * t233 + t210 * t286;
t250 = t234 * t282;
t251 = t234 * t283;
t280 = t229 * t230;
t190 = t220 * t251 + t221 * t287 - t250 * t280;
t232 = sin(qJ(6));
t236 = cos(qJ(6));
t125 = -t168 * t232 + t190 * t236;
t126 = t168 * t236 + t190 * t232;
t92 = Icges(7,5) * t126 + Icges(7,6) * t125 + Icges(7,3) * t167;
t94 = Icges(7,4) * t126 + Icges(7,2) * t125 + Icges(7,6) * t167;
t96 = Icges(7,1) * t126 + Icges(7,4) * t125 + Icges(7,5) * t167;
t28 = t125 * t94 + t126 * t96 + t167 * t92;
t170 = t191 * t233 + t211 * t286;
t192 = t223 * t287 + (t222 * t283 + t228 * t252) * t234;
t127 = -t170 * t232 + t192 * t236;
t128 = t170 * t236 + t192 * t232;
t93 = Icges(7,5) * t128 + Icges(7,6) * t127 + Icges(7,3) * t169;
t95 = Icges(7,4) * t128 + Icges(7,2) * t127 + Icges(7,6) * t169;
t97 = Icges(7,1) * t128 + Icges(7,4) * t127 + Icges(7,5) * t169;
t29 = t125 * t95 + t126 * t97 + t167 * t93;
t194 = t208 * t233 + t219 * t286;
t209 = t231 * t250 + (t235 * t287 + t237 * t251) * t229;
t164 = -t194 * t232 + t209 * t236;
t165 = t194 * t236 + t209 * t232;
t111 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t193;
t114 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t193;
t117 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t193;
t47 = t111 * t167 + t114 * t125 + t117 * t126;
t1 = t167 * t28 + t169 * t29 + t193 * t47;
t293 = t1 / 0.2e1;
t30 = t127 * t94 + t128 * t96 + t169 * t92;
t31 = t127 * t95 + t128 * t97 + t169 * t93;
t48 = t111 * t169 + t114 * t127 + t117 * t128;
t2 = t167 * t30 + t169 * t31 + t193 * t48;
t292 = t2 / 0.2e1;
t38 = t164 * t94 + t165 * t96 + t193 * t92;
t39 = t164 * t95 + t165 * t97 + t193 * t93;
t54 = t111 * t193 + t114 * t164 + t117 * t165;
t9 = t167 * t38 + t169 * t39 + t193 * t54;
t291 = t9 / 0.2e1;
t290 = t167 / 0.2e1;
t289 = t169 / 0.2e1;
t288 = t193 / 0.2e1;
t98 = rSges(7,1) * t126 + rSges(7,2) * t125 + rSges(7,3) * t167;
t285 = pkin(5) * t168 + pkin(11) * t167 + t98;
t99 = rSges(7,1) * t128 + rSges(7,2) * t127 + rSges(7,3) * t169;
t284 = pkin(5) * t170 + pkin(11) * t169 + t99;
t281 = t228 * t229;
t120 = rSges(7,1) * t165 + rSges(7,2) * t164 + rSges(7,3) * t193;
t276 = pkin(5) * t194 + pkin(11) * t193 + t120;
t146 = rSges(5,1) * t210 - rSges(5,2) * t190 + rSges(5,3) * t189;
t156 = pkin(3) * t190 + qJ(4) * t189;
t275 = -t146 - t156;
t149 = t211 * t156;
t171 = pkin(4) * t210 + pkin(10) * t190;
t274 = t171 * t211 + t149;
t157 = pkin(3) * t192 + qJ(4) * t191;
t151 = t219 * t157;
t172 = pkin(4) * t211 + pkin(10) * t192;
t273 = t172 * t219 + t151;
t198 = pkin(2) * t223 + pkin(9) * t211;
t195 = t231 * t198;
t272 = t157 * t231 + t195;
t271 = -t156 - t171;
t270 = -t157 - t172;
t184 = pkin(3) * t209 + qJ(4) * t208;
t161 = t210 * t184;
t196 = pkin(4) * t219 + pkin(10) * t209;
t269 = t196 * t210 + t161;
t180 = rSges(5,1) * t219 - rSges(5,2) * t209 + rSges(5,3) * t208;
t268 = -t180 - t184;
t267 = -t184 - t196;
t197 = t221 * pkin(2) + pkin(9) * t210;
t266 = t197 * t281 + t198 * t280;
t112 = Icges(6,5) * t168 - Icges(6,6) * t167 + Icges(6,3) * t190;
t115 = Icges(6,4) * t168 - Icges(6,2) * t167 + Icges(6,6) * t190;
t118 = Icges(6,1) * t168 - Icges(6,4) * t167 + Icges(6,5) * t190;
t55 = t112 * t190 - t115 * t167 + t118 * t168;
t113 = Icges(6,5) * t170 - Icges(6,6) * t169 + Icges(6,3) * t192;
t116 = Icges(6,4) * t170 - Icges(6,2) * t169 + Icges(6,6) * t192;
t119 = Icges(6,1) * t170 - Icges(6,4) * t169 + Icges(6,5) * t192;
t56 = t113 * t190 - t116 * t167 + t119 * t168;
t141 = Icges(6,5) * t194 - Icges(6,6) * t193 + Icges(6,3) * t209;
t142 = Icges(6,4) * t194 - Icges(6,2) * t193 + Icges(6,6) * t209;
t143 = Icges(6,1) * t194 - Icges(6,4) * t193 + Icges(6,5) * t209;
t67 = t141 * t190 - t142 * t167 + t143 * t168;
t13 = t190 * t55 + t192 * t56 + t209 * t67;
t3 = t190 * t28 + t192 * t29 + t209 * t47;
t265 = t3 / 0.2e1 + t13 / 0.2e1;
t57 = t112 * t192 - t115 * t169 + t118 * t170;
t58 = t113 * t192 - t116 * t169 + t119 * t170;
t68 = t141 * t192 - t142 * t169 + t143 * t170;
t14 = t190 * t57 + t192 * t58 + t209 * t68;
t4 = t190 * t30 + t192 * t31 + t209 * t48;
t264 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t210 * t55 + t211 * t56 + t219 * t67;
t5 = t210 * t28 + t211 * t29 + t219 * t47;
t263 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t210 * t57 + t211 * t58 + t219 * t68;
t6 = t210 * t30 + t211 * t31 + t219 * t48;
t262 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t67 * t231 + (t228 * t56 - t230 * t55) * t229;
t7 = t47 * t231 + (t228 * t29 - t230 * t28) * t229;
t261 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t68 * t231 + (t228 * t58 - t230 * t57) * t229;
t8 = t48 * t231 + (t228 * t31 - t230 * t30) * t229;
t260 = t8 / 0.2e1 + t18 / 0.2e1;
t10 = t190 * t38 + t192 * t39 + t209 * t54;
t59 = t112 * t209 - t115 * t193 + t118 * t194;
t60 = t113 * t209 - t116 * t193 + t119 * t194;
t80 = t141 * t209 - t142 * t193 + t143 * t194;
t19 = t190 * t59 + t192 * t60 + t209 * t80;
t259 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = t210 * t38 + t211 * t39 + t219 * t54;
t20 = t210 * t59 + t211 * t60 + t219 * t80;
t258 = t11 / 0.2e1 + t20 / 0.2e1;
t12 = t54 * t231 + (t228 * t39 - t230 * t38) * t229;
t21 = t80 * t231 + (t228 * t60 - t230 * t59) * t229;
t257 = t12 / 0.2e1 + t21 / 0.2e1;
t121 = rSges(6,1) * t168 - rSges(6,2) * t167 + rSges(6,3) * t190;
t256 = -t121 + t271;
t148 = rSges(6,1) * t194 - rSges(6,2) * t193 + rSges(6,3) * t209;
t255 = -t148 + t267;
t254 = t172 * t231 + t272;
t179 = rSges(4,1) * t209 - rSges(4,2) * t208 + rSges(4,3) * t219;
t212 = pkin(2) * t279 + pkin(9) * t219;
t249 = (-t179 - t212) * t229;
t247 = t271 - t285;
t246 = t267 - t276;
t245 = t156 * t281 + t157 * t280 + t266;
t244 = (-t212 + t268) * t229;
t241 = (-t212 + t255) * t229;
t239 = t171 * t281 + t172 * t280 + t245;
t238 = (-t212 + t246) * t229;
t218 = t231 * rSges(3,3) + (rSges(3,1) * t235 + rSges(3,2) * t237) * t229;
t217 = Icges(3,5) * t231 + (Icges(3,1) * t235 + Icges(3,4) * t237) * t229;
t216 = Icges(3,6) * t231 + (Icges(3,4) * t235 + Icges(3,2) * t237) * t229;
t215 = Icges(3,3) * t231 + (Icges(3,5) * t235 + Icges(3,6) * t237) * t229;
t206 = rSges(3,1) * t223 + rSges(3,2) * t222 + rSges(3,3) * t281;
t205 = rSges(3,1) * t221 + rSges(3,2) * t220 - rSges(3,3) * t280;
t204 = Icges(3,1) * t223 + Icges(3,4) * t222 + Icges(3,5) * t281;
t203 = Icges(3,1) * t221 + Icges(3,4) * t220 - Icges(3,5) * t280;
t202 = Icges(3,4) * t223 + Icges(3,2) * t222 + Icges(3,6) * t281;
t201 = Icges(3,4) * t221 + Icges(3,2) * t220 - Icges(3,6) * t280;
t200 = Icges(3,5) * t223 + Icges(3,6) * t222 + Icges(3,3) * t281;
t199 = Icges(3,5) * t221 + Icges(3,6) * t220 - Icges(3,3) * t280;
t183 = -t205 * t231 - t218 * t280;
t182 = t206 * t231 - t218 * t281;
t178 = Icges(4,1) * t209 - Icges(4,4) * t208 + Icges(4,5) * t219;
t177 = Icges(5,1) * t219 - Icges(5,4) * t209 + Icges(5,5) * t208;
t176 = Icges(4,4) * t209 - Icges(4,2) * t208 + Icges(4,6) * t219;
t175 = Icges(5,4) * t219 - Icges(5,2) * t209 + Icges(5,6) * t208;
t174 = Icges(4,5) * t209 - Icges(4,6) * t208 + Icges(4,3) * t219;
t173 = Icges(5,5) * t219 - Icges(5,6) * t209 + Icges(5,3) * t208;
t160 = (t205 * t228 + t206 * t230) * t229;
t147 = rSges(5,1) * t211 - rSges(5,2) * t192 + rSges(5,3) * t191;
t145 = rSges(4,1) * t192 - rSges(4,2) * t191 + rSges(4,3) * t211;
t144 = rSges(4,1) * t190 - rSges(4,2) * t189 + rSges(4,3) * t210;
t140 = Icges(4,1) * t192 - Icges(4,4) * t191 + Icges(4,5) * t211;
t139 = Icges(4,1) * t190 - Icges(4,4) * t189 + Icges(4,5) * t210;
t138 = Icges(5,1) * t211 - Icges(5,4) * t192 + Icges(5,5) * t191;
t137 = Icges(5,1) * t210 - Icges(5,4) * t190 + Icges(5,5) * t189;
t136 = Icges(4,4) * t192 - Icges(4,2) * t191 + Icges(4,6) * t211;
t135 = Icges(4,4) * t190 - Icges(4,2) * t189 + Icges(4,6) * t210;
t134 = Icges(5,4) * t211 - Icges(5,2) * t192 + Icges(5,6) * t191;
t133 = Icges(5,4) * t210 - Icges(5,2) * t190 + Icges(5,6) * t189;
t132 = Icges(4,5) * t192 - Icges(4,6) * t191 + Icges(4,3) * t211;
t131 = Icges(4,5) * t190 - Icges(4,6) * t189 + Icges(4,3) * t210;
t130 = Icges(5,5) * t211 - Icges(5,6) * t192 + Icges(5,3) * t191;
t129 = Icges(5,5) * t210 - Icges(5,6) * t190 + Icges(5,3) * t189;
t122 = rSges(6,1) * t170 - rSges(6,2) * t169 + rSges(6,3) * t192;
t110 = t145 * t219 - t179 * t211;
t109 = -t144 * t219 + t179 * t210;
t108 = (-t144 - t197) * t231 + t230 * t249;
t107 = t231 * t145 + t228 * t249 + t195;
t106 = t173 * t208 - t175 * t209 + t177 * t219;
t105 = t174 * t219 - t176 * t208 + t178 * t209;
t104 = t144 * t211 - t145 * t210;
t103 = t173 * t191 - t175 * t192 + t177 * t211;
t102 = t173 * t189 - t175 * t190 + t177 * t210;
t101 = t174 * t211 - t176 * t191 + t178 * t192;
t100 = t174 * t210 - t176 * t189 + t178 * t190;
t91 = (t144 * t228 + t145 * t230) * t229 + t266;
t90 = t122 * t209 - t148 * t192;
t89 = -t121 * t209 + t148 * t190;
t88 = t219 * t147 + t211 * t268 + t151;
t87 = t210 * t180 + t219 * t275 + t161;
t86 = (-t197 + t275) * t231 + t230 * t244;
t85 = t231 * t147 + t228 * t244 + t272;
t84 = t130 * t208 - t134 * t209 + t138 * t219;
t83 = t129 * t208 - t133 * t209 + t137 * t219;
t82 = t132 * t219 - t136 * t208 + t140 * t209;
t81 = t131 * t219 - t135 * t208 + t139 * t209;
t79 = t130 * t191 - t134 * t192 + t138 * t211;
t78 = t129 * t191 - t133 * t192 + t137 * t211;
t77 = t130 * t189 - t134 * t190 + t138 * t210;
t76 = t129 * t189 - t133 * t190 + t137 * t210;
t75 = t132 * t211 - t136 * t191 + t140 * t192;
t74 = t131 * t211 - t135 * t191 + t139 * t192;
t73 = t132 * t210 - t136 * t189 + t140 * t190;
t72 = t131 * t210 - t135 * t189 + t139 * t190;
t71 = t121 * t192 - t122 * t190;
t70 = t211 * t146 + t149 + (-t147 - t157) * t210;
t69 = (t146 * t228 + t147 * t230) * t229 + t245;
t66 = -t120 * t169 + t193 * t99;
t65 = t120 * t167 - t193 * t98;
t64 = (-t197 + t256) * t231 + t230 * t241;
t63 = t231 * t122 + t228 * t241 + t254;
t62 = t219 * t122 + t211 * t255 + t273;
t61 = t210 * t148 + t219 * t256 + t269;
t53 = -t167 * t99 + t169 * t98;
t52 = (t121 * t228 + t122 * t230) * t229 + t239;
t51 = t211 * t121 + (-t122 + t270) * t210 + t274;
t50 = -t192 * t276 + t209 * t284;
t49 = t190 * t276 - t209 * t285;
t46 = t106 * t231 + (t228 * t84 - t230 * t83) * t229;
t45 = t105 * t231 + (t228 * t82 - t230 * t81) * t229;
t44 = (-t197 + t247) * t231 + t230 * t238;
t43 = t228 * t238 + t231 * t284 + t254;
t42 = -t190 * t284 + t192 * t285;
t41 = t106 * t219 + t210 * t83 + t211 * t84;
t40 = t105 * t219 + t210 * t81 + t211 * t82;
t37 = t211 * t246 + t219 * t284 + t273;
t36 = t210 * t276 + t219 * t247 + t269;
t35 = t103 * t231 + (t228 * t79 - t230 * t78) * t229;
t34 = t102 * t231 + (t228 * t77 - t230 * t76) * t229;
t33 = t101 * t231 + (t228 * t75 - t230 * t74) * t229;
t32 = t100 * t231 + (t228 * t73 - t230 * t72) * t229;
t27 = t103 * t219 + t210 * t78 + t211 * t79;
t26 = t102 * t219 + t210 * t76 + t211 * t77;
t25 = t101 * t219 + t210 * t74 + t211 * t75;
t24 = t100 * t219 + t210 * t72 + t211 * t73;
t23 = (t228 * t285 + t230 * t284) * t229 + t239;
t22 = t285 * t211 + (t270 - t284) * t210 + t274;
t123 = [m(3) + m(4) + m(2) + t294; m(3) * t160 + m(4) * t91 + m(5) * t69 + m(6) * t52 + m(7) * t23; m(7) * (t23 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t52 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t69 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t107 ^ 2 + t108 ^ 2 + t91 ^ 2) + m(3) * (t160 ^ 2 + t182 ^ 2 + t183 ^ 2) + (t8 + t18 + t35 + t33 + (t200 * t281 + t202 * t222 + t204 * t223) * t281) * t281 + (-t7 - t17 - t32 - t34 + (-t199 * t280 + t201 * t220 + t203 * t221) * t280 + (-t199 * t281 + t200 * t280 - t201 * t222 - t202 * t220 - t203 * t223 - t204 * t221) * t281) * t280 + (t12 + (t215 * t281 + t222 * t216 + t223 * t217) * t281 - (-t215 * t280 + t220 * t216 + t221 * t217) * t280 + t21 + t46 + t45 + ((t202 * t237 + t204 * t235) * t228 - (t201 * t237 + t203 * t235) * t230) * t229 ^ 2 + ((-t199 * t230 + t200 * t228 + t216 * t237 + t217 * t235) * t229 + t231 * t215) * t231) * t231; m(4) * t104 + m(5) * t70 + m(6) * t51 + m(7) * t22; (t41 / 0.2e1 + t40 / 0.2e1 + t258) * t231 + (t46 / 0.2e1 + t45 / 0.2e1 + t257) * t219 + (t35 / 0.2e1 + t33 / 0.2e1 + t260) * t211 + (t34 / 0.2e1 + t32 / 0.2e1 + t261) * t210 + m(7) * (t22 * t23 + t36 * t44 + t37 * t43) + m(6) * (t51 * t52 + t61 * t64 + t62 * t63) + m(5) * (t69 * t70 + t85 * t88 + t86 * t87) + m(4) * (t104 * t91 + t107 * t110 + t108 * t109) + ((-t24 / 0.2e1 - t26 / 0.2e1 - t263) * t230 + (t25 / 0.2e1 + t27 / 0.2e1 + t262) * t228) * t229; (t11 + t20 + t41 + t40) * t219 + (t6 + t16 + t27 + t25) * t211 + (t5 + t15 + t26 + t24) * t210 + m(7) * (t22 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t51 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t70 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t104 ^ 2 + t109 ^ 2 + t110 ^ 2); t208 * t294; m(7) * (t189 * t43 + t191 * t44 + t208 * t23) + m(6) * (t189 * t63 + t191 * t64 + t208 * t52) + m(5) * (t189 * t85 + t191 * t86 + t208 * t69); m(7) * (t189 * t37 + t191 * t36 + t208 * t22) + m(6) * (t189 * t62 + t191 * t61 + t208 * t51) + m(5) * (t189 * t88 + t191 * t87 + t208 * t70); (t189 ^ 2 + t191 ^ 2 + t208 ^ 2) * t294; m(6) * t71 + m(7) * t42; t259 * t231 + t257 * t209 + t260 * t192 + t261 * t190 + m(7) * (t23 * t42 + t43 * t50 + t44 * t49) + m(6) * (t52 * t71 + t63 * t90 + t64 * t89) + (t228 * t264 - t230 * t265) * t229; t259 * t219 + t264 * t211 + t265 * t210 + t258 * t209 + t262 * t192 + t263 * t190 + m(7) * (t22 * t42 + t36 * t49 + t37 * t50) + m(6) * (t51 * t71 + t61 * t89 + t62 * t90); m(6) * (t189 * t90 + t191 * t89 + t208 * t71) + m(7) * (t189 * t50 + t191 * t49 + t208 * t42); (t10 + t19) * t209 + (t4 + t14) * t192 + (t3 + t13) * t190 + m(7) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t71 ^ 2 + t89 ^ 2 + t90 ^ 2); m(7) * t53; m(7) * (t23 * t53 + t43 * t66 + t44 * t65) + t7 * t290 + t8 * t289 + t12 * t288 + t231 * t291 + (-t230 * t1 / 0.2e1 + t228 * t292) * t229; t210 * t293 + m(7) * (t22 * t53 + t36 * t65 + t37 * t66) + t219 * t291 + t5 * t290 + t6 * t289 + t211 * t292 + t11 * t288; m(7) * (t189 * t66 + t191 * t65 + t208 * t53); t190 * t293 + m(7) * (t42 * t53 + t49 * t65 + t50 * t66) + t192 * t292 + t3 * t290 + t4 * t289 + t10 * t288 + t209 * t291; t169 * t2 + t167 * t1 + t193 * t9 + m(7) * (t53 ^ 2 + t65 ^ 2 + t66 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t123(1) t123(2) t123(4) t123(7) t123(11) t123(16); t123(2) t123(3) t123(5) t123(8) t123(12) t123(17); t123(4) t123(5) t123(6) t123(9) t123(13) t123(18); t123(7) t123(8) t123(9) t123(10) t123(14) t123(19); t123(11) t123(12) t123(13) t123(14) t123(15) t123(20); t123(16) t123(17) t123(18) t123(19) t123(20) t123(21);];
Mq  = res;
