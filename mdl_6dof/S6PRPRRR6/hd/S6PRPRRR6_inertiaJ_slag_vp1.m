% Calculate joint inertia matrix for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:15
% EndTime: 2019-03-08 20:45:25
% DurationCPUTime: 5.53s
% Computational Cost: add. (18839->517), mult. (41222->745), div. (0->0), fcn. (52572->12), ass. (0->235)
t285 = Icges(3,1) + Icges(4,2);
t283 = Icges(3,4) + Icges(4,6);
t282 = Icges(3,5) - Icges(4,4);
t284 = Icges(3,2) + Icges(4,3);
t281 = Icges(3,6) - Icges(4,5);
t280 = Icges(3,3) + Icges(4,1);
t220 = cos(pkin(6));
t219 = cos(pkin(11));
t225 = cos(qJ(2));
t248 = t225 * t219;
t217 = sin(pkin(11));
t223 = sin(qJ(2));
t251 = t223 * t217;
t203 = -t220 * t248 + t251;
t249 = t225 * t217;
t250 = t223 * t219;
t204 = t220 * t250 + t249;
t218 = sin(pkin(6));
t255 = t218 * t219;
t279 = t284 * t203 - t283 * t204 + t281 * t255;
t278 = t283 * t203 - t285 * t204 + t282 * t255;
t205 = t220 * t249 + t250;
t206 = -t220 * t251 + t248;
t256 = t217 * t218;
t277 = t284 * t205 - t283 * t206 - t281 * t256;
t276 = -t283 * t205 + t285 * t206 + t282 * t256;
t275 = t281 * t203 - t282 * t204 + t280 * t255;
t274 = -t281 * t205 + t282 * t206 + t280 * t256;
t273 = t280 * t220 + (t282 * t223 + t281 * t225) * t218;
t272 = t281 * t220 + (t283 * t223 + t284 * t225) * t218;
t271 = t282 * t220 + (t285 * t223 + t283 * t225) * t218;
t222 = sin(qJ(4));
t254 = t218 * t222;
t263 = cos(qJ(4));
t186 = -t205 * t263 + t217 * t254;
t269 = t186 / 0.2e1;
t188 = t203 * t263 + t219 * t254;
t268 = -t188 / 0.2e1;
t267 = t204 / 0.2e1;
t266 = t206 / 0.2e1;
t236 = t218 * t263;
t207 = t220 * t222 + t225 * t236;
t265 = t207 / 0.2e1;
t264 = t220 / 0.2e1;
t224 = cos(qJ(5));
t262 = pkin(5) * t224;
t187 = t205 * t222 + t217 * t236;
t216 = qJ(5) + qJ(6);
t214 = sin(t216);
t215 = cos(t216);
t145 = -t187 * t214 + t206 * t215;
t146 = t187 * t215 + t206 * t214;
t103 = rSges(7,1) * t146 + rSges(7,2) * t145 + rSges(7,3) * t186;
t221 = sin(qJ(5));
t257 = t206 * t221;
t95 = pkin(5) * t257 + pkin(10) * t186 + t187 * t262;
t260 = t103 + t95;
t189 = t203 * t222 - t219 * t236;
t147 = -t189 * t214 + t204 * t215;
t148 = t189 * t215 + t204 * t214;
t104 = rSges(7,1) * t148 + rSges(7,2) * t147 - rSges(7,3) * t188;
t258 = t204 * t221;
t96 = pkin(5) * t258 - pkin(10) * t188 + t189 * t262;
t259 = t104 + t96;
t253 = t218 * t223;
t252 = t218 * t225;
t208 = t220 * t263 - t222 * t252;
t181 = -t208 * t214 + t215 * t253;
t182 = t208 * t215 + t214 * t253;
t122 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t207;
t240 = t221 * t253;
t123 = pkin(5) * t240 + pkin(10) * t207 + t208 * t262;
t247 = t122 + t123;
t178 = pkin(2) * t204 + qJ(3) * t203;
t179 = pkin(2) * t206 + qJ(3) * t205;
t246 = t178 * t256 + t179 * t255;
t177 = t220 * t179;
t192 = pkin(3) * t256 + pkin(8) * t206;
t245 = t220 * t192 + t177;
t193 = -pkin(3) * t255 + pkin(8) * t204;
t244 = -t178 - t193;
t209 = (pkin(2) * t223 - qJ(3) * t225) * t218;
t243 = -pkin(3) * t220 - pkin(8) * t253 - t209;
t101 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t186;
t97 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t186;
t99 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t186;
t58 = t101 * t182 + t181 * t99 + t207 * t97;
t100 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t188;
t102 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t188;
t98 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t188;
t59 = t100 * t181 + t102 * t182 + t207 * t98;
t119 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t207;
t120 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t207;
t121 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t207;
t72 = t119 * t207 + t120 * t181 + t121 * t182;
t26 = t186 * t58 - t188 * t59 + t207 * t72;
t47 = t101 * t146 + t145 * t99 + t186 * t97;
t48 = t100 * t145 + t102 * t146 + t186 * t98;
t65 = t119 * t186 + t120 * t145 + t121 * t146;
t7 = t186 * t47 - t188 * t48 + t207 * t65;
t49 = t101 * t148 + t147 * t99 - t188 * t97;
t50 = t100 * t147 + t102 * t148 - t188 * t98;
t66 = -t119 * t188 + t120 * t147 + t121 * t148;
t8 = t186 * t49 - t188 * t50 + t207 * t66;
t242 = t186 * t7 - t188 * t8 + t207 * t26;
t241 = -m(4) - m(5) - m(6) - m(7);
t143 = t187 * pkin(4) + t186 * pkin(9);
t239 = t220 * t143 + t245;
t144 = t189 * pkin(4) - t188 * pkin(9);
t238 = -t144 + t244;
t180 = t208 * pkin(4) + t207 * pkin(9);
t237 = -t180 + t243;
t235 = t253 / 0.2e1;
t234 = (-rSges(4,1) * t220 - (-rSges(4,2) * t223 - rSges(4,3) * t225) * t218 - t209) * t218;
t233 = t192 * t255 + t193 * t256 + t246;
t173 = rSges(5,1) * t208 - rSges(5,2) * t207 + rSges(5,3) * t253;
t232 = (-t173 + t243) * t218;
t13 = t204 * t48 + t206 * t47 + t253 * t65;
t14 = t204 * t50 + t206 * t49 + t253 * t66;
t28 = t204 * t59 + t206 * t58 + t253 * t72;
t231 = t13 * t269 + t14 * t268 + t26 * t235 + t28 * t265 + t7 * t266 + t8 * t267;
t15 = t220 * t65 + (t217 * t47 - t219 * t48) * t218;
t16 = t220 * t66 + (t217 * t49 - t219 * t50) * t218;
t30 = t220 * t72 + (t217 * t58 - t219 * t59) * t218;
t230 = t15 * t269 + t16 * t268 + t26 * t264 + t30 * t265 + t7 * t256 / 0.2e1 - t8 * t255 / 0.2e1;
t190 = -t208 * t221 + t224 * t253;
t191 = t208 * t224 + t240;
t136 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t207;
t229 = (-t136 + t237) * t218;
t228 = t143 * t255 + t144 * t256 + t233;
t227 = (t237 - t247) * t218;
t200 = rSges(3,3) * t220 + (rSges(3,1) * t223 + rSges(3,2) * t225) * t218;
t172 = Icges(5,1) * t208 - Icges(5,4) * t207 + Icges(5,5) * t253;
t171 = Icges(5,4) * t208 - Icges(5,2) * t207 + Icges(5,6) * t253;
t170 = Icges(5,5) * t208 - Icges(5,6) * t207 + Icges(5,3) * t253;
t169 = rSges(3,1) * t206 - rSges(3,2) * t205 + rSges(3,3) * t256;
t168 = rSges(3,1) * t204 - rSges(3,2) * t203 - rSges(3,3) * t255;
t167 = -rSges(4,1) * t255 - rSges(4,2) * t204 + rSges(4,3) * t203;
t166 = rSges(4,1) * t256 - rSges(4,2) * t206 + rSges(4,3) * t205;
t153 = t204 * t180;
t152 = t189 * t224 + t258;
t151 = -t189 * t221 + t204 * t224;
t150 = t187 * t224 + t257;
t149 = -t187 * t221 + t206 * t224;
t141 = t143 * t253;
t138 = -t168 * t220 - t200 * t255;
t137 = t169 * t220 - t200 * t256;
t135 = t206 * t144;
t134 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t207;
t133 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t207;
t132 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t207;
t131 = rSges(5,1) * t189 + rSges(5,2) * t188 + rSges(5,3) * t204;
t130 = rSges(5,1) * t187 - rSges(5,2) * t186 + rSges(5,3) * t206;
t129 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t204;
t128 = Icges(5,1) * t187 - Icges(5,4) * t186 + Icges(5,5) * t206;
t127 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t204;
t126 = Icges(5,4) * t187 - Icges(5,2) * t186 + Icges(5,6) * t206;
t125 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t204;
t124 = Icges(5,5) * t187 - Icges(5,6) * t186 + Icges(5,3) * t206;
t118 = (t168 * t217 + t169 * t219) * t218;
t117 = t188 * t122;
t116 = (-t167 - t178) * t220 + t219 * t234;
t115 = t166 * t220 + t217 * t234 + t177;
t114 = rSges(6,1) * t152 + rSges(6,2) * t151 - rSges(6,3) * t188;
t113 = rSges(6,1) * t150 + rSges(6,2) * t149 + rSges(6,3) * t186;
t112 = Icges(6,1) * t152 + Icges(6,4) * t151 - Icges(6,5) * t188;
t111 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t186;
t110 = Icges(6,4) * t152 + Icges(6,2) * t151 - Icges(6,6) * t188;
t109 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t186;
t108 = Icges(6,5) * t152 + Icges(6,6) * t151 - Icges(6,3) * t188;
t107 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t186;
t106 = t130 * t253 - t173 * t206;
t105 = -t131 * t253 + t173 * t204;
t94 = t207 * t103;
t93 = t170 * t253 - t171 * t207 + t172 * t208;
t92 = t186 * t104;
t91 = (t166 * t219 + t167 * t217) * t218 + t246;
t90 = -t130 * t204 + t131 * t206;
t89 = t170 * t204 + t171 * t188 + t172 * t189;
t88 = t170 * t206 - t171 * t186 + t172 * t187;
t87 = (-t131 + t244) * t220 + t219 * t232;
t86 = t130 * t220 + t217 * t232 + t245;
t85 = -t114 * t207 - t136 * t188;
t84 = t113 * t207 - t136 * t186;
t83 = -t104 * t207 - t117;
t82 = -t122 * t186 + t94;
t81 = t125 * t253 - t127 * t207 + t129 * t208;
t80 = t124 * t253 - t126 * t207 + t128 * t208;
t79 = t132 * t207 + t133 * t190 + t134 * t191;
t78 = t125 * t204 + t127 * t188 + t129 * t189;
t77 = t124 * t204 + t126 * t188 + t128 * t189;
t76 = t125 * t206 - t127 * t186 + t129 * t187;
t75 = t124 * t206 - t126 * t186 + t128 * t187;
t74 = (t130 * t219 + t131 * t217) * t218 + t233;
t73 = t113 * t188 + t114 * t186;
t71 = t103 * t188 + t92;
t70 = t113 * t253 + t141 + (-t136 - t180) * t206;
t69 = t136 * t204 + t153 + (-t114 - t144) * t253;
t68 = -t132 * t188 + t133 * t151 + t134 * t152;
t67 = t132 * t186 + t133 * t149 + t134 * t150;
t64 = (-t114 + t238) * t220 + t219 * t229;
t63 = t113 * t220 + t217 * t229 + t239;
t62 = t114 * t206 + t135 + (-t113 - t143) * t204;
t61 = t108 * t207 + t110 * t190 + t112 * t191;
t60 = t107 * t207 + t109 * t190 + t111 * t191;
t57 = -t108 * t188 + t110 * t151 + t112 * t152;
t56 = -t107 * t188 + t109 * t151 + t111 * t152;
t55 = t108 * t186 + t110 * t149 + t112 * t150;
t54 = t107 * t186 + t109 * t149 + t111 * t150;
t53 = -t123 * t188 - t207 * t259 - t117;
t52 = -t186 * t247 + t207 * t95 + t94;
t51 = (t113 * t219 + t114 * t217) * t218 + t228;
t46 = t141 + t260 * t253 + (-t180 - t247) * t206;
t45 = t153 + t247 * t204 + (-t144 - t259) * t253;
t44 = (t238 - t259) * t220 + t219 * t227;
t43 = t217 * t227 + t220 * t260 + t239;
t42 = t186 * t96 + t188 * t260 + t92;
t41 = t220 * t93 + (t217 * t80 - t219 * t81) * t218;
t40 = t204 * t81 + t206 * t80 + t253 * t93;
t39 = t135 + t259 * t206 + (-t143 - t260) * t204;
t38 = t220 * t89 + (t217 * t77 - t219 * t78) * t218;
t37 = t220 * t88 + (t217 * t75 - t219 * t76) * t218;
t36 = t204 * t78 + t206 * t77 + t253 * t89;
t35 = t204 * t76 + t206 * t75 + t253 * t88;
t34 = (t217 * t259 + t219 * t260) * t218 + t228;
t33 = t220 * t79 + (t217 * t60 - t219 * t61) * t218;
t32 = t204 * t61 + t206 * t60 + t253 * t79;
t31 = t186 * t60 - t188 * t61 + t207 * t79;
t22 = t220 * t68 + (t217 * t56 - t219 * t57) * t218;
t21 = t220 * t67 + (t217 * t54 - t219 * t55) * t218;
t20 = t204 * t57 + t206 * t56 + t253 * t68;
t19 = t204 * t55 + t206 * t54 + t253 * t67;
t18 = t186 * t56 - t188 * t57 + t207 * t68;
t17 = t186 * t54 - t188 * t55 + t207 * t67;
t1 = [m(2) + m(3) - t241; m(3) * t118 + m(4) * t91 + m(5) * t74 + m(6) * t51 + m(7) * t34; m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t51 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t74 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t115 ^ 2 + t116 ^ 2 + t91 ^ 2) + m(3) * (t118 ^ 2 + t137 ^ 2 + t138 ^ 2) + (t30 + t33 + t41 + t273 * t220 ^ 2 + ((t271 * t223 + t272 * t225) * t220 + (t275 * t220 + (t278 * t223 + t279 * t225) * t218) * t219 + (t274 * t220 + (t276 * t223 - t277 * t225) * t218) * t217) * t218) * t220 + (t15 + t21 + t37 + (t277 * t205 + t276 * t206 + t274 * t256) * t256 + (-t272 * t205 + t271 * t206 + t273 * t256) * t220) * t256 + (-t16 - t22 - t38 + (t279 * t203 - t278 * t204 + t275 * t255) * t255 + (t272 * t203 - t271 * t204 + t273 * t255) * t220 + (-t277 * t203 - t276 * t204 - t279 * t205 + t278 * t206 + t274 * t255 + t275 * t256) * t256) * t255; t241 * t252; m(7) * (t203 * t43 + t205 * t44 - t252 * t34) + m(6) * (t203 * t63 + t205 * t64 - t252 * t51) + m(5) * (t203 * t86 + t205 * t87 - t252 * t74) + m(4) * (t115 * t203 + t116 * t205 - t252 * t91); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t218 ^ 2 * t225 ^ 2 + t203 ^ 2 + t205 ^ 2); m(5) * t90 + m(6) * t62 + m(7) * t39; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1) * t220 + (t15 / 0.2e1 + t21 / 0.2e1 + t37 / 0.2e1) * t206 + (t16 / 0.2e1 + t22 / 0.2e1 + t38 / 0.2e1) * t204 + m(7) * (t34 * t39 + t43 * t46 + t44 * t45) + m(6) * (t51 * t62 + t63 * t70 + t64 * t69) + m(5) * (t105 * t87 + t106 * t86 + t74 * t90) + ((t30 / 0.2e1 + t33 / 0.2e1 + t41 / 0.2e1) * t223 + (-t14 / 0.2e1 - t20 / 0.2e1 - t36 / 0.2e1) * t219 + (t13 / 0.2e1 + t19 / 0.2e1 + t35 / 0.2e1) * t217) * t218; m(5) * (t105 * t205 + t106 * t203 - t252 * t90) + m(6) * (t203 * t70 + t205 * t69 - t252 * t62) + m(7) * (t203 * t46 + t205 * t45 - t252 * t39); (t28 + t32 + t40) * t253 + (t13 + t19 + t35) * t206 + (t14 + t20 + t36) * t204 + m(7) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t62 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t105 ^ 2 + t106 ^ 2 + t90 ^ 2); m(6) * t73 + m(7) * t42; t21 * t269 + t31 * t264 + t33 * t265 + t22 * t268 + (t217 * t17 / 0.2e1 - t219 * t18 / 0.2e1) * t218 + m(7) * (t34 * t42 + t43 * t52 + t44 * t53) + m(6) * (t51 * t73 + t63 * t84 + t64 * t85) + t230; m(6) * (t203 * t84 + t205 * t85 - t252 * t73) + m(7) * (t203 * t52 + t205 * t53 - t252 * t42); t20 * t268 + t19 * t269 + t32 * t265 + t17 * t266 + t18 * t267 + t31 * t235 + m(7) * (t39 * t42 + t45 * t53 + t46 * t52) + m(6) * (t62 * t73 + t69 * t85 + t70 * t84) + t231; t186 * t17 - t188 * t18 + t207 * t31 + m(7) * (t42 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t73 ^ 2 + t84 ^ 2 + t85 ^ 2) + t242; m(7) * t71; m(7) * (t34 * t71 + t43 * t82 + t44 * t83) + t230; m(7) * (t203 * t82 + t205 * t83 - t252 * t71); m(7) * (t39 * t71 + t45 * t83 + t46 * t82) + t231; m(7) * (t42 * t71 + t52 * t82 + t53 * t83) + t242; m(7) * (t71 ^ 2 + t82 ^ 2 + t83 ^ 2) + t242;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
