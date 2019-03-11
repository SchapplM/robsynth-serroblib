% Calculate joint inertia matrix for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:34
% EndTime: 2019-03-08 20:58:42
% DurationCPUTime: 3.96s
% Computational Cost: add. (20192->551), mult. (32026->807), div. (0->0), fcn. (40199->14), ass. (0->261)
t223 = sin(pkin(6));
t288 = t223 ^ 2;
t251 = m(6) / 0.2e1 + m(7) / 0.2e1;
t287 = 0.2e1 * t251;
t226 = cos(pkin(6));
t222 = sin(pkin(10));
t232 = cos(qJ(2));
t266 = t232 * t222;
t225 = cos(pkin(10));
t230 = sin(qJ(2));
t269 = t225 * t230;
t206 = t226 * t269 + t266;
t252 = qJ(3) + pkin(11);
t218 = sin(t252);
t243 = cos(t252);
t236 = t223 * t243;
t191 = t206 * t218 + t225 * t236;
t265 = t232 * t225;
t267 = t230 * t222;
t208 = -t226 * t267 + t265;
t193 = t208 * t218 - t222 * t236;
t272 = t223 * t230;
t202 = t218 * t272 - t226 * t243;
t275 = t222 * t223;
t194 = t208 * t243 + t218 * t275;
t207 = t226 * t266 + t269;
t220 = pkin(12) + qJ(6);
t217 = sin(t220);
t219 = cos(t220);
t155 = -t194 * t217 + t207 * t219;
t156 = t194 * t219 + t207 * t217;
t274 = t223 * t225;
t192 = t206 * t243 - t218 * t274;
t205 = -t226 * t265 + t267;
t153 = -t192 * t217 + t205 * t219;
t154 = t192 * t219 + t205 * t217;
t91 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t191;
t93 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t191;
t95 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t191;
t38 = t155 * t93 + t156 * t95 + t193 * t91;
t92 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t193;
t94 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t193;
t96 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t193;
t39 = t155 * t94 + t156 * t96 + t193 * t92;
t203 = t226 * t218 + t230 * t236;
t270 = t223 * t232;
t185 = -t203 * t217 - t219 * t270;
t186 = t203 * t219 - t217 * t270;
t112 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t202;
t113 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t202;
t114 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t202;
t54 = t112 * t193 + t113 * t155 + t114 * t156;
t2 = t191 * t38 + t193 * t39 + t202 * t54;
t286 = t2 / 0.2e1;
t285 = t191 / 0.2e1;
t284 = t193 / 0.2e1;
t283 = t202 / 0.2e1;
t231 = cos(qJ(3));
t282 = pkin(3) * t231;
t224 = cos(pkin(12));
t281 = pkin(5) * t224;
t221 = sin(pkin(12));
t277 = t205 * t221;
t97 = rSges(7,1) * t154 + rSges(7,2) * t153 + rSges(7,3) * t191;
t279 = pkin(5) * t277 + pkin(9) * t191 + t192 * t281 + t97;
t276 = t207 * t221;
t98 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t193;
t278 = pkin(5) * t276 + pkin(9) * t193 + t194 * t281 + t98;
t229 = sin(qJ(3));
t273 = t223 * t229;
t271 = t223 * t231;
t268 = t226 * t229;
t116 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t202;
t250 = t221 * t270;
t263 = -pkin(5) * t250 + pkin(9) * t202 + t203 * t281 + t116;
t248 = t225 * t273;
t135 = -pkin(3) * t248 + qJ(4) * t205 + t206 * t282;
t111 = t207 * t135;
t151 = pkin(4) * t192 + qJ(5) * t191;
t262 = t207 * t151 + t111;
t179 = pkin(3) * t268 + (-qJ(4) * t232 + t230 * t282) * t223;
t261 = t135 * t270 + t205 * t179;
t249 = t222 * t273;
t136 = pkin(3) * t249 + qJ(4) * t207 + t208 * t282;
t188 = pkin(2) * t208 + pkin(8) * t207;
t184 = t226 * t188;
t260 = t226 * t136 + t184;
t133 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t207;
t259 = -t133 - t136;
t187 = pkin(2) * t206 + pkin(8) * t205;
t258 = -t135 - t187;
t152 = pkin(4) * t194 + qJ(5) * t193;
t257 = -t136 - t152;
t166 = t203 * rSges(5,1) - t202 * rSges(5,2) - rSges(5,3) * t270;
t256 = -t166 - t179;
t167 = t203 * pkin(4) + t202 * qJ(5);
t255 = -t167 - t179;
t254 = t187 * t275 + t188 * t274;
t253 = -m(5) - m(6) - m(7);
t161 = -t194 * t221 + t207 * t224;
t162 = t194 * t224 + t276;
t106 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t193;
t247 = -t106 + t257;
t189 = -t203 * t221 - t224 * t270;
t190 = t203 * t224 - t250;
t123 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t202;
t246 = -t123 + t255;
t245 = t226 * t152 + t260;
t244 = -t151 + t258;
t209 = t226 * t231 - t229 * t272;
t210 = t230 * t271 + t268;
t180 = t210 * rSges(4,1) + t209 * rSges(4,2) - rSges(4,3) * t270;
t211 = (pkin(2) * t230 - pkin(8) * t232) * t223;
t242 = (-t180 - t211) * t223;
t241 = t257 - t278;
t240 = t255 - t263;
t239 = t135 * t275 + t136 * t274 + t254;
t238 = t151 * t270 + t205 * t167 + t261;
t237 = (-t211 + t256) * t223;
t235 = (-t211 + t246) * t223;
t234 = t151 * t275 + t152 * t274 + t239;
t233 = (-t211 + t240) * t223;
t204 = t226 * rSges(3,3) + (rSges(3,1) * t230 + rSges(3,2) * t232) * t223;
t201 = Icges(3,5) * t226 + (Icges(3,1) * t230 + Icges(3,4) * t232) * t223;
t200 = Icges(3,6) * t226 + (Icges(3,4) * t230 + Icges(3,2) * t232) * t223;
t199 = Icges(3,3) * t226 + (Icges(3,5) * t230 + Icges(3,6) * t232) * t223;
t198 = t208 * t231 + t249;
t197 = -t208 * t229 + t222 * t271;
t196 = t206 * t231 - t248;
t195 = -t206 * t229 - t225 * t271;
t178 = Icges(4,1) * t210 + Icges(4,4) * t209 - Icges(4,5) * t270;
t177 = Icges(4,4) * t210 + Icges(4,2) * t209 - Icges(4,6) * t270;
t176 = Icges(4,5) * t210 + Icges(4,6) * t209 - Icges(4,3) * t270;
t175 = rSges(3,1) * t208 - rSges(3,2) * t207 + rSges(3,3) * t275;
t174 = rSges(3,1) * t206 - rSges(3,2) * t205 - rSges(3,3) * t274;
t173 = Icges(3,1) * t208 - Icges(3,4) * t207 + Icges(3,5) * t275;
t172 = Icges(3,1) * t206 - Icges(3,4) * t205 - Icges(3,5) * t274;
t171 = Icges(3,4) * t208 - Icges(3,2) * t207 + Icges(3,6) * t275;
t170 = Icges(3,4) * t206 - Icges(3,2) * t205 - Icges(3,6) * t274;
t169 = Icges(3,5) * t208 - Icges(3,6) * t207 + Icges(3,3) * t275;
t168 = Icges(3,5) * t206 - Icges(3,6) * t205 - Icges(3,3) * t274;
t165 = Icges(5,1) * t203 - Icges(5,4) * t202 - Icges(5,5) * t270;
t164 = Icges(5,4) * t203 - Icges(5,2) * t202 - Icges(5,6) * t270;
t163 = Icges(5,5) * t203 - Icges(5,6) * t202 - Icges(5,3) * t270;
t160 = t192 * t224 + t277;
t159 = -t192 * t221 + t205 * t224;
t150 = -t174 * t226 - t204 * t274;
t149 = t175 * t226 - t204 * t275;
t144 = rSges(4,1) * t198 + rSges(4,2) * t197 + rSges(4,3) * t207;
t143 = rSges(4,1) * t196 + rSges(4,2) * t195 + rSges(4,3) * t205;
t142 = Icges(4,1) * t198 + Icges(4,4) * t197 + Icges(4,5) * t207;
t141 = Icges(4,1) * t196 + Icges(4,4) * t195 + Icges(4,5) * t205;
t140 = Icges(4,4) * t198 + Icges(4,2) * t197 + Icges(4,6) * t207;
t139 = Icges(4,4) * t196 + Icges(4,2) * t195 + Icges(4,6) * t205;
t138 = Icges(4,5) * t198 + Icges(4,6) * t197 + Icges(4,3) * t207;
t137 = Icges(4,5) * t196 + Icges(4,6) * t195 + Icges(4,3) * t205;
t132 = rSges(5,1) * t192 - rSges(5,2) * t191 + rSges(5,3) * t205;
t131 = Icges(5,1) * t194 - Icges(5,4) * t193 + Icges(5,5) * t207;
t130 = Icges(5,1) * t192 - Icges(5,4) * t191 + Icges(5,5) * t205;
t129 = Icges(5,4) * t194 - Icges(5,2) * t193 + Icges(5,6) * t207;
t128 = Icges(5,4) * t192 - Icges(5,2) * t191 + Icges(5,6) * t205;
t127 = Icges(5,5) * t194 - Icges(5,6) * t193 + Icges(5,3) * t207;
t126 = Icges(5,5) * t192 - Icges(5,6) * t191 + Icges(5,3) * t205;
t122 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t202;
t121 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t202;
t120 = Icges(6,5) * t190 + Icges(6,6) * t189 + Icges(6,3) * t202;
t115 = (t174 * t222 + t175 * t225) * t223;
t108 = -t144 * t270 - t207 * t180;
t107 = t143 * t270 + t205 * t180;
t105 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t191;
t104 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t193;
t103 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t191;
t102 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t193;
t101 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t191;
t100 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t193;
t99 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t191;
t88 = -t176 * t270 + t209 * t177 + t210 * t178;
t87 = t143 * t207 - t144 * t205;
t86 = (-t143 - t187) * t226 + t225 * t242;
t85 = t144 * t226 + t222 * t242 + t184;
t84 = -t163 * t270 - t202 * t164 + t203 * t165;
t83 = t176 * t207 + t177 * t197 + t178 * t198;
t82 = t176 * t205 + t177 * t195 + t178 * t196;
t81 = t163 * t207 - t164 * t193 + t165 * t194;
t80 = t163 * t205 - t164 * t191 + t165 * t192;
t79 = (t143 * t222 + t144 * t225) * t223 + t254;
t78 = -t138 * t270 + t209 * t140 + t210 * t142;
t77 = -t137 * t270 + t209 * t139 + t210 * t141;
t76 = -t127 * t270 - t202 * t129 + t203 * t131;
t75 = -t126 * t270 - t202 * t128 + t203 * t130;
t74 = t207 * t256 + t259 * t270;
t73 = t132 * t270 + t205 * t166 + t261;
t72 = -t116 * t193 + t202 * t98;
t71 = t116 * t191 - t202 * t97;
t70 = t138 * t207 + t140 * t197 + t142 * t198;
t69 = t137 * t207 + t139 * t197 + t141 * t198;
t68 = t138 * t205 + t140 * t195 + t142 * t196;
t67 = t137 * t205 + t139 * t195 + t141 * t196;
t66 = (-t132 + t258) * t226 + t225 * t237;
t65 = t133 * t226 + t222 * t237 + t260;
t64 = t127 * t207 - t129 * t193 + t131 * t194;
t63 = t126 * t207 - t128 * t193 + t130 * t194;
t62 = t127 * t205 - t129 * t191 + t131 * t192;
t61 = t126 * t205 - t128 * t191 + t130 * t192;
t60 = t120 * t202 + t121 * t189 + t122 * t190;
t59 = -t191 * t98 + t193 * t97;
t58 = t112 * t202 + t113 * t185 + t114 * t186;
t57 = t132 * t207 + t205 * t259 + t111;
t56 = t120 * t193 + t121 * t161 + t122 * t162;
t55 = t120 * t191 + t121 * t159 + t122 * t160;
t53 = t112 * t191 + t113 * t153 + t114 * t154;
t52 = (t132 * t222 + t133 * t225) * t223 + t239;
t51 = t100 * t202 + t102 * t189 + t104 * t190;
t50 = t101 * t189 + t103 * t190 + t202 * t99;
t49 = t207 * t246 + t247 * t270;
t48 = t105 * t270 + t205 * t123 + t238;
t47 = t185 * t94 + t186 * t96 + t202 * t92;
t46 = t185 * t93 + t186 * t95 + t202 * t91;
t45 = (-t105 + t244) * t226 + t225 * t235;
t44 = t106 * t226 + t222 * t235 + t245;
t43 = t100 * t193 + t102 * t161 + t104 * t162;
t42 = t101 * t161 + t103 * t162 + t193 * t99;
t41 = t100 * t191 + t102 * t159 + t104 * t160;
t40 = t101 * t159 + t103 * t160 + t191 * t99;
t37 = t153 * t94 + t154 * t96 + t191 * t92;
t36 = t153 * t93 + t154 * t95 + t191 * t91;
t35 = t105 * t207 + t205 * t247 + t262;
t34 = t226 * t88 + (t222 * t78 - t225 * t77) * t223;
t33 = (t105 * t222 + t106 * t225) * t223 + t234;
t32 = t77 * t205 + t78 * t207 - t270 * t88;
t31 = t226 * t84 + (t222 * t76 - t225 * t75) * t223;
t30 = t226 * t83 + (t222 * t70 - t225 * t69) * t223;
t29 = t226 * t82 + (t222 * t68 - t225 * t67) * t223;
t28 = t207 * t240 + t241 * t270;
t27 = t205 * t263 + t270 * t279 + t238;
t26 = t75 * t205 + t76 * t207 - t270 * t84;
t25 = t69 * t205 + t70 * t207 - t270 * t83;
t24 = t67 * t205 + t68 * t207 - t270 * t82;
t23 = (t244 - t279) * t226 + t225 * t233;
t22 = t222 * t233 + t226 * t278 + t245;
t21 = t226 * t81 + (t222 * t64 - t225 * t63) * t223;
t20 = t226 * t80 + (t222 * t62 - t225 * t61) * t223;
t19 = t63 * t205 + t64 * t207 - t270 * t81;
t18 = t61 * t205 + t62 * t207 - t270 * t80;
t17 = t205 * t241 + t207 * t279 + t262;
t16 = (t222 * t279 + t225 * t278) * t223 + t234;
t15 = t226 * t60 + (t222 * t51 - t225 * t50) * t223;
t14 = t50 * t205 + t51 * t207 - t270 * t60;
t13 = t226 * t58 + (t222 * t47 - t225 * t46) * t223;
t12 = t46 * t205 + t47 * t207 - t270 * t58;
t11 = t191 * t46 + t193 * t47 + t202 * t58;
t10 = t226 * t56 + (t222 * t43 - t225 * t42) * t223;
t9 = t226 * t55 + (t222 * t41 - t225 * t40) * t223;
t8 = t42 * t205 + t43 * t207 - t270 * t56;
t7 = t40 * t205 + t41 * t207 - t270 * t55;
t6 = t226 * t54 + (t222 * t39 - t225 * t38) * t223;
t5 = t226 * t53 + (t222 * t37 - t225 * t36) * t223;
t4 = t38 * t205 + t39 * t207 - t270 * t54;
t3 = t36 * t205 + t37 * t207 - t270 * t53;
t1 = t191 * t36 + t193 * t37 + t202 * t53;
t89 = [m(2) + m(3) + m(4) - t253; m(3) * t115 + m(4) * t79 + m(5) * t52 + m(6) * t33 + m(7) * t16; m(7) * (t16 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(6) * (t33 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t52 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t79 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(3) * (t115 ^ 2 + t149 ^ 2 + t150 ^ 2) + (t6 + t30 + t21 + t10 + (t169 * t275 - t171 * t207 + t173 * t208) * t275) * t275 + (-t5 - t9 - t29 - t20 + (-t168 * t274 - t170 * t205 + t172 * t206) * t274 + (-t168 * t275 + t169 * t274 + t170 * t207 + t171 * t205 - t172 * t208 - t173 * t206) * t275) * t274 + (-(-t199 * t274 - t200 * t205 + t201 * t206) * t274 + (t199 * t275 - t207 * t200 + t208 * t201) * t275 + t13 + t15 + t34 + t31 + ((t171 * t232 + t173 * t230) * t222 - (t170 * t232 + t172 * t230) * t225) * t288 + ((-t168 * t225 + t169 * t222 + t200 * t232 + t201 * t230) * t223 + t226 * t199) * t226) * t226; m(4) * t87 + m(5) * t57 + m(6) * t35 + m(7) * t17; (t12 / 0.2e1 + t14 / 0.2e1 + t26 / 0.2e1 + t32 / 0.2e1) * t226 + (t6 / 0.2e1 + t10 / 0.2e1 + t21 / 0.2e1 + t30 / 0.2e1) * t207 + (t5 / 0.2e1 + t9 / 0.2e1 + t20 / 0.2e1 + t29 / 0.2e1) * t205 + m(7) * (t16 * t17 + t22 * t28 + t23 * t27) + m(6) * (t33 * t35 + t44 * t49 + t45 * t48) + m(5) * (t52 * t57 + t65 * t74 + t66 * t73) + m(4) * (t107 * t86 + t108 * t85 + t79 * t87) + ((-t13 / 0.2e1 - t15 / 0.2e1 - t31 / 0.2e1 - t34 / 0.2e1) * t232 + (-t3 / 0.2e1 - t7 / 0.2e1 - t18 / 0.2e1 - t24 / 0.2e1) * t225 + (t4 / 0.2e1 + t8 / 0.2e1 + t19 / 0.2e1 + t25 / 0.2e1) * t222) * t223; (-t12 - t14 - t26 - t32) * t270 + (t4 + t19 + t25 + t8) * t207 + (t3 + t24 + t18 + t7) * t205 + m(7) * (t17 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t35 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t57 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t107 ^ 2 + t108 ^ 2 + t87 ^ 2); t253 * t270; m(7) * (-t16 * t270 + t205 * t22 + t207 * t23) + m(6) * (t205 * t44 + t207 * t45 - t270 * t33) + m(5) * (t205 * t65 + t207 * t66 - t270 * t52); m(7) * (-t17 * t270 + t205 * t28 + t207 * t27) + m(6) * (t205 * t49 + t207 * t48 - t270 * t35) + m(5) * (t205 * t74 + t207 * t73 - t270 * t57); 0.2e1 * (m(5) / 0.2e1 + t251) * (t288 * t232 ^ 2 + t205 ^ 2 + t207 ^ 2); t202 * t287; m(7) * (t16 * t202 + t191 * t22 + t193 * t23) + m(6) * (t191 * t44 + t193 * t45 + t202 * t33); m(7) * (t17 * t202 + t191 * t28 + t193 * t27) + m(6) * (t191 * t49 + t193 * t48 + t202 * t35); (t191 * t205 + t193 * t207 - t202 * t270) * t287; (t191 ^ 2 + t193 ^ 2 + t202 ^ 2) * t287; m(7) * t59; m(7) * (t16 * t59 + t22 * t72 + t23 * t71) + t226 * t11 / 0.2e1 + t6 * t284 + t13 * t283 + t5 * t285 + (t222 * t286 - t225 * t1 / 0.2e1) * t223; m(7) * (t17 * t59 + t27 * t71 + t28 * t72) + t3 * t285 + t4 * t284 + t207 * t286 + t12 * t283 + t205 * t1 / 0.2e1 - t11 * t270 / 0.2e1; m(7) * (t72 * t205 + t71 * t207 - t270 * t59); m(7) * (t191 * t72 + t193 * t71 + t202 * t59); t193 * t2 + t191 * t1 + t202 * t11 + m(7) * (t59 ^ 2 + t71 ^ 2 + t72 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t89(1) t89(2) t89(4) t89(7) t89(11) t89(16); t89(2) t89(3) t89(5) t89(8) t89(12) t89(17); t89(4) t89(5) t89(6) t89(9) t89(13) t89(18); t89(7) t89(8) t89(9) t89(10) t89(14) t89(19); t89(11) t89(12) t89(13) t89(14) t89(15) t89(20); t89(16) t89(17) t89(18) t89(19) t89(20) t89(21);];
Mq  = res;
