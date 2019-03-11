% Calculate joint inertia matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:26
% EndTime: 2019-03-08 21:28:38
% DurationCPUTime: 4.86s
% Computational Cost: add. (22059->559), mult. (38283->808), div. (0->0), fcn. (48492->12), ass. (0->254)
t281 = rSges(7,1) + pkin(5);
t280 = rSges(7,3) + qJ(6);
t218 = sin(pkin(6));
t279 = t218 ^ 2;
t278 = cos(qJ(5));
t225 = cos(qJ(3));
t277 = pkin(3) * t225;
t217 = sin(pkin(10));
t275 = t217 * t218;
t219 = cos(pkin(10));
t274 = t218 * t219;
t223 = sin(qJ(3));
t273 = t218 * t223;
t224 = sin(qJ(2));
t272 = t218 * t224;
t271 = t218 * t225;
t226 = cos(qJ(2));
t270 = t218 * t226;
t220 = cos(pkin(6));
t269 = t220 * t223;
t268 = t220 * t224;
t267 = t220 * t226;
t206 = t217 * t226 + t219 * t268;
t253 = qJ(3) + pkin(11);
t216 = sin(t253);
t237 = cos(t253);
t188 = t206 * t237 - t216 * t274;
t205 = t217 * t224 - t219 * t267;
t222 = sin(qJ(5));
t159 = t188 * t222 - t205 * t278;
t160 = t188 * t278 + t205 * t222;
t230 = t218 * t237;
t187 = t206 * t216 + t219 * t230;
t266 = rSges(7,2) * t187 + t280 * t159 + t160 * t281;
t208 = -t217 * t268 + t219 * t226;
t190 = t208 * t237 + t216 * t275;
t207 = t217 * t267 + t219 * t224;
t161 = t190 * t222 - t207 * t278;
t162 = t190 * t278 + t207 * t222;
t189 = t208 * t216 - t217 * t230;
t265 = rSges(7,2) * t189 + t280 * t161 + t162 * t281;
t242 = t219 * t273;
t137 = -pkin(3) * t242 + qJ(4) * t205 + t206 * t277;
t115 = t207 * t137;
t154 = pkin(4) * t188 + pkin(9) * t187;
t264 = t207 * t154 + t115;
t179 = pkin(3) * t269 + (-qJ(4) * t226 + t224 * t277) * t218;
t263 = t137 * t270 + t205 * t179;
t201 = t220 * t216 + t224 * t230;
t191 = t201 * t222 + t270 * t278;
t192 = t201 * t278 - t222 * t270;
t200 = t216 * t272 - t220 * t237;
t262 = rSges(7,2) * t200 + t280 * t191 + t192 * t281;
t243 = t217 * t273;
t138 = pkin(3) * t243 + qJ(4) * t207 + t208 * t277;
t186 = pkin(2) * t208 + pkin(8) * t207;
t184 = t220 * t186;
t261 = t220 * t138 + t184;
t136 = rSges(5,1) * t190 - rSges(5,2) * t189 + rSges(5,3) * t207;
t260 = -t136 - t138;
t185 = pkin(2) * t206 + pkin(8) * t205;
t259 = -t137 - t185;
t155 = pkin(4) * t190 + pkin(9) * t189;
t258 = -t138 - t155;
t166 = t201 * rSges(5,1) - t200 * rSges(5,2) - rSges(5,3) * t270;
t257 = -t166 - t179;
t173 = pkin(4) * t201 + pkin(9) * t200;
t256 = -t173 - t179;
t255 = t185 * t275 + t186 * t274;
t254 = -m(5) - m(6) - m(7);
t103 = Icges(7,1) * t160 + Icges(7,4) * t187 + Icges(7,5) * t159;
t95 = Icges(7,5) * t160 + Icges(7,6) * t187 + Icges(7,3) * t159;
t99 = Icges(7,4) * t160 + Icges(7,2) * t187 + Icges(7,6) * t159;
t40 = t103 * t160 + t159 * t95 + t187 * t99;
t100 = Icges(7,4) * t162 + Icges(7,2) * t189 + Icges(7,6) * t161;
t104 = Icges(7,1) * t162 + Icges(7,4) * t189 + Icges(7,5) * t161;
t96 = Icges(7,5) * t162 + Icges(7,6) * t189 + Icges(7,3) * t161;
t41 = t100 * t187 + t104 * t160 + t159 * t96;
t120 = Icges(7,5) * t192 + Icges(7,6) * t200 + Icges(7,3) * t191;
t122 = Icges(7,4) * t192 + Icges(7,2) * t200 + Icges(7,6) * t191;
t124 = Icges(7,1) * t192 + Icges(7,4) * t200 + Icges(7,5) * t191;
t59 = t120 * t159 + t122 * t187 + t124 * t160;
t1 = t187 * t40 + t189 * t41 + t200 * t59;
t101 = Icges(6,4) * t160 - Icges(6,2) * t159 + Icges(6,6) * t187;
t105 = Icges(6,1) * t160 - Icges(6,4) * t159 + Icges(6,5) * t187;
t97 = Icges(6,5) * t160 - Icges(6,6) * t159 + Icges(6,3) * t187;
t42 = -t101 * t159 + t105 * t160 + t187 * t97;
t102 = Icges(6,4) * t162 - Icges(6,2) * t161 + Icges(6,6) * t189;
t106 = Icges(6,1) * t162 - Icges(6,4) * t161 + Icges(6,5) * t189;
t98 = Icges(6,5) * t162 - Icges(6,6) * t161 + Icges(6,3) * t189;
t43 = -t102 * t159 + t106 * t160 + t187 * t98;
t121 = Icges(6,5) * t192 - Icges(6,6) * t191 + Icges(6,3) * t200;
t123 = Icges(6,4) * t192 - Icges(6,2) * t191 + Icges(6,6) * t200;
t125 = Icges(6,1) * t192 - Icges(6,4) * t191 + Icges(6,5) * t200;
t60 = t121 * t187 - t123 * t159 + t125 * t160;
t2 = t187 * t42 + t189 * t43 + t200 * t60;
t252 = t2 / 0.2e1 + t1 / 0.2e1;
t44 = t103 * t162 + t161 * t95 + t189 * t99;
t45 = t100 * t189 + t104 * t162 + t161 * t96;
t61 = t120 * t161 + t122 * t189 + t124 * t162;
t3 = t187 * t44 + t189 * t45 + t200 * t61;
t46 = -t101 * t161 + t105 * t162 + t189 * t97;
t47 = -t102 * t161 + t106 * t162 + t189 * t98;
t62 = t121 * t189 - t123 * t161 + t125 * t162;
t4 = t187 * t46 + t189 * t47 + t200 * t62;
t251 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t40 * t205 + t41 * t207 - t270 * t59;
t6 = t42 * t205 + t43 * t207 - t270 * t60;
t250 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t44 * t205 + t45 * t207 - t270 * t61;
t8 = t46 * t205 + t47 * t207 - t270 * t62;
t249 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t220 * t60 + (t217 * t43 - t219 * t42) * t218;
t9 = t220 * t59 + (t217 * t41 - t219 * t40) * t218;
t248 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t220 * t61 + (t217 * t45 - t219 * t44) * t218;
t12 = t220 * t62 + (t217 * t47 - t219 * t46) * t218;
t247 = t12 / 0.2e1 + t11 / 0.2e1;
t52 = t103 * t192 + t191 * t95 + t200 * t99;
t53 = t100 * t200 + t104 * t192 + t191 * t96;
t65 = t120 * t191 + t122 * t200 + t124 * t192;
t13 = t187 * t52 + t189 * t53 + t200 * t65;
t54 = -t101 * t191 + t105 * t192 + t200 * t97;
t55 = -t102 * t191 + t106 * t192 + t200 * t98;
t66 = t121 * t200 - t123 * t191 + t125 * t192;
t14 = t187 * t54 + t189 * t55 + t200 * t66;
t246 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t52 * t205 + t53 * t207 - t270 * t65;
t16 = t54 * t205 + t55 * t207 - t270 * t66;
t245 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t220 * t65 + (t217 * t53 - t219 * t52) * t218;
t18 = t220 * t66 + (t217 * t55 - t219 * t54) * t218;
t244 = t18 / 0.2e1 + t17 / 0.2e1;
t110 = rSges(6,1) * t162 - rSges(6,2) * t161 + rSges(6,3) * t189;
t241 = -t110 + t258;
t127 = rSges(6,1) * t192 - rSges(6,2) * t191 + rSges(6,3) * t200;
t240 = -t127 + t256;
t239 = t220 * t155 + t261;
t238 = -t154 + t259;
t209 = t220 * t225 - t223 * t272;
t210 = t224 * t271 + t269;
t180 = t210 * rSges(4,1) + t209 * rSges(4,2) - rSges(4,3) * t270;
t211 = (pkin(2) * t224 - pkin(8) * t226) * t218;
t236 = (-t180 - t211) * t218;
t235 = t258 - t265;
t234 = t137 * t275 + t138 * t274 + t255;
t233 = t154 * t270 + t205 * t173 + t263;
t232 = t256 - t262;
t231 = (-t211 + t257) * t218;
t229 = (-t211 + t240) * t218;
t228 = t154 * t275 + t155 * t274 + t234;
t227 = (-t211 + t232) * t218;
t202 = t220 * rSges(3,3) + (rSges(3,1) * t224 + rSges(3,2) * t226) * t218;
t199 = Icges(3,5) * t220 + (Icges(3,1) * t224 + Icges(3,4) * t226) * t218;
t198 = Icges(3,6) * t220 + (Icges(3,4) * t224 + Icges(3,2) * t226) * t218;
t197 = Icges(3,3) * t220 + (Icges(3,5) * t224 + Icges(3,6) * t226) * t218;
t196 = t208 * t225 + t243;
t195 = -t208 * t223 + t217 * t271;
t194 = t206 * t225 - t242;
t193 = -t206 * t223 - t219 * t271;
t178 = Icges(4,1) * t210 + Icges(4,4) * t209 - Icges(4,5) * t270;
t177 = Icges(4,4) * t210 + Icges(4,2) * t209 - Icges(4,6) * t270;
t176 = Icges(4,5) * t210 + Icges(4,6) * t209 - Icges(4,3) * t270;
t175 = rSges(3,1) * t208 - rSges(3,2) * t207 + rSges(3,3) * t275;
t174 = rSges(3,1) * t206 - rSges(3,2) * t205 - rSges(3,3) * t274;
t172 = Icges(3,1) * t208 - Icges(3,4) * t207 + Icges(3,5) * t275;
t171 = Icges(3,1) * t206 - Icges(3,4) * t205 - Icges(3,5) * t274;
t170 = Icges(3,4) * t208 - Icges(3,2) * t207 + Icges(3,6) * t275;
t169 = Icges(3,4) * t206 - Icges(3,2) * t205 - Icges(3,6) * t274;
t168 = Icges(3,5) * t208 - Icges(3,6) * t207 + Icges(3,3) * t275;
t167 = Icges(3,5) * t206 - Icges(3,6) * t205 - Icges(3,3) * t274;
t165 = Icges(5,1) * t201 - Icges(5,4) * t200 - Icges(5,5) * t270;
t164 = Icges(5,4) * t201 - Icges(5,2) * t200 - Icges(5,6) * t270;
t163 = Icges(5,5) * t201 - Icges(5,6) * t200 - Icges(5,3) * t270;
t153 = -t174 * t220 - t202 * t274;
t152 = t175 * t220 - t202 * t275;
t147 = rSges(4,1) * t196 + rSges(4,2) * t195 + rSges(4,3) * t207;
t146 = rSges(4,1) * t194 + rSges(4,2) * t193 + rSges(4,3) * t205;
t145 = Icges(4,1) * t196 + Icges(4,4) * t195 + Icges(4,5) * t207;
t144 = Icges(4,1) * t194 + Icges(4,4) * t193 + Icges(4,5) * t205;
t143 = Icges(4,4) * t196 + Icges(4,2) * t195 + Icges(4,6) * t207;
t142 = Icges(4,4) * t194 + Icges(4,2) * t193 + Icges(4,6) * t205;
t141 = Icges(4,5) * t196 + Icges(4,6) * t195 + Icges(4,3) * t207;
t140 = Icges(4,5) * t194 + Icges(4,6) * t193 + Icges(4,3) * t205;
t135 = rSges(5,1) * t188 - rSges(5,2) * t187 + rSges(5,3) * t205;
t134 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t207;
t133 = Icges(5,1) * t188 - Icges(5,4) * t187 + Icges(5,5) * t205;
t132 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t207;
t131 = Icges(5,4) * t188 - Icges(5,2) * t187 + Icges(5,6) * t205;
t130 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t207;
t129 = Icges(5,5) * t188 - Icges(5,6) * t187 + Icges(5,3) * t205;
t116 = (t174 * t217 + t175 * t219) * t218;
t112 = -t147 * t270 - t207 * t180;
t111 = t146 * t270 + t205 * t180;
t108 = rSges(6,1) * t160 - rSges(6,2) * t159 + rSges(6,3) * t187;
t94 = -t176 * t270 + t209 * t177 + t210 * t178;
t93 = t146 * t207 - t147 * t205;
t92 = (-t146 - t185) * t220 + t219 * t236;
t91 = t147 * t220 + t217 * t236 + t184;
t90 = -t163 * t270 - t200 * t164 + t201 * t165;
t89 = t176 * t207 + t177 * t195 + t178 * t196;
t88 = t176 * t205 + t177 * t193 + t178 * t194;
t87 = t163 * t207 - t164 * t189 + t165 * t190;
t86 = t163 * t205 - t164 * t187 + t165 * t188;
t85 = (t146 * t217 + t147 * t219) * t218 + t255;
t84 = -t141 * t270 + t209 * t143 + t210 * t145;
t83 = -t140 * t270 + t209 * t142 + t210 * t144;
t82 = t110 * t200 - t127 * t189;
t81 = -t108 * t200 + t127 * t187;
t80 = -t130 * t270 - t200 * t132 + t201 * t134;
t79 = -t129 * t270 - t200 * t131 + t201 * t133;
t78 = t207 * t257 + t260 * t270;
t77 = t135 * t270 + t205 * t166 + t263;
t76 = t141 * t207 + t143 * t195 + t145 * t196;
t75 = t140 * t207 + t142 * t195 + t144 * t196;
t74 = t141 * t205 + t143 * t193 + t145 * t194;
t73 = t140 * t205 + t142 * t193 + t144 * t194;
t72 = (-t135 + t259) * t220 + t219 * t231;
t71 = t136 * t220 + t217 * t231 + t261;
t70 = t130 * t207 - t132 * t189 + t134 * t190;
t69 = t129 * t207 - t131 * t189 + t133 * t190;
t68 = t130 * t205 - t132 * t187 + t134 * t188;
t67 = t129 * t205 - t131 * t187 + t133 * t188;
t64 = t108 * t189 - t110 * t187;
t63 = t135 * t207 + t205 * t260 + t115;
t58 = (t135 * t217 + t136 * t219) * t218 + t234;
t57 = -t189 * t262 + t200 * t265;
t56 = t187 * t262 - t200 * t266;
t51 = t207 * t240 + t241 * t270;
t50 = t108 * t270 + t205 * t127 + t233;
t49 = (-t108 + t238) * t220 + t219 * t229;
t48 = t110 * t220 + t217 * t229 + t239;
t39 = -t187 * t265 + t189 * t266;
t38 = t108 * t207 + t205 * t241 + t264;
t37 = (t108 * t217 + t110 * t219) * t218 + t228;
t36 = t220 * t94 + (t217 * t84 - t219 * t83) * t218;
t35 = t207 * t232 + t235 * t270;
t34 = t205 * t262 + t266 * t270 + t233;
t33 = t83 * t205 + t84 * t207 - t270 * t94;
t32 = (t238 - t266) * t220 + t219 * t227;
t31 = t217 * t227 + t220 * t265 + t239;
t30 = t220 * t90 + (t217 * t80 - t219 * t79) * t218;
t29 = t220 * t89 + (t217 * t76 - t219 * t75) * t218;
t28 = t220 * t88 + (t217 * t74 - t219 * t73) * t218;
t27 = t79 * t205 + t80 * t207 - t270 * t90;
t26 = t75 * t205 + t76 * t207 - t270 * t89;
t25 = t73 * t205 + t74 * t207 - t270 * t88;
t24 = t220 * t87 + (t217 * t70 - t219 * t69) * t218;
t23 = t220 * t86 + (t217 * t68 - t219 * t67) * t218;
t22 = t69 * t205 + t70 * t207 - t270 * t87;
t21 = t67 * t205 + t68 * t207 - t270 * t86;
t20 = t205 * t235 + t207 * t266 + t264;
t19 = (t217 * t266 + t219 * t265) * t218 + t228;
t107 = [m(2) + m(3) + m(4) - t254; m(3) * t116 + m(4) * t85 + m(5) * t58 + m(6) * t37 + m(7) * t19; m(7) * (t19 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t37 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t85 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(3) * (t116 ^ 2 + t152 ^ 2 + t153 ^ 2) + (t12 + t11 + t29 + t24 + (t168 * t275 - t170 * t207 + t172 * t208) * t275) * t275 + (-t9 - t10 - t23 - t28 + (-t167 * t274 - t169 * t205 + t171 * t206) * t274 + (-t167 * t275 + t168 * t274 + t169 * t207 + t170 * t205 - t171 * t208 - t172 * t206) * t275) * t274 + (-(-t197 * t274 - t198 * t205 + t199 * t206) * t274 + (t197 * t275 - t198 * t207 + t199 * t208) * t275 + t18 + t17 + t36 + t30 + ((t170 * t226 + t172 * t224) * t217 - (t169 * t226 + t171 * t224) * t219) * t279 + ((-t167 * t219 + t168 * t217 + t198 * t226 + t199 * t224) * t218 + t220 * t197) * t220) * t220; m(4) * t93 + m(5) * t63 + m(6) * t38 + m(7) * t20; (t27 / 0.2e1 + t33 / 0.2e1 + t245) * t220 + (t24 / 0.2e1 + t29 / 0.2e1 + t247) * t207 + (t23 / 0.2e1 + t28 / 0.2e1 + t248) * t205 + m(7) * (t19 * t20 + t31 * t35 + t32 * t34) + m(6) * (t37 * t38 + t48 * t51 + t49 * t50) + m(5) * (t58 * t63 + t71 * t78 + t72 * t77) + m(4) * (t111 * t92 + t112 * t91 + t85 * t93) + ((-t30 / 0.2e1 - t36 / 0.2e1 - t244) * t226 + (-t21 / 0.2e1 - t25 / 0.2e1 - t250) * t219 + (t22 / 0.2e1 + t26 / 0.2e1 + t249) * t217) * t218; (-t15 - t16 - t27 - t33) * t270 + (t8 + t7 + t26 + t22) * t207 + (t6 + t5 + t25 + t21) * t205 + m(7) * (t20 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t38 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t63 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2 + t93 ^ 2); t254 * t270; m(7) * (-t19 * t270 + t205 * t31 + t207 * t32) + m(6) * (t205 * t48 + t207 * t49 - t270 * t37) + m(5) * (t205 * t71 + t207 * t72 - t270 * t58); m(7) * (-t20 * t270 + t205 * t35 + t207 * t34) + m(6) * (t205 * t51 + t207 * t50 - t270 * t38) + m(5) * (t205 * t78 + t207 * t77 - t270 * t63); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t226 ^ 2 * t279 + t205 ^ 2 + t207 ^ 2); m(6) * t64 + m(7) * t39; t246 * t220 + t244 * t200 + t247 * t189 + t248 * t187 + m(7) * (t19 * t39 + t31 * t57 + t32 * t56) + m(6) * (t37 * t64 + t48 * t82 + t49 * t81) + (t217 * t251 - t219 * t252) * t218; -t246 * t270 + t251 * t207 + t252 * t205 + t245 * t200 + t249 * t189 + t250 * t187 + m(7) * (t20 * t39 + t34 * t56 + t35 * t57) + m(6) * (t38 * t64 + t50 * t81 + t51 * t82); m(6) * (t82 * t205 + t81 * t207 - t270 * t64) + m(7) * (t57 * t205 + t56 * t207 - t270 * t39); (t14 + t13) * t200 + (t3 + t4) * t189 + (t1 + t2) * t187 + m(7) * (t39 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t64 ^ 2 + t81 ^ 2 + t82 ^ 2); m(7) * t191; m(7) * (t159 * t31 + t161 * t32 + t19 * t191); m(7) * (t159 * t35 + t161 * t34 + t191 * t20); m(7) * (t159 * t205 + t161 * t207 - t191 * t270); m(7) * (t159 * t57 + t161 * t56 + t191 * t39); m(7) * (t159 ^ 2 + t161 ^ 2 + t191 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t107(1) t107(2) t107(4) t107(7) t107(11) t107(16); t107(2) t107(3) t107(5) t107(8) t107(12) t107(17); t107(4) t107(5) t107(6) t107(9) t107(13) t107(18); t107(7) t107(8) t107(9) t107(10) t107(14) t107(19); t107(11) t107(12) t107(13) t107(14) t107(15) t107(20); t107(16) t107(17) t107(18) t107(19) t107(20) t107(21);];
Mq  = res;
