% Calculate joint inertia matrix for
% S6PRRPRP1
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:41
% EndTime: 2019-03-08 21:22:52
% DurationCPUTime: 4.88s
% Computational Cost: add. (22963->564), mult. (39055->813), div. (0->0), fcn. (49304->12), ass. (0->259)
t287 = rSges(7,3) + qJ(6);
t217 = sin(pkin(6));
t286 = t217 ^ 2;
t226 = cos(qJ(3));
t285 = pkin(3) * t226;
t225 = cos(qJ(5));
t284 = pkin(5) * t225;
t219 = cos(pkin(6));
t216 = sin(pkin(10));
t227 = cos(qJ(2));
t268 = t227 * t216;
t218 = cos(pkin(10));
t224 = sin(qJ(2));
t271 = t218 * t224;
t204 = t219 * t271 + t268;
t255 = qJ(3) + pkin(11);
t215 = sin(t255);
t238 = cos(t255);
t276 = t217 * t218;
t188 = t204 * t238 - t215 * t276;
t267 = t227 * t218;
t269 = t224 * t216;
t203 = -t219 * t267 + t269;
t222 = sin(qJ(5));
t159 = -t188 * t222 + t203 * t225;
t279 = t203 * t222;
t160 = t188 * t225 + t279;
t231 = t217 * t238;
t187 = t204 * t215 + t218 * t231;
t281 = rSges(7,1) * t160 + rSges(7,2) * t159 + pkin(5) * t279 + t287 * t187 + t188 * t284;
t206 = -t219 * t269 + t267;
t277 = t216 * t217;
t190 = t206 * t238 + t215 * t277;
t205 = t219 * t268 + t271;
t161 = -t190 * t222 + t205 * t225;
t278 = t205 * t222;
t162 = t190 * t225 + t278;
t189 = t206 * t215 - t216 * t231;
t280 = rSges(7,1) * t162 + rSges(7,2) * t161 + pkin(5) * t278 + t287 * t189 + t190 * t284;
t223 = sin(qJ(3));
t275 = t217 * t223;
t274 = t217 * t224;
t273 = t217 * t226;
t272 = t217 * t227;
t270 = t219 * t223;
t244 = t218 * t275;
t138 = -pkin(3) * t244 + qJ(4) * t203 + t204 * t285;
t115 = t205 * t138;
t155 = pkin(4) * t188 + pkin(9) * t187;
t266 = t205 * t155 + t115;
t201 = t219 * t215 + t224 * t231;
t191 = -t201 * t222 - t225 * t272;
t243 = t222 * t272;
t192 = t201 * t225 - t243;
t200 = t215 * t274 - t219 * t238;
t265 = rSges(7,1) * t192 + rSges(7,2) * t191 - pkin(5) * t243 + t287 * t200 + t201 * t284;
t179 = pkin(3) * t270 + (-qJ(4) * t227 + t224 * t285) * t217;
t264 = t138 * t272 + t203 * t179;
t245 = t216 * t275;
t139 = pkin(3) * t245 + qJ(4) * t205 + t206 * t285;
t186 = pkin(2) * t206 + pkin(8) * t205;
t184 = t219 * t186;
t263 = t219 * t139 + t184;
t137 = rSges(5,1) * t190 - rSges(5,2) * t189 + rSges(5,3) * t205;
t262 = -t137 - t139;
t185 = pkin(2) * t204 + pkin(8) * t203;
t261 = -t138 - t185;
t156 = pkin(4) * t190 + pkin(9) * t189;
t260 = -t139 - t156;
t166 = t201 * rSges(5,1) - t200 * rSges(5,2) - rSges(5,3) * t272;
t259 = -t166 - t179;
t173 = t201 * pkin(4) + t200 * pkin(9);
t258 = -t173 - t179;
t257 = t185 * t277 + t186 * t276;
t256 = -m(5) - m(6) - m(7);
t101 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t187;
t105 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t187;
t97 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t187;
t42 = t101 * t159 + t105 * t160 + t187 * t97;
t102 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t189;
t106 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t189;
t98 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t189;
t43 = t102 * t159 + t106 * t160 + t187 * t98;
t121 = Icges(7,5) * t192 + Icges(7,6) * t191 + Icges(7,3) * t200;
t123 = Icges(7,4) * t192 + Icges(7,2) * t191 + Icges(7,6) * t200;
t125 = Icges(7,1) * t192 + Icges(7,4) * t191 + Icges(7,5) * t200;
t59 = t121 * t187 + t123 * t159 + t125 * t160;
t1 = t187 * t42 + t189 * t43 + t200 * t59;
t103 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t187;
t107 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t187;
t99 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t187;
t44 = t103 * t159 + t107 * t160 + t187 * t99;
t100 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t189;
t104 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t189;
t108 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t189;
t45 = t100 * t187 + t104 * t159 + t108 * t160;
t122 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t200;
t124 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t200;
t126 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t200;
t60 = t122 * t187 + t124 * t159 + t126 * t160;
t2 = t187 * t44 + t189 * t45 + t200 * t60;
t254 = t1 / 0.2e1 + t2 / 0.2e1;
t46 = t101 * t161 + t105 * t162 + t189 * t97;
t47 = t102 * t161 + t106 * t162 + t189 * t98;
t61 = t121 * t189 + t123 * t161 + t125 * t162;
t3 = t187 * t46 + t189 * t47 + t200 * t61;
t48 = t103 * t161 + t107 * t162 + t189 * t99;
t49 = t100 * t189 + t104 * t161 + t108 * t162;
t62 = t122 * t189 + t124 * t161 + t126 * t162;
t4 = t187 * t48 + t189 * t49 + t200 * t62;
t253 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t42 * t203 + t43 * t205 - t272 * t59;
t6 = t44 * t203 + t45 * t205 - t272 * t60;
t252 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t46 * t203 + t47 * t205 - t272 * t61;
t8 = t48 * t203 + t49 * t205 - t272 * t62;
t251 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t219 * t60 + (t216 * t45 - t218 * t44) * t217;
t9 = t219 * t59 + (t216 * t43 - t218 * t42) * t217;
t250 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t219 * t61 + (t216 * t47 - t218 * t46) * t217;
t12 = t219 * t62 + (t216 * t49 - t218 * t48) * t217;
t249 = t12 / 0.2e1 + t11 / 0.2e1;
t54 = t101 * t191 + t105 * t192 + t200 * t97;
t55 = t102 * t191 + t106 * t192 + t200 * t98;
t65 = t121 * t200 + t123 * t191 + t125 * t192;
t13 = t187 * t54 + t189 * t55 + t200 * t65;
t56 = t103 * t191 + t107 * t192 + t200 * t99;
t57 = t100 * t200 + t104 * t191 + t108 * t192;
t66 = t122 * t200 + t124 * t191 + t126 * t192;
t14 = t187 * t56 + t189 * t57 + t200 * t66;
t248 = -t14 / 0.2e1 - t13 / 0.2e1;
t15 = t54 * t203 + t55 * t205 - t272 * t65;
t16 = t56 * t203 + t57 * t205 - t272 * t66;
t247 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t219 * t65 + (t216 * t55 - t218 * t54) * t217;
t18 = t219 * t66 + (t216 * t57 - t218 * t56) * t217;
t246 = t17 / 0.2e1 + t18 / 0.2e1;
t112 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t189;
t242 = -t112 + t260;
t128 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t200;
t241 = -t128 + t258;
t240 = t219 * t156 + t263;
t239 = -t155 + t261;
t207 = t219 * t226 - t223 * t274;
t208 = t224 * t273 + t270;
t180 = t208 * rSges(4,1) + t207 * rSges(4,2) - rSges(4,3) * t272;
t209 = (pkin(2) * t224 - pkin(8) * t227) * t217;
t237 = (-t180 - t209) * t217;
t236 = t260 - t280;
t235 = t258 - t265;
t234 = t138 * t277 + t139 * t276 + t257;
t233 = t155 * t272 + t203 * t173 + t264;
t232 = (-t209 + t259) * t217;
t230 = (-t209 + t241) * t217;
t229 = t155 * t277 + t156 * t276 + t234;
t228 = (-t209 + t235) * t217;
t202 = t219 * rSges(3,3) + (rSges(3,1) * t224 + rSges(3,2) * t227) * t217;
t199 = Icges(3,5) * t219 + (Icges(3,1) * t224 + Icges(3,4) * t227) * t217;
t198 = Icges(3,6) * t219 + (Icges(3,4) * t224 + Icges(3,2) * t227) * t217;
t197 = Icges(3,3) * t219 + (Icges(3,5) * t224 + Icges(3,6) * t227) * t217;
t196 = t206 * t226 + t245;
t195 = -t206 * t223 + t216 * t273;
t194 = t204 * t226 - t244;
t193 = -t204 * t223 - t218 * t273;
t178 = Icges(4,1) * t208 + Icges(4,4) * t207 - Icges(4,5) * t272;
t177 = Icges(4,4) * t208 + Icges(4,2) * t207 - Icges(4,6) * t272;
t176 = Icges(4,5) * t208 + Icges(4,6) * t207 - Icges(4,3) * t272;
t175 = rSges(3,1) * t206 - rSges(3,2) * t205 + rSges(3,3) * t277;
t174 = rSges(3,1) * t204 - rSges(3,2) * t203 - rSges(3,3) * t276;
t172 = Icges(3,1) * t206 - Icges(3,4) * t205 + Icges(3,5) * t277;
t171 = Icges(3,1) * t204 - Icges(3,4) * t203 - Icges(3,5) * t276;
t170 = Icges(3,4) * t206 - Icges(3,2) * t205 + Icges(3,6) * t277;
t169 = Icges(3,4) * t204 - Icges(3,2) * t203 - Icges(3,6) * t276;
t168 = Icges(3,5) * t206 - Icges(3,6) * t205 + Icges(3,3) * t277;
t167 = Icges(3,5) * t204 - Icges(3,6) * t203 - Icges(3,3) * t276;
t165 = Icges(5,1) * t201 - Icges(5,4) * t200 - Icges(5,5) * t272;
t164 = Icges(5,4) * t201 - Icges(5,2) * t200 - Icges(5,6) * t272;
t163 = Icges(5,5) * t201 - Icges(5,6) * t200 - Icges(5,3) * t272;
t154 = -t174 * t219 - t202 * t276;
t153 = t175 * t219 - t202 * t277;
t148 = rSges(4,1) * t196 + rSges(4,2) * t195 + rSges(4,3) * t205;
t147 = rSges(4,1) * t194 + rSges(4,2) * t193 + rSges(4,3) * t203;
t146 = Icges(4,1) * t196 + Icges(4,4) * t195 + Icges(4,5) * t205;
t145 = Icges(4,1) * t194 + Icges(4,4) * t193 + Icges(4,5) * t203;
t144 = Icges(4,4) * t196 + Icges(4,2) * t195 + Icges(4,6) * t205;
t143 = Icges(4,4) * t194 + Icges(4,2) * t193 + Icges(4,6) * t203;
t142 = Icges(4,5) * t196 + Icges(4,6) * t195 + Icges(4,3) * t205;
t141 = Icges(4,5) * t194 + Icges(4,6) * t193 + Icges(4,3) * t203;
t136 = rSges(5,1) * t188 - rSges(5,2) * t187 + rSges(5,3) * t203;
t135 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t205;
t134 = Icges(5,1) * t188 - Icges(5,4) * t187 + Icges(5,5) * t203;
t133 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t205;
t132 = Icges(5,4) * t188 - Icges(5,2) * t187 + Icges(5,6) * t203;
t131 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t205;
t130 = Icges(5,5) * t188 - Icges(5,6) * t187 + Icges(5,3) * t203;
t117 = (t174 * t216 + t175 * t218) * t217;
t114 = -t148 * t272 - t205 * t180;
t113 = t147 * t272 + t203 * t180;
t110 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t187;
t94 = -t176 * t272 + t207 * t177 + t208 * t178;
t93 = t147 * t205 - t148 * t203;
t92 = (-t147 - t185) * t219 + t218 * t237;
t91 = t148 * t219 + t216 * t237 + t184;
t90 = -t163 * t272 - t200 * t164 + t201 * t165;
t89 = t176 * t205 + t177 * t195 + t178 * t196;
t88 = t176 * t203 + t177 * t193 + t178 * t194;
t87 = t163 * t205 - t164 * t189 + t165 * t190;
t86 = t163 * t203 - t164 * t187 + t165 * t188;
t85 = (t147 * t216 + t148 * t218) * t217 + t257;
t84 = -t142 * t272 + t207 * t144 + t208 * t146;
t83 = -t141 * t272 + t207 * t143 + t208 * t145;
t82 = t112 * t200 - t128 * t189;
t81 = -t110 * t200 + t128 * t187;
t80 = -t131 * t272 - t200 * t133 + t201 * t135;
t79 = -t130 * t272 - t200 * t132 + t201 * t134;
t78 = t205 * t259 + t262 * t272;
t77 = t136 * t272 + t203 * t166 + t264;
t76 = t142 * t205 + t144 * t195 + t146 * t196;
t75 = t141 * t205 + t143 * t195 + t145 * t196;
t74 = t142 * t203 + t144 * t193 + t146 * t194;
t73 = t141 * t203 + t143 * t193 + t145 * t194;
t72 = (-t136 + t261) * t219 + t218 * t232;
t71 = t137 * t219 + t216 * t232 + t263;
t70 = t131 * t205 - t133 * t189 + t135 * t190;
t69 = t130 * t205 - t132 * t189 + t134 * t190;
t68 = t131 * t203 - t133 * t187 + t135 * t188;
t67 = t130 * t203 - t132 * t187 + t134 * t188;
t64 = t110 * t189 - t112 * t187;
t63 = t136 * t205 + t203 * t262 + t115;
t58 = (t136 * t216 + t137 * t218) * t217 + t234;
t53 = t205 * t241 + t242 * t272;
t52 = t110 * t272 + t203 * t128 + t233;
t51 = (-t110 + t239) * t219 + t218 * t230;
t50 = t112 * t219 + t216 * t230 + t240;
t41 = -t189 * t265 + t200 * t280;
t40 = t187 * t265 - t200 * t281;
t39 = t110 * t205 + t203 * t242 + t266;
t38 = (t110 * t216 + t112 * t218) * t217 + t229;
t37 = t219 * t94 + (t216 * t84 - t218 * t83) * t217;
t36 = t83 * t203 + t84 * t205 - t272 * t94;
t35 = -t187 * t280 + t189 * t281;
t34 = t205 * t235 + t236 * t272;
t33 = t265 * t203 + t272 * t281 + t233;
t32 = t219 * t90 + (t216 * t80 - t218 * t79) * t217;
t31 = t219 * t89 + (t216 * t76 - t218 * t75) * t217;
t30 = t219 * t88 + (t216 * t74 - t218 * t73) * t217;
t29 = (t239 - t281) * t219 + t218 * t228;
t28 = t216 * t228 + t219 * t280 + t240;
t27 = t79 * t203 + t80 * t205 - t272 * t90;
t26 = t75 * t203 + t76 * t205 - t272 * t89;
t25 = t73 * t203 + t74 * t205 - t272 * t88;
t24 = t219 * t87 + (t216 * t70 - t218 * t69) * t217;
t23 = t219 * t86 + (t216 * t68 - t218 * t67) * t217;
t22 = t69 * t203 + t70 * t205 - t272 * t87;
t21 = t67 * t203 + t68 * t205 - t272 * t86;
t20 = t203 * t236 + t205 * t281 + t266;
t19 = (t216 * t281 + t218 * t280) * t217 + t229;
t95 = [m(2) + m(3) + m(4) - t256; m(3) * t117 + m(4) * t85 + m(5) * t58 + m(6) * t38 + m(7) * t19; m(7) * (t19 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t38 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t85 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(3) * (t117 ^ 2 + t153 ^ 2 + t154 ^ 2) + (t12 + t11 + t31 + t24 + (t168 * t277 - t170 * t205 + t172 * t206) * t277) * t277 + (-t9 - t10 - t30 - t23 + (-t167 * t276 - t169 * t203 + t171 * t204) * t276 + (-t167 * t277 + t168 * t276 + t169 * t205 + t170 * t203 - t171 * t206 - t172 * t204) * t277) * t276 + (-(-t197 * t276 - t198 * t203 + t199 * t204) * t276 + (t197 * t277 - t198 * t205 + t199 * t206) * t277 + t17 + t18 + t37 + t32 + ((t170 * t227 + t172 * t224) * t216 - (t169 * t227 + t171 * t224) * t218) * t286 + ((-t167 * t218 + t168 * t216 + t198 * t227 + t199 * t224) * t217 + t219 * t197) * t219) * t219; m(4) * t93 + m(5) * t63 + m(6) * t39 + m(7) * t20; (t27 / 0.2e1 + t36 / 0.2e1 + t247) * t219 + (t24 / 0.2e1 + t31 / 0.2e1 + t249) * t205 + (t23 / 0.2e1 + t30 / 0.2e1 + t250) * t203 + m(7) * (t19 * t20 + t28 * t34 + t29 * t33) + m(6) * (t38 * t39 + t50 * t53 + t51 * t52) + m(5) * (t58 * t63 + t71 * t78 + t72 * t77) + m(4) * (t113 * t92 + t114 * t91 + t85 * t93) + ((-t32 / 0.2e1 - t37 / 0.2e1 - t246) * t227 + (-t21 / 0.2e1 - t25 / 0.2e1 - t252) * t218 + (t22 / 0.2e1 + t26 / 0.2e1 + t251) * t216) * t217; (-t15 - t16 - t27 - t36) * t272 + (t8 + t7 + t26 + t22) * t205 + (t5 + t6 + t21 + t25) * t203 + m(7) * (t20 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t39 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t63 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2 + t93 ^ 2); t256 * t272; m(7) * (-t19 * t272 + t203 * t28 + t205 * t29) + m(6) * (t203 * t50 + t205 * t51 - t272 * t38) + m(5) * (t203 * t71 + t205 * t72 - t272 * t58); m(7) * (-t20 * t272 + t203 * t34 + t205 * t33) + m(6) * (t203 * t53 + t205 * t52 - t272 * t39) + m(5) * (t203 * t78 + t205 * t77 - t272 * t63); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t286 * t227 ^ 2 + t203 ^ 2 + t205 ^ 2); m(6) * t64 + m(7) * t35; -t248 * t219 + t246 * t200 + t249 * t189 + t250 * t187 + m(7) * (t19 * t35 + t28 * t41 + t29 * t40) + m(6) * (t38 * t64 + t50 * t82 + t51 * t81) + (t216 * t253 - t218 * t254) * t217; t248 * t272 + t253 * t205 + t254 * t203 + t247 * t200 + t251 * t189 + t252 * t187 + m(7) * (t20 * t35 + t33 * t40 + t34 * t41) + m(6) * (t39 * t64 + t52 * t81 + t53 * t82); m(6) * (t82 * t203 + t81 * t205 - t272 * t64) + m(7) * (t41 * t203 + t40 * t205 - t272 * t35); (t14 + t13) * t200 + (t4 + t3) * t189 + (t2 + t1) * t187 + m(7) * (t35 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t64 ^ 2 + t81 ^ 2 + t82 ^ 2); m(7) * t200; m(7) * (t187 * t28 + t189 * t29 + t19 * t200); m(7) * (t187 * t34 + t189 * t33 + t20 * t200); m(7) * (t187 * t203 + t189 * t205 - t200 * t272); m(7) * (t187 * t41 + t189 * t40 + t200 * t35); m(7) * (t187 ^ 2 + t189 ^ 2 + t200 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t95(1) t95(2) t95(4) t95(7) t95(11) t95(16); t95(2) t95(3) t95(5) t95(8) t95(12) t95(17); t95(4) t95(5) t95(6) t95(9) t95(13) t95(18); t95(7) t95(8) t95(9) t95(10) t95(14) t95(19); t95(11) t95(12) t95(13) t95(14) t95(15) t95(20); t95(16) t95(17) t95(18) t95(19) t95(20) t95(21);];
Mq  = res;
