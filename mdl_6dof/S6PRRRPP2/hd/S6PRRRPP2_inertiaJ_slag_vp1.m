% Calculate joint inertia matrix for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:58
% EndTime: 2019-03-08 22:49:08
% DurationCPUTime: 5.07s
% Computational Cost: add. (18410->578), mult. (47614->810), div. (0->0), fcn. (61223->10), ass. (0->249)
t275 = m(6) + m(7);
t277 = rSges(7,1) + pkin(5);
t276 = -rSges(7,3) - qJ(6);
t274 = cos(qJ(3));
t273 = cos(qJ(4));
t221 = sin(pkin(10));
t222 = sin(pkin(6));
t272 = t221 * t222;
t223 = cos(pkin(10));
t271 = t222 * t223;
t226 = sin(qJ(3));
t270 = t222 * t226;
t228 = cos(qJ(2));
t269 = t222 * t228;
t224 = cos(pkin(6));
t227 = sin(qJ(2));
t268 = t224 * t227;
t267 = t224 * t228;
t211 = t221 * t228 + t223 * t268;
t198 = t211 * t274 - t223 * t270;
t210 = t221 * t227 - t223 * t267;
t225 = sin(qJ(4));
t172 = t198 * t225 - t210 * t273;
t173 = t198 * t273 + t210 * t225;
t247 = t222 * t274;
t196 = t211 * t226 + t223 * t247;
t266 = rSges(7,2) * t172 + t277 * t173 + t276 * t196;
t213 = -t221 * t268 + t223 * t228;
t201 = t213 * t274 + t221 * t270;
t212 = t221 * t267 + t223 * t227;
t174 = t201 * t225 - t212 * t273;
t175 = t201 * t273 + t212 * t225;
t199 = t213 * t226 - t221 * t247;
t265 = rSges(7,2) * t174 + t277 * t175 + t276 * t199;
t125 = rSges(6,1) * t175 + rSges(6,2) * t199 + rSges(6,3) * t174;
t136 = pkin(4) * t175 + qJ(5) * t174;
t264 = -t125 - t136;
t126 = rSges(5,1) * t175 - rSges(5,2) * t174 + rSges(5,3) * t199;
t170 = pkin(3) * t201 + pkin(9) * t199;
t263 = -t126 - t170;
t135 = pkin(4) * t173 + qJ(5) * t172;
t169 = pkin(3) * t198 + pkin(9) * t196;
t158 = t212 * t169;
t262 = t212 * t135 + t158;
t216 = t224 * t226 + t227 * t247;
t202 = t216 * t225 + t273 * t269;
t203 = t216 * t273 - t225 * t269;
t214 = -t224 * t274 + t227 * t270;
t261 = rSges(7,2) * t202 + t277 * t203 + t276 * t214;
t160 = rSges(6,1) * t203 + rSges(6,2) * t214 + rSges(6,3) * t202;
t171 = pkin(4) * t203 + qJ(5) * t202;
t260 = -t160 - t171;
t161 = rSges(5,1) * t203 - rSges(5,2) * t202 + rSges(5,3) * t214;
t195 = pkin(3) * t216 + pkin(9) * t214;
t259 = -t161 - t195;
t258 = t169 * t269 + t210 * t195;
t194 = pkin(2) * t213 + pkin(8) * t212;
t192 = t224 * t194;
t257 = t224 * t170 + t192;
t193 = pkin(2) * t211 + pkin(8) * t210;
t256 = -t169 - t193;
t255 = t193 * t272 + t194 * t271;
t253 = -t136 - t265;
t252 = -t170 + t264;
t251 = t224 * t136 + t257;
t250 = -t135 + t256;
t249 = -t171 - t261;
t248 = -t195 + t260;
t189 = t216 * rSges(4,1) - t214 * rSges(4,2) - rSges(4,3) * t269;
t217 = (pkin(2) * t227 - pkin(8) * t228) * t222;
t246 = (-t189 - t217) * t222;
t103 = Icges(7,5) * t173 + Icges(7,6) * t172 - Icges(7,3) * t196;
t109 = Icges(7,4) * t173 + Icges(7,2) * t172 - Icges(7,6) * t196;
t115 = Icges(7,1) * t173 + Icges(7,4) * t172 - Icges(7,5) * t196;
t48 = -t103 * t196 + t109 * t172 + t115 * t173;
t104 = Icges(7,5) * t175 + Icges(7,6) * t174 - Icges(7,3) * t199;
t110 = Icges(7,4) * t175 + Icges(7,2) * t174 - Icges(7,6) * t199;
t116 = Icges(7,1) * t175 + Icges(7,4) * t174 - Icges(7,5) * t199;
t49 = -t104 * t196 + t110 * t172 + t116 * t173;
t149 = Icges(7,5) * t203 + Icges(7,6) * t202 - Icges(7,3) * t214;
t152 = Icges(7,4) * t203 + Icges(7,2) * t202 - Icges(7,6) * t214;
t155 = Icges(7,1) * t203 + Icges(7,4) * t202 - Icges(7,5) * t214;
t74 = -t149 * t196 + t152 * t172 + t155 * t173;
t1 = t196 * t48 + t199 * t49 + t214 * t74;
t105 = Icges(6,5) * t173 + Icges(6,6) * t196 + Icges(6,3) * t172;
t111 = Icges(6,4) * t173 + Icges(6,2) * t196 + Icges(6,6) * t172;
t117 = Icges(6,1) * t173 + Icges(6,4) * t196 + Icges(6,5) * t172;
t50 = t105 * t172 + t111 * t196 + t117 * t173;
t106 = Icges(6,5) * t175 + Icges(6,6) * t199 + Icges(6,3) * t174;
t112 = Icges(6,4) * t175 + Icges(6,2) * t199 + Icges(6,6) * t174;
t118 = Icges(6,1) * t175 + Icges(6,4) * t199 + Icges(6,5) * t174;
t51 = t106 * t172 + t112 * t196 + t118 * t173;
t150 = Icges(6,5) * t203 + Icges(6,6) * t214 + Icges(6,3) * t202;
t153 = Icges(6,4) * t203 + Icges(6,2) * t214 + Icges(6,6) * t202;
t156 = Icges(6,1) * t203 + Icges(6,4) * t214 + Icges(6,5) * t202;
t75 = t150 * t172 + t153 * t196 + t156 * t173;
t2 = t196 * t50 + t199 * t51 + t214 * t75;
t107 = Icges(5,5) * t173 - Icges(5,6) * t172 + Icges(5,3) * t196;
t113 = Icges(5,4) * t173 - Icges(5,2) * t172 + Icges(5,6) * t196;
t119 = Icges(5,1) * t173 - Icges(5,4) * t172 + Icges(5,5) * t196;
t52 = t107 * t196 - t113 * t172 + t119 * t173;
t108 = Icges(5,5) * t175 - Icges(5,6) * t174 + Icges(5,3) * t199;
t114 = Icges(5,4) * t175 - Icges(5,2) * t174 + Icges(5,6) * t199;
t120 = Icges(5,1) * t175 - Icges(5,4) * t174 + Icges(5,5) * t199;
t53 = t108 * t196 - t114 * t172 + t120 * t173;
t151 = Icges(5,5) * t203 - Icges(5,6) * t202 + Icges(5,3) * t214;
t154 = Icges(5,4) * t203 - Icges(5,2) * t202 + Icges(5,6) * t214;
t157 = Icges(5,1) * t203 - Icges(5,4) * t202 + Icges(5,5) * t214;
t76 = t151 * t196 - t154 * t172 + t157 * t173;
t3 = t196 * t52 + t199 * t53 + t214 * t76;
t245 = t3 / 0.2e1 + t2 / 0.2e1 + t1 / 0.2e1;
t54 = -t103 * t199 + t109 * t174 + t115 * t175;
t55 = -t104 * t199 + t110 * t174 + t116 * t175;
t77 = -t149 * t199 + t152 * t174 + t155 * t175;
t4 = t196 * t54 + t199 * t55 + t214 * t77;
t56 = t105 * t174 + t111 * t199 + t117 * t175;
t57 = t106 * t174 + t112 * t199 + t118 * t175;
t78 = t150 * t174 + t153 * t199 + t156 * t175;
t5 = t196 * t56 + t199 * t57 + t214 * t78;
t58 = t107 * t199 - t113 * t174 + t119 * t175;
t59 = t108 * t199 - t114 * t174 + t120 * t175;
t79 = t151 * t199 - t154 * t174 + t157 * t175;
t6 = t196 * t58 + t199 * t59 + t214 * t79;
t244 = t6 / 0.2e1 + t5 / 0.2e1 + t4 / 0.2e1;
t7 = t48 * t210 + t49 * t212 - t74 * t269;
t8 = t50 * t210 + t51 * t212 - t75 * t269;
t9 = t52 * t210 + t53 * t212 - t76 * t269;
t243 = t8 / 0.2e1 + t7 / 0.2e1 + t9 / 0.2e1;
t242 = -t170 + t253;
t241 = t135 * t269 + t210 * t171 + t258;
t240 = -t195 + t249;
t239 = t169 * t272 + t170 * t271 + t255;
t10 = t54 * t210 + t55 * t212 - t77 * t269;
t11 = t56 * t210 + t57 * t212 - t78 * t269;
t12 = t58 * t210 + t59 * t212 - t79 * t269;
t238 = t10 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1;
t13 = t74 * t224 + (t221 * t49 - t223 * t48) * t222;
t14 = t75 * t224 + (t221 * t51 - t223 * t50) * t222;
t15 = t76 * t224 + (t221 * t53 - t223 * t52) * t222;
t237 = t14 / 0.2e1 + t13 / 0.2e1 + t15 / 0.2e1;
t16 = t77 * t224 + (t221 * t55 - t223 * t54) * t222;
t17 = t224 * t78 + (t221 * t57 - t223 * t56) * t222;
t18 = t224 * t79 + (t221 * t59 - t223 * t58) * t222;
t236 = t18 / 0.2e1 + t17 / 0.2e1 + t16 / 0.2e1;
t63 = -t103 * t214 + t109 * t202 + t115 * t203;
t64 = -t104 * t214 + t110 * t202 + t116 * t203;
t87 = -t149 * t214 + t152 * t202 + t155 * t203;
t19 = t196 * t63 + t199 * t64 + t214 * t87;
t65 = t105 * t202 + t111 * t214 + t117 * t203;
t66 = t106 * t202 + t112 * t214 + t118 * t203;
t88 = t150 * t202 + t153 * t214 + t156 * t203;
t20 = t196 * t65 + t199 * t66 + t214 * t88;
t67 = t107 * t214 - t113 * t202 + t119 * t203;
t68 = t108 * t214 - t114 * t202 + t120 * t203;
t89 = t151 * t214 - t154 * t202 + t157 * t203;
t21 = t196 * t67 + t199 * t68 + t214 * t89;
t235 = t20 / 0.2e1 + t21 / 0.2e1 + t19 / 0.2e1;
t22 = t63 * t210 + t64 * t212 - t87 * t269;
t23 = t65 * t210 + t66 * t212 - t88 * t269;
t24 = t67 * t210 + t68 * t212 - t89 * t269;
t234 = t24 / 0.2e1 + t23 / 0.2e1 + t22 / 0.2e1;
t25 = t224 * t87 + (t221 * t64 - t223 * t63) * t222;
t26 = t224 * t88 + (t221 * t66 - t223 * t65) * t222;
t27 = t224 * t89 + (t221 * t68 - t223 * t67) * t222;
t233 = t26 / 0.2e1 + t25 / 0.2e1 + t27 / 0.2e1;
t232 = (-t217 + t259) * t222;
t231 = (-t217 + t248) * t222;
t230 = t135 * t272 + t136 * t271 + t239;
t229 = (-t217 + t240) * t222;
t207 = t224 * rSges(3,3) + (rSges(3,1) * t227 + rSges(3,2) * t228) * t222;
t206 = Icges(3,5) * t224 + (Icges(3,1) * t227 + Icges(3,4) * t228) * t222;
t205 = Icges(3,6) * t224 + (Icges(3,4) * t227 + Icges(3,2) * t228) * t222;
t204 = Icges(3,3) * t224 + (Icges(3,5) * t227 + Icges(3,6) * t228) * t222;
t188 = Icges(4,1) * t216 - Icges(4,4) * t214 - Icges(4,5) * t269;
t187 = Icges(4,4) * t216 - Icges(4,2) * t214 - Icges(4,6) * t269;
t186 = Icges(4,5) * t216 - Icges(4,6) * t214 - Icges(4,3) * t269;
t185 = rSges(3,1) * t213 - rSges(3,2) * t212 + rSges(3,3) * t272;
t184 = rSges(3,1) * t211 - rSges(3,2) * t210 - rSges(3,3) * t271;
t183 = Icges(3,1) * t213 - Icges(3,4) * t212 + Icges(3,5) * t272;
t182 = Icges(3,1) * t211 - Icges(3,4) * t210 - Icges(3,5) * t271;
t181 = Icges(3,4) * t213 - Icges(3,2) * t212 + Icges(3,6) * t272;
t180 = Icges(3,4) * t211 - Icges(3,2) * t210 - Icges(3,6) * t271;
t179 = Icges(3,5) * t213 - Icges(3,6) * t212 + Icges(3,3) * t272;
t178 = Icges(3,5) * t211 - Icges(3,6) * t210 - Icges(3,3) * t271;
t164 = -t184 * t224 - t207 * t271;
t163 = t185 * t224 - t207 * t272;
t148 = rSges(4,1) * t201 - rSges(4,2) * t199 + rSges(4,3) * t212;
t147 = rSges(4,1) * t198 - rSges(4,2) * t196 + rSges(4,3) * t210;
t146 = Icges(4,1) * t201 - Icges(4,4) * t199 + Icges(4,5) * t212;
t145 = Icges(4,1) * t198 - Icges(4,4) * t196 + Icges(4,5) * t210;
t144 = Icges(4,4) * t201 - Icges(4,2) * t199 + Icges(4,6) * t212;
t143 = Icges(4,4) * t198 - Icges(4,2) * t196 + Icges(4,6) * t210;
t142 = Icges(4,5) * t201 - Icges(4,6) * t199 + Icges(4,3) * t212;
t141 = Icges(4,5) * t198 - Icges(4,6) * t196 + Icges(4,3) * t210;
t138 = t196 * t171;
t137 = (t184 * t221 + t185 * t223) * t222;
t130 = t214 * t136;
t127 = t199 * t135;
t123 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t196;
t122 = rSges(6,1) * t173 + rSges(6,2) * t196 + rSges(6,3) * t172;
t102 = -t148 * t269 - t212 * t189;
t101 = t147 * t269 + t210 * t189;
t100 = -t186 * t269 - t214 * t187 + t216 * t188;
t99 = t147 * t212 - t148 * t210;
t98 = (-t147 - t193) * t224 + t223 * t246;
t97 = t148 * t224 + t221 * t246 + t192;
t96 = t186 * t212 - t187 * t199 + t188 * t201;
t95 = t186 * t210 - t187 * t196 + t188 * t198;
t94 = (t147 * t221 + t148 * t223) * t222 + t255;
t93 = t126 * t214 - t161 * t199;
t92 = -t123 * t214 + t161 * t196;
t91 = -t142 * t269 - t214 * t144 + t216 * t146;
t90 = -t141 * t269 - t214 * t143 + t216 * t145;
t86 = t142 * t212 - t144 * t199 + t146 * t201;
t85 = t141 * t212 - t143 * t199 + t145 * t201;
t84 = t142 * t210 - t144 * t196 + t146 * t198;
t83 = t141 * t210 - t143 * t196 + t145 * t198;
t82 = t123 * t199 - t126 * t196;
t81 = t259 * t212 + t263 * t269;
t80 = t123 * t269 + t210 * t161 + t258;
t73 = (-t123 + t256) * t224 + t223 * t232;
t72 = t126 * t224 + t221 * t232 + t257;
t71 = t123 * t212 + t263 * t210 + t158;
t70 = t125 * t214 + t260 * t199 + t130;
t69 = t160 * t196 + t138 + (-t122 - t135) * t214;
t62 = (t123 * t221 + t126 * t223) * t222 + t239;
t61 = t248 * t212 + t252 * t269;
t60 = t122 * t269 + t210 * t160 + t241;
t47 = (-t122 + t250) * t224 + t223 * t231;
t46 = t125 * t224 + t221 * t231 + t251;
t45 = t122 * t199 + t264 * t196 + t127;
t44 = t249 * t199 + t265 * t214 + t130;
t43 = t138 + t261 * t196 + (-t135 - t266) * t214;
t42 = t240 * t212 + t242 * t269;
t41 = t261 * t210 + t266 * t269 + t241;
t40 = (t250 - t266) * t224 + t223 * t229;
t39 = t221 * t229 + t265 * t224 + t251;
t38 = t122 * t212 + t252 * t210 + t262;
t37 = (t122 * t221 + t125 * t223) * t222 + t230;
t36 = t100 * t224 + (t221 * t91 - t223 * t90) * t222;
t35 = -t100 * t269 + t90 * t210 + t91 * t212;
t34 = t253 * t196 + t266 * t199 + t127;
t33 = t224 * t96 + (t221 * t86 - t223 * t85) * t222;
t32 = t224 * t95 + (t221 * t84 - t223 * t83) * t222;
t31 = t85 * t210 + t86 * t212 - t96 * t269;
t30 = t83 * t210 + t84 * t212 - t95 * t269;
t29 = t242 * t210 + t266 * t212 + t262;
t28 = (t266 * t221 + t265 * t223) * t222 + t230;
t121 = [m(2) + m(3) + m(4) + m(5) + t275; m(3) * t137 + m(4) * t94 + m(5) * t62 + m(6) * t37 + m(7) * t28; m(7) * (t28 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t37 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t62 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(4) * (t94 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(3) * (t137 ^ 2 + t163 ^ 2 + t164 ^ 2) + (t18 + t17 + t16 + t33 + (t179 * t272 - t181 * t212 + t183 * t213) * t272) * t272 + (-t15 - t14 - t13 - t32 + (-t178 * t271 - t180 * t210 + t182 * t211) * t271 + (-t178 * t272 + t179 * t271 + t180 * t212 + t181 * t210 - t182 * t213 - t183 * t211) * t272) * t271 + ((t204 * t272 - t212 * t205 + t213 * t206) * t272 - (-t204 * t271 - t210 * t205 + t211 * t206) * t271 + t26 + t25 + t27 + t36 + ((t181 * t228 + t183 * t227) * t221 - (t180 * t228 + t182 * t227) * t223) * t222 ^ 2 + ((-t178 * t223 + t179 * t221 + t205 * t228 + t206 * t227) * t222 + t224 * t204) * t224) * t224; m(4) * t99 + m(5) * t71 + m(6) * t38 + m(7) * t29; (t35 / 0.2e1 + t234) * t224 + (t33 / 0.2e1 + t236) * t212 + (t32 / 0.2e1 + t237) * t210 + m(7) * (t28 * t29 + t39 * t42 + t40 * t41) + m(6) * (t37 * t38 + t46 * t61 + t47 * t60) + m(5) * (t62 * t71 + t72 * t81 + t73 * t80) + m(4) * (t101 * t98 + t102 * t97 + t94 * t99) + ((-t36 / 0.2e1 - t233) * t228 + (-t30 / 0.2e1 - t243) * t223 + (t31 / 0.2e1 + t238) * t221) * t222; (-t22 - t23 - t24 - t35) * t269 + (t12 + t10 + t11 + t31) * t212 + (t9 + t8 + t7 + t30) * t210 + m(7) * (t29 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t38 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t71 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2 + t99 ^ 2); m(5) * t82 + m(6) * t45 + m(7) * t34; t235 * t224 + t233 * t214 + t236 * t199 + t237 * t196 + m(7) * (t28 * t34 + t39 * t44 + t40 * t43) + m(6) * (t37 * t45 + t46 * t70 + t47 * t69) + m(5) * (t62 * t82 + t72 * t93 + t73 * t92) + (t244 * t221 - t245 * t223) * t222; -t235 * t269 + t234 * t214 + t244 * t212 + t245 * t210 + t238 * t199 + t243 * t196 + m(7) * (t29 * t34 + t41 * t43 + t42 * t44) + m(6) * (t38 * t45 + t60 * t69 + t61 * t70) + m(5) * (t71 * t82 + t80 * t92 + t81 * t93); (t21 + t20 + t19) * t214 + (t6 + t5 + t4) * t199 + (t3 + t2 + t1) * t196 + m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t45 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t82 ^ 2 + t92 ^ 2 + t93 ^ 2); t202 * t275; m(7) * (t172 * t39 + t174 * t40 + t202 * t28) + m(6) * (t172 * t46 + t174 * t47 + t202 * t37); m(7) * (t172 * t42 + t174 * t41 + t202 * t29) + m(6) * (t172 * t61 + t174 * t60 + t202 * t38); m(7) * (t172 * t44 + t174 * t43 + t202 * t34) + m(6) * (t172 * t70 + t174 * t69 + t202 * t45); (t172 ^ 2 + t174 ^ 2 + t202 ^ 2) * t275; -m(7) * t214; m(7) * (-t196 * t39 - t199 * t40 - t214 * t28); m(7) * (-t196 * t42 - t199 * t41 - t214 * t29); m(7) * (-t196 * t44 - t199 * t43 - t214 * t34); m(7) * (-t172 * t196 - t174 * t199 - t202 * t214); m(7) * (t196 ^ 2 + t199 ^ 2 + t214 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t121(1) t121(2) t121(4) t121(7) t121(11) t121(16); t121(2) t121(3) t121(5) t121(8) t121(12) t121(17); t121(4) t121(5) t121(6) t121(9) t121(13) t121(18); t121(7) t121(8) t121(9) t121(10) t121(14) t121(19); t121(11) t121(12) t121(13) t121(14) t121(15) t121(20); t121(16) t121(17) t121(18) t121(19) t121(20) t121(21);];
Mq  = res;
