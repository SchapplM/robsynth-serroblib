% Calculate joint inertia matrix for
% S6PRRRPP3
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:28
% EndTime: 2019-03-08 22:54:36
% DurationCPUTime: 4.68s
% Computational Cost: add. (18436->577), mult. (47682->810), div. (0->0), fcn. (61317->10), ass. (0->249)
t273 = m(6) + m(7);
t275 = rSges(7,1) + pkin(5);
t274 = rSges(7,3) + qJ(6);
t272 = cos(qJ(3));
t271 = cos(qJ(4));
t270 = sin(qJ(4));
t219 = sin(pkin(10));
t220 = sin(pkin(6));
t269 = t219 * t220;
t221 = cos(pkin(10));
t268 = t220 * t221;
t223 = sin(qJ(3));
t267 = t220 * t223;
t225 = cos(qJ(2));
t266 = t220 * t225;
t222 = cos(pkin(6));
t224 = sin(qJ(2));
t265 = t222 * t224;
t264 = t222 * t225;
t209 = t219 * t225 + t221 * t265;
t197 = t209 * t272 - t221 * t267;
t208 = t219 * t224 - t221 * t264;
t172 = t197 * t270 - t208 * t271;
t173 = t197 * t271 + t208 * t270;
t244 = t220 * t272;
t196 = t209 * t223 + t221 * t244;
t263 = rSges(7,2) * t172 + t173 * t274 + t196 * t275;
t211 = -t219 * t265 + t221 * t225;
t199 = t211 * t272 + t219 * t267;
t210 = t219 * t264 + t221 * t224;
t174 = t199 * t270 - t210 * t271;
t175 = t199 * t271 + t210 * t270;
t198 = t211 * t223 - t219 * t244;
t262 = rSges(7,2) * t174 + t175 * t274 + t198 * t275;
t124 = rSges(6,1) * t198 - rSges(6,2) * t175 + rSges(6,3) * t174;
t136 = pkin(4) * t175 + qJ(5) * t174;
t261 = -t124 - t136;
t126 = rSges(5,1) * t175 - rSges(5,2) * t174 + rSges(5,3) * t198;
t170 = pkin(3) * t199 + pkin(9) * t198;
t260 = -t126 - t170;
t135 = pkin(4) * t173 + qJ(5) * t172;
t169 = pkin(3) * t197 + pkin(9) * t196;
t158 = t210 * t169;
t259 = t210 * t135 + t158;
t213 = t222 * t223 + t224 * t244;
t200 = t213 * t270 + t271 * t266;
t201 = t213 * t271 - t270 * t266;
t212 = -t222 * t272 + t224 * t267;
t258 = rSges(7,2) * t200 + t201 * t274 + t212 * t275;
t160 = rSges(6,1) * t212 - rSges(6,2) * t201 + rSges(6,3) * t200;
t171 = pkin(4) * t201 + qJ(5) * t200;
t257 = -t160 - t171;
t161 = rSges(5,1) * t201 - rSges(5,2) * t200 + rSges(5,3) * t212;
t195 = pkin(3) * t213 + pkin(9) * t212;
t256 = -t161 - t195;
t255 = t169 * t266 + t208 * t195;
t194 = pkin(2) * t211 + pkin(8) * t210;
t192 = t222 * t194;
t254 = t222 * t170 + t192;
t193 = pkin(2) * t209 + pkin(8) * t208;
t253 = -t169 - t193;
t252 = t193 * t269 + t194 * t268;
t250 = -t136 - t262;
t249 = -t170 + t261;
t248 = t222 * t136 + t254;
t247 = -t135 + t253;
t246 = -t171 - t258;
t245 = -t195 + t257;
t189 = t213 * rSges(4,1) - t212 * rSges(4,2) - rSges(4,3) * t266;
t214 = (pkin(2) * t224 - pkin(8) * t225) * t220;
t243 = (-t189 - t214) * t220;
t103 = Icges(7,5) * t196 + Icges(7,6) * t172 + Icges(7,3) * t173;
t109 = Icges(7,4) * t196 + Icges(7,2) * t172 + Icges(7,6) * t173;
t115 = Icges(7,1) * t196 + Icges(7,4) * t172 + Icges(7,5) * t173;
t48 = t103 * t173 + t109 * t172 + t115 * t196;
t104 = Icges(7,5) * t198 + Icges(7,6) * t174 + Icges(7,3) * t175;
t110 = Icges(7,4) * t198 + Icges(7,2) * t174 + Icges(7,6) * t175;
t116 = Icges(7,1) * t198 + Icges(7,4) * t174 + Icges(7,5) * t175;
t49 = t104 * t173 + t110 * t172 + t116 * t196;
t149 = Icges(7,5) * t212 + Icges(7,6) * t200 + Icges(7,3) * t201;
t152 = Icges(7,4) * t212 + Icges(7,2) * t200 + Icges(7,6) * t201;
t155 = Icges(7,1) * t212 + Icges(7,4) * t200 + Icges(7,5) * t201;
t74 = t149 * t173 + t152 * t172 + t155 * t196;
t1 = t196 * t48 + t198 * t49 + t212 * t74;
t105 = Icges(6,5) * t196 - Icges(6,6) * t173 + Icges(6,3) * t172;
t111 = Icges(6,4) * t196 - Icges(6,2) * t173 + Icges(6,6) * t172;
t117 = Icges(6,1) * t196 - Icges(6,4) * t173 + Icges(6,5) * t172;
t50 = t105 * t172 - t111 * t173 + t117 * t196;
t106 = Icges(6,5) * t198 - Icges(6,6) * t175 + Icges(6,3) * t174;
t112 = Icges(6,4) * t198 - Icges(6,2) * t175 + Icges(6,6) * t174;
t118 = Icges(6,1) * t198 - Icges(6,4) * t175 + Icges(6,5) * t174;
t51 = t106 * t172 - t112 * t173 + t118 * t196;
t150 = Icges(6,5) * t212 - Icges(6,6) * t201 + Icges(6,3) * t200;
t153 = Icges(6,4) * t212 - Icges(6,2) * t201 + Icges(6,6) * t200;
t156 = Icges(6,1) * t212 - Icges(6,4) * t201 + Icges(6,5) * t200;
t75 = t150 * t172 - t153 * t173 + t156 * t196;
t2 = t196 * t50 + t198 * t51 + t212 * t75;
t107 = Icges(5,5) * t173 - Icges(5,6) * t172 + Icges(5,3) * t196;
t113 = Icges(5,4) * t173 - Icges(5,2) * t172 + Icges(5,6) * t196;
t119 = Icges(5,1) * t173 - Icges(5,4) * t172 + Icges(5,5) * t196;
t56 = t107 * t196 - t113 * t172 + t119 * t173;
t108 = Icges(5,5) * t175 - Icges(5,6) * t174 + Icges(5,3) * t198;
t114 = Icges(5,4) * t175 - Icges(5,2) * t174 + Icges(5,6) * t198;
t120 = Icges(5,1) * t175 - Icges(5,4) * t174 + Icges(5,5) * t198;
t57 = t108 * t196 - t114 * t172 + t120 * t173;
t151 = Icges(5,5) * t201 - Icges(5,6) * t200 + Icges(5,3) * t212;
t154 = Icges(5,4) * t201 - Icges(5,2) * t200 + Icges(5,6) * t212;
t157 = Icges(5,1) * t201 - Icges(5,4) * t200 + Icges(5,5) * t212;
t78 = t151 * t196 - t154 * t172 + t157 * t173;
t5 = t196 * t56 + t198 * t57 + t212 * t78;
t242 = t2 / 0.2e1 + t5 / 0.2e1 + t1 / 0.2e1;
t52 = t103 * t175 + t109 * t174 + t115 * t198;
t53 = t104 * t175 + t110 * t174 + t116 * t198;
t76 = t149 * t175 + t152 * t174 + t155 * t198;
t3 = t196 * t52 + t198 * t53 + t212 * t76;
t54 = t105 * t174 - t111 * t175 + t117 * t198;
t55 = t106 * t174 - t112 * t175 + t118 * t198;
t77 = t150 * t174 - t153 * t175 + t156 * t198;
t4 = t196 * t54 + t198 * t55 + t212 * t77;
t58 = t107 * t198 - t113 * t174 + t119 * t175;
t59 = t108 * t198 - t114 * t174 + t120 * t175;
t79 = t151 * t198 - t154 * t174 + t157 * t175;
t6 = t196 * t58 + t198 * t59 + t212 * t79;
t241 = t3 / 0.2e1 + t6 / 0.2e1 + t4 / 0.2e1;
t11 = t56 * t208 + t57 * t210 - t78 * t266;
t7 = t48 * t208 + t49 * t210 - t74 * t266;
t8 = t50 * t208 + t51 * t210 - t75 * t266;
t240 = t7 / 0.2e1 + t11 / 0.2e1 + t8 / 0.2e1;
t239 = -t170 + t250;
t238 = t135 * t266 + t208 * t171 + t255;
t237 = -t195 + t246;
t236 = t169 * t269 + t170 * t268 + t252;
t10 = t54 * t208 + t55 * t210 - t77 * t266;
t12 = t58 * t208 + t59 * t210 - t79 * t266;
t9 = t52 * t208 + t53 * t210 - t76 * t266;
t235 = t12 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1;
t15 = t222 * t76 + (t219 * t53 - t221 * t52) * t220;
t16 = t222 * t77 + (t219 * t55 - t221 * t54) * t220;
t18 = t222 * t79 + (t219 * t59 - t221 * t58) * t220;
t234 = t15 / 0.2e1 + t18 / 0.2e1 + t16 / 0.2e1;
t13 = t222 * t74 + (t219 * t49 - t221 * t48) * t220;
t14 = t222 * t75 + (t219 * t51 - t221 * t50) * t220;
t17 = t222 * t78 + (t219 * t57 - t221 * t56) * t220;
t233 = t17 / 0.2e1 + t14 / 0.2e1 + t13 / 0.2e1;
t63 = t103 * t201 + t109 * t200 + t115 * t212;
t64 = t104 * t201 + t110 * t200 + t116 * t212;
t87 = t149 * t201 + t152 * t200 + t155 * t212;
t19 = t196 * t63 + t198 * t64 + t212 * t87;
t65 = t105 * t200 - t111 * t201 + t117 * t212;
t66 = t106 * t200 - t112 * t201 + t118 * t212;
t88 = t150 * t200 - t153 * t201 + t156 * t212;
t20 = t196 * t65 + t198 * t66 + t212 * t88;
t67 = t107 * t212 - t113 * t200 + t119 * t201;
t68 = t108 * t212 - t114 * t200 + t120 * t201;
t89 = t151 * t212 - t154 * t200 + t157 * t201;
t21 = t196 * t67 + t198 * t68 + t212 * t89;
t232 = t19 / 0.2e1 + t21 / 0.2e1 + t20 / 0.2e1;
t22 = t63 * t208 + t64 * t210 - t87 * t266;
t23 = t65 * t208 + t66 * t210 - t88 * t266;
t24 = t67 * t208 + t68 * t210 - t89 * t266;
t231 = t22 / 0.2e1 + t24 / 0.2e1 + t23 / 0.2e1;
t25 = t222 * t87 + (t219 * t64 - t221 * t63) * t220;
t26 = t222 * t88 + (t219 * t66 - t221 * t65) * t220;
t27 = t222 * t89 + (t219 * t68 - t221 * t67) * t220;
t230 = t27 / 0.2e1 + t26 / 0.2e1 + t25 / 0.2e1;
t229 = (-t214 + t256) * t220;
t228 = (-t214 + t245) * t220;
t227 = t135 * t269 + t136 * t268 + t236;
t226 = (-t214 + t237) * t220;
t205 = t222 * rSges(3,3) + (rSges(3,1) * t224 + rSges(3,2) * t225) * t220;
t204 = Icges(3,5) * t222 + (Icges(3,1) * t224 + Icges(3,4) * t225) * t220;
t203 = Icges(3,6) * t222 + (Icges(3,4) * t224 + Icges(3,2) * t225) * t220;
t202 = Icges(3,3) * t222 + (Icges(3,5) * t224 + Icges(3,6) * t225) * t220;
t188 = Icges(4,1) * t213 - Icges(4,4) * t212 - Icges(4,5) * t266;
t187 = Icges(4,4) * t213 - Icges(4,2) * t212 - Icges(4,6) * t266;
t186 = Icges(4,5) * t213 - Icges(4,6) * t212 - Icges(4,3) * t266;
t185 = rSges(3,1) * t211 - rSges(3,2) * t210 + rSges(3,3) * t269;
t184 = rSges(3,1) * t209 - rSges(3,2) * t208 - rSges(3,3) * t268;
t183 = Icges(3,1) * t211 - Icges(3,4) * t210 + Icges(3,5) * t269;
t182 = Icges(3,1) * t209 - Icges(3,4) * t208 - Icges(3,5) * t268;
t181 = Icges(3,4) * t211 - Icges(3,2) * t210 + Icges(3,6) * t269;
t180 = Icges(3,4) * t209 - Icges(3,2) * t208 - Icges(3,6) * t268;
t179 = Icges(3,5) * t211 - Icges(3,6) * t210 + Icges(3,3) * t269;
t178 = Icges(3,5) * t209 - Icges(3,6) * t208 - Icges(3,3) * t268;
t164 = -t184 * t222 - t205 * t268;
t163 = t185 * t222 - t205 * t269;
t148 = rSges(4,1) * t199 - rSges(4,2) * t198 + rSges(4,3) * t210;
t147 = rSges(4,1) * t197 - rSges(4,2) * t196 + rSges(4,3) * t208;
t146 = Icges(4,1) * t199 - Icges(4,4) * t198 + Icges(4,5) * t210;
t145 = Icges(4,1) * t197 - Icges(4,4) * t196 + Icges(4,5) * t208;
t144 = Icges(4,4) * t199 - Icges(4,2) * t198 + Icges(4,6) * t210;
t143 = Icges(4,4) * t197 - Icges(4,2) * t196 + Icges(4,6) * t208;
t142 = Icges(4,5) * t199 - Icges(4,6) * t198 + Icges(4,3) * t210;
t141 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t208;
t138 = t196 * t171;
t137 = (t184 * t219 + t185 * t221) * t220;
t130 = t212 * t136;
t127 = t198 * t135;
t125 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t196;
t122 = rSges(6,1) * t196 - rSges(6,2) * t173 + rSges(6,3) * t172;
t102 = -t148 * t266 - t210 * t189;
t101 = t147 * t266 + t208 * t189;
t100 = -t186 * t266 - t212 * t187 + t213 * t188;
t99 = t147 * t210 - t148 * t208;
t98 = (-t147 - t193) * t222 + t221 * t243;
t97 = t148 * t222 + t219 * t243 + t192;
t96 = t186 * t210 - t187 * t198 + t188 * t199;
t95 = t186 * t208 - t187 * t196 + t188 * t197;
t94 = (t147 * t219 + t148 * t221) * t220 + t252;
t93 = t126 * t212 - t161 * t198;
t92 = -t125 * t212 + t161 * t196;
t91 = -t142 * t266 - t212 * t144 + t213 * t146;
t90 = -t141 * t266 - t212 * t143 + t213 * t145;
t86 = t142 * t210 - t144 * t198 + t146 * t199;
t85 = t141 * t210 - t143 * t198 + t145 * t199;
t84 = t142 * t208 - t144 * t196 + t146 * t197;
t83 = t141 * t208 - t143 * t196 + t145 * t197;
t82 = t125 * t198 - t126 * t196;
t81 = t256 * t210 + t260 * t266;
t80 = t125 * t266 + t208 * t161 + t255;
t73 = (-t125 + t253) * t222 + t221 * t229;
t72 = t126 * t222 + t219 * t229 + t254;
t71 = t125 * t210 + t260 * t208 + t158;
t70 = t124 * t212 + t257 * t198 + t130;
t69 = t160 * t196 + t138 + (-t122 - t135) * t212;
t62 = (t125 * t219 + t126 * t221) * t220 + t236;
t61 = t245 * t210 + t249 * t266;
t60 = t122 * t266 + t208 * t160 + t238;
t47 = (-t122 + t247) * t222 + t221 * t228;
t46 = t124 * t222 + t219 * t228 + t248;
t45 = t122 * t198 + t261 * t196 + t127;
t44 = t246 * t198 + t262 * t212 + t130;
t43 = t138 + t258 * t196 + (-t135 - t263) * t212;
t42 = t237 * t210 + t239 * t266;
t41 = t258 * t208 + t263 * t266 + t238;
t40 = (t247 - t263) * t222 + t221 * t226;
t39 = t219 * t226 + t262 * t222 + t248;
t38 = t122 * t210 + t249 * t208 + t259;
t37 = (t122 * t219 + t124 * t221) * t220 + t227;
t36 = t100 * t222 + (t219 * t91 - t221 * t90) * t220;
t35 = -t100 * t266 + t90 * t208 + t91 * t210;
t34 = t250 * t196 + t263 * t198 + t127;
t33 = t222 * t96 + (t219 * t86 - t221 * t85) * t220;
t32 = t222 * t95 + (t219 * t84 - t221 * t83) * t220;
t31 = t85 * t208 + t86 * t210 - t96 * t266;
t30 = t83 * t208 + t84 * t210 - t95 * t266;
t29 = t239 * t208 + t263 * t210 + t259;
t28 = (t263 * t219 + t262 * t221) * t220 + t227;
t121 = [m(2) + m(3) + m(4) + m(5) + t273; m(3) * t137 + m(4) * t94 + m(5) * t62 + m(6) * t37 + m(7) * t28; m(7) * (t28 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t37 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t62 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(4) * (t94 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(3) * (t137 ^ 2 + t163 ^ 2 + t164 ^ 2) + (t15 + t18 + t16 + t33 + (t179 * t269 - t181 * t210 + t183 * t211) * t269) * t269 + (-t13 - t17 - t14 - t32 + (-t178 * t268 - t180 * t208 + t182 * t209) * t268 + (-t178 * t269 + t179 * t268 + t180 * t210 + t181 * t208 - t182 * t211 - t183 * t209) * t269) * t268 + ((t202 * t269 - t210 * t203 + t211 * t204) * t269 - (-t202 * t268 - t208 * t203 + t209 * t204) * t268 + t25 + t27 + t26 + t36 + ((t181 * t225 + t183 * t224) * t219 - (t180 * t225 + t182 * t224) * t221) * t220 ^ 2 + ((-t178 * t221 + t179 * t219 + t203 * t225 + t204 * t224) * t220 + t222 * t202) * t222) * t222; m(4) * t99 + m(5) * t71 + m(6) * t38 + m(7) * t29; (t35 / 0.2e1 + t231) * t222 + (t33 / 0.2e1 + t234) * t210 + (t32 / 0.2e1 + t233) * t208 + m(7) * (t28 * t29 + t39 * t42 + t40 * t41) + m(6) * (t37 * t38 + t46 * t61 + t47 * t60) + m(5) * (t62 * t71 + t72 * t81 + t73 * t80) + m(4) * (t101 * t98 + t102 * t97 + t94 * t99) + ((-t36 / 0.2e1 - t230) * t225 + (-t30 / 0.2e1 - t240) * t221 + (t31 / 0.2e1 + t235) * t219) * t220; (-t22 - t23 - t24 - t35) * t266 + (t12 + t10 + t9 + t31) * t210 + (t8 + t11 + t7 + t30) * t208 + m(7) * (t29 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t38 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t71 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2 + t99 ^ 2); m(5) * t82 + m(6) * t45 + m(7) * t34; t232 * t222 + t230 * t212 + t234 * t198 + t233 * t196 + m(7) * (t28 * t34 + t39 * t44 + t40 * t43) + m(6) * (t37 * t45 + t46 * t70 + t47 * t69) + m(5) * (t62 * t82 + t72 * t93 + t73 * t92) + (t241 * t219 - t242 * t221) * t220; -t232 * t266 + t231 * t212 + t241 * t210 + t242 * t208 + t235 * t198 + t240 * t196 + m(7) * (t29 * t34 + t41 * t43 + t42 * t44) + m(6) * (t38 * t45 + t60 * t69 + t61 * t70) + m(5) * (t71 * t82 + t80 * t92 + t81 * t93); (t19 + t20 + t21) * t212 + (t6 + t3 + t4) * t198 + (t2 + t1 + t5) * t196 + m(7) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t45 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t82 ^ 2 + t92 ^ 2 + t93 ^ 2); t200 * t273; m(7) * (t172 * t39 + t174 * t40 + t200 * t28) + m(6) * (t172 * t46 + t174 * t47 + t200 * t37); m(7) * (t172 * t42 + t174 * t41 + t200 * t29) + m(6) * (t172 * t61 + t174 * t60 + t200 * t38); m(7) * (t172 * t44 + t174 * t43 + t200 * t34) + m(6) * (t172 * t70 + t174 * t69 + t200 * t45); (t172 ^ 2 + t174 ^ 2 + t200 ^ 2) * t273; m(7) * t201; m(7) * (t173 * t39 + t175 * t40 + t201 * t28); m(7) * (t173 * t42 + t175 * t41 + t201 * t29); m(7) * (t173 * t44 + t175 * t43 + t201 * t34); m(7) * (t172 * t173 + t174 * t175 + t200 * t201); m(7) * (t173 ^ 2 + t175 ^ 2 + t201 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t121(1) t121(2) t121(4) t121(7) t121(11) t121(16); t121(2) t121(3) t121(5) t121(8) t121(12) t121(17); t121(4) t121(5) t121(6) t121(9) t121(13) t121(18); t121(7) t121(8) t121(9) t121(10) t121(14) t121(19); t121(11) t121(12) t121(13) t121(14) t121(15) t121(20); t121(16) t121(17) t121(18) t121(19) t121(20) t121(21);];
Mq  = res;
