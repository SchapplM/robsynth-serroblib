% Calculate joint inertia matrix for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:46
% EndTime: 2019-03-08 18:59:57
% DurationCPUTime: 5.61s
% Computational Cost: add. (42937->474), mult. (101871->684), div. (0->0), fcn. (134256->16), ass. (0->231)
t219 = sin(pkin(12));
t221 = cos(pkin(12));
t222 = cos(pkin(6));
t264 = sin(pkin(13));
t240 = t222 * t264;
t266 = cos(pkin(13));
t211 = t219 * t266 + t221 * t240;
t225 = sin(qJ(3));
t241 = t222 * t266;
t233 = t219 * t264 - t221 * t241;
t267 = cos(pkin(7));
t230 = t233 * t267;
t220 = sin(pkin(6));
t265 = sin(pkin(7));
t270 = cos(qJ(3));
t238 = t270 * t265;
t234 = t220 * t238;
t194 = t211 * t225 + t221 * t234 + t230 * t270;
t212 = -t219 * t240 + t221 * t266;
t232 = t219 * t241 + t221 * t264;
t229 = t232 * t267;
t196 = t212 * t225 - t219 * t234 + t229 * t270;
t236 = t267 * t266;
t202 = -t222 * t238 + (t225 * t264 - t236 * t270) * t220;
t242 = t220 * t265;
t195 = t211 * t270 + (-t221 * t242 - t230) * t225;
t243 = t220 * t267;
t204 = -t221 * t243 + t233 * t265;
t260 = qJ(4) + qJ(5);
t217 = sin(t260);
t244 = cos(t260);
t177 = t195 * t244 + t204 * t217;
t223 = sin(qJ(6));
t226 = cos(qJ(6));
t150 = -t177 * t223 + t194 * t226;
t151 = t177 * t226 + t194 * t223;
t176 = t195 * t217 - t204 * t244;
t102 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t176;
t104 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t176;
t106 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t176;
t52 = t102 * t176 + t104 * t150 + t106 * t151;
t197 = t212 * t270 + (t219 * t242 - t229) * t225;
t205 = t219 * t243 + t232 * t265;
t179 = t197 * t244 + t205 * t217;
t152 = -t179 * t223 + t196 * t226;
t153 = t179 * t226 + t196 * t223;
t178 = t197 * t217 - t205 * t244;
t103 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t178;
t105 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t178;
t107 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t178;
t53 = t103 * t176 + t105 * t150 + t107 * t151;
t203 = t222 * t265 * t225 + (t225 * t236 + t264 * t270) * t220;
t210 = t222 * t267 - t242 * t266;
t193 = t203 * t244 + t210 * t217;
t180 = -t193 * t223 + t202 * t226;
t181 = t193 * t226 + t202 * t223;
t192 = t203 * t217 - t210 * t244;
t122 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t192;
t123 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t192;
t124 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t192;
t61 = t122 * t176 + t123 * t150 + t124 * t151;
t11 = t194 * t52 + t196 * t53 + t202 * t61;
t126 = Icges(6,5) * t177 - Icges(6,6) * t176 + Icges(6,3) * t194;
t128 = Icges(6,4) * t177 - Icges(6,2) * t176 + Icges(6,6) * t194;
t130 = Icges(6,1) * t177 - Icges(6,4) * t176 + Icges(6,5) * t194;
t69 = t126 * t194 - t128 * t176 + t130 * t177;
t127 = Icges(6,5) * t179 - Icges(6,6) * t178 + Icges(6,3) * t196;
t129 = Icges(6,4) * t179 - Icges(6,2) * t178 + Icges(6,6) * t196;
t131 = Icges(6,1) * t179 - Icges(6,4) * t178 + Icges(6,5) * t196;
t70 = t127 * t194 - t129 * t176 + t131 * t177;
t155 = Icges(6,5) * t193 - Icges(6,6) * t192 + Icges(6,3) * t202;
t156 = Icges(6,4) * t193 - Icges(6,2) * t192 + Icges(6,6) * t202;
t157 = Icges(6,1) * t193 - Icges(6,4) * t192 + Icges(6,5) * t202;
t86 = t155 * t194 - t156 * t176 + t157 * t177;
t285 = t194 * t69 + t196 * t70 + t202 * t86 + t11;
t54 = t102 * t178 + t104 * t152 + t106 * t153;
t55 = t103 * t178 + t105 * t152 + t107 * t153;
t62 = t122 * t178 + t123 * t152 + t124 * t153;
t12 = t194 * t54 + t196 * t55 + t202 * t62;
t71 = t126 * t196 - t128 * t178 + t130 * t179;
t72 = t127 * t196 - t129 * t178 + t131 * t179;
t87 = t155 * t196 - t156 * t178 + t157 * t179;
t284 = t194 * t71 + t196 * t72 + t202 * t87 + t12;
t15 = t204 * t52 + t205 * t53 + t210 * t61;
t283 = t204 * t69 + t205 * t70 + t210 * t86 + t15;
t16 = t204 * t54 + t205 * t55 + t210 * t62;
t282 = t204 * t71 + t205 * t72 + t210 * t87 + t16;
t56 = t102 * t192 + t104 * t180 + t106 * t181;
t57 = t103 * t192 + t105 * t180 + t107 * t181;
t68 = t122 * t192 + t123 * t180 + t124 * t181;
t22 = t194 * t56 + t196 * t57 + t202 * t68;
t79 = t126 * t202 - t128 * t192 + t130 * t193;
t80 = t127 * t202 - t129 * t192 + t131 * t193;
t93 = t155 * t202 - t156 * t192 + t157 * t193;
t281 = t194 * t79 + t196 * t80 + t202 * t93 + t22;
t24 = t204 * t56 + t205 * t57 + t210 * t68;
t280 = t204 * t79 + t205 * t80 + t210 * t93 + t24;
t108 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t176;
t257 = pkin(5) * t177 + pkin(11) * t176 + t108;
t109 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t178;
t256 = pkin(5) * t179 + pkin(11) * t178 + t109;
t125 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t192;
t250 = pkin(5) * t193 + pkin(11) * t192 + t125;
t279 = t176 / 0.2e1;
t278 = t178 / 0.2e1;
t277 = t192 / 0.2e1;
t276 = t194 / 0.2e1;
t275 = t196 / 0.2e1;
t274 = t202 / 0.2e1;
t273 = t204 / 0.2e1;
t272 = t205 / 0.2e1;
t271 = t210 / 0.2e1;
t227 = cos(qJ(4));
t269 = pkin(4) * t227;
t224 = sin(qJ(4));
t263 = t204 * t224;
t262 = t205 * t224;
t261 = t210 * t224;
t259 = t257 * t196;
t258 = t256 * t202;
t255 = t250 * t194;
t120 = pkin(4) * t263 + pkin(10) * t194 + t195 * t269;
t174 = t195 * pkin(3) + t194 * pkin(9);
t171 = t205 * t174;
t254 = t205 * t120 + t171;
t121 = pkin(4) * t262 + pkin(10) * t196 + t197 * t269;
t175 = t197 * pkin(3) + t196 * pkin(9);
t172 = t210 * t175;
t253 = t210 * t121 + t172;
t132 = rSges(6,1) * t177 - rSges(6,2) * t176 + rSges(6,3) * t194;
t252 = -t120 - t132;
t133 = rSges(6,1) * t179 - rSges(6,2) * t178 + rSges(6,3) * t196;
t251 = -t121 - t133;
t154 = pkin(4) * t261 + pkin(10) * t202 + t203 * t269;
t191 = t203 * pkin(3) + t202 * pkin(9);
t182 = t204 * t191;
t249 = t204 * t154 + t182;
t158 = rSges(6,1) * t193 - rSges(6,2) * t192 + rSges(6,3) * t202;
t248 = -t154 - t158;
t247 = -t120 - t257;
t246 = -t121 - t256;
t245 = -t154 - t250;
t239 = m(3) + m(4) + m(5) + m(6) + m(7);
t18 = t176 * t56 + t178 * t57 + t192 * t68;
t3 = t176 * t52 + t178 * t53 + t192 * t61;
t4 = t176 * t54 + t178 * t55 + t192 * t62;
t237 = t11 * t279 + t12 * t278 + t18 * t274 + t22 * t277 + t4 * t275 + t3 * t276;
t235 = t285 * t194 + t284 * t196 + t281 * t202;
t231 = t281 * t271 + t284 * t272 + t285 * t273 + t280 * t274 + t282 * t275 + t283 * t276;
t199 = t203 * t227 + t261;
t198 = -t203 * t224 + t210 * t227;
t190 = rSges(4,1) * t203 - rSges(4,2) * t202 + rSges(4,3) * t210;
t189 = Icges(4,1) * t203 - Icges(4,4) * t202 + Icges(4,5) * t210;
t188 = Icges(4,4) * t203 - Icges(4,2) * t202 + Icges(4,6) * t210;
t187 = Icges(4,5) * t203 - Icges(4,6) * t202 + Icges(4,3) * t210;
t186 = t197 * t227 + t262;
t185 = -t197 * t224 + t205 * t227;
t184 = t195 * t227 + t263;
t183 = -t195 * t224 + t204 * t227;
t170 = rSges(5,1) * t199 + rSges(5,2) * t198 + rSges(5,3) * t202;
t169 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t202;
t168 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t202;
t167 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t202;
t166 = rSges(4,1) * t197 - rSges(4,2) * t196 + rSges(4,3) * t205;
t165 = rSges(4,1) * t195 - rSges(4,2) * t194 + rSges(4,3) * t204;
t164 = Icges(4,1) * t197 - Icges(4,4) * t196 + Icges(4,5) * t205;
t163 = Icges(4,1) * t195 - Icges(4,4) * t194 + Icges(4,5) * t204;
t162 = Icges(4,4) * t197 - Icges(4,2) * t196 + Icges(4,6) * t205;
t161 = Icges(4,4) * t195 - Icges(4,2) * t194 + Icges(4,6) * t204;
t160 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t205;
t159 = Icges(4,5) * t195 - Icges(4,6) * t194 + Icges(4,3) * t204;
t145 = t194 * t158;
t144 = t194 * t154;
t142 = rSges(5,1) * t186 + rSges(5,2) * t185 + rSges(5,3) * t196;
t141 = rSges(5,1) * t184 + rSges(5,2) * t183 + rSges(5,3) * t194;
t140 = Icges(5,1) * t186 + Icges(5,4) * t185 + Icges(5,5) * t196;
t139 = Icges(5,1) * t184 + Icges(5,4) * t183 + Icges(5,5) * t194;
t138 = Icges(5,4) * t186 + Icges(5,2) * t185 + Icges(5,6) * t196;
t137 = Icges(5,4) * t184 + Icges(5,2) * t183 + Icges(5,6) * t194;
t136 = Icges(5,5) * t186 + Icges(5,6) * t185 + Icges(5,3) * t196;
t135 = Icges(5,5) * t184 + Icges(5,6) * t183 + Icges(5,3) * t194;
t119 = t166 * t210 - t190 * t205;
t118 = -t165 * t210 + t190 * t204;
t117 = t202 * t133;
t114 = t202 * t121;
t113 = t196 * t132;
t111 = t196 * t120;
t110 = t165 * t205 - t166 * t204;
t99 = t142 * t202 - t170 * t196;
t98 = -t141 * t202 + t170 * t194;
t97 = -t158 * t196 + t117;
t96 = -t132 * t202 + t145;
t95 = t167 * t202 + t168 * t198 + t169 * t199;
t94 = t141 * t196 - t142 * t194;
t92 = -t133 * t194 + t113;
t91 = t210 * t142 + t172 + (-t170 - t191) * t205;
t90 = t204 * t170 + t182 + (-t141 - t174) * t210;
t89 = t167 * t196 + t168 * t185 + t169 * t186;
t88 = t167 * t194 + t168 * t183 + t169 * t184;
t85 = t109 * t192 - t125 * t178;
t84 = -t108 * t192 + t125 * t176;
t83 = t205 * t141 + t171 + (-t142 - t175) * t204;
t82 = t136 * t202 + t138 * t198 + t140 * t199;
t81 = t135 * t202 + t137 * t198 + t139 * t199;
t78 = t136 * t196 + t138 * t185 + t140 * t186;
t77 = t135 * t196 + t137 * t185 + t139 * t186;
t76 = t136 * t194 + t138 * t183 + t140 * t184;
t75 = t135 * t194 + t137 * t183 + t139 * t184;
t74 = t196 * t248 + t114 + t117;
t73 = t202 * t252 + t144 + t145;
t67 = t108 * t178 - t109 * t176;
t66 = t210 * t133 + (-t191 + t248) * t205 + t253;
t65 = t204 * t158 + (-t174 + t252) * t210 + t249;
t64 = -t196 * t250 + t258;
t63 = -t202 * t257 + t255;
t60 = t194 * t251 + t111 + t113;
t59 = -t194 * t256 + t259;
t58 = t205 * t132 + (-t175 + t251) * t204 + t254;
t51 = t196 * t245 + t114 + t258;
t50 = t202 * t247 + t144 + t255;
t49 = t256 * t210 + (-t191 + t245) * t205 + t253;
t48 = t250 * t204 + (-t174 + t247) * t210 + t249;
t47 = t194 * t246 + t111 + t259;
t46 = t257 * t205 + (-t175 + t246) * t204 + t254;
t45 = t204 * t81 + t205 * t82 + t210 * t95;
t44 = t194 * t81 + t196 * t82 + t202 * t95;
t38 = t204 * t77 + t205 * t78 + t210 * t89;
t37 = t204 * t75 + t205 * t76 + t210 * t88;
t36 = t194 * t77 + t196 * t78 + t202 * t89;
t35 = t194 * t75 + t196 * t76 + t202 * t88;
t1 = [m(2) + t239; t239 * t222; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1 + m(3) / 0.2e1) * (t222 ^ 2 + (t219 ^ 2 + t221 ^ 2) * t220 ^ 2); m(4) * t110 + m(5) * t83 + m(6) * t58 + m(7) * t46; m(4) * (t110 * t222 + (t118 * t219 - t119 * t221) * t220) + m(5) * (t83 * t222 + (t219 * t90 - t221 * t91) * t220) + m(6) * (t58 * t222 + (t219 * t65 - t221 * t66) * t220) + m(7) * (t46 * t222 + (t219 * t48 - t221 * t49) * t220); (t45 + (t187 * t210 - t188 * t202 + t189 * t203) * t210 + t280) * t210 + (t38 + (t160 * t205 - t162 * t196 + t164 * t197) * t205 + (t160 * t210 - t162 * t202 + t164 * t203 + t187 * t205 - t188 * t196 + t189 * t197) * t210 + t282) * t205 + (t37 + (t159 * t204 - t161 * t194 + t163 * t195) * t204 + (t159 * t210 - t161 * t202 + t163 * t203 + t187 * t204 - t188 * t194 + t189 * t195) * t210 + (t159 * t205 + t160 * t204 - t161 * t196 - t162 * t194 + t163 * t197 + t164 * t195) * t205 + t283) * t204 + m(7) * (t46 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t58 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t83 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(4) * (t110 ^ 2 + t118 ^ 2 + t119 ^ 2); m(5) * t94 + m(6) * t60 + m(7) * t47; m(5) * (t94 * t222 + (t219 * t98 - t221 * t99) * t220) + m(6) * (t60 * t222 + (t219 * t73 - t221 * t74) * t220) + m(7) * (t47 * t222 + (t219 * t50 - t221 * t51) * t220); t35 * t273 + m(7) * (t46 * t47 + t48 * t50 + t49 * t51) + m(6) * (t58 * t60 + t65 * t73 + t66 * t74) + m(5) * (t83 * t94 + t90 * t98 + t91 * t99) + t38 * t275 + t45 * t274 + t36 * t272 + t44 * t271 + t37 * t276 + t231; t194 * t35 + t196 * t36 + t202 * t44 + m(7) * (t47 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t60 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t94 ^ 2 + t98 ^ 2 + t99 ^ 2) + t235; m(6) * t92 + m(7) * t59; m(6) * (t92 * t222 + (t219 * t96 - t221 * t97) * t220) + m(7) * (t59 * t222 + (t219 * t63 - t221 * t64) * t220); m(7) * (t46 * t59 + t48 * t63 + t49 * t64) + m(6) * (t58 * t92 + t65 * t96 + t66 * t97) + t231; m(7) * (t47 * t59 + t50 * t63 + t51 * t64) + m(6) * (t60 * t92 + t73 * t96 + t74 * t97) + t235; m(7) * (t59 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t92 ^ 2 + t96 ^ 2 + t97 ^ 2) + t235; m(7) * t67; m(7) * (t67 * t222 + (t219 * t84 - t221 * t85) * t220); m(7) * (t46 * t67 + t48 * t84 + t49 * t85) + t18 * t271 + t3 * t273 + t24 * t277 + t15 * t279 + t4 * t272 + t16 * t278; m(7) * (t47 * t67 + t50 * t84 + t51 * t85) + t237; m(7) * (t59 * t67 + t63 * t84 + t64 * t85) + t237; t178 * t4 + t176 * t3 + t192 * t18 + m(7) * (t67 ^ 2 + t84 ^ 2 + t85 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
