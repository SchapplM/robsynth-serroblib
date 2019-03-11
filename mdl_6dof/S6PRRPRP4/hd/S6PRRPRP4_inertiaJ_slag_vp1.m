% Calculate joint inertia matrix for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:18
% EndTime: 2019-03-08 21:40:26
% DurationCPUTime: 4.19s
% Computational Cost: add. (14850->552), mult. (37801->791), div. (0->0), fcn. (47967->10), ass. (0->241)
t268 = rSges(7,3) + qJ(6);
t267 = m(5) + m(6) + m(7);
t266 = cos(qJ(3));
t265 = sin(qJ(3));
t216 = cos(qJ(5));
t264 = pkin(5) * t216;
t209 = sin(pkin(10));
t211 = cos(pkin(10));
t217 = cos(qJ(2));
t212 = cos(pkin(6));
t215 = sin(qJ(2));
t256 = t212 * t215;
t198 = t209 * t217 + t211 * t256;
t210 = sin(pkin(6));
t229 = t210 * t266;
t186 = t198 * t265 + t211 * t229;
t214 = sin(qJ(5));
t262 = t186 * t214;
t200 = -t209 * t256 + t211 * t217;
t188 = t200 * t265 - t209 * t229;
t261 = t188 * t214;
t228 = t210 * t265;
t201 = -t212 * t266 + t215 * t228;
t260 = t201 * t214;
t259 = t209 * t210;
t258 = t210 * t211;
t257 = t210 * t217;
t255 = t212 * t217;
t197 = t209 * t215 - t211 * t255;
t156 = t186 * t216 - t197 * t214;
t157 = t197 * t216 + t262;
t187 = t198 * t266 - t211 * t228;
t254 = rSges(7,1) * t157 + rSges(7,2) * t156 + pkin(5) * t262 + t187 * t268 + t264 * t197;
t199 = t209 * t255 + t211 * t215;
t158 = t188 * t216 - t199 * t214;
t159 = t199 * t216 + t261;
t189 = t200 * t266 + t209 * t228;
t253 = rSges(7,1) * t159 + rSges(7,2) * t158 + pkin(5) * t261 + t189 * t268 + t264 * t199;
t129 = rSges(5,1) * t199 - rSges(5,2) * t189 + rSges(5,3) * t188;
t151 = pkin(3) * t189 + qJ(4) * t188;
t252 = -t129 - t151;
t150 = pkin(3) * t187 + qJ(4) * t186;
t132 = t199 * t150;
t161 = t197 * pkin(4) + t187 * pkin(9);
t251 = t199 * t161 + t132;
t190 = t201 * t216 + t214 * t257;
t191 = -t216 * t257 + t260;
t202 = t212 * t265 + t215 * t229;
t250 = rSges(7,1) * t191 + rSges(7,2) * t190 + pkin(5) * t260 + t202 * t268 - t264 * t257;
t185 = pkin(3) * t202 + qJ(4) * t201;
t249 = t150 * t257 + t197 * t185;
t184 = pkin(2) * t200 + pkin(8) * t199;
t182 = t212 * t184;
t248 = t212 * t151 + t182;
t183 = pkin(2) * t198 + pkin(8) * t197;
t247 = -t150 - t183;
t162 = t199 * pkin(4) + t189 * pkin(9);
t246 = -t151 - t162;
t177 = -rSges(5,1) * t257 - t202 * rSges(5,2) + t201 * rSges(5,3);
t245 = -t177 - t185;
t244 = t183 * t259 + t184 * t258;
t192 = -pkin(4) * t257 + t202 * pkin(9);
t243 = -t185 - t192;
t101 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t187;
t105 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t187;
t97 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t187;
t40 = t101 * t156 + t105 * t157 + t187 * t97;
t102 = Icges(7,4) * t159 + Icges(7,2) * t158 + Icges(7,6) * t189;
t106 = Icges(7,1) * t159 + Icges(7,4) * t158 + Icges(7,5) * t189;
t98 = Icges(7,5) * t159 + Icges(7,6) * t158 + Icges(7,3) * t189;
t41 = t102 * t156 + t106 * t157 + t187 * t98;
t133 = Icges(7,5) * t191 + Icges(7,6) * t190 + Icges(7,3) * t202;
t135 = Icges(7,4) * t191 + Icges(7,2) * t190 + Icges(7,6) * t202;
t137 = Icges(7,1) * t191 + Icges(7,4) * t190 + Icges(7,5) * t202;
t59 = t133 * t187 + t135 * t156 + t137 * t157;
t1 = t187 * t40 + t189 * t41 + t202 * t59;
t103 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t187;
t107 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t187;
t99 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t187;
t42 = t103 * t156 + t107 * t157 + t187 * t99;
t100 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t189;
t104 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t189;
t108 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t189;
t43 = t100 * t187 + t104 * t156 + t108 * t157;
t134 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t202;
t136 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t202;
t138 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t202;
t60 = t134 * t187 + t136 * t156 + t138 * t157;
t2 = t187 * t42 + t189 * t43 + t202 * t60;
t242 = -t2 / 0.2e1 - t1 / 0.2e1;
t44 = t101 * t158 + t105 * t159 + t189 * t97;
t45 = t102 * t158 + t106 * t159 + t189 * t98;
t61 = t133 * t189 + t135 * t158 + t137 * t159;
t3 = t187 * t44 + t189 * t45 + t202 * t61;
t46 = t103 * t158 + t107 * t159 + t189 * t99;
t47 = t100 * t189 + t104 * t158 + t108 * t159;
t62 = t134 * t189 + t136 * t158 + t138 * t159;
t4 = t187 * t46 + t189 * t47 + t202 * t62;
t241 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t40 * t197 + t41 * t199 - t59 * t257;
t6 = t42 * t197 + t43 * t199 - t60 * t257;
t240 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t44 * t197 + t45 * t199 - t61 * t257;
t8 = t46 * t197 + t47 * t199 - t62 * t257;
t239 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t60 * t212 + (t209 * t43 - t211 * t42) * t210;
t9 = t59 * t212 + (t209 * t41 - t211 * t40) * t210;
t238 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t61 * t212 + (t209 * t45 - t211 * t44) * t210;
t12 = t62 * t212 + (t209 * t47 - t211 * t46) * t210;
t237 = t11 / 0.2e1 + t12 / 0.2e1;
t50 = t101 * t190 + t105 * t191 + t202 * t97;
t51 = t102 * t190 + t106 * t191 + t202 * t98;
t73 = t133 * t202 + t135 * t190 + t137 * t191;
t13 = t187 * t50 + t189 * t51 + t202 * t73;
t52 = t103 * t190 + t107 * t191 + t202 * t99;
t53 = t100 * t202 + t104 * t190 + t108 * t191;
t74 = t134 * t202 + t136 * t190 + t138 * t191;
t14 = t187 * t52 + t189 * t53 + t202 * t74;
t236 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t50 * t197 + t51 * t199 - t73 * t257;
t16 = t52 * t197 + t53 * t199 - t74 * t257;
t235 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t73 * t212 + (t209 * t51 - t211 * t50) * t210;
t18 = t74 * t212 + (t209 * t53 - t211 * t52) * t210;
t234 = t17 / 0.2e1 + t18 / 0.2e1;
t112 = rSges(6,1) * t159 + rSges(6,2) * t158 + rSges(6,3) * t189;
t233 = -t112 + t246;
t140 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t202;
t232 = -t140 + t243;
t231 = t212 * t162 + t248;
t230 = -t161 + t247;
t178 = t202 * rSges(4,1) - t201 * rSges(4,2) - rSges(4,3) * t257;
t203 = (pkin(2) * t215 - pkin(8) * t217) * t210;
t227 = (-t178 - t203) * t210;
t225 = t246 - t253;
t224 = t243 - t250;
t223 = t150 * t259 + t151 * t258 + t244;
t222 = t161 * t257 + t197 * t192 + t249;
t221 = (-t203 + t245) * t210;
t220 = (-t203 + t232) * t210;
t219 = t161 * t259 + t162 * t258 + t223;
t218 = (-t203 + t224) * t210;
t196 = t212 * rSges(3,3) + (rSges(3,1) * t215 + rSges(3,2) * t217) * t210;
t195 = Icges(3,5) * t212 + (Icges(3,1) * t215 + Icges(3,4) * t217) * t210;
t194 = Icges(3,6) * t212 + (Icges(3,4) * t215 + Icges(3,2) * t217) * t210;
t193 = Icges(3,3) * t212 + (Icges(3,5) * t215 + Icges(3,6) * t217) * t210;
t176 = Icges(4,1) * t202 - Icges(4,4) * t201 - Icges(4,5) * t257;
t175 = Icges(4,4) * t202 - Icges(4,2) * t201 - Icges(4,6) * t257;
t174 = Icges(4,5) * t202 - Icges(4,6) * t201 - Icges(4,3) * t257;
t173 = -Icges(5,1) * t257 - Icges(5,4) * t202 + Icges(5,5) * t201;
t172 = -Icges(5,4) * t257 - Icges(5,2) * t202 + Icges(5,6) * t201;
t171 = -Icges(5,5) * t257 - Icges(5,6) * t202 + Icges(5,3) * t201;
t170 = rSges(3,1) * t200 - rSges(3,2) * t199 + rSges(3,3) * t259;
t169 = rSges(3,1) * t198 - rSges(3,2) * t197 - rSges(3,3) * t258;
t168 = Icges(3,1) * t200 - Icges(3,4) * t199 + Icges(3,5) * t259;
t167 = Icges(3,1) * t198 - Icges(3,4) * t197 - Icges(3,5) * t258;
t166 = Icges(3,4) * t200 - Icges(3,2) * t199 + Icges(3,6) * t259;
t165 = Icges(3,4) * t198 - Icges(3,2) * t197 - Icges(3,6) * t258;
t164 = Icges(3,5) * t200 - Icges(3,6) * t199 + Icges(3,3) * t259;
t163 = Icges(3,5) * t198 - Icges(3,6) * t197 - Icges(3,3) * t258;
t144 = -t169 * t212 - t196 * t258;
t143 = t170 * t212 - t196 * t259;
t131 = rSges(4,1) * t189 - rSges(4,2) * t188 + rSges(4,3) * t199;
t130 = rSges(4,1) * t187 - rSges(4,2) * t186 + rSges(4,3) * t197;
t128 = rSges(5,1) * t197 - rSges(5,2) * t187 + rSges(5,3) * t186;
t127 = Icges(4,1) * t189 - Icges(4,4) * t188 + Icges(4,5) * t199;
t126 = Icges(4,1) * t187 - Icges(4,4) * t186 + Icges(4,5) * t197;
t125 = Icges(5,1) * t199 - Icges(5,4) * t189 + Icges(5,5) * t188;
t124 = Icges(5,1) * t197 - Icges(5,4) * t187 + Icges(5,5) * t186;
t123 = Icges(4,4) * t189 - Icges(4,2) * t188 + Icges(4,6) * t199;
t122 = Icges(4,4) * t187 - Icges(4,2) * t186 + Icges(4,6) * t197;
t121 = Icges(5,4) * t199 - Icges(5,2) * t189 + Icges(5,6) * t188;
t120 = Icges(5,4) * t197 - Icges(5,2) * t187 + Icges(5,6) * t186;
t119 = Icges(4,5) * t189 - Icges(4,6) * t188 + Icges(4,3) * t199;
t118 = Icges(4,5) * t187 - Icges(4,6) * t186 + Icges(4,3) * t197;
t117 = Icges(5,5) * t199 - Icges(5,6) * t189 + Icges(5,3) * t188;
t116 = Icges(5,5) * t197 - Icges(5,6) * t187 + Icges(5,3) * t186;
t115 = (t169 * t209 + t170 * t211) * t210;
t110 = rSges(6,1) * t157 + rSges(6,2) * t156 + rSges(6,3) * t187;
t96 = -t131 * t257 - t199 * t178;
t95 = t130 * t257 + t197 * t178;
t94 = -t174 * t257 - t201 * t175 + t202 * t176;
t93 = t201 * t171 - t202 * t172 - t173 * t257;
t92 = t130 * t199 - t131 * t197;
t91 = (-t130 - t183) * t212 + t211 * t227;
t90 = t212 * t131 + t209 * t227 + t182;
t89 = t174 * t199 - t175 * t188 + t176 * t189;
t88 = t174 * t197 - t175 * t186 + t176 * t187;
t87 = t171 * t188 - t172 * t189 + t173 * t199;
t86 = t171 * t186 - t172 * t187 + t173 * t197;
t85 = (t130 * t209 + t131 * t211) * t210 + t244;
t84 = t112 * t202 - t140 * t189;
t83 = -t110 * t202 + t140 * t187;
t82 = t245 * t199 + t252 * t257;
t81 = t128 * t257 + t197 * t177 + t249;
t80 = -t119 * t257 - t201 * t123 + t202 * t127;
t79 = -t118 * t257 - t201 * t122 + t202 * t126;
t78 = t201 * t117 - t202 * t121 - t125 * t257;
t77 = t201 * t116 - t202 * t120 - t124 * t257;
t76 = (-t128 + t247) * t212 + t211 * t221;
t75 = t212 * t129 + t209 * t221 + t248;
t72 = t119 * t199 - t123 * t188 + t127 * t189;
t71 = t118 * t199 - t122 * t188 + t126 * t189;
t70 = t119 * t197 - t123 * t186 + t127 * t187;
t69 = t118 * t197 - t122 * t186 + t126 * t187;
t68 = t117 * t188 - t121 * t189 + t125 * t199;
t67 = t116 * t188 - t120 * t189 + t124 * t199;
t66 = t117 * t186 - t121 * t187 + t125 * t197;
t65 = t116 * t186 - t120 * t187 + t124 * t197;
t64 = t110 * t189 - t112 * t187;
t63 = t199 * t128 + t252 * t197 + t132;
t58 = (t128 * t209 + t129 * t211) * t210 + t223;
t57 = t199 * t232 + t233 * t257;
t56 = t110 * t257 + t197 * t140 + t222;
t55 = (-t110 + t230) * t212 + t211 * t220;
t54 = t212 * t112 + t209 * t220 + t231;
t49 = -t250 * t189 + t253 * t202;
t48 = t250 * t187 - t254 * t202;
t39 = t199 * t110 + t197 * t233 + t251;
t38 = (t110 * t209 + t112 * t211) * t210 + t219;
t37 = -t253 * t187 + t254 * t189;
t36 = t199 * t224 + t225 * t257;
t35 = t250 * t197 + t254 * t257 + t222;
t34 = (t230 - t254) * t212 + t211 * t218;
t33 = t209 * t218 + t253 * t212 + t231;
t32 = t94 * t212 + (t209 * t80 - t211 * t79) * t210;
t31 = t93 * t212 + (t209 * t78 - t211 * t77) * t210;
t30 = t79 * t197 + t80 * t199 - t94 * t257;
t29 = t77 * t197 + t78 * t199 - t93 * t257;
t28 = t89 * t212 + (t209 * t72 - t211 * t71) * t210;
t27 = t88 * t212 + (t209 * t70 - t211 * t69) * t210;
t26 = t87 * t212 + (t209 * t68 - t211 * t67) * t210;
t25 = t86 * t212 + (t209 * t66 - t211 * t65) * t210;
t24 = t71 * t197 + t72 * t199 - t89 * t257;
t23 = t69 * t197 + t70 * t199 - t88 * t257;
t22 = t67 * t197 + t68 * t199 - t87 * t257;
t21 = t65 * t197 + t66 * t199 - t86 * t257;
t20 = t225 * t197 + t254 * t199 + t251;
t19 = (t254 * t209 + t253 * t211) * t210 + t219;
t109 = [m(2) + m(3) + m(4) + t267; m(3) * t115 + m(4) * t85 + m(5) * t58 + m(6) * t38 + m(7) * t19; m(7) * (t19 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t38 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t58 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t85 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(3) * (t115 ^ 2 + t143 ^ 2 + t144 ^ 2) + (t11 + t12 + t28 + t26 + (t164 * t259 - t166 * t199 + t168 * t200) * t259) * t259 + (-t9 - t10 - t27 - t25 + (-t163 * t258 - t165 * t197 + t167 * t198) * t258 + (-t163 * t259 + t164 * t258 + t165 * t199 + t166 * t197 - t167 * t200 - t168 * t198) * t259) * t258 + ((t193 * t259 - t199 * t194 + t200 * t195) * t259 - (-t193 * t258 - t197 * t194 + t198 * t195) * t258 + t17 + t18 + t32 + t31 + ((t166 * t217 + t168 * t215) * t209 - (t165 * t217 + t167 * t215) * t211) * t210 ^ 2 + ((-t163 * t211 + t164 * t209 + t194 * t217 + t195 * t215) * t210 + t212 * t193) * t212) * t212; m(4) * t92 + m(5) * t63 + m(6) * t39 + m(7) * t20; (t29 / 0.2e1 + t30 / 0.2e1 + t235) * t212 + (t26 / 0.2e1 + t28 / 0.2e1 + t237) * t199 + (t25 / 0.2e1 + t27 / 0.2e1 + t238) * t197 + m(7) * (t19 * t20 + t33 * t36 + t34 * t35) + m(6) * (t38 * t39 + t54 * t57 + t55 * t56) + m(5) * (t58 * t63 + t75 * t82 + t76 * t81) + m(4) * (t85 * t92 + t90 * t96 + t91 * t95) + ((-t31 / 0.2e1 - t32 / 0.2e1 - t234) * t217 + (-t21 / 0.2e1 - t23 / 0.2e1 - t240) * t211 + (t22 / 0.2e1 + t24 / 0.2e1 + t239) * t209) * t210; (-t15 - t16 - t29 - t30) * t257 + (t8 + t7 + t24 + t22) * t199 + (t5 + t6 + t23 + t21) * t197 + m(7) * (t20 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t39 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t63 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2); t201 * t267; m(7) * (t186 * t33 + t188 * t34 + t19 * t201) + m(6) * (t186 * t54 + t188 * t55 + t201 * t38) + m(5) * (t186 * t75 + t188 * t76 + t201 * t58); m(7) * (t186 * t36 + t188 * t35 + t20 * t201) + m(6) * (t186 * t57 + t188 * t56 + t201 * t39) + m(5) * (t186 * t82 + t188 * t81 + t201 * t63); (t186 ^ 2 + t188 ^ 2 + t201 ^ 2) * t267; m(6) * t64 + m(7) * t37; t236 * t212 + t234 * t202 + t237 * t189 + t238 * t187 + m(7) * (t19 * t37 + t33 * t49 + t34 * t48) + m(6) * (t38 * t64 + t54 * t84 + t55 * t83) + (t241 * t209 + t242 * t211) * t210; -t236 * t257 + t235 * t202 + t241 * t199 - t242 * t197 + t239 * t189 + t240 * t187 + m(7) * (t20 * t37 + t35 * t48 + t36 * t49) + m(6) * (t39 * t64 + t56 * t83 + t57 * t84); m(6) * (t186 * t84 + t188 * t83 + t201 * t64) + m(7) * (t186 * t49 + t188 * t48 + t201 * t37); (t14 + t13) * t202 + (t4 + t3) * t189 + (t2 + t1) * t187 + m(7) * (t37 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t64 ^ 2 + t83 ^ 2 + t84 ^ 2); m(7) * t202; m(7) * (t187 * t33 + t189 * t34 + t19 * t202); m(7) * (t187 * t36 + t189 * t35 + t20 * t202); m(7) * (t186 * t187 + t188 * t189 + t201 * t202); m(7) * (t187 * t49 + t189 * t48 + t202 * t37); m(7) * (t187 ^ 2 + t189 ^ 2 + t202 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t109(1) t109(2) t109(4) t109(7) t109(11) t109(16); t109(2) t109(3) t109(5) t109(8) t109(12) t109(17); t109(4) t109(5) t109(6) t109(9) t109(13) t109(18); t109(7) t109(8) t109(9) t109(10) t109(14) t109(19); t109(11) t109(12) t109(13) t109(14) t109(15) t109(20); t109(16) t109(17) t109(18) t109(19) t109(20) t109(21);];
Mq  = res;
