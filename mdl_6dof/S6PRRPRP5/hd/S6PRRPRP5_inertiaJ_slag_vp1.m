% Calculate joint inertia matrix for
% S6PRRPRP5
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:35
% EndTime: 2019-03-08 21:45:44
% DurationCPUTime: 4.39s
% Computational Cost: add. (14426->548), mult. (37287->784), div. (0->0), fcn. (47509->10), ass. (0->238)
t265 = rSges(7,1) + pkin(5);
t264 = rSges(7,3) + qJ(6);
t263 = m(5) + m(6) + m(7);
t262 = cos(qJ(3));
t261 = cos(qJ(5));
t260 = sin(qJ(3));
t211 = sin(pkin(10));
t212 = sin(pkin(6));
t259 = t211 * t212;
t213 = cos(pkin(10));
t258 = t212 * t213;
t217 = cos(qJ(2));
t257 = t212 * t217;
t214 = cos(pkin(6));
t216 = sin(qJ(2));
t256 = t214 * t216;
t255 = t214 * t217;
t201 = t211 * t217 + t213 * t256;
t229 = t212 * t262;
t188 = t201 * t260 + t213 * t229;
t200 = t211 * t216 - t213 * t255;
t215 = sin(qJ(5));
t156 = -t188 * t261 + t200 * t215;
t157 = t188 * t215 + t200 * t261;
t228 = t212 * t260;
t189 = t201 * t262 - t213 * t228;
t254 = t189 * rSges(7,2) + t264 * t156 + t157 * t265;
t203 = -t211 * t256 + t213 * t217;
t190 = t203 * t260 - t211 * t229;
t202 = t211 * t255 + t213 * t216;
t158 = -t190 * t261 + t202 * t215;
t159 = t190 * t215 + t202 * t261;
t191 = t203 * t262 + t211 * t228;
t253 = t191 * rSges(7,2) + t264 * t158 + t159 * t265;
t129 = t202 * rSges(5,1) - t191 * rSges(5,2) + t190 * rSges(5,3);
t150 = t191 * pkin(3) + t190 * qJ(4);
t252 = -t129 - t150;
t149 = t189 * pkin(3) + t188 * qJ(4);
t132 = t202 * t149;
t161 = t200 * pkin(4) + t189 * pkin(9);
t251 = t202 * t161 + t132;
t204 = -t214 * t262 + t216 * t228;
t192 = t204 * t261 + t215 * t257;
t194 = t204 * t215 - t257 * t261;
t205 = t214 * t260 + t216 * t229;
t250 = t205 * rSges(7,2) - t264 * t192 + t194 * t265;
t185 = t205 * pkin(3) + t204 * qJ(4);
t249 = t149 * t257 + t200 * t185;
t184 = t203 * pkin(2) + t202 * pkin(8);
t182 = t214 * t184;
t248 = t214 * t150 + t182;
t183 = t201 * pkin(2) + t200 * pkin(8);
t247 = -t149 - t183;
t162 = t202 * pkin(4) + t191 * pkin(9);
t246 = -t150 - t162;
t177 = -rSges(5,1) * t257 - t205 * rSges(5,2) + t204 * rSges(5,3);
t245 = -t177 - t185;
t244 = t183 * t259 + t184 * t258;
t195 = -pkin(4) * t257 + t205 * pkin(9);
t243 = -t185 - t195;
t101 = Icges(7,4) * t157 + Icges(7,2) * t189 + Icges(7,6) * t156;
t105 = Icges(7,1) * t157 + Icges(7,4) * t189 + Icges(7,5) * t156;
t97 = Icges(7,5) * t157 + Icges(7,6) * t189 + Icges(7,3) * t156;
t40 = t189 * t101 + t157 * t105 + t156 * t97;
t102 = Icges(7,4) * t159 + Icges(7,2) * t191 + Icges(7,6) * t158;
t106 = Icges(7,1) * t159 + Icges(7,4) * t191 + Icges(7,5) * t158;
t98 = Icges(7,5) * t159 + Icges(7,6) * t191 + Icges(7,3) * t158;
t41 = t189 * t102 + t157 * t106 + t156 * t98;
t133 = Icges(7,5) * t194 + Icges(7,6) * t205 - Icges(7,3) * t192;
t135 = Icges(7,4) * t194 + Icges(7,2) * t205 - Icges(7,6) * t192;
t137 = Icges(7,1) * t194 + Icges(7,4) * t205 - Icges(7,5) * t192;
t59 = t156 * t133 + t189 * t135 + t157 * t137;
t1 = t40 * t189 + t41 * t191 + t59 * t205;
t103 = Icges(6,4) * t157 - Icges(6,2) * t156 + Icges(6,6) * t189;
t107 = Icges(6,1) * t157 - Icges(6,4) * t156 + Icges(6,5) * t189;
t99 = Icges(6,5) * t157 - Icges(6,6) * t156 + Icges(6,3) * t189;
t42 = -t156 * t103 + t157 * t107 + t189 * t99;
t100 = Icges(6,5) * t159 - Icges(6,6) * t158 + Icges(6,3) * t191;
t104 = Icges(6,4) * t159 - Icges(6,2) * t158 + Icges(6,6) * t191;
t108 = Icges(6,1) * t159 - Icges(6,4) * t158 + Icges(6,5) * t191;
t43 = t189 * t100 - t156 * t104 + t157 * t108;
t134 = Icges(6,5) * t194 + Icges(6,6) * t192 + Icges(6,3) * t205;
t136 = Icges(6,4) * t194 + Icges(6,2) * t192 + Icges(6,6) * t205;
t138 = Icges(6,1) * t194 + Icges(6,4) * t192 + Icges(6,5) * t205;
t60 = t189 * t134 - t156 * t136 + t157 * t138;
t2 = t42 * t189 + t43 * t191 + t60 * t205;
t242 = t2 / 0.2e1 + t1 / 0.2e1;
t44 = t191 * t101 + t159 * t105 + t158 * t97;
t45 = t191 * t102 + t159 * t106 + t158 * t98;
t61 = t158 * t133 + t191 * t135 + t159 * t137;
t3 = t44 * t189 + t45 * t191 + t61 * t205;
t46 = -t158 * t103 + t159 * t107 + t191 * t99;
t47 = t191 * t100 - t158 * t104 + t159 * t108;
t62 = t191 * t134 - t158 * t136 + t159 * t138;
t4 = t46 * t189 + t47 * t191 + t62 * t205;
t241 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t40 * t200 + t41 * t202 - t257 * t59;
t6 = t42 * t200 + t43 * t202 - t257 * t60;
t240 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t44 * t200 + t45 * t202 - t257 * t61;
t8 = t46 * t200 + t47 * t202 - t257 * t62;
t239 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t60 * t214 + (t211 * t43 - t213 * t42) * t212;
t9 = t59 * t214 + (t211 * t41 - t213 * t40) * t212;
t238 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t61 * t214 + (t211 * t45 - t213 * t44) * t212;
t12 = t62 * t214 + (t211 * t47 - t213 * t46) * t212;
t237 = t12 / 0.2e1 + t11 / 0.2e1;
t48 = t205 * t101 + t194 * t105 - t192 * t97;
t49 = t205 * t102 + t194 * t106 - t192 * t98;
t73 = -t192 * t133 + t205 * t135 + t194 * t137;
t13 = t48 * t189 + t49 * t191 + t73 * t205;
t50 = t192 * t103 + t194 * t107 + t205 * t99;
t51 = t205 * t100 + t192 * t104 + t194 * t108;
t74 = t205 * t134 + t192 * t136 + t194 * t138;
t14 = t50 * t189 + t51 * t191 + t74 * t205;
t236 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t48 * t200 + t49 * t202 - t257 * t73;
t16 = t50 * t200 + t51 * t202 - t257 * t74;
t235 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t73 * t214 + (t211 * t49 - t213 * t48) * t212;
t18 = t74 * t214 + (t211 * t51 - t213 * t50) * t212;
t234 = t18 / 0.2e1 + t17 / 0.2e1;
t112 = t159 * rSges(6,1) - t158 * rSges(6,2) + t191 * rSges(6,3);
t233 = -t112 + t246;
t140 = t194 * rSges(6,1) + t192 * rSges(6,2) + t205 * rSges(6,3);
t232 = -t140 + t243;
t231 = t214 * t162 + t248;
t230 = -t161 + t247;
t178 = t205 * rSges(4,1) - t204 * rSges(4,2) - rSges(4,3) * t257;
t206 = (pkin(2) * t216 - pkin(8) * t217) * t212;
t227 = (-t178 - t206) * t212;
t225 = t246 - t253;
t224 = t243 - t250;
t223 = t149 * t259 + t150 * t258 + t244;
t222 = t161 * t257 + t200 * t195 + t249;
t221 = (-t206 + t245) * t212;
t220 = (-t206 + t232) * t212;
t219 = t161 * t259 + t162 * t258 + t223;
t218 = (-t206 + t224) * t212;
t199 = t214 * rSges(3,3) + (rSges(3,1) * t216 + rSges(3,2) * t217) * t212;
t198 = Icges(3,5) * t214 + (Icges(3,1) * t216 + Icges(3,4) * t217) * t212;
t197 = Icges(3,6) * t214 + (Icges(3,4) * t216 + Icges(3,2) * t217) * t212;
t196 = Icges(3,3) * t214 + (Icges(3,5) * t216 + Icges(3,6) * t217) * t212;
t176 = Icges(4,1) * t205 - Icges(4,4) * t204 - Icges(4,5) * t257;
t175 = Icges(4,4) * t205 - Icges(4,2) * t204 - Icges(4,6) * t257;
t174 = Icges(4,5) * t205 - Icges(4,6) * t204 - Icges(4,3) * t257;
t173 = -Icges(5,1) * t257 - Icges(5,4) * t205 + Icges(5,5) * t204;
t172 = -Icges(5,4) * t257 - Icges(5,2) * t205 + Icges(5,6) * t204;
t171 = -Icges(5,5) * t257 - Icges(5,6) * t205 + Icges(5,3) * t204;
t170 = t203 * rSges(3,1) - t202 * rSges(3,2) + rSges(3,3) * t259;
t169 = t201 * rSges(3,1) - t200 * rSges(3,2) - rSges(3,3) * t258;
t168 = Icges(3,1) * t203 - Icges(3,4) * t202 + Icges(3,5) * t259;
t167 = Icges(3,1) * t201 - Icges(3,4) * t200 - Icges(3,5) * t258;
t166 = Icges(3,4) * t203 - Icges(3,2) * t202 + Icges(3,6) * t259;
t165 = Icges(3,4) * t201 - Icges(3,2) * t200 - Icges(3,6) * t258;
t164 = Icges(3,5) * t203 - Icges(3,6) * t202 + Icges(3,3) * t259;
t163 = Icges(3,5) * t201 - Icges(3,6) * t200 - Icges(3,3) * t258;
t143 = -t214 * t169 - t199 * t258;
t142 = t214 * t170 - t199 * t259;
t131 = t191 * rSges(4,1) - t190 * rSges(4,2) + t202 * rSges(4,3);
t130 = t189 * rSges(4,1) - t188 * rSges(4,2) + t200 * rSges(4,3);
t128 = t200 * rSges(5,1) - t189 * rSges(5,2) + t188 * rSges(5,3);
t127 = Icges(4,1) * t191 - Icges(4,4) * t190 + Icges(4,5) * t202;
t126 = Icges(4,1) * t189 - Icges(4,4) * t188 + Icges(4,5) * t200;
t125 = Icges(5,1) * t202 - Icges(5,4) * t191 + Icges(5,5) * t190;
t124 = Icges(5,1) * t200 - Icges(5,4) * t189 + Icges(5,5) * t188;
t123 = Icges(4,4) * t191 - Icges(4,2) * t190 + Icges(4,6) * t202;
t122 = Icges(4,4) * t189 - Icges(4,2) * t188 + Icges(4,6) * t200;
t121 = Icges(5,4) * t202 - Icges(5,2) * t191 + Icges(5,6) * t190;
t120 = Icges(5,4) * t200 - Icges(5,2) * t189 + Icges(5,6) * t188;
t119 = Icges(4,5) * t191 - Icges(4,6) * t190 + Icges(4,3) * t202;
t118 = Icges(4,5) * t189 - Icges(4,6) * t188 + Icges(4,3) * t200;
t117 = Icges(5,5) * t202 - Icges(5,6) * t191 + Icges(5,3) * t190;
t116 = Icges(5,5) * t200 - Icges(5,6) * t189 + Icges(5,3) * t188;
t115 = (t169 * t211 + t170 * t213) * t212;
t110 = t157 * rSges(6,1) - t156 * rSges(6,2) + t189 * rSges(6,3);
t96 = -t131 * t257 - t202 * t178;
t95 = t130 * t257 + t200 * t178;
t94 = -t174 * t257 - t204 * t175 + t205 * t176;
t93 = t204 * t171 - t205 * t172 - t173 * t257;
t92 = t202 * t130 - t200 * t131;
t91 = (-t130 - t183) * t214 + t213 * t227;
t90 = t214 * t131 + t211 * t227 + t182;
t89 = t202 * t174 - t190 * t175 + t191 * t176;
t88 = t200 * t174 - t188 * t175 + t189 * t176;
t87 = t190 * t171 - t191 * t172 + t202 * t173;
t86 = t188 * t171 - t189 * t172 + t200 * t173;
t85 = (t130 * t211 + t131 * t213) * t212 + t244;
t84 = t205 * t112 - t191 * t140;
t83 = -t205 * t110 + t189 * t140;
t82 = t202 * t245 + t252 * t257;
t81 = t128 * t257 + t200 * t177 + t249;
t80 = -t119 * t257 - t204 * t123 + t205 * t127;
t79 = -t118 * t257 - t204 * t122 + t205 * t126;
t78 = t204 * t117 - t205 * t121 - t125 * t257;
t77 = t204 * t116 - t205 * t120 - t124 * t257;
t76 = (-t128 + t247) * t214 + t213 * t221;
t75 = t214 * t129 + t211 * t221 + t248;
t72 = t202 * t119 - t190 * t123 + t191 * t127;
t71 = t202 * t118 - t190 * t122 + t191 * t126;
t70 = t200 * t119 - t188 * t123 + t189 * t127;
t69 = t200 * t118 - t188 * t122 + t189 * t126;
t68 = t190 * t117 - t191 * t121 + t202 * t125;
t67 = t190 * t116 - t191 * t120 + t202 * t124;
t66 = t188 * t117 - t189 * t121 + t200 * t125;
t65 = t188 * t116 - t189 * t120 + t200 * t124;
t64 = t191 * t110 - t189 * t112;
t63 = t202 * t128 + t200 * t252 + t132;
t58 = (t128 * t211 + t129 * t213) * t212 + t223;
t57 = t202 * t232 + t233 * t257;
t56 = t110 * t257 + t200 * t140 + t222;
t55 = -t191 * t250 + t205 * t253;
t54 = t189 * t250 - t205 * t254;
t53 = (-t110 + t230) * t214 + t213 * t220;
t52 = t214 * t112 + t211 * t220 + t231;
t39 = t202 * t110 + t200 * t233 + t251;
t38 = -t189 * t253 + t191 * t254;
t37 = (t110 * t211 + t112 * t213) * t212 + t219;
t36 = t202 * t224 + t225 * t257;
t35 = t200 * t250 + t254 * t257 + t222;
t34 = (t230 - t254) * t214 + t213 * t218;
t33 = t211 * t218 + t214 * t253 + t231;
t32 = t94 * t214 + (t211 * t80 - t213 * t79) * t212;
t31 = t93 * t214 + (t211 * t78 - t213 * t77) * t212;
t30 = t79 * t200 + t80 * t202 - t257 * t94;
t29 = t77 * t200 + t78 * t202 - t257 * t93;
t28 = t200 * t225 + t202 * t254 + t251;
t27 = (t211 * t254 + t213 * t253) * t212 + t219;
t26 = t89 * t214 + (t211 * t72 - t213 * t71) * t212;
t25 = t88 * t214 + (t211 * t70 - t213 * t69) * t212;
t24 = t87 * t214 + (t211 * t68 - t213 * t67) * t212;
t23 = t86 * t214 + (t211 * t66 - t213 * t65) * t212;
t22 = t71 * t200 + t72 * t202 - t257 * t89;
t21 = t69 * t200 + t70 * t202 - t257 * t88;
t20 = t67 * t200 + t68 * t202 - t257 * t87;
t19 = t65 * t200 + t66 * t202 - t257 * t86;
t109 = [m(2) + m(3) + m(4) + t263; m(3) * t115 + m(4) * t85 + m(5) * t58 + m(6) * t37 + m(7) * t27; m(7) * (t27 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t37 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t58 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t85 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(3) * (t115 ^ 2 + t142 ^ 2 + t143 ^ 2) + (t12 + t11 + t26 + t24 + (t164 * t259 - t166 * t202 + t168 * t203) * t259) * t259 + (-t9 - t10 - t25 - t23 + (-t163 * t258 - t165 * t200 + t167 * t201) * t258 + (-t163 * t259 + t164 * t258 + t165 * t202 + t200 * t166 - t167 * t203 - t201 * t168) * t259) * t258 + (-(-t196 * t258 - t200 * t197 + t201 * t198) * t258 + (t196 * t259 - t202 * t197 + t203 * t198) * t259 + t18 + t17 + t32 + t31 + ((t166 * t217 + t168 * t216) * t211 - (t165 * t217 + t167 * t216) * t213) * t212 ^ 2 + ((-t163 * t213 + t164 * t211 + t197 * t217 + t198 * t216) * t212 + t214 * t196) * t214) * t214; m(4) * t92 + m(5) * t63 + m(6) * t39 + m(7) * t28; (t29 / 0.2e1 + t30 / 0.2e1 + t235) * t214 + (t24 / 0.2e1 + t26 / 0.2e1 + t237) * t202 + (t23 / 0.2e1 + t25 / 0.2e1 + t238) * t200 + m(7) * (t27 * t28 + t33 * t36 + t34 * t35) + m(6) * (t37 * t39 + t52 * t57 + t53 * t56) + m(5) * (t58 * t63 + t75 * t82 + t76 * t81) + m(4) * (t85 * t92 + t90 * t96 + t91 * t95) + ((-t31 / 0.2e1 - t32 / 0.2e1 - t234) * t217 + (-t21 / 0.2e1 - t19 / 0.2e1 - t240) * t213 + (t20 / 0.2e1 + t22 / 0.2e1 + t239) * t211) * t212; (-t15 - t16 - t29 - t30) * t257 + (t8 + t7 + t22 + t20) * t202 + (t6 + t5 + t21 + t19) * t200 + m(7) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t39 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t63 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2); t204 * t263; m(7) * (t188 * t33 + t190 * t34 + t204 * t27) + m(6) * (t188 * t52 + t190 * t53 + t204 * t37) + m(5) * (t188 * t75 + t190 * t76 + t204 * t58); m(7) * (t188 * t36 + t190 * t35 + t204 * t28) + m(6) * (t188 * t57 + t190 * t56 + t204 * t39) + m(5) * (t188 * t82 + t190 * t81 + t204 * t63); (t188 ^ 2 + t190 ^ 2 + t204 ^ 2) * t263; m(6) * t64 + m(7) * t38; t236 * t214 + t234 * t205 + t237 * t191 + t238 * t189 + m(7) * (t27 * t38 + t33 * t55 + t34 * t54) + m(6) * (t37 * t64 + t52 * t84 + t53 * t83) + (t211 * t241 - t213 * t242) * t212; -t236 * t257 + t235 * t205 + t241 * t202 + t242 * t200 + t239 * t191 + t240 * t189 + m(7) * (t28 * t38 + t35 * t54 + t36 * t55) + m(6) * (t39 * t64 + t56 * t83 + t57 * t84); m(6) * (t188 * t84 + t190 * t83 + t204 * t64) + m(7) * (t188 * t55 + t190 * t54 + t204 * t38); (t13 + t14) * t205 + (t4 + t3) * t191 + (t2 + t1) * t189 + m(7) * (t38 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t64 ^ 2 + t83 ^ 2 + t84 ^ 2); -m(7) * t192; m(7) * (t156 * t33 + t158 * t34 - t192 * t27); m(7) * (t156 * t36 + t158 * t35 - t192 * t28); m(7) * (t156 * t188 + t158 * t190 - t192 * t204); m(7) * (t156 * t55 + t158 * t54 - t192 * t38); m(7) * (t156 ^ 2 + t158 ^ 2 + t192 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t109(1) t109(2) t109(4) t109(7) t109(11) t109(16); t109(2) t109(3) t109(5) t109(8) t109(12) t109(17); t109(4) t109(5) t109(6) t109(9) t109(13) t109(18); t109(7) t109(8) t109(9) t109(10) t109(14) t109(19); t109(11) t109(12) t109(13) t109(14) t109(15) t109(20); t109(16) t109(17) t109(18) t109(19) t109(20) t109(21);];
Mq  = res;
