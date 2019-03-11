% Calculate joint inertia matrix for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:19
% EndTime: 2019-03-08 19:51:27
% DurationCPUTime: 4.52s
% Computational Cost: add. (9122->466), mult. (23575->683), div. (0->0), fcn. (29515->10), ass. (0->202)
t242 = Icges(3,1) + Icges(4,2);
t240 = Icges(3,4) + Icges(4,6);
t239 = Icges(3,5) - Icges(4,4);
t241 = Icges(3,2) + Icges(4,3);
t238 = Icges(3,6) - Icges(4,5);
t237 = Icges(3,3) + Icges(4,1);
t185 = sin(pkin(10));
t187 = cos(pkin(10));
t190 = sin(qJ(2));
t188 = cos(pkin(6));
t192 = cos(qJ(2));
t211 = t188 * t192;
t174 = t185 * t190 - t187 * t211;
t212 = t188 * t190;
t175 = t185 * t192 + t187 * t212;
t186 = sin(pkin(6));
t215 = t186 * t187;
t236 = t241 * t174 - t240 * t175 + t238 * t215;
t235 = t240 * t174 - t242 * t175 + t239 * t215;
t176 = t185 * t211 + t187 * t190;
t177 = -t185 * t212 + t187 * t192;
t216 = t185 * t186;
t234 = t241 * t176 - t240 * t177 - t238 * t216;
t233 = -t240 * t176 + t242 * t177 + t239 * t216;
t205 = m(6) / 0.2e1 + m(7) / 0.2e1;
t232 = 0.2e1 * t205;
t231 = t174 * t238 - t175 * t239 + t215 * t237;
t230 = -t176 * t238 + t177 * t239 + t216 * t237;
t229 = t237 * t188 + (t190 * t239 + t192 * t238) * t186;
t228 = t238 * t188 + (t240 * t190 + t241 * t192) * t186;
t227 = t239 * t188 + (t242 * t190 + t240 * t192) * t186;
t220 = cos(qJ(4));
t200 = t186 * t220;
t219 = sin(qJ(4));
t155 = t176 * t219 + t185 * t200;
t158 = -t174 * t219 + t187 * t200;
t199 = t186 * t219;
t179 = t188 * t220 - t192 * t199;
t154 = -t176 * t220 + t185 * t199;
t189 = sin(qJ(6));
t191 = cos(qJ(6));
t115 = t154 * t191 - t177 * t189;
t116 = t154 * t189 + t177 * t191;
t156 = t174 * t220 + t187 * t199;
t113 = -t156 * t191 - t175 * t189;
t114 = -t156 * t189 + t175 * t191;
t71 = Icges(7,5) * t114 + Icges(7,6) * t113 - Icges(7,3) * t158;
t73 = Icges(7,4) * t114 + Icges(7,2) * t113 - Icges(7,6) * t158;
t75 = Icges(7,1) * t114 + Icges(7,4) * t113 - Icges(7,5) * t158;
t26 = t115 * t73 + t116 * t75 + t155 * t71;
t72 = Icges(7,5) * t116 + Icges(7,6) * t115 + Icges(7,3) * t155;
t74 = Icges(7,4) * t116 + Icges(7,2) * t115 + Icges(7,6) * t155;
t76 = Icges(7,1) * t116 + Icges(7,4) * t115 + Icges(7,5) * t155;
t27 = t115 * t74 + t116 * t76 + t155 * t72;
t178 = t188 * t219 + t192 * t200;
t214 = t186 * t190;
t159 = t178 * t191 - t189 * t214;
t160 = t178 * t189 + t191 * t214;
t100 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t179;
t101 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t179;
t102 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t179;
t36 = t100 * t155 + t101 * t115 + t102 * t116;
t2 = t155 * t27 - t158 * t26 + t179 * t36;
t225 = t2 / 0.2e1;
t30 = t159 * t73 + t160 * t75 + t179 * t71;
t31 = t159 * t74 + t160 * t76 + t179 * t72;
t50 = t100 * t179 + t101 * t159 + t102 * t160;
t7 = t155 * t31 - t158 * t30 + t179 * t50;
t224 = t7 / 0.2e1;
t223 = t155 / 0.2e1;
t222 = -t158 / 0.2e1;
t221 = t179 / 0.2e1;
t77 = rSges(7,1) * t114 + rSges(7,2) * t113 - rSges(7,3) * t158;
t218 = pkin(5) * t175 - pkin(9) * t158 + t77;
t78 = rSges(7,1) * t116 + rSges(7,2) * t115 + rSges(7,3) * t155;
t217 = pkin(5) * t177 + pkin(9) * t155 + t78;
t213 = t186 * t192;
t103 = rSges(7,1) * t160 + rSges(7,2) * t159 + rSges(7,3) * t179;
t210 = pkin(5) * t214 + pkin(9) * t179 + t103;
t148 = pkin(2) * t175 + qJ(3) * t174;
t149 = pkin(2) * t177 + qJ(3) * t176;
t209 = t148 * t216 + t149 * t215;
t147 = t188 * t149;
t161 = pkin(3) * t216 + pkin(8) * t177;
t208 = t188 * t161 + t147;
t162 = -pkin(3) * t215 + pkin(8) * t175;
t207 = -t148 - t162;
t180 = (pkin(2) * t190 - qJ(3) * t192) * t186;
t206 = -pkin(3) * t188 - pkin(8) * t214 - t180;
t204 = -m(4) - m(5) - m(6) - m(7);
t111 = pkin(4) * t155 + qJ(5) * t154;
t203 = t188 * t111 + t208;
t112 = -pkin(4) * t158 - qJ(5) * t156;
t202 = -t112 + t207;
t150 = pkin(4) * t179 + qJ(5) * t178;
t201 = -t150 + t206;
t198 = (-t188 * rSges(4,1) - (-rSges(4,2) * t190 - rSges(4,3) * t192) * t186 - t180) * t186;
t197 = t161 * t215 + t162 * t216 + t209;
t143 = rSges(5,1) * t179 - rSges(5,2) * t178 + rSges(5,3) * t214;
t196 = (-t143 + t206) * t186;
t142 = rSges(6,1) * t214 - rSges(6,2) * t179 + rSges(6,3) * t178;
t195 = (-t142 + t201) * t186;
t194 = t111 * t215 + t112 * t216 + t197;
t193 = (t201 - t210) * t186;
t170 = t188 * rSges(3,3) + (rSges(3,1) * t190 + rSges(3,2) * t192) * t186;
t141 = Icges(5,1) * t179 - Icges(5,4) * t178 + Icges(5,5) * t214;
t140 = Icges(5,4) * t179 - Icges(5,2) * t178 + Icges(5,6) * t214;
t139 = Icges(5,5) * t179 - Icges(5,6) * t178 + Icges(5,3) * t214;
t138 = Icges(6,1) * t214 - Icges(6,4) * t179 + Icges(6,5) * t178;
t137 = Icges(6,4) * t214 - Icges(6,2) * t179 + Icges(6,6) * t178;
t136 = Icges(6,5) * t214 - Icges(6,6) * t179 + Icges(6,3) * t178;
t135 = rSges(3,1) * t177 - rSges(3,2) * t176 + rSges(3,3) * t216;
t134 = rSges(3,1) * t175 - rSges(3,2) * t174 - rSges(3,3) * t215;
t133 = -rSges(4,1) * t215 - rSges(4,2) * t175 + rSges(4,3) * t174;
t132 = rSges(4,1) * t216 - rSges(4,2) * t177 + rSges(4,3) * t176;
t117 = t175 * t150;
t109 = t111 * t214;
t106 = -t134 * t188 - t170 * t215;
t105 = t135 * t188 - t170 * t216;
t99 = t177 * t112;
t98 = rSges(6,1) * t177 - rSges(6,2) * t155 + rSges(6,3) * t154;
t97 = rSges(6,1) * t175 + rSges(6,2) * t158 - rSges(6,3) * t156;
t96 = -rSges(5,1) * t158 + rSges(5,2) * t156 + rSges(5,3) * t175;
t95 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t177;
t94 = -Icges(5,1) * t158 + Icges(5,4) * t156 + Icges(5,5) * t175;
t93 = Icges(5,1) * t155 - Icges(5,4) * t154 + Icges(5,5) * t177;
t92 = Icges(6,1) * t177 - Icges(6,4) * t155 + Icges(6,5) * t154;
t91 = Icges(6,1) * t175 + Icges(6,4) * t158 - Icges(6,5) * t156;
t90 = -Icges(5,4) * t158 + Icges(5,2) * t156 + Icges(5,6) * t175;
t89 = Icges(5,4) * t155 - Icges(5,2) * t154 + Icges(5,6) * t177;
t88 = Icges(6,4) * t177 - Icges(6,2) * t155 + Icges(6,6) * t154;
t87 = Icges(6,4) * t175 + Icges(6,2) * t158 - Icges(6,6) * t156;
t86 = -Icges(5,5) * t158 + Icges(5,6) * t156 + Icges(5,3) * t175;
t85 = Icges(5,5) * t155 - Icges(5,6) * t154 + Icges(5,3) * t177;
t84 = Icges(6,5) * t177 - Icges(6,6) * t155 + Icges(6,3) * t154;
t83 = Icges(6,5) * t175 + Icges(6,6) * t158 - Icges(6,3) * t156;
t82 = (t134 * t185 + t135 * t187) * t186;
t80 = (-t133 - t148) * t188 + t187 * t198;
t79 = t132 * t188 + t185 * t198 + t147;
t70 = -t143 * t177 + t214 * t95;
t69 = t143 * t175 - t214 * t96;
t68 = t139 * t214 - t140 * t178 + t141 * t179;
t67 = t136 * t178 - t137 * t179 + t138 * t214;
t66 = (t132 * t187 + t133 * t185) * t186 + t209;
t65 = -t175 * t95 + t177 * t96;
t64 = t136 * t154 - t137 * t155 + t138 * t177;
t63 = -t136 * t156 + t137 * t158 + t138 * t175;
t62 = t139 * t175 + t140 * t156 - t141 * t158;
t61 = t139 * t177 - t140 * t154 + t141 * t155;
t60 = (-t96 + t207) * t188 + t187 * t196;
t59 = t185 * t196 + t188 * t95 + t208;
t58 = -t103 * t158 - t179 * t77;
t57 = -t103 * t155 + t179 * t78;
t56 = t98 * t214 + t109 + (-t142 - t150) * t177;
t55 = t142 * t175 + t117 + (-t112 - t97) * t214;
t54 = -t178 * t90 + t179 * t94 + t214 * t86;
t53 = -t178 * t89 + t179 * t93 + t214 * t85;
t52 = t178 * t84 - t179 * t88 + t214 * t92;
t51 = t178 * t83 - t179 * t87 + t214 * t91;
t49 = t154 * t84 - t155 * t88 + t177 * t92;
t48 = t154 * t83 - t155 * t87 + t177 * t91;
t47 = -t156 * t84 + t158 * t88 + t175 * t92;
t46 = -t156 * t83 + t158 * t87 + t175 * t91;
t45 = t156 * t90 - t158 * t94 + t175 * t86;
t44 = t156 * t89 - t158 * t93 + t175 * t85;
t43 = -t154 * t90 + t155 * t94 + t177 * t86;
t42 = -t154 * t89 + t155 * t93 + t177 * t85;
t41 = (t185 * t96 + t187 * t95) * t186 + t197;
t40 = t155 * t77 + t158 * t78;
t39 = (-t97 + t202) * t188 + t187 * t195;
t38 = t185 * t195 + t188 * t98 + t203;
t37 = t177 * t97 + t99 + (-t111 - t98) * t175;
t35 = -t100 * t158 + t101 * t113 + t102 * t114;
t34 = t109 + t217 * t214 + (-t150 - t210) * t177;
t33 = t117 + t210 * t175 + (-t112 - t218) * t214;
t32 = (t185 * t97 + t187 * t98) * t186 + t194;
t29 = (t202 - t218) * t188 + t187 * t193;
t28 = t185 * t193 + t188 * t217 + t203;
t25 = t113 * t74 + t114 * t76 - t158 * t72;
t24 = t113 * t73 + t114 * t75 - t158 * t71;
t23 = t99 + t218 * t177 + (-t111 - t217) * t175;
t22 = (t185 * t218 + t187 * t217) * t186 + t194;
t21 = t188 * t68 + (t185 * t53 - t187 * t54) * t186;
t20 = t188 * t67 + (t185 * t52 - t187 * t51) * t186;
t19 = t175 * t54 + t177 * t53 + t214 * t68;
t18 = t175 * t51 + t177 * t52 + t214 * t67;
t17 = t188 * t64 + (t185 * t49 - t187 * t48) * t186;
t16 = t188 * t63 + (t185 * t47 - t187 * t46) * t186;
t15 = t188 * t62 + (t185 * t44 - t187 * t45) * t186;
t14 = t188 * t61 + (t185 * t42 - t187 * t43) * t186;
t13 = t175 * t48 + t177 * t49 + t214 * t64;
t12 = t175 * t46 + t177 * t47 + t214 * t63;
t11 = t175 * t45 + t177 * t44 + t214 * t62;
t10 = t175 * t43 + t177 * t42 + t214 * t61;
t9 = t188 * t50 + (t185 * t31 - t187 * t30) * t186;
t8 = t175 * t30 + t177 * t31 + t214 * t50;
t6 = t188 * t36 + (t185 * t27 - t187 * t26) * t186;
t5 = t188 * t35 + (t185 * t25 - t187 * t24) * t186;
t4 = t175 * t26 + t177 * t27 + t214 * t36;
t3 = t175 * t24 + t177 * t25 + t214 * t35;
t1 = t155 * t25 - t158 * t24 + t179 * t35;
t81 = [m(2) + m(3) - t204; m(3) * t82 + m(4) * t66 + m(5) * t41 + m(6) * t32 + m(7) * t22; m(7) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t41 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t66 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(3) * (t105 ^ 2 + t106 ^ 2 + t82 ^ 2) + (t20 + t21 + t9 + t229 * t188 ^ 2 + ((t227 * t190 + t228 * t192) * t188 + (t231 * t188 + (t235 * t190 + t236 * t192) * t186) * t187 + (t230 * t188 + (t233 * t190 - t234 * t192) * t186) * t185) * t186) * t188 + (t6 + t17 + t14 + (t234 * t176 + t233 * t177 + t230 * t216) * t216 + (-t228 * t176 + t227 * t177 + t229 * t216) * t188) * t216 + (-t5 - t15 - t16 + (t236 * t174 - t235 * t175 + t231 * t215) * t215 + (t228 * t174 - t227 * t175 + t229 * t215) * t188 + (-t234 * t174 - t233 * t175 - t236 * t176 + t235 * t177 + t230 * t215 + t231 * t216) * t216) * t215; t204 * t213; m(7) * (t174 * t28 + t176 * t29 - t213 * t22) + m(6) * (t174 * t38 + t176 * t39 - t213 * t32) + m(5) * (t174 * t59 + t176 * t60 - t213 * t41) + m(4) * (t174 * t79 + t176 * t80 - t213 * t66); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t205) * (t186 ^ 2 * t192 ^ 2 + t174 ^ 2 + t176 ^ 2); m(5) * t65 + m(6) * t37 + m(7) * t23; (t8 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t188 + (t6 / 0.2e1 + t17 / 0.2e1 + t14 / 0.2e1) * t177 + (t5 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1) * t175 + m(7) * (t22 * t23 + t28 * t34 + t29 * t33) + m(6) * (t32 * t37 + t38 * t56 + t39 * t55) + m(5) * (t41 * t65 + t59 * t70 + t60 * t69) + ((t9 / 0.2e1 + t21 / 0.2e1 + t20 / 0.2e1) * t190 + (-t3 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1) * t187 + (t4 / 0.2e1 + t10 / 0.2e1 + t13 / 0.2e1) * t185) * t186; m(5) * (t70 * t174 + t69 * t176 - t213 * t65) + m(6) * (t56 * t174 + t55 * t176 - t213 * t37) + m(7) * (t34 * t174 + t33 * t176 - t213 * t23); (t18 + t19 + t8) * t214 + (t4 + t13 + t10) * t177 + (t3 + t12 + t11) * t175 + m(7) * (t23 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t37 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t65 ^ 2 + t69 ^ 2 + t70 ^ 2); t178 * t232; m(7) * (t154 * t29 - t156 * t28 + t178 * t22) + m(6) * (t154 * t39 - t156 * t38 + t178 * t32); (t154 * t176 - t156 * t174 - t178 * t213) * t232; m(7) * (t154 * t33 - t156 * t34 + t178 * t23) + m(6) * (t154 * t55 - t156 * t56 + t178 * t37); (t154 ^ 2 + t156 ^ 2 + t178 ^ 2) * t232; m(7) * t40; m(7) * (t22 * t40 + t28 * t57 + t29 * t58) + t6 * t223 + t188 * t224 + t9 * t221 + t5 * t222 + (t185 * t225 - t187 * t1 / 0.2e1) * t186; m(7) * (t57 * t174 + t58 * t176 - t213 * t40); t8 * t221 + t177 * t225 + t175 * t1 / 0.2e1 + t3 * t222 + t4 * t223 + m(7) * (t23 * t40 + t33 * t58 + t34 * t57) + t214 * t224; m(7) * (t154 * t58 - t156 * t57 + t178 * t40); t155 * t2 - t158 * t1 + t179 * t7 + m(7) * (t40 ^ 2 + t57 ^ 2 + t58 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t81(1) t81(2) t81(4) t81(7) t81(11) t81(16); t81(2) t81(3) t81(5) t81(8) t81(12) t81(17); t81(4) t81(5) t81(6) t81(9) t81(13) t81(18); t81(7) t81(8) t81(9) t81(10) t81(14) t81(19); t81(11) t81(12) t81(13) t81(14) t81(15) t81(20); t81(16) t81(17) t81(18) t81(19) t81(20) t81(21);];
Mq  = res;
