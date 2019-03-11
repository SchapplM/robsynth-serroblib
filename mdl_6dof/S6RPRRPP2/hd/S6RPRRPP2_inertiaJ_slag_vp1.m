% Calculate joint inertia matrix for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:35
% EndTime: 2019-03-09 04:31:39
% DurationCPUTime: 2.66s
% Computational Cost: add. (5998->383), mult. (7317->542), div. (0->0), fcn. (7975->8), ass. (0->180)
t167 = qJ(1) + pkin(9);
t164 = sin(t167);
t165 = cos(t167);
t172 = cos(qJ(4));
t169 = sin(qJ(4));
t173 = cos(qJ(3));
t207 = t169 * t173;
t126 = -t164 * t172 + t165 * t207;
t205 = t172 * t173;
t127 = t164 * t169 + t165 * t205;
t170 = sin(qJ(3));
t210 = t165 * t170;
t238 = rSges(7,3) + qJ(6);
t239 = rSges(7,1) + pkin(5);
t220 = t126 * rSges(7,2) + t239 * t127 - t210 * t238;
t241 = Icges(4,5) * t170;
t240 = t241 / 0.2e1;
t130 = -Icges(5,3) * t173 + (Icges(5,5) * t172 - Icges(5,6) * t169) * t170;
t132 = -Icges(6,2) * t173 + (Icges(6,4) * t172 + Icges(6,6) * t169) * t170;
t237 = -t130 - t132;
t124 = t164 * t207 + t165 * t172;
t125 = t164 * t205 - t165 * t169;
t211 = t164 * t170;
t63 = Icges(7,5) * t125 + Icges(7,6) * t124 - Icges(7,3) * t211;
t69 = Icges(7,4) * t125 + Icges(7,2) * t124 - Icges(7,6) * t211;
t75 = Icges(7,1) * t125 + Icges(7,4) * t124 - Icges(7,5) * t211;
t17 = t124 * t69 + t125 * t75 - t63 * t211;
t64 = Icges(7,5) * t127 + Icges(7,6) * t126 - Icges(7,3) * t210;
t70 = Icges(7,4) * t127 + Icges(7,2) * t126 - Icges(7,6) * t210;
t76 = Icges(7,1) * t127 + Icges(7,4) * t126 - Icges(7,5) * t210;
t18 = t124 * t70 + t125 * t76 - t64 * t211;
t65 = Icges(6,5) * t125 + Icges(6,6) * t211 + Icges(6,3) * t124;
t71 = Icges(6,4) * t125 + Icges(6,2) * t211 + Icges(6,6) * t124;
t77 = Icges(6,1) * t125 + Icges(6,4) * t211 + Icges(6,5) * t124;
t19 = t124 * t65 + t125 * t77 + t71 * t211;
t66 = Icges(6,5) * t127 + Icges(6,6) * t210 + Icges(6,3) * t126;
t72 = Icges(6,4) * t127 + Icges(6,2) * t210 + Icges(6,6) * t126;
t78 = Icges(6,1) * t127 + Icges(6,4) * t210 + Icges(6,5) * t126;
t20 = t124 * t66 + t125 * t78 + t72 * t211;
t67 = Icges(5,5) * t125 - Icges(5,6) * t124 + Icges(5,3) * t211;
t73 = Icges(5,4) * t125 - Icges(5,2) * t124 + Icges(5,6) * t211;
t79 = Icges(5,1) * t125 - Icges(5,4) * t124 + Icges(5,5) * t211;
t21 = -t124 * t73 + t125 * t79 + t67 * t211;
t68 = Icges(5,5) * t127 - Icges(5,6) * t126 + Icges(5,3) * t210;
t74 = Icges(5,4) * t127 - Icges(5,2) * t126 + Icges(5,6) * t210;
t80 = Icges(5,1) * t127 - Icges(5,4) * t126 + Icges(5,5) * t210;
t22 = -t124 * t74 + t125 * t80 + t68 * t211;
t128 = Icges(7,3) * t173 + (Icges(7,5) * t172 + Icges(7,6) * t169) * t170;
t131 = Icges(7,6) * t173 + (Icges(7,4) * t172 + Icges(7,2) * t169) * t170;
t134 = Icges(7,5) * t173 + (Icges(7,1) * t172 + Icges(7,4) * t169) * t170;
t45 = t124 * t131 + t125 * t134 - t128 * t211;
t129 = -Icges(6,6) * t173 + (Icges(6,5) * t172 + Icges(6,3) * t169) * t170;
t135 = -Icges(6,4) * t173 + (Icges(6,1) * t172 + Icges(6,5) * t169) * t170;
t46 = t124 * t129 + t125 * t135 + t132 * t211;
t133 = -Icges(5,6) * t173 + (Icges(5,4) * t172 - Icges(5,2) * t169) * t170;
t136 = -Icges(5,5) * t173 + (Icges(5,1) * t172 - Icges(5,4) * t169) * t170;
t47 = -t124 * t133 + t125 * t136 + t130 * t211;
t235 = (-t45 - t46 - t47) * t173 + ((t18 + t20 + t22) * t165 + (t17 + t19 + t21) * t164) * t170;
t23 = t126 * t69 + t127 * t75 - t63 * t210;
t24 = t126 * t70 + t127 * t76 - t64 * t210;
t25 = t126 * t65 + t127 * t77 + t71 * t210;
t26 = t126 * t66 + t127 * t78 + t72 * t210;
t27 = -t126 * t73 + t127 * t79 + t67 * t210;
t28 = -t126 * t74 + t127 * t80 + t68 * t210;
t48 = t126 * t131 + t127 * t134 - t128 * t210;
t49 = t126 * t129 + t127 * t135 + t132 * t210;
t50 = -t126 * t133 + t127 * t136 + t130 * t210;
t234 = (-t48 - t49 - t50) * t173 + ((t24 + t26 + t28) * t165 + (t23 + t25 + t27) * t164) * t170;
t31 = t173 * t63 + (t169 * t69 + t172 * t75) * t170;
t33 = -t173 * t71 + (t169 * t65 + t172 * t77) * t170;
t35 = -t173 * t67 + (-t169 * t73 + t172 * t79) * t170;
t233 = -t31 - t33 - t35;
t32 = t173 * t64 + (t169 * t70 + t172 * t76) * t170;
t34 = -t173 * t72 + (t169 * t66 + t172 * t78) * t170;
t36 = -t173 * t68 + (-t169 * t74 + t172 * t80) * t170;
t232 = t32 + t34 + t36;
t206 = t170 * t172;
t208 = t169 * t170;
t231 = t173 * t128 + (t129 + t131) * t208 + (t134 + t135 + t136) * t206;
t162 = t164 ^ 2;
t163 = t165 ^ 2;
t230 = t173 ^ 2;
t229 = m(6) + m(7);
t228 = t164 / 0.2e1;
t226 = -t173 / 0.2e1;
t148 = t170 * rSges(4,1) + t173 * rSges(4,2);
t225 = m(4) * t148;
t224 = m(7) * t170;
t223 = pkin(3) * t173;
t171 = sin(qJ(1));
t222 = t171 * pkin(1);
t218 = t124 * rSges(7,2);
t221 = t125 * t239 - t238 * t211 + t218;
t85 = t127 * rSges(6,1) + rSges(6,2) * t210 + t126 * rSges(6,3);
t93 = t127 * pkin(4) + t126 * qJ(5);
t219 = -t85 - t93;
t217 = t124 * rSges(6,3);
t216 = t165 * rSges(4,3);
t142 = (pkin(4) * t172 + qJ(5) * t169) * t170;
t115 = t124 * qJ(5);
t92 = t125 * pkin(4) + t115;
t214 = t142 * t211 + t173 * t92;
t212 = Icges(4,4) * t173;
t209 = t165 * t173;
t198 = pkin(3) * t209 + pkin(8) * t210;
t203 = t162 * (pkin(8) * t170 + t223) + t165 * t198;
t201 = (rSges(7,1) * t172 + rSges(7,2) * t169) * t170 + pkin(5) * t206 + t238 * t173;
t138 = -t173 * rSges(6,2) + (rSges(6,1) * t172 + rSges(6,3) * t169) * t170;
t200 = -t138 - t142;
t139 = -t173 * rSges(5,3) + (rSges(5,1) * t172 - rSges(5,2) * t169) * t170;
t154 = t170 * pkin(3) - t173 * pkin(8);
t199 = -t139 - t154;
t197 = t162 + t163;
t196 = -t133 * t208 + t173 * t237 + t231;
t195 = -t93 - t220;
t86 = t127 * rSges(5,1) - t126 * rSges(5,2) + rSges(5,3) * t210;
t194 = -t142 - t201;
t193 = -t154 + t200;
t174 = cos(qJ(1));
t166 = t174 * pkin(1);
t192 = t165 * pkin(2) + t164 * pkin(7) + t166;
t191 = -pkin(2) - t223;
t190 = t165 * pkin(7) - t222;
t189 = t164 * t92 + t165 * t93 + t203;
t188 = -t154 + t194;
t187 = -t115 + t190;
t186 = rSges(4,1) * t173 - rSges(4,2) * t170;
t185 = -t125 * rSges(5,1) + t124 * rSges(5,2);
t184 = t192 + t198;
t182 = -Icges(4,2) * t170 + t212;
t181 = Icges(4,5) * t173 - Icges(4,6) * t170;
t178 = rSges(4,1) * t209 - rSges(4,2) * t210 + t164 * rSges(4,3);
t177 = t31 / 0.2e1 + t35 / 0.2e1 + t33 / 0.2e1 + t47 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1;
t176 = t32 / 0.2e1 + t34 / 0.2e1 + t50 / 0.2e1 + t49 / 0.2e1 + t48 / 0.2e1 + t36 / 0.2e1;
t175 = t184 + t93;
t168 = t170 ^ 2;
t150 = t174 * rSges(2,1) - t171 * rSges(2,2);
t149 = -t171 * rSges(2,1) - t174 * rSges(2,2);
t145 = Icges(4,6) * t173 + t241;
t141 = t165 * rSges(3,1) - t164 * rSges(3,2) + t166;
t140 = -t164 * rSges(3,1) - t165 * rSges(3,2) - t222;
t101 = Icges(4,3) * t164 + t181 * t165;
t100 = -Icges(4,3) * t165 + t181 * t164;
t97 = t199 * t165;
t96 = t199 * t164;
t95 = t178 + t192;
t94 = t216 + (-pkin(2) - t186) * t164 + t190;
t87 = t92 * t210;
t83 = rSges(5,3) * t211 - t185;
t82 = t125 * rSges(6,1) + rSges(6,2) * t211 + t217;
t62 = t193 * t165;
t61 = t193 * t164;
t60 = t165 * t178 + (t186 * t164 - t216) * t164;
t56 = t188 * t165;
t55 = t188 * t164;
t54 = -t139 * t210 - t173 * t86;
t53 = t139 * t211 + t173 * t83;
t52 = t184 + t86;
t51 = ((-rSges(5,3) - pkin(8)) * t170 + t191) * t164 + t185 + t190;
t44 = (-t164 * t86 + t165 * t83) * t170;
t43 = t175 + t85;
t42 = -t217 + (-rSges(6,1) - pkin(4)) * t125 + ((-rSges(6,2) - pkin(8)) * t170 + t191) * t164 + t187;
t41 = t219 * t173 + t200 * t210;
t40 = t138 * t211 + t173 * t82 + t214;
t39 = t175 + t220;
t38 = -t218 + (-pkin(4) - t239) * t125 + ((-pkin(8) + t238) * t170 + t191) * t164 + t187;
t37 = t164 * t83 + t165 * t86 + t203;
t30 = t195 * t173 + t194 * t210;
t29 = t221 * t173 + t201 * t211 + t214;
t16 = t87 + (t219 * t164 + t165 * t82) * t170;
t15 = t164 * t82 + t165 * t85 + t189;
t14 = t87 + (t195 * t164 + t221 * t165) * t170;
t13 = t221 * t164 + t220 * t165 + t189;
t12 = t28 * t164 - t27 * t165;
t11 = t26 * t164 - t25 * t165;
t10 = t24 * t164 - t23 * t165;
t9 = t22 * t164 - t21 * t165;
t8 = t20 * t164 - t19 * t165;
t7 = t18 * t164 - t17 * t165;
t1 = [Icges(2,3) + Icges(3,3) + (Icges(4,1) * t170 - t169 * t133 + t212) * t170 + (Icges(4,4) * t170 + Icges(4,2) * t173 + t237) * t173 + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t51 ^ 2 + t52 ^ 2) + m(4) * (t94 ^ 2 + t95 ^ 2) + m(3) * (t140 ^ 2 + t141 ^ 2) + m(2) * (t149 ^ 2 + t150 ^ 2) + t231; 0; m(3) + m(4) + m(5) + t229; m(7) * (t56 * t38 + t55 * t39) + m(6) * (t62 * t42 + t61 * t43) + m(5) * (t97 * t51 + t96 * t52) + (t182 * t164 * t226 - t94 * t225 - t177 + (-Icges(4,6) * t226 + t240 + t145 / 0.2e1) * t165) * t165 + (-t95 * t225 + t173 * (Icges(4,6) * t164 + t182 * t165) / 0.2e1 + t164 * t240 + t145 * t228 + t176) * t164; m(4) * t60 + m(5) * t37 + m(6) * t15 + m(7) * t13; m(7) * (t13 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t15 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t37 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(4) * (t197 * t148 ^ 2 + t60 ^ 2) + (-t163 * t100 - t7 - t8 - t9) * t165 + (t162 * t101 + t10 + t11 + t12 + (-t164 * t100 + t165 * t101) * t165) * t164; -t196 * t173 + m(7) * (t29 * t38 + t30 * t39) + m(6) * (t40 * t42 + t41 * t43) + m(5) * (t53 * t51 + t54 * t52) + (t177 * t164 + t176 * t165) * t170; m(5) * t44 + m(6) * t16 + m(7) * t14; m(7) * (t14 * t13 + t29 * t56 + t30 * t55) + m(6) * (t16 * t15 + t40 * t62 + t41 * t61) + m(5) * (t44 * t37 + t53 * t97 + t54 * t96) + ((t12 / 0.2e1 + t11 / 0.2e1 + t10 / 0.2e1) * t165 + (t7 / 0.2e1 + t9 / 0.2e1 + t8 / 0.2e1) * t164) * t170 + t234 * t228 - t235 * t165 / 0.2e1 + (t232 * t164 + t233 * t165) * t226; m(6) * (t16 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(7) * (t14 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(5) * (t44 ^ 2 + t53 ^ 2 + t54 ^ 2) + t196 * t230 + ((-t232 * t173 + t234) * t165 + (t233 * t173 + t235) * t164) * t170; m(7) * (t124 * t39 + t126 * t38) + m(6) * (t124 * t43 + t126 * t42); t229 * t208; m(7) * (t124 * t55 + t126 * t56 + t13 * t208) + m(6) * (t124 * t61 + t126 * t62 + t15 * t208); m(6) * (t124 * t41 + t126 * t40 + t16 * t208) + m(7) * (t124 * t30 + t126 * t29 + t14 * t208); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t168 * t169 ^ 2 + t124 ^ 2 + t126 ^ 2); (-t164 * t39 - t165 * t38) * t224; m(7) * t173; m(7) * (t173 * t13 + (-t164 * t55 - t165 * t56) * t170); m(7) * (t173 * t14 + (-t164 * t30 - t165 * t29) * t170); (-t124 * t164 - t126 * t165 + t207) * t224; m(7) * (t197 * t168 + t230);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
