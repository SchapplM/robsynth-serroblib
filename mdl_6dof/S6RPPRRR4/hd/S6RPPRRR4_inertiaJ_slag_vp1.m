% Calculate joint inertia matrix for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:30
% EndTime: 2019-03-09 02:25:34
% DurationCPUTime: 1.85s
% Computational Cost: add. (6443->340), mult. (11281->499), div. (0->0), fcn. (14214->10), ass. (0->172)
t155 = sin(qJ(4));
t222 = -t155 / 0.2e1;
t158 = cos(qJ(4));
t221 = -t158 / 0.2e1;
t156 = sin(qJ(1));
t159 = cos(qJ(1));
t206 = sin(pkin(10));
t207 = cos(pkin(10));
t131 = -t156 * t206 - t159 * t207;
t220 = t131 ^ 2;
t132 = -t156 * t207 + t159 * t206;
t219 = t132 ^ 2;
t218 = -t131 / 0.2e1;
t217 = t132 / 0.2e1;
t216 = Icges(5,5) * t222 + Icges(5,6) * t221;
t215 = t158 / 0.2e1;
t138 = -rSges(5,1) * t155 - rSges(5,2) * t158;
t214 = m(5) * t138;
t213 = pkin(4) * t158;
t157 = cos(qJ(5));
t146 = pkin(5) * t157 + pkin(4);
t212 = pkin(4) - t146;
t160 = -pkin(9) - pkin(8);
t211 = pkin(8) + t160;
t202 = t131 * t155;
t153 = qJ(5) + qJ(6);
t148 = cos(t153);
t147 = sin(t153);
t197 = t147 * t158;
t95 = t131 * t197 + t132 * t148;
t195 = t148 * t158;
t96 = -t131 * t195 + t132 * t147;
t64 = t96 * rSges(7,1) + t95 * rSges(7,2) - rSges(7,3) * t202;
t198 = t146 * t158;
t154 = sin(qJ(5));
t200 = t132 * t154;
t164 = pkin(5) * t200 - t131 * t198 + t160 * t202;
t201 = t131 * t158;
t188 = pkin(4) * t201 + pkin(8) * t202;
t75 = t164 + t188;
t210 = t64 + t75;
t209 = t131 * rSges(5,3);
t208 = rSges(7,3) - t160;
t205 = Icges(5,4) * t155;
t204 = Icges(5,4) * t158;
t203 = t131 * t154;
t199 = t132 * t155;
t111 = Icges(7,5) * t158 + (-Icges(7,1) * t148 + Icges(7,4) * t147) * t155;
t196 = t148 * t111;
t194 = t154 * t158;
t118 = Icges(6,5) * t158 + (-Icges(6,1) * t157 + Icges(6,4) * t154) * t155;
t193 = t157 * t118;
t192 = t157 * t158;
t109 = Icges(7,3) * t158 + (-Icges(7,5) * t148 + Icges(7,6) * t147) * t155;
t110 = Icges(7,6) * t158 + (-Icges(7,4) * t148 + Icges(7,2) * t147) * t155;
t191 = t155 * t147 * t110 + t158 * t109;
t116 = Icges(6,3) * t158 + (-Icges(6,5) * t157 + Icges(6,6) * t154) * t155;
t117 = Icges(6,6) * t158 + (-Icges(6,4) * t157 + Icges(6,2) * t154) * t155;
t190 = t155 * t154 * t117 + t158 * t116;
t108 = t212 * t155 - t211 * t158;
t113 = t158 * rSges(7,3) + (-rSges(7,1) * t148 + rSges(7,2) * t147) * t155;
t189 = t108 + t113;
t187 = t159 * pkin(1) + t156 * qJ(2);
t93 = -t131 * t148 + t132 * t197;
t94 = -t131 * t147 - t132 * t195;
t58 = Icges(7,5) * t94 + Icges(7,6) * t93 - Icges(7,3) * t199;
t60 = Icges(7,4) * t94 + Icges(7,2) * t93 - Icges(7,6) * t199;
t62 = Icges(7,1) * t94 + Icges(7,4) * t93 - Icges(7,5) * t199;
t30 = t158 * t58 + (t147 * t60 - t148 * t62) * t155;
t59 = Icges(7,5) * t96 + Icges(7,6) * t95 - Icges(7,3) * t202;
t61 = Icges(7,4) * t96 + Icges(7,2) * t95 - Icges(7,6) * t202;
t63 = Icges(7,1) * t96 + Icges(7,4) * t95 - Icges(7,5) * t202;
t31 = t158 * t59 + (t147 * t61 - t148 * t63) * t155;
t20 = -t58 * t199 + t60 * t93 + t62 * t94;
t21 = -t59 * t199 + t61 * t93 + t63 * t94;
t46 = -t109 * t199 + t110 * t93 + t111 * t94;
t5 = t158 * t46 + (-t131 * t21 - t132 * t20) * t155;
t22 = -t58 * t202 + t60 * t95 + t62 * t96;
t23 = -t59 * t202 + t61 * t95 + t63 * t96;
t47 = -t109 * t202 + t110 * t95 + t111 * t96;
t6 = t158 * t47 + (-t131 * t23 - t132 * t22) * t155;
t78 = (-t155 * t196 + t191) * t158;
t186 = -t5 * t199 - t6 * t202 + t158 * (t78 + (-t131 * t31 - t132 * t30) * t155);
t185 = -pkin(8) + t208;
t184 = pkin(5) * t203;
t102 = t131 * t194 + t132 * t157;
t103 = -t131 * t192 + t200;
t74 = t103 * rSges(6,1) + t102 * rSges(6,2) - rSges(6,3) * t202;
t100 = -t131 * t157 + t132 * t194;
t101 = -t132 * t192 - t203;
t67 = Icges(6,5) * t101 + Icges(6,6) * t100 - Icges(6,3) * t199;
t69 = Icges(6,4) * t101 + Icges(6,2) * t100 - Icges(6,6) * t199;
t71 = Icges(6,1) * t101 + Icges(6,4) * t100 - Icges(6,5) * t199;
t33 = t158 * t67 + (t154 * t69 - t157 * t71) * t155;
t48 = t100 * t117 + t101 * t118 - t116 * t199;
t183 = -t48 / 0.2e1 - t33 / 0.2e1;
t68 = Icges(6,5) * t103 + Icges(6,6) * t102 - Icges(6,3) * t202;
t70 = Icges(6,4) * t103 + Icges(6,2) * t102 - Icges(6,6) * t202;
t72 = Icges(6,1) * t103 + Icges(6,4) * t102 - Icges(6,5) * t202;
t34 = t158 * t68 + (t154 * t70 - t157 * t72) * t155;
t49 = t102 * t117 + t103 * t118 - t116 * t202;
t182 = -t49 / 0.2e1 - t34 / 0.2e1;
t181 = t159 * pkin(2) + t187;
t180 = -t202 / 0.2e1;
t179 = -t199 / 0.2e1;
t178 = t212 * t158;
t11 = -t131 * t20 + t132 * t21;
t12 = -t131 * t22 + t132 * t23;
t177 = t11 * t179 + t12 * t180 + t6 * t217 + t5 * t218 + (-t30 * t131 + t31 * t132) * t215;
t176 = t78 + (t31 + t47) * t180 + (t30 + t46) * t179;
t150 = t159 * qJ(2);
t175 = t150 + (-pkin(1) - pkin(2)) * t156;
t174 = t94 * rSges(7,1) + t93 * rSges(7,2);
t173 = -rSges(5,1) * t158 + rSges(5,2) * t155;
t172 = -rSges(6,1) * t101 - rSges(6,2) * t100;
t169 = -Icges(5,1) * t158 + t205;
t168 = Icges(5,2) * t155 - t204;
t167 = -Icges(5,5) * t158 + Icges(5,6) * t155;
t166 = -rSges(5,1) * t201 + rSges(5,2) * t202 + t132 * rSges(5,3);
t165 = t131 * pkin(7) + t175;
t163 = -t131 * pkin(3) + t132 * pkin(7) + t181;
t162 = -t174 + t184;
t161 = -rSges(7,3) * t199 + t174;
t141 = -t155 * pkin(4) + t158 * pkin(8);
t140 = rSges(2,1) * t159 - rSges(2,2) * t156;
t139 = -rSges(2,1) * t156 - rSges(2,2) * t159;
t121 = rSges(3,1) * t159 + rSges(3,3) * t156 + t187;
t120 = rSges(3,3) * t159 + t150 + (-rSges(3,1) - pkin(1)) * t156;
t119 = rSges(6,3) * t158 + (-rSges(6,1) * t157 + rSges(6,2) * t154) * t155;
t112 = t131 * t141;
t104 = (-pkin(8) * t155 - t213) * t132;
t97 = t113 * t199;
t90 = -rSges(4,1) * t131 - rSges(4,2) * t132 + t181;
t89 = rSges(4,1) * t132 - rSges(4,2) * t131 + t175;
t88 = t131 * t188;
t83 = Icges(5,3) * t132 + t167 * t131;
t82 = -Icges(5,3) * t131 + t167 * t132;
t81 = -t119 * t131 - t112;
t80 = (-t119 - t141) * t132;
t79 = (-t155 * t193 + t190) * t158;
t77 = t163 + t166;
t76 = t209 + (pkin(3) - t173) * t132 + t165;
t73 = -rSges(6,3) * t199 - t172;
t66 = -t189 * t131 - t112;
t65 = (-t141 - t189) * t132;
t57 = t158 * t64;
t55 = t161 * t202;
t54 = t131 * t166 + (t173 * t132 - t209) * t132;
t53 = t119 * t202 + t158 * t74;
t52 = -t119 * t199 - t158 * t73;
t51 = t113 * t202 + t57;
t50 = -t158 * t161 - t97;
t45 = t163 + t74 - t188;
t44 = (t213 + pkin(3) + (rSges(6,3) + pkin(8)) * t155) * t132 + t165 + t172;
t41 = t163 + t164 + t64;
t40 = (t208 * t155 + pkin(3) + t198) * t132 + t162 + t165;
t39 = (-t131 * t73 + t132 * t74) * t155;
t37 = -t64 * t199 + t55;
t36 = t158 * t75 + t189 * t202 + t57;
t35 = -t97 + t162 * t158 + (-t158 ^ 2 * t212 + (t185 * t158 - t108) * t155) * t132;
t32 = t131 * t74 - t88 + (t104 + t73) * t132;
t27 = t102 * t70 + t103 * t72 - t68 * t202;
t26 = t102 * t69 + t103 * t71 - t67 * t202;
t25 = t100 * t70 + t101 * t72 - t68 * t199;
t24 = t100 * t69 + t101 * t71 - t67 * t199;
t18 = t55 + (-t184 + (t211 * t155 + t178) * t132) * t202 - t210 * t199;
t17 = -t88 + t210 * t131 + (t132 * t178 - t185 * t199 + t104 - t162) * t132;
t15 = -t131 * t26 + t132 * t27;
t14 = -t131 * t24 + t132 * t25;
t8 = t158 * t49 + (-t131 * t27 - t132 * t26) * t155;
t7 = t158 * t48 + (-t131 * t25 - t132 * t24) * t155;
t1 = [-t158 * (-Icges(5,2) * t158 - t205) + Icges(3,2) + Icges(2,3) + Icges(4,3) + (Icges(5,1) * t155 - t193 - t196 + t204) * t155 + m(7) * (t40 ^ 2 + t41 ^ 2) + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t76 ^ 2 + t77 ^ 2) + m(4) * (t89 ^ 2 + t90 ^ 2) + m(3) * (t120 ^ 2 + t121 ^ 2) + m(2) * (t139 ^ 2 + t140 ^ 2) + t190 + t191; m(7) * (t156 * t40 - t159 * t41) + m(6) * (t156 * t44 - t159 * t45) + m(5) * (t156 * t76 - t159 * t77) + m(4) * (t156 * t89 - t159 * t90) + m(3) * (t120 * t156 - t121 * t159); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t156 ^ 2 + t159 ^ 2); 0; 0; m(4) + m(5) + m(6) + m(7); m(7) * (t66 * t40 + t65 * t41) + m(6) * (t44 * t81 + t45 * t80) + (t47 / 0.2e1 + t31 / 0.2e1 - t77 * t214 + (Icges(5,5) * t132 + t169 * t131) * t222 + (Icges(5,6) * t132 + t168 * t131) * t221 + t132 * t216 - t182) * t132 + (-t46 / 0.2e1 - t30 / 0.2e1 - t76 * t214 + t155 * (-Icges(5,5) * t131 + t169 * t132) / 0.2e1 + (-Icges(5,6) * t131 + t168 * t132) * t215 + t131 * t216 + t183) * t131; m(6) * (t156 * t81 - t159 * t80) + m(7) * (t156 * t66 - t159 * t65) + (-t131 * t156 + t132 * t159) * t214; -m(5) * t54 - m(6) * t32 - m(7) * t17; m(7) * (t17 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t32 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t54 ^ 2 + (t219 + t220) * t138 ^ 2) + (t219 * t83 + t12 + t15) * t132 + (-t220 * t82 - t11 - t14 + (t131 * t83 - t132 * t82) * t132) * t131; t79 + m(7) * (t35 * t40 + t36 * t41) + m(6) * (t44 * t52 + t45 * t53) + (t182 * t131 + t183 * t132) * t155 + t176; m(6) * (t156 * t52 - t159 * t53) + m(7) * (t156 * t35 - t159 * t36); -m(6) * t39 + m(7) * t18; t8 * t217 + (-t33 * t131 + t34 * t132) * t215 + t7 * t218 + (-t132 * t14 / 0.2e1 + t15 * t218) * t155 + m(7) * (-t17 * t18 + t35 * t66 + t36 * t65) + m(6) * (t39 * t32 + t52 * t81 + t53 * t80) + t177; t158 * t79 + m(7) * (t18 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t39 ^ 2 + t52 ^ 2 + t53 ^ 2) + (-t132 * t7 - t131 * t8 + t158 * (-t131 * t34 - t132 * t33)) * t155 + t186; m(7) * (t40 * t50 + t41 * t51) + t176; m(7) * (t156 * t50 - t159 * t51); m(7) * t37; m(7) * (-t17 * t37 + t50 * t66 + t51 * t65) + t177; m(7) * (t18 * t37 + t35 * t50 + t36 * t51) + t186; m(7) * (t37 ^ 2 + t50 ^ 2 + t51 ^ 2) + t186;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
