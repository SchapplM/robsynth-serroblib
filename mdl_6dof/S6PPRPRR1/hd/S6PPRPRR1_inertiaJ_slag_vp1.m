% Calculate joint inertia matrix for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:07
% EndTime: 2019-03-08 18:42:17
% DurationCPUTime: 5.06s
% Computational Cost: add. (30736->438), mult. (83628->653), div. (0->0), fcn. (111203->16), ass. (0->198)
t221 = Icges(4,3) + Icges(5,3);
t186 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t220 = 0.2e1 * t186;
t198 = sin(pkin(13));
t199 = sin(pkin(7));
t182 = t199 * t198;
t200 = cos(pkin(13));
t183 = t200 * t199;
t212 = sin(qJ(3));
t214 = cos(qJ(3));
t150 = t214 * t182 + t212 * t183;
t170 = cos(pkin(7));
t177 = t214 * t198 + t212 * t200;
t152 = t177 * t170;
t165 = sin(pkin(12));
t166 = sin(pkin(11));
t168 = cos(pkin(12));
t169 = cos(pkin(11));
t171 = cos(pkin(6));
t193 = t169 * t171;
t157 = -t165 * t166 + t168 * t193;
t158 = t165 * t193 + t166 * t168;
t162 = -t212 * t198 + t214 * t200;
t167 = sin(pkin(6));
t195 = t167 * t169;
t125 = -t150 * t195 + t152 * t157 + t158 * t162;
t173 = sin(qJ(5));
t194 = t167 * t170;
t178 = t157 * t199 + t169 * t194;
t213 = cos(qJ(5));
t118 = t125 * t173 + t178 * t213;
t196 = t166 * t171;
t159 = -t165 * t169 - t168 * t196;
t160 = -t165 * t196 + t168 * t169;
t197 = t166 * t167;
t127 = t150 * t197 + t152 * t159 + t160 * t162;
t179 = t159 * t199 - t166 * t194;
t120 = t127 * t173 + t179 * t213;
t137 = t150 * t171 + (t152 * t168 + t162 * t165) * t167;
t156 = -t167 * t168 * t199 + t171 * t170;
t128 = t137 * t173 - t156 * t213;
t119 = t125 * t213 - t173 * t178;
t151 = t162 * t170;
t176 = -t212 * t182 + t214 * t183;
t175 = t167 * t176;
t124 = t151 * t157 - t158 * t177 - t169 * t175;
t172 = sin(qJ(6));
t174 = cos(qJ(6));
t81 = -t119 * t172 - t124 * t174;
t82 = t119 * t174 - t124 * t172;
t52 = Icges(7,5) * t82 + Icges(7,6) * t81 + Icges(7,3) * t118;
t54 = Icges(7,4) * t82 + Icges(7,2) * t81 + Icges(7,6) * t118;
t56 = Icges(7,1) * t82 + Icges(7,4) * t81 + Icges(7,5) * t118;
t17 = t118 * t52 + t54 * t81 + t56 * t82;
t121 = t127 * t213 - t173 * t179;
t126 = t159 * t151 - t160 * t177 + t166 * t175;
t83 = -t121 * t172 - t126 * t174;
t84 = t121 * t174 - t126 * t172;
t53 = Icges(7,5) * t84 + Icges(7,6) * t83 + Icges(7,3) * t120;
t55 = Icges(7,4) * t84 + Icges(7,2) * t83 + Icges(7,6) * t120;
t57 = Icges(7,1) * t84 + Icges(7,4) * t83 + Icges(7,5) * t120;
t18 = t118 * t53 + t55 * t81 + t57 * t82;
t129 = t137 * t213 + t156 * t173;
t136 = t171 * t176 + (t168 * t151 - t165 * t177) * t167;
t100 = t129 * t174 - t136 * t172;
t99 = -t129 * t172 - t136 * t174;
t60 = Icges(7,5) * t100 + Icges(7,6) * t99 + Icges(7,3) * t128;
t61 = Icges(7,4) * t100 + Icges(7,2) * t99 + Icges(7,6) * t128;
t62 = Icges(7,1) * t100 + Icges(7,4) * t99 + Icges(7,5) * t128;
t26 = t118 * t60 + t61 * t81 + t62 * t82;
t1 = t118 * t17 + t120 * t18 + t128 * t26;
t219 = -t1 / 0.2e1;
t19 = t120 * t52 + t54 * t83 + t56 * t84;
t20 = t120 * t53 + t55 * t83 + t57 * t84;
t27 = t120 * t60 + t61 * t83 + t62 * t84;
t2 = t118 * t19 + t120 * t20 + t128 * t27;
t218 = -t2 / 0.2e1;
t217 = t118 / 0.2e1;
t216 = t120 / 0.2e1;
t215 = t128 / 0.2e1;
t211 = t214 * pkin(3);
t58 = rSges(7,1) * t82 + rSges(7,2) * t81 + rSges(7,3) * t118;
t210 = pkin(5) * t119 + pkin(10) * t118 + t58;
t59 = rSges(7,1) * t84 + rSges(7,2) * t83 + rSges(7,3) * t120;
t209 = pkin(5) * t121 + pkin(10) * t120 + t59;
t63 = rSges(7,1) * t100 + rSges(7,2) * t99 + rSges(7,3) * t128;
t208 = pkin(5) * t129 + pkin(10) * t128 + t63;
t184 = t199 * t212;
t154 = pkin(3) * t184 + (pkin(8) + qJ(4)) * t170;
t188 = t199 * pkin(8);
t189 = t170 * t212;
t155 = pkin(3) * t189 - t199 * qJ(4) - t188;
t116 = t178 * pkin(8) - t154 * t195 + t157 * t155 + t211 * t158;
t101 = t179 * t116;
t95 = pkin(4) * t125 - pkin(9) * t124;
t207 = -t179 * t95 - t101;
t117 = t179 * pkin(8) + t154 * t197 + t159 * t155 + t211 * t160;
t106 = t156 * t117;
t96 = pkin(4) * t127 - pkin(9) * t126;
t206 = t156 * t96 + t106;
t185 = t214 * t199;
t181 = t167 * t185;
t190 = t170 * t214;
t138 = t157 * t190 - t158 * t212 - t169 * t181;
t180 = t167 * t184;
t139 = t157 * t189 + t158 * t214 - t169 * t180;
t205 = Icges(4,5) * t139 + Icges(5,5) * t125 + Icges(4,6) * t138 + Icges(5,6) * t124 - t178 * t221;
t140 = t159 * t190 - t160 * t212 + t166 * t181;
t141 = t159 * t189 + t160 * t214 + t166 * t180;
t204 = Icges(4,5) * t141 + Icges(5,5) * t127 + Icges(4,6) * t140 + Icges(5,6) * t126 - t179 * t221;
t203 = -t116 - t95;
t202 = -t117 - t96;
t115 = pkin(4) * t137 - pkin(9) * t136;
t134 = (-t170 * pkin(8) + t154) * t171 + ((t188 + t155) * t168 + t211 * t165) * t167;
t122 = t178 * t134;
t201 = -t115 * t178 - t122;
t144 = t171 * t185 + (-t212 * t165 + t168 * t190) * t167;
t145 = t171 * t184 + (t214 * t165 + t168 * t189) * t167;
t192 = Icges(4,5) * t145 + Icges(5,5) * t137 + Icges(4,6) * t144 + Icges(5,6) * t136 + t156 * t221;
t191 = -t115 - t134;
t187 = m(3) + m(4) + m(5) + m(6) + m(7);
t133 = rSges(4,1) * t145 + rSges(4,2) * t144 + rSges(4,3) * t156;
t132 = Icges(4,1) * t145 + Icges(4,4) * t144 + Icges(4,5) * t156;
t131 = Icges(4,4) * t145 + Icges(4,2) * t144 + Icges(4,6) * t156;
t114 = rSges(4,1) * t141 + rSges(4,2) * t140 - rSges(4,3) * t179;
t113 = rSges(4,1) * t139 + rSges(4,2) * t138 - rSges(4,3) * t178;
t112 = Icges(4,1) * t141 + Icges(4,4) * t140 - Icges(4,5) * t179;
t111 = Icges(4,1) * t139 + Icges(4,4) * t138 - Icges(4,5) * t178;
t110 = Icges(4,4) * t141 + Icges(4,2) * t140 - Icges(4,6) * t179;
t109 = Icges(4,4) * t139 + Icges(4,2) * t138 - Icges(4,6) * t178;
t105 = rSges(5,1) * t137 + rSges(5,2) * t136 + rSges(5,3) * t156;
t104 = Icges(5,1) * t137 + Icges(5,4) * t136 + Icges(5,5) * t156;
t103 = Icges(5,4) * t137 + Icges(5,2) * t136 + Icges(5,6) * t156;
t92 = rSges(5,1) * t127 + rSges(5,2) * t126 - rSges(5,3) * t179;
t91 = rSges(5,1) * t125 + rSges(5,2) * t124 - rSges(5,3) * t178;
t90 = Icges(5,1) * t127 + Icges(5,4) * t126 - Icges(5,5) * t179;
t89 = Icges(5,1) * t125 + Icges(5,4) * t124 - Icges(5,5) * t178;
t88 = Icges(5,4) * t127 + Icges(5,2) * t126 - Icges(5,6) * t179;
t87 = Icges(5,4) * t125 + Icges(5,2) * t124 - Icges(5,6) * t178;
t80 = rSges(6,1) * t129 - rSges(6,2) * t128 - rSges(6,3) * t136;
t79 = Icges(6,1) * t129 - Icges(6,4) * t128 - Icges(6,5) * t136;
t78 = Icges(6,4) * t129 - Icges(6,2) * t128 - Icges(6,6) * t136;
t77 = Icges(6,5) * t129 - Icges(6,6) * t128 - Icges(6,3) * t136;
t74 = t114 * t156 + t133 * t179;
t73 = -t113 * t156 - t133 * t178;
t72 = -t113 * t179 + t114 * t178;
t71 = rSges(6,1) * t121 - rSges(6,2) * t120 - rSges(6,3) * t126;
t70 = rSges(6,1) * t119 - rSges(6,2) * t118 - rSges(6,3) * t124;
t69 = Icges(6,1) * t121 - Icges(6,4) * t120 - Icges(6,5) * t126;
t68 = Icges(6,1) * t119 - Icges(6,4) * t118 - Icges(6,5) * t124;
t67 = Icges(6,4) * t121 - Icges(6,2) * t120 - Icges(6,6) * t126;
t66 = Icges(6,4) * t119 - Icges(6,2) * t118 - Icges(6,6) * t124;
t65 = Icges(6,5) * t121 - Icges(6,6) * t120 - Icges(6,3) * t126;
t64 = Icges(6,5) * t119 - Icges(6,6) * t118 - Icges(6,3) * t124;
t51 = t156 * t92 + t106 - (-t105 - t134) * t179;
t50 = -t105 * t178 - t122 + (-t116 - t91) * t156;
t49 = t126 * t80 - t136 * t71;
t48 = -t124 * t80 + t136 * t70;
t47 = -t179 * t91 - t101 - (-t117 - t92) * t178;
t46 = t124 * t71 - t126 * t70;
t45 = -t128 * t78 + t129 * t79 - t136 * t77;
t44 = -t120 * t78 + t121 * t79 - t126 * t77;
t43 = -t118 * t78 + t119 * t79 - t124 * t77;
t42 = -t120 * t63 + t128 * t59;
t41 = t118 * t63 - t128 * t58;
t40 = t156 * t71 - (-t80 + t191) * t179 + t206;
t39 = -t178 * t80 + (-t70 + t203) * t156 + t201;
t38 = -t128 * t67 + t129 * t69 - t136 * t65;
t37 = -t128 * t66 + t129 * t68 - t136 * t64;
t36 = -t120 * t67 + t121 * t69 - t126 * t65;
t35 = -t120 * t66 + t121 * t68 - t126 * t64;
t34 = -t118 * t67 + t119 * t69 - t124 * t65;
t33 = -t118 * t66 + t119 * t68 - t124 * t64;
t32 = -t118 * t59 + t120 * t58;
t31 = -t179 * t70 - (-t71 + t202) * t178 + t207;
t30 = t100 * t62 + t128 * t60 + t61 * t99;
t29 = t208 * t126 - t209 * t136;
t28 = -t208 * t124 + t210 * t136;
t25 = t209 * t124 - t210 * t126;
t24 = t209 * t156 - (t191 - t208) * t179 + t206;
t23 = -t208 * t178 + (t203 - t210) * t156 + t201;
t22 = t100 * t57 + t128 * t53 + t55 * t99;
t21 = t100 * t56 + t128 * t52 + t54 * t99;
t16 = -t210 * t179 - (t202 - t209) * t178 + t207;
t15 = t156 * t45 - t178 * t37 - t179 * t38;
t14 = -t124 * t37 - t126 * t38 - t136 * t45;
t13 = t156 * t44 - t178 * t35 - t179 * t36;
t12 = t156 * t43 - t178 * t33 - t179 * t34;
t11 = -t124 * t35 - t126 * t36 - t136 * t44;
t10 = -t124 * t33 - t126 * t34 - t136 * t43;
t9 = t156 * t30 - t178 * t21 - t179 * t22;
t8 = -t124 * t21 - t126 * t22 - t136 * t30;
t7 = t118 * t21 + t120 * t22 + t128 * t30;
t6 = t156 * t27 - t178 * t19 - t179 * t20;
t5 = t156 * t26 - t17 * t178 - t179 * t18;
t4 = -t124 * t19 - t126 * t20 - t136 * t27;
t3 = -t124 * t17 - t126 * t18 - t136 * t26;
t75 = [m(2) + t187; t187 * t171; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t186) * (t171 ^ 2 + (t166 ^ 2 + t169 ^ 2) * t167 ^ 2); m(4) * t72 + m(5) * t47 + m(6) * t31 + m(7) * t16; m(4) * (t171 * t72 + (t166 * t73 - t169 * t74) * t167) + m(5) * (t171 * t47 + (t166 * t50 - t169 * t51) * t167) + m(6) * (t171 * t31 + (t166 * t39 - t169 * t40) * t167) + m(7) * (t16 * t171 + (t166 * t23 - t169 * t24) * t167); m(7) * (t16 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t31 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t47 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(4) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) + (t9 + t15 + (t103 * t136 + t104 * t137 + t131 * t144 + t132 * t145 + t192 * t156) * t156) * t156 - (t6 + t13 - (t110 * t140 + t112 * t141 + t126 * t88 + t127 * t90 - t179 * t204) * t179 + (t103 * t126 + t104 * t127 + t110 * t144 + t112 * t145 + t131 * t140 + t132 * t141 + t136 * t88 + t137 * t90 + t204 * t156 - t179 * t192) * t156) * t179 - (t5 + t12 - (t109 * t138 + t111 * t139 + t124 * t87 + t125 * t89 - t178 * t205) * t178 + (t103 * t124 + t104 * t125 + t109 * t144 + t111 * t145 + t131 * t138 + t132 * t139 + t136 * t87 + t137 * t89 + t205 * t156 - t178 * t192) * t156 - (t109 * t140 + t110 * t138 + t111 * t141 + t112 * t139 + t124 * t88 + t125 * t90 + t126 * t87 + t127 * t89 - t178 * t204 - t179 * t205) * t179) * t178; t156 * t220; (t156 * t171 + (-t166 * t179 + t169 * t178) * t167) * t220; m(7) * (t156 * t16 - t178 * t24 - t179 * t23) + m(6) * (t156 * t31 - t178 * t40 - t179 * t39) + m(5) * (t156 * t47 - t178 * t51 - t179 * t50); (t156 ^ 2 + t178 ^ 2 + t179 ^ 2) * t220; m(6) * t46 + m(7) * t25; m(6) * (t46 * t171 + (t166 * t48 - t169 * t49) * t167) + m(7) * (t25 * t171 + (t166 * t28 - t169 * t29) * t167); (t8 / 0.2e1 + t14 / 0.2e1) * t156 - (t4 / 0.2e1 + t11 / 0.2e1) * t179 - (t3 / 0.2e1 + t10 / 0.2e1) * t178 + (-t9 / 0.2e1 - t15 / 0.2e1) * t136 + (-t6 / 0.2e1 - t13 / 0.2e1) * t126 + (-t5 / 0.2e1 - t12 / 0.2e1) * t124 + m(7) * (t16 * t25 + t23 * t28 + t24 * t29) + m(6) * (t31 * t46 + t39 * t48 + t40 * t49); m(6) * (t156 * t46 - t178 * t49 - t179 * t48) + m(7) * (t156 * t25 - t178 * t29 - t179 * t28); (-t8 - t14) * t136 + (-t4 - t11) * t126 + (-t3 - t10) * t124 + m(7) * (t25 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t46 ^ 2 + t48 ^ 2 + t49 ^ 2); m(7) * t32; m(7) * (t171 * t32 + (t166 * t41 - t169 * t42) * t167); t178 * t219 + t179 * t218 + m(7) * (t16 * t32 + t23 * t41 + t24 * t42) + t9 * t215 + t156 * t7 / 0.2e1 + t5 * t217 + t6 * t216; m(7) * (t156 * t32 - t178 * t42 - t179 * t41); m(7) * (t25 * t32 + t28 * t41 + t29 * t42) + t8 * t215 + t124 * t219 + t4 * t216 + t3 * t217 - t136 * t7 / 0.2e1 + t126 * t218; t120 * t2 + t118 * t1 + t128 * t7 + m(7) * (t32 ^ 2 + t41 ^ 2 + t42 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t75(1) t75(2) t75(4) t75(7) t75(11) t75(16); t75(2) t75(3) t75(5) t75(8) t75(12) t75(17); t75(4) t75(5) t75(6) t75(9) t75(13) t75(18); t75(7) t75(8) t75(9) t75(10) t75(14) t75(19); t75(11) t75(12) t75(13) t75(14) t75(15) t75(20); t75(16) t75(17) t75(18) t75(19) t75(20) t75(21);];
Mq  = res;
