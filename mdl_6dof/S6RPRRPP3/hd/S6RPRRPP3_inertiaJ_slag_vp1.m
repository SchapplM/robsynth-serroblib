% Calculate joint inertia matrix for
% S6RPRRPP3
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:09
% EndTime: 2019-03-09 04:35:14
% DurationCPUTime: 2.40s
% Computational Cost: add. (6018->384), mult. (7341->545), div. (0->0), fcn. (8009->8), ass. (0->178)
t168 = sin(qJ(3));
t236 = Icges(4,5) * t168;
t235 = t236 / 0.2e1;
t234 = rSges(7,1) + pkin(5);
t233 = rSges(7,3) + qJ(6);
t167 = sin(qJ(4));
t170 = cos(qJ(4));
t171 = cos(qJ(3));
t125 = -Icges(5,3) * t171 + (Icges(5,5) * t170 - Icges(5,6) * t167) * t168;
t132 = -Icges(7,1) * t171 + (Icges(7,4) * t167 + Icges(7,5) * t170) * t168;
t133 = -Icges(6,1) * t171 + (-Icges(6,4) * t170 + Icges(6,5) * t167) * t168;
t232 = t125 + t132 + t133;
t165 = qJ(1) + pkin(9);
t162 = sin(t165);
t163 = cos(t165);
t204 = t167 * t171;
t121 = t162 * t204 + t163 * t170;
t201 = t170 * t171;
t122 = t162 * t201 - t163 * t167;
t208 = t162 * t168;
t63 = Icges(7,5) * t208 + Icges(7,6) * t121 + Icges(7,3) * t122;
t67 = Icges(7,4) * t208 + Icges(7,2) * t121 + Icges(7,6) * t122;
t71 = Icges(7,1) * t208 + Icges(7,4) * t121 + Icges(7,5) * t122;
t17 = t121 * t67 + t122 * t63 + t71 * t208;
t123 = -t162 * t170 + t163 * t204;
t124 = t162 * t167 + t163 * t201;
t207 = t163 * t168;
t64 = Icges(7,5) * t207 + Icges(7,6) * t123 + Icges(7,3) * t124;
t68 = Icges(7,4) * t207 + Icges(7,2) * t123 + Icges(7,6) * t124;
t72 = Icges(7,1) * t207 + Icges(7,4) * t123 + Icges(7,5) * t124;
t18 = t121 * t68 + t122 * t64 + t72 * t208;
t65 = Icges(6,5) * t208 - Icges(6,6) * t122 + Icges(6,3) * t121;
t69 = Icges(6,4) * t208 - Icges(6,2) * t122 + Icges(6,6) * t121;
t73 = Icges(6,1) * t208 - Icges(6,4) * t122 + Icges(6,5) * t121;
t19 = t121 * t65 - t122 * t69 + t73 * t208;
t66 = Icges(6,5) * t207 - Icges(6,6) * t124 + Icges(6,3) * t123;
t70 = Icges(6,4) * t207 - Icges(6,2) * t124 + Icges(6,6) * t123;
t74 = Icges(6,1) * t207 - Icges(6,4) * t124 + Icges(6,5) * t123;
t20 = t121 * t66 - t122 * t70 + t74 * t208;
t75 = Icges(5,5) * t122 - Icges(5,6) * t121 + Icges(5,3) * t208;
t77 = Icges(5,4) * t122 - Icges(5,2) * t121 + Icges(5,6) * t208;
t79 = Icges(5,1) * t122 - Icges(5,4) * t121 + Icges(5,5) * t208;
t25 = -t121 * t77 + t122 * t79 + t75 * t208;
t76 = Icges(5,5) * t124 - Icges(5,6) * t123 + Icges(5,3) * t207;
t78 = Icges(5,4) * t124 - Icges(5,2) * t123 + Icges(5,6) * t207;
t80 = Icges(5,1) * t124 - Icges(5,4) * t123 + Icges(5,5) * t207;
t26 = -t121 * t78 + t122 * t80 + t76 * t208;
t128 = -Icges(7,5) * t171 + (Icges(7,6) * t167 + Icges(7,3) * t170) * t168;
t130 = -Icges(7,4) * t171 + (Icges(7,2) * t167 + Icges(7,6) * t170) * t168;
t45 = t121 * t130 + t122 * t128 + t132 * t208;
t129 = -Icges(6,5) * t171 + (-Icges(6,6) * t170 + Icges(6,3) * t167) * t168;
t131 = -Icges(6,4) * t171 + (-Icges(6,2) * t170 + Icges(6,6) * t167) * t168;
t46 = t121 * t129 - t122 * t131 + t133 * t208;
t126 = -Icges(5,6) * t171 + (Icges(5,4) * t170 - Icges(5,2) * t167) * t168;
t127 = -Icges(5,5) * t171 + (Icges(5,1) * t170 - Icges(5,4) * t167) * t168;
t49 = -t121 * t126 + t122 * t127 + t125 * t208;
t231 = (-t45 - t46 - t49) * t171 + ((t18 + t20 + t26) * t163 + (t17 + t19 + t25) * t162) * t168;
t21 = t123 * t67 + t124 * t63 + t71 * t207;
t22 = t123 * t68 + t124 * t64 + t72 * t207;
t23 = t123 * t65 - t124 * t69 + t73 * t207;
t24 = t123 * t66 - t124 * t70 + t74 * t207;
t27 = -t123 * t77 + t124 * t79 + t75 * t207;
t28 = -t123 * t78 + t124 * t80 + t76 * t207;
t47 = t123 * t130 + t124 * t128 + t132 * t207;
t48 = t123 * t129 - t124 * t131 + t133 * t207;
t50 = -t123 * t126 + t124 * t127 + t125 * t207;
t230 = (-t47 - t48 - t50) * t171 + ((t22 + t24 + t28) * t163 + (t21 + t23 + t27) * t162) * t168;
t31 = -t171 * t75 + (-t167 * t77 + t170 * t79) * t168;
t33 = -t171 * t71 + (t167 * t67 + t170 * t63) * t168;
t35 = -t171 * t73 + (t167 * t65 - t170 * t69) * t168;
t229 = -t31 - t33 - t35;
t32 = -t171 * t76 + (-t167 * t78 + t170 * t80) * t168;
t34 = -t171 * t72 + (t167 * t68 + t170 * t64) * t168;
t36 = -t171 * t74 + (t167 * t66 - t170 * t70) * t168;
t228 = t32 + t34 + t36;
t203 = t168 * t170;
t205 = t167 * t168;
t227 = (t129 + t130) * t205 + (t127 + t128) * t203;
t226 = t162 ^ 2;
t225 = t163 ^ 2;
t224 = m(6) + m(7);
t223 = t162 / 0.2e1;
t221 = -t171 / 0.2e1;
t145 = t168 * rSges(4,1) + t171 * rSges(4,2);
t220 = m(4) * t145;
t219 = pkin(3) * t171;
t169 = sin(qJ(1));
t218 = t169 * pkin(1);
t214 = t121 * rSges(7,2);
t217 = t233 * t122 + t234 * t208 + t214;
t216 = t123 * rSges(7,2) + t233 * t124 + t234 * t207;
t84 = rSges(6,1) * t207 - t124 * rSges(6,2) + t123 * rSges(6,3);
t93 = t124 * pkin(4) + t123 * qJ(5);
t215 = -t84 - t93;
t213 = t121 * rSges(6,3);
t212 = t163 * rSges(4,3);
t139 = (pkin(4) * t170 + qJ(5) * t167) * t168;
t113 = t121 * qJ(5);
t92 = t122 * pkin(4) + t113;
t211 = t139 * t208 + t171 * t92;
t209 = Icges(4,4) * t171;
t206 = t163 * t171;
t202 = t170 * t131;
t195 = pkin(3) * t206 + pkin(8) * t207;
t199 = t226 * (pkin(8) * t168 + t219) + t163 * t195;
t134 = -t171 * rSges(5,3) + (rSges(5,1) * t170 - rSges(5,2) * t167) * t168;
t152 = t168 * pkin(3) - t171 * pkin(8);
t198 = -t134 - t152;
t197 = (rSges(7,2) * t167 + rSges(7,3) * t170) * t168 + qJ(6) * t203 - t234 * t171;
t136 = -t171 * rSges(6,1) + (-rSges(6,2) * t170 + rSges(6,3) * t167) * t168;
t196 = -t136 - t139;
t194 = t126 * t205 + t168 * t202 + t232 * t171 - t227;
t193 = -t93 - t216;
t86 = t124 * rSges(5,1) - t123 * rSges(5,2) + rSges(5,3) * t207;
t192 = -t139 - t197;
t191 = -t152 + t196;
t172 = cos(qJ(1));
t164 = t172 * pkin(1);
t190 = t163 * pkin(2) + t162 * pkin(7) + t164;
t189 = -pkin(2) - t219;
t188 = t163 * pkin(7) - t218;
t187 = t162 * t92 + t163 * t93 + t199;
t186 = -t152 + t192;
t185 = -t113 + t188;
t184 = rSges(4,1) * t171 - rSges(4,2) * t168;
t183 = -t122 * rSges(5,1) + t121 * rSges(5,2);
t182 = t190 + t195;
t180 = -Icges(4,2) * t168 + t209;
t179 = Icges(4,5) * t171 - Icges(4,6) * t168;
t176 = rSges(4,1) * t206 - rSges(4,2) * t207 + t162 * rSges(4,3);
t175 = t35 / 0.2e1 + t49 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1 + t33 / 0.2e1 + t31 / 0.2e1;
t174 = t47 / 0.2e1 + t34 / 0.2e1 + t32 / 0.2e1 + t48 / 0.2e1 + t36 / 0.2e1 + t50 / 0.2e1;
t173 = t182 + t93;
t166 = t168 ^ 2;
t147 = t172 * rSges(2,1) - t169 * rSges(2,2);
t146 = -t169 * rSges(2,1) - t172 * rSges(2,2);
t142 = Icges(4,6) * t171 + t236;
t138 = t163 * rSges(3,1) - t162 * rSges(3,2) + t164;
t137 = -t162 * rSges(3,1) - t163 * rSges(3,2) - t218;
t101 = Icges(4,3) * t162 + t179 * t163;
t100 = -Icges(4,3) * t163 + t179 * t162;
t97 = t198 * t163;
t96 = t198 * t162;
t95 = t176 + t190;
t94 = t212 + (-pkin(2) - t184) * t162 + t188;
t87 = t92 * t207;
t85 = rSges(5,3) * t208 - t183;
t82 = rSges(6,1) * t208 - t122 * rSges(6,2) + t213;
t62 = t191 * t163;
t61 = t191 * t162;
t60 = t163 * t176 + (t184 * t162 - t212) * t162;
t56 = t186 * t163;
t55 = t186 * t162;
t54 = -t134 * t207 - t171 * t86;
t53 = t134 * t208 + t171 * t85;
t52 = t182 + t86;
t51 = ((-rSges(5,3) - pkin(8)) * t168 + t189) * t162 + t183 + t188;
t44 = (-t162 * t86 + t163 * t85) * t168;
t43 = t173 + t84;
t42 = -t213 + (rSges(6,2) - pkin(4)) * t122 + ((-rSges(6,1) - pkin(8)) * t168 + t189) * t162 + t185;
t41 = t215 * t171 + t196 * t207;
t40 = t136 * t208 + t171 * t82 + t211;
t39 = t173 + t216;
t38 = -t214 + (-pkin(4) - t233) * t122 + ((-pkin(8) - t234) * t168 + t189) * t162 + t185;
t37 = t162 * t85 + t163 * t86 + t199;
t30 = t193 * t171 + t192 * t207;
t29 = t217 * t171 + t197 * t208 + t211;
t16 = t87 + (t215 * t162 + t163 * t82) * t168;
t15 = t162 * t82 + t163 * t84 + t187;
t14 = t87 + (t193 * t162 + t217 * t163) * t168;
t13 = t217 * t162 + t216 * t163 + t187;
t12 = t28 * t162 - t27 * t163;
t11 = t26 * t162 - t25 * t163;
t10 = t24 * t162 - t23 * t163;
t9 = t22 * t162 - t21 * t163;
t8 = t20 * t162 - t19 * t163;
t7 = t18 * t162 - t17 * t163;
t1 = [Icges(2,3) + Icges(3,3) + (Icges(4,1) * t168 - t167 * t126 - t202 + t209) * t168 + (Icges(4,4) * t168 + Icges(4,2) * t171 - t232) * t171 + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t51 ^ 2 + t52 ^ 2) + m(4) * (t94 ^ 2 + t95 ^ 2) + m(3) * (t137 ^ 2 + t138 ^ 2) + m(2) * (t146 ^ 2 + t147 ^ 2) + t227; 0; m(3) + m(4) + m(5) + t224; m(7) * (t56 * t38 + t55 * t39) + m(6) * (t62 * t42 + t61 * t43) + m(5) * (t97 * t51 + t96 * t52) + (t180 * t162 * t221 - t94 * t220 - t175 + (-Icges(4,6) * t221 + t235 + t142 / 0.2e1) * t163) * t163 + (-t95 * t220 + t171 * (Icges(4,6) * t162 + t180 * t163) / 0.2e1 + t162 * t235 + t142 * t223 + t174) * t162; m(4) * t60 + m(5) * t37 + m(6) * t15 + m(7) * t13; m(7) * (t13 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t15 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t37 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(4) * (t60 ^ 2 + (t225 + t226) * t145 ^ 2) + (-t225 * t100 - t11 - t7 - t8) * t163 + (t226 * t101 + t10 + t12 + t9 + (-t162 * t100 + t163 * t101) * t163) * t162; t194 * t171 + m(7) * (t29 * t38 + t30 * t39) + m(6) * (t40 * t42 + t41 * t43) + m(5) * (t53 * t51 + t54 * t52) + (t175 * t162 + t174 * t163) * t168; m(5) * t44 + m(6) * t16 + m(7) * t14; m(7) * (t14 * t13 + t29 * t56 + t30 * t55) + m(6) * (t16 * t15 + t40 * t62 + t41 * t61) + m(5) * (t44 * t37 + t53 * t97 + t54 * t96) + ((t9 / 0.2e1 + t12 / 0.2e1 + t10 / 0.2e1) * t163 + (t11 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t162) * t168 + t230 * t223 - t231 * t163 / 0.2e1 + (t228 * t162 + t229 * t163) * t221; m(7) * (t14 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t16 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t44 ^ 2 + t53 ^ 2 + t54 ^ 2) - t194 * t171 ^ 2 + ((-t228 * t171 + t230) * t163 + (t229 * t171 + t231) * t162) * t168; m(7) * (t121 * t39 + t123 * t38) + m(6) * (t121 * t43 + t123 * t42); t224 * t205; m(7) * (t121 * t55 + t123 * t56 + t13 * t205) + m(6) * (t121 * t61 + t123 * t62 + t15 * t205); m(7) * (t121 * t30 + t123 * t29 + t14 * t205) + m(6) * (t121 * t41 + t123 * t40 + t16 * t205); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t166 * t167 ^ 2 + t121 ^ 2 + t123 ^ 2); m(7) * (t122 * t39 + t124 * t38); m(7) * t203; m(7) * (t122 * t55 + t124 * t56 + t13 * t203); m(7) * (t122 * t30 + t124 * t29 + t14 * t203); m(7) * (t166 * t170 * t167 + t122 * t121 + t124 * t123); m(7) * (t166 * t170 ^ 2 + t122 ^ 2 + t124 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
