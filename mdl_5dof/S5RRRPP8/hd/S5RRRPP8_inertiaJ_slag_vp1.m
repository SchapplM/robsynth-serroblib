% Calculate joint inertia matrix for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:31
% EndTime: 2019-12-31 21:07:38
% DurationCPUTime: 2.34s
% Computational Cost: add. (2846->365), mult. (6963->526), div. (0->0), fcn. (7639->6), ass. (0->170)
t162 = sin(qJ(2));
t227 = Icges(3,5) * t162;
t226 = t227 / 0.2e1;
t225 = rSges(6,1) + pkin(4);
t224 = rSges(6,3) + qJ(5);
t161 = sin(qJ(3));
t164 = cos(qJ(3));
t165 = cos(qJ(2));
t105 = -Icges(4,3) * t165 + (Icges(4,5) * t164 - Icges(4,6) * t161) * t162;
t118 = -Icges(6,1) * t165 + (Icges(6,4) * t161 + Icges(6,5) * t164) * t162;
t119 = -Icges(5,1) * t165 + (-Icges(5,4) * t164 + Icges(5,5) * t161) * t162;
t223 = -t105 - t118 - t119;
t163 = sin(qJ(1));
t166 = cos(qJ(1));
t196 = t163 * t165;
t133 = t161 * t196 + t164 * t166;
t193 = t166 * t161;
t134 = t164 * t196 - t193;
t199 = t162 * t163;
t63 = Icges(6,5) * t199 + Icges(6,6) * t133 + Icges(6,3) * t134;
t67 = Icges(6,4) * t199 + Icges(6,2) * t133 + Icges(6,6) * t134;
t71 = Icges(6,1) * t199 + Icges(6,4) * t133 + Icges(6,5) * t134;
t19 = t133 * t67 + t134 * t63 + t71 * t199;
t135 = -t163 * t164 + t165 * t193;
t194 = t165 * t166;
t136 = t163 * t161 + t164 * t194;
t197 = t162 * t166;
t64 = Icges(6,5) * t197 + Icges(6,6) * t135 + Icges(6,3) * t136;
t68 = Icges(6,4) * t197 + Icges(6,2) * t135 + Icges(6,6) * t136;
t72 = Icges(6,1) * t197 + Icges(6,4) * t135 + Icges(6,5) * t136;
t20 = t133 * t68 + t134 * t64 + t72 * t199;
t65 = Icges(5,5) * t199 - Icges(5,6) * t134 + Icges(5,3) * t133;
t69 = Icges(5,4) * t199 - Icges(5,2) * t134 + Icges(5,6) * t133;
t73 = Icges(5,1) * t199 - Icges(5,4) * t134 + Icges(5,5) * t133;
t21 = t133 * t65 - t134 * t69 + t73 * t199;
t66 = Icges(5,5) * t197 - Icges(5,6) * t136 + Icges(5,3) * t135;
t70 = Icges(5,4) * t197 - Icges(5,2) * t136 + Icges(5,6) * t135;
t74 = Icges(5,1) * t197 - Icges(5,4) * t136 + Icges(5,5) * t135;
t22 = t133 * t66 - t134 * t70 + t74 * t199;
t75 = Icges(4,5) * t134 - Icges(4,6) * t133 + Icges(4,3) * t199;
t77 = Icges(4,4) * t134 - Icges(4,2) * t133 + Icges(4,6) * t199;
t79 = Icges(4,1) * t134 - Icges(4,4) * t133 + Icges(4,5) * t199;
t27 = -t133 * t77 + t134 * t79 + t75 * t199;
t76 = Icges(4,5) * t136 - Icges(4,6) * t135 + Icges(4,3) * t197;
t78 = Icges(4,4) * t136 - Icges(4,2) * t135 + Icges(4,6) * t197;
t80 = Icges(4,1) * t136 - Icges(4,4) * t135 + Icges(4,5) * t197;
t28 = -t133 * t78 + t134 * t80 + t76 * t199;
t114 = -Icges(6,5) * t165 + (Icges(6,6) * t161 + Icges(6,3) * t164) * t162;
t116 = -Icges(6,4) * t165 + (Icges(6,2) * t161 + Icges(6,6) * t164) * t162;
t44 = t114 * t134 + t116 * t133 + t118 * t199;
t115 = -Icges(5,5) * t165 + (-Icges(5,6) * t164 + Icges(5,3) * t161) * t162;
t117 = -Icges(5,4) * t165 + (-Icges(5,2) * t164 + Icges(5,6) * t161) * t162;
t45 = t115 * t133 - t117 * t134 + t119 * t199;
t108 = -Icges(4,6) * t165 + (Icges(4,4) * t164 - Icges(4,2) * t161) * t162;
t111 = -Icges(4,5) * t165 + (Icges(4,1) * t164 - Icges(4,4) * t161) * t162;
t48 = t105 * t199 - t108 * t133 + t111 * t134;
t222 = (-t44 - t45 - t48) * t165 + ((t20 + t22 + t28) * t166 + (t19 + t21 + t27) * t163) * t162;
t23 = t135 * t67 + t136 * t63 + t71 * t197;
t24 = t135 * t68 + t136 * t64 + t72 * t197;
t25 = t135 * t65 - t136 * t69 + t73 * t197;
t26 = t135 * t66 - t136 * t70 + t74 * t197;
t29 = -t135 * t77 + t136 * t79 + t75 * t197;
t30 = -t135 * t78 + t136 * t80 + t76 * t197;
t46 = t136 * t114 + t135 * t116 + t118 * t197;
t47 = t135 * t115 - t136 * t117 + t119 * t197;
t49 = t105 * t197 - t135 * t108 + t136 * t111;
t221 = (-t46 - t47 - t49) * t165 + ((t24 + t26 + t30) * t166 + (t23 + t25 + t29) * t163) * t162;
t31 = -t165 * t75 + (-t161 * t77 + t164 * t79) * t162;
t33 = -t165 * t71 + (t161 * t67 + t164 * t63) * t162;
t35 = -t165 * t73 + (t161 * t65 - t164 * t69) * t162;
t220 = -t31 - t33 - t35;
t32 = -t165 * t76 + (-t161 * t78 + t164 * t80) * t162;
t34 = -t165 * t72 + (t161 * t68 + t164 * t64) * t162;
t36 = -t165 * t74 + (t161 * t66 - t164 * t70) * t162;
t219 = t32 + t34 + t36;
t198 = t162 * t164;
t200 = t161 * t162;
t218 = (t115 + t116) * t200 + (t111 + t114) * t198;
t217 = t163 ^ 2;
t216 = t166 ^ 2;
t215 = t163 / 0.2e1;
t214 = -t165 / 0.2e1;
t143 = rSges(3,1) * t162 + rSges(3,2) * t165;
t212 = m(3) * t143;
t211 = pkin(2) * t165;
t84 = rSges(5,1) * t197 - t136 * rSges(5,2) + t135 * rSges(5,3);
t95 = t136 * pkin(3) + qJ(4) * t135;
t210 = -t84 - t95;
t209 = rSges(6,2) * t133;
t208 = rSges(5,3) * t133;
t207 = t166 * rSges(3,3);
t205 = t224 * t134 + t225 * t199 + t209;
t204 = t135 * rSges(6,2) + t224 * t136 + t225 * t197;
t137 = (pkin(3) * t164 + qJ(4) * t161) * t162;
t125 = t133 * qJ(4);
t94 = pkin(3) * t134 + t125;
t203 = t137 * t199 + t165 * t94;
t201 = Icges(3,4) * t165;
t195 = t164 * t117;
t120 = -rSges(4,3) * t165 + (rSges(4,1) * t164 - rSges(4,2) * t161) * t162;
t146 = pkin(2) * t162 - pkin(7) * t165;
t192 = -t120 - t146;
t191 = (rSges(6,2) * t161 + rSges(6,3) * t164) * t162 + qJ(5) * t198 - t225 * t165;
t122 = -rSges(5,1) * t165 + (-rSges(5,2) * t164 + rSges(5,3) * t161) * t162;
t190 = -t122 - t137;
t188 = pkin(2) * t194 + pkin(7) * t197;
t189 = t217 * (pkin(7) * t162 + t211) + t166 * t188;
t187 = t166 * pkin(1) + t163 * pkin(6);
t158 = t166 * pkin(6);
t186 = t158 - t125;
t185 = -t108 * t200 - t162 * t195 + t223 * t165 + t218;
t184 = -t95 - t204;
t183 = -t137 - t191;
t182 = -t146 + t190;
t86 = t136 * rSges(4,1) - t135 * rSges(4,2) + rSges(4,3) * t197;
t181 = -pkin(1) - t211;
t180 = t163 * t94 + t166 * t95 + t189;
t179 = -t146 + t183;
t178 = t187 + t188;
t177 = rSges(3,1) * t165 - rSges(3,2) * t162;
t176 = -rSges(4,1) * t134 + rSges(4,2) * t133;
t174 = -Icges(3,2) * t162 + t201;
t173 = Icges(3,5) * t165 - Icges(3,6) * t162;
t170 = rSges(3,1) * t194 - rSges(3,2) * t197 + t163 * rSges(3,3);
t169 = t178 + t95;
t168 = t31 / 0.2e1 + t35 / 0.2e1 + t33 / 0.2e1 + t48 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1;
t167 = t32 / 0.2e1 + t36 / 0.2e1 + t47 / 0.2e1 + t46 / 0.2e1 + t34 / 0.2e1 + t49 / 0.2e1;
t160 = t162 ^ 2;
t145 = rSges(2,1) * t166 - t163 * rSges(2,2);
t144 = -t163 * rSges(2,1) - rSges(2,2) * t166;
t140 = Icges(3,6) * t165 + t227;
t107 = Icges(3,3) * t163 + t173 * t166;
t106 = -Icges(3,3) * t166 + t173 * t163;
t97 = t170 + t187;
t96 = t207 + t158 + (-pkin(1) - t177) * t163;
t93 = t192 * t166;
t92 = t192 * t163;
t87 = t94 * t197;
t85 = rSges(4,3) * t199 - t176;
t82 = rSges(5,1) * t199 - rSges(5,2) * t134 + t208;
t62 = t166 * t170 + (t177 * t163 - t207) * t163;
t61 = t182 * t166;
t60 = t182 * t163;
t59 = t178 + t86;
t58 = t158 + ((-rSges(4,3) - pkin(7)) * t162 + t181) * t163 + t176;
t57 = t179 * t166;
t56 = t179 * t163;
t55 = -t120 * t197 - t165 * t86;
t54 = t120 * t199 + t165 * t85;
t50 = (-t163 * t86 + t166 * t85) * t162;
t43 = t169 + t84;
t42 = -t208 + (rSges(5,2) - pkin(3)) * t134 + ((-rSges(5,1) - pkin(7)) * t162 + t181) * t163 + t186;
t41 = t169 + t204;
t40 = -t209 + (-pkin(3) - t224) * t134 + ((-pkin(7) - t225) * t162 + t181) * t163 + t186;
t39 = t163 * t85 + t166 * t86 + t189;
t38 = t210 * t165 + t190 * t197;
t37 = t122 * t199 + t165 * t82 + t203;
t18 = t184 * t165 + t183 * t197;
t17 = t205 * t165 + t191 * t199 + t203;
t16 = t87 + (t210 * t163 + t166 * t82) * t162;
t15 = t163 * t82 + t166 * t84 + t180;
t14 = t87 + (t184 * t163 + t205 * t166) * t162;
t13 = t205 * t163 + t204 * t166 + t180;
t12 = t30 * t163 - t166 * t29;
t11 = t28 * t163 - t166 * t27;
t10 = t26 * t163 - t166 * t25;
t9 = t24 * t163 - t166 * t23;
t8 = t22 * t163 - t166 * t21;
t7 = t20 * t163 - t166 * t19;
t1 = [Icges(2,3) + (Icges(3,1) * t162 - t161 * t108 - t195 + t201) * t162 + (Icges(3,4) * t162 + Icges(3,2) * t165 + t223) * t165 + m(6) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2) + m(3) * (t96 ^ 2 + t97 ^ 2) + m(2) * (t144 ^ 2 + t145 ^ 2) + t218; m(6) * (t40 * t57 + t41 * t56) + m(5) * (t42 * t61 + t43 * t60) + m(4) * (t58 * t93 + t59 * t92) + (t174 * t163 * t214 - t96 * t212 - t168 + (-Icges(3,6) * t214 + t226 + t140 / 0.2e1) * t166) * t166 + (-t97 * t212 + (Icges(3,6) * t163 + t174 * t166) * t165 / 0.2e1 + t163 * t226 + t140 * t215 + t167) * t163; m(6) * (t13 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t39 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(3) * (t62 ^ 2 + (t216 + t217) * t143 ^ 2) + (-t216 * t106 - t11 - t7 - t8) * t166 + (t217 * t107 + t10 + t12 + t9 + (-t163 * t106 + t166 * t107) * t166) * t163; -t185 * t165 + m(6) * (t17 * t40 + t18 * t41) + m(5) * (t37 * t42 + t38 * t43) + m(4) * (t54 * t58 + t55 * t59) + (t168 * t163 + t167 * t166) * t162; m(6) * (t13 * t14 + t17 * t57 + t18 * t56) + m(5) * (t15 * t16 + t37 * t61 + t38 * t60) + m(4) * (t39 * t50 + t54 * t93 + t55 * t92) + ((t9 / 0.2e1 + t10 / 0.2e1 + t12 / 0.2e1) * t166 + (t7 / 0.2e1 + t11 / 0.2e1 + t8 / 0.2e1) * t163) * t162 + t221 * t215 + (t219 * t163 + t220 * t166) * t214 - t222 * t166 / 0.2e1; m(6) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(5) * (t16 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(4) * (t50 ^ 2 + t54 ^ 2 + t55 ^ 2) + t185 * t165 ^ 2 + ((-t219 * t165 + t221) * t166 + (t220 * t165 + t222) * t163) * t162; m(6) * (t133 * t41 + t135 * t40) + m(5) * (t133 * t43 + t135 * t42); m(6) * (t13 * t200 + t133 * t56 + t135 * t57) + m(5) * (t133 * t60 + t135 * t61 + t15 * t200); m(6) * (t133 * t18 + t135 * t17 + t14 * t200) + m(5) * (t133 * t38 + t135 * t37 + t16 * t200); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t160 * t161 ^ 2 + t133 ^ 2 + t135 ^ 2); m(6) * (t134 * t41 + t136 * t40); m(6) * (t13 * t198 + t134 * t56 + t136 * t57); m(6) * (t134 * t18 + t136 * t17 + t14 * t198); m(6) * (t160 * t161 * t164 + t133 * t134 + t135 * t136); m(6) * (t160 * t164 ^ 2 + t134 ^ 2 + t136 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
