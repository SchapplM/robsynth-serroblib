% Calculate joint inertia matrix for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:44
% EndTime: 2019-12-31 19:09:49
% DurationCPUTime: 1.95s
% Computational Cost: add. (5953->311), mult. (6007->463), div. (0->0), fcn. (6488->10), ass. (0->164)
t139 = pkin(9) + qJ(3);
t134 = sin(t139);
t205 = Icges(4,5) * t134;
t204 = t205 / 0.2e1;
t148 = cos(qJ(4));
t133 = t148 * pkin(4) + pkin(3);
t150 = -pkin(8) - pkin(7);
t146 = sin(qJ(4));
t147 = sin(qJ(1));
t179 = t147 * t146;
t135 = cos(t139);
t149 = cos(qJ(1));
t182 = t135 * t149;
t183 = t134 * t149;
t142 = qJ(4) + qJ(5);
t136 = sin(t142);
t177 = t149 * t136;
t137 = cos(t142);
t180 = t147 * t137;
t107 = -t135 * t177 + t180;
t176 = t149 * t137;
t181 = t147 * t136;
t108 = t135 * t176 + t181;
t65 = t108 * rSges(6,1) + t107 * rSges(6,2) + rSges(6,3) * t183;
t203 = pkin(4) * t179 + t133 * t182 - t150 * t183 + t65;
t140 = t147 ^ 2;
t141 = t149 ^ 2;
t184 = t134 * t147;
t105 = -t135 * t181 - t176;
t106 = t135 * t180 - t177;
t58 = Icges(6,5) * t106 + Icges(6,6) * t105 + Icges(6,3) * t184;
t60 = Icges(6,4) * t106 + Icges(6,2) * t105 + Icges(6,6) * t184;
t62 = Icges(6,1) * t106 + Icges(6,4) * t105 + Icges(6,5) * t184;
t19 = t105 * t60 + t106 * t62 + t58 * t184;
t59 = Icges(6,5) * t108 + Icges(6,6) * t107 + Icges(6,3) * t183;
t61 = Icges(6,4) * t108 + Icges(6,2) * t107 + Icges(6,6) * t183;
t63 = Icges(6,1) * t108 + Icges(6,4) * t107 + Icges(6,5) * t183;
t20 = t105 * t61 + t106 * t63 + t59 * t184;
t85 = -Icges(6,3) * t135 + (Icges(6,5) * t137 - Icges(6,6) * t136) * t134;
t86 = -Icges(6,6) * t135 + (Icges(6,4) * t137 - Icges(6,2) * t136) * t134;
t87 = -Icges(6,5) * t135 + (Icges(6,1) * t137 - Icges(6,4) * t136) * t134;
t38 = t105 * t86 + t106 * t87 + t85 * t184;
t5 = -t38 * t135 + (t147 * t19 + t149 * t20) * t134;
t21 = t107 * t60 + t108 * t62 + t58 * t183;
t22 = t107 * t61 + t108 * t63 + t59 * t183;
t39 = t107 * t86 + t108 * t87 + t85 * t183;
t6 = -t39 * t135 + (t147 * t21 + t149 * t22) * t134;
t202 = t6 * t183 + t5 * t184;
t201 = -t135 / 0.2e1;
t200 = t147 / 0.2e1;
t199 = -t149 / 0.2e1;
t198 = t149 / 0.2e1;
t119 = t134 * rSges(4,1) + t135 * rSges(4,2);
t197 = m(4) * t119;
t196 = pkin(3) * t135;
t195 = -pkin(3) + t133;
t194 = pkin(7) + t150;
t172 = pkin(3) * t182 + pkin(7) * t183;
t193 = -t172 + t203;
t159 = -t106 * rSges(6,1) - t105 * rSges(6,2);
t64 = rSges(6,3) * t184 - t159;
t88 = -t135 * rSges(6,3) + (rSges(6,1) * t137 - rSges(6,2) * t136) * t134;
t48 = t135 * t64 + t88 * t184;
t84 = t195 * t134 + t194 * t135;
t192 = -t84 - t88;
t191 = t136 * t86;
t92 = -Icges(5,6) * t135 + (Icges(5,4) * t148 - Icges(5,2) * t146) * t134;
t190 = t146 * t92;
t79 = t134 * t137 * t87;
t43 = -t134 * t191 - t135 * t85 + t79;
t189 = t43 * t135;
t188 = rSges(3,3) + qJ(2);
t121 = t134 * pkin(3) - t135 * pkin(7);
t94 = -t135 * rSges(5,3) + (rSges(5,1) * t148 - rSges(5,2) * t146) * t134;
t187 = -t121 - t94;
t185 = Icges(4,4) * t135;
t178 = t147 * t148;
t175 = t149 * t146;
t174 = t149 * t148;
t173 = t140 * (pkin(7) * t134 + t196) + t149 * t172;
t171 = t140 + t141;
t170 = -t121 + t192;
t111 = -t135 * t179 - t174;
t112 = t135 * t178 - t175;
t69 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t184;
t71 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t184;
t73 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t184;
t33 = -t135 * t69 + (-t146 * t71 + t148 * t73) * t134;
t91 = -Icges(5,3) * t135 + (Icges(5,5) * t148 - Icges(5,6) * t146) * t134;
t93 = -Icges(5,5) * t135 + (Icges(5,1) * t148 - Icges(5,4) * t146) * t134;
t41 = t111 * t92 + t112 * t93 + t91 * t184;
t169 = t41 / 0.2e1 + t33 / 0.2e1;
t113 = -t135 * t175 + t178;
t114 = t135 * t174 + t179;
t70 = Icges(5,5) * t114 + Icges(5,6) * t113 + Icges(5,3) * t183;
t72 = Icges(5,4) * t114 + Icges(5,2) * t113 + Icges(5,6) * t183;
t74 = Icges(5,1) * t114 + Icges(5,4) * t113 + Icges(5,5) * t183;
t34 = -t135 * t70 + (-t146 * t72 + t148 * t74) * t134;
t42 = t113 * t92 + t114 * t93 + t91 * t183;
t168 = t42 / 0.2e1 + t34 / 0.2e1;
t76 = t114 * rSges(5,1) + t113 * rSges(5,2) + rSges(5,3) * t183;
t167 = t184 / 0.2e1;
t166 = t183 / 0.2e1;
t27 = -t135 * t58 + (-t136 * t60 + t137 * t62) * t134;
t28 = -t135 * t59 + (-t136 * t61 + t137 * t63) * t134;
t165 = (t27 + t38) * t167 + (t28 + t39) * t166;
t144 = cos(pkin(9));
t132 = t144 * pkin(2) + pkin(1);
t145 = -pkin(6) - qJ(2);
t164 = t149 * t132 - t147 * t145;
t7 = -t189 + (t147 * t27 + t149 * t28) * t134;
t163 = -t135 * t7 + t202;
t12 = t20 * t147 - t19 * t149;
t13 = t22 * t147 - t21 * t149;
t162 = t12 * t167 + t13 * t166 + t5 * t199 + t6 * t200 + (t28 * t147 - t27 * t149) * t201;
t161 = rSges(4,1) * t135 - rSges(4,2) * t134;
t160 = -t112 * rSges(5,1) - t111 * rSges(5,2);
t155 = -Icges(4,2) * t134 + t185;
t154 = Icges(4,5) * t135 - Icges(4,6) * t134;
t153 = rSges(4,1) * t182 - rSges(4,2) * t183 + t147 * rSges(4,3);
t143 = sin(pkin(9));
t152 = rSges(3,1) * t144 - rSges(3,2) * t143 + pkin(1);
t123 = t149 * rSges(2,1) - t147 * rSges(2,2);
t122 = -t147 * rSges(2,1) - t149 * rSges(2,2);
t116 = Icges(4,6) * t135 + t205;
t96 = Icges(4,3) * t147 + t154 * t149;
t95 = -Icges(4,3) * t149 + t154 * t147;
t90 = t188 * t147 + t152 * t149;
t89 = -t152 * t147 + t188 * t149;
t83 = t134 * t148 * t93;
t81 = t153 + t164;
t80 = (rSges(4,3) - t145) * t149 + (-t132 - t161) * t147;
t78 = t187 * t149;
t77 = t187 * t147;
t75 = rSges(5,3) * t184 - t160;
t67 = -pkin(4) * t175 + (-t194 * t134 + t195 * t135) * t147;
t66 = t149 * t153 + (-t149 * rSges(4,3) + t161 * t147) * t147;
t56 = t64 * t183;
t55 = t164 + t76 + t172;
t54 = -t149 * t145 + (-t196 - t132 + (-rSges(5,3) - pkin(7)) * t134) * t147 + t160;
t53 = t170 * t149;
t52 = t170 * t147;
t51 = -t135 * t76 - t94 * t183;
t50 = t135 * t75 + t94 * t184;
t49 = -t135 * t65 - t88 * t183;
t47 = -t134 * t190 - t135 * t91 + t83;
t46 = t164 + t203;
t45 = (pkin(4) * t146 - t145) * t149 + (-t133 * t135 - t132 + (-rSges(6,3) + t150) * t134) * t147 + t159;
t44 = (-t147 * t76 + t149 * t75) * t134;
t40 = -t65 * t184 + t56;
t37 = t147 * t75 + t149 * t76 + t173;
t32 = t113 * t72 + t114 * t74 + t70 * t183;
t31 = t113 * t71 + t114 * t73 + t69 * t183;
t30 = t111 * t72 + t112 * t74 + t70 * t184;
t29 = t111 * t71 + t112 * t73 + t69 * t184;
t26 = -t193 * t135 + t192 * t183;
t25 = t135 * t67 + t84 * t184 + t48;
t18 = t56 + (-t193 * t147 + t149 * t67) * t134;
t17 = t193 * t149 + (t64 + t67) * t147 + t173;
t16 = t32 * t147 - t31 * t149;
t15 = t30 * t147 - t29 * t149;
t9 = -t42 * t135 + (t147 * t31 + t149 * t32) * t134;
t8 = -t41 * t135 + (t147 * t29 + t149 * t30) * t134;
t1 = [Icges(3,2) * t144 ^ 2 + Icges(2,3) + t79 + t83 + (Icges(3,1) * t143 + 0.2e1 * Icges(3,4) * t144) * t143 + (Icges(4,4) * t134 + Icges(4,2) * t135 - t85 - t91) * t135 + (Icges(4,1) * t134 + t185 - t190 - t191) * t134 + m(6) * (t45 ^ 2 + t46 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2) + m(4) * (t80 ^ 2 + t81 ^ 2) + m(3) * (t89 ^ 2 + t90 ^ 2) + m(2) * (t122 ^ 2 + t123 ^ 2); m(6) * (t147 * t45 - t149 * t46) + m(5) * (t147 * t54 - t149 * t55) + m(4) * (t147 * t80 - t149 * t81) + m(3) * (t147 * t89 - t149 * t90); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t171; m(6) * (t45 * t53 + t46 * t52) + m(5) * (t54 * t78 + t55 * t77) + (-t27 / 0.2e1 - t38 / 0.2e1 + t149 * t204 + (-Icges(4,6) * t149 + t155 * t147) * t201 - t80 * t197 + t116 * t198 - t169) * t149 + (t39 / 0.2e1 + t28 / 0.2e1 + t147 * t204 + t135 * (Icges(4,6) * t147 + t155 * t149) / 0.2e1 - t81 * t197 + t116 * t200 + t168) * t147; m(5) * (t78 * t147 - t77 * t149) + m(6) * (t53 * t147 - t52 * t149); m(6) * (t17 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t37 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t171 * t119 ^ 2 + t66 ^ 2) + (t140 * t96 + t13 + t16) * t147 + (-t141 * t95 - t12 - t15 + (-t147 * t95 + t149 * t96) * t147) * t149; (-t43 - t47) * t135 + m(6) * (t25 * t45 + t26 * t46) + m(5) * (t50 * t54 + t51 * t55) + (t169 * t147 + t168 * t149) * t134 + t165; m(5) * (t50 * t147 - t51 * t149) + m(6) * (t25 * t147 - t26 * t149); (t34 * t147 - t33 * t149) * t201 + t8 * t199 + t9 * t200 + (t15 * t200 + t16 * t198) * t134 + m(6) * (t17 * t18 + t25 * t53 + t26 * t52) + m(5) * (t37 * t44 + t50 * t78 + t51 * t77) + t162; (t47 * t135 - t7) * t135 + m(6) * (t18 ^ 2 + t25 ^ 2 + t26 ^ 2) + m(5) * (t44 ^ 2 + t50 ^ 2 + t51 ^ 2) + (t149 * t9 + t147 * t8 - t135 * (t147 * t33 + t149 * t34)) * t134 + t202; -t189 + m(6) * (t45 * t48 + t46 * t49) + t165; m(6) * (t48 * t147 - t49 * t149); m(6) * (t17 * t40 + t48 * t53 + t49 * t52) + t162; m(6) * (t18 * t40 + t25 * t48 + t26 * t49) + t163; m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + t163;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
