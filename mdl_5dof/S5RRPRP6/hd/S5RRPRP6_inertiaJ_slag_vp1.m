% Calculate joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:56:52
% DurationCPUTime: 2.03s
% Computational Cost: add. (3775->321), mult. (4796->470), div. (0->0), fcn. (5061->8), ass. (0->155)
t137 = -qJ(5) - pkin(7);
t221 = rSges(6,3) - t137;
t139 = sin(qJ(4));
t134 = qJ(2) + pkin(8);
t129 = sin(t134);
t130 = cos(t134);
t142 = cos(qJ(4));
t74 = -Icges(6,6) * t130 + (Icges(6,4) * t142 - Icges(6,2) * t139) * t129;
t75 = -Icges(5,6) * t130 + (Icges(5,4) * t142 - Icges(5,2) * t139) * t129;
t220 = (t74 + t75) * t139;
t219 = Icges(3,3) + Icges(4,3);
t140 = sin(qJ(2));
t143 = cos(qJ(2));
t218 = Icges(3,5) * t143 + Icges(4,5) * t130 - Icges(3,6) * t140 - Icges(4,6) * t129;
t72 = -Icges(6,3) * t130 + (Icges(6,5) * t142 - Icges(6,6) * t139) * t129;
t73 = -Icges(5,3) * t130 + (Icges(5,5) * t142 - Icges(5,6) * t139) * t129;
t217 = t72 + t73;
t144 = cos(qJ(1));
t173 = t142 * t144;
t141 = sin(qJ(1));
t175 = t141 * t139;
t100 = -t130 * t175 - t173;
t174 = t141 * t142;
t176 = t139 * t144;
t101 = t130 * t174 - t176;
t181 = t129 * t141;
t46 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t181;
t50 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t181;
t54 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t181;
t12 = t100 * t50 + t101 * t54 + t46 * t181;
t102 = -t130 * t176 + t174;
t103 = t130 * t173 + t175;
t179 = t129 * t144;
t47 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t179;
t51 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t179;
t55 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t179;
t13 = t100 * t51 + t101 * t55 + t47 * t181;
t48 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t181;
t52 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t181;
t56 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t181;
t14 = t100 * t52 + t101 * t56 + t48 * t181;
t49 = Icges(5,5) * t103 + Icges(5,6) * t102 + Icges(5,3) * t179;
t53 = Icges(5,4) * t103 + Icges(5,2) * t102 + Icges(5,6) * t179;
t57 = Icges(5,1) * t103 + Icges(5,4) * t102 + Icges(5,5) * t179;
t15 = t100 * t53 + t101 * t57 + t49 * t181;
t76 = -Icges(6,5) * t130 + (Icges(6,1) * t142 - Icges(6,4) * t139) * t129;
t26 = t100 * t74 + t101 * t76 + t72 * t181;
t77 = -Icges(5,5) * t130 + (Icges(5,1) * t142 - Icges(5,4) * t139) * t129;
t27 = t100 * t75 + t101 * t77 + t73 * t181;
t215 = (-t26 - t27) * t130 + ((t13 + t15) * t144 + (t12 + t14) * t141) * t129;
t16 = t102 * t50 + t103 * t54 + t46 * t179;
t17 = t102 * t51 + t103 * t55 + t47 * t179;
t18 = t102 * t52 + t103 * t56 + t48 * t179;
t19 = t102 * t53 + t103 * t57 + t49 * t179;
t28 = t102 * t74 + t103 * t76 + t72 * t179;
t29 = t102 * t75 + t103 * t77 + t73 * t179;
t214 = (-t28 - t29) * t130 + ((t17 + t19) * t144 + (t16 + t18) * t141) * t129;
t213 = t129 / 0.2e1;
t212 = t130 / 0.2e1;
t211 = t140 / 0.2e1;
t210 = t143 / 0.2e1;
t22 = -t130 * t46 + (-t139 * t50 + t142 * t54) * t129;
t24 = -t130 * t48 + (-t139 * t52 + t142 * t56) * t129;
t209 = -t22 - t24;
t23 = -t130 * t47 + (-t139 * t51 + t142 * t55) * t129;
t25 = -t130 * t49 + (-t139 * t53 + t142 * t57) * t129;
t208 = t23 + t25;
t207 = (t76 + t77) * t129 * t142;
t206 = -t218 * t141 + t219 * t144;
t205 = t219 * t141 + t218 * t144;
t126 = pkin(4) * t142 + pkin(3);
t178 = t130 * t144;
t204 = t103 * rSges(6,1) + t102 * rSges(6,2) + pkin(4) * t175 + t126 * t178 + t221 * t179;
t203 = t130 ^ 2;
t135 = t141 ^ 2;
t136 = t144 ^ 2;
t202 = -t130 / 0.2e1;
t114 = rSges(3,1) * t140 + rSges(3,2) * t143;
t199 = m(3) * t114;
t198 = pkin(2) * t140;
t197 = pkin(3) * t130;
t196 = -pkin(3) + t126;
t195 = pkin(7) + t137;
t194 = t129 * t220 + t217 * t130 - t207;
t160 = -t101 * rSges(6,1) - t100 * rSges(6,2);
t193 = -pkin(4) * t176 + (-t195 * t129 + t196 * t130) * t141 + rSges(6,3) * t181 - t160;
t172 = pkin(3) * t178 + pkin(7) * t179;
t192 = -t172 + t204;
t191 = (t195 - rSges(6,3)) * t130 + (rSges(6,1) * t142 - rSges(6,2) * t139 + t196) * t129;
t127 = pkin(2) * t143 + pkin(1);
t122 = t144 * t127;
t133 = t144 * pkin(6);
t138 = -qJ(3) - pkin(6);
t177 = t138 * t144;
t190 = t141 * (t177 + t133 + (-pkin(1) + t127) * t141) + t144 * (-pkin(1) * t144 + t122 + (-pkin(6) - t138) * t141);
t189 = rSges(3,1) * t143;
t188 = rSges(3,2) * t140;
t187 = t144 * rSges(3,3);
t186 = Icges(3,4) * t140;
t185 = Icges(3,4) * t143;
t184 = Icges(4,4) * t129;
t183 = Icges(4,4) * t130;
t171 = t141 * rSges(3,3) + t144 * t189;
t170 = t135 + t136;
t61 = t103 * rSges(5,1) + t102 * rSges(5,2) + rSges(5,3) * t179;
t169 = Icges(3,5) * t211 + Icges(4,5) * t213 + Icges(3,6) * t210 + Icges(4,6) * t212;
t168 = -rSges(4,1) * t129 - rSges(4,2) * t130 - t198;
t167 = -pkin(3) * t129 + pkin(7) * t130 - t198;
t166 = t135 * (pkin(7) * t129 + t197) + t144 * t172 + t190;
t165 = -t141 * t138 + t122;
t79 = -rSges(5,3) * t130 + (rSges(5,1) * t142 - rSges(5,2) * t139) * t129;
t164 = t167 - t79;
t163 = -t188 + t189;
t162 = rSges(4,1) * t130 - rSges(4,2) * t129;
t161 = -t101 * rSges(5,1) - t100 * rSges(5,2);
t155 = t167 - t191;
t154 = Icges(3,1) * t143 - t186;
t153 = Icges(4,1) * t130 - t184;
t152 = -Icges(3,2) * t140 + t185;
t151 = -Icges(4,2) * t129 + t183;
t148 = rSges(4,1) * t178 - rSges(4,2) * t179 + t141 * rSges(4,3);
t146 = t24 / 0.2e1 + t22 / 0.2e1 + t27 / 0.2e1 + t26 / 0.2e1;
t145 = t29 / 0.2e1 + t28 / 0.2e1 + t25 / 0.2e1 + t23 / 0.2e1;
t116 = rSges(2,1) * t144 - t141 * rSges(2,2);
t115 = -t141 * rSges(2,1) - rSges(2,2) * t144;
t81 = t168 * t144;
t80 = t168 * t141;
t71 = t141 * pkin(6) + (pkin(1) - t188) * t144 + t171;
t70 = t187 + t133 + (-pkin(1) - t163) * t141;
t64 = t148 + t165;
t63 = (rSges(4,3) - t138) * t144 + (-t127 - t162) * t141;
t62 = t144 * (-t144 * t188 + t171) + (t141 * t163 - t187) * t141;
t59 = rSges(5,3) * t181 - t161;
t43 = t164 * t144;
t42 = t164 * t141;
t41 = t165 + t61 + t172;
t40 = -t177 + (-t197 - t127 + (-rSges(5,3) - pkin(7)) * t129) * t141 + t161;
t39 = -t130 * t61 - t79 * t179;
t38 = t130 * t59 + t79 * t181;
t37 = t155 * t144;
t36 = t155 * t141;
t35 = t165 + t204;
t34 = (pkin(4) * t139 - t138) * t144 + (-t126 * t130 - t221 * t129 - t127) * t141 + t160;
t31 = t144 * t148 + (-t144 * rSges(4,3) + t141 * t162) * t141 + t190;
t30 = (-t141 * t61 + t144 * t59) * t129;
t21 = -t192 * t130 - t191 * t179;
t20 = t193 * t130 + t191 * t181;
t11 = t141 * t59 + t144 * t61 + t166;
t10 = (-t192 * t141 + t193 * t144) * t129;
t9 = t193 * t141 + t192 * t144 + t166;
t8 = t19 * t141 - t144 * t18;
t7 = t17 * t141 - t144 * t16;
t6 = -t14 * t144 + t15 * t141;
t5 = -t12 * t144 + t13 * t141;
t1 = [t143 * (Icges(3,2) * t143 + t186) + t140 * (Icges(3,1) * t140 + t185) + Icges(2,3) + (Icges(4,2) * t130 + t184 - t217) * t130 + (Icges(4,1) * t129 + t183 - t220) * t129 + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2) + m(4) * (t63 ^ 2 + t64 ^ 2) + m(3) * (t70 ^ 2 + t71 ^ 2) + m(2) * (t115 ^ 2 + t116 ^ 2) + t207; m(6) * (t34 * t37 + t35 * t36) + m(5) * (t40 * t43 + t41 * t42) + m(4) * (t63 * t81 + t64 * t80) + (-t70 * t199 - t140 * (-Icges(3,5) * t144 + t154 * t141) / 0.2e1 - t143 * (-Icges(3,6) * t144 + t152 * t141) / 0.2e1 - t129 * (-Icges(4,5) * t144 + t153 * t141) / 0.2e1 + (-Icges(4,6) * t144 + t151 * t141) * t202 + t169 * t144 - t146) * t144 + (t169 * t141 - t71 * t199 + (Icges(3,6) * t141 + t152 * t144) * t210 + (Icges(3,5) * t141 + t154 * t144) * t211 + (Icges(4,6) * t141 + t151 * t144) * t212 + (Icges(4,5) * t141 + t153 * t144) * t213 + t145) * t141; m(5) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2 + t9 ^ 2) + m(4) * (t31 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(3) * (t170 * t114 ^ 2 + t62 ^ 2) + (t206 * t136 - t5 - t6) * t144 + (t7 + t8 + t205 * t135 + (t206 * t141 + t205 * t144) * t144) * t141; m(6) * (t141 * t34 - t144 * t35) + m(5) * (t141 * t40 - t144 * t41) + m(4) * (t141 * t63 - t144 * t64); m(5) * (t141 * t43 - t144 * t42) + m(6) * (t141 * t37 - t144 * t36) + m(4) * (t141 * t81 - t144 * t80); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t170; t194 * t130 + m(6) * (t20 * t34 + t21 * t35) + m(5) * (t38 * t40 + t39 * t41) + (t141 * t146 + t144 * t145) * t129; m(5) * (t11 * t30 + t38 * t43 + t39 * t42) + m(6) * (t10 * t9 + t20 * t37 + t21 * t36) + ((t8 / 0.2e1 + t7 / 0.2e1) * t144 + (t5 / 0.2e1 + t6 / 0.2e1) * t141) * t129 + (t208 * t141 + t209 * t144) * t202 + t214 * t141 / 0.2e1 - t215 * t144 / 0.2e1; m(5) * (t38 * t141 - t144 * t39) + m(6) * (t20 * t141 - t144 * t21); m(6) * (t10 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(5) * (t30 ^ 2 + t38 ^ 2 + t39 ^ 2) - t194 * t203 + (t214 * t144 + t215 * t141 + (t209 * t141 - t208 * t144) * t130) * t129; m(6) * (t141 * t35 + t144 * t34) * t129; m(6) * (-t130 * t9 + (t141 * t36 + t144 * t37) * t129); 0; m(6) * (-t130 * t10 + (t141 * t21 + t144 * t20) * t129); m(6) * (t170 * t129 ^ 2 + t203);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
