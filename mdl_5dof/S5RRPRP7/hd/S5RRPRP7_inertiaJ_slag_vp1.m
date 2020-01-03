% Calculate joint inertia matrix for
% S5RRPRP7
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:56
% EndTime: 2019-12-31 20:00:01
% DurationCPUTime: 1.93s
% Computational Cost: add. (3588->312), mult. (4790->465), div. (0->0), fcn. (5128->8), ass. (0->152)
t216 = Icges(3,3) + Icges(4,3);
t133 = qJ(2) + pkin(8);
t128 = sin(t133);
t129 = cos(t133);
t138 = sin(qJ(2));
t141 = cos(qJ(2));
t215 = Icges(3,5) * t141 + Icges(4,5) * t129 - Icges(3,6) * t138 - Icges(4,6) * t128;
t140 = cos(qJ(4));
t142 = cos(qJ(1));
t171 = t140 * t142;
t137 = sin(qJ(4));
t139 = sin(qJ(1));
t173 = t139 * t137;
t102 = t129 * t173 + t171;
t172 = t139 * t140;
t174 = t137 * t142;
t103 = t129 * t172 - t174;
t211 = rSges(6,3) + qJ(5);
t213 = rSges(6,1) + pkin(4);
t214 = -t211 * t102 - t213 * t103;
t73 = -Icges(5,3) * t129 + (Icges(5,5) * t140 - Icges(5,6) * t137) * t128;
t74 = -Icges(6,2) * t129 + (Icges(6,4) * t140 + Icges(6,6) * t137) * t128;
t212 = -t73 - t74;
t179 = t128 * t139;
t44 = Icges(6,5) * t103 + Icges(6,6) * t179 + Icges(6,3) * t102;
t48 = Icges(6,4) * t103 + Icges(6,2) * t179 + Icges(6,6) * t102;
t52 = Icges(6,1) * t103 + Icges(6,4) * t179 + Icges(6,5) * t102;
t12 = t102 * t44 + t103 * t52 + t48 * t179;
t104 = t129 * t174 - t172;
t105 = t129 * t171 + t173;
t177 = t128 * t142;
t45 = Icges(6,5) * t105 + Icges(6,6) * t177 + Icges(6,3) * t104;
t49 = Icges(6,4) * t105 + Icges(6,2) * t177 + Icges(6,6) * t104;
t53 = Icges(6,1) * t105 + Icges(6,4) * t177 + Icges(6,5) * t104;
t13 = t102 * t45 + t103 * t53 + t49 * t179;
t46 = Icges(5,5) * t103 - Icges(5,6) * t102 + Icges(5,3) * t179;
t50 = Icges(5,4) * t103 - Icges(5,2) * t102 + Icges(5,6) * t179;
t54 = Icges(5,1) * t103 - Icges(5,4) * t102 + Icges(5,5) * t179;
t14 = -t102 * t50 + t103 * t54 + t46 * t179;
t47 = Icges(5,5) * t105 - Icges(5,6) * t104 + Icges(5,3) * t177;
t51 = Icges(5,4) * t105 - Icges(5,2) * t104 + Icges(5,6) * t177;
t55 = Icges(5,1) * t105 - Icges(5,4) * t104 + Icges(5,5) * t177;
t15 = -t102 * t51 + t103 * t55 + t47 * t179;
t72 = -Icges(6,6) * t129 + (Icges(6,5) * t140 + Icges(6,3) * t137) * t128;
t76 = -Icges(6,4) * t129 + (Icges(6,1) * t140 + Icges(6,5) * t137) * t128;
t26 = t102 * t72 + t103 * t76 + t74 * t179;
t75 = -Icges(5,6) * t129 + (Icges(5,4) * t140 - Icges(5,2) * t137) * t128;
t77 = -Icges(5,5) * t129 + (Icges(5,1) * t140 - Icges(5,4) * t137) * t128;
t27 = -t102 * t75 + t103 * t77 + t73 * t179;
t210 = (-t26 - t27) * t129 + ((t13 + t15) * t142 + (t12 + t14) * t139) * t128;
t16 = t104 * t44 + t105 * t52 + t48 * t177;
t17 = t104 * t45 + t105 * t53 + t49 * t177;
t18 = -t104 * t50 + t105 * t54 + t46 * t177;
t19 = -t104 * t51 + t105 * t55 + t47 * t177;
t28 = t104 * t72 + t105 * t76 + t74 * t177;
t29 = -t104 * t75 + t105 * t77 + t73 * t177;
t209 = (-t28 - t29) * t129 + ((t17 + t19) * t142 + (t16 + t18) * t139) * t128;
t208 = t128 / 0.2e1;
t207 = t129 / 0.2e1;
t206 = t138 / 0.2e1;
t205 = t141 / 0.2e1;
t20 = -t129 * t48 + (t137 * t44 + t140 * t52) * t128;
t22 = -t129 * t46 + (-t137 * t50 + t140 * t54) * t128;
t204 = -t20 - t22;
t21 = -t129 * t49 + (t137 * t45 + t140 * t53) * t128;
t23 = -t129 * t47 + (-t137 * t51 + t140 * t55) * t128;
t203 = t21 + t23;
t202 = -t215 * t139 + t216 * t142;
t201 = t216 * t139 + t215 * t142;
t180 = t128 * t137;
t200 = t72 * t180 + (t76 + t77) * t128 * t140;
t134 = t139 ^ 2;
t135 = t142 ^ 2;
t199 = -t129 / 0.2e1;
t115 = rSges(3,1) * t138 + rSges(3,2) * t141;
t196 = m(3) * t115;
t195 = pkin(2) * t138;
t194 = pkin(3) * t129;
t193 = t212 * t129 - t75 * t180 + t200;
t192 = rSges(6,2) * t179 - t214;
t191 = rSges(6,2) * t177 + t211 * t104 + t213 * t105;
t127 = pkin(2) * t141 + pkin(1);
t123 = t142 * t127;
t132 = t142 * pkin(6);
t136 = -qJ(3) - pkin(6);
t175 = t136 * t142;
t189 = t139 * (t175 + t132 + (-pkin(1) + t127) * t139) + t142 * (-pkin(1) * t142 + t123 + (-pkin(6) - t136) * t139);
t188 = -rSges(6,2) * t129 + (t211 * t137 + t213 * t140) * t128;
t187 = rSges(3,1) * t141;
t186 = rSges(3,2) * t138;
t185 = t142 * rSges(3,3);
t184 = Icges(3,4) * t138;
t183 = Icges(3,4) * t141;
t182 = Icges(4,4) * t128;
t181 = Icges(4,4) * t129;
t176 = t129 * t142;
t170 = pkin(3) * t176 + pkin(7) * t177;
t169 = t139 * rSges(3,3) + t142 * t187;
t168 = t134 + t135;
t59 = t105 * rSges(5,1) - t104 * rSges(5,2) + rSges(5,3) * t177;
t167 = Icges(3,5) * t206 + Icges(4,5) * t208 + Icges(3,6) * t205 + Icges(4,6) * t207;
t166 = -rSges(4,1) * t128 - rSges(4,2) * t129 - t195;
t165 = -pkin(3) * t128 + pkin(7) * t129 - t195;
t164 = -t127 - t194;
t163 = t134 * (pkin(7) * t128 + t194) + t142 * t170 + t189;
t162 = -t139 * t136 + t123;
t79 = -rSges(5,3) * t129 + (rSges(5,1) * t140 - rSges(5,2) * t137) * t128;
t161 = t165 - t79;
t160 = -t186 + t187;
t159 = rSges(4,1) * t129 - rSges(4,2) * t128;
t158 = -t103 * rSges(5,1) + t102 * rSges(5,2);
t153 = t165 - t188;
t152 = Icges(3,1) * t141 - t184;
t151 = Icges(4,1) * t129 - t182;
t150 = -Icges(3,2) * t138 + t183;
t149 = -Icges(4,2) * t128 + t181;
t146 = rSges(4,1) * t176 - rSges(4,2) * t177 + t139 * rSges(4,3);
t145 = t162 + t170;
t144 = t27 / 0.2e1 + t26 / 0.2e1 + t20 / 0.2e1 + t22 / 0.2e1;
t143 = t28 / 0.2e1 + t23 / 0.2e1 + t21 / 0.2e1 + t29 / 0.2e1;
t117 = rSges(2,1) * t142 - t139 * rSges(2,2);
t116 = -t139 * rSges(2,1) - rSges(2,2) * t142;
t81 = t166 * t142;
t80 = t166 * t139;
t71 = t139 * pkin(6) + (pkin(1) - t186) * t142 + t169;
t70 = t185 + t132 + (-pkin(1) - t160) * t139;
t64 = t146 + t162;
t63 = (rSges(4,3) - t136) * t142 + (-t127 - t159) * t139;
t60 = t142 * (-t142 * t186 + t169) + (t160 * t139 - t185) * t139;
t57 = rSges(5,3) * t179 - t158;
t43 = t161 * t142;
t42 = t161 * t139;
t41 = t153 * t142;
t40 = t153 * t139;
t39 = t145 + t59;
t38 = -t175 + ((-rSges(5,3) - pkin(7)) * t128 + t164) * t139 + t158;
t37 = -t129 * t59 - t79 * t177;
t36 = t129 * t57 + t79 * t179;
t33 = t142 * t146 + (-t142 * rSges(4,3) + t159 * t139) * t139 + t189;
t32 = (-t139 * t59 + t142 * t57) * t128;
t31 = t145 + t191;
t30 = -t175 + ((-rSges(6,2) - pkin(7)) * t128 + t164) * t139 + t214;
t25 = -t191 * t129 - t188 * t177;
t24 = t192 * t129 + t188 * t179;
t11 = t139 * t57 + t142 * t59 + t163;
t10 = (-t191 * t139 + t192 * t142) * t128;
t9 = t192 * t139 + t191 * t142 + t163;
t8 = t19 * t139 - t142 * t18;
t7 = t17 * t139 - t142 * t16;
t6 = t15 * t139 - t14 * t142;
t5 = -t12 * t142 + t13 * t139;
t1 = [t141 * (Icges(3,2) * t141 + t184) + t138 * (Icges(3,1) * t138 + t183) + Icges(2,3) + (Icges(4,1) * t128 - t137 * t75 + t181) * t128 + (Icges(4,2) * t129 + t182 + t212) * t129 + m(6) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t38 ^ 2 + t39 ^ 2) + m(4) * (t63 ^ 2 + t64 ^ 2) + m(3) * (t70 ^ 2 + t71 ^ 2) + m(2) * (t116 ^ 2 + t117 ^ 2) + t200; m(6) * (t30 * t41 + t31 * t40) + m(5) * (t38 * t43 + t39 * t42) + m(4) * (t63 * t81 + t64 * t80) + (-t138 * (-Icges(3,5) * t142 + t152 * t139) / 0.2e1 - t141 * (-Icges(3,6) * t142 + t150 * t139) / 0.2e1 - t128 * (-Icges(4,5) * t142 + t151 * t139) / 0.2e1 + (-Icges(4,6) * t142 + t149 * t139) * t199 - t70 * t196 + t167 * t142 - t144) * t142 + (t167 * t139 - t71 * t196 + (Icges(3,6) * t139 + t150 * t142) * t205 + (Icges(3,5) * t139 + t152 * t142) * t206 + (Icges(4,6) * t139 + t149 * t142) * t207 + (Icges(4,5) * t139 + t151 * t142) * t208 + t143) * t139; m(6) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(4) * (t33 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(3) * (t168 * t115 ^ 2 + t60 ^ 2) + (t202 * t135 - t5 - t6) * t142 + (t7 + t8 + t201 * t134 + (t202 * t139 + t201 * t142) * t142) * t139; m(6) * (t139 * t30 - t142 * t31) + m(5) * (t139 * t38 - t142 * t39) + m(4) * (t139 * t63 - t142 * t64); m(6) * (t139 * t41 - t142 * t40) + m(5) * (t139 * t43 - t142 * t42) + m(4) * (t139 * t81 - t142 * t80); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t168; -t193 * t129 + m(6) * (t24 * t30 + t25 * t31) + m(5) * (t36 * t38 + t37 * t39) + (t144 * t139 + t143 * t142) * t128; m(6) * (t10 * t9 + t24 * t41 + t25 * t40) + m(5) * (t11 * t32 + t36 * t43 + t37 * t42) + ((t7 / 0.2e1 + t8 / 0.2e1) * t142 + (t5 / 0.2e1 + t6 / 0.2e1) * t139) * t128 + (t203 * t139 + t204 * t142) * t199 + t209 * t139 / 0.2e1 - t210 * t142 / 0.2e1; m(5) * (t36 * t139 - t142 * t37) + m(6) * (t24 * t139 - t142 * t25); m(6) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t32 ^ 2 + t36 ^ 2 + t37 ^ 2) + t193 * t129 ^ 2 + (t209 * t142 + t210 * t139 + (t204 * t139 - t203 * t142) * t129) * t128; m(6) * (t102 * t31 + t104 * t30); m(6) * (t102 * t40 + t104 * t41 + t9 * t180); m(6) * (-t102 * t142 + t104 * t139); m(6) * (t10 * t180 + t102 * t25 + t104 * t24); m(6) * (t128 ^ 2 * t137 ^ 2 + t102 ^ 2 + t104 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
