% Calculate joint inertia matrix for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:42
% EndTime: 2019-03-09 02:12:46
% DurationCPUTime: 1.85s
% Computational Cost: add. (3939->309), mult. (4952->443), div. (0->0), fcn. (5175->8), ass. (0->146)
t129 = pkin(9) + qJ(4);
t124 = cos(t129);
t202 = Icges(5,5) * t124;
t123 = sin(t129);
t201 = Icges(5,6) * t123;
t134 = -qJ(6) - pkin(8);
t200 = -rSges(7,3) + t134;
t199 = t202 / 0.2e1 - t201 / 0.2e1;
t138 = cos(qJ(5));
t136 = sin(qJ(5));
t71 = Icges(7,3) * t123 + (Icges(7,5) * t138 - Icges(7,6) * t136) * t124;
t72 = Icges(6,3) * t123 + (Icges(6,5) * t138 - Icges(6,6) * t136) * t124;
t75 = Icges(7,5) * t123 + (Icges(7,1) * t138 - Icges(7,4) * t136) * t124;
t76 = Icges(6,5) * t123 + (Icges(6,1) * t138 - Icges(6,4) * t136) * t124;
t198 = (t75 + t76) * t124 * t138 + (t71 + t72) * t123;
t73 = Icges(7,6) * t123 + (Icges(7,4) * t138 - Icges(7,2) * t136) * t124;
t74 = Icges(6,6) * t123 + (Icges(6,4) * t138 - Icges(6,2) * t136) * t124;
t197 = (-t73 - t74) * t136;
t137 = sin(qJ(1));
t168 = t124 * t137;
t139 = cos(qJ(1));
t162 = t138 * t139;
t164 = t137 * t136;
t93 = -t123 * t164 + t162;
t163 = t137 * t138;
t165 = t136 * t139;
t94 = t123 * t163 + t165;
t44 = Icges(7,5) * t94 + Icges(7,6) * t93 - Icges(7,3) * t168;
t48 = Icges(7,4) * t94 + Icges(7,2) * t93 - Icges(7,6) * t168;
t52 = Icges(7,1) * t94 + Icges(7,4) * t93 - Icges(7,5) * t168;
t11 = -t44 * t168 + t48 * t93 + t52 * t94;
t166 = t124 * t139;
t95 = t123 * t165 + t163;
t96 = -t123 * t162 + t164;
t45 = Icges(7,5) * t96 + Icges(7,6) * t95 + Icges(7,3) * t166;
t49 = Icges(7,4) * t96 + Icges(7,2) * t95 + Icges(7,6) * t166;
t53 = Icges(7,1) * t96 + Icges(7,4) * t95 + Icges(7,5) * t166;
t12 = -t45 * t168 + t49 * t93 + t53 * t94;
t46 = Icges(6,5) * t94 + Icges(6,6) * t93 - Icges(6,3) * t168;
t50 = Icges(6,4) * t94 + Icges(6,2) * t93 - Icges(6,6) * t168;
t54 = Icges(6,1) * t94 + Icges(6,4) * t93 - Icges(6,5) * t168;
t13 = -t46 * t168 + t50 * t93 + t54 * t94;
t47 = Icges(6,5) * t96 + Icges(6,6) * t95 + Icges(6,3) * t166;
t51 = Icges(6,4) * t96 + Icges(6,2) * t95 + Icges(6,6) * t166;
t55 = Icges(6,1) * t96 + Icges(6,4) * t95 + Icges(6,5) * t166;
t14 = -t47 * t168 + t51 * t93 + t55 * t94;
t26 = -t71 * t168 + t73 * t93 + t75 * t94;
t27 = -t72 * t168 + t74 * t93 + t76 * t94;
t196 = ((t12 + t14) * t139 + (-t11 - t13) * t137) * t124 + (t26 + t27) * t123;
t15 = t44 * t166 + t95 * t48 + t96 * t52;
t16 = t45 * t166 + t95 * t49 + t96 * t53;
t17 = t46 * t166 + t95 * t50 + t96 * t54;
t18 = t47 * t166 + t95 * t51 + t96 * t55;
t28 = t71 * t166 + t95 * t73 + t96 * t75;
t29 = t72 * t166 + t95 * t74 + t96 * t76;
t195 = ((t16 + t18) * t139 + (-t15 - t17) * t137) * t124 + (t28 + t29) * t123;
t21 = t123 * t44 + (-t136 * t48 + t138 * t52) * t124;
t23 = t123 * t46 + (-t136 * t50 + t138 * t54) * t124;
t194 = t21 + t23;
t22 = t123 * t45 + (-t136 * t49 + t138 * t53) * t124;
t24 = t123 * t47 + (-t136 * t51 + t138 * t55) * t124;
t193 = t22 + t24;
t121 = pkin(5) * t138 + pkin(4);
t170 = t123 * t137;
t192 = t94 * rSges(7,1) + t93 * rSges(7,2) + pkin(5) * t165 + t121 * t170 + t168 * t200;
t191 = (rSges(5,1) * t123 + rSges(5,2) * t124) * t139;
t130 = t137 ^ 2;
t131 = t139 ^ 2;
t187 = t137 / 0.2e1;
t186 = t139 / 0.2e1;
t102 = rSges(5,1) * t124 - rSges(5,2) * t123;
t185 = m(5) * t102;
t116 = t130 + t131;
t107 = m(5) * t116;
t184 = m(7) * t124;
t132 = sin(pkin(9));
t183 = pkin(3) * t132;
t182 = -pkin(8) - t134;
t181 = (t124 * t197 + t198) * t123;
t114 = pkin(4) * t170;
t90 = -pkin(8) * t168 + t114;
t180 = -t90 + t192;
t115 = t139 * t123 * pkin(4);
t153 = -t96 * rSges(7,1) - t95 * rSges(7,2);
t171 = t121 * t123;
t179 = pkin(5) * t164 + t115 + (t182 * t124 - t171) * t139 + rSges(7,3) * t166 - t153;
t176 = (rSges(7,1) * t138 - rSges(7,2) * t136 - pkin(4) + t121) * t124 + (t182 + rSges(7,3)) * t123;
t175 = t94 * rSges(6,1) + t93 * rSges(6,2);
t174 = rSges(4,3) + qJ(3);
t161 = t139 * pkin(1) + t137 * qJ(2);
t159 = rSges(5,1) * t170 + rSges(5,2) * t168 + t139 * rSges(5,3);
t126 = t139 * qJ(2);
t135 = -pkin(7) - qJ(3);
t158 = t137 * t135 + t139 * t183 + t126;
t157 = t124 * (-rSges(6,3) - pkin(8));
t156 = t176 * t137;
t155 = t107 + (m(4) + m(6) + m(7)) * t116;
t154 = -t96 * rSges(6,1) - t95 * rSges(6,2);
t133 = cos(pkin(9));
t152 = rSges(4,1) * t132 + rSges(4,2) * t133;
t19 = t180 * t123 + t124 * t156;
t20 = -t179 * t123 + t176 * t166;
t148 = t20 * t137 - t139 * t19;
t33 = (-pkin(5) * t136 - pkin(1)) * t137 + (t124 * t200 + t171) * t139 + t153 + t158;
t142 = -t135 * t139 + t137 * t183 + t161;
t34 = t142 + t192;
t147 = t137 * t33 - t139 * t34;
t104 = pkin(4) * t124 + pkin(8) * t123;
t97 = t137 * t104;
t39 = t97 + t156;
t40 = (-t104 - t176) * t139;
t146 = t39 * t137 - t139 * t40;
t143 = Icges(5,5) * t123 + Icges(5,6) * t124;
t141 = t22 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1 + t24 / 0.2e1;
t140 = -t27 / 0.2e1 - t26 / 0.2e1 - t23 / 0.2e1 - t21 / 0.2e1;
t111 = rSges(2,1) * t139 - t137 * rSges(2,2);
t110 = -t137 * rSges(2,1) - rSges(2,2) * t139;
t99 = -t201 + t202;
t92 = -rSges(3,2) * t139 + t137 * rSges(3,3) + t161;
t91 = rSges(3,3) * t139 + t126 + (rSges(3,2) - pkin(1)) * t137;
t85 = t139 * (pkin(8) * t166 - t115);
t80 = Icges(5,3) * t137 - t143 * t139;
t79 = Icges(5,3) * t139 + t143 * t137;
t78 = rSges(6,3) * t123 + (rSges(6,1) * t138 - rSges(6,2) * t136) * t124;
t68 = t152 * t137 + t174 * t139 + t161;
t67 = t126 + t152 * t139 + (-pkin(1) - t174) * t137;
t63 = t142 + t159;
t62 = t191 + (-rSges(5,3) - pkin(1)) * t137 + t158;
t61 = (-t104 - t78) * t139;
t60 = t137 * t78 + t97;
t59 = rSges(6,3) * t166 - t154;
t57 = -rSges(6,3) * t168 + t175;
t41 = -t137 * t159 + (t137 * rSges(5,3) - t191) * t139;
t38 = t137 * t157 + t114 + t142 + t175;
t37 = -t137 * pkin(1) + t139 * t157 + t115 + t154 + t158;
t36 = -t123 * t59 + t78 * t166;
t35 = t123 * t57 + t78 * t168;
t30 = (-t137 * t59 - t139 * t57) * t124;
t25 = t139 * t59 + t85 + (-t57 - t90) * t137;
t10 = (-t179 * t137 - t180 * t139) * t124;
t9 = t85 + t179 * t139 + (-t90 - t180) * t137;
t8 = t18 * t137 + t139 * t17;
t7 = t16 * t137 + t139 * t15;
t6 = t13 * t139 + t14 * t137;
t5 = t11 * t139 + t12 * t137;
t1 = [Icges(4,1) * t133 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t133 + Icges(4,2) * t132) * t132 + (Icges(5,1) * t124 + t197) * t124 + m(7) * (t33 ^ 2 + t34 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t91 ^ 2 + t92 ^ 2) + m(4) * (t67 ^ 2 + t68 ^ 2) + m(2) * (t110 ^ 2 + t111 ^ 2) + t198 + (-0.2e1 * Icges(5,4) * t124 + Icges(5,2) * t123) * t123; m(7) * t147 + m(6) * (t137 * t37 - t139 * t38) + m(5) * (t137 * t62 - t139 * t63) + m(3) * (t137 * t91 - t139 * t92) + m(4) * (t137 * t67 - t139 * t68); m(3) * t116 + t155; m(7) * (t137 * t34 + t139 * t33) + m(6) * (t137 * t38 + t139 * t37) + m(5) * (t137 * t63 + t139 * t62) + m(4) * (t137 * t68 + t139 * t67); 0; t155; m(7) * (t33 * t39 + t34 * t40) + m(6) * (t37 * t60 + t38 * t61) + (t139 * t199 - t63 * t185 + t99 * t186 - t140) * t139 + (t137 * t199 + t62 * t185 + t99 * t187 + t141) * t137; m(6) * (t60 * t137 - t139 * t61) + m(7) * t146 + t102 * t107; m(6) * (t61 * t137 + t139 * t60) + m(7) * (t40 * t137 + t139 * t39); m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t116 * t102 ^ 2 + t41 ^ 2) + (t130 * t80 + t7 + t8) * t137 + (t131 * t79 + t5 + t6 + (t137 * t79 + t139 * t80) * t137) * t139; m(7) * (t19 * t34 + t20 * t33) + m(6) * (t35 * t38 + t36 * t37) + (t140 * t137 + t141 * t139) * t124 + t181; m(6) * (t36 * t137 - t139 * t35) + m(7) * t148; m(6) * (t35 * t137 + t139 * t36) + m(7) * (t19 * t137 + t139 * t20); m(7) * (t10 * t9 + t19 * t40 + t20 * t39) + m(6) * (t25 * t30 + t35 * t61 + t36 * t60) + ((t8 / 0.2e1 + t7 / 0.2e1) * t139 + (-t6 / 0.2e1 - t5 / 0.2e1) * t137) * t124 + (t193 * t137 + t194 * t139) * t123 / 0.2e1 + t195 * t187 + t196 * t186; t181 * t123 + m(7) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(6) * (t30 ^ 2 + t35 ^ 2 + t36 ^ 2) + (t195 * t139 - t196 * t137 + (-t194 * t137 + t193 * t139) * t123) * t124; -t147 * t184; -t116 * t184; 0; m(7) * (t123 * t9 - t146 * t124); m(7) * (t123 * t10 - t148 * t124); m(7) * (t116 * t124 ^ 2 + t123 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
