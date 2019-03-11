% Calculate joint inertia matrix for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:16
% EndTime: 2019-03-09 02:28:18
% DurationCPUTime: 1.46s
% Computational Cost: add. (3239->285), mult. (4419->421), div. (0->0), fcn. (4549->8), ass. (0->150)
t132 = qJ(4) + qJ(5);
t124 = sin(t132);
t125 = cos(t132);
t133 = sin(qJ(6));
t136 = cos(qJ(6));
t63 = t124 * rSges(7,3) + (rSges(7,1) * t136 - rSges(7,2) * t133) * t125;
t198 = pkin(5) * t125 + pkin(9) * t124 + t63;
t135 = sin(qJ(1));
t130 = t135 ^ 2;
t138 = cos(qJ(1));
t178 = Icges(6,4) * t124;
t146 = Icges(6,2) * t125 + t178;
t69 = -Icges(6,6) * t135 + t146 * t138;
t177 = Icges(6,4) * t125;
t148 = Icges(6,1) * t124 + t177;
t71 = -Icges(6,5) * t135 + t148 * t138;
t155 = -t124 * t71 - t125 * t69;
t68 = Icges(6,6) * t138 + t146 * t135;
t70 = Icges(6,5) * t138 + t148 * t135;
t156 = t124 * t70 + t125 * t68;
t144 = Icges(6,5) * t124 + Icges(6,6) * t125;
t66 = Icges(6,3) * t138 + t144 * t135;
t67 = -Icges(6,3) * t135 + t144 * t138;
t174 = t125 * t138;
t175 = t125 * t135;
t169 = t136 * t138;
t173 = t133 * t135;
t90 = -t124 * t173 + t169;
t170 = t135 * t136;
t172 = t133 * t138;
t91 = t124 * t170 + t172;
t43 = Icges(7,5) * t91 + Icges(7,6) * t90 - Icges(7,3) * t175;
t45 = Icges(7,4) * t91 + Icges(7,2) * t90 - Icges(7,6) * t175;
t47 = Icges(7,1) * t91 + Icges(7,4) * t90 - Icges(7,5) * t175;
t92 = -t124 * t172 - t170;
t93 = t124 * t169 - t173;
t14 = -t43 * t174 + t45 * t92 + t47 * t93;
t44 = Icges(7,5) * t93 + Icges(7,6) * t92 - Icges(7,3) * t174;
t46 = Icges(7,4) * t93 + Icges(7,2) * t92 - Icges(7,6) * t174;
t48 = Icges(7,1) * t93 + Icges(7,4) * t92 - Icges(7,5) * t174;
t15 = -t44 * t174 + t46 * t92 + t48 * t93;
t9 = -t135 * t15 + t138 * t14;
t197 = -t130 * t67 - (t156 * t138 + (t155 - t66) * t135) * t138 - t9;
t131 = t138 ^ 2;
t196 = -t135 / 0.2e1;
t195 = t138 / 0.2e1;
t194 = -rSges(5,3) - pkin(7);
t193 = rSges(7,3) + pkin(9);
t12 = -t43 * t175 + t45 * t90 + t47 * t91;
t13 = -t44 * t175 + t46 * t90 + t48 * t91;
t8 = t12 * t138 - t13 * t135;
t192 = (t131 * t66 + t8 + (t155 * t135 + (t156 - t67) * t138) * t135) * t138;
t116 = t130 + t131;
t104 = m(5) * t116;
t103 = m(6) * t116;
t191 = pkin(4) * t135;
t190 = pkin(5) * t124;
t189 = -pkin(1) - qJ(3);
t159 = -rSges(7,1) * t91 - rSges(7,2) * t90;
t49 = -rSges(7,3) * t175 - t159;
t188 = -t49 - (-pkin(9) * t125 + t190) * t135;
t176 = t124 * t138;
t115 = pkin(5) * t176;
t185 = t93 * rSges(7,1) + t92 * rSges(7,2);
t50 = -rSges(7,3) * t174 + t185;
t187 = pkin(9) * t174 - t115 - t50;
t60 = Icges(7,3) * t124 + (Icges(7,5) * t136 - Icges(7,6) * t133) * t125;
t62 = Icges(7,5) * t124 + (Icges(7,1) * t136 - Icges(7,4) * t133) * t125;
t186 = t125 * t136 * t62 + t124 * t60;
t52 = t198 * t135;
t53 = t198 * t138;
t61 = Icges(7,6) * t124 + (Icges(7,4) * t136 - Icges(7,2) * t133) * t125;
t184 = t133 * t61;
t183 = t138 * rSges(6,3);
t20 = t124 * t43 + (-t133 * t45 + t136 * t47) * t125;
t182 = t20 * t138;
t21 = t124 * t44 + (-t133 * t46 + t136 * t48) * t125;
t181 = t21 * t135;
t134 = sin(qJ(4));
t180 = Icges(5,4) * t134;
t137 = cos(qJ(4));
t179 = Icges(5,4) * t137;
t171 = t134 * t138;
t168 = t137 * t138;
t167 = rSges(5,1) * t171 + rSges(5,2) * t168;
t139 = -pkin(8) - pkin(7);
t166 = pkin(4) * t171 + t135 * t139;
t123 = t138 * t139;
t128 = t138 * qJ(2);
t165 = t123 + t128;
t164 = t138 * pkin(1) + t135 * qJ(2);
t163 = t138 * qJ(3) + t164;
t25 = -t60 * t175 + t61 * t90 + t62 * t91;
t3 = t124 * t25 + (-t12 * t135 - t13 * t138) * t125;
t26 = -t60 * t174 + t61 * t92 + t62 * t93;
t4 = t124 * t26 + (-t135 * t14 - t138 * t15) * t125;
t162 = t3 * t195 + t4 * t196 + t124 * (-t181 + t182) / 0.2e1 - t8 * t175 / 0.2e1 - t9 * t174 / 0.2e1;
t161 = t103 + t104 + (m(4) + m(7)) * t116;
t160 = -pkin(4) * t134 + t189;
t73 = rSges(6,1) * t176 + rSges(6,2) * t174 - rSges(6,3) * t135;
t158 = -rSges(5,1) * t134 - rSges(5,2) * t137;
t157 = rSges(6,1) * t124 + rSges(6,2) * t125;
t98 = -Icges(6,2) * t124 + t177;
t99 = Icges(6,1) * t125 - t178;
t154 = t124 * t99 + t125 * t98;
t100 = rSges(6,1) * t125 - rSges(6,2) * t124;
t119 = t137 * t191;
t64 = t100 * t135 + t119;
t121 = pkin(4) * t168;
t65 = t100 * t138 + t121;
t151 = t135 * t64 + t138 * t65;
t150 = t163 + t166;
t149 = Icges(5,1) * t134 + t179;
t147 = Icges(5,2) * t137 + t180;
t145 = Icges(5,5) * t134 + Icges(5,6) * t137;
t143 = t135 * t197 + t192;
t54 = t128 + t194 * t138 + (t158 + t189) * t135;
t55 = t194 * t135 + t163 + t167;
t142 = m(5) * (t135 * t55 + t138 * t54);
t41 = -t183 + (-t157 + t160) * t135 + t165;
t42 = t150 + t73;
t141 = m(6) * (t135 * t42 + t138 * t41);
t97 = Icges(6,5) * t125 - Icges(6,6) * t124;
t140 = t182 / 0.2e1 - t181 / 0.2e1 + (-t124 * t69 + t125 * t71 - t135 * t97 + t154 * t138 + t26) * t196 + (-t124 * t68 + t125 * t70 + t154 * t135 + t138 * t97 + t25) * t195;
t112 = rSges(2,1) * t138 - rSges(2,2) * t135;
t111 = rSges(5,1) * t137 - rSges(5,2) * t134;
t110 = -rSges(2,1) * t135 - rSges(2,2) * t138;
t89 = -pkin(7) * t138 + t134 * t191 - t123;
t88 = pkin(7) * t135 + t166;
t85 = -rSges(3,2) * t138 + rSges(3,3) * t135 + t164;
t84 = rSges(3,3) * t138 + t128 + (rSges(3,2) - pkin(1)) * t135;
t77 = -Icges(5,3) * t135 + t145 * t138;
t76 = Icges(5,3) * t138 + t145 * t135;
t75 = rSges(4,2) * t135 + rSges(4,3) * t138 + t163;
t74 = rSges(4,2) * t138 + t128 + (-rSges(4,3) + t189) * t135;
t72 = t157 * t135 + t183;
t51 = t158 * t130 - t138 * t167;
t40 = t121 + t53;
t39 = t119 + t52;
t38 = -t135 * t72 - t138 * t73;
t33 = t124 * t50 + t63 * t174;
t32 = -t124 * t49 - t63 * t175;
t31 = -t193 * t174 + t115 + t150 + t185;
t30 = (t193 * t125 + t160 - t190) * t135 + t159 + t165;
t29 = (-t73 - t88) * t138 + (-t72 - t89) * t135;
t28 = (-t125 * t184 + t186) * t124;
t27 = (t135 * t50 - t138 * t49) * t125;
t22 = t188 * t135 + t187 * t138;
t17 = (-t88 + t187) * t138 + (-t89 + t188) * t135;
t1 = [-t134 * (-Icges(5,2) * t134 + t179) + t137 * (Icges(5,1) * t137 - t180) - t124 * t98 + Icges(3,1) + Icges(4,1) + Icges(2,3) + (t99 - t184) * t125 + m(7) * (t30 ^ 2 + t31 ^ 2) + m(6) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2) + m(4) * (t74 ^ 2 + t75 ^ 2) + m(3) * (t84 ^ 2 + t85 ^ 2) + m(2) * (t110 ^ 2 + t112 ^ 2) + t186; m(7) * (t135 * t30 - t138 * t31) + m(6) * (t135 * t41 - t138 * t42) + m(5) * (t135 * t54 - t138 * t55) + m(4) * (t135 * t74 - t138 * t75) + m(3) * (t135 * t84 - t138 * t85); m(3) * t116 + t161; m(7) * (t135 * t31 + t138 * t30) + t141 + t142 + m(4) * (t135 * t75 + t138 * t74); 0; t161; (-t134 * (Icges(5,6) * t138 + t147 * t135) + t137 * (Icges(5,5) * t138 + t149 * t135)) * t195 + (-t134 * (-Icges(5,6) * t135 + t147 * t138) + t137 * (-Icges(5,5) * t135 + t149 * t138)) * t196 + m(7) * (t30 * t40 + t31 * t39) + m(6) * (t41 * t65 + t42 * t64) + t111 * t142 + (t131 / 0.2e1 + t130 / 0.2e1) * (Icges(5,5) * t137 - Icges(5,6) * t134) + t140; m(6) * (t135 * t65 - t138 * t64) + m(7) * (t135 * t40 - t138 * t39); m(6) * t151 + m(7) * (t135 * t39 + t138 * t40) + t111 * t104; m(7) * (t17 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t29 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t116 * t111 ^ 2 + t51 ^ 2) + t138 * t131 * t76 + t192 + (-t130 * t77 + (t135 * t76 - t138 * t77) * t138 + t197) * t135; m(7) * (t30 * t53 + t31 * t52) + t100 * t141 + t140; m(7) * (t135 * t53 - t138 * t52); m(7) * (t135 * t52 + t138 * t53) + t100 * t103; m(7) * (t17 * t22 + t39 * t52 + t40 * t53) + m(6) * (t151 * t100 + t29 * t38) + t143; m(6) * (t116 * t100 ^ 2 + t38 ^ 2) + m(7) * (t22 ^ 2 + t52 ^ 2 + t53 ^ 2) + t143; m(7) * (t30 * t32 + t31 * t33) + t28 + ((-t26 / 0.2e1 - t21 / 0.2e1) * t138 + (-t25 / 0.2e1 - t20 / 0.2e1) * t135) * t125; m(7) * (t135 * t32 - t138 * t33); m(7) * (t135 * t33 + t138 * t32); m(7) * (t17 * t27 + t32 * t40 + t33 * t39) + t162; m(7) * (t22 * t27 + t32 * t53 + t33 * t52) + t162; t124 * t28 + m(7) * (t27 ^ 2 + t32 ^ 2 + t33 ^ 2) + (-t138 * t4 - t135 * t3 + t124 * (-t135 * t20 - t138 * t21)) * t125;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
