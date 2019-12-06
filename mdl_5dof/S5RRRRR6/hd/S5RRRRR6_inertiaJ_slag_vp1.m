% Calculate joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:03
% EndTime: 2019-12-05 19:00:06
% DurationCPUTime: 1.12s
% Computational Cost: add. (4193->236), mult. (2794->333), div. (0->0), fcn. (2510->10), ass. (0->133)
t124 = qJ(1) + qJ(2);
t117 = sin(t124);
t119 = cos(t124);
t189 = t117 * t119;
t114 = t117 ^ 2;
t115 = t119 ^ 2;
t129 = -pkin(8) - pkin(7);
t123 = qJ(3) + qJ(4);
t120 = qJ(5) + t123;
t111 = sin(t120);
t112 = cos(t120);
t133 = Icges(6,5) * t112 - Icges(6,6) * t111;
t45 = Icges(6,3) * t119 - t133 * t117;
t46 = Icges(6,3) * t117 + t133 * t119;
t188 = t117 * (t114 * t46 + t45 * t189) + t119 * (t115 * t45 + t46 * t189);
t125 = sin(qJ(3));
t127 = cos(qJ(3));
t98 = t125 * rSges(4,1) + t127 * rSges(4,2);
t187 = m(4) * t98;
t116 = sin(t123);
t118 = cos(t123);
t83 = t116 * rSges(5,1) + t118 * rSges(5,2);
t186 = m(5) * t83;
t78 = t111 * rSges(6,1) + t112 * rSges(6,2);
t185 = m(6) * t78;
t184 = t117 / 0.2e1;
t183 = t119 / 0.2e1;
t182 = pkin(3) * t125;
t126 = sin(qJ(1));
t181 = t126 * pkin(1);
t128 = cos(qJ(1));
t180 = t128 * pkin(1);
t113 = t127 * pkin(3) + pkin(2);
t179 = pkin(2) - t113;
t122 = -pkin(9) + t129;
t104 = t117 * t122;
t105 = t117 * t129;
t90 = pkin(4) * t118 + t113;
t170 = t113 - t90;
t173 = rSges(6,2) * t111;
t152 = -t117 * rSges(6,3) + t119 * t173;
t174 = rSges(6,1) * t112;
t40 = t119 * (t119 * t174 - t152);
t178 = t119 * (-t170 * t119 - t104 + t105) + t40;
t172 = t119 * rSges(6,3) + t117 * t173;
t51 = -t117 * t174 + t172;
t177 = -(-t122 + t129) * t119 - t170 * t117 - t51;
t163 = t116 * t117;
t42 = pkin(4) * t163 + t117 * t78;
t176 = rSges(4,1) * t127;
t175 = rSges(5,1) * t118;
t171 = rSges(5,2) * t163 + t119 * rSges(5,3);
t169 = Icges(4,4) * t125;
t168 = Icges(4,4) * t127;
t167 = Icges(5,4) * t116;
t166 = Icges(5,4) * t118;
t165 = Icges(6,4) * t111;
t164 = Icges(6,4) * t112;
t162 = t117 * t125;
t161 = t119 * t129;
t160 = rSges(4,2) * t162 + t119 * rSges(4,3);
t159 = t114 + t115;
t134 = Icges(5,5) * t118 - Icges(5,6) * t116;
t52 = Icges(5,3) * t119 - t134 * t117;
t53 = Icges(5,3) * t117 + t134 * t119;
t158 = t117 * (t114 * t53 + t52 * t189) + t119 * (t115 * t52 + t53 * t189) + t188;
t136 = -Icges(6,2) * t111 + t164;
t139 = Icges(6,1) * t112 - t165;
t76 = Icges(6,2) * t112 + t165;
t77 = Icges(6,1) * t111 + t164;
t148 = t111 * t76 - t112 * t77;
t75 = Icges(6,5) * t111 + Icges(6,6) * t112;
t157 = (t111 * (Icges(6,5) * t117 + t139 * t119) + t112 * (Icges(6,6) * t117 + t136 * t119) + t117 * t75 - t148 * t119) * t184 + (t111 * (Icges(6,5) * t119 - t139 * t117) + t112 * (Icges(6,6) * t119 - t136 * t117) + t148 * t117 + t119 * t75) * t183;
t156 = -pkin(2) - t176;
t155 = -pkin(4) * t116 - t78;
t154 = -t90 - t174;
t153 = t119 * t116 * rSges(5,2) - t117 * rSges(5,3);
t85 = -t119 * rSges(3,1) + t117 * rSges(3,2);
t151 = -t113 - t175;
t84 = -t117 * rSges(3,1) - t119 * rSges(3,2);
t81 = Icges(5,2) * t118 + t167;
t82 = Icges(5,1) * t116 + t166;
t145 = t116 * t81 - t118 * t82;
t96 = Icges(4,2) * t127 + t169;
t97 = Icges(4,1) * t125 + t168;
t142 = t125 * t96 - t127 * t97;
t141 = Icges(4,1) * t127 - t169;
t140 = Icges(5,1) * t118 - t167;
t138 = -Icges(4,2) * t125 + t168;
t137 = -Icges(5,2) * t116 + t166;
t135 = Icges(4,5) * t127 - Icges(4,6) * t125;
t132 = t111 * t77 + t112 * t76 + t116 * t82 + t118 * t81 + t125 * t97 + t127 * t96 + Icges(3,3);
t80 = Icges(5,5) * t116 + Icges(5,6) * t118;
t131 = t157 + (t116 * (Icges(5,5) * t117 + t140 * t119) + t117 * t80 + t118 * (Icges(5,6) * t117 + t137 * t119) - t145 * t119) * t184 + (t116 * (Icges(5,5) * t119 - t140 * t117) + t145 * t117 + t118 * (Icges(5,6) * t119 - t137 * t117) + t119 * t80) * t183;
t110 = t119 * pkin(7);
t37 = t156 * t117 + t110 + t160;
t26 = t154 * t119 + t104 + t152;
t31 = t151 * t119 + t105 + t153;
t25 = t154 * t117 - t119 * t122 + t172;
t102 = t119 * t125 * rSges(4,2);
t38 = t102 + t156 * t119 + (-rSges(4,3) - pkin(7)) * t117;
t30 = t151 * t117 - t161 + t171;
t95 = Icges(4,5) * t125 + Icges(4,6) * t127;
t130 = t131 + (t117 * t95 - t142 * t119 + t125 * (Icges(4,5) * t117 + t141 * t119) + t127 * (Icges(4,6) * t117 + t138 * t119)) * t184 + (t142 * t117 + t119 * t95 + t125 * (Icges(4,5) * t119 - t141 * t117) + t127 * (Icges(4,6) * t119 - t138 * t117)) * t183;
t103 = pkin(3) * t162;
t100 = -t128 * rSges(2,1) + t126 * rSges(2,2);
t99 = -t126 * rSges(2,1) - t128 * rSges(2,2);
t73 = t85 - t180;
t72 = t84 - t181;
t62 = Icges(4,3) * t117 + t135 * t119;
t61 = Icges(4,3) * t119 - t135 * t117;
t60 = (-t83 - t182) * t119;
t59 = t117 * t83 + t103;
t58 = -t117 * t175 + t171;
t44 = t179 * t117 - t110 - t161;
t43 = t155 * t119;
t41 = t119 * (t119 * t175 - t153);
t39 = t119 * (-t117 * pkin(7) - t179 * t119 - t105);
t36 = (t155 - t182) * t119;
t35 = t103 + t42;
t34 = t38 - t180;
t33 = t37 - t181;
t28 = t31 - t180;
t27 = t30 - t181;
t22 = t119 * (t117 * rSges(4,3) + t119 * t176 - t102) - t117 * (-t117 * t176 + t160);
t21 = t26 - t180;
t20 = t25 - t181;
t17 = -t117 * t58 + t41;
t14 = -t117 * t51 + t40;
t7 = t39 + t41 + (-t44 - t58) * t117;
t6 = t177 * t117 + t178;
t3 = t39 + (-t44 + t177) * t117 + t178;
t1 = [Icges(2,3) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(2) * (t100 ^ 2 + t99 ^ 2) + m(3) * (t72 ^ 2 + t73 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + t132; m(6) * (t25 * t20 + t26 * t21) + m(5) * (t30 * t27 + t31 * t28) + m(3) * (t84 * t72 + t85 * t73) + m(4) * (t37 * t33 + t38 * t34) + t132; m(6) * (t25 ^ 2 + t26 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t84 ^ 2 + t85 ^ 2) + t132; m(6) * (t36 * t20 + t35 * t21) + m(5) * (t60 * t27 + t59 * t28) + (t117 * t34 - t119 * t33) * t187 + t130; m(6) * (t36 * t25 + t35 * t26) + m(5) * (t60 * t30 + t59 * t31) + (t117 * t38 - t119 * t37) * t187 + t130; m(6) * (t3 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t59 ^ 2 + t60 ^ 2 + t7 ^ 2) + t117 * (t114 * t62 + t189 * t61) + t119 * (t115 * t61 + t189 * t62) + m(4) * (t159 * t98 ^ 2 + t22 ^ 2) + t158; m(6) * (t43 * t20 + t42 * t21) + (t117 * t28 - t119 * t27) * t186 + t131; m(6) * (t43 * t25 + t42 * t26) + (t117 * t31 - t119 * t30) * t186 + t131; m(6) * (t6 * t3 + t42 * t35 + t43 * t36) + m(5) * (t17 * t7 + (t117 * t59 - t119 * t60) * t83) + t158; m(5) * (t159 * t83 ^ 2 + t17 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2 + t6 ^ 2) + t158; (t117 * t21 - t119 * t20) * t185 + t157; (t117 * t26 - t119 * t25) * t185 + t157; m(6) * (t14 * t3 + (t117 * t35 - t119 * t36) * t78) + t188; m(6) * (t14 * t6 + (t117 * t42 - t119 * t43) * t78) + t188; m(6) * (t159 * t78 ^ 2 + t14 ^ 2) + t188;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
