% Calculate Gravitation load on the joints for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:16
% EndTime: 2019-03-09 23:29:21
% DurationCPUTime: 2.38s
% Computational Cost: add. (1446->274), mult. (3447->412), div. (0->0), fcn. (4347->16), ass. (0->109)
t136 = cos(pkin(7));
t161 = cos(qJ(3));
t115 = t136 * t161;
t100 = cos(qJ(1));
t93 = sin(pkin(6));
t138 = t100 * t93;
t92 = sin(pkin(7));
t133 = t92 * t138;
t137 = cos(pkin(6));
t162 = cos(qJ(2));
t117 = t137 * t162;
t159 = sin(qJ(2));
t160 = sin(qJ(1));
t70 = -t100 * t117 + t159 * t160;
t116 = t137 * t159;
t71 = t100 * t116 + t160 * t162;
t97 = sin(qJ(3));
t26 = t115 * t70 + t133 * t161 + t71 * t97;
t127 = t93 * t136;
t180 = t100 * t127 - t70 * t92;
t126 = t97 * t136;
t27 = -t126 * t70 - t97 * t133 + t161 * t71;
t91 = qJ(4) + pkin(13);
t88 = sin(t91);
t89 = cos(t91);
t5 = -t180 * t88 + t27 * t89;
t95 = sin(qJ(6));
t98 = cos(qJ(6));
t187 = -t26 * t98 + t5 * t95;
t186 = -t26 * t95 - t5 * t98;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t185 = t180 * t99 + t27 * t96;
t155 = t180 * t96;
t184 = -t27 * t99 + t155;
t107 = t100 * t159 + t117 * t160;
t54 = t107 * t92 + t160 * t127;
t128 = t92 * t137;
t51 = t97 * t128 + (t126 * t162 + t159 * t161) * t93;
t131 = t93 * t162;
t69 = -t131 * t92 + t136 * t137;
t179 = -t51 * t96 + t69 * t99;
t130 = t93 * t160;
t176 = t107 * t136 - t92 * t130;
t72 = t100 * t162 - t116 * t160;
t31 = t72 * t161 - t176 * t97;
t10 = -t31 * t96 + t54 * t99;
t30 = t176 * t161 + t72 * t97;
t129 = t93 * t159;
t50 = -t115 * t131 - t128 * t161 + t129 * t97;
t173 = g(1) * t30 + g(2) * t26 + g(3) * t50;
t177 = -t180 * t89 - t27 * t88;
t163 = pkin(12) + rSges(7,3);
t175 = g(1) * t72 + g(2) * t71;
t172 = pkin(5) * t89;
t166 = t92 * pkin(10);
t164 = pkin(11) + rSges(5,3);
t153 = t54 * t96;
t149 = t88 * t92;
t148 = t89 * t92;
t147 = t89 * t95;
t146 = t89 * t98;
t145 = t92 * t96;
t144 = t92 * t99;
t87 = pkin(4) * t99 + pkin(3);
t94 = -qJ(5) - pkin(11);
t143 = -t26 * t87 - t27 * t94;
t142 = -t30 * t87 - t31 * t94;
t141 = -t50 * t87 - t51 * t94;
t125 = t92 * t129;
t140 = pkin(2) * t131 + pkin(10) * t125;
t139 = t100 * pkin(1) + pkin(9) * t130;
t38 = t115 * t71 - t70 * t97;
t39 = -t126 * t71 - t161 * t70;
t65 = t70 * pkin(2);
t135 = -t38 * t94 + t39 * t87 - t65;
t40 = -t107 * t97 + t115 * t72;
t41 = -t107 * t161 - t126 * t72;
t67 = t107 * pkin(2);
t134 = -t40 * t94 + t41 * t87 - t67;
t123 = t185 * pkin(4);
t122 = t10 * pkin(4);
t121 = t179 * pkin(4);
t119 = -pkin(1) * t160 + pkin(9) * t138;
t112 = t96 * t125;
t114 = t136 * t159;
t62 = (t114 * t161 + t162 * t97) * t93;
t63 = (-t114 * t97 + t161 * t162) * t93;
t118 = pkin(4) * t112 - t62 * t94 + t63 * t87 + t140;
t113 = -rSges(6,1) * t89 + rSges(6,2) * t88;
t108 = -t71 * pkin(2) + t180 * pkin(10) + t119;
t106 = pkin(4) * t155 + t26 * t94 - t27 * t87 + t108;
t104 = t175 * (pkin(4) * t96 + pkin(10)) * t92;
t102 = t72 * pkin(2) + t54 * pkin(10) + t139;
t101 = pkin(4) * t153 - t30 * t94 + t31 * t87 + t102;
t37 = t125 * t88 + t63 * t89;
t36 = -t125 * t89 + t63 * t88;
t21 = t51 * t89 + t69 * t88;
t20 = -t51 * t88 + t69 * t89;
t15 = t149 * t72 + t41 * t89;
t14 = -t148 * t72 + t41 * t88;
t13 = t149 * t71 + t39 * t89;
t12 = -t148 * t71 + t39 * t88;
t11 = t31 * t99 + t153;
t9 = t31 * t89 + t54 * t88;
t8 = t31 * t88 - t54 * t89;
t2 = t30 * t95 + t9 * t98;
t1 = t30 * t98 - t9 * t95;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t160 - t100 * rSges(2,2)) + g(2) * (t100 * rSges(2,1) - rSges(2,2) * t160)) - m(3) * (g(1) * (-t71 * rSges(3,1) + t70 * rSges(3,2) + rSges(3,3) * t138 + t119) + g(2) * (t72 * rSges(3,1) - rSges(3,2) * t107 + rSges(3,3) * t130 + t139)) - m(4) * (g(1) * (-rSges(4,1) * t27 + rSges(4,2) * t26 + rSges(4,3) * t180 + t108) + g(2) * (t31 * rSges(4,1) - t30 * rSges(4,2) + t54 * rSges(4,3) + t102)) - m(5) * (g(1) * (t184 * rSges(5,1) + t185 * rSges(5,2) - t27 * pkin(3) - t164 * t26 + t108) + g(2) * (t11 * rSges(5,1) + t10 * rSges(5,2) + t31 * pkin(3) + t164 * t30 + t102)) - m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t177 - rSges(6,3) * t26 + t106) + g(2) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t30 * rSges(6,3) + t101)) - m(7) * (g(1) * (t186 * rSges(7,1) + t187 * rSges(7,2) - t5 * pkin(5) + t163 * t177 + t106) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t9 * pkin(5) + t163 * t8 + t101)) -m(3) * (g(1) * (-rSges(3,1) * t107 - t72 * rSges(3,2)) + g(2) * (-rSges(3,1) * t70 - rSges(3,2) * t71) + g(3) * (rSges(3,1) * t162 - rSges(3,2) * t159) * t93) - m(4) * (g(1) * (rSges(4,1) * t41 - rSges(4,2) * t40 - t67) + g(2) * (rSges(4,1) * t39 - rSges(4,2) * t38 - t65) + g(3) * (t63 * rSges(4,1) - t62 * rSges(4,2) + t140) + (rSges(4,3) * g(3) * t129 + t175 * (rSges(4,3) + pkin(10))) * t92) - m(5) * (g(1) * (t41 * pkin(3) - t67 + t72 * t166 + (t145 * t72 + t41 * t99) * rSges(5,1) + (t144 * t72 - t41 * t96) * rSges(5,2) + t164 * t40) + g(2) * (t39 * pkin(3) - t65 + t71 * t166 + (t145 * t71 + t39 * t99) * rSges(5,1) + (t144 * t71 - t39 * t96) * rSges(5,2) + t164 * t38) + g(3) * (t63 * pkin(3) + (t63 * t99 + t112) * rSges(5,1) + (t125 * t99 - t63 * t96) * rSges(5,2) + t164 * t62 + t140)) - m(6) * (g(1) * (rSges(6,1) * t15 - rSges(6,2) * t14 + rSges(6,3) * t40 + t134) + g(2) * (rSges(6,1) * t13 - rSges(6,2) * t12 + rSges(6,3) * t38 + t135) + g(3) * (rSges(6,1) * t37 - rSges(6,2) * t36 + rSges(6,3) * t62 + t118) + t104) - m(7) * (g(1) * (t15 * pkin(5) + (t15 * t98 + t40 * t95) * rSges(7,1) + (-t15 * t95 + t40 * t98) * rSges(7,2) + t134 + t163 * t14) + g(2) * (t13 * pkin(5) + (t13 * t98 + t38 * t95) * rSges(7,1) + (-t13 * t95 + t38 * t98) * rSges(7,2) + t135 + t163 * t12) + g(3) * (t37 * pkin(5) + (t37 * t98 + t62 * t95) * rSges(7,1) + (-t37 * t95 + t62 * t98) * rSges(7,2) + t163 * t36 + t118) + t104) -m(4) * (g(1) * (-rSges(4,1) * t30 - rSges(4,2) * t31) + g(2) * (-rSges(4,1) * t26 - rSges(4,2) * t27) + g(3) * (-rSges(4,1) * t50 - rSges(4,2) * t51)) - m(5) * ((g(1) * t31 + g(2) * t27 + g(3) * t51) * t164 + t173 * (-t99 * rSges(5,1) + t96 * rSges(5,2) - pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t31 + t113 * t30 + t142) + g(2) * (rSges(6,3) * t27 + t113 * t26 + t143) + g(3) * (rSges(6,3) * t51 + t113 * t50 + t141)) + (-g(1) * (-t30 * t172 + (-t146 * t30 + t31 * t95) * rSges(7,1) + (t147 * t30 + t31 * t98) * rSges(7,2) + t142) - g(2) * (-t26 * t172 + (-t146 * t26 + t27 * t95) * rSges(7,1) + (t147 * t26 + t27 * t98) * rSges(7,2) + t143) - g(3) * (-t50 * t172 + (-t146 * t50 + t51 * t95) * rSges(7,1) + (t147 * t50 + t51 * t98) * rSges(7,2) + t141) + t173 * t88 * t163) * m(7), -m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (-rSges(5,1) * t185 + t184 * rSges(5,2)) + g(3) * (t179 * rSges(5,1) + (-t51 * t99 - t69 * t96) * rSges(5,2))) - m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9 + t122) + g(2) * (rSges(6,1) * t177 - rSges(6,2) * t5 - t123) + g(3) * (rSges(6,1) * t20 - rSges(6,2) * t21 + t121)) + (-g(1) * (t163 * t9 + t122) - g(2) * (t163 * t5 - t123) - g(3) * (t163 * t21 + t121) - (-g(1) * t8 + g(2) * t177 + g(3) * t20) * (rSges(7,1) * t98 - rSges(7,2) * t95 + pkin(5))) * m(7) (-m(6) - m(7)) * t173, -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t187 * rSges(7,1) + t186 * rSges(7,2)) + g(3) * ((-t21 * t95 + t50 * t98) * rSges(7,1) + (-t21 * t98 - t50 * t95) * rSges(7,2)))];
taug  = t3(:);
