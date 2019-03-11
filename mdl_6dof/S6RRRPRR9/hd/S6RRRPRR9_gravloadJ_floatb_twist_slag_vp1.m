% Calculate Gravitation load on the joints for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:37
% EndTime: 2019-03-09 18:56:39
% DurationCPUTime: 2.46s
% Computational Cost: add. (1582->281), mult. (4165->429), div. (0->0), fcn. (5341->16), ass. (0->120)
t100 = cos(qJ(3));
t138 = cos(pkin(13));
t89 = sin(pkin(13));
t95 = sin(qJ(3));
t108 = t100 * t89 + t138 * t95;
t102 = cos(qJ(1));
t91 = sin(pkin(6));
t140 = t102 * t91;
t74 = -t100 * t138 + t95 * t89;
t90 = sin(pkin(7));
t63 = t74 * t90;
t92 = cos(pkin(7));
t65 = t74 * t92;
t101 = cos(qJ(2));
t139 = cos(pkin(6));
t123 = t102 * t139;
t96 = sin(qJ(2));
t97 = sin(qJ(1));
t70 = -t101 * t123 + t96 * t97;
t71 = t97 * t101 + t123 * t96;
t22 = t108 * t71 - t140 * t63 - t70 * t65;
t64 = t108 * t90;
t66 = t108 * t92;
t104 = t140 * t64 + t70 * t66 + t71 * t74;
t51 = t140 * t92 - t70 * t90;
t94 = sin(qJ(5));
t99 = cos(qJ(5));
t6 = t104 * t99 + t51 * t94;
t93 = sin(qJ(6));
t98 = cos(qJ(6));
t176 = t22 * t98 + t6 * t93;
t175 = -t22 * t93 + t6 * t98;
t131 = t90 * t140;
t144 = t100 * t92;
t136 = pkin(3) * t144;
t159 = t71 * t95;
t103 = (-t100 * t131 - t159) * pkin(3) - t70 * t136;
t172 = -t22 * pkin(4) + t103;
t137 = t100 * t101;
t150 = t95 * t96;
t171 = t92 * t137 - t150;
t170 = t104 * t94 - t51 * t99;
t161 = rSges(7,3) + pkin(12);
t169 = g(1) * t161;
t31 = (t101 * t66 - t74 * t96) * t91 + t139 * t64;
t168 = t100 * (t70 * t92 + t131) + t159;
t167 = pkin(3) * t90;
t165 = g(3) * t91;
t164 = t99 * pkin(5);
t162 = -rSges(6,3) - pkin(11);
t160 = pkin(3) * t100;
t126 = t97 * t139;
t73 = t101 * t102 - t126 * t96;
t158 = t73 * t95;
t157 = t90 * rSges(5,3);
t156 = t90 * t94;
t155 = t90 * t99;
t154 = t91 * t96;
t153 = t91 * t97;
t152 = t92 * t95;
t151 = t93 * t99;
t149 = t98 * t99;
t148 = pkin(10) + qJ(4);
t68 = pkin(3) * t152 - t148 * t90;
t87 = pkin(2) + t160;
t147 = -t71 * t68 - t70 * t87;
t72 = -t101 * t126 - t102 * t96;
t146 = -t73 * t68 + t72 * t87;
t145 = t102 * pkin(1) + pkin(9) * t153;
t143 = t100 * t96;
t142 = t101 * t91;
t141 = t101 * t95;
t135 = t90 * t154;
t134 = t90 * t153;
t35 = -t66 * t71 + t70 * t74;
t133 = t35 * pkin(4) + t147;
t37 = -t66 * t73 - t72 * t74;
t132 = t37 * pkin(4) + t146;
t130 = g(2) * t161;
t129 = g(3) * t161;
t127 = -t97 * pkin(1) + pkin(9) * t140;
t124 = t100 * t139;
t67 = t148 * t92 + t167 * t95;
t121 = t67 * t153 + t72 * t68 + t73 * t87 + t145;
t120 = rSges(6,1) * t99 - rSges(6,2) * t94;
t118 = -pkin(3) * t158 + t134 * t160 + t72 * t136;
t44 = (-t101 * t74 - t66 * t96) * t91;
t78 = t87 * t142;
t117 = t44 * pkin(4) - t68 * t154 + t78;
t26 = t153 * t64 + t66 * t72 - t73 * t74;
t116 = t26 * pkin(4) + t121;
t114 = t153 * t92 - t72 * t90;
t113 = t72 * t92 + t134;
t111 = pkin(3) * t171 * t91 + t124 * t167;
t25 = -t108 * t73 - t153 * t63 - t65 * t72;
t110 = t25 * pkin(4) + t118;
t109 = t67 * t140 + t70 * t68 - t71 * t87 + t127;
t107 = -t100 * t71 + t95 * t131 + t152 * t70;
t30 = -t139 * t63 + (-t101 * t65 - t108 * t96) * t91;
t106 = t30 * pkin(4) + t111;
t105 = pkin(4) * t104 + t109;
t69 = t139 * t92 - t142 * t90;
t43 = -t108 * t142 + t154 * t65;
t41 = t100 * t73 + t113 * t95;
t40 = t100 * t113 - t158;
t39 = t135 * t94 + t44 * t99;
t38 = -t135 * t99 + t44 * t94;
t36 = -t108 * t72 + t65 * t73;
t34 = t108 * t70 + t65 * t71;
t18 = t156 * t73 + t37 * t99;
t17 = -t155 * t73 + t37 * t94;
t16 = t156 * t71 + t35 * t99;
t15 = -t155 * t71 + t35 * t94;
t14 = t31 * t99 + t69 * t94;
t13 = -t31 * t94 + t69 * t99;
t8 = t114 * t94 + t26 * t99;
t7 = -t114 * t99 + t26 * t94;
t2 = -t25 * t93 + t8 * t98;
t1 = -t25 * t98 - t8 * t93;
t3 = [-m(2) * (g(1) * (-t97 * rSges(2,1) - rSges(2,2) * t102) + g(2) * (rSges(2,1) * t102 - t97 * rSges(2,2))) - m(3) * (g(1) * (-t71 * rSges(3,1) + t70 * rSges(3,2) + rSges(3,3) * t140 + t127) + g(2) * (rSges(3,1) * t73 + rSges(3,2) * t72 + rSges(3,3) * t153 + t145)) - m(4) * (g(1) * (t107 * rSges(4,1) + t168 * rSges(4,2) + t51 * rSges(4,3) - t71 * pkin(2) + t127) + g(2) * (t41 * rSges(4,1) + t40 * rSges(4,2) + rSges(4,3) * t114 + t73 * pkin(2) + t145) + (g(1) * t51 + g(2) * t114) * pkin(10)) - m(5) * (g(1) * (rSges(5,1) * t104 + rSges(5,2) * t22 + rSges(5,3) * t51 + t109) + g(2) * (t26 * rSges(5,1) + t25 * rSges(5,2) + rSges(5,3) * t114 + t121)) - m(6) * (g(1) * (rSges(6,1) * t6 - rSges(6,2) * t170 + t162 * t22 + t105) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t162 * t25 + t116)) - m(7) * (g(1) * (t175 * rSges(7,1) - t176 * rSges(7,2) + t6 * pkin(5) - t22 * pkin(11) + t161 * t170 + t105) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t8 - pkin(11) * t25 + t161 * t7 + t116)) -m(3) * (g(1) * (rSges(3,1) * t72 - rSges(3,2) * t73) + g(2) * (-rSges(3,1) * t70 - rSges(3,2) * t71) + (rSges(3,1) * t101 - rSges(3,2) * t96) * t165) - m(4) * (g(1) * (t72 * pkin(2) + (t100 * t72 - t152 * t73) * rSges(4,1) + (-t144 * t73 - t72 * t95) * rSges(4,2)) + g(2) * (-t70 * pkin(2) + (-t100 * t70 - t152 * t71) * rSges(4,1) + (-t144 * t71 + t70 * t95) * rSges(4,2)) + (t101 * pkin(2) + (-t150 * t92 + t137) * rSges(4,1) + (-t143 * t92 - t141) * rSges(4,2)) * t165 + (g(1) * t73 + g(2) * t71 + t165 * t96) * t90 * (rSges(4,3) + pkin(10))) - m(5) * (g(1) * (rSges(5,1) * t37 + rSges(5,2) * t36 + t157 * t73 + t146) + g(2) * (rSges(5,1) * t35 + rSges(5,2) * t34 + t157 * t71 + t147) + g(3) * (rSges(5,1) * t44 + rSges(5,2) * t43 + t78 + (-t68 + t157) * t154)) - m(6) * (g(1) * (rSges(6,1) * t18 - rSges(6,2) * t17 + t162 * t36 + t132) + g(2) * (rSges(6,1) * t16 - rSges(6,2) * t15 + t162 * t34 + t133) + g(3) * (rSges(6,1) * t39 - rSges(6,2) * t38 + t162 * t43 + t117)) - m(7) * (g(1) * (t18 * pkin(5) - t36 * pkin(11) + (t18 * t98 - t36 * t93) * rSges(7,1) + (-t18 * t93 - t36 * t98) * rSges(7,2) + t161 * t17 + t132) + g(2) * (t16 * pkin(5) - t34 * pkin(11) + (t16 * t98 - t34 * t93) * rSges(7,1) + (-t16 * t93 - t34 * t98) * rSges(7,2) + t161 * t15 + t133) + g(3) * (t39 * pkin(5) - t43 * pkin(11) + (t39 * t98 - t43 * t93) * rSges(7,1) + (-t39 * t93 - t43 * t98) * rSges(7,2) + t161 * t38 + t117)) -m(4) * (g(1) * (rSges(4,1) * t40 - rSges(4,2) * t41) + g(2) * (-t168 * rSges(4,1) + t107 * rSges(4,2)) + g(3) * ((-rSges(4,2) * t139 * t95 + rSges(4,1) * t124) * t90 + (t171 * rSges(4,1) + (-t141 * t92 - t143) * rSges(4,2)) * t91)) - m(5) * (g(1) * (rSges(5,1) * t25 - rSges(5,2) * t26 + t118) + g(2) * (-rSges(5,1) * t22 + rSges(5,2) * t104 + t103) + g(3) * (rSges(5,1) * t30 - rSges(5,2) * t31 + t111)) - m(6) * (g(1) * (t120 * t25 - t162 * t26 + t110) + g(2) * (t104 * t162 - t120 * t22 + t172) + g(3) * (t120 * t30 - t162 * t31 + t106)) - m(7) * (g(1) * (t25 * t164 + t26 * pkin(11) + (t149 * t25 + t26 * t93) * rSges(7,1) + (-t151 * t25 + t26 * t98) * rSges(7,2) + t110) + g(2) * (-t22 * t164 - t104 * pkin(11) + (-t104 * t93 - t149 * t22) * rSges(7,1) + (-t104 * t98 + t151 * t22) * rSges(7,2) + t172) + g(3) * (t30 * t164 + t31 * pkin(11) + (t149 * t30 + t31 * t93) * rSges(7,1) + (-t151 * t30 + t31 * t98) * rSges(7,2) + t106) + (t129 * t30 - t130 * t22 + t25 * t169) * t94) (-m(5) - m(6) - m(7)) * (g(1) * t114 - g(2) * t51 + g(3) * t69) -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t170 + rSges(6,2) * t6) + g(3) * (rSges(6,1) * t13 - rSges(6,2) * t14)) - m(7) * (t8 * t169 + t14 * t129 - t6 * t130 + (-g(1) * t7 + g(2) * t170 + g(3) * t13) * (t98 * rSges(7,1) - t93 * rSges(7,2) + pkin(5))) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t176 * rSges(7,1) + t175 * rSges(7,2)) + g(3) * ((-t14 * t93 - t30 * t98) * rSges(7,1) + (-t14 * t98 + t30 * t93) * rSges(7,2)))];
taug  = t3(:);
