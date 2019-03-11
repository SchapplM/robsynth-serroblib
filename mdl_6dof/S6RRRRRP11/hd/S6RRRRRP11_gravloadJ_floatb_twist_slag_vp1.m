% Calculate Gravitation load on the joints for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:45
% EndTime: 2019-03-10 02:32:53
% DurationCPUTime: 2.73s
% Computational Cost: add. (1585->269), mult. (4261->393), div. (0->0), fcn. (5432->14), ass. (0->117)
t158 = cos(qJ(3));
t146 = cos(pkin(6));
t159 = cos(qJ(2));
t131 = t146 * t159;
t156 = sin(qJ(2));
t157 = sin(qJ(1));
t160 = cos(qJ(1));
t108 = -t131 * t160 + t157 * t156;
t143 = sin(pkin(7));
t144 = sin(pkin(6));
t122 = t144 * t143;
t145 = cos(pkin(7));
t177 = t108 * t145 + t160 * t122;
t130 = t146 * t156;
t70 = t130 * t160 + t157 * t159;
t90 = sin(qJ(3));
t40 = -t158 * t70 + t177 * t90;
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t123 = t145 * t144;
t98 = t108 * t143 - t160 * t123;
t20 = t40 * t92 - t89 * t98;
t37 = t177 * t158 + t70 * t90;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t1 = t20 * t88 + t37 * t91;
t180 = t20 * t91 - t37 * t88;
t19 = t40 * t89 + t92 * t98;
t103 = t131 * t157 + t156 * t160;
t93 = t103 * t143 + t157 * t123;
t175 = t103 * t145 - t157 * t122;
t174 = t159 * t123 + t143 * t146;
t71 = -t130 * t157 + t159 * t160;
t42 = t71 * t158 - t175 * t90;
t22 = t42 * t92 + t89 * t93;
t102 = -t122 * t159 + t145 * t146;
t125 = t144 * t156;
t56 = t158 * t125 + t174 * t90;
t36 = t102 * t89 + t56 * t92;
t173 = g(1) * t22 - g(2) * t20 + g(3) * t36;
t21 = t42 * t89 - t92 * t93;
t35 = -t102 * t92 + t56 * t89;
t172 = g(1) * t21 - g(2) * t19 + g(3) * t35;
t41 = t158 * t175 + t71 * t90;
t55 = t125 * t90 - t158 * t174;
t171 = t89 * (-g(1) * t41 - g(2) * t37 - g(3) * t55);
t170 = pkin(4) * t92;
t162 = pkin(11) + rSges(5,3);
t161 = pkin(12) + rSges(6,3);
t155 = t40 * t88;
t154 = t42 * t88;
t153 = t56 * t88;
t85 = pkin(5) * t91 + pkin(4);
t152 = t85 * t92;
t151 = t88 * t92;
t150 = t91 * t92;
t112 = t156 * t122;
t127 = t159 * t144;
t149 = pkin(2) * t127 + pkin(10) * t112;
t126 = t144 * t157;
t148 = t160 * pkin(1) + pkin(9) * t126;
t147 = qJ(6) + pkin(12) + rSges(7,3);
t115 = t156 * t123;
t64 = -t115 * t90 + t127 * t158;
t142 = t64 * pkin(3) + t149;
t141 = pkin(5) * t88 + pkin(11);
t31 = t37 * pkin(3);
t140 = -pkin(11) * t40 - t31;
t33 = t41 * pkin(3);
t139 = t42 * pkin(11) - t33;
t54 = t55 * pkin(3);
t138 = t56 * pkin(11) - t54;
t137 = t143 * pkin(10);
t136 = t89 * t143;
t135 = t90 * t145;
t134 = t92 * t143;
t128 = t160 * t144;
t133 = -pkin(1) * t157 + pkin(9) * t128;
t129 = t145 * t158;
t124 = -rSges(5,1) * t92 + rSges(5,2) * t89;
t5 = -t22 * t88 + t41 * t91;
t15 = -t36 * t88 + t55 * t91;
t46 = -t108 * t158 - t135 * t70;
t66 = t108 * pkin(2);
t118 = t46 * pkin(3) + t137 * t70 - t66;
t48 = -t103 * t158 - t135 * t71;
t68 = t103 * pkin(2);
t117 = t48 * pkin(3) + t137 * t71 - t68;
t109 = rSges(4,3) * t143 + t137;
t99 = -t70 * pkin(2) - pkin(10) * t98 + t133;
t96 = t40 * pkin(3) + t99;
t95 = t71 * pkin(2) + pkin(10) * t93 + t148;
t94 = t42 * pkin(3) + t95;
t63 = t115 * t158 + t127 * t90;
t50 = t112 * t89 + t64 * t92;
t49 = -t112 * t92 + t64 * t89;
t47 = -t103 * t90 + t129 * t71;
t45 = -t108 * t90 + t129 * t70;
t30 = t50 * t91 + t63 * t88;
t29 = -t50 * t88 + t63 * t91;
t28 = t136 * t71 + t48 * t92;
t27 = -t134 * t71 + t48 * t89;
t26 = t136 * t70 + t46 * t92;
t25 = -t134 * t70 + t46 * t89;
t24 = -t150 * t55 + t153;
t23 = t151 * t55 + t56 * t91;
t16 = -t36 * t91 - t55 * t88;
t14 = t28 * t91 + t47 * t88;
t13 = -t28 * t88 + t47 * t91;
t12 = t26 * t91 + t45 * t88;
t11 = -t26 * t88 + t45 * t91;
t10 = -t150 * t41 + t154;
t9 = t151 * t41 + t42 * t91;
t8 = -t150 * t37 - t155;
t7 = t151 * t37 - t40 * t91;
t6 = t22 * t91 + t41 * t88;
t2 = [-m(2) * (g(1) * (-rSges(2,1) * t157 - rSges(2,2) * t160) + g(2) * (rSges(2,1) * t160 - rSges(2,2) * t157)) - m(3) * (g(1) * (-t70 * rSges(3,1) + rSges(3,2) * t108 + rSges(3,3) * t128 + t133) + g(2) * (t71 * rSges(3,1) - rSges(3,2) * t103 + rSges(3,3) * t126 + t148)) - m(4) * (g(1) * (t40 * rSges(4,1) + rSges(4,2) * t37 - rSges(4,3) * t98 + t99) + g(2) * (t42 * rSges(4,1) - t41 * rSges(4,2) + rSges(4,3) * t93 + t95)) - m(5) * (g(1) * (t20 * rSges(5,1) - t19 * rSges(5,2) - t162 * t37 + t96) + g(2) * (t22 * rSges(5,1) - t21 * rSges(5,2) + t162 * t41 + t94)) - m(6) * (g(1) * (rSges(6,1) * t180 - rSges(6,2) * t1 + t20 * pkin(4) - pkin(11) * t37 + t161 * t19 + t96) + g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t22 * pkin(4) + t41 * pkin(11) + t161 * t21 + t94)) - m(7) * (g(1) * (rSges(7,1) * t180 - rSges(7,2) * t1 - t141 * t37 + t147 * t19 + t20 * t85 + t96) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t141 * t41 + t147 * t21 + t22 * t85 + t94)) -m(3) * (g(1) * (-rSges(3,1) * t103 - t71 * rSges(3,2)) + g(2) * (-rSges(3,1) * t108 - t70 * rSges(3,2)) + g(3) * (rSges(3,1) * t127 - rSges(3,2) * t125)) - m(4) * (g(1) * (t48 * rSges(4,1) - t47 * rSges(4,2) + t109 * t71 - t68) + g(2) * (t46 * rSges(4,1) - t45 * rSges(4,2) + t109 * t70 - t66) + g(3) * (t64 * rSges(4,1) - t63 * rSges(4,2) + rSges(4,3) * t112 + t149)) - m(5) * (g(1) * (t28 * rSges(5,1) - t27 * rSges(5,2) + t162 * t47 + t117) + g(2) * (t26 * rSges(5,1) - t25 * rSges(5,2) + t162 * t45 + t118) + g(3) * (rSges(5,1) * t50 - rSges(5,2) * t49 + t162 * t63 + t142)) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t28 * pkin(4) + t47 * pkin(11) + t161 * t27 + t117) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t26 * pkin(4) + t45 * pkin(11) + t161 * t25 + t118) + g(3) * (rSges(6,1) * t30 + rSges(6,2) * t29 + pkin(4) * t50 + pkin(11) * t63 + t161 * t49 + t142)) - m(7) * (g(1) * (t14 * rSges(7,1) + t13 * rSges(7,2) + t141 * t47 + t147 * t27 + t28 * t85 + t117) + g(2) * (t12 * rSges(7,1) + t11 * rSges(7,2) + t141 * t45 + t147 * t25 + t26 * t85 + t118) + g(3) * (rSges(7,1) * t30 + rSges(7,2) * t29 + t141 * t63 + t147 * t49 + t50 * t85 + t142)) -m(4) * (g(1) * (-rSges(4,1) * t41 - rSges(4,2) * t42) + g(2) * (-rSges(4,1) * t37 + rSges(4,2) * t40) + g(3) * (-rSges(4,1) * t55 - rSges(4,2) * t56)) - m(5) * (g(1) * (t124 * t41 + t162 * t42 - t33) + g(2) * (t124 * t37 - t162 * t40 - t31) + g(3) * (t124 * t55 + t162 * t56 - t54)) + (-g(1) * (t10 * rSges(7,1) + t9 * rSges(7,2) + pkin(5) * t154 - t152 * t41 + t139) - g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) - pkin(5) * t155 - t152 * t37 + t140) - g(3) * (t24 * rSges(7,1) + t23 * rSges(7,2) + pkin(5) * t153 - t152 * t55 + t138) - t147 * t171) * m(7) + (-g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) - t170 * t41 + t139) - g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) - t170 * t37 + t140) - g(3) * (t24 * rSges(6,1) + t23 * rSges(6,2) - t170 * t55 + t138) - t161 * t171) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t21 - rSges(5,2) * t22) + g(2) * (rSges(5,1) * t19 + rSges(5,2) * t20) + g(3) * (-rSges(5,1) * t35 - rSges(5,2) * t36)) - m(6) * (t173 * t161 + t172 * (-rSges(6,1) * t91 + rSges(6,2) * t88 - pkin(4))) - m(7) * (t173 * t147 + t172 * (-rSges(7,1) * t91 + rSges(7,2) * t88 - t85)) -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t180) + g(3) * (rSges(6,1) * t15 + rSges(6,2) * t16)) + (-g(1) * (rSges(7,1) * t5 - rSges(7,2) * t6) - g(2) * (rSges(7,1) * t1 + rSges(7,2) * t180) - g(3) * (rSges(7,1) * t15 + rSges(7,2) * t16) - (g(1) * t5 + g(2) * t1 + g(3) * t15) * pkin(5)) * m(7), -m(7) * t172];
taug  = t2(:);
