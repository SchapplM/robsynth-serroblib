% Calculate Gravitation load on the joints for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:32
% EndTime: 2019-03-10 02:00:36
% DurationCPUTime: 1.60s
% Computational Cost: add. (876->213), mult. (1743->303), div. (0->0), fcn. (2087->12), ass. (0->88)
t113 = sin(qJ(1));
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t114 = cos(qJ(1));
t96 = cos(pkin(6));
t87 = t96 * t114;
t42 = t113 * t71 + t68 * t87;
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t65 = sin(pkin(6));
t92 = t65 * t114;
t28 = t42 * t70 - t67 * t92;
t41 = t113 * t68 - t71 * t87;
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t136 = t28 * t66 - t41 * t69;
t135 = -t28 * t69 - t41 * t66;
t64 = qJ(4) + qJ(5);
t59 = sin(t64);
t60 = cos(t64);
t11 = t28 * t59 - t41 * t60;
t134 = -t28 * t60 - t41 * t59;
t86 = t96 * t113;
t43 = t114 * t68 + t71 * t86;
t133 = g(1) * t43 + g(2) * t41;
t72 = -pkin(11) - pkin(10);
t100 = rSges(7,3) + qJ(6) - t72;
t61 = t69 * pkin(4);
t47 = pkin(5) * t60 + t61;
t45 = pkin(3) + t47;
t132 = t100 * t67 + t45 * t70;
t101 = rSges(6,3) - t72;
t58 = t61 + pkin(3);
t131 = t101 * t67 + t58 * t70;
t44 = t114 * t71 - t68 * t86;
t91 = t65 * t113;
t32 = t44 * t70 + t67 * t91;
t104 = t65 * t68;
t40 = t70 * t104 + t96 * t67;
t130 = g(1) * t32 + g(2) * t28 + g(3) * t40;
t27 = t42 * t67 + t70 * t92;
t31 = t44 * t67 - t70 * t91;
t39 = t67 * t104 - t96 * t70;
t129 = g(1) * t31 + g(2) * t27 + g(3) * t39;
t115 = pkin(10) + rSges(5,3);
t81 = rSges(5,1) * t69 - rSges(5,2) * t66 + pkin(3);
t128 = t115 * t67 + t81 * t70;
t127 = -rSges(7,1) * t11 + rSges(7,2) * t134;
t13 = -t32 * t59 + t43 * t60;
t14 = t32 * t60 + t43 * t59;
t126 = t13 * rSges(7,1) - t14 * rSges(7,2);
t125 = pkin(4) * t66;
t118 = g(3) * t65;
t117 = rSges(4,3) + pkin(9);
t46 = pkin(5) * t59 + t125;
t116 = pkin(9) + t46;
t106 = t59 * t70;
t105 = t60 * t70;
t103 = t65 * t71;
t102 = t70 * t71;
t25 = -t60 * t103 - t40 * t59;
t26 = t59 * t103 - t40 * t60;
t99 = t25 * rSges(7,1) + t26 * rSges(7,2);
t98 = t114 * pkin(1) + pkin(8) * t91;
t97 = pkin(2) * t103 + pkin(9) * t104;
t95 = t44 * pkin(2) + t98;
t94 = pkin(9) + t125;
t93 = g(3) * t97;
t35 = t41 * pkin(2);
t90 = t42 * pkin(9) - t35;
t37 = t43 * pkin(2);
t89 = t44 * pkin(9) - t37;
t88 = -t113 * pkin(1) + pkin(8) * t92;
t85 = rSges(4,1) * t70 - rSges(4,2) * t67;
t84 = rSges(5,1) * t66 + rSges(5,2) * t69;
t15 = -t32 * t66 + t43 * t69;
t82 = -t42 * pkin(2) + t88;
t80 = pkin(9) + t84;
t79 = -t69 * t103 - t40 * t66;
t73 = m(6) * (g(1) * (t13 * rSges(6,1) - t14 * rSges(6,2)) + g(2) * (-rSges(6,1) * t11 + rSges(6,2) * t134) + g(3) * (t25 * rSges(6,1) + t26 * rSges(6,2)));
t34 = (t60 * t102 + t59 * t68) * t65;
t33 = (-t59 * t102 + t60 * t68) * t65;
t20 = -t43 * t105 + t44 * t59;
t19 = t43 * t106 + t44 * t60;
t18 = -t41 * t105 + t42 * t59;
t17 = t41 * t106 + t42 * t60;
t16 = t32 * t69 + t43 * t66;
t1 = [-m(2) * (g(1) * (-t113 * rSges(2,1) - t114 * rSges(2,2)) + g(2) * (t114 * rSges(2,1) - t113 * rSges(2,2))) - m(3) * (g(1) * (-t42 * rSges(3,1) + t41 * rSges(3,2) + rSges(3,3) * t92 + t88) + g(2) * (t44 * rSges(3,1) - t43 * rSges(3,2) + rSges(3,3) * t91 + t98)) - m(4) * (g(1) * (-rSges(4,1) * t28 + rSges(4,2) * t27 - t117 * t41 + t82) + g(2) * (rSges(4,1) * t32 - rSges(4,2) * t31 + t117 * t43 + t95)) - m(5) * (g(1) * (t135 * rSges(5,1) + rSges(5,2) * t136 - t28 * pkin(3) - t41 * pkin(9) - t115 * t27 + t82) + g(2) * (rSges(5,1) * t16 + rSges(5,2) * t15 + pkin(3) * t32 + t43 * pkin(9) + t115 * t31 + t95)) - m(6) * (g(1) * (rSges(6,1) * t134 + t11 * rSges(6,2) - t101 * t27 - t28 * t58 - t94 * t41 + t82) + g(2) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t101 * t31 + t32 * t58 + t94 * t43 + t95)) - m(7) * (g(1) * (rSges(7,1) * t134 + t11 * rSges(7,2) - t100 * t27 - t116 * t41 - t28 * t45 + t82) + g(2) * (rSges(7,1) * t14 + rSges(7,2) * t13 + t100 * t31 + t116 * t43 + t32 * t45 + t95)) -m(3) * (g(1) * (-rSges(3,1) * t43 - rSges(3,2) * t44) + g(2) * (-rSges(3,1) * t41 - rSges(3,2) * t42) + (rSges(3,1) * t71 - rSges(3,2) * t68) * t118) - m(4) * (g(1) * (t117 * t44 - t85 * t43 - t37) + g(2) * (t117 * t42 - t85 * t41 - t35) + t93 + (rSges(4,3) * t68 + t85 * t71) * t118) - m(5) * (t93 + (t128 * t71 + t84 * t68) * t118 - t133 * t128 + (t80 * t42 - t35) * g(2) + (t80 * t44 - t37) * g(1)) - m(6) * (g(1) * (t20 * rSges(6,1) + t19 * rSges(6,2) + t44 * t125 + t89) + g(2) * (t18 * rSges(6,1) + t17 * rSges(6,2) + t42 * t125 + t90) + g(3) * (t34 * rSges(6,1) + t33 * rSges(6,2) + t97) + (t68 * t125 + t131 * t71) * t118 - t133 * t131) - m(7) * (g(1) * (rSges(7,1) * t20 + rSges(7,2) * t19 + t44 * t46 + t89) + g(2) * (rSges(7,1) * t18 + rSges(7,2) * t17 + t42 * t46 + t90) + g(3) * (rSges(7,1) * t34 + rSges(7,2) * t33 + t97) + (t132 * t71 + t46 * t68) * t118 - t133 * t132) -m(4) * (g(1) * (-rSges(4,1) * t31 - rSges(4,2) * t32) + g(2) * (-rSges(4,1) * t27 - rSges(4,2) * t28) + g(3) * (-rSges(4,1) * t39 - rSges(4,2) * t40)) - m(5) * (t130 * t115 - t129 * t81) - m(6) * (t129 * (-rSges(6,1) * t60 + rSges(6,2) * t59 - t58) + t130 * t101) - m(7) * (t129 * (-rSges(7,1) * t60 + rSges(7,2) * t59 - t45) + t130 * t100) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (-rSges(5,1) * t136 + rSges(5,2) * t135) + g(3) * (t79 * rSges(5,1) + (t66 * t103 - t40 * t69) * rSges(5,2))) - t73 - m(7) * (g(1) * (-t32 * t46 + t43 * t47 + t126) + g(2) * (-t28 * t46 + t41 * t47 + t127) + g(3) * (-t47 * t103 - t40 * t46 + t99)) - m(6) * (g(1) * t15 - g(2) * t136 + g(3) * t79) * pkin(4), -t73 + (-g(1) * t126 - g(2) * t127 - g(3) * t99 - (g(1) * t13 - g(2) * t11 + g(3) * t25) * pkin(5)) * m(7), -m(7) * t129];
taug  = t1(:);
