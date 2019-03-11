% Calculate Gravitation load on the joints for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:30
% EndTime: 2019-03-09 06:34:34
% DurationCPUTime: 1.75s
% Computational Cost: add. (1219->185), mult. (3260->269), div. (0->0), fcn. (4170->14), ass. (0->84)
t118 = cos(qJ(4));
t119 = cos(qJ(3));
t110 = cos(pkin(7));
t120 = cos(qJ(1));
t106 = sin(pkin(12));
t117 = sin(qJ(1));
t109 = cos(pkin(12));
t111 = cos(pkin(6));
t93 = t111 * t109;
t79 = t117 * t106 - t120 * t93;
t107 = sin(pkin(7));
t108 = sin(pkin(6));
t90 = t108 * t107;
t138 = t79 * t110 + t120 * t90;
t92 = t111 * t106;
t46 = t117 * t109 + t120 * t92;
t62 = sin(qJ(3));
t30 = -t46 * t119 + t138 * t62;
t61 = sin(qJ(4));
t91 = t110 * t108;
t68 = -t79 * t107 + t120 * t91;
t16 = t30 * t118 + t68 * t61;
t27 = t138 * t119 + t46 * t62;
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t1 = t16 * t60 + t27 * t63;
t141 = t16 * t63 - t27 * t60;
t15 = -t68 * t118 + t30 * t61;
t74 = t120 * t106 + t117 * t93;
t64 = t74 * t107 + t117 * t91;
t135 = t74 * t110 - t117 * t90;
t134 = t107 * t111 + t109 * t91;
t47 = t120 * t109 - t117 * t92;
t32 = t47 * t119 - t135 * t62;
t18 = t32 * t118 + t64 * t61;
t89 = t108 * t106;
t38 = t119 * t89 + t134 * t62;
t73 = -t109 * t90 + t111 * t110;
t26 = t38 * t118 + t73 * t61;
t132 = g(1) * t18 - g(2) * t16 + g(3) * t26;
t17 = -t64 * t118 + t32 * t61;
t25 = -t73 * t118 + t38 * t61;
t131 = g(1) * t17 - g(2) * t15 + g(3) * t25;
t31 = t135 * t119 + t47 * t62;
t37 = -t134 * t119 + t62 * t89;
t130 = t61 * (-g(1) * t31 - g(2) * t27 - g(3) * t37);
t122 = pkin(10) + rSges(5,3);
t121 = pkin(11) + rSges(6,3);
t116 = t30 * t60;
t115 = t32 * t60;
t114 = t38 * t60;
t94 = t108 * t117;
t113 = t120 * pkin(1) + qJ(2) * t94;
t112 = qJ(6) + pkin(11) + rSges(7,3);
t105 = t118 * pkin(4);
t104 = pkin(5) * t60 + pkin(10);
t103 = t60 * t118;
t102 = t63 * t118;
t57 = pkin(5) * t63 + pkin(4);
t101 = t118 * t57;
t21 = t27 * pkin(3);
t100 = -pkin(10) * t30 - t21;
t23 = t31 * pkin(3);
t99 = t32 * pkin(10) - t23;
t36 = t37 * pkin(3);
t98 = t38 * pkin(10) - t36;
t95 = t120 * t108;
t97 = -t117 * pkin(1) + qJ(2) * t95;
t5 = -t18 * t60 + t31 * t63;
t11 = -t26 * t60 + t37 * t63;
t85 = -t118 * rSges(5,1) + rSges(5,2) * t61;
t70 = -t46 * pkin(2) + pkin(9) * t68 + t97;
t67 = t30 * pkin(3) + t70;
t66 = t47 * pkin(2) + pkin(9) * t64 + t113;
t65 = t32 * pkin(3) + t66;
t20 = -t37 * t102 + t114;
t19 = t37 * t103 + t38 * t63;
t12 = -t26 * t63 - t37 * t60;
t10 = -t31 * t102 + t115;
t9 = t31 * t103 + t32 * t63;
t8 = -t27 * t102 - t116;
t7 = t27 * t103 - t30 * t63;
t6 = t18 * t63 + t31 * t60;
t2 = [-m(2) * (g(1) * (-t117 * rSges(2,1) - t120 * rSges(2,2)) + g(2) * (t120 * rSges(2,1) - t117 * rSges(2,2))) - m(3) * (g(1) * (-t46 * rSges(3,1) + t79 * rSges(3,2) + rSges(3,3) * t95 + t97) + g(2) * (t47 * rSges(3,1) - t74 * rSges(3,2) + rSges(3,3) * t94 + t113)) - m(4) * (g(1) * (t30 * rSges(4,1) + rSges(4,2) * t27 + t68 * rSges(4,3) + t70) + g(2) * (t32 * rSges(4,1) - t31 * rSges(4,2) + t64 * rSges(4,3) + t66)) - m(5) * (g(1) * (t16 * rSges(5,1) - t15 * rSges(5,2) - t122 * t27 + t67) + g(2) * (t18 * rSges(5,1) - t17 * rSges(5,2) + t122 * t31 + t65)) - m(6) * (g(1) * (rSges(6,1) * t141 - rSges(6,2) * t1 + t16 * pkin(4) - pkin(10) * t27 + t121 * t15 + t67) + g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t18 * pkin(4) + t31 * pkin(10) + t121 * t17 + t65)) - m(7) * (g(1) * (rSges(7,1) * t141 - rSges(7,2) * t1 - t104 * t27 + t112 * t15 + t16 * t57 + t67) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t104 * t31 + t112 * t17 + t18 * t57 + t65)) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t94 - g(2) * t95 + g(3) * t111) -m(4) * (g(1) * (-rSges(4,1) * t31 - rSges(4,2) * t32) + g(2) * (-rSges(4,1) * t27 + rSges(4,2) * t30) + g(3) * (-rSges(4,1) * t37 - rSges(4,2) * t38)) - m(5) * (g(1) * (t122 * t32 + t85 * t31 - t23) + g(2) * (-t122 * t30 + t85 * t27 - t21) + g(3) * (t122 * t38 + t85 * t37 - t36)) + (-g(1) * (t10 * rSges(7,1) + t9 * rSges(7,2) + pkin(5) * t115 - t31 * t101 + t99) - g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) - pkin(5) * t116 - t27 * t101 + t100) - g(3) * (t20 * rSges(7,1) + t19 * rSges(7,2) + pkin(5) * t114 - t37 * t101 + t98) - t112 * t130) * m(7) + (-g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) - t31 * t105 + t99) - g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) - t27 * t105 + t100) - g(3) * (t20 * rSges(6,1) + t19 * rSges(6,2) - t37 * t105 + t98) - t121 * t130) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (rSges(5,1) * t15 + rSges(5,2) * t16) + g(3) * (-rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (t131 * (-rSges(6,1) * t63 + rSges(6,2) * t60 - pkin(4)) + t132 * t121) - m(7) * (t131 * (-rSges(7,1) * t63 + rSges(7,2) * t60 - t57) + t132 * t112) -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t141) + g(3) * (rSges(6,1) * t11 + rSges(6,2) * t12)) + (-g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) - g(2) * (t1 * rSges(7,1) + rSges(7,2) * t141) - g(3) * (t11 * rSges(7,1) + t12 * rSges(7,2)) - (g(1) * t5 + g(2) * t1 + g(3) * t11) * pkin(5)) * m(7), -m(7) * t131];
taug  = t2(:);
