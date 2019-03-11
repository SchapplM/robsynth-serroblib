% Calculate Gravitation load on the joints for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:52
% EndTime: 2019-03-09 22:38:55
% DurationCPUTime: 1.24s
% Computational Cost: add. (673->189), mult. (904->248), div. (0->0), fcn. (958->10), ass. (0->85)
t57 = sin(qJ(1));
t60 = cos(qJ(2));
t104 = t57 * t60;
t53 = qJ(3) + qJ(4);
t48 = sin(t53);
t49 = cos(t53);
t61 = cos(qJ(1));
t27 = t104 * t48 + t49 * t61;
t102 = t61 * t48;
t28 = t104 * t49 - t102;
t54 = sin(qJ(6));
t58 = cos(qJ(6));
t126 = t27 * t58 - t28 * t54;
t74 = t27 * t54 + t28 * t58;
t131 = t126 * rSges(7,1) - t74 * rSges(7,2);
t29 = t102 * t60 - t57 * t49;
t103 = t60 * t61;
t30 = t103 * t49 + t48 * t57;
t6 = t29 * t58 - t30 * t54;
t7 = t29 * t54 + t30 * t58;
t130 = rSges(7,1) * t6 - rSges(7,2) * t7;
t55 = sin(qJ(3));
t59 = cos(qJ(3));
t33 = -t55 * t103 + t57 * t59;
t56 = sin(qJ(2));
t72 = t48 * t54 + t49 * t58;
t73 = t48 * t58 - t49 * t54;
t129 = (rSges(7,1) * t73 - rSges(7,2) * t72) * t56;
t47 = pkin(3) * t59 + pkin(2);
t38 = t60 * t47;
t88 = -pkin(1) - t38;
t128 = -t27 * rSges(6,1) + t28 * rSges(6,3);
t127 = -t29 * rSges(6,1) + t30 * rSges(6,3);
t116 = g(2) * t57;
t125 = g(1) * t61 + t116;
t124 = -t27 * pkin(5) - t131;
t123 = -t29 * pkin(5) - t130;
t112 = rSges(4,3) + pkin(8);
t122 = t60 * pkin(2) + t112 * t56;
t121 = -pkin(4) - pkin(5);
t119 = pkin(3) * t55;
t118 = g(1) * t57;
t115 = g(3) * t56;
t113 = -rSges(6,1) - pkin(4);
t111 = rSges(7,3) + pkin(10);
t110 = rSges(3,2) * t56;
t108 = t49 * t56;
t107 = t55 * t57;
t106 = t55 * t61;
t62 = -pkin(9) - pkin(8);
t101 = rSges(6,2) - t62;
t100 = rSges(5,3) - t62;
t99 = -t27 * rSges(5,1) - t28 * rSges(5,2);
t98 = -t29 * rSges(5,1) - t30 * rSges(5,2);
t36 = qJ(5) * t108;
t97 = rSges(6,3) * t108 + t36;
t96 = t61 * pkin(1) + t57 * pkin(7);
t95 = rSges(6,3) + qJ(5);
t94 = t36 - t129;
t93 = -t62 - t111;
t51 = t61 * pkin(7);
t91 = t57 * t56 * t62 + pkin(3) * t106 + t51;
t90 = t121 * t48;
t89 = t113 * t48;
t87 = t101 * t61;
t86 = t100 * t61;
t84 = -t27 * pkin(4) + qJ(5) * t28;
t83 = -t29 * pkin(4) + qJ(5) * t30;
t82 = pkin(3) * t107 + t47 * t103 + t96;
t81 = g(3) * (t38 + (t49 * pkin(4) + t48 * qJ(5)) * t60);
t80 = t93 * t61;
t79 = t30 * pkin(4) + t82;
t78 = rSges(3,1) * t60 - t110;
t76 = rSges(5,1) * t49 - rSges(5,2) * t48;
t75 = -rSges(5,1) * t48 - rSges(5,2) * t49;
t71 = t33 * pkin(3);
t70 = rSges(4,1) * t59 - rSges(4,2) * t55 + pkin(2);
t31 = t104 * t55 + t59 * t61;
t68 = -t28 * pkin(4) - t27 * qJ(5) + t91;
t67 = t31 * pkin(3);
t66 = t71 + t83;
t64 = -t67 + t84;
t34 = t103 * t59 + t107;
t32 = -t104 * t59 + t106;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t57 - rSges(2,2) * t61) + g(2) * (rSges(2,1) * t61 - rSges(2,2) * t57)) - m(3) * (g(1) * (rSges(3,3) * t61 + t51) + g(2) * (rSges(3,1) * t103 - t110 * t61 + t96) + (g(1) * (-pkin(1) - t78) + g(2) * rSges(3,3)) * t57) - m(4) * (g(1) * (rSges(4,1) * t32 + rSges(4,2) * t31 + t51) + (-pkin(1) - t122) * t118 + (rSges(4,1) * t34 + rSges(4,2) * t33 + t122 * t61 + t96) * g(2)) - m(5) * (g(1) * (-rSges(5,1) * t28 + rSges(5,2) * t27 + t91) + g(2) * (t30 * rSges(5,1) - t29 * rSges(5,2) + t56 * t86 + t82) + (-t56 * rSges(5,3) + t88) * t118) - m(6) * (g(1) * (-rSges(6,1) * t28 - rSges(6,3) * t27 + t68) + g(2) * (t30 * rSges(6,1) + t29 * t95 + t56 * t87 + t79) + (-t56 * rSges(6,2) + t88) * t118) - m(7) * (g(1) * (-t74 * rSges(7,1) - rSges(7,2) * t126 - t28 * pkin(5) + t88 * t57 + t68) + g(2) * (t7 * rSges(7,1) + t6 * rSges(7,2) + t30 * pkin(5) + t29 * qJ(5) + t79) + (g(2) * t80 + t111 * t118) * t56) -m(3) * (g(3) * t78 + t125 * (-rSges(3,1) * t56 - rSges(3,2) * t60)) - m(4) * ((g(3) * t70 + t125 * t112) * t60 + (g(3) * t112 - t125 * t70) * t56) - m(5) * (g(3) * t38 + (g(1) * t86 + g(3) * t76 + t100 * t116) * t60 + (g(3) * t100 + t125 * (-t47 - t76)) * t56) - m(6) * (t81 + (g(3) * (rSges(6,1) * t49 + rSges(6,3) * t48) + g(1) * t87 + t101 * t116) * t60 + (g(3) * t101 + t125 * (t113 * t49 - t48 * t95 - t47)) * t56) - m(7) * (t81 + (g(3) * (rSges(7,1) * t72 + rSges(7,2) * t73 + t49 * pkin(5)) + g(1) * t80 + t93 * t116) * t60 + (g(3) * t93 + t125 * (-t47 + (-rSges(7,1) * t54 - rSges(7,2) * t58 - qJ(5)) * t48 + (-rSges(7,1) * t58 + rSges(7,2) * t54 + t121) * t49)) * t56) -m(4) * (g(1) * (rSges(4,1) * t33 - rSges(4,2) * t34) + g(2) * (-rSges(4,1) * t31 + rSges(4,2) * t32) + (-rSges(4,1) * t55 - rSges(4,2) * t59) * t115) - m(5) * (g(1) * (t71 + t98) + g(2) * (-t67 + t99) + (t75 - t119) * t115) - m(6) * (g(1) * (t66 + t127) + g(2) * (t64 + t128) + g(3) * t97 + (t89 - t119) * t115) - m(7) * (g(1) * (t66 + t123) + g(2) * (t64 + t124) + g(3) * t94 + (t90 - t119) * t115) -m(5) * (g(1) * t98 + g(2) * t99 + t75 * t115) - m(6) * (g(1) * (t83 + t127) + g(2) * (t84 + t128) + g(3) * (t56 * t89 + t97)) - m(7) * (g(1) * (t83 + t123) + g(2) * (t84 + t124) + g(3) * (t56 * t90 + t94)) (-m(6) - m(7)) * (g(1) * t29 + g(2) * t27 + t115 * t48) -m(7) * (g(1) * t130 + g(2) * t131 + g(3) * t129)];
taug  = t1(:);
