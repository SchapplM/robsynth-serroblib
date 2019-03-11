% Calculate Gravitation load on the joints for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:53
% EndTime: 2019-03-09 21:41:56
% DurationCPUTime: 1.38s
% Computational Cost: add. (838->216), mult. (2060->305), div. (0->0), fcn. (2527->10), ass. (0->89)
t75 = sin(qJ(4));
t78 = cos(qJ(4));
t136 = pkin(4) * t78 + qJ(5) * t75;
t128 = cos(qJ(1));
t74 = sin(pkin(6));
t102 = t74 * t128;
t127 = sin(qJ(1));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t110 = cos(pkin(6));
t94 = t110 * t128;
t57 = t127 * t80 + t77 * t94;
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t30 = -t76 * t102 + t57 * t79;
t56 = t127 * t77 - t80 * t94;
t6 = t30 * t75 - t56 * t78;
t7 = t30 * t78 + t56 * t75;
t113 = rSges(7,2) + qJ(5);
t111 = rSges(7,3) + qJ(6);
t135 = pkin(3) * t79;
t133 = g(3) * t74;
t132 = rSges(6,1) + pkin(10);
t131 = rSges(7,1) + pkin(5);
t130 = rSges(4,3) + pkin(9);
t129 = rSges(5,3) + pkin(10);
t126 = rSges(7,2) * t75;
t124 = t56 * t76;
t93 = t110 * t127;
t58 = t128 * t77 + t80 * t93;
t122 = t58 * t76;
t121 = t74 * t77;
t120 = t74 * t80;
t119 = t75 * t79;
t118 = t78 * t79;
t117 = t79 * t80;
t101 = t74 * t127;
t116 = t128 * pkin(1) + pkin(8) * t101;
t115 = pkin(2) * t120 + pkin(9) * t121;
t112 = rSges(6,3) + qJ(5);
t109 = pkin(10) + t131;
t108 = t76 * t120;
t107 = t75 * t120;
t100 = -t79 * t102 - t57 * t76;
t23 = t100 * pkin(3);
t106 = t136 * t100 + t23;
t59 = t128 * t80 - t77 * t93;
t33 = -t101 * t79 + t59 * t76;
t25 = t33 * pkin(3);
t105 = -t136 * t33 - t25;
t54 = t110 * t79 - t121 * t76;
t49 = t54 * pkin(3);
t104 = t136 * t54 + t49;
t103 = t59 * pkin(2) + t116;
t99 = t74 * pkin(3) * t117 + pkin(10) * t108 + t115;
t98 = g(1) * t109;
t97 = g(2) * t109;
t96 = -pkin(1) * t127 + pkin(8) * t102;
t37 = (t117 * t78 + t75 * t77) * t74;
t95 = t37 * pkin(4) + t99;
t92 = rSges(4,1) * t79 - rSges(4,2) * t76;
t91 = rSges(5,1) * t78 - rSges(5,2) * t75;
t90 = rSges(6,2) * t78 - rSges(6,3) * t75;
t50 = t56 * pkin(2);
t89 = t57 * pkin(9) - pkin(10) * t124 - t56 * t135 - t50;
t52 = t58 * pkin(2);
t88 = t59 * pkin(9) - pkin(10) * t122 - t58 * t135 - t52;
t87 = -t57 * pkin(2) + t96;
t15 = -t118 * t56 + t57 * t75;
t86 = t15 * pkin(4) + t89;
t17 = -t118 * t58 + t59 * t75;
t85 = t17 * pkin(4) + t88;
t34 = t101 * t76 + t59 * t79;
t84 = t34 * pkin(3) + pkin(9) * t58 + t103;
t11 = t34 * t78 + t58 * t75;
t83 = t11 * pkin(4) + t84;
t82 = -pkin(3) * t30 - t56 * pkin(9) + t87;
t81 = -pkin(4) * t7 + t82;
t55 = t110 * t76 + t121 * t79;
t36 = t107 * t79 - t121 * t78;
t28 = t55 * t78 - t107;
t27 = t120 * t78 + t55 * t75;
t22 = t27 * pkin(4);
t16 = -t119 * t58 - t59 * t78;
t14 = -t119 * t56 - t57 * t78;
t10 = t34 * t75 - t58 * t78;
t4 = t10 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t127 - rSges(2,2) * t128) + g(2) * (rSges(2,1) * t128 - rSges(2,2) * t127)) - m(3) * (g(1) * (-t57 * rSges(3,1) + t56 * rSges(3,2) + rSges(3,3) * t102 + t96) + g(2) * (t59 * rSges(3,1) - t58 * rSges(3,2) + rSges(3,3) * t101 + t116)) - m(4) * (g(1) * (-rSges(4,1) * t30 - rSges(4,2) * t100 - t130 * t56 + t87) + g(2) * (rSges(4,1) * t34 - rSges(4,2) * t33 + t130 * t58 + t103)) - m(5) * (g(1) * (-rSges(5,1) * t7 + rSges(5,2) * t6 + t100 * t129 + t82) + g(2) * (rSges(5,1) * t11 - rSges(5,2) * t10 + t129 * t33 + t84)) - m(6) * (g(1) * (rSges(6,2) * t7 + t100 * t132 - t112 * t6 + t81) + g(2) * (-rSges(6,2) * t11 + t10 * t112 + t132 * t33 + t83)) - m(7) * (g(1) * (-t111 * t7 - t113 * t6 + t81) + g(2) * (t113 * t10 + t111 * t11 + t83) + t33 * t97 + t100 * t98) -m(3) * (g(1) * (-rSges(3,1) * t58 - rSges(3,2) * t59) + g(2) * (-rSges(3,1) * t56 - rSges(3,2) * t57) + (rSges(3,1) * t80 - rSges(3,2) * t77) * t133) - m(4) * (g(1) * (t130 * t59 - t58 * t92 - t52) + g(2) * (t130 * t57 - t56 * t92 - t50) + g(3) * t115 + (rSges(4,3) * t77 + t80 * t92) * t133) - m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t16 - rSges(5,3) * t122 + t88) + g(2) * (rSges(5,1) * t15 - rSges(5,2) * t14 - rSges(5,3) * t124 + t89) + g(3) * (t37 * rSges(5,1) - t36 * rSges(5,2) + rSges(5,3) * t108 + t99)) - m(6) * (g(1) * (-rSges(6,1) * t122 - rSges(6,2) * t17 + t112 * t16 + t85) + g(2) * (-rSges(6,1) * t124 - rSges(6,2) * t15 + t112 * t14 + t86) + g(3) * (rSges(6,1) * t108 - t37 * rSges(6,2) + t112 * t36 + t95)) - m(7) * (g(1) * (t111 * t17 + t113 * t16 + t85) + g(2) * (t111 * t15 + t113 * t14 + t86) + g(3) * (t111 * t37 + t113 * t36 + t95) + (-g(1) * t58 - g(2) * t56 + g(3) * t120) * t76 * t131) -m(4) * (g(1) * (-rSges(4,1) * t33 - rSges(4,2) * t34) + g(2) * (rSges(4,1) * t100 - rSges(4,2) * t30) + g(3) * (rSges(4,1) * t54 - rSges(4,2) * t55)) - m(5) * (g(1) * (t129 * t34 - t33 * t91 - t25) + g(2) * (t100 * t91 + t129 * t30 + t23) + g(3) * (t129 * t55 + t54 * t91 + t49)) - m(6) * (g(1) * (t132 * t34 + t33 * t90 + t105) + g(2) * (-t100 * t90 + t132 * t30 + t106) + g(3) * (t132 * t55 - t54 * t90 + t104)) + (-g(1) * (-t126 * t33 + t105) - g(2) * (t100 * t126 + t106) - g(3) * (t126 * t54 + t104) - g(3) * t109 * t55 - t34 * t98 - t30 * t97 - (-g(1) * t33 + g(2) * t100 + g(3) * t54) * t78 * t111) * m(7), -m(5) * (g(1) * (-rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7) + g(3) * (-rSges(5,1) * t27 - rSges(5,2) * t28)) - m(6) * (g(1) * (rSges(6,2) * t10 + t11 * t112 - t4) + g(2) * (rSges(6,2) * t6 + t112 * t7 - t2) + g(3) * (rSges(6,2) * t27 + t112 * t28 - t22)) - m(7) * (g(1) * (-t10 * t111 + t11 * t113 - t4) + g(2) * (-t111 * t6 + t113 * t7 - t2) + g(3) * (-t111 * t27 + t113 * t28 - t22)) (-m(6) - m(7)) * (g(1) * t10 + g(2) * t6 + g(3) * t27) -m(7) * (g(1) * t11 + g(2) * t7 + g(3) * t28)];
taug  = t1(:);
