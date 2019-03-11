% Calculate Gravitation load on the joints for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:58:00
% EndTime: 2019-03-09 11:58:04
% DurationCPUTime: 1.36s
% Computational Cost: add. (889->197), mult. (2151->288), div. (0->0), fcn. (2705->12), ass. (0->92)
t58 = sin(pkin(11));
t63 = sin(qJ(2));
t67 = cos(qJ(2));
t99 = cos(pkin(11));
t45 = -t67 * t58 - t63 * t99;
t64 = sin(qJ(1));
t68 = cos(qJ(1));
t59 = cos(pkin(6));
t74 = -t63 * t58 + t67 * t99;
t71 = t59 * t74;
t29 = t45 * t68 - t64 * t71;
t107 = t64 * t67;
t111 = t63 * t68;
t42 = -t59 * t107 - t111;
t73 = t42 * pkin(2);
t128 = t29 * pkin(3) + t73;
t104 = t67 * t68;
t108 = t64 * t63;
t40 = -t59 * t104 + t108;
t101 = t45 * t59;
t25 = t68 * t101 - t64 * t74;
t62 = sin(qJ(4));
t66 = cos(qJ(4));
t98 = sin(pkin(6));
t85 = t68 * t98;
t14 = -t25 * t66 - t62 * t85;
t26 = t64 * t45 + t68 * t71;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t3 = t14 * t61 + t26 * t65;
t127 = -t14 * t65 + t26 * t61;
t78 = t99 * t98;
t86 = t67 * t98;
t38 = t58 * t86 + t63 * t78;
t32 = t38 * t66 + t59 * t62;
t126 = g(3) * t32;
t30 = t64 * t101 + t68 * t74;
t13 = -t25 * t62 + t66 * t85;
t87 = t64 * t98;
t17 = t30 * t62 - t66 * t87;
t31 = t38 * t62 - t59 * t66;
t125 = g(1) * t17 + g(2) * t13 + g(3) * t31;
t124 = pkin(4) * t66;
t88 = t63 * t98;
t37 = t58 * t88 - t67 * t78;
t121 = g(3) * t37;
t120 = rSges(5,3) + pkin(9);
t119 = rSges(6,3) + pkin(10);
t118 = t25 * t61;
t115 = t30 * t61;
t114 = t38 * t61;
t56 = pkin(5) * t65 + pkin(4);
t113 = t56 * t66;
t112 = t61 * t66;
t57 = pkin(2) * t67 + pkin(1);
t109 = t64 * t57;
t106 = t65 * t66;
t102 = rSges(7,3) + qJ(6) + pkin(10);
t53 = pkin(2) * t86;
t100 = -t37 * pkin(3) + t53;
t95 = g(1) * t119;
t94 = g(2) * t119;
t93 = pkin(5) * t61 + pkin(9);
t92 = g(1) * t102;
t91 = g(2) * t102;
t90 = t98 * pkin(8);
t84 = t40 * pkin(2);
t83 = t38 * pkin(9) + t100;
t39 = t59 * t63 * pkin(2) - t98 * qJ(3) - t90;
t50 = t68 * t57;
t82 = t30 * pkin(3) - t64 * t39 + t50;
t81 = rSges(5,1) * t66 - rSges(5,2) * t62;
t18 = t30 * t66 + t62 * t87;
t5 = -t18 * t61 - t29 * t65;
t11 = -t32 * t61 + t37 * t65;
t80 = t98 * rSges(4,3) - t39;
t79 = t26 * pkin(3) - t84;
t75 = pkin(3) * t25 - t68 * t39 - t109;
t72 = -t25 * pkin(9) + t79;
t70 = t98 * rSges(3,3) + t90;
t69 = pkin(9) * t30 + t128;
t43 = -t59 * t108 + t104;
t41 = -t59 * t111 - t107;
t20 = -t37 * t106 + t114;
t19 = t37 * t112 + t38 * t65;
t12 = -t32 * t65 - t37 * t61;
t10 = t29 * t106 + t115;
t9 = -t29 * t112 + t30 * t65;
t8 = t26 * t106 - t118;
t7 = -t26 * t112 - t25 * t65;
t6 = t18 * t65 - t29 * t61;
t1 = [-m(2) * (g(1) * (-t64 * rSges(2,1) - rSges(2,2) * t68) + g(2) * (rSges(2,1) * t68 - t64 * rSges(2,2))) - m(3) * (g(1) * (t41 * rSges(3,1) + t40 * rSges(3,2) - t64 * pkin(1) + t70 * t68) + g(2) * (t43 * rSges(3,1) + t42 * rSges(3,2) + t68 * pkin(1) + t70 * t64)) - m(4) * (g(1) * (rSges(4,1) * t25 - t26 * rSges(4,2) + t80 * t68 - t109) + g(2) * (t30 * rSges(4,1) + t29 * rSges(4,2) + t80 * t64 + t50)) - m(5) * (g(1) * (-rSges(5,1) * t14 + rSges(5,2) * t13 + t120 * t26 + t75) + g(2) * (rSges(5,1) * t18 - rSges(5,2) * t17 - t120 * t29 + t82)) - m(6) * (g(1) * (rSges(6,1) * t127 + t3 * rSges(6,2) - pkin(4) * t14 + t26 * pkin(9) - t119 * t13 + t75) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 + pkin(4) * t18 - pkin(9) * t29 + t119 * t17 + t82)) - m(7) * (g(1) * (rSges(7,1) * t127 + t3 * rSges(7,2) - t102 * t13 - t14 * t56 + t93 * t26 + t75) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t102 * t17 + t18 * t56 - t93 * t29 + t82)) -m(3) * (g(1) * (rSges(3,1) * t42 - t43 * rSges(3,2)) + g(2) * (-t40 * rSges(3,1) + t41 * rSges(3,2)) + g(3) * (rSges(3,1) * t86 - rSges(3,2) * t88)) - m(4) * (g(1) * (t29 * rSges(4,1) - rSges(4,2) * t30 + t73) + g(2) * (t26 * rSges(4,1) + t25 * rSges(4,2) - t84) + g(3) * (-t37 * rSges(4,1) - t38 * rSges(4,2) + t53)) - m(5) * (g(1) * (t120 * t30 + t81 * t29 + t128) + g(2) * (-t120 * t25 + t81 * t26 + t79) + g(3) * (t120 * t38 - t81 * t37 + t100)) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t29 * t124 + t69) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t26 * t124 + t72) + g(3) * (t20 * rSges(6,1) + t19 * rSges(6,2) - t37 * t124 + t83) + (-t119 * t121 + t26 * t94 + t29 * t95) * t62) - m(7) * (g(1) * (t10 * rSges(7,1) + t9 * rSges(7,2) + pkin(5) * t115 + t29 * t113 + t69) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) - pkin(5) * t118 + t26 * t113 + t72) + g(3) * (t20 * rSges(7,1) + t19 * rSges(7,2) + pkin(5) * t114 - t37 * t113 + t83) + (-t102 * t121 + t26 * t91 + t29 * t92) * t62) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t87 - g(2) * t85 + g(3) * t59) -m(5) * (g(1) * (-rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t13 - rSges(5,2) * t14) + g(3) * (-rSges(5,1) * t31 - rSges(5,2) * t32)) - m(6) * (t119 * t126 + t14 * t94 + t18 * t95 + t125 * (-rSges(6,1) * t65 + rSges(6,2) * t61 - pkin(4))) - m(7) * (t102 * t126 + t14 * t91 + t18 * t92 + t125 * (-rSges(7,1) * t65 + rSges(7,2) * t61 - t56)) -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t3 + rSges(6,2) * t127) + g(3) * (rSges(6,1) * t11 + rSges(6,2) * t12)) + (-g(1) * (rSges(7,1) * t5 - rSges(7,2) * t6) - g(2) * (-rSges(7,1) * t3 + rSges(7,2) * t127) - g(3) * (rSges(7,1) * t11 + rSges(7,2) * t12) - (g(1) * t5 - g(2) * t3 + g(3) * t11) * pkin(5)) * m(7), -m(7) * t125];
taug  = t1(:);
