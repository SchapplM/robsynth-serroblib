% Calculate Gravitation load on the joints for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:53
% EndTime: 2019-03-08 22:22:57
% DurationCPUTime: 1.19s
% Computational Cost: add. (989->198), mult. (2346->311), div. (0->0), fcn. (2947->16), ass. (0->91)
t115 = cos(qJ(3));
t101 = cos(pkin(7));
t63 = sin(pkin(7));
t114 = sin(qJ(2));
t116 = cos(qJ(2));
t100 = cos(pkin(12));
t102 = cos(pkin(6));
t86 = t102 * t100;
t98 = sin(pkin(12));
t71 = t114 * t98 - t116 * t86;
t99 = sin(pkin(6));
t83 = t100 * t99;
t131 = t71 * t101 + t63 * t83;
t47 = t114 * t86 + t116 * t98;
t67 = sin(qJ(3));
t16 = t115 * t131 + t47 * t67;
t85 = t102 * t98;
t72 = t100 * t114 + t116 * t85;
t82 = t99 * t98;
t130 = t72 * t101 - t63 * t82;
t48 = t100 * t116 - t114 * t85;
t18 = t115 * t130 + t48 * t67;
t84 = t101 * t99;
t129 = t102 * t63 + t116 * t84;
t88 = t99 * t114;
t32 = -t115 * t129 + t67 * t88;
t126 = g(1) * t18 + g(2) * t16 + g(3) * t32;
t117 = rSges(7,3) + pkin(11);
t128 = g(1) * t48 + g(2) * t47;
t61 = pkin(13) + qJ(5);
t60 = cos(t61);
t120 = t60 * pkin(5);
t119 = t63 * pkin(9);
t59 = sin(t61);
t113 = t59 * t63;
t112 = t60 * t63;
t66 = sin(qJ(6));
t111 = t60 * t66;
t68 = cos(qJ(6));
t110 = t60 * t68;
t62 = sin(pkin(13));
t109 = t62 * t63;
t64 = cos(pkin(13));
t108 = t63 * t64;
t17 = t115 * t47 - t131 * t67;
t58 = pkin(4) * t64 + pkin(3);
t65 = -pkin(10) - qJ(4);
t107 = -t16 * t58 - t17 * t65;
t19 = t48 * t115 - t130 * t67;
t106 = -t18 * t58 - t19 * t65;
t33 = t115 * t88 + t129 * t67;
t105 = -t32 * t58 - t33 * t65;
t79 = t63 * t88;
t89 = t116 * t99;
t104 = pkin(2) * t89 + pkin(9) * t79;
t103 = qJ(4) + rSges(5,3);
t97 = -m(5) - m(6) - m(7);
t90 = t101 * t115;
t24 = t47 * t90 - t67 * t71;
t93 = t67 * t101;
t25 = -t115 * t71 - t47 * t93;
t44 = t71 * pkin(2);
t96 = -t24 * t65 + t25 * t58 - t44;
t26 = t48 * t90 - t67 * t72;
t27 = -t115 * t72 - t48 * t93;
t45 = t72 * pkin(2);
t95 = -t26 * t65 + t27 * t58 - t45;
t74 = t114 * t84;
t42 = t115 * t74 + t67 * t89;
t43 = t115 * t89 - t67 * t74;
t76 = t62 * t79;
t91 = pkin(4) * t76 - t42 * t65 + t43 * t58 + t104;
t87 = -rSges(6,1) * t60 + rSges(6,2) * t59;
t80 = t68 * rSges(7,1) - t66 * rSges(7,2) + pkin(5);
t73 = t128 * (pkin(4) * t62 + pkin(9)) * t63;
t46 = t101 * t102 - t63 * t89;
t35 = t101 * t82 + t63 * t72;
t34 = -t101 * t83 + t63 * t71;
t29 = t43 * t60 + t59 * t79;
t28 = t43 * t59 - t60 * t79;
t13 = t33 * t60 + t46 * t59;
t12 = -t33 * t59 + t46 * t60;
t9 = t113 * t48 + t27 * t60;
t8 = -t112 * t48 + t27 * t59;
t7 = t113 * t47 + t25 * t60;
t6 = -t112 * t47 + t25 * t59;
t5 = t19 * t60 + t35 * t59;
t4 = -t19 * t59 + t35 * t60;
t3 = t17 * t60 + t34 * t59;
t2 = -t17 * t59 + t34 * t60;
t1 = [(-m(2) - m(3) - m(4) + t97) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t72 - t48 * rSges(3,2)) + g(2) * (-rSges(3,1) * t71 - t47 * rSges(3,2)) + g(3) * (rSges(3,1) * t89 - rSges(3,2) * t88)) - m(4) * (g(1) * (rSges(4,1) * t27 - rSges(4,2) * t26 - t45) + g(2) * (rSges(4,1) * t25 - rSges(4,2) * t24 - t44) + g(3) * (t43 * rSges(4,1) - t42 * rSges(4,2) + t104) + (rSges(4,3) * g(3) * t88 + t128 * (rSges(4,3) + pkin(9))) * t63) - m(5) * (g(1) * (t27 * pkin(3) - t45 + t48 * t119 + (t109 * t48 + t27 * t64) * rSges(5,1) + (t108 * t48 - t27 * t62) * rSges(5,2) + t103 * t26) + g(2) * (t25 * pkin(3) - t44 + t47 * t119 + (t109 * t47 + t25 * t64) * rSges(5,1) + (t108 * t47 - t25 * t62) * rSges(5,2) + t103 * t24) + g(3) * (t43 * pkin(3) + (t43 * t64 + t76) * rSges(5,1) + (-t43 * t62 + t64 * t79) * rSges(5,2) + t103 * t42 + t104)) - m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t8 + rSges(6,3) * t26 + t95) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + rSges(6,3) * t24 + t96) + g(3) * (rSges(6,1) * t29 - rSges(6,2) * t28 + t42 * rSges(6,3) + t91) + t73) - m(7) * (g(1) * (t9 * pkin(5) + (t26 * t66 + t68 * t9) * rSges(7,1) + (t26 * t68 - t9 * t66) * rSges(7,2) + t95 + t117 * t8) + g(2) * (t7 * pkin(5) + (t24 * t66 + t68 * t7) * rSges(7,1) + (t24 * t68 - t7 * t66) * rSges(7,2) + t96 + t117 * t6) + g(3) * (t29 * pkin(5) + (t29 * t68 + t42 * t66) * rSges(7,1) + (-t29 * t66 + t42 * t68) * rSges(7,2) + t117 * t28 + t91) + t73) -m(4) * (g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(3) * (-rSges(4,1) * t32 - rSges(4,2) * t33)) - m(5) * (t126 * (-t64 * rSges(5,1) + t62 * rSges(5,2) - pkin(3)) + (g(1) * t19 + g(2) * t17 + g(3) * t33) * t103) - m(6) * (g(1) * (rSges(6,3) * t19 + t18 * t87 + t106) + g(2) * (rSges(6,3) * t17 + t16 * t87 + t107) + g(3) * (rSges(6,3) * t33 + t32 * t87 + t105)) + (-g(1) * (-t18 * t120 + (-t110 * t18 + t19 * t66) * rSges(7,1) + (t111 * t18 + t19 * t68) * rSges(7,2) + t106) - g(2) * (-t16 * t120 + (-t110 * t16 + t17 * t66) * rSges(7,1) + (t111 * t16 + t17 * t68) * rSges(7,2) + t107) - g(3) * (-t32 * t120 + (-t110 * t32 + t33 * t66) * rSges(7,1) + (t111 * t32 + t33 * t68) * rSges(7,2) + t105) + t126 * t59 * t117) * m(7), t97 * t126, -m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(3) * (rSges(6,1) * t12 - rSges(6,2) * t13)) - m(7) * (g(1) * (t117 * t5 + t4 * t80) + (t117 * t13 + t80 * t12) * g(3) + (t117 * t3 + t80 * t2) * g(2)) -m(7) * (g(1) * ((t18 * t68 - t5 * t66) * rSges(7,1) + (-t18 * t66 - t5 * t68) * rSges(7,2)) + g(2) * ((t16 * t68 - t3 * t66) * rSges(7,1) + (-t16 * t66 - t3 * t68) * rSges(7,2)) + g(3) * ((-t13 * t66 + t32 * t68) * rSges(7,1) + (-t13 * t68 - t32 * t66) * rSges(7,2)))];
taug  = t1(:);
