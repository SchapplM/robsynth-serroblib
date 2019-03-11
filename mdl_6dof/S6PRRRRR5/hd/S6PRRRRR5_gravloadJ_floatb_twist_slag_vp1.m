% Calculate Gravitation load on the joints for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:56
% EndTime: 2019-03-09 01:03:59
% DurationCPUTime: 1.39s
% Computational Cost: add. (1276->224), mult. (3337->345), div. (0->0), fcn. (4249->16), ass. (0->102)
t106 = cos(pkin(13));
t107 = cos(pkin(7));
t103 = sin(pkin(13));
t119 = sin(qJ(2));
t121 = cos(qJ(2));
t108 = cos(pkin(6));
t87 = t108 * t106;
t68 = t103 * t119 - t121 * t87;
t104 = sin(pkin(7));
t105 = sin(pkin(6));
t84 = t105 * t104;
t137 = t106 * t84 + t68 * t107;
t86 = t108 * t103;
t69 = t106 * t119 + t121 * t86;
t83 = t105 * t103;
t136 = -t104 * t83 + t69 * t107;
t85 = t107 * t105;
t135 = t104 * t108 + t121 * t85;
t120 = cos(qJ(3));
t62 = sin(qJ(3));
t91 = t105 * t119;
t33 = t120 * t91 + t135 * t62;
t45 = t108 * t107 - t121 * t84;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t21 = -t33 * t61 + t45 * t64;
t46 = t103 * t121 + t119 * t87;
t18 = t46 * t120 - t137 * t62;
t34 = t68 * t104 - t106 * t85;
t7 = -t18 * t61 + t34 * t64;
t47 = t106 * t121 - t119 * t86;
t20 = t47 * t120 - t136 * t62;
t35 = t69 * t104 + t107 * t83;
t9 = -t20 * t61 + t35 * t64;
t134 = g(1) * t9 + g(2) * t7 + g(3) * t21;
t10 = t20 * t64 + t35 * t61;
t22 = t33 * t64 + t45 * t61;
t8 = t18 * t64 + t34 * t61;
t133 = g(1) * t10 + g(2) * t8 + g(3) * t22;
t17 = t137 * t120 + t46 * t62;
t19 = t136 * t120 + t47 * t62;
t32 = -t135 * t120 + t62 * t91;
t132 = (-g(1) * t19 - g(2) * t17 - g(3) * t32) * t61;
t124 = t64 * pkin(4);
t123 = rSges(5,3) + pkin(10);
t122 = rSges(6,3) + pkin(11);
t60 = sin(qJ(5));
t118 = t18 * t60;
t117 = t20 * t60;
t116 = t33 * t60;
t59 = qJ(5) + qJ(6);
t57 = sin(t59);
t115 = t57 * t64;
t58 = cos(t59);
t114 = t58 * t64;
t113 = t60 * t64;
t63 = cos(qJ(5));
t112 = t63 * t64;
t56 = pkin(5) * t63 + pkin(4);
t111 = t64 * t56;
t110 = rSges(7,3) + pkin(12) + pkin(11);
t75 = t119 * t84;
t92 = t121 * t105;
t109 = pkin(2) * t92 + pkin(9) * t75;
t76 = t119 * t85;
t42 = t120 * t92 - t62 * t76;
t102 = t42 * pkin(3) + t109;
t15 = t17 * pkin(3);
t101 = t18 * pkin(10) - t15;
t16 = t19 * pkin(3);
t100 = t20 * pkin(10) - t16;
t31 = t32 * pkin(3);
t99 = t33 * pkin(10) - t31;
t98 = t104 * pkin(9);
t97 = t61 * t104;
t96 = t62 * t107;
t95 = t64 * t104;
t94 = t17 * t63 - t60 * t8;
t93 = t107 * t120;
t90 = -rSges(5,1) * t64 + rSges(5,2) * t61;
t89 = -t10 * t60 + t19 * t63;
t88 = -t22 * t60 + t32 * t63;
t80 = t58 * rSges(7,1) - t57 * rSges(7,2) + t56;
t26 = -t68 * t120 - t46 * t96;
t43 = t68 * pkin(2);
t79 = t26 * pkin(3) + t46 * t98 - t43;
t28 = -t69 * t120 - t47 * t96;
t44 = t69 * pkin(2);
t78 = t28 * pkin(3) + t47 * t98 - t44;
t72 = t57 * rSges(7,1) + t58 * rSges(7,2) + t60 * pkin(5) + pkin(10);
t71 = t104 * rSges(4,3) + t98;
t70 = m(7) * (g(1) * ((-t10 * t57 + t19 * t58) * rSges(7,1) + (-t10 * t58 - t19 * t57) * rSges(7,2)) + g(2) * ((t17 * t58 - t57 * t8) * rSges(7,1) + (-t17 * t57 - t58 * t8) * rSges(7,2)) + g(3) * ((-t22 * t57 + t32 * t58) * rSges(7,1) + (-t22 * t58 - t32 * t57) * rSges(7,2)));
t41 = t120 * t76 + t62 * t92;
t30 = t42 * t64 + t61 * t75;
t29 = t42 * t61 - t64 * t75;
t27 = t47 * t93 - t69 * t62;
t25 = t46 * t93 - t68 * t62;
t14 = t28 * t64 + t47 * t97;
t13 = t28 * t61 - t47 * t95;
t12 = t26 * t64 + t46 * t97;
t11 = t26 * t61 - t46 * t95;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-t69 * rSges(3,1) - t47 * rSges(3,2)) + g(2) * (-t68 * rSges(3,1) - t46 * rSges(3,2)) + g(3) * (rSges(3,1) * t92 - rSges(3,2) * t91)) - m(4) * (g(1) * (t28 * rSges(4,1) - t27 * rSges(4,2) + t71 * t47 - t44) + g(2) * (t26 * rSges(4,1) - t25 * rSges(4,2) + t71 * t46 - t43) + g(3) * (t42 * rSges(4,1) - t41 * rSges(4,2) + rSges(4,3) * t75 + t109)) - m(5) * (g(1) * (t14 * rSges(5,1) - t13 * rSges(5,2) + t123 * t27 + t78) + g(2) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t123 * t25 + t79) + g(3) * (rSges(5,1) * t30 - rSges(5,2) * t29 + t123 * t41 + t102)) - m(6) * (g(1) * (t14 * pkin(4) + t27 * pkin(10) + (t14 * t63 + t27 * t60) * rSges(6,1) + (-t14 * t60 + t27 * t63) * rSges(6,2) + t122 * t13 + t78) + g(2) * (t12 * pkin(4) + t25 * pkin(10) + (t12 * t63 + t25 * t60) * rSges(6,1) + (-t12 * t60 + t25 * t63) * rSges(6,2) + t122 * t11 + t79) + g(3) * (t30 * pkin(4) + t41 * pkin(10) + (t30 * t63 + t41 * t60) * rSges(6,1) + (-t30 * t60 + t41 * t63) * rSges(6,2) + t122 * t29 + t102)) - m(7) * (g(1) * (t110 * t13 + t80 * t14 + t72 * t27 + t78) + g(2) * (t110 * t11 + t80 * t12 + t72 * t25 + t79) + g(3) * (t110 * t29 + t80 * t30 + t72 * t41 + t102)) -m(4) * (g(1) * (-rSges(4,1) * t19 - rSges(4,2) * t20) + g(2) * (-rSges(4,1) * t17 - rSges(4,2) * t18) + g(3) * (-rSges(4,1) * t32 - rSges(4,2) * t33)) - m(5) * (g(1) * (t123 * t20 + t90 * t19 - t16) + g(2) * (t123 * t18 + t90 * t17 - t15) + g(3) * (t123 * t33 + t90 * t32 - t31)) + (-g(1) * (-t19 * t111 + pkin(5) * t117 + (-t19 * t114 + t20 * t57) * rSges(7,1) + (t19 * t115 + t20 * t58) * rSges(7,2) + t100) - g(2) * (-t17 * t111 + pkin(5) * t118 + (-t17 * t114 + t18 * t57) * rSges(7,1) + (t17 * t115 + t18 * t58) * rSges(7,2) + t101) - g(3) * (-t32 * t111 + pkin(5) * t116 + (-t32 * t114 + t33 * t57) * rSges(7,1) + (t32 * t115 + t33 * t58) * rSges(7,2) + t99) - t110 * t132) * m(7) + (-g(1) * (-t19 * t124 + (-t19 * t112 + t117) * rSges(6,1) + (t19 * t113 + t20 * t63) * rSges(6,2) + t100) - g(2) * (-t17 * t124 + (-t17 * t112 + t118) * rSges(6,1) + (t17 * t113 + t18 * t63) * rSges(6,2) + t101) - g(3) * (-t32 * t124 + (-t32 * t112 + t116) * rSges(6,1) + (t32 * t113 + t33 * t63) * rSges(6,2) + t99) - t122 * t132) * m(6), -m(5) * (g(1) * (rSges(5,1) * t9 - rSges(5,2) * t10) + g(2) * (rSges(5,1) * t7 - rSges(5,2) * t8) + g(3) * (rSges(5,1) * t21 - rSges(5,2) * t22)) - m(6) * (t134 * (t63 * rSges(6,1) - t60 * rSges(6,2) + pkin(4)) + t133 * t122) - m(7) * (t133 * t110 + t134 * t80) -m(6) * (g(1) * (t89 * rSges(6,1) + (-t10 * t63 - t19 * t60) * rSges(6,2)) + g(2) * (t94 * rSges(6,1) + (-t17 * t60 - t63 * t8) * rSges(6,2)) + g(3) * (t88 * rSges(6,1) + (-t22 * t63 - t32 * t60) * rSges(6,2))) - t70 - m(7) * (g(1) * t89 + g(2) * t94 + g(3) * t88) * pkin(5), -t70];
taug  = t1(:);
