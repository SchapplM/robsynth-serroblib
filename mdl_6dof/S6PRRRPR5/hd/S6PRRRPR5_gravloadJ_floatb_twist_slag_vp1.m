% Calculate Gravitation load on the joints for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:03
% EndTime: 2019-03-08 23:23:07
% DurationCPUTime: 1.37s
% Computational Cost: add. (1060->217), mult. (2541->337), div. (0->0), fcn. (3195->16), ass. (0->96)
t125 = cos(qJ(3));
t107 = cos(pkin(6));
t126 = cos(qJ(2));
t65 = sin(pkin(7));
t104 = sin(pkin(6));
t106 = cos(pkin(7));
t87 = t106 * t104;
t141 = t107 * t65 + t126 * t87;
t69 = sin(qJ(3));
t124 = sin(qJ(2));
t91 = t104 * t124;
t35 = t125 * t91 + t141 * t69;
t92 = t126 * t104;
t49 = t107 * t106 - t65 * t92;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t146 = -t35 * t68 + t49 * t71;
t105 = cos(pkin(12));
t103 = sin(pkin(12));
t88 = t107 * t103;
t75 = t105 * t124 + t126 * t88;
t85 = t104 * t103;
t142 = t75 * t106 - t65 * t85;
t51 = t105 * t126 - t124 * t88;
t19 = t51 * t125 - t142 * t69;
t37 = t106 * t85 + t75 * t65;
t145 = -t19 * t68 + t37 * t71;
t89 = t107 * t105;
t74 = t103 * t124 - t126 * t89;
t86 = t105 * t104;
t143 = t74 * t106 + t65 * t86;
t50 = t103 * t126 + t124 * t89;
t17 = t50 * t125 - t143 * t69;
t36 = -t106 * t86 + t74 * t65;
t144 = -t17 * t68 + t36 * t71;
t16 = t143 * t125 + t50 * t69;
t18 = t142 * t125 + t51 * t69;
t34 = -t141 * t125 + t69 * t91;
t138 = g(1) * t18 + g(2) * t16 + g(3) * t34;
t128 = rSges(7,3) + pkin(11);
t140 = g(1) * t51 + g(2) * t50;
t137 = -m(6) - m(7);
t64 = qJ(4) + pkin(13);
t63 = cos(t64);
t131 = t63 * pkin(5);
t130 = t65 * pkin(9);
t127 = pkin(10) + rSges(5,3);
t62 = sin(t64);
t117 = t62 * t65;
t116 = t63 * t65;
t67 = sin(qJ(6));
t115 = t63 * t67;
t70 = cos(qJ(6));
t114 = t63 * t70;
t113 = t65 * t68;
t112 = t65 * t71;
t61 = pkin(4) * t71 + pkin(3);
t66 = -qJ(5) - pkin(10);
t111 = -t16 * t61 - t17 * t66;
t110 = -t18 * t61 - t19 * t66;
t109 = -t34 * t61 - t35 * t66;
t82 = t65 * t91;
t108 = pkin(2) * t92 + pkin(9) * t82;
t93 = t106 * t125;
t24 = t50 * t93 - t74 * t69;
t99 = t69 * t106;
t25 = -t74 * t125 - t50 * t99;
t47 = t74 * pkin(2);
t102 = -t24 * t66 + t25 * t61 - t47;
t26 = t51 * t93 - t75 * t69;
t27 = -t75 * t125 - t51 * t99;
t48 = t75 * pkin(2);
t101 = -t26 * t66 + t27 * t61 - t48;
t97 = t144 * pkin(4);
t96 = t145 * pkin(4);
t95 = t146 * pkin(4);
t77 = t124 * t87;
t44 = t125 * t77 + t69 * t92;
t45 = t125 * t92 - t69 * t77;
t79 = t68 * t82;
t94 = pkin(4) * t79 - t44 * t66 + t45 * t61 + t108;
t90 = -rSges(6,1) * t63 + rSges(6,2) * t62;
t76 = t140 * (pkin(4) * t68 + pkin(9)) * t65;
t29 = t45 * t63 + t62 * t82;
t28 = t45 * t62 - t63 * t82;
t13 = t35 * t63 + t49 * t62;
t12 = -t35 * t62 + t49 * t63;
t9 = t51 * t117 + t27 * t63;
t8 = -t51 * t116 + t27 * t62;
t7 = t50 * t117 + t25 * t63;
t6 = -t50 * t116 + t25 * t62;
t5 = t19 * t63 + t37 * t62;
t4 = -t19 * t62 + t37 * t63;
t3 = t17 * t63 + t36 * t62;
t2 = -t17 * t62 + t36 * t63;
t1 = [(-m(2) - m(3) - m(4) - m(5) + t137) * g(3), -m(3) * (g(1) * (-t75 * rSges(3,1) - t51 * rSges(3,2)) + g(2) * (-t74 * rSges(3,1) - t50 * rSges(3,2)) + g(3) * (rSges(3,1) * t92 - rSges(3,2) * t91)) - m(4) * (g(1) * (rSges(4,1) * t27 - rSges(4,2) * t26 - t48) + g(2) * (rSges(4,1) * t25 - rSges(4,2) * t24 - t47) + g(3) * (t45 * rSges(4,1) - t44 * rSges(4,2) + t108) + (g(3) * rSges(4,3) * t91 + t140 * (rSges(4,3) + pkin(9))) * t65) - m(5) * (g(1) * (t27 * pkin(3) - t48 + t51 * t130 + (t51 * t113 + t27 * t71) * rSges(5,1) + (t51 * t112 - t27 * t68) * rSges(5,2) + t127 * t26) + g(2) * (t25 * pkin(3) - t47 + t50 * t130 + (t50 * t113 + t25 * t71) * rSges(5,1) + (t50 * t112 - t25 * t68) * rSges(5,2) + t127 * t24) + g(3) * (t45 * pkin(3) + (t45 * t71 + t79) * rSges(5,1) + (-t45 * t68 + t71 * t82) * rSges(5,2) + t127 * t44 + t108)) - m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t8 + rSges(6,3) * t26 + t101) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + rSges(6,3) * t24 + t102) + g(3) * (rSges(6,1) * t29 - t28 * rSges(6,2) + rSges(6,3) * t44 + t94) + t76) - m(7) * (g(1) * (t9 * pkin(5) + (t26 * t67 + t70 * t9) * rSges(7,1) + (t26 * t70 - t67 * t9) * rSges(7,2) + t101 + t128 * t8) + g(2) * (t7 * pkin(5) + (t24 * t67 + t7 * t70) * rSges(7,1) + (t24 * t70 - t67 * t7) * rSges(7,2) + t102 + t128 * t6) + g(3) * (t29 * pkin(5) + (t29 * t70 + t44 * t67) * rSges(7,1) + (-t29 * t67 + t44 * t70) * rSges(7,2) + t128 * t28 + t94) + t76) -m(4) * (g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(3) * (-rSges(4,1) * t34 - rSges(4,2) * t35)) - m(5) * (t138 * (-t71 * rSges(5,1) + t68 * rSges(5,2) - pkin(3)) + (g(1) * t19 + g(2) * t17 + g(3) * t35) * t127) - m(6) * (g(1) * (rSges(6,3) * t19 + t90 * t18 + t110) + g(2) * (rSges(6,3) * t17 + t90 * t16 + t111) + g(3) * (rSges(6,3) * t35 + t90 * t34 + t109)) + (-g(1) * (-t18 * t131 + (-t18 * t114 + t19 * t67) * rSges(7,1) + (t18 * t115 + t19 * t70) * rSges(7,2) + t110) - g(2) * (-t16 * t131 + (-t16 * t114 + t17 * t67) * rSges(7,1) + (t16 * t115 + t17 * t70) * rSges(7,2) + t111) - g(3) * (-t34 * t131 + (-t34 * t114 + t35 * t67) * rSges(7,1) + (t34 * t115 + t35 * t70) * rSges(7,2) + t109) + t138 * t62 * t128) * m(7), -m(5) * (g(1) * (t145 * rSges(5,1) + (-t19 * t71 - t37 * t68) * rSges(5,2)) + g(2) * (t144 * rSges(5,1) + (-t17 * t71 - t36 * t68) * rSges(5,2)) + g(3) * (t146 * rSges(5,1) + (-t35 * t71 - t49 * t68) * rSges(5,2))) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t5 + t96) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t3 + t97) + g(3) * (rSges(6,1) * t12 - rSges(6,2) * t13 + t95)) + (-g(1) * (t128 * t5 + t96) - g(2) * (t128 * t3 + t97) - g(3) * (t128 * t13 + t95) - (g(1) * t4 + g(2) * t2 + g(3) * t12) * (t70 * rSges(7,1) - t67 * rSges(7,2) + pkin(5))) * m(7), t137 * t138, -m(7) * (g(1) * ((t18 * t70 - t5 * t67) * rSges(7,1) + (-t18 * t67 - t5 * t70) * rSges(7,2)) + g(2) * ((t16 * t70 - t3 * t67) * rSges(7,1) + (-t16 * t67 - t3 * t70) * rSges(7,2)) + g(3) * ((-t13 * t67 + t34 * t70) * rSges(7,1) + (-t13 * t70 - t34 * t67) * rSges(7,2)))];
taug  = t1(:);
