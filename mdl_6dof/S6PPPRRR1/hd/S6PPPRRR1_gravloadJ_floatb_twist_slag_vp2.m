% Calculate Gravitation load on the joints for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:48:39
% EndTime: 2018-11-23 14:48:40
% DurationCPUTime: 0.76s
% Computational Cost: add. (4268->103), mult. (4321->151), div. (0->0), fcn. (4281->30), ass. (0->90)
t120 = m(6) + m(7);
t44 = sin(qJ(6));
t47 = cos(qJ(6));
t125 = m(7) * pkin(5) + mrSges(7,1) * t47 - mrSges(7,2) * t44 + mrSges(6,1);
t118 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t124 = pkin(4) * t120 - t118 * t45 + t125 * t48 + mrSges(5,1);
t98 = pkin(8) + qJ(4);
t83 = sin(t98) / 0.2e1;
t99 = pkin(8) - qJ(4);
t90 = sin(t99);
t117 = t83 - t90 / 0.2e1;
t108 = cos(pkin(6));
t94 = pkin(7) + pkin(14);
t76 = sin(t94) / 0.2e1;
t95 = pkin(7) - pkin(14);
t86 = sin(t95);
t37 = t76 - t86 / 0.2e1;
t78 = cos(t95) / 0.2e1;
t88 = cos(t94);
t38 = t78 - t88 / 0.2e1;
t42 = cos(pkin(14));
t96 = pkin(6) + pkin(13);
t77 = sin(t96) / 0.2e1;
t97 = pkin(6) - pkin(13);
t87 = sin(t97);
t63 = t77 + t87 / 0.2e1;
t79 = cos(t97) / 0.2e1;
t89 = cos(t96);
t67 = t79 - t89 / 0.2e1;
t30 = t108 * t38 + t63 * t37 + t67 * t42;
t49 = cos(qJ(4));
t100 = sin(pkin(14));
t62 = t76 + t86 / 0.2e1;
t65 = t78 + t88 / 0.2e1;
t52 = t67 * t100 - t108 * t62 - t63 * t65;
t123 = -t117 * t52 + t30 * t49;
t101 = sin(pkin(13));
t102 = sin(pkin(12));
t106 = cos(pkin(12));
t66 = t79 + t89 / 0.2e1;
t57 = -t106 * t101 - t102 * t66;
t105 = cos(pkin(13));
t64 = t77 - t87 / 0.2e1;
t58 = -t102 * t64 + t106 * t105;
t104 = sin(pkin(6));
t80 = t104 * t102;
t27 = t57 * t37 + t38 * t80 + t58 * t42;
t61 = t62 * t104;
t51 = t58 * t100 - t102 * t61 - t57 * t65;
t122 = -t117 * t51 + t27 * t49;
t55 = -t102 * t101 + t106 * t66;
t56 = t102 * t105 + t106 * t64;
t81 = t106 * t104;
t26 = t55 * t37 - t38 * t81 + t56 * t42;
t50 = t56 * t100 + t106 * t61 - t55 * t65;
t121 = -t117 * t50 + t26 * t49;
t114 = m(4) + m(5) + t120;
t113 = m(3) + t114;
t112 = -t44 * mrSges(7,1) - t47 * mrSges(7,2) - t120 * pkin(10) + mrSges(5,2) - mrSges(6,3);
t107 = cos(pkin(7));
t103 = sin(pkin(7));
t92 = cos(t99);
t91 = cos(t98);
t85 = t92 / 0.2e1;
t84 = t91 / 0.2e1;
t75 = t85 - t91 / 0.2e1;
t71 = t85 + t84;
t70 = t84 - t92 / 0.2e1;
t68 = t83 + t90 / 0.2e1;
t59 = t63 * t103 - t108 * t107;
t54 = t57 * t103 - t107 * t80;
t53 = t55 * t103 + t107 * t81;
t46 = sin(qJ(4));
t43 = cos(pkin(8));
t41 = sin(pkin(8));
t25 = t52 * t41 - t59 * t43;
t20 = t51 * t41 - t54 * t43;
t19 = t50 * t41 - t53 * t43;
t17 = -t59 * t75 + t123;
t16 = t30 * t46 + t52 * t71 + t59 * t68;
t13 = -t54 * t75 + t122;
t12 = t27 * t46 + t51 * t71 + t54 * t68;
t10 = -t53 * t75 + t121;
t9 = t26 * t46 + t50 * t71 + t53 * t68;
t6 = t17 * t48 + t25 * t45;
t4 = t13 * t48 + t20 * t45;
t2 = t10 * t48 + t19 * t45;
t1 = [(-m(2) - t113) * g(3) (-g(1) * t80 + g(2) * t81 - g(3) * t108) * t113 (g(1) * t54 + g(2) * t53 + g(3) * t59) * t114 (t112 * (t59 * t70 + t123) + t124 * t16) * g(3) + (t112 * (t53 * t70 + t121) + t124 * t9) * g(2) + (t112 * (t54 * t70 + t122) + t124 * t12) * g(1) (t118 * t6 - t125 * (-t17 * t45 + t25 * t48)) * g(3) + (t118 * t2 - t125 * (-t10 * t45 + t19 * t48)) * g(2) + (t118 * t4 - t125 * (-t13 * t45 + t20 * t48)) * g(1), -g(1) * ((t12 * t47 - t4 * t44) * mrSges(7,1) + (-t12 * t44 - t4 * t47) * mrSges(7,2)) - g(2) * ((-t2 * t44 + t47 * t9) * mrSges(7,1) + (-t2 * t47 - t44 * t9) * mrSges(7,2)) - g(3) * ((t16 * t47 - t44 * t6) * mrSges(7,1) + (-t16 * t44 - t47 * t6) * mrSges(7,2))];
taug  = t1(:);
