% Calculate Gravitation load on the joints for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2018-11-23 15:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:25:18
% EndTime: 2018-11-23 15:25:19
% DurationCPUTime: 1.10s
% Computational Cost: add. (1700->135), mult. (2120->183), div. (0->0), fcn. (2182->16), ass. (0->87)
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t134 = pkin(4) * t70 + qJ(5) * t66;
t65 = sin(qJ(6));
t69 = cos(qJ(6));
t120 = -t65 * mrSges(7,1) - t69 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t119 = -m(7) * pkin(5) - t69 * mrSges(7,1) + t65 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t133 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t132 = m(6) + m(7);
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t131 = t71 * mrSges(4,1) - mrSges(4,2) * t67 + mrSges(3,1);
t130 = mrSges(3,2) - mrSges(4,3);
t107 = cos(pkin(11));
t63 = sin(pkin(11));
t72 = cos(qJ(2));
t111 = t63 * t72;
t104 = pkin(6) + qJ(2);
t92 = sin(t104);
t60 = t92 / 0.2e1;
t105 = pkin(6) - qJ(2);
t93 = sin(t105);
t82 = t60 - t93 / 0.2e1;
t43 = t107 * t82 + t111;
t106 = sin(pkin(6));
t83 = t107 * t106;
t22 = -t43 * t67 - t71 * t83;
t128 = t134 * t22;
t95 = t107 * t72;
t46 = -t63 * t82 + t95;
t96 = t63 * t106;
t24 = -t46 * t67 + t71 * t96;
t127 = t134 * t24;
t90 = cos(t104) / 0.2e1;
t94 = cos(t105);
t58 = t90 - t94 / 0.2e1;
t64 = cos(pkin(6));
t48 = t58 * t67 + t64 * t71;
t126 = t134 * t48;
t125 = pkin(4) * t132 - t119;
t124 = t119 * t70 + t120 * t66 - mrSges(4,1);
t123 = mrSges(4,2) - m(7) * (pkin(9) - pkin(10)) + t133;
t122 = m(7) * pkin(10) + t133;
t116 = pkin(3) * t71;
t68 = sin(qJ(2));
t75 = t94 / 0.2e1 + t90;
t42 = -t107 * t75 + t63 * t68;
t114 = t42 * t67;
t45 = t107 * t68 + t63 * t75;
t113 = t45 * t67;
t61 = t93 / 0.2e1;
t57 = t60 + t61;
t112 = t57 * t67;
t110 = t66 * t71;
t109 = t70 * t71;
t81 = t61 - t92 / 0.2e1;
t44 = -t107 * t81 + t111;
t102 = -t42 * pkin(2) + pkin(8) * t44;
t47 = t63 * t81 + t95;
t101 = -t45 * pkin(2) + pkin(8) * t47;
t100 = t57 * pkin(2) - pkin(8) * t58;
t20 = t22 * pkin(3);
t23 = t43 * t71 - t67 * t83;
t99 = pkin(9) * t23 + t20;
t21 = t24 * pkin(3);
t25 = t46 * t71 + t67 * t96;
t98 = pkin(9) * t25 + t21;
t41 = t48 * pkin(3);
t49 = -t58 * t71 + t64 * t67;
t97 = pkin(9) * t49 + t41;
t86 = -pkin(9) * t114 - t42 * t116 + t102;
t85 = -pkin(9) * t113 - t45 * t116 + t101;
t84 = pkin(9) * t112 + t57 * t116 + t100;
t74 = -qJ(5) * t132 + t120;
t28 = t109 * t57 - t58 * t66;
t27 = t110 * t57 + t58 * t70;
t15 = t49 * t70 - t57 * t66;
t14 = t49 * t66 + t57 * t70;
t12 = -t109 * t45 + t47 * t66;
t11 = -t110 * t45 - t47 * t70;
t10 = -t109 * t42 + t44 * t66;
t9 = -t110 * t42 - t44 * t70;
t6 = t25 * t70 + t45 * t66;
t5 = t25 * t66 - t45 * t70;
t4 = t23 * t70 + t42 * t66;
t3 = t23 * t66 - t42 * t70;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t132) * g(3) (-m(4) * t100 - m(5) * t84 - t132 * (t28 * pkin(4) + qJ(5) * t27 + t84) - t130 * t58 - t131 * t57 + t119 * t28 + t120 * t27 + t122 * t112) * g(3) + (-m(4) * t102 - m(5) * t86 - t132 * (t10 * pkin(4) + qJ(5) * t9 + t86) + t130 * t44 + t131 * t42 + t120 * t9 + t119 * t10 - t122 * t114) * g(2) + (-m(4) * t101 - m(5) * t85 - t132 * (t12 * pkin(4) + qJ(5) * t11 + t85) + t130 * t47 + t131 * t45 + t119 * t12 + t120 * t11 - t122 * t113) * g(1) (-m(5) * t97 - m(6) * (t97 + t126) - m(7) * (t41 + t126) + t123 * t49 + t124 * t48) * g(3) + (-m(5) * t99 - m(6) * (t99 + t128) - m(7) * (t20 + t128) + t123 * t23 + t124 * t22) * g(2) + (-m(5) * t98 - m(6) * (t98 + t127) - m(7) * (t21 + t127) + t123 * t25 + t124 * t24) * g(1) (t125 * t14 + t15 * t74) * g(3) + (t125 * t3 + t4 * t74) * g(2) + (t125 * t5 + t6 * t74) * g(1), t132 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) -g(1) * ((t5 * t69 - t6 * t65) * mrSges(7,1) + (-t5 * t65 - t6 * t69) * mrSges(7,2)) - g(2) * ((t3 * t69 - t4 * t65) * mrSges(7,1) + (-t3 * t65 - t4 * t69) * mrSges(7,2)) - g(3) * ((t14 * t69 - t15 * t65) * mrSges(7,1) + (-t14 * t65 - t15 * t69) * mrSges(7,2))];
taug  = t1(:);
