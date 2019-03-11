% Calculate Gravitation load on the joints for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:41
% EndTime: 2019-03-10 01:00:44
% DurationCPUTime: 1.02s
% Computational Cost: add. (699->126), mult. (645->133), div. (0->0), fcn. (574->10), ass. (0->74)
t134 = -mrSges(6,1) - mrSges(7,1);
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t133 = -t43 * mrSges(7,3) + t134 * t46;
t132 = -mrSges(6,3) - mrSges(7,2);
t125 = pkin(5) * t46 + qJ(6) * t43;
t42 = qJ(2) + qJ(3);
t39 = qJ(4) + t42;
t34 = sin(t39);
t131 = (-m(7) * (-pkin(4) - t125) - t133) * t34;
t130 = t132 * t34;
t37 = sin(t42);
t38 = cos(t42);
t35 = cos(t39);
t65 = mrSges(5,1) * t34 + mrSges(5,2) * t35;
t129 = mrSges(4,1) * t37 + mrSges(4,2) * t38 + t65;
t48 = cos(qJ(1));
t93 = t43 * mrSges(6,2);
t82 = t34 * t93;
t96 = t35 * t48;
t127 = t132 * t96 - t48 * t82;
t45 = sin(qJ(1));
t97 = t35 * t45;
t126 = t132 * t97 - t45 * t82;
t124 = g(1) * t48 + g(2) * t45;
t123 = m(6) + m(7);
t122 = t35 * mrSges(5,1) - t34 * mrSges(5,2);
t85 = t35 * pkin(4) + t34 * pkin(10);
t75 = t38 * mrSges(4,1) - t37 * mrSges(4,2);
t106 = pkin(4) * t34;
t107 = pkin(3) * t37;
t121 = m(7) * t107 - m(6) * (-t106 - t107) + t131;
t120 = t131 * t48 + t127;
t119 = t131 * t45 + t126;
t117 = -t122 + t130 + (t93 + t133) * t35;
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t44 = sin(qJ(2));
t69 = t47 * mrSges(3,1) - t44 * mrSges(3,2);
t115 = -mrSges(2,1) - m(3) * pkin(1) - t69 - m(4) * (t40 + pkin(1)) - t75 - t122;
t49 = -pkin(8) - pkin(7);
t114 = -m(3) * pkin(7) + m(4) * t49 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t113 = -t75 + t117;
t112 = m(7) * pkin(5) - t134;
t111 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t20 = pkin(10) * t97;
t110 = m(7) * t20;
t23 = pkin(10) * t96;
t109 = m(7) * t23;
t108 = pkin(2) * t44;
t31 = pkin(3) * t38;
t102 = g(3) * t34;
t98 = t34 * t48;
t91 = t43 * t45;
t90 = t45 * t46;
t87 = t46 * t48;
t86 = t48 * t43;
t84 = t31 + t40;
t79 = t31 + t85;
t10 = pkin(1) + t84;
t41 = -pkin(9) + t49;
t76 = t48 * t10 - t41 * t45;
t73 = t125 * t35 + t85;
t72 = -t45 * t106 + t20;
t71 = -pkin(4) * t98 + t23;
t70 = t31 + t73;
t11 = -t107 - t108;
t7 = t48 * t11;
t6 = t45 * t11;
t4 = t35 * t87 + t91;
t3 = t35 * t86 - t90;
t2 = t35 * t90 - t86;
t1 = t35 * t91 + t87;
t5 = [(-m(5) * t76 + t132 * t98 - t123 * (pkin(4) * t96 + pkin(10) * t98 + t76) - t112 * t4 - t111 * t3 + t115 * t48 + t114 * t45) * g(2) + (t112 * t2 + t111 * t1 + (m(5) * t10 - t123 * (-t10 - t85) - t115 - t130) * t45 + (t114 + (m(5) + t123) * t41) * t48) * g(1) (-m(6) * (t6 + t72) - m(7) * (t20 + t6) + t119) * g(2) + (-m(6) * (t7 + t71) - m(7) * (t23 + t7) + t120) * g(1) + (-t69 - m(4) * t40 - m(5) * t84 - m(6) * (t40 + t79) - m(7) * (t40 + t70) + t113) * g(3) + t124 * (m(4) * t108 - m(5) * t11 + mrSges(3,1) * t44 + mrSges(3,2) * t47 + t129) (-m(6) * t20 + t121 * t45 - t110 + t126) * g(2) + (-m(6) * t23 + t121 * t48 - t109 + t127) * g(1) + (-m(5) * t31 - m(6) * t79 - m(7) * t70 + t113) * g(3) + t124 * (m(5) * t107 + t129) t124 * t65 + (-m(6) * t72 - t110 + t119) * g(2) + (-m(6) * t71 - t109 + t120) * g(1) + (-m(6) * t85 - m(7) * t73 + t117) * g(3) (-t111 * t46 + t112 * t43) * t102 + (t1 * t112 - t111 * t2) * g(2) + (-t111 * t4 + t112 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - t43 * t102) * m(7)];
taug  = t5(:);
