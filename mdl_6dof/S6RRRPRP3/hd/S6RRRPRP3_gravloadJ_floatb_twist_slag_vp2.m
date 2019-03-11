% Calculate Gravitation load on the joints for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:18
% EndTime: 2019-03-09 16:39:21
% DurationCPUTime: 1.00s
% Computational Cost: add. (585->123), mult. (620->142), div. (0->0), fcn. (559->10), ass. (0->69)
t128 = -mrSges(6,3) - mrSges(7,2);
t132 = mrSges(5,3) - t128;
t130 = -mrSges(6,1) - mrSges(7,1);
t112 = m(7) * pkin(5) - t130;
t129 = mrSges(6,2) - mrSges(7,3);
t44 = pkin(10) + qJ(5);
t39 = sin(t44);
t45 = qJ(2) + qJ(3);
t41 = sin(t45);
t100 = t39 * t41;
t46 = sin(pkin(10));
t91 = t46 * mrSges(5,2);
t127 = -mrSges(6,2) * t100 - t41 * t91;
t47 = cos(pkin(10));
t90 = t47 * mrSges(5,1);
t125 = t90 - t91;
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t119 = g(1) * t52 + g(2) * t50;
t124 = t132 * t41;
t36 = pkin(4) * t47 + pkin(3);
t40 = cos(t44);
t85 = qJ(6) * t39;
t123 = (t90 - m(7) * (-t36 - t85) + t39 * mrSges(7,3) + t112 * t40) * t41;
t122 = -m(6) - m(7);
t42 = cos(t45);
t121 = t42 * mrSges(4,1) - t41 * mrSges(4,2);
t49 = sin(qJ(2));
t109 = pkin(2) * t49;
t48 = -pkin(9) - qJ(4);
t94 = t42 * t48;
t68 = -t36 * t41 - t94;
t118 = -m(7) * (-t94 - t109) - m(6) * (t68 - t109) - m(5) * (-pkin(3) * t41 - t109) + t123;
t117 = -m(6) * t68 + m(7) * t94 + t123;
t92 = t42 * t52;
t116 = t127 * t52 - t132 * t92;
t93 = t42 * t50;
t115 = t127 * t50 - t132 * t93;
t53 = -pkin(8) - pkin(7);
t114 = -m(3) * pkin(7) + m(5) * t53 - t46 * mrSges(5,1) - t47 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t98 = t40 * t42;
t113 = t130 * t98 - t121 - t124 + (t129 * t39 - t125) * t42;
t111 = m(7) * qJ(6) - t129;
t51 = cos(qJ(2));
t72 = t51 * mrSges(3,1) - t49 * mrSges(3,2);
t110 = -(m(5) * pkin(3) + t125) * t42 - mrSges(2,1) - m(3) * pkin(1) - t72 - t121;
t108 = pkin(4) * t46;
t43 = t51 * pkin(2);
t96 = t41 * t50;
t95 = t41 * t52;
t11 = t42 * t36;
t89 = t50 * t40;
t88 = t52 * t39;
t31 = t41 * qJ(4);
t87 = t42 * pkin(3) + t31;
t86 = qJ(4) * t42;
t76 = -t41 * t48 + t11;
t38 = t43 + pkin(1);
t20 = t52 * t38;
t75 = -t50 * t53 + t20;
t69 = mrSges(4,1) * t41 + mrSges(4,2) * t42;
t65 = pkin(5) * t98 + t42 * t85 + t76;
t19 = t52 * t86;
t18 = t50 * t86;
t4 = t39 * t50 + t40 * t92;
t3 = t42 * t88 - t89;
t2 = t42 * t89 - t88;
t1 = t39 * t93 + t40 * t52;
t5 = [(-m(4) * t75 - m(5) * t20 + t128 * t95 + t122 * (t50 * t108 + t36 * t92 - t48 * t95 + t75) - t112 * t4 - t111 * t3 + (-(m(5) * qJ(4) + mrSges(5,3)) * t41 + t110) * t52 + t114 * t50) * g(2) + (t122 * t48 * t96 + t112 * t2 + t111 * t1 + (m(4) * t38 - m(5) * (-t38 - t31) + t122 * (-t38 - t11) - t110 + t124) * t50 + (t122 * (t108 - t53) + m(4) * t53 + t114) * t52) * g(1) (-m(5) * t18 + t118 * t50 + t115) * g(2) + (-m(5) * t19 + t118 * t52 + t116) * g(1) + (-t72 - m(4) * t43 - m(5) * (t43 + t87) - m(6) * (t43 + t76) - m(7) * (t43 + t65) + t113) * g(3) + t119 * (m(4) * t109 + mrSges(3,1) * t49 + mrSges(3,2) * t51 + t69) t119 * t69 + (-m(5) * (-pkin(3) * t96 + t18) + t117 * t50 + t115) * g(2) + (-m(5) * (-pkin(3) * t95 + t19) + t117 * t52 + t116) * g(1) + (-m(5) * t87 - m(6) * t76 - m(7) * t65 + t113) * g(3) (t42 * g(3) - t119 * t41) * (m(5) - t122) (-t111 * t40 + t112 * t39) * g(3) * t41 + (t1 * t112 - t111 * t2) * g(2) + (-t111 * t4 + t112 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t100) * m(7)];
taug  = t5(:);
