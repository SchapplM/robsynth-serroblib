% Calculate Gravitation load on the joints for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:51:59
% EndTime: 2018-11-23 17:52:00
% DurationCPUTime: 1.02s
% Computational Cost: add. (618->140), mult. (761->159), div. (0->0), fcn. (761->10), ass. (0->72)
t123 = -mrSges(4,1) - mrSges(5,1);
t112 = m(6) + m(7);
t47 = qJ(2) + qJ(3);
t44 = sin(t47);
t45 = cos(t47);
t122 = t123 * t45 + (mrSges(4,2) - mrSges(5,3)) * t44;
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t121 = g(1) * t53 + g(2) * t50;
t48 = sin(qJ(6));
t51 = cos(qJ(6));
t119 = mrSges(7,1) * t51 - mrSges(7,2) * t48;
t99 = sin(qJ(5));
t81 = t53 * t99;
t100 = cos(qJ(5));
t82 = t53 * t100;
t15 = -t44 * t81 - t45 * t82;
t16 = -t44 * t82 + t45 * t81;
t101 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t94 = t45 * t53;
t118 = -mrSges(5,3) * t94 - t15 * mrSges(7,3) - t119 * t16 - t101;
t23 = t44 * t100 - t45 * t99;
t13 = t23 * t50;
t22 = t45 * t100 + t44 * t99;
t14 = t22 * t50;
t106 = -t13 * mrSges(6,1) + t14 * mrSges(6,2);
t117 = -t50 * t45 * mrSges(5,3) + t14 * mrSges(7,3) + t119 * t13 - t106;
t85 = m(7) * pkin(10) + mrSges(7,3);
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t73 = t52 * mrSges(3,1) - t49 * mrSges(3,2);
t116 = -m(3) * pkin(1) - mrSges(2,1) + t122 - t73;
t93 = t22 * mrSges(6,1) + t23 * mrSges(6,2);
t115 = t23 * mrSges(7,3) - t119 * t22 + t122 - t93;
t114 = -mrSges(6,2) + t85;
t54 = -pkin(8) - pkin(7);
t113 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3) - t112 * (-pkin(9) - t54);
t90 = qJ(4) * t45;
t28 = t50 * t90;
t105 = pkin(2) * t49;
t107 = pkin(3) + pkin(4);
t64 = -t107 * t44 - t105;
t111 = t64 * t50 + t28;
t37 = t44 * qJ(4);
t110 = pkin(3) * t94 + t53 * t37;
t30 = t53 * t90;
t109 = t64 * t53 + t30;
t42 = t45 * pkin(3);
t46 = t52 * pkin(2);
t97 = mrSges(4,2) * t45;
t91 = t42 + t37;
t89 = m(5) + t112;
t41 = t45 * pkin(4);
t88 = t41 + t91;
t87 = t46 + t91;
t80 = -t13 * pkin(5) - pkin(10) * t14;
t79 = t16 * pkin(5) + pkin(10) * t15;
t43 = t46 + pkin(1);
t31 = t53 * t43;
t77 = -t50 * t54 + t31;
t76 = -t43 - t37;
t75 = pkin(4) * t94 + t110 + t31;
t69 = -t14 * t51 - t48 * t53;
t68 = t14 * t48 - t51 * t53;
t63 = t22 * pkin(5) - pkin(10) * t23 + t88;
t61 = -m(7) * pkin(5) - t119;
t57 = (-t107 * t45 + t76) * t50;
t56 = m(5) * (-pkin(3) * t44 - t105) - t44 * mrSges(5,1);
t55 = t97 + (m(5) * pkin(3) + t112 * t107 - t123) * t44;
t2 = -t15 * t51 - t48 * t50;
t1 = t15 * t48 - t50 * t51;
t3 = [(-m(4) * t77 - m(5) * (t77 + t110) - m(6) * t75 + t15 * mrSges(6,1) - m(7) * (-t15 * pkin(5) + t75) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t114 * t16 + t116 * t53 + t113 * t50) * g(2) + (-m(6) * t57 + t14 * mrSges(6,1) - t69 * mrSges(7,1) - t68 * mrSges(7,2) - (-t14 * pkin(5) + t57) * m(7) - t114 * t13 + (m(4) * t43 - m(5) * (t76 - t42) - t116) * t50 + ((m(4) + m(5)) * t54 + t113) * t53) * g(1) (-m(5) * t28 - t56 * t50 - t111 * m(6) - (t80 + t111) * m(7) + t117) * g(2) + (-m(5) * t30 - t56 * t53 - t109 * m(6) - (t79 + t109) * m(7) + t118) * g(1) + (-t73 - m(4) * t46 - m(5) * t87 - m(6) * (t41 + t87) - m(7) * (t46 + t63) + t115) * g(3) + t121 * (m(4) * t105 + mrSges(3,1) * t49 + mrSges(4,1) * t44 + mrSges(3,2) * t52 + t97) (-m(5) * t91 - m(6) * t88 - m(7) * t63 + t115) * g(3) + (-m(7) * t80 - t89 * t28 + t55 * t50 + t117) * g(2) + (-m(7) * t79 - t89 * t30 + t55 * t53 + t118) * g(1) (t45 * g(3) - t121 * t44) * t89 (-t61 * t22 - t85 * t23 + t93) * g(3) + (t61 * t13 - t85 * t14 + t106) * g(2) + (t85 * t15 - t61 * t16 + t101) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t68 * mrSges(7,1) + t69 * mrSges(7,2)) - g(3) * (-t48 * mrSges(7,1) - t51 * mrSges(7,2)) * t23];
taug  = t3(:);
