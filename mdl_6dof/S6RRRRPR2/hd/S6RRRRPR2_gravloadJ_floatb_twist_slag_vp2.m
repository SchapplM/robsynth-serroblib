% Calculate Gravitation load on the joints for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:24
% EndTime: 2019-03-09 21:57:26
% DurationCPUTime: 0.89s
% Computational Cost: add. (652->131), mult. (553->136), div. (0->0), fcn. (471->12), ass. (0->76)
t40 = pkin(11) + qJ(6);
t34 = sin(t40);
t101 = t34 * mrSges(7,2);
t42 = qJ(2) + qJ(3);
t38 = qJ(4) + t42;
t31 = sin(t38);
t32 = cos(t38);
t43 = sin(pkin(11));
t94 = t43 * mrSges(6,2);
t138 = t32 * (-mrSges(6,3) - mrSges(7,3)) + (-t101 - t94) * t31;
t44 = cos(pkin(11));
t106 = mrSges(6,1) * t44;
t35 = cos(t40);
t98 = t35 * mrSges(7,1);
t137 = (t106 + t98) * t31;
t29 = pkin(5) * t44 + pkin(4);
t45 = -pkin(10) - qJ(5);
t67 = -t29 * t31 - t32 * t45;
t134 = -m(7) * t67 + t137;
t122 = t32 * t29 - t31 * t45;
t133 = m(7) * t122;
t36 = sin(t42);
t37 = cos(t42);
t69 = mrSges(5,1) * t31 + mrSges(5,2) * t32;
t132 = mrSges(4,1) * t36 + mrSges(4,2) * t37 + t69;
t131 = t106 - t94;
t49 = cos(qJ(1));
t130 = t138 * t49;
t47 = sin(qJ(1));
t129 = t138 * t47;
t128 = -t32 * mrSges(5,1) + (mrSges(5,2) - mrSges(7,3)) * t31;
t120 = g(1) * t49 + g(2) * t47;
t80 = t37 * mrSges(4,1) - t36 * mrSges(4,2);
t126 = t134 * t49 + t130;
t125 = t134 * t47 + t129;
t112 = pkin(4) * t31;
t113 = pkin(3) * t36;
t124 = -m(7) * (t67 - t113) - m(6) * (-t112 - t113) + t137;
t25 = t31 * mrSges(6,3);
t119 = -t25 + t128 + (t101 - t98 - t131) * t32;
t117 = -t80 + t119;
t50 = -pkin(8) - pkin(7);
t41 = -pkin(9) + t50;
t116 = -m(3) * pkin(7) + m(4) * t50 + m(6) * t41 - t43 * mrSges(6,1) - t44 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t48 = cos(qJ(2));
t39 = t48 * pkin(2);
t46 = sin(qJ(2));
t73 = t48 * mrSges(3,1) - t46 * mrSges(3,2);
t115 = -(m(6) * pkin(4) + t131) * t32 - mrSges(2,1) - m(4) * (t39 + pkin(1)) - t80 - m(3) * pkin(1) - t73 + t128;
t114 = pkin(2) * t46;
t30 = pkin(3) * t37;
t111 = pkin(5) * t43;
t100 = t34 * t47;
t99 = t34 * t49;
t97 = t35 * t47;
t96 = t35 * t49;
t23 = t31 * qJ(5);
t91 = t32 * pkin(4) + t23;
t90 = t30 + t39;
t89 = qJ(5) * t32;
t84 = t30 + t91;
t17 = t47 * t89;
t77 = -t112 * t47 + t17;
t18 = t49 * t89;
t76 = -t112 * t49 + t18;
t75 = t30 + t122;
t15 = -t113 - t114;
t14 = pkin(1) + t90;
t9 = t49 * t15;
t8 = t47 * t15;
t5 = t49 * t14;
t4 = t32 * t96 + t100;
t3 = -t32 * t99 + t97;
t2 = -t32 * t97 + t99;
t1 = t100 * t32 + t96;
t6 = [(-m(6) * t5 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + (-m(5) - m(7)) * (-t41 * t47 + t5) + (-m(7) * t111 + t116) * t47 + (-(m(6) * qJ(5) + mrSges(6,3)) * t31 - t133 + t115) * t49) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + (m(5) * t41 - m(7) * (-t41 + t111) + t116) * t49 + (m(5) * t14 - m(6) * (-t14 - t23) + t25 - m(7) * (-t14 - t122) - t115) * t47) * g(1) (-m(6) * (t77 + t8) - m(7) * t8 + t125) * g(2) + (-m(6) * (t76 + t9) - m(7) * t9 + t126) * g(1) + (-t73 - m(4) * t39 - m(5) * t90 - m(6) * (t39 + t84) - m(7) * (t39 + t75) + t117) * g(3) + t120 * (m(4) * t114 - m(5) * t15 + mrSges(3,1) * t46 + mrSges(3,2) * t48 + t132) (-m(6) * t17 + t124 * t47 + t129) * g(2) + (-m(6) * t18 + t124 * t49 + t130) * g(1) + (-m(5) * t30 - m(6) * t84 - m(7) * t75 + t117) * g(3) + t120 * (m(5) * t113 + t132) t120 * t69 + (-m(6) * t77 + t125) * g(2) + (-m(6) * t76 + t126) * g(1) + (-m(6) * t91 + t119 - t133) * g(3) (g(3) * t32 - t120 * t31) * (m(6) + m(7)) -g(1) * (mrSges(7,1) * t3 - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - g(3) * (-mrSges(7,1) * t34 - mrSges(7,2) * t35) * t31];
taug  = t6(:);
