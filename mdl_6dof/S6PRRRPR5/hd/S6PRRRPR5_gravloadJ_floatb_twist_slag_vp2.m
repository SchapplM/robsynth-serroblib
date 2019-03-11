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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:13
% EndTime: 2019-03-08 23:23:16
% DurationCPUTime: 1.36s
% Computational Cost: add. (1059->140), mult. (2559->218), div. (0->0), fcn. (3195->16), ass. (0->73)
t140 = m(6) + m(7);
t145 = t140 * pkin(4);
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t144 = m(7) * pkin(5) + mrSges(7,1) * t69 - mrSges(7,2) * t66 + mrSges(6,1);
t136 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t129 = cos(qJ(3));
t114 = cos(pkin(6));
t130 = cos(qJ(2));
t64 = sin(pkin(7));
t111 = sin(pkin(6));
t113 = cos(pkin(7));
t91 = t113 * t111;
t135 = t114 * t64 + t130 * t91;
t68 = sin(qJ(3));
t128 = sin(qJ(2));
t96 = t111 * t128;
t34 = t129 * t96 + t135 * t68;
t97 = t130 * t111;
t48 = t113 * t114 - t64 * t97;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t143 = -t34 * t67 + t48 * t70;
t112 = cos(pkin(12));
t110 = sin(pkin(12));
t92 = t114 * t110;
t76 = t112 * t128 + t130 * t92;
t89 = t111 * t110;
t138 = t113 * t76 - t64 * t89;
t50 = t112 * t130 - t128 * t92;
t18 = t129 * t50 - t138 * t68;
t36 = t113 * t89 + t64 * t76;
t142 = -t18 * t67 + t36 * t70;
t93 = t114 * t112;
t75 = t110 * t128 - t130 * t93;
t90 = t112 * t111;
t139 = t113 * t75 + t64 * t90;
t49 = t110 * t130 + t128 * t93;
t16 = t49 * t129 - t139 * t68;
t35 = -t113 * t90 + t64 * t75;
t141 = -t16 * t67 + t35 * t70;
t79 = -m(5) * pkin(3) - t70 * mrSges(5,1) + t67 * mrSges(5,2) - mrSges(4,1);
t109 = -m(4) - m(5) - t140;
t137 = -pkin(2) * t109 + mrSges(3,1);
t63 = qJ(4) + pkin(13);
t61 = sin(t63);
t62 = cos(t63);
t134 = -t136 * t61 + t144 * t62 - t79;
t77 = -m(5) * pkin(10) - mrSges(7,1) * t66 - mrSges(7,2) * t69 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t121 = t61 * t64;
t120 = t62 * t64;
t88 = t64 * t96;
t115 = pkin(2) * t97 + pkin(9) * t88;
t106 = t68 * t113;
t98 = t113 * t129;
t85 = t67 * t88;
t81 = t128 * t91;
t71 = mrSges(3,2) + (-t70 * mrSges(5,2) - mrSges(4,3) + (-mrSges(5,1) - t145) * t67 + t109 * pkin(9)) * t64;
t65 = -qJ(5) - pkin(10);
t60 = pkin(4) * t70 + pkin(3);
t44 = t129 * t97 - t68 * t81;
t43 = t129 * t81 + t68 * t97;
t33 = -t129 * t135 + t68 * t96;
t26 = -t106 * t50 - t129 * t76;
t25 = t50 * t98 - t68 * t76;
t24 = -t106 * t49 - t129 * t75;
t23 = t49 * t98 - t68 * t75;
t17 = t129 * t138 + t50 * t68;
t15 = t129 * t139 + t49 * t68;
t12 = t34 * t62 + t48 * t61;
t4 = t18 * t62 + t36 * t61;
t2 = t16 * t62 + t35 * t61;
t1 = [(-m(2) - m(3) + t109) * g(3) (-mrSges(3,1) * t97 + mrSges(3,2) * t96 - m(4) * t115 - t44 * mrSges(4,1) - mrSges(4,3) * t88 - m(5) * (pkin(3) * t44 + t115) - (t44 * t70 + t85) * mrSges(5,1) - (-t44 * t67 + t70 * t88) * mrSges(5,2) + t136 * (t44 * t61 - t62 * t88) - t144 * (t44 * t62 + t61 * t88) + t77 * t43 - t140 * (pkin(4) * t85 - t43 * t65 + t44 * t60 + t115)) * g(3) + (t136 * (-t120 * t49 + t24 * t61) - t144 * (t121 * t49 + t24 * t62) + t79 * t24 + t77 * t23 + t71 * t49 + t137 * t75 - t140 * (-t23 * t65 + t24 * t60)) * g(2) + (t136 * (-t120 * t50 + t26 * t61) - t144 * (t121 * t50 + t26 * t62) + t79 * t26 + t77 * t25 + t71 * t50 + t137 * t76 - t140 * (-t25 * t65 + t26 * t60)) * g(1) (-t140 * (-t33 * t60 - t34 * t65) + t77 * t34 + t134 * t33) * g(3) + (-t140 * (-t15 * t60 - t16 * t65) + t77 * t16 + t134 * t15) * g(2) + (-t140 * (-t17 * t60 - t18 * t65) + t77 * t18 + t134 * t17) * g(1) (-t143 * mrSges(5,1) - (-t34 * t70 - t48 * t67) * mrSges(5,2) + t136 * t12 - t144 * (-t34 * t61 + t48 * t62)) * g(3) + (-t141 * mrSges(5,1) - (-t16 * t70 - t35 * t67) * mrSges(5,2) + t136 * t2 - t144 * (-t16 * t61 + t35 * t62)) * g(2) + (-t142 * mrSges(5,1) - (-t18 * t70 - t36 * t67) * mrSges(5,2) + t136 * t4 - t144 * (-t18 * t61 + t36 * t62)) * g(1) + (-g(1) * t142 - g(2) * t141 - g(3) * t143) * t145, t140 * (-g(1) * t17 - g(2) * t15 - g(3) * t33) -g(1) * ((t17 * t69 - t4 * t66) * mrSges(7,1) + (-t17 * t66 - t4 * t69) * mrSges(7,2)) - g(2) * ((t15 * t69 - t2 * t66) * mrSges(7,1) + (-t15 * t66 - t2 * t69) * mrSges(7,2)) - g(3) * ((-t12 * t66 + t33 * t69) * mrSges(7,1) + (-t12 * t69 - t33 * t66) * mrSges(7,2))];
taug  = t1(:);
