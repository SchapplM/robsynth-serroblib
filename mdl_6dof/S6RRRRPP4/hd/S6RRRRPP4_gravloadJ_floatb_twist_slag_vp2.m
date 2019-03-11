% Calculate Gravitation load on the joints for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:30
% EndTime: 2019-03-09 20:58:33
% DurationCPUTime: 1.14s
% Computational Cost: add. (638->149), mult. (717->171), div. (0->0), fcn. (681->10), ass. (0->79)
t140 = mrSges(6,1) + mrSges(7,1);
t139 = mrSges(6,2) - mrSges(7,3);
t132 = m(6) + m(7);
t141 = mrSges(4,3) + mrSges(5,3);
t138 = -mrSges(6,3) - mrSges(7,2);
t62 = -pkin(9) - pkin(8);
t137 = -m(4) * pkin(8) + m(5) * t62 - t141;
t55 = qJ(3) + qJ(4);
t47 = pkin(10) + t55;
t43 = sin(t47);
t44 = cos(t47);
t59 = cos(qJ(3));
t51 = t59 * pkin(3);
t46 = t51 + pkin(2);
t48 = sin(t55);
t49 = cos(t55);
t56 = sin(qJ(3));
t136 = m(4) * pkin(2) + m(5) * t46 + t59 * mrSges(4,1) + t49 * mrSges(5,1) - t56 * mrSges(4,2) - t48 * mrSges(5,2) - t139 * t43 + t140 * t44;
t135 = m(7) * qJ(6) + mrSges(7,3);
t58 = sin(qJ(1));
t60 = cos(qJ(2));
t61 = cos(qJ(1));
t97 = t60 * t61;
t21 = -t48 * t97 + t49 * t58;
t57 = sin(qJ(2));
t134 = t138 * t57;
t133 = g(1) * t61 + g(2) * t58;
t120 = m(5) * pkin(3);
t130 = mrSges(2,2) - mrSges(3,3);
t129 = mrSges(5,1) * t48 + mrSges(6,1) * t43 + mrSges(5,2) * t49 + mrSges(6,2) * t44;
t128 = mrSges(4,1) + t120;
t127 = -m(3) - m(4) - m(5);
t96 = t61 * t43;
t13 = -t58 * t44 + t60 * t96;
t14 = t43 * t58 + t44 * t97;
t22 = t48 * t58 + t49 * t97;
t126 = -t21 * mrSges(5,1) + t22 * mrSges(5,2) + t13 * t140 + t139 * t14;
t98 = t58 * t60;
t11 = t43 * t98 + t44 * t61;
t12 = t44 * t98 - t96;
t19 = t48 * t98 + t49 * t61;
t20 = t48 * t61 - t49 * t98;
t125 = t19 * mrSges(5,1) - t20 * mrSges(5,2) + t11 * t140 + t139 * t12;
t79 = t60 * mrSges(3,1) - t57 * mrSges(3,2);
t124 = t141 * t57 + mrSges(2,1) + t79;
t123 = m(7) * pkin(5) + t140;
t122 = -mrSges(6,2) + t135;
t115 = pkin(3) * t56;
t114 = pkin(4) * t48;
t113 = pkin(5) * t43;
t110 = g(3) * t57;
t109 = t43 * mrSges(7,1);
t106 = t56 * t58;
t105 = t56 * t61;
t54 = -qJ(5) + t62;
t100 = t57 * t54;
t99 = t57 * t61;
t36 = pkin(4) * t49 + t51;
t34 = pkin(2) + t36;
t30 = t60 * t34;
t35 = t114 + t115;
t33 = t61 * t35;
t93 = -t60 * t33 + t58 * t36;
t92 = t61 * pkin(1) + t58 * pkin(7);
t87 = t135 * t44 * t57;
t85 = -t11 * pkin(5) + qJ(6) * t12;
t84 = -t35 * t98 - t36 * t61;
t82 = -t13 * pkin(5) + qJ(6) * t14;
t81 = pkin(2) * t60 + pkin(8) * t57;
t73 = pkin(5) * t44 + qJ(6) * t43;
t72 = t60 * t46 - t57 * t62;
t71 = t21 * pkin(4);
t28 = -t56 * t97 + t58 * t59;
t26 = t56 * t98 + t59 * t61;
t68 = t19 * pkin(4);
t52 = t61 * pkin(7);
t29 = t59 * t97 + t106;
t27 = -t59 * t98 + t105;
t1 = [(-t106 * t120 - t29 * mrSges(4,1) - t22 * mrSges(5,1) - t28 * mrSges(4,2) - t21 * mrSges(5,2) + t138 * t99 + t127 * t92 - t132 * (t34 * t97 + t58 * t35 - t54 * t99 + t92) + t130 * t58 - t123 * t14 - t122 * t13 + (-m(4) * t81 - m(5) * t72 - t124) * t61) * g(2) + (-t105 * t120 - t27 * mrSges(4,1) - t20 * mrSges(5,1) - t26 * mrSges(4,2) - t19 * mrSges(5,2) - t132 * (t58 * t100 + t33 + t52) + t130 * t61 + t127 * t52 + t123 * t12 + t122 * t11 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t81) - m(5) * (-pkin(1) - t72) - t132 * (-pkin(1) - t30) + t124 - t134) * t58) * g(1) (-t79 - t132 * (t30 - t100) + t134) * g(3) + ((-m(7) * t73 - t136) * g(3) + t133 * (t132 * t54 + mrSges(3,2) + t137 + t138)) * t60 + (t137 * g(3) + t133 * (mrSges(3,1) + m(6) * t34 - m(7) * (-t34 - t73) + t136)) * t57, -g(3) * ((m(7) * (-t35 - t113) - t109) * t57 + t87) + (m(5) * t115 + m(6) * t35 + mrSges(4,1) * t56 + mrSges(4,2) * t59 + t129) * t110 + (-t27 * mrSges(4,2) - m(6) * t84 - m(7) * (t84 + t85) + t128 * t26 + t125) * g(2) + (t29 * mrSges(4,2) - m(6) * t93 - m(7) * (t82 + t93) - t128 * t28 + t126) * g(1), -g(3) * ((m(7) * (-t113 - t114) - t109) * t57 + t87) + (m(6) * t114 + t129) * t110 + (m(6) * t68 - m(7) * (-t68 + t85) + t125) * g(2) + (-m(6) * t71 - m(7) * (t71 + t82) + t126) * g(1) (t60 * g(3) - t57 * t133) * t132 (-g(1) * t13 - g(2) * t11 - t43 * t110) * m(7)];
taug  = t1(:);
