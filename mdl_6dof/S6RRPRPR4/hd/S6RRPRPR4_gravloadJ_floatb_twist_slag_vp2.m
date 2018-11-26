% Calculate Gravitation load on the joints for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:01:52
% EndTime: 2018-11-23 17:01:54
% DurationCPUTime: 1.13s
% Computational Cost: add. (1542->153), mult. (1244->199), div. (0->0), fcn. (1146->22), ass. (0->81)
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t135 = m(7) * pkin(5) + mrSges(7,1) * t79 - mrSges(7,2) * t75 + mrSges(6,1);
t134 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t137 = m(6) + m(7);
t72 = sin(pkin(6));
t78 = sin(qJ(1));
t123 = t72 * t78;
t76 = sin(qJ(4));
t80 = cos(qJ(4));
t69 = qJ(2) + pkin(11);
t108 = pkin(6) + t69;
t90 = sin(t108) / 0.2e1;
t109 = pkin(6) - t69;
t96 = sin(t109);
t35 = t90 - t96 / 0.2e1;
t65 = cos(t69);
t82 = cos(qJ(1));
t99 = -t35 * t78 + t65 * t82;
t9 = t123 * t80 - t76 * t99;
t91 = cos(t109) / 0.2e1;
t97 = cos(t108);
t37 = t91 - t97 / 0.2e1;
t73 = cos(pkin(6));
t139 = -t37 * t76 + t73 * t80;
t138 = -m(4) - m(5);
t136 = m(5) * pkin(3) + t80 * mrSges(5,1) - t76 * mrSges(5,2) + mrSges(4,1);
t98 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t133 = mrSges(7,1) * t75 + mrSges(7,2) * t79 - t98;
t70 = pkin(6) + qJ(2);
t57 = cos(t70) / 0.2e1;
t71 = pkin(6) - qJ(2);
t67 = cos(t71);
t40 = t67 / 0.2e1 + t57;
t100 = t35 * t82 + t65 * t78;
t122 = t72 * t82;
t93 = -t100 * t76 - t122 * t80;
t68 = qJ(4) + pkin(12);
t60 = sin(t68);
t64 = cos(t68);
t132 = -t134 * t60 + t135 * t64 + t136;
t129 = sin(t70) / 0.2e1;
t63 = sin(t71);
t127 = pkin(2) * t63;
t77 = sin(qJ(2));
t120 = t77 * t82;
t119 = t78 * t77;
t52 = pkin(2) * t129;
t118 = t127 / 0.2e1 + t52;
t117 = t76 * t123;
t48 = t76 * t122;
t113 = t137 - t138;
t4 = t100 * t64 - t122 * t60;
t3 = -t100 * t60 - t122 * t64;
t61 = sin(t69);
t84 = t97 / 0.2e1 + t91;
t24 = t61 * t82 + t78 * t84;
t81 = cos(qJ(2));
t59 = pkin(2) * t81 + pkin(1);
t47 = t82 * t59;
t58 = pkin(4) * t80 + pkin(3);
t74 = -qJ(5) - pkin(9);
t111 = pkin(4) * t117 - t24 * t74 + t58 * t99 + t47;
t107 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t38 = t40 * pkin(2);
t105 = -pkin(2) * t119 + t38 * t82;
t95 = m(3) * pkin(1) + mrSges(3,1) * t81 + mrSges(2,1);
t94 = -pkin(2) * t120 - t78 * t38;
t39 = t129 - t63 / 0.2e1;
t85 = mrSges(3,1) * t39 + t113 * (-t127 / 0.2e1 + t52 - t72 * (pkin(8) + qJ(3))) + mrSges(2,2);
t36 = t96 / 0.2e1 + t90;
t28 = -t40 * t78 - t120;
t27 = -t40 * t82 + t119;
t21 = t61 * t78 - t82 * t84;
t20 = t37 * t64 + t60 * t73;
t10 = t80 * t99 + t117;
t8 = t123 * t60 + t64 * t99;
t7 = -t123 * t64 + t60 * t99;
t2 = t24 * t75 + t79 * t8;
t1 = t24 * t79 - t75 * t8;
t5 = [(-t28 * mrSges(3,2) - m(4) * t47 - t99 * mrSges(4,1) - m(5) * (pkin(3) * t99 + t47) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t111 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t111) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t134 * t7 - t95 * t82 + (t107 * t72 + t85) * t78 + t98 * t24) * g(2) + (-t27 * mrSges(3,2) - t48 * mrSges(5,1) + t134 * t3 + t136 * t100 + (t113 * t59 + t95) * t78 + ((-mrSges(5,2) * t80 + t107) * t72 + t85) * t82 + t135 * t4 + t133 * t21 + t137 * (-pkin(4) * t48 + t100 * t58 - t21 * t74)) * g(1) (-(t129 + t63 / 0.2e1) * mrSges(3,1) - (t57 - t67 / 0.2e1) * mrSges(3,2) + t138 * t118 - t137 * (t36 * t58 - t37 * t74 + t118) - t133 * t37 - t132 * t36) * g(3) + (t27 * mrSges(3,1) - (-t39 * t82 - t78 * t81) * mrSges(3,2) - t137 * (-t100 * t74 - t21 * t58 + t105) + t138 * t105 - t133 * t100 + t132 * t21) * g(2) + (-t28 * mrSges(3,1) - (t39 * t78 - t81 * t82) * mrSges(3,2) + t138 * t94 - t137 * (-t24 * t58 - t74 * t99 + t94) - t133 * t99 + t132 * t24) * g(1) (-g(3) * t73 + (-g(1) * t78 + g(2) * t82) * t72) * t113 (-t139 * mrSges(5,1) - (-t37 * t80 - t73 * t76) * mrSges(5,2) + t134 * t20 - t135 * (-t37 * t60 + t64 * t73)) * g(3) + (-t93 * mrSges(5,1) - (-t100 * t80 + t48) * mrSges(5,2) - t135 * t3 + t134 * t4) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t134 * t8 + t135 * t7) * g(1) + (-g(1) * t9 - g(2) * t93 - g(3) * t139) * t137 * pkin(4), t137 * (-g(1) * t24 - g(2) * t21 + g(3) * t36) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t21 * t79 - t4 * t75) * mrSges(7,1) + (-t21 * t75 - t4 * t79) * mrSges(7,2)) - g(3) * ((-t20 * t75 - t36 * t79) * mrSges(7,1) + (-t20 * t79 + t36 * t75) * mrSges(7,2))];
taug  = t5(:);
