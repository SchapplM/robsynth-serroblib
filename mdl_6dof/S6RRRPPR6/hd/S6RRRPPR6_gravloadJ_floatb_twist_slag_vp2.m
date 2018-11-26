% Calculate Gravitation load on the joints for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:36:41
% EndTime: 2018-11-23 17:36:42
% DurationCPUTime: 1.23s
% Computational Cost: add. (1411->146), mult. (1512->191), div. (0->0), fcn. (1466->16), ass. (0->73)
t134 = m(6) + m(7);
t133 = mrSges(5,2) - mrSges(6,3);
t85 = -t134 * qJ(5) + t133;
t66 = sin(qJ(6));
t70 = cos(qJ(6));
t91 = -t66 * mrSges(7,1) - t70 * mrSges(7,2);
t131 = t85 + t91;
t138 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t128 = cos(qJ(1));
t63 = sin(pkin(6));
t106 = t63 * t128;
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t110 = pkin(6) + qJ(2);
t93 = sin(t110) / 0.2e1;
t111 = pkin(6) - qJ(2);
t99 = sin(t111);
t43 = t93 - t99 / 0.2e1;
t69 = sin(qJ(1));
t72 = cos(qJ(2));
t83 = t128 * t43 + t69 * t72;
t77 = t71 * t106 + t83 * t67;
t74 = t77 * pkin(3);
t118 = t63 * t69;
t82 = t128 * t72 - t69 * t43;
t13 = t71 * t118 - t82 * t67;
t100 = cos(t111);
t94 = cos(t110) / 0.2e1;
t44 = t94 - t100 / 0.2e1;
t64 = cos(pkin(6));
t135 = t44 * t67 + t64 * t71;
t132 = m(7) * pkin(10) + t138;
t75 = -m(4) * pkin(9) - m(7) * pkin(5) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t129 = t70 * mrSges(7,1) - t66 * mrSges(7,2) - t75;
t62 = qJ(3) + pkin(11);
t59 = sin(t62);
t60 = cos(t62);
t130 = m(4) * pkin(2) + t71 * mrSges(4,1) - t67 * mrSges(4,2) + mrSges(3,1) + t138 * t60 + (-t91 - t133) * t59;
t68 = sin(qJ(2));
t73 = t100 / 0.2e1 + t94;
t30 = -t128 * t73 + t68 * t69;
t127 = t30 * t60;
t33 = t128 * t68 + t69 * t73;
t125 = t33 * t60;
t42 = t93 + t99 / 0.2e1;
t123 = t42 * t60;
t58 = pkin(3) * t71 + pkin(2);
t65 = -qJ(4) - pkin(9);
t116 = -t30 * t58 - t65 * t83;
t115 = -t33 * t58 - t65 * t82;
t114 = t42 * t58 + t44 * t65;
t113 = t128 * pkin(1) + pkin(8) * t118;
t112 = qJ(5) * t59;
t109 = t67 * t118;
t105 = -pkin(1) * t69 + pkin(8) * t106;
t8 = -t59 * t106 + t60 * t83;
t50 = t67 * t106;
t104 = -t71 * t83 + t50;
t103 = -pkin(4) * t127 - t30 * t112 + t116;
t102 = -pkin(4) * t125 - t33 * t112 + t115;
t101 = pkin(4) * t123 + t42 * t112 + t114;
t97 = t13 * pkin(3);
t96 = t135 * pkin(3);
t95 = pkin(3) * t109 - t33 * t65 + t58 * t82 + t113;
t84 = pkin(3) * t50 + t30 * t65 - t58 * t83 + t105;
t7 = t106 * t60 + t59 * t83;
t24 = -t44 * t59 - t64 * t60;
t14 = t71 * t82 + t109;
t12 = t118 * t59 + t60 * t82;
t11 = -t118 * t60 + t59 * t82;
t2 = t11 * t66 + t33 * t70;
t1 = t11 * t70 - t33 * t66;
t3 = [(-t128 * mrSges(2,1) - m(3) * t113 - t82 * mrSges(3,1) - m(4) * (pkin(2) * t82 + t113) - t14 * mrSges(4,1) - t13 * mrSges(4,2) - m(5) * t95 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t63 + mrSges(2,2)) * t69 + t85 * t11 - t132 * t12 + t75 * t33 - t134 * (t12 * pkin(4) + t95)) * g(2) + (t69 * mrSges(2,1) + t128 * mrSges(2,2) - m(3) * t105 + t83 * mrSges(3,1) - mrSges(3,3) * t106 - m(4) * (-pkin(2) * t83 + t105) - t104 * mrSges(4,1) - t77 * mrSges(4,2) - m(5) * t84 + t132 * t8 - t131 * t7 + t129 * t30 + t134 * (pkin(4) * t8 - t84)) * g(1) (-m(5) * t114 - m(6) * t101 - m(7) * (pkin(10) * t123 + t101) + t129 * t44 - t130 * t42) * g(3) + (-m(5) * t116 - m(6) * t103 - m(7) * (-pkin(10) * t127 + t103) - t129 * t83 + t130 * t30) * g(2) + (-m(5) * t115 - m(6) * t102 - m(7) * (-pkin(10) * t125 + t102) - t129 * t82 + t130 * t33) * g(1) (-t135 * mrSges(4,1) - (t44 * t71 - t64 * t67) * mrSges(4,2) - m(5) * t96 - t134 * (-t24 * pkin(4) + t96) + t131 * (-t44 * t60 + t59 * t64) + t132 * t24) * g(3) + (m(5) * t74 + t77 * mrSges(4,1) - t104 * mrSges(4,2) + t131 * t8 + t132 * t7 + t134 * (t7 * pkin(4) + t74)) * g(2) + (-m(5) * t97 - t13 * mrSges(4,1) + t14 * mrSges(4,2) - t134 * (-t11 * pkin(4) + t97) + t131 * t12 + t132 * t11) * g(1) (m(5) + t134) * (-g(1) * t33 - g(2) * t30 + g(3) * t42) t134 * (-g(1) * t11 - g(2) * t7 - g(3) * t24) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t30 * t66 + t7 * t70) * mrSges(7,1) + (-t30 * t70 - t66 * t7) * mrSges(7,2)) - g(3) * ((t24 * t70 + t42 * t66) * mrSges(7,1) + (-t24 * t66 + t42 * t70) * mrSges(7,2))];
taug  = t3(:);
