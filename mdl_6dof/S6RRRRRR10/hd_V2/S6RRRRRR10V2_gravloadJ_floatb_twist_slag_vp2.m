% Calculate Gravitation load on the joints for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:38
% EndTime: 2019-04-11 14:41:43
% DurationCPUTime: 1.00s
% Computational Cost: add. (696->134), mult. (1006->184), div. (0->0), fcn. (1076->12), ass. (0->73)
t137 = -mrSges(7,3) + mrSges(6,2);
t69 = sin(qJ(6));
t74 = cos(qJ(6));
t87 = mrSges(7,1) * t74 - mrSges(7,2) * t69 + mrSges(6,1);
t95 = -m(7) * pkin(6) + t137;
t70 = sin(qJ(5));
t75 = cos(qJ(5));
t126 = -t95 * t70 + t87 * t75 + mrSges(5,1);
t71 = sin(qJ(4));
t76 = cos(qJ(4));
t135 = t76 * mrSges(5,1) + t71 * mrSges(6,3);
t68 = qJ(2) + qJ(3);
t65 = sin(t68);
t73 = sin(qJ(1));
t117 = t65 * t73;
t78 = cos(qJ(1));
t101 = t78 * t71;
t105 = t73 * t76;
t66 = cos(t68);
t40 = t105 * t66 - t101;
t10 = t117 * t70 + t40 * t75;
t102 = t76 * t78;
t106 = t73 * t71;
t39 = t106 * t66 + t102;
t134 = t10 * t69 - t39 * t74;
t133 = -t10 * t74 - t39 * t69;
t111 = t70 * t76;
t88 = t111 * t65 + t66 * t75;
t132 = t88 * t95;
t122 = -m(5) - m(6);
t130 = mrSges(5,2) - mrSges(6,3);
t100 = t66 * pkin(3) + t65 * pkin(5);
t129 = -t66 * mrSges(4,1) + t65 * mrSges(4,2);
t128 = t122 - m(7);
t72 = sin(qJ(2));
t107 = t72 * mrSges(3,2);
t77 = cos(qJ(2));
t67 = t77 * pkin(2);
t64 = t67 + pkin(1);
t127 = m(3) * pkin(1) + m(4) * t64 + t77 * mrSges(3,1) + mrSges(2,1) - t107 - t129;
t125 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t124 = -mrSges(7,1) * t69 - mrSges(7,2) * t74 + t130;
t116 = t65 * t75;
t115 = t65 * t78;
t114 = t66 * t73;
t113 = t66 * t78;
t112 = t69 * t71;
t109 = t71 * t74;
t104 = t75 * t76;
t99 = t65 * t106;
t98 = t65 * t109;
t97 = t65 * t101;
t34 = t111 * t66 - t116;
t94 = pkin(6) * t34 + t100;
t92 = -pkin(2) * t72 - pkin(3) * t65;
t11 = t116 * t73 - t40 * t70;
t33 = t104 * t65 - t66 * t70;
t26 = t33 * t73;
t52 = pkin(5) * t114;
t84 = t26 * mrSges(6,1) - (t26 * t69 - t73 * t98) * mrSges(7,2) - (-t26 * t74 - t69 * t99) * mrSges(7,1) - mrSges(5,2) * t99 - mrSges(5,3) * t114 - m(7) * t52 - t73 * t132;
t28 = t33 * t78;
t54 = pkin(5) * t113;
t83 = t28 * mrSges(6,1) - mrSges(5,2) * t97 - (t28 * t69 - t74 * t97) * mrSges(7,2) - mrSges(5,3) * t113 - (-t28 * t74 - t69 * t97) * mrSges(7,1) - m(7) * t54 - t78 * t132;
t81 = mrSges(4,2) * t66 + (m(7) * pkin(3) + mrSges(4,1) + t135) * t65;
t80 = -t65 * mrSges(5,3) + t129 + t137 * t34 - t87 * (t104 * t66 + t65 * t70) + (-t112 * mrSges(7,1) + mrSges(5,2) * t71 - t109 * mrSges(7,2) - t135) * t66;
t79 = mrSges(3,2) * t77 + (mrSges(3,1) + (m(4) + m(7)) * pkin(2)) * t72 + t81;
t42 = t102 * t66 + t106;
t41 = t101 * t66 - t105;
t14 = t115 * t70 + t42 * t75;
t13 = -t115 * t75 + t42 * t70;
t2 = t14 * t74 + t41 * t69;
t1 = -t14 * t69 + t41 * t74;
t3 = [(-t42 * mrSges(5,1) - t14 * mrSges(6,1) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(5,3) * t115 + t130 * t41 + t128 * (pkin(3) * t113 + pkin(5) * t115 + t78 * t64) + t95 * t13 - t127 * t78 + t125 * t73) * g(2) + (t40 * mrSges(5,1) + mrSges(5,3) * t117 + t10 * mrSges(6,1) - t133 * mrSges(7,1) - t134 * mrSges(7,2) - t130 * t39 + t95 * t11 + t125 * t78 + (t127 + t128 * (-t64 - t100)) * t73) * g(1) (t80 + t122 * (t67 + t100) + t107 - m(7) * (t67 + t94) + (-m(4) * pkin(2) - mrSges(3,1)) * t77) * g(3) + (t122 * (t73 * t92 + t52) + t79 * t73 + t84) * g(2) + (t122 * (t78 * t92 + t54) + t79 * t78 + t83) * g(1) (-m(7) * t94 + t122 * t100 + t80) * g(3) + (t122 * (-pkin(3) * t117 + t52) + t81 * t73 + t84) * g(2) + (t122 * (-pkin(3) * t115 + t54) + t81 * t78 + t83) * g(1) (t124 * t76 + t126 * t71) * g(3) * t65 + (t124 * t40 + t126 * t39) * g(2) + (t124 * t42 + t126 * t41) * g(1) (t33 * t95 + t87 * t88) * g(3) + (t10 * t95 - t11 * t87) * g(2) + (t13 * t87 + t14 * t95) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t134 * mrSges(7,1) + t133 * mrSges(7,2)) - g(3) * ((-t33 * t69 + t98) * mrSges(7,1) + (-t112 * t65 - t33 * t74) * mrSges(7,2))];
taug  = t3(:);
