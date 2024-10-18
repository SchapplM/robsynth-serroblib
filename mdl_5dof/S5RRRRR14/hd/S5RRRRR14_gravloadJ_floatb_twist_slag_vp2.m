% Calculate Gravitation load on the joints for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:24
% EndTime: 2024-09-27 18:42:24
% DurationCPUTime: 0.45s
% Computational Cost: add. (723->117), mult. (507->150), div. (0->0), fcn. (450->22), ass. (0->79)
t124 = -m(5) * pkin(3) - mrSges(4,1);
t123 = -mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t76 = qJ(3) + qJ(4);
t68 = pkin(5) + t76;
t101 = sin(t68) / 0.2e1;
t63 = qJ(5) + t68;
t102 = sin(t63) / 0.2e1;
t55 = cos(t63) / 0.2e1;
t98 = pkin(5) - t76;
t97 = -qJ(5) + t98;
t59 = cos(t97);
t89 = sin(t97);
t105 = (t102 + t89 / 0.2e1) * mrSges(6,1) + (t55 - t59 / 0.2e1) * mrSges(6,2);
t58 = cos(t68) / 0.2e1;
t61 = cos(t98);
t96 = sin(t98);
t122 = -(t101 + t96 / 0.2e1) * mrSges(5,1) - (t58 - t61 / 0.2e1) * mrSges(5,2) - t105;
t36 = t59 / 0.2e1 + t55;
t73 = qJ(5) + t76;
t64 = sin(t73);
t77 = qJ(1) + qJ(2);
t70 = sin(t77);
t72 = cos(t77);
t10 = -t70 * t36 - t72 * t64;
t35 = t102 - t89 / 0.2e1;
t65 = cos(t73);
t94 = t70 * t35 - t72 * t65;
t116 = t10 * mrSges(6,1) + t94 * mrSges(6,2);
t69 = sin(t76);
t112 = t72 * t69;
t40 = t61 / 0.2e1 + t58;
t20 = -t70 * t40 - t112;
t39 = t101 - t96 / 0.2e1;
t71 = cos(t76);
t92 = t70 * t39 - t72 * t71;
t121 = -t20 * mrSges(5,1) - t92 * mrSges(5,2) - t116;
t9 = -t72 * t36 + t70 * t64;
t95 = -t72 * t35 - t70 * t65;
t117 = -t9 * mrSges(6,1) + t95 * mrSges(6,2);
t114 = t70 * t69;
t19 = -t72 * t40 + t114;
t93 = -t72 * t39 - t70 * t71;
t120 = t19 * mrSges(5,1) - t93 * mrSges(5,2) - t117;
t118 = pkin(8) + pkin(9);
t78 = sin(pkin(5));
t115 = m(6) * t78;
t84 = cos(qJ(3));
t74 = t84 * pkin(3);
t67 = t74 + pkin(2);
t113 = t70 * t78;
t79 = cos(pkin(5));
t81 = sin(qJ(3));
t110 = t79 * t81;
t109 = t79 * t84;
t80 = sin(qJ(4));
t108 = t80 * t81;
t103 = t72 * pkin(2) + pkin(8) * t113;
t83 = cos(qJ(4));
t66 = t83 * pkin(4) + pkin(3);
t21 = -t78 * (pkin(10) + t118) + (t84 * t80 * pkin(4) + t81 * t66) * t79;
t42 = pkin(4) * t71 + t67;
t100 = -t70 * t21 + t72 * t42;
t41 = pkin(3) * t110 - t78 * t118;
t99 = -t70 * t41 + t72 * t67;
t91 = -pkin(4) * t108 + t66 * t84;
t90 = t72 * t109 - t70 * t81;
t30 = -t70 * t109 - t72 * t81;
t88 = pkin(4) * (t83 * t84 - t108);
t31 = -t70 * t110 + t72 * t84;
t87 = -t72 * mrSges(3,1) - t31 * mrSges(4,1) + t92 * mrSges(5,1) + t94 * mrSges(6,1) + t70 * mrSges(3,2) - t30 * mrSges(4,2) - t20 * mrSges(5,2) - t10 * mrSges(6,2) + t113 * t123;
t29 = -t72 * t110 - t70 * t84;
t86 = -t29 * mrSges(4,1) - t93 * mrSges(5,1) - t95 * mrSges(6,1) + t90 * mrSges(4,2) - t19 * mrSges(5,2) - t9 * mrSges(6,2) + (m(4) * pkin(2) + m(5) * t67 + m(6) * t42 + mrSges(3,1)) * t70 + (m(5) * t41 + m(6) * t21 + mrSges(3,2) + (-m(4) * pkin(8) + t123) * t78) * t72;
t85 = cos(qJ(1));
t82 = sin(qJ(1));
t75 = t85 * pkin(1);
t43 = -t81 * pkin(3) - pkin(4) * t69;
t27 = t79 * t88;
t22 = t91 * t79;
t1 = [(-m(4) * (t75 + t103) - m(5) * (t75 + t99) - m(6) * (t100 + t75) + t82 * mrSges(2,2) + t87 + (-m(3) * pkin(1) - mrSges(2,1)) * t85) * g(2) + (t85 * mrSges(2,2) + t86 + (mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1)) * t82) * g(1), (-m(4) * t103 - m(5) * t99 - m(6) * t100 + t87) * g(2) + t86 * g(1), (-t91 * t115 + (-m(5) * t74 - mrSges(4,1) * t84 + mrSges(4,2) * t81) * t78 + t122) * g(3) + (-t29 * mrSges(4,2) - m(6) * (t72 * t22 + t70 * t43) + t120 + t124 * t90) * g(2) + (t31 * mrSges(4,2) - m(6) * (-t70 * t22 + t72 * t43) + t124 * t30 + t121) * g(1), (-t88 * t115 + t122) * g(3) + (-m(6) * (-pkin(4) * t114 + t72 * t27) + t120) * g(2) + (-m(6) * (-pkin(4) * t112 - t70 * t27) + t121) * g(1), -g(1) * t116 - g(2) * t117 - g(3) * t105];
taug = t1(:);
