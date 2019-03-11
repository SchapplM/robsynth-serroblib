% Calculate Gravitation load on the joints for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:31
% EndTime: 2019-03-09 13:50:34
% DurationCPUTime: 1.08s
% Computational Cost: add. (527->162), mult. (828->183), div. (0->0), fcn. (847->10), ass. (0->92)
t135 = mrSges(6,2) - mrSges(7,3);
t62 = cos(qJ(2));
t101 = qJ(3) * t62;
t59 = sin(qJ(1));
t43 = t59 * t101;
t61 = cos(qJ(4));
t48 = pkin(4) * t61 + pkin(3);
t114 = -pkin(2) - t48;
t58 = sin(qJ(2));
t95 = t58 * t114;
t134 = t59 * t95 + t43;
t133 = -m(4) - m(5);
t131 = (mrSges(3,1) + mrSges(4,1)) * t62 + (-mrSges(3,2) + mrSges(4,3)) * t58;
t56 = sin(qJ(6));
t60 = cos(qJ(6));
t130 = t60 * mrSges(7,1) - t56 * mrSges(7,2);
t129 = -mrSges(2,1) - t131;
t128 = pkin(10) * m(7) - t135;
t127 = pkin(8) * m(5) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t100 = qJ(4) + qJ(5);
t89 = sin(t100);
t90 = cos(t100);
t27 = t58 * t90 - t62 * t89;
t63 = cos(qJ(1));
t54 = t63 * pkin(7);
t57 = sin(qJ(4));
t64 = -pkin(9) - pkin(8);
t126 = t63 * t64 + t54 + (-pkin(1) + t114 * t62 + (-pkin(4) * t57 - qJ(3)) * t58) * t59;
t125 = g(1) * t63 + g(2) * t59;
t13 = t27 * t59;
t26 = t58 * t89 + t62 * t90;
t14 = t26 * t59;
t124 = t13 * mrSges(6,1) - t135 * t14;
t82 = t63 * t89;
t83 = t63 * t90;
t15 = -t58 * t82 - t62 * t83;
t16 = -t58 * t83 + t62 * t82;
t123 = -t16 * mrSges(6,1) + t135 * t15;
t122 = -t26 * mrSges(6,1) - t135 * t27;
t121 = -pkin(2) - pkin(3);
t120 = pkin(4) * t59;
t119 = pkin(10) * t14;
t118 = pkin(10) * t15;
t117 = pkin(10) * t27;
t53 = t62 * pkin(2);
t112 = t57 * t58;
t111 = t57 * t62;
t109 = t58 * t61;
t108 = t58 * t63;
t106 = t61 * t62;
t105 = t62 * t63;
t51 = t58 * qJ(3);
t103 = t53 + t51;
t102 = t63 * pkin(1) + t59 * pkin(7);
t46 = pkin(4) * t112;
t99 = t57 * t105;
t98 = t61 * t108;
t97 = t58 * t121;
t36 = t111 * t120;
t94 = -t36 + t119;
t39 = pkin(4) * t99;
t93 = -t39 - t118;
t92 = -t26 * pkin(5) + t117;
t91 = -pkin(1) - t51;
t88 = t62 * t48 + t103 + t46;
t87 = pkin(2) * t105 + t63 * t51 + t102;
t86 = -pkin(4) * t106 - t46;
t73 = -t109 + t111;
t22 = t73 * t59;
t72 = t106 + t112;
t23 = t72 * t59;
t78 = t22 * mrSges(5,1) + t23 * mrSges(5,2);
t24 = -t98 + t99;
t25 = t72 * t63;
t77 = t24 * mrSges(5,1) + t25 * mrSges(5,2);
t76 = t72 * mrSges(5,1) - t73 * mrSges(5,2);
t75 = -t14 * t60 - t56 * t63;
t74 = t14 * t56 - t60 * t63;
t71 = t48 * t105 + t63 * t46 + t59 * t64 + t87;
t70 = -t130 * t13 - t124;
t69 = t130 * t16 - t123;
t68 = m(7) * pkin(5) + t130;
t67 = t130 * t26 - t122;
t66 = t62 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t58;
t44 = t63 * t101;
t38 = pkin(4) * t98;
t35 = t109 * t120;
t12 = t16 * pkin(5);
t11 = t13 * pkin(5);
t2 = -t15 * t60 - t56 * t59;
t1 = t15 * t56 - t59 * t60;
t3 = [(-m(3) * t102 - m(4) * t87 - m(5) * (pkin(3) * t105 + t87) - t25 * mrSges(5,1) + t24 * mrSges(5,2) - m(6) * t71 + t15 * mrSges(6,1) - m(7) * (-pkin(5) * t15 + t71) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t128 * t16 + t129 * t63 + t127 * t59) * g(2) + (t23 * mrSges(5,1) - t22 * mrSges(5,2) - m(6) * t126 + t14 * mrSges(6,1) - t75 * mrSges(7,1) - t74 * mrSges(7,2) - (-pkin(5) * t14 + t126) * m(7) - t128 * t13 + (-m(3) + t133) * t54 + (m(3) * pkin(1) - m(4) * (t91 - t53) - m(5) * (t121 * t62 + t91) - t129) * t59 + t127 * t63) * g(1), t125 * (mrSges(3,1) * t58 + mrSges(3,2) * t62) + (-m(4) * t43 - t59 * t66 - m(5) * (t59 * t97 + t43) - t78 - m(6) * (t36 + t134) - m(7) * (-t94 + t134) + t68 * t13 + t124) * g(2) + (-m(4) * t44 - t63 * t66 - m(5) * (t63 * t97 + t44) - t77 - m(6) * (t63 * t95 + t39 + t44) - m(7) * (t114 * t108 + t44 - t93) - t68 * t16 + t123) * g(1) + (-m(4) * t103 - m(5) * (pkin(3) * t62 + t103) - t76 - m(6) * t88 - m(7) * (t88 - t117) - t68 * t26 + t122 - t131) * g(3) (t62 * g(3) - t125 * t58) * (m(6) + m(7) - t133) (-m(6) * t86 - m(7) * (t86 + t92) + t67 + t76) * g(3) + (-m(6) * (-t36 + t35) - m(7) * (t11 + t35 + t94) + t70 + t78) * g(2) + (-m(6) * (-t39 + t38) - m(7) * (-t12 + t38 + t93) + t69 + t77) * g(1) (-m(7) * t92 + t67) * g(3) + (-m(7) * (t11 + t119) + t70) * g(2) + (-m(7) * (-t12 - t118) + t69) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t74 * mrSges(7,1) + t75 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t56 - mrSges(7,2) * t60) * t27];
taug  = t3(:);
