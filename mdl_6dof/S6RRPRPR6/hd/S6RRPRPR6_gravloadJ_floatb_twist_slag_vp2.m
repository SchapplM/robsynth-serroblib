% Calculate Gravitation load on the joints for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:03:42
% EndTime: 2018-11-23 17:03:44
% DurationCPUTime: 1.15s
% Computational Cost: add. (1604->150), mult. (1372->188), div. (0->0), fcn. (1290->20), ass. (0->96)
t139 = mrSges(5,2) - mrSges(6,3);
t140 = m(6) + m(7);
t87 = -qJ(5) * t140 + t139;
t68 = sin(qJ(6));
t72 = cos(qJ(6));
t96 = -t68 * mrSges(7,1) - t72 * mrSges(7,2);
t130 = t87 + t96;
t144 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t131 = m(7) * pkin(10) + t144;
t141 = pkin(4) * t140 + t131;
t69 = sin(qJ(4));
t106 = qJ(5) * t69;
t63 = qJ(2) + pkin(11);
t57 = sin(t63);
t71 = sin(qJ(1));
t75 = cos(qJ(1));
t100 = pkin(6) - t63;
t85 = cos(t100) / 0.2e1;
t99 = pkin(6) + t63;
t90 = cos(t99);
t77 = t90 / 0.2e1 + t85;
t21 = t57 * t71 - t75 * t77;
t73 = cos(qJ(4));
t120 = t21 * t73;
t138 = -pkin(4) * t120 - t21 * t106;
t24 = t75 * t57 + t71 * t77;
t119 = t24 * t73;
t137 = -pkin(4) * t119 - t24 * t106;
t51 = sin(t99);
t127 = t51 / 0.2e1;
t89 = sin(t100);
t84 = t89 / 0.2e1;
t39 = t84 + t127;
t118 = t39 * t73;
t136 = pkin(4) * t118 + t39 * t106;
t64 = pkin(6) + qJ(2);
t54 = cos(t64) / 0.2e1;
t65 = pkin(6) - qJ(2);
t62 = cos(t65);
t43 = t62 / 0.2e1 + t54;
t135 = mrSges(4,2) - mrSges(6,1) - mrSges(5,3);
t105 = m(5) + t140;
t133 = -m(7) * pkin(5) - pkin(9) * t105 + t135;
t132 = t72 * mrSges(7,1) - t68 * mrSges(7,2);
t129 = mrSges(4,1) + t144 * t73 + (-t96 - t139) * t69;
t128 = -m(7) * (pkin(5) + pkin(9)) - t132 + t135;
t126 = sin(t64) / 0.2e1;
t59 = sin(t65);
t122 = pkin(2) * t59;
t66 = sin(pkin(6));
t117 = t66 * t71;
t116 = t66 * t75;
t70 = sin(qJ(2));
t113 = t70 * t75;
t60 = cos(t63);
t112 = t71 * t60;
t111 = t71 * t70;
t109 = t75 * t60;
t104 = t127 - t89 / 0.2e1;
t25 = -t104 * t71 + t109;
t74 = cos(qJ(2));
t56 = pkin(2) * t74 + pkin(1);
t44 = t75 * t56;
t108 = t25 * pkin(3) + t44;
t49 = pkin(2) * t126;
t107 = t122 / 0.2e1 + t49;
t102 = t39 * pkin(3) + t107;
t101 = m(4) + t105;
t22 = t104 * t75 + t112;
t8 = -t69 * t116 + t22 * t73;
t41 = t43 * pkin(2);
t98 = -pkin(2) * t111 + t75 * t41;
t92 = -t21 * pkin(3) + t98;
t40 = t85 - t90 / 0.2e1;
t91 = pkin(9) * t40 + t102;
t88 = m(3) * pkin(1) + t74 * mrSges(3,1) + mrSges(2,1);
t86 = -pkin(2) * t113 - t71 * t41;
t7 = t116 * t73 + t22 * t69;
t83 = -t24 * pkin(3) + t86;
t79 = -t51 / 0.2e1 + t84;
t23 = -t75 * t79 + t112;
t80 = pkin(9) * t23 + t92;
t26 = t71 * t79 + t109;
t78 = t26 * pkin(9) + t83;
t42 = t126 - t59 / 0.2e1;
t76 = t42 * mrSges(3,1) + mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t66 + t101 * (-t122 / 0.2e1 + t49 - t66 * (pkin(8) + qJ(3)));
t67 = cos(pkin(6));
t31 = -t71 * t43 - t113;
t30 = -t43 * t75 + t111;
t28 = t40 * t69 - t67 * t73;
t18 = t22 * pkin(3);
t12 = t117 * t69 + t25 * t73;
t11 = -t117 * t73 + t25 * t69;
t2 = t11 * t68 + t24 * t72;
t1 = t11 * t72 - t24 * t68;
t3 = [(-t31 * mrSges(3,2) - m(4) * t44 - t25 * mrSges(4,1) - m(5) * t108 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t87 * t11 - t88 * t75 + t76 * t71 - t131 * t12 + t133 * t24 - t140 * (t12 * pkin(4) + t108)) * g(2) + (-t30 * mrSges(3,2) + t22 * mrSges(4,1) + m(5) * t18 + (t101 * t56 + t88) * t71 + t76 * t75 + t131 * t8 - t130 * t7 + (t132 - t133) * t21 + t140 * (pkin(4) * t8 + t18)) * g(1) (-(t126 + t59 / 0.2e1) * mrSges(3,1) - (t54 - t62 / 0.2e1) * mrSges(3,2) - m(4) * t107 - m(5) * t91 - m(6) * (t91 + t136) - m(7) * (pkin(10) * t118 + t102 + t136) + t128 * t40 - t129 * t39) * g(3) + (t30 * mrSges(3,1) - (-t42 * t75 - t71 * t74) * mrSges(3,2) - m(4) * t98 - m(5) * t80 - m(6) * (t80 + t138) - m(7) * (-pkin(10) * t120 + t138 + t92) + t128 * t23 + t129 * t21) * g(2) + (-t31 * mrSges(3,1) - (t71 * t42 - t74 * t75) * mrSges(3,2) - m(4) * t86 - m(5) * t78 - m(6) * (t78 + t137) - m(7) * (-pkin(10) * t119 + t137 + t83) + t128 * t26 + t129 * t24) * g(1) (-t67 * g(3) + (-g(1) * t71 + g(2) * t75) * t66) * t101 (t130 * (t40 * t73 + t67 * t69) + t141 * t28) * g(3) + (t130 * t8 + t141 * t7) * g(2) + (t11 * t141 + t130 * t12) * g(1), t140 * (-g(1) * t11 - g(2) * t7 - g(3) * t28) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t21 * t68 + t7 * t72) * mrSges(7,1) + (-t21 * t72 - t68 * t7) * mrSges(7,2)) - g(3) * ((t28 * t72 + t39 * t68) * mrSges(7,1) + (-t28 * t68 + t39 * t72) * mrSges(7,2))];
taug  = t3(:);
