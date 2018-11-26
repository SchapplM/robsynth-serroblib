% Calculate Gravitation load on the joints for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:19:37
% EndTime: 2018-11-23 17:19:38
% DurationCPUTime: 1.05s
% Computational Cost: add. (1292->127), mult. (1567->167), div. (0->0), fcn. (1529->14), ass. (0->73)
t94 = m(5) + m(6) + m(7);
t90 = m(4) + t94;
t114 = qJ(3) * t90 - mrSges(3,2) + mrSges(4,3);
t122 = mrSges(6,1) + mrSges(7,1);
t119 = mrSges(6,2) + mrSges(7,2);
t59 = cos(qJ(5));
t50 = pkin(5) * t59 + pkin(4);
t121 = m(6) * pkin(4) + m(7) * t50 + mrSges(5,1);
t120 = m(6) * pkin(10) - m(7) * (-qJ(6) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t110 = m(7) * pkin(5);
t118 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t55 = sin(qJ(5));
t117 = -t119 * t55 + t122 * t59 + t121;
t116 = -t110 - t122;
t115 = pkin(9) * t94 + t55 * t110 + t118;
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t100 = t58 * t61;
t62 = cos(qJ(1));
t92 = pkin(6) + qJ(2);
t80 = sin(t92);
t76 = t80 / 0.2e1;
t93 = pkin(6) - qJ(2);
t81 = sin(t93);
t72 = t76 - t81 / 0.2e1;
t31 = t62 * t72 + t100;
t52 = sin(pkin(6));
t103 = t52 * t62;
t57 = sin(qJ(2));
t78 = cos(t93) / 0.2e1;
t82 = cos(t92);
t65 = t78 + t82 / 0.2e1;
t30 = t57 * t58 - t62 * t65;
t56 = sin(qJ(4));
t60 = cos(qJ(4));
t71 = t60 * t103 - t30 * t56;
t113 = t31 * t59 + t55 * t71;
t112 = -t31 * t55 + t59 * t71;
t111 = t120 * t60 - t121 * t56 - t114;
t109 = t30 * t55;
t33 = t62 * t57 + t58 * t65;
t106 = t33 * t55;
t77 = t81 / 0.2e1;
t41 = t76 + t77;
t105 = t41 * t55;
t104 = t52 * t58;
t102 = t55 * t56;
t101 = t56 * t59;
t99 = t62 * t61;
t95 = t62 * pkin(1) + pkin(8) * t104;
t34 = -t58 * t72 + t99;
t91 = t34 * pkin(2) + t95;
t87 = -pkin(1) * t58 + pkin(8) * t103;
t83 = pkin(3) * t104 + t91;
t79 = -t31 * pkin(2) + t87;
t75 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t52;
t14 = t104 * t60 + t33 * t56;
t1 = -t14 * t55 + t34 * t59;
t74 = pkin(3) * t103 + t79;
t73 = t77 - t80 / 0.2e1;
t15 = t56 * t103 + t30 * t60;
t53 = cos(pkin(6));
t42 = t78 - t82 / 0.2e1;
t40 = t41 * pkin(2);
t35 = t58 * t73 + t99;
t32 = -t62 * t73 + t100;
t29 = -t41 * t56 + t53 * t60;
t28 = t41 * t60 + t53 * t56;
t26 = t33 * pkin(2);
t24 = t30 * pkin(2);
t13 = t104 * t56 - t33 * t60;
t2 = t14 * t59 + t34 * t55;
t3 = [(-t62 * mrSges(2,1) - m(3) * t95 - m(4) * t91 - m(5) * t83 - t14 * mrSges(5,1) - m(6) * (pkin(4) * t14 + t83) - m(7) * (t14 * t50 + t83) - t114 * t33 - t122 * t2 - t119 * t1 + t75 * t58 - t120 * t13 - t115 * t34) * g(2) + (t58 * mrSges(2,1) - m(3) * t87 - m(4) * t79 - m(5) * t74 - t71 * mrSges(5,1) - m(6) * (pkin(4) * t71 + t74) - m(7) * (t50 * t71 + t74) - t122 * t112 + t114 * t30 + t119 * t113 + t75 * t62 - t120 * t15 + t115 * t31) * g(1) (-t105 * t110 - m(4) * t40 - t94 * (pkin(9) * t41 + t40) - t122 * (t101 * t42 + t105) - t119 * (-t102 * t42 + t41 * t59) - t118 * t41 + t111 * t42) * g(3) + (t109 * t110 + m(4) * t24 - t94 * (-pkin(9) * t30 - t24) - t122 * (t101 * t32 - t109) - t119 * (-t102 * t32 - t30 * t59) + t118 * t30 + t111 * t32) * g(2) + (t106 * t110 + m(4) * t26 - t119 * (-t102 * t35 - t33 * t59) - t94 * (-pkin(9) * t33 - t26) - t122 * (t101 * t35 - t106) + t118 * t33 + t111 * t35) * g(1) (-g(1) * t33 - g(2) * t30 + g(3) * t41) * t90 (t117 * t28 - t120 * t29) * g(3) + (-t117 * t15 + t120 * t71) * g(2) + (t117 * t13 - t120 * t14) * g(1) (-t119 * (-t29 * t59 - t42 * t55) + t116 * (-t29 * t55 + t42 * t59)) * g(3) + (-t119 * t112 + t116 * t113) * g(2) + (t116 * t1 + t119 * t2) * g(1) (-g(1) * t13 + g(2) * t15 - g(3) * t28) * m(7)];
taug  = t3(:);
