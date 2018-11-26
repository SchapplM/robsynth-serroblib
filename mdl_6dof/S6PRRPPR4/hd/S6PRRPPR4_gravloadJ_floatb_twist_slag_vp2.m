% Calculate Gravitation load on the joints for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:59
% EndTime: 2018-11-23 15:10:01
% DurationCPUTime: 1.17s
% Computational Cost: add. (1348->124), mult. (1690->168), div. (0->0), fcn. (1706->16), ass. (0->82)
t136 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2) + mrSges(7,3);
t58 = sin(pkin(11));
t59 = cos(pkin(11));
t134 = -pkin(4) * t59 - qJ(5) * t58;
t112 = cos(qJ(3));
t61 = sin(qJ(3));
t133 = -mrSges(4,1) * t112 - mrSges(3,1) + (m(7) * pkin(9) + t136) * t61;
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t132 = t60 * mrSges(7,1) + t63 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t131 = m(7) * pkin(5) + mrSges(7,1) * t63 - mrSges(7,2) * t60 + mrSges(5,1) + mrSges(6,1);
t130 = t112 * pkin(3) + qJ(4) * t61;
t127 = m(6) + m(7);
t126 = mrSges(3,2) - mrSges(4,3);
t104 = cos(pkin(10));
t100 = pkin(6) + qJ(2);
t86 = sin(t100);
t82 = t86 / 0.2e1;
t101 = pkin(6) - qJ(2);
t87 = sin(t101);
t74 = t82 - t87 / 0.2e1;
t102 = sin(pkin(10));
t64 = cos(qJ(2));
t92 = t102 * t64;
t38 = t104 * t74 + t92;
t103 = sin(pkin(6));
t76 = t104 * t103;
t20 = t112 * t76 + t38 * t61;
t124 = t134 * t20;
t93 = t104 * t64;
t41 = -t102 * t74 + t93;
t75 = t103 * t102;
t22 = -t112 * t75 + t41 * t61;
t123 = t134 * t22;
t105 = cos(pkin(6));
t84 = cos(t100) / 0.2e1;
t88 = cos(t101);
t52 = t84 - t88 / 0.2e1;
t43 = -t105 * t112 - t52 * t61;
t122 = t134 * t43;
t121 = m(5) + t127;
t119 = t131 * t59 + t132 * t58 + mrSges(4,1);
t118 = -m(7) * (-pkin(9) + qJ(4)) + t136;
t98 = t58 * t112;
t97 = t59 * t112;
t62 = sin(qJ(2));
t67 = t88 / 0.2e1 + t84;
t37 = t102 * t62 - t104 * t67;
t83 = t87 / 0.2e1;
t66 = t83 - t86 / 0.2e1;
t39 = -t104 * t66 + t92;
t96 = -t37 * pkin(2) + t39 * pkin(8);
t40 = t102 * t67 + t104 * t62;
t42 = t102 * t66 + t93;
t95 = -t40 * pkin(2) + t42 * pkin(8);
t51 = t82 + t83;
t94 = t51 * pkin(2) - t52 * pkin(8);
t17 = t20 * pkin(3);
t21 = t112 * t38 - t61 * t76;
t91 = qJ(4) * t21 - t17;
t18 = t22 * pkin(3);
t23 = t112 * t41 + t61 * t75;
t90 = qJ(4) * t23 - t18;
t36 = t43 * pkin(3);
t44 = t105 * t61 - t112 * t52;
t89 = qJ(4) * t44 - t36;
t79 = -t130 * t37 + t96;
t78 = -t130 * t40 + t95;
t77 = t130 * t51 + t94;
t25 = t51 * t97 - t52 * t58;
t24 = t51 * t98 + t52 * t59;
t12 = t44 * t59 - t51 * t58;
t11 = t44 * t58 + t51 * t59;
t10 = -t40 * t97 + t42 * t58;
t9 = -t40 * t98 - t42 * t59;
t8 = -t37 * t97 + t39 * t58;
t7 = -t37 * t98 - t39 * t59;
t4 = t23 * t59 + t40 * t58;
t3 = t23 * t58 - t40 * t59;
t2 = t21 * t59 + t37 * t58;
t1 = t21 * t58 - t37 * t59;
t5 = [(-m(2) - m(3) - m(4) - t121) * g(3) (-m(4) * t94 - m(5) * t77 - t127 * (t25 * pkin(4) + qJ(5) * t24 + t77) - t126 * t52 - t131 * t25 - t132 * t24 + t133 * t51) * g(3) + (-m(4) * t96 - m(5) * t79 - t127 * (t8 * pkin(4) + t7 * qJ(5) + t79) - t131 * t8 - t132 * t7 + t126 * t39 - t133 * t37) * g(2) + (-m(4) * t95 - m(5) * t78 - t127 * (t10 * pkin(4) + t9 * qJ(5) + t78) + t126 * t42 - t132 * t9 - t131 * t10 - t133 * t40) * g(1) (-m(5) * t89 - m(6) * (t89 + t122) - m(7) * (-t36 + t122) + t118 * t44 + t119 * t43) * g(3) + (-m(5) * t91 - m(6) * (t91 + t124) - m(7) * (-t17 + t124) + t118 * t21 + t119 * t20) * g(2) + (-m(5) * t90 - m(6) * (t90 + t123) - m(7) * (-t18 + t123) + t118 * t23 + t119 * t22) * g(1), t121 * (-g(1) * t22 - g(2) * t20 - g(3) * t43) t127 * (-g(1) * t3 - g(2) * t1 - g(3) * t11) -g(1) * ((t3 * t63 - t4 * t60) * mrSges(7,1) + (-t3 * t60 - t4 * t63) * mrSges(7,2)) - g(2) * ((t1 * t63 - t2 * t60) * mrSges(7,1) + (-t1 * t60 - t2 * t63) * mrSges(7,2)) - g(3) * ((t11 * t63 - t12 * t60) * mrSges(7,1) + (-t11 * t60 - t12 * t63) * mrSges(7,2))];
taug  = t5(:);
