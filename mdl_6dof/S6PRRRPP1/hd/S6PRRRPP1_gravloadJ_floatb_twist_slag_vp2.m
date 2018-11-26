% Calculate Gravitation load on the joints for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2018-11-23 15:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:20:08
% EndTime: 2018-11-23 15:20:09
% DurationCPUTime: 1.08s
% Computational Cost: add. (1378->106), mult. (1575->142), div. (0->0), fcn. (1561->16), ass. (0->63)
t136 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t61 = qJ(4) + pkin(11);
t59 = sin(t61);
t60 = cos(t61);
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t91 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t125 = m(5) * pkin(3) + t66 * mrSges(5,1) - t63 * mrSges(5,2) + t136 * t59 + t60 * t91 + mrSges(4,1);
t127 = m(6) + m(7);
t139 = t127 * pkin(4);
t124 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t58 = pkin(4) * t66 + pkin(3);
t62 = -qJ(5) - pkin(9);
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t137 = t127 * (-t58 * t67 + t62 * t64) + t124 * t64 - mrSges(3,1) - t125 * t67;
t135 = -m(4) - t127;
t106 = cos(pkin(6));
t101 = pkin(6) + qJ(2);
t86 = cos(t101) / 0.2e1;
t102 = pkin(6) - qJ(2);
t94 = cos(t102);
t52 = t86 - t94 / 0.2e1;
t42 = t106 * t64 - t52 * t67;
t92 = sin(t101);
t84 = t92 / 0.2e1;
t93 = sin(t102);
t85 = t93 / 0.2e1;
t51 = t84 + t85;
t133 = t42 * t63 + t51 * t66;
t103 = sin(pkin(10));
t79 = t84 - t93 / 0.2e1;
t105 = cos(pkin(10));
t68 = cos(qJ(2));
t96 = t105 * t68;
t39 = -t103 * t79 + t96;
t104 = sin(pkin(6));
t80 = t104 * t103;
t20 = t39 * t67 + t64 * t80;
t65 = sin(qJ(2));
t72 = t94 / 0.2e1 + t86;
t38 = t103 * t72 + t105 * t65;
t132 = -t20 * t63 + t38 * t66;
t95 = t103 * t68;
t36 = t105 * t79 + t95;
t81 = t105 * t104;
t18 = t36 * t67 - t64 * t81;
t35 = t103 * t65 - t105 * t72;
t131 = -t18 * t63 + t35 * t66;
t129 = m(5) * pkin(8) + t66 * mrSges(5,2) + t91 * t59 - t136 * t60 - mrSges(3,2) + mrSges(4,3) + (mrSges(5,1) + t139) * t63;
t71 = t85 - t92 / 0.2e1;
t50 = t51 * pkin(2);
t41 = -t106 * t67 - t52 * t64;
t40 = t103 * t71 + t96;
t37 = -t105 * t71 + t95;
t34 = t38 * pkin(2);
t33 = t35 * pkin(2);
t19 = t39 * t64 - t67 * t80;
t17 = t36 * t64 + t67 * t81;
t9 = t42 * t59 + t51 * t60;
t3 = t20 * t59 - t38 * t60;
t1 = t18 * t59 - t35 * t60;
t2 = [(-m(2) - m(3) - m(5) + t135) * g(3) (-m(5) * t50 + t135 * (-pkin(8) * t52 + t50) + t129 * t52 + t137 * t51) * g(3) + (m(5) * t33 + t135 * (pkin(8) * t37 - t33) - t129 * t37 - t137 * t35) * g(2) + (m(5) * t34 + t135 * (pkin(8) * t40 - t34) - t129 * t40 - t137 * t38) * g(1) (-t127 * (-t41 * t58 - t42 * t62) + t124 * t42 + t125 * t41) * g(3) + (-t127 * (-t17 * t58 - t18 * t62) + t124 * t18 + t125 * t17) * g(2) + (-t127 * (-t19 * t58 - t20 * t62) + t124 * t20 + t125 * t19) * g(1) (t133 * mrSges(5,1) - (-t42 * t66 + t51 * t63) * mrSges(5,2) + t91 * t9 - t136 * (t42 * t60 - t51 * t59)) * g(3) + (-t131 * mrSges(5,1) - (-t18 * t66 - t35 * t63) * mrSges(5,2) - t136 * (t18 * t60 + t35 * t59) + t91 * t1) * g(2) + (-t132 * mrSges(5,1) - (-t20 * t66 - t38 * t63) * mrSges(5,2) - t136 * (t20 * t60 + t38 * t59) + t91 * t3) * g(1) + (-t132 * g(1) - t131 * g(2) + t133 * g(3)) * t139, t127 * (-g(1) * t19 - g(2) * t17 - g(3) * t41) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
