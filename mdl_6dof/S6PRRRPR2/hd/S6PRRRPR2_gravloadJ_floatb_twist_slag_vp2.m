% Calculate Gravitation load on the joints for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:22:51
% EndTime: 2018-11-23 15:22:52
% DurationCPUTime: 1.03s
% Computational Cost: add. (1299->119), mult. (1314->153), div. (0->0), fcn. (1278->18), ass. (0->64)
t65 = pkin(12) + qJ(6);
t61 = sin(t65);
t62 = cos(t65);
t67 = sin(pkin(12));
t69 = cos(pkin(12));
t154 = t69 * mrSges(6,1) + t62 * mrSges(7,1) - t67 * mrSges(6,2) - t61 * mrSges(7,2) + mrSges(5,1);
t150 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t59 = pkin(5) * t69 + pkin(4);
t66 = qJ(3) + qJ(4);
t63 = sin(t66);
t64 = cos(t66);
t70 = -pkin(10) - qJ(5);
t71 = sin(qJ(3));
t73 = cos(qJ(3));
t134 = m(4) * pkin(2) + t73 * mrSges(4,1) - t71 * mrSges(4,2) + mrSges(3,1) + (m(6) * pkin(4) + m(7) * t59 + t154) * t64 + (m(6) * qJ(5) - m(7) * t70 - t150) * t63;
t109 = sin(pkin(6));
t68 = sin(pkin(11));
t103 = t68 * t109;
t110 = cos(pkin(11));
t107 = pkin(6) + qJ(2);
t92 = sin(t107) / 0.2e1;
t108 = pkin(6) - qJ(2);
t97 = sin(t108);
t49 = t92 - t97 / 0.2e1;
t74 = cos(qJ(2));
t86 = t110 * t74 - t68 * t49;
t148 = t73 * t103 - t86 * t71;
t111 = cos(pkin(6));
t93 = cos(t107) / 0.2e1;
t98 = cos(t108);
t50 = t93 - t98 / 0.2e1;
t147 = t111 * t73 + t50 * t71;
t87 = t110 * t49 + t68 * t74;
t88 = t110 * t109;
t19 = t63 * t87 + t64 * t88;
t20 = -t63 * t88 + t64 * t87;
t146 = t150 * t20 + t154 * t19;
t21 = -t103 * t64 + t63 * t86;
t22 = t103 * t63 + t64 * t86;
t145 = t150 * t22 + t154 * t21;
t34 = -t111 * t64 - t50 * t63;
t35 = t111 * t63 - t50 * t64;
t144 = t150 * t35 + t154 * t34;
t143 = m(6) + m(7);
t139 = -m(5) - t143;
t133 = -m(4) * pkin(8) - t61 * mrSges(7,1) - t69 * mrSges(6,2) - t62 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t67;
t132 = -t19 * t59 - t20 * t70;
t131 = -t21 * t59 - t22 * t70;
t118 = -t34 * t59 - t35 * t70;
t101 = -t19 * pkin(4) + t20 * qJ(5);
t100 = -t21 * pkin(4) + qJ(5) * t22;
t99 = -t34 * pkin(4) + qJ(5) * t35;
t96 = t148 * pkin(3);
t95 = t147 * pkin(3);
t80 = -t71 * t87 - t73 * t88;
t79 = t98 / 0.2e1 + t93;
t78 = t80 * pkin(3);
t75 = -pkin(9) - pkin(8);
t72 = sin(qJ(2));
t60 = pkin(3) * t73 + pkin(2);
t48 = t92 + t97 / 0.2e1;
t41 = t110 * t72 + t68 * t79;
t38 = -t110 * t79 + t68 * t72;
t1 = [(-m(2) - m(3) - m(4) + t139) * g(3) (t139 * (t48 * t60 + t50 * t75) - t133 * t50 - t134 * t48) * g(3) + (t139 * (-t38 * t60 - t75 * t87) + t133 * t87 + t134 * t38) * g(2) + (t139 * (-t41 * t60 - t75 * t86) + t133 * t86 + t134 * t41) * g(1) (-t147 * mrSges(4,1) - (-t111 * t71 + t50 * t73) * mrSges(4,2) - m(5) * t95 - m(6) * (t95 + t99) - m(7) * (t95 + t118) + t144) * g(3) + (-t80 * mrSges(4,1) - (t71 * t88 - t73 * t87) * mrSges(4,2) - m(5) * t78 - m(6) * (t101 + t78) - m(7) * (t78 + t132) + t146) * g(2) + (-t148 * mrSges(4,1) - (-t103 * t71 - t73 * t86) * mrSges(4,2) - m(5) * t96 - m(6) * (t100 + t96) - m(7) * (t96 + t131) + t145) * g(1) (-m(6) * t99 - m(7) * t118 + t144) * g(3) + (-m(6) * t101 - m(7) * t132 + t146) * g(2) + (-m(6) * t100 - m(7) * t131 + t145) * g(1), t143 * (-g(1) * t21 - g(2) * t19 - g(3) * t34) -g(1) * ((-t22 * t61 + t41 * t62) * mrSges(7,1) + (-t22 * t62 - t41 * t61) * mrSges(7,2)) - g(2) * ((-t20 * t61 + t38 * t62) * mrSges(7,1) + (-t20 * t62 - t38 * t61) * mrSges(7,2)) - g(3) * ((-t35 * t61 - t48 * t62) * mrSges(7,1) + (-t35 * t62 + t48 * t61) * mrSges(7,2))];
taug  = t1(:);
