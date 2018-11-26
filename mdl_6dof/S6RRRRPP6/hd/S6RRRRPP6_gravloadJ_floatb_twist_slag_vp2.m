% Calculate Gravitation load on the joints for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:07:56
% EndTime: 2018-11-23 18:07:58
% DurationCPUTime: 1.24s
% Computational Cost: add. (574->146), mult. (796->162), div. (0->0), fcn. (782->8), ass. (0->74)
t135 = -mrSges(3,2) + mrSges(5,3) + mrSges(6,1);
t132 = -mrSges(7,2) - mrSges(6,3);
t121 = mrSges(5,2) + t132;
t133 = mrSges(5,1) - mrSges(6,2);
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t120 = t50 * mrSges(3,1) + t135 * t47;
t134 = t47 * mrSges(4,3) + mrSges(2,1) + t120;
t130 = mrSges(7,3) + t133;
t45 = qJ(3) + qJ(4);
t40 = sin(t45);
t41 = cos(t45);
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t129 = m(4) * pkin(2) + t49 * mrSges(4,1) - t46 * mrSges(4,2) - t121 * t40 + t133 * t41;
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t128 = g(1) * t51 + g(2) * t48;
t91 = t51 * t50;
t23 = -t46 * t91 + t48 * t49;
t126 = -m(3) - m(4);
t125 = m(6) + m(7);
t124 = mrSges(2,2) - mrSges(3,3);
t88 = qJ(5) * t40;
t92 = t50 * t41;
t123 = pkin(4) * t92 + t50 * t88;
t102 = t41 * t47;
t103 = t40 * t47;
t26 = qJ(5) * t102;
t122 = -m(7) * t26 - mrSges(6,2) * t103 + t132 * t102;
t19 = t40 * t91 - t48 * t41;
t20 = t40 * t48 + t41 * t91;
t119 = t121 * t20 + t130 * t19;
t93 = t48 * t50;
t17 = t40 * t93 + t41 * t51;
t18 = -t51 * t40 + t48 * t92;
t118 = t121 * t18 + t130 * t17;
t117 = m(7) * qJ(6) + mrSges(7,3);
t115 = t117 + t133;
t111 = pkin(3) * t46;
t108 = g(3) * t47;
t101 = t46 * t48;
t100 = t46 * t51;
t52 = -pkin(9) - pkin(8);
t95 = t47 * t52;
t39 = pkin(3) * t49 + pkin(2);
t28 = t50 * t39;
t89 = t51 * pkin(1) + t48 * pkin(7);
t87 = qJ(6) * t17;
t86 = qJ(6) * t19;
t84 = t51 * t95;
t43 = t51 * pkin(7);
t83 = pkin(3) * t100 + t48 * t95 + t43;
t82 = m(4) * pkin(8) + mrSges(4,3);
t79 = t28 - t95;
t78 = -t17 * pkin(4) + qJ(5) * t18;
t77 = -t19 * pkin(4) + qJ(5) * t20;
t76 = -t39 - t88;
t75 = pkin(3) * t101 + t39 * t91 + t89;
t74 = m(7) * (pkin(5) - t52) + mrSges(7,1);
t73 = pkin(2) * t50 + pkin(8) * t47;
t72 = m(7) * (-pkin(4) - qJ(6)) - mrSges(7,3);
t68 = -mrSges(5,1) * t40 - mrSges(5,2) * t41;
t66 = t23 * pkin(3);
t21 = t46 * t93 + t49 * t51;
t65 = t74 * t47;
t64 = t72 * t40;
t61 = t21 * pkin(3);
t60 = t20 * pkin(4) + t19 * qJ(5) + t75;
t58 = t66 + t77;
t57 = -t61 + t78;
t24 = t49 * t91 + t101;
t22 = -t49 * t93 + t100;
t1 = [(-t24 * mrSges(4,1) - t23 * mrSges(4,2) - m(5) * (t75 - t84) - m(6) * (t60 - t84) - m(7) * t60 + t126 * t89 + t124 * t48 - t115 * t20 + t121 * t19 + (-m(4) * t73 - t134 - t65) * t51) * g(2) + (-m(5) * t83 - t22 * mrSges(4,1) - t21 * mrSges(4,2) - t125 * (-t18 * pkin(4) - qJ(5) * t17 + t83) + t124 * t51 + t126 * t43 + t115 * t18 - t121 * t17 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t73) + (m(7) * pkin(5) + mrSges(7,1)) * t47 + (-m(5) - t125) * (-pkin(1) - t28) + t134) * t48) * g(1) (-m(5) * t79 - m(6) * (t79 + t123) - m(7) * (t28 + t123) - t65 - t120) * g(3) + ((-t117 * t41 - t129) * g(3) + t128 * (-t74 - t82 + (m(5) + m(6)) * t52 - t135)) * t50 + (-t82 * g(3) + t128 * (mrSges(3,1) + m(5) * t39 - m(6) * (-pkin(4) * t41 + t76) - m(7) * t76 - t41 * t72 + t129)) * t47 (m(5) * t111 + mrSges(4,1) * t46 + mrSges(4,2) * t49 - t68) * t108 + (-m(6) * (t26 + (-pkin(4) * t40 - t111) * t47) - (-m(7) * t111 + t64) * t47 + t122) * g(3) + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + m(5) * t61 - m(6) * t57 - m(7) * (t57 - t87) + t118) * g(2) + (-t23 * mrSges(4,1) + t24 * mrSges(4,2) - m(5) * t66 - m(6) * t58 - m(7) * (t58 - t86) + t119) * g(1), -t68 * t108 + (-m(6) * (-pkin(4) * t103 + t26) - t47 * t64 + t122) * g(3) + (-m(6) * t78 - m(7) * (t78 - t87) + t118) * g(2) + (-m(6) * t77 - m(7) * (t77 - t86) + t119) * g(1), t125 * (-g(1) * t19 - g(2) * t17 - g(3) * t103) (-g(1) * t20 - g(2) * t18 - g(3) * t102) * m(7)];
taug  = t1(:);
