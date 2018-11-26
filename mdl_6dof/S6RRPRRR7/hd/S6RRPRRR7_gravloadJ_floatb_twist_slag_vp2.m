% Calculate Gravitation load on the joints for
% S6RRPRRR7
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:25:02
% EndTime: 2018-11-23 17:25:03
% DurationCPUTime: 1.20s
% Computational Cost: add. (431->118), mult. (891->147), div. (0->0), fcn. (940->10), ass. (0->63)
t113 = mrSges(5,2) - m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) - mrSges(6,3) - mrSges(7,3);
t92 = -m(5) - m(6) - m(7);
t38 = cos(qJ(5));
t25 = pkin(5) * t38 + pkin(4);
t34 = qJ(5) + qJ(6);
t26 = sin(t34);
t27 = cos(t34);
t35 = sin(qJ(5));
t112 = m(6) * pkin(4) + m(7) * t25 + t38 * mrSges(6,1) + t27 * mrSges(7,1) - t35 * mrSges(6,2) - t26 * mrSges(7,2) + mrSges(5,1);
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t80 = sin(qJ(4));
t81 = cos(qJ(4));
t13 = t36 * t80 + t39 * t81;
t37 = sin(qJ(1));
t10 = t13 * t37;
t14 = t36 * t81 - t39 * t80;
t9 = t14 * t37;
t111 = t113 * t10 - t112 * t9;
t110 = t112 * t13 + t113 * t14;
t40 = cos(qJ(1));
t65 = t40 * t80;
t66 = t40 * t81;
t11 = -t36 * t65 - t39 * t66;
t12 = -t36 * t66 + t39 * t65;
t109 = -t113 * t11 + t112 * t12;
t88 = -pkin(2) - pkin(3);
t106 = t36 * t88;
t105 = (mrSges(3,1) + mrSges(4,1)) * t39 + (-mrSges(3,2) + mrSges(4,3)) * t36;
t104 = -mrSges(2,1) - t105;
t100 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t97 = m(7) * pkin(5) + mrSges(6,1);
t32 = t40 * pkin(7);
t28 = t36 * qJ(3);
t62 = -pkin(1) - t28;
t96 = (t88 * t39 + t62) * t37 - t40 * pkin(8) + t32;
t93 = g(1) * t40 + g(2) * t37;
t51 = t10 * t26 - t27 * t40;
t52 = -t10 * t27 - t26 * t40;
t87 = -t51 * mrSges(7,1) + t52 * mrSges(7,2);
t5 = t11 * t26 - t27 * t37;
t6 = -t11 * t27 - t26 * t37;
t86 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t85 = pkin(5) * t35;
t82 = g(3) * t14;
t31 = t39 * pkin(2);
t79 = t35 * t40;
t76 = t39 * t40;
t75 = t31 + t28;
t74 = t40 * pkin(1) + t37 * pkin(7);
t73 = qJ(3) * t39;
t61 = pkin(2) * t76 + t40 * t28 + t74;
t60 = pkin(3) * t76 + t61;
t53 = -mrSges(7,1) * t26 - mrSges(7,2) * t27;
t50 = -t10 * t38 - t79;
t49 = t10 * t35 - t38 * t40;
t7 = t11 * t35 - t37 * t38;
t48 = -t37 * pkin(8) + t60;
t43 = t39 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t36;
t20 = t40 * t73;
t19 = t37 * t73;
t8 = -t11 * t38 - t35 * t37;
t1 = [(-m(3) * t74 - m(4) * t61 - m(5) * t48 + t11 * mrSges(5,1) - m(6) * (-pkin(4) * t11 + t48) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (-t11 * t25 + t60) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t104 * t40 + t113 * t12 + (-m(7) * (-pkin(8) - t85) + t100) * t37) * g(2) + (-t96 * m(5) + t10 * mrSges(5,1) - t50 * mrSges(6,1) - t49 * mrSges(6,2) - (-pkin(4) * t10 + t96) * m(6) - t52 * mrSges(7,1) - t51 * mrSges(7,2) - (-pkin(5) * t79 - t10 * t25 + t96) * m(7) + (-m(3) - m(4)) * t32 + t113 * t9 + (m(3) * pkin(1) - m(4) * (t62 - t31) - t104) * t37 + t100 * t40) * g(1), t93 * (mrSges(3,1) * t36 + mrSges(3,2) * t39) + (-m(4) * t19 - t43 * t37 - t111 + t92 * (t106 * t37 + t19)) * g(2) + (-m(4) * t20 - t43 * t40 - t109 + t92 * (t106 * t40 + t20)) * g(1) + (-m(4) * t75 + t92 * (t39 * pkin(3) + t75) - t105 - t110) * g(3) (t39 * g(3) - t93 * t36) * (m(4) - t92) t109 * g(1) + t111 * g(2) + t110 * g(3) (m(7) * t85 + t35 * mrSges(6,1) + t38 * mrSges(6,2) - t53) * t82 + (-t50 * mrSges(6,2) + t97 * t49 - t87) * g(2) + (t8 * mrSges(6,2) - t97 * t7 - t86) * g(1), -g(1) * t86 - g(2) * t87 - t53 * t82];
taug  = t1(:);
