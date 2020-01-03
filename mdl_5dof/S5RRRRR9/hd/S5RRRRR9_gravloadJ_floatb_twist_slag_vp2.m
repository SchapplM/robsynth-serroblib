% Calculate Gravitation load on the joints for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:44
% EndTime: 2019-12-31 22:27:46
% DurationCPUTime: 0.66s
% Computational Cost: add. (369->110), mult. (469->132), div. (0->0), fcn. (439->10), ass. (0->64)
t95 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t36 = qJ(3) + qJ(4);
t29 = cos(t36);
t40 = cos(qJ(3));
t32 = t40 * pkin(3);
t23 = pkin(4) * t29 + t32;
t21 = pkin(2) + t23;
t30 = qJ(5) + t36;
t25 = sin(t30);
t26 = cos(t30);
t27 = t32 + pkin(2);
t28 = sin(t36);
t37 = sin(qJ(3));
t94 = -m(4) * pkin(2) - m(5) * t27 - m(6) * t21 - t40 * mrSges(4,1) - t29 * mrSges(5,1) - t26 * mrSges(6,1) + t37 * mrSges(4,2) + t28 * mrSges(5,2) + t25 * mrSges(6,2);
t43 = -pkin(8) - pkin(7);
t35 = -pkin(9) + t43;
t93 = -m(4) * pkin(7) + m(5) * t43 + m(6) * t35 - t95;
t83 = m(5) * pkin(3);
t75 = t37 * pkin(3);
t79 = pkin(4) * t28;
t22 = t75 + t79;
t92 = m(6) * t22;
t91 = mrSges(4,1) + t83;
t90 = mrSges(2,2) - mrSges(3,3);
t52 = -mrSges(6,1) * t25 - mrSges(6,2) * t26;
t89 = mrSges(5,1) * t28 + mrSges(5,2) * t29 - t52;
t39 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t65 = t42 * t28;
t15 = t39 * t29 - t41 * t65;
t64 = t42 * t29;
t16 = t39 * t28 + t41 * t64;
t67 = t42 * t25;
t7 = t39 * t26 - t41 * t67;
t66 = t42 * t26;
t8 = t39 * t25 + t41 * t66;
t80 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t88 = -t15 * mrSges(5,1) + t16 * mrSges(5,2) - t80;
t69 = t39 * t41;
t13 = t28 * t69 + t64;
t14 = -t29 * t69 + t65;
t5 = t25 * t69 + t66;
t6 = -t26 * t69 + t67;
t81 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t87 = t13 * mrSges(5,1) - t14 * mrSges(5,2) - t81;
t86 = -m(3) - m(4) - m(5) - m(6);
t38 = sin(qJ(2));
t55 = t41 * mrSges(3,1) - t38 * mrSges(3,2);
t84 = t95 * t38 + mrSges(2,1) + t55;
t82 = m(6) * pkin(4);
t76 = g(3) * t38;
t70 = t39 * t37;
t68 = t42 * t22;
t63 = t42 * t37;
t62 = t42 * t40;
t56 = t41 * pkin(2) + t38 * pkin(7);
t51 = t41 * t21 - t38 * t35;
t50 = t41 * t27 - t38 * t43;
t19 = t39 * t40 - t41 * t63;
t17 = t37 * t69 + t62;
t20 = t41 * t62 + t70;
t18 = -t40 * t69 + t63;
t1 = [(-t70 * t83 - t20 * mrSges(4,1) - t16 * mrSges(5,1) - t8 * mrSges(6,1) - t19 * mrSges(4,2) - t15 * mrSges(5,2) - t7 * mrSges(6,2) + t86 * (t42 * pkin(1) + t39 * pkin(6)) + (t90 - t92) * t39 + (-m(4) * t56 - m(5) * t50 - m(6) * t51 - t84) * t42) * g(2) + (-t63 * t83 - m(6) * t68 - t18 * mrSges(4,1) - t14 * mrSges(5,1) - t6 * mrSges(6,1) - t17 * mrSges(4,2) - t13 * mrSges(5,2) - t5 * mrSges(6,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - t56) - m(5) * (-pkin(1) - t50) - m(6) * (-pkin(1) - t51) + t84) * t39 + (t86 * pkin(6) + t90) * t42) * g(1), (t93 * t38 + t94 * t41 - t55) * g(3) + (g(1) * t42 + g(2) * t39) * ((mrSges(3,2) + t93) * t41 + (mrSges(3,1) - t94) * t38), (m(5) * t75 + mrSges(4,1) * t37 + mrSges(4,2) * t40 + t89 + t92) * t76 + (-t18 * mrSges(4,2) - m(6) * (-t22 * t69 - t42 * t23) + t91 * t17 + t87) * g(2) + (t20 * mrSges(4,2) - m(6) * (t39 * t23 - t41 * t68) - t91 * t19 + t88) * g(1), (m(6) * t79 + t89) * t76 + (t13 * t82 + t87) * g(2) + (-t15 * t82 + t88) * g(1), -g(1) * t80 - g(2) * t81 - t52 * t76];
taug = t1(:);
