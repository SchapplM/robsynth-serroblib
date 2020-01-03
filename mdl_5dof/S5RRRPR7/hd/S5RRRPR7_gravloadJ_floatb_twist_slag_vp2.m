% Calculate Gravitation load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:16:09
% DurationCPUTime: 0.66s
% Computational Cost: add. (326->96), mult. (376->104), div. (0->0), fcn. (325->10), ass. (0->55)
t108 = m(5) * qJ(4) + mrSges(5,3);
t32 = qJ(2) + qJ(3);
t28 = sin(t32);
t29 = cos(t32);
t33 = sin(pkin(9));
t72 = t33 * mrSges(5,2);
t31 = pkin(9) + qJ(5);
t26 = sin(t31);
t78 = t26 * mrSges(6,2);
t107 = (-t72 - t78) * t28 + (-mrSges(6,3) - t108) * t29;
t34 = cos(pkin(9));
t23 = t34 * pkin(4) + pkin(3);
t35 = -pkin(8) - qJ(4);
t97 = t29 * t23 - t28 * t35;
t102 = m(6) * t97;
t101 = -t29 * mrSges(4,1) + (mrSges(4,2) - mrSges(6,3)) * t28;
t71 = t34 * mrSges(5,1);
t27 = cos(t31);
t77 = t27 * mrSges(6,1);
t100 = (t71 + t77) * t28;
t99 = t71 - t72;
t50 = -t23 * t28 - t29 * t35;
t86 = pkin(3) * t28;
t36 = sin(qJ(2));
t87 = pkin(2) * t36;
t95 = -m(6) * (t50 - t87) - m(5) * (-t86 - t87) + t100;
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t94 = g(1) * t39 + g(2) * t37;
t21 = t28 * mrSges(5,3);
t93 = -t21 + t101 + (-t77 + t78 - t99) * t29;
t40 = -pkin(7) - pkin(6);
t92 = -m(3) * pkin(6) + m(5) * t40 - t33 * mrSges(5,1) - t34 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t91 = t107 * t37;
t90 = t107 * t39;
t89 = m(5) * t86 - m(6) * t50 + t100;
t38 = cos(qJ(2));
t55 = t38 * mrSges(3,1) - t36 * mrSges(3,2);
t88 = -(m(5) * pkin(3) + t99) * t29 - mrSges(2,1) - m(3) * pkin(1) - t55 + t101;
t85 = pkin(4) * t33;
t30 = t38 * pkin(2);
t70 = t37 * t26;
t69 = t37 * t27;
t68 = t39 * t26;
t67 = t39 * t27;
t19 = t28 * qJ(4);
t66 = t29 * pkin(3) + t19;
t52 = mrSges(4,1) * t28 + mrSges(4,2) * t29;
t25 = t30 + pkin(1);
t14 = t39 * t25;
t4 = t29 * t67 + t70;
t3 = -t29 * t68 + t69;
t2 = -t29 * t69 + t68;
t1 = t29 * t70 + t67;
t5 = [(-m(5) * t14 - t4 * mrSges(6,1) - t3 * mrSges(6,2) + (-m(4) - m(6)) * (-t37 * t40 + t14) + (-m(6) * t85 + t92) * t37 + (-t108 * t28 - t102 + t88) * t39) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (m(4) * t40 - m(6) * (-t40 + t85) + t92) * t39 + (m(4) * t25 - m(5) * (-t25 - t19) + t21 - m(6) * (-t25 - t97) - t88) * t37) * g(1), (t95 * t37 + t91) * g(2) + (t95 * t39 + t90) * g(1) + (-t55 - m(4) * t30 - m(5) * (t30 + t66) - m(6) * (t30 + t97) + t93) * g(3) + (m(4) * t87 + mrSges(3,1) * t36 + mrSges(3,2) * t38 + t52) * t94, t94 * t52 + (t89 * t37 + t91) * g(2) + (t89 * t39 + t90) * g(1) + (-m(5) * t66 - t102 + t93) * g(3), (t29 * g(3) - t28 * t94) * (m(5) + m(6)), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t26 - mrSges(6,2) * t27) * t28];
taug = t5(:);
