% Calculate Gravitation load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:53
% EndTime: 2019-03-08 20:58:55
% DurationCPUTime: 0.98s
% Computational Cost: add. (540->90), mult. (909->128), div. (0->0), fcn. (1023->14), ass. (0->44)
t29 = pkin(12) + qJ(6);
t25 = sin(t29);
t27 = cos(t29);
t31 = sin(pkin(12));
t33 = cos(pkin(12));
t80 = mrSges(5,1) + m(7) * (pkin(5) * t33 + pkin(4)) + t27 * mrSges(7,1) - t25 * mrSges(7,2) + m(6) * pkin(4) + t33 * mrSges(6,1) - t31 * mrSges(6,2);
t77 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3);
t30 = qJ(3) + pkin(11);
t26 = sin(t30);
t28 = cos(t30);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t79 = m(4) * pkin(2) + t38 * mrSges(4,1) - t36 * mrSges(4,2) - t77 * t26 + t80 * t28 + mrSges(3,1);
t88 = m(6) + m(7);
t81 = m(5) + t88;
t37 = sin(qJ(2));
t39 = cos(qJ(2));
t66 = sin(pkin(10));
t68 = cos(pkin(6));
t50 = t68 * t66;
t67 = cos(pkin(10));
t14 = -t37 * t50 + t67 * t39;
t32 = sin(pkin(6));
t60 = t32 * t66;
t91 = -t14 * t36 + t38 * t60;
t70 = t32 * t37;
t90 = -t36 * t70 + t68 * t38;
t78 = -m(4) * pkin(8) - t25 * mrSges(7,1) - t33 * mrSges(6,2) - t27 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t31;
t51 = t68 * t67;
t12 = t37 * t51 + t66 * t39;
t61 = t32 * t67;
t83 = -t12 * t36 - t38 * t61;
t69 = t32 * t39;
t34 = -qJ(4) - pkin(8);
t24 = pkin(3) * t38 + pkin(2);
t13 = t67 * t37 + t39 * t50;
t11 = t66 * t37 - t39 * t51;
t8 = t68 * t26 + t28 * t70;
t7 = t26 * t70 - t68 * t28;
t4 = t14 * t28 + t26 * t60;
t3 = t14 * t26 - t28 * t60;
t2 = t12 * t28 - t26 * t61;
t1 = t12 * t26 + t28 * t61;
t5 = [(-m(2) - m(3) - m(4) - t81) * g(3) (-t81 * (-t11 * t24 - t12 * t34) + t78 * t12 + t79 * t11) * g(2) + (-t81 * (-t13 * t24 - t14 * t34) + t78 * t14 + t79 * t13) * g(1) + (-t81 * t24 * t69 + (-t79 * t39 + (t34 * t81 + t78) * t37) * t32) * g(3) (-t90 * mrSges(4,1) - (-t68 * t36 - t38 * t70) * mrSges(4,2) + t77 * t8 + t80 * t7) * g(3) + (-(-t12 * t38 + t36 * t61) * mrSges(4,2) - mrSges(4,1) * t83 + t77 * t2 + t80 * t1) * g(2) + (-t91 * mrSges(4,1) - (-t14 * t38 - t36 * t60) * mrSges(4,2) + t77 * t4 + t80 * t3) * g(1) + (-g(1) * t91 - t83 * g(2) - g(3) * t90) * t81 * pkin(3), t81 * (-g(1) * t13 - g(2) * t11 + g(3) * t69) t88 * (-g(1) * t3 - g(2) * t1 - g(3) * t7) -g(1) * ((t13 * t27 - t25 * t4) * mrSges(7,1) + (-t13 * t25 - t27 * t4) * mrSges(7,2)) - g(2) * ((t11 * t27 - t2 * t25) * mrSges(7,1) + (-t11 * t25 - t2 * t27) * mrSges(7,2)) - g(3) * ((-t8 * t25 - t27 * t69) * mrSges(7,1) + (t25 * t69 - t8 * t27) * mrSges(7,2))];
taug  = t5(:);
