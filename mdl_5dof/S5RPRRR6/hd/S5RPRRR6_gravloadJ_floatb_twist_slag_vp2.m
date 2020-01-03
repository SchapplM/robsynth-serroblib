% Calculate Gravitation load on the joints for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:00:58
% DurationCPUTime: 0.40s
% Computational Cost: add. (273->76), mult. (251->90), div. (0->0), fcn. (210->10), ass. (0->42)
t25 = qJ(3) + qJ(4);
t20 = sin(t25);
t21 = cos(t25);
t26 = sin(qJ(5));
t50 = t26 * mrSges(6,2);
t72 = -t20 * t50 + t21 * (-m(6) * pkin(8) - mrSges(6,3));
t65 = t21 * pkin(4) + t20 * pkin(8);
t71 = m(6) * t65;
t70 = -t21 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t20;
t67 = m(3) + m(4);
t66 = -m(5) - m(6);
t29 = cos(qJ(5));
t49 = t29 * mrSges(6,1);
t64 = -(t49 - t50) * t21 + t70;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t40 = t30 * mrSges(4,1) - t27 * mrSges(4,2);
t62 = m(4) * pkin(2) + mrSges(3,1) + t40 - t70;
t61 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t60 = pkin(3) * t27;
t28 = sin(qJ(1));
t57 = t28 * pkin(1);
t22 = t30 * pkin(3);
t31 = cos(qJ(1));
t23 = t31 * pkin(1);
t56 = mrSges(5,2) * t21;
t52 = t21 * t26;
t51 = t21 * t29;
t24 = qJ(1) + pkin(9);
t18 = sin(t24);
t45 = t72 * t18;
t19 = cos(t24);
t44 = t72 * t19;
t34 = m(6) * (-pkin(4) * t20 - t60) - t20 * t49;
t33 = t56 + (m(6) * pkin(4) + mrSges(5,1) + t49) * t20;
t32 = -pkin(7) - pkin(6);
t17 = t22 + pkin(2);
t4 = t18 * t26 + t19 * t51;
t3 = t18 * t29 - t19 * t52;
t2 = -t18 * t51 + t19 * t26;
t1 = t18 * t52 + t19 * t29;
t5 = [(-t31 * mrSges(2,1) - t4 * mrSges(6,1) + t28 * mrSges(2,2) - t3 * mrSges(6,2) + t66 * (t19 * t17 - t18 * t32 + t23) - t67 * t23 + t61 * t18 + (-t62 - t71) * t19) * g(2) + (t28 * mrSges(2,1) - t2 * mrSges(6,1) + t31 * mrSges(2,2) - t1 * mrSges(6,2) + t67 * t57 + t66 * (-t19 * t32 - t57) + t61 * t19 + (m(5) * t17 - m(6) * (-t17 - t65) + t62) * t18) * g(1), (t66 - t67) * g(3), -g(1) * (t34 * t19 - t44) - g(2) * (t34 * t18 - t45) + (-t40 - m(5) * t22 - m(6) * (t22 + t65) + t64) * g(3) + (m(5) * t60 + mrSges(4,1) * t27 + mrSges(5,1) * t20 + mrSges(4,2) * t30 + t56) * (g(1) * t19 + g(2) * t18), (t64 - t71) * g(3) + (t33 * t18 + t45) * g(2) + (t33 * t19 + t44) * g(1), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t26 - mrSges(6,2) * t29) * t20];
taug = t5(:);
