% Calculate Gravitation load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:07
% EndTime: 2018-11-23 15:37:08
% DurationCPUTime: 0.29s
% Computational Cost: add. (210->61), mult. (218->74), div. (0->0), fcn. (173->8), ass. (0->30)
t42 = -m(6) - m(7);
t12 = qJ(1) + pkin(9);
t10 = cos(t12);
t9 = sin(t12);
t41 = g(1) * t10 + g(2) * t9;
t31 = -m(5) + t42;
t40 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t23 = t14 * mrSges(6,1) + t17 * mrSges(6,2);
t39 = m(7) * (pkin(5) * t14 - pkin(8) * t17) + mrSges(3,1) - mrSges(4,2) + mrSges(5,3) + t23 - t17 * mrSges(7,3);
t15 = sin(qJ(1));
t37 = pkin(1) * t15;
t18 = cos(qJ(1));
t11 = t18 * pkin(1);
t13 = sin(qJ(6));
t35 = t14 * t13;
t16 = cos(qJ(6));
t34 = t14 * t16;
t30 = t10 * pkin(2) + t9 * qJ(3) + t11;
t29 = m(4) - t31;
t28 = m(7) * pkin(8) + mrSges(7,3);
t27 = t10 * qJ(3) - t37;
t26 = t10 * qJ(4) + t30;
t20 = m(7) * pkin(5) + t16 * mrSges(7,1) - t13 * mrSges(7,2);
t4 = t10 * t34 - t13 * t9;
t3 = -t10 * t35 - t16 * t9;
t2 = -t10 * t13 - t9 * t34;
t1 = -t10 * t16 + t9 * t35;
t5 = [(-m(3) * t11 - m(4) * t30 - m(5) * t26 - mrSges(2,1) * t18 - t4 * mrSges(7,1) + t15 * mrSges(2,2) - t3 * mrSges(7,2) + t42 * (-pkin(7) * t9 + t26) + t40 * t9 - t39 * t10) * g(2) + (m(3) * t37 + t15 * mrSges(2,1) - t2 * mrSges(7,1) + mrSges(2,2) * t18 - t1 * mrSges(7,2) + (-m(4) - m(5)) * t27 + t42 * (-pkin(7) * t10 + t27) + t40 * t10 + (m(4) * pkin(2) + t31 * (-pkin(2) - qJ(4)) + t39) * t9) * g(1) (-m(3) - t29) * g(3) (-g(1) * t9 + g(2) * t10) * t29, t41 * t31 (t20 * t14 - t28 * t17 + t23) * g(3) + ((-mrSges(6,1) - t20) * t17 + (mrSges(6,2) - t28) * t14) * t41, -g(1) * (mrSges(7,1) * t3 - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - g(3) * (-mrSges(7,1) * t13 - mrSges(7,2) * t16) * t17];
taug  = t5(:);
