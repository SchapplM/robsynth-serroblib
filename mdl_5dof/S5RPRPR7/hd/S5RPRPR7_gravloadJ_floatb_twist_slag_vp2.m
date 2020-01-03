% Calculate Gravitation load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:54
% EndTime: 2019-12-31 18:18:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (225->63), mult. (209->73), div. (0->0), fcn. (173->10), ass. (0->35)
t47 = m(5) + m(6);
t50 = -mrSges(5,2) + mrSges(6,3);
t17 = sin(qJ(3));
t20 = cos(qJ(3));
t13 = qJ(3) + pkin(9);
t7 = sin(t13);
t9 = cos(t13);
t49 = -t20 * mrSges(4,1) - t9 * mrSges(5,1) + t17 * mrSges(4,2) - t50 * t7;
t48 = m(3) + m(4);
t44 = m(4) * pkin(2) + mrSges(3,1) - t49;
t43 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t41 = t7 * pkin(7);
t18 = sin(qJ(1));
t38 = t18 * pkin(1);
t11 = t20 * pkin(3);
t21 = cos(qJ(1));
t12 = t21 * pkin(1);
t16 = sin(qJ(5));
t14 = qJ(1) + pkin(8);
t8 = sin(t14);
t36 = t8 * t16;
t19 = cos(qJ(5));
t35 = t8 * t19;
t10 = cos(t14);
t34 = t10 * t16;
t33 = t10 * t19;
t31 = t9 * pkin(4) + t41;
t24 = m(6) * pkin(4) + t19 * mrSges(6,1) - t16 * mrSges(6,2);
t15 = -qJ(4) - pkin(6);
t6 = t11 + pkin(2);
t4 = t9 * t33 + t36;
t3 = -t9 * t34 + t35;
t2 = -t9 * t35 + t34;
t1 = t9 * t36 + t33;
t5 = [(-t21 * mrSges(2,1) - t4 * mrSges(6,1) + t18 * mrSges(2,2) - t3 * mrSges(6,2) - t47 * (t10 * t6 - t8 * t15 + t12) - t48 * t12 + t43 * t8 + (-m(6) * t31 - t44) * t10) * g(2) + (t18 * mrSges(2,1) - t2 * mrSges(6,1) + t21 * mrSges(2,2) - t1 * mrSges(6,2) + t48 * t38 - t47 * (-t10 * t15 - t38) + t43 * t10 + (m(5) * t6 - m(6) * (-t31 - t6) + t44) * t8) * g(1), (-t47 - t48) * g(3), (-m(5) * t11 - m(6) * (t11 + t41) - t24 * t9 + t49) * g(3) + (g(1) * t10 + g(2) * t8) * (mrSges(4,2) * t20 + (-m(6) * pkin(7) - t50) * t9 + (mrSges(5,1) + t24) * t7 + (t47 * pkin(3) + mrSges(4,1)) * t17), t47 * (-g(1) * t8 + g(2) * t10), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t16 - mrSges(6,2) * t19) * t7];
taug = t5(:);
