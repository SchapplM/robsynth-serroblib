% Calculate Gravitation load on the joints for
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:06
% EndTime: 2019-12-05 18:16:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (252->46), mult. (169->46), div. (0->0), fcn. (125->10), ass. (0->33)
t22 = qJ(4) + qJ(5);
t20 = cos(t22);
t12 = t20 * mrSges(6,1);
t25 = cos(qJ(4));
t23 = sin(qJ(4));
t42 = mrSges(5,2) * t23;
t49 = -t42 + m(5) * pkin(3) + m(6) * (pkin(4) * t25 + pkin(3)) + t25 * mrSges(5,1) + mrSges(4,1) + t12;
t48 = -m(5) * pkin(7) - mrSges(5,3) - mrSges(6,3) + mrSges(4,2) + m(6) * (-pkin(8) - pkin(7));
t35 = m(6) * pkin(4) + mrSges(5,1);
t45 = -mrSges(5,2) * t25 - t23 * t35;
t21 = qJ(1) + pkin(9);
t18 = qJ(3) + t21;
t13 = sin(t18);
t19 = sin(t22);
t38 = t13 * t19;
t39 = mrSges(6,2) * t20;
t44 = mrSges(6,1) * t38 + t13 * t39;
t14 = cos(t18);
t43 = g(3) * t14;
t40 = mrSges(6,2) * t19;
t37 = m(4) + m(5) + m(6);
t36 = m(3) + t37;
t34 = t12 - t40;
t33 = -mrSges(6,1) * t19 - t39;
t32 = t37 * pkin(2) + mrSges(3,1);
t31 = pkin(1) * t36 + mrSges(2,1);
t29 = (-t40 + t49) * t14 - t48 * t13;
t28 = -mrSges(6,2) * t38 + t49 * t13 + t48 * t14;
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t17 = cos(t21);
t16 = sin(t21);
t1 = [(mrSges(2,2) * t26 + mrSges(3,2) * t17 + t16 * t32 + t24 * t31 + t28) * g(3) + (-mrSges(2,2) * t24 - mrSges(3,2) * t16 + t17 * t32 + t26 * t31 + t29) * g(2), -t36 * g(1), g(2) * t29 + g(3) * t28, (t45 * t13 - t44) * g(2) + (-t25 * t35 - t34 + t42) * g(1) + (-t33 - t45) * t43, -g(1) * t34 - g(2) * t44 - t33 * t43];
taug = t1(:);
