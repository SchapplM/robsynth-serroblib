% Calculate Gravitation load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:13
% EndTime: 2022-01-20 09:12:15
% DurationCPUTime: 0.36s
% Computational Cost: add. (177->57), mult. (173->71), div. (0->0), fcn. (148->10), ass. (0->26)
t36 = m(5) + m(6);
t30 = m(4) + t36;
t16 = sin(pkin(9));
t18 = cos(pkin(9));
t38 = -t18 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t16;
t17 = sin(pkin(8));
t19 = cos(pkin(8));
t37 = -mrSges(3,1) + (-m(5) * pkin(3) - t18 * mrSges(5,1) + t16 * mrSges(5,2) - mrSges(4,1)) * t19 + (mrSges(4,2) - mrSges(6,3)) * t17;
t21 = sin(qJ(1));
t35 = pkin(1) * t21;
t22 = cos(qJ(1));
t13 = t22 * pkin(1);
t15 = qJ(1) + pkin(7);
t10 = sin(t15);
t33 = t10 * t19;
t12 = cos(t15);
t32 = t12 * t19;
t26 = -t17 * (-pkin(6) - qJ(4)) + t19 * (pkin(4) * t18 + pkin(3));
t14 = pkin(9) + qJ(5);
t11 = cos(t14);
t9 = sin(t14);
t4 = t10 * t9 + t11 * t32;
t3 = t10 * t11 - t32 * t9;
t2 = -t11 * t33 + t12 * t9;
t1 = t11 * t12 + t33 * t9;
t5 = [(-m(3) * t13 - mrSges(2,1) * t22 - t4 * mrSges(6,1) + t21 * mrSges(2,2) - t3 * mrSges(6,2) - t30 * (t12 * pkin(2) + t10 * qJ(3) + t13) + t38 * t10 + (-(m(5) * qJ(4) + mrSges(5,3)) * t17 - m(6) * t26 + t37) * t12) * g(2) + (m(3) * t35 + t21 * mrSges(2,1) - t2 * mrSges(6,1) + mrSges(2,2) * t22 - t1 * mrSges(6,2) - t30 * (t12 * qJ(3) - t35) + t38 * t12 + (m(4) * pkin(2) - m(5) * (-qJ(4) * t17 - pkin(2)) + t17 * mrSges(5,3) - m(6) * (-pkin(2) - t26) - t37) * t10) * g(1), (-m(3) - t30) * g(3), (-g(1) * t10 + g(2) * t12) * t30, (g(3) * t19 + t17 * (-g(1) * t12 - g(2) * t10)) * t36, -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t9 - mrSges(6,2) * t11) * t17];
taug = t5(:);
