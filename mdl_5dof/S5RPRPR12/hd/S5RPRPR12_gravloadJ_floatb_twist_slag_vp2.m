% Calculate Gravitation load on the joints for
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:26
% EndTime: 2019-12-31 18:29:27
% DurationCPUTime: 0.36s
% Computational Cost: add. (235->74), mult. (271->81), div. (0->0), fcn. (233->10), ass. (0->38)
t20 = sin(qJ(1));
t21 = cos(qJ(1));
t51 = g(1) * t21 + g(2) * t20;
t46 = m(5) + m(6);
t13 = pkin(8) + qJ(3);
t11 = cos(t13);
t14 = sin(pkin(9));
t16 = cos(pkin(9));
t26 = m(5) * pkin(3) + t16 * mrSges(5,1) - t14 * mrSges(5,2);
t50 = t26 * t11;
t17 = cos(pkin(8));
t9 = sin(t13);
t29 = t11 * mrSges(4,1) - t9 * mrSges(4,2);
t48 = mrSges(2,1) + m(3) * pkin(1) + t17 * mrSges(3,1) - sin(pkin(8)) * mrSges(3,2) + t29 + t9 * mrSges(6,3);
t19 = -pkin(6) - qJ(2);
t47 = -m(3) * qJ(2) + m(5) * t19 - t14 * mrSges(5,1) - t16 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t45 = pkin(4) * t14;
t12 = pkin(9) + qJ(5);
t8 = sin(t12);
t42 = t20 * t8;
t41 = t21 * t8;
t10 = cos(t12);
t39 = t20 * t10;
t38 = t21 * t10;
t18 = -pkin(7) - qJ(4);
t36 = -m(6) * t18 + mrSges(6,3);
t33 = m(5) * qJ(4) + mrSges(5,3);
t6 = t16 * pkin(4) + pkin(3);
t31 = t11 * t6 - t9 * t18;
t28 = m(6) * t6 + t10 * mrSges(6,1) - t8 * mrSges(6,2);
t22 = t33 * t9 + t50;
t7 = t17 * pkin(2) + pkin(1);
t5 = t21 * t7;
t4 = t11 * t38 + t42;
t3 = -t11 * t41 + t39;
t2 = -t11 * t39 + t41;
t1 = t11 * t42 + t38;
t15 = [(-m(5) * t5 - t4 * mrSges(6,1) - t3 * mrSges(6,2) + (-m(4) - m(6)) * (-t20 * t19 + t5) + (-m(6) * t45 + t47) * t20 + (-m(6) * t31 - t22 - t48) * t21) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (m(4) * t19 - m(6) * (-t19 + t45) + t47) * t21 + (m(4) * t7 - m(5) * (-t9 * qJ(4) - t7) + t9 * mrSges(5,3) + t50 - m(6) * (-t31 - t7) + t48) * t20) * g(1), (-g(1) * t20 + g(2) * t21) * (m(3) + m(4) + t46), (-t22 - t29) * g(3) + (-t36 * g(3) + t51 * (mrSges(4,1) + t26 + t28)) * t9 + (-t28 * g(3) + t51 * (mrSges(4,2) - t33 - t36)) * t11, (g(3) * t11 - t9 * t51) * t46, -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t8 - mrSges(6,2) * t10) * t9];
taug = t15(:);
