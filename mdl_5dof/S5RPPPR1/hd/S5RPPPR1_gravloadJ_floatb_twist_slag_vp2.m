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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:05
% EndTime: 2020-01-03 11:20:07
% DurationCPUTime: 0.30s
% Computational Cost: add. (177->55), mult. (173->67), div. (0->0), fcn. (148->10), ass. (0->26)
t36 = m(5) + m(6);
t31 = m(4) + t36;
t17 = sin(pkin(9));
t18 = sin(pkin(8));
t19 = cos(pkin(9));
t20 = cos(pkin(8));
t38 = -mrSges(3,1) + (-mrSges(4,1) - m(5) * pkin(3) - t19 * mrSges(5,1) + t17 * mrSges(5,2) - m(6) * (pkin(4) * t19 + pkin(3))) * t20 + (mrSges(4,2) - m(5) * qJ(4) - mrSges(5,3) - m(6) * (pkin(6) + qJ(4)) - mrSges(6,3)) * t18;
t37 = -t19 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t17;
t22 = sin(qJ(1));
t13 = t22 * pkin(1);
t23 = cos(qJ(1));
t14 = t23 * pkin(1);
t16 = qJ(1) + pkin(7);
t10 = sin(t16);
t34 = t10 * pkin(2) + t13;
t33 = t10 * t20;
t12 = cos(t16);
t32 = t12 * t20;
t15 = pkin(9) + qJ(5);
t11 = cos(t15);
t9 = sin(t15);
t4 = t10 * t9 + t11 * t32;
t3 = -t10 * t11 + t9 * t32;
t2 = t11 * t33 - t12 * t9;
t1 = -t11 * t12 - t9 * t33;
t5 = [(-m(3) * t13 - m(5) * t34 - t22 * mrSges(2,1) - t2 * mrSges(6,1) - mrSges(2,2) * t23 - t1 * mrSges(6,2) + (-m(4) - m(6)) * (-qJ(3) * t12 + t34) + (m(5) * qJ(3) - t37) * t12 + t38 * t10) * g(3) + (-m(3) * t14 - mrSges(2,1) * t23 - t4 * mrSges(6,1) + t22 * mrSges(2,2) + t3 * mrSges(6,2) - t31 * (t12 * pkin(2) + t10 * qJ(3) + t14) + t37 * t10 + t38 * t12) * g(2), (-m(3) - t31) * g(1), (g(2) * t12 + g(3) * t10) * t31, (g(1) * t20 + t18 * (-g(2) * t10 + g(3) * t12)) * t36, -g(2) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(3) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(1) * (-mrSges(6,1) * t9 - mrSges(6,2) * t11) * t18];
taug = t5(:);
