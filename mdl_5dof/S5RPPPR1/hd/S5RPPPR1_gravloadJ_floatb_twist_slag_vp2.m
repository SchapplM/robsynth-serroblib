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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:28:45
% EndTime: 2019-12-05 17:28:47
% DurationCPUTime: 0.32s
% Computational Cost: add. (177->55), mult. (173->65), div. (0->0), fcn. (148->10), ass. (0->25)
t33 = m(5) + m(6);
t36 = m(4) + t33;
t13 = sin(pkin(9));
t14 = sin(pkin(8));
t15 = cos(pkin(9));
t16 = cos(pkin(8));
t35 = mrSges(3,1) + (mrSges(4,1) + m(5) * pkin(3) + t15 * mrSges(5,1) - t13 * mrSges(5,2) + m(6) * (pkin(4) * t15 + pkin(3))) * t16 + t36 * pkin(2) + (-mrSges(4,2) + m(5) * qJ(4) + mrSges(5,3) - m(6) * (-pkin(6) - qJ(4)) + mrSges(6,3)) * t14;
t34 = -t15 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t13;
t18 = sin(qJ(1));
t32 = pkin(1) * t18;
t19 = cos(qJ(1));
t31 = pkin(1) * t19;
t12 = qJ(1) + pkin(7);
t8 = sin(t12);
t29 = t16 * t8;
t10 = cos(t12);
t28 = t10 * t16;
t11 = pkin(9) + qJ(5);
t9 = cos(t11);
t7 = sin(t11);
t4 = -t28 * t9 - t7 * t8;
t3 = t28 * t7 - t8 * t9;
t2 = -t10 * t7 + t29 * t9;
t1 = t10 * t9 + t29 * t7;
t5 = [(m(3) * t32 + t18 * mrSges(2,1) + t2 * mrSges(6,1) + mrSges(2,2) * t19 - t1 * mrSges(6,2) - t36 * (t10 * qJ(3) - t32) + t34 * t10 + t35 * t8) * g(3) + (mrSges(2,1) * t19 - t4 * mrSges(6,1) - t18 * mrSges(2,2) - t3 * mrSges(6,2) + (m(3) + m(5)) * t31 + (-m(4) - m(6)) * (-t8 * qJ(3) - t31) + (m(5) * qJ(3) - t34) * t8 + t35 * t10) * g(2), (-m(3) - t36) * g(1), -(g(2) * t10 + g(3) * t8) * t36, (g(1) * t16 + t14 * (g(2) * t8 - g(3) * t10)) * t33, -g(2) * (mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(1) * (-mrSges(6,1) * t7 - mrSges(6,2) * t9) * t14];
taug = t5(:);
