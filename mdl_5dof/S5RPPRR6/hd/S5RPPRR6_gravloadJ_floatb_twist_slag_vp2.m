% Calculate Gravitation load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (209->60), mult. (187->70), div. (0->0), fcn. (155->10), ass. (0->33)
t43 = m(3) + m(4);
t42 = -m(5) - m(6);
t15 = cos(pkin(9));
t12 = pkin(9) + qJ(4);
t7 = sin(t12);
t9 = cos(t12);
t26 = t9 * mrSges(5,1) - t7 * mrSges(5,2);
t40 = mrSges(3,1) + m(4) * pkin(2) + t15 * mrSges(4,1) - sin(pkin(9)) * mrSges(4,2) + t26 + t7 * mrSges(6,3);
t39 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t18 = sin(qJ(1));
t37 = pkin(1) * t18;
t17 = sin(qJ(5));
t13 = qJ(1) + pkin(8);
t8 = sin(t13);
t35 = t17 * t8;
t19 = cos(qJ(5));
t34 = t19 * t8;
t20 = cos(qJ(1));
t11 = t20 * pkin(1);
t10 = cos(t13);
t32 = t10 * t17;
t31 = t10 * t19;
t30 = m(4) - t42;
t29 = m(6) * pkin(7) + mrSges(6,3);
t27 = pkin(4) * t9 + pkin(7) * t7;
t22 = m(6) * pkin(4) + t19 * mrSges(6,1) - t17 * mrSges(6,2);
t16 = -pkin(6) - qJ(3);
t6 = pkin(3) * t15 + pkin(2);
t4 = t9 * t31 + t35;
t3 = -t9 * t32 + t34;
t2 = -t9 * t34 + t32;
t1 = t9 * t35 + t31;
t5 = [(-mrSges(2,1) * t20 - t4 * mrSges(6,1) + t18 * mrSges(2,2) - t3 * mrSges(6,2) + t42 * (t10 * t6 - t16 * t8 + t11) - t43 * t11 + t39 * t8 + (-m(6) * t27 - t40) * t10) * g(2) + (t18 * mrSges(2,1) - t2 * mrSges(6,1) + mrSges(2,2) * t20 - t1 * mrSges(6,2) + t43 * t37 + t42 * (-t10 * t16 - t37) + t39 * t10 + (m(5) * t6 - m(6) * (-t27 - t6) + t40) * t8) * g(1), (-m(3) - t30) * g(3), (-g(1) * t8 + g(2) * t10) * t30, (-t22 * t9 - t29 * t7 - t26) * g(3) + ((mrSges(5,2) - t29) * t9 + (mrSges(5,1) + t22) * t7) * (g(1) * t10 + g(2) * t8), -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t17 - mrSges(6,2) * t19) * t7];
taug = t5(:);
