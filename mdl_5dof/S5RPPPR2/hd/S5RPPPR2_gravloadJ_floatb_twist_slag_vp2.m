% Calculate Gravitation load on the joints for
% S5RPPPR2
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:57
% EndTime: 2019-12-05 17:31:00
% DurationCPUTime: 0.50s
% Computational Cost: add. (166->70), mult. (381->99), div. (0->0), fcn. (411->10), ass. (0->38)
t53 = mrSges(2,2) - mrSges(3,3);
t52 = -mrSges(4,2) + mrSges(5,3);
t19 = sin(pkin(7));
t22 = cos(pkin(7));
t33 = -pkin(2) * t22 - qJ(3) * t19 - pkin(1);
t51 = m(3) * pkin(1) - m(4) * t33 + t22 * mrSges(3,1) + mrSges(2,1) + (-mrSges(3,2) + mrSges(4,3)) * t19;
t18 = sin(pkin(8));
t21 = cos(pkin(8));
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t38 = t26 * t22;
t13 = t18 * t38 - t24 * t21;
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t14 = t24 * t18 + t21 * t38;
t17 = sin(pkin(9));
t20 = cos(pkin(9));
t39 = t26 * t19;
t5 = t14 * t20 + t17 * t39;
t50 = -t13 * t25 + t5 * t23;
t49 = -t13 * t23 - t5 * t25;
t48 = -pkin(6) * m(6) + mrSges(5,2) - mrSges(6,3);
t47 = m(5) + m(6);
t37 = t24 * qJ(2);
t46 = -t14 * pkin(3) - t13 * qJ(4) + t33 * t26 - t37;
t40 = t24 * t22;
t11 = t18 * t40 + t26 * t21;
t12 = t26 * t18 - t21 * t40;
t16 = t26 * qJ(2);
t45 = t12 * pkin(3) - t11 * qJ(4) + t33 * t24 + t16;
t42 = t18 * t19;
t41 = t24 * t19;
t36 = m(4) + t47;
t10 = t19 * t21 * t20 - t22 * t17;
t4 = t12 * t20 - t17 * t41;
t2 = -t11 * t23 + t4 * t25;
t1 = -t11 * t25 - t4 * t23;
t3 = [(-t12 * mrSges(4,1) - t45 * m(5) - t4 * mrSges(5,1) - t2 * mrSges(6,1) - t1 * mrSges(6,2) - (t4 * pkin(4) + t45) * m(6) + t48 * (t12 * t17 + t20 * t41) + t53 * t26 + (-m(3) - m(4)) * t16 + t52 * t11 + t51 * t24) * g(3) + (m(4) * t37 + t14 * mrSges(4,1) - t46 * m(5) + t5 * mrSges(5,1) - t49 * mrSges(6,1) - t50 * mrSges(6,2) - (-pkin(4) * t5 + t46) * m(6) + t48 * (-t14 * t17 + t20 * t39) + (m(3) * qJ(2) - t53) * t24 + t52 * t13 + t51 * t26) * g(2), (g(2) * t26 + g(3) * t24) * (-m(3) - t36), (g(1) * t22 + (g(2) * t24 - g(3) * t26) * t19) * t36, t47 * (-g(1) * t42 + g(2) * t11 - g(3) * t13), -g(1) * ((-t10 * t23 + t25 * t42) * mrSges(6,1) + (-t10 * t25 - t23 * t42) * mrSges(6,2)) - g(2) * (t1 * mrSges(6,1) - t2 * mrSges(6,2)) - g(3) * (-t50 * mrSges(6,1) + t49 * mrSges(6,2))];
taug = t3(:);
