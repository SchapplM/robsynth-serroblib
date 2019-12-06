% Calculate Gravitation load on the joints for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:22
% EndTime: 2019-12-05 15:16:24
% DurationCPUTime: 0.29s
% Computational Cost: add. (160->51), mult. (312->79), div. (0->0), fcn. (328->10), ass. (0->30)
t49 = -m(6) * pkin(4) - mrSges(5,1);
t15 = qJ(4) + qJ(5);
t13 = sin(t15);
t14 = cos(t15);
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t48 = mrSges(4,1) + m(6) * (t22 * pkin(4) + pkin(3)) + t14 * mrSges(6,1) - t13 * mrSges(6,2) + m(5) * pkin(3) + t22 * mrSges(5,1) - t20 * mrSges(5,2);
t47 = mrSges(4,2) + m(6) * (-pkin(7) - pkin(6)) - mrSges(6,3) - m(5) * pkin(6) - mrSges(5,3);
t16 = sin(pkin(9));
t17 = sin(pkin(8));
t41 = t16 * t17;
t18 = cos(pkin(9));
t19 = cos(pkin(8));
t21 = sin(qJ(3));
t34 = t19 * t21;
t23 = cos(qJ(3));
t35 = t17 * t23;
t8 = t18 * t35 - t34;
t45 = (-t8 * t13 + t14 * t41) * mrSges(6,1) + (-t13 * t41 - t8 * t14) * mrSges(6,2);
t33 = t19 * t23;
t36 = t17 * t21;
t10 = t18 * t33 + t36;
t40 = t16 * t19;
t44 = (-t10 * t13 + t14 * t40) * mrSges(6,1) + (-t10 * t14 - t13 * t40) * mrSges(6,2);
t37 = t16 * t23;
t43 = (-t13 * t37 - t18 * t14) * mrSges(6,1) + (t18 * t13 - t14 * t37) * mrSges(6,2);
t39 = t16 * t20;
t38 = t16 * t22;
t32 = m(3) + m(4) + m(5) + m(6);
t1 = [(-m(2) - t32) * g(3), (-g(1) * t17 + g(2) * t19) * t32, (t48 * t21 + t47 * t23) * g(3) * t16 + (t47 * t8 - t48 * (-t18 * t36 - t33)) * g(2) + (-t48 * (-t18 * t34 + t35) + t47 * t10) * g(1), (-(t18 * t20 - t22 * t37) * mrSges(5,2) - t43 + t49 * (-t18 * t22 - t20 * t37)) * g(3) + (-(-t17 * t39 - t8 * t22) * mrSges(5,2) - t45 + t49 * (t17 * t38 - t8 * t20)) * g(2) + (-(-t10 * t22 - t19 * t39) * mrSges(5,2) - t44 + t49 * (-t10 * t20 + t19 * t38)) * g(1), -g(1) * t44 - g(2) * t45 - g(3) * t43];
taug = t1(:);
