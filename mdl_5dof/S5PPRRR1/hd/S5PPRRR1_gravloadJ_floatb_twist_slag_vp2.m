% Calculate Gravitation load on the joints for
% S5PPRRR1
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:36
% EndTime: 2019-12-05 15:12:37
% DurationCPUTime: 0.23s
% Computational Cost: add. (180->45), mult. (170->63), div. (0->0), fcn. (134->8), ass. (0->29)
t21 = cos(qJ(5));
t33 = t21 * mrSges(6,1);
t51 = -mrSges(5,1) - t33;
t17 = pkin(9) + qJ(3);
t16 = qJ(4) + t17;
t12 = sin(t16);
t13 = cos(t16);
t20 = sin(qJ(5));
t34 = t20 * mrSges(6,2);
t50 = -t12 * t34 + t13 * (-m(6) * pkin(7) - mrSges(6,3));
t47 = (t34 + t51) * t13 + (mrSges(5,2) - mrSges(6,3)) * t12;
t14 = sin(t17);
t45 = pkin(3) * t14;
t15 = cos(t17);
t44 = pkin(3) * t15;
t41 = mrSges(5,2) * t13;
t18 = sin(pkin(8));
t38 = t18 * t20;
t37 = t18 * t21;
t19 = cos(pkin(8));
t36 = t19 * t20;
t35 = t19 * t21;
t32 = t13 * pkin(4) + t12 * pkin(7);
t30 = m(3) + m(4) + m(5) + m(6);
t28 = t50 * t18;
t27 = t50 * t19;
t23 = m(6) * (-pkin(4) * t12 - t45) - t12 * t33;
t22 = t41 + (m(6) * pkin(4) - t51) * t12;
t1 = [(-m(2) - t30) * g(3), (-g(1) * t18 + g(2) * t19) * t30, -g(1) * (t23 * t19 - t27) - g(2) * (t23 * t18 - t28) + (-mrSges(4,1) * t15 + mrSges(4,2) * t14 - m(5) * t44 - m(6) * (t32 + t44) + t47) * g(3) + (m(5) * t45 + mrSges(4,1) * t14 + mrSges(5,1) * t12 + mrSges(4,2) * t15 + t41) * (g(1) * t19 + g(2) * t18), (-m(6) * t32 + t47) * g(3) + (t22 * t18 + t28) * g(2) + (t22 * t19 + t27) * g(1), -g(1) * ((-t13 * t36 + t37) * mrSges(6,1) + (-t13 * t35 - t38) * mrSges(6,2)) - g(2) * ((-t13 * t38 - t35) * mrSges(6,1) + (-t13 * t37 + t36) * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t20 - mrSges(6,2) * t21) * t12];
taug = t1(:);
