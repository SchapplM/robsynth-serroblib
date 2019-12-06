% Calculate potential energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:44
% EndTime: 2019-12-05 14:57:45
% DurationCPUTime: 0.56s
% Computational Cost: add. (180->67), mult. (228->63), div. (0->0), fcn. (216->10), ass. (0->30)
t19 = sin(pkin(9));
t22 = cos(pkin(9));
t55 = -m(4) * pkin(2) - t22 * mrSges(4,1) + t19 * mrSges(4,2) - mrSges(3,1);
t26 = sin(qJ(5));
t27 = cos(qJ(5));
t54 = -m(4) * qJ(3) - t26 * mrSges(6,1) - t27 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t51 = -m(1) - m(2);
t50 = m(3) + m(4);
t49 = -m(5) - m(6);
t48 = -t19 * mrSges(4,1) - t22 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t46 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t45 = -m(6) * pkin(4) - mrSges(6,1) * t27 + mrSges(6,2) * t26 - mrSges(5,1);
t20 = sin(pkin(8));
t23 = cos(pkin(8));
t44 = t54 * t20 + t55 * t23 - mrSges(2,1);
t43 = pkin(3) * t19;
t25 = -pkin(5) - qJ(3);
t41 = t20 * t25;
t21 = sin(pkin(7));
t38 = t21 * t23;
t24 = cos(pkin(7));
t37 = t23 * t24;
t36 = t21 * pkin(1) + r_base(2);
t17 = qJ(1) + r_base(3);
t35 = t24 * pkin(1) + t21 * qJ(2) + r_base(1);
t18 = pkin(9) + qJ(4);
t13 = cos(t18);
t12 = sin(t18);
t11 = pkin(3) * t22 + pkin(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t49 * (t20 * t11 + t23 * t25 + t17) + (-m(2) - t50) * t17 - t54 * t23 + (t46 * t12 + t45 * t13 + t55) * t20) * g(3) + (-mrSges(1,2) + t51 * r_base(2) - t50 * t36 + t49 * (-t21 * t41 + t11 * t38 + (-qJ(2) - t43) * t24 + t36) + t46 * (t12 * t38 + t13 * t24) + (t50 * qJ(2) - t48) * t24 + t45 * (-t12 * t24 + t13 * t38) + t44 * t21) * g(2) + (-mrSges(1,1) + t51 * r_base(1) - t50 * t35 + t49 * (t11 * t37 + t21 * t43 - t24 * t41 + t35) + t46 * (t12 * t37 - t21 * t13) + t45 * (t12 * t21 + t13 * t37) + t48 * t21 + t44 * t24) * g(1);
U = t1;
