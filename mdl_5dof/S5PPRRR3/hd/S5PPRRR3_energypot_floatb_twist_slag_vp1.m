% Calculate potential energy for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:17
% EndTime: 2019-12-05 15:16:17
% DurationCPUTime: 0.47s
% Computational Cost: add. (168->95), mult. (254->114), div. (0->0), fcn. (266->10), ass. (0->35)
t45 = pkin(7) + pkin(6) + rSges(6,3);
t44 = pkin(6) + rSges(5,3);
t18 = sin(pkin(9));
t19 = sin(pkin(8));
t43 = t18 * t19;
t21 = cos(pkin(8));
t42 = t18 * t21;
t22 = sin(qJ(4));
t41 = t18 * t22;
t24 = cos(qJ(4));
t40 = t18 * t24;
t20 = cos(pkin(9));
t39 = t19 * t20;
t23 = sin(qJ(3));
t38 = t19 * t23;
t25 = cos(qJ(3));
t37 = t19 * t25;
t36 = t21 * t23;
t35 = t21 * t25;
t34 = t19 * pkin(1) + r_base(2);
t33 = qJ(1) + r_base(3);
t32 = t21 * pkin(1) + t19 * qJ(2) + r_base(1);
t31 = t18 * pkin(2) + t33;
t30 = t21 * t20 * pkin(2) + pkin(5) * t42 + t32;
t17 = qJ(4) + qJ(5);
t12 = sin(t17);
t13 = cos(t17);
t29 = t13 * rSges(6,1) - t12 * rSges(6,2) + t24 * pkin(4) + pkin(3);
t28 = pkin(2) * t39 + pkin(5) * t43 - t21 * qJ(2) + t34;
t27 = t12 * rSges(6,1) + t13 * rSges(6,2) + t22 * pkin(4);
t4 = t20 * t35 + t38;
t3 = t20 * t36 - t37;
t2 = t20 * t37 - t36;
t1 = t20 * t38 + t35;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t19 * rSges(2,2) + r_base(1)) + g(2) * (t19 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t19 * rSges(3,3) + t32) + g(2) * (rSges(3,1) * t39 - rSges(3,2) * t43 + t34) + g(3) * (t18 * rSges(3,1) + t20 * rSges(3,2) + t33) + (g(1) * (rSges(3,1) * t20 - rSges(3,2) * t18) + g(2) * (-rSges(3,3) - qJ(2))) * t21) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + rSges(4,3) * t42 + t30) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + rSges(4,3) * t43 + t28) + g(3) * ((-rSges(4,3) - pkin(5)) * t20 + (rSges(4,1) * t25 - rSges(4,2) * t23) * t18 + t31)) - m(5) * (g(1) * (t4 * pkin(3) + (t21 * t41 + t4 * t24) * rSges(5,1) + (t21 * t40 - t4 * t22) * rSges(5,2) + t44 * t3 + t30) + g(2) * (t2 * pkin(3) + (t19 * t41 + t2 * t24) * rSges(5,1) + (t19 * t40 - t2 * t22) * rSges(5,2) + t44 * t1 + t28) + g(3) * ((-t22 * rSges(5,1) - t24 * rSges(5,2) - pkin(5)) * t20 + (t44 * t23 + (t24 * rSges(5,1) - t22 * rSges(5,2) + pkin(3)) * t25) * t18 + t31)) - m(6) * (g(1) * (t45 * t3 + t30) + g(2) * (t45 * t1 + t28) + (g(1) * t21 + g(2) * t19) * t27 * t18 + (g(1) * t4 + g(2) * t2) * t29 + (t31 + (-pkin(5) - t27) * t20 + (t45 * t23 + t29 * t25) * t18) * g(3));
U = t5;
