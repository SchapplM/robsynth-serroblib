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
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:44
% EndTime: 2019-12-05 14:57:45
% DurationCPUTime: 0.53s
% Computational Cost: add. (180->91), mult. (216->105), div. (0->0), fcn. (216->10), ass. (0->35)
t50 = pkin(6) + rSges(6,3);
t25 = sin(qJ(5));
t26 = cos(qJ(5));
t49 = -t25 * rSges(6,1) - t26 * rSges(6,2);
t20 = sin(pkin(7));
t23 = cos(pkin(7));
t48 = g(1) * t23 + g(2) * t20;
t47 = rSges(4,3) + qJ(3);
t19 = sin(pkin(8));
t44 = rSges(3,2) * t19;
t18 = sin(pkin(9));
t43 = t20 * t18;
t22 = cos(pkin(8));
t42 = t20 * t22;
t41 = t23 * t18;
t40 = t23 * t22;
t35 = t20 * pkin(1) + r_base(2);
t34 = qJ(1) + r_base(3);
t33 = t23 * pkin(1) + t20 * qJ(2) + r_base(1);
t21 = cos(pkin(9));
t11 = pkin(3) * t21 + pkin(2);
t24 = -pkin(5) - qJ(3);
t32 = t19 * t11 + t22 * t24 + t34;
t31 = pkin(3) * t43 + t11 * t40 + t33;
t30 = -t23 * qJ(2) + t35;
t29 = rSges(6,1) * t26 - rSges(6,2) * t25 + pkin(4);
t27 = -pkin(3) * t41 + t11 * t42 + t30;
t17 = pkin(9) + qJ(4);
t13 = cos(t17);
t12 = sin(t17);
t4 = t12 * t20 + t13 * t40;
t3 = t12 * t40 - t20 * t13;
t2 = -t12 * t23 + t13 * t42;
t1 = t12 * t42 + t13 * t23;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t23 - rSges(2,2) * t20 + r_base(1)) + g(2) * (rSges(2,1) * t20 + rSges(2,2) * t23 + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t20 * rSges(3,3) + t33) + g(2) * (rSges(3,1) * t42 - t20 * t44 + t35) + g(3) * (rSges(3,1) * t19 + rSges(3,2) * t22 + t34) + (g(1) * (rSges(3,1) * t22 - t44) + g(2) * (-rSges(3,3) - qJ(2))) * t23) - m(4) * (g(1) * (pkin(2) * t40 + (t21 * t40 + t43) * rSges(4,1) + (-t18 * t40 + t20 * t21) * rSges(4,2) + t33) + g(2) * (pkin(2) * t42 + (t21 * t42 - t41) * rSges(4,1) + (-t18 * t42 - t21 * t23) * rSges(4,2) + t30) + g(3) * (-t47 * t22 + t34) + (g(3) * (rSges(4,1) * t21 - rSges(4,2) * t18 + pkin(2)) + t48 * t47) * t19) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t31) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t27) + g(3) * (-t22 * rSges(5,3) + t32) + (g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12) + t48 * (rSges(5,3) - t24)) * t19) - m(6) * (g(3) * (t49 * t22 + t32) + (g(3) * (t12 * t50 + t29 * t13) + t48 * (-t24 - t49)) * t19 + (t1 * t50 + t29 * t2 + t27) * g(2) + (t29 * t4 + t3 * t50 + t31) * g(1));
U = t5;
