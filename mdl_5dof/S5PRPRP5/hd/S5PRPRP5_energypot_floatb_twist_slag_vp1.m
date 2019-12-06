% Calculate potential energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:29
% EndTime: 2019-12-05 15:37:29
% DurationCPUTime: 0.42s
% Computational Cost: add. (167->89), mult. (195->101), div. (0->0), fcn. (187->8), ass. (0->32)
t45 = rSges(6,1) + pkin(4);
t19 = sin(pkin(7));
t21 = cos(pkin(7));
t44 = g(1) * t21 + g(2) * t19;
t43 = rSges(4,3) + qJ(3);
t42 = rSges(6,3) + qJ(5);
t23 = sin(qJ(2));
t39 = rSges(3,2) * t23;
t18 = sin(pkin(8));
t38 = t19 * t18;
t24 = cos(qJ(2));
t37 = t19 * t24;
t36 = t21 * t18;
t35 = t21 * t24;
t31 = t19 * pkin(1) + r_base(2);
t30 = qJ(1) + r_base(3);
t29 = t21 * pkin(1) + t19 * pkin(5) + r_base(1);
t20 = cos(pkin(8));
t10 = t20 * pkin(3) + pkin(2);
t22 = -pkin(6) - qJ(3);
t28 = t23 * t10 + t24 * t22 + t30;
t27 = -t21 * pkin(5) + t31;
t26 = pkin(3) * t38 + t10 * t35 + t29;
t25 = -pkin(3) * t36 + t10 * t37 + t27;
t17 = pkin(8) + qJ(4);
t13 = cos(t17);
t12 = sin(t17);
t4 = t19 * t12 + t13 * t35;
t3 = t12 * t35 - t19 * t13;
t2 = -t21 * t12 + t13 * t37;
t1 = t12 * t37 + t21 * t13;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t19 * rSges(2,2) + r_base(1)) + g(2) * (t19 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t19 * rSges(3,3) + t29) + g(2) * (rSges(3,1) * t37 - t19 * t39 + t31) + g(3) * (t23 * rSges(3,1) + t24 * rSges(3,2) + t30) + (g(1) * (rSges(3,1) * t24 - t39) + g(2) * (-rSges(3,3) - pkin(5))) * t21) - m(4) * (g(1) * (pkin(2) * t35 + (t20 * t35 + t38) * rSges(4,1) + (-t18 * t35 + t19 * t20) * rSges(4,2) + t29) + g(2) * (pkin(2) * t37 + (t20 * t37 - t36) * rSges(4,1) + (-t18 * t37 - t21 * t20) * rSges(4,2) + t27) + g(3) * (-t43 * t24 + t30) + (g(3) * (rSges(4,1) * t20 - rSges(4,2) * t18 + pkin(2)) + t44 * t43) * t23) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t26) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t25) + g(3) * (-t24 * rSges(5,3) + t28) + (g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12) + t44 * (rSges(5,3) - t22)) * t23) - m(6) * (g(1) * (t42 * t3 + t45 * t4 + t26) + g(2) * (t42 * t1 + t45 * t2 + t25) + g(3) * (-t24 * rSges(6,2) + t28) + (g(3) * (t42 * t12 + t45 * t13) + t44 * (rSges(6,2) - t22)) * t23);
U = t5;
