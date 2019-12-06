% Calculate potential energy for
% S5PPRPR3
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:40
% EndTime: 2019-12-05 15:04:41
% DurationCPUTime: 0.55s
% Computational Cost: add. (180->91), mult. (216->105), div. (0->0), fcn. (216->10), ass. (0->37)
t52 = rSges(4,3) + pkin(5);
t51 = pkin(6) + rSges(6,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t50 = -t23 * rSges(6,1) - t25 * rSges(6,2);
t19 = sin(pkin(7));
t21 = cos(pkin(7));
t49 = g(1) * t21 + g(2) * t19;
t18 = sin(pkin(8));
t45 = rSges(3,2) * t18;
t20 = cos(pkin(8));
t44 = t19 * t20;
t24 = sin(qJ(3));
t43 = t19 * t24;
t26 = cos(qJ(3));
t42 = t19 * t26;
t41 = t21 * t20;
t40 = t21 * t24;
t39 = t21 * t26;
t35 = t19 * pkin(1) + r_base(2);
t34 = qJ(1) + r_base(3);
t33 = t21 * pkin(1) + t19 * qJ(2) + r_base(1);
t11 = t26 * pkin(3) + pkin(2);
t22 = -qJ(4) - pkin(5);
t32 = t18 * t11 + t20 * t22 + t34;
t31 = pkin(3) * t43 + t11 * t41 + t33;
t30 = -t21 * qJ(2) + t35;
t29 = t25 * rSges(6,1) - t23 * rSges(6,2) + pkin(4);
t27 = -pkin(3) * t40 + t11 * t44 + t30;
t17 = qJ(3) + pkin(9);
t13 = cos(t17);
t12 = sin(t17);
t4 = t19 * t12 + t13 * t41;
t3 = t12 * t41 - t19 * t13;
t2 = -t21 * t12 + t13 * t44;
t1 = t12 * t44 + t21 * t13;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t19 * rSges(2,2) + r_base(1)) + g(2) * (t19 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t19 * rSges(3,3) + t33) + g(2) * (rSges(3,1) * t44 - t19 * t45 + t35) + g(3) * (t18 * rSges(3,1) + t20 * rSges(3,2) + t34) + (g(1) * (rSges(3,1) * t20 - t45) + g(2) * (-rSges(3,3) - qJ(2))) * t21) - m(4) * (g(1) * (pkin(2) * t41 + (t20 * t39 + t43) * rSges(4,1) + (-t20 * t40 + t42) * rSges(4,2) + t33) + g(2) * (pkin(2) * t44 + (t20 * t42 - t40) * rSges(4,1) + (-t20 * t43 - t39) * rSges(4,2) + t30) + g(3) * (-t52 * t20 + t34) + (g(3) * (rSges(4,1) * t26 - rSges(4,2) * t24 + pkin(2)) + t49 * t52) * t18) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t31) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t27) + g(3) * (-t20 * rSges(5,3) + t32) + (g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12) + t49 * (rSges(5,3) - t22)) * t18) - m(6) * (g(3) * (t50 * t20 + t32) + (g(3) * (t51 * t12 + t29 * t13) + t49 * (-t22 - t50)) * t18 + (t51 * t1 + t29 * t2 + t27) * g(2) + (t29 * t4 + t51 * t3 + t31) * g(1));
U = t5;
