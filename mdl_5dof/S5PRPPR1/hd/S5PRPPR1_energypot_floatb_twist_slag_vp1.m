% Calculate potential energy for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:41
% EndTime: 2019-12-05 15:21:42
% DurationCPUTime: 0.44s
% Computational Cost: add. (175->84), mult. (141->96), div. (0->0), fcn. (125->10), ass. (0->31)
t42 = rSges(6,3) + pkin(6) + qJ(4);
t14 = sin(pkin(8));
t17 = cos(pkin(8));
t41 = rSges(4,1) * t17 - rSges(4,2) * t14;
t12 = pkin(7) + qJ(2);
t6 = sin(t12);
t8 = cos(t12);
t40 = g(1) * t8 + g(2) * t6;
t39 = rSges(5,3) + qJ(4);
t13 = sin(pkin(9));
t36 = t6 * t13;
t35 = t6 * t17;
t34 = t8 * t13;
t33 = t8 * t17;
t30 = t13 * t17;
t16 = cos(pkin(9));
t29 = t16 * t17;
t15 = sin(pkin(7));
t26 = t15 * pkin(1) + r_base(2);
t18 = cos(pkin(7));
t25 = t18 * pkin(1) + r_base(1);
t24 = qJ(1) + r_base(3);
t23 = t6 * pkin(2) + t26;
t22 = pkin(5) + t24;
t21 = t8 * pkin(2) + t6 * qJ(3) + t25;
t20 = -t8 * qJ(3) + t23;
t11 = pkin(9) + qJ(5);
t7 = cos(t11);
t5 = sin(t11);
t4 = t16 * pkin(4) + pkin(3);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t18 * rSges(2,1) - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + t18 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t24)) - m(3) * (g(1) * (rSges(3,1) * t8 - rSges(3,2) * t6 + t25) + g(2) * (rSges(3,1) * t6 + rSges(3,2) * t8 + t26) + g(3) * (rSges(3,3) + t22)) - m(4) * (g(1) * (t6 * rSges(4,3) + t21) + g(2) * (t41 * t6 + t23) + g(3) * (rSges(4,1) * t14 + rSges(4,2) * t17 + t22) + (g(1) * t41 + g(2) * (-rSges(4,3) - qJ(3))) * t8) - m(5) * (g(1) * (pkin(3) * t33 + (t29 * t8 + t36) * rSges(5,1) + (t16 * t6 - t30 * t8) * rSges(5,2) + t21) + g(2) * (pkin(3) * t35 + (t29 * t6 - t34) * rSges(5,1) + (-t16 * t8 - t30 * t6) * rSges(5,2) + t20) + g(3) * (-t39 * t17 + t22) + (g(3) * (rSges(5,1) * t16 - rSges(5,2) * t13 + pkin(3)) + t40 * t39) * t14) - m(6) * (g(1) * (t4 * t33 + pkin(4) * t36 + (t33 * t7 + t5 * t6) * rSges(6,1) + (-t33 * t5 + t6 * t7) * rSges(6,2) + t21) + g(2) * (t4 * t35 - pkin(4) * t34 + (t35 * t7 - t5 * t8) * rSges(6,1) + (-t35 * t5 - t7 * t8) * rSges(6,2) + t20) + g(3) * (-t42 * t17 + t22) + (g(3) * (rSges(6,1) * t7 - rSges(6,2) * t5 + t4) + t40 * t42) * t14);
U = t1;
