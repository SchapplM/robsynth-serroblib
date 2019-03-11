% Calculate potential energy for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:16
% EndTime: 2019-03-09 12:50:16
% DurationCPUTime: 0.52s
% Computational Cost: add. (204->109), mult. (252->122), div. (0->0), fcn. (240->8), ass. (0->39)
t55 = rSges(7,1) + pkin(5);
t54 = rSges(5,3) + pkin(8);
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t53 = g(1) * t26 + g(2) * t23;
t52 = rSges(7,3) + qJ(6);
t22 = sin(qJ(2));
t48 = t22 * t23;
t47 = t22 * t26;
t24 = cos(qJ(4));
t46 = t23 * t24;
t25 = cos(qJ(2));
t45 = t23 * t25;
t44 = t24 * t26;
t41 = qJ(3) * t22;
t40 = pkin(6) + r_base(3);
t39 = t23 * pkin(1) + r_base(2);
t21 = sin(qJ(4));
t38 = t21 * t48;
t37 = t21 * t47;
t36 = t22 * pkin(2) + t40;
t35 = -pkin(4) * t21 - qJ(3);
t34 = t26 * pkin(1) + t23 * pkin(7) + r_base(1);
t33 = pkin(2) * t45 + t23 * t41 + t39;
t32 = t34 + (pkin(2) * t25 + t41) * t26;
t27 = -pkin(9) - pkin(8);
t31 = -t22 * t27 + t36;
t30 = -t26 * pkin(7) + t33;
t13 = pkin(4) * t24 + pkin(3);
t29 = pkin(4) * t37 + t23 * t13 + t32;
t28 = pkin(4) * t38 - t26 * t13 + t30;
t20 = qJ(4) + qJ(5);
t15 = cos(t20);
t14 = sin(t20);
t4 = t14 * t48 - t15 * t26;
t3 = t14 * t26 + t15 * t48;
t2 = t14 * t47 + t15 * t23;
t1 = t14 * t23 - t15 * t47;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t26 - rSges(2,2) * t23 + r_base(1)) + g(2) * (rSges(2,1) * t23 + rSges(2,2) * t26 + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(1) * (rSges(3,3) * t23 + t34) + g(2) * (rSges(3,1) * t45 - rSges(3,2) * t48 + t39) + g(3) * (rSges(3,1) * t22 + rSges(3,2) * t25 + t40) + (g(1) * (rSges(3,1) * t25 - rSges(3,2) * t22) + g(2) * (-rSges(3,3) - pkin(7))) * t26) - m(4) * (g(1) * (rSges(4,1) * t23 + t32) + g(2) * (-rSges(4,2) * t45 + rSges(4,3) * t48 + t33) + g(3) * (-rSges(4,2) * t22 + (-rSges(4,3) - qJ(3)) * t25 + t36) + (g(1) * (-rSges(4,2) * t25 + rSges(4,3) * t22) + g(2) * (-rSges(4,1) - pkin(7))) * t26) - m(5) * (g(1) * (t23 * pkin(3) + (t37 + t46) * rSges(5,1) + (-t21 * t23 + t22 * t44) * rSges(5,2) + t32) + g(2) * (-t26 * pkin(3) + (t38 - t44) * rSges(5,1) + (t21 * t26 + t22 * t46) * rSges(5,2) + t30) + g(3) * (t54 * t22 + t36) + (g(3) * (-rSges(5,1) * t21 - rSges(5,2) * t24 - qJ(3)) + t53 * t54) * t25) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t29) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t28) + g(3) * (t22 * rSges(6,3) + t31) + (g(3) * (-rSges(6,1) * t14 - rSges(6,2) * t15 + t35) + t53 * (rSges(6,3) - t27)) * t25) - m(7) * (g(1) * (t52 * t1 + t55 * t2 + t29) + g(2) * (-t52 * t3 + t55 * t4 + t28) + g(3) * (t22 * rSges(7,2) + t31) + (g(3) * (-t55 * t14 + t52 * t15 + t35) + t53 * (rSges(7,2) - t27)) * t25);
U  = t5;
