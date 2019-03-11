% Calculate potential energy for
% S6RRPRRP11
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:18
% EndTime: 2019-03-09 12:45:18
% DurationCPUTime: 0.51s
% Computational Cost: add. (197->113), mult. (237->127), div. (0->0), fcn. (221->8), ass. (0->43)
t57 = rSges(5,3) + pkin(8);
t27 = -pkin(9) - pkin(8);
t56 = rSges(6,3) - t27;
t55 = rSges(7,3) + qJ(6) - t27;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t54 = g(1) * t26 + g(2) * t23;
t21 = sin(qJ(4));
t51 = t21 * pkin(4);
t24 = cos(qJ(4));
t11 = t24 * pkin(4) + pkin(3);
t22 = sin(qJ(2));
t49 = t22 * t23;
t20 = qJ(4) + qJ(5);
t12 = sin(t20);
t48 = t23 * t12;
t13 = cos(t20);
t47 = t23 * t13;
t46 = t23 * t21;
t45 = t23 * t24;
t25 = cos(qJ(2));
t44 = t23 * t25;
t43 = t26 * t12;
t42 = t26 * t13;
t41 = t26 * t21;
t40 = t26 * t24;
t37 = qJ(3) * t22;
t36 = pkin(6) + r_base(3);
t35 = t23 * pkin(1) + r_base(2);
t34 = t22 * t46;
t33 = t22 * t41;
t32 = t22 * pkin(2) + t36;
t31 = t26 * pkin(1) + t23 * pkin(7) + r_base(1);
t30 = pkin(2) * t44 + t23 * t37 + t35;
t29 = t31 + (pkin(2) * t25 + t37) * t26;
t28 = -t26 * pkin(7) + t30;
t6 = pkin(5) * t12 + t51;
t5 = pkin(5) * t13 + t11;
t4 = t22 * t48 - t42;
t3 = t22 * t47 + t43;
t2 = t22 * t43 + t47;
t1 = t22 * t42 - t48;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t23 * rSges(3,3) + t31) + g(2) * (rSges(3,1) * t44 - rSges(3,2) * t49 + t35) + g(3) * (t22 * rSges(3,1) + t25 * rSges(3,2) + t36) + (g(1) * (rSges(3,1) * t25 - rSges(3,2) * t22) + g(2) * (-rSges(3,3) - pkin(7))) * t26) - m(4) * (g(1) * (t23 * rSges(4,1) + t29) + g(2) * (-rSges(4,2) * t44 + rSges(4,3) * t49 + t30) + g(3) * (-t22 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t25 + t32) + (g(1) * (-rSges(4,2) * t25 + rSges(4,3) * t22) + g(2) * (-rSges(4,1) - pkin(7))) * t26) - m(5) * (g(1) * (t23 * pkin(3) + (t33 + t45) * rSges(5,1) + (t22 * t40 - t46) * rSges(5,2) + t29) + g(2) * (-t26 * pkin(3) + (t34 - t40) * rSges(5,1) + (t22 * t45 + t41) * rSges(5,2) + t28) + g(3) * (t57 * t22 + t32) + (g(3) * (-rSges(5,1) * t21 - rSges(5,2) * t24 - qJ(3)) + t54 * t57) * t25) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + pkin(4) * t33 + t23 * t11 + t29) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t34 - t26 * t11 + t28) + g(3) * (t56 * t22 + t32) + (g(3) * (-rSges(6,1) * t12 - rSges(6,2) * t13 - qJ(3) - t51) + t54 * t56) * t25) - m(7) * (g(1) * (t26 * t22 * t6 + t2 * rSges(7,1) + t1 * rSges(7,2) + t23 * t5 + t29) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) - t26 * t5 + t6 * t49 + t28) + g(3) * (t55 * t22 + t32) + (g(3) * (-rSges(7,1) * t12 - rSges(7,2) * t13 - qJ(3) - t6) + t54 * t55) * t25);
U  = t7;
