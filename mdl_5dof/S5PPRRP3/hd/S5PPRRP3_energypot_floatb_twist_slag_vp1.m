% Calculate potential energy for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:29
% EndTime: 2019-12-05 15:10:29
% DurationCPUTime: 0.28s
% Computational Cost: add. (163->91), mult. (279->108), div. (0->0), fcn. (301->8), ass. (0->43)
t57 = rSges(6,1) + pkin(4);
t56 = rSges(6,2) + pkin(6);
t55 = rSges(5,3) + pkin(6);
t26 = sin(pkin(8));
t27 = sin(pkin(7));
t54 = t26 * t27;
t29 = cos(pkin(7));
t53 = t26 * t29;
t30 = sin(qJ(4));
t52 = t26 * t30;
t31 = sin(qJ(3));
t51 = t26 * t31;
t32 = cos(qJ(4));
t50 = t26 * t32;
t33 = cos(qJ(3));
t49 = t26 * t33;
t28 = cos(pkin(8));
t48 = t27 * t28;
t47 = t27 * t31;
t46 = t27 * t33;
t45 = t29 * t31;
t44 = t29 * t33;
t43 = rSges(6,3) + qJ(5);
t42 = t27 * pkin(1) + r_base(2);
t41 = qJ(1) + r_base(3);
t40 = t29 * pkin(1) + t27 * qJ(2) + r_base(1);
t39 = t26 * pkin(2) + t41;
t38 = t29 * t28 * pkin(2) + pkin(5) * t53 + t40;
t10 = t28 * t44 + t47;
t37 = t10 * pkin(3) + t38;
t36 = pkin(2) * t48 + pkin(5) * t54 - t29 * qJ(2) + t42;
t35 = pkin(3) * t49 - t28 * pkin(5) + pkin(6) * t51 + t39;
t8 = t28 * t46 - t45;
t34 = t8 * pkin(3) + t36;
t12 = -t28 * t30 + t32 * t49;
t11 = t28 * t32 + t30 * t49;
t9 = t28 * t45 - t46;
t7 = t28 * t47 + t44;
t4 = t10 * t32 + t29 * t52;
t3 = t10 * t30 - t29 * t50;
t2 = t27 * t52 + t8 * t32;
t1 = -t27 * t50 + t8 * t30;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t27 * rSges(3,3) + t40) + g(2) * (rSges(3,1) * t48 - rSges(3,2) * t54 + t42) + g(3) * (t26 * rSges(3,1) + t28 * rSges(3,2) + t41) + (g(1) * (rSges(3,1) * t28 - rSges(3,2) * t26) + g(2) * (-rSges(3,3) - qJ(2))) * t29) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + rSges(4,3) * t53 + t38) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + rSges(4,3) * t54 + t36) + g(3) * ((-rSges(4,3) - pkin(5)) * t28 + (rSges(4,1) * t33 - rSges(4,2) * t31) * t26 + t39)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t55 * t9 + t37) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t55 * t7 + t34) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) + rSges(5,3) * t51 + t35)) - m(6) * (g(1) * (t43 * t3 + t57 * t4 + t56 * t9 + t37) + g(2) * (t43 * t1 + t57 * t2 + t56 * t7 + t34) + g(3) * (rSges(6,2) * t51 + t43 * t11 + t57 * t12 + t35));
U = t5;
