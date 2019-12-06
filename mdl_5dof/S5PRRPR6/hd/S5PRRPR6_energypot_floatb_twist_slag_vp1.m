% Calculate potential energy for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:29:59
% EndTime: 2019-12-05 16:30:00
% DurationCPUTime: 0.35s
% Computational Cost: add. (237->103), mult. (435->127), div. (0->0), fcn. (511->12), ass. (0->43)
t27 = sin(pkin(5));
t54 = pkin(6) * t27;
t53 = rSges(4,3) + pkin(7);
t32 = sin(qJ(3));
t52 = t27 * t32;
t33 = sin(qJ(2));
t51 = t27 * t33;
t34 = cos(qJ(3));
t50 = t27 * t34;
t35 = cos(qJ(2));
t49 = t27 * t35;
t30 = cos(pkin(5));
t48 = t30 * t33;
t47 = t30 * t35;
t46 = pkin(8) + qJ(4) + rSges(6,3);
t45 = qJ(4) + rSges(5,3);
t26 = sin(pkin(9));
t44 = t26 * pkin(1) + r_base(2);
t43 = qJ(1) + r_base(3);
t29 = cos(pkin(9));
t42 = t29 * pkin(1) + t26 * t54 + r_base(1);
t41 = t30 * pkin(6) + t43;
t10 = -t26 * t48 + t29 * t35;
t40 = t10 * pkin(2) + t42;
t39 = pkin(2) * t51 + t41;
t24 = pkin(10) + qJ(5);
t19 = sin(t24);
t20 = cos(t24);
t28 = cos(pkin(10));
t38 = t20 * rSges(6,1) - t19 * rSges(6,2) + t28 * pkin(4) + pkin(3);
t8 = t26 * t35 + t29 * t48;
t37 = t8 * pkin(2) - t29 * t54 + t44;
t25 = sin(pkin(10));
t36 = t19 * rSges(6,1) + t20 * rSges(6,2) + t25 * pkin(4) + pkin(7);
t12 = t30 * t32 + t33 * t50;
t11 = -t30 * t34 + t32 * t51;
t9 = t26 * t47 + t29 * t33;
t7 = t26 * t33 - t29 * t47;
t4 = t10 * t34 + t26 * t52;
t3 = t10 * t32 - t26 * t50;
t2 = -t29 * t52 + t8 * t34;
t1 = t29 * t50 + t8 * t32;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t43)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t42) + g(2) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t44) + g(3) * (t30 * rSges(3,3) + t41) + (g(1) * rSges(3,3) * t26 + g(3) * (rSges(3,1) * t33 + rSges(3,2) * t35) + g(2) * (-rSges(3,3) - pkin(6)) * t29) * t27) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t53 * t9 + t40) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t53 * t7 + t37) + g(3) * (t12 * rSges(4,1) - t11 * rSges(4,2) - t53 * t49 + t39)) - m(5) * (g(1) * (t4 * pkin(3) + t9 * pkin(7) + (t9 * t25 + t4 * t28) * rSges(5,1) + (-t4 * t25 + t9 * t28) * rSges(5,2) + t45 * t3 + t40) + g(2) * (t2 * pkin(3) + t7 * pkin(7) + (t2 * t28 + t7 * t25) * rSges(5,1) + (-t2 * t25 + t7 * t28) * rSges(5,2) + t45 * t1 + t37) + g(3) * (t12 * pkin(3) - pkin(7) * t49 + (t12 * t28 - t25 * t49) * rSges(5,1) + (-t12 * t25 - t28 * t49) * rSges(5,2) + t45 * t11 + t39)) - m(6) * (g(1) * (t46 * t3 + t36 * t9 + t38 * t4 + t40) + g(2) * (t46 * t1 + t38 * t2 + t36 * t7 + t37) + g(3) * (t46 * t11 + t38 * t12 - t36 * t49 + t39));
U = t5;
