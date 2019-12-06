% Calculate potential energy for
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:53
% EndTime: 2019-12-05 16:34:54
% DurationCPUTime: 0.34s
% Computational Cost: add. (261->108), mult. (536->139), div. (0->0), fcn. (651->12), ass. (0->49)
t34 = sin(pkin(5));
t63 = pkin(6) * t34;
t62 = rSges(4,3) + pkin(7);
t61 = pkin(8) + rSges(6,3);
t39 = sin(qJ(3));
t60 = t34 * t39;
t40 = sin(qJ(2));
t59 = t34 * t40;
t42 = cos(qJ(3));
t58 = t34 * t42;
t43 = cos(qJ(2));
t57 = t34 * t43;
t37 = cos(pkin(5));
t56 = t37 * t40;
t55 = t37 * t43;
t54 = rSges(5,3) + qJ(4);
t33 = sin(pkin(9));
t53 = t33 * pkin(1) + r_base(2);
t52 = qJ(1) + r_base(3);
t36 = cos(pkin(9));
t51 = t36 * pkin(1) + t33 * t63 + r_base(1);
t50 = t37 * pkin(6) + t52;
t21 = -t33 * t56 + t36 * t43;
t49 = t21 * pkin(2) + t51;
t48 = pkin(2) * t59 + t50;
t19 = t33 * t43 + t36 * t56;
t47 = t19 * pkin(2) - t36 * t63 + t53;
t12 = t21 * t42 + t33 * t60;
t20 = t33 * t55 + t36 * t40;
t46 = t12 * pkin(3) + t20 * pkin(7) + t49;
t23 = t37 * t39 + t40 * t58;
t45 = t23 * pkin(3) - pkin(7) * t57 + t48;
t10 = t19 * t42 - t36 * t60;
t18 = t33 * t40 - t36 * t55;
t44 = t10 * pkin(3) + t18 * pkin(7) + t47;
t41 = cos(qJ(5));
t38 = sin(qJ(5));
t35 = cos(pkin(10));
t32 = sin(pkin(10));
t22 = -t37 * t42 + t39 * t59;
t11 = t21 * t39 - t33 * t58;
t9 = t19 * t39 + t36 * t58;
t8 = t23 * t35 - t32 * t57;
t7 = t23 * t32 + t35 * t57;
t4 = t12 * t35 + t20 * t32;
t3 = t12 * t32 - t20 * t35;
t2 = t10 * t35 + t18 * t32;
t1 = t10 * t32 - t18 * t35;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t36 * rSges(2,1) - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + t36 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t52)) - m(3) * (g(1) * (t21 * rSges(3,1) - t20 * rSges(3,2) + t51) + g(2) * (t19 * rSges(3,1) - t18 * rSges(3,2) + t53) + g(3) * (t37 * rSges(3,3) + t50) + (g(1) * rSges(3,3) * t33 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t43) + g(2) * (-rSges(3,3) - pkin(6)) * t36) * t34) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t62 * t20 + t49) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t62 * t18 + t47) + g(3) * (t23 * rSges(4,1) - t22 * rSges(4,2) - t62 * t57 + t48)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t54 * t11 + t46) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t54 * t9 + t44) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t54 * t22 + t45)) - m(6) * (g(1) * (t4 * pkin(4) + t11 * qJ(4) + (t11 * t38 + t4 * t41) * rSges(6,1) + (t11 * t41 - t4 * t38) * rSges(6,2) + t61 * t3 + t46) + g(2) * (t2 * pkin(4) + t9 * qJ(4) + (t2 * t41 + t9 * t38) * rSges(6,1) + (-t2 * t38 + t9 * t41) * rSges(6,2) + t61 * t1 + t44) + g(3) * (t8 * pkin(4) + t22 * qJ(4) + (t22 * t38 + t8 * t41) * rSges(6,1) + (t22 * t41 - t8 * t38) * rSges(6,2) + t61 * t7 + t45));
U = t5;
