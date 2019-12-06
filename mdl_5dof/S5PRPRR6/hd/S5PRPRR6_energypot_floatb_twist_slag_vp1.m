% Calculate potential energy for
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:21
% EndTime: 2019-12-05 15:56:21
% DurationCPUTime: 0.38s
% Computational Cost: add. (245->111), mult. (371->140), div. (0->0), fcn. (422->12), ass. (0->43)
t30 = sin(pkin(10));
t57 = pkin(3) * t30;
t56 = pkin(8) + rSges(6,3);
t31 = sin(pkin(9));
t32 = sin(pkin(5));
t55 = t31 * t32;
t34 = cos(pkin(9));
t54 = t32 * t34;
t38 = sin(qJ(2));
t53 = t32 * t38;
t40 = cos(qJ(2));
t52 = t32 * t40;
t35 = cos(pkin(5));
t51 = t35 * t38;
t50 = t35 * t40;
t49 = qJ(3) + rSges(4,3);
t48 = t31 * pkin(1) + r_base(2);
t47 = t30 * t55;
t46 = qJ(1) + r_base(3);
t45 = t34 * pkin(1) + pkin(6) * t55 + r_base(1);
t44 = t35 * pkin(6) + t46;
t13 = t31 * t50 + t34 * t38;
t14 = -t31 * t51 + t34 * t40;
t33 = cos(pkin(10));
t23 = t33 * pkin(3) + pkin(2);
t36 = -pkin(7) - qJ(3);
t43 = pkin(3) * t47 - t13 * t36 + t14 * t23 + t45;
t42 = t23 * t53 + t35 * t57 + t36 * t52 + t44;
t11 = t31 * t38 - t34 * t50;
t12 = t31 * t40 + t34 * t51;
t41 = t12 * t23 + (-pkin(6) - t57) * t54 - t11 * t36 + t48;
t39 = cos(qJ(5));
t37 = sin(qJ(5));
t29 = pkin(10) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t35 * t24 + t25 * t53;
t7 = t24 * t53 - t35 * t25;
t4 = t14 * t25 + t24 * t55;
t3 = t14 * t24 - t25 * t55;
t2 = t12 * t25 - t24 * t54;
t1 = t12 * t24 + t25 * t54;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t31 * rSges(2,2) + r_base(1)) + g(2) * (t31 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t46)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t45) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t48) + g(3) * (t35 * rSges(3,3) + t44) + (g(1) * rSges(3,3) * t31 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t40) + g(2) * (-rSges(3,3) - pkin(6)) * t34) * t32) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t33 + t47) * rSges(4,1) + (-t14 * t30 + t33 * t55) * rSges(4,2) + t49 * t13 + t45) + g(2) * (t12 * pkin(2) - pkin(6) * t54 + (t12 * t33 - t30 * t54) * rSges(4,1) + (-t12 * t30 - t33 * t54) * rSges(4,2) + t49 * t11 + t48) + g(3) * ((t30 * rSges(4,1) + t33 * rSges(4,2)) * t35 + (-t49 * t40 + (t33 * rSges(4,1) - t30 * rSges(4,2) + pkin(2)) * t38) * t32 + t44)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t43) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t41) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t52 + t42)) - m(6) * (g(1) * (t4 * pkin(4) + (t13 * t37 + t4 * t39) * rSges(6,1) + (t13 * t39 - t4 * t37) * rSges(6,2) + t56 * t3 + t43) + g(2) * (t2 * pkin(4) + (t11 * t37 + t2 * t39) * rSges(6,1) + (t11 * t39 - t2 * t37) * rSges(6,2) + t56 * t1 + t41) + g(3) * (t8 * pkin(4) + (-t37 * t52 + t8 * t39) * rSges(6,1) + (-t8 * t37 - t39 * t52) * rSges(6,2) + t56 * t7 + t42));
U = t5;
