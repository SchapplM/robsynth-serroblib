% Calculate potential energy for
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:42
% EndTime: 2019-12-05 17:14:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (245->111), mult. (371->142), div. (0->0), fcn. (422->12), ass. (0->45)
t35 = sin(qJ(3));
t59 = pkin(3) * t35;
t58 = pkin(7) + rSges(4,3);
t57 = pkin(9) + rSges(6,3);
t30 = sin(pkin(10));
t31 = sin(pkin(5));
t56 = t30 * t31;
t32 = cos(pkin(10));
t55 = t31 * t32;
t54 = t31 * t35;
t36 = sin(qJ(2));
t53 = t31 * t36;
t38 = cos(qJ(3));
t52 = t31 * t38;
t39 = cos(qJ(2));
t51 = t31 * t39;
t33 = cos(pkin(5));
t50 = t33 * t36;
t49 = t33 * t39;
t48 = t30 * pkin(1) + r_base(2);
t47 = t30 * t54;
t46 = qJ(1) + r_base(3);
t45 = t32 * pkin(1) + pkin(6) * t56 + r_base(1);
t44 = t33 * pkin(6) + t46;
t13 = t30 * t49 + t32 * t36;
t14 = -t30 * t50 + t32 * t39;
t23 = t38 * pkin(3) + pkin(2);
t40 = -pkin(8) - pkin(7);
t43 = pkin(3) * t47 - t13 * t40 + t14 * t23 + t45;
t42 = t23 * t53 + t33 * t59 + t40 * t51 + t44;
t11 = t30 * t36 - t32 * t49;
t12 = t30 * t39 + t32 * t50;
t41 = t12 * t23 + (-pkin(6) - t59) * t55 - t11 * t40 + t48;
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t29 = qJ(3) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t33 * t24 + t25 * t53;
t7 = t24 * t53 - t33 * t25;
t4 = t14 * t25 + t24 * t56;
t3 = t14 * t24 - t25 * t56;
t2 = t12 * t25 - t24 * t55;
t1 = t12 * t24 + t25 * t55;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t32 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t32 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t46)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t45) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t48) + g(3) * (t33 * rSges(3,3) + t44) + (g(1) * rSges(3,3) * t30 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t39) + g(2) * (-rSges(3,3) - pkin(6)) * t32) * t31) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t38 + t47) * rSges(4,1) + (-t14 * t35 + t30 * t52) * rSges(4,2) + t58 * t13 + t45) + g(2) * (t12 * pkin(2) - pkin(6) * t55 + (t12 * t38 - t32 * t54) * rSges(4,1) + (-t12 * t35 - t32 * t52) * rSges(4,2) + t58 * t11 + t48) + g(3) * ((t35 * rSges(4,1) + t38 * rSges(4,2)) * t33 + (-t58 * t39 + (t38 * rSges(4,1) - t35 * rSges(4,2) + pkin(2)) * t36) * t31 + t44)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t43) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t41) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t51 + t42)) - m(6) * (g(1) * (t4 * pkin(4) + (t13 * t34 + t4 * t37) * rSges(6,1) + (t13 * t37 - t4 * t34) * rSges(6,2) + t57 * t3 + t43) + g(2) * (t2 * pkin(4) + (t11 * t34 + t2 * t37) * rSges(6,1) + (t11 * t37 - t2 * t34) * rSges(6,2) + t57 * t1 + t41) + g(3) * (t8 * pkin(4) + (-t34 * t51 + t8 * t37) * rSges(6,1) + (-t8 * t34 - t37 * t51) * rSges(6,2) + t57 * t7 + t42));
U = t5;
