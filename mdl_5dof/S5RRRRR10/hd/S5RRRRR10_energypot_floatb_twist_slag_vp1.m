% Calculate potential energy for
% S5RRRRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR10_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:53
% EndTime: 2019-12-31 22:31:54
% DurationCPUTime: 0.38s
% Computational Cost: add. (245->111), mult. (371->138), div. (0->0), fcn. (422->12), ass. (0->45)
t33 = sin(qJ(3));
t59 = pkin(3) * t33;
t58 = pkin(8) + rSges(4,3);
t57 = pkin(10) + rSges(6,3);
t30 = sin(pkin(5));
t34 = sin(qJ(2));
t56 = t30 * t34;
t35 = sin(qJ(1));
t55 = t30 * t35;
t38 = cos(qJ(2));
t54 = t30 * t38;
t39 = cos(qJ(1));
t53 = t30 * t39;
t52 = t35 * t34;
t51 = t35 * t38;
t50 = t39 * t34;
t49 = t39 * t38;
t48 = pkin(6) + r_base(3);
t47 = t35 * pkin(1) + r_base(2);
t46 = t33 * t55;
t31 = cos(pkin(5));
t45 = t31 * pkin(7) + t48;
t44 = t39 * pkin(1) + pkin(7) * t55 + r_base(1);
t37 = cos(qJ(3));
t23 = t37 * pkin(3) + pkin(2);
t40 = -pkin(9) - pkin(8);
t43 = t23 * t56 + t31 * t59 + t40 * t54 + t45;
t13 = t31 * t51 + t50;
t14 = -t31 * t52 + t49;
t42 = pkin(3) * t46 - t13 * t40 + t14 * t23 + t44;
t11 = -t31 * t49 + t52;
t12 = t31 * t50 + t51;
t41 = t12 * t23 + (-pkin(7) - t59) * t53 - t11 * t40 + t47;
t36 = cos(qJ(5));
t32 = sin(qJ(5));
t29 = qJ(3) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t31 * t24 + t25 * t56;
t7 = t24 * t56 - t31 * t25;
t4 = t14 * t25 + t24 * t55;
t3 = t14 * t24 - t25 * t55;
t2 = t12 * t25 - t24 * t53;
t1 = t12 * t24 + t25 * t53;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t39 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t39 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t48)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t44) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t47) + g(3) * (t31 * rSges(3,3) + t45) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t34 + rSges(3,2) * t38) + g(2) * (-rSges(3,3) - pkin(7)) * t39) * t30) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t37 + t46) * rSges(4,1) + (-t14 * t33 + t37 * t55) * rSges(4,2) + t58 * t13 + t44) + g(2) * (t12 * pkin(2) - pkin(7) * t53 + (t12 * t37 - t33 * t53) * rSges(4,1) + (-t12 * t33 - t37 * t53) * rSges(4,2) + t58 * t11 + t47) + g(3) * ((t33 * rSges(4,1) + t37 * rSges(4,2)) * t31 + (-t58 * t38 + (t37 * rSges(4,1) - t33 * rSges(4,2) + pkin(2)) * t34) * t30 + t45)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t42) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t41) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t54 + t43)) - m(6) * (g(1) * (t4 * pkin(4) + (t13 * t32 + t4 * t36) * rSges(6,1) + (t13 * t36 - t4 * t32) * rSges(6,2) + t57 * t3 + t42) + g(2) * (t2 * pkin(4) + (t11 * t32 + t2 * t36) * rSges(6,1) + (t11 * t36 - t2 * t32) * rSges(6,2) + t57 * t1 + t41) + g(3) * (t8 * pkin(4) + (-t32 * t54 + t8 * t36) * rSges(6,1) + (-t8 * t32 - t36 * t54) * rSges(6,2) + t57 * t7 + t43));
U = t5;
