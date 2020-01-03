% Calculate potential energy for
% S5RRPRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR14_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR14_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:41
% EndTime: 2019-12-31 20:35:41
% DurationCPUTime: 0.38s
% Computational Cost: add. (245->111), mult. (371->138), div. (0->0), fcn. (422->12), ass. (0->45)
t30 = sin(pkin(10));
t59 = pkin(3) * t30;
t58 = pkin(9) + rSges(6,3);
t31 = sin(pkin(5));
t36 = sin(qJ(2));
t57 = t31 * t36;
t37 = sin(qJ(1));
t56 = t31 * t37;
t39 = cos(qJ(2));
t55 = t31 * t39;
t40 = cos(qJ(1));
t54 = t31 * t40;
t53 = t37 * t36;
t52 = t37 * t39;
t51 = t40 * t36;
t50 = t40 * t39;
t49 = qJ(3) + rSges(4,3);
t48 = pkin(6) + r_base(3);
t47 = t37 * pkin(1) + r_base(2);
t46 = t30 * t56;
t33 = cos(pkin(5));
t45 = t33 * pkin(7) + t48;
t44 = t40 * pkin(1) + pkin(7) * t56 + r_base(1);
t32 = cos(pkin(10));
t23 = t32 * pkin(3) + pkin(2);
t34 = -pkin(8) - qJ(3);
t43 = t23 * t57 + t33 * t59 + t34 * t55 + t45;
t13 = t33 * t52 + t51;
t14 = -t33 * t53 + t50;
t42 = pkin(3) * t46 - t13 * t34 + t14 * t23 + t44;
t11 = -t33 * t50 + t53;
t12 = t33 * t51 + t52;
t41 = t12 * t23 + (-pkin(7) - t59) * t54 - t11 * t34 + t47;
t38 = cos(qJ(5));
t35 = sin(qJ(5));
t29 = pkin(10) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t8 = t33 * t24 + t25 * t57;
t7 = t24 * t57 - t33 * t25;
t4 = t14 * t25 + t24 * t56;
t3 = t14 * t24 - t25 * t56;
t2 = t12 * t25 - t24 * t54;
t1 = t12 * t24 + t25 * t54;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t48)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t44) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t47) + g(3) * (t33 * rSges(3,3) + t45) + (g(1) * rSges(3,3) * t37 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t39) + g(2) * (-rSges(3,3) - pkin(7)) * t40) * t31) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t32 + t46) * rSges(4,1) + (-t14 * t30 + t32 * t56) * rSges(4,2) + t49 * t13 + t44) + g(2) * (t12 * pkin(2) - pkin(7) * t54 + (t12 * t32 - t30 * t54) * rSges(4,1) + (-t12 * t30 - t32 * t54) * rSges(4,2) + t49 * t11 + t47) + g(3) * ((t30 * rSges(4,1) + t32 * rSges(4,2)) * t33 + (-t49 * t39 + (t32 * rSges(4,1) - t30 * rSges(4,2) + pkin(2)) * t36) * t31 + t45)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t42) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t41) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t55 + t43)) - m(6) * (g(1) * (t4 * pkin(4) + (t13 * t35 + t4 * t38) * rSges(6,1) + (t13 * t38 - t4 * t35) * rSges(6,2) + t58 * t3 + t42) + g(2) * (t2 * pkin(4) + (t11 * t35 + t2 * t38) * rSges(6,1) + (t11 * t38 - t2 * t35) * rSges(6,2) + t58 * t1 + t41) + g(3) * (t8 * pkin(4) + (-t35 * t55 + t8 * t38) * rSges(6,1) + (-t8 * t35 - t38 * t55) * rSges(6,2) + t58 * t7 + t43));
U = t5;
