% Calculate potential energy for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:12
% EndTime: 2019-12-31 21:43:13
% DurationCPUTime: 0.30s
% Computational Cost: add. (214->93), mult. (413->108), div. (0->0), fcn. (481->10), ass. (0->46)
t59 = rSges(5,1) + pkin(8);
t58 = rSges(4,3) + pkin(8);
t57 = pkin(9) + rSges(6,3);
t25 = sin(pkin(5));
t29 = sin(qJ(2));
t56 = t25 * t29;
t30 = sin(qJ(1));
t55 = t25 * t30;
t32 = cos(qJ(3));
t54 = t25 * t32;
t33 = cos(qJ(2));
t53 = t25 * t33;
t34 = cos(qJ(1));
t52 = t25 * t34;
t51 = t30 * t29;
t50 = t30 * t33;
t49 = t34 * t29;
t48 = t34 * t33;
t47 = rSges(5,3) + qJ(4);
t46 = pkin(6) + r_base(3);
t45 = t30 * pkin(1) + r_base(2);
t26 = cos(pkin(5));
t44 = t26 * pkin(7) + t46;
t43 = t34 * pkin(1) + pkin(7) * t55 + r_base(1);
t42 = pkin(2) * t56 + t44;
t15 = -t26 * t51 + t48;
t41 = t15 * pkin(2) + t43;
t28 = sin(qJ(3));
t11 = t26 * t28 + t29 * t54;
t40 = t11 * pkin(3) + t42;
t6 = t15 * t32 + t28 * t55;
t39 = t6 * pkin(3) + t41;
t27 = sin(qJ(5));
t31 = cos(qJ(5));
t38 = t27 * rSges(6,1) + t31 * rSges(6,2) + qJ(4);
t37 = t31 * rSges(6,1) - t27 * rSges(6,2) + pkin(4) + pkin(8);
t13 = t26 * t49 + t50;
t36 = t13 * pkin(2) - pkin(7) * t52 + t45;
t4 = t13 * t32 - t28 * t52;
t35 = t4 * pkin(3) + t36;
t14 = t26 * t50 + t49;
t12 = -t26 * t48 + t51;
t10 = -t26 * t32 + t28 * t56;
t5 = t15 * t28 - t30 * t54;
t3 = t13 * t28 + t32 * t52;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t46)) - m(3) * (g(1) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t43) + g(2) * (t13 * rSges(3,1) - t12 * rSges(3,2) + t45) + g(3) * (t26 * rSges(3,3) + t44) + (g(1) * rSges(3,3) * t30 + g(3) * (rSges(3,1) * t29 + rSges(3,2) * t33) + g(2) * (-rSges(3,3) - pkin(7)) * t34) * t25) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t58 * t14 + t41) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t58 * t12 + t36) + g(3) * (t11 * rSges(4,1) - t10 * rSges(4,2) - t58 * t53 + t42)) - m(5) * (g(1) * (-t6 * rSges(5,2) + t59 * t14 + t47 * t5 + t39) + g(2) * (-t4 * rSges(5,2) + t59 * t12 + t47 * t3 + t35) + g(3) * (-t11 * rSges(5,2) + t47 * t10 - t59 * t53 + t40)) - m(6) * (g(1) * (t37 * t14 + t38 * t5 + t57 * t6 + t39) + g(2) * (t37 * t12 + t38 * t3 + t57 * t4 + t35) + g(3) * (t38 * t10 + t57 * t11 - t37 * t53 + t40));
U = t1;
