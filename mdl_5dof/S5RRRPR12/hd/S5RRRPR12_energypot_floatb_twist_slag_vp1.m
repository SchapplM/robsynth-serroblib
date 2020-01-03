% Calculate potential energy for
% S5RRRPR12
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:58
% EndTime: 2019-12-31 21:36:58
% DurationCPUTime: 0.34s
% Computational Cost: add. (237->103), mult. (435->125), div. (0->0), fcn. (511->12), ass. (0->45)
t56 = rSges(4,3) + pkin(8);
t26 = sin(pkin(5));
t31 = sin(qJ(2));
t55 = t26 * t31;
t32 = sin(qJ(1));
t54 = t26 * t32;
t33 = cos(qJ(3));
t53 = t26 * t33;
t34 = cos(qJ(2));
t52 = t26 * t34;
t35 = cos(qJ(1));
t51 = t26 * t35;
t50 = t32 * t31;
t49 = t32 * t34;
t48 = t35 * t31;
t47 = t35 * t34;
t46 = pkin(9) + qJ(4) + rSges(6,3);
t45 = qJ(4) + rSges(5,3);
t44 = pkin(6) + r_base(3);
t43 = t32 * pkin(1) + r_base(2);
t28 = cos(pkin(5));
t42 = t28 * pkin(7) + t44;
t41 = t35 * pkin(1) + pkin(7) * t54 + r_base(1);
t40 = pkin(2) * t55 + t42;
t12 = -t28 * t50 + t47;
t39 = t12 * pkin(2) + t41;
t24 = pkin(10) + qJ(5);
t19 = sin(t24);
t20 = cos(t24);
t27 = cos(pkin(10));
t38 = t20 * rSges(6,1) - t19 * rSges(6,2) + t27 * pkin(4) + pkin(3);
t10 = t28 * t48 + t49;
t37 = t10 * pkin(2) - pkin(7) * t51 + t43;
t25 = sin(pkin(10));
t36 = t19 * rSges(6,1) + t20 * rSges(6,2) + t25 * pkin(4) + pkin(8);
t30 = sin(qJ(3));
t11 = t28 * t49 + t48;
t9 = -t28 * t47 + t50;
t8 = t28 * t30 + t31 * t53;
t7 = -t28 * t33 + t30 * t55;
t4 = t12 * t33 + t30 * t54;
t3 = t12 * t30 - t32 * t53;
t2 = t10 * t33 - t30 * t51;
t1 = t10 * t30 + t33 * t51;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t35 * rSges(2,1) - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + t35 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t41) + g(2) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t43) + g(3) * (t28 * rSges(3,3) + t42) + (g(1) * rSges(3,3) * t32 + g(3) * (rSges(3,1) * t31 + rSges(3,2) * t34) + g(2) * (-rSges(3,3) - pkin(7)) * t35) * t26) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t56 * t11 + t39) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t56 * t9 + t37) + g(3) * (t8 * rSges(4,1) - t7 * rSges(4,2) - t56 * t52 + t40)) - m(5) * (g(1) * (t4 * pkin(3) + t11 * pkin(8) + (t11 * t25 + t4 * t27) * rSges(5,1) + (t11 * t27 - t4 * t25) * rSges(5,2) + t45 * t3 + t39) + g(2) * (t2 * pkin(3) + t9 * pkin(8) + (t2 * t27 + t9 * t25) * rSges(5,1) + (-t2 * t25 + t9 * t27) * rSges(5,2) + t45 * t1 + t37) + g(3) * (t8 * pkin(3) - pkin(8) * t52 + (-t25 * t52 + t8 * t27) * rSges(5,1) + (-t8 * t25 - t27 * t52) * rSges(5,2) + t45 * t7 + t40)) - m(6) * (g(1) * (t36 * t11 + t46 * t3 + t38 * t4 + t39) + g(2) * (t46 * t1 + t38 * t2 + t36 * t9 + t37) + g(3) * (-t36 * t52 + t38 * t8 + t46 * t7 + t40));
U = t5;
