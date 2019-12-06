% Calculate potential energy for
% S5PRRRR9
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:18:57
% EndTime: 2019-12-05 17:18:57
% DurationCPUTime: 0.35s
% Computational Cost: add. (237->103), mult. (435->127), div. (0->0), fcn. (511->12), ass. (0->43)
t26 = sin(pkin(5));
t54 = pkin(6) * t26;
t53 = rSges(4,3) + pkin(7);
t52 = pkin(8) + rSges(5,3);
t30 = sin(qJ(3));
t51 = t26 * t30;
t31 = sin(qJ(2));
t50 = t26 * t31;
t33 = cos(qJ(3));
t49 = t26 * t33;
t34 = cos(qJ(2));
t48 = t26 * t34;
t28 = cos(pkin(5));
t47 = t28 * t31;
t46 = t28 * t34;
t45 = pkin(9) + pkin(8) + rSges(6,3);
t25 = sin(pkin(10));
t44 = t25 * pkin(1) + r_base(2);
t43 = qJ(1) + r_base(3);
t27 = cos(pkin(10));
t42 = t27 * pkin(1) + t25 * t54 + r_base(1);
t41 = t28 * pkin(6) + t43;
t10 = -t25 * t47 + t27 * t34;
t40 = t10 * pkin(2) + t42;
t39 = pkin(2) * t50 + t41;
t24 = qJ(4) + qJ(5);
t19 = sin(t24);
t20 = cos(t24);
t32 = cos(qJ(4));
t38 = t20 * rSges(6,1) - t19 * rSges(6,2) + t32 * pkin(4) + pkin(3);
t8 = t25 * t34 + t27 * t47;
t37 = t8 * pkin(2) - t27 * t54 + t44;
t29 = sin(qJ(4));
t36 = t19 * rSges(6,1) + t20 * rSges(6,2) + t29 * pkin(4) + pkin(7);
t12 = t28 * t30 + t31 * t49;
t11 = -t28 * t33 + t30 * t50;
t9 = t25 * t46 + t27 * t31;
t7 = t25 * t31 - t27 * t46;
t4 = t10 * t33 + t25 * t51;
t3 = t10 * t30 - t25 * t49;
t2 = -t27 * t51 + t8 * t33;
t1 = t27 * t49 + t8 * t30;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t27 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t27 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t43)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t42) + g(2) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t44) + g(3) * (t28 * rSges(3,3) + t41) + (g(1) * rSges(3,3) * t25 + g(3) * (rSges(3,1) * t31 + rSges(3,2) * t34) + g(2) * (-rSges(3,3) - pkin(6)) * t27) * t26) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t53 * t9 + t40) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t53 * t7 + t37) + g(3) * (t12 * rSges(4,1) - t11 * rSges(4,2) - t53 * t48 + t39)) - m(5) * (g(1) * (t4 * pkin(3) + t9 * pkin(7) + (t9 * t29 + t4 * t32) * rSges(5,1) + (-t4 * t29 + t9 * t32) * rSges(5,2) + t52 * t3 + t40) + g(2) * (t2 * pkin(3) + t7 * pkin(7) + (t2 * t32 + t7 * t29) * rSges(5,1) + (-t2 * t29 + t7 * t32) * rSges(5,2) + t52 * t1 + t37) + g(3) * (t12 * pkin(3) - pkin(7) * t48 + (t12 * t32 - t29 * t48) * rSges(5,1) + (-t12 * t29 - t32 * t48) * rSges(5,2) + t52 * t11 + t39)) - m(6) * (g(1) * (t45 * t3 + t36 * t9 + t38 * t4 + t40) + g(2) * (t45 * t1 + t38 * t2 + t36 * t7 + t37) + g(3) * (t45 * t11 + t38 * t12 - t36 * t48 + t39));
U = t5;
