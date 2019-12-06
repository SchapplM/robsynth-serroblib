% Calculate potential energy for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:19
% EndTime: 2019-12-05 16:02:20
% DurationCPUTime: 0.41s
% Computational Cost: add. (193->105), mult. (358->133), div. (0->0), fcn. (405->10), ass. (0->43)
t24 = sin(pkin(9));
t55 = g(1) * t24;
t26 = cos(pkin(9));
t54 = g(2) * t26;
t53 = pkin(8) + rSges(6,3);
t25 = sin(pkin(5));
t52 = t24 * t25;
t29 = sin(qJ(4));
t51 = t25 * t29;
t30 = sin(qJ(2));
t50 = t25 * t30;
t32 = cos(qJ(4));
t49 = t25 * t32;
t33 = cos(qJ(2));
t48 = t25 * t33;
t27 = cos(pkin(5));
t47 = t27 * t30;
t46 = t27 * t33;
t45 = t33 * qJ(3);
t44 = t24 * pkin(1) + r_base(2);
t43 = qJ(1) + r_base(3);
t42 = (-pkin(3) - pkin(6)) * t26;
t41 = t26 * pkin(1) + pkin(6) * t52 + r_base(1);
t40 = t27 * pkin(6) + t43;
t39 = pkin(2) * t50 + t40;
t8 = t24 * t30 - t26 * t46;
t9 = t24 * t33 + t26 * t47;
t38 = t9 * pkin(2) + t8 * qJ(3) + t44;
t37 = t27 * pkin(3) + pkin(7) * t50 + t39;
t10 = t24 * t46 + t26 * t30;
t11 = -t24 * t47 + t26 * t33;
t36 = t11 * pkin(2) + t10 * qJ(3) + t41;
t35 = pkin(3) * t52 + t36;
t34 = t9 * pkin(7) + t38;
t31 = cos(qJ(5));
t28 = sin(qJ(5));
t13 = t27 * t32 - t29 * t48;
t12 = t27 * t29 + t32 * t48;
t4 = -t26 * t49 + t8 * t29;
t3 = t26 * t51 + t8 * t32;
t2 = t10 * t29 + t24 * t49;
t1 = -t10 * t32 + t24 * t51;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t43)) - m(3) * (g(1) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t41) + g(2) * (t9 * rSges(3,1) - t8 * rSges(3,2) + t44) + g(3) * (t27 * rSges(3,3) + t40) + (rSges(3,3) * t55 + g(3) * (rSges(3,1) * t30 + rSges(3,2) * t33) + (-rSges(3,3) - pkin(6)) * t54) * t25) - m(4) * (g(1) * (-t11 * rSges(4,2) + t10 * rSges(4,3) + t36) + g(2) * (-t9 * rSges(4,2) + t8 * rSges(4,3) + t38) + g(3) * (t27 * rSges(4,1) + t39) + (rSges(4,1) * t55 + g(3) * (-rSges(4,2) * t30 - rSges(4,3) * t33 - t45) + (-rSges(4,1) - pkin(6)) * t54) * t25) - m(5) * (g(1) * (t2 * rSges(5,1) - t1 * rSges(5,2) + (rSges(5,3) + pkin(7)) * t11 + t35) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t9 * rSges(5,3) + t34) + g(3) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t37) + (g(3) * (rSges(5,3) * t30 - t45) + g(2) * t42) * t25) - m(6) * (g(1) * (t2 * pkin(4) + t11 * pkin(7) + (t11 * t28 + t2 * t31) * rSges(6,1) + (t11 * t31 - t2 * t28) * rSges(6,2) + t53 * t1 + t35) + g(2) * (t4 * pkin(4) + (t9 * t28 + t4 * t31) * rSges(6,1) + (-t4 * t28 + t9 * t31) * rSges(6,2) - t53 * t3 + t25 * t42 + t34) + g(3) * (t13 * pkin(4) - t25 * t45 + (t13 * t31 + t28 * t50) * rSges(6,1) + (-t13 * t28 + t31 * t50) * rSges(6,2) + t53 * t12 + t37));
U = t5;
