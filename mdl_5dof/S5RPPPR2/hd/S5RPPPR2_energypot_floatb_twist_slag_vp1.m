% Calculate potential energy for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:53
% EndTime: 2019-12-05 17:30:53
% DurationCPUTime: 0.47s
% Computational Cost: add. (176->102), mult. (313->124), div. (0->0), fcn. (348->10), ass. (0->44)
t54 = pkin(6) + rSges(6,3);
t30 = sin(qJ(1));
t55 = g(2) * t30;
t24 = sin(pkin(8));
t25 = sin(pkin(7));
t53 = t24 * t25;
t27 = cos(pkin(8));
t52 = t25 * t27;
t32 = cos(qJ(1));
t51 = t25 * t32;
t28 = cos(pkin(7));
t50 = t28 * t32;
t49 = t30 * t24;
t48 = t30 * t25;
t47 = t30 * t27;
t46 = t32 * t24;
t45 = t32 * t27;
t44 = t25 * qJ(3);
t43 = -rSges(4,3) - qJ(3);
t42 = pkin(5) + r_base(1);
t41 = t32 * qJ(2) + r_base(2);
t40 = -t28 * pkin(2) - pkin(1);
t39 = t25 * pkin(2) + t42;
t38 = t32 * pkin(1) + t30 * qJ(2) + r_base(3);
t37 = pkin(2) * t50 + t32 * t44 + t38;
t10 = -t28 * t47 + t46;
t9 = t28 * t49 + t45;
t36 = t10 * pkin(3) - t9 * qJ(4) + t41;
t12 = t28 * t45 + t49;
t35 = t12 * pkin(3) + t37;
t34 = pkin(3) * t52 - t28 * qJ(3) + qJ(4) * t53 + t39;
t33 = (t40 - t44) * t55;
t31 = cos(qJ(5));
t29 = sin(qJ(5));
t26 = cos(pkin(9));
t23 = sin(pkin(9));
t11 = t28 * t46 - t47;
t8 = -t28 * t23 + t26 * t52;
t7 = t23 * t52 + t28 * t26;
t4 = t12 * t26 + t23 * t51;
t3 = t12 * t23 - t26 * t51;
t2 = t10 * t26 - t23 * t48;
t1 = t10 * t23 + t26 * t48;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t42) + g(2) * (-t30 * rSges(2,1) - t32 * rSges(2,2) + r_base(2)) + g(3) * (t32 * rSges(2,1) - t30 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (t25 * rSges(3,1) + t28 * rSges(3,2) + t42) + g(2) * (t32 * rSges(3,3) + t41) + g(3) * (rSges(3,1) * t50 - rSges(3,2) * t51 + t38) + (g(2) * (-rSges(3,1) * t28 + rSges(3,2) * t25 - pkin(1)) + g(3) * rSges(3,3)) * t30) - m(4) * (g(1) * (t43 * t28 + (rSges(4,1) * t27 - rSges(4,2) * t24) * t25 + t39) + g(2) * (t10 * rSges(4,1) + t9 * rSges(4,2) + t41) + g(3) * (t12 * rSges(4,1) - t11 * rSges(4,2) + rSges(4,3) * t51 + t37) + (t43 * t25 + t40) * t55) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) + rSges(5,3) * t53 + t34) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) - t9 * rSges(5,3) + t36) + g(3) * (t4 * rSges(5,1) - t3 * rSges(5,2) + (rSges(5,3) + qJ(4)) * t11 + t35) + t33) - m(6) * (g(1) * (t8 * pkin(4) + (t29 * t53 + t8 * t31) * rSges(6,1) + (-t8 * t29 + t31 * t53) * rSges(6,2) + t54 * t7 + t34) + g(2) * (t2 * pkin(4) + (t2 * t31 - t9 * t29) * rSges(6,1) + (-t2 * t29 - t9 * t31) * rSges(6,2) + t36 + t54 * t1) + g(3) * (t4 * pkin(4) + t11 * qJ(4) + (t11 * t29 + t4 * t31) * rSges(6,1) + (t11 * t31 - t4 * t29) * rSges(6,2) + t54 * t3 + t35) + t33);
U = t5;
