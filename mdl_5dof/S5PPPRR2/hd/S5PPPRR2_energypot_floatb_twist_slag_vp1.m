% Calculate potential energy for
% S5PPPRR2
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:22
% EndTime: 2019-12-05 14:59:22
% DurationCPUTime: 0.34s
% Computational Cost: add. (176->100), mult. (313->125), div. (0->0), fcn. (348->10), ass. (0->40)
t54 = rSges(5,3) + pkin(5);
t53 = pkin(6) + rSges(6,3);
t26 = sin(pkin(9));
t27 = sin(pkin(8));
t52 = t26 * t27;
t28 = sin(pkin(7));
t51 = t27 * t28;
t33 = sin(qJ(4));
t50 = t27 * t33;
t35 = cos(qJ(4));
t49 = t27 * t35;
t30 = cos(pkin(8));
t48 = t28 * t30;
t31 = cos(pkin(7));
t47 = t31 * t26;
t29 = cos(pkin(9));
t46 = t31 * t29;
t45 = qJ(3) * t27;
t44 = t28 * pkin(1) + r_base(2);
t43 = qJ(1) + r_base(3);
t42 = t31 * pkin(1) + t28 * qJ(2) + r_base(1);
t41 = t27 * pkin(2) + t43;
t40 = t42 + (pkin(2) * t30 + t45) * t31;
t10 = t28 * t26 + t30 * t46;
t39 = t10 * pkin(3) + t40;
t38 = pkin(2) * t48 - t31 * qJ(2) + t28 * t45 + t44;
t8 = t29 * t48 - t47;
t37 = t8 * pkin(3) + t38;
t36 = t27 * t29 * pkin(3) + pkin(5) * t52 - t30 * qJ(3) + t41;
t34 = cos(qJ(5));
t32 = sin(qJ(5));
t12 = t29 * t49 - t30 * t33;
t11 = t29 * t50 + t30 * t35;
t9 = -t28 * t29 + t30 * t47;
t7 = t26 * t48 + t46;
t4 = t10 * t35 + t31 * t50;
t3 = t10 * t33 - t31 * t49;
t2 = t28 * t50 + t8 * t35;
t1 = -t28 * t49 + t8 * t33;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t31 * rSges(2,1) - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + t31 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t43)) - m(3) * (g(1) * (t28 * rSges(3,3) + t42) + g(2) * (rSges(3,1) * t48 - rSges(3,2) * t51 + t44) + g(3) * (t27 * rSges(3,1) + t30 * rSges(3,2) + t43) + (g(1) * (rSges(3,1) * t30 - rSges(3,2) * t27) + g(2) * (-rSges(3,3) - qJ(2))) * t31) - m(4) * (g(1) * (t31 * t27 * rSges(4,3) + t10 * rSges(4,1) - t9 * rSges(4,2) + t40) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + rSges(4,3) * t51 + t38) + g(3) * ((-rSges(4,3) - qJ(3)) * t30 + (rSges(4,1) * t29 - rSges(4,2) * t26) * t27 + t41)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t54 * t9 + t39) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t54 * t7 + t37) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) + rSges(5,3) * t52 + t36)) - m(6) * (g(1) * (t4 * pkin(4) + t9 * pkin(5) + (t9 * t32 + t4 * t34) * rSges(6,1) + (-t4 * t32 + t9 * t34) * rSges(6,2) + t53 * t3 + t39) + g(2) * (t2 * pkin(4) + t7 * pkin(5) + (t2 * t34 + t7 * t32) * rSges(6,1) + (-t2 * t32 + t7 * t34) * rSges(6,2) + t53 * t1 + t37) + g(3) * (t12 * pkin(4) + (t12 * t34 + t32 * t52) * rSges(6,1) + (-t12 * t32 + t34 * t52) * rSges(6,2) + t53 * t11 + t36));
U = t5;
