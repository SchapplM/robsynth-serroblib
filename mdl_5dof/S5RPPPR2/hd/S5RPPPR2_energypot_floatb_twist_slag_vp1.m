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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:15
% EndTime: 2020-01-03 11:22:15
% DurationCPUTime: 0.49s
% Computational Cost: add. (176->101), mult. (313->123), div. (0->0), fcn. (348->10), ass. (0->44)
t23 = sin(pkin(7));
t26 = cos(pkin(7));
t55 = t26 * pkin(2) + t23 * qJ(3);
t54 = pkin(6) + rSges(6,3);
t52 = g(2) * qJ(2);
t22 = sin(pkin(8));
t51 = t22 * t23;
t25 = cos(pkin(8));
t50 = t23 * t25;
t28 = sin(qJ(1));
t49 = t23 * t28;
t48 = t28 * t22;
t47 = t28 * t25;
t30 = cos(qJ(1));
t46 = t30 * t22;
t45 = t30 * t23;
t44 = t30 * t25;
t42 = -rSges(3,3) - qJ(2);
t41 = pkin(5) + r_base(1);
t40 = t28 * pkin(1) + r_base(2);
t39 = t23 * pkin(2) + t41;
t38 = -t28 * qJ(2) + r_base(3);
t37 = t55 * t28 + t40;
t36 = rSges(3,1) * t26 - rSges(3,2) * t23;
t35 = -pkin(1) - t55;
t11 = t26 * t46 - t47;
t12 = -t26 * t44 - t48;
t34 = t12 * pkin(3) - t11 * qJ(4) + t38;
t10 = t26 * t47 - t46;
t9 = t26 * t48 + t44;
t33 = t10 * pkin(3) + t9 * qJ(4) + t37;
t32 = pkin(3) * t50 - t26 * qJ(3) + qJ(4) * t51 + t39;
t31 = (g(3) * t35 - t52) * t30;
t29 = cos(qJ(5));
t27 = sin(qJ(5));
t24 = cos(pkin(9));
t21 = sin(pkin(9));
t8 = -t26 * t21 + t24 * t50;
t7 = t21 * t50 + t26 * t24;
t4 = t12 * t24 - t21 * t45;
t3 = t12 * t21 + t24 * t45;
t2 = t10 * t24 + t21 * t49;
t1 = t10 * t21 - t24 * t49;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t41) + g(2) * (t28 * rSges(2,1) + t30 * rSges(2,2) + r_base(2)) + g(3) * (-t30 * rSges(2,1) + t28 * rSges(2,2) + r_base(3))) - m(3) * (g(1) * (t23 * rSges(3,1) + t26 * rSges(3,2) + t41) + g(2) * t40 + g(3) * r_base(3) + (g(2) * t36 + g(3) * t42) * t28 + (g(2) * t42 + g(3) * (-pkin(1) - t36)) * t30) - m(4) * (g(1) * ((-rSges(4,3) - qJ(3)) * t26 + (rSges(4,1) * t25 - rSges(4,2) * t22) * t23 + t39) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + rSges(4,3) * t49 + t37) + g(3) * (t12 * rSges(4,1) + t11 * rSges(4,2) + t38) + (-t52 + g(3) * (-rSges(4,3) * t23 + t35)) * t30) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) + rSges(5,3) * t51 + t32) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t9 * rSges(5,3) + t33) + g(3) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t11 * rSges(5,3) + t34) + t31) - m(6) * (g(1) * (t8 * pkin(4) + (t27 * t51 + t8 * t29) * rSges(6,1) + (-t8 * t27 + t29 * t51) * rSges(6,2) + t54 * t7 + t32) + g(2) * (t2 * pkin(4) + (t2 * t29 + t9 * t27) * rSges(6,1) + (-t2 * t27 + t9 * t29) * rSges(6,2) + t33 + t54 * t1) + g(3) * (t4 * pkin(4) + (-t11 * t27 + t4 * t29) * rSges(6,1) + (-t11 * t29 - t4 * t27) * rSges(6,2) + t34 + t54 * t3) + t31);
U = t5;
