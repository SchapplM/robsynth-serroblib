% Calculate potential energy for
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR14_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:17
% EndTime: 2024-09-27 18:42:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (260->98), mult. (178->105), div. (0->0), fcn. (154->22), ass. (0->41)
t42 = pkin(8) + pkin(9);
t57 = pkin(10) + t42 + rSges(6,3);
t56 = t42 + rSges(5,3);
t55 = rSges(4,3) + pkin(8);
t40 = cos(qJ(3));
t21 = t40 * pkin(3) + pkin(2);
t36 = cos(pkin(5));
t38 = sin(qJ(3));
t53 = t36 * t38;
t52 = t36 * t40;
t33 = qJ(3) + qJ(4);
t51 = pkin(6) + r_base(3);
t39 = sin(qJ(1));
t50 = t39 * pkin(1) + r_base(2);
t41 = cos(qJ(1));
t49 = t41 * pkin(1) + r_base(1);
t48 = pkin(7) + t51;
t23 = pkin(5) - t33;
t22 = pkin(5) + t33;
t47 = t40 * sin(qJ(4)) * pkin(4) + t38 * (cos(qJ(4)) * pkin(4) + pkin(3));
t26 = cos(t33);
t28 = qJ(5) + t33;
t46 = cos(t28) * rSges(6,1) - sin(t28) * rSges(6,2) + pkin(4) * t26 + t21;
t45 = t26 * rSges(5,1) - sin(t33) * rSges(5,2) + t21;
t10 = sin(t22) / 0.2e1;
t11 = cos(t23) / 0.2e1;
t14 = sin(t23);
t15 = cos(t22);
t35 = sin(pkin(5));
t44 = (t10 - t14 / 0.2e1) * rSges(5,1) + (t11 + t15 / 0.2e1) * rSges(5,2) + pkin(3) * t53 - t56 * t35;
t17 = -qJ(5) + t23;
t12 = sin(t17);
t16 = qJ(5) + t22;
t13 = cos(t16);
t8 = sin(t16) / 0.2e1;
t9 = cos(t17) / 0.2e1;
t43 = (t8 - t12 / 0.2e1) * rSges(6,1) + (t9 + t13 / 0.2e1) * rSges(6,2) + t47 * t36 - t57 * t35;
t34 = qJ(1) + qJ(2);
t27 = cos(t34);
t25 = sin(t34);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t41 * rSges(2,1) - t39 * rSges(2,2) + r_base(1)) + g(2) * (t39 * rSges(2,1) + t41 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (t27 * rSges(3,1) - t25 * rSges(3,2) + t49) + g(2) * (t25 * rSges(3,1) + t27 * rSges(3,2) + t50) + g(3) * (rSges(3,3) + t48)) - m(4) * (g(1) * (t27 * pkin(2) + (-t25 * t53 + t27 * t40) * rSges(4,1) + (-t25 * t52 - t27 * t38) * rSges(4,2) + t49) + g(2) * (t25 * pkin(2) + (t25 * t40 + t27 * t53) * rSges(4,1) + (-t25 * t38 + t27 * t52) * rSges(4,2) + t50) + g(3) * (t55 * t36 + t48) + (g(3) * (rSges(4,1) * t38 + rSges(4,2) * t40) + (g(1) * t25 - g(2) * t27) * t55) * t35) - m(5) * (g(1) * (-t44 * t25 + t45 * t27 + t49) + g(2) * (t45 * t25 + t44 * t27 + t50) + g(3) * (t35 * t38 * pkin(3) + (t11 - t15 / 0.2e1) * rSges(5,1) + (t10 + t14 / 0.2e1) * rSges(5,2) + t56 * t36 + t48)) - m(6) * (g(1) * (-t43 * t25 + t46 * t27 + t49) + g(2) * (t46 * t25 + t43 * t27 + t50) + g(3) * ((t9 - t13 / 0.2e1) * rSges(6,1) + (t8 + t12 / 0.2e1) * rSges(6,2) + t57 * t36 + t47 * t35 + t48));
U = t1;
