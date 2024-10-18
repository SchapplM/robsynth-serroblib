% Calculate potential energy for
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:07
% EndTime: 2024-09-27 21:45:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (260->98), mult. (178->105), div. (0->0), fcn. (154->22), ass. (0->41)
t42 = pkin(7) + pkin(8);
t57 = pkin(9) + t42 + rSges(6,3);
t56 = t42 + rSges(5,3);
t55 = rSges(4,3) + pkin(7);
t41 = cos(qJ(3));
t21 = t41 * pkin(3) + pkin(2);
t38 = cos(pkin(5));
t40 = sin(qJ(3));
t53 = t38 * t40;
t52 = t38 * t41;
t34 = qJ(3) + qJ(4);
t35 = sin(pkin(10));
t51 = t35 * pkin(1) + r_base(2);
t37 = cos(pkin(10));
t50 = t37 * pkin(1) + r_base(1);
t49 = qJ(1) + r_base(3);
t25 = pkin(5) - t34;
t24 = pkin(5) + t34;
t48 = pkin(6) + t49;
t47 = t41 * sin(qJ(4)) * pkin(4) + t40 * (cos(qJ(4)) * pkin(4) + pkin(3));
t27 = cos(t34);
t30 = qJ(5) + t34;
t46 = cos(t30) * rSges(6,1) - sin(t30) * rSges(6,2) + pkin(4) * t27 + t21;
t45 = t27 * rSges(5,1) - sin(t34) * rSges(5,2) + t21;
t10 = sin(t24) / 0.2e1;
t11 = cos(t25) / 0.2e1;
t14 = sin(t25);
t15 = cos(t24);
t36 = sin(pkin(5));
t44 = (t10 - t14 / 0.2e1) * rSges(5,1) + (t11 + t15 / 0.2e1) * rSges(5,2) + pkin(3) * t53 - t56 * t36;
t17 = -qJ(5) + t25;
t12 = sin(t17);
t16 = qJ(5) + t24;
t13 = cos(t16);
t8 = sin(t16) / 0.2e1;
t9 = cos(t17) / 0.2e1;
t43 = (t8 - t12 / 0.2e1) * rSges(6,1) + (t9 + t13 / 0.2e1) * rSges(6,2) + t47 * t38 - t57 * t36;
t32 = pkin(10) + qJ(2);
t23 = cos(t32);
t22 = sin(t32);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t49)) - m(3) * (g(1) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t50) + g(2) * (t22 * rSges(3,1) + t23 * rSges(3,2) + t51) + g(3) * (rSges(3,3) + t48)) - m(4) * (g(1) * (t23 * pkin(2) + (-t22 * t53 + t23 * t41) * rSges(4,1) + (-t22 * t52 - t23 * t40) * rSges(4,2) + t50) + g(2) * (t22 * pkin(2) + (t22 * t41 + t23 * t53) * rSges(4,1) + (-t22 * t40 + t23 * t52) * rSges(4,2) + t51) + g(3) * (t55 * t38 + t48) + (g(3) * (rSges(4,1) * t40 + rSges(4,2) * t41) + (g(1) * t22 - g(2) * t23) * t55) * t36) - m(5) * (g(1) * (-t44 * t22 + t45 * t23 + t50) + g(2) * (t45 * t22 + t44 * t23 + t51) + g(3) * (t36 * t40 * pkin(3) + (t11 - t15 / 0.2e1) * rSges(5,1) + (t10 + t14 / 0.2e1) * rSges(5,2) + t56 * t38 + t48)) - m(6) * (g(1) * (-t43 * t22 + t46 * t23 + t50) + g(2) * (t46 * t22 + t43 * t23 + t51) + g(3) * ((t9 - t13 / 0.2e1) * rSges(6,1) + (t8 + t12 / 0.2e1) * rSges(6,2) + t57 * t38 + t47 * t36 + t48));
U = t1;
