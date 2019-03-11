% Calculate potential energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:26
% EndTime: 2019-03-09 15:56:27
% DurationCPUTime: 0.58s
% Computational Cost: add. (227->122), mult. (362->146), div. (0->0), fcn. (388->10), ass. (0->41)
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t56 = g(1) * t30 + g(2) * t27;
t29 = cos(qJ(2));
t53 = g(3) * t29;
t22 = sin(pkin(10));
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t47 = t27 * t29;
t3 = t25 * t47 + t28 * t30;
t52 = t22 * t3;
t46 = t29 * t30;
t5 = t25 * t46 - t27 * t28;
t51 = t22 * t5;
t50 = t22 * t25;
t26 = sin(qJ(2));
t49 = t26 * t27;
t48 = t26 * t30;
t45 = -rSges(7,3) - pkin(9) - qJ(5);
t44 = rSges(5,3) + qJ(4);
t43 = -rSges(6,3) - qJ(5);
t42 = pkin(6) + r_base(3);
t41 = t27 * pkin(1) + r_base(2);
t40 = t26 * pkin(2) + t42;
t39 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t38 = t40 + (pkin(3) * t28 + qJ(4) * t25) * t26;
t37 = pkin(2) * t46 + pkin(8) * t48 + t39;
t6 = t27 * t25 + t28 * t46;
t36 = t6 * pkin(3) + t37;
t35 = g(3) * t38;
t34 = pkin(2) * t47 - t30 * pkin(7) + pkin(8) * t49 + t41;
t4 = -t25 * t30 + t28 * t47;
t33 = t4 * pkin(3) + t34;
t32 = t5 * qJ(4) + t36;
t31 = t3 * qJ(4) + t33;
t23 = cos(pkin(10));
t21 = pkin(10) + qJ(6);
t16 = cos(t21);
t15 = sin(t21);
t13 = pkin(5) * t23 + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t42)) - m(3) * (g(1) * (t27 * rSges(3,3) + t39) + g(2) * (rSges(3,1) * t47 - rSges(3,2) * t49 + t41) + g(3) * (rSges(3,1) * t26 + rSges(3,2) * t29 + t42) + (g(1) * (rSges(3,1) * t29 - rSges(3,2) * t26) + g(2) * (-rSges(3,3) - pkin(7))) * t30) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + rSges(4,3) * t48 + t37) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + rSges(4,3) * t49 + t34) + g(3) * ((-rSges(4,3) - pkin(8)) * t29 + (rSges(4,1) * t28 - rSges(4,2) * t25) * t26 + t40)) - m(5) * (g(1) * (t6 * rSges(5,1) + rSges(5,2) * t48 + t44 * t5 + t36) + g(2) * (t4 * rSges(5,1) + rSges(5,2) * t49 + t44 * t3 + t33) + g(3) * ((-rSges(5,2) - pkin(8)) * t29 + (rSges(5,1) * t28 + rSges(5,3) * t25) * t26 + t38)) - m(6) * (g(1) * (t6 * pkin(4) + (t23 * t6 + t51) * rSges(6,1) + (-t22 * t6 + t23 * t5) * rSges(6,2) + t32) + g(2) * (t4 * pkin(4) + (t23 * t4 + t52) * rSges(6,1) + (-t22 * t4 + t23 * t3) * rSges(6,2) + t31) + t35 + (-pkin(8) - t43) * t53 + (g(3) * (t28 * pkin(4) + (t23 * t28 + t50) * rSges(6,1) + (-t22 * t28 + t23 * t25) * rSges(6,2)) + t56 * t43) * t26) - m(7) * (g(1) * (t6 * t13 + pkin(5) * t51 + (t15 * t5 + t16 * t6) * rSges(7,1) + (-t15 * t6 + t16 * t5) * rSges(7,2) + t32) + g(2) * (t4 * t13 + pkin(5) * t52 + (t15 * t3 + t16 * t4) * rSges(7,1) + (-t15 * t4 + t16 * t3) * rSges(7,2) + t31) + t35 + (-pkin(8) - t45) * t53 + (g(3) * (t28 * t13 + pkin(5) * t50 + (t15 * t25 + t16 * t28) * rSges(7,1) + (-t15 * t28 + t16 * t25) * rSges(7,2)) + t56 * t45) * t26);
U  = t1;
