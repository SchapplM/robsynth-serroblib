% Calculate potential energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:24
% EndTime: 2019-03-08 21:12:25
% DurationCPUTime: 0.49s
% Computational Cost: add. (382->126), mult. (807->156), div. (0->0), fcn. (997->12), ass. (0->53)
t37 = sin(pkin(6));
t71 = pkin(7) * t37;
t70 = rSges(4,3) + pkin(8);
t69 = cos(qJ(3));
t41 = sin(qJ(3));
t68 = t37 * t41;
t42 = sin(qJ(2));
t67 = t37 * t42;
t44 = cos(qJ(2));
t66 = t37 * t44;
t65 = rSges(6,2) + qJ(4);
t64 = rSges(5,3) + qJ(4);
t63 = rSges(6,3) + qJ(5);
t62 = cos(pkin(6));
t36 = sin(pkin(10));
t61 = t36 * pkin(1) + r_base(2);
t60 = qJ(1) + r_base(3);
t58 = t37 * t69;
t57 = t42 * t62;
t56 = t44 * t62;
t39 = cos(pkin(10));
t55 = t39 * pkin(1) + t36 * t71 + r_base(1);
t54 = t62 * pkin(7) + t60;
t24 = -t36 * t57 + t39 * t44;
t53 = t24 * pkin(2) + t55;
t52 = pkin(2) * t67 + t54;
t22 = t36 * t44 + t39 * t57;
t51 = t22 * pkin(2) - t39 * t71 + t61;
t15 = t24 * t69 + t36 * t68;
t23 = t36 * t56 + t39 * t42;
t50 = t15 * pkin(3) + t23 * pkin(8) + t53;
t35 = sin(pkin(11));
t38 = cos(pkin(11));
t6 = t15 * t38 + t23 * t35;
t49 = t6 * pkin(4) + t50;
t26 = t62 * t41 + t42 * t58;
t48 = t26 * pkin(3) - pkin(8) * t66 + t52;
t11 = t26 * t38 - t35 * t66;
t47 = t11 * pkin(4) + t48;
t13 = t22 * t69 - t39 * t68;
t21 = t36 * t42 - t39 * t56;
t46 = t13 * pkin(3) + t21 * pkin(8) + t51;
t4 = t13 * t38 + t21 * t35;
t45 = t4 * pkin(4) + t46;
t43 = cos(qJ(6));
t40 = sin(qJ(6));
t25 = t41 * t67 - t62 * t69;
t14 = t24 * t41 - t36 * t58;
t12 = t22 * t41 + t39 * t58;
t10 = t26 * t35 + t38 * t66;
t5 = t15 * t35 - t23 * t38;
t3 = t13 * t35 - t21 * t38;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t39 * rSges(2,1) - t36 * rSges(2,2) + r_base(1)) + g(2) * (t36 * rSges(2,1) + t39 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t60)) - m(3) * (g(1) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t55) + g(2) * (t22 * rSges(3,1) - t21 * rSges(3,2) + t61) + g(3) * (t62 * rSges(3,3) + t54) + (g(1) * rSges(3,3) * t36 + g(3) * (rSges(3,1) * t42 + rSges(3,2) * t44) + g(2) * (-rSges(3,3) - pkin(7)) * t39) * t37) - m(4) * (g(1) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t70 * t23 + t53) + g(2) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t70 * t21 + t51) + g(3) * (t26 * rSges(4,1) - t25 * rSges(4,2) - t70 * t66 + t52)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t64 * t14 + t50) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t64 * t12 + t46) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t64 * t25 + t48)) - m(6) * (g(1) * (t6 * rSges(6,1) + t65 * t14 + t63 * t5 + t49) + g(2) * (t4 * rSges(6,1) + t65 * t12 + t63 * t3 + t45) + g(3) * (t11 * rSges(6,1) + t63 * t10 + t65 * t25 + t47)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t5 * t40 + t6 * t43) * rSges(7,1) + (-t6 * t40 + t5 * t43) * rSges(7,2) + t49) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t3 * t40 + t4 * t43) * rSges(7,1) + (t3 * t43 - t4 * t40) * rSges(7,2) + t45) + g(3) * (t11 * pkin(5) + t10 * qJ(5) + (t10 * t40 + t11 * t43) * rSges(7,1) + (t10 * t43 - t11 * t40) * rSges(7,2) + t47) + (g(1) * t14 + g(2) * t12 + g(3) * t25) * (-pkin(9) + qJ(4) - rSges(7,3)));
U  = t1;
