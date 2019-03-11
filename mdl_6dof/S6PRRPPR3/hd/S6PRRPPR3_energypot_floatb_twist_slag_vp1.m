% Calculate potential energy for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:03
% EndTime: 2019-03-08 21:08:03
% DurationCPUTime: 0.42s
% Computational Cost: add. (299->113), mult. (590->127), div. (0->0), fcn. (697->10), ass. (0->51)
t68 = rSges(6,1) + qJ(4);
t30 = sin(pkin(6));
t67 = pkin(7) * t30;
t66 = rSges(5,2) + pkin(8);
t65 = rSges(4,3) + pkin(8);
t64 = rSges(6,3) - pkin(8);
t63 = pkin(9) + rSges(7,3);
t62 = cos(qJ(3));
t33 = sin(qJ(3));
t61 = t30 * t33;
t34 = sin(qJ(2));
t60 = t30 * t34;
t36 = cos(qJ(2));
t59 = t30 * t36;
t58 = rSges(5,3) + qJ(4);
t57 = cos(pkin(6));
t29 = sin(pkin(10));
t56 = t29 * pkin(1) + r_base(2);
t55 = qJ(1) + r_base(3);
t54 = -qJ(5) - t64;
t53 = t30 * t62;
t52 = t34 * t57;
t51 = t36 * t57;
t31 = cos(pkin(10));
t50 = t31 * pkin(1) + t29 * t67 + r_base(1);
t49 = t57 * pkin(7) + t55;
t16 = -t29 * t52 + t31 * t36;
t48 = t16 * pkin(2) + t50;
t47 = pkin(2) * t60 + t49;
t8 = t16 * t62 + t29 * t61;
t46 = t8 * pkin(3) + t48;
t18 = t57 * t33 + t34 * t53;
t45 = t18 * pkin(3) + t47;
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t44 = rSges(7,1) * t32 + rSges(7,2) * t35 - pkin(8);
t43 = t8 * pkin(4) + t46;
t14 = t29 * t36 + t31 * t52;
t42 = t14 * pkin(2) - t31 * t67 + t56;
t41 = rSges(7,1) * t35 - rSges(7,2) * t32 + pkin(5) + qJ(4);
t40 = -qJ(5) - t44;
t39 = t18 * pkin(4) + qJ(5) * t59 + t45;
t6 = t14 * t62 - t31 * t61;
t38 = t6 * pkin(3) + t42;
t37 = t6 * pkin(4) + t38;
t17 = t33 * t60 - t57 * t62;
t15 = t29 * t51 + t31 * t34;
t13 = t29 * t34 - t31 * t51;
t7 = t16 * t33 - t29 * t53;
t5 = t14 * t33 + t31 * t53;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t31 - rSges(2,2) * t29 + r_base(1)) + g(2) * (rSges(2,1) * t29 + rSges(2,2) * t31 + r_base(2)) + g(3) * (rSges(2,3) + t55)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t50) + g(2) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t56) + g(3) * (t57 * rSges(3,3) + t49) + (g(1) * rSges(3,3) * t29 + g(3) * (rSges(3,1) * t34 + rSges(3,2) * t36) + g(2) * (-rSges(3,3) - pkin(7)) * t31) * t30) - m(4) * (g(1) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t65 * t15 + t48) + g(2) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t65 * t13 + t42) + g(3) * (t18 * rSges(4,1) - t17 * rSges(4,2) - t65 * t59 + t47)) - m(5) * (g(1) * (t8 * rSges(5,1) + t66 * t15 + t58 * t7 + t46) + g(2) * (t6 * rSges(5,1) + t66 * t13 + t58 * t5 + t38) + g(3) * (t18 * rSges(5,1) + t58 * t17 - t66 * t59 + t45)) - m(6) * (g(3) * (-t18 * rSges(6,2) + t17 * t68 + t64 * t59 + t39) + (-t6 * rSges(6,2) + t54 * t13 + t5 * t68 + t37) * g(2) + (-t8 * rSges(6,2) + t54 * t15 + t68 * t7 + t43) * g(1)) - m(7) * (g(1) * (t40 * t15 + t41 * t7 + t63 * t8 + t43) + g(2) * (t40 * t13 + t41 * t5 + t63 * t6 + t37) + g(3) * (t41 * t17 + t63 * t18 + t44 * t59 + t39));
U  = t1;
