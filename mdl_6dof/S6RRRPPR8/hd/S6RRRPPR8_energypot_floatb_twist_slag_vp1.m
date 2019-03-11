% Calculate potential energy for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:29
% EndTime: 2019-03-09 16:03:29
% DurationCPUTime: 0.42s
% Computational Cost: add. (299->113), mult. (590->127), div. (0->0), fcn. (697->10), ass. (0->51)
t68 = rSges(6,1) + qJ(4);
t67 = rSges(5,2) + pkin(9);
t66 = rSges(4,3) + pkin(9);
t65 = rSges(6,3) - pkin(9);
t64 = pkin(10) + rSges(7,3);
t63 = cos(qJ(3));
t29 = sin(pkin(6));
t32 = sin(qJ(2));
t62 = t29 * t32;
t33 = sin(qJ(1));
t61 = t29 * t33;
t35 = cos(qJ(2));
t60 = t29 * t35;
t36 = cos(qJ(1));
t59 = t29 * t36;
t58 = rSges(5,3) + qJ(4);
t57 = cos(pkin(6));
t56 = pkin(7) + r_base(3);
t55 = t33 * pkin(1) + r_base(2);
t54 = -qJ(5) - t65;
t53 = t29 * t63;
t52 = t57 * pkin(8) + t56;
t51 = t33 * t57;
t50 = t36 * t57;
t49 = t36 * pkin(1) + pkin(8) * t61 + r_base(1);
t48 = pkin(2) * t62 + t52;
t18 = -t32 * t51 + t36 * t35;
t47 = t18 * pkin(2) + t49;
t31 = sin(qJ(3));
t14 = t57 * t31 + t32 * t53;
t46 = t14 * pkin(3) + t48;
t8 = t18 * t63 + t31 * t61;
t45 = t8 * pkin(3) + t47;
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t44 = rSges(7,1) * t30 + rSges(7,2) * t34 - pkin(9);
t43 = t8 * pkin(4) + t45;
t42 = rSges(7,1) * t34 - rSges(7,2) * t30 + pkin(5) + qJ(4);
t41 = -qJ(5) - t44;
t16 = t32 * t50 + t33 * t35;
t40 = t16 * pkin(2) - pkin(8) * t59 + t55;
t39 = t14 * pkin(4) + qJ(5) * t60 + t46;
t6 = t16 * t63 - t31 * t59;
t38 = t6 * pkin(3) + t40;
t37 = t6 * pkin(4) + t38;
t17 = t36 * t32 + t35 * t51;
t15 = t32 * t33 - t35 * t50;
t13 = t31 * t62 - t57 * t63;
t7 = t18 * t31 - t33 * t53;
t5 = t16 * t31 + t36 * t53;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t36 - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + rSges(2,2) * t36 + r_base(2)) + g(3) * (rSges(2,3) + t56)) - m(3) * (g(1) * (rSges(3,1) * t18 - rSges(3,2) * t17 + t49) + g(2) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t55) + g(3) * (t57 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t33 + g(3) * (rSges(3,1) * t32 + rSges(3,2) * t35) + g(2) * (-rSges(3,3) - pkin(8)) * t36) * t29) - m(4) * (g(1) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t66 * t17 + t47) + g(2) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t66 * t15 + t40) + g(3) * (t14 * rSges(4,1) - t13 * rSges(4,2) - t66 * t60 + t48)) - m(5) * (g(1) * (t8 * rSges(5,1) + t67 * t17 + t58 * t7 + t45) + g(2) * (t6 * rSges(5,1) + t67 * t15 + t58 * t5 + t38) + g(3) * (t14 * rSges(5,1) + t58 * t13 - t67 * t60 + t46)) - m(6) * (g(3) * (-t14 * rSges(6,2) + t13 * t68 + t65 * t60 + t39) + (-t6 * rSges(6,2) + t54 * t15 + t5 * t68 + t37) * g(2) + (-t8 * rSges(6,2) + t54 * t17 + t68 * t7 + t43) * g(1)) - m(7) * (g(1) * (t41 * t17 + t42 * t7 + t64 * t8 + t43) + g(2) * (t41 * t15 + t42 * t5 + t64 * t6 + t37) + g(3) * (t42 * t13 + t64 * t14 + t44 * t60 + t39));
U  = t1;
