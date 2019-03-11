% Calculate potential energy for
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:31
% EndTime: 2019-03-09 16:11:31
% DurationCPUTime: 0.49s
% Computational Cost: add. (382->126), mult. (807->156), div. (0->0), fcn. (997->12), ass. (0->53)
t71 = rSges(4,3) + pkin(9);
t70 = cos(qJ(3));
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t69 = t36 * t40;
t41 = sin(qJ(1));
t68 = t36 * t41;
t43 = cos(qJ(2));
t67 = t36 * t43;
t44 = cos(qJ(1));
t66 = t36 * t44;
t65 = rSges(6,2) + qJ(4);
t64 = rSges(5,3) + qJ(4);
t63 = rSges(6,3) + qJ(5);
t62 = cos(pkin(6));
t61 = pkin(7) + r_base(3);
t60 = t41 * pkin(1) + r_base(2);
t58 = t36 * t70;
t57 = t62 * pkin(8) + t61;
t56 = t41 * t62;
t55 = t44 * t62;
t54 = t44 * pkin(1) + pkin(8) * t68 + r_base(1);
t53 = pkin(2) * t69 + t57;
t26 = -t40 * t56 + t44 * t43;
t52 = t26 * pkin(2) + t54;
t24 = t40 * t55 + t41 * t43;
t51 = t24 * pkin(2) - pkin(8) * t66 + t60;
t39 = sin(qJ(3));
t15 = t26 * t70 + t39 * t68;
t25 = t44 * t40 + t43 * t56;
t50 = t15 * pkin(3) + t25 * pkin(9) + t52;
t35 = sin(pkin(11));
t37 = cos(pkin(11));
t6 = t15 * t37 + t25 * t35;
t49 = t6 * pkin(4) + t50;
t22 = t62 * t39 + t40 * t58;
t48 = t22 * pkin(3) - pkin(9) * t67 + t53;
t11 = t22 * t37 - t35 * t67;
t47 = t11 * pkin(4) + t48;
t13 = t24 * t70 - t39 * t66;
t23 = t41 * t40 - t43 * t55;
t46 = t13 * pkin(3) + t23 * pkin(9) + t51;
t4 = t13 * t37 + t23 * t35;
t45 = t4 * pkin(4) + t46;
t42 = cos(qJ(6));
t38 = sin(qJ(6));
t21 = t39 * t69 - t62 * t70;
t14 = t26 * t39 - t41 * t58;
t12 = t24 * t39 + t44 * t58;
t10 = t22 * t35 + t37 * t67;
t5 = t15 * t35 - t25 * t37;
t3 = t13 * t35 - t23 * t37;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t44 * rSges(2,1) - t41 * rSges(2,2) + r_base(1)) + g(2) * (t41 * rSges(2,1) + t44 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t61)) - m(3) * (g(1) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t54) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t60) + g(3) * (t62 * rSges(3,3) + t57) + (g(1) * rSges(3,3) * t41 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t43) + g(2) * (-rSges(3,3) - pkin(8)) * t44) * t36) - m(4) * (g(1) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t71 * t25 + t52) + g(2) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t71 * t23 + t51) + g(3) * (t22 * rSges(4,1) - t21 * rSges(4,2) - t71 * t67 + t53)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t64 * t14 + t50) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t64 * t12 + t46) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t64 * t21 + t48)) - m(6) * (g(1) * (t6 * rSges(6,1) + t65 * t14 + t63 * t5 + t49) + g(2) * (t4 * rSges(6,1) + t65 * t12 + t63 * t3 + t45) + g(3) * (t11 * rSges(6,1) + t63 * t10 + t65 * t21 + t47)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t5 * t38 + t6 * t42) * rSges(7,1) + (-t6 * t38 + t5 * t42) * rSges(7,2) + t49) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t3 * t38 + t4 * t42) * rSges(7,1) + (t3 * t42 - t4 * t38) * rSges(7,2) + t45) + g(3) * (t11 * pkin(5) + t10 * qJ(5) + (t10 * t38 + t11 * t42) * rSges(7,1) + (t10 * t42 - t11 * t38) * rSges(7,2) + t47) + (g(1) * t14 + g(2) * t12 + g(3) * t21) * (-pkin(10) + qJ(4) - rSges(7,3)));
U  = t1;
