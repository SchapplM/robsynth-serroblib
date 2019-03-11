% Calculate potential energy for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:21
% EndTime: 2019-03-08 21:45:22
% DurationCPUTime: 0.40s
% Computational Cost: add. (322->116), mult. (650->137), div. (0->0), fcn. (780->10), ass. (0->55)
t72 = pkin(4) + pkin(8);
t36 = sin(pkin(6));
t71 = pkin(7) * t36;
t70 = rSges(5,1) + pkin(8);
t69 = rSges(7,1) + pkin(5);
t68 = rSges(7,2) + pkin(9);
t67 = rSges(4,3) + pkin(8);
t66 = rSges(6,3) + pkin(9);
t65 = cos(qJ(3));
t39 = sin(qJ(3));
t64 = t36 * t39;
t40 = sin(qJ(2));
t63 = t36 * t40;
t42 = cos(qJ(2));
t62 = t36 * t42;
t61 = rSges(5,3) + qJ(4);
t60 = rSges(7,3) + qJ(6);
t59 = cos(pkin(6));
t35 = sin(pkin(10));
t58 = pkin(1) * t35 + r_base(2);
t57 = qJ(1) + r_base(3);
t56 = t36 * t65;
t55 = t40 * t59;
t54 = t42 * t59;
t37 = cos(pkin(10));
t53 = pkin(1) * t37 + t35 * t71 + r_base(1);
t52 = pkin(7) * t59 + t57;
t23 = -t35 * t55 + t37 * t42;
t51 = pkin(2) * t23 + t53;
t50 = pkin(2) * t63 + t52;
t12 = t23 * t65 + t35 * t64;
t49 = pkin(3) * t12 + t51;
t25 = t39 * t59 + t40 * t56;
t48 = pkin(3) * t25 + t50;
t21 = t35 * t42 + t37 * t55;
t47 = pkin(2) * t21 - t37 * t71 + t58;
t10 = t21 * t65 - t37 * t64;
t46 = pkin(3) * t10 + t47;
t11 = t23 * t39 - t35 * t56;
t22 = t35 * t54 + t37 * t40;
t45 = t11 * qJ(4) + t22 * t72 + t49;
t24 = t39 * t63 - t59 * t65;
t44 = t24 * qJ(4) - t62 * t72 + t48;
t20 = t35 * t40 - t37 * t54;
t9 = t21 * t39 + t37 * t56;
t43 = t9 * qJ(4) + t20 * t72 + t46;
t41 = cos(qJ(5));
t38 = sin(qJ(5));
t14 = t24 * t38 - t41 * t62;
t13 = t24 * t41 + t38 * t62;
t4 = t11 * t38 + t22 * t41;
t3 = -t11 * t41 + t22 * t38;
t2 = t20 * t41 + t38 * t9;
t1 = t20 * t38 - t41 * t9;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t37 - rSges(2,2) * t35 + r_base(1)) + g(2) * (rSges(2,1) * t35 + rSges(2,2) * t37 + r_base(2)) + g(3) * (rSges(2,3) + t57)) - m(3) * (g(1) * (rSges(3,1) * t23 - rSges(3,2) * t22 + t53) + g(2) * (t21 * rSges(3,1) - t20 * rSges(3,2) + t58) + g(3) * (t59 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t42) + g(2) * (-rSges(3,3) - pkin(7)) * t37) * t36) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t22 * t67 + t51) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t20 * t67 + t47) + g(3) * (t25 * rSges(4,1) - t24 * rSges(4,2) - t62 * t67 + t50)) - m(5) * (g(1) * (-t12 * rSges(5,2) + t11 * t61 + t22 * t70 + t49) + g(2) * (-t10 * rSges(5,2) + t20 * t70 + t61 * t9 + t46) + g(3) * (-t25 * rSges(5,2) + t24 * t61 - t62 * t70 + t48)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t12 * t66 + t45) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t10 * t66 + t43) + g(3) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t25 * t66 + t44)) - m(7) * (g(1) * (t12 * t68 + t3 * t60 + t4 * t69 + t45) + g(2) * (t1 * t60 + t10 * t68 + t2 * t69 + t43) + g(3) * (-t13 * t60 + t14 * t69 + t25 * t68 + t44));
U  = t5;
