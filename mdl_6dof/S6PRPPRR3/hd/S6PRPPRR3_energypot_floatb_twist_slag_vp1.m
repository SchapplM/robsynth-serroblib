% Calculate potential energy for
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:20
% EndTime: 2019-03-08 19:21:21
% DurationCPUTime: 0.49s
% Computational Cost: add. (333->130), mult. (677->165), div. (0->0), fcn. (818->12), ass. (0->54)
t37 = sin(pkin(6));
t71 = pkin(7) * t37;
t36 = sin(pkin(10));
t70 = g(1) * t36;
t39 = cos(pkin(10));
t69 = g(2) * t39;
t68 = rSges(6,3) + pkin(8);
t67 = pkin(9) + rSges(7,3);
t42 = sin(qJ(5));
t66 = t37 * t42;
t43 = sin(qJ(2));
t65 = t37 * t43;
t45 = cos(qJ(5));
t64 = t37 * t45;
t40 = cos(pkin(6));
t63 = t40 * t43;
t46 = cos(qJ(2));
t62 = t40 * t46;
t61 = qJ(3) * t46;
t60 = qJ(4) * t37;
t59 = t36 * pkin(1) + r_base(2);
t58 = qJ(1) + r_base(3);
t57 = t39 * pkin(1) + t36 * t71 + r_base(1);
t56 = t40 * pkin(7) + t58;
t55 = pkin(2) * t65 + t56;
t22 = t36 * t43 - t39 * t62;
t23 = t36 * t46 + t39 * t63;
t54 = t23 * pkin(2) + t22 * qJ(3) + t59;
t24 = t36 * t62 + t39 * t43;
t25 = -t36 * t63 + t39 * t46;
t53 = t25 * pkin(2) + t24 * qJ(3) + t57;
t52 = t23 * pkin(3) + t39 * t60 + t54;
t51 = t25 * pkin(3) + t53;
t50 = pkin(3) * t65 - t40 * qJ(4) + t55;
t35 = sin(pkin(11));
t38 = cos(pkin(11));
t8 = t22 * t35 + t23 * t38;
t49 = t8 * pkin(4) - t39 * t71 + t52;
t10 = t24 * t35 + t25 * t38;
t48 = t10 * pkin(4) - t36 * t60 + t51;
t17 = (-t35 * t46 + t38 * t43) * t37;
t47 = t17 * pkin(4) - t37 * t61 + t50;
t44 = cos(qJ(6));
t41 = sin(qJ(6));
t16 = (t35 * t43 + t38 * t46) * t37;
t12 = t17 * t45 - t40 * t42;
t11 = t17 * t42 + t40 * t45;
t9 = -t24 * t38 + t25 * t35;
t7 = -t22 * t38 + t23 * t35;
t4 = t10 * t45 - t36 * t66;
t3 = t10 * t42 + t36 * t64;
t2 = t39 * t66 + t45 * t8;
t1 = -t39 * t64 + t42 * t8;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t39 - rSges(2,2) * t36 + r_base(1)) + g(2) * (rSges(2,1) * t36 + rSges(2,2) * t39 + r_base(2)) + g(3) * (rSges(2,3) + t58)) - m(3) * (g(1) * (rSges(3,1) * t25 - rSges(3,2) * t24 + t57) + g(2) * (rSges(3,1) * t23 - rSges(3,2) * t22 + t59) + g(3) * (t40 * rSges(3,3) + t56) + (rSges(3,3) * t70 + g(3) * (rSges(3,1) * t43 + rSges(3,2) * t46) + (-rSges(3,3) - pkin(7)) * t69) * t37) - m(4) * (g(1) * (rSges(4,1) * t25 + rSges(4,3) * t24 + t53) + g(2) * (rSges(4,1) * t23 + rSges(4,3) * t22 + t54) + g(3) * (t40 * rSges(4,2) + t55) + (rSges(4,2) * t70 + g(3) * (rSges(4,1) * t43 - rSges(4,3) * t46 - t61) + (-rSges(4,2) - pkin(7)) * t69) * t37) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t51) + g(2) * (rSges(5,1) * t8 - rSges(5,2) * t7 + t52) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) - t40 * rSges(5,3) + t50) + (-g(3) * t61 + (rSges(5,3) - pkin(7)) * t69 + (-rSges(5,3) - qJ(4)) * t70) * t37) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t68 * t9 + t48) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t68 * t7 + t49) + g(3) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t68 * t16 + t47)) - m(7) * (g(1) * (t4 * pkin(5) + t9 * pkin(8) + (t4 * t44 + t41 * t9) * rSges(7,1) + (-t4 * t41 + t44 * t9) * rSges(7,2) + t67 * t3 + t48) + g(2) * (t2 * pkin(5) + t7 * pkin(8) + (t2 * t44 + t41 * t7) * rSges(7,1) + (-t2 * t41 + t44 * t7) * rSges(7,2) + t67 * t1 + t49) + g(3) * (t12 * pkin(5) + t16 * pkin(8) + (t12 * t44 + t16 * t41) * rSges(7,1) + (-t12 * t41 + t16 * t44) * rSges(7,2) + t67 * t11 + t47));
U  = t5;
