% Calculate potential energy for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:20
% EndTime: 2019-03-08 19:34:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (393->119), mult. (818->144), div. (0->0), fcn. (1011->12), ass. (0->52)
t35 = sin(pkin(11));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t61 = cos(pkin(11));
t24 = -t42 * t35 + t45 * t61;
t71 = rSges(6,1) + pkin(8);
t69 = rSges(5,3) + pkin(8);
t68 = pkin(9) + rSges(7,3);
t36 = sin(pkin(10));
t37 = sin(pkin(6));
t67 = t36 * t37;
t38 = cos(pkin(10));
t66 = t38 * t37;
t39 = cos(pkin(6));
t65 = t39 * t42;
t64 = t39 * t45;
t62 = rSges(6,3) + qJ(5);
t32 = t45 * pkin(2) + pkin(1);
t60 = t38 * t32 + r_base(1);
t59 = qJ(1) + r_base(3);
t22 = pkin(2) * t65 + (-pkin(7) - qJ(3)) * t37;
t57 = t38 * t22 + t36 * t32 + r_base(2);
t56 = t39 * pkin(7) + t59;
t23 = -t45 * t35 - t42 * t61;
t21 = t23 * t39;
t10 = -t38 * t21 + t36 * t24;
t55 = t10 * pkin(3) + t57;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t4 = t10 * t44 - t41 * t66;
t54 = t4 * pkin(4) + t55;
t53 = t37 * t42 * pkin(2) + t39 * qJ(3) + t56;
t12 = t36 * t21 + t38 * t24;
t52 = t12 * pkin(3) - t36 * t22 + t60;
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t51 = t40 * rSges(7,1) + t43 * rSges(7,2) + qJ(5);
t20 = t23 * t37;
t50 = -t20 * pkin(3) + t53;
t49 = t43 * rSges(7,1) - t40 * rSges(7,2) + pkin(5) + pkin(8);
t6 = t12 * t44 + t41 * t67;
t48 = t6 * pkin(4) + t52;
t15 = -t20 * t44 + t39 * t41;
t47 = t15 * pkin(4) + t50;
t46 = t24 * t39;
t19 = t24 * t37;
t14 = -t20 * t41 - t39 * t44;
t11 = t38 * t23 - t36 * t46;
t9 = t36 * t23 + t38 * t46;
t5 = t12 * t41 - t44 * t67;
t3 = t10 * t41 + t44 * t66;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t36 * rSges(2,2) + r_base(1)) + g(2) * (t36 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t59)) - m(3) * (g(1) * (t38 * pkin(1) + r_base(1) + (-t36 * t65 + t38 * t45) * rSges(3,1) + (-t36 * t64 - t38 * t42) * rSges(3,2)) + g(2) * (t36 * pkin(1) + r_base(2) + (t36 * t45 + t38 * t65) * rSges(3,1) + (-t36 * t42 + t38 * t64) * rSges(3,2)) + g(3) * (t39 * rSges(3,3) + t56) + (g(3) * (rSges(3,1) * t42 + rSges(3,2) * t45) + (g(1) * t36 - g(2) * t38) * (rSges(3,3) + pkin(7))) * t37) - m(4) * (g(1) * (t12 * rSges(4,1) + t11 * rSges(4,2) + (rSges(4,3) * t37 - t22) * t36 + t60) + g(2) * (t10 * rSges(4,1) + t9 * rSges(4,2) - rSges(4,3) * t66 + t57) + g(3) * (-t20 * rSges(4,1) + t19 * rSges(4,2) + t39 * rSges(4,3) + t53)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) - t11 * t69 + t52) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t69 * t9 + t55) + g(3) * (t15 * rSges(5,1) - t14 * rSges(5,2) - t19 * t69 + t50)) - m(6) * (g(1) * (-t6 * rSges(6,2) - t11 * t71 + t62 * t5 + t48) + g(2) * (-t4 * rSges(6,2) + t62 * t3 - t71 * t9 + t54) + g(3) * (-t15 * rSges(6,2) + t62 * t14 - t19 * t71 + t47)) - m(7) * (g(1) * (-t11 * t49 + t5 * t51 + t6 * t68 + t48) + g(2) * (t3 * t51 + t4 * t68 - t49 * t9 + t54) + g(3) * (t14 * t51 + t15 * t68 - t19 * t49 + t47));
U  = t1;
