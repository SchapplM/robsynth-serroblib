% Calculate potential energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:25
% EndTime: 2019-03-08 19:29:26
% DurationCPUTime: 0.60s
% Computational Cost: add. (421->129), mult. (853->161), div. (0->0), fcn. (1059->14), ass. (0->51)
t36 = sin(pkin(11));
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t59 = cos(pkin(11));
t21 = -t44 * t36 + t46 * t59;
t67 = rSges(5,3) + pkin(8);
t37 = sin(pkin(10));
t38 = sin(pkin(6));
t66 = t37 * t38;
t40 = cos(pkin(10));
t65 = t40 * t38;
t41 = cos(pkin(6));
t64 = t41 * t44;
t63 = t41 * t46;
t61 = pkin(9) + qJ(5) + rSges(7,3);
t60 = qJ(5) + rSges(6,3);
t29 = t46 * pkin(2) + pkin(1);
t58 = t40 * t29 + r_base(1);
t57 = qJ(1) + r_base(3);
t19 = pkin(2) * t64 + (-pkin(7) - qJ(3)) * t38;
t55 = t40 * t19 + t37 * t29 + r_base(2);
t54 = t41 * pkin(7) + t57;
t20 = -t46 * t36 - t44 * t59;
t18 = t20 * t41;
t8 = -t40 * t18 + t37 * t21;
t53 = t8 * pkin(3) + t55;
t52 = t38 * t44 * pkin(2) + t41 * qJ(3) + t54;
t34 = pkin(12) + qJ(6);
t30 = sin(t34);
t31 = cos(t34);
t39 = cos(pkin(12));
t51 = t31 * rSges(7,1) - t30 * rSges(7,2) + t39 * pkin(5) + pkin(4);
t10 = t37 * t18 + t40 * t21;
t50 = t10 * pkin(3) - t37 * t19 + t58;
t17 = t20 * t38;
t49 = -t17 * pkin(3) + t52;
t35 = sin(pkin(12));
t48 = t30 * rSges(7,1) + t31 * rSges(7,2) + t35 * pkin(5) + pkin(8);
t47 = t21 * t41;
t45 = cos(qJ(4));
t43 = sin(qJ(4));
t16 = t21 * t38;
t12 = -t17 * t45 + t41 * t43;
t11 = -t17 * t43 - t41 * t45;
t9 = t40 * t20 - t37 * t47;
t7 = t37 * t20 + t40 * t47;
t4 = t10 * t45 + t43 * t66;
t3 = t10 * t43 - t45 * t66;
t2 = -t43 * t65 + t8 * t45;
t1 = t8 * t43 + t45 * t65;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t57)) - m(3) * (g(1) * (t40 * pkin(1) + r_base(1) + (-t37 * t64 + t40 * t46) * rSges(3,1) + (-t37 * t63 - t40 * t44) * rSges(3,2)) + g(2) * (t37 * pkin(1) + r_base(2) + (t37 * t46 + t40 * t64) * rSges(3,1) + (-t37 * t44 + t40 * t63) * rSges(3,2)) + g(3) * (t41 * rSges(3,3) + t54) + (g(3) * (rSges(3,1) * t44 + rSges(3,2) * t46) + (g(1) * t37 - g(2) * t40) * (rSges(3,3) + pkin(7))) * t38) - m(4) * (g(1) * (t10 * rSges(4,1) + t9 * rSges(4,2) + (rSges(4,3) * t38 - t19) * t37 + t58) + g(2) * (t8 * rSges(4,1) + t7 * rSges(4,2) - rSges(4,3) * t65 + t55) + g(3) * (-t17 * rSges(4,1) + t16 * rSges(4,2) + t41 * rSges(4,3) + t52)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t67 * t9 + t50) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) - t67 * t7 + t53) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) - t67 * t16 + t49)) - m(6) * (g(1) * (t4 * pkin(4) - t9 * pkin(8) + (-t9 * t35 + t4 * t39) * rSges(6,1) + (-t4 * t35 - t9 * t39) * rSges(6,2) + t60 * t3 + t50) + g(2) * (t2 * pkin(4) - t7 * pkin(8) + (t2 * t39 - t7 * t35) * rSges(6,1) + (-t2 * t35 - t7 * t39) * rSges(6,2) + t60 * t1 + t53) + g(3) * (t12 * pkin(4) - t16 * pkin(8) + (t12 * t39 - t16 * t35) * rSges(6,1) + (-t12 * t35 - t16 * t39) * rSges(6,2) + t60 * t11 + t49)) - m(7) * (g(1) * (t61 * t3 + t51 * t4 - t48 * t9 + t50) + g(2) * (t61 * t1 + t51 * t2 - t48 * t7 + t53) + g(3) * (t61 * t11 + t51 * t12 - t48 * t16 + t49));
U  = t5;
