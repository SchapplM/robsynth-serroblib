% Calculate potential energy for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:04
% EndTime: 2019-03-09 09:01:05
% DurationCPUTime: 0.57s
% Computational Cost: add. (357->129), mult. (724->160), div. (0->0), fcn. (881->12), ass. (0->52)
t33 = sin(pkin(11));
t35 = cos(pkin(11));
t39 = sin(qJ(2));
t43 = cos(qJ(2));
t23 = -t39 * t33 + t43 * t35;
t71 = pkin(2) * t39;
t69 = rSges(6,3) + pkin(9);
t68 = pkin(10) + rSges(7,3);
t34 = sin(pkin(6));
t40 = sin(qJ(1));
t66 = t40 * t34;
t65 = t40 * t39;
t64 = t40 * t43;
t44 = cos(qJ(1));
t62 = t44 * t34;
t61 = t44 * t39;
t60 = t44 * t43;
t59 = rSges(5,3) + qJ(4);
t58 = pkin(7) + r_base(3);
t29 = t43 * pkin(2) + pkin(1);
t57 = t44 * t29 + r_base(1);
t36 = cos(pkin(6));
t56 = t36 * pkin(8) + t58;
t52 = t43 * t33 + t39 * t35;
t49 = t52 * t36;
t11 = t44 * t23 - t40 * t49;
t55 = t11 * pkin(3) + t57;
t21 = t36 * t71 + (-pkin(8) - qJ(3)) * t34;
t54 = t44 * t21 + t40 * t29 + r_base(2);
t9 = t40 * t23 + t44 * t49;
t53 = t9 * pkin(3) + t54;
t51 = t36 * qJ(3) + t34 * t71 + t56;
t20 = t52 * t34;
t50 = t20 * pkin(3) + t51;
t48 = t23 * t36;
t10 = -t40 * t48 - t44 * t52;
t47 = pkin(4) * t66 - t10 * qJ(4) - t40 * t21 + t55;
t19 = t23 * t34;
t46 = t36 * pkin(4) - t19 * qJ(4) + t50;
t8 = -t40 * t52 + t44 * t48;
t45 = -pkin(4) * t62 - t8 * qJ(4) + t53;
t42 = cos(qJ(5));
t41 = cos(qJ(6));
t38 = sin(qJ(5));
t37 = sin(qJ(6));
t13 = -t19 * t38 + t36 * t42;
t12 = t19 * t42 + t36 * t38;
t4 = -t8 * t38 - t42 * t62;
t3 = t38 * t62 - t8 * t42;
t2 = -t10 * t38 + t42 * t66;
t1 = t10 * t42 + t38 * t66;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t44 * rSges(2,1) - t40 * rSges(2,2) + r_base(1)) + g(2) * (t40 * rSges(2,1) + t44 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t58)) - m(3) * (g(1) * (t44 * pkin(1) + r_base(1) + (-t36 * t65 + t60) * rSges(3,1) + (-t36 * t64 - t61) * rSges(3,2)) + g(2) * (t40 * pkin(1) + r_base(2) + (t36 * t61 + t64) * rSges(3,1) + (t36 * t60 - t65) * rSges(3,2)) + g(3) * (t36 * rSges(3,3) + t56) + (g(3) * (rSges(3,1) * t39 + rSges(3,2) * t43) + (g(1) * t40 - g(2) * t44) * (rSges(3,3) + pkin(8))) * t34) - m(4) * (g(1) * (t11 * rSges(4,1) + t10 * rSges(4,2) + (rSges(4,3) * t34 - t21) * t40 + t57) + g(2) * (t9 * rSges(4,1) + t8 * rSges(4,2) - rSges(4,3) * t62 + t54) + g(3) * (t20 * rSges(4,1) + t19 * rSges(4,2) + t36 * rSges(4,3) + t51)) - m(5) * (g(1) * (-t11 * rSges(5,2) + (rSges(5,1) * t34 - t21) * t40 - t59 * t10 + t55) + g(2) * (-rSges(5,1) * t62 - t9 * rSges(5,2) - t59 * t8 + t53) + g(3) * (t36 * rSges(5,1) - t20 * rSges(5,2) - t59 * t19 + t50)) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t69 * t11 + t47) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t69 * t9 + t45) + g(3) * (t13 * rSges(6,1) - t12 * rSges(6,2) + t20 * t69 + t46)) - m(7) * (g(1) * (t2 * pkin(5) + t11 * pkin(9) + (t11 * t37 + t2 * t41) * rSges(7,1) + (t11 * t41 - t2 * t37) * rSges(7,2) + t68 * t1 + t47) + g(2) * (t4 * pkin(5) + t9 * pkin(9) + (t9 * t37 + t4 * t41) * rSges(7,1) + (-t4 * t37 + t9 * t41) * rSges(7,2) - t68 * t3 + t45) + g(3) * (t13 * pkin(5) + t20 * pkin(9) + (t13 * t41 + t20 * t37) * rSges(7,1) + (-t13 * t37 + t20 * t41) * rSges(7,2) + t68 * t12 + t46));
U  = t5;
