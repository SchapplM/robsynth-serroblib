% Calculate potential energy for
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:11
% EndTime: 2019-03-08 22:16:12
% DurationCPUTime: 0.45s
% Computational Cost: add. (356->124), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->52)
t32 = sin(pkin(12));
t66 = t32 * pkin(4) + pkin(8);
t34 = sin(pkin(6));
t65 = pkin(7) * t34;
t63 = rSges(4,3) + pkin(8);
t35 = cos(pkin(12));
t22 = t35 * pkin(4) + pkin(3);
t62 = cos(qJ(3));
t38 = sin(qJ(3));
t61 = t34 * t38;
t39 = sin(qJ(2));
t60 = t34 * t39;
t40 = cos(qJ(2));
t59 = t34 * t40;
t37 = -pkin(9) - qJ(4);
t58 = pkin(10) - t37 + rSges(7,3);
t57 = -t37 + rSges(6,3);
t56 = qJ(4) + rSges(5,3);
t55 = cos(pkin(6));
t31 = pkin(12) + qJ(5);
t33 = sin(pkin(11));
t54 = t33 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t52 = t34 * t62;
t51 = t39 * t55;
t50 = t40 * t55;
t36 = cos(pkin(11));
t49 = t36 * pkin(1) + t33 * t65 + r_base(1);
t48 = t55 * pkin(7) + t53;
t10 = -t33 * t51 + t36 * t40;
t47 = t10 * pkin(2) + t49;
t46 = pkin(2) * t60 + t48;
t23 = sin(t31);
t24 = cos(t31);
t45 = t24 * rSges(6,1) - t23 * rSges(6,2) + t22;
t25 = qJ(6) + t31;
t20 = sin(t25);
t21 = cos(t25);
t44 = t21 * rSges(7,1) - t20 * rSges(7,2) + pkin(5) * t24 + t22;
t43 = t20 * rSges(7,1) + t21 * rSges(7,2) + pkin(5) * t23 + t66;
t8 = t33 * t40 + t36 * t51;
t42 = t8 * pkin(2) - t36 * t65 + t54;
t41 = t23 * rSges(6,1) + t24 * rSges(6,2) + t66;
t12 = t55 * t38 + t39 * t52;
t11 = t38 * t60 - t55 * t62;
t9 = t33 * t50 + t36 * t39;
t7 = t33 * t39 - t36 * t50;
t4 = t10 * t62 + t33 * t61;
t3 = t10 * t38 - t33 * t52;
t2 = -t36 * t61 + t8 * t62;
t1 = t36 * t52 + t8 * t38;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t36 * rSges(2,1) - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + t36 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t49) + g(2) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t54) + g(3) * (t55 * rSges(3,3) + t48) + (g(1) * rSges(3,3) * t33 + g(3) * (rSges(3,1) * t39 + rSges(3,2) * t40) + g(2) * (-rSges(3,3) - pkin(7)) * t36) * t34) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t63 * t9 + t47) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t63 * t7 + t42) + g(3) * (t12 * rSges(4,1) - t11 * rSges(4,2) - t63 * t59 + t46)) - m(5) * (g(1) * (t4 * pkin(3) + t9 * pkin(8) + (t9 * t32 + t4 * t35) * rSges(5,1) + (-t4 * t32 + t9 * t35) * rSges(5,2) + t56 * t3 + t47) + g(2) * (t2 * pkin(3) + t7 * pkin(8) + (t2 * t35 + t7 * t32) * rSges(5,1) + (-t2 * t32 + t7 * t35) * rSges(5,2) + t56 * t1 + t42) + g(3) * (t12 * pkin(3) - pkin(8) * t59 + (t12 * t35 - t32 * t59) * rSges(5,1) + (-t12 * t32 - t35 * t59) * rSges(5,2) + t56 * t11 + t46)) - m(6) * (g(1) * (t57 * t3 + t45 * t4 + t41 * t9 + t47) + g(2) * (t57 * t1 + t45 * t2 + t41 * t7 + t42) + g(3) * (t57 * t11 + t45 * t12 - t41 * t59 + t46)) - m(7) * (g(1) * (t58 * t3 + t44 * t4 + t43 * t9 + t47) + g(2) * (t58 * t1 + t44 * t2 + t43 * t7 + t42) + g(3) * (t58 * t11 + t44 * t12 - t43 * t59 + t46));
U  = t5;
