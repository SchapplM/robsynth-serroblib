% Calculate potential energy for
% S6PRRPRP4
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:10
% EndTime: 2019-03-08 21:40:10
% DurationCPUTime: 0.38s
% Computational Cost: add. (310->122), mult. (612->145), div. (0->0), fcn. (727->10), ass. (0->52)
t67 = pkin(4) + pkin(8);
t33 = sin(pkin(6));
t66 = pkin(7) * t33;
t65 = rSges(5,1) + pkin(8);
t64 = rSges(4,3) + pkin(8);
t63 = rSges(6,3) + pkin(9);
t39 = cos(qJ(5));
t62 = t39 * pkin(5) + t67;
t61 = cos(qJ(3));
t37 = sin(qJ(3));
t60 = t33 * t37;
t38 = sin(qJ(2));
t59 = t33 * t38;
t40 = cos(qJ(2));
t58 = t33 * t40;
t57 = rSges(7,3) + qJ(6) + pkin(9);
t56 = rSges(5,3) + qJ(4);
t55 = cos(pkin(6));
t32 = sin(pkin(10));
t54 = t32 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t52 = t33 * t61;
t36 = sin(qJ(5));
t51 = pkin(5) * t36 + qJ(4);
t50 = t38 * t55;
t49 = t40 * t55;
t34 = cos(pkin(10));
t48 = t34 * pkin(1) + t32 * t66 + r_base(1);
t47 = t55 * pkin(7) + t53;
t19 = -t32 * t50 + t34 * t40;
t46 = t19 * pkin(2) + t48;
t45 = pkin(2) * t59 + t47;
t10 = t19 * t61 + t32 * t60;
t44 = t10 * pkin(3) + t46;
t21 = t55 * t37 + t38 * t52;
t43 = t21 * pkin(3) + t45;
t17 = t32 * t40 + t34 * t50;
t42 = t17 * pkin(2) - t34 * t66 + t54;
t8 = t17 * t61 - t34 * t60;
t41 = t8 * pkin(3) + t42;
t20 = t37 * t59 - t55 * t61;
t18 = t32 * t49 + t34 * t38;
t16 = t32 * t38 - t34 * t49;
t12 = t20 * t36 - t39 * t58;
t11 = t20 * t39 + t36 * t58;
t9 = t19 * t37 - t32 * t52;
t7 = t17 * t37 + t34 * t52;
t4 = t18 * t39 + t9 * t36;
t3 = -t18 * t36 + t9 * t39;
t2 = t16 * t39 + t7 * t36;
t1 = -t16 * t36 + t7 * t39;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t19 * rSges(3,1) - t18 * rSges(3,2) + t48) + g(2) * (t17 * rSges(3,1) - t16 * rSges(3,2) + t54) + g(3) * (t55 * rSges(3,3) + t47) + (g(1) * rSges(3,3) * t32 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t40) + g(2) * (-rSges(3,3) - pkin(7)) * t34) * t33) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t64 * t18 + t46) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t64 * t16 + t42) + g(3) * (t21 * rSges(4,1) - t20 * rSges(4,2) - t64 * t58 + t45)) - m(5) * (g(1) * (-t10 * rSges(5,2) + t65 * t18 + t56 * t9 + t44) + g(2) * (-t8 * rSges(5,2) + t65 * t16 + t56 * t7 + t41) + g(3) * (-t21 * rSges(5,2) + t56 * t20 - t65 * t58 + t43)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t9 * qJ(4) + t63 * t10 + t67 * t18 + t44) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t7 * qJ(4) + t67 * t16 + t63 * t8 + t41) + g(3) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t20 * qJ(4) + t63 * t21 - t67 * t58 + t43)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t57 * t10 + t62 * t18 + t51 * t9 + t44) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t62 * t16 + t51 * t7 + t57 * t8 + t41) + g(3) * (t12 * rSges(7,1) + t11 * rSges(7,2) + t51 * t20 + t57 * t21 - t62 * t58 + t43));
U  = t5;
