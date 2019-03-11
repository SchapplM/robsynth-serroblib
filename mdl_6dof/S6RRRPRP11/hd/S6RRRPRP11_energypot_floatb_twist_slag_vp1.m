% Calculate potential energy for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:28
% EndTime: 2019-03-09 17:41:28
% DurationCPUTime: 0.38s
% Computational Cost: add. (310->122), mult. (612->145), div. (0->0), fcn. (727->10), ass. (0->52)
t67 = pkin(4) + pkin(9);
t66 = rSges(5,1) + pkin(9);
t65 = rSges(4,3) + pkin(9);
t64 = rSges(6,3) + pkin(10);
t38 = cos(qJ(5));
t63 = t38 * pkin(5) + t67;
t62 = cos(qJ(3));
t32 = sin(pkin(6));
t36 = sin(qJ(2));
t61 = t32 * t36;
t37 = sin(qJ(1));
t60 = t32 * t37;
t39 = cos(qJ(2));
t59 = t32 * t39;
t40 = cos(qJ(1));
t58 = t32 * t40;
t57 = rSges(7,3) + qJ(6) + pkin(10);
t56 = rSges(5,3) + qJ(4);
t55 = cos(pkin(6));
t54 = pkin(7) + r_base(3);
t53 = t37 * pkin(1) + r_base(2);
t52 = t32 * t62;
t51 = t55 * pkin(8) + t54;
t34 = sin(qJ(5));
t50 = pkin(5) * t34 + qJ(4);
t49 = t37 * t55;
t48 = t40 * t55;
t47 = t40 * pkin(1) + pkin(8) * t60 + r_base(1);
t46 = pkin(2) * t61 + t51;
t21 = -t36 * t49 + t40 * t39;
t45 = t21 * pkin(2) + t47;
t35 = sin(qJ(3));
t17 = t55 * t35 + t36 * t52;
t44 = t17 * pkin(3) + t46;
t12 = t21 * t62 + t35 * t60;
t43 = t12 * pkin(3) + t45;
t19 = t36 * t48 + t37 * t39;
t42 = t19 * pkin(2) - pkin(8) * t58 + t53;
t10 = t19 * t62 - t35 * t58;
t41 = t10 * pkin(3) + t42;
t20 = t40 * t36 + t39 * t49;
t18 = t37 * t36 - t39 * t48;
t16 = t35 * t61 - t55 * t62;
t11 = t21 * t35 - t37 * t52;
t9 = t19 * t35 + t40 * t52;
t8 = t16 * t34 - t38 * t59;
t7 = t16 * t38 + t34 * t59;
t4 = t11 * t34 + t20 * t38;
t3 = t11 * t38 - t20 * t34;
t2 = t18 * t38 + t9 * t34;
t1 = -t18 * t34 + t9 * t38;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t21 * rSges(3,1) - t20 * rSges(3,2) + t47) + g(2) * (t19 * rSges(3,1) - t18 * rSges(3,2) + t53) + g(3) * (t55 * rSges(3,3) + t51) + (g(1) * rSges(3,3) * t37 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t39) + g(2) * (-rSges(3,3) - pkin(8)) * t40) * t32) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t65 * t20 + t45) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t65 * t18 + t42) + g(3) * (t17 * rSges(4,1) - t16 * rSges(4,2) - t65 * t59 + t46)) - m(5) * (g(1) * (-t12 * rSges(5,2) + t56 * t11 + t66 * t20 + t43) + g(2) * (-t10 * rSges(5,2) + t66 * t18 + t56 * t9 + t41) + g(3) * (-t17 * rSges(5,2) + t56 * t16 - t66 * t59 + t44)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t11 * qJ(4) + t64 * t12 + t67 * t20 + t43) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t9 * qJ(4) + t64 * t10 + t67 * t18 + t41) + g(3) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t16 * qJ(4) + t64 * t17 - t67 * t59 + t44)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t50 * t11 + t57 * t12 + t63 * t20 + t43) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t57 * t10 + t63 * t18 + t50 * t9 + t41) + g(3) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t50 * t16 + t57 * t17 - t63 * t59 + t44));
U  = t5;
