% Calculate potential energy for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:08
% EndTime: 2019-03-08 22:30:08
% DurationCPUTime: 0.38s
% Computational Cost: add. (322->113), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->51)
t66 = pkin(4) + pkin(8);
t30 = sin(pkin(6));
t65 = pkin(7) * t30;
t64 = rSges(5,1) + pkin(8);
t63 = rSges(4,3) + pkin(8);
t62 = pkin(9) + rSges(6,3);
t61 = cos(qJ(3));
t33 = sin(qJ(3));
t60 = t30 * t33;
t34 = sin(qJ(2));
t59 = t30 * t34;
t36 = cos(qJ(2));
t58 = t30 * t36;
t57 = pkin(10) + pkin(9) + rSges(7,3);
t56 = rSges(5,3) + qJ(4);
t55 = cos(pkin(6));
t29 = sin(pkin(11));
t54 = t29 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t52 = t30 * t61;
t51 = t34 * t55;
t50 = t36 * t55;
t31 = cos(pkin(11));
t49 = t31 * pkin(1) + t29 * t65 + r_base(1);
t48 = t55 * pkin(7) + t53;
t13 = -t29 * t51 + t31 * t36;
t47 = t13 * pkin(2) + t49;
t46 = pkin(2) * t59 + t48;
t6 = t13 * t61 + t29 * t60;
t45 = t6 * pkin(3) + t47;
t15 = t55 * t33 + t34 * t52;
t44 = t15 * pkin(3) + t46;
t32 = sin(qJ(5));
t35 = cos(qJ(5));
t43 = t32 * rSges(6,1) + t35 * rSges(6,2) + qJ(4);
t42 = t35 * rSges(6,1) - t32 * rSges(6,2) + t66;
t28 = qJ(5) + qJ(6);
t23 = sin(t28);
t24 = cos(t28);
t41 = t24 * rSges(7,1) - t23 * rSges(7,2) + t35 * pkin(5) + t66;
t11 = t29 * t36 + t31 * t51;
t40 = t11 * pkin(2) - t31 * t65 + t54;
t4 = t11 * t61 - t31 * t60;
t39 = t4 * pkin(3) + t40;
t38 = t23 * rSges(7,1) + t24 * rSges(7,2) + t32 * pkin(5) + qJ(4);
t14 = t33 * t59 - t55 * t61;
t12 = t29 * t50 + t31 * t34;
t10 = t29 * t34 - t31 * t50;
t5 = t13 * t33 - t29 * t52;
t3 = t11 * t33 + t31 * t52;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t31 * rSges(2,1) - t29 * rSges(2,2) + r_base(1)) + g(2) * (t29 * rSges(2,1) + t31 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t13 * rSges(3,1) - t12 * rSges(3,2) + t49) + g(2) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t54) + g(3) * (t55 * rSges(3,3) + t48) + (g(1) * rSges(3,3) * t29 + g(3) * (rSges(3,1) * t34 + rSges(3,2) * t36) + g(2) * (-rSges(3,3) - pkin(7)) * t31) * t30) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t63 * t12 + t47) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t63 * t10 + t40) + g(3) * (t15 * rSges(4,1) - t14 * rSges(4,2) - t63 * t58 + t46)) - m(5) * (g(1) * (-t6 * rSges(5,2) + t64 * t12 + t56 * t5 + t45) + g(2) * (-t4 * rSges(5,2) + t64 * t10 + t56 * t3 + t39) + g(3) * (-t15 * rSges(5,2) + t56 * t14 - t64 * t58 + t44)) - m(6) * (g(1) * (t42 * t12 + t43 * t5 + t62 * t6 + t45) + g(2) * (t42 * t10 + t43 * t3 + t62 * t4 + t39) + g(3) * (t43 * t14 + t62 * t15 - t42 * t58 + t44)) - m(7) * (g(1) * (t41 * t12 + t38 * t5 + t57 * t6 + t45) + g(2) * (t41 * t10 + t38 * t3 + t57 * t4 + t39) + g(3) * (t38 * t14 + t57 * t15 - t41 * t58 + t44));
U  = t1;
