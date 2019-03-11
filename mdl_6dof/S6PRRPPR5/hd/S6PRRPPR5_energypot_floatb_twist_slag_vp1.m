% Calculate potential energy for
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:04
% EndTime: 2019-03-08 21:18:05
% DurationCPUTime: 0.38s
% Computational Cost: add. (322->113), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->51)
t66 = pkin(4) + pkin(8);
t31 = sin(pkin(6));
t65 = pkin(7) * t31;
t64 = rSges(5,1) + pkin(8);
t63 = rSges(4,3) + pkin(8);
t62 = cos(qJ(3));
t35 = sin(qJ(3));
t61 = t31 * t35;
t36 = sin(qJ(2));
t60 = t31 * t36;
t37 = cos(qJ(2));
t59 = t31 * t37;
t58 = pkin(9) + qJ(5) + rSges(7,3);
t57 = rSges(5,3) + qJ(4);
t56 = qJ(5) + rSges(6,3);
t55 = cos(pkin(6));
t30 = sin(pkin(10));
t54 = t30 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t52 = t31 * t62;
t51 = t36 * t55;
t50 = t37 * t55;
t33 = cos(pkin(10));
t49 = t33 * pkin(1) + t30 * t65 + r_base(1);
t48 = t55 * pkin(7) + t53;
t13 = -t30 * t51 + t33 * t37;
t47 = t13 * pkin(2) + t49;
t46 = pkin(2) * t60 + t48;
t6 = t13 * t62 + t30 * t61;
t45 = t6 * pkin(3) + t47;
t15 = t55 * t35 + t36 * t52;
t44 = t15 * pkin(3) + t46;
t29 = sin(pkin(11));
t32 = cos(pkin(11));
t43 = t29 * rSges(6,1) + t32 * rSges(6,2) + qJ(4);
t42 = t32 * rSges(6,1) - t29 * rSges(6,2) + t66;
t28 = pkin(11) + qJ(6);
t23 = sin(t28);
t24 = cos(t28);
t41 = t24 * rSges(7,1) - t23 * rSges(7,2) + t32 * pkin(5) + t66;
t11 = t30 * t37 + t33 * t51;
t40 = t11 * pkin(2) - t33 * t65 + t54;
t4 = t11 * t62 - t33 * t61;
t39 = t4 * pkin(3) + t40;
t38 = t23 * rSges(7,1) + t24 * rSges(7,2) + t29 * pkin(5) + qJ(4);
t14 = t35 * t60 - t55 * t62;
t12 = t30 * t50 + t33 * t36;
t10 = t30 * t36 - t33 * t50;
t5 = t13 * t35 - t30 * t52;
t3 = t11 * t35 + t33 * t52;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t33 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t33 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t13 * rSges(3,1) - t12 * rSges(3,2) + t49) + g(2) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t54) + g(3) * (t55 * rSges(3,3) + t48) + (g(1) * rSges(3,3) * t30 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t37) + g(2) * (-rSges(3,3) - pkin(7)) * t33) * t31) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t63 * t12 + t47) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t63 * t10 + t40) + g(3) * (t15 * rSges(4,1) - t14 * rSges(4,2) - t63 * t59 + t46)) - m(5) * (g(1) * (-t6 * rSges(5,2) + t64 * t12 + t57 * t5 + t45) + g(2) * (-t4 * rSges(5,2) + t64 * t10 + t57 * t3 + t39) + g(3) * (-t15 * rSges(5,2) + t57 * t14 - t64 * t59 + t44)) - m(6) * (g(1) * (t42 * t12 + t43 * t5 + t56 * t6 + t45) + g(2) * (t42 * t10 + t43 * t3 + t56 * t4 + t39) + g(3) * (t43 * t14 + t56 * t15 - t42 * t59 + t44)) - m(7) * (g(1) * (t41 * t12 + t38 * t5 + t58 * t6 + t45) + g(2) * (t41 * t10 + t38 * t3 + t58 * t4 + t39) + g(3) * (t38 * t14 + t58 * t15 - t41 * t59 + t44));
U  = t1;
