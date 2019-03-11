% Calculate potential energy for
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR14_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR14_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR14_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:16
% EndTime: 2019-03-09 20:12:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (322->113), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->51)
t66 = pkin(4) + pkin(9);
t65 = rSges(5,1) + pkin(9);
t64 = rSges(4,3) + pkin(9);
t63 = pkin(10) + rSges(6,3);
t62 = cos(qJ(3));
t29 = sin(pkin(6));
t32 = sin(qJ(2));
t61 = t29 * t32;
t33 = sin(qJ(1));
t60 = t29 * t33;
t35 = cos(qJ(2));
t59 = t29 * t35;
t36 = cos(qJ(1));
t58 = t29 * t36;
t57 = pkin(11) + pkin(10) + rSges(7,3);
t56 = rSges(5,3) + qJ(4);
t55 = cos(pkin(6));
t54 = pkin(7) + r_base(3);
t53 = t33 * pkin(1) + r_base(2);
t52 = t29 * t62;
t51 = t55 * pkin(8) + t54;
t50 = t33 * t55;
t49 = t36 * t55;
t48 = t36 * pkin(1) + pkin(8) * t60 + r_base(1);
t47 = pkin(2) * t61 + t51;
t15 = -t32 * t50 + t36 * t35;
t46 = t15 * pkin(2) + t48;
t31 = sin(qJ(3));
t11 = t55 * t31 + t32 * t52;
t45 = t11 * pkin(3) + t47;
t6 = t15 * t62 + t31 * t60;
t44 = t6 * pkin(3) + t46;
t30 = sin(qJ(5));
t34 = cos(qJ(5));
t43 = t30 * rSges(6,1) + t34 * rSges(6,2) + qJ(4);
t42 = t34 * rSges(6,1) - t30 * rSges(6,2) + t66;
t28 = qJ(5) + qJ(6);
t23 = sin(t28);
t24 = cos(t28);
t41 = t24 * rSges(7,1) - t23 * rSges(7,2) + t34 * pkin(5) + t66;
t13 = t32 * t49 + t33 * t35;
t40 = t13 * pkin(2) - pkin(8) * t58 + t53;
t4 = t13 * t62 - t31 * t58;
t39 = t4 * pkin(3) + t40;
t38 = t23 * rSges(7,1) + t24 * rSges(7,2) + t30 * pkin(5) + qJ(4);
t14 = t36 * t32 + t35 * t50;
t12 = t33 * t32 - t35 * t49;
t10 = t31 * t61 - t55 * t62;
t5 = t15 * t31 - t33 * t52;
t3 = t13 * t31 + t36 * t52;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t36 * rSges(2,1) - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + t36 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t48) + g(2) * (t13 * rSges(3,1) - t12 * rSges(3,2) + t53) + g(3) * (t55 * rSges(3,3) + t51) + (g(1) * rSges(3,3) * t33 + g(3) * (rSges(3,1) * t32 + rSges(3,2) * t35) + g(2) * (-rSges(3,3) - pkin(8)) * t36) * t29) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t64 * t14 + t46) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t64 * t12 + t40) + g(3) * (t11 * rSges(4,1) - t10 * rSges(4,2) - t64 * t59 + t47)) - m(5) * (g(1) * (-t6 * rSges(5,2) + t65 * t14 + t56 * t5 + t44) + g(2) * (-t4 * rSges(5,2) + t65 * t12 + t56 * t3 + t39) + g(3) * (-t11 * rSges(5,2) + t56 * t10 - t65 * t59 + t45)) - m(6) * (g(1) * (t42 * t14 + t43 * t5 + t63 * t6 + t44) + g(2) * (t42 * t12 + t43 * t3 + t63 * t4 + t39) + g(3) * (t43 * t10 + t63 * t11 - t42 * t59 + t45)) - m(7) * (g(1) * (t41 * t14 + t38 * t5 + t57 * t6 + t44) + g(2) * (t41 * t12 + t38 * t3 + t57 * t4 + t39) + g(3) * (t38 * t10 + t57 * t11 - t41 * t59 + t45));
U  = t1;
