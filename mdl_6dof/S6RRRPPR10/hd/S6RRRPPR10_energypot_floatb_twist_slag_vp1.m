% Calculate potential energy for
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:22
% EndTime: 2019-03-09 16:21:22
% DurationCPUTime: 0.38s
% Computational Cost: add. (322->113), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->51)
t66 = pkin(4) + pkin(9);
t65 = rSges(5,1) + pkin(9);
t64 = rSges(4,3) + pkin(9);
t63 = cos(qJ(3));
t30 = sin(pkin(6));
t34 = sin(qJ(2));
t62 = t30 * t34;
t35 = sin(qJ(1));
t61 = t30 * t35;
t36 = cos(qJ(2));
t60 = t30 * t36;
t37 = cos(qJ(1));
t59 = t30 * t37;
t58 = pkin(10) + qJ(5) + rSges(7,3);
t57 = rSges(5,3) + qJ(4);
t56 = qJ(5) + rSges(6,3);
t55 = cos(pkin(6));
t54 = pkin(7) + r_base(3);
t53 = t35 * pkin(1) + r_base(2);
t52 = t30 * t63;
t51 = t55 * pkin(8) + t54;
t50 = t35 * t55;
t49 = t37 * t55;
t48 = t37 * pkin(1) + pkin(8) * t61 + r_base(1);
t47 = pkin(2) * t62 + t51;
t15 = -t34 * t50 + t37 * t36;
t46 = t15 * pkin(2) + t48;
t33 = sin(qJ(3));
t11 = t55 * t33 + t34 * t52;
t45 = t11 * pkin(3) + t47;
t6 = t15 * t63 + t33 * t61;
t44 = t6 * pkin(3) + t46;
t29 = sin(pkin(11));
t31 = cos(pkin(11));
t43 = t29 * rSges(6,1) + t31 * rSges(6,2) + qJ(4);
t42 = t31 * rSges(6,1) - t29 * rSges(6,2) + t66;
t28 = pkin(11) + qJ(6);
t23 = sin(t28);
t24 = cos(t28);
t41 = t24 * rSges(7,1) - t23 * rSges(7,2) + t31 * pkin(5) + t66;
t13 = t34 * t49 + t35 * t36;
t40 = t13 * pkin(2) - pkin(8) * t59 + t53;
t4 = t13 * t63 - t33 * t59;
t39 = t4 * pkin(3) + t40;
t38 = t23 * rSges(7,1) + t24 * rSges(7,2) + t29 * pkin(5) + qJ(4);
t14 = t37 * t34 + t36 * t50;
t12 = t35 * t34 - t36 * t49;
t10 = t33 * t62 - t55 * t63;
t5 = t15 * t33 - t35 * t52;
t3 = t13 * t33 + t37 * t52;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t48) + g(2) * (t13 * rSges(3,1) - t12 * rSges(3,2) + t53) + g(3) * (t55 * rSges(3,3) + t51) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t34 + rSges(3,2) * t36) + g(2) * (-rSges(3,3) - pkin(8)) * t37) * t30) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t64 * t14 + t46) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t64 * t12 + t40) + g(3) * (t11 * rSges(4,1) - t10 * rSges(4,2) - t64 * t60 + t47)) - m(5) * (g(1) * (-t6 * rSges(5,2) + t65 * t14 + t57 * t5 + t44) + g(2) * (-t4 * rSges(5,2) + t65 * t12 + t57 * t3 + t39) + g(3) * (-t11 * rSges(5,2) + t57 * t10 - t65 * t60 + t45)) - m(6) * (g(1) * (t42 * t14 + t43 * t5 + t56 * t6 + t44) + g(2) * (t42 * t12 + t43 * t3 + t56 * t4 + t39) + g(3) * (t43 * t10 + t56 * t11 - t42 * t60 + t45)) - m(7) * (g(1) * (t41 * t14 + t38 * t5 + t58 * t6 + t44) + g(2) * (t41 * t12 + t38 * t3 + t58 * t4 + t39) + g(3) * (t38 * t10 + t58 * t11 - t41 * t60 + t45));
U  = t1;
