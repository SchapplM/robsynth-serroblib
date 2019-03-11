% Calculate potential energy for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:13
% EndTime: 2019-03-09 09:28:14
% DurationCPUTime: 0.58s
% Computational Cost: add. (249->124), mult. (459->142), div. (0->0), fcn. (516->10), ass. (0->55)
t71 = pkin(10) + rSges(7,3);
t70 = -pkin(3) - pkin(8);
t33 = sin(qJ(1));
t69 = g(1) * t33;
t37 = cos(qJ(1));
t68 = g(2) * t37;
t36 = cos(qJ(2));
t67 = g(3) * t36;
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t66 = t28 * t32;
t65 = t28 * t33;
t35 = cos(qJ(5));
t64 = t28 * t35;
t63 = t28 * t37;
t62 = t32 * t33;
t61 = t32 * t37;
t60 = t33 * t36;
t59 = t36 * t37;
t58 = qJ(3) * t36;
t29 = cos(pkin(6));
t10 = -t29 * t59 + t62;
t57 = t10 * qJ(3);
t12 = t29 * t60 + t61;
t56 = t12 * qJ(3);
t55 = rSges(6,3) - qJ(3);
t54 = pkin(7) + r_base(3);
t53 = t33 * pkin(1) + r_base(2);
t52 = -pkin(9) - t55;
t51 = t29 * pkin(8) + t54;
t11 = t29 * t61 + t60;
t50 = t11 * pkin(2) + t53;
t49 = t37 * pkin(1) + pkin(8) * t65 + r_base(1);
t48 = pkin(2) * t66 + t51;
t13 = -t29 * t62 + t59;
t47 = t13 * pkin(2) + t49;
t46 = (-pkin(4) + t70) * t68;
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t45 = t34 * rSges(7,1) - t30 * rSges(7,2) + pkin(5);
t44 = rSges(7,1) * t30 + rSges(7,2) * t34 - qJ(3);
t43 = t11 * qJ(4) + t50;
t42 = t29 * pkin(3) + qJ(4) * t66 + t48;
t41 = pkin(3) * t65 + t13 * qJ(4) + t47;
t40 = t28 * t36 * pkin(9) + t29 * pkin(4) + t42;
t39 = t43 + t57;
t38 = pkin(4) * t65 + t41;
t31 = sin(qJ(5));
t9 = t29 * t35 + t31 * t66;
t8 = t29 * t31 - t32 * t64;
t4 = t11 * t31 - t35 * t63;
t3 = t11 * t35 + t31 * t63;
t2 = t13 * t31 + t33 * t64;
t1 = -t13 * t35 + t31 * t65;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t37 - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + rSges(2,2) * t37 + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (rSges(3,1) * t13 - rSges(3,2) * t12 + t49) + g(2) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t53) + g(3) * (rSges(3,3) * t29 + t51) + (rSges(3,3) * t69 + g(3) * (rSges(3,1) * t32 + rSges(3,2) * t36) + (-rSges(3,3) - pkin(8)) * t68) * t28) - m(4) * (g(1) * (-rSges(4,2) * t13 + rSges(4,3) * t12 + t47 + t56) + g(2) * (-t11 * rSges(4,2) + t10 * rSges(4,3) + t50 + t57) + g(3) * (rSges(4,1) * t29 + t48) + (rSges(4,1) * t69 + g(3) * (-rSges(4,2) * t32 - rSges(4,3) * t36 - t58) + (-rSges(4,1) - pkin(8)) * t68) * t28) - m(5) * (g(1) * (rSges(5,2) * t12 + rSges(5,3) * t13 + t41 + t56) + g(2) * (t10 * rSges(5,2) + t11 * rSges(5,3) + t39) + g(3) * (rSges(5,1) * t29 + t42) + (rSges(5,1) * t69 + g(3) * (-rSges(5,2) * t36 + rSges(5,3) * t32 - t58) + (-rSges(5,1) + t70) * t68) * t28) - m(6) * (g(3) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t40) + (t55 * t67 + t46) * t28 + (t4 * rSges(6,1) + t3 * rSges(6,2) + t52 * t10 + t43) * g(2) + (rSges(6,1) * t2 - rSges(6,2) * t1 + t52 * t12 + t38) * g(1)) - m(7) * (g(1) * (t45 * t2 + (-pkin(9) - t44) * t12 + t71 * t1 + t38) + g(2) * (t4 * pkin(5) - t10 * pkin(9) + (-t10 * t30 + t34 * t4) * rSges(7,1) + (-t10 * t34 - t30 * t4) * rSges(7,2) + t39 - t71 * t3) + (t44 * t67 + t46) * t28 + (t45 * t9 + t71 * t8 + t40) * g(3));
U  = t5;
