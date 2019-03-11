% Calculate potential energy for
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:50
% EndTime: 2019-03-09 14:44:51
% DurationCPUTime: 0.71s
% Computational Cost: add. (293->131), mult. (536->157), div. (0->0), fcn. (622->12), ass. (0->53)
t66 = -pkin(11) - pkin(10) - rSges(7,3);
t33 = sin(qJ(1));
t65 = g(1) * t33;
t37 = cos(qJ(1));
t64 = g(2) * t37;
t63 = pkin(10) + rSges(6,3);
t29 = cos(pkin(6));
t36 = cos(qJ(2));
t55 = t33 * t36;
t32 = sin(qJ(2));
t56 = t32 * t37;
t11 = t29 * t56 + t55;
t30 = sin(qJ(5));
t62 = t11 * t30;
t28 = sin(pkin(6));
t61 = t28 * t32;
t60 = t28 * t33;
t59 = t28 * t36;
t58 = t28 * t37;
t57 = t32 * t33;
t54 = t36 * t37;
t53 = qJ(3) * t36;
t52 = pkin(7) + r_base(3);
t51 = pkin(1) * t33 + r_base(2);
t50 = (-pkin(3) - pkin(8)) * t37;
t49 = pkin(8) * t29 + t52;
t48 = pkin(1) * t37 + pkin(8) * t60 + r_base(1);
t47 = g(2) * t50;
t46 = pkin(2) * t61 + t49;
t34 = cos(qJ(5));
t20 = pkin(5) * t34 + pkin(4);
t27 = qJ(5) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t45 = rSges(7,1) * t22 - rSges(7,2) * t21 + t20;
t10 = -t29 * t54 + t57;
t44 = pkin(2) * t11 + t10 * qJ(3) + t51;
t43 = pkin(3) * t29 + pkin(9) * t61 + t46;
t12 = t29 * t55 + t56;
t13 = -t29 * t57 + t54;
t42 = pkin(2) * t13 + t12 * qJ(3) + t48;
t41 = rSges(7,1) * t21 + rSges(7,2) * t22 + pkin(5) * t30;
t40 = pkin(3) * t60 + t42;
t39 = t11 * pkin(9) + t44;
t35 = cos(qJ(4));
t31 = sin(qJ(4));
t9 = t29 * t35 - t31 * t59;
t8 = t29 * t31 + t35 * t59;
t4 = t10 * t31 - t35 * t58;
t3 = t10 * t35 + t31 * t58;
t2 = t12 * t31 + t35 * t60;
t1 = -t12 * t35 + t31 * t60;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t37 - rSges(2,2) * t33 + r_base(1)) + g(2) * (rSges(2,1) * t33 + rSges(2,2) * t37 + r_base(2)) + g(3) * (rSges(2,3) + t52)) - m(3) * (g(1) * (rSges(3,1) * t13 - rSges(3,2) * t12 + t48) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t10 + t51) + g(3) * (rSges(3,3) * t29 + t49) + (rSges(3,3) * t65 + g(3) * (rSges(3,1) * t32 + rSges(3,2) * t36) + (-rSges(3,3) - pkin(8)) * t64) * t28) - m(4) * (g(1) * (-rSges(4,2) * t13 + rSges(4,3) * t12 + t42) + g(2) * (-rSges(4,2) * t11 + rSges(4,3) * t10 + t44) + g(3) * (rSges(4,1) * t29 + t46) + (rSges(4,1) * t65 + g(3) * (-rSges(4,2) * t32 - rSges(4,3) * t36 - t53) + (-rSges(4,1) - pkin(8)) * t64) * t28) - m(5) * (g(1) * (rSges(5,1) * t2 - rSges(5,2) * t1 + (rSges(5,3) + pkin(9)) * t13 + t40) + g(2) * (rSges(5,1) * t4 + rSges(5,2) * t3 + rSges(5,3) * t11 + t39) + g(3) * (rSges(5,1) * t9 - t8 * rSges(5,2) + t43) + (g(3) * (rSges(5,3) * t32 - t53) + t47) * t28) - m(6) * (g(1) * (t2 * pkin(4) + t13 * pkin(9) + (t13 * t30 + t2 * t34) * rSges(6,1) + (t13 * t34 - t2 * t30) * rSges(6,2) + t63 * t1 + t40) + g(2) * (t4 * pkin(4) + (t34 * t4 + t62) * rSges(6,1) + (t11 * t34 - t30 * t4) * rSges(6,2) - t63 * t3 + t28 * t50 + t39) + g(3) * (t9 * pkin(4) - t28 * t53 + (t30 * t61 + t34 * t9) * rSges(6,1) + (-t30 * t9 + t34 * t61) * rSges(6,2) + t63 * t8 + t43)) - m(7) * (g(1) * (t45 * t2 + (pkin(9) + t41) * t13 - t66 * t1 + t40) + g(2) * (t4 * t20 + pkin(5) * t62 + (t11 * t21 + t22 * t4) * rSges(7,1) + (t11 * t22 - t21 * t4) * rSges(7,2) + t39 + t66 * t3) + t47 * t28 + (t43 + t45 * t9 + (t32 * t41 - t53) * t28 - t66 * t8) * g(3));
U  = t5;
