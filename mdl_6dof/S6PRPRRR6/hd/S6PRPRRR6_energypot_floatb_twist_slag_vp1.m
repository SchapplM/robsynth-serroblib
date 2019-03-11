% Calculate potential energy for
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:04
% EndTime: 2019-03-08 20:45:05
% DurationCPUTime: 0.67s
% Computational Cost: add. (293->131), mult. (536->160), div. (0->0), fcn. (622->12), ass. (0->52)
t65 = -pkin(10) - pkin(9) - rSges(7,3);
t28 = sin(pkin(11));
t64 = g(1) * t28;
t30 = cos(pkin(11));
t63 = g(2) * t30;
t32 = sin(qJ(5));
t37 = cos(qJ(2));
t31 = cos(pkin(6));
t34 = sin(qJ(2));
t55 = t31 * t34;
t9 = t28 * t37 + t30 * t55;
t62 = t32 * t9;
t61 = pkin(9) + rSges(6,3);
t29 = sin(pkin(6));
t60 = t28 * t29;
t33 = sin(qJ(4));
t59 = t29 * t33;
t58 = t29 * t34;
t36 = cos(qJ(4));
t57 = t29 * t36;
t56 = t29 * t37;
t54 = t31 * t37;
t53 = t37 * qJ(3);
t52 = t28 * pkin(1) + r_base(2);
t51 = qJ(1) + r_base(3);
t50 = (-pkin(3) - pkin(7)) * t30;
t49 = t30 * pkin(1) + pkin(7) * t60 + r_base(1);
t48 = t31 * pkin(7) + t51;
t47 = g(2) * t50;
t46 = pkin(2) * t58 + t48;
t35 = cos(qJ(5));
t20 = pkin(5) * t35 + pkin(4);
t27 = qJ(5) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t45 = rSges(7,1) * t22 - rSges(7,2) * t21 + t20;
t8 = t28 * t34 - t30 * t54;
t44 = t9 * pkin(2) + t8 * qJ(3) + t52;
t43 = t31 * pkin(3) + pkin(8) * t58 + t46;
t10 = t28 * t54 + t30 * t34;
t11 = -t28 * t55 + t30 * t37;
t42 = t11 * pkin(2) + t10 * qJ(3) + t49;
t41 = t21 * rSges(7,1) + t22 * rSges(7,2) + t32 * pkin(5);
t40 = pkin(3) * t60 + t42;
t39 = t9 * pkin(8) + t44;
t13 = t31 * t36 - t33 * t56;
t12 = t31 * t33 + t36 * t56;
t4 = -t30 * t57 + t33 * t8;
t3 = t30 * t59 + t36 * t8;
t2 = t10 * t33 + t28 * t57;
t1 = -t10 * t36 + t28 * t59;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t30 - rSges(2,2) * t28 + r_base(1)) + g(2) * (rSges(2,1) * t28 + rSges(2,2) * t30 + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (rSges(3,1) * t11 - rSges(3,2) * t10 + t49) + g(2) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t52) + g(3) * (rSges(3,3) * t31 + t48) + (rSges(3,3) * t64 + g(3) * (rSges(3,1) * t34 + rSges(3,2) * t37) + (-rSges(3,3) - pkin(7)) * t63) * t29) - m(4) * (g(1) * (-rSges(4,2) * t11 + rSges(4,3) * t10 + t42) + g(2) * (-rSges(4,2) * t9 + rSges(4,3) * t8 + t44) + g(3) * (rSges(4,1) * t31 + t46) + (rSges(4,1) * t64 + g(3) * (-rSges(4,2) * t34 - rSges(4,3) * t37 - t53) + (-rSges(4,1) - pkin(7)) * t63) * t29) - m(5) * (g(1) * (rSges(5,1) * t2 - rSges(5,2) * t1 + (rSges(5,3) + pkin(8)) * t11 + t40) + g(2) * (rSges(5,1) * t4 + rSges(5,2) * t3 + rSges(5,3) * t9 + t39) + g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12 + t43) + (g(3) * (rSges(5,3) * t34 - t53) + t47) * t29) - m(6) * (g(1) * (t2 * pkin(4) + t11 * pkin(8) + (t11 * t32 + t2 * t35) * rSges(6,1) + (t11 * t35 - t2 * t32) * rSges(6,2) + t61 * t1 + t40) + g(2) * (t4 * pkin(4) + (t35 * t4 + t62) * rSges(6,1) + (-t32 * t4 + t35 * t9) * rSges(6,2) - t61 * t3 + t29 * t50 + t39) + g(3) * (t13 * pkin(4) - t29 * t53 + (t13 * t35 + t32 * t58) * rSges(6,1) + (-t13 * t32 + t35 * t58) * rSges(6,2) + t61 * t12 + t43)) - m(7) * (g(1) * (t45 * t2 + (pkin(8) + t41) * t11 - t65 * t1 + t40) + g(2) * (t4 * t20 + pkin(5) * t62 + (t21 * t9 + t22 * t4) * rSges(7,1) + (-t21 * t4 + t22 * t9) * rSges(7,2) + t39 + t65 * t3) + t47 * t29 + (t43 + t45 * t13 + (t34 * t41 - t53) * t29 - t65 * t12) * g(3));
U  = t5;
