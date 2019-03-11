% Calculate potential energy for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:33
% EndTime: 2019-03-08 19:46:34
% DurationCPUTime: 0.68s
% Computational Cost: add. (293->131), mult. (536->160), div. (0->0), fcn. (622->12), ass. (0->52)
t65 = -pkin(9) - qJ(5) - rSges(7,3);
t29 = sin(pkin(10));
t64 = g(1) * t29;
t32 = cos(pkin(10));
t63 = g(2) * t32;
t28 = sin(pkin(11));
t38 = cos(qJ(2));
t33 = cos(pkin(6));
t36 = sin(qJ(2));
t56 = t33 * t36;
t9 = t29 * t38 + t32 * t56;
t62 = t28 * t9;
t30 = sin(pkin(6));
t61 = t29 * t30;
t35 = sin(qJ(4));
t60 = t30 * t35;
t59 = t30 * t36;
t37 = cos(qJ(4));
t58 = t30 * t37;
t57 = t30 * t38;
t55 = t33 * t38;
t54 = t38 * qJ(3);
t53 = qJ(5) + rSges(6,3);
t52 = t29 * pkin(1) + r_base(2);
t51 = qJ(1) + r_base(3);
t50 = (-pkin(3) - pkin(7)) * t32;
t49 = t32 * pkin(1) + pkin(7) * t61 + r_base(1);
t48 = t33 * pkin(7) + t51;
t47 = g(2) * t50;
t46 = pkin(2) * t59 + t48;
t31 = cos(pkin(11));
t20 = pkin(5) * t31 + pkin(4);
t27 = pkin(11) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t45 = rSges(7,1) * t22 - rSges(7,2) * t21 + t20;
t8 = t29 * t36 - t32 * t55;
t44 = t9 * pkin(2) + t8 * qJ(3) + t52;
t43 = t33 * pkin(3) + pkin(8) * t59 + t46;
t10 = t29 * t55 + t32 * t36;
t11 = -t29 * t56 + t32 * t38;
t42 = t11 * pkin(2) + t10 * qJ(3) + t49;
t41 = t21 * rSges(7,1) + t22 * rSges(7,2) + t28 * pkin(5);
t40 = pkin(3) * t61 + t42;
t39 = t9 * pkin(8) + t44;
t13 = t33 * t37 - t35 * t57;
t12 = t33 * t35 + t37 * t57;
t4 = -t32 * t58 + t35 * t8;
t3 = t32 * t60 + t37 * t8;
t2 = t10 * t35 + t29 * t58;
t1 = -t10 * t37 + t29 * t60;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t32 - rSges(2,2) * t29 + r_base(1)) + g(2) * (rSges(2,1) * t29 + rSges(2,2) * t32 + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (rSges(3,1) * t11 - rSges(3,2) * t10 + t49) + g(2) * (rSges(3,1) * t9 - rSges(3,2) * t8 + t52) + g(3) * (t33 * rSges(3,3) + t48) + (rSges(3,3) * t64 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t38) + (-rSges(3,3) - pkin(7)) * t63) * t30) - m(4) * (g(1) * (-rSges(4,2) * t11 + rSges(4,3) * t10 + t42) + g(2) * (-rSges(4,2) * t9 + rSges(4,3) * t8 + t44) + g(3) * (t33 * rSges(4,1) + t46) + (rSges(4,1) * t64 + g(3) * (-rSges(4,2) * t36 - rSges(4,3) * t38 - t54) + (-rSges(4,1) - pkin(7)) * t63) * t30) - m(5) * (g(1) * (rSges(5,1) * t2 - rSges(5,2) * t1 + (rSges(5,3) + pkin(8)) * t11 + t40) + g(2) * (rSges(5,1) * t4 + rSges(5,2) * t3 + rSges(5,3) * t9 + t39) + g(3) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t43) + (g(3) * (rSges(5,3) * t36 - t54) + t47) * t30) - m(6) * (g(1) * (t2 * pkin(4) + t11 * pkin(8) + (t11 * t28 + t2 * t31) * rSges(6,1) + (t11 * t31 - t2 * t28) * rSges(6,2) + t53 * t1 + t40) + g(2) * (t4 * pkin(4) + (t31 * t4 + t62) * rSges(6,1) + (-t28 * t4 + t31 * t9) * rSges(6,2) + t30 * t50 - t53 * t3 + t39) + g(3) * (t13 * pkin(4) - t30 * t54 + (t13 * t31 + t28 * t59) * rSges(6,1) + (-t13 * t28 + t31 * t59) * rSges(6,2) + t53 * t12 + t43)) - m(7) * (g(1) * (t45 * t2 + (pkin(8) + t41) * t11 - t65 * t1 + t40) + g(2) * (t4 * t20 + pkin(5) * t62 + (t21 * t9 + t22 * t4) * rSges(7,1) + (-t21 * t4 + t22 * t9) * rSges(7,2) + t39 + t65 * t3) + t47 * t30 + (t43 + t45 * t13 + (t36 * t41 - t54) * t30 - t65 * t12) * g(3));
U  = t5;
