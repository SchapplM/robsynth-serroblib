% Calculate potential energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:35:49
% EndTime: 2019-03-08 20:35:49
% DurationCPUTime: 0.48s
% Computational Cost: add. (376->138), mult. (545->172), div. (0->0), fcn. (633->14), ass. (0->51)
t34 = sin(pkin(12));
t66 = pkin(3) * t34;
t65 = pkin(9) + rSges(6,3);
t35 = sin(pkin(11));
t38 = cos(pkin(11));
t42 = sin(qJ(2));
t39 = cos(pkin(6));
t44 = cos(qJ(2));
t57 = t39 * t44;
t11 = t35 * t42 - t38 * t57;
t41 = sin(qJ(5));
t64 = t11 * t41;
t13 = t35 * t57 + t38 * t42;
t63 = t13 * t41;
t36 = sin(pkin(6));
t62 = t35 * t36;
t61 = t36 * t38;
t60 = t36 * t42;
t59 = t36 * t44;
t58 = t39 * t42;
t56 = pkin(10) + pkin(9) + rSges(7,3);
t55 = qJ(3) + rSges(4,3);
t54 = t35 * pkin(1) + r_base(2);
t53 = t34 * t62;
t52 = t41 * t59;
t51 = qJ(1) + r_base(3);
t50 = t38 * pkin(1) + pkin(7) * t62 + r_base(1);
t49 = t39 * pkin(7) + t51;
t14 = -t35 * t58 + t38 * t44;
t37 = cos(pkin(12));
t23 = t37 * pkin(3) + pkin(2);
t40 = -pkin(8) - qJ(3);
t48 = pkin(3) * t53 - t13 * t40 + t14 * t23 + t50;
t47 = t23 * t60 + t39 * t66 + t40 * t59 + t49;
t12 = t35 * t44 + t38 * t58;
t46 = t12 * t23 + (-pkin(7) - t66) * t61 - t11 * t40 + t54;
t43 = cos(qJ(5));
t33 = qJ(5) + qJ(6);
t32 = pkin(12) + qJ(4);
t28 = cos(t33);
t27 = sin(t33);
t26 = cos(t32);
t25 = sin(t32);
t24 = t43 * pkin(5) + pkin(4);
t8 = t39 * t25 + t26 * t60;
t7 = t25 * t60 - t39 * t26;
t4 = t14 * t26 + t25 * t62;
t3 = t14 * t25 - t26 * t62;
t2 = t12 * t26 - t25 * t61;
t1 = t12 * t25 + t26 * t61;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t50) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t54) + g(3) * (t39 * rSges(3,3) + t49) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t42 + rSges(3,2) * t44) + g(2) * (-rSges(3,3) - pkin(7)) * t38) * t36) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t37 + t53) * rSges(4,1) + (-t14 * t34 + t37 * t62) * rSges(4,2) + t55 * t13 + t50) + g(2) * (t12 * pkin(2) - pkin(7) * t61 + (t12 * t37 - t34 * t61) * rSges(4,1) + (-t12 * t34 - t37 * t61) * rSges(4,2) + t55 * t11 + t54) + g(3) * ((t34 * rSges(4,1) + t37 * rSges(4,2)) * t39 + (-t55 * t44 + (t37 * rSges(4,1) - t34 * rSges(4,2) + pkin(2)) * t42) * t36 + t49)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t48) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t46) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t59 + t47)) - m(6) * (g(1) * (t4 * pkin(4) + (t4 * t43 + t63) * rSges(6,1) + (t13 * t43 - t4 * t41) * rSges(6,2) + t65 * t3 + t48) + g(2) * (t2 * pkin(4) + (t2 * t43 + t64) * rSges(6,1) + (t11 * t43 - t2 * t41) * rSges(6,2) + t65 * t1 + t46) + g(3) * (t8 * pkin(4) + (t8 * t43 - t52) * rSges(6,1) + (-t8 * t41 - t43 * t59) * rSges(6,2) + t65 * t7 + t47)) - m(7) * (g(1) * (t4 * t24 + pkin(5) * t63 + (t13 * t27 + t4 * t28) * rSges(7,1) + (t13 * t28 - t4 * t27) * rSges(7,2) + t56 * t3 + t48) + g(2) * (t2 * t24 + pkin(5) * t64 + (t11 * t27 + t2 * t28) * rSges(7,1) + (t11 * t28 - t2 * t27) * rSges(7,2) + t56 * t1 + t46) + g(3) * (t8 * t24 - pkin(5) * t52 + (-t27 * t59 + t8 * t28) * rSges(7,1) + (-t8 * t27 - t28 * t59) * rSges(7,2) + t56 * t7 + t47));
U  = t5;
