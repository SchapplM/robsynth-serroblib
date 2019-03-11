% Calculate potential energy for
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:38:56
% EndTime: 2019-03-08 19:38:56
% DurationCPUTime: 0.48s
% Computational Cost: add. (376->138), mult. (545->172), div. (0->0), fcn. (633->14), ass. (0->51)
t35 = sin(pkin(11));
t66 = pkin(3) * t35;
t36 = sin(pkin(10));
t40 = cos(pkin(10));
t44 = sin(qJ(2));
t41 = cos(pkin(6));
t45 = cos(qJ(2));
t58 = t41 * t45;
t11 = t36 * t44 - t40 * t58;
t34 = sin(pkin(12));
t65 = t11 * t34;
t13 = t36 * t58 + t40 * t44;
t64 = t13 * t34;
t37 = sin(pkin(6));
t63 = t36 * t37;
t62 = t37 * t40;
t61 = t37 * t44;
t60 = t37 * t45;
t59 = t41 * t44;
t57 = pkin(9) + qJ(5) + rSges(7,3);
t56 = qJ(3) + rSges(4,3);
t55 = qJ(5) + rSges(6,3);
t54 = t36 * pkin(1) + r_base(2);
t53 = t35 * t63;
t52 = t34 * t60;
t51 = qJ(1) + r_base(3);
t50 = t40 * pkin(1) + pkin(7) * t63 + r_base(1);
t49 = t41 * pkin(7) + t51;
t14 = -t36 * t59 + t40 * t45;
t39 = cos(pkin(11));
t24 = t39 * pkin(3) + pkin(2);
t43 = -pkin(8) - qJ(3);
t48 = pkin(3) * t53 - t13 * t43 + t14 * t24 + t50;
t47 = t24 * t61 + t41 * t66 + t43 * t60 + t49;
t12 = t36 * t45 + t40 * t59;
t46 = t12 * t24 + (-pkin(7) - t66) * t62 - t11 * t43 + t54;
t38 = cos(pkin(12));
t33 = pkin(11) + qJ(4);
t32 = pkin(12) + qJ(6);
t28 = cos(t33);
t27 = cos(t32);
t26 = sin(t33);
t25 = sin(t32);
t23 = t38 * pkin(5) + pkin(4);
t8 = t41 * t26 + t28 * t61;
t7 = t26 * t61 - t41 * t28;
t4 = t14 * t28 + t26 * t63;
t3 = t14 * t26 - t28 * t63;
t2 = t12 * t28 - t26 * t62;
t1 = t12 * t26 + t28 * t62;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t36 * rSges(2,2) + r_base(1)) + g(2) * (t36 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t50) + g(2) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t54) + g(3) * (t41 * rSges(3,3) + t49) + (g(1) * rSges(3,3) * t36 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t45) + g(2) * (-rSges(3,3) - pkin(7)) * t40) * t37) - m(4) * (g(1) * (t14 * pkin(2) + (t14 * t39 + t53) * rSges(4,1) + (-t14 * t35 + t39 * t63) * rSges(4,2) + t56 * t13 + t50) + g(2) * (t12 * pkin(2) - pkin(7) * t62 + (t12 * t39 - t35 * t62) * rSges(4,1) + (-t12 * t35 - t39 * t62) * rSges(4,2) + t56 * t11 + t54) + g(3) * ((t35 * rSges(4,1) + t39 * rSges(4,2)) * t41 + (-t56 * t45 + (t39 * rSges(4,1) - t35 * rSges(4,2) + pkin(2)) * t44) * t37 + t49)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t13 * rSges(5,3) + t48) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t11 * rSges(5,3) + t46) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) - rSges(5,3) * t60 + t47)) - m(6) * (g(1) * (t4 * pkin(4) + (t4 * t38 + t64) * rSges(6,1) + (t13 * t38 - t4 * t34) * rSges(6,2) + t55 * t3 + t48) + g(2) * (t2 * pkin(4) + (t2 * t38 + t65) * rSges(6,1) + (t11 * t38 - t2 * t34) * rSges(6,2) + t55 * t1 + t46) + g(3) * (t8 * pkin(4) + (t8 * t38 - t52) * rSges(6,1) + (-t8 * t34 - t38 * t60) * rSges(6,2) + t55 * t7 + t47)) - m(7) * (g(1) * (t4 * t23 + pkin(5) * t64 + (t13 * t25 + t4 * t27) * rSges(7,1) + (t13 * t27 - t4 * t25) * rSges(7,2) + t57 * t3 + t48) + g(2) * (t2 * t23 + pkin(5) * t65 + (t11 * t25 + t2 * t27) * rSges(7,1) + (t11 * t27 - t2 * t25) * rSges(7,2) + t57 * t1 + t46) + g(3) * (t8 * t23 - pkin(5) * t52 + (-t25 * t60 + t8 * t27) * rSges(7,1) + (-t8 * t25 - t27 * t60) * rSges(7,2) + t57 * t7 + t47));
U  = t5;
