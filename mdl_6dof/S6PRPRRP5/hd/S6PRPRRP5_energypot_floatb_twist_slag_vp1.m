% Calculate potential energy for
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:56
% EndTime: 2019-03-08 20:13:57
% DurationCPUTime: 0.53s
% Computational Cost: add. (281->128), mult. (536->157), div. (0->0), fcn. (622->10), ass. (0->54)
t67 = rSges(7,3) + qJ(6) + pkin(9);
t31 = sin(pkin(10));
t66 = g(1) * t31;
t33 = cos(pkin(10));
t65 = g(2) * t33;
t64 = rSges(6,3) + pkin(9);
t41 = cos(qJ(2));
t34 = cos(pkin(6));
t38 = sin(qJ(2));
t57 = t34 * t38;
t15 = t31 * t41 + t33 * t57;
t36 = sin(qJ(5));
t63 = t15 * t36;
t32 = sin(pkin(6));
t62 = t31 * t32;
t37 = sin(qJ(4));
t61 = t32 * t37;
t60 = t32 * t38;
t40 = cos(qJ(4));
t59 = t32 * t40;
t58 = t32 * t41;
t56 = t34 * t41;
t55 = t36 * t38;
t54 = qJ(3) * t41;
t53 = t31 * pkin(1) + r_base(2);
t52 = qJ(1) + r_base(3);
t51 = (-pkin(3) - pkin(7)) * t33;
t50 = t33 * pkin(1) + pkin(7) * t62 + r_base(1);
t49 = t34 * pkin(7) + t52;
t48 = g(2) * t51;
t47 = pkin(2) * t60 + t49;
t14 = t31 * t38 - t33 * t56;
t46 = t15 * pkin(2) + qJ(3) * t14 + t53;
t45 = t34 * pkin(3) + pkin(8) * t60 + t47;
t16 = t31 * t56 + t33 * t38;
t17 = -t31 * t57 + t33 * t41;
t44 = t17 * pkin(2) + qJ(3) * t16 + t50;
t43 = pkin(3) * t62 + t44;
t42 = pkin(8) * t15 + t46;
t39 = cos(qJ(5));
t26 = pkin(5) * t39 + pkin(4);
t19 = t34 * t40 - t37 * t58;
t18 = t34 * t37 + t40 * t58;
t10 = t19 * t39 + t32 * t55;
t9 = -t19 * t36 + t39 * t60;
t8 = t14 * t37 - t33 * t59;
t7 = t14 * t40 + t33 * t61;
t6 = t16 * t37 + t31 * t59;
t5 = -t16 * t40 + t31 * t61;
t4 = t39 * t8 + t63;
t3 = t15 * t39 - t36 * t8;
t2 = t17 * t36 + t39 * t6;
t1 = t17 * t39 - t36 * t6;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t33 - rSges(2,2) * t31 + r_base(1)) + g(2) * (rSges(2,1) * t31 + rSges(2,2) * t33 + r_base(2)) + g(3) * (rSges(2,3) + t52)) - m(3) * (g(1) * (rSges(3,1) * t17 - rSges(3,2) * t16 + t50) + g(2) * (rSges(3,1) * t15 - rSges(3,2) * t14 + t53) + g(3) * (t34 * rSges(3,3) + t49) + (rSges(3,3) * t66 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t41) + (-rSges(3,3) - pkin(7)) * t65) * t32) - m(4) * (g(1) * (-rSges(4,2) * t17 + rSges(4,3) * t16 + t44) + g(2) * (-rSges(4,2) * t15 + rSges(4,3) * t14 + t46) + g(3) * (t34 * rSges(4,1) + t47) + (rSges(4,1) * t66 + g(3) * (-rSges(4,2) * t38 - rSges(4,3) * t41 - t54) + (-rSges(4,1) - pkin(7)) * t65) * t32) - m(5) * (g(1) * (rSges(5,1) * t6 - rSges(5,2) * t5 + (rSges(5,3) + pkin(8)) * t17 + t43) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t7 + rSges(5,3) * t15 + t42) + g(3) * (t19 * rSges(5,1) - t18 * rSges(5,2) + t45) + (g(3) * (rSges(5,3) * t38 - t54) + t48) * t32) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1 + pkin(4) * t6 + pkin(8) * t17 + t64 * t5 + t43) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t3 + pkin(4) * t8 + t32 * t51 - t64 * t7 + t42) + g(3) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t19 * pkin(4) + t64 * t18 - t32 * t54 + t45)) - m(7) * (g(1) * (rSges(7,1) * t2 + rSges(7,2) * t1 + t26 * t6 + t67 * t5 + (pkin(5) * t36 + pkin(8)) * t17 + t43) + g(2) * (rSges(7,1) * t4 + rSges(7,2) * t3 + pkin(5) * t63 + t26 * t8 - t67 * t7 + t42) + g(3) * (t10 * rSges(7,1) + t9 * rSges(7,2) + t67 * t18 + t19 * t26 + t45) + (g(3) * (pkin(5) * t55 - t54) + t48) * t32);
U  = t11;
