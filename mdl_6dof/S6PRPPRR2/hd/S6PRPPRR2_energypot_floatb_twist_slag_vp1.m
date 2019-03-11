% Calculate potential energy for
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:35
% EndTime: 2019-03-08 19:17:35
% DurationCPUTime: 0.57s
% Computational Cost: add. (357->129), mult. (724->162), div. (0->0), fcn. (881->12), ass. (0->49)
t33 = sin(pkin(11));
t36 = cos(pkin(11));
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t23 = -t41 * t33 + t44 * t36;
t67 = rSges(6,3) + pkin(8);
t66 = pkin(9) + rSges(7,3);
t34 = sin(pkin(10));
t35 = sin(pkin(6));
t65 = t34 * t35;
t37 = cos(pkin(10));
t64 = t37 * t35;
t38 = cos(pkin(6));
t63 = t38 * t41;
t62 = t38 * t44;
t59 = rSges(5,3) + qJ(4);
t29 = t44 * pkin(2) + pkin(1);
t58 = t37 * t29 + r_base(1);
t57 = qJ(1) + r_base(3);
t52 = t44 * t33 + t41 * t36;
t50 = t52 * t38;
t11 = t37 * t23 - t34 * t50;
t56 = t11 * pkin(3) + t58;
t21 = pkin(2) * t63 + (-pkin(7) - qJ(3)) * t35;
t55 = t37 * t21 + t34 * t29 + r_base(2);
t54 = t38 * pkin(7) + t57;
t9 = t34 * t23 + t37 * t50;
t53 = t9 * pkin(3) + t55;
t51 = t35 * t41 * pkin(2) + t38 * qJ(3) + t54;
t49 = t23 * t38;
t20 = t52 * t35;
t48 = t20 * pkin(3) + t51;
t10 = -t34 * t49 - t37 * t52;
t47 = pkin(4) * t65 - t10 * qJ(4) - t34 * t21 + t56;
t8 = -t34 * t52 + t37 * t49;
t46 = -pkin(4) * t64 - t8 * qJ(4) + t53;
t19 = t23 * t35;
t45 = t38 * pkin(4) - t19 * qJ(4) + t48;
t43 = cos(qJ(5));
t42 = cos(qJ(6));
t40 = sin(qJ(5));
t39 = sin(qJ(6));
t13 = -t19 * t40 + t38 * t43;
t12 = t19 * t43 + t38 * t40;
t4 = -t8 * t40 - t43 * t64;
t3 = t40 * t64 - t8 * t43;
t2 = -t10 * t40 + t43 * t65;
t1 = t10 * t43 + t40 * t65;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t37 * rSges(2,1) - t34 * rSges(2,2) + r_base(1)) + g(2) * (t34 * rSges(2,1) + t37 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t57)) - m(3) * (g(1) * (t37 * pkin(1) + r_base(1) + (-t34 * t63 + t37 * t44) * rSges(3,1) + (-t34 * t62 - t37 * t41) * rSges(3,2)) + g(2) * (t34 * pkin(1) + r_base(2) + (t34 * t44 + t37 * t63) * rSges(3,1) + (-t34 * t41 + t37 * t62) * rSges(3,2)) + g(3) * (t38 * rSges(3,3) + t54) + (g(3) * (rSges(3,1) * t41 + rSges(3,2) * t44) + (g(1) * t34 - g(2) * t37) * (rSges(3,3) + pkin(7))) * t35) - m(4) * (g(1) * (t11 * rSges(4,1) + t10 * rSges(4,2) + (rSges(4,3) * t35 - t21) * t34 + t58) + g(2) * (t9 * rSges(4,1) + t8 * rSges(4,2) - rSges(4,3) * t64 + t55) + g(3) * (t20 * rSges(4,1) + t19 * rSges(4,2) + t38 * rSges(4,3) + t51)) - m(5) * (g(1) * (-t11 * rSges(5,2) + (rSges(5,1) * t35 - t21) * t34 - t59 * t10 + t56) + g(2) * (-rSges(5,1) * t64 - t9 * rSges(5,2) - t59 * t8 + t53) + g(3) * (t38 * rSges(5,1) - t20 * rSges(5,2) - t59 * t19 + t48)) - m(6) * (g(1) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t67 * t11 + t47) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t67 * t9 + t46) + g(3) * (t13 * rSges(6,1) - t12 * rSges(6,2) + t67 * t20 + t45)) - m(7) * (g(1) * (t2 * pkin(5) + t11 * pkin(8) + (t11 * t39 + t2 * t42) * rSges(7,1) + (t11 * t42 - t2 * t39) * rSges(7,2) + t66 * t1 + t47) + g(2) * (t4 * pkin(5) + t9 * pkin(8) + (t9 * t39 + t4 * t42) * rSges(7,1) + (-t4 * t39 + t9 * t42) * rSges(7,2) - t66 * t3 + t46) + g(3) * (t13 * pkin(5) + t20 * pkin(8) + (t13 * t42 + t20 * t39) * rSges(7,1) + (-t13 * t39 + t20 * t42) * rSges(7,2) + t66 * t12 + t45));
U  = t5;
