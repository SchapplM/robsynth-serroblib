% Calculate potential energy for
% S6PRPRRR5
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:47
% EndTime: 2019-03-08 20:40:48
% DurationCPUTime: 0.67s
% Computational Cost: add. (305->129), mult. (482->152), div. (0->0), fcn. (547->12), ass. (0->52)
t67 = pkin(8) + rSges(5,3);
t66 = pkin(10) + rSges(7,3);
t36 = sin(qJ(4));
t65 = pkin(4) * t36;
t31 = sin(pkin(11));
t64 = g(1) * t31;
t33 = cos(pkin(11));
t63 = g(2) * t33;
t32 = sin(pkin(6));
t62 = t31 * t32;
t61 = t32 * t33;
t40 = cos(qJ(2));
t60 = t32 * t40;
t34 = cos(pkin(6));
t37 = sin(qJ(2));
t59 = t34 * t37;
t58 = t34 * t40;
t57 = t31 * pkin(1) + r_base(2);
t56 = qJ(1) + r_base(3);
t15 = t31 * t40 + t33 * t59;
t55 = t15 * pkin(2) + t57;
t54 = t33 * pkin(1) + pkin(7) * t62 + r_base(1);
t53 = t34 * pkin(7) + t56;
t39 = cos(qJ(4));
t24 = pkin(4) * t39 + pkin(3);
t52 = (-pkin(7) - t24) * t63;
t17 = -t31 * t59 + t33 * t40;
t51 = t17 * pkin(2) + t54;
t50 = t32 * t37 * pkin(2) + t53;
t49 = (-qJ(3) - t65) * t40;
t48 = t34 * t24 + t50;
t47 = rSges(5,1) * t39 - rSges(5,2) * t36 + pkin(3);
t46 = rSges(5,1) * t36 + rSges(5,2) * t39 + qJ(3);
t14 = t31 * t37 - t33 * t58;
t45 = t14 * qJ(3) + t55;
t16 = t31 * t58 + t33 * t37;
t44 = t16 * qJ(3) + t51;
t41 = -pkin(9) - pkin(8);
t43 = t14 * t65 - t15 * t41 + t45;
t42 = t16 * t65 - t17 * t41 + t24 * t62 + t44;
t38 = cos(qJ(6));
t35 = sin(qJ(6));
t30 = qJ(4) + qJ(5);
t26 = cos(t30);
t25 = sin(t30);
t9 = -t25 * t60 + t26 * t34;
t8 = t25 * t34 + t26 * t60;
t4 = t14 * t25 - t26 * t61;
t3 = t14 * t26 + t25 * t61;
t2 = t16 * t25 + t26 * t62;
t1 = -t16 * t26 + t25 * t62;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t33 - rSges(2,2) * t31 + r_base(1)) + g(2) * (rSges(2,1) * t31 + rSges(2,2) * t33 + r_base(2)) + g(3) * (rSges(2,3) + t56)) - m(3) * (g(1) * (rSges(3,1) * t17 - rSges(3,2) * t16 + t54) + g(2) * (rSges(3,1) * t15 - rSges(3,2) * t14 + t57) + g(3) * (rSges(3,3) * t34 + t53) + (rSges(3,3) * t64 + g(3) * (rSges(3,1) * t37 + rSges(3,2) * t40) + (-rSges(3,3) - pkin(7)) * t63) * t32) - m(4) * (g(1) * (-rSges(4,2) * t17 + rSges(4,3) * t16 + t44) + g(2) * (-rSges(4,2) * t15 + rSges(4,3) * t14 + t45) + g(3) * (rSges(4,1) * t34 + t50) + (rSges(4,1) * t64 + g(3) * (-rSges(4,2) * t37 + (-rSges(4,3) - qJ(3)) * t40) + (-rSges(4,1) - pkin(7)) * t63) * t32) - m(5) * ((t47 * t64 + (-pkin(7) - t47) * t63) * t32 + (t50 + t47 * t34 + (t67 * t37 - t46 * t40) * t32) * g(3) + (t46 * t14 + t67 * t15 + t55) * g(2) + (t46 * t16 + t67 * t17 + t51) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t17 + t42) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t3 + rSges(6,3) * t15 + t43) + g(3) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t48) + (g(3) * (t49 + (rSges(6,3) - t41) * t37) + t52) * t32) - m(7) * (g(1) * (t2 * pkin(5) + (t17 * t35 + t2 * t38) * rSges(7,1) + (t17 * t38 - t2 * t35) * rSges(7,2) + t66 * t1 + t42) + g(2) * (t4 * pkin(5) + (t15 * t35 + t38 * t4) * rSges(7,1) + (t15 * t38 - t35 * t4) * rSges(7,2) + t43 - t66 * t3) + t52 * t32 + (t48 + (t38 * rSges(7,1) - t35 * rSges(7,2) + pkin(5)) * t9 + (t49 + (t35 * rSges(7,1) + t38 * rSges(7,2) - t41) * t37) * t32 + t66 * t8) * g(3));
U  = t5;
