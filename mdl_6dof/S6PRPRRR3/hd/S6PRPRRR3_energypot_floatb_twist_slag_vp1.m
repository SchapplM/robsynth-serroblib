% Calculate potential energy for
% S6PRPRRR3
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:34
% EndTime: 2019-03-08 20:31:35
% DurationCPUTime: 0.65s
% Computational Cost: add. (378->133), mult. (471->161), div. (0->0), fcn. (534->14), ass. (0->53)
t43 = -pkin(8) - qJ(3);
t68 = -t43 + rSges(5,3);
t37 = sin(pkin(12));
t67 = pkin(3) * t37;
t38 = sin(pkin(11));
t66 = g(1) * t38;
t41 = cos(pkin(11));
t65 = g(2) * t41;
t64 = pkin(10) + rSges(7,3);
t40 = cos(pkin(12));
t27 = t40 * pkin(3) + pkin(2);
t39 = sin(pkin(6));
t63 = t38 * t39;
t62 = t39 * t41;
t45 = sin(qJ(2));
t61 = t39 * t45;
t47 = cos(qJ(2));
t60 = t39 * t47;
t42 = cos(pkin(6));
t59 = t42 * t45;
t58 = t42 * t47;
t57 = qJ(3) + rSges(4,3);
t36 = pkin(12) + qJ(4);
t56 = t38 * pkin(1) + r_base(2);
t55 = qJ(1) + r_base(3);
t54 = t41 * pkin(1) + pkin(7) * t63 + r_base(1);
t53 = t42 * pkin(7) + t55;
t28 = sin(t36);
t29 = cos(t36);
t52 = rSges(5,1) * t29 - rSges(5,2) * t28 + t27;
t15 = t38 * t58 + t41 * t45;
t16 = -t38 * t59 + t41 * t47;
t19 = pkin(4) * t29 + t27;
t20 = pkin(4) * t28 + t67;
t35 = -pkin(9) + t43;
t51 = -t15 * t35 + t16 * t19 + t20 * t63 + t54;
t50 = t19 * t61 + t42 * t20 + t35 * t60 + t53;
t49 = rSges(5,1) * t28 + rSges(5,2) * t29 + t67;
t13 = t38 * t45 - t41 * t58;
t14 = t38 * t47 + t41 * t59;
t48 = t14 * t19 + (-pkin(7) - t20) * t62 - t13 * t35 + t56;
t46 = cos(qJ(6));
t44 = sin(qJ(6));
t30 = qJ(5) + t36;
t26 = cos(t30);
t25 = sin(t30);
t8 = t25 * t42 + t26 * t61;
t7 = t25 * t61 - t42 * t26;
t4 = t16 * t26 + t25 * t63;
t3 = t16 * t25 - t26 * t63;
t2 = t14 * t26 - t25 * t62;
t1 = t14 * t25 + t26 * t62;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t41 - rSges(2,2) * t38 + r_base(1)) + g(2) * (rSges(2,1) * t38 + rSges(2,2) * t41 + r_base(2)) + g(3) * (rSges(2,3) + t55)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t54) + g(2) * (rSges(3,1) * t14 - rSges(3,2) * t13 + t56) + g(3) * (t42 * rSges(3,3) + t53) + (rSges(3,3) * t66 + g(3) * (rSges(3,1) * t45 + rSges(3,2) * t47) + (-rSges(3,3) - pkin(7)) * t65) * t39) - m(4) * (g(1) * (t16 * pkin(2) + (t16 * t40 + t37 * t63) * rSges(4,1) + (-t16 * t37 + t40 * t63) * rSges(4,2) + t57 * t15 + t54) + g(2) * (t14 * pkin(2) - pkin(7) * t62 + (t14 * t40 - t37 * t62) * rSges(4,1) + (-t14 * t37 - t40 * t62) * rSges(4,2) + t57 * t13 + t56) + g(3) * ((rSges(4,1) * t37 + rSges(4,2) * t40) * t42 + (-t57 * t47 + (rSges(4,1) * t40 - rSges(4,2) * t37 + pkin(2)) * t45) * t39 + t53)) - m(5) * ((t49 * t66 + (-pkin(7) - t49) * t65) * t39 + (t53 + t49 * t42 + (t52 * t45 - t68 * t47) * t39) * g(3) + (t68 * t13 + t52 * t14 + t56) * g(2) + (t68 * t15 + t52 * t16 + t54) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t15 + t51) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t13 + t48) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) - rSges(6,3) * t60 + t50)) - m(7) * (g(1) * (t4 * pkin(5) + (t15 * t44 + t4 * t46) * rSges(7,1) + (t15 * t46 - t4 * t44) * rSges(7,2) + t64 * t3 + t51) + g(2) * (t2 * pkin(5) + (t13 * t44 + t2 * t46) * rSges(7,1) + (t13 * t46 - t2 * t44) * rSges(7,2) + t64 * t1 + t48) + g(3) * (t8 * pkin(5) + (-t44 * t60 + t8 * t46) * rSges(7,1) + (-t8 * t44 - t46 * t60) * rSges(7,2) + t64 * t7 + t50));
U  = t5;
