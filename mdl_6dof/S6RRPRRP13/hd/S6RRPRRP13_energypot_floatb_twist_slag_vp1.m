% Calculate potential energy for
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP13_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:28
% EndTime: 2019-03-09 12:55:28
% DurationCPUTime: 0.56s
% Computational Cost: add. (281->128), mult. (536->154), div. (0->0), fcn. (622->10), ass. (0->55)
t68 = rSges(7,3) + qJ(6) + pkin(10);
t37 = sin(qJ(1));
t67 = g(1) * t37;
t41 = cos(qJ(1));
t66 = g(2) * t41;
t65 = rSges(6,3) + pkin(10);
t32 = cos(pkin(6));
t40 = cos(qJ(2));
t56 = t37 * t40;
t36 = sin(qJ(2));
t57 = t36 * t41;
t17 = t32 * t57 + t56;
t34 = sin(qJ(5));
t64 = t17 * t34;
t31 = sin(pkin(6));
t63 = t31 * t36;
t62 = t31 * t37;
t61 = t31 * t40;
t60 = t31 * t41;
t59 = t34 * t36;
t58 = t36 * t37;
t55 = t40 * t41;
t54 = qJ(3) * t40;
t53 = pkin(7) + r_base(3);
t52 = t37 * pkin(1) + r_base(2);
t51 = (-pkin(3) - pkin(8)) * t41;
t50 = t32 * pkin(8) + t53;
t49 = t41 * pkin(1) + pkin(8) * t62 + r_base(1);
t48 = g(2) * t51;
t47 = pkin(2) * t63 + t50;
t46 = t32 * pkin(3) + pkin(9) * t63 + t47;
t16 = -t32 * t55 + t58;
t45 = t17 * pkin(2) + t16 * qJ(3) + t52;
t18 = t32 * t56 + t57;
t19 = -t32 * t58 + t55;
t44 = t19 * pkin(2) + qJ(3) * t18 + t49;
t43 = pkin(3) * t62 + t44;
t42 = t17 * pkin(9) + t45;
t39 = cos(qJ(4));
t38 = cos(qJ(5));
t35 = sin(qJ(4));
t26 = pkin(5) * t38 + pkin(4);
t15 = t32 * t39 - t35 * t61;
t14 = t32 * t35 + t39 * t61;
t10 = t16 * t35 - t39 * t60;
t9 = t16 * t39 + t35 * t60;
t8 = t18 * t35 + t39 * t62;
t7 = -t18 * t39 + t35 * t62;
t6 = t15 * t38 + t31 * t59;
t5 = -t15 * t34 + t38 * t63;
t4 = t10 * t38 + t64;
t3 = -t10 * t34 + t17 * t38;
t2 = t19 * t34 + t38 * t8;
t1 = t19 * t38 - t34 * t8;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t41 - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + rSges(2,2) * t41 + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (rSges(3,1) * t19 - rSges(3,2) * t18 + t49) + g(2) * (t17 * rSges(3,1) - t16 * rSges(3,2) + t52) + g(3) * (rSges(3,3) * t32 + t50) + (rSges(3,3) * t67 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t40) + (-rSges(3,3) - pkin(8)) * t66) * t31) - m(4) * (g(1) * (-rSges(4,2) * t19 + rSges(4,3) * t18 + t44) + g(2) * (-t17 * rSges(4,2) + t16 * rSges(4,3) + t45) + g(3) * (rSges(4,1) * t32 + t47) + (rSges(4,1) * t67 + g(3) * (-rSges(4,2) * t36 - rSges(4,3) * t40 - t54) + (-rSges(4,1) - pkin(8)) * t66) * t31) - m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t7 + (rSges(5,3) + pkin(9)) * t19 + t43) + g(2) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t17 * rSges(5,3) + t42) + g(3) * (rSges(5,1) * t15 - rSges(5,2) * t14 + t46) + (g(3) * (rSges(5,3) * t36 - t54) + t48) * t31) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1 + pkin(4) * t8 + pkin(9) * t19 + t65 * t7 + t43) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * pkin(4) + t31 * t51 - t65 * t9 + t42) + g(3) * (rSges(6,1) * t6 + rSges(6,2) * t5 + pkin(4) * t15 + t65 * t14 - t31 * t54 + t46)) - m(7) * (g(1) * (rSges(7,1) * t2 + rSges(7,2) * t1 + t26 * t8 + t68 * t7 + (pkin(5) * t34 + pkin(9)) * t19 + t43) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t64 + t10 * t26 - t68 * t9 + t42) + g(3) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t68 * t14 + t15 * t26 + t46) + (g(3) * (pkin(5) * t59 - t54) + t48) * t31);
U  = t11;
