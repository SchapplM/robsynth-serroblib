% Calculate potential energy for
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:03
% EndTime: 2019-03-09 09:17:04
% DurationCPUTime: 0.53s
% Computational Cost: add. (253->127), mult. (469->153), div. (0->0), fcn. (530->10), ass. (0->47)
t33 = sin(qJ(1));
t63 = g(1) * t33;
t37 = cos(qJ(1));
t62 = g(2) * t37;
t61 = pkin(10) + rSges(7,3);
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t60 = t28 * t32;
t59 = t28 * t33;
t36 = cos(qJ(2));
t58 = t28 * t36;
t57 = t28 * t37;
t56 = t32 * t33;
t55 = t32 * t37;
t54 = t33 * t36;
t53 = t36 * t37;
t52 = qJ(3) * t36;
t51 = qJ(4) * t33;
t50 = pkin(7) + r_base(3);
t49 = t33 * pkin(1) + r_base(2);
t29 = cos(pkin(6));
t48 = t29 * pkin(8) + t50;
t47 = t37 * pkin(1) + pkin(8) * t59 + r_base(1);
t46 = pkin(2) * t60 + t48;
t13 = -t29 * t53 + t56;
t14 = t29 * t55 + t54;
t45 = t14 * pkin(2) + t13 * qJ(3) + t49;
t15 = t29 * t54 + t55;
t16 = -t29 * t56 + t53;
t44 = t16 * pkin(2) + t15 * qJ(3) + t47;
t43 = t14 * pkin(3) + qJ(4) * t57 + t45;
t42 = pkin(3) * t60 - t29 * qJ(4) + t46;
t41 = t16 * pkin(3) + t44;
t40 = pkin(9) * t60 + t42;
t39 = t13 * pkin(4) + t14 * pkin(9) + t43;
t38 = t15 * pkin(4) + t16 * pkin(9) + t41;
t35 = cos(qJ(5));
t34 = cos(qJ(6));
t31 = sin(qJ(5));
t30 = sin(qJ(6));
t12 = -t29 * t31 - t35 * t58;
t11 = -t29 * t35 + t31 * t58;
t4 = t15 * t35 - t31 * t59;
t3 = t15 * t31 + t35 * t59;
t2 = t13 * t35 + t31 * t57;
t1 = t13 * t31 - t35 * t57;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t37 - t33 * rSges(2,2) + r_base(1)) + g(2) * (t33 * rSges(2,1) + rSges(2,2) * t37 + r_base(2)) + g(3) * (rSges(2,3) + t50)) - m(3) * (g(1) * (rSges(3,1) * t16 - rSges(3,2) * t15 + t47) + g(2) * (t14 * rSges(3,1) - t13 * rSges(3,2) + t49) + g(3) * (rSges(3,3) * t29 + t48) + (rSges(3,3) * t63 + g(3) * (rSges(3,1) * t32 + rSges(3,2) * t36) + (-rSges(3,3) - pkin(8)) * t62) * t28) - m(4) * (g(1) * (rSges(4,1) * t16 + rSges(4,3) * t15 + t44) + g(2) * (t14 * rSges(4,1) + t13 * rSges(4,3) + t45) + g(3) * (rSges(4,2) * t29 + t46) + (rSges(4,2) * t63 + g(3) * (rSges(4,1) * t32 - rSges(4,3) * t36 - t52) + (-rSges(4,2) - pkin(8)) * t62) * t28) - m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16 + t41) + g(2) * (t13 * rSges(5,1) - t14 * rSges(5,2) + t43) + g(3) * (-rSges(5,3) * t29 + t42) + (g(3) * (-rSges(5,2) * t32 + (-rSges(5,1) - qJ(3)) * t36) + (rSges(5,3) - pkin(8)) * t62 + (-rSges(5,3) - qJ(4)) * t63) * t28) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t16 + t38) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t14 * rSges(6,3) + t39) + g(3) * (rSges(6,1) * t12 + rSges(6,2) * t11 + t40) + (-g(1) * t51 - pkin(8) * t62 + g(3) * (rSges(6,3) * t32 - pkin(4) * t36 - t52)) * t28) - m(7) * (g(1) * (t4 * pkin(5) - t28 * t51 + (t16 * t30 + t34 * t4) * rSges(7,1) + (t16 * t34 - t30 * t4) * rSges(7,2) + t61 * t3 + t38) + g(2) * (t2 * pkin(5) - pkin(8) * t57 + (t14 * t30 + t2 * t34) * rSges(7,1) + (t14 * t34 - t2 * t30) * rSges(7,2) + t61 * t1 + t39) + g(3) * ((t34 * rSges(7,1) - t30 * rSges(7,2) + pkin(5)) * t12 + ((-pkin(4) - qJ(3)) * t36 + (rSges(7,1) * t30 + rSges(7,2) * t34) * t32) * t28 - t61 * t11 + t40));
U  = t5;
