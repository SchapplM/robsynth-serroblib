% Calculate potential energy for
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:05
% EndTime: 2019-03-08 19:51:05
% DurationCPUTime: 0.58s
% Computational Cost: add. (270->123), mult. (514->144), div. (0->0), fcn. (592->10), ass. (0->49)
t66 = pkin(9) + rSges(7,3);
t65 = rSges(6,3) + qJ(5);
t29 = sin(pkin(10));
t64 = g(1) * t29;
t31 = cos(pkin(10));
t63 = g(2) * t31;
t30 = sin(pkin(6));
t62 = t29 * t30;
t34 = sin(qJ(4));
t61 = t30 * t34;
t35 = sin(qJ(2));
t60 = t30 * t35;
t37 = cos(qJ(4));
t59 = t30 * t37;
t38 = cos(qJ(2));
t58 = t30 * t38;
t32 = cos(pkin(6));
t57 = t32 * t35;
t56 = t32 * t38;
t55 = t38 * qJ(3);
t54 = t29 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t52 = t31 * pkin(1) + pkin(7) * t62 + r_base(1);
t51 = t32 * pkin(7) + t53;
t50 = (-pkin(3) - pkin(7)) * t63;
t49 = pkin(2) * t60 + t51;
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t48 = rSges(7,1) * t36 - rSges(7,2) * t33 + pkin(5);
t47 = t33 * rSges(7,1) + t36 * rSges(7,2) + qJ(5);
t12 = t29 * t35 - t31 * t56;
t13 = t29 * t38 + t31 * t57;
t46 = t13 * pkin(2) + t12 * qJ(3) + t54;
t45 = t32 * pkin(3) + pkin(8) * t60 + t49;
t14 = t29 * t56 + t31 * t35;
t15 = -t29 * t57 + t31 * t38;
t44 = t15 * pkin(2) + t14 * qJ(3) + t52;
t17 = t32 * t37 - t34 * t58;
t43 = t17 * pkin(4) + t45;
t42 = pkin(3) * t62 + t44;
t41 = t13 * pkin(8) + t46;
t4 = t14 * t34 + t29 * t59;
t40 = t4 * pkin(4) + t42;
t5 = t12 * t37 + t31 * t61;
t6 = -t12 * t34 + t31 * t59;
t39 = -t6 * pkin(4) - t5 * qJ(5) + t41;
t16 = t32 * t34 + t37 * t58;
t3 = -t14 * t37 + t29 * t61;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t31 - rSges(2,2) * t29 + r_base(1)) + g(2) * (rSges(2,1) * t29 + rSges(2,2) * t31 + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (rSges(3,1) * t15 - rSges(3,2) * t14 + t52) + g(2) * (rSges(3,1) * t13 - rSges(3,2) * t12 + t54) + g(3) * (t32 * rSges(3,3) + t51) + (rSges(3,3) * t64 + g(3) * (rSges(3,1) * t35 + rSges(3,2) * t38) + (-rSges(3,3) - pkin(7)) * t63) * t30) - m(4) * (g(1) * (-rSges(4,2) * t15 + rSges(4,3) * t14 + t44) + g(2) * (-rSges(4,2) * t13 + rSges(4,3) * t12 + t46) + g(3) * (t32 * rSges(4,1) + t49) + (rSges(4,1) * t64 + g(3) * (-rSges(4,2) * t35 - rSges(4,3) * t38 - t55) + (-rSges(4,1) - pkin(7)) * t63) * t30) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + (rSges(5,3) + pkin(8)) * t15 + t42) + g(2) * (-rSges(5,1) * t6 + rSges(5,2) * t5 + rSges(5,3) * t13 + t41) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) + t45) + (g(3) * (rSges(5,3) * t35 - t55) + t50) * t30) - m(6) * (g(1) * (-rSges(6,2) * t4 + t65 * t3 + (rSges(6,1) + pkin(8)) * t15 + t40) + g(2) * (rSges(6,1) * t13 + rSges(6,2) * t6 - rSges(6,3) * t5 + t39) + g(3) * (-t17 * rSges(6,2) + t16 * t65 + t43) + (g(3) * (rSges(6,1) * t35 - t55) + t50) * t30) - m(7) * (g(1) * (t66 * t4 + t47 * t3 + (pkin(8) + t48) * t15 + t40) + g(2) * (t13 * pkin(5) + (t13 * t36 - t5 * t33) * rSges(7,1) + (-t13 * t33 - t36 * t5) * rSges(7,2) + t39 - t66 * t6) + t50 * t30 + (t43 + t47 * t16 + (t48 * t35 - t55) * t30 + t66 * t17) * g(3));
U  = t1;
