% Calculate potential energy for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR14_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:27
% EndTime: 2019-03-09 11:32:28
% DurationCPUTime: 0.63s
% Computational Cost: add. (270->123), mult. (514->141), div. (0->0), fcn. (592->10), ass. (0->50)
t67 = pkin(10) + rSges(7,3);
t66 = rSges(6,3) + qJ(5);
t34 = sin(qJ(1));
t65 = g(1) * t34;
t38 = cos(qJ(1));
t64 = g(2) * t38;
t29 = sin(pkin(6));
t33 = sin(qJ(2));
t63 = t29 * t33;
t62 = t29 * t34;
t37 = cos(qJ(2));
t61 = t29 * t37;
t60 = t29 * t38;
t59 = t33 * t34;
t58 = t33 * t38;
t57 = t34 * t37;
t56 = t37 * t38;
t55 = qJ(3) * t37;
t54 = pkin(7) + r_base(3);
t53 = t34 * pkin(1) + r_base(2);
t30 = cos(pkin(6));
t52 = t30 * pkin(8) + t54;
t51 = t38 * pkin(1) + pkin(8) * t62 + r_base(1);
t50 = (-pkin(3) - pkin(8)) * t64;
t49 = pkin(2) * t63 + t52;
t31 = sin(qJ(6));
t35 = cos(qJ(6));
t48 = rSges(7,1) * t35 - rSges(7,2) * t31 + pkin(5);
t47 = t31 * rSges(7,1) + t35 * rSges(7,2) + qJ(5);
t46 = t30 * pkin(3) + pkin(9) * t63 + t49;
t14 = -t30 * t56 + t59;
t15 = t30 * t58 + t57;
t45 = t15 * pkin(2) + t14 * qJ(3) + t53;
t32 = sin(qJ(4));
t36 = cos(qJ(4));
t13 = t30 * t36 - t32 * t61;
t44 = t13 * pkin(4) + t46;
t16 = t30 * t57 + t58;
t17 = -t30 * t59 + t56;
t43 = t17 * pkin(2) + t16 * qJ(3) + t51;
t42 = pkin(3) * t62 + t43;
t41 = t15 * pkin(9) + t45;
t4 = t16 * t32 + t36 * t62;
t40 = t4 * pkin(4) + t42;
t5 = t14 * t36 + t32 * t60;
t6 = -t14 * t32 + t36 * t60;
t39 = -t6 * pkin(4) - t5 * qJ(5) + t41;
t12 = t30 * t32 + t36 * t61;
t3 = -t16 * t36 + t32 * t62;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t38 - t34 * rSges(2,2) + r_base(1)) + g(2) * (t34 * rSges(2,1) + rSges(2,2) * t38 + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (rSges(3,1) * t17 - rSges(3,2) * t16 + t51) + g(2) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t53) + g(3) * (rSges(3,3) * t30 + t52) + (rSges(3,3) * t65 + g(3) * (rSges(3,1) * t33 + rSges(3,2) * t37) + (-rSges(3,3) - pkin(8)) * t64) * t29) - m(4) * (g(1) * (-rSges(4,2) * t17 + rSges(4,3) * t16 + t43) + g(2) * (-t15 * rSges(4,2) + t14 * rSges(4,3) + t45) + g(3) * (rSges(4,1) * t30 + t49) + (rSges(4,1) * t65 + g(3) * (-rSges(4,2) * t33 - rSges(4,3) * t37 - t55) + (-rSges(4,1) - pkin(8)) * t64) * t29) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + (rSges(5,3) + pkin(9)) * t17 + t42) + g(2) * (-t6 * rSges(5,1) + t5 * rSges(5,2) + t15 * rSges(5,3) + t41) + g(3) * (rSges(5,1) * t13 - rSges(5,2) * t12 + t46) + (g(3) * (rSges(5,3) * t33 - t55) + t50) * t29) - m(6) * (g(1) * (-rSges(6,2) * t4 + t66 * t3 + (rSges(6,1) + pkin(9)) * t17 + t40) + g(2) * (t15 * rSges(6,1) + t6 * rSges(6,2) - t5 * rSges(6,3) + t39) + g(3) * (-rSges(6,2) * t13 + t66 * t12 + t44) + (g(3) * (rSges(6,1) * t33 - t55) + t50) * t29) - m(7) * (g(1) * (t67 * t4 + t47 * t3 + (pkin(9) + t48) * t17 + t40) + g(2) * (t15 * pkin(5) + (t15 * t35 - t31 * t5) * rSges(7,1) + (-t15 * t31 - t35 * t5) * rSges(7,2) + t39 - t67 * t6) + t50 * t29 + (t44 + t47 * t12 + (t33 * t48 - t55) * t29 + t67 * t13) * g(3));
U  = t1;
