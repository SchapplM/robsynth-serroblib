% Calculate potential energy for
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:17
% EndTime: 2019-03-09 11:04:17
% DurationCPUTime: 0.42s
% Computational Cost: add. (353->130), mult. (523->157), div. (0->0), fcn. (603->12), ass. (0->49)
t34 = sin(pkin(11));
t67 = pkin(3) * t34;
t66 = pkin(10) + rSges(7,3);
t35 = sin(pkin(6));
t40 = sin(qJ(2));
t65 = t35 * t40;
t41 = sin(qJ(1));
t64 = t35 * t41;
t43 = cos(qJ(2));
t63 = t35 * t43;
t44 = cos(qJ(1));
t62 = t35 * t44;
t61 = t41 * t40;
t60 = t41 * t43;
t59 = t44 * t40;
t58 = t44 * t43;
t57 = rSges(6,3) + qJ(5);
t56 = qJ(3) + rSges(4,3);
t55 = pkin(7) + r_base(3);
t54 = t41 * pkin(1) + r_base(2);
t53 = t34 * t64;
t37 = cos(pkin(6));
t52 = t37 * pkin(8) + t55;
t51 = t44 * pkin(1) + pkin(8) * t64 + r_base(1);
t36 = cos(pkin(11));
t27 = t36 * pkin(3) + pkin(2);
t38 = -pkin(9) - qJ(3);
t50 = t27 * t65 + t37 * t67 + t38 * t63 + t52;
t16 = t37 * t60 + t59;
t17 = -t37 * t61 + t58;
t49 = pkin(3) * t53 - t16 * t38 + t17 * t27 + t51;
t33 = pkin(11) + qJ(4);
t28 = sin(t33);
t29 = cos(t33);
t11 = t37 * t28 + t29 * t65;
t48 = t11 * pkin(4) + t50;
t6 = t17 * t29 + t28 * t64;
t47 = t6 * pkin(4) + t49;
t14 = -t37 * t58 + t61;
t15 = t37 * t59 + t60;
t46 = t15 * t27 + (-pkin(8) - t67) * t62 - t14 * t38 + t54;
t4 = t15 * t29 - t28 * t62;
t45 = t4 * pkin(4) + t46;
t42 = cos(qJ(6));
t39 = sin(qJ(6));
t10 = t28 * t65 - t37 * t29;
t5 = t17 * t28 - t29 * t64;
t3 = t15 * t28 + t29 * t62;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t44 * rSges(2,1) - t41 * rSges(2,2) + r_base(1)) + g(2) * (t41 * rSges(2,1) + t44 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t55)) - m(3) * (g(1) * (t17 * rSges(3,1) - t16 * rSges(3,2) + t51) + g(2) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t54) + g(3) * (t37 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t41 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t43) + g(2) * (-rSges(3,3) - pkin(8)) * t44) * t35) - m(4) * (g(1) * (t17 * pkin(2) + (t17 * t36 + t53) * rSges(4,1) + (-t17 * t34 + t36 * t64) * rSges(4,2) + t56 * t16 + t51) + g(2) * (t15 * pkin(2) - pkin(8) * t62 + (t15 * t36 - t34 * t62) * rSges(4,1) + (-t15 * t34 - t36 * t62) * rSges(4,2) + t56 * t14 + t54) + g(3) * ((t34 * rSges(4,1) + t36 * rSges(4,2)) * t37 + (-t56 * t43 + (t36 * rSges(4,1) - t34 * rSges(4,2) + pkin(2)) * t40) * t35 + t52)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t16 * rSges(5,3) + t49) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t14 * rSges(5,3) + t46) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) - rSges(5,3) * t63 + t50)) - m(6) * (g(1) * (t16 * rSges(6,1) - t6 * rSges(6,2) + t57 * t5 + t47) + g(2) * (t14 * rSges(6,1) - t4 * rSges(6,2) + t57 * t3 + t45) + g(3) * (-rSges(6,1) * t63 - t11 * rSges(6,2) + t57 * t10 + t48)) - m(7) * (g(1) * (t16 * pkin(5) + t5 * qJ(5) + (t16 * t42 + t5 * t39) * rSges(7,1) + (-t16 * t39 + t5 * t42) * rSges(7,2) + t66 * t6 + t47) + g(2) * (t14 * pkin(5) + t3 * qJ(5) + (t14 * t42 + t3 * t39) * rSges(7,1) + (-t14 * t39 + t3 * t42) * rSges(7,2) + t66 * t4 + t45) + g(3) * (-pkin(5) * t63 + t10 * qJ(5) + (t10 * t39 - t42 * t63) * rSges(7,1) + (t10 * t42 + t39 * t63) * rSges(7,2) + t66 * t11 + t48));
U  = t1;
