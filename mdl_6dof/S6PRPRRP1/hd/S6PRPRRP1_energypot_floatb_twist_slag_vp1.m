% Calculate potential energy for
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:07
% EndTime: 2019-03-08 19:55:07
% DurationCPUTime: 0.55s
% Computational Cost: add. (409->130), mult. (853->163), div. (0->0), fcn. (1059->12), ass. (0->54)
t38 = sin(pkin(11));
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t61 = cos(pkin(11));
t27 = -t46 * t38 + t49 * t61;
t69 = rSges(5,3) + pkin(8);
t68 = rSges(6,3) + pkin(9);
t39 = sin(pkin(10));
t40 = sin(pkin(6));
t67 = t39 * t40;
t41 = cos(pkin(10));
t66 = t41 * t40;
t42 = cos(pkin(6));
t65 = t42 * t46;
t64 = t42 * t49;
t62 = rSges(7,3) + qJ(6) + pkin(9);
t35 = t49 * pkin(2) + pkin(1);
t60 = t41 * t35 + r_base(1);
t59 = qJ(1) + r_base(3);
t44 = sin(qJ(5));
t58 = pkin(5) * t44 + pkin(8);
t25 = pkin(2) * t65 + (-pkin(7) - qJ(3)) * t40;
t56 = t41 * t25 + t39 * t35 + r_base(2);
t55 = t42 * pkin(7) + t59;
t26 = -t49 * t38 - t46 * t61;
t24 = t26 * t42;
t14 = -t41 * t24 + t39 * t27;
t54 = t14 * pkin(3) + t56;
t53 = t40 * t46 * pkin(2) + t42 * qJ(3) + t55;
t16 = t39 * t24 + t41 * t27;
t52 = t16 * pkin(3) - t39 * t25 + t60;
t23 = t26 * t40;
t51 = -t23 * pkin(3) + t53;
t50 = t27 * t42;
t48 = cos(qJ(4));
t47 = cos(qJ(5));
t45 = sin(qJ(4));
t34 = t47 * pkin(5) + pkin(4);
t22 = t27 * t40;
t18 = -t23 * t48 + t42 * t45;
t17 = -t23 * t45 - t42 * t48;
t15 = t41 * t26 - t39 * t50;
t13 = t39 * t26 + t41 * t50;
t10 = t16 * t48 + t45 * t67;
t9 = t16 * t45 - t48 * t67;
t8 = t14 * t48 - t45 * t66;
t7 = t14 * t45 + t48 * t66;
t6 = t18 * t47 - t22 * t44;
t5 = -t18 * t44 - t22 * t47;
t4 = t10 * t47 - t15 * t44;
t3 = -t10 * t44 - t15 * t47;
t2 = -t13 * t44 + t8 * t47;
t1 = -t13 * t47 - t8 * t44;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t41 * rSges(2,1) - t39 * rSges(2,2) + r_base(1)) + g(2) * (t39 * rSges(2,1) + t41 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t59)) - m(3) * (g(1) * (t41 * pkin(1) + r_base(1) + (-t39 * t65 + t41 * t49) * rSges(3,1) + (-t39 * t64 - t41 * t46) * rSges(3,2)) + g(2) * (t39 * pkin(1) + r_base(2) + (t39 * t49 + t41 * t65) * rSges(3,1) + (-t39 * t46 + t41 * t64) * rSges(3,2)) + g(3) * (t42 * rSges(3,3) + t55) + (g(3) * (rSges(3,1) * t46 + rSges(3,2) * t49) + (g(1) * t39 - g(2) * t41) * (rSges(3,3) + pkin(7))) * t40) - m(4) * (g(1) * (t16 * rSges(4,1) + t15 * rSges(4,2) + (rSges(4,3) * t40 - t25) * t39 + t60) + g(2) * (t14 * rSges(4,1) + t13 * rSges(4,2) - rSges(4,3) * t66 + t56) + g(3) * (-t23 * rSges(4,1) + t22 * rSges(4,2) + t42 * rSges(4,3) + t53)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) - t15 * t69 + t52) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) - t13 * t69 + t54) + g(3) * (t18 * rSges(5,1) - t17 * rSges(5,2) - t22 * t69 + t51)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * pkin(4) - t15 * pkin(8) + t68 * t9 + t52) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) - t13 * pkin(8) + t68 * t7 + t54) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t18 * pkin(4) - t22 * pkin(8) + t17 * t68 + t51)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t10 * t34 - t15 * t58 + t62 * t9 + t52) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t13 * t58 + t8 * t34 + t62 * t7 + t54) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t17 * t62 + t18 * t34 - t22 * t58 + t51));
U  = t11;
