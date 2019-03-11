% Calculate potential energy for
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:25
% EndTime: 2019-03-09 12:36:26
% DurationCPUTime: 0.38s
% Computational Cost: add. (391->126), mult. (591->151), div. (0->0), fcn. (697->12), ass. (0->57)
t41 = sin(pkin(11));
t76 = pkin(3) * t41;
t75 = rSges(7,1) + pkin(5);
t74 = rSges(7,2) + pkin(10);
t73 = rSges(6,3) + pkin(10);
t42 = sin(pkin(6));
t47 = sin(qJ(2));
t72 = t42 * t47;
t48 = sin(qJ(1));
t71 = t42 * t48;
t50 = cos(qJ(2));
t70 = t42 * t50;
t51 = cos(qJ(1));
t69 = t42 * t51;
t68 = t48 * t47;
t67 = t48 * t50;
t66 = t51 * t47;
t65 = t51 * t50;
t64 = rSges(7,3) + qJ(6);
t63 = qJ(3) + rSges(4,3);
t62 = pkin(7) + r_base(3);
t61 = t48 * pkin(1) + r_base(2);
t60 = t41 * t71;
t44 = cos(pkin(6));
t59 = t44 * pkin(8) + t62;
t58 = t51 * pkin(1) + pkin(8) * t71 + r_base(1);
t43 = cos(pkin(11));
t34 = t43 * pkin(3) + pkin(2);
t45 = -pkin(9) - qJ(3);
t57 = t34 * t72 + t44 * t76 + t45 * t70 + t59;
t24 = t44 * t67 + t66;
t25 = -t44 * t68 + t65;
t56 = pkin(3) * t60 - t24 * t45 + t25 * t34 + t58;
t40 = pkin(11) + qJ(4);
t35 = sin(t40);
t36 = cos(t40);
t17 = t44 * t35 + t36 * t72;
t55 = t17 * pkin(4) + t57;
t12 = t25 * t36 + t35 * t71;
t54 = t12 * pkin(4) + t56;
t22 = -t44 * t65 + t68;
t23 = t44 * t66 + t67;
t53 = t23 * t34 + (-pkin(8) - t76) * t69 - t22 * t45 + t61;
t10 = t23 * t36 - t35 * t69;
t52 = t10 * pkin(4) + t53;
t49 = cos(qJ(5));
t46 = sin(qJ(5));
t16 = t35 * t72 - t44 * t36;
t11 = t25 * t35 - t36 * t71;
t9 = t23 * t35 + t36 * t69;
t8 = t17 * t49 - t46 * t70;
t7 = t17 * t46 + t49 * t70;
t4 = t12 * t49 + t24 * t46;
t3 = t12 * t46 - t24 * t49;
t2 = t10 * t49 + t22 * t46;
t1 = t10 * t46 - t22 * t49;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t51 * rSges(2,1) - t48 * rSges(2,2) + r_base(1)) + g(2) * (t48 * rSges(2,1) + t51 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t62)) - m(3) * (g(1) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t58) + g(2) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t61) + g(3) * (t44 * rSges(3,3) + t59) + (g(1) * rSges(3,3) * t48 + g(3) * (rSges(3,1) * t47 + rSges(3,2) * t50) + g(2) * (-rSges(3,3) - pkin(8)) * t51) * t42) - m(4) * (g(1) * (t25 * pkin(2) + (t25 * t43 + t60) * rSges(4,1) + (-t25 * t41 + t43 * t71) * rSges(4,2) + t63 * t24 + t58) + g(2) * (t23 * pkin(2) - pkin(8) * t69 + (t23 * t43 - t41 * t69) * rSges(4,1) + (-t23 * t41 - t43 * t69) * rSges(4,2) + t63 * t22 + t61) + g(3) * ((t41 * rSges(4,1) + t43 * rSges(4,2)) * t44 + (-t63 * t50 + (t43 * rSges(4,1) - t41 * rSges(4,2) + pkin(2)) * t47) * t42 + t59)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t24 * rSges(5,3) + t56) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t22 * rSges(5,3) + t53) + g(3) * (t17 * rSges(5,1) - t16 * rSges(5,2) - rSges(5,3) * t70 + t57)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t11 * t73 + t54) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t73 * t9 + t52) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t73 * t16 + t55)) - m(7) * (g(1) * (t11 * t74 + t3 * t64 + t4 * t75 + t54) + g(2) * (t1 * t64 + t2 * t75 + t74 * t9 + t52) + g(3) * (t16 * t74 + t64 * t7 + t75 * t8 + t55));
U  = t5;
