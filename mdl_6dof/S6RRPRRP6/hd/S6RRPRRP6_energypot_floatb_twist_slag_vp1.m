% Calculate potential energy for
% S6RRPRRP6
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:39
% EndTime: 2019-03-09 12:07:39
% DurationCPUTime: 0.51s
% Computational Cost: add. (439->124), mult. (938->153), div. (0->0), fcn. (1177->12), ass. (0->60)
t43 = sin(pkin(11));
t45 = cos(pkin(11));
t49 = sin(qJ(2));
t53 = cos(qJ(2));
t33 = -t49 * t43 + t53 * t45;
t82 = pkin(2) * t49;
t81 = rSges(7,1) + pkin(5);
t80 = rSges(7,2) + pkin(10);
t78 = rSges(5,3) + pkin(9);
t77 = rSges(6,3) + pkin(10);
t44 = sin(pkin(6));
t50 = sin(qJ(1));
t75 = t50 * t44;
t74 = t50 * t49;
t73 = t50 * t53;
t54 = cos(qJ(1));
t71 = t54 * t44;
t70 = t54 * t49;
t69 = t54 * t53;
t68 = rSges(7,3) + qJ(6);
t67 = pkin(7) + r_base(3);
t40 = t53 * pkin(2) + pkin(1);
t66 = t54 * t40 + r_base(1);
t46 = cos(pkin(6));
t65 = t46 * pkin(8) + t67;
t31 = t46 * t82 + (-pkin(8) - qJ(3)) * t44;
t64 = t54 * t31 + t50 * t40 + r_base(2);
t62 = t53 * t43 + t49 * t45;
t30 = t62 * t46;
t18 = t54 * t30 + t50 * t33;
t63 = t18 * pkin(3) + t64;
t61 = t46 * qJ(3) + t44 * t82 + t65;
t20 = -t50 * t30 + t54 * t33;
t60 = t20 * pkin(3) - t50 * t31 + t66;
t29 = t62 * t44;
t59 = t29 * pkin(3) + t61;
t58 = t33 * t46;
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t10 = t18 * t52 - t48 * t71;
t17 = -t50 * t62 + t54 * t58;
t57 = t10 * pkin(4) - t17 * pkin(9) + t63;
t12 = t20 * t52 + t48 * t75;
t19 = -t50 * t58 - t54 * t62;
t56 = t12 * pkin(4) - t19 * pkin(9) + t60;
t23 = t29 * t52 + t46 * t48;
t28 = t33 * t44;
t55 = t23 * pkin(4) - t28 * pkin(9) + t59;
t51 = cos(qJ(5));
t47 = sin(qJ(5));
t22 = t29 * t48 - t46 * t52;
t11 = t20 * t48 - t52 * t75;
t9 = t18 * t48 + t52 * t71;
t6 = t23 * t51 - t28 * t47;
t5 = t23 * t47 + t28 * t51;
t4 = t12 * t51 - t19 * t47;
t3 = t12 * t47 + t19 * t51;
t2 = t10 * t51 - t17 * t47;
t1 = t10 * t47 + t17 * t51;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t54 * rSges(2,1) - t50 * rSges(2,2) + r_base(1)) + g(2) * (t50 * rSges(2,1) + t54 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t67)) - m(3) * (g(1) * (t54 * pkin(1) + r_base(1) + (-t46 * t74 + t69) * rSges(3,1) + (-t46 * t73 - t70) * rSges(3,2)) + g(2) * (t50 * pkin(1) + r_base(2) + (t46 * t70 + t73) * rSges(3,1) + (t46 * t69 - t74) * rSges(3,2)) + g(3) * (t46 * rSges(3,3) + t65) + (g(3) * (rSges(3,1) * t49 + rSges(3,2) * t53) + (g(1) * t50 - g(2) * t54) * (rSges(3,3) + pkin(8))) * t44) - m(4) * (g(1) * (t20 * rSges(4,1) + t19 * rSges(4,2) + (rSges(4,3) * t44 - t31) * t50 + t66) + g(2) * (t18 * rSges(4,1) + t17 * rSges(4,2) - rSges(4,3) * t71 + t64) + g(3) * (t29 * rSges(4,1) + t28 * rSges(4,2) + t46 * rSges(4,3) + t61)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) - t19 * t78 + t60) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) - t17 * t78 + t63) + g(3) * (t23 * rSges(5,1) - t22 * rSges(5,2) - t28 * t78 + t59)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t11 * t77 + t56) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t77 * t9 + t57) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t22 * t77 + t55)) - m(7) * (g(1) * (t11 * t80 + t3 * t68 + t4 * t81 + t56) + g(2) * (t1 * t68 + t2 * t81 + t80 * t9 + t57) + g(3) * (t22 * t80 + t5 * t68 + t6 * t81 + t55));
U  = t7;
