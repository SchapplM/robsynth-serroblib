% Calculate potential energy for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:54
% EndTime: 2019-03-08 22:35:54
% DurationCPUTime: 0.50s
% Computational Cost: add. (498->135), mult. (1189->177), div. (0->0), fcn. (1483->14), ass. (0->62)
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t49 = cos(pkin(6));
t46 = sin(pkin(6));
t56 = cos(qJ(2));
t80 = t46 * t56;
t28 = -t45 * t80 + t49 * t48;
t44 = sin(pkin(12));
t47 = cos(pkin(12));
t53 = sin(qJ(2));
t76 = t49 * t56;
t31 = -t44 * t76 - t47 * t53;
t82 = t46 * t48;
t22 = -t31 * t45 + t44 * t82;
t89 = rSges(6,3) + pkin(10);
t88 = pkin(11) + rSges(7,3);
t87 = cos(qJ(3));
t85 = t44 * t46;
t52 = sin(qJ(3));
t84 = t45 * t52;
t83 = t46 * t47;
t81 = t46 * t53;
t79 = t48 * t52;
t77 = t49 * t53;
t75 = rSges(5,3) + qJ(4);
t74 = t44 * pkin(1) + r_base(2);
t71 = qJ(1) + r_base(3);
t70 = t45 * t87;
t69 = t48 * t87;
t68 = t47 * pkin(1) + pkin(8) * t85 + r_base(1);
t67 = t49 * pkin(8) + t71;
t66 = t46 * t70;
t29 = -t44 * t53 + t47 * t76;
t21 = -t29 * t45 - t47 * t82;
t32 = -t44 * t77 + t47 * t56;
t65 = t32 * pkin(2) + t22 * pkin(9) + t68;
t12 = t32 * t87 + (t31 * t48 + t45 * t85) * t52;
t64 = t12 * pkin(3) + t65;
t63 = pkin(2) * t81 + t28 * pkin(9) + t67;
t20 = t49 * t84 + (t87 * t53 + t56 * t79) * t46;
t62 = t20 * pkin(3) + t63;
t11 = -t31 * t69 + t32 * t52 - t44 * t66;
t61 = t22 * pkin(4) + t11 * qJ(4) + t64;
t19 = -t49 * t70 + t52 * t81 - t69 * t80;
t60 = t28 * pkin(4) + t19 * qJ(4) + t62;
t30 = t44 * t56 + t47 * t77;
t59 = t30 * pkin(2) - pkin(8) * t83 + t21 * pkin(9) + t74;
t10 = t29 * t79 + t30 * t87 - t83 * t84;
t58 = t10 * pkin(3) + t59;
t9 = -t29 * t69 + t30 * t52 + t47 * t66;
t57 = t21 * pkin(4) + t9 * qJ(4) + t58;
t55 = cos(qJ(5));
t54 = cos(qJ(6));
t51 = sin(qJ(5));
t50 = sin(qJ(6));
t14 = t19 * t51 + t28 * t55;
t13 = -t19 * t55 + t28 * t51;
t4 = t11 * t51 + t22 * t55;
t3 = -t11 * t55 + t22 * t51;
t2 = t21 * t55 + t9 * t51;
t1 = t21 * t51 - t9 * t55;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t47 * rSges(2,1) - t44 * rSges(2,2) + r_base(1)) + g(2) * (t44 * rSges(2,1) + t47 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t71)) - m(3) * (g(1) * (t32 * rSges(3,1) + t31 * rSges(3,2) + t68) + g(2) * (t30 * rSges(3,1) + t29 * rSges(3,2) + t74) + g(3) * (t49 * rSges(3,3) + t67) + (g(1) * rSges(3,3) * t44 + g(3) * (rSges(3,1) * t53 + rSges(3,2) * t56) + g(2) * (-rSges(3,3) - pkin(8)) * t47) * t46) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t22 * rSges(4,3) + t65) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t21 * rSges(4,3) + t59) + g(3) * (t20 * rSges(4,1) - t19 * rSges(4,2) + t28 * rSges(4,3) + t63)) - m(5) * (g(1) * (t22 * rSges(5,1) - t12 * rSges(5,2) + t75 * t11 + t64) + g(2) * (t21 * rSges(5,1) - t10 * rSges(5,2) + t75 * t9 + t58) + g(3) * (t28 * rSges(5,1) - t20 * rSges(5,2) + t75 * t19 + t62)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t89 * t12 + t61) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t89 * t10 + t57) + g(3) * (t14 * rSges(6,1) - t13 * rSges(6,2) + t89 * t20 + t60)) - m(7) * (g(1) * (t4 * pkin(5) + t12 * pkin(10) + (t12 * t50 + t4 * t54) * rSges(7,1) + (t12 * t54 - t4 * t50) * rSges(7,2) + t88 * t3 + t61) + g(2) * (t2 * pkin(5) + t10 * pkin(10) + (t10 * t50 + t2 * t54) * rSges(7,1) + (t10 * t54 - t2 * t50) * rSges(7,2) + t88 * t1 + t57) + g(3) * (t14 * pkin(5) + t20 * pkin(10) + (t14 * t54 + t20 * t50) * rSges(7,1) + (-t14 * t50 + t20 * t54) * rSges(7,2) + t88 * t13 + t60));
U  = t5;
