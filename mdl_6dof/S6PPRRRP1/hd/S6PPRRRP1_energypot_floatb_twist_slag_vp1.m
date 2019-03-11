% Calculate potential energy for
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:13
% EndTime: 2019-03-08 18:52:14
% DurationCPUTime: 0.51s
% Computational Cost: add. (572->137), mult. (1388->180), div. (0->0), fcn. (1753->14), ass. (0->65)
t46 = sin(pkin(7));
t50 = cos(pkin(7));
t51 = cos(pkin(6));
t47 = sin(pkin(6));
t48 = cos(pkin(12));
t81 = t47 * t48;
t63 = -t46 * t81 + t50 * t51;
t44 = sin(pkin(12));
t49 = cos(pkin(11));
t45 = sin(pkin(11));
t83 = t45 * t51;
t31 = -t44 * t49 - t48 * t83;
t80 = t47 * t50;
t64 = -t31 * t46 + t45 * t80;
t89 = rSges(5,3) + pkin(9);
t88 = rSges(6,3) + pkin(10);
t87 = cos(qJ(3));
t86 = cos(qJ(4));
t84 = t44 * t47;
t82 = t46 * t47;
t79 = t49 * t51;
t77 = rSges(7,3) + qJ(6) + pkin(10);
t76 = t47 * qJ(2);
t75 = pkin(1) * t45 + r_base(2);
t72 = qJ(1) + r_base(3);
t53 = sin(qJ(5));
t71 = pkin(5) * t53 + pkin(9);
t70 = t46 * t87;
t69 = t50 * t87;
t68 = pkin(1) * t49 + t45 * t76 + r_base(1);
t67 = qJ(2) * t51 + t72;
t66 = t47 * t70;
t29 = -t44 * t45 + t48 * t79;
t65 = -t29 * t46 - t49 * t80;
t32 = -t44 * t83 + t48 * t49;
t62 = t32 * pkin(2) + pkin(8) * t64 + t68;
t55 = sin(qJ(3));
t16 = t32 * t87 + (t31 * t50 + t45 * t82) * t55;
t61 = pkin(3) * t16 + t62;
t60 = pkin(2) * t84 + pkin(8) * t63 + t67;
t23 = t51 * t46 * t55 + (t48 * t50 * t55 + t44 * t87) * t47;
t59 = pkin(3) * t23 + t60;
t30 = t44 * t79 + t45 * t48;
t58 = pkin(2) * t30 + pkin(8) * t65 - t49 * t76 + t75;
t14 = t30 * t87 + (t29 * t50 - t49 * t82) * t55;
t57 = pkin(3) * t14 + t58;
t56 = cos(qJ(5));
t54 = sin(qJ(4));
t40 = pkin(5) * t56 + pkin(4);
t22 = -t51 * t70 + t55 * t84 - t69 * t81;
t18 = t23 * t86 + t54 * t63;
t17 = t23 * t54 - t63 * t86;
t15 = -t31 * t69 + t32 * t55 - t45 * t66;
t13 = -t29 * t69 + t30 * t55 + t49 * t66;
t10 = t18 * t56 + t22 * t53;
t9 = -t18 * t53 + t22 * t56;
t8 = t16 * t86 + t54 * t64;
t7 = t16 * t54 - t64 * t86;
t6 = t14 * t86 + t54 * t65;
t5 = t14 * t54 - t65 * t86;
t4 = t15 * t53 + t56 * t8;
t3 = t15 * t56 - t53 * t8;
t2 = t13 * t53 + t56 * t6;
t1 = t13 * t56 - t53 * t6;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t49 - rSges(2,2) * t45 + r_base(1)) + g(2) * (rSges(2,1) * t45 + rSges(2,2) * t49 + r_base(2)) + g(3) * (rSges(2,3) + t72)) - m(3) * (g(1) * (rSges(3,1) * t32 + rSges(3,2) * t31 + t68) + g(2) * (t30 * rSges(3,1) + t29 * rSges(3,2) + t75) + g(3) * (t51 * rSges(3,3) + t67) + (g(1) * rSges(3,3) * t45 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t48) + g(2) * (-rSges(3,3) - qJ(2)) * t49) * t47) - m(4) * (g(1) * (t16 * rSges(4,1) - t15 * rSges(4,2) + rSges(4,3) * t64 + t62) + g(2) * (t14 * rSges(4,1) - t13 * rSges(4,2) + rSges(4,3) * t65 + t58) + g(3) * (t23 * rSges(4,1) - t22 * rSges(4,2) + rSges(4,3) * t63 + t60)) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t15 * t89 + t61) + g(2) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t13 * t89 + t57) + g(3) * (t18 * rSges(5,1) - t17 * rSges(5,2) + t22 * t89 + t59)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t8 * pkin(4) + t15 * pkin(9) + t7 * t88 + t61) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t6 * pkin(4) + t13 * pkin(9) + t5 * t88 + t57) + g(3) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t18 * pkin(4) + t22 * pkin(9) + t17 * t88 + t59)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t15 * t71 + t8 * t40 + t7 * t77 + t61) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t13 * t71 + t6 * t40 + t5 * t77 + t57) + g(3) * (t10 * rSges(7,1) + t9 * rSges(7,2) + t17 * t77 + t18 * t40 + t22 * t71 + t59));
U  = t11;
