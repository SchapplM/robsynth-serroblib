% Calculate potential energy for
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:04
% EndTime: 2019-03-09 00:20:05
% DurationCPUTime: 0.51s
% Computational Cost: add. (572->137), mult. (1388->180), div. (0->0), fcn. (1753->14), ass. (0->65)
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t49 = cos(pkin(6));
t46 = sin(pkin(6));
t56 = cos(qJ(2));
t80 = t46 * t56;
t63 = -t45 * t80 + t48 * t49;
t44 = sin(pkin(12));
t47 = cos(pkin(12));
t54 = sin(qJ(2));
t77 = t49 * t56;
t31 = -t44 * t77 - t47 * t54;
t82 = t46 * t48;
t64 = -t31 * t45 + t44 * t82;
t89 = rSges(5,3) + pkin(10);
t88 = rSges(6,3) + pkin(11);
t87 = cos(qJ(3));
t86 = cos(qJ(4));
t84 = t44 * t46;
t83 = t46 * t47;
t81 = t46 * t54;
t78 = t49 * t54;
t76 = rSges(7,3) + qJ(6) + pkin(11);
t75 = pkin(1) * t44 + r_base(2);
t72 = qJ(1) + r_base(3);
t51 = sin(qJ(5));
t71 = pkin(5) * t51 + pkin(10);
t70 = t45 * t87;
t69 = t48 * t87;
t68 = pkin(1) * t47 + pkin(8) * t84 + r_base(1);
t67 = pkin(8) * t49 + t72;
t66 = t46 * t70;
t29 = -t44 * t54 + t47 * t77;
t65 = -t29 * t45 - t47 * t82;
t32 = -t44 * t78 + t47 * t56;
t62 = t32 * pkin(2) + pkin(9) * t64 + t68;
t53 = sin(qJ(3));
t16 = t32 * t87 + (t31 * t48 + t45 * t84) * t53;
t61 = pkin(3) * t16 + t62;
t60 = pkin(2) * t81 + pkin(9) * t63 + t67;
t23 = t49 * t45 * t53 + (t48 * t53 * t56 + t54 * t87) * t46;
t59 = pkin(3) * t23 + t60;
t30 = t44 * t56 + t47 * t78;
t58 = pkin(2) * t30 - pkin(8) * t83 + pkin(9) * t65 + t75;
t14 = t30 * t87 + (t29 * t48 - t45 * t83) * t53;
t57 = pkin(3) * t14 + t58;
t55 = cos(qJ(5));
t52 = sin(qJ(4));
t40 = pkin(5) * t55 + pkin(4);
t22 = -t49 * t70 + t53 * t81 - t69 * t80;
t18 = t23 * t86 + t52 * t63;
t17 = t23 * t52 - t63 * t86;
t15 = -t31 * t69 + t32 * t53 - t44 * t66;
t13 = -t29 * t69 + t30 * t53 + t47 * t66;
t10 = t16 * t86 + t52 * t64;
t9 = t16 * t52 - t64 * t86;
t8 = t14 * t86 + t52 * t65;
t7 = t14 * t52 - t65 * t86;
t6 = t18 * t55 + t22 * t51;
t5 = -t18 * t51 + t22 * t55;
t4 = t10 * t55 + t15 * t51;
t3 = -t10 * t51 + t15 * t55;
t2 = t13 * t51 + t55 * t8;
t1 = t13 * t55 - t51 * t8;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t47 - rSges(2,2) * t44 + r_base(1)) + g(2) * (rSges(2,1) * t44 + rSges(2,2) * t47 + r_base(2)) + g(3) * (rSges(2,3) + t72)) - m(3) * (g(1) * (rSges(3,1) * t32 + rSges(3,2) * t31 + t68) + g(2) * (t30 * rSges(3,1) + t29 * rSges(3,2) + t75) + g(3) * (t49 * rSges(3,3) + t67) + (g(1) * rSges(3,3) * t44 + g(3) * (rSges(3,1) * t54 + rSges(3,2) * t56) + g(2) * (-rSges(3,3) - pkin(8)) * t47) * t46) - m(4) * (g(1) * (t16 * rSges(4,1) - t15 * rSges(4,2) + rSges(4,3) * t64 + t62) + g(2) * (t14 * rSges(4,1) - t13 * rSges(4,2) + rSges(4,3) * t65 + t58) + g(3) * (t23 * rSges(4,1) - t22 * rSges(4,2) + rSges(4,3) * t63 + t60)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t15 * t89 + t61) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t13 * t89 + t57) + g(3) * (t18 * rSges(5,1) - t17 * rSges(5,2) + t22 * t89 + t59)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * pkin(4) + t15 * pkin(10) + t88 * t9 + t61) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) + t13 * pkin(10) + t7 * t88 + t57) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t18 * pkin(4) + t22 * pkin(10) + t17 * t88 + t59)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t10 * t40 + t15 * t71 + t76 * t9 + t61) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t13 * t71 + t8 * t40 + t7 * t76 + t57) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t17 * t76 + t18 * t40 + t22 * t71 + t59));
U  = t11;
