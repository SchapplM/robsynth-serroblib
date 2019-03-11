% Calculate potential energy for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:27
% EndTime: 2019-03-09 00:28:27
% DurationCPUTime: 0.49s
% Computational Cost: add. (621->131), mult. (1532->172), div. (0->0), fcn. (1949->14), ass. (0->68)
t50 = sin(pkin(7));
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t51 = sin(pkin(6));
t60 = cos(qJ(2));
t86 = t51 * t60;
t70 = -t50 * t86 + t53 * t54;
t49 = sin(pkin(12));
t52 = cos(pkin(12));
t58 = sin(qJ(2));
t83 = t54 * t60;
t37 = -t49 * t83 - t52 * t58;
t88 = t51 * t53;
t71 = -t37 * t50 + t49 * t88;
t97 = rSges(7,1) + pkin(5);
t96 = rSges(7,2) + pkin(11);
t95 = rSges(5,3) + pkin(10);
t94 = rSges(6,3) + pkin(11);
t93 = cos(qJ(3));
t92 = cos(qJ(4));
t90 = t49 * t51;
t89 = t51 * t52;
t87 = t51 * t58;
t84 = t54 * t58;
t82 = rSges(7,3) + qJ(6);
t81 = pkin(1) * t49 + r_base(2);
t78 = qJ(1) + r_base(3);
t77 = t50 * t93;
t76 = t53 * t93;
t75 = pkin(1) * t52 + pkin(8) * t90 + r_base(1);
t74 = pkin(8) * t54 + t78;
t73 = t51 * t77;
t35 = -t49 * t58 + t52 * t83;
t72 = -t35 * t50 - t52 * t88;
t38 = -t49 * t84 + t52 * t60;
t69 = t38 * pkin(2) + pkin(9) * t71 + t75;
t57 = sin(qJ(3));
t21 = t38 * t93 + (t37 * t53 + t50 * t90) * t57;
t68 = pkin(3) * t21 + t69;
t67 = pkin(2) * t87 + pkin(9) * t70 + t74;
t29 = t54 * t50 * t57 + (t53 * t57 * t60 + t58 * t93) * t51;
t66 = pkin(3) * t29 + t67;
t56 = sin(qJ(4));
t12 = t21 * t92 + t56 * t71;
t20 = -t37 * t76 + t38 * t57 - t49 * t73;
t65 = pkin(4) * t12 + t20 * pkin(10) + t68;
t23 = t29 * t92 + t56 * t70;
t28 = -t54 * t77 + t57 * t87 - t76 * t86;
t64 = pkin(4) * t23 + t28 * pkin(10) + t66;
t36 = t49 * t60 + t52 * t84;
t63 = pkin(2) * t36 - pkin(8) * t89 + pkin(9) * t72 + t81;
t19 = t36 * t93 + (t35 * t53 - t50 * t89) * t57;
t62 = pkin(3) * t19 + t63;
t10 = t19 * t92 + t56 * t72;
t18 = -t35 * t76 + t36 * t57 + t52 * t73;
t61 = pkin(4) * t10 + t18 * pkin(10) + t62;
t59 = cos(qJ(5));
t55 = sin(qJ(5));
t22 = t29 * t56 - t70 * t92;
t11 = t21 * t56 - t71 * t92;
t9 = t19 * t56 - t72 * t92;
t8 = t23 * t59 + t28 * t55;
t7 = t23 * t55 - t28 * t59;
t4 = t12 * t59 + t20 * t55;
t3 = t12 * t55 - t20 * t59;
t2 = t10 * t59 + t18 * t55;
t1 = t10 * t55 - t18 * t59;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t52 - rSges(2,2) * t49 + r_base(1)) + g(2) * (rSges(2,1) * t49 + rSges(2,2) * t52 + r_base(2)) + g(3) * (rSges(2,3) + t78)) - m(3) * (g(1) * (rSges(3,1) * t38 + rSges(3,2) * t37 + t75) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + t81) + g(3) * (t54 * rSges(3,3) + t74) + (g(1) * rSges(3,3) * t49 + g(3) * (rSges(3,1) * t58 + rSges(3,2) * t60) + g(2) * (-rSges(3,3) - pkin(8)) * t52) * t51) - m(4) * (g(1) * (t21 * rSges(4,1) - t20 * rSges(4,2) + rSges(4,3) * t71 + t69) + g(2) * (t19 * rSges(4,1) - t18 * rSges(4,2) + rSges(4,3) * t72 + t63) + g(3) * (t29 * rSges(4,1) - t28 * rSges(4,2) + rSges(4,3) * t70 + t67)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t20 * t95 + t68) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t18 * t95 + t62) + g(3) * (t23 * rSges(5,1) - t22 * rSges(5,2) + t28 * t95 + t66)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t11 * t94 + t65) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t9 * t94 + t61) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t22 * t94 + t64)) - m(7) * (g(1) * (t11 * t96 + t3 * t82 + t4 * t97 + t65) + g(2) * (t1 * t82 + t2 * t97 + t9 * t96 + t61) + g(3) * (t22 * t96 + t7 * t82 + t8 * t97 + t64));
U  = t5;
