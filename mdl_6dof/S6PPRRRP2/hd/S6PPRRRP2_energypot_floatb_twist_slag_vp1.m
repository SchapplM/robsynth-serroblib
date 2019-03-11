% Calculate potential energy for
% S6PPRRRP2
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:55:55
% EndTime: 2019-03-08 18:55:56
% DurationCPUTime: 0.50s
% Computational Cost: add. (621->131), mult. (1532->172), div. (0->0), fcn. (1949->14), ass. (0->68)
t51 = sin(pkin(7));
t55 = cos(pkin(7));
t56 = cos(pkin(6));
t52 = sin(pkin(6));
t53 = cos(pkin(12));
t87 = t52 * t53;
t70 = -t51 * t87 + t55 * t56;
t49 = sin(pkin(12));
t54 = cos(pkin(11));
t50 = sin(pkin(11));
t89 = t50 * t56;
t37 = -t49 * t54 - t53 * t89;
t86 = t52 * t55;
t71 = -t37 * t51 + t50 * t86;
t97 = rSges(7,1) + pkin(5);
t96 = rSges(7,2) + pkin(10);
t95 = rSges(5,3) + pkin(9);
t94 = rSges(6,3) + pkin(10);
t93 = cos(qJ(3));
t92 = cos(qJ(4));
t90 = t49 * t52;
t88 = t51 * t52;
t85 = t54 * t56;
t83 = t52 * qJ(2);
t82 = rSges(7,3) + qJ(6);
t81 = pkin(1) * t50 + r_base(2);
t78 = qJ(1) + r_base(3);
t77 = t51 * t93;
t76 = t55 * t93;
t75 = pkin(1) * t54 + t50 * t83 + r_base(1);
t74 = qJ(2) * t56 + t78;
t73 = t52 * t77;
t35 = -t49 * t50 + t53 * t85;
t72 = -t35 * t51 - t54 * t86;
t38 = -t49 * t89 + t53 * t54;
t69 = t38 * pkin(2) + pkin(8) * t71 + t75;
t59 = sin(qJ(3));
t21 = t38 * t93 + (t37 * t55 + t50 * t88) * t59;
t68 = pkin(3) * t21 + t69;
t67 = pkin(2) * t90 + pkin(8) * t70 + t74;
t29 = t56 * t51 * t59 + (t53 * t55 * t59 + t49 * t93) * t52;
t66 = pkin(3) * t29 + t67;
t58 = sin(qJ(4));
t10 = t21 * t92 + t58 * t71;
t20 = -t37 * t76 + t38 * t59 - t50 * t73;
t65 = pkin(4) * t10 + t20 * pkin(9) + t68;
t23 = t29 * t92 + t58 * t70;
t28 = -t56 * t77 + t59 * t90 - t76 * t87;
t64 = pkin(4) * t23 + t28 * pkin(9) + t66;
t36 = t49 * t85 + t50 * t53;
t63 = pkin(2) * t36 + pkin(8) * t72 - t54 * t83 + t81;
t19 = t36 * t93 + (t35 * t55 - t54 * t88) * t59;
t62 = pkin(3) * t19 + t63;
t18 = -t35 * t76 + t36 * t59 + t54 * t73;
t8 = t19 * t92 + t58 * t72;
t61 = pkin(4) * t8 + t18 * pkin(9) + t62;
t60 = cos(qJ(5));
t57 = sin(qJ(5));
t22 = t29 * t58 - t70 * t92;
t12 = t23 * t60 + t28 * t57;
t11 = t23 * t57 - t28 * t60;
t9 = t21 * t58 - t71 * t92;
t7 = t19 * t58 - t72 * t92;
t4 = t10 * t60 + t20 * t57;
t3 = t10 * t57 - t20 * t60;
t2 = t18 * t57 + t60 * t8;
t1 = -t18 * t60 + t57 * t8;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t54 - rSges(2,2) * t50 + r_base(1)) + g(2) * (rSges(2,1) * t50 + rSges(2,2) * t54 + r_base(2)) + g(3) * (rSges(2,3) + t78)) - m(3) * (g(1) * (rSges(3,1) * t38 + rSges(3,2) * t37 + t75) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + t81) + g(3) * (t56 * rSges(3,3) + t74) + (g(1) * rSges(3,3) * t50 + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t53) + g(2) * (-rSges(3,3) - qJ(2)) * t54) * t52) - m(4) * (g(1) * (t21 * rSges(4,1) - t20 * rSges(4,2) + rSges(4,3) * t71 + t69) + g(2) * (t19 * rSges(4,1) - t18 * rSges(4,2) + rSges(4,3) * t72 + t63) + g(3) * (t29 * rSges(4,1) - t28 * rSges(4,2) + rSges(4,3) * t70 + t67)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t20 * t95 + t68) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t18 * t95 + t62) + g(3) * (t23 * rSges(5,1) - t22 * rSges(5,2) + t28 * t95 + t66)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t9 * t94 + t65) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t7 * t94 + t61) + g(3) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t22 * t94 + t64)) - m(7) * (g(1) * (t3 * t82 + t4 * t97 + t9 * t96 + t65) + g(2) * (t1 * t82 + t2 * t97 + t7 * t96 + t61) + g(3) * (t11 * t82 + t12 * t97 + t22 * t96 + t64));
U  = t5;
