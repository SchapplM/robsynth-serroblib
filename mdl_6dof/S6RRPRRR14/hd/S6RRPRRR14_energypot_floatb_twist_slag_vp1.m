% Calculate potential energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:10
% EndTime: 2019-03-09 14:57:11
% DurationCPUTime: 0.62s
% Computational Cost: add. (972->147), mult. (2565->199), div. (0->0), fcn. (3324->18), ass. (0->75)
t58 = sin(pkin(7));
t59 = sin(pkin(6));
t62 = cos(pkin(7));
t63 = cos(pkin(6));
t71 = cos(qJ(2));
t44 = -t58 * t59 * t71 + t62 * t63;
t68 = sin(qJ(1));
t93 = t68 * t71;
t67 = sin(qJ(2));
t72 = cos(qJ(1));
t95 = t67 * t72;
t47 = -t63 * t93 - t95;
t99 = t59 * t68;
t39 = -t47 * t58 + t62 * t99;
t101 = t58 * t63;
t56 = sin(pkin(14));
t60 = cos(pkin(14));
t96 = t62 * t71;
t36 = t60 * t101 + (-t56 * t67 + t60 * t96) * t59;
t57 = sin(pkin(8));
t61 = cos(pkin(8));
t23 = -t36 * t57 + t44 * t61;
t92 = t71 * t72;
t94 = t68 * t67;
t48 = -t63 * t94 + t92;
t82 = t47 * t62 + t58 * t99;
t28 = -t48 * t56 + t82 * t60;
t19 = -t28 * t57 + t39 * t61;
t46 = t63 * t95 + t93;
t45 = t63 * t92 - t94;
t98 = t59 * t72;
t83 = t45 * t62 - t58 * t98;
t26 = -t46 * t56 + t83 * t60;
t38 = -t45 * t58 - t62 * t98;
t18 = -t26 * t57 + t38 * t61;
t111 = rSges(6,3) + pkin(12);
t110 = pkin(13) + rSges(7,3);
t109 = cos(qJ(4));
t100 = t59 * t67;
t91 = pkin(9) + r_base(3);
t90 = t68 * pkin(1) + r_base(2);
t87 = t57 * t109;
t86 = t61 * t109;
t85 = t63 * pkin(10) + t91;
t84 = t72 * pkin(1) + pkin(10) * t99 + r_base(1);
t81 = t48 * pkin(2) + t39 * qJ(3) + t84;
t80 = pkin(2) * t100 + t44 * qJ(3) + t85;
t29 = t48 * t60 + t82 * t56;
t79 = t29 * pkin(3) + t19 * pkin(11) + t81;
t37 = t60 * t100 + (t59 * t96 + t101) * t56;
t78 = t37 * pkin(3) + t23 * pkin(11) + t80;
t66 = sin(qJ(4));
t12 = t29 * t109 + (t28 * t61 + t39 * t57) * t66;
t77 = t12 * pkin(4) + t79;
t76 = t46 * pkin(2) - pkin(10) * t98 + t38 * qJ(3) + t90;
t15 = t37 * t109 + (t36 * t61 + t44 * t57) * t66;
t75 = t15 * pkin(4) + t78;
t27 = t46 * t60 + t83 * t56;
t74 = t27 * pkin(3) + t18 * pkin(11) + t76;
t10 = t27 * t109 + (t26 * t61 + t38 * t57) * t66;
t73 = t10 * pkin(4) + t74;
t70 = cos(qJ(5));
t69 = cos(qJ(6));
t65 = sin(qJ(5));
t64 = sin(qJ(6));
t14 = -t36 * t86 + t37 * t66 - t44 * t87;
t11 = -t28 * t86 + t29 * t66 - t39 * t87;
t9 = -t26 * t86 + t27 * t66 - t38 * t87;
t6 = t15 * t70 + t23 * t65;
t5 = t15 * t65 - t23 * t70;
t4 = t12 * t70 + t19 * t65;
t3 = t12 * t65 - t19 * t70;
t2 = t10 * t70 + t18 * t65;
t1 = t10 * t65 - t18 * t70;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t72 - t68 * rSges(2,2) + r_base(1)) + g(2) * (t68 * rSges(2,1) + rSges(2,2) * t72 + r_base(2)) + g(3) * (rSges(2,3) + t91)) - m(3) * (g(1) * (rSges(3,1) * t48 + rSges(3,2) * t47 + t84) + g(2) * (t46 * rSges(3,1) + t45 * rSges(3,2) + t90) + g(3) * (rSges(3,3) * t63 + t85) + (g(1) * rSges(3,3) * t68 + g(3) * (rSges(3,1) * t67 + rSges(3,2) * t71) + g(2) * (-rSges(3,3) - pkin(10)) * t72) * t59) - m(4) * (g(1) * (rSges(4,1) * t29 + rSges(4,2) * t28 + rSges(4,3) * t39 + t81) + g(2) * (t27 * rSges(4,1) + t26 * rSges(4,2) + t38 * rSges(4,3) + t76) + g(3) * (rSges(4,1) * t37 + rSges(4,2) * t36 + rSges(4,3) * t44 + t80)) - m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t11 + rSges(5,3) * t19 + t79) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t18 * rSges(5,3) + t74) + g(3) * (rSges(5,1) * t15 - rSges(5,2) * t14 + rSges(5,3) * t23 + t78)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t111 * t11 + t77) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t111 * t9 + t73) + g(3) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t111 * t14 + t75)) - m(7) * (g(1) * (t77 + (t11 * t64 + t4 * t69) * rSges(7,1) + (t11 * t69 - t4 * t64) * rSges(7,2) + t110 * t3 + t11 * pkin(12) + t4 * pkin(5)) + g(2) * (t2 * pkin(5) + t9 * pkin(12) + (t2 * t69 + t64 * t9) * rSges(7,1) + (-t2 * t64 + t69 * t9) * rSges(7,2) + t110 * t1 + t73) + g(3) * (t110 * t5 + (t14 * t64 + t6 * t69) * rSges(7,1) + (t14 * t69 - t6 * t64) * rSges(7,2) + t14 * pkin(12) + t6 * pkin(5) + t75));
U  = t7;
