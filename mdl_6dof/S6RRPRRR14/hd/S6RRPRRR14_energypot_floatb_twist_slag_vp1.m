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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-12-10 18:09:45
% EndTime: 2018-12-10 18:09:46
% DurationCPUTime: 0.78s
% Computational Cost: add. (3330->169), mult. (3351->206), div. (0->0), fcn. (3324->30), ass. (0->97)
t74 = pkin(7) + pkin(14);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - pkin(14);
t67 = sin(t75);
t50 = t63 + t67 / 0.2e1;
t64 = cos(t75) / 0.2e1;
t68 = cos(t74);
t52 = t64 + t68 / 0.2e1;
t76 = pkin(6) + qJ(2);
t65 = sin(t76) / 0.2e1;
t77 = pkin(6) - qJ(2);
t69 = sin(t77);
t55 = t65 + t69 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t70 = cos(t76);
t59 = t66 - t70 / 0.2e1;
t78 = sin(pkin(14));
t85 = cos(pkin(6));
t32 = t50 * t85 + t52 * t55 - t59 * t78;
t84 = cos(pkin(7));
t119 = t84 * t85;
t80 = sin(pkin(7));
t44 = -t55 * t80 + t119;
t79 = sin(pkin(8));
t83 = cos(pkin(8));
t23 = -t32 * t79 + t44 * t83;
t81 = sin(pkin(6));
t90 = sin(qJ(1));
t121 = t81 * t90;
t58 = t66 + t70 / 0.2e1;
t89 = sin(qJ(2));
t95 = cos(qJ(1));
t47 = -t90 * t58 - t89 * t95;
t56 = t65 - t69 / 0.2e1;
t94 = cos(qJ(2));
t48 = -t90 * t56 + t94 * t95;
t28 = t121 * t50 + t47 * t52 - t48 * t78;
t113 = t84 * t121;
t39 = -t47 * t80 + t113;
t19 = -t28 * t79 + t39 * t83;
t120 = t81 * t95;
t45 = t58 * t95 - t90 * t89;
t46 = t56 * t95 + t90 * t94;
t26 = -t120 * t50 + t45 * t52 - t46 * t78;
t38 = -t120 * t84 - t45 * t80;
t18 = -t26 * t79 + t38 * t83;
t129 = rSges(6,3) + pkin(12);
t128 = pkin(13) + rSges(7,3);
t118 = qJ(3) * t80;
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t115 = pkin(9) + r_base(3);
t114 = t90 * pkin(1) + r_base(2);
t112 = t85 * pkin(10) + t115;
t111 = t95 * pkin(1) + pkin(10) * t121 + r_base(1);
t110 = cos(t116);
t109 = sin(t117);
t108 = cos(t117) / 0.2e1;
t107 = sin(t116) / 0.2e1;
t106 = t59 * pkin(2) + qJ(3) * t119 - t118 * t55 + t112;
t105 = t48 * pkin(2) + qJ(3) * t113 - t118 * t47 + t111;
t104 = t108 + t110 / 0.2e1;
t103 = t107 + t109 / 0.2e1;
t51 = t63 - t67 / 0.2e1;
t53 = t64 - t68 / 0.2e1;
t82 = cos(pkin(14));
t33 = t51 * t55 + t53 * t85 + t59 * t82;
t102 = t33 * pkin(3) + t23 * pkin(11) + t106;
t29 = t121 * t53 + t47 * t51 + t48 * t82;
t101 = t29 * pkin(3) + t19 * pkin(11) + t105;
t54 = t107 - t109 / 0.2e1;
t57 = t108 - t110 / 0.2e1;
t93 = cos(qJ(4));
t15 = t32 * t54 + t33 * t93 + t44 * t57;
t100 = t15 * pkin(4) + t102;
t12 = t28 * t54 + t29 * t93 + t39 * t57;
t99 = t12 * pkin(4) + t101;
t98 = t46 * pkin(2) - pkin(10) * t120 + qJ(3) * t38 + t114;
t27 = -t120 * t53 + t45 * t51 + t46 * t82;
t97 = t27 * pkin(3) + t18 * pkin(11) + t98;
t10 = t26 * t54 + t27 * t93 + t38 * t57;
t96 = t10 * pkin(4) + t97;
t92 = cos(qJ(5));
t91 = cos(qJ(6));
t88 = sin(qJ(4));
t87 = sin(qJ(5));
t86 = sin(qJ(6));
t14 = -t103 * t44 - t104 * t32 + t33 * t88;
t11 = -t103 * t39 - t104 * t28 + t29 * t88;
t9 = -t103 * t38 - t104 * t26 + t27 * t88;
t6 = t15 * t92 + t23 * t87;
t5 = t15 * t87 - t23 * t92;
t4 = t12 * t92 + t19 * t87;
t3 = t12 * t87 - t19 * t92;
t2 = t10 * t92 + t18 * t87;
t1 = t10 * t87 - t18 * t92;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t95 - t90 * rSges(2,2) + r_base(1)) + g(2) * (t90 * rSges(2,1) + rSges(2,2) * t95 + r_base(2)) + g(3) * (rSges(2,3) + t115)) - m(3) * (g(1) * (rSges(3,1) * t48 + rSges(3,2) * t47 + rSges(3,3) * t121 + t111) + g(2) * (t46 * rSges(3,1) + t45 * rSges(3,2) + (-rSges(3,3) - pkin(10)) * t120 + t114) + g(3) * (rSges(3,1) * t59 + rSges(3,2) * t55 + rSges(3,3) * t85 + t112)) - m(4) * (g(1) * (rSges(4,1) * t29 + rSges(4,2) * t28 + rSges(4,3) * t39 + t105) + g(2) * (t27 * rSges(4,1) + t26 * rSges(4,2) + t38 * rSges(4,3) + t98) + g(3) * (rSges(4,1) * t33 + rSges(4,2) * t32 + rSges(4,3) * t44 + t106)) - m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t11 + rSges(5,3) * t19 + t101) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t18 * rSges(5,3) + t97) + g(3) * (rSges(5,1) * t15 - rSges(5,2) * t14 + rSges(5,3) * t23 + t102)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t11 * t129 + t99) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t129 * t9 + t96) + g(3) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t129 * t14 + t100)) - m(7) * (g(1) * (t99 + (t11 * t86 + t4 * t91) * rSges(7,1) + (t11 * t91 - t4 * t86) * rSges(7,2) + t11 * pkin(12) + t4 * pkin(5) + t128 * t3) + g(2) * (t2 * pkin(5) + t9 * pkin(12) + (t2 * t91 + t86 * t9) * rSges(7,1) + (-t2 * t86 + t9 * t91) * rSges(7,2) + t128 * t1 + t96) + g(3) * (t100 + (t14 * t86 + t6 * t91) * rSges(7,1) + (t14 * t91 - t6 * t86) * rSges(7,2) + t14 * pkin(12) + t6 * pkin(5) + t128 * t5));
U  = t7;
