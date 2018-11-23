% Calculate potential energy for
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:06
% EndTime: 2018-11-23 10:29:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (3330->167), mult. (3351->203), div. (0->0), fcn. (3324->30), ass. (0->94)
t80 = sin(pkin(6));
t89 = sin(qJ(1));
t120 = t80 * t89;
t77 = pkin(6) - qJ(2);
t66 = cos(t77) / 0.2e1;
t76 = pkin(6) + qJ(2);
t70 = cos(t76);
t58 = t66 + t70 / 0.2e1;
t88 = sin(qJ(2));
t95 = cos(qJ(1));
t47 = -t89 * t58 - t88 * t95;
t79 = sin(pkin(7));
t82 = cos(pkin(7));
t39 = t82 * t120 - t47 * t79;
t64 = sin(t76) / 0.2e1;
t68 = sin(t77);
t53 = t64 + t68 / 0.2e1;
t83 = cos(pkin(6));
t44 = -t53 * t79 + t82 * t83;
t74 = pkin(7) + qJ(3);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - qJ(3);
t67 = sin(t75);
t51 = t63 + t67 / 0.2e1;
t65 = cos(t75) / 0.2e1;
t69 = cos(t74);
t56 = t65 + t69 / 0.2e1;
t59 = t66 - t70 / 0.2e1;
t87 = sin(qJ(3));
t32 = t51 * t83 + t53 * t56 - t59 * t87;
t78 = sin(pkin(8));
t81 = cos(pkin(8));
t23 = -t32 * t78 + t44 * t81;
t54 = t64 - t68 / 0.2e1;
t94 = cos(qJ(2));
t48 = -t89 * t54 + t94 * t95;
t28 = t120 * t51 + t47 * t56 - t48 * t87;
t19 = -t28 * t78 + t39 * t81;
t119 = t80 * t95;
t45 = t58 * t95 - t89 * t88;
t46 = t54 * t95 + t89 * t94;
t26 = -t119 * t51 + t45 * t56 - t46 * t87;
t38 = -t119 * t82 - t45 * t79;
t18 = -t26 * t78 + t38 * t81;
t130 = rSges(6,3) + pkin(13);
t129 = pkin(14) + rSges(7,3);
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t115 = pkin(9) + r_base(3);
t114 = t89 * pkin(1) + r_base(2);
t112 = t83 * pkin(10) + t115;
t111 = t95 * pkin(1) + pkin(10) * t120 + r_base(1);
t110 = cos(t116);
t109 = sin(t117);
t108 = cos(t117) / 0.2e1;
t107 = sin(t116) / 0.2e1;
t106 = t59 * pkin(2) + t44 * pkin(11) + t112;
t105 = t48 * pkin(2) + t39 * pkin(11) + t111;
t104 = t108 + t110 / 0.2e1;
t103 = t107 + t109 / 0.2e1;
t52 = t63 - t67 / 0.2e1;
t57 = t65 - t69 / 0.2e1;
t93 = cos(qJ(3));
t33 = t52 * t53 + t57 * t83 + t59 * t93;
t102 = t33 * pkin(3) + t23 * pkin(12) + t106;
t29 = t120 * t57 + t47 * t52 + t48 * t93;
t101 = t29 * pkin(3) + t19 * pkin(12) + t105;
t50 = t107 - t109 / 0.2e1;
t55 = t108 - t110 / 0.2e1;
t92 = cos(qJ(4));
t15 = t32 * t50 + t33 * t92 + t44 * t55;
t100 = t15 * pkin(4) + t102;
t12 = t28 * t50 + t29 * t92 + t39 * t55;
t99 = t12 * pkin(4) + t101;
t98 = t46 * pkin(2) - pkin(10) * t119 + pkin(11) * t38 + t114;
t27 = -t119 * t57 + t45 * t52 + t46 * t93;
t97 = t27 * pkin(3) + t18 * pkin(12) + t98;
t10 = t26 * t50 + t27 * t92 + t38 * t55;
t96 = t10 * pkin(4) + t97;
t91 = cos(qJ(5));
t90 = cos(qJ(6));
t86 = sin(qJ(4));
t85 = sin(qJ(5));
t84 = sin(qJ(6));
t14 = -t103 * t44 - t104 * t32 + t33 * t86;
t11 = -t103 * t39 - t104 * t28 + t29 * t86;
t9 = -t103 * t38 - t104 * t26 + t27 * t86;
t6 = t15 * t91 + t23 * t85;
t5 = t15 * t85 - t23 * t91;
t4 = t12 * t91 + t19 * t85;
t3 = t12 * t85 - t19 * t91;
t2 = t10 * t91 + t18 * t85;
t1 = t10 * t85 - t18 * t91;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t95 - t89 * rSges(2,2) + r_base(1)) + g(2) * (t89 * rSges(2,1) + rSges(2,2) * t95 + r_base(2)) + g(3) * (rSges(2,3) + t115)) - m(3) * (g(1) * (rSges(3,1) * t48 + rSges(3,2) * t47 + rSges(3,3) * t120 + t111) + g(2) * (t46 * rSges(3,1) + t45 * rSges(3,2) + (-rSges(3,3) - pkin(10)) * t119 + t114) + g(3) * (rSges(3,1) * t59 + rSges(3,2) * t53 + rSges(3,3) * t83 + t112)) - m(4) * (g(1) * (rSges(4,1) * t29 + rSges(4,2) * t28 + rSges(4,3) * t39 + t105) + g(2) * (t27 * rSges(4,1) + t26 * rSges(4,2) + t38 * rSges(4,3) + t98) + g(3) * (rSges(4,1) * t33 + rSges(4,2) * t32 + rSges(4,3) * t44 + t106)) - m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t11 + rSges(5,3) * t19 + t101) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t18 * rSges(5,3) + t97) + g(3) * (rSges(5,1) * t15 - rSges(5,2) * t14 + rSges(5,3) * t23 + t102)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t11 * t130 + t99) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t130 * t9 + t96) + g(3) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t130 * t14 + t100)) - m(7) * (g(1) * (t99 + (t11 * t84 + t4 * t90) * rSges(7,1) + (t11 * t90 - t4 * t84) * rSges(7,2) + t11 * pkin(13) + t4 * pkin(5) + t129 * t3) + g(2) * (t2 * pkin(5) + t9 * pkin(13) + (t2 * t90 + t84 * t9) * rSges(7,1) + (-t2 * t84 + t9 * t90) * rSges(7,2) + t129 * t1 + t96) + g(3) * (t100 + (t14 * t84 + t6 * t90) * rSges(7,1) + (t14 * t90 - t6 * t84) * rSges(7,2) + t14 * pkin(13) + t6 * pkin(5) + t129 * t5));
U  = t7;
