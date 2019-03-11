% Calculate Gravitation load on the joints for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:34
% EndTime: 2019-03-09 00:37:37
% DurationCPUTime: 0.85s
% Computational Cost: add. (833->150), mult. (1063->216), div. (0->0), fcn. (1220->14), ass. (0->74)
t122 = rSges(7,3) + pkin(11);
t61 = sin(pkin(6));
t64 = sin(qJ(2));
t104 = t61 * t64;
t59 = qJ(3) + qJ(4);
t54 = sin(t59);
t55 = cos(t59);
t93 = cos(pkin(6));
t121 = -t54 * t104 + t93 * t55;
t60 = sin(pkin(12));
t105 = t60 * t61;
t67 = cos(qJ(2));
t89 = t60 * t93;
t92 = cos(pkin(12));
t40 = -t64 * t89 + t67 * t92;
t120 = t55 * t105 - t40 * t54;
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t119 = rSges(7,1) * t65 - rSges(7,2) * t62 + pkin(5);
t118 = t62 * rSges(7,1) + t65 * rSges(7,2);
t56 = qJ(5) + t59;
t51 = sin(t56);
t52 = cos(t56);
t117 = t52 * t119 + t122 * t51;
t68 = -pkin(9) - pkin(8);
t83 = t93 * t92;
t38 = t60 * t67 + t64 * t83;
t88 = t61 * t92;
t11 = -t38 * t51 - t52 * t88;
t12 = t38 * t52 - t51 * t88;
t116 = t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t105 * t52 - t40 * t51;
t14 = t105 * t51 + t40 * t52;
t115 = t13 * rSges(6,1) - t14 * rSges(6,2);
t37 = t60 * t64 - t67 * t83;
t114 = g(2) * t37;
t113 = g(2) * t38;
t102 = t61 * t67;
t66 = cos(qJ(3));
t57 = t66 * pkin(3);
t46 = pkin(4) * t55 + t57;
t44 = pkin(2) + t46;
t112 = g(3) * t44 * t102;
t111 = g(3) * t61;
t110 = rSges(4,3) + pkin(8);
t103 = t61 * t66;
t99 = rSges(5,3) - t68;
t58 = -pkin(10) + t68;
t98 = -t37 * t44 - t38 * t58;
t39 = t64 * t92 + t67 * t89;
t97 = -t39 * t44 - t40 * t58;
t63 = sin(qJ(3));
t45 = -pkin(3) * t63 - pkin(4) * t54;
t96 = t46 * t105 + t40 * t45;
t28 = -t104 * t51 + t52 * t93;
t29 = t104 * t52 + t51 * t93;
t95 = t28 * rSges(6,1) - t29 * rSges(6,2);
t94 = t45 * t104 + t93 * t46;
t86 = t120 * pkin(4);
t85 = rSges(6,1) * t52 - rSges(6,2) * t51;
t84 = t121 * pkin(4);
t82 = rSges(4,1) * t66 - rSges(4,2) * t63 + pkin(2);
t80 = t103 * t60 - t40 * t63;
t79 = rSges(5,1) * t55 - rSges(5,2) * t54 + pkin(2) + t57;
t78 = t38 * t45 - t46 * t88;
t77 = -t38 * t54 - t55 * t88;
t76 = -t38 * t63 - t66 * t88;
t75 = -t104 * t63 + t66 * t93;
t74 = t119 * t11 + t122 * t12;
t73 = t119 * t13 + t122 * t14;
t72 = t119 * t28 + t122 * t29;
t71 = t77 * pkin(4);
t70 = g(1) * (t120 * rSges(5,1) + (-t105 * t54 - t40 * t55) * rSges(5,2)) + g(2) * (t77 * rSges(5,1) + (-t38 * t55 + t54 * t88) * rSges(5,2)) + g(3) * (t121 * rSges(5,1) + (-t104 * t55 - t54 * t93) * rSges(5,2));
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t39 - rSges(3,2) * t40) + g(2) * (-rSges(3,1) * t37 - rSges(3,2) * t38) + (rSges(3,1) * t67 - rSges(3,2) * t64) * t111) - m(4) * (g(1) * (t110 * t40 - t39 * t82) + t110 * t113 - t82 * t114 + (t110 * t64 + t67 * t82) * t111) - m(5) * (g(1) * (-t39 * t79 + t40 * t99) + t99 * t113 - t79 * t114 + (t64 * t99 + t67 * t79) * t111) - m(6) * (g(1) * (rSges(6,3) * t40 - t39 * t85 + t97) + g(2) * (rSges(6,3) * t38 - t37 * t85 + t98) + t112 + (t85 * t67 + (rSges(6,3) - t58) * t64) * t111) - m(7) * (g(2) * (t118 * t38 + t98) + t112 - t117 * t114 + ((-t58 + t118) * t64 + t117 * t67) * t111 + (-t117 * t39 + t118 * t40 + t97) * g(1)) -m(4) * (g(1) * (t80 * rSges(4,1) + (-t105 * t63 - t40 * t66) * rSges(4,2)) + g(2) * (t76 * rSges(4,1) + (-t38 * t66 + t63 * t88) * rSges(4,2)) + g(3) * (t75 * rSges(4,1) + (-t103 * t64 - t63 * t93) * rSges(4,2))) - m(5) * ((g(1) * t80 + g(2) * t76 + g(3) * t75) * pkin(3) + t70) - m(6) * (g(1) * (t96 + t115) + g(2) * (t78 + t116) + g(3) * (t94 + t95)) - m(7) * (g(1) * (t73 + t96) + g(2) * (t74 + t78) + g(3) * (t72 + t94)) -m(5) * t70 - m(6) * (g(1) * (t86 + t115) + g(2) * (t71 + t116) + g(3) * (t84 + t95)) - m(7) * (g(1) * (t73 + t86) + g(2) * (t71 + t74) + g(3) * (t72 + t84)) -m(6) * (g(1) * t115 + g(2) * t116 + g(3) * t95) - m(7) * (g(1) * t73 + g(2) * t74 + g(3) * t72) -m(7) * (g(1) * ((-t14 * t62 + t39 * t65) * rSges(7,1) + (-t14 * t65 - t39 * t62) * rSges(7,2)) + g(2) * ((-t12 * t62 + t37 * t65) * rSges(7,1) + (-t12 * t65 - t37 * t62) * rSges(7,2)) + g(3) * ((-t102 * t65 - t29 * t62) * rSges(7,1) + (t102 * t62 - t29 * t65) * rSges(7,2)))];
taug  = t1(:);
