% Calculate Gravitation load on the joints for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:41
% EndTime: 2019-03-09 13:16:44
% DurationCPUTime: 0.81s
% Computational Cost: add. (558->143), mult. (475->194), div. (0->0), fcn. (429->12), ass. (0->82)
t115 = rSges(6,3) + pkin(9);
t49 = qJ(2) + pkin(11);
t44 = qJ(4) + t49;
t38 = sin(t44);
t39 = cos(t44);
t55 = cos(qJ(5));
t40 = t55 * pkin(5) + pkin(4);
t114 = t38 * rSges(7,3) + t39 * t40;
t113 = t39 * rSges(5,1) - t38 * rSges(5,2);
t42 = sin(t49);
t43 = cos(t49);
t56 = cos(qJ(2));
t47 = t56 * pkin(2);
t112 = t43 * rSges(4,1) - t42 * rSges(4,2) + t47;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t111 = g(1) * t57 + g(2) * t54;
t63 = t39 * pkin(4) + t115 * t38;
t50 = qJ(5) + qJ(6);
t46 = cos(t50);
t88 = t57 * t46;
t45 = sin(t50);
t93 = t54 * t45;
t5 = t39 * t93 + t88;
t89 = t57 * t45;
t92 = t54 * t46;
t6 = -t39 * t92 + t89;
t110 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t39 * t89 + t92;
t8 = t39 * t88 + t93;
t109 = t7 * rSges(7,1) - t8 * rSges(7,2);
t52 = sin(qJ(5));
t108 = pkin(5) * t52;
t105 = g(3) * t38;
t53 = sin(qJ(2));
t104 = t53 * pkin(2);
t103 = rSges(3,3) + pkin(7);
t102 = rSges(6,1) * t55;
t101 = rSges(7,1) * t46;
t100 = rSges(6,2) * t52;
t99 = rSges(7,2) * t45;
t58 = -pkin(10) - pkin(9);
t97 = t38 * t58;
t96 = t39 * t54;
t95 = t39 * t57;
t94 = t39 * t58;
t91 = t54 * t52;
t90 = t54 * t55;
t87 = t57 * t52;
t86 = t57 * t55;
t51 = -qJ(3) - pkin(7);
t85 = rSges(4,3) - t51;
t48 = -pkin(8) + t51;
t84 = rSges(5,3) - t48;
t79 = t38 * t99;
t83 = rSges(7,3) * t96 + t54 * t79;
t82 = rSges(7,3) * t95 + t57 * t79;
t81 = pkin(3) * t43 + t47;
t80 = t38 * t100;
t78 = t115 * t96 + t54 * t80;
t77 = t115 * t95 + t57 * t80;
t76 = -pkin(4) - t102;
t75 = -t48 + t108;
t74 = -t40 - t101;
t72 = t56 * rSges(3,1) - t53 * rSges(3,2);
t68 = -rSges(5,1) * t38 - rSges(5,2) * t39;
t67 = -rSges(7,1) * t45 - rSges(7,2) * t46;
t66 = pkin(1) + t72;
t11 = -t39 * t87 + t90;
t9 = t39 * t91 + t86;
t65 = pkin(1) + t112;
t64 = t114 + (t101 - t99) * t39;
t62 = t63 + (-t100 + t102) * t39;
t60 = -t97 + t114;
t23 = -pkin(3) * t42 - t104;
t22 = pkin(1) + t81;
t17 = t57 * t23;
t16 = t54 * t23;
t13 = t57 * t22;
t12 = t39 * t86 + t91;
t10 = -t39 * t90 + t87;
t1 = [-m(2) * (g(1) * (-t54 * rSges(2,1) - t57 * rSges(2,2)) + g(2) * (t57 * rSges(2,1) - t54 * rSges(2,2))) - m(3) * ((g(1) * t103 + g(2) * t66) * t57 + (-g(1) * t66 + g(2) * t103) * t54) - m(4) * ((g(1) * t85 + g(2) * t65) * t57 + (-g(1) * t65 + g(2) * t85) * t54) - m(5) * (g(2) * t13 + (g(1) * t84 + g(2) * t113) * t57 + (g(1) * (-t22 - t113) + g(2) * t84) * t54) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2)) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t13) + (-g(1) * t48 + g(2) * t63) * t57 + (g(1) * (-t22 - t63) - g(2) * t48) * t54) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2)) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t13) + (g(1) * t75 + g(2) * t60) * t57 + (g(1) * (-t22 - t60) + g(2) * t75) * t54) -m(3) * (g(3) * t72 + t111 * (-rSges(3,1) * t53 - rSges(3,2) * t56)) - m(4) * (g(3) * t112 + t111 * (-rSges(4,1) * t42 - rSges(4,2) * t43 - t104)) - m(5) * (g(1) * (t68 * t57 + t17) + g(2) * (t68 * t54 + t16) + g(3) * (t113 + t81)) - m(6) * (g(1) * (t17 + t77) + g(2) * (t16 + t78) + g(3) * (t62 + t81) + t111 * t38 * t76) - m(7) * (g(1) * (-t57 * t94 + t17 + t82) + g(2) * (-t54 * t94 + t16 + t83) + g(3) * (t64 + t81) + (-g(3) * t58 + t111 * t74) * t38) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t54 - g(2) * t57) -m(5) * g(3) * t113 - m(6) * (g(1) * t77 + g(2) * t78 + g(3) * t62) - m(7) * (g(1) * t82 + g(2) * t83 + g(3) * (t64 - t97)) + t111 * ((m(5) * rSges(5,2) + m(7) * t58) * t39 + (m(5) * rSges(5,1) - m(6) * t76 - m(7) * t74) * t38) -m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2))) - m(7) * (g(1) * (t11 * pkin(5) + t109) + g(2) * (-t9 * pkin(5) + t110)) + (-m(6) * (-rSges(6,1) * t52 - rSges(6,2) * t55) - m(7) * (t67 - t108)) * t105, -m(7) * (g(1) * t109 + g(2) * t110 + t67 * t105)];
taug  = t1(:);
