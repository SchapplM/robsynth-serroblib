% Calculate Gravitation load on the joints for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:19
% EndTime: 2019-03-08 19:07:22
% DurationCPUTime: 1.01s
% Computational Cost: add. (1747->176), mult. (4990->290), div. (0->0), fcn. (6532->18), ass. (0->95)
t101 = cos(pkin(8));
t51 = sin(pkin(8));
t108 = sin(qJ(3));
t110 = cos(qJ(3));
t100 = cos(pkin(13));
t102 = cos(pkin(7));
t103 = cos(pkin(6));
t95 = sin(pkin(14));
t96 = sin(pkin(13));
t79 = t96 * t95;
t99 = cos(pkin(14));
t86 = t100 * t99;
t71 = t103 * t86 - t79;
t97 = sin(pkin(7));
t98 = sin(pkin(6));
t83 = t98 * t97;
t118 = -t100 * t83 + t71 * t102;
t80 = t96 * t99;
t84 = t100 * t95;
t72 = t103 * t84 + t80;
t59 = t108 * t72 - t110 * t118;
t85 = t100 * t98;
t65 = -t102 * t85 - t71 * t97;
t121 = t59 * t101 - t65 * t51;
t73 = -t103 * t80 - t84;
t82 = t98 * t96;
t117 = t73 * t102 + t97 * t82;
t74 = -t103 * t79 + t86;
t60 = t108 * t74 - t110 * t117;
t66 = t102 * t82 - t73 * t97;
t120 = t60 * t101 - t66 * t51;
t116 = t102 * t99 * t98 + t103 * t97;
t81 = t98 * t95;
t64 = t108 * t81 - t110 * t116;
t70 = t102 * t103 - t83 * t99;
t119 = t64 * t101 - t70 * t51;
t115 = pkin(10) * t51;
t56 = cos(qJ(5));
t114 = t56 * pkin(5);
t112 = rSges(6,3) + pkin(11);
t111 = rSges(7,3) + pkin(12);
t109 = cos(qJ(4));
t53 = sin(qJ(5));
t107 = t51 * t53;
t106 = t51 * t56;
t52 = sin(qJ(6));
t105 = t52 * t56;
t55 = cos(qJ(6));
t104 = t55 * t56;
t54 = sin(qJ(4));
t94 = t54 * t101;
t93 = -m(3) - m(4) - m(5) - m(6) - m(7);
t92 = t101 * t109;
t91 = -rSges(6,1) * t56 + rSges(6,2) * t53;
t39 = t108 * t118 + t110 * t72;
t20 = -t109 * t59 - t39 * t94;
t37 = t59 * pkin(3);
t90 = t20 * pkin(4) + t115 * t39 - t37;
t40 = t108 * t117 + t110 * t74;
t22 = -t109 * t60 - t40 * t94;
t38 = t60 * pkin(3);
t89 = t22 * pkin(4) + t115 * t40 - t38;
t47 = t108 * t116 + t110 * t81;
t32 = -t109 * t64 - t47 * t94;
t46 = t64 * pkin(3);
t88 = t32 * pkin(4) + t115 * t47 - t46;
t78 = t55 * rSges(7,1) - t52 * rSges(7,2) + pkin(5);
t41 = t101 * t70 + t51 * t64;
t31 = t47 * t92 - t54 * t64;
t29 = t101 * t66 + t51 * t60;
t28 = t101 * t65 + t51 * t59;
t27 = t47 * t109 - t119 * t54;
t26 = t109 * t119 + t47 * t54;
t25 = t26 * pkin(4);
t24 = t107 * t47 + t32 * t56;
t23 = -t106 * t47 + t32 * t53;
t21 = t40 * t92 - t54 * t60;
t19 = t39 * t92 - t54 * t59;
t16 = t27 * t56 + t41 * t53;
t15 = -t27 * t53 + t41 * t56;
t14 = t40 * t109 - t120 * t54;
t13 = t109 * t120 + t40 * t54;
t12 = t39 * t109 - t121 * t54;
t11 = t109 * t121 + t39 * t54;
t10 = t13 * pkin(4);
t9 = t11 * pkin(4);
t8 = t107 * t40 + t22 * t56;
t7 = -t106 * t40 + t22 * t53;
t6 = t107 * t39 + t20 * t56;
t5 = -t106 * t39 + t20 * t53;
t4 = t14 * t56 + t29 * t53;
t3 = -t14 * t53 + t29 * t56;
t2 = t12 * t56 + t28 * t53;
t1 = -t12 * t53 + t28 * t56;
t17 = [(-m(2) + t93) * g(3), t93 * (g(1) * t82 - g(2) * t85 + g(3) * t103) -m(4) * (g(1) * (-rSges(4,1) * t60 - t40 * rSges(4,2)) + g(2) * (-rSges(4,1) * t59 - t39 * rSges(4,2)) + g(3) * (-rSges(4,1) * t64 - t47 * rSges(4,2))) - m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t112 * t21 + t89) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + t112 * t19 + t90) + g(3) * (rSges(6,1) * t24 - rSges(6,2) * t23 + t112 * t31 + t88)) - m(7) * (g(1) * (t8 * pkin(5) + t21 * pkin(11) + (t21 * t52 + t55 * t8) * rSges(7,1) + (t21 * t55 - t52 * t8) * rSges(7,2) + t111 * t7 + t89) + g(2) * (t6 * pkin(5) + t19 * pkin(11) + (t19 * t52 + t55 * t6) * rSges(7,1) + (t19 * t55 - t52 * t6) * rSges(7,2) + t111 * t5 + t90) + g(3) * (t24 * pkin(5) + t31 * pkin(11) + (t24 * t55 + t31 * t52) * rSges(7,1) + (-t24 * t52 + t31 * t55) * rSges(7,2) + t111 * t23 + t88)) + (-g(1) * (rSges(5,1) * t22 - rSges(5,2) * t21 - t38) - g(2) * (rSges(5,1) * t20 - rSges(5,2) * t19 - t37) - g(3) * (rSges(5,1) * t32 - rSges(5,2) * t31 - t46) - (g(1) * t40 + g(2) * t39 + g(3) * t47) * t51 * (rSges(5,3) + pkin(10))) * m(5), -m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t14) + g(2) * (-rSges(5,1) * t11 - rSges(5,2) * t12) + g(3) * (-rSges(5,1) * t26 - rSges(5,2) * t27)) - m(6) * (g(1) * (t112 * t14 + t13 * t91 - t10) + g(2) * (t11 * t91 + t112 * t12 - t9) + g(3) * (t112 * t27 + t26 * t91 - t25)) + (-g(1) * (-t13 * t114 - t10 + t14 * pkin(11) + (-t104 * t13 + t14 * t52) * rSges(7,1) + (t105 * t13 + t14 * t55) * rSges(7,2)) - g(2) * (-t11 * t114 - t9 + t12 * pkin(11) + (-t104 * t11 + t12 * t52) * rSges(7,1) + (t105 * t11 + t12 * t55) * rSges(7,2)) - g(3) * (-t26 * t114 - t25 + t27 * pkin(11) + (-t104 * t26 + t27 * t52) * rSges(7,1) + (t105 * t26 + t27 * t55) * rSges(7,2)) - (-g(1) * t13 - g(2) * t11 - g(3) * t26) * t53 * t111) * m(7), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (rSges(6,1) * t15 - rSges(6,2) * t16)) - m(7) * (g(1) * (t111 * t4 + t3 * t78) + (t111 * t16 + t78 * t15) * g(3) + (t78 * t1 + t111 * t2) * g(2)) -m(7) * (g(1) * ((t13 * t55 - t4 * t52) * rSges(7,1) + (-t13 * t52 - t4 * t55) * rSges(7,2)) + g(2) * ((t11 * t55 - t2 * t52) * rSges(7,1) + (-t11 * t52 - t2 * t55) * rSges(7,2)) + g(3) * ((-t16 * t52 + t26 * t55) * rSges(7,1) + (-t16 * t55 - t26 * t52) * rSges(7,2)))];
taug  = t17(:);
