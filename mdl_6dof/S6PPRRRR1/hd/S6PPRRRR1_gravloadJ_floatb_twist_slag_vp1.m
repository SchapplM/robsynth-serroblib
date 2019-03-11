% Calculate Gravitation load on the joints for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:44
% EndTime: 2019-03-08 18:59:45
% DurationCPUTime: 0.67s
% Computational Cost: add. (864->124), mult. (2021->201), div. (0->0), fcn. (2566->16), ass. (0->73)
t122 = rSges(7,3) + pkin(11);
t106 = cos(qJ(3));
t86 = sin(pkin(7));
t87 = sin(pkin(6));
t88 = cos(pkin(13));
t90 = cos(pkin(7));
t91 = cos(pkin(6));
t115 = t88 * t90 * t87 + t91 * t86;
t55 = sin(qJ(3));
t84 = sin(pkin(13));
t72 = t87 * t84;
t35 = t106 * t72 + t115 * t55;
t74 = t87 * t86;
t42 = -t88 * t74 + t91 * t90;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t121 = -t35 * t54 + t42 * t57;
t85 = sin(pkin(12));
t71 = t85 * t88;
t89 = cos(pkin(12));
t75 = t89 * t84;
t62 = t91 * t71 + t75;
t73 = t87 * t85;
t116 = t62 * t90 - t86 * t73;
t70 = t85 * t84;
t77 = t89 * t88;
t44 = -t91 * t70 + t77;
t29 = t44 * t106 - t116 * t55;
t37 = t62 * t86 + t90 * t73;
t120 = -t29 * t54 + t37 * t57;
t61 = -t91 * t77 + t70;
t117 = t61 * t90 + t89 * t74;
t43 = t91 * t75 + t71;
t27 = t43 * t106 - t117 * t55;
t76 = t89 * t87;
t36 = t61 * t86 - t90 * t76;
t119 = -t27 * t54 + t36 * t57;
t53 = sin(qJ(6));
t56 = cos(qJ(6));
t118 = t56 * rSges(7,1) - t53 * rSges(7,2) + pkin(5);
t26 = t117 * t106 + t43 * t55;
t28 = t116 * t106 + t44 * t55;
t34 = -t115 * t106 + t55 * t72;
t114 = g(1) * t28 + g(2) * t26 + g(3) * t34;
t52 = qJ(4) + qJ(5);
t50 = sin(t52);
t51 = cos(t52);
t11 = -t27 * t50 + t36 * t51;
t12 = t27 * t51 + t36 * t50;
t113 = t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = -t29 * t50 + t37 * t51;
t14 = t29 * t51 + t37 * t50;
t112 = t13 * rSges(6,1) - t14 * rSges(6,2);
t109 = t51 * pkin(5);
t99 = t51 * t53;
t98 = t51 * t56;
t49 = t57 * pkin(4) + pkin(3);
t58 = -pkin(10) - pkin(9);
t95 = -t26 * t49 - t27 * t58;
t94 = -t28 * t49 - t29 * t58;
t22 = -t35 * t50 + t42 * t51;
t23 = t35 * t51 + t42 * t50;
t93 = t22 * rSges(6,1) - t23 * rSges(6,2);
t92 = -t34 * t49 - t35 * t58;
t83 = -m(3) - m(4) - m(5) - m(6) - m(7);
t82 = t119 * pkin(4);
t81 = t120 * pkin(4);
t80 = t121 * pkin(4);
t79 = -rSges(6,1) * t51 + rSges(6,2) * t50;
t65 = t118 * t11 + t122 * t12;
t64 = t118 * t13 + t122 * t14;
t63 = t118 * t22 + t122 * t23;
t1 = [(-m(2) + t83) * g(3), t83 * (g(1) * t73 - g(2) * t76 + g(3) * t91) -m(4) * (g(1) * (-t28 * rSges(4,1) - t29 * rSges(4,2)) + g(2) * (-t26 * rSges(4,1) - t27 * rSges(4,2)) + g(3) * (-t34 * rSges(4,1) - t35 * rSges(4,2))) - m(5) * (t114 * (-t57 * rSges(5,1) + t54 * rSges(5,2) - pkin(3)) + (g(1) * t29 + g(2) * t27 + g(3) * t35) * (rSges(5,3) + pkin(9))) - m(6) * (g(1) * (t29 * rSges(6,3) + t28 * t79 + t94) + g(2) * (t27 * rSges(6,3) + t26 * t79 + t95) + g(3) * (t35 * rSges(6,3) + t34 * t79 + t92)) + (-g(1) * (-t28 * t109 + (-t28 * t98 + t29 * t53) * rSges(7,1) + (t28 * t99 + t29 * t56) * rSges(7,2) + t94) - g(2) * (-t26 * t109 + (-t26 * t98 + t27 * t53) * rSges(7,1) + (t26 * t99 + t27 * t56) * rSges(7,2) + t95) - g(3) * (-t34 * t109 + (-t34 * t98 + t35 * t53) * rSges(7,1) + (t34 * t99 + t35 * t56) * rSges(7,2) + t92) + t114 * t50 * t122) * m(7), -m(5) * (g(1) * (t120 * rSges(5,1) + (-t29 * t57 - t37 * t54) * rSges(5,2)) + g(2) * (t119 * rSges(5,1) + (-t27 * t57 - t36 * t54) * rSges(5,2)) + g(3) * (t121 * rSges(5,1) + (-t35 * t57 - t42 * t54) * rSges(5,2))) - m(6) * (g(1) * (t81 + t112) + g(2) * (t82 + t113) + g(3) * (t80 + t93)) - m(7) * (g(1) * (t64 + t81) + g(2) * (t65 + t82) + g(3) * (t63 + t80)) -m(6) * (g(1) * t112 + g(2) * t113 + g(3) * t93) - m(7) * (g(1) * t64 + g(2) * t65 + g(3) * t63) -m(7) * (g(1) * ((-t14 * t53 + t28 * t56) * rSges(7,1) + (-t14 * t56 - t28 * t53) * rSges(7,2)) + g(2) * ((-t12 * t53 + t26 * t56) * rSges(7,1) + (-t12 * t56 - t26 * t53) * rSges(7,2)) + g(3) * ((-t23 * t53 + t34 * t56) * rSges(7,1) + (-t23 * t56 - t34 * t53) * rSges(7,2)))];
taug  = t1(:);
