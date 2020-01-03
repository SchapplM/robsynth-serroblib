% Calculate Gravitation load on the joints for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:00
% EndTime: 2019-12-31 22:32:03
% DurationCPUTime: 0.90s
% Computational Cost: add. (564->148), mult. (982->217), div. (0->0), fcn. (1143->12), ass. (0->69)
t108 = pkin(10) + rSges(6,3);
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t90 = cos(pkin(5));
t58 = sin(pkin(5));
t61 = sin(qJ(2));
t98 = t58 * t61;
t119 = -t60 * t98 + t90 * t64;
t106 = cos(qJ(1));
t65 = cos(qJ(2));
t62 = sin(qJ(1));
t82 = t62 * t90;
t39 = t106 * t65 - t61 * t82;
t96 = t58 * t64;
t19 = -t39 * t60 + t62 * t96;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t118 = rSges(6,1) * t63 - rSges(6,2) * t59 + pkin(4);
t78 = t90 * t106;
t37 = t61 * t78 + t62 * t65;
t57 = qJ(3) + qJ(4);
t54 = sin(t57);
t55 = cos(t57);
t86 = t58 * t106;
t14 = t37 * t55 - t54 * t86;
t36 = t62 * t61 - t65 * t78;
t117 = t14 * t59 - t36 * t63;
t116 = -t14 * t63 - t36 * t59;
t115 = rSges(6,1) * t59 + rSges(6,2) * t63;
t114 = t108 * t54 + t118 * t55;
t84 = -t37 * t54 - t55 * t86;
t113 = rSges(5,1) * t84 - t14 * rSges(5,2);
t112 = g(2) * t36;
t53 = pkin(3) * t64 + pkin(2);
t95 = t58 * t65;
t111 = g(3) * t53 * t95;
t110 = g(3) * t58;
t109 = -pkin(8) - rSges(4,3);
t97 = t58 * t62;
t17 = t39 * t54 - t55 * t97;
t18 = t39 * t55 + t54 * t97;
t107 = -t17 * rSges(5,1) - t18 * rSges(5,2);
t66 = -pkin(9) - pkin(8);
t94 = -t36 * t53 - t37 * t66;
t38 = t106 * t61 + t65 * t82;
t93 = -t38 * t53 - t39 * t66;
t30 = -t54 * t98 + t90 * t55;
t31 = t90 * t54 + t55 * t98;
t92 = t30 * rSges(5,1) - t31 * rSges(5,2);
t91 = t106 * pkin(1) + pkin(7) * t97;
t89 = t60 * t97;
t85 = -t62 * pkin(1) + pkin(7) * t86;
t48 = t60 * t86;
t83 = -t37 * t64 + t48;
t80 = t19 * pkin(3);
t79 = pkin(3) * t89 - t38 * t66 + t39 * t53 + t91;
t77 = rSges(5,1) * t55 - rSges(5,2) * t54;
t76 = t119 * pkin(3);
t75 = rSges(4,1) * t64 - rSges(4,2) * t60 + pkin(2);
t73 = pkin(3) * t48 + t36 * t66 - t37 * t53 + t85;
t72 = t37 * t60 + t64 * t86;
t71 = t108 * t14 + t118 * t84;
t70 = t108 * t18 - t118 * t17;
t69 = t108 * t31 + t118 * t30;
t68 = t72 * pkin(3);
t20 = t39 * t64 + t89;
t2 = t18 * t63 + t38 * t59;
t1 = -t18 * t59 + t38 * t63;
t3 = [-m(2) * (g(1) * (-t62 * rSges(2,1) - t106 * rSges(2,2)) + g(2) * (t106 * rSges(2,1) - t62 * rSges(2,2))) - m(3) * (g(1) * (-t37 * rSges(3,1) + t36 * rSges(3,2) + rSges(3,3) * t86 + t85) + g(2) * (t39 * rSges(3,1) - t38 * rSges(3,2) + rSges(3,3) * t97 + t91)) - m(4) * (g(1) * (t83 * rSges(4,1) + t72 * rSges(4,2) - t37 * pkin(2) + t109 * t36 + t85) + g(2) * (t20 * rSges(4,1) + t19 * rSges(4,2) + t39 * pkin(2) - t109 * t38 + t91)) - m(5) * (g(1) * (-rSges(5,1) * t14 - rSges(5,2) * t84 - t36 * rSges(5,3) + t73) + g(2) * (t18 * rSges(5,1) - t17 * rSges(5,2) + t38 * rSges(5,3) + t79)) - m(6) * (g(1) * (t116 * rSges(6,1) + t117 * rSges(6,2) - t14 * pkin(4) + t108 * t84 + t73) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t18 * pkin(4) + t108 * t17 + t79)), -m(3) * (g(1) * (-rSges(3,1) * t38 - rSges(3,2) * t39) + g(2) * (-rSges(3,1) * t36 - rSges(3,2) * t37) + (rSges(3,1) * t65 - rSges(3,2) * t61) * t110) - m(4) * (g(1) * (-t109 * t39 - t75 * t38) - g(2) * t109 * t37 - t75 * t112 + (-t109 * t61 + t75 * t65) * t110) - m(5) * (g(1) * (t39 * rSges(5,3) - t77 * t38 + t93) + g(2) * (t37 * rSges(5,3) - t77 * t36 + t94) + t111 + (t77 * t65 + (rSges(5,3) - t66) * t61) * t110) - m(6) * (g(2) * (t115 * t37 + t94) + t111 - t114 * t112 + ((-t66 + t115) * t61 + t114 * t65) * t110 + (-t114 * t38 + t115 * t39 + t93) * g(1)), -m(4) * (g(1) * (rSges(4,1) * t19 - rSges(4,2) * t20) + g(2) * (-t72 * rSges(4,1) + t83 * rSges(4,2)) + g(3) * (t119 * rSges(4,1) + (-t90 * t60 - t61 * t96) * rSges(4,2))) - m(5) * (g(1) * (t80 + t107) + g(2) * (-t68 + t113) + g(3) * (t76 + t92)) - m(6) * (g(1) * (t70 + t80) + g(2) * (-t68 + t71) + g(3) * (t69 + t76)), -m(5) * (g(1) * t107 + g(2) * t113 + g(3) * t92) - m(6) * (g(1) * t70 + g(2) * t71 + g(3) * t69), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (-t117 * rSges(6,1) + t116 * rSges(6,2)) + g(3) * ((-t31 * t59 - t63 * t95) * rSges(6,1) + (-t31 * t63 + t59 * t95) * rSges(6,2)))];
taug = t3(:);
