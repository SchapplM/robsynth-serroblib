% Calculate Gravitation load on the joints for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:29:02
% EndTime: 2019-03-10 04:29:06
% DurationCPUTime: 1.44s
% Computational Cost: add. (968->196), mult. (1790->277), div. (0->0), fcn. (2146->14), ass. (0->85)
t104 = sin(qJ(1));
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t105 = cos(qJ(1));
t95 = cos(pkin(6));
t88 = t95 * t105;
t32 = t104 * t64 + t61 * t88;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t58 = sin(pkin(6));
t92 = t58 * t105;
t20 = t32 * t63 - t60 * t92;
t31 = t104 * t61 - t64 * t88;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t128 = t20 * t59 - t31 * t62;
t127 = -t20 * t62 - t31 * t59;
t87 = t95 * t104;
t34 = t105 * t64 - t61 * t87;
t126 = g(1) * t34 + g(2) * t32;
t33 = t105 * t61 + t64 * t87;
t125 = g(1) * t33 + g(2) * t31;
t91 = t58 * t104;
t24 = t34 * t63 + t60 * t91;
t101 = t58 * t61;
t30 = t63 * t101 + t95 * t60;
t124 = g(1) * t24 + g(2) * t20 + g(3) * t30;
t23 = t34 * t60 - t63 * t91;
t29 = -t60 * t101 + t95 * t63;
t90 = -t32 * t60 - t63 * t92;
t123 = -g(1) * t23 + g(2) * t90 + g(3) * t29;
t57 = qJ(4) + qJ(5);
t52 = cos(t57);
t54 = t62 * pkin(4);
t37 = pkin(5) * t52 + t54;
t35 = pkin(3) + t37;
t53 = qJ(6) + t57;
t48 = sin(t53);
t49 = cos(t53);
t76 = rSges(7,1) * t49 - rSges(7,2) * t48 + t35;
t65 = -pkin(11) - pkin(10);
t97 = pkin(12) - t65 + rSges(7,3);
t122 = t97 * t60 + t76 * t63;
t50 = t54 + pkin(3);
t51 = sin(t57);
t77 = rSges(6,1) * t52 - rSges(6,2) * t51 + t50;
t96 = -t65 + rSges(6,3);
t121 = t96 * t60 + t77 * t63;
t106 = pkin(10) + rSges(5,3);
t81 = rSges(5,1) * t62 - rSges(5,2) * t59 + pkin(3);
t120 = t106 * t60 + t81 * t63;
t119 = (-t20 * t48 + t31 * t49) * rSges(7,1) + (-t20 * t49 - t31 * t48) * rSges(7,2);
t5 = -t24 * t48 + t33 * t49;
t6 = t24 * t49 + t33 * t48;
t118 = t5 * rSges(7,1) - t6 * rSges(7,2);
t109 = g(3) * t58;
t108 = t59 * pkin(4);
t107 = rSges(4,3) + pkin(9);
t100 = t58 * t64;
t99 = (-t49 * t100 - t30 * t48) * rSges(7,1) + (t48 * t100 - t30 * t49) * rSges(7,2);
t98 = t105 * pkin(1) + pkin(8) * t91;
t94 = t34 * pkin(2) + t98;
t93 = g(3) * (pkin(2) * t100 + pkin(9) * t101);
t89 = -t104 * pkin(1) + pkin(8) * t92;
t86 = rSges(4,1) * t63 - rSges(4,2) * t60;
t85 = t59 * rSges(5,1) + t62 * rSges(5,2);
t84 = -t20 * t51 + t31 * t52;
t11 = -t24 * t51 + t33 * t52;
t13 = -t24 * t59 + t33 * t62;
t82 = -t32 * pkin(2) + t89;
t79 = -t52 * t100 - t30 * t51;
t78 = -t62 * t100 - t30 * t59;
t36 = pkin(5) * t51 + t108;
t75 = t48 * rSges(7,1) + t49 * rSges(7,2) + t36;
t74 = pkin(9) + t75;
t73 = t51 * rSges(6,1) + t52 * rSges(6,2) + t108;
t72 = pkin(9) + t73;
t25 = t31 * pkin(2);
t27 = t33 * pkin(2);
t71 = -g(1) * t27 - g(2) * t25 + t93;
t67 = m(7) * (g(1) * t118 + g(2) * t119 + g(3) * t99);
t12 = t24 * t52 + t33 * t51;
t66 = m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (t84 * rSges(6,1) + (-t20 * t52 - t31 * t51) * rSges(6,2)) + g(3) * (t79 * rSges(6,1) + (t51 * t100 - t30 * t52) * rSges(6,2)));
t14 = t24 * t62 + t33 * t59;
t1 = [-m(2) * (g(1) * (-t104 * rSges(2,1) - t105 * rSges(2,2)) + g(2) * (t105 * rSges(2,1) - t104 * rSges(2,2))) - m(3) * (g(1) * (-t32 * rSges(3,1) + t31 * rSges(3,2) + rSges(3,3) * t92 + t89) + g(2) * (t34 * rSges(3,1) - t33 * rSges(3,2) + rSges(3,3) * t91 + t98)) - m(4) * (g(1) * (-rSges(4,1) * t20 - rSges(4,2) * t90 - t107 * t31 + t82) + g(2) * (rSges(4,1) * t24 - rSges(4,2) * t23 + t107 * t33 + t94)) - m(5) * (g(1) * (t127 * rSges(5,1) + t128 * rSges(5,2) - t20 * pkin(3) - t31 * pkin(9) + t106 * t90 + t82) + g(2) * (rSges(5,1) * t14 + rSges(5,2) * t13 + pkin(3) * t24 + t33 * pkin(9) + t106 * t23 + t94)) - m(6) * (g(1) * (-t20 * t77 - t72 * t31 + t90 * t96 + t82) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t24 * t50 + (pkin(9) + t108) * t33 + t96 * t23 + t94)) - m(7) * (g(1) * (-t20 * t76 - t74 * t31 + t90 * t97 + t82) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t24 * t35 + (pkin(9) + t36) * t33 + t97 * t23 + t94)) -m(3) * (g(1) * (-rSges(3,1) * t33 - rSges(3,2) * t34) + g(2) * (-rSges(3,1) * t31 - rSges(3,2) * t32) + (rSges(3,1) * t64 - rSges(3,2) * t61) * t109) - m(4) * (g(1) * (t107 * t34 - t86 * t33 - t27) + g(2) * (t107 * t32 - t86 * t31 - t25) + t93 + (rSges(4,3) * t61 + t86 * t64) * t109) - m(5) * ((t120 * t64 + t85 * t61) * t109 + t71 + t126 * (pkin(9) + t85) - t125 * t120) - m(6) * ((t121 * t64 + t73 * t61) * t109 + t71 + t126 * t72 - t125 * t121) - m(7) * ((t122 * t64 + t75 * t61) * t109 + t71 + t126 * t74 - t125 * t122) -m(4) * (g(1) * (-rSges(4,1) * t23 - rSges(4,2) * t24) + g(2) * (rSges(4,1) * t90 - rSges(4,2) * t20) + g(3) * (rSges(4,1) * t29 - rSges(4,2) * t30)) - m(5) * (t124 * t106 + t123 * t81) - m(6) * (t123 * t77 + t124 * t96) - m(7) * (t123 * t76 + t124 * t97) -m(5) * (g(1) * (rSges(5,1) * t13 - rSges(5,2) * t14) + g(2) * (-rSges(5,1) * t128 + t127 * rSges(5,2)) + g(3) * (t78 * rSges(5,1) + (t59 * t100 - t30 * t62) * rSges(5,2))) - t66 - m(7) * (g(1) * (-t24 * t36 + t33 * t37 + t118) + g(2) * (-t20 * t36 + t31 * t37 + t119) + g(3) * (-t37 * t100 - t30 * t36 + t99)) - m(6) * (g(1) * t13 - g(2) * t128 + g(3) * t78) * pkin(4), -t66 - t67 - m(7) * (g(1) * t11 + g(2) * t84 + g(3) * t79) * pkin(5), -t67];
taug  = t1(:);
