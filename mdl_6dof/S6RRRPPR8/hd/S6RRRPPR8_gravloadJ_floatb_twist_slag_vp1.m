% Calculate Gravitation load on the joints for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:37
% EndTime: 2019-03-09 16:03:39
% DurationCPUTime: 1.06s
% Computational Cost: add. (595->194), mult. (1414->267), div. (0->0), fcn. (1677->10), ass. (0->78)
t104 = pkin(10) + rSges(7,3);
t90 = rSges(6,1) + qJ(4);
t109 = -m(6) - m(7);
t100 = sin(qJ(1));
t101 = cos(qJ(2));
t53 = sin(qJ(2));
t102 = cos(qJ(1));
t87 = cos(pkin(6));
t66 = t87 * t102;
t34 = t100 * t101 + t53 * t66;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t50 = sin(pkin(6));
t82 = t50 * t102;
t13 = t34 * t52 + t55 * t82;
t6 = t13 * pkin(3);
t108 = -t13 * pkin(4) - t6;
t107 = g(3) * t50;
t106 = rSges(5,2) + pkin(9);
t105 = rSges(4,3) + pkin(9);
t65 = t87 * t100;
t36 = t102 * t101 - t53 * t65;
t80 = t50 * t100;
t17 = t36 * t52 - t55 * t80;
t10 = t17 * pkin(3);
t103 = -t17 * pkin(4) - t10;
t99 = rSges(6,1) * t52;
t33 = t100 * t53 - t101 * t66;
t98 = t33 * t55;
t35 = t101 * t65 + t102 * t53;
t97 = t35 * t55;
t96 = t50 * t53;
t95 = qJ(5) - pkin(9);
t31 = t52 * t96 - t87 * t55;
t26 = t31 * pkin(3);
t94 = -t31 * pkin(4) - t26;
t93 = t102 * pkin(1) + pkin(8) * t80;
t81 = t50 * t101;
t92 = pkin(2) * t81 + pkin(9) * t96;
t91 = qJ(4) * t52;
t89 = rSges(5,3) + qJ(4);
t88 = -rSges(6,3) - qJ(5);
t86 = -pkin(9) - t88;
t27 = t33 * pkin(2);
t85 = -pkin(3) * t98 - t33 * t91 - t27;
t29 = t35 * pkin(2);
t84 = -pkin(3) * t97 - t35 * t91 - t29;
t83 = t36 * pkin(2) + t93;
t79 = t52 * t101;
t54 = cos(qJ(6));
t78 = t54 * t101;
t77 = t55 * t101;
t14 = t34 * t55 - t52 * t82;
t18 = t36 * t55 + t52 * t80;
t76 = t18 * pkin(3) + t83;
t75 = -pkin(4) * t98 + t85;
t74 = -pkin(4) * t97 + t84;
t70 = t50 * t77;
t71 = t50 * t79;
t73 = pkin(3) * t70 + qJ(4) * t71 + t92;
t72 = g(2) * t86;
t69 = -t100 * pkin(1) + pkin(8) * t82;
t68 = t18 * pkin(4) + t76;
t67 = pkin(4) * t70 + t73;
t64 = -rSges(4,1) * t55 + rSges(4,2) * t52;
t63 = -rSges(5,1) * t55 - rSges(5,3) * t52;
t62 = -t34 * pkin(2) + t69;
t51 = sin(qJ(6));
t61 = -t54 * rSges(7,1) + t51 * rSges(7,2) - pkin(5);
t60 = -pkin(3) * t14 + t62;
t59 = qJ(4) - t61;
t58 = t51 * rSges(7,1) + t54 * rSges(7,2) + t95;
t57 = -pkin(4) * t14 + t60;
t56 = -t104 * t55 + t61 * t52;
t32 = t87 * t52 + t55 * t96;
t3 = t17 * t54 - t35 * t51;
t2 = -t17 * t51 - t35 * t54;
t1 = [-m(2) * (g(1) * (-t100 * rSges(2,1) - t102 * rSges(2,2)) + g(2) * (t102 * rSges(2,1) - t100 * rSges(2,2))) - m(3) * (g(1) * (-t34 * rSges(3,1) + t33 * rSges(3,2) + rSges(3,3) * t82 + t69) + g(2) * (t36 * rSges(3,1) - t35 * rSges(3,2) + rSges(3,3) * t80 + t93)) - m(4) * (g(1) * (-rSges(4,1) * t14 + rSges(4,2) * t13 - t105 * t33 + t62) + g(2) * (rSges(4,1) * t18 - rSges(4,2) * t17 + t105 * t35 + t83)) - m(5) * (g(1) * (-rSges(5,1) * t14 - t106 * t33 - t13 * t89 + t60) + g(2) * (rSges(5,1) * t18 + t106 * t35 + t89 * t17 + t76)) - m(6) * (g(2) * (-rSges(6,2) * t18 + t90 * t17 + t68) - t35 * t72 + (rSges(6,2) * t14 - t13 * t90 + t86 * t33 + t57) * g(1)) - m(7) * (g(1) * (-t104 * t14 - t13 * t59 + t58 * t33 + t57) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 - t95 * t35 + t104 * t18 + (pkin(5) + qJ(4)) * t17 + t68)) -m(3) * (g(1) * (-rSges(3,1) * t35 - rSges(3,2) * t36) + g(2) * (-rSges(3,1) * t33 - rSges(3,2) * t34) + (t101 * rSges(3,1) - rSges(3,2) * t53) * t107) - m(4) * (g(1) * (t105 * t36 + t64 * t35 - t29) + g(2) * (t105 * t34 + t64 * t33 - t27) + g(3) * t92 + (rSges(4,1) * t77 - rSges(4,2) * t79 + rSges(4,3) * t53) * t107) - m(5) * (g(1) * (t106 * t36 + t63 * t35 + t84) + g(2) * (t106 * t34 + t63 * t33 + t85) + g(3) * t73 + (rSges(5,1) * t77 + rSges(5,2) * t53 + rSges(5,3) * t79) * t107) - m(6) * (g(2) * (rSges(6,2) * t98 - t33 * t99 + t75) + g(3) * t67 - t34 * t72 + (rSges(6,1) * t79 - rSges(6,2) * t77 + t88 * t53) * t107 + (rSges(6,2) * t97 - t35 * t99 - t86 * t36 + t74) * g(1)) - m(7) * (g(1) * (t56 * t35 - t58 * t36 + t74) + g(2) * (t56 * t33 - t58 * t34 + t75) + g(3) * (pkin(5) * t71 - qJ(5) * t96 + t67 + t104 * t70 + ((-t51 * t53 + t52 * t78) * rSges(7,1) + (-t51 * t79 - t53 * t54) * rSges(7,2)) * t50)) -m(4) * (g(1) * (-rSges(4,1) * t17 - rSges(4,2) * t18) + g(2) * (-rSges(4,1) * t13 - rSges(4,2) * t14) + g(3) * (-rSges(4,1) * t31 - rSges(4,2) * t32)) - m(5) * (g(1) * (-rSges(5,1) * t17 + t89 * t18 - t10) + g(2) * (-rSges(5,1) * t13 + t89 * t14 - t6) + g(3) * (-rSges(5,1) * t31 + t89 * t32 - t26)) - m(6) * (g(1) * (rSges(6,2) * t17 + t90 * t18 + t103) + g(2) * (rSges(6,2) * t13 + t90 * t14 + t108) + g(3) * (rSges(6,2) * t31 + t90 * t32 + t94)) + (-g(1) * (-t104 * t17 + t103) - g(2) * (-t104 * t13 + t108) - g(3) * (-t104 * t31 + t94) - (g(1) * t18 + g(2) * t14 + g(3) * t32) * t59) * m(7) (-m(5) + t109) * (g(1) * t17 + g(2) * t13 + g(3) * t31) t109 * (-g(1) * t35 - g(2) * t33 + g(3) * t81) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((-t13 * t51 - t33 * t54) * rSges(7,1) + (-t13 * t54 + t33 * t51) * rSges(7,2)) + g(3) * ((-t31 * t51 + t50 * t78) * rSges(7,1) + (-t31 * t54 - t51 * t81) * rSges(7,2)))];
taug  = t1(:);
