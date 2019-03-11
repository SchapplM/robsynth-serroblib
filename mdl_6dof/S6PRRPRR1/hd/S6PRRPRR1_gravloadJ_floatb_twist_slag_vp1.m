% Calculate Gravitation load on the joints for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:10
% EndTime: 2019-03-08 21:51:12
% DurationCPUTime: 0.85s
% Computational Cost: add. (652->138), mult. (874->206), div. (0->0), fcn. (997->14), ass. (0->70)
t108 = rSges(7,3) + pkin(10);
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t107 = rSges(7,1) * t59 - rSges(7,2) * t56 + pkin(5);
t54 = sin(pkin(6));
t61 = cos(qJ(2));
t88 = t54 * t61;
t106 = g(3) * t88;
t105 = t56 * rSges(7,1) + t59 * rSges(7,2);
t58 = sin(qJ(2));
t53 = sin(pkin(11));
t79 = cos(pkin(6));
t76 = t53 * t79;
t78 = cos(pkin(11));
t34 = t58 * t78 + t61 * t76;
t73 = t79 * t78;
t32 = t53 * t58 - t61 * t73;
t99 = g(2) * t32;
t104 = g(1) * t34 + t99;
t52 = qJ(3) + pkin(12);
t49 = qJ(5) + t52;
t44 = sin(t49);
t45 = cos(t49);
t103 = t45 * t107 + t108 * t44;
t33 = t53 * t61 + t58 * t73;
t75 = t54 * t78;
t12 = -t33 * t44 - t45 * t75;
t13 = t33 * t45 - t44 * t75;
t102 = t12 * rSges(6,1) - t13 * rSges(6,2);
t35 = -t58 * t76 + t61 * t78;
t91 = t53 * t54;
t14 = -t35 * t44 + t45 * t91;
t15 = t35 * t45 + t44 * t91;
t101 = t14 * rSges(6,1) - t15 * rSges(6,2);
t98 = g(2) * t33;
t48 = cos(t52);
t60 = cos(qJ(3));
t50 = t60 * pkin(3);
t40 = pkin(4) * t48 + t50;
t38 = pkin(2) + t40;
t97 = t38 * t106;
t96 = g(3) * t54;
t95 = rSges(4,3) + pkin(8);
t90 = t54 * t58;
t89 = t54 * t60;
t55 = -qJ(4) - pkin(8);
t85 = rSges(5,3) - t55;
t51 = -pkin(9) + t55;
t84 = -t32 * t38 - t33 * t51;
t83 = -t34 * t38 - t35 * t51;
t47 = sin(t52);
t57 = sin(qJ(3));
t39 = -pkin(3) * t57 - pkin(4) * t47;
t82 = t35 * t39 + t40 * t91;
t25 = -t44 * t90 + t45 * t79;
t26 = t44 * t79 + t45 * t90;
t81 = t25 * rSges(6,1) - t26 * rSges(6,2);
t80 = t39 * t90 + t79 * t40;
t77 = -m(5) - m(6) - m(7);
t74 = rSges(6,1) * t45 - rSges(6,2) * t44;
t72 = rSges(4,1) * t60 - rSges(4,2) * t57 + pkin(2);
t70 = -t35 * t57 + t53 * t89;
t69 = rSges(5,1) * t48 - rSges(5,2) * t47 + pkin(2) + t50;
t68 = t33 * t39 - t40 * t75;
t67 = -t33 * t57 - t60 * t75;
t66 = -t57 * t90 + t60 * t79;
t65 = t107 * t12 + t108 * t13;
t64 = t107 * t14 + t108 * t15;
t63 = t107 * t25 + t108 * t26;
t1 = [(-m(2) - m(3) - m(4) + t77) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + g(2) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + (rSges(3,1) * t61 - rSges(3,2) * t58) * t96) - m(4) * (g(1) * (-t34 * t72 + t35 * t95) + t95 * t98 - t72 * t99 + (t58 * t95 + t61 * t72) * t96) - m(5) * (g(1) * (-t34 * t69 + t35 * t85) + t85 * t98 - t69 * t99 + (t58 * t85 + t61 * t69) * t96) - m(6) * (g(1) * (rSges(6,3) * t35 - t34 * t74 + t83) + g(2) * (rSges(6,3) * t33 - t74 * t32 + t84) + t97 + (t74 * t61 + (rSges(6,3) - t51) * t58) * t96) - m(7) * (g(1) * (t105 * t35 + t83) + g(2) * (t105 * t33 + t84) + t97 + ((-t51 + t105) * t58 + t103 * t61) * t96 - t104 * t103) -m(4) * (g(1) * (t70 * rSges(4,1) + (-t35 * t60 - t57 * t91) * rSges(4,2)) + g(2) * (t67 * rSges(4,1) + (-t33 * t60 + t57 * t75) * rSges(4,2)) + g(3) * (t66 * rSges(4,1) + (-t57 * t79 - t58 * t89) * rSges(4,2))) - m(5) * (g(1) * ((-t35 * t47 + t48 * t91) * rSges(5,1) + (-t35 * t48 - t47 * t91) * rSges(5,2)) + g(2) * ((-t33 * t47 - t48 * t75) * rSges(5,1) + (-t33 * t48 + t47 * t75) * rSges(5,2)) + g(3) * ((-t47 * t90 + t48 * t79) * rSges(5,1) + (-t47 * t79 - t48 * t90) * rSges(5,2)) + (g(1) * t70 + g(2) * t67 + g(3) * t66) * pkin(3)) - m(6) * (g(1) * (t82 + t101) + g(2) * (t68 + t102) + g(3) * (t80 + t81)) - m(7) * (g(1) * (t64 + t82) + g(2) * (t65 + t68) + g(3) * (t63 + t80)) t77 * (t104 - t106) -m(6) * (g(1) * t101 + g(2) * t102 + g(3) * t81) - m(7) * (g(1) * t64 + g(2) * t65 + g(3) * t63) -m(7) * (g(1) * ((-t15 * t56 + t34 * t59) * rSges(7,1) + (-t15 * t59 - t34 * t56) * rSges(7,2)) + g(2) * ((-t13 * t56 + t32 * t59) * rSges(7,1) + (-t13 * t59 - t32 * t56) * rSges(7,2)) + g(3) * ((-t26 * t56 - t59 * t88) * rSges(7,1) + (-t26 * t59 + t56 * t88) * rSges(7,2)))];
taug  = t1(:);
