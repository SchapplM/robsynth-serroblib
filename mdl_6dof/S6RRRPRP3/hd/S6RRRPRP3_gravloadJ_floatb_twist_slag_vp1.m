% Calculate Gravitation load on the joints for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:08
% EndTime: 2019-03-09 16:39:11
% DurationCPUTime: 1.11s
% Computational Cost: add. (586->148), mult. (595->189), div. (0->0), fcn. (559->10), ass. (0->72)
t46 = qJ(2) + qJ(3);
t42 = sin(t46);
t43 = cos(t46);
t48 = cos(pkin(10));
t77 = -rSges(5,1) * t48 - pkin(3);
t47 = sin(pkin(10));
t96 = rSges(5,2) * t47;
t60 = (rSges(5,3) + qJ(4)) * t42 + (-t77 - t96) * t43;
t115 = rSges(7,1) + pkin(5);
t116 = qJ(4) * t43 + t42 * t96;
t86 = rSges(7,3) + qJ(6);
t112 = t43 * rSges(4,1) - rSges(4,2) * t42;
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t111 = g(1) * t53 + g(2) * t51;
t110 = t111 * t42;
t50 = sin(qJ(2));
t109 = pkin(2) * t50;
t108 = pkin(4) * t47;
t52 = cos(qJ(2));
t44 = t52 * pkin(2);
t39 = t44 + pkin(1);
t21 = t53 * t39;
t106 = g(2) * t21;
t54 = -pkin(8) - pkin(7);
t104 = g(2) * t54;
t49 = -pkin(9) - qJ(4);
t103 = g(3) * t49;
t101 = rSges(3,3) + pkin(7);
t45 = pkin(10) + qJ(5);
t40 = sin(t45);
t95 = t40 * t42;
t81 = rSges(6,2) * t95;
t91 = t51 * t43;
t100 = rSges(6,3) * t91 + t51 * t81;
t89 = t53 * t43;
t99 = rSges(6,3) * t89 + t53 * t81;
t94 = t40 * t43;
t41 = cos(t45);
t93 = t41 * t43;
t35 = t42 * rSges(7,2);
t33 = t42 * rSges(6,3);
t92 = t42 * t49;
t37 = pkin(4) * t48 + pkin(3);
t12 = t43 * t37;
t90 = t53 * t42;
t88 = rSges(4,3) - t54;
t85 = t51 * t109;
t84 = t53 * t109;
t83 = rSges(5,3) * t91 + t116 * t51;
t80 = t49 * t91;
t79 = t49 * t89;
t78 = rSges(5,3) * t89 + t116 * t53;
t76 = -rSges(6,1) * t41 - t37;
t74 = -t39 - t12;
t73 = rSges(7,2) * t91 - t80;
t72 = rSges(7,2) * t89 - t79;
t71 = t51 * t92 + (t108 - t54) * t53;
t70 = rSges(3,1) * t52 - rSges(3,2) * t50;
t67 = -rSges(4,1) * t42 - rSges(4,2) * t43;
t66 = t115 * t93 + t86 * t94 + t12 + t35;
t65 = pkin(1) + t70;
t63 = t47 * rSges(5,1) + t48 * rSges(5,2) - t54;
t62 = t51 * t108 + t37 * t89 - t49 * t90 + t21;
t61 = rSges(6,1) * t93 - rSges(6,2) * t94 + t12 + t33;
t56 = t77 * t110;
t55 = (-t103 + t111 * (-t115 * t41 - t86 * t40 - t37)) * t42;
t5 = t40 * t51 + t41 * t89;
t4 = t40 * t89 - t51 * t41;
t3 = -t53 * t40 + t41 * t91;
t2 = t40 * t91 + t41 * t53;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t51 - rSges(2,2) * t53) + g(2) * (rSges(2,1) * t53 - rSges(2,2) * t51)) - m(3) * ((g(1) * t101 + g(2) * t65) * t53 + (-g(1) * t65 + g(2) * t101) * t51) - m(4) * (t106 + (g(1) * t88 + g(2) * t112) * t53 + (g(1) * (-t39 - t112) + g(2) * t88) * t51) - m(5) * (t106 + (g(1) * t63 + t60 * g(2)) * t53 + (g(2) * t63 + (-t39 - t60) * g(1)) * t51) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t71) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + rSges(6,3) * t90 + t62) + (g(1) * (t74 - t33) - t104) * t51) - m(7) * (g(1) * (-t115 * t3 - t86 * t2 + t71) + g(2) * (rSges(7,2) * t90 + t115 * t5 + t86 * t4 + t62) + (g(1) * (t74 - t35) - t104) * t51) -m(3) * (g(3) * t70 + t111 * (-rSges(3,1) * t50 - rSges(3,2) * t52)) - m(4) * (g(3) * (t44 + t112) + t111 * (t67 - t109)) - m(5) * (g(1) * (t78 - t84) + g(2) * (t83 - t85) + g(3) * (t44 + t60) + t56) - m(6) * (g(1) * t99 + g(2) * t100 + g(3) * (t44 + t61 - t92) + t111 * (t76 * t42 - t43 * t49 - t109)) - m(7) * (g(1) * (t72 - t84) + g(2) * (t73 - t85) + g(3) * (t44 + t66) + t55) -m(4) * (g(3) * t112 + t111 * t67) - m(5) * (g(1) * t78 + g(2) * t83 + g(3) * t60 + t56) - m(6) * (g(1) * (-t79 + t99) + g(2) * (-t80 + t100) + g(3) * t61 + (t111 * t76 - t103) * t42) - m(7) * (g(1) * t72 + g(2) * t73 + g(3) * t66 + t55) (-m(5) - m(6) - m(7)) * (-g(3) * t43 + t110) -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3)) - m(7) * (g(1) * (-t115 * t4 + t86 * t5) + g(2) * (-t115 * t2 + t86 * t3)) + (-m(6) * (-rSges(6,1) * t40 - rSges(6,2) * t41) - m(7) * (-t115 * t40 + t86 * t41)) * g(3) * t42, -m(7) * (g(1) * t4 + g(2) * t2 + g(3) * t95)];
taug  = t1(:);
