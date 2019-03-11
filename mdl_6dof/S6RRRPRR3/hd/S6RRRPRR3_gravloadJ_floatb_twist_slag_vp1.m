% Calculate Gravitation load on the joints for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:56
% EndTime: 2019-03-09 18:10:59
% DurationCPUTime: 0.93s
% Computational Cost: add. (619->144), mult. (739->175), div. (0->0), fcn. (761->10), ass. (0->75)
t49 = sin(qJ(6));
t52 = cos(qJ(6));
t112 = t52 * rSges(7,1) - t49 * rSges(7,2) + pkin(5);
t48 = qJ(2) + qJ(3);
t45 = sin(t48);
t46 = cos(t48);
t96 = sin(qJ(5));
t97 = cos(qJ(5));
t23 = t45 * t96 + t46 * t97;
t24 = t45 * t97 - t46 * t96;
t99 = rSges(7,3) + pkin(10);
t118 = -t112 * t23 + t24 * t99;
t54 = cos(qJ(1));
t79 = t54 * t96;
t80 = t54 * t97;
t16 = -t45 * t79 - t46 * t80;
t17 = -t45 * t80 + t46 * t79;
t117 = -t112 * t17 - t16 * t99;
t51 = sin(qJ(1));
t14 = t24 * t51;
t15 = t23 * t51;
t116 = t112 * t14 + t15 * t99;
t111 = t14 * rSges(6,1) - t15 * rSges(6,2);
t55 = -pkin(8) - pkin(7);
t98 = -pkin(9) - t55;
t110 = t17 * rSges(6,1) - t16 * rSges(6,2);
t109 = -t23 * rSges(6,1) - t24 * rSges(6,2);
t38 = t45 * qJ(4);
t89 = t46 * pkin(3) + t38;
t68 = t46 * rSges(5,1) + t45 * rSges(5,3);
t108 = t46 * rSges(4,1) - rSges(4,2) * t45;
t103 = g(1) * t54;
t107 = g(2) * t51 + t103;
t106 = t107 * t45;
t50 = sin(qJ(2));
t104 = pkin(2) * t50;
t42 = t46 * pkin(4);
t100 = rSges(3,3) + pkin(7);
t94 = t46 * t54;
t91 = rSges(5,2) - t55;
t90 = rSges(4,3) - t55;
t88 = qJ(4) * t46;
t87 = t51 * t104;
t86 = t54 * t104;
t85 = -rSges(6,3) + t98;
t53 = cos(qJ(2));
t47 = t53 * pkin(2);
t44 = t47 + pkin(1);
t32 = t54 * t44;
t84 = pkin(3) * t94 + t54 * t38 + t32;
t83 = t42 + t89;
t77 = pkin(4) * t94 + t84;
t76 = t89 + t68;
t29 = t51 * t88;
t75 = t29 - t87;
t31 = t54 * t88;
t74 = t31 - t86;
t73 = t83 - t109;
t72 = rSges(3,1) * t53 - rSges(3,2) * t50;
t69 = -rSges(4,1) * t45 - rSges(4,2) * t46;
t67 = -t15 * t52 - t49 * t54;
t66 = t15 * t49 - t52 * t54;
t65 = pkin(1) + t72;
t63 = -t44 - t89;
t61 = -t116 + t29;
t60 = -t117 + t31;
t59 = g(1) * (t63 - t42);
t58 = (-pkin(3) - pkin(4)) * t106;
t57 = (-rSges(5,1) - pkin(3)) * t106;
t56 = -t118 + t83;
t34 = rSges(5,3) * t94;
t33 = t51 * t46 * rSges(5,3);
t2 = -t16 * t52 - t49 * t51;
t1 = t16 * t49 - t51 * t52;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t51 - rSges(2,2) * t54) + g(2) * (rSges(2,1) * t54 - rSges(2,2) * t51)) - m(3) * ((g(1) * t100 + g(2) * t65) * t54 + (-g(1) * t65 + g(2) * t100) * t51) - m(4) * (g(2) * t32 + (g(1) * t90 + g(2) * t108) * t54 + (g(1) * (-t44 - t108) + g(2) * t90) * t51) - m(5) * (g(2) * t84 + (g(1) * t91 + g(2) * t68) * t54 + (g(1) * (t63 - t68) + g(2) * t91) * t51) - m(6) * (g(1) * (-t15 * rSges(6,1) - t14 * rSges(6,2)) + g(2) * (-t16 * rSges(6,1) - t17 * rSges(6,2) + t77) + t85 * t103 + (g(2) * t85 + t59) * t51) - m(7) * (g(1) * (rSges(7,1) * t67 + rSges(7,2) * t66 - t15 * pkin(5) + t99 * t14 + t98 * t54) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t16 * pkin(5) + t99 * t17 + t77) + (g(2) * t98 + t59) * t51) -m(3) * (g(3) * t72 + t107 * (-rSges(3,1) * t50 - rSges(3,2) * t53)) - m(4) * (g(3) * (t47 + t108) + t107 * (t69 - t104)) - m(5) * (g(1) * (t34 + t74) + g(2) * (t33 + t75) + g(3) * (t47 + t76) + t57) - m(6) * (g(1) * (t74 + t110) + g(2) * (t75 - t111) + g(3) * (t47 + t73) + t58) - m(7) * (g(1) * (t60 - t86) + g(2) * (t61 - t87) + g(3) * (t47 + t56) + t58) -m(4) * (g(3) * t108 + t107 * t69) - m(5) * (g(1) * (t31 + t34) + g(2) * (t29 + t33) + g(3) * t76 + t57) - m(6) * (g(1) * (t31 + t110) + g(2) * (t29 - t111) + g(3) * t73 + t58) - m(7) * (g(1) * t60 + g(2) * t61 + g(3) * t56 + t58) (-m(5) - m(6) - m(7)) * (-g(3) * t46 + t106) -m(6) * (-g(1) * t110 + g(2) * t111 + g(3) * t109) - m(7) * (t117 * g(1) + t116 * g(2) + g(3) * t118) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-rSges(7,1) * t66 + rSges(7,2) * t67) + g(3) * (-t49 * rSges(7,1) - t52 * rSges(7,2)) * t24)];
taug  = t3(:);
