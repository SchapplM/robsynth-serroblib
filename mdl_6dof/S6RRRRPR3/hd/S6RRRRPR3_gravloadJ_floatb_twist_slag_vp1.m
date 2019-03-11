% Calculate Gravitation load on the joints for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:02:01
% EndTime: 2019-03-09 22:02:03
% DurationCPUTime: 0.78s
% Computational Cost: add. (593->139), mult. (494->174), div. (0->0), fcn. (430->10), ass. (0->81)
t42 = sin(qJ(6));
t45 = cos(qJ(6));
t111 = rSges(7,1) * t42 + rSges(7,2) * t45;
t100 = rSges(7,3) + pkin(10);
t41 = qJ(2) + qJ(3);
t38 = qJ(4) + t41;
t34 = cos(t38);
t110 = t100 * t34;
t109 = t111 * t34;
t33 = sin(t38);
t62 = -rSges(6,2) * t34 + t33 * rSges(6,3);
t108 = t34 * rSges(5,1) - rSges(5,2) * t33;
t36 = sin(t41);
t37 = cos(t41);
t73 = t37 * rSges(4,1) - rSges(4,2) * t36;
t44 = sin(qJ(1));
t103 = g(2) * t44;
t47 = cos(qJ(1));
t107 = g(1) * t47 + t103;
t48 = -pkin(8) - pkin(7);
t43 = sin(qJ(2));
t106 = pkin(2) * t43;
t105 = pkin(3) * t36;
t102 = g(3) * t34;
t30 = t34 * pkin(4);
t101 = rSges(3,3) + pkin(7);
t40 = -pkin(9) + t48;
t99 = pkin(5) - t40;
t46 = cos(qJ(2));
t39 = t46 * pkin(2);
t35 = t39 + pkin(1);
t93 = t33 * t44;
t92 = t33 * t47;
t91 = t34 * t47;
t90 = t42 * t44;
t89 = t42 * t47;
t88 = t44 * t45;
t87 = t45 * t47;
t86 = rSges(6,1) - t40;
t85 = rSges(4,3) - t48;
t84 = rSges(5,3) - t40;
t25 = t33 * qJ(5);
t83 = t25 + t30;
t82 = qJ(5) * t34;
t81 = -pkin(4) - t100;
t15 = t44 * t82;
t80 = t109 * t44 + t15;
t32 = pkin(3) * t37;
t13 = t32 + t35;
t6 = t47 * t13;
t79 = pkin(4) * t91 + t47 * t25 + t6;
t17 = t47 * t82;
t76 = t109 * t47 + t17;
t75 = t44 * t34 * rSges(6,3) + rSges(6,2) * t93 + t15;
t74 = rSges(6,2) * t92 + rSges(6,3) * t91 + t17;
t71 = -t13 - t25;
t70 = g(1) * t81;
t69 = t32 + t108;
t68 = -pkin(4) * t33 - t105;
t67 = rSges(3,1) * t46 - rSges(3,2) * t43;
t65 = -rSges(4,1) * t36 - rSges(4,2) * t37;
t63 = -rSges(5,1) * t33 - rSges(5,2) * t34;
t61 = t83 + t62;
t60 = t111 * t33 + t110 + t83;
t59 = pkin(1) + t67;
t58 = t35 + t73;
t57 = -pkin(4) * t93 + t75;
t56 = -pkin(4) * t92 + t74;
t55 = t32 + t61;
t54 = t63 * t44;
t53 = t63 * t47;
t52 = t32 + t60;
t49 = (t81 * t103 + t47 * t70) * t33;
t14 = -t105 - t106;
t8 = t47 * t14;
t7 = t44 * t14;
t5 = -t33 * t90 + t87;
t4 = t33 * t88 + t89;
t3 = t33 * t89 + t88;
t2 = t33 * t87 - t90;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t44 - rSges(2,2) * t47) + g(2) * (rSges(2,1) * t47 - rSges(2,2) * t44)) - m(3) * ((g(1) * t101 + g(2) * t59) * t47 + (-g(1) * t59 + g(2) * t101) * t44) - m(4) * ((g(1) * t85 + g(2) * t58) * t47 + (-g(1) * t58 + g(2) * t85) * t44) - m(5) * (g(2) * t6 + (g(1) * t84 + g(2) * t108) * t47 + (g(1) * (-t13 - t108) + g(2) * t84) * t44) - m(6) * (g(2) * t79 + (g(1) * t86 + g(2) * t62) * t47 + (g(1) * (-t62 + t71 - t30) + g(2) * t86) * t44) - m(7) * (g(1) * (rSges(7,1) * t5 - rSges(7,2) * t4) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + t79) + (g(1) * t99 + g(2) * t110) * t47 + (g(1) * t71 + g(2) * t99 + t34 * t70) * t44) -m(3) * (g(3) * t67 + t107 * (-rSges(3,1) * t43 - rSges(3,2) * t46)) - m(4) * (g(3) * (t39 + t73) + t107 * (t65 - t106)) - m(5) * (g(1) * (t8 + t53) + g(2) * (t7 + t54) + g(3) * (t39 + t69)) - m(6) * (g(1) * (t56 + t8) + g(2) * (t57 + t7) + g(3) * (t39 + t55)) - m(7) * (g(1) * (t8 + t76) + g(2) * (t7 + t80) + g(3) * (t39 + t52) + t49) -m(4) * (g(3) * t73 + t107 * t65) - m(5) * (g(3) * t69 + t107 * (t63 - t105)) - m(6) * (g(1) * (t68 * t47 + t74) + g(2) * (t68 * t44 + t75) + g(3) * t55) - m(7) * (g(1) * (-t47 * t105 + t76) + g(2) * (-t44 * t105 + t80) + g(3) * t52 + t49) -m(5) * (g(1) * t53 + g(2) * t54 + g(3) * t108) - m(6) * (g(1) * t56 + g(2) * t57 + g(3) * t61) - m(7) * (g(1) * t76 + g(2) * t80 + g(3) * t60 + t49) (-m(6) - m(7)) * (t107 * t33 - t102) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (rSges(7,1) * t4 + rSges(7,2) * t5) + (-rSges(7,1) * t45 + rSges(7,2) * t42) * t102)];
taug  = t1(:);
