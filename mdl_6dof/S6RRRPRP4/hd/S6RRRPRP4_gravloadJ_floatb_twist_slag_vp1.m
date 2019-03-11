% Calculate Gravitation load on the joints for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:36
% EndTime: 2019-03-09 16:43:38
% DurationCPUTime: 0.85s
% Computational Cost: add. (473->148), mult. (571->190), div. (0->0), fcn. (531->8), ass. (0->76)
t98 = rSges(7,1) + pkin(5);
t45 = sin(qJ(1));
t42 = qJ(2) + qJ(3);
t38 = cos(t42);
t43 = sin(qJ(5));
t89 = t43 * t45;
t77 = t38 * t89;
t46 = cos(qJ(5));
t91 = t38 * t46;
t78 = rSges(6,2) * t91;
t110 = rSges(6,1) * t77 + t45 * t78;
t109 = t98 * t77;
t48 = cos(qJ(1));
t88 = t43 * t48;
t76 = t38 * t88;
t108 = t98 * t76;
t37 = sin(t42);
t107 = t37 * t46;
t106 = rSges(6,1) * t76 + t48 * t78;
t58 = -rSges(5,2) * t38 + t37 * rSges(5,3);
t105 = t38 * rSges(4,1) - rSges(4,2) * t37;
t81 = rSges(7,3) + qJ(6);
t100 = g(2) * t45;
t104 = g(1) * t48 + t100;
t103 = -pkin(3) - pkin(9);
t44 = sin(qJ(2));
t102 = pkin(2) * t44;
t99 = g(3) * t38;
t35 = t38 * pkin(3);
t97 = rSges(3,3) + pkin(7);
t94 = t37 * t43;
t93 = t37 * t45;
t92 = t37 * t48;
t90 = t38 * t48;
t87 = t45 * t46;
t86 = t46 * t48;
t49 = -pkin(8) - pkin(7);
t85 = rSges(5,1) - t49;
t84 = rSges(4,3) - t49;
t29 = t37 * qJ(4);
t83 = t29 + t35;
t82 = qJ(4) * t38;
t80 = -rSges(7,2) + t103;
t79 = -rSges(6,3) + t103;
t15 = t45 * t82;
t75 = t45 * t38 * rSges(5,3) + rSges(5,2) * t93 + t15;
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t36 = t40 + pkin(1);
t18 = t48 * t36;
t74 = pkin(3) * t90 + t48 * t29 + t18;
t17 = t48 * t82;
t73 = rSges(5,2) * t92 + rSges(5,3) * t90 + t17;
t72 = t38 * pkin(9) + t83;
t70 = (pkin(4) - t49) * t48;
t69 = -t36 - t29;
t68 = g(1) * t80;
t67 = g(1) * t79;
t66 = -t102 * t45 + t15;
t65 = -t102 * t48 + t17;
t64 = t45 * pkin(4) + pkin(9) * t90 + t74;
t63 = -pkin(3) * t37 - t102;
t62 = rSges(3,1) * t47 - rSges(3,2) * t44;
t59 = -rSges(4,1) * t37 - rSges(4,2) * t38;
t57 = t83 + t58;
t56 = t38 * rSges(7,2) + t98 * t94 + t72;
t55 = rSges(6,1) * t94 + rSges(6,2) * t107 + t38 * rSges(6,3) + t72;
t54 = pkin(1) + t62;
t52 = g(1) * t69 - g(2) * t49;
t51 = (t100 * t79 + t48 * t67) * t37;
t50 = (t100 * t80 + t48 * t68) * t37 + (-g(3) * t107 - t104 * t91) * t81;
t5 = -t37 * t89 + t86;
t4 = t37 * t87 + t88;
t3 = t37 * t88 + t87;
t2 = -t37 * t86 + t89;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t45 - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - rSges(2,2) * t45)) - m(3) * ((g(1) * t97 + g(2) * t54) * t48 + (-g(1) * t54 + g(2) * t97) * t45) - m(4) * (g(2) * t18 + (g(1) * t84 + g(2) * t105) * t48 + (g(1) * (-t36 - t105) + g(2) * t84) * t45) - m(5) * (g(2) * t74 + (g(1) * t85 + g(2) * t58) * t48 + (g(1) * (-t58 + t69 - t35) + g(2) * t85) * t45) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t70) + g(2) * (t3 * rSges(6,1) - t2 * rSges(6,2) + rSges(6,3) * t90 + t64) + (t38 * t67 + t52) * t45) - m(7) * (g(1) * (t81 * t4 + t98 * t5 + t70) + g(2) * (rSges(7,2) * t90 + t81 * t2 + t98 * t3 + t64) + (t38 * t68 + t52) * t45) -m(3) * (g(3) * t62 + t104 * (-rSges(3,1) * t44 - rSges(3,2) * t47)) - m(4) * (g(3) * (t40 + t105) + t104 * (t59 - t102)) - m(5) * (g(1) * (t48 * t63 + t73) + g(2) * (t45 * t63 + t75) + g(3) * (t40 + t57)) - m(6) * (g(1) * (t65 + t106) + g(2) * (t66 + t110) + g(3) * (t40 + t55) + t51) - m(7) * (g(1) * (t65 + t108) + g(2) * (t66 + t109) + g(3) * (t40 + t56) + t50) -m(4) * (g(3) * t105 + t104 * t59) - m(5) * (g(1) * (-pkin(3) * t92 + t73) + g(2) * (-pkin(3) * t93 + t75) + g(3) * t57) - m(6) * (g(1) * (t17 + t106) + g(2) * (t15 + t110) + g(3) * t55 + t51) - m(7) * (g(1) * (t17 + t108) + g(2) * (t15 + t109) + g(3) * t56 + t50) (-m(5) - m(6) - m(7)) * (t104 * t37 - t99) -m(6) * (g(1) * (-rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5)) - m(7) * (g(1) * (-t2 * t98 + t3 * t81) + g(2) * (t4 * t98 - t5 * t81)) + (-m(6) * (-rSges(6,1) * t46 + rSges(6,2) * t43) - m(7) * (-t81 * t43 - t98 * t46)) * t99, -m(7) * (g(1) * t2 - g(2) * t4 + g(3) * t91)];
taug  = t1(:);
