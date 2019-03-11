% Calculate Gravitation load on the joints for
% S6RRRPRR10
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:39
% EndTime: 2019-03-09 19:15:44
% DurationCPUTime: 1.50s
% Computational Cost: add. (507->206), mult. (1027->270), div. (0->0), fcn. (1126->10), ass. (0->76)
t47 = sin(qJ(1));
t45 = sin(qJ(3));
t51 = cos(qJ(1));
t86 = t51 * t45;
t49 = cos(qJ(3));
t50 = cos(qJ(2));
t88 = t49 * t50;
t23 = t47 * t88 - t86;
t48 = cos(qJ(5));
t89 = t47 * t50;
t22 = t45 * t89 + t49 * t51;
t44 = sin(qJ(5));
t97 = t22 * t44;
t65 = t23 * t48 + t97;
t66 = t22 * t48 - t23 * t44;
t109 = rSges(6,1) * t66 - rSges(6,2) * t65;
t46 = sin(qJ(2));
t43 = qJ(5) + qJ(6);
t36 = sin(t43);
t37 = cos(t43);
t63 = t36 * t45 + t37 * t49;
t64 = t36 * t49 - t37 * t45;
t107 = (-rSges(7,1) * t64 - rSges(7,2) * t63) * t46;
t14 = t22 * t37;
t106 = -t23 * t36 + t14;
t101 = g(2) * t47;
t105 = g(1) * t51 + t101;
t104 = -pkin(3) - pkin(4);
t103 = g(1) * t47;
t100 = g(3) * t46;
t40 = t50 * pkin(2);
t99 = -rSges(5,1) - pkin(3);
t98 = -pkin(9) - rSges(6,3);
t94 = t44 * t45;
t93 = t44 * t49;
t92 = t45 * t46;
t91 = t45 * t50;
t90 = t46 * t51;
t87 = t50 * t51;
t85 = t46 * pkin(8) + t40;
t84 = t51 * pkin(1) + t47 * pkin(7);
t83 = -pkin(10) - pkin(9) - rSges(7,3);
t82 = rSges(5,3) + qJ(4);
t81 = -pkin(1) - t40;
t80 = t98 * t51;
t79 = t83 * t51;
t35 = pkin(5) * t48 + pkin(4);
t78 = t36 * rSges(7,2) - t35;
t77 = pkin(5) * t44 + qJ(4);
t75 = pkin(3) * t88 + qJ(4) * t91 + t85;
t74 = pkin(2) * t87 + pkin(8) * t90 + t84;
t41 = t51 * pkin(7);
t73 = -t23 * pkin(3) - t22 * qJ(4) + t41;
t25 = t45 * t47 + t49 * t87;
t72 = t25 * pkin(3) + t74;
t24 = -t47 * t49 + t50 * t86;
t8 = t24 * t48 - t25 * t44;
t9 = t24 * t44 + t25 * t48;
t71 = rSges(6,1) * t8 - rSges(6,2) * t9;
t70 = rSges(3,1) * t50 - rSges(3,2) * t46;
t61 = t48 * t49 + t94;
t62 = -t45 * t48 + t93;
t68 = (-rSges(6,1) * t62 - rSges(6,2) * t61) * t46;
t67 = -t22 * t36 - t23 * t37;
t59 = -t37 * rSges(7,1) + t78;
t58 = rSges(7,1) * t36 + t37 * rSges(7,2) + t77;
t28 = pkin(8) * t89;
t31 = pkin(8) * t87;
t56 = g(1) * t31 + g(2) * t28 + g(3) * t75;
t5 = t24 * t37 - t25 * t36;
t6 = t24 * t36 + t25 * t37;
t55 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (t106 * rSges(7,1) + t67 * rSges(7,2)) + g(3) * t107);
t26 = t46 * t49 * qJ(4);
t20 = t24 * pkin(3);
t18 = t22 * pkin(3);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t47 - rSges(2,2) * t51) + g(2) * (rSges(2,1) * t51 - rSges(2,2) * t47)) - m(3) * (g(1) * (rSges(3,3) * t51 + t41) + g(2) * (rSges(3,1) * t87 - rSges(3,2) * t90 + t84) + (g(1) * (-pkin(1) - t70) + g(2) * rSges(3,3)) * t47) - m(4) * (g(1) * (-rSges(4,1) * t23 + rSges(4,2) * t22 + t41) + g(2) * (rSges(4,1) * t25 - rSges(4,2) * t24 + rSges(4,3) * t90 + t74) + ((-rSges(4,3) - pkin(8)) * t46 + t81) * t103) - m(5) * (g(1) * (-rSges(5,1) * t23 - rSges(5,3) * t22 + t73) + g(2) * (rSges(5,1) * t25 + rSges(5,2) * t90 + t24 * t82 + t72) + ((-rSges(5,2) - pkin(8)) * t46 + t81) * t103) - m(6) * (g(1) * (-rSges(6,1) * t65 - rSges(6,2) * t66 - t23 * pkin(4) + t73) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + pkin(4) * t25 + t24 * qJ(4) + t46 * t80 + t72) + ((-pkin(8) - t98) * t46 + t81) * t103) - m(7) * (g(1) * (t67 * rSges(7,1) - t106 * rSges(7,2) - pkin(5) * t97 - t23 * t35 + t73) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t24 * t77 + t25 * t35 + t46 * t79 + t72) + ((-pkin(8) - t83) * t46 + t81) * t103) -m(3) * (g(3) * t70 + t105 * (-rSges(3,1) * t46 - rSges(3,2) * t50)) - m(4) * (g(1) * (rSges(4,3) * t87 + t31) + g(2) * (rSges(4,3) * t89 + t28) + g(3) * (rSges(4,1) * t88 - rSges(4,2) * t91 + t85) + (g(3) * rSges(4,3) + t105 * (-rSges(4,1) * t49 + rSges(4,2) * t45 - pkin(2))) * t46) - m(5) * (g(1) * (rSges(5,2) * t87 + t31) + g(2) * (rSges(5,2) * t89 + t28) + g(3) * (rSges(5,1) * t88 + rSges(5,3) * t91 + t75) + (g(3) * rSges(5,2) + t105 * (-t82 * t45 + t99 * t49 - pkin(2))) * t46) - m(6) * ((g(3) * (rSges(6,1) * t61 - rSges(6,2) * t62 + t49 * pkin(4)) + g(1) * t80 + t98 * t101) * t50 + (g(3) * t98 + t105 * (-pkin(2) + (-rSges(6,1) * t44 - rSges(6,2) * t48 - qJ(4)) * t45 + (-t48 * rSges(6,1) + t44 * rSges(6,2) + t104) * t49)) * t46 + t56) - m(7) * ((g(3) * (rSges(7,1) * t63 - rSges(7,2) * t64 + pkin(5) * t94 + t49 * t35) + g(1) * t79 + t83 * t101) * t50 + (g(3) * t83 + t105 * (-pkin(2) + (-pkin(3) + t59) * t49 - t58 * t45)) * t46 + t56) -m(4) * (g(1) * (-rSges(4,1) * t24 - rSges(4,2) * t25) + g(2) * (-rSges(4,1) * t22 - rSges(4,2) * t23) + (-rSges(4,1) * t45 - rSges(4,2) * t49) * t100) - m(5) * (g(1) * (-rSges(5,1) * t24 + t25 * t82 - t20) + g(2) * (-rSges(5,1) * t22 + t23 * t82 - t18) + g(3) * t26 + (rSges(5,3) * t49 + t99 * t45) * t100) - m(6) * (g(1) * (-t24 * pkin(4) + t25 * qJ(4) - t20 - t71) + g(2) * (-t22 * pkin(4) + t23 * qJ(4) - t109 - t18) + g(3) * (t104 * t92 + t26 - t68)) - m(7) * (g(1) * (t24 * t59 + t25 * t58 - t20) + g(2) * (-t14 * rSges(7,1) + t22 * t78 + t23 * t58 - t18) + g(3) * (t26 - t107) + (pkin(5) * t93 + (-pkin(3) - t35) * t45) * t100) (-m(5) - m(6) - m(7)) * (g(1) * t24 + g(2) * t22 + g(3) * t92) -m(6) * (g(1) * t71 + g(2) * t109 + g(3) * t68) - t55 - m(7) * (g(1) * t8 + g(2) * t66 - t62 * t100) * pkin(5), -t55];
taug  = t1(:);
