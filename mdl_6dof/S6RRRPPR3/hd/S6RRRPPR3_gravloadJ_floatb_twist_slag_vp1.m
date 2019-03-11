% Calculate Gravitation load on the joints for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:25
% EndTime: 2019-03-09 15:28:27
% DurationCPUTime: 0.73s
% Computational Cost: add. (415->138), mult. (470->164), div. (0->0), fcn. (415->8), ass. (0->78)
t107 = rSges(7,3) + pkin(9);
t38 = qJ(2) + qJ(3);
t36 = cos(t38);
t39 = sin(qJ(6));
t86 = rSges(7,2) * t39;
t42 = cos(qJ(6));
t90 = rSges(7,1) * t42;
t106 = (-t86 + t90) * t36;
t41 = sin(qJ(1));
t85 = t36 * t41;
t35 = sin(t38);
t88 = rSges(6,2) * t35;
t105 = rSges(6,1) * t85 + t41 * t88;
t44 = cos(qJ(1));
t84 = t36 * t44;
t104 = rSges(6,1) * t84 + t44 * t88;
t26 = t35 * rSges(6,1);
t103 = -rSges(6,2) * t36 + t26;
t55 = t36 * rSges(5,1) + t35 * rSges(5,3);
t102 = t36 * rSges(4,1) - rSges(4,2) * t35;
t94 = g(2) * t41;
t101 = g(1) * t44 + t94;
t30 = t35 * pkin(5);
t100 = t107 * t36 + t30;
t99 = t101 * t35;
t98 = -m(6) - m(7);
t97 = -pkin(3) - pkin(4);
t40 = sin(qJ(2));
t96 = pkin(2) * t40;
t93 = g(3) * t36;
t33 = t36 * pkin(3);
t91 = rSges(3,3) + pkin(7);
t83 = t39 * t41;
t82 = t39 * t44;
t81 = t41 * t42;
t80 = t42 * t44;
t45 = -pkin(8) - pkin(7);
t79 = rSges(5,2) - t45;
t78 = rSges(4,3) - t45;
t24 = t35 * qJ(4);
t77 = t24 + t33;
t76 = qJ(4) * t36;
t75 = -qJ(5) - t45;
t74 = t41 * t96;
t73 = t44 * t96;
t43 = cos(qJ(2));
t37 = t43 * pkin(2);
t34 = t37 + pkin(1);
t12 = t44 * t34;
t70 = pkin(3) * t84 + t44 * t24 + t12;
t69 = t36 * pkin(4) + t77;
t68 = -rSges(6,3) + t75;
t67 = t97 - t107;
t65 = -t34 - t24;
t64 = pkin(4) * t84 + t70;
t63 = t77 + t55;
t9 = t41 * t76;
t62 = t9 - t74;
t11 = t44 * t76;
t61 = t11 - t73;
t60 = g(1) * t67;
t59 = rSges(3,1) * t43 - rSges(3,2) * t40;
t56 = -rSges(4,1) * t35 - rSges(4,2) * t36;
t54 = pkin(1) + t59;
t53 = t69 + t103;
t52 = t35 * t90 + t100 + t69;
t51 = pkin(5) * t85 + t106 * t41 + t9;
t50 = pkin(5) * t84 + t106 * t44 + t11;
t48 = t97 * t99;
t47 = (-rSges(5,1) - pkin(3)) * t99;
t46 = (-g(3) * t86 + t44 * t60 + t67 * t94) * t35;
t18 = rSges(5,3) * t84;
t15 = rSges(5,3) * t85;
t5 = t35 * t80 - t83;
t4 = -t35 * t82 - t81;
t3 = -t35 * t81 - t82;
t2 = t35 * t83 - t80;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t41 - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - rSges(2,2) * t41)) - m(3) * ((g(1) * t91 + g(2) * t54) * t44 + (-g(1) * t54 + g(2) * t91) * t41) - m(4) * (g(2) * t12 + (g(1) * t78 + g(2) * t102) * t44 + (g(1) * (-t34 - t102) + g(2) * t78) * t41) - m(5) * (g(2) * t70 + (g(1) * t79 + g(2) * t55) * t44 + (g(1) * (-t55 + t65 - t33) + g(2) * t79) * t41) - m(6) * (g(2) * t64 + (g(1) * t68 + g(2) * t103) * t44 + (g(2) * t68 + (t65 - t26 + (rSges(6,2) + t97) * t36) * g(1)) * t41) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2)) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t64) + (g(1) * t75 + g(2) * t100) * t44 + (g(1) * (t65 - t30) + g(2) * t75 + t36 * t60) * t41) -m(3) * (g(3) * t59 + t101 * (-rSges(3,1) * t40 - rSges(3,2) * t43)) - m(4) * (g(3) * (t37 + t102) + t101 * (t56 - t96)) - m(5) * (g(1) * (t18 + t61) + g(2) * (t15 + t62) + g(3) * (t37 + t63) + t47) - m(6) * (g(1) * (t61 + t104) + g(2) * (t62 + t105) + g(3) * (t37 + t53) + t48) - m(7) * (g(1) * (t50 - t73) + g(2) * (t51 - t74) + g(3) * (t37 + t52) + t46) -m(4) * (g(3) * t102 + t101 * t56) - m(5) * (g(1) * (t11 + t18) + g(2) * (t15 + t9) + g(3) * t63 + t47) - m(6) * (g(1) * (t11 + t104) + g(2) * (t9 + t105) + g(3) * t53 + t48) - m(7) * (g(1) * t50 + g(2) * t51 + g(3) * t52 + t46) (-m(5) + t98) * (-t93 + t99) t98 * (-g(1) * t41 + g(2) * t44) -m(7) * (g(1) * (rSges(7,1) * t4 - rSges(7,2) * t5) + g(2) * (-rSges(7,1) * t2 + rSges(7,2) * t3) + (rSges(7,1) * t39 + rSges(7,2) * t42) * t93)];
taug  = t1(:);
