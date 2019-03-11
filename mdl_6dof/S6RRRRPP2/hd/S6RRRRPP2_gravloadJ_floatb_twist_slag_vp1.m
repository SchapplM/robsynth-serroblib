% Calculate Gravitation load on the joints for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:57
% EndTime: 2019-03-09 20:49:00
% DurationCPUTime: 0.93s
% Computational Cost: add. (561->164), mult. (710->195), div. (0->0), fcn. (689->8), ass. (0->75)
t105 = rSges(7,1) + pkin(5);
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t101 = g(1) * t47 + g(2) * t44;
t41 = qJ(2) + qJ(3);
t38 = sin(t41);
t100 = t101 * t38;
t106 = g(3) * t38;
t42 = sin(qJ(4));
t90 = t38 * t42;
t73 = rSges(5,2) * t90;
t39 = cos(t41);
t87 = t39 * t44;
t104 = rSges(5,3) * t87 + t44 * t73;
t85 = t39 * t47;
t103 = rSges(5,3) * t85 + t47 * t73;
t102 = t39 * rSges(4,1) - rSges(4,2) * t38;
t75 = -rSges(7,3) - qJ(6);
t43 = sin(qJ(2));
t99 = pkin(2) * t43;
t48 = -pkin(8) - pkin(7);
t96 = g(2) * t48;
t45 = cos(qJ(4));
t95 = t45 * qJ(5) * t106;
t33 = t38 * pkin(9);
t34 = t39 * pkin(3);
t94 = -rSges(6,1) - pkin(4);
t92 = rSges(3,3) + pkin(7);
t31 = t38 * rSges(6,2);
t30 = t38 * rSges(5,3);
t89 = t38 * t47;
t88 = t39 * t42;
t86 = t39 * t45;
t84 = t42 * t44;
t83 = t44 * t45;
t82 = t45 * t47;
t81 = t47 * t42;
t80 = t47 * t48;
t79 = rSges(4,3) - t48;
t78 = t33 + t34;
t77 = rSges(7,2) + qJ(5);
t76 = rSges(6,3) + qJ(5);
t74 = -pkin(4) - t105;
t46 = cos(qJ(2));
t40 = t46 * pkin(2);
t37 = t40 + pkin(1);
t14 = t47 * t37;
t72 = pkin(3) * t85 + pkin(9) * t89 + t14;
t70 = -t37 - t34;
t68 = pkin(4) * t86 + qJ(5) * t88 + t78;
t24 = pkin(9) * t87;
t67 = -t44 * t99 + t24;
t28 = pkin(9) * t85;
t66 = -t47 * t99 + t28;
t65 = rSges(3,1) * t46 - rSges(3,2) * t43;
t62 = -rSges(4,1) * t38 - rSges(4,2) * t39;
t61 = t70 - t33;
t60 = pkin(1) + t65;
t6 = t39 * t84 + t82;
t7 = t39 * t83 - t81;
t59 = -t7 * pkin(4) - t6 * qJ(5) - t80;
t8 = t39 * t81 - t83;
t9 = t39 * t82 + t84;
t58 = t9 * pkin(4) + t8 * qJ(5) + t72;
t57 = rSges(7,2) * t88 + t105 * t86 + t68;
t56 = rSges(6,1) * t86 + rSges(6,3) * t88 + t31 + t68;
t55 = rSges(5,1) * t86 - rSges(5,2) * t88 + t30 + t78;
t51 = (-rSges(5,1) * t45 - pkin(3)) * t100;
t50 = (-t76 * t42 + t94 * t45 - pkin(3)) * t100;
t49 = (-t77 * t42 + t74 * t45 - pkin(3)) * t100 + (t101 * t39 + t106) * t75;
t23 = rSges(6,2) * t85;
t16 = rSges(6,2) * t87;
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t44 - rSges(2,2) * t47) + g(2) * (rSges(2,1) * t47 - rSges(2,2) * t44)) - m(3) * ((g(1) * t92 + g(2) * t60) * t47 + (-g(1) * t60 + g(2) * t92) * t44) - m(4) * (g(2) * t14 + (g(1) * t79 + g(2) * t102) * t47 + (g(1) * (-t37 - t102) + g(2) * t79) * t44) - m(5) * (g(1) * (-t7 * rSges(5,1) + t6 * rSges(5,2) - t80) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + rSges(5,3) * t89 + t72) + (g(1) * (t61 - t30) - t96) * t44) - m(6) * (g(1) * (-t7 * rSges(6,1) - t6 * rSges(6,3) + t59) + g(2) * (t9 * rSges(6,1) + rSges(6,2) * t89 + t8 * rSges(6,3) + t58) + (g(1) * (t61 - t31) - t96) * t44) - m(7) * (g(1) * (-t6 * rSges(7,2) - t105 * t7 + t59) + g(2) * (t8 * rSges(7,2) + t105 * t9 + t75 * t89 + t58) + (-t96 + (t70 + (-pkin(9) - t75) * t38) * g(1)) * t44) -m(3) * (g(3) * t65 + t101 * (-rSges(3,1) * t43 - rSges(3,2) * t46)) - m(4) * (g(3) * (t40 + t102) + t101 * (t62 - t99)) - m(5) * (g(1) * (t66 + t103) + g(2) * (t67 + t104) + g(3) * (t40 + t55) + t51) - m(6) * (g(1) * (t23 + t66) + g(2) * (t16 + t67) + g(3) * (t40 + t56) + t50) - m(7) * (g(1) * t66 + g(2) * t67 + g(3) * (t40 + t57) + t49) -m(4) * (g(3) * t102 + t101 * t62) - m(5) * (g(1) * (t28 + t103) + g(2) * (t24 + t104) + g(3) * t55 + t51) - m(6) * (g(1) * (t23 + t28) + g(2) * (t16 + t24) + g(3) * t56 + t50) - m(7) * (g(1) * t28 + g(2) * t24 + g(3) * t57 + t49) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (-rSges(6,1) * t8 + t76 * t9 - t4) + g(2) * (-rSges(6,1) * t6 + t76 * t7 - t2) + t95) - m(7) * (g(1) * (-t105 * t8 + t77 * t9 - t4) + g(2) * (-t105 * t6 + t77 * t7 - t2) + t95) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t45 + (m(5) * rSges(5,1) - m(6) * t94 - m(7) * t74) * t42) * t106 (-m(6) - m(7)) * (g(1) * t8 + g(2) * t6 + g(3) * t90) -m(7) * (g(3) * t39 - t100)];
taug  = t1(:);
