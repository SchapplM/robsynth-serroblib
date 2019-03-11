% Calculate Gravitation load on the joints for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:35
% EndTime: 2019-03-09 17:41:38
% DurationCPUTime: 1.46s
% Computational Cost: add. (654->204), mult. (1558->286), div. (0->0), fcn. (1860->10), ass. (0->78)
t57 = sin(qJ(3));
t61 = cos(qJ(3));
t117 = pkin(3) * t61 + qJ(4) * t57;
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t62 = cos(qJ(2));
t83 = cos(pkin(6));
t99 = cos(qJ(1));
t71 = t83 * t99;
t38 = t58 * t71 + t59 * t62;
t54 = sin(pkin(6));
t78 = t54 * t99;
t18 = t38 * t57 + t61 * t78;
t37 = t58 * t59 - t62 * t71;
t56 = sin(qJ(5));
t60 = cos(qJ(5));
t116 = -t18 * t56 - t37 * t60;
t2 = t18 * t60 - t37 * t56;
t101 = rSges(6,3) + pkin(10);
t88 = rSges(7,3) + qJ(6) + pkin(10);
t75 = t59 * t83;
t39 = t99 * t58 + t62 * t75;
t115 = g(1) * t39 + g(2) * t37;
t19 = t38 * t61 - t57 * t78;
t40 = -t58 * t75 + t99 * t62;
t95 = t54 * t59;
t23 = t40 * t61 + t57 * t95;
t94 = t54 * t61;
t36 = t83 * t57 + t58 * t94;
t114 = g(1) * t23 + g(2) * t19 + g(3) * t36;
t92 = t56 * t57;
t113 = pkin(5) * t92 + t88 * t61;
t112 = m(7) * t114;
t111 = pkin(4) + pkin(9);
t104 = g(3) * t54;
t103 = rSges(5,1) + pkin(9);
t102 = rSges(4,3) + pkin(9);
t52 = pkin(5) * t60 + pkin(4);
t100 = pkin(9) + t52;
t96 = t54 * t58;
t93 = t54 * t62;
t91 = t56 * t62;
t90 = t57 * t60;
t89 = t60 * t62;
t87 = t99 * pkin(1) + pkin(8) * t95;
t86 = pkin(2) * t93 + pkin(9) * t96;
t84 = rSges(5,3) + qJ(4);
t31 = t37 * pkin(2);
t81 = -t117 * t37 - t31;
t33 = t39 * pkin(2);
t80 = -t117 * t39 - t33;
t79 = t40 * pkin(2) + t87;
t77 = -t59 * pkin(1) + pkin(8) * t78;
t76 = pkin(5) * t56 + qJ(4);
t74 = t23 * pkin(3) + t79;
t73 = t117 * t93 + t86;
t72 = -t38 * pkin(2) + t77;
t70 = rSges(4,1) * t61 - rSges(4,2) * t57;
t69 = rSges(5,2) * t61 - rSges(5,3) * t57;
t22 = t40 * t57 - t59 * t94;
t6 = t22 * t60 - t39 * t56;
t68 = -pkin(3) * t19 + t72;
t67 = pkin(9) * t38 + t81;
t66 = pkin(9) * t40 + t80;
t35 = t57 * t96 - t83 * t61;
t16 = t35 * t60 + t54 * t91;
t30 = t35 * pkin(3);
t25 = (t57 * t91 + t58 * t60) * t54;
t24 = (-t56 * t58 + t57 * t89) * t54;
t17 = -t35 * t56 + t54 * t89;
t14 = t22 * pkin(3);
t12 = t18 * pkin(3);
t11 = -t39 * t92 + t40 * t60;
t10 = -t39 * t90 - t40 * t56;
t9 = -t37 * t92 + t38 * t60;
t8 = -t37 * t90 - t38 * t56;
t7 = t22 * t56 + t39 * t60;
t1 = [-m(2) * (g(1) * (-t59 * rSges(2,1) - t99 * rSges(2,2)) + g(2) * (t99 * rSges(2,1) - t59 * rSges(2,2))) - m(3) * (g(1) * (-t38 * rSges(3,1) + t37 * rSges(3,2) + rSges(3,3) * t78 + t77) + g(2) * (rSges(3,1) * t40 - rSges(3,2) * t39 + rSges(3,3) * t95 + t87)) - m(4) * (g(1) * (-rSges(4,1) * t19 + rSges(4,2) * t18 - t102 * t37 + t72) + g(2) * (rSges(4,1) * t23 - rSges(4,2) * t22 + t102 * t39 + t79)) - m(5) * (g(1) * (rSges(5,2) * t19 - t103 * t37 - t18 * t84 + t68) + g(2) * (-rSges(5,2) * t23 + t103 * t39 + t84 * t22 + t74)) - m(6) * (g(1) * (rSges(6,1) * t116 - rSges(6,2) * t2 - qJ(4) * t18 - t101 * t19 - t111 * t37 + t68) + g(2) * (rSges(6,1) * t7 + rSges(6,2) * t6 + qJ(4) * t22 + t101 * t23 + t111 * t39 + t74)) - m(7) * (g(1) * (rSges(7,1) * t116 - rSges(7,2) * t2 - t100 * t37 - t18 * t76 - t19 * t88 + t68) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 + t100 * t39 + t76 * t22 + t88 * t23 + t74)) -m(3) * (g(1) * (-rSges(3,1) * t39 - rSges(3,2) * t40) + g(2) * (-rSges(3,1) * t37 - rSges(3,2) * t38) + (rSges(3,1) * t62 - rSges(3,2) * t58) * t104) - m(4) * (g(1) * (t102 * t40 - t70 * t39 - t33) + g(2) * (t102 * t38 - t70 * t37 - t31) + g(3) * t86 + (rSges(4,3) * t58 + t70 * t62) * t104) - m(5) * (g(1) * (t103 * t40 + t69 * t39 + t80) + g(2) * (t103 * t38 + t69 * t37 + t81) + g(3) * t73 + (rSges(5,1) * t58 - t69 * t62) * t104) - m(6) * (g(1) * (rSges(6,1) * t11 + rSges(6,2) * t10 + pkin(4) * t40 + t66) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + pkin(4) * t38 + t67) + g(3) * (t25 * rSges(6,1) + t24 * rSges(6,2) + pkin(4) * t96 + t73) + (g(3) * t93 - t115) * t61 * t101) - m(7) * (g(1) * (rSges(7,1) * t11 + rSges(7,2) * t10 + t40 * t52 + t66) + g(2) * (rSges(7,1) * t9 + rSges(7,2) * t8 + t38 * t52 + t67) + g(3) * (t25 * rSges(7,1) + t24 * rSges(7,2) + t73) + (t113 * t62 + t52 * t58) * t104 - t115 * t113) -m(4) * (g(1) * (-rSges(4,1) * t22 - rSges(4,2) * t23) + g(2) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(3) * (-rSges(4,1) * t35 - rSges(4,2) * t36)) - m(5) * (g(1) * (rSges(5,2) * t22 + t84 * t23 - t14) + g(2) * (rSges(5,2) * t18 + t84 * t19 - t12) + g(3) * (rSges(5,2) * t35 + t84 * t36 - t30)) - m(7) * (g(1) * (-t22 * t88 - t14) + g(2) * (-t18 * t88 - t12) + g(3) * (-t35 * t88 - t30)) - (rSges(7,2) * t60 + qJ(4) + (rSges(7,1) + pkin(5)) * t56) * t112 + (-g(1) * (-t101 * t22 - t14) - g(2) * (-t101 * t18 - t12) - g(3) * (-t101 * t35 - t30) - t114 * (rSges(6,1) * t56 + rSges(6,2) * t60 + qJ(4))) * m(6) (-m(5) - m(6) - m(7)) * (g(1) * t22 + g(2) * t18 + g(3) * t35) -m(6) * (g(1) * (rSges(6,1) * t6 - rSges(6,2) * t7) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t116) + g(3) * (rSges(6,1) * t16 + rSges(6,2) * t17)) + (-g(1) * (rSges(7,1) * t6 - rSges(7,2) * t7) - g(2) * (rSges(7,1) * t2 + rSges(7,2) * t116) - g(3) * (t16 * rSges(7,1) + t17 * rSges(7,2)) - (g(1) * t6 + g(2) * t2 + g(3) * t16) * pkin(5)) * m(7), -t112];
taug  = t1(:);
