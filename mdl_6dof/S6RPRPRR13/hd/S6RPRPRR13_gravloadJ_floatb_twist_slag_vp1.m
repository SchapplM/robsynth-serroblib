% Calculate Gravitation load on the joints for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:44
% EndTime: 2019-03-09 04:21:46
% DurationCPUTime: 1.17s
% Computational Cost: add. (910->153), mult. (2425->227), div. (0->0), fcn. (3073->14), ass. (0->70)
t55 = cos(pkin(12));
t52 = sin(pkin(12));
t60 = sin(qJ(1));
t92 = t60 * t52;
t56 = cos(pkin(6));
t63 = cos(qJ(1));
t93 = t56 * t63;
t39 = -t55 * t93 + t92;
t90 = t60 * t55;
t40 = t52 * t93 + t90;
t59 = sin(qJ(3));
t86 = cos(pkin(7));
t81 = t59 * t86;
t54 = sin(pkin(6));
t94 = t54 * t63;
t53 = sin(pkin(7));
t96 = t53 * t59;
t98 = cos(qJ(3));
t19 = -t39 * t81 + t40 * t98 - t94 * t96;
t57 = sin(qJ(6));
t61 = cos(qJ(6));
t82 = t54 * t86;
t110 = -t39 * t53 + t63 * t82;
t77 = t86 * t98;
t84 = t53 * t98;
t80 = t54 * t84;
t18 = t39 * t77 + t40 * t59 + t63 * t80;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t7 = -t110 * t62 + t18 * t58;
t114 = -t19 * t61 + t57 * t7;
t113 = -t19 * t57 - t61 * t7;
t74 = t52 * t63 + t56 * t90;
t31 = t74 * t53 + t60 * t82;
t109 = t110 * t58 + t18 * t62;
t100 = rSges(6,3) + pkin(10);
t41 = t55 * t63 - t56 * t92;
t70 = t74 * t86;
t91 = t60 * t54;
t23 = t41 * t98 + (t53 * t91 - t70) * t59;
t26 = t56 * t96 + (t52 * t98 + t55 * t81) * t54;
t108 = g(1) * t23 + g(2) * t19 + g(3) * t26;
t22 = t41 * t59 - t60 * t80 + t70 * t98;
t95 = t54 * t55;
t25 = t52 * t54 * t59 - t56 * t84 - t77 * t95;
t107 = g(1) * t22 + g(2) * t18 + g(3) * t25;
t99 = rSges(7,3) + pkin(11);
t88 = qJ(2) * t54;
t89 = t63 * pkin(1) + t60 * t88;
t87 = rSges(5,3) + qJ(4);
t85 = -m(5) - m(6) - m(7);
t83 = -t60 * pkin(1) + t63 * t88;
t76 = rSges(7,1) * t61 - rSges(7,2) * t57 + pkin(5);
t71 = -t40 * pkin(2) + t110 * pkin(9) + t83;
t69 = -pkin(3) * t19 + t71;
t68 = t41 * pkin(2) + t31 * pkin(9) + t89;
t66 = pkin(4) * t110 - qJ(4) * t18 + t69;
t65 = t23 * pkin(3) + t68;
t64 = t31 * pkin(4) + t22 * qJ(4) + t65;
t38 = -t53 * t95 + t56 * t86;
t24 = t25 * pkin(3);
t17 = t25 * t58 + t38 * t62;
t16 = t25 * t62 - t38 * t58;
t14 = t22 * pkin(3);
t12 = t18 * pkin(3);
t9 = t22 * t58 + t31 * t62;
t8 = -t22 * t62 + t31 * t58;
t2 = t23 * t57 + t61 * t9;
t1 = t23 * t61 - t57 * t9;
t3 = [-m(2) * (g(1) * (-t60 * rSges(2,1) - rSges(2,2) * t63) + g(2) * (rSges(2,1) * t63 - t60 * rSges(2,2))) - m(3) * (g(1) * (-t40 * rSges(3,1) + t39 * rSges(3,2) + rSges(3,3) * t94 + t83) + g(2) * (t41 * rSges(3,1) - rSges(3,2) * t74 + rSges(3,3) * t91 + t89)) - m(4) * (g(1) * (-rSges(4,1) * t19 + rSges(4,2) * t18 + rSges(4,3) * t110 + t71) + g(2) * (t23 * rSges(4,1) - t22 * rSges(4,2) + t31 * rSges(4,3) + t68)) - m(5) * (g(1) * (rSges(5,1) * t110 + rSges(5,2) * t19 - t18 * t87 + t69) + g(2) * (t31 * rSges(5,1) - t23 * rSges(5,2) + t22 * t87 + t65)) - m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t109 - t100 * t19 + t66) + g(2) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t100 * t23 + t64)) - m(7) * (g(1) * (t113 * rSges(7,1) + t114 * rSges(7,2) - t7 * pkin(5) - t19 * pkin(10) + t99 * t109 + t66) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t9 * pkin(5) + t23 * pkin(10) + t8 * t99 + t64)) (-m(3) - m(4) + t85) * (g(3) * t56 + (g(1) * t60 - g(2) * t63) * t54) -m(4) * (g(1) * (-rSges(4,1) * t22 - rSges(4,2) * t23) + g(2) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(3) * (-rSges(4,1) * t25 - rSges(4,2) * t26)) - m(5) * (g(1) * (rSges(5,2) * t22 + t23 * t87 - t14) + g(2) * (rSges(5,2) * t18 + t19 * t87 - t12) + g(3) * (rSges(5,2) * t25 + t26 * t87 - t24)) + (-g(1) * (-t100 * t22 - t14) - g(2) * (-t100 * t18 - t12) - g(3) * (-t100 * t25 - t24) - t108 * (rSges(6,1) * t58 + rSges(6,2) * t62 + qJ(4))) * m(6) + (g(1) * t14 + g(2) * t12 + g(3) * t24 - t107 * (-rSges(7,1) * t57 - rSges(7,2) * t61 - pkin(10)) - t108 * (t58 * t76 - t62 * t99 + qJ(4))) * m(7), t85 * t107, -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t109 - rSges(6,2) * t7) + g(3) * (rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (-t76 * t8 + t9 * t99) + (t76 * t16 + t99 * t17) * g(3) + (t109 * t76 + t99 * t7) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t114 * rSges(7,1) + t113 * rSges(7,2)) + g(3) * ((-t17 * t57 + t26 * t61) * rSges(7,1) + (-t17 * t61 - t26 * t57) * rSges(7,2)))];
taug  = t3(:);
