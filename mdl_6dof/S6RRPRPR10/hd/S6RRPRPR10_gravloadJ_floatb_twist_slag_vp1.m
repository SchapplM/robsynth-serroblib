% Calculate Gravitation load on the joints for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:31
% EndTime: 2019-03-09 11:04:33
% DurationCPUTime: 1.09s
% Computational Cost: add. (693->184), mult. (1167->257), div. (0->0), fcn. (1363->12), ass. (0->68)
t50 = pkin(11) + qJ(4);
t47 = sin(t50);
t48 = cos(t50);
t103 = pkin(4) * t48 + qJ(5) * t47;
t56 = sin(qJ(2));
t57 = sin(qJ(1));
t59 = cos(qJ(2));
t78 = cos(pkin(6));
t92 = cos(qJ(1));
t68 = t78 * t92;
t30 = t56 * t57 - t59 * t68;
t55 = sin(qJ(6));
t58 = cos(qJ(6));
t31 = t56 * t68 + t57 * t59;
t52 = sin(pkin(6));
t76 = t52 * t92;
t9 = t31 * t47 + t48 * t76;
t102 = -t30 * t58 - t55 * t9;
t101 = -t30 * t55 + t58 * t9;
t93 = pkin(10) + rSges(7,3);
t100 = t93 * t48;
t74 = t57 * t78;
t32 = t92 * t56 + t59 * t74;
t95 = g(2) * t30;
t99 = -g(1) * t32 - t95;
t98 = -m(6) - m(7);
t94 = g(3) * t52;
t89 = t47 * t55;
t88 = t47 * t58;
t87 = t52 * t56;
t86 = t52 * t57;
t85 = t52 * t59;
t53 = cos(pkin(11));
t46 = pkin(3) * t53 + pkin(2);
t54 = -pkin(9) - qJ(3);
t84 = -t30 * t46 - t31 * t54;
t33 = -t56 * t74 + t92 * t59;
t83 = -t32 * t46 - t33 * t54;
t82 = t92 * pkin(1) + pkin(8) * t86;
t80 = rSges(6,3) + qJ(5);
t79 = qJ(3) + rSges(4,3);
t51 = sin(pkin(11));
t77 = t51 * t86;
t75 = -t57 * pkin(1) + pkin(8) * t76;
t10 = t31 * t48 - t47 * t76;
t73 = -t103 * t30 + t84;
t72 = -t103 * t32 + t83;
t35 = t46 * t85;
t71 = g(3) * (t103 * t85 + t35);
t70 = t51 * t76;
t69 = pkin(3) * t77 - t32 * t54 + t33 * t46 + t82;
t67 = rSges(5,1) * t48 - rSges(5,2) * t47;
t66 = rSges(7,1) * t55 + rSges(7,2) * t58;
t65 = rSges(6,2) * t48 - rSges(6,3) * t47;
t14 = t33 * t48 + t47 * t86;
t64 = t14 * pkin(4) + t69;
t63 = rSges(4,1) * t53 - rSges(4,2) * t51 + pkin(2);
t61 = pkin(3) * t70 + t30 * t54 - t31 * t46 + t75;
t60 = -pkin(4) * t10 + t61;
t25 = t78 * t47 + t48 * t87;
t24 = t47 * t87 - t78 * t48;
t23 = t24 * pkin(4);
t13 = t33 * t47 - t48 * t86;
t7 = t13 * pkin(4);
t5 = t9 * pkin(4);
t3 = t13 * t55 + t32 * t58;
t2 = t13 * t58 - t32 * t55;
t1 = [-m(2) * (g(1) * (-t57 * rSges(2,1) - t92 * rSges(2,2)) + g(2) * (t92 * rSges(2,1) - t57 * rSges(2,2))) - m(3) * (g(1) * (-t31 * rSges(3,1) + t30 * rSges(3,2) + rSges(3,3) * t76 + t75) + g(2) * (rSges(3,1) * t33 - rSges(3,2) * t32 + rSges(3,3) * t86 + t82)) - m(4) * (g(1) * (-t31 * pkin(2) + (-t31 * t53 + t70) * rSges(4,1) + (t31 * t51 + t53 * t76) * rSges(4,2) - t79 * t30 + t75) + g(2) * (t33 * pkin(2) + (t33 * t53 + t77) * rSges(4,1) + (-t33 * t51 + t53 * t86) * rSges(4,2) + t79 * t32 + t82)) - m(5) * (g(1) * (-rSges(5,1) * t10 + rSges(5,2) * t9 - rSges(5,3) * t30 + t61) + g(2) * (rSges(5,1) * t14 - rSges(5,2) * t13 + rSges(5,3) * t32 + t69)) - m(6) * (g(1) * (-rSges(6,1) * t30 + rSges(6,2) * t10 - t80 * t9 + t60) + g(2) * (rSges(6,1) * t32 - rSges(6,2) * t14 + t80 * t13 + t64)) - m(7) * (g(1) * (t102 * rSges(7,1) - t101 * rSges(7,2) - t30 * pkin(5) - t9 * qJ(5) - t93 * t10 + t60) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t32 + qJ(5) * t13 + t93 * t14 + t64)) -m(3) * (g(1) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + g(2) * (-rSges(3,1) * t30 - rSges(3,2) * t31) + (rSges(3,1) * t59 - rSges(3,2) * t56) * t94) - m(4) * (g(1) * (-t63 * t32 + t79 * t33) + g(2) * t79 * t31 - t63 * t95 + (t79 * t56 + t63 * t59) * t94) - m(5) * (g(1) * (rSges(5,3) * t33 - t67 * t32 + t83) + g(2) * (rSges(5,3) * t31 - t67 * t30 + t84) + g(3) * t35 + (t67 * t59 + (rSges(5,3) - t54) * t56) * t94) - m(6) * (g(1) * (rSges(6,1) * t33 + t65 * t32 + t72) + g(2) * (rSges(6,1) * t31 + t65 * t30 + t73) + t71 + (-t65 * t59 + (rSges(6,1) - t54) * t56) * t94) - m(7) * (g(1) * (t33 * pkin(5) + (-t32 * t89 + t33 * t58) * rSges(7,1) + (-t32 * t88 - t33 * t55) * rSges(7,2) + t72) + g(2) * (t31 * pkin(5) + (-t30 * t89 + t31 * t58) * rSges(7,1) + (-t30 * t88 - t31 * t55) * rSges(7,2) + t73) + t71 + t99 * t100 + ((rSges(7,1) * t58 - rSges(7,2) * t55 + pkin(5) - t54) * t56 + (t66 * t47 + t100) * t59) * t94) (-m(4) - m(5) + t98) * (-g(3) * t85 - t99) -m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t14) + g(2) * (-rSges(5,1) * t9 - rSges(5,2) * t10) + g(3) * (-rSges(5,1) * t24 - rSges(5,2) * t25)) - m(6) * (g(1) * (rSges(6,2) * t13 + t80 * t14 - t7) + g(2) * (rSges(6,2) * t9 + t80 * t10 - t5) + g(3) * (rSges(6,2) * t24 + t80 * t25 - t23)) + (-g(1) * (-t93 * t13 - t7) - g(2) * (-t93 * t9 - t5) - g(3) * (-t93 * t24 - t23) - (g(1) * t14 + g(2) * t10 + g(3) * t25) * (qJ(5) + t66)) * m(7), t98 * (g(1) * t13 + g(2) * t9 + g(3) * t24) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (t101 * rSges(7,1) + t102 * rSges(7,2)) + g(3) * ((t24 * t58 + t55 * t85) * rSges(7,1) + (-t24 * t55 + t58 * t85) * rSges(7,2)))];
taug  = t1(:);
