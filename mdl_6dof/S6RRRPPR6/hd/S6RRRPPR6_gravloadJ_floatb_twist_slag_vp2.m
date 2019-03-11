% Calculate Gravitation load on the joints for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:47
% EndTime: 2019-03-09 15:47:51
% DurationCPUTime: 1.45s
% Computational Cost: add. (727->143), mult. (1284->196), div. (0->0), fcn. (1466->12), ass. (0->71)
t131 = m(6) + m(7);
t56 = sin(qJ(6));
t60 = cos(qJ(6));
t135 = t56 * mrSges(7,1) + t60 * mrSges(7,2);
t128 = mrSges(5,2) - mrSges(6,3);
t73 = -t131 * qJ(5) + t128;
t121 = t73 - t135;
t64 = -m(4) * pkin(9) - m(7) * pkin(5) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t139 = t60 * mrSges(7,1) - t56 * mrSges(7,2) - t64;
t137 = -mrSges(5,1) + mrSges(6,2);
t53 = qJ(3) + pkin(11);
t50 = sin(t53);
t51 = cos(t53);
t57 = sin(qJ(3));
t61 = cos(qJ(3));
t136 = m(4) * pkin(2) + t61 * mrSges(4,1) - t57 * mrSges(4,2) - t128 * t50 - t137 * t51 + mrSges(3,1);
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t62 = cos(qJ(2));
t116 = cos(qJ(1));
t97 = cos(pkin(6));
t81 = t97 * t116;
t31 = t58 * t81 + t59 * t62;
t54 = sin(pkin(6));
t91 = t54 * t116;
t67 = t31 * t57 + t61 * t91;
t63 = t67 * pkin(3);
t107 = t54 * t58;
t132 = -t57 * t107 + t97 * t61;
t105 = t54 * t61;
t88 = t59 * t97;
t33 = t116 * t62 - t58 * t88;
t13 = t59 * t105 - t33 * t57;
t126 = m(5) + t131;
t125 = -m(7) * pkin(10) - mrSges(7,3);
t122 = -t125 - t137;
t119 = t51 * mrSges(7,3) + t135 * t50 + t136;
t30 = t58 * t59 - t62 * t81;
t115 = t30 * t51;
t32 = t116 * t58 + t62 * t88;
t113 = t32 * t51;
t106 = t54 * t59;
t104 = t54 * t62;
t103 = t56 * t62;
t102 = t60 * t62;
t49 = pkin(3) * t61 + pkin(2);
t55 = -qJ(4) - pkin(9);
t101 = -t30 * t49 - t31 * t55;
t100 = -t32 * t49 - t33 * t55;
t99 = t116 * pkin(1) + pkin(8) * t106;
t98 = qJ(5) * t50;
t96 = t57 * t106;
t90 = -pkin(1) * t59 + pkin(8) * t91;
t8 = t31 * t51 - t50 * t91;
t43 = t57 * t91;
t89 = -t31 * t61 + t43;
t86 = -pkin(4) * t115 - t30 * t98 + t101;
t85 = -pkin(4) * t113 - t32 * t98 + t100;
t83 = t13 * pkin(3);
t82 = pkin(3) * t96 - t32 * t55 + t33 * t49 + t99;
t75 = t132 * pkin(3);
t71 = pkin(3) * t43 + t30 * t55 - t31 * t49 + t90;
t7 = t31 * t50 + t51 * t91;
t36 = t49 * t104;
t24 = t107 * t50 - t51 * t97;
t14 = t33 * t61 + t96;
t12 = t106 * t50 + t33 * t51;
t11 = -t106 * t51 + t33 * t50;
t2 = t11 * t56 + t32 * t60;
t1 = t11 * t60 - t32 * t56;
t3 = [(-t116 * mrSges(2,1) - m(3) * t99 - t33 * mrSges(3,1) - m(4) * (pkin(2) * t33 + t99) - t14 * mrSges(4,1) - t13 * mrSges(4,2) - m(5) * t82 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t54 + mrSges(2,2)) * t59 + t73 * t11 - t122 * t12 + t64 * t32 - t131 * (t12 * pkin(4) + t82)) * g(2) + (t59 * mrSges(2,1) + t116 * mrSges(2,2) - m(3) * t90 + t31 * mrSges(3,1) - mrSges(3,3) * t91 - m(4) * (-pkin(2) * t31 + t90) - t89 * mrSges(4,1) - t67 * mrSges(4,2) - m(5) * t71 + t122 * t8 - t121 * t7 + t139 * t30 + t131 * (pkin(4) * t8 - t71)) * g(1) (-m(5) * t101 - m(6) * t86 - m(7) * (-pkin(10) * t115 + t86) - t139 * t31 + t119 * t30) * g(2) + (-m(5) * t100 - m(6) * t85 - m(7) * (-pkin(10) * t113 + t85) - t139 * t33 + t119 * t32) * g(1) + (-m(5) * t36 - t131 * (t36 + (pkin(4) * t51 + t98) * t104) + ((-t103 * mrSges(7,1) - t102 * mrSges(7,2)) * t50 + (t126 * t55 - t139) * t58 + (t125 * t51 - t136) * t62) * t54) * g(3) (-t132 * mrSges(4,1) - (-t105 * t58 - t57 * t97) * mrSges(4,2) - m(5) * t75 - t131 * (-t24 * pkin(4) + t75) + t121 * (t107 * t51 + t50 * t97) + t122 * t24) * g(3) + (m(5) * t63 + mrSges(4,1) * t67 - mrSges(4,2) * t89 + t121 * t8 + t122 * t7 + t131 * (t7 * pkin(4) + t63)) * g(2) + (-m(5) * t83 - t13 * mrSges(4,1) + t14 * mrSges(4,2) - t131 * (-t11 * pkin(4) + t83) + t121 * t12 + t122 * t11) * g(1), t126 * (-g(1) * t32 - g(2) * t30 + g(3) * t104) t131 * (-g(1) * t11 - g(2) * t7 - g(3) * t24) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t30 * t56 + t60 * t7) * mrSges(7,1) + (-t30 * t60 - t56 * t7) * mrSges(7,2)) - g(3) * ((t103 * t54 + t24 * t60) * mrSges(7,1) + (t102 * t54 - t24 * t56) * mrSges(7,2))];
taug  = t3(:);
