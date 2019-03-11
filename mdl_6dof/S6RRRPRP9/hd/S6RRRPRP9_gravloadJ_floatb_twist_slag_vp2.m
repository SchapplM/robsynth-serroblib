% Calculate Gravitation load on the joints for
% S6RRRPRP9
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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:53
% EndTime: 2019-03-09 17:21:57
% DurationCPUTime: 1.27s
% Computational Cost: add. (473->119), mult. (1124->148), div. (0->0), fcn. (1232->8), ass. (0->56)
t122 = mrSges(4,1) + mrSges(5,1);
t120 = mrSges(4,2) - mrSges(5,3);
t54 = sin(qJ(5));
t55 = sin(qJ(3));
t59 = cos(qJ(3));
t87 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t58 = cos(qJ(5));
t99 = t58 * t59;
t137 = t87 * (t54 * t55 + t99) + t122 * t59 - t120 * t55;
t136 = m(6) + m(7);
t134 = mrSges(4,3) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t56 = sin(qJ(2));
t108 = t55 * t56;
t102 = t56 * t59;
t107 = t55 * t58;
t19 = t54 * t102 - t56 * t107;
t86 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t132 = t87 * t19 + t86 * (-t54 * t108 - t56 * t99);
t57 = sin(qJ(1));
t60 = cos(qJ(2));
t61 = cos(qJ(1));
t96 = t61 * t55;
t32 = -t57 * t59 + t60 * t96;
t97 = t60 * t61;
t33 = t57 * t55 + t59 * t97;
t10 = t32 * t54 + t33 * t58;
t77 = -t32 * t58 + t33 * t54;
t131 = -t86 * t10 + t77 * t87;
t100 = t57 * t60;
t30 = t55 * t100 + t59 * t61;
t98 = t59 * t60;
t31 = t57 * t98 - t96;
t119 = -t30 * t58 + t31 * t54;
t78 = t30 * t54 + t31 * t58;
t129 = t87 * t119 - t86 * t78;
t121 = mrSges(2,2) - mrSges(3,3);
t95 = t60 * pkin(2) + t56 * pkin(8);
t118 = m(5) + t136;
t117 = -m(4) - t118;
t82 = t60 * mrSges(3,1) - t56 * mrSges(3,2);
t114 = t134 * t56 + t82;
t111 = -pkin(3) - pkin(4);
t110 = pkin(9) * t56;
t101 = t56 * t61;
t94 = t61 * pkin(1) + t57 * pkin(7);
t93 = qJ(4) * t55;
t92 = -pkin(2) - t93;
t91 = -t30 * pkin(3) + qJ(4) * t31;
t90 = -t32 * pkin(3) + qJ(4) * t33;
t89 = pkin(3) * t98 + t60 * t93 + t95;
t88 = pkin(2) * t97 + pkin(8) * t101 + t94;
t52 = t61 * pkin(7);
t84 = -t31 * pkin(3) - qJ(4) * t30 + t52;
t69 = t33 * pkin(3) + t32 * qJ(4) + t88;
t38 = qJ(4) * t102;
t1 = [(-m(3) * t94 - m(4) * t88 - m(5) * t69 - t136 * (t33 * pkin(4) - pkin(9) * t101 + t69) - t86 * t77 + (-mrSges(2,1) - t82) * t61 + t121 * t57 - t122 * t33 + t120 * t32 - t87 * t10 - t134 * t101) * g(2) + (-m(5) * t84 - t136 * (-t31 * pkin(4) + t57 * t110 + t84) + t121 * t61 + t87 * t78 + (-m(3) - m(4)) * t52 + t86 * t119 + t122 * t31 - t120 * t30 + (m(3) * pkin(1) + mrSges(2,1) + t117 * (-pkin(1) - t95) + t114) * t57) * g(1) (-m(4) * t95 - m(5) * t89 - t136 * (pkin(4) * t98 - t110 + t89) - t86 * t54 * t98 - t114 + (t86 * t107 - t137) * t60) * g(3) + (t97 * g(1) + t100 * g(2)) * pkin(8) * t117 + (t61 * g(1) + t57 * g(2)) * ((t136 * pkin(9) + mrSges(3,2) - t134) * t60 + ((t54 * t59 - t107) * t86 - t136 * (t111 * t59 + t92) + mrSges(3,1) + m(4) * pkin(2) - m(5) * (-pkin(3) * t59 + t92) + t137) * t56) (-m(5) * t38 - t136 * (t111 * t108 + t38) + (t120 * t59 + (m(5) * pkin(3) + t122) * t55) * t56 - t132) * g(3) + (-m(5) * t91 - t136 * (-t30 * pkin(4) + t91) + t120 * t31 + t122 * t30 - t129) * g(2) + (-m(5) * t90 - t136 * (-t32 * pkin(4) + t90) + t120 * t33 + t122 * t32 - t131) * g(1), t118 * (-g(1) * t32 - g(2) * t30 - g(3) * t108) t131 * g(1) + t129 * g(2) + t132 * g(3) (-g(1) * t77 - g(2) * t119 - g(3) * t19) * m(7)];
taug  = t1(:);
