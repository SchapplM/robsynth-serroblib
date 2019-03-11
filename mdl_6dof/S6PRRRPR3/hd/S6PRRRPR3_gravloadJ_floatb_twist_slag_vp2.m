% Calculate Gravitation load on the joints for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:03
% EndTime: 2019-03-08 23:11:05
% DurationCPUTime: 0.95s
% Computational Cost: add. (646->124), mult. (1076->166), div. (0->0), fcn. (1230->12), ass. (0->71)
t141 = -mrSges(5,1) + mrSges(6,2);
t140 = mrSges(5,2) - mrSges(6,3);
t59 = sin(qJ(6));
t62 = cos(qJ(6));
t134 = mrSges(7,1) * t59 + mrSges(7,2) * t62;
t123 = -m(4) * pkin(8) - m(7) * pkin(5) - t62 * mrSges(7,1) + t59 * mrSges(7,2) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t56 = qJ(3) + qJ(4);
t54 = sin(t56);
t55 = cos(t56);
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t139 = -m(4) * pkin(2) - t63 * mrSges(4,1) + t60 * mrSges(4,2) + t140 * t54 + t141 * t55 - mrSges(3,1);
t138 = mrSges(7,3) - t141;
t137 = -t134 + t140;
t58 = sin(pkin(6));
t61 = sin(qJ(2));
t106 = t58 * t61;
t94 = cos(pkin(6));
t136 = -t60 * t106 + t94 * t63;
t105 = t58 * t63;
t64 = cos(qJ(2));
t57 = sin(pkin(11));
t85 = t57 * t94;
t93 = cos(pkin(11));
t42 = -t61 * t85 + t64 * t93;
t135 = t57 * t105 - t42 * t60;
t133 = m(6) + m(7);
t129 = m(5) + t133;
t35 = t106 * t54 - t55 * t94;
t36 = t106 * t55 + t54 * t94;
t128 = t137 * t36 + t138 * t35;
t107 = t57 * t58;
t19 = -t107 * t55 + t42 * t54;
t20 = t107 * t54 + t42 * t55;
t127 = t137 * t20 + t138 * t19;
t73 = t94 * t93;
t40 = t57 * t64 + t61 * t73;
t84 = t58 * t93;
t17 = t40 * t54 + t55 * t84;
t18 = t40 * t55 - t54 * t84;
t126 = t137 * t18 + t138 * t17;
t124 = t55 * mrSges(7,3) + t134 * t54 - t139;
t120 = pkin(10) * t19;
t119 = pkin(10) * t35;
t117 = t17 * pkin(10);
t39 = t57 * t61 - t64 * t73;
t114 = t39 * t55;
t41 = t61 * t93 + t64 * t85;
t113 = t41 * t55;
t104 = t58 * t64;
t103 = t59 * t64;
t102 = t62 * t64;
t53 = pkin(3) * t63 + pkin(2);
t65 = -pkin(9) - pkin(8);
t99 = -t39 * t53 - t40 * t65;
t98 = -t41 * t53 - t42 * t65;
t95 = qJ(5) * t54;
t82 = -t17 * pkin(4) + t18 * qJ(5);
t81 = -t19 * pkin(4) + qJ(5) * t20;
t80 = -t35 * pkin(4) + qJ(5) * t36;
t79 = -pkin(4) * t114 - t39 * t95 + t99;
t78 = -pkin(4) * t113 - t41 * t95 + t98;
t77 = t135 * pkin(3);
t74 = t136 * pkin(3);
t70 = -t40 * t60 - t63 * t84;
t69 = t77 + t81;
t68 = t74 + t80;
t67 = t70 * pkin(3);
t66 = t67 + t82;
t44 = t53 * t104;
t1 = [(-m(2) - m(3) - m(4) - t129) * g(3) (-m(5) * t99 - m(6) * t79 - m(7) * (-pkin(10) * t114 + t79) + t123 * t40 + t124 * t39) * g(2) + (-m(5) * t98 - m(6) * t78 - m(7) * (-pkin(10) * t113 + t78) + t123 * t42 + t124 * t41) * g(1) + (-m(5) * t44 - t133 * (t44 + (pkin(4) * t55 + t95) * t104) + ((-t103 * mrSges(7,1) - t102 * mrSges(7,2)) * t54 + ((-m(7) * pkin(10) - mrSges(7,3)) * t55 + t139) * t64 + (t129 * t65 + t123) * t61) * t58) * g(3) (-t136 * mrSges(4,1) - (-t105 * t61 - t60 * t94) * mrSges(4,2) - m(5) * t74 - m(6) * t68 - m(7) * (t68 - t119) + t128) * g(3) + (-t70 * mrSges(4,1) - (-t40 * t63 + t60 * t84) * mrSges(4,2) - m(5) * t67 - m(6) * t66 - m(7) * (t66 - t117) + t126) * g(2) + (-t135 * mrSges(4,1) - (-t107 * t60 - t42 * t63) * mrSges(4,2) - m(5) * t77 - m(6) * t69 - m(7) * (t69 - t120) + t127) * g(1) (-m(6) * t80 - m(7) * (t80 - t119) + t128) * g(3) + (-m(6) * t82 - m(7) * (t82 - t117) + t126) * g(2) + (-m(6) * t81 - m(7) * (t81 - t120) + t127) * g(1), t133 * (-g(1) * t19 - g(2) * t17 - g(3) * t35) -g(1) * ((t19 * t62 - t41 * t59) * mrSges(7,1) + (-t19 * t59 - t41 * t62) * mrSges(7,2)) - g(2) * ((t17 * t62 - t39 * t59) * mrSges(7,1) + (-t17 * t59 - t39 * t62) * mrSges(7,2)) - g(3) * ((t103 * t58 + t35 * t62) * mrSges(7,1) + (t102 * t58 - t35 * t59) * mrSges(7,2))];
taug  = t1(:);
