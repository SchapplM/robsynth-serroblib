% Calculate Gravitation load on the joints for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:34
% EndTime: 2019-03-09 00:37:37
% DurationCPUTime: 0.96s
% Computational Cost: add. (832->126), mult. (1078->170), div. (0->0), fcn. (1220->14), ass. (0->64)
t137 = mrSges(6,2) - mrSges(7,3);
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t136 = -t65 * mrSges(7,1) + t62 * mrSges(7,2) - mrSges(6,1);
t61 = sin(pkin(6));
t64 = sin(qJ(2));
t112 = t61 * t64;
t59 = qJ(3) + qJ(4);
t54 = sin(t59);
t55 = cos(t59);
t99 = cos(pkin(6));
t135 = -t54 * t112 + t99 * t55;
t60 = sin(pkin(12));
t113 = t60 * t61;
t67 = cos(qJ(2));
t89 = t60 * t99;
t98 = cos(pkin(12));
t40 = -t64 * t89 + t98 * t67;
t134 = t55 * t113 - t40 * t54;
t68 = -pkin(9) - pkin(8);
t119 = -m(4) * pkin(8) + m(5) * t68 - t62 * mrSges(7,1) - t65 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t56 = qJ(5) + t59;
t51 = sin(t56);
t52 = cos(t56);
t66 = cos(qJ(3));
t57 = t66 * pkin(3);
t63 = sin(qJ(3));
t132 = -mrSges(3,1) - m(5) * (t57 + pkin(2)) - t55 * mrSges(5,1) + t54 * mrSges(5,2) - m(4) * pkin(2) - t66 * mrSges(4,1) + t63 * mrSges(4,2) + (-m(7) * pkin(5) + t136) * t52 + (-m(7) * pkin(11) + t137) * t51;
t131 = m(6) + m(7);
t129 = -m(5) * pkin(3) - mrSges(4,1);
t28 = -t51 * t112 + t99 * t52;
t29 = t52 * t112 + t99 * t51;
t128 = t136 * t28 + t137 * t29;
t13 = t52 * t113 - t40 * t51;
t14 = t51 * t113 + t40 * t52;
t127 = t136 * t13 + t137 * t14;
t80 = t99 * t98;
t38 = t60 * t67 + t64 * t80;
t88 = t61 * t98;
t11 = -t38 * t51 - t52 * t88;
t12 = t38 * t52 - t51 * t88;
t126 = t136 * t11 + t137 * t12;
t123 = -t135 * mrSges(5,1) - (-t55 * t112 - t99 * t54) * mrSges(5,2) + t128;
t122 = -t134 * mrSges(5,1) - (-t54 * t113 - t40 * t55) * mrSges(5,2) + t127;
t73 = -t38 * t54 - t55 * t88;
t121 = -t73 * mrSges(5,1) - (-t38 * t55 + t54 * t88) * mrSges(5,2) + t126;
t111 = t61 * t66;
t110 = t61 * t67;
t45 = -pkin(3) * t63 - pkin(4) * t54;
t46 = pkin(4) * t55 + t57;
t103 = t46 * t113 + t40 * t45;
t100 = t45 * t112 + t99 * t46;
t93 = t11 * pkin(5) + t12 * pkin(11);
t91 = t13 * pkin(5) + pkin(11) * t14;
t90 = t28 * pkin(5) + pkin(11) * t29;
t86 = t134 * pkin(4);
t81 = t135 * pkin(4);
t77 = t38 * t45 - t46 * t88;
t70 = t73 * pkin(4);
t58 = -pkin(10) + t68;
t44 = pkin(2) + t46;
t39 = t98 * t64 + t67 * t89;
t37 = t60 * t64 - t67 * t80;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t131) * g(3) (-t131 * (-t37 * t44 - t38 * t58) + t119 * t38 - t132 * t37) * g(2) + (-t131 * (-t39 * t44 - t40 * t58) + t119 * t40 - t132 * t39) * g(1) + (-t131 * t44 * t110 + (t132 * t67 + (t131 * t58 + t119) * t64) * t61) * g(3) (-(-t64 * t111 - t99 * t63) * mrSges(4,2) - m(6) * t100 - m(7) * (t90 + t100) + t129 * (-t63 * t112 + t99 * t66) + t123) * g(3) + (-(-t38 * t66 + t63 * t88) * mrSges(4,2) - m(6) * t77 - m(7) * (t77 + t93) + t129 * (-t38 * t63 - t66 * t88) + t121) * g(2) + (-(-t63 * t113 - t40 * t66) * mrSges(4,2) - m(6) * t103 - m(7) * (t91 + t103) + t129 * (t60 * t111 - t40 * t63) + t122) * g(1) (-m(6) * t81 - m(7) * (t81 + t90) + t123) * g(3) + (-m(6) * t70 - m(7) * (t70 + t93) + t121) * g(2) + (-m(6) * t86 - m(7) * (t86 + t91) + t122) * g(1) (-m(7) * t90 + t128) * g(3) + (-m(7) * t93 + t126) * g(2) + (-m(7) * t91 + t127) * g(1), -g(1) * ((-t14 * t62 + t39 * t65) * mrSges(7,1) + (-t14 * t65 - t39 * t62) * mrSges(7,2)) - g(2) * ((-t12 * t62 + t37 * t65) * mrSges(7,1) + (-t12 * t65 - t37 * t62) * mrSges(7,2)) - g(3) * ((-t65 * t110 - t29 * t62) * mrSges(7,1) + (t62 * t110 - t29 * t65) * mrSges(7,2))];
taug  = t1(:);
