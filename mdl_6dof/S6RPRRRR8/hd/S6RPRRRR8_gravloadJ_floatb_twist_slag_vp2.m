% Calculate Gravitation load on the joints for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:46
% EndTime: 2019-03-09 07:19:48
% DurationCPUTime: 0.74s
% Computational Cost: add. (389->115), mult. (467->136), div. (0->0), fcn. (407->10), ass. (0->73)
t42 = qJ(5) + qJ(6);
t35 = sin(t42);
t44 = sin(qJ(5));
t122 = mrSges(6,2) * t44 + mrSges(7,2) * t35;
t121 = -mrSges(7,3) - mrSges(6,3);
t43 = qJ(3) + qJ(4);
t38 = cos(t43);
t120 = t122 * t38;
t47 = cos(qJ(5));
t33 = t47 * pkin(5) + pkin(4);
t36 = sin(t43);
t50 = -pkin(10) - pkin(9);
t61 = t36 * t33 + t38 * t50;
t69 = -t36 * pkin(4) + t38 * pkin(9);
t119 = -m(6) * t69 + m(7) * t61;
t118 = t36 * mrSges(5,1) + (mrSges(5,2) + t121) * t38;
t105 = m(7) * pkin(5);
t117 = -m(3) - m(4);
t114 = mrSges(6,1) + t105;
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t113 = -g(1) * t46 + g(2) * t49;
t112 = -m(5) - m(6) - m(7);
t95 = mrSges(6,1) * t47;
t72 = t38 * t95;
t37 = cos(t42);
t94 = mrSges(7,1) * t37;
t74 = t38 * t94;
t85 = t46 * t38;
t89 = t36 * t46;
t111 = -mrSges(5,1) * t85 + t121 * t89 + (-t74 - t72 + t120) * t46;
t88 = t36 * t50;
t93 = mrSges(5,2) * t36;
t110 = (-m(7) * t88 - t120 - t93) * t49;
t109 = t118 + (t94 + t95 - t122) * t36;
t108 = mrSges(5,1) * t38 - t121 * t36 + t72;
t107 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t64 = t45 * mrSges(4,1) + t48 * mrSges(4,2);
t106 = mrSges(2,2) - mrSges(3,3) - t64 - t118 - t119;
t81 = t49 * t37;
t87 = t46 * t35;
t5 = -t36 * t87 + t81;
t82 = t49 * t35;
t86 = t46 * t37;
t6 = t36 * t86 + t82;
t104 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t36 * t82 + t86;
t8 = t36 * t81 - t87;
t103 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t102 = pkin(3) * t48;
t99 = g(3) * t38;
t97 = t45 * pkin(3);
t84 = t46 * t44;
t83 = t46 * t47;
t80 = t49 * t44;
t79 = t49 * t47;
t78 = pkin(4) * t85 + pkin(9) * t89;
t77 = t49 * pkin(1) + t46 * qJ(2);
t76 = m(5) * t102;
t40 = t49 * qJ(2);
t70 = -t46 * pkin(1) + t40;
t67 = t33 * t85 - t46 * t88;
t66 = -pkin(4) * t38 - pkin(9) * t36;
t62 = -mrSges(7,1) * t35 - mrSges(7,2) * t37;
t11 = t36 * t80 + t83;
t9 = -t36 * t84 + t79;
t51 = -pkin(8) - pkin(7);
t28 = t46 * t102;
t12 = t36 * t79 - t84;
t10 = t36 * t83 + t80;
t1 = [(-t80 * t105 - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t117 * t77 + t112 * (t46 * t97 - t49 * t51 + t77) + (-m(4) * pkin(7) - t107) * t49 + t106 * t46) * g(2) + (t84 * t105 - m(3) * t70 - m(4) * t40 - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + t112 * (t46 * t51 + t49 * t97 + t70) + (-m(4) * (-pkin(1) - pkin(7)) + t107) * t46 + t106 * t49) * g(1), t113 * (-t112 - t117) t113 * (mrSges(4,1) * t48 - mrSges(4,2) * t45) + ((t76 - m(6) * (t66 - t102) - m(7) * (-t33 * t38 - t102) + t74 + t108) * t49 + t110) * g(2) + (-(t76 - t93) * t46 - m(6) * (t28 + t78) - m(7) * (t28 + t67) + t111) * g(1) + (t64 + m(5) * t97 - m(6) * (t69 - t97) - m(7) * (-t61 - t97) + t109) * g(3) (t109 + t119) * g(3) + ((-m(6) * t66 - (-m(7) * t33 - t94) * t38 + t108) * t49 + t110) * g(2) + (-m(6) * t78 - m(7) * t67 + mrSges(5,2) * t89 + t111) * g(1) (mrSges(6,2) * t47 + t114 * t44 - t62) * t99 + (-t12 * mrSges(6,2) - t114 * t11 - t103) * g(2) + (t10 * mrSges(6,2) - t114 * t9 - t104) * g(1), -g(1) * t104 - g(2) * t103 - t62 * t99];
taug  = t1(:);
