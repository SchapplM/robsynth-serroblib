% Calculate Gravitation load on the joints for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:18
% EndTime: 2019-03-09 18:06:20
% DurationCPUTime: 0.88s
% Computational Cost: add. (585->128), mult. (520->143), div. (0->0), fcn. (452->12), ass. (0->79)
t49 = qJ(5) + qJ(6);
t43 = sin(t49);
t51 = sin(qJ(5));
t142 = -t51 * mrSges(6,2) - t43 * mrSges(7,2);
t140 = -mrSges(6,3) - mrSges(7,3);
t45 = cos(t49);
t54 = cos(qJ(5));
t141 = mrSges(6,1) * t54 + t45 * mrSges(7,1);
t50 = qJ(2) + qJ(3);
t42 = pkin(11) + t50;
t37 = sin(t42);
t139 = t142 * t37;
t121 = m(7) * pkin(5);
t38 = cos(t42);
t44 = sin(t50);
t46 = cos(t50);
t138 = mrSges(4,1) * t44 + mrSges(5,1) * t37 + mrSges(4,2) * t46 + mrSges(5,2) * t38;
t137 = t141 * t37;
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t136 = g(1) * t56 + g(2) * t53;
t135 = -t46 * mrSges(4,1) - t38 * mrSges(5,1) + t44 * mrSges(4,2) + (mrSges(5,2) + t140) * t37;
t134 = t51 * t121;
t77 = t38 * pkin(4) + t37 * pkin(9);
t116 = pkin(4) * t37;
t117 = pkin(3) * t44;
t40 = pkin(5) * t54 + pkin(4);
t57 = -pkin(10) - pkin(9);
t69 = -t37 * t40 - t38 * t57;
t133 = -m(7) * (t69 - t117) - m(6) * (-t116 - t117) + t137;
t107 = t38 * t56;
t131 = t107 * t140 + t139 * t56;
t108 = t38 * t53;
t130 = t108 * t140 + t139 * t53;
t129 = mrSges(6,1) + t121;
t70 = -t37 * t57 + t38 * t40;
t127 = m(5) + m(6) + m(7);
t58 = -pkin(8) - pkin(7);
t125 = -m(3) * pkin(7) + m(4) * t58 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t124 = t135 + (-t141 - t142) * t38;
t123 = m(6) * t116 - m(7) * t69 + t137;
t55 = cos(qJ(2));
t47 = t55 * pkin(2);
t52 = sin(qJ(2));
t76 = t55 * mrSges(3,1) - t52 * mrSges(3,2);
t122 = mrSges(2,1) + m(4) * (t47 + pkin(1)) + m(3) * pkin(1) + t76 - t135;
t100 = t45 * t56;
t105 = t43 * t53;
t5 = t105 * t38 + t100;
t101 = t45 * t53;
t104 = t43 * t56;
t6 = -t101 * t38 + t104;
t120 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t7 = -t104 * t38 + t101;
t8 = t100 * t38 + t105;
t119 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t118 = pkin(2) * t52;
t39 = pkin(3) * t46;
t112 = g(3) * t37;
t98 = t51 * t53;
t97 = t51 * t56;
t96 = t53 * t54;
t95 = t54 * t56;
t90 = t39 + t47;
t84 = t39 + t77;
t78 = t39 + t70;
t71 = -mrSges(7,1) * t43 - mrSges(7,2) * t45;
t11 = -t38 * t97 + t96;
t9 = t38 * t98 + t95;
t48 = -qJ(4) + t58;
t30 = pkin(9) * t107;
t29 = pkin(9) * t108;
t23 = -t117 - t118;
t22 = pkin(1) + t90;
t17 = t56 * t23;
t16 = t53 * t23;
t12 = t38 * t95 + t98;
t10 = -t38 * t96 + t97;
t1 = [(-t98 * t121 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t127 * (t56 * t22 - t48 * t53) + t125 * t53 + (-m(6) * t77 - m(7) * t70 - t122) * t56) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t127 * t48 + t125 - t134) * t56 + (m(5) * t22 - m(6) * (-t22 - t77) - m(7) * (-t22 - t70) + t122) * t53) * g(1) (-m(6) * (t16 + t29) - m(7) * t16 + t123 * t53 + t130) * g(2) + (-m(6) * (t17 + t30) - m(7) * t17 + t123 * t56 + t131) * g(1) + (-t76 - m(4) * t47 - m(5) * t90 - m(6) * (t47 + t84) - m(7) * (t47 + t78) + t124) * g(3) + t136 * (m(4) * t118 - m(5) * t23 + mrSges(3,1) * t52 + mrSges(3,2) * t55 + t138) (-m(6) * t29 + t133 * t53 + t130) * g(2) + (-m(6) * t30 + t133 * t56 + t131) * g(1) + (-m(5) * t39 - m(6) * t84 - m(7) * t78 + t124) * g(3) + t136 * (m(5) * t117 + t138) (-g(1) * t53 + g(2) * t56) * t127 (mrSges(6,1) * t51 + mrSges(6,2) * t54 + t134 - t71) * t112 + (-t10 * mrSges(6,2) + t129 * t9 - t120) * g(2) + (t12 * mrSges(6,2) - t129 * t11 - t119) * g(1), -g(1) * t119 - g(2) * t120 - t112 * t71];
taug  = t1(:);
