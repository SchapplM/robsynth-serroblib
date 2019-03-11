% Calculate Gravitation load on the joints for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:52
% EndTime: 2019-03-09 22:38:56
% DurationCPUTime: 1.43s
% Computational Cost: add. (672->160), mult. (926->183), div. (0->0), fcn. (958->10), ass. (0->81)
t148 = -mrSges(3,2) + mrSges(5,3) + mrSges(6,2);
t134 = mrSges(5,1) + mrSges(6,1);
t132 = mrSges(5,2) - mrSges(6,3);
t55 = sin(qJ(2));
t59 = cos(qJ(2));
t129 = t59 * mrSges(3,1) + t148 * t55;
t147 = t55 * mrSges(4,3) + mrSges(2,1) + t129;
t52 = qJ(3) + qJ(4);
t47 = sin(t52);
t48 = cos(t52);
t145 = pkin(4) * t48 + qJ(5) * t47;
t54 = sin(qJ(3));
t58 = cos(qJ(3));
t53 = sin(qJ(6));
t57 = cos(qJ(6));
t76 = t47 * t53 + t48 * t57;
t77 = t47 * t57 - t48 * t53;
t144 = m(4) * pkin(2) + t58 * mrSges(4,1) + t76 * mrSges(7,1) - t54 * mrSges(4,2) + t77 * mrSges(7,2) - t132 * t47 + t134 * t48;
t56 = sin(qJ(1));
t60 = cos(qJ(1));
t143 = g(1) * t60 + g(2) * t56;
t112 = t48 * t55;
t122 = (-mrSges(7,1) * t77 + mrSges(7,2) * t76) * t55;
t141 = -mrSges(6,3) * t112 - t122;
t102 = t59 * t60;
t32 = -t54 * t102 + t56 * t58;
t103 = t56 * t59;
t26 = t103 * t47 + t48 * t60;
t101 = t60 * t47;
t27 = t103 * t48 - t101;
t78 = t26 * t53 + t27 * t57;
t89 = -t26 * t57 + t27 * t53;
t124 = t89 * mrSges(7,1) + t78 * mrSges(7,2);
t28 = t101 * t59 - t56 * t48;
t29 = t102 * t48 + t47 * t56;
t5 = t28 * t57 - t29 * t53;
t6 = t28 * t53 + t29 * t57;
t131 = mrSges(7,1) * t5 - mrSges(7,2) * t6;
t140 = t132 * t29 + t134 * t28 + t131;
t139 = t132 * t27 + t134 * t26 - t124;
t138 = -m(3) - m(4);
t137 = m(6) + m(7);
t133 = mrSges(2,2) - mrSges(3,3);
t130 = t145 * t59;
t128 = m(7) * pkin(5) + t134;
t125 = -pkin(4) - pkin(5);
t120 = pkin(3) * t54;
t118 = pkin(5) * t48;
t115 = g(3) * t55;
t114 = mrSges(5,2) * t48;
t111 = t54 * t56;
t110 = t54 * t60;
t61 = -pkin(9) - pkin(8);
t105 = t55 * t61;
t46 = pkin(3) * t58 + pkin(2);
t37 = t59 * t46;
t96 = t60 * pkin(1) + t56 * pkin(7);
t93 = t60 * t105;
t50 = t60 * pkin(7);
t92 = pkin(3) * t110 + t56 * t105 + t50;
t91 = m(4) * pkin(8) + mrSges(4,3);
t88 = t37 - t105;
t87 = -t26 * pkin(4) + qJ(5) * t27;
t86 = -t28 * pkin(4) + qJ(5) * t29;
t85 = pkin(3) * t111 + t46 * t102 + t96;
t84 = m(7) * (-pkin(10) - t61) - mrSges(7,3);
t83 = pkin(2) * t59 + pkin(8) * t55;
t75 = t32 * pkin(3);
t30 = t103 * t54 + t58 * t60;
t74 = -t46 - t145;
t73 = t84 * t55;
t70 = t30 * pkin(3);
t69 = t29 * pkin(4) + t28 * qJ(5) + t85;
t67 = t75 + t86;
t66 = -t70 + t87;
t35 = qJ(5) * t112;
t33 = t102 * t58 + t111;
t31 = -t103 * t58 + t110;
t23 = t28 * pkin(5);
t20 = t26 * pkin(5);
t1 = [(-t33 * mrSges(4,1) - t32 * mrSges(4,2) - m(5) * (t85 - t93) - m(6) * (t69 - t93) - m(7) * t69 - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t138 * t96 + t133 * t56 - t128 * t29 + t132 * t28 + (-m(4) * t83 - t147 - t73) * t60) * g(2) + (-m(5) * t92 - t31 * mrSges(4,1) + t78 * mrSges(7,1) - t30 * mrSges(4,2) - t89 * mrSges(7,2) - t137 * (-t27 * pkin(4) - qJ(5) * t26 + t92) + t133 * t60 + t138 * t50 + t128 * t27 - t132 * t26 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t83) + (-m(7) * pkin(10) - mrSges(7,3)) * t55 + (-m(5) - t137) * (-pkin(1) - t37) + t147) * t56) * g(1) (-m(5) * t88 - m(6) * (t88 + t130) - m(7) * (t37 + t130) - t73 - t129) * g(3) + ((-m(7) * t118 - t144) * g(3) + t143 * (-t84 - t91 + (m(5) + m(6)) * t61 - t148)) * t59 + (-t91 * g(3) + t143 * (mrSges(3,1) + m(5) * t46 - m(6) * t74 - m(7) * (t74 - t118) + t144)) * t55 (m(5) * t120 + mrSges(4,1) * t54 + mrSges(5,1) * t47 + mrSges(4,2) * t58 + t114) * t115 + (-m(6) * t35 - (m(6) * (-pkin(4) * t47 - t120) - t47 * mrSges(6,1)) * t55 - (t35 + (t125 * t47 - t120) * t55) * m(7) + t141) * g(3) + (t30 * mrSges(4,1) - t31 * mrSges(4,2) + m(5) * t70 - m(6) * t66 - m(7) * (-t20 + t66) + t139) * g(2) + (-t32 * mrSges(4,1) + t33 * mrSges(4,2) - m(5) * t75 - m(6) * t67 - m(7) * (-t23 + t67) + t140) * g(1) (-t137 * t35 + (t114 + (m(6) * pkin(4) - m(7) * t125 + t134) * t47) * t55 + t141) * g(3) + (-m(6) * t87 - m(7) * (-t20 + t87) + t139) * g(2) + (-m(6) * t86 - m(7) * (-t23 + t86) + t140) * g(1), t137 * (-g(1) * t28 - g(2) * t26 - t115 * t47) -g(1) * t131 + g(2) * t124 + g(3) * t122];
taug  = t1(:);
