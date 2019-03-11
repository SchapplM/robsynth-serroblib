% Calculate Gravitation load on the joints for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:54:00
% EndTime: 2019-03-09 20:54:03
% DurationCPUTime: 1.11s
% Computational Cost: add. (562->130), mult. (740->150), div. (0->0), fcn. (696->8), ass. (0->66)
t125 = -mrSges(6,3) - mrSges(7,2);
t134 = mrSges(6,1) + mrSges(5,3);
t48 = sin(qJ(4));
t133 = t125 * t48;
t126 = m(6) + m(7);
t132 = m(5) + t126;
t131 = -mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t47 = qJ(2) + qJ(3);
t44 = sin(t47);
t51 = cos(qJ(4));
t76 = m(7) * (-pkin(4) - qJ(6)) - mrSges(7,3);
t91 = qJ(5) * t48;
t81 = -pkin(3) - t91;
t124 = (-m(6) * (-pkin(4) * t51 + t81) - m(7) * t81 - t51 * t76 - t133) * t44;
t130 = m(5) * pkin(3) * t44 + t124;
t129 = -m(7) * pkin(5) - mrSges(7,1) - t134;
t128 = t134 * t44;
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t127 = g(1) * t53 + g(2) * t50;
t45 = cos(t47);
t122 = t45 * mrSges(4,1) - t44 * mrSges(4,2);
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t73 = t52 * mrSges(3,1) - t49 * mrSges(3,2);
t121 = -m(3) * pkin(1) - mrSges(2,1) - t122 - t73;
t120 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t119 = mrSges(5,2) + t125;
t106 = t44 * t51;
t113 = pkin(2) * t49;
t118 = mrSges(5,1) * t106 + t113 * t132 + t130;
t102 = t45 * t53;
t31 = pkin(9) * t102;
t96 = t51 * t53;
t88 = t44 * t96;
t101 = t48 * mrSges(5,2);
t90 = t44 * t101;
t117 = -mrSges(6,2) * t88 + t129 * t102 - t126 * t31 - t53 * t90;
t103 = t45 * t51;
t116 = -t44 * mrSges(7,1) - t122 - t128 + (t101 + t133) * t45 + t131 * t103;
t97 = t50 * t51;
t89 = t44 * t97;
t115 = -mrSges(6,2) * t89 + (-t90 + (-pkin(9) * t132 + t129) * t45) * t50;
t114 = m(7) * qJ(6) - t131;
t109 = g(3) * t44;
t38 = t44 * pkin(9);
t40 = t45 * pkin(3);
t46 = t52 * pkin(2);
t105 = t44 * t53;
t98 = t50 * t48;
t95 = t53 * t48;
t54 = -pkin(8) - pkin(7);
t94 = t53 * t54;
t92 = t40 + t38;
t43 = t46 + pkin(1);
t85 = -t43 - t40;
t80 = t53 * t43 - t50 * t54;
t79 = pkin(4) * t103 + t45 * t91 + t92;
t70 = mrSges(4,1) * t44 + mrSges(4,2) * t45;
t67 = pkin(3) * t102 + pkin(9) * t105 + t80;
t66 = t44 * pkin(5) + qJ(6) * t103 + t79;
t8 = t45 * t96 + t98;
t7 = t45 * t95 - t97;
t6 = t45 * t97 - t95;
t5 = t45 * t98 + t96;
t1 = [(-m(4) * t80 - m(5) * t67 - t126 * (t8 * pkin(4) + t7 * qJ(5) + t67) - t114 * t8 + t119 * t7 + t121 * t53 + t120 * t50 + t129 * t105) * g(2) + (m(5) * t94 - t126 * (-t6 * pkin(4) - t5 * qJ(5) - t94) + t114 * t6 + (m(4) * t54 + t120) * t53 - t119 * t5 + (m(4) * t43 - m(7) * t85 - (m(7) * (-pkin(5) - pkin(9)) - mrSges(7,1)) * t44 + (-m(5) - m(6)) * (t85 - t38) - t121 + t128) * t50) * g(1) (t118 * t50 + t115) * g(2) + (-m(5) * t31 + t118 * t53 + t117) * g(1) + (-t73 - m(4) * t46 - m(5) * (t46 + t92) - m(6) * (t46 + t79) - m(7) * (t46 + t66) + t116) * g(3) + t127 * (m(4) * t113 + mrSges(3,1) * t49 + mrSges(3,2) * t52 + t70) t127 * t70 + (mrSges(5,1) * t89 + t130 * t50 + t115) * g(2) + (-m(5) * (-pkin(3) * t105 + t31) + mrSges(5,1) * t88 + t124 * t53 + t117) * g(1) + (-m(5) * t92 - m(6) * t79 - m(7) * t66 + t116) * g(3) -(-mrSges(5,1) * t48 - mrSges(5,2) * t51) * t109 + ((t125 * t51 + (m(6) * pkin(4) - mrSges(6,2) - t76) * t48) * t44 - t126 * qJ(5) * t106) * g(3) + (-t126 * (-t5 * pkin(4) + qJ(5) * t6) + t119 * t6 + t114 * t5) * g(2) + (-t126 * (-t7 * pkin(4) + qJ(5) * t8) + t119 * t8 + t114 * t7) * g(1), t126 * (-g(1) * t7 - g(2) * t5 - t109 * t48) (-g(1) * t8 - g(2) * t6 - g(3) * t106) * m(7)];
taug  = t1(:);
