% Calculate Gravitation load on the joints for
% S6RRRPRP8
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:27
% EndTime: 2019-03-09 17:14:31
% DurationCPUTime: 1.36s
% Computational Cost: add. (433->136), mult. (1008->163), div. (0->0), fcn. (1071->8), ass. (0->74)
t124 = mrSges(4,1) + mrSges(5,1);
t121 = mrSges(4,2) - mrSges(5,3);
t139 = -mrSges(6,3) - mrSges(7,3);
t44 = -qJ(6) - pkin(9);
t104 = m(6) * pkin(9) - m(7) * t44 - t139;
t136 = mrSges(5,2) + mrSges(4,3);
t138 = t104 - t136;
t120 = mrSges(6,2) + mrSges(7,2);
t46 = sin(qJ(3));
t50 = cos(qJ(3));
t49 = cos(qJ(5));
t92 = t46 * t49;
t45 = sin(qJ(5));
t94 = t45 * t50;
t63 = -t92 + t94;
t137 = -t120 * t63 - t121 * t46 + t124 * t50;
t123 = mrSges(6,1) + mrSges(7,1);
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t117 = g(1) * t52 + g(2) * t48;
t47 = sin(qJ(2));
t132 = t117 * t47;
t51 = cos(qJ(2));
t84 = t52 * t46;
t26 = -t48 * t50 + t51 * t84;
t85 = t51 * t52;
t27 = t48 * t46 + t50 * t85;
t64 = t26 * t45 + t27 * t49;
t130 = t120 * t64;
t128 = m(6) * pkin(4);
t103 = m(7) * pkin(5);
t127 = -m(5) - m(6);
t126 = -m(6) - m(7);
t95 = t45 * t46;
t62 = t49 * t50 + t95;
t125 = t62 * t47;
t122 = mrSges(2,2) - mrSges(3,3);
t89 = t47 * t50;
t119 = -t123 * (t45 * t89 - t47 * t92) - t120 * t125;
t87 = t48 * t51;
t24 = t46 * t87 + t50 * t52;
t86 = t50 * t51;
t25 = t48 * t86 - t84;
t118 = t24 * t49 - t25 * t45;
t116 = m(5) - t126;
t70 = t51 * mrSges(3,1) - t47 * mrSges(3,2);
t115 = t136 * t47 + t70;
t114 = -t103 - t123;
t97 = t24 * t45;
t1 = t25 * t49 + t97;
t112 = t120 * t1;
t78 = pkin(5) * t45 + qJ(4);
t107 = -m(7) * t78 + t121;
t38 = pkin(5) * t49 + pkin(4);
t105 = m(7) * t38 + t124 + t128;
t102 = -pkin(3) - pkin(4);
t101 = m(7) * t47;
t39 = t47 * pkin(8);
t41 = t51 * pkin(2);
t98 = -pkin(3) - t38;
t93 = t46 * t47;
t88 = t47 * t52;
t83 = t41 + t39;
t82 = t52 * pkin(1) + t48 * pkin(7);
t81 = qJ(4) * t46;
t80 = -pkin(1) - t41;
t79 = -pkin(2) - t81;
t73 = pkin(2) * t85 + pkin(8) * t88 + t82;
t71 = t27 * pkin(3) + t73;
t9 = t26 * t49 - t27 * t45;
t42 = t52 * pkin(7);
t22 = t26 * pkin(3);
t20 = t24 * pkin(3);
t2 = [(-m(3) * t82 - m(4) * t73 - m(7) * t71 - t120 * t9 + t127 * (t26 * qJ(4) + t71) + (-mrSges(2,1) - t70) * t52 + t122 * t48 - t105 * t27 + t107 * t26 - t123 * t64 + t138 * t88) * g(2) + (t97 * t103 - t116 * (-t25 * pkin(3) - qJ(4) * t24 + t42) + t123 * t1 + t122 * t52 + t120 * t118 + (-m(3) - m(4)) * t42 + t105 * t25 - t121 * t24 + (m(3) * pkin(1) + mrSges(2,1) + t126 * t80 + (-m(4) - m(5)) * (t80 - t39) + (-m(6) * (-pkin(8) + pkin(9)) - m(7) * (-pkin(8) - t44) + t139) * t47 + t115) * t48) * g(1), mrSges(3,1) * t132 + (-t86 * t128 - m(4) * t83 - t116 * (pkin(3) * t86 + t83) + t104 * t47 - t115) * g(3) + (t85 * g(1) + t87 * g(2)) * pkin(8) * (-m(4) - t116) + (t117 * mrSges(3,2) + (-t116 * t81 - m(7) * (pkin(5) * t95 + t38 * t50) - t123 * t62 - t137) * g(3)) * t51 + t117 * (-(-t78 * t46 + t98 * t50 - pkin(2)) * t101 + t138 * t51 + t125 * t123 + (-m(5) * (-pkin(3) * t50 + t79) - m(6) * (t102 * t50 + t79) + m(4) * pkin(2) + t137) * t47) (-m(6) * t102 * t93 - (pkin(5) * t94 + t98 * t46) * t101 + (t121 * t50 + (m(5) * pkin(3) + t124) * t46) * t47 - t116 * qJ(4) * t89 + t119) * g(3) + (m(7) * t20 + t127 * (qJ(4) * t25 - t20) + t107 * t25 + t105 * t24 + t123 * t118 - t112) * g(2) + (m(7) * t22 + t123 * t9 + t127 * (qJ(4) * t27 - t22) - t130 + t107 * t27 + t105 * t26) * g(1), t116 * (-g(1) * t26 - g(2) * t24 - g(3) * t93) (t63 * pkin(5) * t101 - t119) * g(3) + (t114 * t118 + t112) * g(2) + (t114 * t9 + t130) * g(1) (-g(3) * t51 + t132) * m(7)];
taug  = t2(:);
