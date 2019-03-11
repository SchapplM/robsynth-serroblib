% Calculate potential energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:47
% EndTime: 2019-03-08 20:13:47
% DurationCPUTime: 0.66s
% Computational Cost: add. (260->101), mult. (545->118), div. (0->0), fcn. (622->10), ass. (0->48)
t144 = -mrSges(3,1) + mrSges(4,2);
t143 = -mrSges(6,1) - mrSges(7,1);
t142 = -mrSges(6,2) - mrSges(7,2);
t141 = mrSges(3,3) + mrSges(4,1);
t140 = -mrSges(4,3) + mrSges(3,2);
t110 = sin(qJ(5));
t139 = -m(7) * (pkin(5) * t110 + pkin(8)) - mrSges(5,3) + t144;
t138 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t107 = cos(pkin(10));
t105 = sin(pkin(10));
t106 = sin(pkin(6));
t136 = t105 * t106;
t137 = t107 * pkin(1) + pkin(7) * t136;
t135 = t106 * t107;
t111 = sin(qJ(4));
t134 = t106 * t111;
t112 = sin(qJ(2));
t133 = t106 * t112;
t114 = cos(qJ(4));
t132 = t106 * t114;
t115 = cos(qJ(2));
t131 = t106 * t115;
t108 = cos(pkin(6));
t130 = t108 * t112;
t129 = t108 * t115;
t128 = t110 * t112;
t127 = t108 * pkin(7) + qJ(1);
t126 = pkin(7) * t135;
t125 = pkin(2) * t133 + t127;
t101 = t105 * pkin(1);
t88 = t105 * t112 - t107 * t129;
t89 = t105 * t115 + t107 * t130;
t123 = t89 * pkin(2) + qJ(3) * t88 + t101;
t122 = t108 * pkin(3) + pkin(8) * t133 + t125;
t90 = t105 * t129 + t107 * t112;
t91 = -t105 * t130 + t107 * t115;
t121 = t91 * pkin(2) + qJ(3) * t90 + t137;
t120 = pkin(3) * t136 + t121;
t119 = pkin(8) * t91 + t120;
t118 = -qJ(3) * t131 + t122;
t117 = (-pkin(3) - pkin(7)) * t135 + t123;
t116 = pkin(8) * t89 + t117;
t113 = cos(qJ(5));
t100 = pkin(5) * t113 + pkin(4);
t93 = t108 * t114 - t111 * t131;
t82 = -t107 * t132 + t111 * t88;
t80 = t105 * t132 + t111 * t90;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t127 - m(4) * t125 - m(5) * t118 - t93 * mrSges(5,1) - mrSges(5,3) * t133 - m(6) * (t93 * pkin(4) + t118) - m(7) * (t93 * t100 + t122) + t143 * (t106 * t128 + t113 * t93) + t142 * (-t110 * t93 + t113 * t133) - t141 * t108 + (-m(7) * pkin(5) * t128 + ((m(4) + m(7)) * qJ(3) - t140) * t115 + t144 * t112) * t106 - t138 * (t108 * t111 + t114 * t131)) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t105 - mrSges(2,2) * t107 - m(3) * (t101 - t126) - m(4) * (t123 - t126) - m(5) * t116 - t82 * mrSges(5,1) - m(6) * (pkin(4) * t82 + t116) - m(7) * (t100 * t82 + t117) + t140 * t88 + t143 * (t110 * t89 + t113 * t82) + t142 * (-t110 * t82 + t113 * t89) + t141 * t135 + t139 * t89 + t138 * (t107 * t134 + t114 * t88)) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t107 + mrSges(2,2) * t105 - m(3) * t137 - m(4) * t121 - m(5) * t119 - t80 * mrSges(5,1) - m(6) * (t80 * pkin(4) + t119) - m(7) * (t100 * t80 + t120) + t140 * t90 + t143 * (t110 * t91 + t113 * t80) + t142 * (-t110 * t80 + t113 * t91) - t141 * t136 + t139 * t91 - t138 * (t105 * t134 - t90 * t114)) * g(1);
U  = t1;
