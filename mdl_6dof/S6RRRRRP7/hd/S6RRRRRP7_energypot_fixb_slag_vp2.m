% Calculate potential energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:23
% EndTime: 2019-03-10 01:34:24
% DurationCPUTime: 0.77s
% Computational Cost: add. (343->112), mult. (554->129), div. (0->0), fcn. (633->12), ass. (0->49)
t147 = -mrSges(6,1) - mrSges(7,1);
t146 = -mrSges(6,2) - mrSges(7,2);
t145 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3);
t144 = -mrSges(5,3) - t145;
t143 = -m(6) * pkin(11) + m(7) * (-qJ(6) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t116 = sin(qJ(3));
t142 = pkin(3) * t116;
t113 = cos(pkin(6));
t141 = t113 * pkin(8) + pkin(7);
t115 = sin(qJ(5));
t121 = cos(qJ(2));
t122 = cos(qJ(1));
t131 = t121 * t122;
t117 = sin(qJ(2));
t118 = sin(qJ(1));
t134 = t117 * t118;
t92 = -t113 * t131 + t134;
t140 = t115 * t92;
t132 = t118 * t121;
t133 = t117 * t122;
t94 = t113 * t132 + t133;
t139 = t115 * t94;
t112 = sin(pkin(6));
t138 = t112 * t117;
t137 = t112 * t118;
t136 = t112 * t121;
t135 = t112 * t122;
t130 = t122 * pkin(1) + pkin(8) * t137;
t129 = t115 * t136;
t128 = t116 * t137;
t109 = t118 * pkin(1);
t127 = -pkin(8) * t135 + t109;
t120 = cos(qJ(3));
t105 = pkin(3) * t120 + pkin(2);
t123 = -pkin(10) - pkin(9);
t95 = -t113 * t134 + t131;
t126 = pkin(3) * t128 + t95 * t105 - t94 * t123 + t130;
t125 = t105 * t138 + t113 * t142 + t123 * t136 + t141;
t93 = t113 * t133 + t132;
t124 = t109 + t93 * t105 - t92 * t123 + (-pkin(8) - t142) * t135;
t119 = cos(qJ(5));
t111 = qJ(3) + qJ(4);
t107 = cos(t111);
t106 = sin(t111);
t104 = pkin(5) * t119 + pkin(4);
t89 = t106 * t113 + t107 * t138;
t85 = t106 * t137 + t107 * t95;
t83 = -t106 * t135 + t107 * t93;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(5) * t125 - t89 * mrSges(5,1) + mrSges(5,3) * t136 - m(6) * (t89 * pkin(4) + t125) - m(7) * (-pkin(5) * t129 + t104 * t89 + t125) + t147 * (t119 * t89 - t129) + t146 * (-t115 * t89 - t119 * t136) + (-m(3) - m(4)) * t141 + (-t116 * mrSges(4,1) - t120 * mrSges(4,2) - mrSges(3,3)) * t113 + (t145 * t121 + (-m(4) * pkin(2) - t120 * mrSges(4,1) + t116 * mrSges(4,2) - mrSges(3,1)) * t117) * t112 + t143 * (t106 * t138 - t113 * t107)) * g(3) + (-mrSges(1,2) - t118 * mrSges(2,1) - t122 * mrSges(2,2) - m(3) * t127 - t93 * mrSges(3,1) + mrSges(3,3) * t135 - m(4) * (pkin(2) * t93 + t127) - (-t116 * t135 + t120 * t93) * mrSges(4,1) - (-t116 * t93 - t120 * t135) * mrSges(4,2) - m(5) * t124 - t83 * mrSges(5,1) - m(6) * (pkin(4) * t83 + t124) - m(7) * (pkin(5) * t140 + t104 * t83 + t124) + t147 * (t119 * t83 + t140) + t146 * (-t115 * t83 + t119 * t92) + t144 * t92 + t143 * (t106 * t93 + t107 * t135)) * g(2) + (-mrSges(1,1) - t122 * mrSges(2,1) + t118 * mrSges(2,2) - m(3) * t130 - t95 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * (pkin(2) * t95 + t130) - (t120 * t95 + t128) * mrSges(4,1) - (-t116 * t95 + t120 * t137) * mrSges(4,2) - m(5) * t126 - t85 * mrSges(5,1) - m(6) * (pkin(4) * t85 + t126) - m(7) * (pkin(5) * t139 + t104 * t85 + t126) + t147 * (t119 * t85 + t139) + t146 * (-t115 * t85 + t119 * t94) + t144 * t94 + t143 * (t106 * t95 - t107 * t137)) * g(1);
U  = t1;
