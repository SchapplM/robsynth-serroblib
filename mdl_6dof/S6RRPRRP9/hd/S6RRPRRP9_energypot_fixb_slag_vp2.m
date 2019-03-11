% Calculate potential energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:25
% EndTime: 2019-03-09 12:27:25
% DurationCPUTime: 0.73s
% Computational Cost: add. (343->112), mult. (554->129), div. (0->0), fcn. (633->12), ass. (0->49)
t147 = -mrSges(6,1) - mrSges(7,1);
t146 = -mrSges(6,2) - mrSges(7,2);
t145 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t144 = -mrSges(5,3) - t145;
t143 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t112 = sin(pkin(11));
t142 = pkin(3) * t112;
t115 = cos(pkin(6));
t141 = t115 * pkin(8) + pkin(7);
t118 = sin(qJ(5));
t122 = cos(qJ(2));
t123 = cos(qJ(1));
t131 = t122 * t123;
t119 = sin(qJ(2));
t120 = sin(qJ(1));
t134 = t119 * t120;
t92 = -t115 * t131 + t134;
t140 = t118 * t92;
t132 = t120 * t122;
t133 = t119 * t123;
t94 = t115 * t132 + t133;
t139 = t118 * t94;
t113 = sin(pkin(6));
t138 = t113 * t119;
t137 = t113 * t120;
t136 = t113 * t122;
t135 = t113 * t123;
t130 = t123 * pkin(1) + pkin(8) * t137;
t129 = t112 * t137;
t128 = t118 * t136;
t109 = t120 * pkin(1);
t127 = -pkin(8) * t135 + t109;
t114 = cos(pkin(11));
t104 = pkin(3) * t114 + pkin(2);
t117 = -pkin(9) - qJ(3);
t126 = t104 * t138 + t115 * t142 + t117 * t136 + t141;
t95 = -t115 * t134 + t131;
t125 = pkin(3) * t129 + t95 * t104 - t94 * t117 + t130;
t93 = t115 * t133 + t132;
t124 = t109 + t93 * t104 - t92 * t117 + (-pkin(8) - t142) * t135;
t121 = cos(qJ(5));
t111 = pkin(11) + qJ(4);
t107 = cos(t111);
t106 = sin(t111);
t105 = pkin(5) * t121 + pkin(4);
t89 = t106 * t115 + t107 * t138;
t85 = t106 * t137 + t107 * t95;
t83 = -t106 * t135 + t93 * t107;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(5) * t126 - t89 * mrSges(5,1) + mrSges(5,3) * t136 - m(6) * (pkin(4) * t89 + t126) - m(7) * (-pkin(5) * t128 + t105 * t89 + t126) + t147 * (t121 * t89 - t128) + t146 * (-t118 * t89 - t121 * t136) + (-m(3) - m(4)) * t141 + (-t112 * mrSges(4,1) - t114 * mrSges(4,2) - mrSges(3,3)) * t115 + (t145 * t122 + (-m(4) * pkin(2) - t114 * mrSges(4,1) + t112 * mrSges(4,2) - mrSges(3,1)) * t119) * t113 + t143 * (t106 * t138 - t115 * t107)) * g(3) + (-mrSges(1,2) - t120 * mrSges(2,1) - t123 * mrSges(2,2) - m(3) * t127 - t93 * mrSges(3,1) + mrSges(3,3) * t135 - m(4) * (t93 * pkin(2) + t127) - (-t112 * t135 + t93 * t114) * mrSges(4,1) - (-t93 * t112 - t114 * t135) * mrSges(4,2) - m(5) * t124 - t83 * mrSges(5,1) - m(6) * (t83 * pkin(4) + t124) - m(7) * (pkin(5) * t140 + t83 * t105 + t124) + t147 * (t121 * t83 + t140) + t146 * (-t118 * t83 + t121 * t92) + t144 * t92 + t143 * (t93 * t106 + t107 * t135)) * g(2) + (-mrSges(1,1) - t123 * mrSges(2,1) + t120 * mrSges(2,2) - m(3) * t130 - t95 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * (pkin(2) * t95 + t130) - (t114 * t95 + t129) * mrSges(4,1) - (-t112 * t95 + t114 * t137) * mrSges(4,2) - m(5) * t125 - t85 * mrSges(5,1) - m(6) * (pkin(4) * t85 + t125) - m(7) * (pkin(5) * t139 + t105 * t85 + t125) + t147 * (t121 * t85 + t139) + t146 * (-t118 * t85 + t121 * t94) + t144 * t94 + t143 * (t106 * t95 - t107 * t137)) * g(1);
U  = t1;
