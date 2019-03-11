% Calculate potential energy for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:18
% EndTime: 2019-03-09 17:41:19
% DurationCPUTime: 0.58s
% Computational Cost: add. (289->96), mult. (621->109), div. (0->0), fcn. (727->10), ass. (0->50)
t142 = pkin(4) + pkin(9);
t148 = -mrSges(6,1) - mrSges(7,1);
t147 = -mrSges(6,2) - mrSges(7,2);
t109 = sin(qJ(5));
t146 = -m(7) * (pkin(5) * t109 + qJ(4)) + mrSges(4,2) - mrSges(5,3);
t113 = cos(qJ(5));
t145 = m(6) * t142 + m(7) * (t113 * pkin(5) + t142) + mrSges(5,1) + mrSges(4,3);
t144 = mrSges(3,2) - t145;
t143 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t111 = sin(qJ(2));
t112 = sin(qJ(1));
t114 = cos(qJ(2));
t115 = cos(qJ(1));
t136 = cos(pkin(6));
t124 = t115 * t136;
t93 = t112 * t111 - t114 * t124;
t141 = t93 * pkin(9);
t125 = t112 * t136;
t95 = t115 * t111 + t114 * t125;
t140 = t95 * pkin(9);
t139 = cos(qJ(3));
t137 = t136 * pkin(8) + pkin(7);
t107 = sin(pkin(6));
t135 = t107 * t111;
t134 = t107 * t114;
t133 = t112 * t107;
t132 = t115 * t107;
t131 = t115 * pkin(1) + pkin(8) * t133;
t130 = pkin(9) * t134;
t129 = pkin(2) * t135 + t137;
t96 = -t111 * t125 + t115 * t114;
t128 = t96 * pkin(2) + t131;
t127 = t107 * t139;
t110 = sin(qJ(3));
t92 = t136 * t110 + t111 * t127;
t123 = t92 * pkin(3) + t129;
t87 = t110 * t133 + t96 * t139;
t122 = t87 * pkin(3) + t128;
t121 = t112 * pkin(1) - pkin(8) * t132;
t94 = t111 * t124 + t112 * t114;
t120 = t94 * pkin(2) + t121;
t85 = -t110 * t132 + t94 * t139;
t119 = t85 * pkin(3) + t120;
t91 = t110 * t135 - t136 * t139;
t118 = t91 * qJ(4) + t123;
t86 = t96 * t110 - t112 * t127;
t117 = t86 * qJ(4) + t122;
t84 = t94 * t110 + t115 * t127;
t116 = t84 * qJ(4) + t119;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t137 - t136 * mrSges(3,3) - (t111 * mrSges(3,1) + t114 * mrSges(3,2)) * t107 - m(4) * (t129 - t130) - m(5) * (t118 - t130) - m(6) * t118 - m(7) * t123 + t146 * t91 + t148 * (t91 * t109 - t113 * t134) + t147 * (t109 * t134 + t91 * t113) + t145 * t134 + t143 * t92) * g(3) + (-mrSges(1,2) - t112 * mrSges(2,1) - t115 * mrSges(2,2) - m(3) * t121 - t94 * mrSges(3,1) + mrSges(3,3) * t132 - m(4) * (t120 + t141) - m(5) * (t116 + t141) - m(6) * t116 - m(7) * t119 + t146 * t84 + t148 * (t84 * t109 + t93 * t113) + t147 * (-t93 * t109 + t84 * t113) + t144 * t93 + t143 * t85) * g(2) + (-mrSges(1,1) - t115 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t131 - t96 * mrSges(3,1) - mrSges(3,3) * t133 - m(4) * (t128 + t140) - m(5) * (t117 + t140) - m(6) * t117 - m(7) * t122 + t146 * t86 + t148 * (t86 * t109 + t95 * t113) + t147 * (-t95 * t109 + t86 * t113) + t144 * t95 + t143 * t87) * g(1);
U  = t1;
