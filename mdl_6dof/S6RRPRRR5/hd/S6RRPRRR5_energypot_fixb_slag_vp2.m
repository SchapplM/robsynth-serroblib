% Calculate potential energy for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:37:56
% EndTime: 2019-03-09 13:37:57
% DurationCPUTime: 0.74s
% Computational Cost: add. (400->93), mult. (862->112), div. (0->0), fcn. (1059->14), ass. (0->49)
t152 = -m(5) - m(6);
t151 = -mrSges(3,3) - mrSges(4,3);
t112 = sin(pkin(12));
t117 = sin(qJ(2));
t121 = cos(qJ(2));
t142 = cos(pkin(12));
t98 = -t117 * t112 + t121 * t142;
t150 = -m(3) * pkin(1) - mrSges(2,1);
t149 = m(3) * pkin(8) - t151;
t148 = -m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t111 = qJ(5) + qJ(6);
t108 = sin(t111);
t109 = cos(t111);
t115 = sin(qJ(5));
t119 = cos(qJ(5));
t147 = m(7) * (pkin(5) * t115 + pkin(9)) + t115 * mrSges(6,1) + t108 * mrSges(7,1) + t119 * mrSges(6,2) + t109 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t146 = -m(6) * pkin(4) - m(7) * (pkin(5) * t119 + pkin(4)) - t119 * mrSges(6,1) - t109 * mrSges(7,1) + t115 * mrSges(6,2) + t108 * mrSges(7,2) - mrSges(5,1);
t145 = pkin(2) * t117;
t114 = cos(pkin(6));
t144 = t114 * pkin(8) + pkin(7);
t106 = pkin(2) * t121 + pkin(1);
t118 = sin(qJ(1));
t122 = cos(qJ(1));
t113 = sin(pkin(6));
t96 = t114 * t145 + (-pkin(8) - qJ(3)) * t113;
t143 = t118 * t106 + t122 * t96;
t141 = t113 * t118;
t139 = t117 * t118;
t138 = t117 * t122;
t137 = t118 * t121;
t136 = t121 * t122;
t135 = t122 * t113;
t97 = -t121 * t112 - t117 * t142;
t95 = t97 * t114;
t85 = t118 * t98 - t122 * t95;
t134 = t85 * pkin(3) + t143;
t131 = t122 * t106 - t118 * t96;
t130 = t114 * qJ(3) + t113 * t145 + t144;
t87 = t118 * t95 + t122 * t98;
t129 = t87 * pkin(3) + t131;
t94 = t97 * t113;
t128 = -t94 * pkin(3) + t130;
t124 = t114 * t98;
t120 = cos(qJ(4));
t116 = sin(qJ(4));
t93 = t98 * t113;
t86 = -t118 * t124 + t122 * t97;
t84 = t118 * t97 + t122 * t124;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t144 - (t117 * mrSges(3,1) + t121 * mrSges(3,2)) * t113 - m(4) * t130 + t94 * mrSges(4,1) - m(7) * t128 + t152 * (-t93 * pkin(9) + t128) + t146 * (t114 * t116 - t120 * t94) + t147 * t93 + t151 * t114 + t148 * (-t114 * t120 - t116 * t94)) * g(3) + (-mrSges(1,2) - t122 * mrSges(2,2) - (t114 * t138 + t137) * mrSges(3,1) - (t114 * t136 - t139) * mrSges(3,2) - m(4) * t143 - t85 * mrSges(4,1) - m(7) * t134 + t152 * (-t84 * pkin(9) + t134) + t146 * (-t116 * t135 + t120 * t85) + t147 * t84 + t150 * t118 + t149 * t135 + t148 * (t116 * t85 + t120 * t135)) * g(2) + (-mrSges(1,1) + t118 * mrSges(2,2) - (-t114 * t139 + t136) * mrSges(3,1) - (-t114 * t137 - t138) * mrSges(3,2) - m(4) * t131 - t87 * mrSges(4,1) - m(7) * t129 + t152 * (-t86 * pkin(9) + t129) + t146 * (t116 * t141 + t120 * t87) + t147 * t86 + t150 * t122 - t149 * t141 + t148 * (t116 * t87 - t120 * t141)) * g(1);
U  = t1;
