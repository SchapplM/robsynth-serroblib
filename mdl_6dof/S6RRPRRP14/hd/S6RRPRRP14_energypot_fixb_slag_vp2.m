% Calculate potential energy for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:39
% EndTime: 2019-03-09 13:03:40
% DurationCPUTime: 0.58s
% Computational Cost: add. (275->92), mult. (591->105), div. (0->0), fcn. (686->10), ass. (0->47)
t149 = -m(6) - m(7);
t148 = -mrSges(3,1) + mrSges(4,2);
t147 = mrSges(3,3) + mrSges(4,1);
t146 = -mrSges(4,3) + mrSges(3,2);
t145 = -mrSges(5,3) + t148;
t144 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t143 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t142 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t113 = cos(pkin(6));
t141 = t113 * pkin(8) + pkin(7);
t112 = sin(pkin(6));
t116 = sin(qJ(2));
t140 = t112 * t116;
t117 = sin(qJ(1));
t139 = t112 * t117;
t120 = cos(qJ(2));
t138 = t112 * t120;
t137 = t116 * t117;
t121 = cos(qJ(1));
t136 = t116 * t121;
t135 = t117 * t120;
t134 = t120 * t121;
t133 = t121 * t112;
t132 = t121 * pkin(1) + pkin(8) * t139;
t131 = pkin(8) * t133;
t130 = pkin(2) * t140 + t141;
t110 = t117 * pkin(1);
t97 = -t113 * t134 + t137;
t98 = t113 * t136 + t135;
t129 = t98 * pkin(2) + t97 * qJ(3) + t110;
t100 = -t113 * t137 + t134;
t99 = t113 * t135 + t136;
t128 = t100 * pkin(2) + qJ(3) * t99 + t132;
t127 = t113 * pkin(3) + pkin(9) * t140 - qJ(3) * t138 + t130;
t126 = pkin(3) * t139 + pkin(9) * t100 + t128;
t125 = t98 * pkin(9) + (-pkin(3) - pkin(8)) * t133 + t129;
t119 = cos(qJ(4));
t118 = cos(qJ(5));
t115 = sin(qJ(4));
t114 = sin(qJ(5));
t96 = t113 * t119 - t115 * t138;
t95 = t113 * t115 + t119 * t138;
t88 = t97 * t115 - t119 * t133;
t87 = t115 * t133 + t97 * t119;
t86 = t115 * t99 + t119 * t139;
t85 = t115 * t139 - t99 * t119;
t1 = (-m(2) * pkin(7) - m(3) * t141 - m(4) * t130 - m(5) * t127 - t96 * mrSges(5,1) - mrSges(5,3) * t140 - mrSges(1,3) - mrSges(2,3) + t149 * (t96 * pkin(4) + pkin(10) * t95 + t127) + t143 * (t114 * t140 + t118 * t96) + t142 * (t114 * t96 - t118 * t140) - t147 * t113 + ((m(4) * qJ(3) - t146) * t120 + t148 * t116) * t112 + t144 * t95) * g(3) + (-mrSges(1,2) - t117 * mrSges(2,1) - t121 * mrSges(2,2) - m(3) * (t110 - t131) - m(4) * (t129 - t131) - m(5) * t125 - t88 * mrSges(5,1) + t146 * t97 + t149 * (t88 * pkin(4) - t87 * pkin(10) + t125) + t143 * (t114 * t98 + t118 * t88) + t142 * (t114 * t88 - t98 * t118) + t147 * t133 + t145 * t98 - t144 * t87) * g(2) + (-m(3) * t132 - m(4) * t128 - m(5) * t126 - t121 * mrSges(2,1) - t86 * mrSges(5,1) + t117 * mrSges(2,2) - mrSges(1,1) + t146 * t99 + t149 * (t86 * pkin(4) + pkin(10) * t85 + t126) + t143 * (t100 * t114 + t118 * t86) + t142 * (-t100 * t118 + t114 * t86) - t147 * t139 + t144 * t85 + t145 * t100) * g(1);
U  = t1;
