% Calculate potential energy for
% S6PRPRRP6
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:02
% EndTime: 2019-03-08 20:18:03
% DurationCPUTime: 0.59s
% Computational Cost: add. (275->92), mult. (591->107), div. (0->0), fcn. (686->10), ass. (0->45)
t147 = -m(6) - m(7);
t146 = -mrSges(3,1) + mrSges(4,2);
t145 = mrSges(3,3) + mrSges(4,1);
t144 = -mrSges(4,3) + mrSges(3,2);
t143 = -mrSges(5,3) + t146;
t142 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t141 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t140 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t112 = sin(pkin(10));
t113 = sin(pkin(6));
t139 = t112 * t113;
t118 = sin(qJ(2));
t138 = t113 * t118;
t121 = cos(qJ(2));
t137 = t113 * t121;
t114 = cos(pkin(10));
t136 = t114 * t113;
t115 = cos(pkin(6));
t135 = t115 * t118;
t134 = t115 * t121;
t133 = t114 * pkin(1) + pkin(7) * t139;
t132 = t115 * pkin(7) + qJ(1);
t131 = pkin(7) * t136;
t130 = pkin(2) * t138 + t132;
t108 = t112 * pkin(1);
t95 = t112 * t118 - t114 * t134;
t96 = t112 * t121 + t114 * t135;
t129 = t96 * pkin(2) + t95 * qJ(3) + t108;
t97 = t112 * t134 + t114 * t118;
t98 = -t112 * t135 + t114 * t121;
t128 = t98 * pkin(2) + t97 * qJ(3) + t133;
t127 = pkin(3) * t139 + t98 * pkin(8) + t128;
t126 = t115 * pkin(3) + pkin(8) * t138 - qJ(3) * t137 + t130;
t125 = t96 * pkin(8) + (-pkin(3) - pkin(7)) * t136 + t129;
t120 = cos(qJ(4));
t119 = cos(qJ(5));
t117 = sin(qJ(4));
t116 = sin(qJ(5));
t100 = t115 * t120 - t117 * t137;
t99 = t115 * t117 + t120 * t137;
t86 = t95 * t117 - t120 * t136;
t85 = t117 * t136 + t95 * t120;
t84 = t97 * t117 + t120 * t139;
t83 = t117 * t139 - t97 * t120;
t1 = (-m(2) * qJ(1) - m(3) * t132 - m(4) * t130 - m(5) * t126 - t100 * mrSges(5,1) - mrSges(5,3) * t138 - mrSges(1,3) - mrSges(2,3) + t147 * (t100 * pkin(4) + t99 * pkin(9) + t126) + t141 * (t100 * t119 + t116 * t138) + t140 * (t100 * t116 - t119 * t138) - t145 * t115 + ((m(4) * qJ(3) - t144) * t121 + t146 * t118) * t113 + t142 * t99) * g(3) + (-mrSges(1,2) - t112 * mrSges(2,1) - t114 * mrSges(2,2) - m(3) * (t108 - t131) - m(4) * (t129 - t131) - m(5) * t125 - t86 * mrSges(5,1) + t144 * t95 + t147 * (t86 * pkin(4) - t85 * pkin(9) + t125) + t141 * (t96 * t116 + t86 * t119) + t140 * (t86 * t116 - t96 * t119) + t145 * t136 + t143 * t96 - t142 * t85) * g(2) + (-m(3) * t133 - m(4) * t128 - m(5) * t127 - t114 * mrSges(2,1) - t84 * mrSges(5,1) + t112 * mrSges(2,2) - mrSges(1,1) + t144 * t97 + t147 * (t84 * pkin(4) + t83 * pkin(9) + t127) + t141 * (t98 * t116 + t84 * t119) + t140 * (t84 * t116 - t98 * t119) - t145 * t139 + t143 * t98 + t142 * t83) * g(1);
U  = t1;
