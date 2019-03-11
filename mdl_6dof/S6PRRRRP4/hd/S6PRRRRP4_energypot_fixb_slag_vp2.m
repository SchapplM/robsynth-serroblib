% Calculate potential energy for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:01
% EndTime: 2019-03-09 00:13:02
% DurationCPUTime: 0.66s
% Computational Cost: add. (346->100), mult. (669->121), div. (0->0), fcn. (793->12), ass. (0->50)
t158 = -m(4) - m(5);
t157 = -m(6) - m(7);
t156 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t155 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t128 = cos(qJ(4));
t154 = -t128 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t153 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t125 = sin(qJ(4));
t152 = -m(5) * pkin(3) - t128 * mrSges(5,1) + t125 * mrSges(5,2) - mrSges(4,1);
t121 = sin(pkin(11));
t123 = cos(pkin(11));
t127 = sin(qJ(2));
t124 = cos(pkin(6));
t130 = cos(qJ(2));
t142 = t124 * t130;
t103 = t121 * t127 - t123 * t142;
t151 = t103 * t125;
t105 = t121 * t142 + t123 * t127;
t150 = t105 * t125;
t122 = sin(pkin(6));
t149 = t121 * t122;
t148 = t122 * t123;
t126 = sin(qJ(3));
t147 = t122 * t126;
t146 = t122 * t127;
t129 = cos(qJ(3));
t145 = t122 * t129;
t144 = t122 * t130;
t143 = t124 * t127;
t141 = t123 * pkin(1) + pkin(7) * t149;
t140 = t124 * pkin(7) + qJ(1);
t139 = pkin(2) * t146 + t140;
t138 = t121 * pkin(1) - pkin(7) * t148;
t106 = -t121 * t143 + t123 * t130;
t137 = t106 * pkin(2) + pkin(8) * t105 + t141;
t136 = -pkin(8) * t144 + t139;
t104 = t121 * t130 + t123 * t143;
t135 = t104 * pkin(2) + pkin(8) * t103 + t138;
t131 = -pkin(10) - pkin(9);
t120 = qJ(4) + qJ(5);
t116 = cos(t120);
t115 = sin(t120);
t114 = pkin(4) * t128 + pkin(3);
t108 = t124 * t126 + t127 * t145;
t107 = -t124 * t129 + t126 * t146;
t94 = t106 * t129 + t121 * t147;
t93 = t106 * t126 - t121 * t145;
t92 = t104 * t129 - t123 * t147;
t91 = t104 * t126 + t123 * t145;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t140 - t124 * mrSges(3,3) - (t127 * mrSges(3,1) + t130 * mrSges(3,2)) * t122 - m(4) * t136 - t108 * mrSges(4,1) + mrSges(4,3) * t144 - m(5) * (pkin(3) * t108 + t136) - (t108 * t128 - t125 * t144) * mrSges(5,1) - (-t108 * t125 - t128 * t144) * mrSges(5,2) + t157 * (-t107 * t131 + t108 * t114 + (-pkin(4) * t125 - pkin(8)) * t144 + t139) + t156 * (t108 * t116 - t115 * t144) + t155 * (t108 * t115 + t116 * t144) + t153 * t107) * g(3) + (-m(3) * t138 - mrSges(2,1) * t121 - t104 * mrSges(3,1) - t151 * mrSges(5,1) - mrSges(2,2) * t123 + mrSges(3,3) * t148 - mrSges(1,2) + t158 * t135 + t152 * t92 + t157 * (pkin(4) * t151 + t92 * t114 - t91 * t131 + t135) + t156 * (t103 * t115 + t116 * t92) + t155 * (-t103 * t116 + t115 * t92) + t154 * t103 + t153 * t91) * g(2) + (-m(3) * t141 - mrSges(2,1) * t123 - t106 * mrSges(3,1) - t150 * mrSges(5,1) + mrSges(2,2) * t121 - mrSges(3,3) * t149 - mrSges(1,1) + t158 * t137 + t152 * t94 + t157 * (pkin(4) * t150 + t94 * t114 - t93 * t131 + t137) + t156 * (t105 * t115 + t116 * t94) + t155 * (-t105 * t116 + t115 * t94) + t154 * t105 + t153 * t93) * g(1);
U  = t1;
