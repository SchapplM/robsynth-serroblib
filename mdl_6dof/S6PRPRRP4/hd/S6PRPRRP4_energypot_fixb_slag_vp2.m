% Calculate potential energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:09
% EndTime: 2019-03-08 20:09:09
% DurationCPUTime: 0.70s
% Computational Cost: add. (370->104), mult. (600->124), div. (0->0), fcn. (697->12), ass. (0->47)
t151 = -m(6) - m(7);
t150 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t149 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t148 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t147 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t146 = -mrSges(5,3) - t148;
t118 = sin(pkin(11));
t145 = pkin(3) * t118;
t119 = sin(pkin(10));
t120 = sin(pkin(6));
t144 = t119 * t120;
t122 = cos(pkin(10));
t143 = t120 * t122;
t126 = sin(qJ(2));
t142 = t120 * t126;
t128 = cos(qJ(2));
t141 = t120 * t128;
t123 = cos(pkin(6));
t140 = t123 * t126;
t139 = t123 * t128;
t138 = t122 * pkin(1) + pkin(7) * t144;
t137 = t123 * pkin(7) + qJ(1);
t136 = t118 * t144;
t114 = t119 * pkin(1);
t135 = -pkin(7) * t143 + t114;
t101 = t119 * t139 + t122 * t126;
t102 = -t119 * t140 + t122 * t128;
t121 = cos(pkin(11));
t111 = pkin(3) * t121 + pkin(2);
t124 = -pkin(8) - qJ(3);
t134 = pkin(3) * t136 - t101 * t124 + t102 * t111 + t138;
t133 = t111 * t142 + t123 * t145 + t124 * t141 + t137;
t100 = t119 * t128 + t122 * t140;
t99 = t119 * t126 - t122 * t139;
t130 = t114 + t100 * t111 - t99 * t124 + (-pkin(7) - t145) * t143;
t127 = cos(qJ(5));
t125 = sin(qJ(5));
t117 = pkin(11) + qJ(4);
t113 = cos(t117);
t112 = sin(t117);
t94 = t112 * t123 + t113 * t142;
t93 = t112 * t142 - t123 * t113;
t87 = t102 * t113 + t112 * t144;
t86 = t102 * t112 - t113 * t144;
t85 = t100 * t113 - t112 * t143;
t84 = t100 * t112 + t113 * t143;
t1 = (-m(2) * qJ(1) - m(5) * t133 - t94 * mrSges(5,1) + mrSges(5,3) * t141 - mrSges(1,3) - mrSges(2,3) + t151 * (t94 * pkin(4) + pkin(9) * t93 + t133) + t149 * (-t125 * t141 + t94 * t127) + t147 * (t94 * t125 + t127 * t141) + (-m(3) - m(4)) * t137 + (-t118 * mrSges(4,1) - t121 * mrSges(4,2) - mrSges(3,3)) * t123 + (t148 * t128 + (-m(4) * pkin(2) - t121 * mrSges(4,1) + t118 * mrSges(4,2) - mrSges(3,1)) * t126) * t120 + t150 * t93) * g(3) + (-mrSges(1,2) - t119 * mrSges(2,1) - t122 * mrSges(2,2) - m(3) * t135 - t100 * mrSges(3,1) + mrSges(3,3) * t143 - m(4) * (pkin(2) * t100 + t135) - (t100 * t121 - t118 * t143) * mrSges(4,1) - (-t100 * t118 - t121 * t143) * mrSges(4,2) - m(5) * t130 - t85 * mrSges(5,1) + t151 * (t85 * pkin(4) + t84 * pkin(9) + t130) + t149 * (t125 * t99 + t127 * t85) + t147 * (t125 * t85 - t99 * t127) + t146 * t99 + t150 * t84) * g(2) + (-mrSges(1,1) - t122 * mrSges(2,1) + t119 * mrSges(2,2) - m(3) * t138 - t102 * mrSges(3,1) - mrSges(3,3) * t144 - m(4) * (pkin(2) * t102 + t138) - (t102 * t121 + t136) * mrSges(4,1) - (-t102 * t118 + t121 * t144) * mrSges(4,2) - m(5) * t134 - t87 * mrSges(5,1) + t151 * (t87 * pkin(4) + pkin(9) * t86 + t134) + t149 * (t101 * t125 + t127 * t87) + t147 * (-t101 * t127 + t125 * t87) + t150 * t86 + t146 * t101) * g(1);
U  = t1;
