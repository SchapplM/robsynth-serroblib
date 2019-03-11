% Calculate potential energy for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:08
% EndTime: 2019-03-08 19:55:08
% DurationCPUTime: 0.74s
% Computational Cost: add. (388->97), mult. (862->116), div. (0->0), fcn. (1059->12), ass. (0->50)
t132 = cos(qJ(2));
t164 = t132 * mrSges(3,2);
t128 = sin(qJ(4));
t155 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t163 = t155 * t128 - mrSges(4,1);
t123 = sin(pkin(6));
t125 = cos(pkin(6));
t131 = cos(qJ(4));
t129 = sin(qJ(2));
t148 = t125 * t129;
t159 = -mrSges(3,3) - mrSges(4,3);
t162 = -t148 * mrSges(3,1) - t125 * t164 - mrSges(2,2) + (m(3) * pkin(7) + t155 * t131 - t159) * t123;
t161 = -mrSges(6,1) - mrSges(7,1);
t160 = -mrSges(6,2) - mrSges(7,2);
t121 = sin(pkin(11));
t153 = cos(pkin(11));
t110 = -t129 * t121 + t132 * t153;
t127 = sin(qJ(5));
t157 = m(7) * (pkin(5) * t127 + pkin(8)) - mrSges(4,2) + mrSges(5,3);
t154 = -m(3) * pkin(1) - mrSges(3,1) * t132 + mrSges(3,2) * t129 - mrSges(2,1);
t150 = t123 * t128;
t108 = pkin(2) * t148 + (-pkin(7) - qJ(3)) * t123;
t118 = pkin(2) * t132 + pkin(1);
t122 = sin(pkin(10));
t124 = cos(pkin(10));
t145 = t124 * t108 + t122 * t118;
t144 = t125 * pkin(7) + qJ(1);
t109 = -t132 * t121 - t129 * t153;
t107 = t109 * t125;
t97 = -t107 * t124 + t110 * t122;
t143 = t97 * pkin(3) + t145;
t140 = -t108 * t122 + t124 * t118;
t139 = t123 * t129 * pkin(2) + t125 * qJ(3) + t144;
t99 = t107 * t122 + t110 * t124;
t138 = t99 * pkin(3) + t140;
t106 = t109 * t123;
t137 = -t106 * pkin(3) + t139;
t133 = t125 * t110;
t96 = t122 * t109 + t124 * t133;
t136 = -pkin(8) * t96 + t143;
t98 = t109 * t124 - t122 * t133;
t135 = -pkin(8) * t98 + t138;
t105 = t110 * t123;
t134 = -t105 * pkin(8) + t137;
t130 = cos(qJ(5));
t117 = pkin(5) * t130 + pkin(4);
t101 = -t106 * t131 + t125 * t128;
t93 = t122 * t150 + t131 * t99;
t91 = -t124 * t150 + t131 * t97;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t144 - (t129 * mrSges(3,1) + t164) * t123 - m(4) * t139 + t106 * mrSges(4,1) - m(5) * t134 - t101 * mrSges(5,1) - m(6) * (pkin(4) * t101 + t134) - m(7) * (t101 * t117 + t137) + t161 * (t101 * t130 - t105 * t127) + t160 * (-t101 * t127 - t105 * t130) + t159 * t125 + t157 * t105 + t155 * (-t106 * t128 - t125 * t131)) * g(3) + (-mrSges(1,2) - m(4) * t145 - m(5) * t136 - t91 * mrSges(5,1) - m(6) * (pkin(4) * t91 + t136) - m(7) * (t117 * t91 + t143) + t157 * t96 + t161 * (-t127 * t96 + t130 * t91) + t160 * (-t127 * t91 - t130 * t96) + t154 * t122 + t163 * t97 + t162 * t124) * g(2) + (-mrSges(1,1) - m(4) * t140 - m(5) * t135 - t93 * mrSges(5,1) - m(6) * (pkin(4) * t93 + t135) - m(7) * (t117 * t93 + t138) + t157 * t98 + t161 * (-t127 * t98 + t130 * t93) + t160 * (-t127 * t93 - t130 * t98) + t154 * t124 + t163 * t99 - t162 * t122) * g(1);
U  = t1;
