% Calculate potential energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:38
% EndTime: 2019-03-08 20:40:39
% DurationCPUTime: 0.85s
% Computational Cost: add. (284->92), mult. (491->99), div. (0->0), fcn. (547->12), ass. (0->44)
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t147 = t114 * mrSges(5,1) - t111 * mrSges(5,2);
t144 = -m(4) - m(5);
t143 = m(6) + m(7);
t142 = mrSges(4,1) + mrSges(3,3);
t141 = -m(5) * pkin(3) - t142;
t140 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t139 = t114 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3);
t107 = sin(pkin(6));
t138 = t147 * t107 - mrSges(2,2);
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t137 = -m(7) * pkin(5) - t113 * mrSges(7,1) + t110 * mrSges(7,2) - mrSges(6,1);
t136 = -m(5) * pkin(8) - t110 * mrSges(7,1) - t113 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t106 = sin(pkin(11));
t108 = cos(pkin(11));
t112 = sin(qJ(2));
t109 = cos(pkin(6));
t115 = cos(qJ(2));
t125 = t109 * t115;
t89 = t106 * t112 - t108 * t125;
t134 = t111 * t89;
t91 = t106 * t125 + t108 * t112;
t133 = t111 * t91;
t131 = t106 * t107;
t132 = t108 * pkin(1) + pkin(7) * t131;
t130 = t107 * t108;
t127 = t107 * t115;
t126 = t109 * t112;
t124 = t109 * pkin(7) + qJ(1);
t123 = pkin(7) * t130;
t122 = t107 * t112 * pkin(2) + t124;
t102 = t106 * pkin(1);
t90 = t106 * t115 + t108 * t126;
t120 = t90 * pkin(2) + qJ(3) * t89 + t102;
t92 = -t106 * t126 + t108 * t115;
t119 = t92 * pkin(2) + qJ(3) * t91 + t132;
t116 = -pkin(9) - pkin(8);
t105 = qJ(4) + qJ(5);
t101 = cos(t105);
t100 = sin(t105);
t99 = pkin(4) * t114 + pkin(3);
t1 = (-m(2) * qJ(1) - m(3) * t124 - mrSges(1,3) - mrSges(2,3) + t137 * (-t100 * t127 + t101 * t109) - t143 * (t109 * t99 + t122) - t140 * (t100 * t109 + t101 * t127) + t144 * t122 + (t141 - t147) * t109 + (((t143 - t144) * qJ(3) + t139 + (t143 * pkin(4) + mrSges(5,1)) * t111) * t115 + (t143 * t116 + t136) * t112) * t107) * g(3) + (-mrSges(1,2) - t106 * mrSges(2,1) - m(3) * (t102 - t123) - m(4) * (t120 - t123) - m(5) * t120 - t134 * mrSges(5,1) + t138 * t108 - t139 * t89 - t143 * (pkin(4) * t134 - t90 * t116 + (-pkin(7) - t99) * t130 + t120) + t140 * (t100 * t130 + t101 * t89) + (-m(5) * (-pkin(3) - pkin(7)) + t142) * t130 + t137 * (t100 * t89 - t101 * t130) + t136 * t90) * g(2) + (-m(3) * t132 - t108 * mrSges(2,1) - t133 * mrSges(5,1) - mrSges(1,1) - t138 * t106 - t139 * t91 - t143 * (pkin(4) * t133 - t92 * t116 + t99 * t131 + t119) - t140 * (t100 * t131 - t91 * t101) + t144 * t119 + t141 * t131 + t137 * (t100 * t91 + t101 * t131) + t136 * t92) * g(1);
U  = t1;
