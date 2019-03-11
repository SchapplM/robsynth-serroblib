% Calculate potential energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:44
% EndTime: 2019-03-08 20:26:45
% DurationCPUTime: 0.73s
% Computational Cost: add. (400->84), mult. (862->94), div. (0->0), fcn. (1059->14), ass. (0->45)
t125 = cos(qJ(2));
t158 = mrSges(3,2) * t125;
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t114 = qJ(5) + qJ(6);
t111 = sin(t114);
t112 = cos(t114);
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t148 = -m(6) * pkin(4) - m(7) * (pkin(5) * t123 + pkin(4)) - t123 * mrSges(6,1) - t112 * mrSges(7,1) + t120 * mrSges(6,2) + t111 * mrSges(7,2) - mrSges(5,1);
t151 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t157 = t121 * t151 + t124 * t148 - mrSges(4,1);
t117 = sin(pkin(6));
t119 = cos(pkin(6));
t122 = sin(qJ(2));
t141 = t119 * t122;
t154 = -mrSges(3,3) - mrSges(4,3);
t156 = -t141 * mrSges(3,1) - t119 * t158 - mrSges(2,2) + (m(3) * pkin(7) - t148 * t121 + t151 * t124 - t154) * t117;
t155 = -m(5) - m(6);
t115 = sin(pkin(12));
t146 = cos(pkin(12));
t101 = -t122 * t115 + t125 * t146;
t150 = -m(3) * pkin(1) - t125 * mrSges(3,1) + t122 * mrSges(3,2) - mrSges(2,1);
t149 = m(7) * (pkin(5) * t120 + pkin(8)) + t120 * mrSges(6,1) + t111 * mrSges(7,1) + t123 * mrSges(6,2) + t112 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t109 = pkin(2) * t125 + pkin(1);
t116 = sin(pkin(11));
t118 = cos(pkin(11));
t99 = pkin(2) * t141 + (-pkin(7) - qJ(3)) * t117;
t147 = t116 * t109 + t118 * t99;
t138 = t119 * pkin(7) + qJ(1);
t100 = -t125 * t115 - t122 * t146;
t98 = t100 * t119;
t88 = t101 * t116 - t118 * t98;
t137 = t88 * pkin(3) + t147;
t134 = t118 * t109 - t116 * t99;
t133 = t117 * t122 * pkin(2) + t119 * qJ(3) + t138;
t90 = t101 * t118 + t116 * t98;
t132 = t90 * pkin(3) + t134;
t97 = t100 * t117;
t131 = -t97 * pkin(3) + t133;
t127 = t119 * t101;
t96 = t101 * t117;
t89 = t100 * t118 - t116 * t127;
t87 = t116 * t100 + t118 * t127;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t138 - (mrSges(3,1) * t122 + t158) * t117 - m(4) * t133 + t97 * mrSges(4,1) - m(7) * t131 + t155 * (-t96 * pkin(8) + t131) + t148 * (t119 * t121 - t124 * t97) + t149 * t96 + t154 * t119 + t151 * (-t119 * t124 - t121 * t97)) * g(3) + (-m(4) * t147 - m(7) * t137 - mrSges(1,2) + t155 * (-t87 * pkin(8) + t137) + t149 * t87 + t150 * t116 + t157 * t88 + t156 * t118) * g(2) + (-m(4) * t134 - m(7) * t132 - mrSges(1,1) + t155 * (-t89 * pkin(8) + t132) + t149 * t89 + t150 * t118 + t157 * t90 - t156 * t116) * g(1);
U  = t1;
