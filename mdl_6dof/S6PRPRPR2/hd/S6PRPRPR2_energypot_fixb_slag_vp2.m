% Calculate potential energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:26
% EndTime: 2019-03-08 19:29:27
% DurationCPUTime: 0.77s
% Computational Cost: add. (400->84), mult. (862->94), div. (0->0), fcn. (1059->14), ass. (0->45)
t126 = cos(qJ(2));
t158 = mrSges(3,2) * t126;
t123 = sin(qJ(4));
t125 = cos(qJ(4));
t114 = pkin(12) + qJ(6);
t110 = sin(t114);
t111 = cos(t114);
t115 = sin(pkin(12));
t119 = cos(pkin(12));
t148 = -m(6) * pkin(4) - m(7) * (pkin(5) * t119 + pkin(4)) - t119 * mrSges(6,1) - t111 * mrSges(7,1) + t115 * mrSges(6,2) + t110 * mrSges(7,2) - mrSges(5,1);
t151 = -m(6) * qJ(5) + m(7) * (-pkin(9) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t157 = t151 * t123 + t148 * t125 - mrSges(4,1);
t118 = sin(pkin(6));
t121 = cos(pkin(6));
t124 = sin(qJ(2));
t141 = t121 * t124;
t154 = -mrSges(3,3) - mrSges(4,3);
t156 = -t141 * mrSges(3,1) - t121 * t158 - mrSges(2,2) + (m(3) * pkin(7) - t148 * t123 + t151 * t125 - t154) * t118;
t155 = -m(5) - m(6);
t116 = sin(pkin(11));
t146 = cos(pkin(11));
t101 = -t124 * t116 + t126 * t146;
t150 = -m(3) * pkin(1) - t126 * mrSges(3,1) + t124 * mrSges(3,2) - mrSges(2,1);
t149 = m(7) * (pkin(5) * t115 + pkin(8)) + t115 * mrSges(6,1) + t110 * mrSges(7,1) + t119 * mrSges(6,2) + t111 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t109 = pkin(2) * t126 + pkin(1);
t117 = sin(pkin(10));
t120 = cos(pkin(10));
t99 = pkin(2) * t141 + (-pkin(7) - qJ(3)) * t118;
t147 = t117 * t109 + t120 * t99;
t138 = t121 * pkin(7) + qJ(1);
t100 = -t126 * t116 - t124 * t146;
t98 = t100 * t121;
t88 = t101 * t117 - t120 * t98;
t137 = t88 * pkin(3) + t147;
t134 = t120 * t109 - t117 * t99;
t133 = t118 * t124 * pkin(2) + t121 * qJ(3) + t138;
t90 = t101 * t120 + t117 * t98;
t132 = t90 * pkin(3) + t134;
t97 = t100 * t118;
t131 = -t97 * pkin(3) + t133;
t127 = t121 * t101;
t96 = t101 * t118;
t89 = t100 * t120 - t117 * t127;
t87 = t117 * t100 + t120 * t127;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t138 - (mrSges(3,1) * t124 + t158) * t118 - m(4) * t133 + t97 * mrSges(4,1) - m(7) * t131 + t155 * (-pkin(8) * t96 + t131) + t148 * (t121 * t123 - t125 * t97) + t149 * t96 + t154 * t121 + t151 * (-t121 * t125 - t123 * t97)) * g(3) + (-m(4) * t147 - m(7) * t137 - mrSges(1,2) + t155 * (-pkin(8) * t87 + t137) + t149 * t87 + t150 * t117 + t157 * t88 + t156 * t120) * g(2) + (-m(4) * t134 - m(7) * t132 - mrSges(1,1) + t155 * (-pkin(8) * t89 + t132) + t149 * t89 + t150 * t120 + t157 * t90 - t156 * t117) * g(1);
U  = t1;
