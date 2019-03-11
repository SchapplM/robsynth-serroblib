% Calculate potential energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:20
% EndTime: 2019-03-08 19:34:21
% DurationCPUTime: 0.67s
% Computational Cost: add. (372->91), mult. (827->104), div. (0->0), fcn. (1011->12), ass. (0->54)
t125 = cos(qJ(2));
t158 = t125 * mrSges(3,2);
t117 = sin(pkin(6));
t119 = cos(pkin(6));
t122 = sin(qJ(2));
t139 = t119 * t122;
t156 = -mrSges(3,3) - mrSges(4,3);
t157 = -t139 * mrSges(3,1) + t117 * (m(3) * pkin(7) - t156) - t119 * t158 - mrSges(2,2);
t115 = sin(pkin(11));
t144 = cos(pkin(11));
t104 = -t122 * t115 + t125 * t144;
t154 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t152 = -t120 * mrSges(7,1) - t123 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t151 = -m(3) * pkin(1) - t125 * mrSges(3,1) + t122 * mrSges(3,2) - mrSges(2,1);
t150 = m(7) * (pkin(5) + pkin(8)) + t123 * mrSges(7,1) - t120 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t103 = -t125 * t115 - t122 * t144;
t116 = sin(pkin(10));
t118 = cos(pkin(10));
t127 = t119 * t104;
t89 = t116 * t103 + t118 * t127;
t148 = pkin(8) * t89;
t91 = t103 * t118 - t116 * t127;
t147 = pkin(8) * t91;
t99 = t104 * t117;
t146 = pkin(8) * t99;
t102 = pkin(2) * t139 + (-pkin(7) - qJ(3)) * t117;
t112 = pkin(2) * t125 + pkin(1);
t145 = t118 * t102 + t116 * t112;
t121 = sin(qJ(4));
t141 = t117 * t121;
t124 = cos(qJ(4));
t140 = t117 * t124;
t136 = t119 * pkin(7) + qJ(1);
t101 = t103 * t119;
t90 = -t101 * t118 + t104 * t116;
t135 = t90 * pkin(3) + t145;
t133 = -t102 * t116 + t118 * t112;
t132 = t117 * t122 * pkin(2) + t119 * qJ(3) + t136;
t92 = t101 * t116 + t104 * t118;
t131 = t92 * pkin(3) + t133;
t100 = t103 * t117;
t130 = -t100 * pkin(3) + t132;
t83 = t118 * t140 + t121 * t90;
t84 = -t118 * t141 + t124 * t90;
t129 = t84 * pkin(4) + t83 * qJ(5) + t135;
t85 = -t116 * t140 + t121 * t92;
t86 = t116 * t141 + t124 * t92;
t128 = t86 * pkin(4) + t85 * qJ(5) + t131;
t94 = -t100 * t121 - t119 * t124;
t95 = -t100 * t124 + t119 * t121;
t126 = t95 * pkin(4) + qJ(5) * t94 + t130;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t136 - (t122 * mrSges(3,1) + t158) * t117 - m(4) * t132 + t100 * mrSges(4,1) - m(5) * (t130 - t146) - m(6) * (t126 - t146) - m(7) * t126 + t156 * t119 + t152 * t94 + t150 * t99 + t154 * t95) * g(3) + (-mrSges(1,2) - m(4) * t145 - t90 * mrSges(4,1) - m(5) * (t135 - t148) - m(6) * (t129 - t148) - m(7) * t129 + t151 * t116 + t152 * t83 + t150 * t89 + t154 * t84 + t157 * t118) * g(2) + (-mrSges(1,1) - m(4) * t133 - t92 * mrSges(4,1) - m(5) * (t131 - t147) - m(6) * (t128 - t147) - m(7) * t128 + t151 * t118 + t152 * t85 + t150 * t91 + t154 * t86 - t157 * t116) * g(1);
U  = t1;
