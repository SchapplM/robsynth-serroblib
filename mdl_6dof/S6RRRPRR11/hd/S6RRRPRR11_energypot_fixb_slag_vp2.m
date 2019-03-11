% Calculate potential energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:44
% EndTime: 2019-03-09 19:24:44
% DurationCPUTime: 0.70s
% Computational Cost: add. (334->94), mult. (746->110), div. (0->0), fcn. (900->12), ass. (0->48)
t159 = -m(6) - m(7);
t119 = sin(qJ(6));
t124 = cos(qJ(6));
t160 = -t119 * mrSges(7,1) - t124 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t158 = -mrSges(4,1) - mrSges(5,1);
t157 = mrSges(4,2) - mrSges(5,3);
t156 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t154 = -m(7) * pkin(5) - t124 * mrSges(7,1) + t119 * mrSges(7,2) - mrSges(6,1);
t153 = mrSges(3,2) - t160 + t159 * (pkin(9) - pkin(10));
t151 = cos(qJ(3));
t122 = sin(qJ(2));
t126 = cos(qJ(2));
t127 = cos(qJ(1));
t123 = sin(qJ(1));
t147 = cos(pkin(6));
t139 = t123 * t147;
t107 = t127 * t122 + t126 * t139;
t150 = pkin(9) * t107;
t138 = t127 * t147;
t105 = t122 * t123 - t126 * t138;
t149 = t105 * pkin(9);
t148 = t147 * pkin(8) + pkin(7);
t118 = sin(pkin(6));
t146 = t118 * t122;
t145 = t118 * t123;
t144 = t118 * t126;
t143 = t118 * t127;
t142 = t127 * pkin(1) + pkin(8) * t145;
t108 = -t122 * t139 + t127 * t126;
t141 = t108 * pkin(2) + t142;
t140 = t118 * t151;
t136 = t123 * pkin(1) - pkin(8) * t143;
t106 = t122 * t138 + t123 * t126;
t135 = t106 * pkin(2) + t136;
t134 = pkin(2) * t146 - pkin(9) * t144 + t148;
t121 = sin(qJ(3));
t96 = t108 * t121 - t123 * t140;
t97 = t108 * t151 + t121 * t145;
t133 = t97 * pkin(3) + qJ(4) * t96 + t141;
t94 = t106 * t121 + t127 * t140;
t95 = t106 * t151 - t121 * t143;
t131 = t95 * pkin(3) + t94 * qJ(4) + t135;
t103 = t121 * t146 - t147 * t151;
t104 = t147 * t121 + t122 * t140;
t129 = t104 * pkin(3) + qJ(4) * t103 + t134;
t125 = cos(qJ(5));
t120 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t148 - t147 * mrSges(3,3) - (mrSges(3,1) * t122 + mrSges(3,2) * t126) * t118 - m(4) * t134 - m(5) * t129 + t159 * (t104 * pkin(4) + pkin(10) * t144 + t129) + t156 * (-t103 * t125 + t104 * t120) + t158 * t104 + t157 * t103 + t154 * (t103 * t120 + t104 * t125) + t160 * t144) * g(3) + (-mrSges(1,2) - t123 * mrSges(2,1) - t127 * mrSges(2,2) - m(3) * t136 - t106 * mrSges(3,1) + mrSges(3,3) * t143 - m(4) * (t135 + t149) - m(5) * (t131 + t149) + t158 * t95 + t157 * t94 + t159 * (t95 * pkin(4) + t131) + t156 * (t120 * t95 - t94 * t125) + t154 * (t120 * t94 + t125 * t95) + t153 * t105) * g(2) + (-mrSges(1,1) - t127 * mrSges(2,1) + t123 * mrSges(2,2) - m(3) * t142 - t108 * mrSges(3,1) - mrSges(3,3) * t145 - m(4) * (t141 + t150) - m(5) * (t133 + t150) + t158 * t97 + t157 * t96 + t159 * (t97 * pkin(4) + t133) + t156 * (t120 * t97 - t96 * t125) + t154 * (t120 * t96 + t125 * t97) + t153 * t107) * g(1);
U  = t1;
