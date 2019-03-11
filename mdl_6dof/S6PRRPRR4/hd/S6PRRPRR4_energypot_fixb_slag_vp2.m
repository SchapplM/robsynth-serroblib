% Calculate potential energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:10
% EndTime: 2019-03-08 22:10:11
% DurationCPUTime: 0.70s
% Computational Cost: add. (334->94), mult. (746->111), div. (0->0), fcn. (900->12), ass. (0->49)
t161 = -m(6) - m(7);
t122 = sin(qJ(6));
t126 = cos(qJ(6));
t162 = -t122 * mrSges(7,1) - t126 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t160 = -mrSges(4,1) - mrSges(5,1);
t159 = mrSges(4,2) - mrSges(5,3);
t158 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t156 = -m(7) * pkin(5) - t126 * mrSges(7,1) + t122 * mrSges(7,2) - mrSges(6,1);
t155 = mrSges(3,2) - t162 + t161 * (pkin(8) - pkin(9));
t153 = cos(qJ(3));
t119 = sin(pkin(11));
t121 = cos(pkin(11));
t125 = sin(qJ(2));
t128 = cos(qJ(2));
t150 = cos(pkin(6));
t139 = t128 * t150;
t104 = t119 * t125 - t121 * t139;
t152 = pkin(8) * t104;
t106 = t119 * t139 + t121 * t125;
t151 = pkin(8) * t106;
t120 = sin(pkin(6));
t149 = t119 * t120;
t148 = t120 * t121;
t124 = sin(qJ(3));
t147 = t120 * t124;
t146 = t120 * t125;
t145 = t120 * t128;
t144 = t121 * pkin(1) + pkin(7) * t149;
t143 = t150 * pkin(7) + qJ(1);
t140 = t125 * t150;
t107 = -t119 * t140 + t121 * t128;
t142 = t107 * pkin(2) + t144;
t141 = t120 * t153;
t137 = t119 * pkin(1) - pkin(7) * t148;
t105 = t119 * t128 + t121 * t140;
t136 = t105 * pkin(2) + t137;
t135 = pkin(2) * t146 - pkin(8) * t145 + t143;
t97 = t107 * t124 - t119 * t141;
t98 = t107 * t153 + t119 * t147;
t134 = t98 * pkin(3) + qJ(4) * t97 + t142;
t95 = t105 * t124 + t121 * t141;
t96 = t105 * t153 - t121 * t147;
t132 = t96 * pkin(3) + qJ(4) * t95 + t136;
t108 = t124 * t146 - t150 * t153;
t109 = t124 * t150 + t125 * t141;
t130 = t109 * pkin(3) + t108 * qJ(4) + t135;
t127 = cos(qJ(5));
t123 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t143 - t150 * mrSges(3,3) - (t125 * mrSges(3,1) + t128 * mrSges(3,2)) * t120 - m(4) * t135 - m(5) * t130 + t161 * (t109 * pkin(4) + pkin(9) * t145 + t130) + t158 * (-t108 * t127 + t109 * t123) + t160 * t109 + t159 * t108 + t156 * (t108 * t123 + t109 * t127) + t162 * t145) * g(3) + (-mrSges(1,2) - t119 * mrSges(2,1) - t121 * mrSges(2,2) - m(3) * t137 - t105 * mrSges(3,1) + mrSges(3,3) * t148 - m(4) * (t136 + t152) - m(5) * (t132 + t152) + t160 * t96 + t159 * t95 + t161 * (t96 * pkin(4) + t132) + t158 * (t123 * t96 - t95 * t127) + t156 * (t123 * t95 + t127 * t96) + t155 * t104) * g(2) + (-mrSges(1,1) - t121 * mrSges(2,1) + t119 * mrSges(2,2) - m(3) * t144 - t107 * mrSges(3,1) - mrSges(3,3) * t149 - m(4) * (t142 + t151) - m(5) * (t134 + t151) + t160 * t98 + t159 * t97 + t161 * (t98 * pkin(4) + t134) + t158 * (t123 * t98 - t97 * t127) + t156 * (t123 * t97 + t127 * t98) + t155 * t106) * g(1);
U  = t1;
