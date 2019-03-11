% Calculate potential energy for
% S6PRPRPR1
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:33
% EndTime: 2019-03-08 19:24:34
% DurationCPUTime: 0.81s
% Computational Cost: add. (383->91), mult. (733->110), div. (0->0), fcn. (880->14), ass. (0->53)
t138 = cos(qJ(4));
t172 = t138 * mrSges(5,2);
t139 = cos(qJ(2));
t171 = t139 * mrSges(3,2);
t170 = -m(4) - m(5);
t169 = -m(6) - m(7);
t168 = mrSges(3,3) + mrSges(4,3);
t127 = sin(pkin(11));
t130 = cos(pkin(11));
t136 = sin(qJ(2));
t167 = t136 * t127 - t130 * t139;
t166 = m(3) * pkin(7) + t168;
t165 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t164 = -m(3) * pkin(1) - mrSges(3,1) * t139 + mrSges(3,2) * t136 - mrSges(2,1);
t135 = sin(qJ(4));
t163 = -m(5) * pkin(3) - t138 * mrSges(5,1) + t135 * mrSges(5,2) - mrSges(4,1);
t134 = sin(qJ(6));
t137 = cos(qJ(6));
t162 = -m(7) * pkin(5) - mrSges(7,1) * t137 + mrSges(7,2) * t134 - mrSges(6,1);
t129 = sin(pkin(6));
t132 = cos(pkin(6));
t153 = t132 * t136;
t161 = mrSges(3,1) * t153 - t129 * t172 + t132 * t171 + mrSges(2,2);
t160 = m(5) * pkin(8) + mrSges(7,1) * t134 + mrSges(7,2) * t137 - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t128 = sin(pkin(10));
t159 = t128 * t129;
t131 = cos(pkin(10));
t158 = t129 * t131;
t157 = t129 * t135;
t154 = t132 * t135;
t109 = pkin(2) * t153 + (-pkin(7) - qJ(3)) * t129;
t121 = pkin(2) * t139 + pkin(1);
t150 = t131 * t109 + t128 * t121;
t149 = t132 * pkin(7) + qJ(1);
t148 = t128 * t157;
t147 = t131 * t157;
t146 = -t109 * t128 + t131 * t121;
t145 = t129 * t136 * pkin(2) + t132 * qJ(3) + t149;
t144 = t127 * t139 + t136 * t130;
t142 = t167 * t132;
t133 = -qJ(5) - pkin(8);
t126 = qJ(4) + pkin(12);
t123 = cos(t126);
t122 = sin(t126);
t120 = pkin(4) * t138 + pkin(3);
t108 = t144 * t132;
t107 = t144 * t129;
t106 = t167 * t129;
t99 = -t108 * t128 - t131 * t167;
t98 = t128 * t142 - t131 * t144;
t97 = t108 * t131 - t128 * t167;
t96 = -t128 * t144 - t131 * t142;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t149 - (t136 * mrSges(3,1) + t171) * t129 - t154 * mrSges(5,1) + t170 * t145 + t169 * (pkin(4) * t154 - t106 * t133 + t107 * t120 + t145) + t163 * t107 + (-t168 - t172) * t132 + t165 * (t107 * t122 - t132 * t123) + t162 * (t107 * t123 + t122 * t132) - t160 * t106) * g(3) + (t147 * mrSges(5,1) - mrSges(1,2) + t170 * t150 + t163 * t97 + t169 * (-pkin(4) * t147 + t97 * t120 + t96 * t133 + t150) + t165 * (t122 * t97 + t123 * t158) - t161 * t131 + t164 * t128 + t166 * t158 + t162 * (-t122 * t158 + t123 * t97) + t160 * t96) * g(2) + (-t148 * mrSges(5,1) - mrSges(1,1) + t170 * t146 + t163 * t99 + t169 * (pkin(4) * t148 + t99 * t120 + t98 * t133 + t146) + t165 * (t122 * t99 - t123 * t159) + t161 * t128 + t164 * t131 - t166 * t159 + t162 * (t122 * t159 + t123 * t99) + t160 * t98) * g(1);
U  = t1;
