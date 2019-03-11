% Calculate potential energy for
% S6PRPRRP2
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:50
% EndTime: 2019-03-08 19:59:51
% DurationCPUTime: 0.67s
% Computational Cost: add. (418->92), mult. (947->115), div. (0->0), fcn. (1177->12), ass. (0->50)
t142 = cos(qJ(2));
t174 = mrSges(3,2) * t142;
t133 = sin(pkin(6));
t136 = cos(pkin(6));
t139 = sin(qJ(2));
t157 = t136 * t139;
t170 = -mrSges(3,3) - mrSges(4,3);
t173 = -t157 * mrSges(3,1) + (m(3) * pkin(7) - t170) * t133 - t136 * t174 - mrSges(2,2);
t172 = -m(6) - m(7);
t171 = mrSges(4,2) - mrSges(5,3);
t131 = sin(pkin(11));
t134 = cos(pkin(11));
t169 = t139 * t131 - t134 * t142;
t168 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t166 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t165 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t163 = -m(3) * pkin(1) - t142 * mrSges(3,1) + t139 * mrSges(3,2) - mrSges(2,1);
t138 = sin(qJ(4));
t160 = t133 * t138;
t141 = cos(qJ(4));
t159 = t133 * t141;
t119 = pkin(2) * t157 + (-pkin(7) - qJ(3)) * t133;
t128 = pkin(2) * t142 + pkin(1);
t132 = sin(pkin(10));
t135 = cos(pkin(10));
t154 = t135 * t119 + t132 * t128;
t153 = t136 * pkin(7) + qJ(1);
t152 = -t119 * t132 + t135 * t128;
t151 = t133 * t139 * pkin(2) + t136 * qJ(3) + t153;
t150 = t131 * t142 + t139 * t134;
t148 = t169 * t136;
t105 = -t132 * t150 - t135 * t148;
t118 = t150 * t136;
t106 = t118 * t135 - t132 * t169;
t149 = t106 * pkin(3) - pkin(8) * t105 + t154;
t107 = t132 * t148 - t135 * t150;
t108 = -t118 * t132 - t135 * t169;
t147 = t108 * pkin(3) - pkin(8) * t107 + t152;
t116 = t169 * t133;
t117 = t150 * t133;
t146 = t117 * pkin(3) + pkin(8) * t116 + t151;
t140 = cos(qJ(5));
t137 = sin(qJ(5));
t111 = t117 * t141 + t136 * t138;
t110 = t117 * t138 - t136 * t141;
t100 = t108 * t141 + t132 * t160;
t99 = t108 * t138 - t132 * t159;
t98 = t106 * t141 - t135 * t160;
t97 = t106 * t138 + t135 * t159;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t153 - (mrSges(3,1) * t139 + t174) * t133 - m(4) * t151 - t117 * mrSges(4,1) - m(5) * t146 - t111 * mrSges(5,1) + t172 * (t111 * pkin(4) + pkin(9) * t110 + t146) + t166 * (t111 * t140 + t116 * t137) + t165 * (t111 * t137 - t116 * t140) + t170 * t136 + t171 * t116 + t168 * t110) * g(3) + (-m(4) * t154 - m(5) * t149 - t106 * mrSges(4,1) - t98 * mrSges(5,1) - mrSges(1,2) + t172 * (t98 * pkin(4) + pkin(9) * t97 + t149) + t166 * (-t105 * t137 + t140 * t98) + t165 * (t105 * t140 + t137 * t98) + t163 * t132 - t171 * t105 + t168 * t97 + t173 * t135) * g(2) + (-m(4) * t152 - m(5) * t147 - t108 * mrSges(4,1) - t100 * mrSges(5,1) - mrSges(1,1) + t172 * (t100 * pkin(4) + pkin(9) * t99 + t147) + t166 * (t100 * t140 - t107 * t137) + t165 * (t100 * t137 + t107 * t140) + t163 * t135 - t171 * t107 + t168 * t99 - t173 * t132) * g(1);
U  = t1;
