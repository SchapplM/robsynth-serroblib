% Calculate potential energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:35
% EndTime: 2019-03-08 19:17:36
% DurationCPUTime: 0.69s
% Computational Cost: add. (336->85), mult. (733->104), div. (0->0), fcn. (881->12), ass. (0->47)
t122 = cos(qJ(2));
t154 = t122 * mrSges(3,2);
t153 = -m(6) - m(7);
t152 = mrSges(4,2) - mrSges(5,3);
t111 = sin(pkin(11));
t114 = cos(pkin(11));
t119 = sin(qJ(2));
t151 = t119 * t111 - t114 * t122;
t150 = -mrSges(3,3) - mrSges(4,3) - mrSges(5,1);
t149 = m(7) * pkin(9) - mrSges(6,2) + mrSges(7,3);
t148 = m(3) * pkin(7) - t150;
t116 = cos(pkin(6));
t137 = t116 * t119;
t147 = t137 * mrSges(3,1) + t116 * t154 + mrSges(2,2);
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t146 = -mrSges(7,1) * t117 - mrSges(7,2) * t120 - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t145 = -m(3) * pkin(1) - t122 * mrSges(3,1) + t119 * mrSges(3,2) - mrSges(2,1);
t144 = -m(7) * pkin(5) - t120 * mrSges(7,1) + t117 * mrSges(7,2) - mrSges(6,1);
t107 = pkin(2) * t122 + pkin(1);
t112 = sin(pkin(10));
t115 = cos(pkin(10));
t113 = sin(pkin(6));
t99 = pkin(2) * t137 + (-pkin(7) - qJ(3)) * t113;
t143 = t112 * t107 + t115 * t99;
t142 = t112 * t113;
t141 = t113 * t115;
t118 = sin(qJ(5));
t140 = t113 * t118;
t121 = cos(qJ(5));
t139 = t113 * t121;
t134 = t116 * pkin(7) + qJ(1);
t133 = t115 * t107 - t112 * t99;
t132 = t113 * t119 * pkin(2) + t116 * qJ(3) + t134;
t131 = t111 * t122 + t119 * t114;
t128 = t151 * t116;
t86 = -t112 * t131 - t115 * t128;
t129 = t131 * t116;
t87 = -t112 * t151 + t115 * t129;
t130 = t87 * pkin(3) - qJ(4) * t86 + t143;
t88 = t112 * t128 - t115 * t131;
t89 = -t112 * t129 - t115 * t151;
t127 = t89 * pkin(3) - qJ(4) * t88 + t133;
t97 = t151 * t113;
t98 = t131 * t113;
t126 = t98 * pkin(3) + qJ(4) * t97 + t132;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t134 - (t119 * mrSges(3,1) + t154) * t113 - m(4) * t132 - m(5) * t126 + t152 * t97 + t153 * (t116 * pkin(4) + pkin(8) * t98 + t126) - t149 * (t116 * t118 - t97 * t121) + t144 * (t116 * t121 + t118 * t97) + t146 * t98 + t150 * t116) * g(3) + (-m(4) * t143 - m(5) * t130 - mrSges(1,2) - t152 * t86 + t153 * (-pkin(4) * t141 + pkin(8) * t87 + t130) + t149 * (t115 * t140 - t121 * t86) - t147 * t115 + t145 * t112 + t144 * (-t115 * t139 - t118 * t86) + t146 * t87 + t148 * t141) * g(2) + (-m(4) * t133 - m(5) * t127 - mrSges(1,1) - t152 * t88 + t153 * (pkin(4) * t142 + pkin(8) * t89 + t127) - t149 * (t112 * t140 + t88 * t121) + t147 * t112 + t145 * t115 + t144 * (t112 * t139 - t118 * t88) + t146 * t89 - t148 * t142) * g(1);
U  = t1;
