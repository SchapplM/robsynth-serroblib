% Calculate potential energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:55
% EndTime: 2019-03-08 19:13:55
% DurationCPUTime: 0.76s
% Computational Cost: add. (383->90), mult. (733->108), div. (0->0), fcn. (880->14), ass. (0->51)
t136 = cos(qJ(2));
t166 = mrSges(3,2) * t136;
t165 = -m(4) - m(5);
t164 = -m(6) - m(7);
t125 = sin(pkin(11));
t129 = cos(pkin(11));
t134 = sin(qJ(2));
t163 = t134 * t125 - t129 * t136;
t162 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t128 = cos(pkin(12));
t161 = -mrSges(5,2) * t128 - mrSges(3,3) - mrSges(4,3);
t131 = cos(pkin(6));
t150 = t131 * t134;
t160 = t150 * mrSges(3,1) + t131 * t166 + mrSges(2,2);
t159 = m(3) * pkin(7) - t161;
t158 = -m(3) * pkin(1) - t136 * mrSges(3,1) + t134 * mrSges(3,2) - mrSges(2,1);
t124 = sin(pkin(12));
t157 = -m(5) * pkin(3) - t128 * mrSges(5,1) + t124 * mrSges(5,2) - mrSges(4,1);
t133 = sin(qJ(6));
t135 = cos(qJ(6));
t156 = -m(7) * pkin(5) - t135 * mrSges(7,1) + t133 * mrSges(7,2) - mrSges(6,1);
t155 = m(5) * qJ(4) + t133 * mrSges(7,1) + t135 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t154 = t124 * t131;
t126 = sin(pkin(10));
t127 = sin(pkin(6));
t153 = t126 * t127;
t130 = cos(pkin(10));
t152 = t127 * t130;
t106 = pkin(2) * t150 + (-pkin(7) - qJ(3)) * t127;
t118 = pkin(2) * t136 + pkin(1);
t147 = t130 * t106 + t126 * t118;
t146 = t131 * pkin(7) + qJ(1);
t145 = t124 * t153;
t144 = t124 * t152;
t143 = -t106 * t126 + t130 * t118;
t142 = t127 * t134 * pkin(2) + t131 * qJ(3) + t146;
t141 = t125 * t136 + t134 * t129;
t139 = t163 * t131;
t132 = -pkin(8) - qJ(4);
t123 = pkin(12) + qJ(5);
t120 = cos(t123);
t119 = sin(t123);
t117 = pkin(4) * t128 + pkin(3);
t105 = t141 * t131;
t104 = t141 * t127;
t103 = t163 * t127;
t96 = -t105 * t126 - t130 * t163;
t95 = t126 * t139 - t130 * t141;
t94 = t105 * t130 - t126 * t163;
t93 = -t126 * t141 - t130 * t139;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t146 - (mrSges(3,1) * t134 + t166) * t127 - t154 * mrSges(5,1) + t164 * (pkin(4) * t154 - t103 * t132 + t104 * t117 + t142) + t162 * (t104 * t119 - t131 * t120) + t165 * t142 + t157 * t104 + t161 * t131 + t156 * (t104 * t120 + t119 * t131) - t155 * t103) * g(3) + (t144 * mrSges(5,1) - mrSges(1,2) + t165 * t147 + t157 * t94 + t164 * (-pkin(4) * t144 + t94 * t117 + t93 * t132 + t147) + t162 * (t119 * t94 + t120 * t152) - t160 * t130 + t158 * t126 + t159 * t152 + t156 * (-t119 * t152 + t120 * t94) + t155 * t93) * g(2) + (-t145 * mrSges(5,1) - mrSges(1,1) + t165 * t143 + t157 * t96 + t164 * (pkin(4) * t145 + t96 * t117 + t95 * t132 + t143) + t162 * (t119 * t96 - t120 * t153) + t160 * t126 + t158 * t130 - t159 * t153 + t156 * (t119 * t153 + t120 * t96) + t155 * t95) * g(1);
U  = t1;
