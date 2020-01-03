% Calculate potential energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:30
% EndTime: 2019-12-31 22:46:30
% DurationCPUTime: 0.59s
% Computational Cost: add. (356->90), mult. (896->121), div. (0->0), fcn. (1108->14), ass. (0->49)
t126 = cos(pkin(5));
t131 = sin(qJ(1));
t134 = cos(qJ(2));
t150 = t131 * t134;
t130 = sin(qJ(2));
t135 = cos(qJ(1));
t152 = t130 * t135;
t111 = -t126 * t150 - t152;
t123 = sin(pkin(6));
t125 = cos(pkin(6));
t124 = sin(pkin(5));
t156 = t124 * t131;
t102 = -t111 * t123 + t125 * t156;
t155 = t124 * t134;
t108 = -t123 * t155 + t125 * t126;
t164 = -m(5) - m(6);
t163 = -m(6) * pkin(11) + mrSges(5,2) - mrSges(6,3);
t127 = sin(qJ(5));
t132 = cos(qJ(5));
t162 = -t127 * mrSges(6,1) - t132 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t161 = -m(6) * pkin(4) - mrSges(6,1) * t132 + mrSges(6,2) * t127 - mrSges(5,1);
t160 = cos(qJ(3));
t159 = t126 * pkin(8) + pkin(7);
t157 = t124 * t130;
t154 = t124 * t135;
t151 = t131 * t130;
t149 = t134 * t135;
t148 = t135 * pkin(1) + pkin(8) * t156;
t145 = t123 * t160;
t144 = t125 * t160;
t143 = t124 * t145;
t142 = t131 * pkin(1) - pkin(8) * t154;
t109 = t126 * t149 - t151;
t101 = -t109 * t123 - t125 * t154;
t112 = -t126 * t151 + t149;
t141 = t112 * pkin(2) + t102 * pkin(9) + t148;
t140 = pkin(2) * t157 + t108 * pkin(9) + t159;
t110 = t126 * t152 + t150;
t137 = t110 * pkin(2) + pkin(9) * t101 + t142;
t133 = cos(qJ(4));
t129 = sin(qJ(3));
t128 = sin(qJ(4));
t100 = t126 * t123 * t129 + (t125 * t129 * t134 + t130 * t160) * t124;
t99 = -t126 * t145 + t129 * t157 - t144 * t155;
t95 = t112 * t160 + (t111 * t125 + t123 * t156) * t129;
t94 = -t111 * t144 + t112 * t129 - t131 * t143;
t93 = t110 * t160 + (t109 * t125 - t123 * t154) * t129;
t92 = -t109 * t144 + t110 * t129 + t135 * t143;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t159 - t126 * mrSges(3,3) - (t130 * mrSges(3,1) + t134 * mrSges(3,2)) * t124 - m(4) * t140 - t100 * mrSges(4,1) - t108 * mrSges(4,3) + t164 * (t100 * pkin(3) + pkin(10) * t99 + t140) + t161 * (t100 * t133 + t108 * t128) + t162 * t99 + t163 * (t100 * t128 - t108 * t133)) * g(3) + (-m(3) * t142 - m(4) * t137 - t131 * mrSges(2,1) - t110 * mrSges(3,1) - t93 * mrSges(4,1) - t135 * mrSges(2,2) - t109 * mrSges(3,2) + mrSges(3,3) * t154 - t101 * mrSges(4,3) - mrSges(1,2) + t164 * (t93 * pkin(3) + t92 * pkin(10) + t137) + t161 * (t101 * t128 + t133 * t93) + t162 * t92 + t163 * (-t101 * t133 + t128 * t93)) * g(2) + (-m(3) * t148 - m(4) * t141 - t135 * mrSges(2,1) - t112 * mrSges(3,1) - t95 * mrSges(4,1) + t131 * mrSges(2,2) - t111 * mrSges(3,2) - mrSges(3,3) * t156 - t102 * mrSges(4,3) - mrSges(1,1) + t164 * (t95 * pkin(3) + pkin(10) * t94 + t141) + t161 * (t102 * t128 + t133 * t95) + t162 * t94 + t163 * (-t102 * t133 + t128 * t95)) * g(1);
U = t1;
