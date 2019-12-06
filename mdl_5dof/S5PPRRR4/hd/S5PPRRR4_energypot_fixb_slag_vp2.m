% Calculate potential energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:31
% EndTime: 2019-12-05 15:18:32
% DurationCPUTime: 0.62s
% Computational Cost: add. (356->90), mult. (896->124), div. (0->0), fcn. (1108->14), ass. (0->48)
t125 = sin(pkin(6));
t129 = cos(pkin(6));
t130 = cos(pkin(5));
t126 = sin(pkin(5));
t127 = cos(pkin(11));
t154 = t126 * t127;
t108 = -t125 * t154 + t129 * t130;
t123 = sin(pkin(11));
t128 = cos(pkin(10));
t124 = sin(pkin(10));
t155 = t124 * t130;
t111 = -t123 * t128 - t127 * t155;
t152 = t126 * t129;
t102 = -t111 * t125 + t124 * t152;
t163 = -m(5) - m(6);
t162 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t161 = -t131 * mrSges(6,1) - t134 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t160 = -m(6) * pkin(4) - t134 * mrSges(6,1) + t131 * mrSges(6,2) - mrSges(5,1);
t159 = cos(qJ(3));
t157 = t123 * t126;
t156 = t124 * t126;
t153 = t126 * t128;
t151 = t128 * t130;
t149 = t130 * qJ(2) + qJ(1);
t148 = t128 * pkin(1) + qJ(2) * t156;
t145 = t125 * t159;
t144 = t129 * t159;
t143 = t126 * t145;
t142 = t124 * pkin(1) - qJ(2) * t153;
t109 = -t123 * t124 + t127 * t151;
t101 = -t109 * t125 - t128 * t152;
t112 = -t123 * t155 + t127 * t128;
t141 = t112 * pkin(2) + t102 * pkin(7) + t148;
t140 = pkin(2) * t157 + t108 * pkin(7) + t149;
t110 = t123 * t151 + t124 * t127;
t137 = t110 * pkin(2) + pkin(7) * t101 + t142;
t135 = cos(qJ(4));
t133 = sin(qJ(3));
t132 = sin(qJ(4));
t100 = t130 * t125 * t133 + (t127 * t129 * t133 + t123 * t159) * t126;
t99 = -t130 * t145 + t133 * t157 - t144 * t154;
t93 = t112 * t159 + (t111 * t129 + t125 * t156) * t133;
t92 = -t111 * t144 + t112 * t133 - t124 * t143;
t91 = t110 * t159 + (t109 * t129 - t125 * t153) * t133;
t90 = -t109 * t144 + t110 * t133 + t128 * t143;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t149 - t130 * mrSges(3,3) - (mrSges(3,1) * t123 + mrSges(3,2) * t127) * t126 - m(4) * t140 - t100 * mrSges(4,1) - t108 * mrSges(4,3) + t163 * (t100 * pkin(3) + pkin(8) * t99 + t140) + t160 * (t100 * t135 + t108 * t132) + t161 * t99 + t162 * (t100 * t132 - t108 * t135)) * g(3) + (-m(3) * t142 - m(4) * t137 - t124 * mrSges(2,1) - t110 * mrSges(3,1) - t91 * mrSges(4,1) - t128 * mrSges(2,2) - t109 * mrSges(3,2) + mrSges(3,3) * t153 - t101 * mrSges(4,3) - mrSges(1,2) + t163 * (t91 * pkin(3) + pkin(8) * t90 + t137) + t160 * (t101 * t132 + t135 * t91) + t161 * t90 + t162 * (-t101 * t135 + t132 * t91)) * g(2) + (-m(3) * t148 - m(4) * t141 - t128 * mrSges(2,1) - t112 * mrSges(3,1) - t93 * mrSges(4,1) + t124 * mrSges(2,2) - t111 * mrSges(3,2) - mrSges(3,3) * t156 - t102 * mrSges(4,3) - mrSges(1,1) + t163 * (t93 * pkin(3) + pkin(8) * t92 + t141) + t160 * (t102 * t132 + t135 * t93) + t161 * t92 + t162 * (-t102 * t135 + t132 * t93)) * g(1);
U = t1;
