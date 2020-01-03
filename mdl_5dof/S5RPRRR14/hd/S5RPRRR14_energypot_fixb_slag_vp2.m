% Calculate potential energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:48
% EndTime: 2019-12-31 19:16:48
% DurationCPUTime: 0.60s
% Computational Cost: add. (356->90), mult. (896->122), div. (0->0), fcn. (1108->14), ass. (0->48)
t123 = sin(pkin(11));
t128 = cos(pkin(5));
t135 = cos(qJ(1));
t126 = cos(pkin(11));
t132 = sin(qJ(1));
t149 = t132 * t126;
t111 = -t123 * t135 - t128 * t149;
t124 = sin(pkin(6));
t127 = cos(pkin(6));
t125 = sin(pkin(5));
t154 = t125 * t132;
t102 = -t111 * t124 + t127 * t154;
t155 = t125 * t126;
t108 = -t124 * t155 + t127 * t128;
t163 = -m(5) - m(6);
t162 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t161 = -t129 * mrSges(6,1) - t133 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t160 = -m(6) * pkin(4) - mrSges(6,1) * t133 + mrSges(6,2) * t129 - mrSges(5,1);
t159 = cos(qJ(3));
t158 = t128 * qJ(2) + pkin(7);
t156 = t123 * t125;
t153 = t125 * t135;
t151 = t128 * t135;
t150 = t132 * t123;
t148 = t135 * pkin(1) + qJ(2) * t154;
t145 = t124 * t159;
t144 = t127 * t159;
t143 = t125 * t145;
t142 = t132 * pkin(1) - qJ(2) * t153;
t109 = t126 * t151 - t150;
t101 = -t109 * t124 - t127 * t153;
t112 = t126 * t135 - t128 * t150;
t141 = t112 * pkin(2) + t102 * pkin(8) + t148;
t140 = pkin(2) * t156 + t108 * pkin(8) + t158;
t110 = t123 * t151 + t149;
t137 = t110 * pkin(2) + pkin(8) * t101 + t142;
t134 = cos(qJ(4));
t131 = sin(qJ(3));
t130 = sin(qJ(4));
t100 = t128 * t124 * t131 + (t126 * t127 * t131 + t123 * t159) * t125;
t99 = -t128 * t145 + t131 * t156 - t144 * t155;
t95 = t112 * t159 + (t111 * t127 + t124 * t154) * t131;
t94 = -t111 * t144 + t112 * t131 - t132 * t143;
t93 = t110 * t159 + (t109 * t127 - t124 * t153) * t131;
t92 = -t109 * t144 + t110 * t131 + t135 * t143;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t158 - t128 * mrSges(3,3) - (t123 * mrSges(3,1) + t126 * mrSges(3,2)) * t125 - m(4) * t140 - t100 * mrSges(4,1) - t108 * mrSges(4,3) + t163 * (t100 * pkin(3) + pkin(9) * t99 + t140) + t160 * (t100 * t134 + t108 * t130) + t161 * t99 + t162 * (t100 * t130 - t108 * t134)) * g(3) + (-m(3) * t142 - m(4) * t137 - t132 * mrSges(2,1) - t110 * mrSges(3,1) - t93 * mrSges(4,1) - t135 * mrSges(2,2) - t109 * mrSges(3,2) + mrSges(3,3) * t153 - t101 * mrSges(4,3) - mrSges(1,2) + t163 * (t93 * pkin(3) + t92 * pkin(9) + t137) + t160 * (t101 * t130 + t134 * t93) + t161 * t92 + t162 * (-t101 * t134 + t130 * t93)) * g(2) + (-m(3) * t148 - m(4) * t141 - t135 * mrSges(2,1) - t112 * mrSges(3,1) - t95 * mrSges(4,1) + t132 * mrSges(2,2) - t111 * mrSges(3,2) - mrSges(3,3) * t154 - t102 * mrSges(4,3) - mrSges(1,1) + t163 * (t95 * pkin(3) + pkin(9) * t94 + t141) + t160 * (t102 * t130 + t134 * t95) + t161 * t94 + t162 * (-t102 * t134 + t130 * t95)) * g(1);
U = t1;
