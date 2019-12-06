% Calculate potential energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:10
% EndTime: 2019-12-05 17:23:10
% DurationCPUTime: 0.61s
% Computational Cost: add. (356->90), mult. (896->124), div. (0->0), fcn. (1108->14), ass. (0->48)
t123 = sin(pkin(6));
t126 = cos(pkin(6));
t127 = cos(pkin(5));
t124 = sin(pkin(5));
t134 = cos(qJ(2));
t152 = t124 * t134;
t107 = -t123 * t152 + t127 * t126;
t122 = sin(pkin(11));
t125 = cos(pkin(11));
t131 = sin(qJ(2));
t149 = t127 * t134;
t110 = -t122 * t149 - t125 * t131;
t154 = t124 * t126;
t101 = -t110 * t123 + t122 * t154;
t162 = -m(5) - m(6);
t161 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t128 = sin(qJ(5));
t132 = cos(qJ(5));
t160 = -t128 * mrSges(6,1) - t132 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t159 = -m(6) * pkin(4) - t132 * mrSges(6,1) + t128 * mrSges(6,2) - mrSges(5,1);
t158 = cos(qJ(3));
t156 = t122 * t124;
t155 = t124 * t125;
t153 = t124 * t131;
t150 = t127 * t131;
t148 = t125 * pkin(1) + pkin(7) * t156;
t147 = t127 * pkin(7) + qJ(1);
t144 = t123 * t158;
t143 = t126 * t158;
t142 = t124 * t144;
t141 = t122 * pkin(1) - pkin(7) * t155;
t108 = -t122 * t131 + t125 * t149;
t100 = -t108 * t123 - t125 * t154;
t111 = -t122 * t150 + t125 * t134;
t140 = t111 * pkin(2) + t101 * pkin(8) + t148;
t139 = pkin(2) * t153 + t107 * pkin(8) + t147;
t109 = t122 * t134 + t125 * t150;
t136 = t109 * pkin(2) + pkin(8) * t100 + t141;
t133 = cos(qJ(4));
t130 = sin(qJ(3));
t129 = sin(qJ(4));
t99 = t127 * t123 * t130 + (t126 * t130 * t134 + t131 * t158) * t124;
t98 = -t127 * t144 + t130 * t153 - t143 * t152;
t92 = t111 * t158 + (t110 * t126 + t123 * t156) * t130;
t91 = -t110 * t143 + t111 * t130 - t122 * t142;
t90 = t109 * t158 + (t108 * t126 - t123 * t155) * t130;
t89 = -t108 * t143 + t109 * t130 + t125 * t142;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t147 - t127 * mrSges(3,3) - (mrSges(3,1) * t131 + mrSges(3,2) * t134) * t124 - m(4) * t139 - t99 * mrSges(4,1) - t107 * mrSges(4,3) + t162 * (t99 * pkin(3) + t98 * pkin(9) + t139) + t159 * (t107 * t129 + t133 * t99) + t160 * t98 + t161 * (-t107 * t133 + t129 * t99)) * g(3) + (-m(3) * t141 - m(4) * t136 - t122 * mrSges(2,1) - t109 * mrSges(3,1) - t90 * mrSges(4,1) - t125 * mrSges(2,2) - t108 * mrSges(3,2) + mrSges(3,3) * t155 - t100 * mrSges(4,3) - mrSges(1,2) + t162 * (t90 * pkin(3) + pkin(9) * t89 + t136) + t159 * (t100 * t129 + t133 * t90) + t160 * t89 + t161 * (-t100 * t133 + t129 * t90)) * g(2) + (-m(3) * t148 - m(4) * t140 - t125 * mrSges(2,1) - t111 * mrSges(3,1) - t92 * mrSges(4,1) + t122 * mrSges(2,2) - t110 * mrSges(3,2) - mrSges(3,3) * t156 - t101 * mrSges(4,3) - mrSges(1,1) + t162 * (t92 * pkin(3) + pkin(9) * t91 + t140) + t159 * (t101 * t129 + t133 * t92) + t160 * t91 + t161 * (-t101 * t133 + t129 * t92)) * g(1);
U = t1;
