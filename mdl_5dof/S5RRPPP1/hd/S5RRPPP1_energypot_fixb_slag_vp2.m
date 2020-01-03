% Calculate potential energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:17
% EndTime: 2019-12-31 19:23:17
% DurationCPUTime: 0.43s
% Computational Cost: add. (168->70), mult. (374->79), div. (0->0), fcn. (406->8), ass. (0->41)
t132 = -m(5) - m(6);
t131 = mrSges(2,2) - mrSges(3,3);
t100 = cos(qJ(2));
t98 = sin(qJ(2));
t130 = -mrSges(3,1) * t100 + mrSges(3,2) * t98 - mrSges(2,1);
t97 = cos(pkin(5));
t115 = t100 * t97;
t94 = sin(pkin(8));
t96 = cos(pkin(8));
t76 = -t96 * t115 + t94 * t98;
t121 = t96 * t98;
t77 = t94 * t115 + t121;
t129 = t77 * pkin(3) + qJ(4) * t76;
t101 = cos(qJ(1));
t99 = sin(qJ(1));
t120 = t99 * t97;
t95 = sin(pkin(5));
t128 = t101 * t95 + t98 * t120;
t127 = mrSges(5,1) + mrSges(6,1) + mrSges(4,3);
t126 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t125 = -m(6) * pkin(4) - t127;
t124 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t123 = t98 * pkin(2) + pkin(6);
t122 = t95 * t99;
t119 = t101 * pkin(1) + t99 * pkin(7);
t118 = qJ(3) * t97;
t116 = t100 * t95;
t114 = t100 * t99;
t112 = t101 * t98;
t111 = t100 * t101;
t109 = t98 * t122;
t108 = t95 * t112;
t107 = pkin(2) * t111 + qJ(3) * t108 + t99 * t118 + t119;
t105 = -qJ(3) * t116 + t123;
t92 = t99 * pkin(1);
t103 = qJ(3) * t109 + pkin(2) * t114 + t92 + (-pkin(7) - t118) * t101;
t74 = t96 * t111 + (-t112 * t97 + t122) * t94;
t73 = -t96 * t122 + (t100 * t94 + t97 * t121) * t101;
t72 = t96 * t114 - t128 * t94;
t71 = t94 * t114 + t128 * t96;
t1 = (-mrSges(1,3) - mrSges(2,3) - t98 * mrSges(3,1) - t100 * mrSges(3,2) - m(4) * t105 - m(5) * (t105 + t129) - m(6) * (t123 + t129) + (-m(2) - m(3)) * pkin(6) + t124 * t77 + t126 * t76 + (-m(6) * (-pkin(4) - qJ(3)) + t127) * t116) * g(3) + (-m(3) * t92 - m(4) * t103 - mrSges(1,2) + t130 * t99 + t132 * (t72 * pkin(3) + t71 * qJ(4) + t103) + (m(3) * pkin(7) - t131) * t101 + t125 * (-t101 * t97 + t109) + t124 * t72 + t126 * t71) * g(2) + (-m(3) * t119 - m(4) * t107 - mrSges(1,1) + t131 * t99 + t132 * (t74 * pkin(3) + qJ(4) * t73 + t107) + t130 * t101 + t125 * (t108 + t120) + t124 * t74 + t126 * t73) * g(1);
U = t1;
