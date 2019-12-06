% Calculate potential energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:31
% EndTime: 2019-12-05 15:49:32
% DurationCPUTime: 0.66s
% Computational Cost: add. (253->74), mult. (556->93), div. (0->0), fcn. (668->12), ass. (0->41)
t115 = cos(qJ(2));
t140 = t115 * mrSges(3,2);
t139 = -m(5) - m(6);
t138 = -mrSges(3,3) - mrSges(4,3);
t104 = sin(pkin(10));
t107 = cos(pkin(10));
t112 = sin(qJ(2));
t137 = t112 * t104 - t107 * t115;
t136 = m(3) * pkin(6) - t138;
t135 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t109 = cos(pkin(5));
t126 = t109 * t112;
t134 = mrSges(3,1) * t126 + t109 * t140 + mrSges(2,2);
t110 = sin(qJ(5));
t113 = cos(qJ(5));
t133 = t110 * mrSges(6,1) + t113 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t132 = -m(3) * pkin(1) - mrSges(3,1) * t115 + mrSges(3,2) * t112 - mrSges(2,1);
t131 = -m(6) * pkin(4) - mrSges(6,1) * t113 + mrSges(6,2) * t110 - mrSges(5,1);
t101 = pkin(2) * t115 + pkin(1);
t105 = sin(pkin(9));
t108 = cos(pkin(9));
t106 = sin(pkin(5));
t92 = pkin(2) * t126 + (-pkin(6) - qJ(3)) * t106;
t130 = t105 * t101 + t108 * t92;
t129 = t105 * t106;
t127 = t108 * t106;
t123 = t109 * pkin(6) + qJ(1);
t122 = t108 * t101 - t105 * t92;
t121 = t106 * t112 * pkin(2) + t109 * qJ(3) + t123;
t119 = t104 * t115 + t107 * t112;
t117 = t137 * t109;
t114 = cos(qJ(4));
t111 = sin(qJ(4));
t91 = t119 * t109;
t90 = t119 * t106;
t89 = t137 * t106;
t83 = -t105 * t91 - t108 * t137;
t82 = t105 * t117 - t108 * t119;
t81 = -t105 * t137 + t108 * t91;
t80 = -t105 * t119 - t108 * t117;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t123 - (mrSges(3,1) * t112 + t140) * t106 - m(4) * t121 - t90 * mrSges(4,1) + t139 * (pkin(3) * t90 + pkin(7) * t89 + t121) + t131 * (t109 * t111 + t114 * t90) - t133 * t89 + t135 * (-t109 * t114 + t111 * t90) + t138 * t109) * g(3) + (-m(4) * t130 - t81 * mrSges(4,1) - mrSges(1,2) + t139 * (pkin(3) * t81 - pkin(7) * t80 + t130) + t131 * (-t111 * t127 + t114 * t81) + t133 * t80 + t135 * (t111 * t81 + t114 * t127) - t134 * t108 + t132 * t105 + t136 * t127) * g(2) + (-m(4) * t122 - t83 * mrSges(4,1) - mrSges(1,1) + t139 * (pkin(3) * t83 - pkin(7) * t82 + t122) + t131 * (t111 * t129 + t114 * t83) + t133 * t82 + t135 * (t111 * t83 - t114 * t129) + t134 * t105 + t132 * t108 - t136 * t129) * g(1);
U = t1;
