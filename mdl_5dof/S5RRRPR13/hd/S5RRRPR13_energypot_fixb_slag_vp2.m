% Calculate potential energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:13
% EndTime: 2019-12-31 21:43:13
% DurationCPUTime: 0.44s
% Computational Cost: add. (196->76), mult. (420->84), div. (0->0), fcn. (481->10), ass. (0->45)
t121 = -m(6) * pkin(9) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t120 = -t89 * mrSges(6,1) - t93 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t119 = m(6) * (pkin(4) + pkin(8)) + t93 * mrSges(6,1) - t89 * mrSges(6,2) + mrSges(5,1) + mrSges(4,3);
t118 = mrSges(3,2) - t119;
t95 = cos(qJ(2));
t96 = cos(qJ(1));
t106 = t96 * t95;
t91 = sin(qJ(2));
t92 = sin(qJ(1));
t110 = t92 * t91;
t88 = cos(pkin(5));
t74 = -t88 * t106 + t110;
t116 = t74 * pkin(8);
t107 = t96 * t91;
t109 = t92 * t95;
t76 = t88 * t109 + t107;
t115 = t76 * pkin(8);
t114 = t88 * pkin(7) + pkin(6);
t87 = sin(pkin(5));
t113 = t87 * t91;
t112 = t87 * t95;
t111 = t92 * t87;
t108 = t96 * t87;
t105 = t96 * pkin(1) + pkin(7) * t111;
t104 = pkin(8) * t112;
t103 = pkin(2) * t113 + t114;
t77 = -t88 * t110 + t106;
t102 = t77 * pkin(2) + t105;
t101 = t92 * pkin(1) - pkin(7) * t108;
t75 = t88 * t107 + t109;
t100 = t75 * pkin(2) + t101;
t90 = sin(qJ(3));
t94 = cos(qJ(3));
t72 = t90 * t113 - t88 * t94;
t73 = t94 * t113 + t88 * t90;
t99 = t73 * pkin(3) + t72 * qJ(4) + t103;
t67 = -t94 * t111 + t77 * t90;
t68 = t90 * t111 + t77 * t94;
t98 = t68 * pkin(3) + t67 * qJ(4) + t102;
t65 = t94 * t108 + t75 * t90;
t66 = -t90 * t108 + t75 * t94;
t97 = t66 * pkin(3) + t65 * qJ(4) + t100;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t114 - t88 * mrSges(3,3) - (t91 * mrSges(3,1) + t95 * mrSges(3,2)) * t87 - m(4) * (t103 - t104) - m(5) * (t99 - t104) - m(6) * t99 + t120 * t72 + t119 * t112 + t121 * t73) * g(3) + (-mrSges(1,2) - t92 * mrSges(2,1) - t96 * mrSges(2,2) - m(3) * t101 - t75 * mrSges(3,1) + mrSges(3,3) * t108 - m(4) * (t100 + t116) - m(5) * (t97 + t116) - m(6) * t97 + t120 * t65 + t118 * t74 + t121 * t66) * g(2) + (-mrSges(1,1) - t96 * mrSges(2,1) + t92 * mrSges(2,2) - m(3) * t105 - t77 * mrSges(3,1) - mrSges(3,3) * t111 - m(4) * (t102 + t115) - m(5) * (t98 + t115) - m(6) * t98 + t120 * t67 + t118 * t76 + t121 * t68) * g(1);
U = t1;
