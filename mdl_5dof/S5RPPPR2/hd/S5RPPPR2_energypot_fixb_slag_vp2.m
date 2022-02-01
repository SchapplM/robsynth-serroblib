% Calculate potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:58:59
% EndTime: 2022-01-23 08:59:00
% DurationCPUTime: 0.59s
% Computational Cost: add. (158->88), mult. (270->102), div. (0->0), fcn. (280->10), ass. (0->37)
t109 = mrSges(4,2) - mrSges(5,3);
t81 = sin(pkin(8));
t114 = t109 * t81;
t86 = sin(qJ(5));
t103 = t81 * t86;
t80 = sin(pkin(9));
t83 = cos(pkin(9));
t76 = t83 * pkin(4) + t80 * pkin(6) + pkin(3);
t84 = cos(pkin(8));
t88 = cos(qJ(5));
t93 = qJ(4) * t84 - qJ(2);
t113 = (m(3) + m(4)) * qJ(2) - m(5) * (-t81 * pkin(3) + t93) - m(6) * (-t76 * t81 + t93) - (t103 * t83 + t88 * t84) * mrSges(6,2) - mrSges(2,2) + mrSges(3,3);
t111 = mrSges(6,1) * t88;
t85 = cos(pkin(7));
t100 = t85 * t83;
t82 = sin(pkin(7));
t67 = t100 * t84 + t82 * t80;
t110 = t67 * t88;
t108 = mrSges(5,2) - mrSges(6,3);
t102 = t81 * t88;
t105 = t82 * qJ(3) + pkin(1);
t104 = t81 * qJ(4) + pkin(2);
t64 = t76 * t84 + t104;
t75 = t84 * pkin(3) + t104;
t92 = t80 * pkin(4) - t83 * pkin(6) + qJ(3);
t106 = -m(4) * (pkin(2) * t85 + t105) - m(5) * (t75 * t85 + t105) - m(6) * (t64 * t85 + t82 * t92 + pkin(1)) - (t102 * t85 - t67 * t86) * mrSges(6,2) - mrSges(2,1) - m(3) * pkin(1) - t85 * mrSges(3,1) + t82 * mrSges(3,2);
t101 = t82 * t84;
t87 = sin(qJ(1));
t99 = t87 * t82;
t98 = t87 * t85;
t89 = cos(qJ(1));
t97 = t89 * t81;
t96 = t89 * t82;
t95 = t89 * t85;
t72 = t87 * t81 + t84 * t95;
t70 = t84 * t98 - t97;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(4) - m(5)) * (-t85 * qJ(3) + pkin(5)) + t108 * (t101 * t80 + t100) + (m(6) * t92 - mrSges(3,2) + mrSges(4,3)) * t85 + (-m(2) - m(3) - m(6)) * pkin(5) + (mrSges(6,2) * t86 - mrSges(5,1) - t111) * (t101 * t83 - t85 * t80) + (-m(4) * pkin(2) - m(5) * t75 - m(6) * t64 - t84 * mrSges(4,1) - t103 * mrSges(6,1) - t102 * mrSges(6,2) - mrSges(3,1) + t114) * t82) * g(3) + (-mrSges(1,2) - t70 * mrSges(4,1) - mrSges(4,3) * t99 - (t70 * t83 + t80 * t99) * mrSges(5,1) + t83 * t97 * t111 + (-mrSges(6,1) * t110 + t106) * t87 + t113 * t89 + (-mrSges(6,1) * t86 + t109) * (t81 * t98 + t89 * t84) + t108 * (t70 * t80 - t83 * t99)) * g(2) + (-mrSges(1,1) - t72 * mrSges(4,1) - mrSges(4,3) * t96 - (t72 * t83 + t80 * t96) * mrSges(5,1) + (-(t103 * t85 + t110) * mrSges(6,1) + t106) * t89 + t95 * t114 + t108 * (t72 * t80 - t83 * t96) + (-(t102 * t83 - t86 * t84) * mrSges(6,1) - t109 * t84 - t113) * t87) * g(1);
U = t1;
