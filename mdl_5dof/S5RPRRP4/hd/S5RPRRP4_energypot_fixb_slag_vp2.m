% Calculate potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:43
% EndTime: 2022-01-23 09:31:43
% DurationCPUTime: 0.46s
% Computational Cost: add. (142->62), mult. (180->57), div. (0->0), fcn. (161->8), ass. (0->29)
t85 = pkin(7) + pkin(6);
t118 = -m(4) * pkin(6) - m(5) * t85 + m(6) * (-qJ(5) - t85) + mrSges(3,2) - mrSges(6,3) - mrSges(4,3) - mrSges(5,3);
t117 = -m(4) - m(5);
t81 = sin(qJ(3));
t104 = t81 * pkin(3);
t111 = -m(3) - m(6);
t78 = qJ(3) + qJ(4);
t70 = sin(t78);
t116 = (m(4) - t111) * qJ(2) + m(5) * (qJ(2) + t104) + m(6) * (pkin(4) * t70 + t104) - mrSges(2,2) + mrSges(3,3);
t83 = cos(qJ(3));
t103 = t83 * pkin(3);
t71 = cos(t78);
t115 = -m(5) * t103 - m(6) * (pkin(4) * t71 + pkin(2) + t103) - mrSges(3,1);
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t112 = t111 * pkin(1) - mrSges(2,1) + t117 * (pkin(2) * t80 + pkin(1)) + t115 * t80 + t118 * t79;
t110 = -mrSges(5,1) - mrSges(6,1);
t109 = mrSges(5,2) + mrSges(6,2);
t82 = sin(qJ(1));
t98 = t82 * t70;
t97 = t82 * t71;
t96 = t82 * t81;
t95 = t82 * t83;
t84 = cos(qJ(1));
t94 = t84 * t70;
t93 = t84 * t71;
t92 = t84 * t81;
t91 = t84 * t83;
t1 = (-mrSges(1,3) - mrSges(2,3) + t117 * (t79 * pkin(2) + pkin(5)) + (-m(2) + t111) * pkin(5) - t118 * t80 + (-t83 * mrSges(4,1) + t81 * mrSges(4,2) + t109 * t70 + t110 * t71 + t115) * t79) * g(3) + (-mrSges(1,2) - (t80 * t95 - t92) * mrSges(4,1) - (-t80 * t96 - t91) * mrSges(4,2) + t110 * (t80 * t97 - t94) - t109 * (-t80 * t98 - t93) + t112 * t82 + t116 * t84) * g(2) + (-mrSges(1,1) - (t80 * t91 + t96) * mrSges(4,1) - (-t80 * t92 + t95) * mrSges(4,2) + t110 * (t80 * t93 + t98) - t109 * (-t80 * t94 + t97) + t112 * t84 - t116 * t82) * g(1);
U = t1;
