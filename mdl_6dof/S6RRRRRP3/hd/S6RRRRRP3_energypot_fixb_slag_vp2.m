% Calculate potential energy for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:45
% EndTime: 2019-03-10 01:05:45
% DurationCPUTime: 0.41s
% Computational Cost: add. (226->67), mult. (221->61), div. (0->0), fcn. (196->10), ass. (0->29)
t78 = cos(qJ(4));
t63 = t78 * pkin(4) + pkin(3);
t73 = qJ(4) + qJ(5);
t68 = cos(t73);
t75 = sin(qJ(4));
t105 = -m(6) * t63 - m(7) * (pkin(5) * t68 + t63) - mrSges(4,1) - m(5) * pkin(3) - t78 * mrSges(5,1) + t75 * mrSges(5,2);
t81 = -pkin(10) - pkin(9);
t104 = m(6) * t81 + m(7) * (-qJ(6) + t81) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(9) - mrSges(5,3);
t103 = -mrSges(6,1) - mrSges(7,1);
t102 = mrSges(6,2) + mrSges(7,2);
t101 = -m(4) - m(6) - m(7);
t100 = -m(5) + t101;
t74 = qJ(2) + qJ(3);
t67 = sin(t74);
t69 = cos(t74);
t76 = sin(qJ(2));
t79 = cos(qJ(2));
t98 = -m(3) * pkin(1) - t79 * mrSges(3,1) + t76 * mrSges(3,2) + t104 * t67 + t105 * t69 - mrSges(2,1);
t66 = sin(t73);
t96 = pkin(4) * t75;
t97 = m(6) * t96 + m(7) * (pkin(5) * t66 + t96) - mrSges(2,2) + mrSges(4,3) + t75 * mrSges(5,1) + t78 * mrSges(5,2) + m(3) * pkin(7) + mrSges(3,3);
t77 = sin(qJ(1));
t94 = t69 * t77;
t80 = cos(qJ(1));
t93 = t69 * t80;
t82 = -pkin(8) - pkin(7);
t64 = pkin(2) * t79 + pkin(1);
t62 = t80 * t64;
t1 = (-mrSges(3,1) * t76 - mrSges(3,2) * t79 - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t100 * (t76 * pkin(2) + pkin(6)) - t104 * t69 + (t102 * t66 + t103 * t68 + t105) * t67) * g(3) + (-mrSges(1,2) + t103 * (-t66 * t80 + t68 * t94) - t102 * (-t66 * t94 - t68 * t80) + t100 * (t77 * t64 + t80 * t82) + t97 * t80 + t98 * t77) * g(2) + (-m(5) * t62 - mrSges(1,1) + t103 * (t66 * t77 + t68 * t93) - t102 * (-t66 * t93 + t68 * t77) + t101 * (-t77 * t82 + t62) + (m(5) * t82 - t97) * t77 + t98 * t80) * g(1);
U  = t1;
