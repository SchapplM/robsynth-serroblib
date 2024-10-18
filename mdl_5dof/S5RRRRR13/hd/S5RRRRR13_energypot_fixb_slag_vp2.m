% Calculate potential energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:17
% EndTime: 2024-09-27 17:30:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (195->59), mult. (135->52), div. (0->0), fcn. (106->16), ass. (0->28)
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t98 = -t79 * mrSges(5,2) + (-m(6) * pkin(4) - mrSges(5,1)) * t77;
t96 = m(6) * (pkin(9) + pkin(10)) + mrSges(5,3) + mrSges(6,3) + m(5) * pkin(9);
t93 = -m(4) - m(5) - m(6);
t73 = qJ(4) + qJ(5);
t92 = -m(5) * pkin(3) - m(6) * (pkin(4) * t79 + pkin(3)) - t79 * mrSges(5,1) - cos(t73) * mrSges(6,1) + t77 * mrSges(5,2) + sin(t73) * mrSges(6,2) - mrSges(4,1);
t64 = pkin(5) + t73;
t55 = sin(t64) / 0.2e1;
t65 = pkin(5) - t73;
t56 = cos(t65) / 0.2e1;
t57 = sin(t65);
t58 = cos(t64);
t75 = sin(pkin(5));
t76 = cos(pkin(5));
t91 = -(t55 - t57 / 0.2e1) * mrSges(6,1) - (t56 + t58 / 0.2e1) * mrSges(6,2) - mrSges(4,2) + t96 * t75 + t98 * t76;
t90 = pkin(7) + pkin(6);
t78 = sin(qJ(1));
t71 = t78 * pkin(1);
t80 = cos(qJ(1));
t72 = t80 * pkin(1);
t74 = qJ(1) + qJ(2);
t70 = qJ(3) + t74;
t69 = cos(t74);
t67 = sin(t74);
t62 = cos(t70);
t61 = sin(t70);
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t90 - mrSges(3,3) - mrSges(4,3) - (t56 - t58 / 0.2e1) * mrSges(6,1) - (t55 + t57 / 0.2e1) * mrSges(6,2) + t98 * t75 + t93 * (pkin(8) + t90) - t96 * t76) * g(3) + (-m(3) * t71 - t78 * mrSges(2,1) - t67 * mrSges(3,1) - t80 * mrSges(2,2) - t69 * mrSges(3,2) - mrSges(1,2) + t93 * (pkin(2) * t67 + t71) + t92 * t61 + t91 * t62) * g(2) + (-m(3) * t72 - t80 * mrSges(2,1) - t69 * mrSges(3,1) + t78 * mrSges(2,2) + t67 * mrSges(3,2) - mrSges(1,1) + t93 * (pkin(2) * t69 + t72) + t92 * t62 - t91 * t61) * g(1);
U = t1;
