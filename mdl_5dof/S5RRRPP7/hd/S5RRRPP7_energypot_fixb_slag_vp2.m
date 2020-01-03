% Calculate potential energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:53
% EndTime: 2019-12-31 21:03:53
% DurationCPUTime: 0.38s
% Computational Cost: add. (121->60), mult. (223->61), div. (0->0), fcn. (214->6), ass. (0->28)
t90 = -m(5) - m(6);
t89 = mrSges(2,2) - mrSges(3,3);
t64 = sin(qJ(3));
t65 = sin(qJ(2));
t67 = cos(qJ(3));
t88 = (pkin(3) * t67 + qJ(4) * t64) * t65;
t87 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t86 = -mrSges(4,3) - mrSges(5,2) + mrSges(6,3);
t68 = cos(qJ(2));
t85 = -t68 * mrSges(3,1) - mrSges(2,1) + (m(6) * qJ(5) + mrSges(3,2)) * t65;
t84 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1) - mrSges(6,1);
t83 = t65 * pkin(2) + pkin(5);
t66 = sin(qJ(1));
t82 = t65 * t66;
t69 = cos(qJ(1));
t81 = t65 * t69;
t80 = t66 * t68;
t79 = t68 * t69;
t78 = t69 * pkin(1) + t66 * pkin(6);
t76 = t66 * pkin(1) - pkin(6) * t69;
t75 = pkin(2) * t79 + pkin(7) * t81 + t78;
t74 = -pkin(7) * t68 + t83;
t72 = pkin(2) * t80 + pkin(7) * t82 + t76;
t52 = t66 * t64 + t67 * t79;
t51 = t64 * t79 - t66 * t67;
t50 = -t64 * t69 + t67 * t80;
t49 = t64 * t80 + t67 * t69;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t74 - m(5) * (t74 + t88) - m(6) * (t83 + t88) + (-m(2) - m(3)) * pkin(5) + (-mrSges(3,2) - m(6) * (-pkin(7) + qJ(5)) - t86) * t68 + (t87 * t64 + t84 * t67 - mrSges(3,1)) * t65) * g(3) + (-m(3) * t76 - m(4) * t72 - mrSges(1,2) + t90 * (t50 * pkin(3) + t49 * qJ(4) + t72) - t89 * t69 + t85 * t66 + t86 * t82 + t84 * t50 + t87 * t49) * g(2) + (-m(3) * t78 - m(4) * t75 - mrSges(1,1) + t90 * (t52 * pkin(3) + t51 * qJ(4) + t75) + t85 * t69 + t89 * t66 + t86 * t81 + t84 * t52 + t87 * t51) * g(1);
U = t1;
