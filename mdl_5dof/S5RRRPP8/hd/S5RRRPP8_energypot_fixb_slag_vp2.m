% Calculate potential energy for
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:26
% EndTime: 2019-12-31 21:07:26
% DurationCPUTime: 0.35s
% Computational Cost: add. (121->60), mult. (223->61), div. (0->0), fcn. (214->6), ass. (0->29)
t92 = -m(5) - m(6);
t67 = sin(qJ(2));
t70 = cos(qJ(2));
t91 = -t70 * mrSges(3,1) + t67 * mrSges(3,2) - mrSges(2,1);
t90 = mrSges(2,2) - mrSges(3,3);
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t89 = (pkin(3) * t69 + qJ(4) * t66) * t67;
t88 = mrSges(5,1) + mrSges(6,1) + mrSges(4,3);
t87 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t86 = -m(6) * pkin(4) - t88;
t85 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t84 = t67 * pkin(2) + pkin(5);
t68 = sin(qJ(1));
t83 = t68 * t67;
t82 = t68 * t70;
t71 = cos(qJ(1));
t81 = t71 * t67;
t80 = t71 * t70;
t79 = t71 * pkin(1) + t68 * pkin(6);
t78 = t68 * pkin(1) - t71 * pkin(6);
t77 = pkin(2) * t80 + pkin(7) * t81 + t79;
t76 = -t70 * pkin(7) + t84;
t74 = pkin(2) * t82 + pkin(7) * t83 + t78;
t53 = t68 * t66 + t69 * t80;
t52 = t66 * t80 - t68 * t69;
t51 = -t71 * t66 + t69 * t82;
t50 = t66 * t82 + t71 * t69;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t76 - m(5) * (t76 + t89) - m(6) * (t84 + t89) + (-m(2) - m(3)) * pkin(5) + (-mrSges(3,2) - m(6) * (-pkin(4) - pkin(7)) + t88) * t70 + (t87 * t66 + t85 * t69 - mrSges(3,1)) * t67) * g(3) + (-m(3) * t78 - m(4) * t74 - mrSges(1,2) + t92 * (t51 * pkin(3) + t50 * qJ(4) + t74) - t90 * t71 + t91 * t68 + t86 * t83 + t85 * t51 + t87 * t50) * g(2) + (-m(3) * t79 - m(4) * t77 - mrSges(1,1) + t92 * (t53 * pkin(3) + t52 * qJ(4) + t77) + t91 * t71 + t90 * t68 + t86 * t81 + t85 * t53 + t87 * t52) * g(1);
U = t1;
