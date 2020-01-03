% Calculate potential energy for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:20
% EndTime: 2019-12-31 21:32:21
% DurationCPUTime: 0.43s
% Computational Cost: add. (131->63), mult. (249->65), div. (0->0), fcn. (250->8), ass. (0->31)
t94 = -m(5) - m(6);
t68 = sin(qJ(2));
t72 = cos(qJ(2));
t93 = -t72 * mrSges(3,1) + t68 * mrSges(3,2) - mrSges(2,1);
t92 = mrSges(2,2) - mrSges(3,3);
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t91 = (pkin(3) * t71 + qJ(4) * t67) * t68;
t90 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t89 = m(6) * pkin(8) - t90;
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t88 = -t66 * mrSges(6,1) - t70 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t87 = -m(6) * pkin(4) - t70 * mrSges(6,1) + t66 * mrSges(6,2) - mrSges(4,1) - mrSges(5,1);
t86 = t68 * pkin(2) + pkin(5);
t69 = sin(qJ(1));
t85 = t69 * t68;
t84 = t69 * t72;
t73 = cos(qJ(1));
t83 = t73 * t68;
t82 = t73 * t72;
t81 = t73 * pkin(1) + t69 * pkin(6);
t80 = t69 * pkin(1) - pkin(6) * t73;
t79 = pkin(2) * t82 + pkin(7) * t83 + t81;
t78 = -pkin(7) * t72 + t86;
t76 = pkin(2) * t84 + pkin(7) * t85 + t80;
t54 = t69 * t67 + t71 * t82;
t53 = t67 * t82 - t69 * t71;
t52 = -t67 * t73 + t71 * t84;
t51 = t67 * t84 + t71 * t73;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t78 - m(5) * (t78 + t91) - m(6) * (t86 + t91) + (-m(2) - m(3)) * pkin(5) + (-mrSges(3,2) - m(6) * (-pkin(7) + pkin(8)) + t90) * t72 + (t88 * t67 + t87 * t71 - mrSges(3,1)) * t68) * g(3) + (-m(3) * t80 - m(4) * t76 - mrSges(1,2) + t94 * (t52 * pkin(3) + t51 * qJ(4) + t76) - t92 * t73 + t93 * t69 + t87 * t52 + t88 * t51 + t89 * t85) * g(2) + (-m(3) * t81 - m(4) * t79 - mrSges(1,1) + t94 * (t54 * pkin(3) + t53 * qJ(4) + t79) + t93 * t73 + t92 * t69 + t87 * t54 + t88 * t53 + t89 * t83) * g(1);
U = t1;
