% Calculate potential energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:17
% EndTime: 2019-12-31 20:12:18
% DurationCPUTime: 0.40s
% Computational Cost: add. (108->59), mult. (189->57), div. (0->0), fcn. (170->6), ass. (0->27)
t98 = -mrSges(3,1) + mrSges(4,2);
t97 = mrSges(3,2) - mrSges(4,3);
t96 = -m(5) - m(6);
t95 = -mrSges(5,3) - mrSges(6,2);
t71 = sin(qJ(1));
t70 = sin(qJ(2));
t82 = qJ(3) * t70;
t73 = cos(qJ(2));
t86 = t71 * t73;
t94 = pkin(2) * t86 + t71 * t82;
t93 = t97 * t70 + t98 * t73 - mrSges(2,1);
t92 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t91 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t90 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t89 = t70 * pkin(2) + pkin(5);
t69 = sin(qJ(4));
t88 = t69 * t71;
t72 = cos(qJ(4));
t87 = t71 * t72;
t74 = cos(qJ(1));
t85 = t74 * t70;
t84 = t74 * t73;
t83 = t74 * pkin(1) + t71 * pkin(6);
t67 = t71 * pkin(1);
t80 = -pkin(6) * t74 + t67;
t79 = pkin(2) * t84 + t74 * t82 + t83;
t1 = (-m(4) * t89 - mrSges(1,3) - mrSges(2,3) + t96 * (t70 * pkin(7) + t89) + (-m(2) - m(3)) * pkin(5) + (-t90 * t72 + t91 * t69 + (m(4) - t96) * qJ(3) - t97) * t73 + (t95 + t98) * t70) * g(3) + (-mrSges(1,2) - m(3) * t80 - m(4) * (t80 + t94) + t95 * t86 + t96 * (pkin(7) * t86 + t67 + (-pkin(3) - pkin(6)) * t74 + t94) - t91 * (t70 * t88 - t72 * t74) + t90 * (t69 * t74 + t70 * t87) - t92 * t74 + t93 * t71) * g(2) + (-m(3) * t83 - m(4) * t79 - mrSges(1,1) + t95 * t84 + t96 * (t71 * pkin(3) + pkin(7) * t84 + t79) - t91 * (t69 * t85 + t87) - t90 * (-t72 * t85 + t88) + t93 * t74 + t92 * t71) * g(1);
U = t1;
