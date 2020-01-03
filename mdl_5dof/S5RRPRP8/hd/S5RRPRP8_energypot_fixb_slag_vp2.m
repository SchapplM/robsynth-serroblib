% Calculate potential energy for
% S5RRPRP8
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:16
% EndTime: 2019-12-31 20:03:16
% DurationCPUTime: 0.42s
% Computational Cost: add. (108->61), mult. (186->59), div. (0->0), fcn. (167->6), ass. (0->31)
t101 = -mrSges(3,1) - mrSges(4,1);
t100 = mrSges(3,2) - mrSges(4,3);
t96 = -m(4) - m(5);
t99 = -m(6) + t96;
t69 = sin(qJ(4));
t70 = sin(qJ(2));
t72 = cos(qJ(4));
t73 = cos(qJ(2));
t58 = -t69 * t73 + t70 * t72;
t87 = t69 * t70;
t75 = t72 * t73 + t87;
t94 = -mrSges(5,2) - mrSges(6,2);
t95 = -mrSges(5,1) - mrSges(6,1);
t98 = t100 * t70 + t101 * t73 + t58 * t94 + t75 * t95 - mrSges(2,1);
t97 = m(5) * pkin(3);
t71 = sin(qJ(1));
t83 = qJ(3) * t70;
t86 = t71 * t73;
t93 = pkin(2) * t86 + t71 * t83;
t89 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3) - mrSges(6,3);
t88 = t70 * pkin(2) + pkin(5);
t74 = cos(qJ(1));
t85 = t73 * t74;
t84 = t74 * pkin(1) + t71 * pkin(6);
t82 = pkin(4) * t87;
t66 = t71 * pkin(1);
t81 = t66 + t93;
t80 = -pkin(6) * t74 + t66;
t68 = -qJ(5) - pkin(7);
t63 = pkin(4) * t72 + pkin(3);
t1 = (-m(6) * t88 - mrSges(1,3) - mrSges(2,3) + t96 * (-t73 * qJ(3) + t88) + (-m(6) * (-pkin(4) * t69 - qJ(3)) - t100) * t73 + (-m(6) * t63 + t101 - t97) * t70 + t95 * t58 - t94 * t75 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) - m(3) * t80 - m(4) * (t80 + t93) - m(5) * (pkin(3) * t86 + t81) - m(6) * (t63 * t86 + t81)) * g(2) + (-t85 * t97 - m(3) * t84 - mrSges(1,1) + t99 * (pkin(2) * t85 + t84)) * g(1) + ((-m(5) * (-pkin(6) + pkin(7)) - m(6) * (-pkin(6) - t68) + t89) * g(2) + (t99 * t83 - m(6) * (t63 * t73 + t82) + t98) * g(1)) * t74 + ((-m(6) * t82 + t98) * g(2) + (m(5) * pkin(7) - m(6) * t68 - t89) * g(1)) * t71;
U = t1;
