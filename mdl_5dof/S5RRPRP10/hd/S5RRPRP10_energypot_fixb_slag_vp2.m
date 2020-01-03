% Calculate potential energy for
% S5RRPRP10
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:19
% EndTime: 2019-12-31 20:09:20
% DurationCPUTime: 0.39s
% Computational Cost: add. (106->60), mult. (179->56), div. (0->0), fcn. (156->6), ass. (0->27)
t64 = sin(qJ(4));
t95 = -m(6) * pkin(4) * t64 + mrSges(3,2) - mrSges(4,3);
t94 = m(6) * (-qJ(5) - pkin(7)) - mrSges(3,1) + mrSges(4,2) - mrSges(6,3);
t93 = -m(4) - m(6);
t92 = mrSges(5,1) + mrSges(6,1);
t91 = mrSges(5,2) + mrSges(6,2);
t66 = sin(qJ(1));
t65 = sin(qJ(2));
t76 = qJ(3) * t65;
t68 = cos(qJ(2));
t81 = t66 * t68;
t90 = pkin(2) * t81 + t66 * t76;
t89 = m(5) - t93;
t88 = -m(5) * pkin(7) - mrSges(5,3);
t87 = t95 * t65 + t94 * t68 - mrSges(2,1);
t67 = cos(qJ(4));
t86 = m(6) * (t67 * pkin(4) + pkin(3)) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t83 = t66 * t64;
t82 = t66 * t67;
t69 = cos(qJ(1));
t80 = t69 * t64;
t79 = t69 * t67;
t78 = t69 * t68;
t77 = t69 * pkin(1) + t66 * pkin(6);
t61 = t66 * pkin(1);
t75 = -t69 * pkin(6) + t61;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) - t89 * (t65 * pkin(2) + pkin(5)) + (t89 * qJ(3) + t92 * t64 + t91 * t67 - t95) * t68 + (t88 + t94) * t65) * g(3) + (-mrSges(1,2) - m(3) * t75 - m(5) * (pkin(7) * t81 + t61 + t90) - mrSges(5,3) * t81 + t93 * (t75 + t90) - t92 * (t65 * t83 - t79) - t91 * (t65 * t82 + t80) + (-m(5) * (-pkin(3) - pkin(6)) + t86) * t69 + t87 * t66) * g(2) + (-m(3) * t77 - mrSges(1,1) + t88 * t78 - t89 * (pkin(2) * t78 + t69 * t76 + t77) - t92 * (t65 * t80 + t82) - t91 * (t65 * t79 - t83) + (-m(5) * pkin(3) - t86) * t66 + t87 * t69) * g(1);
U = t1;
