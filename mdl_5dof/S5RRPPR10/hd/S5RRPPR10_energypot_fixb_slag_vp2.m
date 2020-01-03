% Calculate potential energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:57
% EndTime: 2019-12-31 19:42:58
% DurationCPUTime: 0.46s
% Computational Cost: add. (131->60), mult. (249->62), div. (0->0), fcn. (250->8), ass. (0->29)
t100 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t99 = -t75 * mrSges(3,1) - mrSges(2,1) + (m(6) * pkin(7) - t100) * t72;
t98 = -m(5) - m(6);
t96 = mrSges(2,2) - mrSges(3,3);
t69 = sin(pkin(8));
t70 = cos(pkin(8));
t95 = (pkin(3) * t70 + qJ(4) * t69) * t72;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t92 = -t71 * mrSges(6,1) - t74 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t91 = -m(6) * pkin(4) - t74 * mrSges(6,1) + t71 * mrSges(6,2) - mrSges(4,1) - mrSges(5,1);
t90 = t72 * pkin(2) + pkin(5);
t73 = sin(qJ(1));
t87 = t73 * t75;
t76 = cos(qJ(1));
t86 = t75 * t76;
t85 = t76 * pkin(1) + t73 * pkin(6);
t84 = qJ(3) * t72;
t83 = t73 * pkin(1) - pkin(6) * t76;
t82 = pkin(2) * t86 + t76 * t84 + t85;
t81 = -qJ(3) * t75 + t90;
t79 = pkin(2) * t87 + t73 * t84 + t83;
t57 = t73 * t69 + t70 * t86;
t56 = t69 * t86 - t73 * t70;
t55 = -t69 * t76 + t70 * t87;
t54 = t69 * t87 + t70 * t76;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t81 - m(5) * (t81 + t95) - m(6) * (t90 + t95) + (-m(2) - m(3)) * pkin(5) + (-m(6) * (pkin(7) - qJ(3)) + t100) * t75 + (t92 * t69 + t91 * t70 - mrSges(3,1)) * t72) * g(3) + (-m(3) * t83 - m(4) * t79 - mrSges(1,2) + t98 * (t55 * pkin(3) + t54 * qJ(4) + t79) - t96 * t76 + t91 * t55 + t92 * t54 + t99 * t73) * g(2) + (-m(3) * t85 - m(4) * t82 - mrSges(1,1) + t98 * (t57 * pkin(3) + t56 * qJ(4) + t82) + t96 * t73 + t91 * t57 + t92 * t56 + t99 * t76) * g(1);
U = t1;
