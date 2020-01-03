% Calculate potential energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:20
% EndTime: 2019-12-31 19:46:20
% DurationCPUTime: 0.43s
% Computational Cost: add. (116->64), mult. (179->60), div. (0->0), fcn. (156->8), ass. (0->27)
t63 = pkin(8) + qJ(5);
t57 = sin(t63);
t58 = cos(t63);
t64 = sin(pkin(8));
t98 = -m(6) * pkin(4) * t64 - t57 * mrSges(6,1) - t58 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t97 = -mrSges(3,1) + mrSges(4,2) + m(6) * (-pkin(7) - qJ(4)) - mrSges(6,3) - m(5) * qJ(4);
t95 = -m(5) - m(6);
t68 = sin(qJ(1));
t67 = sin(qJ(2));
t80 = qJ(3) * t67;
t69 = cos(qJ(2));
t83 = t68 * t69;
t94 = pkin(2) * t83 + t68 * t80;
t93 = m(4) - t95;
t90 = t58 * mrSges(6,1) - t57 * mrSges(6,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t89 = t98 * t67 + t97 * t69 - mrSges(2,1);
t70 = cos(qJ(1));
t86 = t67 * t70;
t85 = t68 * t64;
t65 = cos(pkin(8));
t84 = t68 * t65;
t82 = t69 * t70;
t81 = t70 * pkin(1) + t68 * pkin(6);
t61 = t68 * pkin(1);
t77 = -pkin(6) * t70 + t61;
t56 = pkin(4) * t65 + pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) - t93 * (t67 * pkin(2) + pkin(5)) + (t64 * mrSges(5,1) + t65 * mrSges(5,2) + t93 * qJ(3) - t98) * t69 + (-mrSges(5,3) + t97) * t67) * g(3) + (-mrSges(1,2) - m(3) * t77 - m(4) * (t77 + t94) - mrSges(5,3) * t83 + t95 * (t61 + t94) + (-t85 * mrSges(5,1) - t84 * mrSges(5,2)) * t67 + (-m(5) * (-pkin(3) - pkin(6)) + t65 * mrSges(5,1) - t64 * mrSges(5,2) - m(6) * (-pkin(6) - t56) + t90) * t70 + t89 * t68) * g(2) + (-mrSges(1,1) - m(3) * t81 - (t64 * t86 + t84) * mrSges(5,1) - (t65 * t86 - t85) * mrSges(5,2) - mrSges(5,3) * t82 - t93 * (pkin(2) * t82 + t70 * t80 + t81) + t89 * t70 + (-m(5) * pkin(3) - m(6) * t56 - t90) * t68) * g(1);
U = t1;
