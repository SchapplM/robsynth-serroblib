% Calculate potential energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:03
% EndTime: 2019-03-09 06:26:04
% DurationCPUTime: 0.42s
% Computational Cost: add. (164->62), mult. (209->55), div. (0->0), fcn. (184->8), ass. (0->24)
t70 = cos(qJ(4));
t56 = t70 * pkin(4) + pkin(3);
t66 = qJ(4) + qJ(5);
t58 = cos(t66);
t67 = sin(qJ(4));
t95 = -m(6) * t56 - m(7) * (pkin(5) * t58 + t56) - mrSges(4,1) - m(5) * pkin(3) - t70 * mrSges(5,1) + t67 * mrSges(5,2);
t73 = -pkin(9) - pkin(8);
t94 = m(6) * t73 + m(7) * (-qJ(6) + t73) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t93 = -mrSges(6,1) - mrSges(7,1);
t92 = mrSges(6,2) + mrSges(7,2);
t91 = -m(4) - m(5) - m(6) - m(7);
t57 = sin(t66);
t86 = pkin(4) * t67;
t89 = -m(6) * t86 - m(7) * (pkin(5) * t57 + t86) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - t67 * mrSges(5,1) - t70 * mrSges(5,2);
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t88 = t95 * t68 - t94 * t71 + mrSges(2,2) - mrSges(3,3);
t69 = sin(qJ(1));
t85 = t68 * t69;
t72 = cos(qJ(1));
t84 = t68 * t72;
t80 = t72 * pkin(1) + t69 * qJ(2);
t61 = t69 * pkin(1);
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t91 * (pkin(2) + pkin(6)) + (t92 * t57 + t93 * t58 + t95) * t71 + t94 * t68) * g(3) + (-m(3) * t61 - mrSges(1,2) + t93 * (t57 * t69 - t58 * t84) - t92 * (t57 * t84 + t58 * t69) + t91 * (t69 * pkin(7) + t61) + t89 * t69 + ((m(3) - t91) * qJ(2) - t88) * t72) * g(2) + (-m(3) * t80 - mrSges(1,1) + t93 * (t57 * t72 + t58 * t85) - t92 * (-t57 * t85 + t58 * t72) + t91 * (t72 * pkin(7) + t80) + t89 * t72 + t88 * t69) * g(1);
U  = t1;
