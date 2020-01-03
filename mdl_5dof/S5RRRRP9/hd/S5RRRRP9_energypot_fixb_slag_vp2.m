% Calculate potential energy for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:42
% EndTime: 2019-12-31 22:03:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (149->60), mult. (202->59), div. (0->0), fcn. (187->8), ass. (0->29)
t69 = sin(qJ(3));
t72 = cos(qJ(3));
t99 = -m(4) * pkin(2) - t72 * mrSges(4,1) + t69 * mrSges(4,2) - mrSges(3,1);
t98 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3);
t97 = m(3) + m(4);
t96 = -m(5) - m(6);
t95 = -mrSges(5,3) - mrSges(6,2);
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t94 = t98 * t70 + t99 * t73 - mrSges(2,1);
t93 = -t69 * mrSges(4,1) - t72 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t91 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t90 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t89 = pkin(3) * t69;
t71 = sin(qJ(1));
t88 = t71 * t70;
t87 = t71 * t73;
t68 = qJ(3) + qJ(4);
t63 = sin(t68);
t74 = cos(qJ(1));
t86 = t74 * t63;
t64 = cos(t68);
t85 = t74 * t64;
t84 = t74 * t70;
t83 = t74 * pkin(1) + t71 * pkin(6);
t75 = -pkin(8) - pkin(7);
t66 = t71 * pkin(1);
t61 = t72 * pkin(3) + pkin(2);
t1 = (-mrSges(1,3) - mrSges(2,3) + t96 * (t70 * t61 + t73 * t75 + pkin(5)) + (-m(2) - t97) * pkin(5) + (-t95 - t98) * t73 + (t90 * t63 + t91 * t64 + t99) * t70) * g(3) + (-mrSges(1,2) + t95 * t88 + t96 * (-t75 * t88 + t61 * t87 + t66 + (-pkin(6) - t89) * t74) - t97 * t66 + t91 * (t64 * t87 - t86) + t90 * (t63 * t87 + t85) + (t97 * pkin(6) - t93) * t74 + t94 * t71) * g(2) + (-mrSges(1,1) + t95 * t84 - t97 * t83 + t96 * (t74 * t73 * t61 + t71 * t89 - t75 * t84 + t83) + t91 * (t71 * t63 + t73 * t85) + t90 * (-t71 * t64 + t73 * t86) + t94 * t74 + t93 * t71) * g(1);
U = t1;
