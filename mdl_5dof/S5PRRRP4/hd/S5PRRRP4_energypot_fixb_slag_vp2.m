% Calculate potential energy for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:32
% EndTime: 2019-12-05 16:45:33
% DurationCPUTime: 0.33s
% Computational Cost: add. (149->58), mult. (174->57), div. (0->0), fcn. (155->8), ass. (0->28)
t96 = -m(5) - m(6);
t95 = -mrSges(5,3) - mrSges(6,2);
t67 = qJ(2) + qJ(3);
t64 = sin(t67);
t65 = cos(t67);
t71 = sin(qJ(2));
t73 = cos(qJ(2));
t94 = -m(3) * pkin(1) - t73 * mrSges(3,1) - t65 * mrSges(4,1) + t71 * mrSges(3,2) + t64 * mrSges(4,2) - mrSges(2,1);
t93 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t92 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t91 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t90 = pkin(3) * t65;
t69 = cos(pkin(8));
t89 = t64 * t69;
t68 = sin(pkin(8));
t88 = t68 * t64;
t70 = sin(qJ(4));
t87 = t68 * t70;
t72 = cos(qJ(4));
t86 = t68 * t72;
t85 = t69 * t70;
t84 = t69 * t72;
t63 = pkin(2) * t73 + pkin(1);
t74 = -pkin(6) - pkin(5);
t83 = t68 * t63 + t69 * t74;
t82 = t71 * pkin(2) + qJ(1);
t80 = t69 * t63 - t68 * t74;
t1 = (-m(4) * t82 - t71 * mrSges(3,1) - t73 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t96 * (t64 * pkin(3) - pkin(7) * t65 + t82) + (-m(2) - m(3)) * qJ(1) + (-mrSges(4,2) - t95) * t65 + (t91 * t70 + t92 * t72 - mrSges(4,1)) * t64) * g(3) + (-m(4) * t83 - mrSges(1,2) + t95 * t88 + t96 * (pkin(7) * t88 + t68 * t90 + t83) + t92 * (t65 * t86 - t85) + t91 * (t65 * t87 + t84) - t93 * t69 + t94 * t68) * g(2) + (-m(4) * t80 - mrSges(1,1) + t95 * t89 + t96 * (pkin(7) * t89 + t69 * t90 + t80) + t92 * (t65 * t84 + t87) + t91 * (t65 * t85 - t86) + t94 * t69 + t93 * t68) * g(1);
U = t1;
