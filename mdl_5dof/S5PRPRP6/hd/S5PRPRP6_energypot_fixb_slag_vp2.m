% Calculate potential energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:16
% EndTime: 2019-12-05 15:40:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (108->55), mult. (189->56), div. (0->0), fcn. (170->6), ass. (0->25)
t97 = mrSges(3,2) - mrSges(4,3);
t96 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t95 = -m(5) - m(6);
t68 = sin(pkin(7));
t71 = sin(qJ(2));
t82 = qJ(3) * t71;
t73 = cos(qJ(2));
t89 = t68 * t73;
t94 = pkin(2) * t89 + t68 * t82;
t93 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t92 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t91 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t90 = t97 * t71 + t96 * t73 - mrSges(2,1);
t69 = cos(pkin(7));
t88 = t69 * t73;
t70 = sin(qJ(4));
t87 = t70 * t71;
t72 = cos(qJ(4));
t86 = t71 * t72;
t83 = t69 * pkin(1) + t68 * pkin(5);
t81 = t71 * pkin(2) + qJ(1);
t64 = t68 * pkin(1);
t79 = -pkin(5) * t69 + t64;
t78 = pkin(2) * t88 + t69 * t82 + t83;
t1 = (-m(4) * t81 - mrSges(1,3) - mrSges(2,3) + t95 * (t71 * pkin(6) + t81) + (-m(2) - m(3)) * qJ(1) + (-t91 * t72 + t92 * t70 + (m(4) - t95) * qJ(3) - t97) * t73 + t96 * t71) * g(3) + (-mrSges(1,2) - m(3) * t79 - m(4) * (t79 + t94) + t95 * (pkin(6) * t89 + t64 + (-pkin(3) - pkin(5)) * t69 + t94) - t92 * (t68 * t87 - t69 * t72) + t91 * (t68 * t86 + t69 * t70) - t93 * t69 + t90 * t68) * g(2) + (-m(3) * t83 - m(4) * t78 - mrSges(1,1) + t95 * (t68 * pkin(3) + pkin(6) * t88 + t78) - t92 * (t68 * t72 + t69 * t87) - t91 * (t68 * t70 - t69 * t86) + t93 * t68 + t90 * t69) * g(1);
U = t1;
