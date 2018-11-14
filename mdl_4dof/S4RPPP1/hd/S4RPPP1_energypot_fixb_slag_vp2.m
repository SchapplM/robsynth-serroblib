% Calculate potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:22
% EndTime: 2018-11-14 13:45:23
% DurationCPUTime: 0.27s
% Computational Cost: add. (192->53), mult. (213->53), div. (0->0), fcn. (177->10), ass. (0->32)
t98 = -m(4) - m(5);
t97 = mrSges(4,1) + mrSges(5,1) + mrSges(3,3);
t96 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t95 = -m(5) * pkin(3) - t97;
t94 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t74 = cos(pkin(4));
t93 = t74 * qJ(2) + pkin(5);
t75 = sin(qJ(1));
t76 = cos(qJ(1));
t72 = sin(pkin(4));
t89 = qJ(2) * t72;
t90 = t76 * pkin(1) + t75 * t89;
t88 = pkin(4) - pkin(6);
t87 = pkin(4) + pkin(6);
t86 = t76 * t89;
t85 = cos(t87);
t84 = sin(t88);
t71 = sin(pkin(6));
t82 = cos(t88) / 0.2e1;
t78 = t82 + t85 / 0.2e1;
t54 = t71 * t75 - t76 * t78;
t73 = cos(pkin(6));
t81 = sin(t87) / 0.2e1;
t77 = t81 - t84 / 0.2e1;
t55 = t75 * t73 + t76 * t77;
t69 = t75 * pkin(1);
t83 = t55 * pkin(2) + t54 * qJ(3) + t69;
t62 = t82 - t85 / 0.2e1;
t61 = t81 + t84 / 0.2e1;
t57 = t73 * t76 - t75 * t77;
t56 = t76 * t71 + t75 * t78;
t1 = (-m(2) * pkin(5) - m(3) * t93 - mrSges(1,3) - mrSges(2,3) + t98 * (t62 * pkin(2) - qJ(3) * t61 + t93) + t95 * t74 + t94 * t62 - t96 * t61) * g(3) + (-mrSges(1,2) - t75 * mrSges(2,1) - m(3) * (t69 - t86) - m(4) * (t83 - t86) - m(5) * t83 + t94 * t55 + t96 * t54 + (-mrSges(2,2) + (-m(5) * (-pkin(3) - qJ(2)) + t97) * t72) * t76) * g(2) + (-m(3) * t90 - t76 * mrSges(2,1) - mrSges(1,1) + t98 * (t57 * pkin(2) + t56 * qJ(3) + t90) + t94 * t57 + t96 * t56 + (t95 * t72 + mrSges(2,2)) * t75) * g(1);
U  = t1;
