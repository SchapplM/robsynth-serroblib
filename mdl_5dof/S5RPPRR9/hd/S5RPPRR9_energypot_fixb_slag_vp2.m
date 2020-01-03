% Calculate potential energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:15
% EndTime: 2019-12-31 18:02:16
% DurationCPUTime: 0.32s
% Computational Cost: add. (114->48), mult. (184->42), div. (0->0), fcn. (189->8), ass. (0->22)
t82 = m(5) + m(6);
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t81 = m(6) * pkin(4) + t58 * mrSges(6,1) - t56 * mrSges(6,2) + mrSges(5,1);
t80 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t78 = -mrSges(2,1) - mrSges(3,1);
t77 = mrSges(2,2) - mrSges(3,3);
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t76 = -t80 * t57 + t81 * t59 + mrSges(4,1);
t74 = t56 * mrSges(6,1) + t58 * mrSges(6,2) + t82 * pkin(6) - mrSges(4,2) + mrSges(5,3);
t73 = cos(qJ(1));
t72 = sin(qJ(1));
t71 = t73 * pkin(1) + t72 * qJ(2);
t70 = cos(pkin(8));
t69 = sin(pkin(8));
t68 = t73 * pkin(2) + t71;
t66 = t72 * pkin(1) - t73 * qJ(2);
t64 = t72 * pkin(2) + t66;
t46 = t73 * t69 - t72 * t70;
t45 = -t72 * t69 - t73 * t70;
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t81 * t57 + t80 * t59 + (-m(4) - t82) * (-qJ(3) + pkin(5)) + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t66 - m(4) * t64 - mrSges(1,2) - t77 * t73 + t78 * t72 - t82 * (-t46 * pkin(3) + t64) + t76 * t46 + t74 * t45) * g(2) + (-m(3) * t71 - m(4) * t68 - mrSges(1,1) + t78 * t73 + t77 * t72 - t82 * (-t45 * pkin(3) + t68) - t74 * t46 + t76 * t45) * g(1);
U = t1;
