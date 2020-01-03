% Calculate potential energy for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:49
% EndTime: 2019-12-31 19:28:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (137->55), mult. (151->47), div. (0->0), fcn. (126->8), ass. (0->24)
t87 = mrSges(4,2) - mrSges(5,3);
t58 = qJ(2) + pkin(8);
t55 = sin(t58);
t56 = cos(t58);
t86 = pkin(3) * t56 + qJ(4) * t55;
t85 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1);
t83 = -m(5) - m(6);
t65 = cos(qJ(1));
t82 = t86 * t65;
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t69 = t55 * t60 + t56 * t63;
t70 = t55 * t63 - t56 * t60;
t81 = -m(3) * pkin(1) - t64 * mrSges(3,1) - t69 * mrSges(6,1) + t61 * mrSges(3,2) - t70 * mrSges(6,2) + t55 * t87 + t56 * t85 - mrSges(2,1);
t80 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3);
t78 = t61 * pkin(2) + pkin(5);
t54 = pkin(2) * t64 + pkin(1);
t51 = t65 * t54;
t59 = -qJ(3) - pkin(6);
t62 = sin(qJ(1));
t74 = -t59 * t62 + t51;
t1 = (-m(4) * t78 - t61 * mrSges(3,1) - t70 * mrSges(6,1) - t64 * mrSges(3,2) + t69 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) + t83 * (pkin(3) * t55 - qJ(4) * t56 + t78) - t87 * t56 + t85 * t55 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + (-m(6) * pkin(7) - t80) * t65 + (-m(4) + t83) * (t62 * t54 + t65 * t59) + (t83 * t86 + t81) * t62) * g(2) + (-mrSges(1,1) - m(4) * t74 - m(5) * (t74 + t82) - m(6) * (t51 + t82) + t81 * t65 + (-m(6) * (-pkin(7) - t59) + t80) * t62) * g(1);
U = t1;
