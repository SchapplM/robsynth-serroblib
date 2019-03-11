% Calculate potential energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:27
% EndTime: 2019-03-09 01:48:27
% DurationCPUTime: 0.33s
% Computational Cost: add. (128->62), mult. (172->45), div. (0->0), fcn. (143->8), ass. (0->24)
t53 = pkin(9) + qJ(6);
t45 = sin(t53);
t46 = cos(t53);
t54 = sin(pkin(9));
t55 = cos(pkin(9));
t83 = -mrSges(5,1) - m(6) * pkin(4) - t55 * mrSges(6,1) + t54 * mrSges(6,2) - m(7) * (pkin(5) * t55 + pkin(4)) - t46 * mrSges(7,1) + t45 * mrSges(7,2);
t82 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3);
t81 = -m(6) - m(7);
t80 = m(5) - t81;
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t77 = t83 * t57 - t82 * t59 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t76 = -t54 * mrSges(6,1) - t45 * mrSges(7,1) - t55 * mrSges(6,2) - t46 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t75 = pkin(2) + pkin(6);
t74 = pkin(5) * t54;
t58 = sin(qJ(1));
t60 = cos(qJ(1));
t73 = t60 * pkin(1) + t58 * qJ(2);
t50 = t58 * pkin(1);
t69 = -t60 * qJ(2) + t50;
t47 = t58 * qJ(3);
t68 = t47 + t69;
t51 = t60 * pkin(7);
t1 = (-m(4) * t75 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t80 * (pkin(3) + t75) + t83 * t59 + t82 * t57) * g(3) + (-mrSges(1,2) - m(3) * t69 - m(4) * t68 - m(5) * (t51 + t68) + t81 * (t47 + t50 + t51) + (m(6) * qJ(2) - m(7) * (-qJ(2) + t74) + t76) * t60 + t77 * t58) * g(2) + (-m(3) * t73 - mrSges(1,1) + (-m(4) - t80) * (t60 * qJ(3) + t73) + t77 * t60 + (m(7) * t74 + t80 * pkin(7) - t76) * t58) * g(1);
U  = t1;
