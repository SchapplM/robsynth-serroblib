% Calculate potential energy for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:21
% EndTime: 2019-12-31 18:29:22
% DurationCPUTime: 0.28s
% Computational Cost: add. (151->49), mult. (161->40), div. (0->0), fcn. (138->10), ass. (0->20)
t53 = pkin(9) + qJ(5);
t48 = sin(t53);
t50 = cos(t53);
t55 = sin(pkin(9));
t57 = cos(pkin(9));
t81 = -mrSges(4,1) - m(5) * pkin(3) - t57 * mrSges(5,1) + t55 * mrSges(5,2) - m(6) * (pkin(4) * t57 + pkin(3)) - t50 * mrSges(6,1) + t48 * mrSges(6,2);
t80 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) + m(6) * (-pkin(7) - qJ(4)) - mrSges(6,3);
t79 = m(4) + m(5) + m(6);
t54 = pkin(8) + qJ(3);
t49 = sin(t54);
t51 = cos(t54);
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t76 = -m(3) * pkin(1) - t58 * mrSges(3,1) + t56 * mrSges(3,2) + t80 * t49 + t81 * t51 - mrSges(2,1);
t75 = m(3) * qJ(2) + t48 * mrSges(6,1) + t57 * mrSges(5,2) + t50 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t55;
t62 = cos(qJ(1));
t61 = sin(qJ(1));
t60 = -pkin(6) - qJ(2);
t46 = pkin(2) * t58 + pkin(1);
t1 = (-t56 * mrSges(3,1) - t58 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) - t79 * (t56 * pkin(2) + pkin(5)) - t80 * t51 + t81 * t49) * g(3) + (-mrSges(1,2) - t79 * (t61 * t46 + t62 * t60) + t75 * t62 + t76 * t61) * g(2) + (-mrSges(1,1) + (t79 * t60 - t75) * t61 + (-t79 * t46 + t76) * t62) * g(1);
U = t1;
