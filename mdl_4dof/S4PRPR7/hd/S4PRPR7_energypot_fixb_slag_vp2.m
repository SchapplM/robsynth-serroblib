% Calculate potential energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:34
% EndTime: 2019-12-31 16:25:34
% DurationCPUTime: 0.26s
% Computational Cost: add. (69->44), mult. (117->41), div. (0->0), fcn. (96->6), ass. (0->20)
t71 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t70 = -t48 * mrSges(5,1) - t50 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t67 = m(4) + m(5);
t46 = sin(pkin(6));
t49 = sin(qJ(2));
t57 = qJ(3) * t49;
t51 = cos(qJ(2));
t63 = t46 * t51;
t66 = pkin(2) * t63 + t46 * t57;
t65 = t50 * mrSges(5,1) - t48 * mrSges(5,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t64 = t70 * t49 + t71 * t51 - mrSges(2,1);
t47 = cos(pkin(6));
t62 = t47 * t51;
t58 = t47 * pkin(1) + t46 * pkin(4);
t43 = t46 * pkin(1);
t55 = -t47 * pkin(4) + t43;
t54 = pkin(2) * t62 + t47 * t57 + t58;
t1 = (-mrSges(1,3) - mrSges(2,3) - t67 * (t49 * pkin(2) + qJ(1)) + (-m(2) - m(3)) * qJ(1) + (t67 * qJ(3) - t70) * t51 + (-m(5) * pkin(5) + t71) * t49) * g(3) + (-mrSges(1,2) - m(3) * t55 - m(4) * (t55 + t66) - m(5) * (pkin(5) * t63 + t43 + t66) + (-m(5) * (-pkin(3) - pkin(4)) + t65) * t47 + t64 * t46) * g(2) + (-mrSges(1,1) - m(3) * t58 - m(4) * t54 - m(5) * (pkin(5) * t62 + t54) + (-m(5) * pkin(3) - t65) * t46 + t64 * t47) * g(1);
U = t1;
