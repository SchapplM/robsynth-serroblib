% Calculate potential energy for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:42
% EndTime: 2019-12-05 15:21:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (157->49), mult. (148->41), div. (0->0), fcn. (125->10), ass. (0->20)
t52 = pkin(9) + qJ(5);
t46 = sin(t52);
t48 = cos(t52);
t54 = sin(pkin(9));
t57 = cos(pkin(9));
t78 = -mrSges(4,1) - m(5) * pkin(3) - t57 * mrSges(5,1) + t54 * mrSges(5,2) - m(6) * (pkin(4) * t57 + pkin(3)) - t48 * mrSges(6,1) + t46 * mrSges(6,2);
t77 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) + m(6) * (-pkin(6) - qJ(4)) - mrSges(6,3);
t76 = m(4) + m(5) + m(6);
t55 = sin(pkin(8));
t58 = cos(pkin(8));
t73 = t77 * t55 + t78 * t58 - mrSges(3,1);
t72 = t46 * mrSges(6,1) + t57 * mrSges(5,2) + t48 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t54;
t56 = sin(pkin(7));
t50 = t56 * pkin(1);
t59 = cos(pkin(7));
t51 = t59 * pkin(1);
t53 = pkin(7) + qJ(2);
t49 = cos(t53);
t47 = sin(t53);
t1 = (-m(2) * qJ(1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - t76) * (pkin(5) + qJ(1)) - t77 * t58 + t78 * t55) * g(3) + (-m(3) * t50 - t56 * mrSges(2,1) - t59 * mrSges(2,2) - mrSges(1,2) - t76 * (t47 * pkin(2) + t50) + (t76 * qJ(3) + t72) * t49 + t73 * t47) * g(2) + (-m(3) * t51 - t59 * mrSges(2,1) + t56 * mrSges(2,2) - mrSges(1,1) - t76 * (t49 * pkin(2) + t47 * qJ(3) + t51) + t73 * t49 - t72 * t47) * g(1);
U = t1;
