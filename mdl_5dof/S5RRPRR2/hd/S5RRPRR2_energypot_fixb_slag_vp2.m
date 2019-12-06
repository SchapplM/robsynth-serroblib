% Calculate potential energy for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:00
% EndTime: 2019-12-05 18:27:00
% DurationCPUTime: 0.21s
% Computational Cost: add. (138->52), mult. (117->41), div. (0->0), fcn. (86->10), ass. (0->22)
t56 = cos(qJ(2));
t44 = t56 * pkin(2) + pkin(1);
t52 = qJ(2) + pkin(9);
t46 = cos(t52);
t36 = pkin(3) * t46 + t44;
t47 = qJ(4) + t52;
t43 = qJ(5) + t47;
t37 = sin(t43);
t38 = cos(t43);
t41 = sin(t47);
t42 = cos(t47);
t45 = sin(t52);
t54 = sin(qJ(2));
t69 = -mrSges(2,1) - m(6) * (pkin(4) * t42 + t36) - t38 * mrSges(6,1) + t37 * mrSges(6,2) - m(5) * t36 - t42 * mrSges(5,1) + t41 * mrSges(5,2) - m(4) * t44 - t46 * mrSges(4,1) + t45 * mrSges(4,2) - m(3) * pkin(1) - t56 * mrSges(3,1) + t54 * mrSges(3,2);
t53 = -qJ(3) - pkin(6);
t51 = -pkin(7) + t53;
t68 = mrSges(2,2) + m(6) * (-pkin(8) + t51) - mrSges(6,3) + m(5) * t51 - mrSges(5,3) + m(4) * t53 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t67 = t54 * pkin(2) + pkin(5);
t66 = pkin(3) * t45 + t67;
t57 = cos(qJ(1));
t55 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) - t54 * mrSges(3,1) - t56 * mrSges(3,2) - m(4) * t67 - t45 * mrSges(4,1) - t46 * mrSges(4,2) - m(5) * t66 - t41 * mrSges(5,1) - t42 * mrSges(5,2) - m(6) * (pkin(4) * t41 + t66) - t37 * mrSges(6,1) - t38 * mrSges(6,2) + (-m(2) - m(3)) * pkin(5)) * g(3) + (t69 * t55 - t68 * t57 - mrSges(1,2)) * g(2) + (t68 * t55 + t69 * t57 - mrSges(1,1)) * g(1);
U = t1;
