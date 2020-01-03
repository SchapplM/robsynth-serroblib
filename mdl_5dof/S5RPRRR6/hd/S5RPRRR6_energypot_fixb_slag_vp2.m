% Calculate potential energy for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:53
% EndTime: 2019-12-31 19:00:53
% DurationCPUTime: 0.27s
% Computational Cost: add. (149->50), mult. (124->42), div. (0->0), fcn. (97->10), ass. (0->24)
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t76 = -m(6) * pkin(4) - t58 * mrSges(6,1) + t55 * mrSges(6,2) - mrSges(5,1);
t75 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t74 = -m(3) - m(4);
t73 = m(5) + m(6);
t53 = qJ(3) + qJ(4);
t47 = sin(t53);
t48 = cos(t53);
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t71 = -m(4) * pkin(2) - t59 * mrSges(4,1) + t56 * mrSges(4,2) + t75 * t47 + t76 * t48 - mrSges(3,1);
t70 = m(4) * pkin(6) + t55 * mrSges(6,1) + t58 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t57 = sin(qJ(1));
t50 = t57 * pkin(1);
t60 = cos(qJ(1));
t51 = t60 * pkin(1);
t54 = qJ(2) + pkin(5);
t61 = -pkin(7) - pkin(6);
t52 = qJ(1) + pkin(9);
t46 = cos(t52);
t45 = sin(t52);
t44 = t59 * pkin(3) + pkin(2);
t1 = (-m(2) * pkin(5) - t56 * mrSges(4,1) - t59 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t73 * (t56 * pkin(3) + t54) + t74 * t54 - t75 * t48 + t76 * t47) * g(3) + (-t57 * mrSges(2,1) - t60 * mrSges(2,2) - mrSges(1,2) - t73 * (t45 * t44 + t46 * t61 + t50) + t74 * t50 + t70 * t46 + t71 * t45) * g(2) + (-t60 * mrSges(2,1) + t57 * mrSges(2,2) - mrSges(1,1) - t73 * (t46 * t44 + t51) + t74 * t51 + t71 * t46 + (t73 * t61 - t70) * t45) * g(1);
U = t1;
