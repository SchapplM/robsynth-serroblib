% Calculate potential energy for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:35
% EndTime: 2019-12-05 15:42:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (135->51), mult. (104->38), div. (0->0), fcn. (73->10), ass. (0->22)
t64 = -m(3) - m(4);
t65 = -mrSges(2,1) + pkin(1) * (-m(5) - m(6) + t64);
t49 = pkin(9) + qJ(4);
t43 = qJ(5) + t49;
t36 = sin(t43);
t37 = cos(t43);
t53 = cos(pkin(9));
t38 = t53 * pkin(3) + pkin(2);
t39 = sin(t49);
t41 = cos(t49);
t51 = sin(pkin(9));
t63 = -mrSges(3,1) - m(4) * pkin(2) - t53 * mrSges(4,1) + t51 * mrSges(4,2) - m(5) * t38 - t41 * mrSges(5,1) + t39 * mrSges(5,2) - m(6) * (pkin(4) * t41 + t38) - t37 * mrSges(6,1) + t36 * mrSges(6,2);
t55 = -pkin(6) - qJ(3);
t61 = m(4) * qJ(3) - m(5) * t55 - m(6) * (-pkin(7) + t55) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t56 = pkin(5) + qJ(1);
t60 = t51 * pkin(3) + t56;
t54 = cos(pkin(8));
t52 = sin(pkin(8));
t50 = pkin(8) + qJ(2);
t42 = cos(t50);
t40 = sin(t50);
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - mrSges(3,3) - t51 * mrSges(4,1) - t53 * mrSges(4,2) - m(5) * t60 - t39 * mrSges(5,1) - t41 * mrSges(5,2) - m(6) * (pkin(4) * t39 + t60) - t36 * mrSges(6,1) - t37 * mrSges(6,2) + t64 * t56) * g(3) + (-t54 * mrSges(2,2) + t63 * t40 + t61 * t42 + t52 * t65 - mrSges(1,2)) * g(2) + (t52 * mrSges(2,2) - t61 * t40 + t63 * t42 + t54 * t65 - mrSges(1,1)) * g(1);
U = t1;
