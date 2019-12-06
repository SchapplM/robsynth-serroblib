% Calculate potential energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:18:59
% EndTime: 2019-12-05 16:18:59
% DurationCPUTime: 0.22s
% Computational Cost: add. (135->51), mult. (104->38), div. (0->0), fcn. (73->10), ass. (0->22)
t64 = -m(3) - m(4);
t65 = -mrSges(2,1) + pkin(1) * (-m(5) - m(6) + t64);
t50 = qJ(3) + pkin(9);
t43 = qJ(5) + t50;
t36 = sin(t43);
t37 = cos(t43);
t56 = cos(qJ(3));
t38 = t56 * pkin(3) + pkin(2);
t40 = sin(t50);
t42 = cos(t50);
t55 = sin(qJ(3));
t63 = -mrSges(3,1) - m(4) * pkin(2) - t56 * mrSges(4,1) + t55 * mrSges(4,2) - m(5) * t38 - t42 * mrSges(5,1) + t40 * mrSges(5,2) - m(6) * (pkin(4) * t42 + t38) - t37 * mrSges(6,1) + t36 * mrSges(6,2);
t53 = -qJ(4) - pkin(6);
t61 = m(4) * pkin(6) - m(5) * t53 - m(6) * (-pkin(7) + t53) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t54 = pkin(5) + qJ(1);
t60 = t55 * pkin(3) + t54;
t52 = cos(pkin(8));
t51 = sin(pkin(8));
t49 = pkin(8) + qJ(2);
t41 = cos(t49);
t39 = sin(t49);
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - mrSges(3,3) - t55 * mrSges(4,1) - t56 * mrSges(4,2) - m(5) * t60 - t40 * mrSges(5,1) - t42 * mrSges(5,2) - m(6) * (pkin(4) * t40 + t60) - t36 * mrSges(6,1) - t37 * mrSges(6,2) + t64 * t54) * g(3) + (-t52 * mrSges(2,2) + t63 * t39 + t61 * t41 + t51 * t65 - mrSges(1,2)) * g(2) + (t51 * mrSges(2,2) - t61 * t39 + t63 * t41 + t52 * t65 - mrSges(1,1)) * g(1);
U = t1;
